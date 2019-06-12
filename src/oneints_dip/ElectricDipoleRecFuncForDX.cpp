//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricDipoleRecFuncForDX.hpp"

namespace ediprecfunc { // ediprecfunc namespace

    void
    compElectricDipoleForDD(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForDD_0_36(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  nOSFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForDD_36_72(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   nOSFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        ediprecfunc::compElectricDipoleForDD_72_108(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    nOSFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 
    }

    void
    compElectricDipoleForDD_0_36(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const int32_t              nOSFactors,
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

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx); 

            auto tdy_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx); 

            auto tdz_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx); 

            auto tdx_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 1); 

            auto tdy_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 1); 

            auto tdz_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 1); 

            auto tdx_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 2); 

            auto tdy_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 2); 

            auto tdz_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 2); 

            auto tdx_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 3); 

            auto tdy_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 3); 

            auto tdz_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 3); 

            auto tdx_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 4); 

            auto tdy_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 4); 

            auto tdz_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 4); 

            auto tdx_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 5); 

            auto tdy_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 5); 

            auto tdz_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 5); 

            auto tdx_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 6); 

            auto tdy_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tdz_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tdx_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 7); 

            auto tdy_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tdz_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tdx_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 8); 

            auto tdy_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tdz_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tdx_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 9); 

            auto tdy_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tdz_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdx_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 10); 

            auto tdy_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tdz_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tdx_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 11); 

            auto tdy_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tdz_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tdx_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx); 

            auto tdy_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx); 

            auto tdz_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx); 

            auto tdx_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 1); 

            auto tdy_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 1); 

            auto tdz_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 1); 

            auto tdx_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 2); 

            auto tdy_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 2); 

            auto tdz_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 2); 

            auto tdx_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 3); 

            auto tdy_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 3); 

            auto tdz_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 3); 

            auto tdx_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 4); 

            auto tdy_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 4); 

            auto tdz_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 4); 

            auto tdx_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 5); 

            auto tdy_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 5); 

            auto tdz_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 5); 

            auto tdx_x_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx); 

            auto tdy_x_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx); 

            auto tdz_x_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx); 

            auto tdx_x_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 1); 

            auto tdy_x_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 1); 

            auto tdz_x_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 1); 

            auto tdx_x_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 2); 

            auto tdy_x_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 2); 

            auto tdz_x_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 2); 

            auto tdx_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 3); 

            auto tdy_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 3); 

            auto tdz_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 3); 

            auto tdx_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 4); 

            auto tdy_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 4); 

            auto tdz_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 4); 

            auto tdx_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 5); 

            auto tdy_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 5); 

            auto tdz_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 5); 

            auto ts_x_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx); 

            auto ts_x_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 1); 

            auto ts_x_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 2); 

            auto ts_x_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 3); 

            auto ts_x_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 4); 

            auto ts_x_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 5); 

            auto ts_y_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 6); 

            auto ts_y_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 7); 

            auto ts_y_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 8); 

            auto ts_y_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 9); 

            auto ts_y_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 10); 

            auto ts_y_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 11); 

            // set up pointers to integrals

            auto tdx_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx); 

            auto tdy_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx); 

            auto tdz_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx); 

            auto tdx_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 1); 

            auto tdy_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 1); 

            auto tdz_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 1); 

            auto tdx_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 2); 

            auto tdy_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 2); 

            auto tdz_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 2); 

            auto tdx_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 3); 

            auto tdy_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 3); 

            auto tdz_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 3); 

            auto tdx_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 4); 

            auto tdy_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 4); 

            auto tdz_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 4); 

            auto tdx_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 5); 

            auto tdy_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 5); 

            auto tdz_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 5); 

            auto tdx_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 6); 

            auto tdy_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 6); 

            auto tdz_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 6); 

            auto tdx_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 7); 

            auto tdy_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 7); 

            auto tdz_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 7); 

            auto tdx_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 8); 

            auto tdy_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 8); 

            auto tdz_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 8); 

            auto tdx_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 9); 

            auto tdy_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 9); 

            auto tdz_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 9); 

            auto tdx_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 10); 

            auto tdy_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 10); 

            auto tdz_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 10); 

            auto tdx_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 11); 

            auto tdy_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 11); 

            auto tdz_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 11); 

            // Batch of Integrals (0,36)

            #pragma omp simd aligned(fx, pa_x, tdx_0_xx_0, tdx_0_xy_0, tdx_0_xz_0, tdx_0_yy_0, tdx_0_yz_0, \
                                     tdx_0_zz_0, tdx_x_x_0, tdx_x_xx_0, tdx_x_xy_0, tdx_x_xz_0, tdx_x_y_0, tdx_x_yy_0, \
                                     tdx_x_yz_0, tdx_x_z_0, tdx_x_zz_0, tdx_xx_xx_0, tdx_xx_xy_0, tdx_xx_xz_0, \
                                     tdx_xx_yy_0, tdx_xx_yz_0, tdx_xx_zz_0, tdx_xy_xx_0, tdx_xy_xy_0, tdx_xy_xz_0, \
                                     tdx_xy_yy_0, tdx_xy_yz_0, tdx_xy_zz_0, tdx_y_x_0, tdx_y_xx_0, tdx_y_xy_0, \
                                     tdx_y_xz_0, tdx_y_y_0, tdx_y_yy_0, tdx_y_yz_0, tdx_y_z_0, tdx_y_zz_0, tdy_0_xx_0, \
                                     tdy_0_xy_0, tdy_0_xz_0, tdy_0_yy_0, tdy_0_yz_0, tdy_0_zz_0, tdy_x_x_0, tdy_x_xx_0, \
                                     tdy_x_xy_0, tdy_x_xz_0, tdy_x_y_0, tdy_x_yy_0, tdy_x_yz_0, tdy_x_z_0, tdy_x_zz_0, \
                                     tdy_xx_xx_0, tdy_xx_xy_0, tdy_xx_xz_0, tdy_xx_yy_0, tdy_xx_yz_0, tdy_xx_zz_0, \
                                     tdy_xy_xx_0, tdy_xy_xy_0, tdy_xy_xz_0, tdy_xy_yy_0, tdy_xy_yz_0, tdy_xy_zz_0, \
                                     tdy_y_x_0, tdy_y_xx_0, tdy_y_xy_0, tdy_y_xz_0, tdy_y_y_0, tdy_y_yy_0, tdy_y_yz_0, \
                                     tdy_y_z_0, tdy_y_zz_0, tdz_0_xx_0, tdz_0_xy_0, tdz_0_xz_0, tdz_0_yy_0, tdz_0_yz_0, \
                                     tdz_0_zz_0, tdz_x_x_0, tdz_x_xx_0, tdz_x_xy_0, tdz_x_xz_0, tdz_x_y_0, tdz_x_yy_0, \
                                     tdz_x_yz_0, tdz_x_z_0, tdz_x_zz_0, tdz_xx_xx_0, tdz_xx_xy_0, tdz_xx_xz_0, \
                                     tdz_xx_yy_0, tdz_xx_yz_0, tdz_xx_zz_0, tdz_xy_xx_0, tdz_xy_xy_0, tdz_xy_xz_0, \
                                     tdz_xy_yy_0, tdz_xy_yz_0, tdz_xy_zz_0, tdz_y_x_0, tdz_y_xx_0, tdz_y_xy_0, \
                                     tdz_y_xz_0, tdz_y_y_0, tdz_y_yy_0, tdz_y_yz_0, tdz_y_z_0, tdz_y_zz_0, ts_x_xx_0, \
                                     ts_x_xy_0, ts_x_xz_0, ts_x_yy_0, ts_x_yz_0, ts_x_zz_0, ts_y_xx_0, ts_y_xy_0, \
                                     ts_y_xz_0, ts_y_yy_0, ts_y_yz_0, ts_y_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xx_xx_0[j] = pa_x[j] * tdx_x_xx_0[j] + 0.5 * fl1_fx * tdx_0_xx_0[j] + fl1_fx * tdx_x_x_0[j] + 0.5 * fl1_fx * ts_x_xx_0[j];

                tdy_xx_xx_0[j] = pa_x[j] * tdy_x_xx_0[j] + 0.5 * fl1_fx * tdy_0_xx_0[j] + fl1_fx * tdy_x_x_0[j];

                tdz_xx_xx_0[j] = pa_x[j] * tdz_x_xx_0[j] + 0.5 * fl1_fx * tdz_0_xx_0[j] + fl1_fx * tdz_x_x_0[j];

                tdx_xx_xy_0[j] = pa_x[j] * tdx_x_xy_0[j] + 0.5 * fl1_fx * tdx_0_xy_0[j] + 0.5 * fl1_fx * tdx_x_y_0[j] + 0.5 * fl1_fx * ts_x_xy_0[j];

                tdy_xx_xy_0[j] = pa_x[j] * tdy_x_xy_0[j] + 0.5 * fl1_fx * tdy_0_xy_0[j] + 0.5 * fl1_fx * tdy_x_y_0[j];

                tdz_xx_xy_0[j] = pa_x[j] * tdz_x_xy_0[j] + 0.5 * fl1_fx * tdz_0_xy_0[j] + 0.5 * fl1_fx * tdz_x_y_0[j];

                tdx_xx_xz_0[j] = pa_x[j] * tdx_x_xz_0[j] + 0.5 * fl1_fx * tdx_0_xz_0[j] + 0.5 * fl1_fx * tdx_x_z_0[j] + 0.5 * fl1_fx * ts_x_xz_0[j];

                tdy_xx_xz_0[j] = pa_x[j] * tdy_x_xz_0[j] + 0.5 * fl1_fx * tdy_0_xz_0[j] + 0.5 * fl1_fx * tdy_x_z_0[j];

                tdz_xx_xz_0[j] = pa_x[j] * tdz_x_xz_0[j] + 0.5 * fl1_fx * tdz_0_xz_0[j] + 0.5 * fl1_fx * tdz_x_z_0[j];

                tdx_xx_yy_0[j] = pa_x[j] * tdx_x_yy_0[j] + 0.5 * fl1_fx * tdx_0_yy_0[j] + 0.5 * fl1_fx * ts_x_yy_0[j];

                tdy_xx_yy_0[j] = pa_x[j] * tdy_x_yy_0[j] + 0.5 * fl1_fx * tdy_0_yy_0[j];

                tdz_xx_yy_0[j] = pa_x[j] * tdz_x_yy_0[j] + 0.5 * fl1_fx * tdz_0_yy_0[j];

                tdx_xx_yz_0[j] = pa_x[j] * tdx_x_yz_0[j] + 0.5 * fl1_fx * tdx_0_yz_0[j] + 0.5 * fl1_fx * ts_x_yz_0[j];

                tdy_xx_yz_0[j] = pa_x[j] * tdy_x_yz_0[j] + 0.5 * fl1_fx * tdy_0_yz_0[j];

                tdz_xx_yz_0[j] = pa_x[j] * tdz_x_yz_0[j] + 0.5 * fl1_fx * tdz_0_yz_0[j];

                tdx_xx_zz_0[j] = pa_x[j] * tdx_x_zz_0[j] + 0.5 * fl1_fx * tdx_0_zz_0[j] + 0.5 * fl1_fx * ts_x_zz_0[j];

                tdy_xx_zz_0[j] = pa_x[j] * tdy_x_zz_0[j] + 0.5 * fl1_fx * tdy_0_zz_0[j];

                tdz_xx_zz_0[j] = pa_x[j] * tdz_x_zz_0[j] + 0.5 * fl1_fx * tdz_0_zz_0[j];

                tdx_xy_xx_0[j] = pa_x[j] * tdx_y_xx_0[j] + fl1_fx * tdx_y_x_0[j] + 0.5 * fl1_fx * ts_y_xx_0[j];

                tdy_xy_xx_0[j] = pa_x[j] * tdy_y_xx_0[j] + fl1_fx * tdy_y_x_0[j];

                tdz_xy_xx_0[j] = pa_x[j] * tdz_y_xx_0[j] + fl1_fx * tdz_y_x_0[j];

                tdx_xy_xy_0[j] = pa_x[j] * tdx_y_xy_0[j] + 0.5 * fl1_fx * tdx_y_y_0[j] + 0.5 * fl1_fx * ts_y_xy_0[j];

                tdy_xy_xy_0[j] = pa_x[j] * tdy_y_xy_0[j] + 0.5 * fl1_fx * tdy_y_y_0[j];

                tdz_xy_xy_0[j] = pa_x[j] * tdz_y_xy_0[j] + 0.5 * fl1_fx * tdz_y_y_0[j];

                tdx_xy_xz_0[j] = pa_x[j] * tdx_y_xz_0[j] + 0.5 * fl1_fx * tdx_y_z_0[j] + 0.5 * fl1_fx * ts_y_xz_0[j];

                tdy_xy_xz_0[j] = pa_x[j] * tdy_y_xz_0[j] + 0.5 * fl1_fx * tdy_y_z_0[j];

                tdz_xy_xz_0[j] = pa_x[j] * tdz_y_xz_0[j] + 0.5 * fl1_fx * tdz_y_z_0[j];

                tdx_xy_yy_0[j] = pa_x[j] * tdx_y_yy_0[j] + 0.5 * fl1_fx * ts_y_yy_0[j];

                tdy_xy_yy_0[j] = pa_x[j] * tdy_y_yy_0[j];

                tdz_xy_yy_0[j] = pa_x[j] * tdz_y_yy_0[j];

                tdx_xy_yz_0[j] = pa_x[j] * tdx_y_yz_0[j] + 0.5 * fl1_fx * ts_y_yz_0[j];

                tdy_xy_yz_0[j] = pa_x[j] * tdy_y_yz_0[j];

                tdz_xy_yz_0[j] = pa_x[j] * tdz_y_yz_0[j];

                tdx_xy_zz_0[j] = pa_x[j] * tdx_y_zz_0[j] + 0.5 * fl1_fx * ts_y_zz_0[j];

                tdy_xy_zz_0[j] = pa_x[j] * tdy_y_zz_0[j];

                tdz_xy_zz_0[j] = pa_x[j] * tdz_y_zz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDD_36_72(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
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

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tdx_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 6); 

            auto tdy_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tdz_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tdx_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 7); 

            auto tdy_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tdz_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tdx_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 8); 

            auto tdy_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tdz_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tdx_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 9); 

            auto tdy_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tdz_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdx_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 10); 

            auto tdy_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tdz_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tdx_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 11); 

            auto tdy_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tdz_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 11); 

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

            auto tdx_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx); 

            auto tdy_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx); 

            auto tdz_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx); 

            auto tdx_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 1); 

            auto tdy_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 1); 

            auto tdz_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 1); 

            auto tdx_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 2); 

            auto tdy_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 2); 

            auto tdz_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 2); 

            auto tdx_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 3); 

            auto tdy_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 3); 

            auto tdz_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 3); 

            auto tdx_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 4); 

            auto tdy_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 4); 

            auto tdz_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 4); 

            auto tdx_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 5); 

            auto tdy_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 5); 

            auto tdz_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 5); 

            auto tdx_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 3); 

            auto tdy_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 3); 

            auto tdz_y_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 3); 

            auto tdx_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 4); 

            auto tdy_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 4); 

            auto tdz_y_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 4); 

            auto tdx_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 5); 

            auto tdy_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 5); 

            auto tdz_y_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 5); 

            auto tdx_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 6); 

            auto tdy_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 6); 

            auto tdz_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 6); 

            auto tdx_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 7); 

            auto tdy_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 7); 

            auto tdz_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 7); 

            auto tdx_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 8); 

            auto tdy_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 8); 

            auto tdz_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 8); 

            auto ts_y_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 6); 

            auto ts_y_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 7); 

            auto ts_y_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 8); 

            auto ts_y_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 9); 

            auto ts_y_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 10); 

            auto ts_y_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 11); 

            auto ts_z_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 12); 

            auto ts_z_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 13); 

            auto ts_z_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 14); 

            auto ts_z_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 15); 

            auto ts_z_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 16); 

            auto ts_z_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 17); 

            // set up pointers to integrals

            auto tdx_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 12); 

            auto tdy_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 12); 

            auto tdz_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 12); 

            auto tdx_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 13); 

            auto tdy_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 13); 

            auto tdz_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 13); 

            auto tdx_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 14); 

            auto tdy_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 14); 

            auto tdz_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 14); 

            auto tdx_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 15); 

            auto tdy_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 15); 

            auto tdz_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 15); 

            auto tdx_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 16); 

            auto tdy_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 16); 

            auto tdz_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 16); 

            auto tdx_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 17); 

            auto tdy_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 17); 

            auto tdz_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 17); 

            auto tdx_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 18); 

            auto tdy_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tdz_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tdx_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 19); 

            auto tdy_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tdz_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tdx_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 20); 

            auto tdy_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tdz_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tdx_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 21); 

            auto tdy_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tdz_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tdx_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 22); 

            auto tdy_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tdz_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tdx_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 23); 

            auto tdy_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tdz_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            // Batch of Integrals (36,72)

            #pragma omp simd aligned(fx, pa_x, pa_y, tdx_0_xx_0, tdx_0_xy_0, tdx_0_xz_0, tdx_0_yy_0, tdx_0_yz_0, \
                                     tdx_0_zz_0, tdx_xz_xx_0, tdx_xz_xy_0, tdx_xz_xz_0, tdx_xz_yy_0, tdx_xz_yz_0, \
                                     tdx_xz_zz_0, tdx_y_x_0, tdx_y_xx_0, tdx_y_xy_0, tdx_y_xz_0, tdx_y_y_0, tdx_y_yy_0, \
                                     tdx_y_yz_0, tdx_y_z_0, tdx_y_zz_0, tdx_yy_xx_0, tdx_yy_xy_0, tdx_yy_xz_0, \
                                     tdx_yy_yy_0, tdx_yy_yz_0, tdx_yy_zz_0, tdx_z_x_0, tdx_z_xx_0, tdx_z_xy_0, \
                                     tdx_z_xz_0, tdx_z_y_0, tdx_z_yy_0, tdx_z_yz_0, tdx_z_z_0, tdx_z_zz_0, tdy_0_xx_0, \
                                     tdy_0_xy_0, tdy_0_xz_0, tdy_0_yy_0, tdy_0_yz_0, tdy_0_zz_0, tdy_xz_xx_0, \
                                     tdy_xz_xy_0, tdy_xz_xz_0, tdy_xz_yy_0, tdy_xz_yz_0, tdy_xz_zz_0, tdy_y_x_0, \
                                     tdy_y_xx_0, tdy_y_xy_0, tdy_y_xz_0, tdy_y_y_0, tdy_y_yy_0, tdy_y_yz_0, tdy_y_z_0, \
                                     tdy_y_zz_0, tdy_yy_xx_0, tdy_yy_xy_0, tdy_yy_xz_0, tdy_yy_yy_0, tdy_yy_yz_0, \
                                     tdy_yy_zz_0, tdy_z_x_0, tdy_z_xx_0, tdy_z_xy_0, tdy_z_xz_0, tdy_z_y_0, tdy_z_yy_0, \
                                     tdy_z_yz_0, tdy_z_z_0, tdy_z_zz_0, tdz_0_xx_0, tdz_0_xy_0, tdz_0_xz_0, tdz_0_yy_0, \
                                     tdz_0_yz_0, tdz_0_zz_0, tdz_xz_xx_0, tdz_xz_xy_0, tdz_xz_xz_0, tdz_xz_yy_0, \
                                     tdz_xz_yz_0, tdz_xz_zz_0, tdz_y_x_0, tdz_y_xx_0, tdz_y_xy_0, tdz_y_xz_0, tdz_y_y_0, \
                                     tdz_y_yy_0, tdz_y_yz_0, tdz_y_z_0, tdz_y_zz_0, tdz_yy_xx_0, tdz_yy_xy_0, \
                                     tdz_yy_xz_0, tdz_yy_yy_0, tdz_yy_yz_0, tdz_yy_zz_0, tdz_z_x_0, tdz_z_xx_0, \
                                     tdz_z_xy_0, tdz_z_xz_0, tdz_z_y_0, tdz_z_yy_0, tdz_z_yz_0, tdz_z_z_0, tdz_z_zz_0, \
                                     ts_y_xx_0, ts_y_xy_0, ts_y_xz_0, ts_y_yy_0, ts_y_yz_0, ts_y_zz_0, ts_z_xx_0, \
                                     ts_z_xy_0, ts_z_xz_0, ts_z_yy_0, ts_z_yz_0, ts_z_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xz_xx_0[j] = pa_x[j] * tdx_z_xx_0[j] + fl1_fx * tdx_z_x_0[j] + 0.5 * fl1_fx * ts_z_xx_0[j];

                tdy_xz_xx_0[j] = pa_x[j] * tdy_z_xx_0[j] + fl1_fx * tdy_z_x_0[j];

                tdz_xz_xx_0[j] = pa_x[j] * tdz_z_xx_0[j] + fl1_fx * tdz_z_x_0[j];

                tdx_xz_xy_0[j] = pa_x[j] * tdx_z_xy_0[j] + 0.5 * fl1_fx * tdx_z_y_0[j] + 0.5 * fl1_fx * ts_z_xy_0[j];

                tdy_xz_xy_0[j] = pa_x[j] * tdy_z_xy_0[j] + 0.5 * fl1_fx * tdy_z_y_0[j];

                tdz_xz_xy_0[j] = pa_x[j] * tdz_z_xy_0[j] + 0.5 * fl1_fx * tdz_z_y_0[j];

                tdx_xz_xz_0[j] = pa_x[j] * tdx_z_xz_0[j] + 0.5 * fl1_fx * tdx_z_z_0[j] + 0.5 * fl1_fx * ts_z_xz_0[j];

                tdy_xz_xz_0[j] = pa_x[j] * tdy_z_xz_0[j] + 0.5 * fl1_fx * tdy_z_z_0[j];

                tdz_xz_xz_0[j] = pa_x[j] * tdz_z_xz_0[j] + 0.5 * fl1_fx * tdz_z_z_0[j];

                tdx_xz_yy_0[j] = pa_x[j] * tdx_z_yy_0[j] + 0.5 * fl1_fx * ts_z_yy_0[j];

                tdy_xz_yy_0[j] = pa_x[j] * tdy_z_yy_0[j];

                tdz_xz_yy_0[j] = pa_x[j] * tdz_z_yy_0[j];

                tdx_xz_yz_0[j] = pa_x[j] * tdx_z_yz_0[j] + 0.5 * fl1_fx * ts_z_yz_0[j];

                tdy_xz_yz_0[j] = pa_x[j] * tdy_z_yz_0[j];

                tdz_xz_yz_0[j] = pa_x[j] * tdz_z_yz_0[j];

                tdx_xz_zz_0[j] = pa_x[j] * tdx_z_zz_0[j] + 0.5 * fl1_fx * ts_z_zz_0[j];

                tdy_xz_zz_0[j] = pa_x[j] * tdy_z_zz_0[j];

                tdz_xz_zz_0[j] = pa_x[j] * tdz_z_zz_0[j];

                tdx_yy_xx_0[j] = pa_y[j] * tdx_y_xx_0[j] + 0.5 * fl1_fx * tdx_0_xx_0[j];

                tdy_yy_xx_0[j] = pa_y[j] * tdy_y_xx_0[j] + 0.5 * fl1_fx * tdy_0_xx_0[j] + 0.5 * fl1_fx * ts_y_xx_0[j];

                tdz_yy_xx_0[j] = pa_y[j] * tdz_y_xx_0[j] + 0.5 * fl1_fx * tdz_0_xx_0[j];

                tdx_yy_xy_0[j] = pa_y[j] * tdx_y_xy_0[j] + 0.5 * fl1_fx * tdx_0_xy_0[j] + 0.5 * fl1_fx * tdx_y_x_0[j];

                tdy_yy_xy_0[j] = pa_y[j] * tdy_y_xy_0[j] + 0.5 * fl1_fx * tdy_0_xy_0[j] + 0.5 * fl1_fx * tdy_y_x_0[j] + 0.5 * fl1_fx * ts_y_xy_0[j];

                tdz_yy_xy_0[j] = pa_y[j] * tdz_y_xy_0[j] + 0.5 * fl1_fx * tdz_0_xy_0[j] + 0.5 * fl1_fx * tdz_y_x_0[j];

                tdx_yy_xz_0[j] = pa_y[j] * tdx_y_xz_0[j] + 0.5 * fl1_fx * tdx_0_xz_0[j];

                tdy_yy_xz_0[j] = pa_y[j] * tdy_y_xz_0[j] + 0.5 * fl1_fx * tdy_0_xz_0[j] + 0.5 * fl1_fx * ts_y_xz_0[j];

                tdz_yy_xz_0[j] = pa_y[j] * tdz_y_xz_0[j] + 0.5 * fl1_fx * tdz_0_xz_0[j];

                tdx_yy_yy_0[j] = pa_y[j] * tdx_y_yy_0[j] + 0.5 * fl1_fx * tdx_0_yy_0[j] + fl1_fx * tdx_y_y_0[j];

                tdy_yy_yy_0[j] = pa_y[j] * tdy_y_yy_0[j] + 0.5 * fl1_fx * tdy_0_yy_0[j] + fl1_fx * tdy_y_y_0[j] + 0.5 * fl1_fx * ts_y_yy_0[j];

                tdz_yy_yy_0[j] = pa_y[j] * tdz_y_yy_0[j] + 0.5 * fl1_fx * tdz_0_yy_0[j] + fl1_fx * tdz_y_y_0[j];

                tdx_yy_yz_0[j] = pa_y[j] * tdx_y_yz_0[j] + 0.5 * fl1_fx * tdx_0_yz_0[j] + 0.5 * fl1_fx * tdx_y_z_0[j];

                tdy_yy_yz_0[j] = pa_y[j] * tdy_y_yz_0[j] + 0.5 * fl1_fx * tdy_0_yz_0[j] + 0.5 * fl1_fx * tdy_y_z_0[j] + 0.5 * fl1_fx * ts_y_yz_0[j];

                tdz_yy_yz_0[j] = pa_y[j] * tdz_y_yz_0[j] + 0.5 * fl1_fx * tdz_0_yz_0[j] + 0.5 * fl1_fx * tdz_y_z_0[j];

                tdx_yy_zz_0[j] = pa_y[j] * tdx_y_zz_0[j] + 0.5 * fl1_fx * tdx_0_zz_0[j];

                tdy_yy_zz_0[j] = pa_y[j] * tdy_y_zz_0[j] + 0.5 * fl1_fx * tdy_0_zz_0[j] + 0.5 * fl1_fx * ts_y_zz_0[j];

                tdz_yy_zz_0[j] = pa_y[j] * tdz_y_zz_0[j] + 0.5 * fl1_fx * tdz_0_zz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDD_72_108(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
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

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto tdx_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx); 

            auto tdy_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx); 

            auto tdz_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx); 

            auto tdx_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 1); 

            auto tdy_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 1); 

            auto tdz_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 1); 

            auto tdx_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 2); 

            auto tdy_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 2); 

            auto tdz_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 2); 

            auto tdx_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 3); 

            auto tdy_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 3); 

            auto tdz_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 3); 

            auto tdx_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 4); 

            auto tdy_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 4); 

            auto tdz_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 4); 

            auto tdx_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 5); 

            auto tdy_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 5); 

            auto tdz_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 5); 

            auto tdx_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 6); 

            auto tdy_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 6); 

            auto tdz_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 6); 

            auto tdx_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 7); 

            auto tdy_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 7); 

            auto tdz_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 7); 

            auto tdx_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 8); 

            auto tdy_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 8); 

            auto tdz_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 8); 

            auto ts_z_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 12); 

            auto ts_z_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 13); 

            auto ts_z_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 14); 

            auto ts_z_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 15); 

            auto ts_z_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 16); 

            auto ts_z_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 17); 

            // set up pointers to integrals

            auto tdx_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 24); 

            auto tdy_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tdz_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tdx_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 25); 

            auto tdy_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tdz_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tdx_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 26); 

            auto tdy_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tdz_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tdx_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 27); 

            auto tdy_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tdz_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tdx_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 28); 

            auto tdy_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tdz_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tdx_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 29); 

            auto tdy_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 29); 

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

            // Batch of Integrals (72,108)

            #pragma omp simd aligned(fx, pa_y, pa_z, tdx_0_xx_0, tdx_0_xy_0, tdx_0_xz_0, tdx_0_yy_0, tdx_0_yz_0, \
                                     tdx_0_zz_0, tdx_yz_xx_0, tdx_yz_xy_0, tdx_yz_xz_0, tdx_yz_yy_0, tdx_yz_yz_0, \
                                     tdx_yz_zz_0, tdx_z_x_0, tdx_z_xx_0, tdx_z_xy_0, tdx_z_xz_0, tdx_z_y_0, tdx_z_yy_0, \
                                     tdx_z_yz_0, tdx_z_z_0, tdx_z_zz_0, tdx_zz_xx_0, tdx_zz_xy_0, tdx_zz_xz_0, \
                                     tdx_zz_yy_0, tdx_zz_yz_0, tdx_zz_zz_0, tdy_0_xx_0, tdy_0_xy_0, tdy_0_xz_0, \
                                     tdy_0_yy_0, tdy_0_yz_0, tdy_0_zz_0, tdy_yz_xx_0, tdy_yz_xy_0, tdy_yz_xz_0, \
                                     tdy_yz_yy_0, tdy_yz_yz_0, tdy_yz_zz_0, tdy_z_x_0, tdy_z_xx_0, tdy_z_xy_0, \
                                     tdy_z_xz_0, tdy_z_y_0, tdy_z_yy_0, tdy_z_yz_0, tdy_z_z_0, tdy_z_zz_0, tdy_zz_xx_0, \
                                     tdy_zz_xy_0, tdy_zz_xz_0, tdy_zz_yy_0, tdy_zz_yz_0, tdy_zz_zz_0, tdz_0_xx_0, \
                                     tdz_0_xy_0, tdz_0_xz_0, tdz_0_yy_0, tdz_0_yz_0, tdz_0_zz_0, tdz_yz_xx_0, \
                                     tdz_yz_xy_0, tdz_yz_xz_0, tdz_yz_yy_0, tdz_yz_yz_0, tdz_yz_zz_0, tdz_z_x_0, \
                                     tdz_z_xx_0, tdz_z_xy_0, tdz_z_xz_0, tdz_z_y_0, tdz_z_yy_0, tdz_z_yz_0, tdz_z_z_0, \
                                     tdz_z_zz_0, tdz_zz_xx_0, tdz_zz_xy_0, tdz_zz_xz_0, tdz_zz_yy_0, tdz_zz_yz_0, \
                                     tdz_zz_zz_0, ts_z_xx_0, ts_z_xy_0, ts_z_xz_0, ts_z_yy_0, ts_z_yz_0, ts_z_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yz_xx_0[j] = pa_y[j] * tdx_z_xx_0[j];

                tdy_yz_xx_0[j] = pa_y[j] * tdy_z_xx_0[j] + 0.5 * fl1_fx * ts_z_xx_0[j];

                tdz_yz_xx_0[j] = pa_y[j] * tdz_z_xx_0[j];

                tdx_yz_xy_0[j] = pa_y[j] * tdx_z_xy_0[j] + 0.5 * fl1_fx * tdx_z_x_0[j];

                tdy_yz_xy_0[j] = pa_y[j] * tdy_z_xy_0[j] + 0.5 * fl1_fx * tdy_z_x_0[j] + 0.5 * fl1_fx * ts_z_xy_0[j];

                tdz_yz_xy_0[j] = pa_y[j] * tdz_z_xy_0[j] + 0.5 * fl1_fx * tdz_z_x_0[j];

                tdx_yz_xz_0[j] = pa_y[j] * tdx_z_xz_0[j];

                tdy_yz_xz_0[j] = pa_y[j] * tdy_z_xz_0[j] + 0.5 * fl1_fx * ts_z_xz_0[j];

                tdz_yz_xz_0[j] = pa_y[j] * tdz_z_xz_0[j];

                tdx_yz_yy_0[j] = pa_y[j] * tdx_z_yy_0[j] + fl1_fx * tdx_z_y_0[j];

                tdy_yz_yy_0[j] = pa_y[j] * tdy_z_yy_0[j] + fl1_fx * tdy_z_y_0[j] + 0.5 * fl1_fx * ts_z_yy_0[j];

                tdz_yz_yy_0[j] = pa_y[j] * tdz_z_yy_0[j] + fl1_fx * tdz_z_y_0[j];

                tdx_yz_yz_0[j] = pa_y[j] * tdx_z_yz_0[j] + 0.5 * fl1_fx * tdx_z_z_0[j];

                tdy_yz_yz_0[j] = pa_y[j] * tdy_z_yz_0[j] + 0.5 * fl1_fx * tdy_z_z_0[j] + 0.5 * fl1_fx * ts_z_yz_0[j];

                tdz_yz_yz_0[j] = pa_y[j] * tdz_z_yz_0[j] + 0.5 * fl1_fx * tdz_z_z_0[j];

                tdx_yz_zz_0[j] = pa_y[j] * tdx_z_zz_0[j];

                tdy_yz_zz_0[j] = pa_y[j] * tdy_z_zz_0[j] + 0.5 * fl1_fx * ts_z_zz_0[j];

                tdz_yz_zz_0[j] = pa_y[j] * tdz_z_zz_0[j];

                tdx_zz_xx_0[j] = pa_z[j] * tdx_z_xx_0[j] + 0.5 * fl1_fx * tdx_0_xx_0[j];

                tdy_zz_xx_0[j] = pa_z[j] * tdy_z_xx_0[j] + 0.5 * fl1_fx * tdy_0_xx_0[j];

                tdz_zz_xx_0[j] = pa_z[j] * tdz_z_xx_0[j] + 0.5 * fl1_fx * tdz_0_xx_0[j] + 0.5 * fl1_fx * ts_z_xx_0[j];

                tdx_zz_xy_0[j] = pa_z[j] * tdx_z_xy_0[j] + 0.5 * fl1_fx * tdx_0_xy_0[j];

                tdy_zz_xy_0[j] = pa_z[j] * tdy_z_xy_0[j] + 0.5 * fl1_fx * tdy_0_xy_0[j];

                tdz_zz_xy_0[j] = pa_z[j] * tdz_z_xy_0[j] + 0.5 * fl1_fx * tdz_0_xy_0[j] + 0.5 * fl1_fx * ts_z_xy_0[j];

                tdx_zz_xz_0[j] = pa_z[j] * tdx_z_xz_0[j] + 0.5 * fl1_fx * tdx_0_xz_0[j] + 0.5 * fl1_fx * tdx_z_x_0[j];

                tdy_zz_xz_0[j] = pa_z[j] * tdy_z_xz_0[j] + 0.5 * fl1_fx * tdy_0_xz_0[j] + 0.5 * fl1_fx * tdy_z_x_0[j];

                tdz_zz_xz_0[j] = pa_z[j] * tdz_z_xz_0[j] + 0.5 * fl1_fx * tdz_0_xz_0[j] + 0.5 * fl1_fx * tdz_z_x_0[j] + 0.5 * fl1_fx * ts_z_xz_0[j];

                tdx_zz_yy_0[j] = pa_z[j] * tdx_z_yy_0[j] + 0.5 * fl1_fx * tdx_0_yy_0[j];

                tdy_zz_yy_0[j] = pa_z[j] * tdy_z_yy_0[j] + 0.5 * fl1_fx * tdy_0_yy_0[j];

                tdz_zz_yy_0[j] = pa_z[j] * tdz_z_yy_0[j] + 0.5 * fl1_fx * tdz_0_yy_0[j] + 0.5 * fl1_fx * ts_z_yy_0[j];

                tdx_zz_yz_0[j] = pa_z[j] * tdx_z_yz_0[j] + 0.5 * fl1_fx * tdx_0_yz_0[j] + 0.5 * fl1_fx * tdx_z_y_0[j];

                tdy_zz_yz_0[j] = pa_z[j] * tdy_z_yz_0[j] + 0.5 * fl1_fx * tdy_0_yz_0[j] + 0.5 * fl1_fx * tdy_z_y_0[j];

                tdz_zz_yz_0[j] = pa_z[j] * tdz_z_yz_0[j] + 0.5 * fl1_fx * tdz_0_yz_0[j] + 0.5 * fl1_fx * tdz_z_y_0[j] + 0.5 * fl1_fx * ts_z_yz_0[j];

                tdx_zz_zz_0[j] = pa_z[j] * tdx_z_zz_0[j] + 0.5 * fl1_fx * tdx_0_zz_0[j] + fl1_fx * tdx_z_z_0[j];

                tdy_zz_zz_0[j] = pa_z[j] * tdy_z_zz_0[j] + 0.5 * fl1_fx * tdy_0_zz_0[j] + fl1_fx * tdy_z_z_0[j];

                tdz_zz_zz_0[j] = pa_z[j] * tdz_z_zz_0[j] + 0.5 * fl1_fx * tdz_0_zz_0[j] + fl1_fx * tdz_z_z_0[j] + 0.5 * fl1_fx * ts_z_zz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDF(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForDF_0_45(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  nOSFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForDF_45_90(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   nOSFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        ediprecfunc::compElectricDipoleForDF_90_135(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    nOSFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        ediprecfunc::compElectricDipoleForDF_135_180(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     nOSFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 
    }

    void
    compElectricDipoleForDF_0_45(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const int32_t              nOSFactors,
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

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx); 

            auto tdy_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx); 

            auto tdz_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx); 

            auto tdx_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 1); 

            auto tdy_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 1); 

            auto tdz_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 1); 

            auto tdx_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 2); 

            auto tdy_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 2); 

            auto tdz_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 2); 

            auto tdx_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 3); 

            auto tdy_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 3); 

            auto tdz_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 3); 

            auto tdx_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 4); 

            auto tdy_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 4); 

            auto tdz_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 4); 

            auto tdx_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 5); 

            auto tdy_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 5); 

            auto tdz_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 5); 

            auto tdx_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 6); 

            auto tdy_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 6); 

            auto tdz_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 6); 

            auto tdx_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 7); 

            auto tdy_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 7); 

            auto tdz_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 7); 

            auto tdx_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 8); 

            auto tdy_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 8); 

            auto tdz_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 8); 

            auto tdx_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 9); 

            auto tdy_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 9); 

            auto tdz_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 9); 

            auto tdx_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 10); 

            auto tdy_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tdz_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tdx_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 11); 

            auto tdy_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tdz_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tdx_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 12); 

            auto tdy_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tdz_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tdx_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 13); 

            auto tdy_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tdz_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tdx_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 14); 

            auto tdy_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tdz_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tdx_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx); 

            auto tdy_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx); 

            auto tdz_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx); 

            auto tdx_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 1); 

            auto tdy_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 1); 

            auto tdz_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 1); 

            auto tdx_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 2); 

            auto tdy_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 2); 

            auto tdz_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 2); 

            auto tdx_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 3); 

            auto tdy_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 3); 

            auto tdz_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 3); 

            auto tdx_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 4); 

            auto tdy_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 4); 

            auto tdz_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 4); 

            auto tdx_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 5); 

            auto tdy_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 5); 

            auto tdz_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 5); 

            auto tdx_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 6); 

            auto tdy_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 6); 

            auto tdz_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 6); 

            auto tdx_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 7); 

            auto tdy_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 7); 

            auto tdz_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 7); 

            auto tdx_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 8); 

            auto tdy_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 8); 

            auto tdz_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 8); 

            auto tdx_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 9); 

            auto tdy_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 9); 

            auto tdz_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 9); 

            auto tdx_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx); 

            auto tdy_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx); 

            auto tdz_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx); 

            auto tdx_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 1); 

            auto tdy_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 1); 

            auto tdz_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 1); 

            auto tdx_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 2); 

            auto tdy_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 2); 

            auto tdz_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 2); 

            auto tdx_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 3); 

            auto tdy_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 3); 

            auto tdz_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 3); 

            auto tdx_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 4); 

            auto tdy_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 4); 

            auto tdz_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 4); 

            auto tdx_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 5); 

            auto tdy_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 5); 

            auto tdz_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 5); 

            auto tdx_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 6); 

            auto tdy_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tdz_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tdx_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 7); 

            auto tdy_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tdz_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tdx_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 8); 

            auto tdy_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tdz_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tdx_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 9); 

            auto tdy_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tdz_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdx_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 10); 

            auto tdy_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tdz_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto ts_x_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx); 

            auto ts_x_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 1); 

            auto ts_x_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 2); 

            auto ts_x_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 3); 

            auto ts_x_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 4); 

            auto ts_x_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 5); 

            auto ts_x_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 6); 

            auto ts_x_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 7); 

            auto ts_x_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 8); 

            auto ts_x_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 9); 

            auto ts_y_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 10); 

            auto ts_y_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 11); 

            auto ts_y_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 12); 

            auto ts_y_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 13); 

            auto ts_y_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 14); 

            // set up pointers to integrals

            auto tdx_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx); 

            auto tdy_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx); 

            auto tdz_xx_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx); 

            auto tdx_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 1); 

            auto tdy_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 1); 

            auto tdz_xx_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 1); 

            auto tdx_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 2); 

            auto tdy_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 2); 

            auto tdz_xx_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 2); 

            auto tdx_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 3); 

            auto tdy_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 3); 

            auto tdz_xx_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 3); 

            auto tdx_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 4); 

            auto tdy_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 4); 

            auto tdz_xx_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 4); 

            auto tdx_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 5); 

            auto tdy_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 5); 

            auto tdz_xx_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 5); 

            auto tdx_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 6); 

            auto tdy_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 6); 

            auto tdz_xx_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 6); 

            auto tdx_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 7); 

            auto tdy_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 7); 

            auto tdz_xx_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 7); 

            auto tdx_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 8); 

            auto tdy_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 8); 

            auto tdz_xx_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 8); 

            auto tdx_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 9); 

            auto tdy_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 9); 

            auto tdz_xx_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 9); 

            auto tdx_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 10); 

            auto tdy_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 10); 

            auto tdz_xy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 10); 

            auto tdx_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 11); 

            auto tdy_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 11); 

            auto tdz_xy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 11); 

            auto tdx_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 12); 

            auto tdy_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 12); 

            auto tdz_xy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 12); 

            auto tdx_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 13); 

            auto tdy_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 13); 

            auto tdz_xy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 13); 

            auto tdx_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 14); 

            auto tdy_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 14); 

            auto tdz_xy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 14); 

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, tdx_0_xxx_0, tdx_0_xxy_0, tdx_0_xxz_0, tdx_0_xyy_0, \
                                     tdx_0_xyz_0, tdx_0_xzz_0, tdx_0_yyy_0, tdx_0_yyz_0, tdx_0_yzz_0, tdx_0_zzz_0, \
                                     tdx_x_xx_0, tdx_x_xxx_0, tdx_x_xxy_0, tdx_x_xxz_0, tdx_x_xy_0, tdx_x_xyy_0, \
                                     tdx_x_xyz_0, tdx_x_xz_0, tdx_x_xzz_0, tdx_x_yy_0, tdx_x_yyy_0, tdx_x_yyz_0, \
                                     tdx_x_yz_0, tdx_x_yzz_0, tdx_x_zz_0, tdx_x_zzz_0, tdx_xx_xxx_0, tdx_xx_xxy_0, \
                                     tdx_xx_xxz_0, tdx_xx_xyy_0, tdx_xx_xyz_0, tdx_xx_xzz_0, tdx_xx_yyy_0, tdx_xx_yyz_0, \
                                     tdx_xx_yzz_0, tdx_xx_zzz_0, tdx_xy_xxx_0, tdx_xy_xxy_0, tdx_xy_xxz_0, tdx_xy_xyy_0, \
                                     tdx_xy_xyz_0, tdx_y_xx_0, tdx_y_xxx_0, tdx_y_xxy_0, tdx_y_xxz_0, tdx_y_xy_0, \
                                     tdx_y_xyy_0, tdx_y_xyz_0, tdx_y_xz_0, tdx_y_yy_0, tdx_y_yz_0, tdy_0_xxx_0, \
                                     tdy_0_xxy_0, tdy_0_xxz_0, tdy_0_xyy_0, tdy_0_xyz_0, tdy_0_xzz_0, tdy_0_yyy_0, \
                                     tdy_0_yyz_0, tdy_0_yzz_0, tdy_0_zzz_0, tdy_x_xx_0, tdy_x_xxx_0, tdy_x_xxy_0, \
                                     tdy_x_xxz_0, tdy_x_xy_0, tdy_x_xyy_0, tdy_x_xyz_0, tdy_x_xz_0, tdy_x_xzz_0, \
                                     tdy_x_yy_0, tdy_x_yyy_0, tdy_x_yyz_0, tdy_x_yz_0, tdy_x_yzz_0, tdy_x_zz_0, \
                                     tdy_x_zzz_0, tdy_xx_xxx_0, tdy_xx_xxy_0, tdy_xx_xxz_0, tdy_xx_xyy_0, tdy_xx_xyz_0, \
                                     tdy_xx_xzz_0, tdy_xx_yyy_0, tdy_xx_yyz_0, tdy_xx_yzz_0, tdy_xx_zzz_0, tdy_xy_xxx_0, \
                                     tdy_xy_xxy_0, tdy_xy_xxz_0, tdy_xy_xyy_0, tdy_xy_xyz_0, tdy_y_xx_0, tdy_y_xxx_0, \
                                     tdy_y_xxy_0, tdy_y_xxz_0, tdy_y_xy_0, tdy_y_xyy_0, tdy_y_xyz_0, tdy_y_xz_0, \
                                     tdy_y_yy_0, tdy_y_yz_0, tdz_0_xxx_0, tdz_0_xxy_0, tdz_0_xxz_0, tdz_0_xyy_0, \
                                     tdz_0_xyz_0, tdz_0_xzz_0, tdz_0_yyy_0, tdz_0_yyz_0, tdz_0_yzz_0, tdz_0_zzz_0, \
                                     tdz_x_xx_0, tdz_x_xxx_0, tdz_x_xxy_0, tdz_x_xxz_0, tdz_x_xy_0, tdz_x_xyy_0, \
                                     tdz_x_xyz_0, tdz_x_xz_0, tdz_x_xzz_0, tdz_x_yy_0, tdz_x_yyy_0, tdz_x_yyz_0, \
                                     tdz_x_yz_0, tdz_x_yzz_0, tdz_x_zz_0, tdz_x_zzz_0, tdz_xx_xxx_0, tdz_xx_xxy_0, \
                                     tdz_xx_xxz_0, tdz_xx_xyy_0, tdz_xx_xyz_0, tdz_xx_xzz_0, tdz_xx_yyy_0, tdz_xx_yyz_0, \
                                     tdz_xx_yzz_0, tdz_xx_zzz_0, tdz_xy_xxx_0, tdz_xy_xxy_0, tdz_xy_xxz_0, tdz_xy_xyy_0, \
                                     tdz_xy_xyz_0, tdz_y_xx_0, tdz_y_xxx_0, tdz_y_xxy_0, tdz_y_xxz_0, tdz_y_xy_0, \
                                     tdz_y_xyy_0, tdz_y_xyz_0, tdz_y_xz_0, tdz_y_yy_0, tdz_y_yz_0, ts_x_xxx_0, \
                                     ts_x_xxy_0, ts_x_xxz_0, ts_x_xyy_0, ts_x_xyz_0, ts_x_xzz_0, ts_x_yyy_0, ts_x_yyz_0, \
                                     ts_x_yzz_0, ts_x_zzz_0, ts_y_xxx_0, ts_y_xxy_0, ts_y_xxz_0, ts_y_xyy_0, ts_y_xyz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xx_xxx_0[j] = pa_x[j] * tdx_x_xxx_0[j] + 0.5 * fl1_fx * tdx_0_xxx_0[j] + 1.5 * fl1_fx * tdx_x_xx_0[j] + 0.5 * fl1_fx * ts_x_xxx_0[j];

                tdy_xx_xxx_0[j] = pa_x[j] * tdy_x_xxx_0[j] + 0.5 * fl1_fx * tdy_0_xxx_0[j] + 1.5 * fl1_fx * tdy_x_xx_0[j];

                tdz_xx_xxx_0[j] = pa_x[j] * tdz_x_xxx_0[j] + 0.5 * fl1_fx * tdz_0_xxx_0[j] + 1.5 * fl1_fx * tdz_x_xx_0[j];

                tdx_xx_xxy_0[j] = pa_x[j] * tdx_x_xxy_0[j] + 0.5 * fl1_fx * tdx_0_xxy_0[j] + fl1_fx * tdx_x_xy_0[j] + 0.5 * fl1_fx * ts_x_xxy_0[j];

                tdy_xx_xxy_0[j] = pa_x[j] * tdy_x_xxy_0[j] + 0.5 * fl1_fx * tdy_0_xxy_0[j] + fl1_fx * tdy_x_xy_0[j];

                tdz_xx_xxy_0[j] = pa_x[j] * tdz_x_xxy_0[j] + 0.5 * fl1_fx * tdz_0_xxy_0[j] + fl1_fx * tdz_x_xy_0[j];

                tdx_xx_xxz_0[j] = pa_x[j] * tdx_x_xxz_0[j] + 0.5 * fl1_fx * tdx_0_xxz_0[j] + fl1_fx * tdx_x_xz_0[j] + 0.5 * fl1_fx * ts_x_xxz_0[j];

                tdy_xx_xxz_0[j] = pa_x[j] * tdy_x_xxz_0[j] + 0.5 * fl1_fx * tdy_0_xxz_0[j] + fl1_fx * tdy_x_xz_0[j];

                tdz_xx_xxz_0[j] = pa_x[j] * tdz_x_xxz_0[j] + 0.5 * fl1_fx * tdz_0_xxz_0[j] + fl1_fx * tdz_x_xz_0[j];

                tdx_xx_xyy_0[j] = pa_x[j] * tdx_x_xyy_0[j] + 0.5 * fl1_fx * tdx_0_xyy_0[j] + 0.5 * fl1_fx * tdx_x_yy_0[j] + 0.5 * fl1_fx * ts_x_xyy_0[j];

                tdy_xx_xyy_0[j] = pa_x[j] * tdy_x_xyy_0[j] + 0.5 * fl1_fx * tdy_0_xyy_0[j] + 0.5 * fl1_fx * tdy_x_yy_0[j];

                tdz_xx_xyy_0[j] = pa_x[j] * tdz_x_xyy_0[j] + 0.5 * fl1_fx * tdz_0_xyy_0[j] + 0.5 * fl1_fx * tdz_x_yy_0[j];

                tdx_xx_xyz_0[j] = pa_x[j] * tdx_x_xyz_0[j] + 0.5 * fl1_fx * tdx_0_xyz_0[j] + 0.5 * fl1_fx * tdx_x_yz_0[j] + 0.5 * fl1_fx * ts_x_xyz_0[j];

                tdy_xx_xyz_0[j] = pa_x[j] * tdy_x_xyz_0[j] + 0.5 * fl1_fx * tdy_0_xyz_0[j] + 0.5 * fl1_fx * tdy_x_yz_0[j];

                tdz_xx_xyz_0[j] = pa_x[j] * tdz_x_xyz_0[j] + 0.5 * fl1_fx * tdz_0_xyz_0[j] + 0.5 * fl1_fx * tdz_x_yz_0[j];

                tdx_xx_xzz_0[j] = pa_x[j] * tdx_x_xzz_0[j] + 0.5 * fl1_fx * tdx_0_xzz_0[j] + 0.5 * fl1_fx * tdx_x_zz_0[j] + 0.5 * fl1_fx * ts_x_xzz_0[j];

                tdy_xx_xzz_0[j] = pa_x[j] * tdy_x_xzz_0[j] + 0.5 * fl1_fx * tdy_0_xzz_0[j] + 0.5 * fl1_fx * tdy_x_zz_0[j];

                tdz_xx_xzz_0[j] = pa_x[j] * tdz_x_xzz_0[j] + 0.5 * fl1_fx * tdz_0_xzz_0[j] + 0.5 * fl1_fx * tdz_x_zz_0[j];

                tdx_xx_yyy_0[j] = pa_x[j] * tdx_x_yyy_0[j] + 0.5 * fl1_fx * tdx_0_yyy_0[j] + 0.5 * fl1_fx * ts_x_yyy_0[j];

                tdy_xx_yyy_0[j] = pa_x[j] * tdy_x_yyy_0[j] + 0.5 * fl1_fx * tdy_0_yyy_0[j];

                tdz_xx_yyy_0[j] = pa_x[j] * tdz_x_yyy_0[j] + 0.5 * fl1_fx * tdz_0_yyy_0[j];

                tdx_xx_yyz_0[j] = pa_x[j] * tdx_x_yyz_0[j] + 0.5 * fl1_fx * tdx_0_yyz_0[j] + 0.5 * fl1_fx * ts_x_yyz_0[j];

                tdy_xx_yyz_0[j] = pa_x[j] * tdy_x_yyz_0[j] + 0.5 * fl1_fx * tdy_0_yyz_0[j];

                tdz_xx_yyz_0[j] = pa_x[j] * tdz_x_yyz_0[j] + 0.5 * fl1_fx * tdz_0_yyz_0[j];

                tdx_xx_yzz_0[j] = pa_x[j] * tdx_x_yzz_0[j] + 0.5 * fl1_fx * tdx_0_yzz_0[j] + 0.5 * fl1_fx * ts_x_yzz_0[j];

                tdy_xx_yzz_0[j] = pa_x[j] * tdy_x_yzz_0[j] + 0.5 * fl1_fx * tdy_0_yzz_0[j];

                tdz_xx_yzz_0[j] = pa_x[j] * tdz_x_yzz_0[j] + 0.5 * fl1_fx * tdz_0_yzz_0[j];

                tdx_xx_zzz_0[j] = pa_x[j] * tdx_x_zzz_0[j] + 0.5 * fl1_fx * tdx_0_zzz_0[j] + 0.5 * fl1_fx * ts_x_zzz_0[j];

                tdy_xx_zzz_0[j] = pa_x[j] * tdy_x_zzz_0[j] + 0.5 * fl1_fx * tdy_0_zzz_0[j];

                tdz_xx_zzz_0[j] = pa_x[j] * tdz_x_zzz_0[j] + 0.5 * fl1_fx * tdz_0_zzz_0[j];

                tdx_xy_xxx_0[j] = pa_x[j] * tdx_y_xxx_0[j] + 1.5 * fl1_fx * tdx_y_xx_0[j] + 0.5 * fl1_fx * ts_y_xxx_0[j];

                tdy_xy_xxx_0[j] = pa_x[j] * tdy_y_xxx_0[j] + 1.5 * fl1_fx * tdy_y_xx_0[j];

                tdz_xy_xxx_0[j] = pa_x[j] * tdz_y_xxx_0[j] + 1.5 * fl1_fx * tdz_y_xx_0[j];

                tdx_xy_xxy_0[j] = pa_x[j] * tdx_y_xxy_0[j] + fl1_fx * tdx_y_xy_0[j] + 0.5 * fl1_fx * ts_y_xxy_0[j];

                tdy_xy_xxy_0[j] = pa_x[j] * tdy_y_xxy_0[j] + fl1_fx * tdy_y_xy_0[j];

                tdz_xy_xxy_0[j] = pa_x[j] * tdz_y_xxy_0[j] + fl1_fx * tdz_y_xy_0[j];

                tdx_xy_xxz_0[j] = pa_x[j] * tdx_y_xxz_0[j] + fl1_fx * tdx_y_xz_0[j] + 0.5 * fl1_fx * ts_y_xxz_0[j];

                tdy_xy_xxz_0[j] = pa_x[j] * tdy_y_xxz_0[j] + fl1_fx * tdy_y_xz_0[j];

                tdz_xy_xxz_0[j] = pa_x[j] * tdz_y_xxz_0[j] + fl1_fx * tdz_y_xz_0[j];

                tdx_xy_xyy_0[j] = pa_x[j] * tdx_y_xyy_0[j] + 0.5 * fl1_fx * tdx_y_yy_0[j] + 0.5 * fl1_fx * ts_y_xyy_0[j];

                tdy_xy_xyy_0[j] = pa_x[j] * tdy_y_xyy_0[j] + 0.5 * fl1_fx * tdy_y_yy_0[j];

                tdz_xy_xyy_0[j] = pa_x[j] * tdz_y_xyy_0[j] + 0.5 * fl1_fx * tdz_y_yy_0[j];

                tdx_xy_xyz_0[j] = pa_x[j] * tdx_y_xyz_0[j] + 0.5 * fl1_fx * tdx_y_yz_0[j] + 0.5 * fl1_fx * ts_y_xyz_0[j];

                tdy_xy_xyz_0[j] = pa_x[j] * tdy_y_xyz_0[j] + 0.5 * fl1_fx * tdy_y_yz_0[j];

                tdz_xy_xyz_0[j] = pa_x[j] * tdz_y_xyz_0[j] + 0.5 * fl1_fx * tdz_y_yz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDF_45_90(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
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

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 15); 

            auto tdy_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tdz_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tdx_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 16); 

            auto tdy_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tdz_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tdx_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 17); 

            auto tdy_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tdz_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tdx_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 18); 

            auto tdy_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tdz_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tdx_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 19); 

            auto tdy_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tdz_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20); 

            auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21); 

            auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22); 

            auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23); 

            auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24); 

            auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24); 

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

            auto tdx_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 11); 

            auto tdy_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tdz_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 11); 

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

            auto ts_y_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 15); 

            auto ts_y_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 16); 

            auto ts_y_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 17); 

            auto ts_y_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 18); 

            auto ts_y_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 19); 

            auto ts_z_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 20); 

            auto ts_z_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 21); 

            auto ts_z_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 22); 

            auto ts_z_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 23); 

            auto ts_z_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 24); 

            auto ts_z_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 25); 

            auto ts_z_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 26); 

            auto ts_z_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 27); 

            auto ts_z_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 28); 

            auto ts_z_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 29); 

            // set up pointers to integrals

            auto tdx_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 15); 

            auto tdy_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 15); 

            auto tdz_xy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 15); 

            auto tdx_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 16); 

            auto tdy_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 16); 

            auto tdz_xy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 16); 

            auto tdx_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 17); 

            auto tdy_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 17); 

            auto tdz_xy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 17); 

            auto tdx_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 18); 

            auto tdy_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 18); 

            auto tdz_xy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 18); 

            auto tdx_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 19); 

            auto tdy_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 19); 

            auto tdz_xy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 19); 

            auto tdx_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 20); 

            auto tdy_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 20); 

            auto tdz_xz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 20); 

            auto tdx_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 21); 

            auto tdy_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 21); 

            auto tdz_xz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 21); 

            auto tdx_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 22); 

            auto tdy_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 22); 

            auto tdz_xz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 22); 

            auto tdx_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 23); 

            auto tdy_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 23); 

            auto tdz_xz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 23); 

            auto tdx_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 24); 

            auto tdy_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 24); 

            auto tdz_xz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 24); 

            auto tdx_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 25); 

            auto tdy_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 25); 

            auto tdz_xz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 25); 

            auto tdx_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 26); 

            auto tdy_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 26); 

            auto tdz_xz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 26); 

            auto tdx_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 27); 

            auto tdy_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 27); 

            auto tdz_xz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 27); 

            auto tdx_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 28); 

            auto tdy_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 28); 

            auto tdz_xz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 28); 

            auto tdx_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 29); 

            auto tdy_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 29); 

            auto tdz_xz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 29); 

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, tdx_xy_xzz_0, tdx_xy_yyy_0, tdx_xy_yyz_0, tdx_xy_yzz_0, \
                                     tdx_xy_zzz_0, tdx_xz_xxx_0, tdx_xz_xxy_0, tdx_xz_xxz_0, tdx_xz_xyy_0, tdx_xz_xyz_0, \
                                     tdx_xz_xzz_0, tdx_xz_yyy_0, tdx_xz_yyz_0, tdx_xz_yzz_0, tdx_xz_zzz_0, tdx_y_xzz_0, \
                                     tdx_y_yyy_0, tdx_y_yyz_0, tdx_y_yzz_0, tdx_y_zz_0, tdx_y_zzz_0, tdx_z_xx_0, \
                                     tdx_z_xxx_0, tdx_z_xxy_0, tdx_z_xxz_0, tdx_z_xy_0, tdx_z_xyy_0, tdx_z_xyz_0, \
                                     tdx_z_xz_0, tdx_z_xzz_0, tdx_z_yy_0, tdx_z_yyy_0, tdx_z_yyz_0, tdx_z_yz_0, \
                                     tdx_z_yzz_0, tdx_z_zz_0, tdx_z_zzz_0, tdy_xy_xzz_0, tdy_xy_yyy_0, tdy_xy_yyz_0, \
                                     tdy_xy_yzz_0, tdy_xy_zzz_0, tdy_xz_xxx_0, tdy_xz_xxy_0, tdy_xz_xxz_0, tdy_xz_xyy_0, \
                                     tdy_xz_xyz_0, tdy_xz_xzz_0, tdy_xz_yyy_0, tdy_xz_yyz_0, tdy_xz_yzz_0, tdy_xz_zzz_0, \
                                     tdy_y_xzz_0, tdy_y_yyy_0, tdy_y_yyz_0, tdy_y_yzz_0, tdy_y_zz_0, tdy_y_zzz_0, \
                                     tdy_z_xx_0, tdy_z_xxx_0, tdy_z_xxy_0, tdy_z_xxz_0, tdy_z_xy_0, tdy_z_xyy_0, \
                                     tdy_z_xyz_0, tdy_z_xz_0, tdy_z_xzz_0, tdy_z_yy_0, tdy_z_yyy_0, tdy_z_yyz_0, \
                                     tdy_z_yz_0, tdy_z_yzz_0, tdy_z_zz_0, tdy_z_zzz_0, tdz_xy_xzz_0, tdz_xy_yyy_0, \
                                     tdz_xy_yyz_0, tdz_xy_yzz_0, tdz_xy_zzz_0, tdz_xz_xxx_0, tdz_xz_xxy_0, tdz_xz_xxz_0, \
                                     tdz_xz_xyy_0, tdz_xz_xyz_0, tdz_xz_xzz_0, tdz_xz_yyy_0, tdz_xz_yyz_0, tdz_xz_yzz_0, \
                                     tdz_xz_zzz_0, tdz_y_xzz_0, tdz_y_yyy_0, tdz_y_yyz_0, tdz_y_yzz_0, tdz_y_zz_0, \
                                     tdz_y_zzz_0, tdz_z_xx_0, tdz_z_xxx_0, tdz_z_xxy_0, tdz_z_xxz_0, tdz_z_xy_0, \
                                     tdz_z_xyy_0, tdz_z_xyz_0, tdz_z_xz_0, tdz_z_xzz_0, tdz_z_yy_0, tdz_z_yyy_0, \
                                     tdz_z_yyz_0, tdz_z_yz_0, tdz_z_yzz_0, tdz_z_zz_0, tdz_z_zzz_0, ts_y_xzz_0, \
                                     ts_y_yyy_0, ts_y_yyz_0, ts_y_yzz_0, ts_y_zzz_0, ts_z_xxx_0, ts_z_xxy_0, ts_z_xxz_0, \
                                     ts_z_xyy_0, ts_z_xyz_0, ts_z_xzz_0, ts_z_yyy_0, ts_z_yyz_0, ts_z_yzz_0, ts_z_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xy_xzz_0[j] = pa_x[j] * tdx_y_xzz_0[j] + 0.5 * fl1_fx * tdx_y_zz_0[j] + 0.5 * fl1_fx * ts_y_xzz_0[j];

                tdy_xy_xzz_0[j] = pa_x[j] * tdy_y_xzz_0[j] + 0.5 * fl1_fx * tdy_y_zz_0[j];

                tdz_xy_xzz_0[j] = pa_x[j] * tdz_y_xzz_0[j] + 0.5 * fl1_fx * tdz_y_zz_0[j];

                tdx_xy_yyy_0[j] = pa_x[j] * tdx_y_yyy_0[j] + 0.5 * fl1_fx * ts_y_yyy_0[j];

                tdy_xy_yyy_0[j] = pa_x[j] * tdy_y_yyy_0[j];

                tdz_xy_yyy_0[j] = pa_x[j] * tdz_y_yyy_0[j];

                tdx_xy_yyz_0[j] = pa_x[j] * tdx_y_yyz_0[j] + 0.5 * fl1_fx * ts_y_yyz_0[j];

                tdy_xy_yyz_0[j] = pa_x[j] * tdy_y_yyz_0[j];

                tdz_xy_yyz_0[j] = pa_x[j] * tdz_y_yyz_0[j];

                tdx_xy_yzz_0[j] = pa_x[j] * tdx_y_yzz_0[j] + 0.5 * fl1_fx * ts_y_yzz_0[j];

                tdy_xy_yzz_0[j] = pa_x[j] * tdy_y_yzz_0[j];

                tdz_xy_yzz_0[j] = pa_x[j] * tdz_y_yzz_0[j];

                tdx_xy_zzz_0[j] = pa_x[j] * tdx_y_zzz_0[j] + 0.5 * fl1_fx * ts_y_zzz_0[j];

                tdy_xy_zzz_0[j] = pa_x[j] * tdy_y_zzz_0[j];

                tdz_xy_zzz_0[j] = pa_x[j] * tdz_y_zzz_0[j];

                tdx_xz_xxx_0[j] = pa_x[j] * tdx_z_xxx_0[j] + 1.5 * fl1_fx * tdx_z_xx_0[j] + 0.5 * fl1_fx * ts_z_xxx_0[j];

                tdy_xz_xxx_0[j] = pa_x[j] * tdy_z_xxx_0[j] + 1.5 * fl1_fx * tdy_z_xx_0[j];

                tdz_xz_xxx_0[j] = pa_x[j] * tdz_z_xxx_0[j] + 1.5 * fl1_fx * tdz_z_xx_0[j];

                tdx_xz_xxy_0[j] = pa_x[j] * tdx_z_xxy_0[j] + fl1_fx * tdx_z_xy_0[j] + 0.5 * fl1_fx * ts_z_xxy_0[j];

                tdy_xz_xxy_0[j] = pa_x[j] * tdy_z_xxy_0[j] + fl1_fx * tdy_z_xy_0[j];

                tdz_xz_xxy_0[j] = pa_x[j] * tdz_z_xxy_0[j] + fl1_fx * tdz_z_xy_0[j];

                tdx_xz_xxz_0[j] = pa_x[j] * tdx_z_xxz_0[j] + fl1_fx * tdx_z_xz_0[j] + 0.5 * fl1_fx * ts_z_xxz_0[j];

                tdy_xz_xxz_0[j] = pa_x[j] * tdy_z_xxz_0[j] + fl1_fx * tdy_z_xz_0[j];

                tdz_xz_xxz_0[j] = pa_x[j] * tdz_z_xxz_0[j] + fl1_fx * tdz_z_xz_0[j];

                tdx_xz_xyy_0[j] = pa_x[j] * tdx_z_xyy_0[j] + 0.5 * fl1_fx * tdx_z_yy_0[j] + 0.5 * fl1_fx * ts_z_xyy_0[j];

                tdy_xz_xyy_0[j] = pa_x[j] * tdy_z_xyy_0[j] + 0.5 * fl1_fx * tdy_z_yy_0[j];

                tdz_xz_xyy_0[j] = pa_x[j] * tdz_z_xyy_0[j] + 0.5 * fl1_fx * tdz_z_yy_0[j];

                tdx_xz_xyz_0[j] = pa_x[j] * tdx_z_xyz_0[j] + 0.5 * fl1_fx * tdx_z_yz_0[j] + 0.5 * fl1_fx * ts_z_xyz_0[j];

                tdy_xz_xyz_0[j] = pa_x[j] * tdy_z_xyz_0[j] + 0.5 * fl1_fx * tdy_z_yz_0[j];

                tdz_xz_xyz_0[j] = pa_x[j] * tdz_z_xyz_0[j] + 0.5 * fl1_fx * tdz_z_yz_0[j];

                tdx_xz_xzz_0[j] = pa_x[j] * tdx_z_xzz_0[j] + 0.5 * fl1_fx * tdx_z_zz_0[j] + 0.5 * fl1_fx * ts_z_xzz_0[j];

                tdy_xz_xzz_0[j] = pa_x[j] * tdy_z_xzz_0[j] + 0.5 * fl1_fx * tdy_z_zz_0[j];

                tdz_xz_xzz_0[j] = pa_x[j] * tdz_z_xzz_0[j] + 0.5 * fl1_fx * tdz_z_zz_0[j];

                tdx_xz_yyy_0[j] = pa_x[j] * tdx_z_yyy_0[j] + 0.5 * fl1_fx * ts_z_yyy_0[j];

                tdy_xz_yyy_0[j] = pa_x[j] * tdy_z_yyy_0[j];

                tdz_xz_yyy_0[j] = pa_x[j] * tdz_z_yyy_0[j];

                tdx_xz_yyz_0[j] = pa_x[j] * tdx_z_yyz_0[j] + 0.5 * fl1_fx * ts_z_yyz_0[j];

                tdy_xz_yyz_0[j] = pa_x[j] * tdy_z_yyz_0[j];

                tdz_xz_yyz_0[j] = pa_x[j] * tdz_z_yyz_0[j];

                tdx_xz_yzz_0[j] = pa_x[j] * tdx_z_yzz_0[j] + 0.5 * fl1_fx * ts_z_yzz_0[j];

                tdy_xz_yzz_0[j] = pa_x[j] * tdy_z_yzz_0[j];

                tdz_xz_yzz_0[j] = pa_x[j] * tdz_z_yzz_0[j];

                tdx_xz_zzz_0[j] = pa_x[j] * tdx_z_zzz_0[j] + 0.5 * fl1_fx * ts_z_zzz_0[j];

                tdy_xz_zzz_0[j] = pa_x[j] * tdy_z_zzz_0[j];

                tdz_xz_zzz_0[j] = pa_x[j] * tdz_z_zzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDF_90_135(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
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

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tdx_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 10); 

            auto tdy_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tdz_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tdx_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 11); 

            auto tdy_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tdz_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tdx_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 12); 

            auto tdy_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tdz_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tdx_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 13); 

            auto tdy_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tdz_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tdx_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 14); 

            auto tdy_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tdz_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tdx_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 15); 

            auto tdy_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tdz_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tdx_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 16); 

            auto tdy_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tdz_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tdx_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 17); 

            auto tdy_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tdz_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tdx_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 18); 

            auto tdy_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tdz_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tdx_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 19); 

            auto tdy_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tdz_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20); 

            auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21); 

            auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22); 

            auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23); 

            auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24); 

            auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tdx_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx); 

            auto tdy_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx); 

            auto tdz_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx); 

            auto tdx_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 1); 

            auto tdy_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 1); 

            auto tdz_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 1); 

            auto tdx_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 2); 

            auto tdy_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 2); 

            auto tdz_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 2); 

            auto tdx_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 3); 

            auto tdy_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 3); 

            auto tdz_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 3); 

            auto tdx_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 4); 

            auto tdy_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 4); 

            auto tdz_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 4); 

            auto tdx_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 5); 

            auto tdy_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 5); 

            auto tdz_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 5); 

            auto tdx_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 6); 

            auto tdy_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 6); 

            auto tdz_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 6); 

            auto tdx_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 7); 

            auto tdy_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 7); 

            auto tdz_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 7); 

            auto tdx_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 8); 

            auto tdy_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 8); 

            auto tdz_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 8); 

            auto tdx_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 9); 

            auto tdy_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 9); 

            auto tdz_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 9); 

            auto tdx_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 6); 

            auto tdy_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tdz_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tdx_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 7); 

            auto tdy_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tdz_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tdx_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 8); 

            auto tdy_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tdz_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tdx_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 9); 

            auto tdy_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tdz_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdx_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 10); 

            auto tdy_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tdz_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tdx_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 11); 

            auto tdy_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tdz_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tdx_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 12); 

            auto tdy_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tdz_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tdx_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 13); 

            auto tdy_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tdz_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tdx_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 14); 

            auto tdy_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tdz_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto ts_y_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 10); 

            auto ts_y_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 11); 

            auto ts_y_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 12); 

            auto ts_y_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 13); 

            auto ts_y_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 14); 

            auto ts_y_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 15); 

            auto ts_y_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 16); 

            auto ts_y_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 17); 

            auto ts_y_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 18); 

            auto ts_y_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 19); 

            auto ts_z_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 20); 

            auto ts_z_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 21); 

            auto ts_z_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 22); 

            auto ts_z_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 23); 

            auto ts_z_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 24); 

            // set up pointers to integrals

            auto tdx_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 30); 

            auto tdy_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 30); 

            auto tdz_yy_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 30); 

            auto tdx_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 31); 

            auto tdy_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 31); 

            auto tdz_yy_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 31); 

            auto tdx_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 32); 

            auto tdy_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 32); 

            auto tdz_yy_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 32); 

            auto tdx_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 33); 

            auto tdy_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 33); 

            auto tdz_yy_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 33); 

            auto tdx_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 34); 

            auto tdy_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 34); 

            auto tdz_yy_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 34); 

            auto tdx_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 35); 

            auto tdy_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 35); 

            auto tdz_yy_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 35); 

            auto tdx_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 36); 

            auto tdy_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 36); 

            auto tdz_yy_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdx_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 37); 

            auto tdy_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 37); 

            auto tdz_yy_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdx_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 38); 

            auto tdy_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 38); 

            auto tdz_yy_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdx_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 39); 

            auto tdy_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 39); 

            auto tdz_yy_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdx_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 40); 

            auto tdy_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 40); 

            auto tdz_yz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdx_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 41); 

            auto tdy_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 41); 

            auto tdz_yz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 41); 

            auto tdx_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 42); 

            auto tdy_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 42); 

            auto tdz_yz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdx_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 43); 

            auto tdy_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 43); 

            auto tdz_yz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdx_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 44); 

            auto tdy_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 44); 

            auto tdz_yz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 44); 

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fx, pa_y, tdx_0_xxx_0, tdx_0_xxy_0, tdx_0_xxz_0, tdx_0_xyy_0, \
                                     tdx_0_xyz_0, tdx_0_xzz_0, tdx_0_yyy_0, tdx_0_yyz_0, tdx_0_yzz_0, tdx_0_zzz_0, \
                                     tdx_y_xx_0, tdx_y_xxx_0, tdx_y_xxy_0, tdx_y_xxz_0, tdx_y_xy_0, tdx_y_xyy_0, \
                                     tdx_y_xyz_0, tdx_y_xz_0, tdx_y_xzz_0, tdx_y_yy_0, tdx_y_yyy_0, tdx_y_yyz_0, \
                                     tdx_y_yz_0, tdx_y_yzz_0, tdx_y_zz_0, tdx_y_zzz_0, tdx_yy_xxx_0, tdx_yy_xxy_0, \
                                     tdx_yy_xxz_0, tdx_yy_xyy_0, tdx_yy_xyz_0, tdx_yy_xzz_0, tdx_yy_yyy_0, tdx_yy_yyz_0, \
                                     tdx_yy_yzz_0, tdx_yy_zzz_0, tdx_yz_xxx_0, tdx_yz_xxy_0, tdx_yz_xxz_0, tdx_yz_xyy_0, \
                                     tdx_yz_xyz_0, tdx_z_xx_0, tdx_z_xxx_0, tdx_z_xxy_0, tdx_z_xxz_0, tdx_z_xy_0, \
                                     tdx_z_xyy_0, tdx_z_xyz_0, tdx_z_xz_0, tdy_0_xxx_0, tdy_0_xxy_0, tdy_0_xxz_0, \
                                     tdy_0_xyy_0, tdy_0_xyz_0, tdy_0_xzz_0, tdy_0_yyy_0, tdy_0_yyz_0, tdy_0_yzz_0, \
                                     tdy_0_zzz_0, tdy_y_xx_0, tdy_y_xxx_0, tdy_y_xxy_0, tdy_y_xxz_0, tdy_y_xy_0, \
                                     tdy_y_xyy_0, tdy_y_xyz_0, tdy_y_xz_0, tdy_y_xzz_0, tdy_y_yy_0, tdy_y_yyy_0, \
                                     tdy_y_yyz_0, tdy_y_yz_0, tdy_y_yzz_0, tdy_y_zz_0, tdy_y_zzz_0, tdy_yy_xxx_0, \
                                     tdy_yy_xxy_0, tdy_yy_xxz_0, tdy_yy_xyy_0, tdy_yy_xyz_0, tdy_yy_xzz_0, tdy_yy_yyy_0, \
                                     tdy_yy_yyz_0, tdy_yy_yzz_0, tdy_yy_zzz_0, tdy_yz_xxx_0, tdy_yz_xxy_0, tdy_yz_xxz_0, \
                                     tdy_yz_xyy_0, tdy_yz_xyz_0, tdy_z_xx_0, tdy_z_xxx_0, tdy_z_xxy_0, tdy_z_xxz_0, \
                                     tdy_z_xy_0, tdy_z_xyy_0, tdy_z_xyz_0, tdy_z_xz_0, tdz_0_xxx_0, tdz_0_xxy_0, \
                                     tdz_0_xxz_0, tdz_0_xyy_0, tdz_0_xyz_0, tdz_0_xzz_0, tdz_0_yyy_0, tdz_0_yyz_0, \
                                     tdz_0_yzz_0, tdz_0_zzz_0, tdz_y_xx_0, tdz_y_xxx_0, tdz_y_xxy_0, tdz_y_xxz_0, \
                                     tdz_y_xy_0, tdz_y_xyy_0, tdz_y_xyz_0, tdz_y_xz_0, tdz_y_xzz_0, tdz_y_yy_0, \
                                     tdz_y_yyy_0, tdz_y_yyz_0, tdz_y_yz_0, tdz_y_yzz_0, tdz_y_zz_0, tdz_y_zzz_0, \
                                     tdz_yy_xxx_0, tdz_yy_xxy_0, tdz_yy_xxz_0, tdz_yy_xyy_0, tdz_yy_xyz_0, tdz_yy_xzz_0, \
                                     tdz_yy_yyy_0, tdz_yy_yyz_0, tdz_yy_yzz_0, tdz_yy_zzz_0, tdz_yz_xxx_0, tdz_yz_xxy_0, \
                                     tdz_yz_xxz_0, tdz_yz_xyy_0, tdz_yz_xyz_0, tdz_z_xx_0, tdz_z_xxx_0, tdz_z_xxy_0, \
                                     tdz_z_xxz_0, tdz_z_xy_0, tdz_z_xyy_0, tdz_z_xyz_0, tdz_z_xz_0, ts_y_xxx_0, \
                                     ts_y_xxy_0, ts_y_xxz_0, ts_y_xyy_0, ts_y_xyz_0, ts_y_xzz_0, ts_y_yyy_0, ts_y_yyz_0, \
                                     ts_y_yzz_0, ts_y_zzz_0, ts_z_xxx_0, ts_z_xxy_0, ts_z_xxz_0, ts_z_xyy_0, ts_z_xyz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yy_xxx_0[j] = pa_y[j] * tdx_y_xxx_0[j] + 0.5 * fl1_fx * tdx_0_xxx_0[j];

                tdy_yy_xxx_0[j] = pa_y[j] * tdy_y_xxx_0[j] + 0.5 * fl1_fx * tdy_0_xxx_0[j] + 0.5 * fl1_fx * ts_y_xxx_0[j];

                tdz_yy_xxx_0[j] = pa_y[j] * tdz_y_xxx_0[j] + 0.5 * fl1_fx * tdz_0_xxx_0[j];

                tdx_yy_xxy_0[j] = pa_y[j] * tdx_y_xxy_0[j] + 0.5 * fl1_fx * tdx_0_xxy_0[j] + 0.5 * fl1_fx * tdx_y_xx_0[j];

                tdy_yy_xxy_0[j] = pa_y[j] * tdy_y_xxy_0[j] + 0.5 * fl1_fx * tdy_0_xxy_0[j] + 0.5 * fl1_fx * tdy_y_xx_0[j] + 0.5 * fl1_fx * ts_y_xxy_0[j];

                tdz_yy_xxy_0[j] = pa_y[j] * tdz_y_xxy_0[j] + 0.5 * fl1_fx * tdz_0_xxy_0[j] + 0.5 * fl1_fx * tdz_y_xx_0[j];

                tdx_yy_xxz_0[j] = pa_y[j] * tdx_y_xxz_0[j] + 0.5 * fl1_fx * tdx_0_xxz_0[j];

                tdy_yy_xxz_0[j] = pa_y[j] * tdy_y_xxz_0[j] + 0.5 * fl1_fx * tdy_0_xxz_0[j] + 0.5 * fl1_fx * ts_y_xxz_0[j];

                tdz_yy_xxz_0[j] = pa_y[j] * tdz_y_xxz_0[j] + 0.5 * fl1_fx * tdz_0_xxz_0[j];

                tdx_yy_xyy_0[j] = pa_y[j] * tdx_y_xyy_0[j] + 0.5 * fl1_fx * tdx_0_xyy_0[j] + fl1_fx * tdx_y_xy_0[j];

                tdy_yy_xyy_0[j] = pa_y[j] * tdy_y_xyy_0[j] + 0.5 * fl1_fx * tdy_0_xyy_0[j] + fl1_fx * tdy_y_xy_0[j] + 0.5 * fl1_fx * ts_y_xyy_0[j];

                tdz_yy_xyy_0[j] = pa_y[j] * tdz_y_xyy_0[j] + 0.5 * fl1_fx * tdz_0_xyy_0[j] + fl1_fx * tdz_y_xy_0[j];

                tdx_yy_xyz_0[j] = pa_y[j] * tdx_y_xyz_0[j] + 0.5 * fl1_fx * tdx_0_xyz_0[j] + 0.5 * fl1_fx * tdx_y_xz_0[j];

                tdy_yy_xyz_0[j] = pa_y[j] * tdy_y_xyz_0[j] + 0.5 * fl1_fx * tdy_0_xyz_0[j] + 0.5 * fl1_fx * tdy_y_xz_0[j] + 0.5 * fl1_fx * ts_y_xyz_0[j];

                tdz_yy_xyz_0[j] = pa_y[j] * tdz_y_xyz_0[j] + 0.5 * fl1_fx * tdz_0_xyz_0[j] + 0.5 * fl1_fx * tdz_y_xz_0[j];

                tdx_yy_xzz_0[j] = pa_y[j] * tdx_y_xzz_0[j] + 0.5 * fl1_fx * tdx_0_xzz_0[j];

                tdy_yy_xzz_0[j] = pa_y[j] * tdy_y_xzz_0[j] + 0.5 * fl1_fx * tdy_0_xzz_0[j] + 0.5 * fl1_fx * ts_y_xzz_0[j];

                tdz_yy_xzz_0[j] = pa_y[j] * tdz_y_xzz_0[j] + 0.5 * fl1_fx * tdz_0_xzz_0[j];

                tdx_yy_yyy_0[j] = pa_y[j] * tdx_y_yyy_0[j] + 0.5 * fl1_fx * tdx_0_yyy_0[j] + 1.5 * fl1_fx * tdx_y_yy_0[j];

                tdy_yy_yyy_0[j] = pa_y[j] * tdy_y_yyy_0[j] + 0.5 * fl1_fx * tdy_0_yyy_0[j] + 1.5 * fl1_fx * tdy_y_yy_0[j] + 0.5 * fl1_fx * ts_y_yyy_0[j];

                tdz_yy_yyy_0[j] = pa_y[j] * tdz_y_yyy_0[j] + 0.5 * fl1_fx * tdz_0_yyy_0[j] + 1.5 * fl1_fx * tdz_y_yy_0[j];

                tdx_yy_yyz_0[j] = pa_y[j] * tdx_y_yyz_0[j] + 0.5 * fl1_fx * tdx_0_yyz_0[j] + fl1_fx * tdx_y_yz_0[j];

                tdy_yy_yyz_0[j] = pa_y[j] * tdy_y_yyz_0[j] + 0.5 * fl1_fx * tdy_0_yyz_0[j] + fl1_fx * tdy_y_yz_0[j] + 0.5 * fl1_fx * ts_y_yyz_0[j];

                tdz_yy_yyz_0[j] = pa_y[j] * tdz_y_yyz_0[j] + 0.5 * fl1_fx * tdz_0_yyz_0[j] + fl1_fx * tdz_y_yz_0[j];

                tdx_yy_yzz_0[j] = pa_y[j] * tdx_y_yzz_0[j] + 0.5 * fl1_fx * tdx_0_yzz_0[j] + 0.5 * fl1_fx * tdx_y_zz_0[j];

                tdy_yy_yzz_0[j] = pa_y[j] * tdy_y_yzz_0[j] + 0.5 * fl1_fx * tdy_0_yzz_0[j] + 0.5 * fl1_fx * tdy_y_zz_0[j] + 0.5 * fl1_fx * ts_y_yzz_0[j];

                tdz_yy_yzz_0[j] = pa_y[j] * tdz_y_yzz_0[j] + 0.5 * fl1_fx * tdz_0_yzz_0[j] + 0.5 * fl1_fx * tdz_y_zz_0[j];

                tdx_yy_zzz_0[j] = pa_y[j] * tdx_y_zzz_0[j] + 0.5 * fl1_fx * tdx_0_zzz_0[j];

                tdy_yy_zzz_0[j] = pa_y[j] * tdy_y_zzz_0[j] + 0.5 * fl1_fx * tdy_0_zzz_0[j] + 0.5 * fl1_fx * ts_y_zzz_0[j];

                tdz_yy_zzz_0[j] = pa_y[j] * tdz_y_zzz_0[j] + 0.5 * fl1_fx * tdz_0_zzz_0[j];

                tdx_yz_xxx_0[j] = pa_y[j] * tdx_z_xxx_0[j];

                tdy_yz_xxx_0[j] = pa_y[j] * tdy_z_xxx_0[j] + 0.5 * fl1_fx * ts_z_xxx_0[j];

                tdz_yz_xxx_0[j] = pa_y[j] * tdz_z_xxx_0[j];

                tdx_yz_xxy_0[j] = pa_y[j] * tdx_z_xxy_0[j] + 0.5 * fl1_fx * tdx_z_xx_0[j];

                tdy_yz_xxy_0[j] = pa_y[j] * tdy_z_xxy_0[j] + 0.5 * fl1_fx * tdy_z_xx_0[j] + 0.5 * fl1_fx * ts_z_xxy_0[j];

                tdz_yz_xxy_0[j] = pa_y[j] * tdz_z_xxy_0[j] + 0.5 * fl1_fx * tdz_z_xx_0[j];

                tdx_yz_xxz_0[j] = pa_y[j] * tdx_z_xxz_0[j];

                tdy_yz_xxz_0[j] = pa_y[j] * tdy_z_xxz_0[j] + 0.5 * fl1_fx * ts_z_xxz_0[j];

                tdz_yz_xxz_0[j] = pa_y[j] * tdz_z_xxz_0[j];

                tdx_yz_xyy_0[j] = pa_y[j] * tdx_z_xyy_0[j] + fl1_fx * tdx_z_xy_0[j];

                tdy_yz_xyy_0[j] = pa_y[j] * tdy_z_xyy_0[j] + fl1_fx * tdy_z_xy_0[j] + 0.5 * fl1_fx * ts_z_xyy_0[j];

                tdz_yz_xyy_0[j] = pa_y[j] * tdz_z_xyy_0[j] + fl1_fx * tdz_z_xy_0[j];

                tdx_yz_xyz_0[j] = pa_y[j] * tdx_z_xyz_0[j] + 0.5 * fl1_fx * tdx_z_xz_0[j];

                tdy_yz_xyz_0[j] = pa_y[j] * tdy_z_xyz_0[j] + 0.5 * fl1_fx * tdy_z_xz_0[j] + 0.5 * fl1_fx * ts_z_xyz_0[j];

                tdz_yz_xyz_0[j] = pa_y[j] * tdz_z_xyz_0[j] + 0.5 * fl1_fx * tdz_z_xz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDF_135_180(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
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

        auto pidx_d_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20); 

            auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21); 

            auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22); 

            auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23); 

            auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24); 

            auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24); 

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

            auto tdx_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx); 

            auto tdy_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx); 

            auto tdz_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx); 

            auto tdx_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 1); 

            auto tdy_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 1); 

            auto tdz_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 1); 

            auto tdx_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 2); 

            auto tdy_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 2); 

            auto tdz_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 2); 

            auto tdx_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 3); 

            auto tdy_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 3); 

            auto tdz_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 3); 

            auto tdx_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 4); 

            auto tdy_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 4); 

            auto tdz_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 4); 

            auto tdx_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 5); 

            auto tdy_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 5); 

            auto tdz_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 5); 

            auto tdx_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 6); 

            auto tdy_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 6); 

            auto tdz_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 6); 

            auto tdx_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 7); 

            auto tdy_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 7); 

            auto tdz_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 7); 

            auto tdx_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 8); 

            auto tdy_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 8); 

            auto tdz_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 8); 

            auto tdx_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 9); 

            auto tdy_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 9); 

            auto tdz_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 9); 

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

            auto ts_z_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 20); 

            auto ts_z_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 21); 

            auto ts_z_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 22); 

            auto ts_z_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 23); 

            auto ts_z_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 24); 

            auto ts_z_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 25); 

            auto ts_z_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 26); 

            auto ts_z_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 27); 

            auto ts_z_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 28); 

            auto ts_z_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 29); 

            // set up pointers to integrals

            auto tdx_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 45); 

            auto tdy_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 45); 

            auto tdz_yz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdx_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 46); 

            auto tdy_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 46); 

            auto tdz_yz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdx_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 47); 

            auto tdy_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 47); 

            auto tdz_yz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdx_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 48); 

            auto tdy_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 48); 

            auto tdz_yz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdx_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 49); 

            auto tdy_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 49); 

            auto tdz_yz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdx_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 50); 

            auto tdy_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_zz_xxx_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 51); 

            auto tdy_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_zz_xxy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 52); 

            auto tdy_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_zz_xxz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 53); 

            auto tdy_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 53); 

            auto tdz_zz_xyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 53); 

            auto tdx_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 54); 

            auto tdy_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 54); 

            auto tdz_zz_xyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 54); 

            auto tdx_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 55); 

            auto tdy_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 55); 

            auto tdz_zz_xzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 55); 

            auto tdx_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 56); 

            auto tdy_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 56); 

            auto tdz_zz_yyy_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 56); 

            auto tdx_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 57); 

            auto tdy_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 57); 

            auto tdz_zz_yyz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 57); 

            auto tdx_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 58); 

            auto tdy_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 58); 

            auto tdz_zz_yzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 58); 

            auto tdx_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * idx + 59); 

            auto tdy_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 60 * bdim + 60 * idx + 59); 

            auto tdz_zz_zzz_0 = primBuffer.data(pidx_d_2_3_m0 + 120 * bdim + 60 * idx + 59); 

            // Batch of Integrals (135,180)

            #pragma omp simd aligned(fx, pa_y, pa_z, tdx_0_xxx_0, tdx_0_xxy_0, tdx_0_xxz_0, tdx_0_xyy_0, \
                                     tdx_0_xyz_0, tdx_0_xzz_0, tdx_0_yyy_0, tdx_0_yyz_0, tdx_0_yzz_0, tdx_0_zzz_0, \
                                     tdx_yz_xzz_0, tdx_yz_yyy_0, tdx_yz_yyz_0, tdx_yz_yzz_0, tdx_yz_zzz_0, tdx_z_xx_0, \
                                     tdx_z_xxx_0, tdx_z_xxy_0, tdx_z_xxz_0, tdx_z_xy_0, tdx_z_xyy_0, tdx_z_xyz_0, \
                                     tdx_z_xz_0, tdx_z_xzz_0, tdx_z_yy_0, tdx_z_yyy_0, tdx_z_yyz_0, tdx_z_yz_0, \
                                     tdx_z_yzz_0, tdx_z_zz_0, tdx_z_zzz_0, tdx_zz_xxx_0, tdx_zz_xxy_0, tdx_zz_xxz_0, \
                                     tdx_zz_xyy_0, tdx_zz_xyz_0, tdx_zz_xzz_0, tdx_zz_yyy_0, tdx_zz_yyz_0, tdx_zz_yzz_0, \
                                     tdx_zz_zzz_0, tdy_0_xxx_0, tdy_0_xxy_0, tdy_0_xxz_0, tdy_0_xyy_0, tdy_0_xyz_0, \
                                     tdy_0_xzz_0, tdy_0_yyy_0, tdy_0_yyz_0, tdy_0_yzz_0, tdy_0_zzz_0, tdy_yz_xzz_0, \
                                     tdy_yz_yyy_0, tdy_yz_yyz_0, tdy_yz_yzz_0, tdy_yz_zzz_0, tdy_z_xx_0, tdy_z_xxx_0, \
                                     tdy_z_xxy_0, tdy_z_xxz_0, tdy_z_xy_0, tdy_z_xyy_0, tdy_z_xyz_0, tdy_z_xz_0, \
                                     tdy_z_xzz_0, tdy_z_yy_0, tdy_z_yyy_0, tdy_z_yyz_0, tdy_z_yz_0, tdy_z_yzz_0, \
                                     tdy_z_zz_0, tdy_z_zzz_0, tdy_zz_xxx_0, tdy_zz_xxy_0, tdy_zz_xxz_0, tdy_zz_xyy_0, \
                                     tdy_zz_xyz_0, tdy_zz_xzz_0, tdy_zz_yyy_0, tdy_zz_yyz_0, tdy_zz_yzz_0, tdy_zz_zzz_0, \
                                     tdz_0_xxx_0, tdz_0_xxy_0, tdz_0_xxz_0, tdz_0_xyy_0, tdz_0_xyz_0, tdz_0_xzz_0, \
                                     tdz_0_yyy_0, tdz_0_yyz_0, tdz_0_yzz_0, tdz_0_zzz_0, tdz_yz_xzz_0, tdz_yz_yyy_0, \
                                     tdz_yz_yyz_0, tdz_yz_yzz_0, tdz_yz_zzz_0, tdz_z_xx_0, tdz_z_xxx_0, tdz_z_xxy_0, \
                                     tdz_z_xxz_0, tdz_z_xy_0, tdz_z_xyy_0, tdz_z_xyz_0, tdz_z_xz_0, tdz_z_xzz_0, \
                                     tdz_z_yy_0, tdz_z_yyy_0, tdz_z_yyz_0, tdz_z_yz_0, tdz_z_yzz_0, tdz_z_zz_0, \
                                     tdz_z_zzz_0, tdz_zz_xxx_0, tdz_zz_xxy_0, tdz_zz_xxz_0, tdz_zz_xyy_0, tdz_zz_xyz_0, \
                                     tdz_zz_xzz_0, tdz_zz_yyy_0, tdz_zz_yyz_0, tdz_zz_yzz_0, tdz_zz_zzz_0, ts_z_xxx_0, \
                                     ts_z_xxy_0, ts_z_xxz_0, ts_z_xyy_0, ts_z_xyz_0, ts_z_xzz_0, ts_z_yyy_0, ts_z_yyz_0, \
                                     ts_z_yzz_0, ts_z_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yz_xzz_0[j] = pa_y[j] * tdx_z_xzz_0[j];

                tdy_yz_xzz_0[j] = pa_y[j] * tdy_z_xzz_0[j] + 0.5 * fl1_fx * ts_z_xzz_0[j];

                tdz_yz_xzz_0[j] = pa_y[j] * tdz_z_xzz_0[j];

                tdx_yz_yyy_0[j] = pa_y[j] * tdx_z_yyy_0[j] + 1.5 * fl1_fx * tdx_z_yy_0[j];

                tdy_yz_yyy_0[j] = pa_y[j] * tdy_z_yyy_0[j] + 1.5 * fl1_fx * tdy_z_yy_0[j] + 0.5 * fl1_fx * ts_z_yyy_0[j];

                tdz_yz_yyy_0[j] = pa_y[j] * tdz_z_yyy_0[j] + 1.5 * fl1_fx * tdz_z_yy_0[j];

                tdx_yz_yyz_0[j] = pa_y[j] * tdx_z_yyz_0[j] + fl1_fx * tdx_z_yz_0[j];

                tdy_yz_yyz_0[j] = pa_y[j] * tdy_z_yyz_0[j] + fl1_fx * tdy_z_yz_0[j] + 0.5 * fl1_fx * ts_z_yyz_0[j];

                tdz_yz_yyz_0[j] = pa_y[j] * tdz_z_yyz_0[j] + fl1_fx * tdz_z_yz_0[j];

                tdx_yz_yzz_0[j] = pa_y[j] * tdx_z_yzz_0[j] + 0.5 * fl1_fx * tdx_z_zz_0[j];

                tdy_yz_yzz_0[j] = pa_y[j] * tdy_z_yzz_0[j] + 0.5 * fl1_fx * tdy_z_zz_0[j] + 0.5 * fl1_fx * ts_z_yzz_0[j];

                tdz_yz_yzz_0[j] = pa_y[j] * tdz_z_yzz_0[j] + 0.5 * fl1_fx * tdz_z_zz_0[j];

                tdx_yz_zzz_0[j] = pa_y[j] * tdx_z_zzz_0[j];

                tdy_yz_zzz_0[j] = pa_y[j] * tdy_z_zzz_0[j] + 0.5 * fl1_fx * ts_z_zzz_0[j];

                tdz_yz_zzz_0[j] = pa_y[j] * tdz_z_zzz_0[j];

                tdx_zz_xxx_0[j] = pa_z[j] * tdx_z_xxx_0[j] + 0.5 * fl1_fx * tdx_0_xxx_0[j];

                tdy_zz_xxx_0[j] = pa_z[j] * tdy_z_xxx_0[j] + 0.5 * fl1_fx * tdy_0_xxx_0[j];

                tdz_zz_xxx_0[j] = pa_z[j] * tdz_z_xxx_0[j] + 0.5 * fl1_fx * tdz_0_xxx_0[j] + 0.5 * fl1_fx * ts_z_xxx_0[j];

                tdx_zz_xxy_0[j] = pa_z[j] * tdx_z_xxy_0[j] + 0.5 * fl1_fx * tdx_0_xxy_0[j];

                tdy_zz_xxy_0[j] = pa_z[j] * tdy_z_xxy_0[j] + 0.5 * fl1_fx * tdy_0_xxy_0[j];

                tdz_zz_xxy_0[j] = pa_z[j] * tdz_z_xxy_0[j] + 0.5 * fl1_fx * tdz_0_xxy_0[j] + 0.5 * fl1_fx * ts_z_xxy_0[j];

                tdx_zz_xxz_0[j] = pa_z[j] * tdx_z_xxz_0[j] + 0.5 * fl1_fx * tdx_0_xxz_0[j] + 0.5 * fl1_fx * tdx_z_xx_0[j];

                tdy_zz_xxz_0[j] = pa_z[j] * tdy_z_xxz_0[j] + 0.5 * fl1_fx * tdy_0_xxz_0[j] + 0.5 * fl1_fx * tdy_z_xx_0[j];

                tdz_zz_xxz_0[j] = pa_z[j] * tdz_z_xxz_0[j] + 0.5 * fl1_fx * tdz_0_xxz_0[j] + 0.5 * fl1_fx * tdz_z_xx_0[j] + 0.5 * fl1_fx * ts_z_xxz_0[j];

                tdx_zz_xyy_0[j] = pa_z[j] * tdx_z_xyy_0[j] + 0.5 * fl1_fx * tdx_0_xyy_0[j];

                tdy_zz_xyy_0[j] = pa_z[j] * tdy_z_xyy_0[j] + 0.5 * fl1_fx * tdy_0_xyy_0[j];

                tdz_zz_xyy_0[j] = pa_z[j] * tdz_z_xyy_0[j] + 0.5 * fl1_fx * tdz_0_xyy_0[j] + 0.5 * fl1_fx * ts_z_xyy_0[j];

                tdx_zz_xyz_0[j] = pa_z[j] * tdx_z_xyz_0[j] + 0.5 * fl1_fx * tdx_0_xyz_0[j] + 0.5 * fl1_fx * tdx_z_xy_0[j];

                tdy_zz_xyz_0[j] = pa_z[j] * tdy_z_xyz_0[j] + 0.5 * fl1_fx * tdy_0_xyz_0[j] + 0.5 * fl1_fx * tdy_z_xy_0[j];

                tdz_zz_xyz_0[j] = pa_z[j] * tdz_z_xyz_0[j] + 0.5 * fl1_fx * tdz_0_xyz_0[j] + 0.5 * fl1_fx * tdz_z_xy_0[j] + 0.5 * fl1_fx * ts_z_xyz_0[j];

                tdx_zz_xzz_0[j] = pa_z[j] * tdx_z_xzz_0[j] + 0.5 * fl1_fx * tdx_0_xzz_0[j] + fl1_fx * tdx_z_xz_0[j];

                tdy_zz_xzz_0[j] = pa_z[j] * tdy_z_xzz_0[j] + 0.5 * fl1_fx * tdy_0_xzz_0[j] + fl1_fx * tdy_z_xz_0[j];

                tdz_zz_xzz_0[j] = pa_z[j] * tdz_z_xzz_0[j] + 0.5 * fl1_fx * tdz_0_xzz_0[j] + fl1_fx * tdz_z_xz_0[j] + 0.5 * fl1_fx * ts_z_xzz_0[j];

                tdx_zz_yyy_0[j] = pa_z[j] * tdx_z_yyy_0[j] + 0.5 * fl1_fx * tdx_0_yyy_0[j];

                tdy_zz_yyy_0[j] = pa_z[j] * tdy_z_yyy_0[j] + 0.5 * fl1_fx * tdy_0_yyy_0[j];

                tdz_zz_yyy_0[j] = pa_z[j] * tdz_z_yyy_0[j] + 0.5 * fl1_fx * tdz_0_yyy_0[j] + 0.5 * fl1_fx * ts_z_yyy_0[j];

                tdx_zz_yyz_0[j] = pa_z[j] * tdx_z_yyz_0[j] + 0.5 * fl1_fx * tdx_0_yyz_0[j] + 0.5 * fl1_fx * tdx_z_yy_0[j];

                tdy_zz_yyz_0[j] = pa_z[j] * tdy_z_yyz_0[j] + 0.5 * fl1_fx * tdy_0_yyz_0[j] + 0.5 * fl1_fx * tdy_z_yy_0[j];

                tdz_zz_yyz_0[j] = pa_z[j] * tdz_z_yyz_0[j] + 0.5 * fl1_fx * tdz_0_yyz_0[j] + 0.5 * fl1_fx * tdz_z_yy_0[j] + 0.5 * fl1_fx * ts_z_yyz_0[j];

                tdx_zz_yzz_0[j] = pa_z[j] * tdx_z_yzz_0[j] + 0.5 * fl1_fx * tdx_0_yzz_0[j] + fl1_fx * tdx_z_yz_0[j];

                tdy_zz_yzz_0[j] = pa_z[j] * tdy_z_yzz_0[j] + 0.5 * fl1_fx * tdy_0_yzz_0[j] + fl1_fx * tdy_z_yz_0[j];

                tdz_zz_yzz_0[j] = pa_z[j] * tdz_z_yzz_0[j] + 0.5 * fl1_fx * tdz_0_yzz_0[j] + fl1_fx * tdz_z_yz_0[j] + 0.5 * fl1_fx * ts_z_yzz_0[j];

                tdx_zz_zzz_0[j] = pa_z[j] * tdx_z_zzz_0[j] + 0.5 * fl1_fx * tdx_0_zzz_0[j] + 1.5 * fl1_fx * tdx_z_zz_0[j];

                tdy_zz_zzz_0[j] = pa_z[j] * tdy_z_zzz_0[j] + 0.5 * fl1_fx * tdy_0_zzz_0[j] + 1.5 * fl1_fx * tdy_z_zz_0[j];

                tdz_zz_zzz_0[j] = pa_z[j] * tdz_z_zzz_0[j] + 0.5 * fl1_fx * tdz_0_zzz_0[j] + 1.5 * fl1_fx * tdz_z_zz_0[j] + 0.5 * fl1_fx * ts_z_zzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForFD(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForFD_0_45(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  nOSFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForFD_45_90(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   nOSFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        ediprecfunc::compElectricDipoleForFD_90_135(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    nOSFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        ediprecfunc::compElectricDipoleForFD_135_180(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     nOSFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 
    }

    void
    compElectricDipoleForFD_0_45(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const int32_t              nOSFactors,
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

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_3_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx); 

            auto tdy_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx); 

            auto tdz_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx); 

            auto tdx_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 1); 

            auto tdy_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 1); 

            auto tdz_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 1); 

            auto tdx_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 2); 

            auto tdy_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 2); 

            auto tdz_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 2); 

            auto tdx_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 3); 

            auto tdy_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 3); 

            auto tdz_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 3); 

            auto tdx_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 4); 

            auto tdy_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 4); 

            auto tdz_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 4); 

            auto tdx_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 5); 

            auto tdy_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 5); 

            auto tdz_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 5); 

            auto tdx_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 6); 

            auto tdy_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 6); 

            auto tdz_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 6); 

            auto tdx_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 7); 

            auto tdy_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 7); 

            auto tdz_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 7); 

            auto tdx_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 8); 

            auto tdy_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 8); 

            auto tdz_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 8); 

            auto tdx_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 9); 

            auto tdy_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 9); 

            auto tdz_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 9); 

            auto tdx_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 10); 

            auto tdy_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 10); 

            auto tdz_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 10); 

            auto tdx_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 11); 

            auto tdy_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 11); 

            auto tdz_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 11); 

            auto tdx_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 12); 

            auto tdy_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 12); 

            auto tdz_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 12); 

            auto tdx_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 13); 

            auto tdy_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 13); 

            auto tdz_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 13); 

            auto tdx_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 14); 

            auto tdy_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 14); 

            auto tdz_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 14); 

            auto tdx_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx); 

            auto tdy_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx); 

            auto tdz_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx); 

            auto tdx_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 1); 

            auto tdy_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 1); 

            auto tdz_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 1); 

            auto tdx_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 2); 

            auto tdy_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 2); 

            auto tdz_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 2); 

            auto tdx_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 3); 

            auto tdy_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 3); 

            auto tdz_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 3); 

            auto tdx_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 4); 

            auto tdy_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 4); 

            auto tdz_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 4); 

            auto tdx_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 5); 

            auto tdy_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 5); 

            auto tdz_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 5); 

            auto tdx_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 6); 

            auto tdy_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tdz_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tdx_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 7); 

            auto tdy_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tdz_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tdx_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 8); 

            auto tdy_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tdz_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tdx_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 9); 

            auto tdy_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tdz_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdx_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 10); 

            auto tdy_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tdz_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tdx_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 11); 

            auto tdy_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tdz_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tdx_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 12); 

            auto tdy_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tdz_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tdx_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 13); 

            auto tdy_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tdz_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tdx_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 14); 

            auto tdy_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tdz_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tdx_xx_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx); 

            auto tdy_xx_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx); 

            auto tdz_xx_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx); 

            auto tdx_xx_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 1); 

            auto tdy_xx_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 1); 

            auto tdz_xx_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 1); 

            auto tdx_xx_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 2); 

            auto tdy_xx_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 2); 

            auto tdz_xx_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 2); 

            auto tdx_xy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 3); 

            auto tdy_xy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 3); 

            auto tdz_xy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 3); 

            auto tdx_xy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 4); 

            auto tdy_xy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 4); 

            auto tdz_xy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 4); 

            auto tdx_xy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 5); 

            auto tdy_xy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 5); 

            auto tdz_xy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 5); 

            auto tdx_xz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 6); 

            auto tdy_xz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 6); 

            auto tdz_xz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 6); 

            auto tdx_xz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 7); 

            auto tdy_xz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 7); 

            auto tdz_xz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 7); 

            auto tdx_xz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 8); 

            auto tdy_xz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 8); 

            auto tdz_xz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 8); 

            auto ts_xx_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx); 

            auto ts_xx_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 1); 

            auto ts_xx_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 2); 

            auto ts_xx_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 3); 

            auto ts_xx_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 4); 

            auto ts_xx_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 5); 

            auto ts_xy_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 6); 

            auto ts_xy_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 7); 

            auto ts_xy_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 8); 

            auto ts_xy_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 9); 

            auto ts_xy_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 10); 

            auto ts_xy_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 11); 

            auto ts_xz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 12); 

            auto ts_xz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 13); 

            auto ts_xz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 14); 

            // set up pointers to integrals

            auto tdx_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx); 

            auto tdy_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx); 

            auto tdz_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx); 

            auto tdx_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 1); 

            auto tdy_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 1); 

            auto tdz_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 1); 

            auto tdx_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 2); 

            auto tdy_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 2); 

            auto tdz_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 2); 

            auto tdx_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 3); 

            auto tdy_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 3); 

            auto tdz_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 3); 

            auto tdx_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 4); 

            auto tdy_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 4); 

            auto tdz_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 4); 

            auto tdx_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 5); 

            auto tdy_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 5); 

            auto tdz_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 5); 

            auto tdx_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 6); 

            auto tdy_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 6); 

            auto tdz_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 6); 

            auto tdx_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 7); 

            auto tdy_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 7); 

            auto tdz_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 7); 

            auto tdx_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 8); 

            auto tdy_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 8); 

            auto tdz_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 8); 

            auto tdx_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 9); 

            auto tdy_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 9); 

            auto tdz_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 9); 

            auto tdx_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 10); 

            auto tdy_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 10); 

            auto tdz_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 10); 

            auto tdx_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 11); 

            auto tdy_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 11); 

            auto tdz_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 11); 

            auto tdx_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 12); 

            auto tdy_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 12); 

            auto tdz_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 12); 

            auto tdx_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 13); 

            auto tdy_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 13); 

            auto tdz_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 13); 

            auto tdx_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 14); 

            auto tdy_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 14); 

            auto tdz_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 14); 

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, tdx_x_xx_0, tdx_x_xy_0, tdx_x_xz_0, tdx_x_yy_0, tdx_x_yz_0, \
                                     tdx_x_zz_0, tdx_xx_x_0, tdx_xx_xx_0, tdx_xx_xy_0, tdx_xx_xz_0, tdx_xx_y_0, \
                                     tdx_xx_yy_0, tdx_xx_yz_0, tdx_xx_z_0, tdx_xx_zz_0, tdx_xxx_xx_0, tdx_xxx_xy_0, \
                                     tdx_xxx_xz_0, tdx_xxx_yy_0, tdx_xxx_yz_0, tdx_xxx_zz_0, tdx_xxy_xx_0, tdx_xxy_xy_0, \
                                     tdx_xxy_xz_0, tdx_xxy_yy_0, tdx_xxy_yz_0, tdx_xxy_zz_0, tdx_xxz_xx_0, tdx_xxz_xy_0, \
                                     tdx_xxz_xz_0, tdx_xy_x_0, tdx_xy_xx_0, tdx_xy_xy_0, tdx_xy_xz_0, tdx_xy_y_0, \
                                     tdx_xy_yy_0, tdx_xy_yz_0, tdx_xy_z_0, tdx_xy_zz_0, tdx_xz_x_0, tdx_xz_xx_0, \
                                     tdx_xz_xy_0, tdx_xz_xz_0, tdx_xz_y_0, tdx_xz_z_0, tdx_y_xx_0, tdx_y_xy_0, \
                                     tdx_y_xz_0, tdx_y_yy_0, tdx_y_yz_0, tdx_y_zz_0, tdx_z_xx_0, tdx_z_xy_0, tdx_z_xz_0, \
                                     tdy_x_xx_0, tdy_x_xy_0, tdy_x_xz_0, tdy_x_yy_0, tdy_x_yz_0, tdy_x_zz_0, tdy_xx_x_0, \
                                     tdy_xx_xx_0, tdy_xx_xy_0, tdy_xx_xz_0, tdy_xx_y_0, tdy_xx_yy_0, tdy_xx_yz_0, \
                                     tdy_xx_z_0, tdy_xx_zz_0, tdy_xxx_xx_0, tdy_xxx_xy_0, tdy_xxx_xz_0, tdy_xxx_yy_0, \
                                     tdy_xxx_yz_0, tdy_xxx_zz_0, tdy_xxy_xx_0, tdy_xxy_xy_0, tdy_xxy_xz_0, tdy_xxy_yy_0, \
                                     tdy_xxy_yz_0, tdy_xxy_zz_0, tdy_xxz_xx_0, tdy_xxz_xy_0, tdy_xxz_xz_0, tdy_xy_x_0, \
                                     tdy_xy_xx_0, tdy_xy_xy_0, tdy_xy_xz_0, tdy_xy_y_0, tdy_xy_yy_0, tdy_xy_yz_0, \
                                     tdy_xy_z_0, tdy_xy_zz_0, tdy_xz_x_0, tdy_xz_xx_0, tdy_xz_xy_0, tdy_xz_xz_0, \
                                     tdy_xz_y_0, tdy_xz_z_0, tdy_y_xx_0, tdy_y_xy_0, tdy_y_xz_0, tdy_y_yy_0, tdy_y_yz_0, \
                                     tdy_y_zz_0, tdy_z_xx_0, tdy_z_xy_0, tdy_z_xz_0, tdz_x_xx_0, tdz_x_xy_0, tdz_x_xz_0, \
                                     tdz_x_yy_0, tdz_x_yz_0, tdz_x_zz_0, tdz_xx_x_0, tdz_xx_xx_0, tdz_xx_xy_0, \
                                     tdz_xx_xz_0, tdz_xx_y_0, tdz_xx_yy_0, tdz_xx_yz_0, tdz_xx_z_0, tdz_xx_zz_0, \
                                     tdz_xxx_xx_0, tdz_xxx_xy_0, tdz_xxx_xz_0, tdz_xxx_yy_0, tdz_xxx_yz_0, tdz_xxx_zz_0, \
                                     tdz_xxy_xx_0, tdz_xxy_xy_0, tdz_xxy_xz_0, tdz_xxy_yy_0, tdz_xxy_yz_0, tdz_xxy_zz_0, \
                                     tdz_xxz_xx_0, tdz_xxz_xy_0, tdz_xxz_xz_0, tdz_xy_x_0, tdz_xy_xx_0, tdz_xy_xy_0, \
                                     tdz_xy_xz_0, tdz_xy_y_0, tdz_xy_yy_0, tdz_xy_yz_0, tdz_xy_z_0, tdz_xy_zz_0, \
                                     tdz_xz_x_0, tdz_xz_xx_0, tdz_xz_xy_0, tdz_xz_xz_0, tdz_xz_y_0, tdz_xz_z_0, \
                                     tdz_y_xx_0, tdz_y_xy_0, tdz_y_xz_0, tdz_y_yy_0, tdz_y_yz_0, tdz_y_zz_0, tdz_z_xx_0, \
                                     tdz_z_xy_0, tdz_z_xz_0, ts_xx_xx_0, ts_xx_xy_0, ts_xx_xz_0, ts_xx_yy_0, ts_xx_yz_0, \
                                     ts_xx_zz_0, ts_xy_xx_0, ts_xy_xy_0, ts_xy_xz_0, ts_xy_yy_0, ts_xy_yz_0, ts_xy_zz_0, \
                                     ts_xz_xx_0, ts_xz_xy_0, ts_xz_xz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxx_xx_0[j] = pa_x[j] * tdx_xx_xx_0[j] + fl1_fx * tdx_x_xx_0[j] + fl1_fx * tdx_xx_x_0[j] + 0.5 * fl1_fx * ts_xx_xx_0[j];

                tdy_xxx_xx_0[j] = pa_x[j] * tdy_xx_xx_0[j] + fl1_fx * tdy_x_xx_0[j] + fl1_fx * tdy_xx_x_0[j];

                tdz_xxx_xx_0[j] = pa_x[j] * tdz_xx_xx_0[j] + fl1_fx * tdz_x_xx_0[j] + fl1_fx * tdz_xx_x_0[j];

                tdx_xxx_xy_0[j] = pa_x[j] * tdx_xx_xy_0[j] + fl1_fx * tdx_x_xy_0[j] + 0.5 * fl1_fx * tdx_xx_y_0[j] + 0.5 * fl1_fx * ts_xx_xy_0[j];

                tdy_xxx_xy_0[j] = pa_x[j] * tdy_xx_xy_0[j] + fl1_fx * tdy_x_xy_0[j] + 0.5 * fl1_fx * tdy_xx_y_0[j];

                tdz_xxx_xy_0[j] = pa_x[j] * tdz_xx_xy_0[j] + fl1_fx * tdz_x_xy_0[j] + 0.5 * fl1_fx * tdz_xx_y_0[j];

                tdx_xxx_xz_0[j] = pa_x[j] * tdx_xx_xz_0[j] + fl1_fx * tdx_x_xz_0[j] + 0.5 * fl1_fx * tdx_xx_z_0[j] + 0.5 * fl1_fx * ts_xx_xz_0[j];

                tdy_xxx_xz_0[j] = pa_x[j] * tdy_xx_xz_0[j] + fl1_fx * tdy_x_xz_0[j] + 0.5 * fl1_fx * tdy_xx_z_0[j];

                tdz_xxx_xz_0[j] = pa_x[j] * tdz_xx_xz_0[j] + fl1_fx * tdz_x_xz_0[j] + 0.5 * fl1_fx * tdz_xx_z_0[j];

                tdx_xxx_yy_0[j] = pa_x[j] * tdx_xx_yy_0[j] + fl1_fx * tdx_x_yy_0[j] + 0.5 * fl1_fx * ts_xx_yy_0[j];

                tdy_xxx_yy_0[j] = pa_x[j] * tdy_xx_yy_0[j] + fl1_fx * tdy_x_yy_0[j];

                tdz_xxx_yy_0[j] = pa_x[j] * tdz_xx_yy_0[j] + fl1_fx * tdz_x_yy_0[j];

                tdx_xxx_yz_0[j] = pa_x[j] * tdx_xx_yz_0[j] + fl1_fx * tdx_x_yz_0[j] + 0.5 * fl1_fx * ts_xx_yz_0[j];

                tdy_xxx_yz_0[j] = pa_x[j] * tdy_xx_yz_0[j] + fl1_fx * tdy_x_yz_0[j];

                tdz_xxx_yz_0[j] = pa_x[j] * tdz_xx_yz_0[j] + fl1_fx * tdz_x_yz_0[j];

                tdx_xxx_zz_0[j] = pa_x[j] * tdx_xx_zz_0[j] + fl1_fx * tdx_x_zz_0[j] + 0.5 * fl1_fx * ts_xx_zz_0[j];

                tdy_xxx_zz_0[j] = pa_x[j] * tdy_xx_zz_0[j] + fl1_fx * tdy_x_zz_0[j];

                tdz_xxx_zz_0[j] = pa_x[j] * tdz_xx_zz_0[j] + fl1_fx * tdz_x_zz_0[j];

                tdx_xxy_xx_0[j] = pa_x[j] * tdx_xy_xx_0[j] + 0.5 * fl1_fx * tdx_y_xx_0[j] + fl1_fx * tdx_xy_x_0[j] + 0.5 * fl1_fx * ts_xy_xx_0[j];

                tdy_xxy_xx_0[j] = pa_x[j] * tdy_xy_xx_0[j] + 0.5 * fl1_fx * tdy_y_xx_0[j] + fl1_fx * tdy_xy_x_0[j];

                tdz_xxy_xx_0[j] = pa_x[j] * tdz_xy_xx_0[j] + 0.5 * fl1_fx * tdz_y_xx_0[j] + fl1_fx * tdz_xy_x_0[j];

                tdx_xxy_xy_0[j] = pa_x[j] * tdx_xy_xy_0[j] + 0.5 * fl1_fx * tdx_y_xy_0[j] + 0.5 * fl1_fx * tdx_xy_y_0[j] + 0.5 * fl1_fx * ts_xy_xy_0[j];

                tdy_xxy_xy_0[j] = pa_x[j] * tdy_xy_xy_0[j] + 0.5 * fl1_fx * tdy_y_xy_0[j] + 0.5 * fl1_fx * tdy_xy_y_0[j];

                tdz_xxy_xy_0[j] = pa_x[j] * tdz_xy_xy_0[j] + 0.5 * fl1_fx * tdz_y_xy_0[j] + 0.5 * fl1_fx * tdz_xy_y_0[j];

                tdx_xxy_xz_0[j] = pa_x[j] * tdx_xy_xz_0[j] + 0.5 * fl1_fx * tdx_y_xz_0[j] + 0.5 * fl1_fx * tdx_xy_z_0[j] + 0.5 * fl1_fx * ts_xy_xz_0[j];

                tdy_xxy_xz_0[j] = pa_x[j] * tdy_xy_xz_0[j] + 0.5 * fl1_fx * tdy_y_xz_0[j] + 0.5 * fl1_fx * tdy_xy_z_0[j];

                tdz_xxy_xz_0[j] = pa_x[j] * tdz_xy_xz_0[j] + 0.5 * fl1_fx * tdz_y_xz_0[j] + 0.5 * fl1_fx * tdz_xy_z_0[j];

                tdx_xxy_yy_0[j] = pa_x[j] * tdx_xy_yy_0[j] + 0.5 * fl1_fx * tdx_y_yy_0[j] + 0.5 * fl1_fx * ts_xy_yy_0[j];

                tdy_xxy_yy_0[j] = pa_x[j] * tdy_xy_yy_0[j] + 0.5 * fl1_fx * tdy_y_yy_0[j];

                tdz_xxy_yy_0[j] = pa_x[j] * tdz_xy_yy_0[j] + 0.5 * fl1_fx * tdz_y_yy_0[j];

                tdx_xxy_yz_0[j] = pa_x[j] * tdx_xy_yz_0[j] + 0.5 * fl1_fx * tdx_y_yz_0[j] + 0.5 * fl1_fx * ts_xy_yz_0[j];

                tdy_xxy_yz_0[j] = pa_x[j] * tdy_xy_yz_0[j] + 0.5 * fl1_fx * tdy_y_yz_0[j];

                tdz_xxy_yz_0[j] = pa_x[j] * tdz_xy_yz_0[j] + 0.5 * fl1_fx * tdz_y_yz_0[j];

                tdx_xxy_zz_0[j] = pa_x[j] * tdx_xy_zz_0[j] + 0.5 * fl1_fx * tdx_y_zz_0[j] + 0.5 * fl1_fx * ts_xy_zz_0[j];

                tdy_xxy_zz_0[j] = pa_x[j] * tdy_xy_zz_0[j] + 0.5 * fl1_fx * tdy_y_zz_0[j];

                tdz_xxy_zz_0[j] = pa_x[j] * tdz_xy_zz_0[j] + 0.5 * fl1_fx * tdz_y_zz_0[j];

                tdx_xxz_xx_0[j] = pa_x[j] * tdx_xz_xx_0[j] + 0.5 * fl1_fx * tdx_z_xx_0[j] + fl1_fx * tdx_xz_x_0[j] + 0.5 * fl1_fx * ts_xz_xx_0[j];

                tdy_xxz_xx_0[j] = pa_x[j] * tdy_xz_xx_0[j] + 0.5 * fl1_fx * tdy_z_xx_0[j] + fl1_fx * tdy_xz_x_0[j];

                tdz_xxz_xx_0[j] = pa_x[j] * tdz_xz_xx_0[j] + 0.5 * fl1_fx * tdz_z_xx_0[j] + fl1_fx * tdz_xz_x_0[j];

                tdx_xxz_xy_0[j] = pa_x[j] * tdx_xz_xy_0[j] + 0.5 * fl1_fx * tdx_z_xy_0[j] + 0.5 * fl1_fx * tdx_xz_y_0[j] + 0.5 * fl1_fx * ts_xz_xy_0[j];

                tdy_xxz_xy_0[j] = pa_x[j] * tdy_xz_xy_0[j] + 0.5 * fl1_fx * tdy_z_xy_0[j] + 0.5 * fl1_fx * tdy_xz_y_0[j];

                tdz_xxz_xy_0[j] = pa_x[j] * tdz_xz_xy_0[j] + 0.5 * fl1_fx * tdz_z_xy_0[j] + 0.5 * fl1_fx * tdz_xz_y_0[j];

                tdx_xxz_xz_0[j] = pa_x[j] * tdx_xz_xz_0[j] + 0.5 * fl1_fx * tdx_z_xz_0[j] + 0.5 * fl1_fx * tdx_xz_z_0[j] + 0.5 * fl1_fx * ts_xz_xz_0[j];

                tdy_xxz_xz_0[j] = pa_x[j] * tdy_xz_xz_0[j] + 0.5 * fl1_fx * tdy_z_xz_0[j] + 0.5 * fl1_fx * tdy_xz_z_0[j];

                tdz_xxz_xz_0[j] = pa_x[j] * tdz_xz_xz_0[j] + 0.5 * fl1_fx * tdz_z_xz_0[j] + 0.5 * fl1_fx * tdz_xz_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForFD_45_90(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
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

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_3_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 15); 

            auto tdy_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 15); 

            auto tdz_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 15); 

            auto tdx_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 16); 

            auto tdy_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 16); 

            auto tdz_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 16); 

            auto tdx_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 17); 

            auto tdy_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 17); 

            auto tdz_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 17); 

            auto tdx_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 18); 

            auto tdy_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tdz_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tdx_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 19); 

            auto tdy_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tdz_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tdx_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 20); 

            auto tdy_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tdz_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tdx_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 21); 

            auto tdy_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tdz_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tdx_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 22); 

            auto tdy_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tdz_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tdx_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 23); 

            auto tdy_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tdz_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tdx_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 24); 

            auto tdy_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tdz_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tdx_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 25); 

            auto tdy_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tdz_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tdx_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 26); 

            auto tdy_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tdz_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tdx_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 27); 

            auto tdy_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tdz_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tdx_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 28); 

            auto tdy_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tdz_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tdx_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 29); 

            auto tdy_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 29); 

            auto tdz_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tdx_z_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 15); 

            auto tdy_z_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tdz_z_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tdx_z_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 16); 

            auto tdy_z_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tdz_z_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tdx_z_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 17); 

            auto tdy_z_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tdz_z_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            auto tdx_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 9); 

            auto tdy_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 9); 

            auto tdz_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdx_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 10); 

            auto tdy_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 10); 

            auto tdz_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 10); 

            auto tdx_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 11); 

            auto tdy_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 11); 

            auto tdz_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 11); 

            auto tdx_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 12); 

            auto tdy_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 12); 

            auto tdz_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 12); 

            auto tdx_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 13); 

            auto tdy_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 13); 

            auto tdz_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 13); 

            auto tdx_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 14); 

            auto tdy_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 14); 

            auto tdz_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 14); 

            auto ts_xz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 15); 

            auto ts_xz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 16); 

            auto ts_xz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 17); 

            auto ts_yy_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 18); 

            auto ts_yy_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 19); 

            auto ts_yy_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 20); 

            auto ts_yy_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 21); 

            auto ts_yy_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 22); 

            auto ts_yy_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 23); 

            auto ts_yz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 24); 

            auto ts_yz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 25); 

            auto ts_yz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 26); 

            auto ts_yz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 27); 

            auto ts_yz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 28); 

            auto ts_yz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 29); 

            // set up pointers to integrals

            auto tdx_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 15); 

            auto tdy_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 15); 

            auto tdz_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 15); 

            auto tdx_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 16); 

            auto tdy_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 16); 

            auto tdz_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 16); 

            auto tdx_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 17); 

            auto tdy_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 17); 

            auto tdz_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 17); 

            auto tdx_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 18); 

            auto tdy_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 18); 

            auto tdz_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 18); 

            auto tdx_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 19); 

            auto tdy_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 19); 

            auto tdz_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 19); 

            auto tdx_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 20); 

            auto tdy_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 20); 

            auto tdz_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 20); 

            auto tdx_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 21); 

            auto tdy_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 21); 

            auto tdz_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 21); 

            auto tdx_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 22); 

            auto tdy_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 22); 

            auto tdz_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 22); 

            auto tdx_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 23); 

            auto tdy_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 23); 

            auto tdz_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 23); 

            auto tdx_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 24); 

            auto tdy_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 24); 

            auto tdz_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 24); 

            auto tdx_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 25); 

            auto tdy_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 25); 

            auto tdz_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 25); 

            auto tdx_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 26); 

            auto tdy_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 26); 

            auto tdz_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 26); 

            auto tdx_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 27); 

            auto tdy_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 27); 

            auto tdz_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 27); 

            auto tdx_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 28); 

            auto tdy_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 28); 

            auto tdz_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 28); 

            auto tdx_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 29); 

            auto tdy_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 29); 

            auto tdz_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 29); 

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, tdx_xxz_yy_0, tdx_xxz_yz_0, tdx_xxz_zz_0, tdx_xyy_xx_0, \
                                     tdx_xyy_xy_0, tdx_xyy_xz_0, tdx_xyy_yy_0, tdx_xyy_yz_0, tdx_xyy_zz_0, tdx_xyz_xx_0, \
                                     tdx_xyz_xy_0, tdx_xyz_xz_0, tdx_xyz_yy_0, tdx_xyz_yz_0, tdx_xyz_zz_0, tdx_xz_yy_0, \
                                     tdx_xz_yz_0, tdx_xz_zz_0, tdx_yy_x_0, tdx_yy_xx_0, tdx_yy_xy_0, tdx_yy_xz_0, \
                                     tdx_yy_y_0, tdx_yy_yy_0, tdx_yy_yz_0, tdx_yy_z_0, tdx_yy_zz_0, tdx_yz_x_0, \
                                     tdx_yz_xx_0, tdx_yz_xy_0, tdx_yz_xz_0, tdx_yz_y_0, tdx_yz_yy_0, tdx_yz_yz_0, \
                                     tdx_yz_z_0, tdx_yz_zz_0, tdx_z_yy_0, tdx_z_yz_0, tdx_z_zz_0, tdy_xxz_yy_0, \
                                     tdy_xxz_yz_0, tdy_xxz_zz_0, tdy_xyy_xx_0, tdy_xyy_xy_0, tdy_xyy_xz_0, tdy_xyy_yy_0, \
                                     tdy_xyy_yz_0, tdy_xyy_zz_0, tdy_xyz_xx_0, tdy_xyz_xy_0, tdy_xyz_xz_0, tdy_xyz_yy_0, \
                                     tdy_xyz_yz_0, tdy_xyz_zz_0, tdy_xz_yy_0, tdy_xz_yz_0, tdy_xz_zz_0, tdy_yy_x_0, \
                                     tdy_yy_xx_0, tdy_yy_xy_0, tdy_yy_xz_0, tdy_yy_y_0, tdy_yy_yy_0, tdy_yy_yz_0, \
                                     tdy_yy_z_0, tdy_yy_zz_0, tdy_yz_x_0, tdy_yz_xx_0, tdy_yz_xy_0, tdy_yz_xz_0, \
                                     tdy_yz_y_0, tdy_yz_yy_0, tdy_yz_yz_0, tdy_yz_z_0, tdy_yz_zz_0, tdy_z_yy_0, \
                                     tdy_z_yz_0, tdy_z_zz_0, tdz_xxz_yy_0, tdz_xxz_yz_0, tdz_xxz_zz_0, tdz_xyy_xx_0, \
                                     tdz_xyy_xy_0, tdz_xyy_xz_0, tdz_xyy_yy_0, tdz_xyy_yz_0, tdz_xyy_zz_0, tdz_xyz_xx_0, \
                                     tdz_xyz_xy_0, tdz_xyz_xz_0, tdz_xyz_yy_0, tdz_xyz_yz_0, tdz_xyz_zz_0, tdz_xz_yy_0, \
                                     tdz_xz_yz_0, tdz_xz_zz_0, tdz_yy_x_0, tdz_yy_xx_0, tdz_yy_xy_0, tdz_yy_xz_0, \
                                     tdz_yy_y_0, tdz_yy_yy_0, tdz_yy_yz_0, tdz_yy_z_0, tdz_yy_zz_0, tdz_yz_x_0, \
                                     tdz_yz_xx_0, tdz_yz_xy_0, tdz_yz_xz_0, tdz_yz_y_0, tdz_yz_yy_0, tdz_yz_yz_0, \
                                     tdz_yz_z_0, tdz_yz_zz_0, tdz_z_yy_0, tdz_z_yz_0, tdz_z_zz_0, ts_xz_yy_0, \
                                     ts_xz_yz_0, ts_xz_zz_0, ts_yy_xx_0, ts_yy_xy_0, ts_yy_xz_0, ts_yy_yy_0, ts_yy_yz_0, \
                                     ts_yy_zz_0, ts_yz_xx_0, ts_yz_xy_0, ts_yz_xz_0, ts_yz_yy_0, ts_yz_yz_0, ts_yz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxz_yy_0[j] = pa_x[j] * tdx_xz_yy_0[j] + 0.5 * fl1_fx * tdx_z_yy_0[j] + 0.5 * fl1_fx * ts_xz_yy_0[j];

                tdy_xxz_yy_0[j] = pa_x[j] * tdy_xz_yy_0[j] + 0.5 * fl1_fx * tdy_z_yy_0[j];

                tdz_xxz_yy_0[j] = pa_x[j] * tdz_xz_yy_0[j] + 0.5 * fl1_fx * tdz_z_yy_0[j];

                tdx_xxz_yz_0[j] = pa_x[j] * tdx_xz_yz_0[j] + 0.5 * fl1_fx * tdx_z_yz_0[j] + 0.5 * fl1_fx * ts_xz_yz_0[j];

                tdy_xxz_yz_0[j] = pa_x[j] * tdy_xz_yz_0[j] + 0.5 * fl1_fx * tdy_z_yz_0[j];

                tdz_xxz_yz_0[j] = pa_x[j] * tdz_xz_yz_0[j] + 0.5 * fl1_fx * tdz_z_yz_0[j];

                tdx_xxz_zz_0[j] = pa_x[j] * tdx_xz_zz_0[j] + 0.5 * fl1_fx * tdx_z_zz_0[j] + 0.5 * fl1_fx * ts_xz_zz_0[j];

                tdy_xxz_zz_0[j] = pa_x[j] * tdy_xz_zz_0[j] + 0.5 * fl1_fx * tdy_z_zz_0[j];

                tdz_xxz_zz_0[j] = pa_x[j] * tdz_xz_zz_0[j] + 0.5 * fl1_fx * tdz_z_zz_0[j];

                tdx_xyy_xx_0[j] = pa_x[j] * tdx_yy_xx_0[j] + fl1_fx * tdx_yy_x_0[j] + 0.5 * fl1_fx * ts_yy_xx_0[j];

                tdy_xyy_xx_0[j] = pa_x[j] * tdy_yy_xx_0[j] + fl1_fx * tdy_yy_x_0[j];

                tdz_xyy_xx_0[j] = pa_x[j] * tdz_yy_xx_0[j] + fl1_fx * tdz_yy_x_0[j];

                tdx_xyy_xy_0[j] = pa_x[j] * tdx_yy_xy_0[j] + 0.5 * fl1_fx * tdx_yy_y_0[j] + 0.5 * fl1_fx * ts_yy_xy_0[j];

                tdy_xyy_xy_0[j] = pa_x[j] * tdy_yy_xy_0[j] + 0.5 * fl1_fx * tdy_yy_y_0[j];

                tdz_xyy_xy_0[j] = pa_x[j] * tdz_yy_xy_0[j] + 0.5 * fl1_fx * tdz_yy_y_0[j];

                tdx_xyy_xz_0[j] = pa_x[j] * tdx_yy_xz_0[j] + 0.5 * fl1_fx * tdx_yy_z_0[j] + 0.5 * fl1_fx * ts_yy_xz_0[j];

                tdy_xyy_xz_0[j] = pa_x[j] * tdy_yy_xz_0[j] + 0.5 * fl1_fx * tdy_yy_z_0[j];

                tdz_xyy_xz_0[j] = pa_x[j] * tdz_yy_xz_0[j] + 0.5 * fl1_fx * tdz_yy_z_0[j];

                tdx_xyy_yy_0[j] = pa_x[j] * tdx_yy_yy_0[j] + 0.5 * fl1_fx * ts_yy_yy_0[j];

                tdy_xyy_yy_0[j] = pa_x[j] * tdy_yy_yy_0[j];

                tdz_xyy_yy_0[j] = pa_x[j] * tdz_yy_yy_0[j];

                tdx_xyy_yz_0[j] = pa_x[j] * tdx_yy_yz_0[j] + 0.5 * fl1_fx * ts_yy_yz_0[j];

                tdy_xyy_yz_0[j] = pa_x[j] * tdy_yy_yz_0[j];

                tdz_xyy_yz_0[j] = pa_x[j] * tdz_yy_yz_0[j];

                tdx_xyy_zz_0[j] = pa_x[j] * tdx_yy_zz_0[j] + 0.5 * fl1_fx * ts_yy_zz_0[j];

                tdy_xyy_zz_0[j] = pa_x[j] * tdy_yy_zz_0[j];

                tdz_xyy_zz_0[j] = pa_x[j] * tdz_yy_zz_0[j];

                tdx_xyz_xx_0[j] = pa_x[j] * tdx_yz_xx_0[j] + fl1_fx * tdx_yz_x_0[j] + 0.5 * fl1_fx * ts_yz_xx_0[j];

                tdy_xyz_xx_0[j] = pa_x[j] * tdy_yz_xx_0[j] + fl1_fx * tdy_yz_x_0[j];

                tdz_xyz_xx_0[j] = pa_x[j] * tdz_yz_xx_0[j] + fl1_fx * tdz_yz_x_0[j];

                tdx_xyz_xy_0[j] = pa_x[j] * tdx_yz_xy_0[j] + 0.5 * fl1_fx * tdx_yz_y_0[j] + 0.5 * fl1_fx * ts_yz_xy_0[j];

                tdy_xyz_xy_0[j] = pa_x[j] * tdy_yz_xy_0[j] + 0.5 * fl1_fx * tdy_yz_y_0[j];

                tdz_xyz_xy_0[j] = pa_x[j] * tdz_yz_xy_0[j] + 0.5 * fl1_fx * tdz_yz_y_0[j];

                tdx_xyz_xz_0[j] = pa_x[j] * tdx_yz_xz_0[j] + 0.5 * fl1_fx * tdx_yz_z_0[j] + 0.5 * fl1_fx * ts_yz_xz_0[j];

                tdy_xyz_xz_0[j] = pa_x[j] * tdy_yz_xz_0[j] + 0.5 * fl1_fx * tdy_yz_z_0[j];

                tdz_xyz_xz_0[j] = pa_x[j] * tdz_yz_xz_0[j] + 0.5 * fl1_fx * tdz_yz_z_0[j];

                tdx_xyz_yy_0[j] = pa_x[j] * tdx_yz_yy_0[j] + 0.5 * fl1_fx * ts_yz_yy_0[j];

                tdy_xyz_yy_0[j] = pa_x[j] * tdy_yz_yy_0[j];

                tdz_xyz_yy_0[j] = pa_x[j] * tdz_yz_yy_0[j];

                tdx_xyz_yz_0[j] = pa_x[j] * tdx_yz_yz_0[j] + 0.5 * fl1_fx * ts_yz_yz_0[j];

                tdy_xyz_yz_0[j] = pa_x[j] * tdy_yz_yz_0[j];

                tdz_xyz_yz_0[j] = pa_x[j] * tdz_yz_yz_0[j];

                tdx_xyz_zz_0[j] = pa_x[j] * tdx_yz_zz_0[j] + 0.5 * fl1_fx * ts_yz_zz_0[j];

                tdy_xyz_zz_0[j] = pa_x[j] * tdy_yz_zz_0[j];

                tdz_xyz_zz_0[j] = pa_x[j] * tdz_yz_zz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForFD_90_135(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
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

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_3_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tdx_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 18); 

            auto tdy_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tdz_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tdx_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 19); 

            auto tdy_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tdz_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tdx_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 20); 

            auto tdy_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tdz_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tdx_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 21); 

            auto tdy_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tdz_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tdx_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 22); 

            auto tdy_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tdz_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tdx_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 23); 

            auto tdy_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tdz_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tdx_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 24); 

            auto tdy_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tdz_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tdx_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 25); 

            auto tdy_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tdz_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tdx_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 26); 

            auto tdy_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tdz_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 26); 

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

            auto tdx_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 6); 

            auto tdy_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tdz_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tdx_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 7); 

            auto tdy_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tdz_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tdx_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 8); 

            auto tdy_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tdz_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tdx_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 9); 

            auto tdy_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tdz_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdx_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 10); 

            auto tdy_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tdz_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tdx_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 11); 

            auto tdy_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tdz_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tdx_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 12); 

            auto tdy_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tdz_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tdx_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 13); 

            auto tdy_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tdz_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tdx_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 14); 

            auto tdy_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tdz_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tdx_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 9); 

            auto tdy_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 9); 

            auto tdz_yy_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdx_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 10); 

            auto tdy_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 10); 

            auto tdz_yy_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 10); 

            auto tdx_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 11); 

            auto tdy_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 11); 

            auto tdz_yy_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 11); 

            auto tdx_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 12); 

            auto tdy_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 12); 

            auto tdz_yz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 12); 

            auto tdx_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 15); 

            auto tdy_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 15); 

            auto tdz_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 15); 

            auto tdx_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 16); 

            auto tdy_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 16); 

            auto tdz_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 16); 

            auto tdx_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 17); 

            auto tdy_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 17); 

            auto tdz_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 17); 

            auto ts_yy_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 18); 

            auto ts_yy_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 19); 

            auto ts_yy_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 20); 

            auto ts_yy_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 21); 

            auto ts_yy_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 22); 

            auto ts_yy_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 23); 

            auto ts_yz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 24); 

            auto ts_yz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 25); 

            auto ts_yz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 26); 

            auto ts_zz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 30); 

            auto ts_zz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 31); 

            auto ts_zz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 32); 

            auto ts_zz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 33); 

            auto ts_zz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 34); 

            auto ts_zz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 35); 

            // set up pointers to integrals

            auto tdx_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 30); 

            auto tdy_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 30); 

            auto tdz_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 30); 

            auto tdx_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 31); 

            auto tdy_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 31); 

            auto tdz_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 31); 

            auto tdx_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 32); 

            auto tdy_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 32); 

            auto tdz_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 32); 

            auto tdx_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 33); 

            auto tdy_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 33); 

            auto tdz_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 33); 

            auto tdx_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 34); 

            auto tdy_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 34); 

            auto tdz_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 34); 

            auto tdx_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 35); 

            auto tdy_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 35); 

            auto tdz_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 35); 

            auto tdx_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 36); 

            auto tdy_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tdz_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdx_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 37); 

            auto tdy_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tdz_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdx_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 38); 

            auto tdy_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tdz_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdx_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 39); 

            auto tdy_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tdz_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdx_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 40); 

            auto tdy_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tdz_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdx_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 41); 

            auto tdy_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tdz_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tdx_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 42); 

            auto tdy_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tdz_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdx_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 43); 

            auto tdy_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tdz_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdx_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 44); 

            auto tdy_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tdz_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fx, pa_x, pa_y, tdx_xzz_xx_0, tdx_xzz_xy_0, tdx_xzz_xz_0, tdx_xzz_yy_0, \
                                     tdx_xzz_yz_0, tdx_xzz_zz_0, tdx_y_xx_0, tdx_y_xy_0, tdx_y_xz_0, tdx_y_yy_0, \
                                     tdx_y_yz_0, tdx_y_zz_0, tdx_yy_x_0, tdx_yy_xx_0, tdx_yy_xy_0, tdx_yy_xz_0, \
                                     tdx_yy_y_0, tdx_yy_yy_0, tdx_yy_yz_0, tdx_yy_z_0, tdx_yy_zz_0, tdx_yyy_xx_0, \
                                     tdx_yyy_xy_0, tdx_yyy_xz_0, tdx_yyy_yy_0, tdx_yyy_yz_0, tdx_yyy_zz_0, tdx_yyz_xx_0, \
                                     tdx_yyz_xy_0, tdx_yyz_xz_0, tdx_yz_x_0, tdx_yz_xx_0, tdx_yz_xy_0, tdx_yz_xz_0, \
                                     tdx_z_xx_0, tdx_z_xy_0, tdx_z_xz_0, tdx_zz_x_0, tdx_zz_xx_0, tdx_zz_xy_0, \
                                     tdx_zz_xz_0, tdx_zz_y_0, tdx_zz_yy_0, tdx_zz_yz_0, tdx_zz_z_0, tdx_zz_zz_0, \
                                     tdy_xzz_xx_0, tdy_xzz_xy_0, tdy_xzz_xz_0, tdy_xzz_yy_0, tdy_xzz_yz_0, tdy_xzz_zz_0, \
                                     tdy_y_xx_0, tdy_y_xy_0, tdy_y_xz_0, tdy_y_yy_0, tdy_y_yz_0, tdy_y_zz_0, tdy_yy_x_0, \
                                     tdy_yy_xx_0, tdy_yy_xy_0, tdy_yy_xz_0, tdy_yy_y_0, tdy_yy_yy_0, tdy_yy_yz_0, \
                                     tdy_yy_z_0, tdy_yy_zz_0, tdy_yyy_xx_0, tdy_yyy_xy_0, tdy_yyy_xz_0, tdy_yyy_yy_0, \
                                     tdy_yyy_yz_0, tdy_yyy_zz_0, tdy_yyz_xx_0, tdy_yyz_xy_0, tdy_yyz_xz_0, tdy_yz_x_0, \
                                     tdy_yz_xx_0, tdy_yz_xy_0, tdy_yz_xz_0, tdy_z_xx_0, tdy_z_xy_0, tdy_z_xz_0, \
                                     tdy_zz_x_0, tdy_zz_xx_0, tdy_zz_xy_0, tdy_zz_xz_0, tdy_zz_y_0, tdy_zz_yy_0, \
                                     tdy_zz_yz_0, tdy_zz_z_0, tdy_zz_zz_0, tdz_xzz_xx_0, tdz_xzz_xy_0, tdz_xzz_xz_0, \
                                     tdz_xzz_yy_0, tdz_xzz_yz_0, tdz_xzz_zz_0, tdz_y_xx_0, tdz_y_xy_0, tdz_y_xz_0, \
                                     tdz_y_yy_0, tdz_y_yz_0, tdz_y_zz_0, tdz_yy_x_0, tdz_yy_xx_0, tdz_yy_xy_0, \
                                     tdz_yy_xz_0, tdz_yy_y_0, tdz_yy_yy_0, tdz_yy_yz_0, tdz_yy_z_0, tdz_yy_zz_0, \
                                     tdz_yyy_xx_0, tdz_yyy_xy_0, tdz_yyy_xz_0, tdz_yyy_yy_0, tdz_yyy_yz_0, tdz_yyy_zz_0, \
                                     tdz_yyz_xx_0, tdz_yyz_xy_0, tdz_yyz_xz_0, tdz_yz_x_0, tdz_yz_xx_0, tdz_yz_xy_0, \
                                     tdz_yz_xz_0, tdz_z_xx_0, tdz_z_xy_0, tdz_z_xz_0, tdz_zz_x_0, tdz_zz_xx_0, \
                                     tdz_zz_xy_0, tdz_zz_xz_0, tdz_zz_y_0, tdz_zz_yy_0, tdz_zz_yz_0, tdz_zz_z_0, \
                                     tdz_zz_zz_0, ts_yy_xx_0, ts_yy_xy_0, ts_yy_xz_0, ts_yy_yy_0, ts_yy_yz_0, ts_yy_zz_0, \
                                     ts_yz_xx_0, ts_yz_xy_0, ts_yz_xz_0, ts_zz_xx_0, ts_zz_xy_0, ts_zz_xz_0, ts_zz_yy_0, \
                                     ts_zz_yz_0, ts_zz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xzz_xx_0[j] = pa_x[j] * tdx_zz_xx_0[j] + fl1_fx * tdx_zz_x_0[j] + 0.5 * fl1_fx * ts_zz_xx_0[j];

                tdy_xzz_xx_0[j] = pa_x[j] * tdy_zz_xx_0[j] + fl1_fx * tdy_zz_x_0[j];

                tdz_xzz_xx_0[j] = pa_x[j] * tdz_zz_xx_0[j] + fl1_fx * tdz_zz_x_0[j];

                tdx_xzz_xy_0[j] = pa_x[j] * tdx_zz_xy_0[j] + 0.5 * fl1_fx * tdx_zz_y_0[j] + 0.5 * fl1_fx * ts_zz_xy_0[j];

                tdy_xzz_xy_0[j] = pa_x[j] * tdy_zz_xy_0[j] + 0.5 * fl1_fx * tdy_zz_y_0[j];

                tdz_xzz_xy_0[j] = pa_x[j] * tdz_zz_xy_0[j] + 0.5 * fl1_fx * tdz_zz_y_0[j];

                tdx_xzz_xz_0[j] = pa_x[j] * tdx_zz_xz_0[j] + 0.5 * fl1_fx * tdx_zz_z_0[j] + 0.5 * fl1_fx * ts_zz_xz_0[j];

                tdy_xzz_xz_0[j] = pa_x[j] * tdy_zz_xz_0[j] + 0.5 * fl1_fx * tdy_zz_z_0[j];

                tdz_xzz_xz_0[j] = pa_x[j] * tdz_zz_xz_0[j] + 0.5 * fl1_fx * tdz_zz_z_0[j];

                tdx_xzz_yy_0[j] = pa_x[j] * tdx_zz_yy_0[j] + 0.5 * fl1_fx * ts_zz_yy_0[j];

                tdy_xzz_yy_0[j] = pa_x[j] * tdy_zz_yy_0[j];

                tdz_xzz_yy_0[j] = pa_x[j] * tdz_zz_yy_0[j];

                tdx_xzz_yz_0[j] = pa_x[j] * tdx_zz_yz_0[j] + 0.5 * fl1_fx * ts_zz_yz_0[j];

                tdy_xzz_yz_0[j] = pa_x[j] * tdy_zz_yz_0[j];

                tdz_xzz_yz_0[j] = pa_x[j] * tdz_zz_yz_0[j];

                tdx_xzz_zz_0[j] = pa_x[j] * tdx_zz_zz_0[j] + 0.5 * fl1_fx * ts_zz_zz_0[j];

                tdy_xzz_zz_0[j] = pa_x[j] * tdy_zz_zz_0[j];

                tdz_xzz_zz_0[j] = pa_x[j] * tdz_zz_zz_0[j];

                tdx_yyy_xx_0[j] = pa_y[j] * tdx_yy_xx_0[j] + fl1_fx * tdx_y_xx_0[j];

                tdy_yyy_xx_0[j] = pa_y[j] * tdy_yy_xx_0[j] + fl1_fx * tdy_y_xx_0[j] + 0.5 * fl1_fx * ts_yy_xx_0[j];

                tdz_yyy_xx_0[j] = pa_y[j] * tdz_yy_xx_0[j] + fl1_fx * tdz_y_xx_0[j];

                tdx_yyy_xy_0[j] = pa_y[j] * tdx_yy_xy_0[j] + fl1_fx * tdx_y_xy_0[j] + 0.5 * fl1_fx * tdx_yy_x_0[j];

                tdy_yyy_xy_0[j] = pa_y[j] * tdy_yy_xy_0[j] + fl1_fx * tdy_y_xy_0[j] + 0.5 * fl1_fx * tdy_yy_x_0[j] + 0.5 * fl1_fx * ts_yy_xy_0[j];

                tdz_yyy_xy_0[j] = pa_y[j] * tdz_yy_xy_0[j] + fl1_fx * tdz_y_xy_0[j] + 0.5 * fl1_fx * tdz_yy_x_0[j];

                tdx_yyy_xz_0[j] = pa_y[j] * tdx_yy_xz_0[j] + fl1_fx * tdx_y_xz_0[j];

                tdy_yyy_xz_0[j] = pa_y[j] * tdy_yy_xz_0[j] + fl1_fx * tdy_y_xz_0[j] + 0.5 * fl1_fx * ts_yy_xz_0[j];

                tdz_yyy_xz_0[j] = pa_y[j] * tdz_yy_xz_0[j] + fl1_fx * tdz_y_xz_0[j];

                tdx_yyy_yy_0[j] = pa_y[j] * tdx_yy_yy_0[j] + fl1_fx * tdx_y_yy_0[j] + fl1_fx * tdx_yy_y_0[j];

                tdy_yyy_yy_0[j] = pa_y[j] * tdy_yy_yy_0[j] + fl1_fx * tdy_y_yy_0[j] + fl1_fx * tdy_yy_y_0[j] + 0.5 * fl1_fx * ts_yy_yy_0[j];

                tdz_yyy_yy_0[j] = pa_y[j] * tdz_yy_yy_0[j] + fl1_fx * tdz_y_yy_0[j] + fl1_fx * tdz_yy_y_0[j];

                tdx_yyy_yz_0[j] = pa_y[j] * tdx_yy_yz_0[j] + fl1_fx * tdx_y_yz_0[j] + 0.5 * fl1_fx * tdx_yy_z_0[j];

                tdy_yyy_yz_0[j] = pa_y[j] * tdy_yy_yz_0[j] + fl1_fx * tdy_y_yz_0[j] + 0.5 * fl1_fx * tdy_yy_z_0[j] + 0.5 * fl1_fx * ts_yy_yz_0[j];

                tdz_yyy_yz_0[j] = pa_y[j] * tdz_yy_yz_0[j] + fl1_fx * tdz_y_yz_0[j] + 0.5 * fl1_fx * tdz_yy_z_0[j];

                tdx_yyy_zz_0[j] = pa_y[j] * tdx_yy_zz_0[j] + fl1_fx * tdx_y_zz_0[j];

                tdy_yyy_zz_0[j] = pa_y[j] * tdy_yy_zz_0[j] + fl1_fx * tdy_y_zz_0[j] + 0.5 * fl1_fx * ts_yy_zz_0[j];

                tdz_yyy_zz_0[j] = pa_y[j] * tdz_yy_zz_0[j] + fl1_fx * tdz_y_zz_0[j];

                tdx_yyz_xx_0[j] = pa_y[j] * tdx_yz_xx_0[j] + 0.5 * fl1_fx * tdx_z_xx_0[j];

                tdy_yyz_xx_0[j] = pa_y[j] * tdy_yz_xx_0[j] + 0.5 * fl1_fx * tdy_z_xx_0[j] + 0.5 * fl1_fx * ts_yz_xx_0[j];

                tdz_yyz_xx_0[j] = pa_y[j] * tdz_yz_xx_0[j] + 0.5 * fl1_fx * tdz_z_xx_0[j];

                tdx_yyz_xy_0[j] = pa_y[j] * tdx_yz_xy_0[j] + 0.5 * fl1_fx * tdx_z_xy_0[j] + 0.5 * fl1_fx * tdx_yz_x_0[j];

                tdy_yyz_xy_0[j] = pa_y[j] * tdy_yz_xy_0[j] + 0.5 * fl1_fx * tdy_z_xy_0[j] + 0.5 * fl1_fx * tdy_yz_x_0[j] + 0.5 * fl1_fx * ts_yz_xy_0[j];

                tdz_yyz_xy_0[j] = pa_y[j] * tdz_yz_xy_0[j] + 0.5 * fl1_fx * tdz_z_xy_0[j] + 0.5 * fl1_fx * tdz_yz_x_0[j];

                tdx_yyz_xz_0[j] = pa_y[j] * tdx_yz_xz_0[j] + 0.5 * fl1_fx * tdx_z_xz_0[j];

                tdy_yyz_xz_0[j] = pa_y[j] * tdy_yz_xz_0[j] + 0.5 * fl1_fx * tdy_z_xz_0[j] + 0.5 * fl1_fx * ts_yz_xz_0[j];

                tdz_yyz_xz_0[j] = pa_y[j] * tdz_yz_xz_0[j] + 0.5 * fl1_fx * tdz_z_xz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForFD_135_180(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
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

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_3_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tdx_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 27); 

            auto tdy_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tdz_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tdx_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 28); 

            auto tdy_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tdz_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tdx_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 29); 

            auto tdy_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 29); 

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

            auto tdx_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 13); 

            auto tdy_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 13); 

            auto tdz_yz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 13); 

            auto tdx_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 14); 

            auto tdy_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 14); 

            auto tdz_yz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 14); 

            auto tdx_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 15); 

            auto tdy_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 15); 

            auto tdz_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 15); 

            auto tdx_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 16); 

            auto tdy_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 16); 

            auto tdz_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 16); 

            auto tdx_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 17); 

            auto tdy_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 17); 

            auto tdz_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 17); 

            auto ts_yz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 27); 

            auto ts_yz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 28); 

            auto ts_yz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 29); 

            auto ts_zz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 30); 

            auto ts_zz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 31); 

            auto ts_zz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 32); 

            auto ts_zz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 33); 

            auto ts_zz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 34); 

            auto ts_zz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 35); 

            // set up pointers to integrals

            auto tdx_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 45); 

            auto tdy_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto tdz_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdx_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 46); 

            auto tdy_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tdz_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdx_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 47); 

            auto tdy_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tdz_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdx_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 48); 

            auto tdy_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tdz_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdx_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 49); 

            auto tdy_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tdz_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdx_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 50); 

            auto tdy_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 51); 

            auto tdy_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 52); 

            auto tdy_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 53); 

            auto tdy_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 53); 

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

            // Batch of Integrals (135,180)

            #pragma omp simd aligned(fx, pa_y, pa_z, tdx_yyz_yy_0, tdx_yyz_yz_0, tdx_yyz_zz_0, tdx_yz_y_0, \
                                     tdx_yz_yy_0, tdx_yz_yz_0, tdx_yz_z_0, tdx_yz_zz_0, tdx_yzz_xx_0, tdx_yzz_xy_0, \
                                     tdx_yzz_xz_0, tdx_yzz_yy_0, tdx_yzz_yz_0, tdx_yzz_zz_0, tdx_z_xx_0, tdx_z_xy_0, \
                                     tdx_z_xz_0, tdx_z_yy_0, tdx_z_yz_0, tdx_z_zz_0, tdx_zz_x_0, tdx_zz_xx_0, \
                                     tdx_zz_xy_0, tdx_zz_xz_0, tdx_zz_y_0, tdx_zz_yy_0, tdx_zz_yz_0, tdx_zz_z_0, \
                                     tdx_zz_zz_0, tdx_zzz_xx_0, tdx_zzz_xy_0, tdx_zzz_xz_0, tdx_zzz_yy_0, tdx_zzz_yz_0, \
                                     tdx_zzz_zz_0, tdy_yyz_yy_0, tdy_yyz_yz_0, tdy_yyz_zz_0, tdy_yz_y_0, tdy_yz_yy_0, \
                                     tdy_yz_yz_0, tdy_yz_z_0, tdy_yz_zz_0, tdy_yzz_xx_0, tdy_yzz_xy_0, tdy_yzz_xz_0, \
                                     tdy_yzz_yy_0, tdy_yzz_yz_0, tdy_yzz_zz_0, tdy_z_xx_0, tdy_z_xy_0, tdy_z_xz_0, \
                                     tdy_z_yy_0, tdy_z_yz_0, tdy_z_zz_0, tdy_zz_x_0, tdy_zz_xx_0, tdy_zz_xy_0, \
                                     tdy_zz_xz_0, tdy_zz_y_0, tdy_zz_yy_0, tdy_zz_yz_0, tdy_zz_z_0, tdy_zz_zz_0, \
                                     tdy_zzz_xx_0, tdy_zzz_xy_0, tdy_zzz_xz_0, tdy_zzz_yy_0, tdy_zzz_yz_0, tdy_zzz_zz_0, \
                                     tdz_yyz_yy_0, tdz_yyz_yz_0, tdz_yyz_zz_0, tdz_yz_y_0, tdz_yz_yy_0, tdz_yz_yz_0, \
                                     tdz_yz_z_0, tdz_yz_zz_0, tdz_yzz_xx_0, tdz_yzz_xy_0, tdz_yzz_xz_0, tdz_yzz_yy_0, \
                                     tdz_yzz_yz_0, tdz_yzz_zz_0, tdz_z_xx_0, tdz_z_xy_0, tdz_z_xz_0, tdz_z_yy_0, \
                                     tdz_z_yz_0, tdz_z_zz_0, tdz_zz_x_0, tdz_zz_xx_0, tdz_zz_xy_0, tdz_zz_xz_0, \
                                     tdz_zz_y_0, tdz_zz_yy_0, tdz_zz_yz_0, tdz_zz_z_0, tdz_zz_zz_0, tdz_zzz_xx_0, \
                                     tdz_zzz_xy_0, tdz_zzz_xz_0, tdz_zzz_yy_0, tdz_zzz_yz_0, tdz_zzz_zz_0, ts_yz_yy_0, \
                                     ts_yz_yz_0, ts_yz_zz_0, ts_zz_xx_0, ts_zz_xy_0, ts_zz_xz_0, ts_zz_yy_0, ts_zz_yz_0, \
                                     ts_zz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yyz_yy_0[j] = pa_y[j] * tdx_yz_yy_0[j] + 0.5 * fl1_fx * tdx_z_yy_0[j] + fl1_fx * tdx_yz_y_0[j];

                tdy_yyz_yy_0[j] = pa_y[j] * tdy_yz_yy_0[j] + 0.5 * fl1_fx * tdy_z_yy_0[j] + fl1_fx * tdy_yz_y_0[j] + 0.5 * fl1_fx * ts_yz_yy_0[j];

                tdz_yyz_yy_0[j] = pa_y[j] * tdz_yz_yy_0[j] + 0.5 * fl1_fx * tdz_z_yy_0[j] + fl1_fx * tdz_yz_y_0[j];

                tdx_yyz_yz_0[j] = pa_y[j] * tdx_yz_yz_0[j] + 0.5 * fl1_fx * tdx_z_yz_0[j] + 0.5 * fl1_fx * tdx_yz_z_0[j];

                tdy_yyz_yz_0[j] = pa_y[j] * tdy_yz_yz_0[j] + 0.5 * fl1_fx * tdy_z_yz_0[j] + 0.5 * fl1_fx * tdy_yz_z_0[j] + 0.5 * fl1_fx * ts_yz_yz_0[j];

                tdz_yyz_yz_0[j] = pa_y[j] * tdz_yz_yz_0[j] + 0.5 * fl1_fx * tdz_z_yz_0[j] + 0.5 * fl1_fx * tdz_yz_z_0[j];

                tdx_yyz_zz_0[j] = pa_y[j] * tdx_yz_zz_0[j] + 0.5 * fl1_fx * tdx_z_zz_0[j];

                tdy_yyz_zz_0[j] = pa_y[j] * tdy_yz_zz_0[j] + 0.5 * fl1_fx * tdy_z_zz_0[j] + 0.5 * fl1_fx * ts_yz_zz_0[j];

                tdz_yyz_zz_0[j] = pa_y[j] * tdz_yz_zz_0[j] + 0.5 * fl1_fx * tdz_z_zz_0[j];

                tdx_yzz_xx_0[j] = pa_y[j] * tdx_zz_xx_0[j];

                tdy_yzz_xx_0[j] = pa_y[j] * tdy_zz_xx_0[j] + 0.5 * fl1_fx * ts_zz_xx_0[j];

                tdz_yzz_xx_0[j] = pa_y[j] * tdz_zz_xx_0[j];

                tdx_yzz_xy_0[j] = pa_y[j] * tdx_zz_xy_0[j] + 0.5 * fl1_fx * tdx_zz_x_0[j];

                tdy_yzz_xy_0[j] = pa_y[j] * tdy_zz_xy_0[j] + 0.5 * fl1_fx * tdy_zz_x_0[j] + 0.5 * fl1_fx * ts_zz_xy_0[j];

                tdz_yzz_xy_0[j] = pa_y[j] * tdz_zz_xy_0[j] + 0.5 * fl1_fx * tdz_zz_x_0[j];

                tdx_yzz_xz_0[j] = pa_y[j] * tdx_zz_xz_0[j];

                tdy_yzz_xz_0[j] = pa_y[j] * tdy_zz_xz_0[j] + 0.5 * fl1_fx * ts_zz_xz_0[j];

                tdz_yzz_xz_0[j] = pa_y[j] * tdz_zz_xz_0[j];

                tdx_yzz_yy_0[j] = pa_y[j] * tdx_zz_yy_0[j] + fl1_fx * tdx_zz_y_0[j];

                tdy_yzz_yy_0[j] = pa_y[j] * tdy_zz_yy_0[j] + fl1_fx * tdy_zz_y_0[j] + 0.5 * fl1_fx * ts_zz_yy_0[j];

                tdz_yzz_yy_0[j] = pa_y[j] * tdz_zz_yy_0[j] + fl1_fx * tdz_zz_y_0[j];

                tdx_yzz_yz_0[j] = pa_y[j] * tdx_zz_yz_0[j] + 0.5 * fl1_fx * tdx_zz_z_0[j];

                tdy_yzz_yz_0[j] = pa_y[j] * tdy_zz_yz_0[j] + 0.5 * fl1_fx * tdy_zz_z_0[j] + 0.5 * fl1_fx * ts_zz_yz_0[j];

                tdz_yzz_yz_0[j] = pa_y[j] * tdz_zz_yz_0[j] + 0.5 * fl1_fx * tdz_zz_z_0[j];

                tdx_yzz_zz_0[j] = pa_y[j] * tdx_zz_zz_0[j];

                tdy_yzz_zz_0[j] = pa_y[j] * tdy_zz_zz_0[j] + 0.5 * fl1_fx * ts_zz_zz_0[j];

                tdz_yzz_zz_0[j] = pa_y[j] * tdz_zz_zz_0[j];

                tdx_zzz_xx_0[j] = pa_z[j] * tdx_zz_xx_0[j] + fl1_fx * tdx_z_xx_0[j];

                tdy_zzz_xx_0[j] = pa_z[j] * tdy_zz_xx_0[j] + fl1_fx * tdy_z_xx_0[j];

                tdz_zzz_xx_0[j] = pa_z[j] * tdz_zz_xx_0[j] + fl1_fx * tdz_z_xx_0[j] + 0.5 * fl1_fx * ts_zz_xx_0[j];

                tdx_zzz_xy_0[j] = pa_z[j] * tdx_zz_xy_0[j] + fl1_fx * tdx_z_xy_0[j];

                tdy_zzz_xy_0[j] = pa_z[j] * tdy_zz_xy_0[j] + fl1_fx * tdy_z_xy_0[j];

                tdz_zzz_xy_0[j] = pa_z[j] * tdz_zz_xy_0[j] + fl1_fx * tdz_z_xy_0[j] + 0.5 * fl1_fx * ts_zz_xy_0[j];

                tdx_zzz_xz_0[j] = pa_z[j] * tdx_zz_xz_0[j] + fl1_fx * tdx_z_xz_0[j] + 0.5 * fl1_fx * tdx_zz_x_0[j];

                tdy_zzz_xz_0[j] = pa_z[j] * tdy_zz_xz_0[j] + fl1_fx * tdy_z_xz_0[j] + 0.5 * fl1_fx * tdy_zz_x_0[j];

                tdz_zzz_xz_0[j] = pa_z[j] * tdz_zz_xz_0[j] + fl1_fx * tdz_z_xz_0[j] + 0.5 * fl1_fx * tdz_zz_x_0[j] + 0.5 * fl1_fx * ts_zz_xz_0[j];

                tdx_zzz_yy_0[j] = pa_z[j] * tdx_zz_yy_0[j] + fl1_fx * tdx_z_yy_0[j];

                tdy_zzz_yy_0[j] = pa_z[j] * tdy_zz_yy_0[j] + fl1_fx * tdy_z_yy_0[j];

                tdz_zzz_yy_0[j] = pa_z[j] * tdz_zz_yy_0[j] + fl1_fx * tdz_z_yy_0[j] + 0.5 * fl1_fx * ts_zz_yy_0[j];

                tdx_zzz_yz_0[j] = pa_z[j] * tdx_zz_yz_0[j] + fl1_fx * tdx_z_yz_0[j] + 0.5 * fl1_fx * tdx_zz_y_0[j];

                tdy_zzz_yz_0[j] = pa_z[j] * tdy_zz_yz_0[j] + fl1_fx * tdy_z_yz_0[j] + 0.5 * fl1_fx * tdy_zz_y_0[j];

                tdz_zzz_yz_0[j] = pa_z[j] * tdz_zz_yz_0[j] + fl1_fx * tdz_z_yz_0[j] + 0.5 * fl1_fx * tdz_zz_y_0[j] + 0.5 * fl1_fx * ts_zz_yz_0[j];

                tdx_zzz_zz_0[j] = pa_z[j] * tdx_zz_zz_0[j] + fl1_fx * tdx_z_zz_0[j] + fl1_fx * tdx_zz_z_0[j];

                tdy_zzz_zz_0[j] = pa_z[j] * tdy_zz_zz_0[j] + fl1_fx * tdy_z_zz_0[j] + fl1_fx * tdy_zz_z_0[j];

                tdz_zzz_zz_0[j] = pa_z[j] * tdz_zz_zz_0[j] + fl1_fx * tdz_z_zz_0[j] + fl1_fx * tdz_zz_z_0[j] + 0.5 * fl1_fx * ts_zz_zz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDG(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForDG_0_45(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  nOSFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForDG_45_90(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   nOSFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        ediprecfunc::compElectricDipoleForDG_90_135(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    nOSFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        ediprecfunc::compElectricDipoleForDG_135_180(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     nOSFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForDG_180_225(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     nOSFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForDG_225_270(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     nOSFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 
    }

    void
    compElectricDipoleForDG_0_45(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const int32_t              nOSFactors,
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

        auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_x_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx); 

            auto tdy_x_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx); 

            auto tdz_x_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx); 

            auto tdx_x_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 1); 

            auto tdy_x_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 1); 

            auto tdz_x_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 1); 

            auto tdx_x_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 2); 

            auto tdy_x_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 2); 

            auto tdz_x_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 2); 

            auto tdx_x_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 3); 

            auto tdy_x_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 3); 

            auto tdz_x_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 3); 

            auto tdx_x_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 4); 

            auto tdy_x_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 4); 

            auto tdz_x_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 4); 

            auto tdx_x_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 5); 

            auto tdy_x_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 5); 

            auto tdz_x_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 5); 

            auto tdx_x_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 6); 

            auto tdy_x_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 6); 

            auto tdz_x_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 6); 

            auto tdx_x_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 7); 

            auto tdy_x_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 7); 

            auto tdz_x_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 7); 

            auto tdx_x_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 8); 

            auto tdy_x_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 8); 

            auto tdz_x_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 8); 

            auto tdx_x_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 9); 

            auto tdy_x_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 9); 

            auto tdz_x_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 9); 

            auto tdx_x_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 10); 

            auto tdy_x_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 10); 

            auto tdz_x_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 10); 

            auto tdx_x_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 11); 

            auto tdy_x_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 11); 

            auto tdz_x_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 11); 

            auto tdx_x_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 12); 

            auto tdy_x_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 12); 

            auto tdz_x_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 12); 

            auto tdx_x_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 13); 

            auto tdy_x_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 13); 

            auto tdz_x_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 13); 

            auto tdx_x_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 14); 

            auto tdy_x_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 14); 

            auto tdz_x_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 14); 

            auto tdx_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx); 

            auto tdy_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx); 

            auto tdz_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx); 

            auto tdx_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 1); 

            auto tdy_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 1); 

            auto tdz_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 1); 

            auto tdx_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 2); 

            auto tdy_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 2); 

            auto tdz_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 2); 

            auto tdx_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 3); 

            auto tdy_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 3); 

            auto tdz_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 3); 

            auto tdx_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 4); 

            auto tdy_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 4); 

            auto tdz_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 4); 

            auto tdx_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 5); 

            auto tdy_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 5); 

            auto tdz_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 5); 

            auto tdx_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 6); 

            auto tdy_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 6); 

            auto tdz_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 6); 

            auto tdx_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 7); 

            auto tdy_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 7); 

            auto tdz_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 7); 

            auto tdx_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 8); 

            auto tdy_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 8); 

            auto tdz_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 8); 

            auto tdx_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 9); 

            auto tdy_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 9); 

            auto tdz_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 9); 

            auto tdx_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 10); 

            auto tdy_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 10); 

            auto tdz_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 10); 

            auto tdx_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 11); 

            auto tdy_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 11); 

            auto tdz_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 11); 

            auto tdx_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 12); 

            auto tdy_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 12); 

            auto tdz_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 12); 

            auto tdx_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 13); 

            auto tdy_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 13); 

            auto tdz_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 13); 

            auto tdx_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 14); 

            auto tdy_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 14); 

            auto tdz_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 14); 

            auto tdx_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx); 

            auto tdy_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx); 

            auto tdz_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx); 

            auto tdx_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 1); 

            auto tdy_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 1); 

            auto tdz_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 1); 

            auto tdx_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 2); 

            auto tdy_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 2); 

            auto tdz_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 2); 

            auto tdx_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 3); 

            auto tdy_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 3); 

            auto tdz_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 3); 

            auto tdx_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 4); 

            auto tdy_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 4); 

            auto tdz_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 4); 

            auto tdx_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 5); 

            auto tdy_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 5); 

            auto tdz_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 5); 

            auto tdx_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 6); 

            auto tdy_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 6); 

            auto tdz_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 6); 

            auto tdx_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 7); 

            auto tdy_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 7); 

            auto tdz_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 7); 

            auto tdx_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 8); 

            auto tdy_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 8); 

            auto tdz_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 8); 

            auto tdx_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 9); 

            auto tdy_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 9); 

            auto tdz_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 9); 

            auto ts_x_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx); 

            auto ts_x_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 1); 

            auto ts_x_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 2); 

            auto ts_x_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 3); 

            auto ts_x_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 4); 

            auto ts_x_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 5); 

            auto ts_x_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 6); 

            auto ts_x_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 7); 

            auto ts_x_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 8); 

            auto ts_x_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 9); 

            auto ts_x_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 10); 

            auto ts_x_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 11); 

            auto ts_x_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 12); 

            auto ts_x_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 13); 

            auto ts_x_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 14); 

            // set up pointers to integrals

            auto tdx_xx_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx); 

            auto tdy_xx_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx); 

            auto tdz_xx_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx); 

            auto tdx_xx_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 1); 

            auto tdy_xx_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 1); 

            auto tdz_xx_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 1); 

            auto tdx_xx_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 2); 

            auto tdy_xx_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 2); 

            auto tdz_xx_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 2); 

            auto tdx_xx_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 3); 

            auto tdy_xx_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 3); 

            auto tdz_xx_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 3); 

            auto tdx_xx_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 4); 

            auto tdy_xx_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 4); 

            auto tdz_xx_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 4); 

            auto tdx_xx_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 5); 

            auto tdy_xx_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 5); 

            auto tdz_xx_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 5); 

            auto tdx_xx_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 6); 

            auto tdy_xx_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 6); 

            auto tdz_xx_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 6); 

            auto tdx_xx_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 7); 

            auto tdy_xx_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 7); 

            auto tdz_xx_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 7); 

            auto tdx_xx_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 8); 

            auto tdy_xx_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 8); 

            auto tdz_xx_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 8); 

            auto tdx_xx_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 9); 

            auto tdy_xx_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 9); 

            auto tdz_xx_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 9); 

            auto tdx_xx_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 10); 

            auto tdy_xx_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 10); 

            auto tdz_xx_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 10); 

            auto tdx_xx_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 11); 

            auto tdy_xx_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 11); 

            auto tdz_xx_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 11); 

            auto tdx_xx_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 12); 

            auto tdy_xx_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 12); 

            auto tdz_xx_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 12); 

            auto tdx_xx_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 13); 

            auto tdy_xx_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 13); 

            auto tdz_xx_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 13); 

            auto tdx_xx_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 14); 

            auto tdy_xx_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 14); 

            auto tdz_xx_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 14); 

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, tdx_0_xxxx_0, tdx_0_xxxy_0, tdx_0_xxxz_0, tdx_0_xxyy_0, \
                                     tdx_0_xxyz_0, tdx_0_xxzz_0, tdx_0_xyyy_0, tdx_0_xyyz_0, tdx_0_xyzz_0, tdx_0_xzzz_0, \
                                     tdx_0_yyyy_0, tdx_0_yyyz_0, tdx_0_yyzz_0, tdx_0_yzzz_0, tdx_0_zzzz_0, tdx_x_xxx_0, \
                                     tdx_x_xxxx_0, tdx_x_xxxy_0, tdx_x_xxxz_0, tdx_x_xxy_0, tdx_x_xxyy_0, tdx_x_xxyz_0, \
                                     tdx_x_xxz_0, tdx_x_xxzz_0, tdx_x_xyy_0, tdx_x_xyyy_0, tdx_x_xyyz_0, tdx_x_xyz_0, \
                                     tdx_x_xyzz_0, tdx_x_xzz_0, tdx_x_xzzz_0, tdx_x_yyy_0, tdx_x_yyyy_0, tdx_x_yyyz_0, \
                                     tdx_x_yyz_0, tdx_x_yyzz_0, tdx_x_yzz_0, tdx_x_yzzz_0, tdx_x_zzz_0, tdx_x_zzzz_0, \
                                     tdx_xx_xxxx_0, tdx_xx_xxxy_0, tdx_xx_xxxz_0, tdx_xx_xxyy_0, tdx_xx_xxyz_0, \
                                     tdx_xx_xxzz_0, tdx_xx_xyyy_0, tdx_xx_xyyz_0, tdx_xx_xyzz_0, tdx_xx_xzzz_0, \
                                     tdx_xx_yyyy_0, tdx_xx_yyyz_0, tdx_xx_yyzz_0, tdx_xx_yzzz_0, tdx_xx_zzzz_0, \
                                     tdy_0_xxxx_0, tdy_0_xxxy_0, tdy_0_xxxz_0, tdy_0_xxyy_0, tdy_0_xxyz_0, tdy_0_xxzz_0, \
                                     tdy_0_xyyy_0, tdy_0_xyyz_0, tdy_0_xyzz_0, tdy_0_xzzz_0, tdy_0_yyyy_0, tdy_0_yyyz_0, \
                                     tdy_0_yyzz_0, tdy_0_yzzz_0, tdy_0_zzzz_0, tdy_x_xxx_0, tdy_x_xxxx_0, tdy_x_xxxy_0, \
                                     tdy_x_xxxz_0, tdy_x_xxy_0, tdy_x_xxyy_0, tdy_x_xxyz_0, tdy_x_xxz_0, tdy_x_xxzz_0, \
                                     tdy_x_xyy_0, tdy_x_xyyy_0, tdy_x_xyyz_0, tdy_x_xyz_0, tdy_x_xyzz_0, tdy_x_xzz_0, \
                                     tdy_x_xzzz_0, tdy_x_yyy_0, tdy_x_yyyy_0, tdy_x_yyyz_0, tdy_x_yyz_0, tdy_x_yyzz_0, \
                                     tdy_x_yzz_0, tdy_x_yzzz_0, tdy_x_zzz_0, tdy_x_zzzz_0, tdy_xx_xxxx_0, \
                                     tdy_xx_xxxy_0, tdy_xx_xxxz_0, tdy_xx_xxyy_0, tdy_xx_xxyz_0, tdy_xx_xxzz_0, \
                                     tdy_xx_xyyy_0, tdy_xx_xyyz_0, tdy_xx_xyzz_0, tdy_xx_xzzz_0, tdy_xx_yyyy_0, \
                                     tdy_xx_yyyz_0, tdy_xx_yyzz_0, tdy_xx_yzzz_0, tdy_xx_zzzz_0, tdz_0_xxxx_0, \
                                     tdz_0_xxxy_0, tdz_0_xxxz_0, tdz_0_xxyy_0, tdz_0_xxyz_0, tdz_0_xxzz_0, tdz_0_xyyy_0, \
                                     tdz_0_xyyz_0, tdz_0_xyzz_0, tdz_0_xzzz_0, tdz_0_yyyy_0, tdz_0_yyyz_0, tdz_0_yyzz_0, \
                                     tdz_0_yzzz_0, tdz_0_zzzz_0, tdz_x_xxx_0, tdz_x_xxxx_0, tdz_x_xxxy_0, tdz_x_xxxz_0, \
                                     tdz_x_xxy_0, tdz_x_xxyy_0, tdz_x_xxyz_0, tdz_x_xxz_0, tdz_x_xxzz_0, tdz_x_xyy_0, \
                                     tdz_x_xyyy_0, tdz_x_xyyz_0, tdz_x_xyz_0, tdz_x_xyzz_0, tdz_x_xzz_0, tdz_x_xzzz_0, \
                                     tdz_x_yyy_0, tdz_x_yyyy_0, tdz_x_yyyz_0, tdz_x_yyz_0, tdz_x_yyzz_0, tdz_x_yzz_0, \
                                     tdz_x_yzzz_0, tdz_x_zzz_0, tdz_x_zzzz_0, tdz_xx_xxxx_0, tdz_xx_xxxy_0, \
                                     tdz_xx_xxxz_0, tdz_xx_xxyy_0, tdz_xx_xxyz_0, tdz_xx_xxzz_0, tdz_xx_xyyy_0, \
                                     tdz_xx_xyyz_0, tdz_xx_xyzz_0, tdz_xx_xzzz_0, tdz_xx_yyyy_0, tdz_xx_yyyz_0, \
                                     tdz_xx_yyzz_0, tdz_xx_yzzz_0, tdz_xx_zzzz_0, ts_x_xxxx_0, ts_x_xxxy_0, ts_x_xxxz_0, \
                                     ts_x_xxyy_0, ts_x_xxyz_0, ts_x_xxzz_0, ts_x_xyyy_0, ts_x_xyyz_0, ts_x_xyzz_0, \
                                     ts_x_xzzz_0, ts_x_yyyy_0, ts_x_yyyz_0, ts_x_yyzz_0, ts_x_yzzz_0, ts_x_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xx_xxxx_0[j] = pa_x[j] * tdx_x_xxxx_0[j] + 0.5 * fl1_fx * tdx_0_xxxx_0[j] + 2.0 * fl1_fx * tdx_x_xxx_0[j] + 0.5 * fl1_fx * ts_x_xxxx_0[j];

                tdy_xx_xxxx_0[j] = pa_x[j] * tdy_x_xxxx_0[j] + 0.5 * fl1_fx * tdy_0_xxxx_0[j] + 2.0 * fl1_fx * tdy_x_xxx_0[j];

                tdz_xx_xxxx_0[j] = pa_x[j] * tdz_x_xxxx_0[j] + 0.5 * fl1_fx * tdz_0_xxxx_0[j] + 2.0 * fl1_fx * tdz_x_xxx_0[j];

                tdx_xx_xxxy_0[j] = pa_x[j] * tdx_x_xxxy_0[j] + 0.5 * fl1_fx * tdx_0_xxxy_0[j] + 1.5 * fl1_fx * tdx_x_xxy_0[j] + 0.5 * fl1_fx * ts_x_xxxy_0[j];

                tdy_xx_xxxy_0[j] = pa_x[j] * tdy_x_xxxy_0[j] + 0.5 * fl1_fx * tdy_0_xxxy_0[j] + 1.5 * fl1_fx * tdy_x_xxy_0[j];

                tdz_xx_xxxy_0[j] = pa_x[j] * tdz_x_xxxy_0[j] + 0.5 * fl1_fx * tdz_0_xxxy_0[j] + 1.5 * fl1_fx * tdz_x_xxy_0[j];

                tdx_xx_xxxz_0[j] = pa_x[j] * tdx_x_xxxz_0[j] + 0.5 * fl1_fx * tdx_0_xxxz_0[j] + 1.5 * fl1_fx * tdx_x_xxz_0[j] + 0.5 * fl1_fx * ts_x_xxxz_0[j];

                tdy_xx_xxxz_0[j] = pa_x[j] * tdy_x_xxxz_0[j] + 0.5 * fl1_fx * tdy_0_xxxz_0[j] + 1.5 * fl1_fx * tdy_x_xxz_0[j];

                tdz_xx_xxxz_0[j] = pa_x[j] * tdz_x_xxxz_0[j] + 0.5 * fl1_fx * tdz_0_xxxz_0[j] + 1.5 * fl1_fx * tdz_x_xxz_0[j];

                tdx_xx_xxyy_0[j] = pa_x[j] * tdx_x_xxyy_0[j] + 0.5 * fl1_fx * tdx_0_xxyy_0[j] + fl1_fx * tdx_x_xyy_0[j] + 0.5 * fl1_fx * ts_x_xxyy_0[j];

                tdy_xx_xxyy_0[j] = pa_x[j] * tdy_x_xxyy_0[j] + 0.5 * fl1_fx * tdy_0_xxyy_0[j] + fl1_fx * tdy_x_xyy_0[j];

                tdz_xx_xxyy_0[j] = pa_x[j] * tdz_x_xxyy_0[j] + 0.5 * fl1_fx * tdz_0_xxyy_0[j] + fl1_fx * tdz_x_xyy_0[j];

                tdx_xx_xxyz_0[j] = pa_x[j] * tdx_x_xxyz_0[j] + 0.5 * fl1_fx * tdx_0_xxyz_0[j] + fl1_fx * tdx_x_xyz_0[j] + 0.5 * fl1_fx * ts_x_xxyz_0[j];

                tdy_xx_xxyz_0[j] = pa_x[j] * tdy_x_xxyz_0[j] + 0.5 * fl1_fx * tdy_0_xxyz_0[j] + fl1_fx * tdy_x_xyz_0[j];

                tdz_xx_xxyz_0[j] = pa_x[j] * tdz_x_xxyz_0[j] + 0.5 * fl1_fx * tdz_0_xxyz_0[j] + fl1_fx * tdz_x_xyz_0[j];

                tdx_xx_xxzz_0[j] = pa_x[j] * tdx_x_xxzz_0[j] + 0.5 * fl1_fx * tdx_0_xxzz_0[j] + fl1_fx * tdx_x_xzz_0[j] + 0.5 * fl1_fx * ts_x_xxzz_0[j];

                tdy_xx_xxzz_0[j] = pa_x[j] * tdy_x_xxzz_0[j] + 0.5 * fl1_fx * tdy_0_xxzz_0[j] + fl1_fx * tdy_x_xzz_0[j];

                tdz_xx_xxzz_0[j] = pa_x[j] * tdz_x_xxzz_0[j] + 0.5 * fl1_fx * tdz_0_xxzz_0[j] + fl1_fx * tdz_x_xzz_0[j];

                tdx_xx_xyyy_0[j] = pa_x[j] * tdx_x_xyyy_0[j] + 0.5 * fl1_fx * tdx_0_xyyy_0[j] + 0.5 * fl1_fx * tdx_x_yyy_0[j] + 0.5 * fl1_fx * ts_x_xyyy_0[j];

                tdy_xx_xyyy_0[j] = pa_x[j] * tdy_x_xyyy_0[j] + 0.5 * fl1_fx * tdy_0_xyyy_0[j] + 0.5 * fl1_fx * tdy_x_yyy_0[j];

                tdz_xx_xyyy_0[j] = pa_x[j] * tdz_x_xyyy_0[j] + 0.5 * fl1_fx * tdz_0_xyyy_0[j] + 0.5 * fl1_fx * tdz_x_yyy_0[j];

                tdx_xx_xyyz_0[j] = pa_x[j] * tdx_x_xyyz_0[j] + 0.5 * fl1_fx * tdx_0_xyyz_0[j] + 0.5 * fl1_fx * tdx_x_yyz_0[j] + 0.5 * fl1_fx * ts_x_xyyz_0[j];

                tdy_xx_xyyz_0[j] = pa_x[j] * tdy_x_xyyz_0[j] + 0.5 * fl1_fx * tdy_0_xyyz_0[j] + 0.5 * fl1_fx * tdy_x_yyz_0[j];

                tdz_xx_xyyz_0[j] = pa_x[j] * tdz_x_xyyz_0[j] + 0.5 * fl1_fx * tdz_0_xyyz_0[j] + 0.5 * fl1_fx * tdz_x_yyz_0[j];

                tdx_xx_xyzz_0[j] = pa_x[j] * tdx_x_xyzz_0[j] + 0.5 * fl1_fx * tdx_0_xyzz_0[j] + 0.5 * fl1_fx * tdx_x_yzz_0[j] + 0.5 * fl1_fx * ts_x_xyzz_0[j];

                tdy_xx_xyzz_0[j] = pa_x[j] * tdy_x_xyzz_0[j] + 0.5 * fl1_fx * tdy_0_xyzz_0[j] + 0.5 * fl1_fx * tdy_x_yzz_0[j];

                tdz_xx_xyzz_0[j] = pa_x[j] * tdz_x_xyzz_0[j] + 0.5 * fl1_fx * tdz_0_xyzz_0[j] + 0.5 * fl1_fx * tdz_x_yzz_0[j];

                tdx_xx_xzzz_0[j] = pa_x[j] * tdx_x_xzzz_0[j] + 0.5 * fl1_fx * tdx_0_xzzz_0[j] + 0.5 * fl1_fx * tdx_x_zzz_0[j] + 0.5 * fl1_fx * ts_x_xzzz_0[j];

                tdy_xx_xzzz_0[j] = pa_x[j] * tdy_x_xzzz_0[j] + 0.5 * fl1_fx * tdy_0_xzzz_0[j] + 0.5 * fl1_fx * tdy_x_zzz_0[j];

                tdz_xx_xzzz_0[j] = pa_x[j] * tdz_x_xzzz_0[j] + 0.5 * fl1_fx * tdz_0_xzzz_0[j] + 0.5 * fl1_fx * tdz_x_zzz_0[j];

                tdx_xx_yyyy_0[j] = pa_x[j] * tdx_x_yyyy_0[j] + 0.5 * fl1_fx * tdx_0_yyyy_0[j] + 0.5 * fl1_fx * ts_x_yyyy_0[j];

                tdy_xx_yyyy_0[j] = pa_x[j] * tdy_x_yyyy_0[j] + 0.5 * fl1_fx * tdy_0_yyyy_0[j];

                tdz_xx_yyyy_0[j] = pa_x[j] * tdz_x_yyyy_0[j] + 0.5 * fl1_fx * tdz_0_yyyy_0[j];

                tdx_xx_yyyz_0[j] = pa_x[j] * tdx_x_yyyz_0[j] + 0.5 * fl1_fx * tdx_0_yyyz_0[j] + 0.5 * fl1_fx * ts_x_yyyz_0[j];

                tdy_xx_yyyz_0[j] = pa_x[j] * tdy_x_yyyz_0[j] + 0.5 * fl1_fx * tdy_0_yyyz_0[j];

                tdz_xx_yyyz_0[j] = pa_x[j] * tdz_x_yyyz_0[j] + 0.5 * fl1_fx * tdz_0_yyyz_0[j];

                tdx_xx_yyzz_0[j] = pa_x[j] * tdx_x_yyzz_0[j] + 0.5 * fl1_fx * tdx_0_yyzz_0[j] + 0.5 * fl1_fx * ts_x_yyzz_0[j];

                tdy_xx_yyzz_0[j] = pa_x[j] * tdy_x_yyzz_0[j] + 0.5 * fl1_fx * tdy_0_yyzz_0[j];

                tdz_xx_yyzz_0[j] = pa_x[j] * tdz_x_yyzz_0[j] + 0.5 * fl1_fx * tdz_0_yyzz_0[j];

                tdx_xx_yzzz_0[j] = pa_x[j] * tdx_x_yzzz_0[j] + 0.5 * fl1_fx * tdx_0_yzzz_0[j] + 0.5 * fl1_fx * ts_x_yzzz_0[j];

                tdy_xx_yzzz_0[j] = pa_x[j] * tdy_x_yzzz_0[j] + 0.5 * fl1_fx * tdy_0_yzzz_0[j];

                tdz_xx_yzzz_0[j] = pa_x[j] * tdz_x_yzzz_0[j] + 0.5 * fl1_fx * tdz_0_yzzz_0[j];

                tdx_xx_zzzz_0[j] = pa_x[j] * tdx_x_zzzz_0[j] + 0.5 * fl1_fx * tdx_0_zzzz_0[j] + 0.5 * fl1_fx * ts_x_zzzz_0[j];

                tdy_xx_zzzz_0[j] = pa_x[j] * tdy_x_zzzz_0[j] + 0.5 * fl1_fx * tdy_0_zzzz_0[j];

                tdz_xx_zzzz_0[j] = pa_x[j] * tdz_x_zzzz_0[j] + 0.5 * fl1_fx * tdz_0_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDG_45_90(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
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

        auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 15); 

            auto tdy_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 15); 

            auto tdz_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 15); 

            auto tdx_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 16); 

            auto tdy_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 16); 

            auto tdz_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 16); 

            auto tdx_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 17); 

            auto tdy_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 17); 

            auto tdz_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 17); 

            auto tdx_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 18); 

            auto tdy_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 18); 

            auto tdz_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 18); 

            auto tdx_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 19); 

            auto tdy_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 19); 

            auto tdz_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 19); 

            auto tdx_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 20); 

            auto tdy_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 20); 

            auto tdz_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 20); 

            auto tdx_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 21); 

            auto tdy_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 21); 

            auto tdz_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 21); 

            auto tdx_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 22); 

            auto tdy_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 22); 

            auto tdz_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 22); 

            auto tdx_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 23); 

            auto tdy_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 23); 

            auto tdz_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 23); 

            auto tdx_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 24); 

            auto tdy_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 24); 

            auto tdz_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 24); 

            auto tdx_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 25); 

            auto tdy_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 25); 

            auto tdz_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 25); 

            auto tdx_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 26); 

            auto tdy_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 26); 

            auto tdz_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 26); 

            auto tdx_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 27); 

            auto tdy_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 27); 

            auto tdz_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 27); 

            auto tdx_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 28); 

            auto tdy_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 28); 

            auto tdz_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 28); 

            auto tdx_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 29); 

            auto tdy_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 29); 

            auto tdz_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 29); 

            auto tdx_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 10); 

            auto tdy_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tdz_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tdx_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 11); 

            auto tdy_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tdz_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tdx_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 12); 

            auto tdy_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tdz_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tdx_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 13); 

            auto tdy_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tdz_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tdx_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 14); 

            auto tdy_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tdz_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tdx_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 15); 

            auto tdy_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tdz_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tdx_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 16); 

            auto tdy_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tdz_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tdx_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 17); 

            auto tdy_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tdz_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tdx_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 18); 

            auto tdy_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tdz_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tdx_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 19); 

            auto tdy_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tdz_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto ts_y_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 15); 

            auto ts_y_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 16); 

            auto ts_y_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 17); 

            auto ts_y_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 18); 

            auto ts_y_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 19); 

            auto ts_y_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 20); 

            auto ts_y_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 21); 

            auto ts_y_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 22); 

            auto ts_y_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 23); 

            auto ts_y_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 24); 

            auto ts_y_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 25); 

            auto ts_y_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 26); 

            auto ts_y_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 27); 

            auto ts_y_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 28); 

            auto ts_y_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 29); 

            // set up pointers to integrals

            auto tdx_xy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 15); 

            auto tdy_xy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 15); 

            auto tdz_xy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 15); 

            auto tdx_xy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 16); 

            auto tdy_xy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 16); 

            auto tdz_xy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 16); 

            auto tdx_xy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 17); 

            auto tdy_xy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 17); 

            auto tdz_xy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 17); 

            auto tdx_xy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 18); 

            auto tdy_xy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 18); 

            auto tdz_xy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 18); 

            auto tdx_xy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 19); 

            auto tdy_xy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 19); 

            auto tdz_xy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 19); 

            auto tdx_xy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 20); 

            auto tdy_xy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 20); 

            auto tdz_xy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 20); 

            auto tdx_xy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 21); 

            auto tdy_xy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 21); 

            auto tdz_xy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 21); 

            auto tdx_xy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 22); 

            auto tdy_xy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 22); 

            auto tdz_xy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 22); 

            auto tdx_xy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 23); 

            auto tdy_xy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 23); 

            auto tdz_xy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 23); 

            auto tdx_xy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 24); 

            auto tdy_xy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 24); 

            auto tdz_xy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 24); 

            auto tdx_xy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 25); 

            auto tdy_xy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 25); 

            auto tdz_xy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 25); 

            auto tdx_xy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 26); 

            auto tdy_xy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 26); 

            auto tdz_xy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 26); 

            auto tdx_xy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 27); 

            auto tdy_xy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 27); 

            auto tdz_xy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 27); 

            auto tdx_xy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 28); 

            auto tdy_xy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 28); 

            auto tdz_xy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 28); 

            auto tdx_xy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 29); 

            auto tdy_xy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 29); 

            auto tdz_xy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 29); 

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, tdx_xy_xxxx_0, tdx_xy_xxxy_0, tdx_xy_xxxz_0, tdx_xy_xxyy_0, \
                                     tdx_xy_xxyz_0, tdx_xy_xxzz_0, tdx_xy_xyyy_0, tdx_xy_xyyz_0, tdx_xy_xyzz_0, \
                                     tdx_xy_xzzz_0, tdx_xy_yyyy_0, tdx_xy_yyyz_0, tdx_xy_yyzz_0, tdx_xy_yzzz_0, \
                                     tdx_xy_zzzz_0, tdx_y_xxx_0, tdx_y_xxxx_0, tdx_y_xxxy_0, tdx_y_xxxz_0, tdx_y_xxy_0, \
                                     tdx_y_xxyy_0, tdx_y_xxyz_0, tdx_y_xxz_0, tdx_y_xxzz_0, tdx_y_xyy_0, tdx_y_xyyy_0, \
                                     tdx_y_xyyz_0, tdx_y_xyz_0, tdx_y_xyzz_0, tdx_y_xzz_0, tdx_y_xzzz_0, tdx_y_yyy_0, \
                                     tdx_y_yyyy_0, tdx_y_yyyz_0, tdx_y_yyz_0, tdx_y_yyzz_0, tdx_y_yzz_0, tdx_y_yzzz_0, \
                                     tdx_y_zzz_0, tdx_y_zzzz_0, tdy_xy_xxxx_0, tdy_xy_xxxy_0, tdy_xy_xxxz_0, \
                                     tdy_xy_xxyy_0, tdy_xy_xxyz_0, tdy_xy_xxzz_0, tdy_xy_xyyy_0, tdy_xy_xyyz_0, \
                                     tdy_xy_xyzz_0, tdy_xy_xzzz_0, tdy_xy_yyyy_0, tdy_xy_yyyz_0, tdy_xy_yyzz_0, \
                                     tdy_xy_yzzz_0, tdy_xy_zzzz_0, tdy_y_xxx_0, tdy_y_xxxx_0, tdy_y_xxxy_0, tdy_y_xxxz_0, \
                                     tdy_y_xxy_0, tdy_y_xxyy_0, tdy_y_xxyz_0, tdy_y_xxz_0, tdy_y_xxzz_0, tdy_y_xyy_0, \
                                     tdy_y_xyyy_0, tdy_y_xyyz_0, tdy_y_xyz_0, tdy_y_xyzz_0, tdy_y_xzz_0, tdy_y_xzzz_0, \
                                     tdy_y_yyy_0, tdy_y_yyyy_0, tdy_y_yyyz_0, tdy_y_yyz_0, tdy_y_yyzz_0, tdy_y_yzz_0, \
                                     tdy_y_yzzz_0, tdy_y_zzz_0, tdy_y_zzzz_0, tdz_xy_xxxx_0, tdz_xy_xxxy_0, \
                                     tdz_xy_xxxz_0, tdz_xy_xxyy_0, tdz_xy_xxyz_0, tdz_xy_xxzz_0, tdz_xy_xyyy_0, \
                                     tdz_xy_xyyz_0, tdz_xy_xyzz_0, tdz_xy_xzzz_0, tdz_xy_yyyy_0, tdz_xy_yyyz_0, \
                                     tdz_xy_yyzz_0, tdz_xy_yzzz_0, tdz_xy_zzzz_0, tdz_y_xxx_0, tdz_y_xxxx_0, \
                                     tdz_y_xxxy_0, tdz_y_xxxz_0, tdz_y_xxy_0, tdz_y_xxyy_0, tdz_y_xxyz_0, tdz_y_xxz_0, \
                                     tdz_y_xxzz_0, tdz_y_xyy_0, tdz_y_xyyy_0, tdz_y_xyyz_0, tdz_y_xyz_0, tdz_y_xyzz_0, \
                                     tdz_y_xzz_0, tdz_y_xzzz_0, tdz_y_yyy_0, tdz_y_yyyy_0, tdz_y_yyyz_0, tdz_y_yyz_0, \
                                     tdz_y_yyzz_0, tdz_y_yzz_0, tdz_y_yzzz_0, tdz_y_zzz_0, tdz_y_zzzz_0, ts_y_xxxx_0, \
                                     ts_y_xxxy_0, ts_y_xxxz_0, ts_y_xxyy_0, ts_y_xxyz_0, ts_y_xxzz_0, ts_y_xyyy_0, \
                                     ts_y_xyyz_0, ts_y_xyzz_0, ts_y_xzzz_0, ts_y_yyyy_0, ts_y_yyyz_0, ts_y_yyzz_0, \
                                     ts_y_yzzz_0, ts_y_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xy_xxxx_0[j] = pa_x[j] * tdx_y_xxxx_0[j] + 2.0 * fl1_fx * tdx_y_xxx_0[j] + 0.5 * fl1_fx * ts_y_xxxx_0[j];

                tdy_xy_xxxx_0[j] = pa_x[j] * tdy_y_xxxx_0[j] + 2.0 * fl1_fx * tdy_y_xxx_0[j];

                tdz_xy_xxxx_0[j] = pa_x[j] * tdz_y_xxxx_0[j] + 2.0 * fl1_fx * tdz_y_xxx_0[j];

                tdx_xy_xxxy_0[j] = pa_x[j] * tdx_y_xxxy_0[j] + 1.5 * fl1_fx * tdx_y_xxy_0[j] + 0.5 * fl1_fx * ts_y_xxxy_0[j];

                tdy_xy_xxxy_0[j] = pa_x[j] * tdy_y_xxxy_0[j] + 1.5 * fl1_fx * tdy_y_xxy_0[j];

                tdz_xy_xxxy_0[j] = pa_x[j] * tdz_y_xxxy_0[j] + 1.5 * fl1_fx * tdz_y_xxy_0[j];

                tdx_xy_xxxz_0[j] = pa_x[j] * tdx_y_xxxz_0[j] + 1.5 * fl1_fx * tdx_y_xxz_0[j] + 0.5 * fl1_fx * ts_y_xxxz_0[j];

                tdy_xy_xxxz_0[j] = pa_x[j] * tdy_y_xxxz_0[j] + 1.5 * fl1_fx * tdy_y_xxz_0[j];

                tdz_xy_xxxz_0[j] = pa_x[j] * tdz_y_xxxz_0[j] + 1.5 * fl1_fx * tdz_y_xxz_0[j];

                tdx_xy_xxyy_0[j] = pa_x[j] * tdx_y_xxyy_0[j] + fl1_fx * tdx_y_xyy_0[j] + 0.5 * fl1_fx * ts_y_xxyy_0[j];

                tdy_xy_xxyy_0[j] = pa_x[j] * tdy_y_xxyy_0[j] + fl1_fx * tdy_y_xyy_0[j];

                tdz_xy_xxyy_0[j] = pa_x[j] * tdz_y_xxyy_0[j] + fl1_fx * tdz_y_xyy_0[j];

                tdx_xy_xxyz_0[j] = pa_x[j] * tdx_y_xxyz_0[j] + fl1_fx * tdx_y_xyz_0[j] + 0.5 * fl1_fx * ts_y_xxyz_0[j];

                tdy_xy_xxyz_0[j] = pa_x[j] * tdy_y_xxyz_0[j] + fl1_fx * tdy_y_xyz_0[j];

                tdz_xy_xxyz_0[j] = pa_x[j] * tdz_y_xxyz_0[j] + fl1_fx * tdz_y_xyz_0[j];

                tdx_xy_xxzz_0[j] = pa_x[j] * tdx_y_xxzz_0[j] + fl1_fx * tdx_y_xzz_0[j] + 0.5 * fl1_fx * ts_y_xxzz_0[j];

                tdy_xy_xxzz_0[j] = pa_x[j] * tdy_y_xxzz_0[j] + fl1_fx * tdy_y_xzz_0[j];

                tdz_xy_xxzz_0[j] = pa_x[j] * tdz_y_xxzz_0[j] + fl1_fx * tdz_y_xzz_0[j];

                tdx_xy_xyyy_0[j] = pa_x[j] * tdx_y_xyyy_0[j] + 0.5 * fl1_fx * tdx_y_yyy_0[j] + 0.5 * fl1_fx * ts_y_xyyy_0[j];

                tdy_xy_xyyy_0[j] = pa_x[j] * tdy_y_xyyy_0[j] + 0.5 * fl1_fx * tdy_y_yyy_0[j];

                tdz_xy_xyyy_0[j] = pa_x[j] * tdz_y_xyyy_0[j] + 0.5 * fl1_fx * tdz_y_yyy_0[j];

                tdx_xy_xyyz_0[j] = pa_x[j] * tdx_y_xyyz_0[j] + 0.5 * fl1_fx * tdx_y_yyz_0[j] + 0.5 * fl1_fx * ts_y_xyyz_0[j];

                tdy_xy_xyyz_0[j] = pa_x[j] * tdy_y_xyyz_0[j] + 0.5 * fl1_fx * tdy_y_yyz_0[j];

                tdz_xy_xyyz_0[j] = pa_x[j] * tdz_y_xyyz_0[j] + 0.5 * fl1_fx * tdz_y_yyz_0[j];

                tdx_xy_xyzz_0[j] = pa_x[j] * tdx_y_xyzz_0[j] + 0.5 * fl1_fx * tdx_y_yzz_0[j] + 0.5 * fl1_fx * ts_y_xyzz_0[j];

                tdy_xy_xyzz_0[j] = pa_x[j] * tdy_y_xyzz_0[j] + 0.5 * fl1_fx * tdy_y_yzz_0[j];

                tdz_xy_xyzz_0[j] = pa_x[j] * tdz_y_xyzz_0[j] + 0.5 * fl1_fx * tdz_y_yzz_0[j];

                tdx_xy_xzzz_0[j] = pa_x[j] * tdx_y_xzzz_0[j] + 0.5 * fl1_fx * tdx_y_zzz_0[j] + 0.5 * fl1_fx * ts_y_xzzz_0[j];

                tdy_xy_xzzz_0[j] = pa_x[j] * tdy_y_xzzz_0[j] + 0.5 * fl1_fx * tdy_y_zzz_0[j];

                tdz_xy_xzzz_0[j] = pa_x[j] * tdz_y_xzzz_0[j] + 0.5 * fl1_fx * tdz_y_zzz_0[j];

                tdx_xy_yyyy_0[j] = pa_x[j] * tdx_y_yyyy_0[j] + 0.5 * fl1_fx * ts_y_yyyy_0[j];

                tdy_xy_yyyy_0[j] = pa_x[j] * tdy_y_yyyy_0[j];

                tdz_xy_yyyy_0[j] = pa_x[j] * tdz_y_yyyy_0[j];

                tdx_xy_yyyz_0[j] = pa_x[j] * tdx_y_yyyz_0[j] + 0.5 * fl1_fx * ts_y_yyyz_0[j];

                tdy_xy_yyyz_0[j] = pa_x[j] * tdy_y_yyyz_0[j];

                tdz_xy_yyyz_0[j] = pa_x[j] * tdz_y_yyyz_0[j];

                tdx_xy_yyzz_0[j] = pa_x[j] * tdx_y_yyzz_0[j] + 0.5 * fl1_fx * ts_y_yyzz_0[j];

                tdy_xy_yyzz_0[j] = pa_x[j] * tdy_y_yyzz_0[j];

                tdz_xy_yyzz_0[j] = pa_x[j] * tdz_y_yyzz_0[j];

                tdx_xy_yzzz_0[j] = pa_x[j] * tdx_y_yzzz_0[j] + 0.5 * fl1_fx * ts_y_yzzz_0[j];

                tdy_xy_yzzz_0[j] = pa_x[j] * tdy_y_yzzz_0[j];

                tdz_xy_yzzz_0[j] = pa_x[j] * tdz_y_yzzz_0[j];

                tdx_xy_zzzz_0[j] = pa_x[j] * tdx_y_zzzz_0[j] + 0.5 * fl1_fx * ts_y_zzzz_0[j];

                tdy_xy_zzzz_0[j] = pa_x[j] * tdy_y_zzzz_0[j];

                tdz_xy_zzzz_0[j] = pa_x[j] * tdz_y_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDG_90_135(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
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

        auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 30); 

            auto tdy_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tdz_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tdx_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 31); 

            auto tdy_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tdz_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tdx_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 32); 

            auto tdy_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tdz_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tdx_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 33); 

            auto tdy_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tdz_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tdx_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 34); 

            auto tdy_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tdz_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tdx_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 35); 

            auto tdy_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tdz_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tdx_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 36); 

            auto tdy_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tdz_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tdx_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 37); 

            auto tdy_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tdz_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tdx_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 38); 

            auto tdy_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tdz_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tdx_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 39); 

            auto tdy_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tdz_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tdx_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 40); 

            auto tdy_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tdz_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tdx_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 41); 

            auto tdy_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tdz_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tdx_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 42); 

            auto tdy_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tdz_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tdx_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 43); 

            auto tdy_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tdz_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tdx_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 44); 

            auto tdy_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            auto tdz_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20); 

            auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21); 

            auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22); 

            auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23); 

            auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24); 

            auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24); 

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

            auto ts_z_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 30); 

            auto ts_z_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 31); 

            auto ts_z_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 32); 

            auto ts_z_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 33); 

            auto ts_z_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 34); 

            auto ts_z_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 35); 

            auto ts_z_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 36); 

            auto ts_z_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 37); 

            auto ts_z_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 38); 

            auto ts_z_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 39); 

            auto ts_z_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 40); 

            auto ts_z_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 41); 

            auto ts_z_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 42); 

            auto ts_z_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 43); 

            auto ts_z_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 44); 

            // set up pointers to integrals

            auto tdx_xz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 30); 

            auto tdy_xz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 30); 

            auto tdz_xz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 30); 

            auto tdx_xz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 31); 

            auto tdy_xz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 31); 

            auto tdz_xz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 31); 

            auto tdx_xz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 32); 

            auto tdy_xz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 32); 

            auto tdz_xz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 32); 

            auto tdx_xz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 33); 

            auto tdy_xz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 33); 

            auto tdz_xz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 33); 

            auto tdx_xz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 34); 

            auto tdy_xz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 34); 

            auto tdz_xz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 34); 

            auto tdx_xz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 35); 

            auto tdy_xz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 35); 

            auto tdz_xz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 35); 

            auto tdx_xz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 36); 

            auto tdy_xz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 36); 

            auto tdz_xz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 36); 

            auto tdx_xz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 37); 

            auto tdy_xz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 37); 

            auto tdz_xz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 37); 

            auto tdx_xz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 38); 

            auto tdy_xz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 38); 

            auto tdz_xz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 38); 

            auto tdx_xz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 39); 

            auto tdy_xz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 39); 

            auto tdz_xz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 39); 

            auto tdx_xz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 40); 

            auto tdy_xz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 40); 

            auto tdz_xz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 40); 

            auto tdx_xz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 41); 

            auto tdy_xz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 41); 

            auto tdz_xz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 41); 

            auto tdx_xz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 42); 

            auto tdy_xz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 42); 

            auto tdz_xz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 42); 

            auto tdx_xz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 43); 

            auto tdy_xz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 43); 

            auto tdz_xz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 43); 

            auto tdx_xz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 44); 

            auto tdy_xz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 44); 

            auto tdz_xz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 44); 

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fx, pa_x, tdx_xz_xxxx_0, tdx_xz_xxxy_0, tdx_xz_xxxz_0, tdx_xz_xxyy_0, \
                                     tdx_xz_xxyz_0, tdx_xz_xxzz_0, tdx_xz_xyyy_0, tdx_xz_xyyz_0, tdx_xz_xyzz_0, \
                                     tdx_xz_xzzz_0, tdx_xz_yyyy_0, tdx_xz_yyyz_0, tdx_xz_yyzz_0, tdx_xz_yzzz_0, \
                                     tdx_xz_zzzz_0, tdx_z_xxx_0, tdx_z_xxxx_0, tdx_z_xxxy_0, tdx_z_xxxz_0, tdx_z_xxy_0, \
                                     tdx_z_xxyy_0, tdx_z_xxyz_0, tdx_z_xxz_0, tdx_z_xxzz_0, tdx_z_xyy_0, tdx_z_xyyy_0, \
                                     tdx_z_xyyz_0, tdx_z_xyz_0, tdx_z_xyzz_0, tdx_z_xzz_0, tdx_z_xzzz_0, tdx_z_yyy_0, \
                                     tdx_z_yyyy_0, tdx_z_yyyz_0, tdx_z_yyz_0, tdx_z_yyzz_0, tdx_z_yzz_0, tdx_z_yzzz_0, \
                                     tdx_z_zzz_0, tdx_z_zzzz_0, tdy_xz_xxxx_0, tdy_xz_xxxy_0, tdy_xz_xxxz_0, \
                                     tdy_xz_xxyy_0, tdy_xz_xxyz_0, tdy_xz_xxzz_0, tdy_xz_xyyy_0, tdy_xz_xyyz_0, \
                                     tdy_xz_xyzz_0, tdy_xz_xzzz_0, tdy_xz_yyyy_0, tdy_xz_yyyz_0, tdy_xz_yyzz_0, \
                                     tdy_xz_yzzz_0, tdy_xz_zzzz_0, tdy_z_xxx_0, tdy_z_xxxx_0, tdy_z_xxxy_0, tdy_z_xxxz_0, \
                                     tdy_z_xxy_0, tdy_z_xxyy_0, tdy_z_xxyz_0, tdy_z_xxz_0, tdy_z_xxzz_0, tdy_z_xyy_0, \
                                     tdy_z_xyyy_0, tdy_z_xyyz_0, tdy_z_xyz_0, tdy_z_xyzz_0, tdy_z_xzz_0, tdy_z_xzzz_0, \
                                     tdy_z_yyy_0, tdy_z_yyyy_0, tdy_z_yyyz_0, tdy_z_yyz_0, tdy_z_yyzz_0, tdy_z_yzz_0, \
                                     tdy_z_yzzz_0, tdy_z_zzz_0, tdy_z_zzzz_0, tdz_xz_xxxx_0, tdz_xz_xxxy_0, \
                                     tdz_xz_xxxz_0, tdz_xz_xxyy_0, tdz_xz_xxyz_0, tdz_xz_xxzz_0, tdz_xz_xyyy_0, \
                                     tdz_xz_xyyz_0, tdz_xz_xyzz_0, tdz_xz_xzzz_0, tdz_xz_yyyy_0, tdz_xz_yyyz_0, \
                                     tdz_xz_yyzz_0, tdz_xz_yzzz_0, tdz_xz_zzzz_0, tdz_z_xxx_0, tdz_z_xxxx_0, \
                                     tdz_z_xxxy_0, tdz_z_xxxz_0, tdz_z_xxy_0, tdz_z_xxyy_0, tdz_z_xxyz_0, tdz_z_xxz_0, \
                                     tdz_z_xxzz_0, tdz_z_xyy_0, tdz_z_xyyy_0, tdz_z_xyyz_0, tdz_z_xyz_0, tdz_z_xyzz_0, \
                                     tdz_z_xzz_0, tdz_z_xzzz_0, tdz_z_yyy_0, tdz_z_yyyy_0, tdz_z_yyyz_0, tdz_z_yyz_0, \
                                     tdz_z_yyzz_0, tdz_z_yzz_0, tdz_z_yzzz_0, tdz_z_zzz_0, tdz_z_zzzz_0, ts_z_xxxx_0, \
                                     ts_z_xxxy_0, ts_z_xxxz_0, ts_z_xxyy_0, ts_z_xxyz_0, ts_z_xxzz_0, ts_z_xyyy_0, \
                                     ts_z_xyyz_0, ts_z_xyzz_0, ts_z_xzzz_0, ts_z_yyyy_0, ts_z_yyyz_0, ts_z_yyzz_0, \
                                     ts_z_yzzz_0, ts_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xz_xxxx_0[j] = pa_x[j] * tdx_z_xxxx_0[j] + 2.0 * fl1_fx * tdx_z_xxx_0[j] + 0.5 * fl1_fx * ts_z_xxxx_0[j];

                tdy_xz_xxxx_0[j] = pa_x[j] * tdy_z_xxxx_0[j] + 2.0 * fl1_fx * tdy_z_xxx_0[j];

                tdz_xz_xxxx_0[j] = pa_x[j] * tdz_z_xxxx_0[j] + 2.0 * fl1_fx * tdz_z_xxx_0[j];

                tdx_xz_xxxy_0[j] = pa_x[j] * tdx_z_xxxy_0[j] + 1.5 * fl1_fx * tdx_z_xxy_0[j] + 0.5 * fl1_fx * ts_z_xxxy_0[j];

                tdy_xz_xxxy_0[j] = pa_x[j] * tdy_z_xxxy_0[j] + 1.5 * fl1_fx * tdy_z_xxy_0[j];

                tdz_xz_xxxy_0[j] = pa_x[j] * tdz_z_xxxy_0[j] + 1.5 * fl1_fx * tdz_z_xxy_0[j];

                tdx_xz_xxxz_0[j] = pa_x[j] * tdx_z_xxxz_0[j] + 1.5 * fl1_fx * tdx_z_xxz_0[j] + 0.5 * fl1_fx * ts_z_xxxz_0[j];

                tdy_xz_xxxz_0[j] = pa_x[j] * tdy_z_xxxz_0[j] + 1.5 * fl1_fx * tdy_z_xxz_0[j];

                tdz_xz_xxxz_0[j] = pa_x[j] * tdz_z_xxxz_0[j] + 1.5 * fl1_fx * tdz_z_xxz_0[j];

                tdx_xz_xxyy_0[j] = pa_x[j] * tdx_z_xxyy_0[j] + fl1_fx * tdx_z_xyy_0[j] + 0.5 * fl1_fx * ts_z_xxyy_0[j];

                tdy_xz_xxyy_0[j] = pa_x[j] * tdy_z_xxyy_0[j] + fl1_fx * tdy_z_xyy_0[j];

                tdz_xz_xxyy_0[j] = pa_x[j] * tdz_z_xxyy_0[j] + fl1_fx * tdz_z_xyy_0[j];

                tdx_xz_xxyz_0[j] = pa_x[j] * tdx_z_xxyz_0[j] + fl1_fx * tdx_z_xyz_0[j] + 0.5 * fl1_fx * ts_z_xxyz_0[j];

                tdy_xz_xxyz_0[j] = pa_x[j] * tdy_z_xxyz_0[j] + fl1_fx * tdy_z_xyz_0[j];

                tdz_xz_xxyz_0[j] = pa_x[j] * tdz_z_xxyz_0[j] + fl1_fx * tdz_z_xyz_0[j];

                tdx_xz_xxzz_0[j] = pa_x[j] * tdx_z_xxzz_0[j] + fl1_fx * tdx_z_xzz_0[j] + 0.5 * fl1_fx * ts_z_xxzz_0[j];

                tdy_xz_xxzz_0[j] = pa_x[j] * tdy_z_xxzz_0[j] + fl1_fx * tdy_z_xzz_0[j];

                tdz_xz_xxzz_0[j] = pa_x[j] * tdz_z_xxzz_0[j] + fl1_fx * tdz_z_xzz_0[j];

                tdx_xz_xyyy_0[j] = pa_x[j] * tdx_z_xyyy_0[j] + 0.5 * fl1_fx * tdx_z_yyy_0[j] + 0.5 * fl1_fx * ts_z_xyyy_0[j];

                tdy_xz_xyyy_0[j] = pa_x[j] * tdy_z_xyyy_0[j] + 0.5 * fl1_fx * tdy_z_yyy_0[j];

                tdz_xz_xyyy_0[j] = pa_x[j] * tdz_z_xyyy_0[j] + 0.5 * fl1_fx * tdz_z_yyy_0[j];

                tdx_xz_xyyz_0[j] = pa_x[j] * tdx_z_xyyz_0[j] + 0.5 * fl1_fx * tdx_z_yyz_0[j] + 0.5 * fl1_fx * ts_z_xyyz_0[j];

                tdy_xz_xyyz_0[j] = pa_x[j] * tdy_z_xyyz_0[j] + 0.5 * fl1_fx * tdy_z_yyz_0[j];

                tdz_xz_xyyz_0[j] = pa_x[j] * tdz_z_xyyz_0[j] + 0.5 * fl1_fx * tdz_z_yyz_0[j];

                tdx_xz_xyzz_0[j] = pa_x[j] * tdx_z_xyzz_0[j] + 0.5 * fl1_fx * tdx_z_yzz_0[j] + 0.5 * fl1_fx * ts_z_xyzz_0[j];

                tdy_xz_xyzz_0[j] = pa_x[j] * tdy_z_xyzz_0[j] + 0.5 * fl1_fx * tdy_z_yzz_0[j];

                tdz_xz_xyzz_0[j] = pa_x[j] * tdz_z_xyzz_0[j] + 0.5 * fl1_fx * tdz_z_yzz_0[j];

                tdx_xz_xzzz_0[j] = pa_x[j] * tdx_z_xzzz_0[j] + 0.5 * fl1_fx * tdx_z_zzz_0[j] + 0.5 * fl1_fx * ts_z_xzzz_0[j];

                tdy_xz_xzzz_0[j] = pa_x[j] * tdy_z_xzzz_0[j] + 0.5 * fl1_fx * tdy_z_zzz_0[j];

                tdz_xz_xzzz_0[j] = pa_x[j] * tdz_z_xzzz_0[j] + 0.5 * fl1_fx * tdz_z_zzz_0[j];

                tdx_xz_yyyy_0[j] = pa_x[j] * tdx_z_yyyy_0[j] + 0.5 * fl1_fx * ts_z_yyyy_0[j];

                tdy_xz_yyyy_0[j] = pa_x[j] * tdy_z_yyyy_0[j];

                tdz_xz_yyyy_0[j] = pa_x[j] * tdz_z_yyyy_0[j];

                tdx_xz_yyyz_0[j] = pa_x[j] * tdx_z_yyyz_0[j] + 0.5 * fl1_fx * ts_z_yyyz_0[j];

                tdy_xz_yyyz_0[j] = pa_x[j] * tdy_z_yyyz_0[j];

                tdz_xz_yyyz_0[j] = pa_x[j] * tdz_z_yyyz_0[j];

                tdx_xz_yyzz_0[j] = pa_x[j] * tdx_z_yyzz_0[j] + 0.5 * fl1_fx * ts_z_yyzz_0[j];

                tdy_xz_yyzz_0[j] = pa_x[j] * tdy_z_yyzz_0[j];

                tdz_xz_yyzz_0[j] = pa_x[j] * tdz_z_yyzz_0[j];

                tdx_xz_yzzz_0[j] = pa_x[j] * tdx_z_yzzz_0[j] + 0.5 * fl1_fx * ts_z_yzzz_0[j];

                tdy_xz_yzzz_0[j] = pa_x[j] * tdy_z_yzzz_0[j];

                tdz_xz_yzzz_0[j] = pa_x[j] * tdz_z_yzzz_0[j];

                tdx_xz_zzzz_0[j] = pa_x[j] * tdx_z_zzzz_0[j] + 0.5 * fl1_fx * ts_z_zzzz_0[j];

                tdy_xz_zzzz_0[j] = pa_x[j] * tdy_z_zzzz_0[j];

                tdz_xz_zzzz_0[j] = pa_x[j] * tdz_z_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDG_135_180(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
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

        auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tdx_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 15); 

            auto tdy_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 15); 

            auto tdz_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 15); 

            auto tdx_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 16); 

            auto tdy_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 16); 

            auto tdz_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 16); 

            auto tdx_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 17); 

            auto tdy_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 17); 

            auto tdz_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 17); 

            auto tdx_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 18); 

            auto tdy_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 18); 

            auto tdz_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 18); 

            auto tdx_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 19); 

            auto tdy_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 19); 

            auto tdz_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 19); 

            auto tdx_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 20); 

            auto tdy_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 20); 

            auto tdz_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 20); 

            auto tdx_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 21); 

            auto tdy_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 21); 

            auto tdz_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 21); 

            auto tdx_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 22); 

            auto tdy_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 22); 

            auto tdz_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 22); 

            auto tdx_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 23); 

            auto tdy_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 23); 

            auto tdz_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 23); 

            auto tdx_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 24); 

            auto tdy_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 24); 

            auto tdz_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 24); 

            auto tdx_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 25); 

            auto tdy_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 25); 

            auto tdz_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 25); 

            auto tdx_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 26); 

            auto tdy_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 26); 

            auto tdz_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 26); 

            auto tdx_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 27); 

            auto tdy_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 27); 

            auto tdz_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 27); 

            auto tdx_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 28); 

            auto tdy_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 28); 

            auto tdz_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 28); 

            auto tdx_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 29); 

            auto tdy_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 29); 

            auto tdz_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 29); 

            auto tdx_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx); 

            auto tdy_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx); 

            auto tdz_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx); 

            auto tdx_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 1); 

            auto tdy_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 1); 

            auto tdz_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 1); 

            auto tdx_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 2); 

            auto tdy_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 2); 

            auto tdz_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 2); 

            auto tdx_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 3); 

            auto tdy_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 3); 

            auto tdz_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 3); 

            auto tdx_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 4); 

            auto tdy_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 4); 

            auto tdz_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 4); 

            auto tdx_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 5); 

            auto tdy_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 5); 

            auto tdz_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 5); 

            auto tdx_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 6); 

            auto tdy_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 6); 

            auto tdz_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 6); 

            auto tdx_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 7); 

            auto tdy_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 7); 

            auto tdz_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 7); 

            auto tdx_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 8); 

            auto tdy_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 8); 

            auto tdz_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 8); 

            auto tdx_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 9); 

            auto tdy_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 9); 

            auto tdz_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 9); 

            auto tdx_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 10); 

            auto tdy_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 10); 

            auto tdz_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 10); 

            auto tdx_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 11); 

            auto tdy_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 11); 

            auto tdz_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 11); 

            auto tdx_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 12); 

            auto tdy_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 12); 

            auto tdz_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 12); 

            auto tdx_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 13); 

            auto tdy_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 13); 

            auto tdz_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 13); 

            auto tdx_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 14); 

            auto tdy_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 14); 

            auto tdz_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 14); 

            auto tdx_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 10); 

            auto tdy_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tdz_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tdx_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 11); 

            auto tdy_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tdz_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tdx_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 12); 

            auto tdy_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tdz_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tdx_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 13); 

            auto tdy_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tdz_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tdx_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 14); 

            auto tdy_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tdz_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tdx_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 15); 

            auto tdy_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tdz_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tdx_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 16); 

            auto tdy_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tdz_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tdx_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 17); 

            auto tdy_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tdz_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tdx_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 18); 

            auto tdy_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tdz_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tdx_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 19); 

            auto tdy_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tdz_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto ts_y_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 15); 

            auto ts_y_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 16); 

            auto ts_y_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 17); 

            auto ts_y_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 18); 

            auto ts_y_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 19); 

            auto ts_y_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 20); 

            auto ts_y_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 21); 

            auto ts_y_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 22); 

            auto ts_y_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 23); 

            auto ts_y_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 24); 

            auto ts_y_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 25); 

            auto ts_y_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 26); 

            auto ts_y_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 27); 

            auto ts_y_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 28); 

            auto ts_y_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 29); 

            // set up pointers to integrals

            auto tdx_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 45); 

            auto tdy_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 45); 

            auto tdz_yy_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 45); 

            auto tdx_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 46); 

            auto tdy_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 46); 

            auto tdz_yy_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 46); 

            auto tdx_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 47); 

            auto tdy_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 47); 

            auto tdz_yy_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 47); 

            auto tdx_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 48); 

            auto tdy_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 48); 

            auto tdz_yy_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 48); 

            auto tdx_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 49); 

            auto tdy_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 49); 

            auto tdz_yy_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 49); 

            auto tdx_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 50); 

            auto tdy_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 50); 

            auto tdz_yy_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 50); 

            auto tdx_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 51); 

            auto tdy_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 51); 

            auto tdz_yy_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 51); 

            auto tdx_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 52); 

            auto tdy_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 52); 

            auto tdz_yy_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 52); 

            auto tdx_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 53); 

            auto tdy_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 53); 

            auto tdz_yy_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 53); 

            auto tdx_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 54); 

            auto tdy_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 54); 

            auto tdz_yy_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 54); 

            auto tdx_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 55); 

            auto tdy_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 55); 

            auto tdz_yy_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 55); 

            auto tdx_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 56); 

            auto tdy_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 56); 

            auto tdz_yy_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 56); 

            auto tdx_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 57); 

            auto tdy_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 57); 

            auto tdz_yy_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 57); 

            auto tdx_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 58); 

            auto tdy_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 58); 

            auto tdz_yy_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 58); 

            auto tdx_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 59); 

            auto tdy_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 59); 

            auto tdz_yy_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 59); 

            // Batch of Integrals (135,180)

            #pragma omp simd aligned(fx, pa_y, tdx_0_xxxx_0, tdx_0_xxxy_0, tdx_0_xxxz_0, tdx_0_xxyy_0, \
                                     tdx_0_xxyz_0, tdx_0_xxzz_0, tdx_0_xyyy_0, tdx_0_xyyz_0, tdx_0_xyzz_0, tdx_0_xzzz_0, \
                                     tdx_0_yyyy_0, tdx_0_yyyz_0, tdx_0_yyzz_0, tdx_0_yzzz_0, tdx_0_zzzz_0, tdx_y_xxx_0, \
                                     tdx_y_xxxx_0, tdx_y_xxxy_0, tdx_y_xxxz_0, tdx_y_xxy_0, tdx_y_xxyy_0, tdx_y_xxyz_0, \
                                     tdx_y_xxz_0, tdx_y_xxzz_0, tdx_y_xyy_0, tdx_y_xyyy_0, tdx_y_xyyz_0, tdx_y_xyz_0, \
                                     tdx_y_xyzz_0, tdx_y_xzz_0, tdx_y_xzzz_0, tdx_y_yyy_0, tdx_y_yyyy_0, tdx_y_yyyz_0, \
                                     tdx_y_yyz_0, tdx_y_yyzz_0, tdx_y_yzz_0, tdx_y_yzzz_0, tdx_y_zzz_0, tdx_y_zzzz_0, \
                                     tdx_yy_xxxx_0, tdx_yy_xxxy_0, tdx_yy_xxxz_0, tdx_yy_xxyy_0, tdx_yy_xxyz_0, \
                                     tdx_yy_xxzz_0, tdx_yy_xyyy_0, tdx_yy_xyyz_0, tdx_yy_xyzz_0, tdx_yy_xzzz_0, \
                                     tdx_yy_yyyy_0, tdx_yy_yyyz_0, tdx_yy_yyzz_0, tdx_yy_yzzz_0, tdx_yy_zzzz_0, \
                                     tdy_0_xxxx_0, tdy_0_xxxy_0, tdy_0_xxxz_0, tdy_0_xxyy_0, tdy_0_xxyz_0, tdy_0_xxzz_0, \
                                     tdy_0_xyyy_0, tdy_0_xyyz_0, tdy_0_xyzz_0, tdy_0_xzzz_0, tdy_0_yyyy_0, tdy_0_yyyz_0, \
                                     tdy_0_yyzz_0, tdy_0_yzzz_0, tdy_0_zzzz_0, tdy_y_xxx_0, tdy_y_xxxx_0, tdy_y_xxxy_0, \
                                     tdy_y_xxxz_0, tdy_y_xxy_0, tdy_y_xxyy_0, tdy_y_xxyz_0, tdy_y_xxz_0, tdy_y_xxzz_0, \
                                     tdy_y_xyy_0, tdy_y_xyyy_0, tdy_y_xyyz_0, tdy_y_xyz_0, tdy_y_xyzz_0, tdy_y_xzz_0, \
                                     tdy_y_xzzz_0, tdy_y_yyy_0, tdy_y_yyyy_0, tdy_y_yyyz_0, tdy_y_yyz_0, tdy_y_yyzz_0, \
                                     tdy_y_yzz_0, tdy_y_yzzz_0, tdy_y_zzz_0, tdy_y_zzzz_0, tdy_yy_xxxx_0, \
                                     tdy_yy_xxxy_0, tdy_yy_xxxz_0, tdy_yy_xxyy_0, tdy_yy_xxyz_0, tdy_yy_xxzz_0, \
                                     tdy_yy_xyyy_0, tdy_yy_xyyz_0, tdy_yy_xyzz_0, tdy_yy_xzzz_0, tdy_yy_yyyy_0, \
                                     tdy_yy_yyyz_0, tdy_yy_yyzz_0, tdy_yy_yzzz_0, tdy_yy_zzzz_0, tdz_0_xxxx_0, \
                                     tdz_0_xxxy_0, tdz_0_xxxz_0, tdz_0_xxyy_0, tdz_0_xxyz_0, tdz_0_xxzz_0, tdz_0_xyyy_0, \
                                     tdz_0_xyyz_0, tdz_0_xyzz_0, tdz_0_xzzz_0, tdz_0_yyyy_0, tdz_0_yyyz_0, tdz_0_yyzz_0, \
                                     tdz_0_yzzz_0, tdz_0_zzzz_0, tdz_y_xxx_0, tdz_y_xxxx_0, tdz_y_xxxy_0, tdz_y_xxxz_0, \
                                     tdz_y_xxy_0, tdz_y_xxyy_0, tdz_y_xxyz_0, tdz_y_xxz_0, tdz_y_xxzz_0, tdz_y_xyy_0, \
                                     tdz_y_xyyy_0, tdz_y_xyyz_0, tdz_y_xyz_0, tdz_y_xyzz_0, tdz_y_xzz_0, tdz_y_xzzz_0, \
                                     tdz_y_yyy_0, tdz_y_yyyy_0, tdz_y_yyyz_0, tdz_y_yyz_0, tdz_y_yyzz_0, tdz_y_yzz_0, \
                                     tdz_y_yzzz_0, tdz_y_zzz_0, tdz_y_zzzz_0, tdz_yy_xxxx_0, tdz_yy_xxxy_0, \
                                     tdz_yy_xxxz_0, tdz_yy_xxyy_0, tdz_yy_xxyz_0, tdz_yy_xxzz_0, tdz_yy_xyyy_0, \
                                     tdz_yy_xyyz_0, tdz_yy_xyzz_0, tdz_yy_xzzz_0, tdz_yy_yyyy_0, tdz_yy_yyyz_0, \
                                     tdz_yy_yyzz_0, tdz_yy_yzzz_0, tdz_yy_zzzz_0, ts_y_xxxx_0, ts_y_xxxy_0, ts_y_xxxz_0, \
                                     ts_y_xxyy_0, ts_y_xxyz_0, ts_y_xxzz_0, ts_y_xyyy_0, ts_y_xyyz_0, ts_y_xyzz_0, \
                                     ts_y_xzzz_0, ts_y_yyyy_0, ts_y_yyyz_0, ts_y_yyzz_0, ts_y_yzzz_0, ts_y_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yy_xxxx_0[j] = pa_y[j] * tdx_y_xxxx_0[j] + 0.5 * fl1_fx * tdx_0_xxxx_0[j];

                tdy_yy_xxxx_0[j] = pa_y[j] * tdy_y_xxxx_0[j] + 0.5 * fl1_fx * tdy_0_xxxx_0[j] + 0.5 * fl1_fx * ts_y_xxxx_0[j];

                tdz_yy_xxxx_0[j] = pa_y[j] * tdz_y_xxxx_0[j] + 0.5 * fl1_fx * tdz_0_xxxx_0[j];

                tdx_yy_xxxy_0[j] = pa_y[j] * tdx_y_xxxy_0[j] + 0.5 * fl1_fx * tdx_0_xxxy_0[j] + 0.5 * fl1_fx * tdx_y_xxx_0[j];

                tdy_yy_xxxy_0[j] = pa_y[j] * tdy_y_xxxy_0[j] + 0.5 * fl1_fx * tdy_0_xxxy_0[j] + 0.5 * fl1_fx * tdy_y_xxx_0[j] + 0.5 * fl1_fx * ts_y_xxxy_0[j];

                tdz_yy_xxxy_0[j] = pa_y[j] * tdz_y_xxxy_0[j] + 0.5 * fl1_fx * tdz_0_xxxy_0[j] + 0.5 * fl1_fx * tdz_y_xxx_0[j];

                tdx_yy_xxxz_0[j] = pa_y[j] * tdx_y_xxxz_0[j] + 0.5 * fl1_fx * tdx_0_xxxz_0[j];

                tdy_yy_xxxz_0[j] = pa_y[j] * tdy_y_xxxz_0[j] + 0.5 * fl1_fx * tdy_0_xxxz_0[j] + 0.5 * fl1_fx * ts_y_xxxz_0[j];

                tdz_yy_xxxz_0[j] = pa_y[j] * tdz_y_xxxz_0[j] + 0.5 * fl1_fx * tdz_0_xxxz_0[j];

                tdx_yy_xxyy_0[j] = pa_y[j] * tdx_y_xxyy_0[j] + 0.5 * fl1_fx * tdx_0_xxyy_0[j] + fl1_fx * tdx_y_xxy_0[j];

                tdy_yy_xxyy_0[j] = pa_y[j] * tdy_y_xxyy_0[j] + 0.5 * fl1_fx * tdy_0_xxyy_0[j] + fl1_fx * tdy_y_xxy_0[j] + 0.5 * fl1_fx * ts_y_xxyy_0[j];

                tdz_yy_xxyy_0[j] = pa_y[j] * tdz_y_xxyy_0[j] + 0.5 * fl1_fx * tdz_0_xxyy_0[j] + fl1_fx * tdz_y_xxy_0[j];

                tdx_yy_xxyz_0[j] = pa_y[j] * tdx_y_xxyz_0[j] + 0.5 * fl1_fx * tdx_0_xxyz_0[j] + 0.5 * fl1_fx * tdx_y_xxz_0[j];

                tdy_yy_xxyz_0[j] = pa_y[j] * tdy_y_xxyz_0[j] + 0.5 * fl1_fx * tdy_0_xxyz_0[j] + 0.5 * fl1_fx * tdy_y_xxz_0[j] + 0.5 * fl1_fx * ts_y_xxyz_0[j];

                tdz_yy_xxyz_0[j] = pa_y[j] * tdz_y_xxyz_0[j] + 0.5 * fl1_fx * tdz_0_xxyz_0[j] + 0.5 * fl1_fx * tdz_y_xxz_0[j];

                tdx_yy_xxzz_0[j] = pa_y[j] * tdx_y_xxzz_0[j] + 0.5 * fl1_fx * tdx_0_xxzz_0[j];

                tdy_yy_xxzz_0[j] = pa_y[j] * tdy_y_xxzz_0[j] + 0.5 * fl1_fx * tdy_0_xxzz_0[j] + 0.5 * fl1_fx * ts_y_xxzz_0[j];

                tdz_yy_xxzz_0[j] = pa_y[j] * tdz_y_xxzz_0[j] + 0.5 * fl1_fx * tdz_0_xxzz_0[j];

                tdx_yy_xyyy_0[j] = pa_y[j] * tdx_y_xyyy_0[j] + 0.5 * fl1_fx * tdx_0_xyyy_0[j] + 1.5 * fl1_fx * tdx_y_xyy_0[j];

                tdy_yy_xyyy_0[j] = pa_y[j] * tdy_y_xyyy_0[j] + 0.5 * fl1_fx * tdy_0_xyyy_0[j] + 1.5 * fl1_fx * tdy_y_xyy_0[j] + 0.5 * fl1_fx * ts_y_xyyy_0[j];

                tdz_yy_xyyy_0[j] = pa_y[j] * tdz_y_xyyy_0[j] + 0.5 * fl1_fx * tdz_0_xyyy_0[j] + 1.5 * fl1_fx * tdz_y_xyy_0[j];

                tdx_yy_xyyz_0[j] = pa_y[j] * tdx_y_xyyz_0[j] + 0.5 * fl1_fx * tdx_0_xyyz_0[j] + fl1_fx * tdx_y_xyz_0[j];

                tdy_yy_xyyz_0[j] = pa_y[j] * tdy_y_xyyz_0[j] + 0.5 * fl1_fx * tdy_0_xyyz_0[j] + fl1_fx * tdy_y_xyz_0[j] + 0.5 * fl1_fx * ts_y_xyyz_0[j];

                tdz_yy_xyyz_0[j] = pa_y[j] * tdz_y_xyyz_0[j] + 0.5 * fl1_fx * tdz_0_xyyz_0[j] + fl1_fx * tdz_y_xyz_0[j];

                tdx_yy_xyzz_0[j] = pa_y[j] * tdx_y_xyzz_0[j] + 0.5 * fl1_fx * tdx_0_xyzz_0[j] + 0.5 * fl1_fx * tdx_y_xzz_0[j];

                tdy_yy_xyzz_0[j] = pa_y[j] * tdy_y_xyzz_0[j] + 0.5 * fl1_fx * tdy_0_xyzz_0[j] + 0.5 * fl1_fx * tdy_y_xzz_0[j] + 0.5 * fl1_fx * ts_y_xyzz_0[j];

                tdz_yy_xyzz_0[j] = pa_y[j] * tdz_y_xyzz_0[j] + 0.5 * fl1_fx * tdz_0_xyzz_0[j] + 0.5 * fl1_fx * tdz_y_xzz_0[j];

                tdx_yy_xzzz_0[j] = pa_y[j] * tdx_y_xzzz_0[j] + 0.5 * fl1_fx * tdx_0_xzzz_0[j];

                tdy_yy_xzzz_0[j] = pa_y[j] * tdy_y_xzzz_0[j] + 0.5 * fl1_fx * tdy_0_xzzz_0[j] + 0.5 * fl1_fx * ts_y_xzzz_0[j];

                tdz_yy_xzzz_0[j] = pa_y[j] * tdz_y_xzzz_0[j] + 0.5 * fl1_fx * tdz_0_xzzz_0[j];

                tdx_yy_yyyy_0[j] = pa_y[j] * tdx_y_yyyy_0[j] + 0.5 * fl1_fx * tdx_0_yyyy_0[j] + 2.0 * fl1_fx * tdx_y_yyy_0[j];

                tdy_yy_yyyy_0[j] = pa_y[j] * tdy_y_yyyy_0[j] + 0.5 * fl1_fx * tdy_0_yyyy_0[j] + 2.0 * fl1_fx * tdy_y_yyy_0[j] + 0.5 * fl1_fx * ts_y_yyyy_0[j];

                tdz_yy_yyyy_0[j] = pa_y[j] * tdz_y_yyyy_0[j] + 0.5 * fl1_fx * tdz_0_yyyy_0[j] + 2.0 * fl1_fx * tdz_y_yyy_0[j];

                tdx_yy_yyyz_0[j] = pa_y[j] * tdx_y_yyyz_0[j] + 0.5 * fl1_fx * tdx_0_yyyz_0[j] + 1.5 * fl1_fx * tdx_y_yyz_0[j];

                tdy_yy_yyyz_0[j] = pa_y[j] * tdy_y_yyyz_0[j] + 0.5 * fl1_fx * tdy_0_yyyz_0[j] + 1.5 * fl1_fx * tdy_y_yyz_0[j] + 0.5 * fl1_fx * ts_y_yyyz_0[j];

                tdz_yy_yyyz_0[j] = pa_y[j] * tdz_y_yyyz_0[j] + 0.5 * fl1_fx * tdz_0_yyyz_0[j] + 1.5 * fl1_fx * tdz_y_yyz_0[j];

                tdx_yy_yyzz_0[j] = pa_y[j] * tdx_y_yyzz_0[j] + 0.5 * fl1_fx * tdx_0_yyzz_0[j] + fl1_fx * tdx_y_yzz_0[j];

                tdy_yy_yyzz_0[j] = pa_y[j] * tdy_y_yyzz_0[j] + 0.5 * fl1_fx * tdy_0_yyzz_0[j] + fl1_fx * tdy_y_yzz_0[j] + 0.5 * fl1_fx * ts_y_yyzz_0[j];

                tdz_yy_yyzz_0[j] = pa_y[j] * tdz_y_yyzz_0[j] + 0.5 * fl1_fx * tdz_0_yyzz_0[j] + fl1_fx * tdz_y_yzz_0[j];

                tdx_yy_yzzz_0[j] = pa_y[j] * tdx_y_yzzz_0[j] + 0.5 * fl1_fx * tdx_0_yzzz_0[j] + 0.5 * fl1_fx * tdx_y_zzz_0[j];

                tdy_yy_yzzz_0[j] = pa_y[j] * tdy_y_yzzz_0[j] + 0.5 * fl1_fx * tdy_0_yzzz_0[j] + 0.5 * fl1_fx * tdy_y_zzz_0[j] + 0.5 * fl1_fx * ts_y_yzzz_0[j];

                tdz_yy_yzzz_0[j] = pa_y[j] * tdz_y_yzzz_0[j] + 0.5 * fl1_fx * tdz_0_yzzz_0[j] + 0.5 * fl1_fx * tdz_y_zzz_0[j];

                tdx_yy_zzzz_0[j] = pa_y[j] * tdx_y_zzzz_0[j] + 0.5 * fl1_fx * tdx_0_zzzz_0[j];

                tdy_yy_zzzz_0[j] = pa_y[j] * tdy_y_zzzz_0[j] + 0.5 * fl1_fx * tdy_0_zzzz_0[j] + 0.5 * fl1_fx * ts_y_zzzz_0[j];

                tdz_yy_zzzz_0[j] = pa_y[j] * tdz_y_zzzz_0[j] + 0.5 * fl1_fx * tdz_0_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDG_180_225(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
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

        auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tdx_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 30); 

            auto tdy_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tdz_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tdx_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 31); 

            auto tdy_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tdz_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tdx_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 32); 

            auto tdy_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tdz_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tdx_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 33); 

            auto tdy_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tdz_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tdx_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 34); 

            auto tdy_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tdz_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tdx_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 35); 

            auto tdy_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tdz_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tdx_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 36); 

            auto tdy_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tdz_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tdx_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 37); 

            auto tdy_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tdz_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tdx_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 38); 

            auto tdy_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tdz_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tdx_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 39); 

            auto tdy_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tdz_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tdx_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 40); 

            auto tdy_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tdz_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tdx_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 41); 

            auto tdy_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tdz_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tdx_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 42); 

            auto tdy_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tdz_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tdx_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 43); 

            auto tdy_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tdz_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tdx_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 44); 

            auto tdy_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            auto tdz_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20); 

            auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21); 

            auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22); 

            auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23); 

            auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24); 

            auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24); 

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

            auto ts_z_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 30); 

            auto ts_z_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 31); 

            auto ts_z_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 32); 

            auto ts_z_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 33); 

            auto ts_z_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 34); 

            auto ts_z_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 35); 

            auto ts_z_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 36); 

            auto ts_z_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 37); 

            auto ts_z_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 38); 

            auto ts_z_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 39); 

            auto ts_z_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 40); 

            auto ts_z_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 41); 

            auto ts_z_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 42); 

            auto ts_z_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 43); 

            auto ts_z_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 44); 

            // set up pointers to integrals

            auto tdx_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 60); 

            auto tdy_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 60); 

            auto tdz_yz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 60); 

            auto tdx_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 61); 

            auto tdy_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 61); 

            auto tdz_yz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 61); 

            auto tdx_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 62); 

            auto tdy_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 62); 

            auto tdz_yz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 62); 

            auto tdx_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 63); 

            auto tdy_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 63); 

            auto tdz_yz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 63); 

            auto tdx_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 64); 

            auto tdy_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 64); 

            auto tdz_yz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 64); 

            auto tdx_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 65); 

            auto tdy_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 65); 

            auto tdz_yz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 65); 

            auto tdx_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 66); 

            auto tdy_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 66); 

            auto tdz_yz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 66); 

            auto tdx_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 67); 

            auto tdy_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 67); 

            auto tdz_yz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 67); 

            auto tdx_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 68); 

            auto tdy_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 68); 

            auto tdz_yz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 68); 

            auto tdx_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 69); 

            auto tdy_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 69); 

            auto tdz_yz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 69); 

            auto tdx_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 70); 

            auto tdy_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 70); 

            auto tdz_yz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 70); 

            auto tdx_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 71); 

            auto tdy_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 71); 

            auto tdz_yz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 71); 

            auto tdx_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 72); 

            auto tdy_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 72); 

            auto tdz_yz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 72); 

            auto tdx_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 73); 

            auto tdy_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 73); 

            auto tdz_yz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 73); 

            auto tdx_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 74); 

            auto tdy_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 74); 

            auto tdz_yz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 74); 

            // Batch of Integrals (180,225)

            #pragma omp simd aligned(fx, pa_y, tdx_yz_xxxx_0, tdx_yz_xxxy_0, tdx_yz_xxxz_0, tdx_yz_xxyy_0, \
                                     tdx_yz_xxyz_0, tdx_yz_xxzz_0, tdx_yz_xyyy_0, tdx_yz_xyyz_0, tdx_yz_xyzz_0, \
                                     tdx_yz_xzzz_0, tdx_yz_yyyy_0, tdx_yz_yyyz_0, tdx_yz_yyzz_0, tdx_yz_yzzz_0, \
                                     tdx_yz_zzzz_0, tdx_z_xxx_0, tdx_z_xxxx_0, tdx_z_xxxy_0, tdx_z_xxxz_0, tdx_z_xxy_0, \
                                     tdx_z_xxyy_0, tdx_z_xxyz_0, tdx_z_xxz_0, tdx_z_xxzz_0, tdx_z_xyy_0, tdx_z_xyyy_0, \
                                     tdx_z_xyyz_0, tdx_z_xyz_0, tdx_z_xyzz_0, tdx_z_xzz_0, tdx_z_xzzz_0, tdx_z_yyy_0, \
                                     tdx_z_yyyy_0, tdx_z_yyyz_0, tdx_z_yyz_0, tdx_z_yyzz_0, tdx_z_yzz_0, tdx_z_yzzz_0, \
                                     tdx_z_zzz_0, tdx_z_zzzz_0, tdy_yz_xxxx_0, tdy_yz_xxxy_0, tdy_yz_xxxz_0, \
                                     tdy_yz_xxyy_0, tdy_yz_xxyz_0, tdy_yz_xxzz_0, tdy_yz_xyyy_0, tdy_yz_xyyz_0, \
                                     tdy_yz_xyzz_0, tdy_yz_xzzz_0, tdy_yz_yyyy_0, tdy_yz_yyyz_0, tdy_yz_yyzz_0, \
                                     tdy_yz_yzzz_0, tdy_yz_zzzz_0, tdy_z_xxx_0, tdy_z_xxxx_0, tdy_z_xxxy_0, tdy_z_xxxz_0, \
                                     tdy_z_xxy_0, tdy_z_xxyy_0, tdy_z_xxyz_0, tdy_z_xxz_0, tdy_z_xxzz_0, tdy_z_xyy_0, \
                                     tdy_z_xyyy_0, tdy_z_xyyz_0, tdy_z_xyz_0, tdy_z_xyzz_0, tdy_z_xzz_0, tdy_z_xzzz_0, \
                                     tdy_z_yyy_0, tdy_z_yyyy_0, tdy_z_yyyz_0, tdy_z_yyz_0, tdy_z_yyzz_0, tdy_z_yzz_0, \
                                     tdy_z_yzzz_0, tdy_z_zzz_0, tdy_z_zzzz_0, tdz_yz_xxxx_0, tdz_yz_xxxy_0, \
                                     tdz_yz_xxxz_0, tdz_yz_xxyy_0, tdz_yz_xxyz_0, tdz_yz_xxzz_0, tdz_yz_xyyy_0, \
                                     tdz_yz_xyyz_0, tdz_yz_xyzz_0, tdz_yz_xzzz_0, tdz_yz_yyyy_0, tdz_yz_yyyz_0, \
                                     tdz_yz_yyzz_0, tdz_yz_yzzz_0, tdz_yz_zzzz_0, tdz_z_xxx_0, tdz_z_xxxx_0, \
                                     tdz_z_xxxy_0, tdz_z_xxxz_0, tdz_z_xxy_0, tdz_z_xxyy_0, tdz_z_xxyz_0, tdz_z_xxz_0, \
                                     tdz_z_xxzz_0, tdz_z_xyy_0, tdz_z_xyyy_0, tdz_z_xyyz_0, tdz_z_xyz_0, tdz_z_xyzz_0, \
                                     tdz_z_xzz_0, tdz_z_xzzz_0, tdz_z_yyy_0, tdz_z_yyyy_0, tdz_z_yyyz_0, tdz_z_yyz_0, \
                                     tdz_z_yyzz_0, tdz_z_yzz_0, tdz_z_yzzz_0, tdz_z_zzz_0, tdz_z_zzzz_0, ts_z_xxxx_0, \
                                     ts_z_xxxy_0, ts_z_xxxz_0, ts_z_xxyy_0, ts_z_xxyz_0, ts_z_xxzz_0, ts_z_xyyy_0, \
                                     ts_z_xyyz_0, ts_z_xyzz_0, ts_z_xzzz_0, ts_z_yyyy_0, ts_z_yyyz_0, ts_z_yyzz_0, \
                                     ts_z_yzzz_0, ts_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yz_xxxx_0[j] = pa_y[j] * tdx_z_xxxx_0[j];

                tdy_yz_xxxx_0[j] = pa_y[j] * tdy_z_xxxx_0[j] + 0.5 * fl1_fx * ts_z_xxxx_0[j];

                tdz_yz_xxxx_0[j] = pa_y[j] * tdz_z_xxxx_0[j];

                tdx_yz_xxxy_0[j] = pa_y[j] * tdx_z_xxxy_0[j] + 0.5 * fl1_fx * tdx_z_xxx_0[j];

                tdy_yz_xxxy_0[j] = pa_y[j] * tdy_z_xxxy_0[j] + 0.5 * fl1_fx * tdy_z_xxx_0[j] + 0.5 * fl1_fx * ts_z_xxxy_0[j];

                tdz_yz_xxxy_0[j] = pa_y[j] * tdz_z_xxxy_0[j] + 0.5 * fl1_fx * tdz_z_xxx_0[j];

                tdx_yz_xxxz_0[j] = pa_y[j] * tdx_z_xxxz_0[j];

                tdy_yz_xxxz_0[j] = pa_y[j] * tdy_z_xxxz_0[j] + 0.5 * fl1_fx * ts_z_xxxz_0[j];

                tdz_yz_xxxz_0[j] = pa_y[j] * tdz_z_xxxz_0[j];

                tdx_yz_xxyy_0[j] = pa_y[j] * tdx_z_xxyy_0[j] + fl1_fx * tdx_z_xxy_0[j];

                tdy_yz_xxyy_0[j] = pa_y[j] * tdy_z_xxyy_0[j] + fl1_fx * tdy_z_xxy_0[j] + 0.5 * fl1_fx * ts_z_xxyy_0[j];

                tdz_yz_xxyy_0[j] = pa_y[j] * tdz_z_xxyy_0[j] + fl1_fx * tdz_z_xxy_0[j];

                tdx_yz_xxyz_0[j] = pa_y[j] * tdx_z_xxyz_0[j] + 0.5 * fl1_fx * tdx_z_xxz_0[j];

                tdy_yz_xxyz_0[j] = pa_y[j] * tdy_z_xxyz_0[j] + 0.5 * fl1_fx * tdy_z_xxz_0[j] + 0.5 * fl1_fx * ts_z_xxyz_0[j];

                tdz_yz_xxyz_0[j] = pa_y[j] * tdz_z_xxyz_0[j] + 0.5 * fl1_fx * tdz_z_xxz_0[j];

                tdx_yz_xxzz_0[j] = pa_y[j] * tdx_z_xxzz_0[j];

                tdy_yz_xxzz_0[j] = pa_y[j] * tdy_z_xxzz_0[j] + 0.5 * fl1_fx * ts_z_xxzz_0[j];

                tdz_yz_xxzz_0[j] = pa_y[j] * tdz_z_xxzz_0[j];

                tdx_yz_xyyy_0[j] = pa_y[j] * tdx_z_xyyy_0[j] + 1.5 * fl1_fx * tdx_z_xyy_0[j];

                tdy_yz_xyyy_0[j] = pa_y[j] * tdy_z_xyyy_0[j] + 1.5 * fl1_fx * tdy_z_xyy_0[j] + 0.5 * fl1_fx * ts_z_xyyy_0[j];

                tdz_yz_xyyy_0[j] = pa_y[j] * tdz_z_xyyy_0[j] + 1.5 * fl1_fx * tdz_z_xyy_0[j];

                tdx_yz_xyyz_0[j] = pa_y[j] * tdx_z_xyyz_0[j] + fl1_fx * tdx_z_xyz_0[j];

                tdy_yz_xyyz_0[j] = pa_y[j] * tdy_z_xyyz_0[j] + fl1_fx * tdy_z_xyz_0[j] + 0.5 * fl1_fx * ts_z_xyyz_0[j];

                tdz_yz_xyyz_0[j] = pa_y[j] * tdz_z_xyyz_0[j] + fl1_fx * tdz_z_xyz_0[j];

                tdx_yz_xyzz_0[j] = pa_y[j] * tdx_z_xyzz_0[j] + 0.5 * fl1_fx * tdx_z_xzz_0[j];

                tdy_yz_xyzz_0[j] = pa_y[j] * tdy_z_xyzz_0[j] + 0.5 * fl1_fx * tdy_z_xzz_0[j] + 0.5 * fl1_fx * ts_z_xyzz_0[j];

                tdz_yz_xyzz_0[j] = pa_y[j] * tdz_z_xyzz_0[j] + 0.5 * fl1_fx * tdz_z_xzz_0[j];

                tdx_yz_xzzz_0[j] = pa_y[j] * tdx_z_xzzz_0[j];

                tdy_yz_xzzz_0[j] = pa_y[j] * tdy_z_xzzz_0[j] + 0.5 * fl1_fx * ts_z_xzzz_0[j];

                tdz_yz_xzzz_0[j] = pa_y[j] * tdz_z_xzzz_0[j];

                tdx_yz_yyyy_0[j] = pa_y[j] * tdx_z_yyyy_0[j] + 2.0 * fl1_fx * tdx_z_yyy_0[j];

                tdy_yz_yyyy_0[j] = pa_y[j] * tdy_z_yyyy_0[j] + 2.0 * fl1_fx * tdy_z_yyy_0[j] + 0.5 * fl1_fx * ts_z_yyyy_0[j];

                tdz_yz_yyyy_0[j] = pa_y[j] * tdz_z_yyyy_0[j] + 2.0 * fl1_fx * tdz_z_yyy_0[j];

                tdx_yz_yyyz_0[j] = pa_y[j] * tdx_z_yyyz_0[j] + 1.5 * fl1_fx * tdx_z_yyz_0[j];

                tdy_yz_yyyz_0[j] = pa_y[j] * tdy_z_yyyz_0[j] + 1.5 * fl1_fx * tdy_z_yyz_0[j] + 0.5 * fl1_fx * ts_z_yyyz_0[j];

                tdz_yz_yyyz_0[j] = pa_y[j] * tdz_z_yyyz_0[j] + 1.5 * fl1_fx * tdz_z_yyz_0[j];

                tdx_yz_yyzz_0[j] = pa_y[j] * tdx_z_yyzz_0[j] + fl1_fx * tdx_z_yzz_0[j];

                tdy_yz_yyzz_0[j] = pa_y[j] * tdy_z_yyzz_0[j] + fl1_fx * tdy_z_yzz_0[j] + 0.5 * fl1_fx * ts_z_yyzz_0[j];

                tdz_yz_yyzz_0[j] = pa_y[j] * tdz_z_yyzz_0[j] + fl1_fx * tdz_z_yzz_0[j];

                tdx_yz_yzzz_0[j] = pa_y[j] * tdx_z_yzzz_0[j] + 0.5 * fl1_fx * tdx_z_zzz_0[j];

                tdy_yz_yzzz_0[j] = pa_y[j] * tdy_z_yzzz_0[j] + 0.5 * fl1_fx * tdy_z_zzz_0[j] + 0.5 * fl1_fx * ts_z_yzzz_0[j];

                tdz_yz_yzzz_0[j] = pa_y[j] * tdz_z_yzzz_0[j] + 0.5 * fl1_fx * tdz_z_zzz_0[j];

                tdx_yz_zzzz_0[j] = pa_y[j] * tdx_z_zzzz_0[j];

                tdy_yz_zzzz_0[j] = pa_y[j] * tdy_z_zzzz_0[j] + 0.5 * fl1_fx * ts_z_zzzz_0[j];

                tdz_yz_zzzz_0[j] = pa_y[j] * tdz_z_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDG_225_270(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
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

        auto pidx_d_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tdx_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 30); 

            auto tdy_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tdz_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tdx_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 31); 

            auto tdy_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tdz_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tdx_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 32); 

            auto tdy_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tdz_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tdx_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 33); 

            auto tdy_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tdz_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tdx_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 34); 

            auto tdy_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tdz_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tdx_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 35); 

            auto tdy_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tdz_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tdx_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 36); 

            auto tdy_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tdz_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tdx_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 37); 

            auto tdy_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tdz_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tdx_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 38); 

            auto tdy_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tdz_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tdx_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 39); 

            auto tdy_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tdz_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tdx_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 40); 

            auto tdy_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tdz_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tdx_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 41); 

            auto tdy_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tdz_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tdx_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 42); 

            auto tdy_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tdz_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tdx_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 43); 

            auto tdy_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tdz_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tdx_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 44); 

            auto tdy_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            auto tdz_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            auto tdx_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx); 

            auto tdy_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx); 

            auto tdz_0_xxxx_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx); 

            auto tdx_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 1); 

            auto tdy_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 1); 

            auto tdz_0_xxxy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 1); 

            auto tdx_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 2); 

            auto tdy_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 2); 

            auto tdz_0_xxxz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 2); 

            auto tdx_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 3); 

            auto tdy_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 3); 

            auto tdz_0_xxyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 3); 

            auto tdx_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 4); 

            auto tdy_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 4); 

            auto tdz_0_xxyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 4); 

            auto tdx_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 5); 

            auto tdy_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 5); 

            auto tdz_0_xxzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 5); 

            auto tdx_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 6); 

            auto tdy_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 6); 

            auto tdz_0_xyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 6); 

            auto tdx_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 7); 

            auto tdy_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 7); 

            auto tdz_0_xyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 7); 

            auto tdx_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 8); 

            auto tdy_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 8); 

            auto tdz_0_xyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 8); 

            auto tdx_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 9); 

            auto tdy_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 9); 

            auto tdz_0_xzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 9); 

            auto tdx_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 10); 

            auto tdy_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 10); 

            auto tdz_0_yyyy_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 10); 

            auto tdx_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 11); 

            auto tdy_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 11); 

            auto tdz_0_yyyz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 11); 

            auto tdx_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 12); 

            auto tdy_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 12); 

            auto tdz_0_yyzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 12); 

            auto tdx_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 13); 

            auto tdy_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 13); 

            auto tdz_0_yzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 13); 

            auto tdx_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * idx + 14); 

            auto tdy_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 15 * bdim + 15 * idx + 14); 

            auto tdz_0_zzzz_0 = primBuffer.data(pidx_d_0_4_m0 + 30 * bdim + 15 * idx + 14); 

            auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20); 

            auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21); 

            auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22); 

            auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23); 

            auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24); 

            auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24); 

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

            auto ts_z_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 30); 

            auto ts_z_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 31); 

            auto ts_z_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 32); 

            auto ts_z_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 33); 

            auto ts_z_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 34); 

            auto ts_z_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 35); 

            auto ts_z_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 36); 

            auto ts_z_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 37); 

            auto ts_z_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 38); 

            auto ts_z_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 39); 

            auto ts_z_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 40); 

            auto ts_z_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 41); 

            auto ts_z_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 42); 

            auto ts_z_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 43); 

            auto ts_z_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 44); 

            // set up pointers to integrals

            auto tdx_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 75); 

            auto tdy_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 75); 

            auto tdz_zz_xxxx_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 75); 

            auto tdx_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 76); 

            auto tdy_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 76); 

            auto tdz_zz_xxxy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 76); 

            auto tdx_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 77); 

            auto tdy_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 77); 

            auto tdz_zz_xxxz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 77); 

            auto tdx_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 78); 

            auto tdy_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 78); 

            auto tdz_zz_xxyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 78); 

            auto tdx_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 79); 

            auto tdy_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 79); 

            auto tdz_zz_xxyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 79); 

            auto tdx_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 80); 

            auto tdy_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 80); 

            auto tdz_zz_xxzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 80); 

            auto tdx_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 81); 

            auto tdy_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 81); 

            auto tdz_zz_xyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 81); 

            auto tdx_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 82); 

            auto tdy_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 82); 

            auto tdz_zz_xyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 82); 

            auto tdx_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 83); 

            auto tdy_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 83); 

            auto tdz_zz_xyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 83); 

            auto tdx_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 84); 

            auto tdy_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 84); 

            auto tdz_zz_xzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 84); 

            auto tdx_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 85); 

            auto tdy_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 85); 

            auto tdz_zz_yyyy_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 85); 

            auto tdx_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 86); 

            auto tdy_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 86); 

            auto tdz_zz_yyyz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 86); 

            auto tdx_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 87); 

            auto tdy_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 87); 

            auto tdz_zz_yyzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 87); 

            auto tdx_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 88); 

            auto tdy_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 88); 

            auto tdz_zz_yzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 88); 

            auto tdx_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * idx + 89); 

            auto tdy_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 90 * bdim + 90 * idx + 89); 

            auto tdz_zz_zzzz_0 = primBuffer.data(pidx_d_2_4_m0 + 180 * bdim + 90 * idx + 89); 

            // Batch of Integrals (225,270)

            #pragma omp simd aligned(fx, pa_z, tdx_0_xxxx_0, tdx_0_xxxy_0, tdx_0_xxxz_0, tdx_0_xxyy_0, \
                                     tdx_0_xxyz_0, tdx_0_xxzz_0, tdx_0_xyyy_0, tdx_0_xyyz_0, tdx_0_xyzz_0, tdx_0_xzzz_0, \
                                     tdx_0_yyyy_0, tdx_0_yyyz_0, tdx_0_yyzz_0, tdx_0_yzzz_0, tdx_0_zzzz_0, tdx_z_xxx_0, \
                                     tdx_z_xxxx_0, tdx_z_xxxy_0, tdx_z_xxxz_0, tdx_z_xxy_0, tdx_z_xxyy_0, tdx_z_xxyz_0, \
                                     tdx_z_xxz_0, tdx_z_xxzz_0, tdx_z_xyy_0, tdx_z_xyyy_0, tdx_z_xyyz_0, tdx_z_xyz_0, \
                                     tdx_z_xyzz_0, tdx_z_xzz_0, tdx_z_xzzz_0, tdx_z_yyy_0, tdx_z_yyyy_0, tdx_z_yyyz_0, \
                                     tdx_z_yyz_0, tdx_z_yyzz_0, tdx_z_yzz_0, tdx_z_yzzz_0, tdx_z_zzz_0, tdx_z_zzzz_0, \
                                     tdx_zz_xxxx_0, tdx_zz_xxxy_0, tdx_zz_xxxz_0, tdx_zz_xxyy_0, tdx_zz_xxyz_0, \
                                     tdx_zz_xxzz_0, tdx_zz_xyyy_0, tdx_zz_xyyz_0, tdx_zz_xyzz_0, tdx_zz_xzzz_0, \
                                     tdx_zz_yyyy_0, tdx_zz_yyyz_0, tdx_zz_yyzz_0, tdx_zz_yzzz_0, tdx_zz_zzzz_0, \
                                     tdy_0_xxxx_0, tdy_0_xxxy_0, tdy_0_xxxz_0, tdy_0_xxyy_0, tdy_0_xxyz_0, tdy_0_xxzz_0, \
                                     tdy_0_xyyy_0, tdy_0_xyyz_0, tdy_0_xyzz_0, tdy_0_xzzz_0, tdy_0_yyyy_0, tdy_0_yyyz_0, \
                                     tdy_0_yyzz_0, tdy_0_yzzz_0, tdy_0_zzzz_0, tdy_z_xxx_0, tdy_z_xxxx_0, tdy_z_xxxy_0, \
                                     tdy_z_xxxz_0, tdy_z_xxy_0, tdy_z_xxyy_0, tdy_z_xxyz_0, tdy_z_xxz_0, tdy_z_xxzz_0, \
                                     tdy_z_xyy_0, tdy_z_xyyy_0, tdy_z_xyyz_0, tdy_z_xyz_0, tdy_z_xyzz_0, tdy_z_xzz_0, \
                                     tdy_z_xzzz_0, tdy_z_yyy_0, tdy_z_yyyy_0, tdy_z_yyyz_0, tdy_z_yyz_0, tdy_z_yyzz_0, \
                                     tdy_z_yzz_0, tdy_z_yzzz_0, tdy_z_zzz_0, tdy_z_zzzz_0, tdy_zz_xxxx_0, \
                                     tdy_zz_xxxy_0, tdy_zz_xxxz_0, tdy_zz_xxyy_0, tdy_zz_xxyz_0, tdy_zz_xxzz_0, \
                                     tdy_zz_xyyy_0, tdy_zz_xyyz_0, tdy_zz_xyzz_0, tdy_zz_xzzz_0, tdy_zz_yyyy_0, \
                                     tdy_zz_yyyz_0, tdy_zz_yyzz_0, tdy_zz_yzzz_0, tdy_zz_zzzz_0, tdz_0_xxxx_0, \
                                     tdz_0_xxxy_0, tdz_0_xxxz_0, tdz_0_xxyy_0, tdz_0_xxyz_0, tdz_0_xxzz_0, tdz_0_xyyy_0, \
                                     tdz_0_xyyz_0, tdz_0_xyzz_0, tdz_0_xzzz_0, tdz_0_yyyy_0, tdz_0_yyyz_0, tdz_0_yyzz_0, \
                                     tdz_0_yzzz_0, tdz_0_zzzz_0, tdz_z_xxx_0, tdz_z_xxxx_0, tdz_z_xxxy_0, tdz_z_xxxz_0, \
                                     tdz_z_xxy_0, tdz_z_xxyy_0, tdz_z_xxyz_0, tdz_z_xxz_0, tdz_z_xxzz_0, tdz_z_xyy_0, \
                                     tdz_z_xyyy_0, tdz_z_xyyz_0, tdz_z_xyz_0, tdz_z_xyzz_0, tdz_z_xzz_0, tdz_z_xzzz_0, \
                                     tdz_z_yyy_0, tdz_z_yyyy_0, tdz_z_yyyz_0, tdz_z_yyz_0, tdz_z_yyzz_0, tdz_z_yzz_0, \
                                     tdz_z_yzzz_0, tdz_z_zzz_0, tdz_z_zzzz_0, tdz_zz_xxxx_0, tdz_zz_xxxy_0, \
                                     tdz_zz_xxxz_0, tdz_zz_xxyy_0, tdz_zz_xxyz_0, tdz_zz_xxzz_0, tdz_zz_xyyy_0, \
                                     tdz_zz_xyyz_0, tdz_zz_xyzz_0, tdz_zz_xzzz_0, tdz_zz_yyyy_0, tdz_zz_yyyz_0, \
                                     tdz_zz_yyzz_0, tdz_zz_yzzz_0, tdz_zz_zzzz_0, ts_z_xxxx_0, ts_z_xxxy_0, ts_z_xxxz_0, \
                                     ts_z_xxyy_0, ts_z_xxyz_0, ts_z_xxzz_0, ts_z_xyyy_0, ts_z_xyyz_0, ts_z_xyzz_0, \
                                     ts_z_xzzz_0, ts_z_yyyy_0, ts_z_yyyz_0, ts_z_yyzz_0, ts_z_yzzz_0, ts_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_zz_xxxx_0[j] = pa_z[j] * tdx_z_xxxx_0[j] + 0.5 * fl1_fx * tdx_0_xxxx_0[j];

                tdy_zz_xxxx_0[j] = pa_z[j] * tdy_z_xxxx_0[j] + 0.5 * fl1_fx * tdy_0_xxxx_0[j];

                tdz_zz_xxxx_0[j] = pa_z[j] * tdz_z_xxxx_0[j] + 0.5 * fl1_fx * tdz_0_xxxx_0[j] + 0.5 * fl1_fx * ts_z_xxxx_0[j];

                tdx_zz_xxxy_0[j] = pa_z[j] * tdx_z_xxxy_0[j] + 0.5 * fl1_fx * tdx_0_xxxy_0[j];

                tdy_zz_xxxy_0[j] = pa_z[j] * tdy_z_xxxy_0[j] + 0.5 * fl1_fx * tdy_0_xxxy_0[j];

                tdz_zz_xxxy_0[j] = pa_z[j] * tdz_z_xxxy_0[j] + 0.5 * fl1_fx * tdz_0_xxxy_0[j] + 0.5 * fl1_fx * ts_z_xxxy_0[j];

                tdx_zz_xxxz_0[j] = pa_z[j] * tdx_z_xxxz_0[j] + 0.5 * fl1_fx * tdx_0_xxxz_0[j] + 0.5 * fl1_fx * tdx_z_xxx_0[j];

                tdy_zz_xxxz_0[j] = pa_z[j] * tdy_z_xxxz_0[j] + 0.5 * fl1_fx * tdy_0_xxxz_0[j] + 0.5 * fl1_fx * tdy_z_xxx_0[j];

                tdz_zz_xxxz_0[j] = pa_z[j] * tdz_z_xxxz_0[j] + 0.5 * fl1_fx * tdz_0_xxxz_0[j] + 0.5 * fl1_fx * tdz_z_xxx_0[j] + 0.5 * fl1_fx * ts_z_xxxz_0[j];

                tdx_zz_xxyy_0[j] = pa_z[j] * tdx_z_xxyy_0[j] + 0.5 * fl1_fx * tdx_0_xxyy_0[j];

                tdy_zz_xxyy_0[j] = pa_z[j] * tdy_z_xxyy_0[j] + 0.5 * fl1_fx * tdy_0_xxyy_0[j];

                tdz_zz_xxyy_0[j] = pa_z[j] * tdz_z_xxyy_0[j] + 0.5 * fl1_fx * tdz_0_xxyy_0[j] + 0.5 * fl1_fx * ts_z_xxyy_0[j];

                tdx_zz_xxyz_0[j] = pa_z[j] * tdx_z_xxyz_0[j] + 0.5 * fl1_fx * tdx_0_xxyz_0[j] + 0.5 * fl1_fx * tdx_z_xxy_0[j];

                tdy_zz_xxyz_0[j] = pa_z[j] * tdy_z_xxyz_0[j] + 0.5 * fl1_fx * tdy_0_xxyz_0[j] + 0.5 * fl1_fx * tdy_z_xxy_0[j];

                tdz_zz_xxyz_0[j] = pa_z[j] * tdz_z_xxyz_0[j] + 0.5 * fl1_fx * tdz_0_xxyz_0[j] + 0.5 * fl1_fx * tdz_z_xxy_0[j] + 0.5 * fl1_fx * ts_z_xxyz_0[j];

                tdx_zz_xxzz_0[j] = pa_z[j] * tdx_z_xxzz_0[j] + 0.5 * fl1_fx * tdx_0_xxzz_0[j] + fl1_fx * tdx_z_xxz_0[j];

                tdy_zz_xxzz_0[j] = pa_z[j] * tdy_z_xxzz_0[j] + 0.5 * fl1_fx * tdy_0_xxzz_0[j] + fl1_fx * tdy_z_xxz_0[j];

                tdz_zz_xxzz_0[j] = pa_z[j] * tdz_z_xxzz_0[j] + 0.5 * fl1_fx * tdz_0_xxzz_0[j] + fl1_fx * tdz_z_xxz_0[j] + 0.5 * fl1_fx * ts_z_xxzz_0[j];

                tdx_zz_xyyy_0[j] = pa_z[j] * tdx_z_xyyy_0[j] + 0.5 * fl1_fx * tdx_0_xyyy_0[j];

                tdy_zz_xyyy_0[j] = pa_z[j] * tdy_z_xyyy_0[j] + 0.5 * fl1_fx * tdy_0_xyyy_0[j];

                tdz_zz_xyyy_0[j] = pa_z[j] * tdz_z_xyyy_0[j] + 0.5 * fl1_fx * tdz_0_xyyy_0[j] + 0.5 * fl1_fx * ts_z_xyyy_0[j];

                tdx_zz_xyyz_0[j] = pa_z[j] * tdx_z_xyyz_0[j] + 0.5 * fl1_fx * tdx_0_xyyz_0[j] + 0.5 * fl1_fx * tdx_z_xyy_0[j];

                tdy_zz_xyyz_0[j] = pa_z[j] * tdy_z_xyyz_0[j] + 0.5 * fl1_fx * tdy_0_xyyz_0[j] + 0.5 * fl1_fx * tdy_z_xyy_0[j];

                tdz_zz_xyyz_0[j] = pa_z[j] * tdz_z_xyyz_0[j] + 0.5 * fl1_fx * tdz_0_xyyz_0[j] + 0.5 * fl1_fx * tdz_z_xyy_0[j] + 0.5 * fl1_fx * ts_z_xyyz_0[j];

                tdx_zz_xyzz_0[j] = pa_z[j] * tdx_z_xyzz_0[j] + 0.5 * fl1_fx * tdx_0_xyzz_0[j] + fl1_fx * tdx_z_xyz_0[j];

                tdy_zz_xyzz_0[j] = pa_z[j] * tdy_z_xyzz_0[j] + 0.5 * fl1_fx * tdy_0_xyzz_0[j] + fl1_fx * tdy_z_xyz_0[j];

                tdz_zz_xyzz_0[j] = pa_z[j] * tdz_z_xyzz_0[j] + 0.5 * fl1_fx * tdz_0_xyzz_0[j] + fl1_fx * tdz_z_xyz_0[j] + 0.5 * fl1_fx * ts_z_xyzz_0[j];

                tdx_zz_xzzz_0[j] = pa_z[j] * tdx_z_xzzz_0[j] + 0.5 * fl1_fx * tdx_0_xzzz_0[j] + 1.5 * fl1_fx * tdx_z_xzz_0[j];

                tdy_zz_xzzz_0[j] = pa_z[j] * tdy_z_xzzz_0[j] + 0.5 * fl1_fx * tdy_0_xzzz_0[j] + 1.5 * fl1_fx * tdy_z_xzz_0[j];

                tdz_zz_xzzz_0[j] = pa_z[j] * tdz_z_xzzz_0[j] + 0.5 * fl1_fx * tdz_0_xzzz_0[j] + 1.5 * fl1_fx * tdz_z_xzz_0[j] + 0.5 * fl1_fx * ts_z_xzzz_0[j];

                tdx_zz_yyyy_0[j] = pa_z[j] * tdx_z_yyyy_0[j] + 0.5 * fl1_fx * tdx_0_yyyy_0[j];

                tdy_zz_yyyy_0[j] = pa_z[j] * tdy_z_yyyy_0[j] + 0.5 * fl1_fx * tdy_0_yyyy_0[j];

                tdz_zz_yyyy_0[j] = pa_z[j] * tdz_z_yyyy_0[j] + 0.5 * fl1_fx * tdz_0_yyyy_0[j] + 0.5 * fl1_fx * ts_z_yyyy_0[j];

                tdx_zz_yyyz_0[j] = pa_z[j] * tdx_z_yyyz_0[j] + 0.5 * fl1_fx * tdx_0_yyyz_0[j] + 0.5 * fl1_fx * tdx_z_yyy_0[j];

                tdy_zz_yyyz_0[j] = pa_z[j] * tdy_z_yyyz_0[j] + 0.5 * fl1_fx * tdy_0_yyyz_0[j] + 0.5 * fl1_fx * tdy_z_yyy_0[j];

                tdz_zz_yyyz_0[j] = pa_z[j] * tdz_z_yyyz_0[j] + 0.5 * fl1_fx * tdz_0_yyyz_0[j] + 0.5 * fl1_fx * tdz_z_yyy_0[j] + 0.5 * fl1_fx * ts_z_yyyz_0[j];

                tdx_zz_yyzz_0[j] = pa_z[j] * tdx_z_yyzz_0[j] + 0.5 * fl1_fx * tdx_0_yyzz_0[j] + fl1_fx * tdx_z_yyz_0[j];

                tdy_zz_yyzz_0[j] = pa_z[j] * tdy_z_yyzz_0[j] + 0.5 * fl1_fx * tdy_0_yyzz_0[j] + fl1_fx * tdy_z_yyz_0[j];

                tdz_zz_yyzz_0[j] = pa_z[j] * tdz_z_yyzz_0[j] + 0.5 * fl1_fx * tdz_0_yyzz_0[j] + fl1_fx * tdz_z_yyz_0[j] + 0.5 * fl1_fx * ts_z_yyzz_0[j];

                tdx_zz_yzzz_0[j] = pa_z[j] * tdx_z_yzzz_0[j] + 0.5 * fl1_fx * tdx_0_yzzz_0[j] + 1.5 * fl1_fx * tdx_z_yzz_0[j];

                tdy_zz_yzzz_0[j] = pa_z[j] * tdy_z_yzzz_0[j] + 0.5 * fl1_fx * tdy_0_yzzz_0[j] + 1.5 * fl1_fx * tdy_z_yzz_0[j];

                tdz_zz_yzzz_0[j] = pa_z[j] * tdz_z_yzzz_0[j] + 0.5 * fl1_fx * tdz_0_yzzz_0[j] + 1.5 * fl1_fx * tdz_z_yzz_0[j] + 0.5 * fl1_fx * ts_z_yzzz_0[j];

                tdx_zz_zzzz_0[j] = pa_z[j] * tdx_z_zzzz_0[j] + 0.5 * fl1_fx * tdx_0_zzzz_0[j] + 2.0 * fl1_fx * tdx_z_zzz_0[j];

                tdy_zz_zzzz_0[j] = pa_z[j] * tdy_z_zzzz_0[j] + 0.5 * fl1_fx * tdy_0_zzzz_0[j] + 2.0 * fl1_fx * tdy_z_zzz_0[j];

                tdz_zz_zzzz_0[j] = pa_z[j] * tdz_z_zzzz_0[j] + 0.5 * fl1_fx * tdz_0_zzzz_0[j] + 2.0 * fl1_fx * tdz_z_zzz_0[j] + 0.5 * fl1_fx * ts_z_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGD(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForGD_0_45(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  nOSFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForGD_45_90(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   nOSFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        ediprecfunc::compElectricDipoleForGD_90_135(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    nOSFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        ediprecfunc::compElectricDipoleForGD_135_180(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     nOSFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForGD_180_225(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     nOSFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        ediprecfunc::compElectricDipoleForGD_225_270(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     nOSFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 
    }

    void
    compElectricDipoleForGD_0_45(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const int32_t              nOSFactors,
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

        auto pidx_d_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx); 

            auto tdy_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx); 

            auto tdz_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx); 

            auto tdx_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 1); 

            auto tdy_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 1); 

            auto tdz_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 1); 

            auto tdx_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 2); 

            auto tdy_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 2); 

            auto tdz_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 2); 

            auto tdx_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 3); 

            auto tdy_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 3); 

            auto tdz_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 3); 

            auto tdx_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 4); 

            auto tdy_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 4); 

            auto tdz_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 4); 

            auto tdx_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 5); 

            auto tdy_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 5); 

            auto tdz_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 5); 

            auto tdx_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 6); 

            auto tdy_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 6); 

            auto tdz_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 6); 

            auto tdx_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 7); 

            auto tdy_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 7); 

            auto tdz_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 7); 

            auto tdx_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 8); 

            auto tdy_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 8); 

            auto tdz_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 8); 

            auto tdx_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 9); 

            auto tdy_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 9); 

            auto tdz_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 9); 

            auto tdx_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 10); 

            auto tdy_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 10); 

            auto tdz_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 10); 

            auto tdx_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 11); 

            auto tdy_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 11); 

            auto tdz_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 11); 

            auto tdx_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 12); 

            auto tdy_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 12); 

            auto tdz_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 12); 

            auto tdx_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 13); 

            auto tdy_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 13); 

            auto tdz_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 13); 

            auto tdx_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 14); 

            auto tdy_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 14); 

            auto tdz_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 14); 

            auto tdx_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx); 

            auto tdy_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx); 

            auto tdz_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx); 

            auto tdx_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 1); 

            auto tdy_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 1); 

            auto tdz_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 1); 

            auto tdx_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 2); 

            auto tdy_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 2); 

            auto tdz_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 2); 

            auto tdx_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 3); 

            auto tdy_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 3); 

            auto tdz_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 3); 

            auto tdx_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 4); 

            auto tdy_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 4); 

            auto tdz_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 4); 

            auto tdx_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 5); 

            auto tdy_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 5); 

            auto tdz_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 5); 

            auto tdx_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 6); 

            auto tdy_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 6); 

            auto tdz_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 6); 

            auto tdx_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 7); 

            auto tdy_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 7); 

            auto tdz_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 7); 

            auto tdx_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 8); 

            auto tdy_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 8); 

            auto tdz_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 8); 

            auto tdx_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 9); 

            auto tdy_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 9); 

            auto tdz_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 9); 

            auto tdx_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 10); 

            auto tdy_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 10); 

            auto tdz_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 10); 

            auto tdx_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 11); 

            auto tdy_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 11); 

            auto tdz_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 11); 

            auto tdx_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 12); 

            auto tdy_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 12); 

            auto tdz_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 12); 

            auto tdx_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 13); 

            auto tdy_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 13); 

            auto tdz_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 13); 

            auto tdx_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 14); 

            auto tdy_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 14); 

            auto tdz_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 14); 

            auto tdx_xxx_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx); 

            auto tdy_xxx_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx); 

            auto tdz_xxx_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx); 

            auto tdx_xxx_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 1); 

            auto tdy_xxx_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 1); 

            auto tdz_xxx_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 1); 

            auto tdx_xxx_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 2); 

            auto tdy_xxx_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 2); 

            auto tdz_xxx_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 2); 

            auto tdx_xxy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 3); 

            auto tdy_xxy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 3); 

            auto tdz_xxy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 3); 

            auto tdx_xxy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 4); 

            auto tdy_xxy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 4); 

            auto tdz_xxy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 4); 

            auto tdx_xxy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 5); 

            auto tdy_xxy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 5); 

            auto tdz_xxy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 5); 

            auto tdx_xxz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 6); 

            auto tdy_xxz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 6); 

            auto tdz_xxz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 6); 

            auto tdx_xxz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 7); 

            auto tdy_xxz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 7); 

            auto tdz_xxz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 7); 

            auto tdx_xxz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 8); 

            auto tdy_xxz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 8); 

            auto tdz_xxz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 8); 

            auto ts_xxx_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx); 

            auto ts_xxx_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 1); 

            auto ts_xxx_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 2); 

            auto ts_xxx_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 3); 

            auto ts_xxx_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 4); 

            auto ts_xxx_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 5); 

            auto ts_xxy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 6); 

            auto ts_xxy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 7); 

            auto ts_xxy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 8); 

            auto ts_xxy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 9); 

            auto ts_xxy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 10); 

            auto ts_xxy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 11); 

            auto ts_xxz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 12); 

            auto ts_xxz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 13); 

            auto ts_xxz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 14); 

            // set up pointers to integrals

            auto tdx_xxxx_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx); 

            auto tdy_xxxx_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx); 

            auto tdz_xxxx_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx); 

            auto tdx_xxxx_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 1); 

            auto tdy_xxxx_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 1); 

            auto tdz_xxxx_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 1); 

            auto tdx_xxxx_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 2); 

            auto tdy_xxxx_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 2); 

            auto tdz_xxxx_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 2); 

            auto tdx_xxxx_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 3); 

            auto tdy_xxxx_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 3); 

            auto tdz_xxxx_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 3); 

            auto tdx_xxxx_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 4); 

            auto tdy_xxxx_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 4); 

            auto tdz_xxxx_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 4); 

            auto tdx_xxxx_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 5); 

            auto tdy_xxxx_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 5); 

            auto tdz_xxxx_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 5); 

            auto tdx_xxxy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 6); 

            auto tdy_xxxy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 6); 

            auto tdz_xxxy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 6); 

            auto tdx_xxxy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 7); 

            auto tdy_xxxy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 7); 

            auto tdz_xxxy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 7); 

            auto tdx_xxxy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 8); 

            auto tdy_xxxy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 8); 

            auto tdz_xxxy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 8); 

            auto tdx_xxxy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 9); 

            auto tdy_xxxy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 9); 

            auto tdz_xxxy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 9); 

            auto tdx_xxxy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 10); 

            auto tdy_xxxy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 10); 

            auto tdz_xxxy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 10); 

            auto tdx_xxxy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 11); 

            auto tdy_xxxy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 11); 

            auto tdz_xxxy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 11); 

            auto tdx_xxxz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 12); 

            auto tdy_xxxz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 12); 

            auto tdz_xxxz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 12); 

            auto tdx_xxxz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 13); 

            auto tdy_xxxz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 13); 

            auto tdz_xxxz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 13); 

            auto tdx_xxxz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 14); 

            auto tdy_xxxz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 14); 

            auto tdz_xxxz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 14); 

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, tdx_xx_xx_0, tdx_xx_xy_0, tdx_xx_xz_0, tdx_xx_yy_0, \
                                     tdx_xx_yz_0, tdx_xx_zz_0, tdx_xxx_x_0, tdx_xxx_xx_0, tdx_xxx_xy_0, tdx_xxx_xz_0, \
                                     tdx_xxx_y_0, tdx_xxx_yy_0, tdx_xxx_yz_0, tdx_xxx_z_0, tdx_xxx_zz_0, tdx_xxxx_xx_0, \
                                     tdx_xxxx_xy_0, tdx_xxxx_xz_0, tdx_xxxx_yy_0, tdx_xxxx_yz_0, tdx_xxxx_zz_0, \
                                     tdx_xxxy_xx_0, tdx_xxxy_xy_0, tdx_xxxy_xz_0, tdx_xxxy_yy_0, tdx_xxxy_yz_0, \
                                     tdx_xxxy_zz_0, tdx_xxxz_xx_0, tdx_xxxz_xy_0, tdx_xxxz_xz_0, tdx_xxy_x_0, \
                                     tdx_xxy_xx_0, tdx_xxy_xy_0, tdx_xxy_xz_0, tdx_xxy_y_0, tdx_xxy_yy_0, tdx_xxy_yz_0, \
                                     tdx_xxy_z_0, tdx_xxy_zz_0, tdx_xxz_x_0, tdx_xxz_xx_0, tdx_xxz_xy_0, tdx_xxz_xz_0, \
                                     tdx_xxz_y_0, tdx_xxz_z_0, tdx_xy_xx_0, tdx_xy_xy_0, tdx_xy_xz_0, tdx_xy_yy_0, \
                                     tdx_xy_yz_0, tdx_xy_zz_0, tdx_xz_xx_0, tdx_xz_xy_0, tdx_xz_xz_0, tdy_xx_xx_0, \
                                     tdy_xx_xy_0, tdy_xx_xz_0, tdy_xx_yy_0, tdy_xx_yz_0, tdy_xx_zz_0, tdy_xxx_x_0, \
                                     tdy_xxx_xx_0, tdy_xxx_xy_0, tdy_xxx_xz_0, tdy_xxx_y_0, tdy_xxx_yy_0, tdy_xxx_yz_0, \
                                     tdy_xxx_z_0, tdy_xxx_zz_0, tdy_xxxx_xx_0, tdy_xxxx_xy_0, tdy_xxxx_xz_0, \
                                     tdy_xxxx_yy_0, tdy_xxxx_yz_0, tdy_xxxx_zz_0, tdy_xxxy_xx_0, tdy_xxxy_xy_0, \
                                     tdy_xxxy_xz_0, tdy_xxxy_yy_0, tdy_xxxy_yz_0, tdy_xxxy_zz_0, tdy_xxxz_xx_0, \
                                     tdy_xxxz_xy_0, tdy_xxxz_xz_0, tdy_xxy_x_0, tdy_xxy_xx_0, tdy_xxy_xy_0, tdy_xxy_xz_0, \
                                     tdy_xxy_y_0, tdy_xxy_yy_0, tdy_xxy_yz_0, tdy_xxy_z_0, tdy_xxy_zz_0, tdy_xxz_x_0, \
                                     tdy_xxz_xx_0, tdy_xxz_xy_0, tdy_xxz_xz_0, tdy_xxz_y_0, tdy_xxz_z_0, tdy_xy_xx_0, \
                                     tdy_xy_xy_0, tdy_xy_xz_0, tdy_xy_yy_0, tdy_xy_yz_0, tdy_xy_zz_0, tdy_xz_xx_0, \
                                     tdy_xz_xy_0, tdy_xz_xz_0, tdz_xx_xx_0, tdz_xx_xy_0, tdz_xx_xz_0, tdz_xx_yy_0, \
                                     tdz_xx_yz_0, tdz_xx_zz_0, tdz_xxx_x_0, tdz_xxx_xx_0, tdz_xxx_xy_0, tdz_xxx_xz_0, \
                                     tdz_xxx_y_0, tdz_xxx_yy_0, tdz_xxx_yz_0, tdz_xxx_z_0, tdz_xxx_zz_0, tdz_xxxx_xx_0, \
                                     tdz_xxxx_xy_0, tdz_xxxx_xz_0, tdz_xxxx_yy_0, tdz_xxxx_yz_0, tdz_xxxx_zz_0, \
                                     tdz_xxxy_xx_0, tdz_xxxy_xy_0, tdz_xxxy_xz_0, tdz_xxxy_yy_0, tdz_xxxy_yz_0, \
                                     tdz_xxxy_zz_0, tdz_xxxz_xx_0, tdz_xxxz_xy_0, tdz_xxxz_xz_0, tdz_xxy_x_0, \
                                     tdz_xxy_xx_0, tdz_xxy_xy_0, tdz_xxy_xz_0, tdz_xxy_y_0, tdz_xxy_yy_0, tdz_xxy_yz_0, \
                                     tdz_xxy_z_0, tdz_xxy_zz_0, tdz_xxz_x_0, tdz_xxz_xx_0, tdz_xxz_xy_0, tdz_xxz_xz_0, \
                                     tdz_xxz_y_0, tdz_xxz_z_0, tdz_xy_xx_0, tdz_xy_xy_0, tdz_xy_xz_0, tdz_xy_yy_0, \
                                     tdz_xy_yz_0, tdz_xy_zz_0, tdz_xz_xx_0, tdz_xz_xy_0, tdz_xz_xz_0, ts_xxx_xx_0, \
                                     ts_xxx_xy_0, ts_xxx_xz_0, ts_xxx_yy_0, ts_xxx_yz_0, ts_xxx_zz_0, ts_xxy_xx_0, \
                                     ts_xxy_xy_0, ts_xxy_xz_0, ts_xxy_yy_0, ts_xxy_yz_0, ts_xxy_zz_0, ts_xxz_xx_0, \
                                     ts_xxz_xy_0, ts_xxz_xz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxxx_xx_0[j] = pa_x[j] * tdx_xxx_xx_0[j] + 1.5 * fl1_fx * tdx_xx_xx_0[j] + fl1_fx * tdx_xxx_x_0[j] + 0.5 * fl1_fx * ts_xxx_xx_0[j];

                tdy_xxxx_xx_0[j] = pa_x[j] * tdy_xxx_xx_0[j] + 1.5 * fl1_fx * tdy_xx_xx_0[j] + fl1_fx * tdy_xxx_x_0[j];

                tdz_xxxx_xx_0[j] = pa_x[j] * tdz_xxx_xx_0[j] + 1.5 * fl1_fx * tdz_xx_xx_0[j] + fl1_fx * tdz_xxx_x_0[j];

                tdx_xxxx_xy_0[j] = pa_x[j] * tdx_xxx_xy_0[j] + 1.5 * fl1_fx * tdx_xx_xy_0[j] + 0.5 * fl1_fx * tdx_xxx_y_0[j] + 0.5 * fl1_fx * ts_xxx_xy_0[j];

                tdy_xxxx_xy_0[j] = pa_x[j] * tdy_xxx_xy_0[j] + 1.5 * fl1_fx * tdy_xx_xy_0[j] + 0.5 * fl1_fx * tdy_xxx_y_0[j];

                tdz_xxxx_xy_0[j] = pa_x[j] * tdz_xxx_xy_0[j] + 1.5 * fl1_fx * tdz_xx_xy_0[j] + 0.5 * fl1_fx * tdz_xxx_y_0[j];

                tdx_xxxx_xz_0[j] = pa_x[j] * tdx_xxx_xz_0[j] + 1.5 * fl1_fx * tdx_xx_xz_0[j] + 0.5 * fl1_fx * tdx_xxx_z_0[j] + 0.5 * fl1_fx * ts_xxx_xz_0[j];

                tdy_xxxx_xz_0[j] = pa_x[j] * tdy_xxx_xz_0[j] + 1.5 * fl1_fx * tdy_xx_xz_0[j] + 0.5 * fl1_fx * tdy_xxx_z_0[j];

                tdz_xxxx_xz_0[j] = pa_x[j] * tdz_xxx_xz_0[j] + 1.5 * fl1_fx * tdz_xx_xz_0[j] + 0.5 * fl1_fx * tdz_xxx_z_0[j];

                tdx_xxxx_yy_0[j] = pa_x[j] * tdx_xxx_yy_0[j] + 1.5 * fl1_fx * tdx_xx_yy_0[j] + 0.5 * fl1_fx * ts_xxx_yy_0[j];

                tdy_xxxx_yy_0[j] = pa_x[j] * tdy_xxx_yy_0[j] + 1.5 * fl1_fx * tdy_xx_yy_0[j];

                tdz_xxxx_yy_0[j] = pa_x[j] * tdz_xxx_yy_0[j] + 1.5 * fl1_fx * tdz_xx_yy_0[j];

                tdx_xxxx_yz_0[j] = pa_x[j] * tdx_xxx_yz_0[j] + 1.5 * fl1_fx * tdx_xx_yz_0[j] + 0.5 * fl1_fx * ts_xxx_yz_0[j];

                tdy_xxxx_yz_0[j] = pa_x[j] * tdy_xxx_yz_0[j] + 1.5 * fl1_fx * tdy_xx_yz_0[j];

                tdz_xxxx_yz_0[j] = pa_x[j] * tdz_xxx_yz_0[j] + 1.5 * fl1_fx * tdz_xx_yz_0[j];

                tdx_xxxx_zz_0[j] = pa_x[j] * tdx_xxx_zz_0[j] + 1.5 * fl1_fx * tdx_xx_zz_0[j] + 0.5 * fl1_fx * ts_xxx_zz_0[j];

                tdy_xxxx_zz_0[j] = pa_x[j] * tdy_xxx_zz_0[j] + 1.5 * fl1_fx * tdy_xx_zz_0[j];

                tdz_xxxx_zz_0[j] = pa_x[j] * tdz_xxx_zz_0[j] + 1.5 * fl1_fx * tdz_xx_zz_0[j];

                tdx_xxxy_xx_0[j] = pa_x[j] * tdx_xxy_xx_0[j] + fl1_fx * tdx_xy_xx_0[j] + fl1_fx * tdx_xxy_x_0[j] + 0.5 * fl1_fx * ts_xxy_xx_0[j];

                tdy_xxxy_xx_0[j] = pa_x[j] * tdy_xxy_xx_0[j] + fl1_fx * tdy_xy_xx_0[j] + fl1_fx * tdy_xxy_x_0[j];

                tdz_xxxy_xx_0[j] = pa_x[j] * tdz_xxy_xx_0[j] + fl1_fx * tdz_xy_xx_0[j] + fl1_fx * tdz_xxy_x_0[j];

                tdx_xxxy_xy_0[j] = pa_x[j] * tdx_xxy_xy_0[j] + fl1_fx * tdx_xy_xy_0[j] + 0.5 * fl1_fx * tdx_xxy_y_0[j] + 0.5 * fl1_fx * ts_xxy_xy_0[j];

                tdy_xxxy_xy_0[j] = pa_x[j] * tdy_xxy_xy_0[j] + fl1_fx * tdy_xy_xy_0[j] + 0.5 * fl1_fx * tdy_xxy_y_0[j];

                tdz_xxxy_xy_0[j] = pa_x[j] * tdz_xxy_xy_0[j] + fl1_fx * tdz_xy_xy_0[j] + 0.5 * fl1_fx * tdz_xxy_y_0[j];

                tdx_xxxy_xz_0[j] = pa_x[j] * tdx_xxy_xz_0[j] + fl1_fx * tdx_xy_xz_0[j] + 0.5 * fl1_fx * tdx_xxy_z_0[j] + 0.5 * fl1_fx * ts_xxy_xz_0[j];

                tdy_xxxy_xz_0[j] = pa_x[j] * tdy_xxy_xz_0[j] + fl1_fx * tdy_xy_xz_0[j] + 0.5 * fl1_fx * tdy_xxy_z_0[j];

                tdz_xxxy_xz_0[j] = pa_x[j] * tdz_xxy_xz_0[j] + fl1_fx * tdz_xy_xz_0[j] + 0.5 * fl1_fx * tdz_xxy_z_0[j];

                tdx_xxxy_yy_0[j] = pa_x[j] * tdx_xxy_yy_0[j] + fl1_fx * tdx_xy_yy_0[j] + 0.5 * fl1_fx * ts_xxy_yy_0[j];

                tdy_xxxy_yy_0[j] = pa_x[j] * tdy_xxy_yy_0[j] + fl1_fx * tdy_xy_yy_0[j];

                tdz_xxxy_yy_0[j] = pa_x[j] * tdz_xxy_yy_0[j] + fl1_fx * tdz_xy_yy_0[j];

                tdx_xxxy_yz_0[j] = pa_x[j] * tdx_xxy_yz_0[j] + fl1_fx * tdx_xy_yz_0[j] + 0.5 * fl1_fx * ts_xxy_yz_0[j];

                tdy_xxxy_yz_0[j] = pa_x[j] * tdy_xxy_yz_0[j] + fl1_fx * tdy_xy_yz_0[j];

                tdz_xxxy_yz_0[j] = pa_x[j] * tdz_xxy_yz_0[j] + fl1_fx * tdz_xy_yz_0[j];

                tdx_xxxy_zz_0[j] = pa_x[j] * tdx_xxy_zz_0[j] + fl1_fx * tdx_xy_zz_0[j] + 0.5 * fl1_fx * ts_xxy_zz_0[j];

                tdy_xxxy_zz_0[j] = pa_x[j] * tdy_xxy_zz_0[j] + fl1_fx * tdy_xy_zz_0[j];

                tdz_xxxy_zz_0[j] = pa_x[j] * tdz_xxy_zz_0[j] + fl1_fx * tdz_xy_zz_0[j];

                tdx_xxxz_xx_0[j] = pa_x[j] * tdx_xxz_xx_0[j] + fl1_fx * tdx_xz_xx_0[j] + fl1_fx * tdx_xxz_x_0[j] + 0.5 * fl1_fx * ts_xxz_xx_0[j];

                tdy_xxxz_xx_0[j] = pa_x[j] * tdy_xxz_xx_0[j] + fl1_fx * tdy_xz_xx_0[j] + fl1_fx * tdy_xxz_x_0[j];

                tdz_xxxz_xx_0[j] = pa_x[j] * tdz_xxz_xx_0[j] + fl1_fx * tdz_xz_xx_0[j] + fl1_fx * tdz_xxz_x_0[j];

                tdx_xxxz_xy_0[j] = pa_x[j] * tdx_xxz_xy_0[j] + fl1_fx * tdx_xz_xy_0[j] + 0.5 * fl1_fx * tdx_xxz_y_0[j] + 0.5 * fl1_fx * ts_xxz_xy_0[j];

                tdy_xxxz_xy_0[j] = pa_x[j] * tdy_xxz_xy_0[j] + fl1_fx * tdy_xz_xy_0[j] + 0.5 * fl1_fx * tdy_xxz_y_0[j];

                tdz_xxxz_xy_0[j] = pa_x[j] * tdz_xxz_xy_0[j] + fl1_fx * tdz_xz_xy_0[j] + 0.5 * fl1_fx * tdz_xxz_y_0[j];

                tdx_xxxz_xz_0[j] = pa_x[j] * tdx_xxz_xz_0[j] + fl1_fx * tdx_xz_xz_0[j] + 0.5 * fl1_fx * tdx_xxz_z_0[j] + 0.5 * fl1_fx * ts_xxz_xz_0[j];

                tdy_xxxz_xz_0[j] = pa_x[j] * tdy_xxz_xz_0[j] + fl1_fx * tdy_xz_xz_0[j] + 0.5 * fl1_fx * tdy_xxz_z_0[j];

                tdz_xxxz_xz_0[j] = pa_x[j] * tdz_xxz_xz_0[j] + fl1_fx * tdz_xz_xz_0[j] + 0.5 * fl1_fx * tdz_xxz_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGD_45_90(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const int32_t              nOSFactors,
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

        auto pidx_d_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 15); 

            auto tdy_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 15); 

            auto tdz_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 15); 

            auto tdx_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 16); 

            auto tdy_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 16); 

            auto tdz_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 16); 

            auto tdx_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 17); 

            auto tdy_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 17); 

            auto tdz_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 17); 

            auto tdx_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 18); 

            auto tdy_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 18); 

            auto tdz_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 18); 

            auto tdx_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 19); 

            auto tdy_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 19); 

            auto tdz_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 19); 

            auto tdx_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 20); 

            auto tdy_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 20); 

            auto tdz_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 20); 

            auto tdx_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 21); 

            auto tdy_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 21); 

            auto tdz_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 21); 

            auto tdx_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 22); 

            auto tdy_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 22); 

            auto tdz_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 22); 

            auto tdx_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 23); 

            auto tdy_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 23); 

            auto tdz_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 23); 

            auto tdx_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 24); 

            auto tdy_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 24); 

            auto tdz_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 24); 

            auto tdx_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 25); 

            auto tdy_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 25); 

            auto tdz_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 25); 

            auto tdx_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 26); 

            auto tdy_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 26); 

            auto tdz_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 26); 

            auto tdx_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 27); 

            auto tdy_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 27); 

            auto tdz_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 27); 

            auto tdx_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 28); 

            auto tdy_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 28); 

            auto tdz_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 28); 

            auto tdx_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 29); 

            auto tdy_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 29); 

            auto tdz_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 29); 

            auto tdx_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 15); 

            auto tdy_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 15); 

            auto tdz_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 15); 

            auto tdx_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 16); 

            auto tdy_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 16); 

            auto tdz_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 16); 

            auto tdx_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 17); 

            auto tdy_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 17); 

            auto tdz_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 17); 

            auto tdx_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 18); 

            auto tdy_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tdz_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tdx_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 19); 

            auto tdy_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tdz_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tdx_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 20); 

            auto tdy_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tdz_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tdx_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 21); 

            auto tdy_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tdz_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tdx_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 22); 

            auto tdy_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tdz_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tdx_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 23); 

            auto tdy_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tdz_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tdx_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 24); 

            auto tdy_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tdz_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tdx_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 25); 

            auto tdy_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tdz_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tdx_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 26); 

            auto tdy_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tdz_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tdx_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 27); 

            auto tdy_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tdz_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tdx_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 28); 

            auto tdy_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tdz_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tdx_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 29); 

            auto tdy_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 29); 

            auto tdz_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tdx_xyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 9); 

            auto tdy_xyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 9); 

            auto tdz_xyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 9); 

            auto tdx_xyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 10); 

            auto tdy_xyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 10); 

            auto tdz_xyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 10); 

            auto tdx_xyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 11); 

            auto tdy_xyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 11); 

            auto tdz_xyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 11); 

            auto tdx_xyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 12); 

            auto tdy_xyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 12); 

            auto tdz_xyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 12); 

            auto tdx_xyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 13); 

            auto tdy_xyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 13); 

            auto tdz_xyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 13); 

            auto tdx_xyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 14); 

            auto tdy_xyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 14); 

            auto tdz_xyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 14); 

            auto ts_xxz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 15); 

            auto ts_xxz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 16); 

            auto ts_xxz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 17); 

            auto ts_xyy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 18); 

            auto ts_xyy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 19); 

            auto ts_xyy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 20); 

            auto ts_xyy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 21); 

            auto ts_xyy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 22); 

            auto ts_xyy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 23); 

            auto ts_xyz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 24); 

            auto ts_xyz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 25); 

            auto ts_xyz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 26); 

            auto ts_xyz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 27); 

            auto ts_xyz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 28); 

            auto ts_xyz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 29); 

            // set up pointers to integrals

            auto tdx_xxxz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 15); 

            auto tdy_xxxz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 15); 

            auto tdz_xxxz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 15); 

            auto tdx_xxxz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 16); 

            auto tdy_xxxz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 16); 

            auto tdz_xxxz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 16); 

            auto tdx_xxxz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 17); 

            auto tdy_xxxz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 17); 

            auto tdz_xxxz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 17); 

            auto tdx_xxyy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 18); 

            auto tdy_xxyy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 18); 

            auto tdz_xxyy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 18); 

            auto tdx_xxyy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 19); 

            auto tdy_xxyy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 19); 

            auto tdz_xxyy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 19); 

            auto tdx_xxyy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 20); 

            auto tdy_xxyy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 20); 

            auto tdz_xxyy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 20); 

            auto tdx_xxyy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 21); 

            auto tdy_xxyy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 21); 

            auto tdz_xxyy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 21); 

            auto tdx_xxyy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 22); 

            auto tdy_xxyy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 22); 

            auto tdz_xxyy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 22); 

            auto tdx_xxyy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 23); 

            auto tdy_xxyy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 23); 

            auto tdz_xxyy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 23); 

            auto tdx_xxyz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 24); 

            auto tdy_xxyz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 24); 

            auto tdz_xxyz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 24); 

            auto tdx_xxyz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 25); 

            auto tdy_xxyz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 25); 

            auto tdz_xxyz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 25); 

            auto tdx_xxyz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 26); 

            auto tdy_xxyz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 26); 

            auto tdz_xxyz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 26); 

            auto tdx_xxyz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 27); 

            auto tdy_xxyz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 27); 

            auto tdz_xxyz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 27); 

            auto tdx_xxyz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 28); 

            auto tdy_xxyz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 28); 

            auto tdz_xxyz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 28); 

            auto tdx_xxyz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 29); 

            auto tdy_xxyz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 29); 

            auto tdz_xxyz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 29); 

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, tdx_xxxz_yy_0, tdx_xxxz_yz_0, tdx_xxxz_zz_0, tdx_xxyy_xx_0, \
                                     tdx_xxyy_xy_0, tdx_xxyy_xz_0, tdx_xxyy_yy_0, tdx_xxyy_yz_0, tdx_xxyy_zz_0, \
                                     tdx_xxyz_xx_0, tdx_xxyz_xy_0, tdx_xxyz_xz_0, tdx_xxyz_yy_0, tdx_xxyz_yz_0, \
                                     tdx_xxyz_zz_0, tdx_xxz_yy_0, tdx_xxz_yz_0, tdx_xxz_zz_0, tdx_xyy_x_0, tdx_xyy_xx_0, \
                                     tdx_xyy_xy_0, tdx_xyy_xz_0, tdx_xyy_y_0, tdx_xyy_yy_0, tdx_xyy_yz_0, tdx_xyy_z_0, \
                                     tdx_xyy_zz_0, tdx_xyz_x_0, tdx_xyz_xx_0, tdx_xyz_xy_0, tdx_xyz_xz_0, tdx_xyz_y_0, \
                                     tdx_xyz_yy_0, tdx_xyz_yz_0, tdx_xyz_z_0, tdx_xyz_zz_0, tdx_xz_yy_0, tdx_xz_yz_0, \
                                     tdx_xz_zz_0, tdx_yy_xx_0, tdx_yy_xy_0, tdx_yy_xz_0, tdx_yy_yy_0, tdx_yy_yz_0, \
                                     tdx_yy_zz_0, tdx_yz_xx_0, tdx_yz_xy_0, tdx_yz_xz_0, tdx_yz_yy_0, tdx_yz_yz_0, \
                                     tdx_yz_zz_0, tdy_xxxz_yy_0, tdy_xxxz_yz_0, tdy_xxxz_zz_0, tdy_xxyy_xx_0, \
                                     tdy_xxyy_xy_0, tdy_xxyy_xz_0, tdy_xxyy_yy_0, tdy_xxyy_yz_0, tdy_xxyy_zz_0, \
                                     tdy_xxyz_xx_0, tdy_xxyz_xy_0, tdy_xxyz_xz_0, tdy_xxyz_yy_0, tdy_xxyz_yz_0, \
                                     tdy_xxyz_zz_0, tdy_xxz_yy_0, tdy_xxz_yz_0, tdy_xxz_zz_0, tdy_xyy_x_0, tdy_xyy_xx_0, \
                                     tdy_xyy_xy_0, tdy_xyy_xz_0, tdy_xyy_y_0, tdy_xyy_yy_0, tdy_xyy_yz_0, tdy_xyy_z_0, \
                                     tdy_xyy_zz_0, tdy_xyz_x_0, tdy_xyz_xx_0, tdy_xyz_xy_0, tdy_xyz_xz_0, tdy_xyz_y_0, \
                                     tdy_xyz_yy_0, tdy_xyz_yz_0, tdy_xyz_z_0, tdy_xyz_zz_0, tdy_xz_yy_0, tdy_xz_yz_0, \
                                     tdy_xz_zz_0, tdy_yy_xx_0, tdy_yy_xy_0, tdy_yy_xz_0, tdy_yy_yy_0, tdy_yy_yz_0, \
                                     tdy_yy_zz_0, tdy_yz_xx_0, tdy_yz_xy_0, tdy_yz_xz_0, tdy_yz_yy_0, tdy_yz_yz_0, \
                                     tdy_yz_zz_0, tdz_xxxz_yy_0, tdz_xxxz_yz_0, tdz_xxxz_zz_0, tdz_xxyy_xx_0, \
                                     tdz_xxyy_xy_0, tdz_xxyy_xz_0, tdz_xxyy_yy_0, tdz_xxyy_yz_0, tdz_xxyy_zz_0, \
                                     tdz_xxyz_xx_0, tdz_xxyz_xy_0, tdz_xxyz_xz_0, tdz_xxyz_yy_0, tdz_xxyz_yz_0, \
                                     tdz_xxyz_zz_0, tdz_xxz_yy_0, tdz_xxz_yz_0, tdz_xxz_zz_0, tdz_xyy_x_0, tdz_xyy_xx_0, \
                                     tdz_xyy_xy_0, tdz_xyy_xz_0, tdz_xyy_y_0, tdz_xyy_yy_0, tdz_xyy_yz_0, tdz_xyy_z_0, \
                                     tdz_xyy_zz_0, tdz_xyz_x_0, tdz_xyz_xx_0, tdz_xyz_xy_0, tdz_xyz_xz_0, tdz_xyz_y_0, \
                                     tdz_xyz_yy_0, tdz_xyz_yz_0, tdz_xyz_z_0, tdz_xyz_zz_0, tdz_xz_yy_0, tdz_xz_yz_0, \
                                     tdz_xz_zz_0, tdz_yy_xx_0, tdz_yy_xy_0, tdz_yy_xz_0, tdz_yy_yy_0, tdz_yy_yz_0, \
                                     tdz_yy_zz_0, tdz_yz_xx_0, tdz_yz_xy_0, tdz_yz_xz_0, tdz_yz_yy_0, tdz_yz_yz_0, \
                                     tdz_yz_zz_0, ts_xxz_yy_0, ts_xxz_yz_0, ts_xxz_zz_0, ts_xyy_xx_0, ts_xyy_xy_0, \
                                     ts_xyy_xz_0, ts_xyy_yy_0, ts_xyy_yz_0, ts_xyy_zz_0, ts_xyz_xx_0, ts_xyz_xy_0, \
                                     ts_xyz_xz_0, ts_xyz_yy_0, ts_xyz_yz_0, ts_xyz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxxz_yy_0[j] = pa_x[j] * tdx_xxz_yy_0[j] + fl1_fx * tdx_xz_yy_0[j] + 0.5 * fl1_fx * ts_xxz_yy_0[j];

                tdy_xxxz_yy_0[j] = pa_x[j] * tdy_xxz_yy_0[j] + fl1_fx * tdy_xz_yy_0[j];

                tdz_xxxz_yy_0[j] = pa_x[j] * tdz_xxz_yy_0[j] + fl1_fx * tdz_xz_yy_0[j];

                tdx_xxxz_yz_0[j] = pa_x[j] * tdx_xxz_yz_0[j] + fl1_fx * tdx_xz_yz_0[j] + 0.5 * fl1_fx * ts_xxz_yz_0[j];

                tdy_xxxz_yz_0[j] = pa_x[j] * tdy_xxz_yz_0[j] + fl1_fx * tdy_xz_yz_0[j];

                tdz_xxxz_yz_0[j] = pa_x[j] * tdz_xxz_yz_0[j] + fl1_fx * tdz_xz_yz_0[j];

                tdx_xxxz_zz_0[j] = pa_x[j] * tdx_xxz_zz_0[j] + fl1_fx * tdx_xz_zz_0[j] + 0.5 * fl1_fx * ts_xxz_zz_0[j];

                tdy_xxxz_zz_0[j] = pa_x[j] * tdy_xxz_zz_0[j] + fl1_fx * tdy_xz_zz_0[j];

                tdz_xxxz_zz_0[j] = pa_x[j] * tdz_xxz_zz_0[j] + fl1_fx * tdz_xz_zz_0[j];

                tdx_xxyy_xx_0[j] = pa_x[j] * tdx_xyy_xx_0[j] + 0.5 * fl1_fx * tdx_yy_xx_0[j] + fl1_fx * tdx_xyy_x_0[j] + 0.5 * fl1_fx * ts_xyy_xx_0[j];

                tdy_xxyy_xx_0[j] = pa_x[j] * tdy_xyy_xx_0[j] + 0.5 * fl1_fx * tdy_yy_xx_0[j] + fl1_fx * tdy_xyy_x_0[j];

                tdz_xxyy_xx_0[j] = pa_x[j] * tdz_xyy_xx_0[j] + 0.5 * fl1_fx * tdz_yy_xx_0[j] + fl1_fx * tdz_xyy_x_0[j];

                tdx_xxyy_xy_0[j] = pa_x[j] * tdx_xyy_xy_0[j] + 0.5 * fl1_fx * tdx_yy_xy_0[j] + 0.5 * fl1_fx * tdx_xyy_y_0[j] + 0.5 * fl1_fx * ts_xyy_xy_0[j];

                tdy_xxyy_xy_0[j] = pa_x[j] * tdy_xyy_xy_0[j] + 0.5 * fl1_fx * tdy_yy_xy_0[j] + 0.5 * fl1_fx * tdy_xyy_y_0[j];

                tdz_xxyy_xy_0[j] = pa_x[j] * tdz_xyy_xy_0[j] + 0.5 * fl1_fx * tdz_yy_xy_0[j] + 0.5 * fl1_fx * tdz_xyy_y_0[j];

                tdx_xxyy_xz_0[j] = pa_x[j] * tdx_xyy_xz_0[j] + 0.5 * fl1_fx * tdx_yy_xz_0[j] + 0.5 * fl1_fx * tdx_xyy_z_0[j] + 0.5 * fl1_fx * ts_xyy_xz_0[j];

                tdy_xxyy_xz_0[j] = pa_x[j] * tdy_xyy_xz_0[j] + 0.5 * fl1_fx * tdy_yy_xz_0[j] + 0.5 * fl1_fx * tdy_xyy_z_0[j];

                tdz_xxyy_xz_0[j] = pa_x[j] * tdz_xyy_xz_0[j] + 0.5 * fl1_fx * tdz_yy_xz_0[j] + 0.5 * fl1_fx * tdz_xyy_z_0[j];

                tdx_xxyy_yy_0[j] = pa_x[j] * tdx_xyy_yy_0[j] + 0.5 * fl1_fx * tdx_yy_yy_0[j] + 0.5 * fl1_fx * ts_xyy_yy_0[j];

                tdy_xxyy_yy_0[j] = pa_x[j] * tdy_xyy_yy_0[j] + 0.5 * fl1_fx * tdy_yy_yy_0[j];

                tdz_xxyy_yy_0[j] = pa_x[j] * tdz_xyy_yy_0[j] + 0.5 * fl1_fx * tdz_yy_yy_0[j];

                tdx_xxyy_yz_0[j] = pa_x[j] * tdx_xyy_yz_0[j] + 0.5 * fl1_fx * tdx_yy_yz_0[j] + 0.5 * fl1_fx * ts_xyy_yz_0[j];

                tdy_xxyy_yz_0[j] = pa_x[j] * tdy_xyy_yz_0[j] + 0.5 * fl1_fx * tdy_yy_yz_0[j];

                tdz_xxyy_yz_0[j] = pa_x[j] * tdz_xyy_yz_0[j] + 0.5 * fl1_fx * tdz_yy_yz_0[j];

                tdx_xxyy_zz_0[j] = pa_x[j] * tdx_xyy_zz_0[j] + 0.5 * fl1_fx * tdx_yy_zz_0[j] + 0.5 * fl1_fx * ts_xyy_zz_0[j];

                tdy_xxyy_zz_0[j] = pa_x[j] * tdy_xyy_zz_0[j] + 0.5 * fl1_fx * tdy_yy_zz_0[j];

                tdz_xxyy_zz_0[j] = pa_x[j] * tdz_xyy_zz_0[j] + 0.5 * fl1_fx * tdz_yy_zz_0[j];

                tdx_xxyz_xx_0[j] = pa_x[j] * tdx_xyz_xx_0[j] + 0.5 * fl1_fx * tdx_yz_xx_0[j] + fl1_fx * tdx_xyz_x_0[j] + 0.5 * fl1_fx * ts_xyz_xx_0[j];

                tdy_xxyz_xx_0[j] = pa_x[j] * tdy_xyz_xx_0[j] + 0.5 * fl1_fx * tdy_yz_xx_0[j] + fl1_fx * tdy_xyz_x_0[j];

                tdz_xxyz_xx_0[j] = pa_x[j] * tdz_xyz_xx_0[j] + 0.5 * fl1_fx * tdz_yz_xx_0[j] + fl1_fx * tdz_xyz_x_0[j];

                tdx_xxyz_xy_0[j] = pa_x[j] * tdx_xyz_xy_0[j] + 0.5 * fl1_fx * tdx_yz_xy_0[j] + 0.5 * fl1_fx * tdx_xyz_y_0[j] + 0.5 * fl1_fx * ts_xyz_xy_0[j];

                tdy_xxyz_xy_0[j] = pa_x[j] * tdy_xyz_xy_0[j] + 0.5 * fl1_fx * tdy_yz_xy_0[j] + 0.5 * fl1_fx * tdy_xyz_y_0[j];

                tdz_xxyz_xy_0[j] = pa_x[j] * tdz_xyz_xy_0[j] + 0.5 * fl1_fx * tdz_yz_xy_0[j] + 0.5 * fl1_fx * tdz_xyz_y_0[j];

                tdx_xxyz_xz_0[j] = pa_x[j] * tdx_xyz_xz_0[j] + 0.5 * fl1_fx * tdx_yz_xz_0[j] + 0.5 * fl1_fx * tdx_xyz_z_0[j] + 0.5 * fl1_fx * ts_xyz_xz_0[j];

                tdy_xxyz_xz_0[j] = pa_x[j] * tdy_xyz_xz_0[j] + 0.5 * fl1_fx * tdy_yz_xz_0[j] + 0.5 * fl1_fx * tdy_xyz_z_0[j];

                tdz_xxyz_xz_0[j] = pa_x[j] * tdz_xyz_xz_0[j] + 0.5 * fl1_fx * tdz_yz_xz_0[j] + 0.5 * fl1_fx * tdz_xyz_z_0[j];

                tdx_xxyz_yy_0[j] = pa_x[j] * tdx_xyz_yy_0[j] + 0.5 * fl1_fx * tdx_yz_yy_0[j] + 0.5 * fl1_fx * ts_xyz_yy_0[j];

                tdy_xxyz_yy_0[j] = pa_x[j] * tdy_xyz_yy_0[j] + 0.5 * fl1_fx * tdy_yz_yy_0[j];

                tdz_xxyz_yy_0[j] = pa_x[j] * tdz_xyz_yy_0[j] + 0.5 * fl1_fx * tdz_yz_yy_0[j];

                tdx_xxyz_yz_0[j] = pa_x[j] * tdx_xyz_yz_0[j] + 0.5 * fl1_fx * tdx_yz_yz_0[j] + 0.5 * fl1_fx * ts_xyz_yz_0[j];

                tdy_xxyz_yz_0[j] = pa_x[j] * tdy_xyz_yz_0[j] + 0.5 * fl1_fx * tdy_yz_yz_0[j];

                tdz_xxyz_yz_0[j] = pa_x[j] * tdz_xyz_yz_0[j] + 0.5 * fl1_fx * tdz_yz_yz_0[j];

                tdx_xxyz_zz_0[j] = pa_x[j] * tdx_xyz_zz_0[j] + 0.5 * fl1_fx * tdx_yz_zz_0[j] + 0.5 * fl1_fx * ts_xyz_zz_0[j];

                tdy_xxyz_zz_0[j] = pa_x[j] * tdy_xyz_zz_0[j] + 0.5 * fl1_fx * tdy_yz_zz_0[j];

                tdz_xxyz_zz_0[j] = pa_x[j] * tdz_xyz_zz_0[j] + 0.5 * fl1_fx * tdz_yz_zz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGD_90_135(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const int32_t              nOSFactors,
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

        auto pidx_d_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 30); 

            auto tdy_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 30); 

            auto tdz_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 30); 

            auto tdx_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 31); 

            auto tdy_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 31); 

            auto tdz_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 31); 

            auto tdx_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 32); 

            auto tdy_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 32); 

            auto tdz_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 32); 

            auto tdx_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 33); 

            auto tdy_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 33); 

            auto tdz_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 33); 

            auto tdx_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 34); 

            auto tdy_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 34); 

            auto tdz_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 34); 

            auto tdx_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 35); 

            auto tdy_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 35); 

            auto tdz_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 35); 

            auto tdx_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 36); 

            auto tdy_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tdz_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdx_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 37); 

            auto tdy_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tdz_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdx_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 38); 

            auto tdy_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tdz_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdx_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 39); 

            auto tdy_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tdz_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdx_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 40); 

            auto tdy_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tdz_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdx_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 41); 

            auto tdy_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tdz_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tdx_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 42); 

            auto tdy_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tdz_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdx_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 43); 

            auto tdy_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tdz_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdx_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 44); 

            auto tdy_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tdz_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 44); 

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

            auto tdx_xzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 15); 

            auto tdy_xzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 15); 

            auto tdz_xzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 15); 

            auto tdx_xzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 16); 

            auto tdy_xzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 16); 

            auto tdz_xzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 16); 

            auto tdx_xzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 17); 

            auto tdy_xzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 17); 

            auto tdz_xzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 17); 

            auto tdx_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 18); 

            auto tdy_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 18); 

            auto tdz_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 18); 

            auto tdx_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 19); 

            auto tdy_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 19); 

            auto tdz_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 19); 

            auto tdx_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 20); 

            auto tdy_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdz_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdx_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 21); 

            auto tdy_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdz_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdx_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 22); 

            auto tdy_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdz_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdx_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 23); 

            auto tdy_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdz_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 23); 

            auto ts_xzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 30); 

            auto ts_xzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 31); 

            auto ts_xzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 32); 

            auto ts_xzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 33); 

            auto ts_xzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 34); 

            auto ts_xzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 35); 

            auto ts_yyy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 36); 

            auto ts_yyy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 37); 

            auto ts_yyy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 38); 

            auto ts_yyy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 39); 

            auto ts_yyy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 40); 

            auto ts_yyy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 41); 

            auto ts_yyz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 42); 

            auto ts_yyz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 43); 

            auto ts_yyz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 44); 

            // set up pointers to integrals

            auto tdx_xxzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 30); 

            auto tdy_xxzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 30); 

            auto tdz_xxzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 30); 

            auto tdx_xxzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 31); 

            auto tdy_xxzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 31); 

            auto tdz_xxzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 31); 

            auto tdx_xxzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 32); 

            auto tdy_xxzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 32); 

            auto tdz_xxzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 32); 

            auto tdx_xxzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 33); 

            auto tdy_xxzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 33); 

            auto tdz_xxzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 33); 

            auto tdx_xxzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 34); 

            auto tdy_xxzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 34); 

            auto tdz_xxzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 34); 

            auto tdx_xxzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 35); 

            auto tdy_xxzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 35); 

            auto tdz_xxzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 35); 

            auto tdx_xyyy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 36); 

            auto tdy_xyyy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 36); 

            auto tdz_xyyy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 36); 

            auto tdx_xyyy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 37); 

            auto tdy_xyyy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 37); 

            auto tdz_xyyy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 37); 

            auto tdx_xyyy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 38); 

            auto tdy_xyyy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 38); 

            auto tdz_xyyy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 38); 

            auto tdx_xyyy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 39); 

            auto tdy_xyyy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 39); 

            auto tdz_xyyy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 39); 

            auto tdx_xyyy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 40); 

            auto tdy_xyyy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 40); 

            auto tdz_xyyy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 40); 

            auto tdx_xyyy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 41); 

            auto tdy_xyyy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 41); 

            auto tdz_xyyy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 41); 

            auto tdx_xyyz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 42); 

            auto tdy_xyyz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 42); 

            auto tdz_xyyz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 42); 

            auto tdx_xyyz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 43); 

            auto tdy_xyyz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 43); 

            auto tdz_xyyz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 43); 

            auto tdx_xyyz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 44); 

            auto tdy_xyyz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 44); 

            auto tdz_xyyz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 44); 

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fx, pa_x, tdx_xxzz_xx_0, tdx_xxzz_xy_0, tdx_xxzz_xz_0, tdx_xxzz_yy_0, \
                                     tdx_xxzz_yz_0, tdx_xxzz_zz_0, tdx_xyyy_xx_0, tdx_xyyy_xy_0, tdx_xyyy_xz_0, \
                                     tdx_xyyy_yy_0, tdx_xyyy_yz_0, tdx_xyyy_zz_0, tdx_xyyz_xx_0, tdx_xyyz_xy_0, \
                                     tdx_xyyz_xz_0, tdx_xzz_x_0, tdx_xzz_xx_0, tdx_xzz_xy_0, tdx_xzz_xz_0, tdx_xzz_y_0, \
                                     tdx_xzz_yy_0, tdx_xzz_yz_0, tdx_xzz_z_0, tdx_xzz_zz_0, tdx_yyy_x_0, tdx_yyy_xx_0, \
                                     tdx_yyy_xy_0, tdx_yyy_xz_0, tdx_yyy_y_0, tdx_yyy_yy_0, tdx_yyy_yz_0, tdx_yyy_z_0, \
                                     tdx_yyy_zz_0, tdx_yyz_x_0, tdx_yyz_xx_0, tdx_yyz_xy_0, tdx_yyz_xz_0, tdx_yyz_y_0, \
                                     tdx_yyz_z_0, tdx_zz_xx_0, tdx_zz_xy_0, tdx_zz_xz_0, tdx_zz_yy_0, tdx_zz_yz_0, \
                                     tdx_zz_zz_0, tdy_xxzz_xx_0, tdy_xxzz_xy_0, tdy_xxzz_xz_0, tdy_xxzz_yy_0, \
                                     tdy_xxzz_yz_0, tdy_xxzz_zz_0, tdy_xyyy_xx_0, tdy_xyyy_xy_0, tdy_xyyy_xz_0, \
                                     tdy_xyyy_yy_0, tdy_xyyy_yz_0, tdy_xyyy_zz_0, tdy_xyyz_xx_0, tdy_xyyz_xy_0, \
                                     tdy_xyyz_xz_0, tdy_xzz_x_0, tdy_xzz_xx_0, tdy_xzz_xy_0, tdy_xzz_xz_0, tdy_xzz_y_0, \
                                     tdy_xzz_yy_0, tdy_xzz_yz_0, tdy_xzz_z_0, tdy_xzz_zz_0, tdy_yyy_x_0, tdy_yyy_xx_0, \
                                     tdy_yyy_xy_0, tdy_yyy_xz_0, tdy_yyy_y_0, tdy_yyy_yy_0, tdy_yyy_yz_0, tdy_yyy_z_0, \
                                     tdy_yyy_zz_0, tdy_yyz_x_0, tdy_yyz_xx_0, tdy_yyz_xy_0, tdy_yyz_xz_0, tdy_yyz_y_0, \
                                     tdy_yyz_z_0, tdy_zz_xx_0, tdy_zz_xy_0, tdy_zz_xz_0, tdy_zz_yy_0, tdy_zz_yz_0, \
                                     tdy_zz_zz_0, tdz_xxzz_xx_0, tdz_xxzz_xy_0, tdz_xxzz_xz_0, tdz_xxzz_yy_0, \
                                     tdz_xxzz_yz_0, tdz_xxzz_zz_0, tdz_xyyy_xx_0, tdz_xyyy_xy_0, tdz_xyyy_xz_0, \
                                     tdz_xyyy_yy_0, tdz_xyyy_yz_0, tdz_xyyy_zz_0, tdz_xyyz_xx_0, tdz_xyyz_xy_0, \
                                     tdz_xyyz_xz_0, tdz_xzz_x_0, tdz_xzz_xx_0, tdz_xzz_xy_0, tdz_xzz_xz_0, tdz_xzz_y_0, \
                                     tdz_xzz_yy_0, tdz_xzz_yz_0, tdz_xzz_z_0, tdz_xzz_zz_0, tdz_yyy_x_0, tdz_yyy_xx_0, \
                                     tdz_yyy_xy_0, tdz_yyy_xz_0, tdz_yyy_y_0, tdz_yyy_yy_0, tdz_yyy_yz_0, tdz_yyy_z_0, \
                                     tdz_yyy_zz_0, tdz_yyz_x_0, tdz_yyz_xx_0, tdz_yyz_xy_0, tdz_yyz_xz_0, tdz_yyz_y_0, \
                                     tdz_yyz_z_0, tdz_zz_xx_0, tdz_zz_xy_0, tdz_zz_xz_0, tdz_zz_yy_0, tdz_zz_yz_0, \
                                     tdz_zz_zz_0, ts_xzz_xx_0, ts_xzz_xy_0, ts_xzz_xz_0, ts_xzz_yy_0, ts_xzz_yz_0, \
                                     ts_xzz_zz_0, ts_yyy_xx_0, ts_yyy_xy_0, ts_yyy_xz_0, ts_yyy_yy_0, ts_yyy_yz_0, \
                                     ts_yyy_zz_0, ts_yyz_xx_0, ts_yyz_xy_0, ts_yyz_xz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxzz_xx_0[j] = pa_x[j] * tdx_xzz_xx_0[j] + 0.5 * fl1_fx * tdx_zz_xx_0[j] + fl1_fx * tdx_xzz_x_0[j] + 0.5 * fl1_fx * ts_xzz_xx_0[j];

                tdy_xxzz_xx_0[j] = pa_x[j] * tdy_xzz_xx_0[j] + 0.5 * fl1_fx * tdy_zz_xx_0[j] + fl1_fx * tdy_xzz_x_0[j];

                tdz_xxzz_xx_0[j] = pa_x[j] * tdz_xzz_xx_0[j] + 0.5 * fl1_fx * tdz_zz_xx_0[j] + fl1_fx * tdz_xzz_x_0[j];

                tdx_xxzz_xy_0[j] = pa_x[j] * tdx_xzz_xy_0[j] + 0.5 * fl1_fx * tdx_zz_xy_0[j] + 0.5 * fl1_fx * tdx_xzz_y_0[j] + 0.5 * fl1_fx * ts_xzz_xy_0[j];

                tdy_xxzz_xy_0[j] = pa_x[j] * tdy_xzz_xy_0[j] + 0.5 * fl1_fx * tdy_zz_xy_0[j] + 0.5 * fl1_fx * tdy_xzz_y_0[j];

                tdz_xxzz_xy_0[j] = pa_x[j] * tdz_xzz_xy_0[j] + 0.5 * fl1_fx * tdz_zz_xy_0[j] + 0.5 * fl1_fx * tdz_xzz_y_0[j];

                tdx_xxzz_xz_0[j] = pa_x[j] * tdx_xzz_xz_0[j] + 0.5 * fl1_fx * tdx_zz_xz_0[j] + 0.5 * fl1_fx * tdx_xzz_z_0[j] + 0.5 * fl1_fx * ts_xzz_xz_0[j];

                tdy_xxzz_xz_0[j] = pa_x[j] * tdy_xzz_xz_0[j] + 0.5 * fl1_fx * tdy_zz_xz_0[j] + 0.5 * fl1_fx * tdy_xzz_z_0[j];

                tdz_xxzz_xz_0[j] = pa_x[j] * tdz_xzz_xz_0[j] + 0.5 * fl1_fx * tdz_zz_xz_0[j] + 0.5 * fl1_fx * tdz_xzz_z_0[j];

                tdx_xxzz_yy_0[j] = pa_x[j] * tdx_xzz_yy_0[j] + 0.5 * fl1_fx * tdx_zz_yy_0[j] + 0.5 * fl1_fx * ts_xzz_yy_0[j];

                tdy_xxzz_yy_0[j] = pa_x[j] * tdy_xzz_yy_0[j] + 0.5 * fl1_fx * tdy_zz_yy_0[j];

                tdz_xxzz_yy_0[j] = pa_x[j] * tdz_xzz_yy_0[j] + 0.5 * fl1_fx * tdz_zz_yy_0[j];

                tdx_xxzz_yz_0[j] = pa_x[j] * tdx_xzz_yz_0[j] + 0.5 * fl1_fx * tdx_zz_yz_0[j] + 0.5 * fl1_fx * ts_xzz_yz_0[j];

                tdy_xxzz_yz_0[j] = pa_x[j] * tdy_xzz_yz_0[j] + 0.5 * fl1_fx * tdy_zz_yz_0[j];

                tdz_xxzz_yz_0[j] = pa_x[j] * tdz_xzz_yz_0[j] + 0.5 * fl1_fx * tdz_zz_yz_0[j];

                tdx_xxzz_zz_0[j] = pa_x[j] * tdx_xzz_zz_0[j] + 0.5 * fl1_fx * tdx_zz_zz_0[j] + 0.5 * fl1_fx * ts_xzz_zz_0[j];

                tdy_xxzz_zz_0[j] = pa_x[j] * tdy_xzz_zz_0[j] + 0.5 * fl1_fx * tdy_zz_zz_0[j];

                tdz_xxzz_zz_0[j] = pa_x[j] * tdz_xzz_zz_0[j] + 0.5 * fl1_fx * tdz_zz_zz_0[j];

                tdx_xyyy_xx_0[j] = pa_x[j] * tdx_yyy_xx_0[j] + fl1_fx * tdx_yyy_x_0[j] + 0.5 * fl1_fx * ts_yyy_xx_0[j];

                tdy_xyyy_xx_0[j] = pa_x[j] * tdy_yyy_xx_0[j] + fl1_fx * tdy_yyy_x_0[j];

                tdz_xyyy_xx_0[j] = pa_x[j] * tdz_yyy_xx_0[j] + fl1_fx * tdz_yyy_x_0[j];

                tdx_xyyy_xy_0[j] = pa_x[j] * tdx_yyy_xy_0[j] + 0.5 * fl1_fx * tdx_yyy_y_0[j] + 0.5 * fl1_fx * ts_yyy_xy_0[j];

                tdy_xyyy_xy_0[j] = pa_x[j] * tdy_yyy_xy_0[j] + 0.5 * fl1_fx * tdy_yyy_y_0[j];

                tdz_xyyy_xy_0[j] = pa_x[j] * tdz_yyy_xy_0[j] + 0.5 * fl1_fx * tdz_yyy_y_0[j];

                tdx_xyyy_xz_0[j] = pa_x[j] * tdx_yyy_xz_0[j] + 0.5 * fl1_fx * tdx_yyy_z_0[j] + 0.5 * fl1_fx * ts_yyy_xz_0[j];

                tdy_xyyy_xz_0[j] = pa_x[j] * tdy_yyy_xz_0[j] + 0.5 * fl1_fx * tdy_yyy_z_0[j];

                tdz_xyyy_xz_0[j] = pa_x[j] * tdz_yyy_xz_0[j] + 0.5 * fl1_fx * tdz_yyy_z_0[j];

                tdx_xyyy_yy_0[j] = pa_x[j] * tdx_yyy_yy_0[j] + 0.5 * fl1_fx * ts_yyy_yy_0[j];

                tdy_xyyy_yy_0[j] = pa_x[j] * tdy_yyy_yy_0[j];

                tdz_xyyy_yy_0[j] = pa_x[j] * tdz_yyy_yy_0[j];

                tdx_xyyy_yz_0[j] = pa_x[j] * tdx_yyy_yz_0[j] + 0.5 * fl1_fx * ts_yyy_yz_0[j];

                tdy_xyyy_yz_0[j] = pa_x[j] * tdy_yyy_yz_0[j];

                tdz_xyyy_yz_0[j] = pa_x[j] * tdz_yyy_yz_0[j];

                tdx_xyyy_zz_0[j] = pa_x[j] * tdx_yyy_zz_0[j] + 0.5 * fl1_fx * ts_yyy_zz_0[j];

                tdy_xyyy_zz_0[j] = pa_x[j] * tdy_yyy_zz_0[j];

                tdz_xyyy_zz_0[j] = pa_x[j] * tdz_yyy_zz_0[j];

                tdx_xyyz_xx_0[j] = pa_x[j] * tdx_yyz_xx_0[j] + fl1_fx * tdx_yyz_x_0[j] + 0.5 * fl1_fx * ts_yyz_xx_0[j];

                tdy_xyyz_xx_0[j] = pa_x[j] * tdy_yyz_xx_0[j] + fl1_fx * tdy_yyz_x_0[j];

                tdz_xyyz_xx_0[j] = pa_x[j] * tdz_yyz_xx_0[j] + fl1_fx * tdz_yyz_x_0[j];

                tdx_xyyz_xy_0[j] = pa_x[j] * tdx_yyz_xy_0[j] + 0.5 * fl1_fx * tdx_yyz_y_0[j] + 0.5 * fl1_fx * ts_yyz_xy_0[j];

                tdy_xyyz_xy_0[j] = pa_x[j] * tdy_yyz_xy_0[j] + 0.5 * fl1_fx * tdy_yyz_y_0[j];

                tdz_xyyz_xy_0[j] = pa_x[j] * tdz_yyz_xy_0[j] + 0.5 * fl1_fx * tdz_yyz_y_0[j];

                tdx_xyyz_xz_0[j] = pa_x[j] * tdx_yyz_xz_0[j] + 0.5 * fl1_fx * tdx_yyz_z_0[j] + 0.5 * fl1_fx * ts_yyz_xz_0[j];

                tdy_xyyz_xz_0[j] = pa_x[j] * tdy_yyz_xz_0[j] + 0.5 * fl1_fx * tdy_yyz_z_0[j];

                tdz_xyyz_xz_0[j] = pa_x[j] * tdz_yyz_xz_0[j] + 0.5 * fl1_fx * tdz_yyz_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGD_135_180(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
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

        auto pidx_d_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tdx_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 45); 

            auto tdy_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto tdz_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdx_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 46); 

            auto tdy_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tdz_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdx_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 47); 

            auto tdy_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tdz_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdx_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 48); 

            auto tdy_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tdz_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdx_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 49); 

            auto tdy_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tdz_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdx_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 50); 

            auto tdy_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 51); 

            auto tdy_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 52); 

            auto tdy_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 53); 

            auto tdy_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 53); 

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

            auto tdx_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 24); 

            auto tdy_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdz_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 24); 

            auto tdx_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 25); 

            auto tdy_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 25); 

            auto tdz_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 25); 

            auto tdx_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 26); 

            auto tdy_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 26); 

            auto tdz_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 26); 

            auto tdx_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 27); 

            auto tdy_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 27); 

            auto tdz_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 27); 

            auto tdx_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 28); 

            auto tdy_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 28); 

            auto tdz_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 28); 

            auto tdx_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 29); 

            auto tdy_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 29); 

            auto tdz_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 29); 

            auto ts_yyz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 45); 

            auto ts_yyz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 46); 

            auto ts_yyz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 47); 

            auto ts_yzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 48); 

            auto ts_yzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 49); 

            auto ts_yzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 50); 

            auto ts_yzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 51); 

            auto ts_yzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 52); 

            auto ts_yzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 53); 

            auto ts_zzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 54); 

            auto ts_zzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 55); 

            auto ts_zzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 56); 

            auto ts_zzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 57); 

            auto ts_zzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 58); 

            auto ts_zzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 59); 

            // set up pointers to integrals

            auto tdx_xyyz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 45); 

            auto tdy_xyyz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 45); 

            auto tdz_xyyz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 45); 

            auto tdx_xyyz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 46); 

            auto tdy_xyyz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 46); 

            auto tdz_xyyz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 46); 

            auto tdx_xyyz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 47); 

            auto tdy_xyyz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 47); 

            auto tdz_xyyz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 47); 

            auto tdx_xyzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 48); 

            auto tdy_xyzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 48); 

            auto tdz_xyzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 48); 

            auto tdx_xyzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 49); 

            auto tdy_xyzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 49); 

            auto tdz_xyzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 49); 

            auto tdx_xyzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 50); 

            auto tdy_xyzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 50); 

            auto tdz_xyzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 50); 

            auto tdx_xyzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 51); 

            auto tdy_xyzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 51); 

            auto tdz_xyzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 51); 

            auto tdx_xyzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 52); 

            auto tdy_xyzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 52); 

            auto tdz_xyzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 52); 

            auto tdx_xyzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 53); 

            auto tdy_xyzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 53); 

            auto tdz_xyzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 53); 

            auto tdx_xzzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 54); 

            auto tdy_xzzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 54); 

            auto tdz_xzzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 54); 

            auto tdx_xzzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 55); 

            auto tdy_xzzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 55); 

            auto tdz_xzzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 55); 

            auto tdx_xzzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 56); 

            auto tdy_xzzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 56); 

            auto tdz_xzzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 56); 

            auto tdx_xzzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 57); 

            auto tdy_xzzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 57); 

            auto tdz_xzzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 57); 

            auto tdx_xzzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 58); 

            auto tdy_xzzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 58); 

            auto tdz_xzzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 58); 

            auto tdx_xzzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 59); 

            auto tdy_xzzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 59); 

            auto tdz_xzzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 59); 

            // Batch of Integrals (135,180)

            #pragma omp simd aligned(fx, pa_x, tdx_xyyz_yy_0, tdx_xyyz_yz_0, tdx_xyyz_zz_0, tdx_xyzz_xx_0, \
                                     tdx_xyzz_xy_0, tdx_xyzz_xz_0, tdx_xyzz_yy_0, tdx_xyzz_yz_0, tdx_xyzz_zz_0, \
                                     tdx_xzzz_xx_0, tdx_xzzz_xy_0, tdx_xzzz_xz_0, tdx_xzzz_yy_0, tdx_xzzz_yz_0, \
                                     tdx_xzzz_zz_0, tdx_yyz_yy_0, tdx_yyz_yz_0, tdx_yyz_zz_0, tdx_yzz_x_0, tdx_yzz_xx_0, \
                                     tdx_yzz_xy_0, tdx_yzz_xz_0, tdx_yzz_y_0, tdx_yzz_yy_0, tdx_yzz_yz_0, tdx_yzz_z_0, \
                                     tdx_yzz_zz_0, tdx_zzz_x_0, tdx_zzz_xx_0, tdx_zzz_xy_0, tdx_zzz_xz_0, tdx_zzz_y_0, \
                                     tdx_zzz_yy_0, tdx_zzz_yz_0, tdx_zzz_z_0, tdx_zzz_zz_0, tdy_xyyz_yy_0, \
                                     tdy_xyyz_yz_0, tdy_xyyz_zz_0, tdy_xyzz_xx_0, tdy_xyzz_xy_0, tdy_xyzz_xz_0, \
                                     tdy_xyzz_yy_0, tdy_xyzz_yz_0, tdy_xyzz_zz_0, tdy_xzzz_xx_0, tdy_xzzz_xy_0, \
                                     tdy_xzzz_xz_0, tdy_xzzz_yy_0, tdy_xzzz_yz_0, tdy_xzzz_zz_0, tdy_yyz_yy_0, \
                                     tdy_yyz_yz_0, tdy_yyz_zz_0, tdy_yzz_x_0, tdy_yzz_xx_0, tdy_yzz_xy_0, tdy_yzz_xz_0, \
                                     tdy_yzz_y_0, tdy_yzz_yy_0, tdy_yzz_yz_0, tdy_yzz_z_0, tdy_yzz_zz_0, tdy_zzz_x_0, \
                                     tdy_zzz_xx_0, tdy_zzz_xy_0, tdy_zzz_xz_0, tdy_zzz_y_0, tdy_zzz_yy_0, tdy_zzz_yz_0, \
                                     tdy_zzz_z_0, tdy_zzz_zz_0, tdz_xyyz_yy_0, tdz_xyyz_yz_0, tdz_xyyz_zz_0, \
                                     tdz_xyzz_xx_0, tdz_xyzz_xy_0, tdz_xyzz_xz_0, tdz_xyzz_yy_0, tdz_xyzz_yz_0, \
                                     tdz_xyzz_zz_0, tdz_xzzz_xx_0, tdz_xzzz_xy_0, tdz_xzzz_xz_0, tdz_xzzz_yy_0, \
                                     tdz_xzzz_yz_0, tdz_xzzz_zz_0, tdz_yyz_yy_0, tdz_yyz_yz_0, tdz_yyz_zz_0, tdz_yzz_x_0, \
                                     tdz_yzz_xx_0, tdz_yzz_xy_0, tdz_yzz_xz_0, tdz_yzz_y_0, tdz_yzz_yy_0, tdz_yzz_yz_0, \
                                     tdz_yzz_z_0, tdz_yzz_zz_0, tdz_zzz_x_0, tdz_zzz_xx_0, tdz_zzz_xy_0, tdz_zzz_xz_0, \
                                     tdz_zzz_y_0, tdz_zzz_yy_0, tdz_zzz_yz_0, tdz_zzz_z_0, tdz_zzz_zz_0, ts_yyz_yy_0, \
                                     ts_yyz_yz_0, ts_yyz_zz_0, ts_yzz_xx_0, ts_yzz_xy_0, ts_yzz_xz_0, ts_yzz_yy_0, \
                                     ts_yzz_yz_0, ts_yzz_zz_0, ts_zzz_xx_0, ts_zzz_xy_0, ts_zzz_xz_0, ts_zzz_yy_0, \
                                     ts_zzz_yz_0, ts_zzz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xyyz_yy_0[j] = pa_x[j] * tdx_yyz_yy_0[j] + 0.5 * fl1_fx * ts_yyz_yy_0[j];

                tdy_xyyz_yy_0[j] = pa_x[j] * tdy_yyz_yy_0[j];

                tdz_xyyz_yy_0[j] = pa_x[j] * tdz_yyz_yy_0[j];

                tdx_xyyz_yz_0[j] = pa_x[j] * tdx_yyz_yz_0[j] + 0.5 * fl1_fx * ts_yyz_yz_0[j];

                tdy_xyyz_yz_0[j] = pa_x[j] * tdy_yyz_yz_0[j];

                tdz_xyyz_yz_0[j] = pa_x[j] * tdz_yyz_yz_0[j];

                tdx_xyyz_zz_0[j] = pa_x[j] * tdx_yyz_zz_0[j] + 0.5 * fl1_fx * ts_yyz_zz_0[j];

                tdy_xyyz_zz_0[j] = pa_x[j] * tdy_yyz_zz_0[j];

                tdz_xyyz_zz_0[j] = pa_x[j] * tdz_yyz_zz_0[j];

                tdx_xyzz_xx_0[j] = pa_x[j] * tdx_yzz_xx_0[j] + fl1_fx * tdx_yzz_x_0[j] + 0.5 * fl1_fx * ts_yzz_xx_0[j];

                tdy_xyzz_xx_0[j] = pa_x[j] * tdy_yzz_xx_0[j] + fl1_fx * tdy_yzz_x_0[j];

                tdz_xyzz_xx_0[j] = pa_x[j] * tdz_yzz_xx_0[j] + fl1_fx * tdz_yzz_x_0[j];

                tdx_xyzz_xy_0[j] = pa_x[j] * tdx_yzz_xy_0[j] + 0.5 * fl1_fx * tdx_yzz_y_0[j] + 0.5 * fl1_fx * ts_yzz_xy_0[j];

                tdy_xyzz_xy_0[j] = pa_x[j] * tdy_yzz_xy_0[j] + 0.5 * fl1_fx * tdy_yzz_y_0[j];

                tdz_xyzz_xy_0[j] = pa_x[j] * tdz_yzz_xy_0[j] + 0.5 * fl1_fx * tdz_yzz_y_0[j];

                tdx_xyzz_xz_0[j] = pa_x[j] * tdx_yzz_xz_0[j] + 0.5 * fl1_fx * tdx_yzz_z_0[j] + 0.5 * fl1_fx * ts_yzz_xz_0[j];

                tdy_xyzz_xz_0[j] = pa_x[j] * tdy_yzz_xz_0[j] + 0.5 * fl1_fx * tdy_yzz_z_0[j];

                tdz_xyzz_xz_0[j] = pa_x[j] * tdz_yzz_xz_0[j] + 0.5 * fl1_fx * tdz_yzz_z_0[j];

                tdx_xyzz_yy_0[j] = pa_x[j] * tdx_yzz_yy_0[j] + 0.5 * fl1_fx * ts_yzz_yy_0[j];

                tdy_xyzz_yy_0[j] = pa_x[j] * tdy_yzz_yy_0[j];

                tdz_xyzz_yy_0[j] = pa_x[j] * tdz_yzz_yy_0[j];

                tdx_xyzz_yz_0[j] = pa_x[j] * tdx_yzz_yz_0[j] + 0.5 * fl1_fx * ts_yzz_yz_0[j];

                tdy_xyzz_yz_0[j] = pa_x[j] * tdy_yzz_yz_0[j];

                tdz_xyzz_yz_0[j] = pa_x[j] * tdz_yzz_yz_0[j];

                tdx_xyzz_zz_0[j] = pa_x[j] * tdx_yzz_zz_0[j] + 0.5 * fl1_fx * ts_yzz_zz_0[j];

                tdy_xyzz_zz_0[j] = pa_x[j] * tdy_yzz_zz_0[j];

                tdz_xyzz_zz_0[j] = pa_x[j] * tdz_yzz_zz_0[j];

                tdx_xzzz_xx_0[j] = pa_x[j] * tdx_zzz_xx_0[j] + fl1_fx * tdx_zzz_x_0[j] + 0.5 * fl1_fx * ts_zzz_xx_0[j];

                tdy_xzzz_xx_0[j] = pa_x[j] * tdy_zzz_xx_0[j] + fl1_fx * tdy_zzz_x_0[j];

                tdz_xzzz_xx_0[j] = pa_x[j] * tdz_zzz_xx_0[j] + fl1_fx * tdz_zzz_x_0[j];

                tdx_xzzz_xy_0[j] = pa_x[j] * tdx_zzz_xy_0[j] + 0.5 * fl1_fx * tdx_zzz_y_0[j] + 0.5 * fl1_fx * ts_zzz_xy_0[j];

                tdy_xzzz_xy_0[j] = pa_x[j] * tdy_zzz_xy_0[j] + 0.5 * fl1_fx * tdy_zzz_y_0[j];

                tdz_xzzz_xy_0[j] = pa_x[j] * tdz_zzz_xy_0[j] + 0.5 * fl1_fx * tdz_zzz_y_0[j];

                tdx_xzzz_xz_0[j] = pa_x[j] * tdx_zzz_xz_0[j] + 0.5 * fl1_fx * tdx_zzz_z_0[j] + 0.5 * fl1_fx * ts_zzz_xz_0[j];

                tdy_xzzz_xz_0[j] = pa_x[j] * tdy_zzz_xz_0[j] + 0.5 * fl1_fx * tdy_zzz_z_0[j];

                tdz_xzzz_xz_0[j] = pa_x[j] * tdz_zzz_xz_0[j] + 0.5 * fl1_fx * tdz_zzz_z_0[j];

                tdx_xzzz_yy_0[j] = pa_x[j] * tdx_zzz_yy_0[j] + 0.5 * fl1_fx * ts_zzz_yy_0[j];

                tdy_xzzz_yy_0[j] = pa_x[j] * tdy_zzz_yy_0[j];

                tdz_xzzz_yy_0[j] = pa_x[j] * tdz_zzz_yy_0[j];

                tdx_xzzz_yz_0[j] = pa_x[j] * tdx_zzz_yz_0[j] + 0.5 * fl1_fx * ts_zzz_yz_0[j];

                tdy_xzzz_yz_0[j] = pa_x[j] * tdy_zzz_yz_0[j];

                tdz_xzzz_yz_0[j] = pa_x[j] * tdz_zzz_yz_0[j];

                tdx_xzzz_zz_0[j] = pa_x[j] * tdx_zzz_zz_0[j] + 0.5 * fl1_fx * ts_zzz_zz_0[j];

                tdy_xzzz_zz_0[j] = pa_x[j] * tdy_zzz_zz_0[j];

                tdz_xzzz_zz_0[j] = pa_x[j] * tdz_zzz_zz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGD_180_225(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
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

        auto pidx_d_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tdx_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 36); 

            auto tdy_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tdz_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdx_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 37); 

            auto tdy_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tdz_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdx_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 38); 

            auto tdy_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tdz_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdx_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 39); 

            auto tdy_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tdz_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdx_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 40); 

            auto tdy_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tdz_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdx_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 41); 

            auto tdy_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tdz_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tdx_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 42); 

            auto tdy_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tdz_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdx_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 43); 

            auto tdy_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tdz_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdx_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 44); 

            auto tdy_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tdz_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            auto tdx_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 45); 

            auto tdy_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto tdz_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdx_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 46); 

            auto tdy_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tdz_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdx_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 47); 

            auto tdy_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tdz_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdx_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 48); 

            auto tdy_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tdz_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdx_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 49); 

            auto tdy_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tdz_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdx_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 50); 

            auto tdy_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 18); 

            auto tdy_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tdz_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tdx_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 19); 

            auto tdy_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tdz_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tdx_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 20); 

            auto tdy_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tdz_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tdx_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 21); 

            auto tdy_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tdz_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tdx_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 22); 

            auto tdy_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tdz_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tdx_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 23); 

            auto tdy_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tdz_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tdx_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 24); 

            auto tdy_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tdz_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tdx_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 25); 

            auto tdy_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tdz_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tdx_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 26); 

            auto tdy_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tdz_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tdx_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 27); 

            auto tdy_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tdz_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tdx_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 28); 

            auto tdy_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tdz_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tdx_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 29); 

            auto tdy_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 29); 

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

            auto tdx_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 18); 

            auto tdy_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 18); 

            auto tdz_yyy_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 18); 

            auto tdx_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 19); 

            auto tdy_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 19); 

            auto tdz_yyy_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 19); 

            auto tdx_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 20); 

            auto tdy_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdz_yyy_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdx_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 21); 

            auto tdy_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdz_yyz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdx_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 22); 

            auto tdy_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdz_yyz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdx_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 23); 

            auto tdy_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdz_yyz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 23); 

            auto tdx_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 24); 

            auto tdy_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdz_yzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 24); 

            auto ts_yyy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 36); 

            auto ts_yyy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 37); 

            auto ts_yyy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 38); 

            auto ts_yyy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 39); 

            auto ts_yyy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 40); 

            auto ts_yyy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 41); 

            auto ts_yyz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 42); 

            auto ts_yyz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 43); 

            auto ts_yyz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 44); 

            auto ts_yyz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 45); 

            auto ts_yyz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 46); 

            auto ts_yyz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 47); 

            auto ts_yzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 48); 

            auto ts_yzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 49); 

            auto ts_yzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 50); 

            // set up pointers to integrals

            auto tdx_yyyy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 60); 

            auto tdy_yyyy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 60); 

            auto tdz_yyyy_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 60); 

            auto tdx_yyyy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 61); 

            auto tdy_yyyy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 61); 

            auto tdz_yyyy_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 61); 

            auto tdx_yyyy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 62); 

            auto tdy_yyyy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 62); 

            auto tdz_yyyy_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 62); 

            auto tdx_yyyy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 63); 

            auto tdy_yyyy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 63); 

            auto tdz_yyyy_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 63); 

            auto tdx_yyyy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 64); 

            auto tdy_yyyy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 64); 

            auto tdz_yyyy_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 64); 

            auto tdx_yyyy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 65); 

            auto tdy_yyyy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 65); 

            auto tdz_yyyy_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 65); 

            auto tdx_yyyz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 66); 

            auto tdy_yyyz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 66); 

            auto tdz_yyyz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 66); 

            auto tdx_yyyz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 67); 

            auto tdy_yyyz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 67); 

            auto tdz_yyyz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 67); 

            auto tdx_yyyz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 68); 

            auto tdy_yyyz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 68); 

            auto tdz_yyyz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 68); 

            auto tdx_yyyz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 69); 

            auto tdy_yyyz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 69); 

            auto tdz_yyyz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 69); 

            auto tdx_yyyz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 70); 

            auto tdy_yyyz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 70); 

            auto tdz_yyyz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 70); 

            auto tdx_yyyz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 71); 

            auto tdy_yyyz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 71); 

            auto tdz_yyyz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 71); 

            auto tdx_yyzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 72); 

            auto tdy_yyzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 72); 

            auto tdz_yyzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 72); 

            auto tdx_yyzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 73); 

            auto tdy_yyzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 73); 

            auto tdz_yyzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 73); 

            auto tdx_yyzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 74); 

            auto tdy_yyzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 74); 

            auto tdz_yyzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 74); 

            // Batch of Integrals (180,225)

            #pragma omp simd aligned(fx, pa_y, tdx_yy_xx_0, tdx_yy_xy_0, tdx_yy_xz_0, tdx_yy_yy_0, \
                                     tdx_yy_yz_0, tdx_yy_zz_0, tdx_yyy_x_0, tdx_yyy_xx_0, tdx_yyy_xy_0, tdx_yyy_xz_0, \
                                     tdx_yyy_y_0, tdx_yyy_yy_0, tdx_yyy_yz_0, tdx_yyy_z_0, tdx_yyy_zz_0, tdx_yyyy_xx_0, \
                                     tdx_yyyy_xy_0, tdx_yyyy_xz_0, tdx_yyyy_yy_0, tdx_yyyy_yz_0, tdx_yyyy_zz_0, \
                                     tdx_yyyz_xx_0, tdx_yyyz_xy_0, tdx_yyyz_xz_0, tdx_yyyz_yy_0, tdx_yyyz_yz_0, \
                                     tdx_yyyz_zz_0, tdx_yyz_x_0, tdx_yyz_xx_0, tdx_yyz_xy_0, tdx_yyz_xz_0, tdx_yyz_y_0, \
                                     tdx_yyz_yy_0, tdx_yyz_yz_0, tdx_yyz_z_0, tdx_yyz_zz_0, tdx_yyzz_xx_0, \
                                     tdx_yyzz_xy_0, tdx_yyzz_xz_0, tdx_yz_xx_0, tdx_yz_xy_0, tdx_yz_xz_0, tdx_yz_yy_0, \
                                     tdx_yz_yz_0, tdx_yz_zz_0, tdx_yzz_x_0, tdx_yzz_xx_0, tdx_yzz_xy_0, tdx_yzz_xz_0, \
                                     tdx_zz_xx_0, tdx_zz_xy_0, tdx_zz_xz_0, tdy_yy_xx_0, tdy_yy_xy_0, tdy_yy_xz_0, \
                                     tdy_yy_yy_0, tdy_yy_yz_0, tdy_yy_zz_0, tdy_yyy_x_0, tdy_yyy_xx_0, tdy_yyy_xy_0, \
                                     tdy_yyy_xz_0, tdy_yyy_y_0, tdy_yyy_yy_0, tdy_yyy_yz_0, tdy_yyy_z_0, tdy_yyy_zz_0, \
                                     tdy_yyyy_xx_0, tdy_yyyy_xy_0, tdy_yyyy_xz_0, tdy_yyyy_yy_0, tdy_yyyy_yz_0, \
                                     tdy_yyyy_zz_0, tdy_yyyz_xx_0, tdy_yyyz_xy_0, tdy_yyyz_xz_0, tdy_yyyz_yy_0, \
                                     tdy_yyyz_yz_0, tdy_yyyz_zz_0, tdy_yyz_x_0, tdy_yyz_xx_0, tdy_yyz_xy_0, tdy_yyz_xz_0, \
                                     tdy_yyz_y_0, tdy_yyz_yy_0, tdy_yyz_yz_0, tdy_yyz_z_0, tdy_yyz_zz_0, tdy_yyzz_xx_0, \
                                     tdy_yyzz_xy_0, tdy_yyzz_xz_0, tdy_yz_xx_0, tdy_yz_xy_0, tdy_yz_xz_0, tdy_yz_yy_0, \
                                     tdy_yz_yz_0, tdy_yz_zz_0, tdy_yzz_x_0, tdy_yzz_xx_0, tdy_yzz_xy_0, tdy_yzz_xz_0, \
                                     tdy_zz_xx_0, tdy_zz_xy_0, tdy_zz_xz_0, tdz_yy_xx_0, tdz_yy_xy_0, tdz_yy_xz_0, \
                                     tdz_yy_yy_0, tdz_yy_yz_0, tdz_yy_zz_0, tdz_yyy_x_0, tdz_yyy_xx_0, tdz_yyy_xy_0, \
                                     tdz_yyy_xz_0, tdz_yyy_y_0, tdz_yyy_yy_0, tdz_yyy_yz_0, tdz_yyy_z_0, tdz_yyy_zz_0, \
                                     tdz_yyyy_xx_0, tdz_yyyy_xy_0, tdz_yyyy_xz_0, tdz_yyyy_yy_0, tdz_yyyy_yz_0, \
                                     tdz_yyyy_zz_0, tdz_yyyz_xx_0, tdz_yyyz_xy_0, tdz_yyyz_xz_0, tdz_yyyz_yy_0, \
                                     tdz_yyyz_yz_0, tdz_yyyz_zz_0, tdz_yyz_x_0, tdz_yyz_xx_0, tdz_yyz_xy_0, tdz_yyz_xz_0, \
                                     tdz_yyz_y_0, tdz_yyz_yy_0, tdz_yyz_yz_0, tdz_yyz_z_0, tdz_yyz_zz_0, tdz_yyzz_xx_0, \
                                     tdz_yyzz_xy_0, tdz_yyzz_xz_0, tdz_yz_xx_0, tdz_yz_xy_0, tdz_yz_xz_0, tdz_yz_yy_0, \
                                     tdz_yz_yz_0, tdz_yz_zz_0, tdz_yzz_x_0, tdz_yzz_xx_0, tdz_yzz_xy_0, tdz_yzz_xz_0, \
                                     tdz_zz_xx_0, tdz_zz_xy_0, tdz_zz_xz_0, ts_yyy_xx_0, ts_yyy_xy_0, ts_yyy_xz_0, \
                                     ts_yyy_yy_0, ts_yyy_yz_0, ts_yyy_zz_0, ts_yyz_xx_0, ts_yyz_xy_0, ts_yyz_xz_0, \
                                     ts_yyz_yy_0, ts_yyz_yz_0, ts_yyz_zz_0, ts_yzz_xx_0, ts_yzz_xy_0, ts_yzz_xz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yyyy_xx_0[j] = pa_y[j] * tdx_yyy_xx_0[j] + 1.5 * fl1_fx * tdx_yy_xx_0[j];

                tdy_yyyy_xx_0[j] = pa_y[j] * tdy_yyy_xx_0[j] + 1.5 * fl1_fx * tdy_yy_xx_0[j] + 0.5 * fl1_fx * ts_yyy_xx_0[j];

                tdz_yyyy_xx_0[j] = pa_y[j] * tdz_yyy_xx_0[j] + 1.5 * fl1_fx * tdz_yy_xx_0[j];

                tdx_yyyy_xy_0[j] = pa_y[j] * tdx_yyy_xy_0[j] + 1.5 * fl1_fx * tdx_yy_xy_0[j] + 0.5 * fl1_fx * tdx_yyy_x_0[j];

                tdy_yyyy_xy_0[j] = pa_y[j] * tdy_yyy_xy_0[j] + 1.5 * fl1_fx * tdy_yy_xy_0[j] + 0.5 * fl1_fx * tdy_yyy_x_0[j] + 0.5 * fl1_fx * ts_yyy_xy_0[j];

                tdz_yyyy_xy_0[j] = pa_y[j] * tdz_yyy_xy_0[j] + 1.5 * fl1_fx * tdz_yy_xy_0[j] + 0.5 * fl1_fx * tdz_yyy_x_0[j];

                tdx_yyyy_xz_0[j] = pa_y[j] * tdx_yyy_xz_0[j] + 1.5 * fl1_fx * tdx_yy_xz_0[j];

                tdy_yyyy_xz_0[j] = pa_y[j] * tdy_yyy_xz_0[j] + 1.5 * fl1_fx * tdy_yy_xz_0[j] + 0.5 * fl1_fx * ts_yyy_xz_0[j];

                tdz_yyyy_xz_0[j] = pa_y[j] * tdz_yyy_xz_0[j] + 1.5 * fl1_fx * tdz_yy_xz_0[j];

                tdx_yyyy_yy_0[j] = pa_y[j] * tdx_yyy_yy_0[j] + 1.5 * fl1_fx * tdx_yy_yy_0[j] + fl1_fx * tdx_yyy_y_0[j];

                tdy_yyyy_yy_0[j] = pa_y[j] * tdy_yyy_yy_0[j] + 1.5 * fl1_fx * tdy_yy_yy_0[j] + fl1_fx * tdy_yyy_y_0[j] + 0.5 * fl1_fx * ts_yyy_yy_0[j];

                tdz_yyyy_yy_0[j] = pa_y[j] * tdz_yyy_yy_0[j] + 1.5 * fl1_fx * tdz_yy_yy_0[j] + fl1_fx * tdz_yyy_y_0[j];

                tdx_yyyy_yz_0[j] = pa_y[j] * tdx_yyy_yz_0[j] + 1.5 * fl1_fx * tdx_yy_yz_0[j] + 0.5 * fl1_fx * tdx_yyy_z_0[j];

                tdy_yyyy_yz_0[j] = pa_y[j] * tdy_yyy_yz_0[j] + 1.5 * fl1_fx * tdy_yy_yz_0[j] + 0.5 * fl1_fx * tdy_yyy_z_0[j] + 0.5 * fl1_fx * ts_yyy_yz_0[j];

                tdz_yyyy_yz_0[j] = pa_y[j] * tdz_yyy_yz_0[j] + 1.5 * fl1_fx * tdz_yy_yz_0[j] + 0.5 * fl1_fx * tdz_yyy_z_0[j];

                tdx_yyyy_zz_0[j] = pa_y[j] * tdx_yyy_zz_0[j] + 1.5 * fl1_fx * tdx_yy_zz_0[j];

                tdy_yyyy_zz_0[j] = pa_y[j] * tdy_yyy_zz_0[j] + 1.5 * fl1_fx * tdy_yy_zz_0[j] + 0.5 * fl1_fx * ts_yyy_zz_0[j];

                tdz_yyyy_zz_0[j] = pa_y[j] * tdz_yyy_zz_0[j] + 1.5 * fl1_fx * tdz_yy_zz_0[j];

                tdx_yyyz_xx_0[j] = pa_y[j] * tdx_yyz_xx_0[j] + fl1_fx * tdx_yz_xx_0[j];

                tdy_yyyz_xx_0[j] = pa_y[j] * tdy_yyz_xx_0[j] + fl1_fx * tdy_yz_xx_0[j] + 0.5 * fl1_fx * ts_yyz_xx_0[j];

                tdz_yyyz_xx_0[j] = pa_y[j] * tdz_yyz_xx_0[j] + fl1_fx * tdz_yz_xx_0[j];

                tdx_yyyz_xy_0[j] = pa_y[j] * tdx_yyz_xy_0[j] + fl1_fx * tdx_yz_xy_0[j] + 0.5 * fl1_fx * tdx_yyz_x_0[j];

                tdy_yyyz_xy_0[j] = pa_y[j] * tdy_yyz_xy_0[j] + fl1_fx * tdy_yz_xy_0[j] + 0.5 * fl1_fx * tdy_yyz_x_0[j] + 0.5 * fl1_fx * ts_yyz_xy_0[j];

                tdz_yyyz_xy_0[j] = pa_y[j] * tdz_yyz_xy_0[j] + fl1_fx * tdz_yz_xy_0[j] + 0.5 * fl1_fx * tdz_yyz_x_0[j];

                tdx_yyyz_xz_0[j] = pa_y[j] * tdx_yyz_xz_0[j] + fl1_fx * tdx_yz_xz_0[j];

                tdy_yyyz_xz_0[j] = pa_y[j] * tdy_yyz_xz_0[j] + fl1_fx * tdy_yz_xz_0[j] + 0.5 * fl1_fx * ts_yyz_xz_0[j];

                tdz_yyyz_xz_0[j] = pa_y[j] * tdz_yyz_xz_0[j] + fl1_fx * tdz_yz_xz_0[j];

                tdx_yyyz_yy_0[j] = pa_y[j] * tdx_yyz_yy_0[j] + fl1_fx * tdx_yz_yy_0[j] + fl1_fx * tdx_yyz_y_0[j];

                tdy_yyyz_yy_0[j] = pa_y[j] * tdy_yyz_yy_0[j] + fl1_fx * tdy_yz_yy_0[j] + fl1_fx * tdy_yyz_y_0[j] + 0.5 * fl1_fx * ts_yyz_yy_0[j];

                tdz_yyyz_yy_0[j] = pa_y[j] * tdz_yyz_yy_0[j] + fl1_fx * tdz_yz_yy_0[j] + fl1_fx * tdz_yyz_y_0[j];

                tdx_yyyz_yz_0[j] = pa_y[j] * tdx_yyz_yz_0[j] + fl1_fx * tdx_yz_yz_0[j] + 0.5 * fl1_fx * tdx_yyz_z_0[j];

                tdy_yyyz_yz_0[j] = pa_y[j] * tdy_yyz_yz_0[j] + fl1_fx * tdy_yz_yz_0[j] + 0.5 * fl1_fx * tdy_yyz_z_0[j] + 0.5 * fl1_fx * ts_yyz_yz_0[j];

                tdz_yyyz_yz_0[j] = pa_y[j] * tdz_yyz_yz_0[j] + fl1_fx * tdz_yz_yz_0[j] + 0.5 * fl1_fx * tdz_yyz_z_0[j];

                tdx_yyyz_zz_0[j] = pa_y[j] * tdx_yyz_zz_0[j] + fl1_fx * tdx_yz_zz_0[j];

                tdy_yyyz_zz_0[j] = pa_y[j] * tdy_yyz_zz_0[j] + fl1_fx * tdy_yz_zz_0[j] + 0.5 * fl1_fx * ts_yyz_zz_0[j];

                tdz_yyyz_zz_0[j] = pa_y[j] * tdz_yyz_zz_0[j] + fl1_fx * tdz_yz_zz_0[j];

                tdx_yyzz_xx_0[j] = pa_y[j] * tdx_yzz_xx_0[j] + 0.5 * fl1_fx * tdx_zz_xx_0[j];

                tdy_yyzz_xx_0[j] = pa_y[j] * tdy_yzz_xx_0[j] + 0.5 * fl1_fx * tdy_zz_xx_0[j] + 0.5 * fl1_fx * ts_yzz_xx_0[j];

                tdz_yyzz_xx_0[j] = pa_y[j] * tdz_yzz_xx_0[j] + 0.5 * fl1_fx * tdz_zz_xx_0[j];

                tdx_yyzz_xy_0[j] = pa_y[j] * tdx_yzz_xy_0[j] + 0.5 * fl1_fx * tdx_zz_xy_0[j] + 0.5 * fl1_fx * tdx_yzz_x_0[j];

                tdy_yyzz_xy_0[j] = pa_y[j] * tdy_yzz_xy_0[j] + 0.5 * fl1_fx * tdy_zz_xy_0[j] + 0.5 * fl1_fx * tdy_yzz_x_0[j] + 0.5 * fl1_fx * ts_yzz_xy_0[j];

                tdz_yyzz_xy_0[j] = pa_y[j] * tdz_yzz_xy_0[j] + 0.5 * fl1_fx * tdz_zz_xy_0[j] + 0.5 * fl1_fx * tdz_yzz_x_0[j];

                tdx_yyzz_xz_0[j] = pa_y[j] * tdx_yzz_xz_0[j] + 0.5 * fl1_fx * tdx_zz_xz_0[j];

                tdy_yyzz_xz_0[j] = pa_y[j] * tdy_yzz_xz_0[j] + 0.5 * fl1_fx * tdy_zz_xz_0[j] + 0.5 * fl1_fx * ts_yzz_xz_0[j];

                tdz_yyzz_xz_0[j] = pa_y[j] * tdz_yzz_xz_0[j] + 0.5 * fl1_fx * tdz_zz_xz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGD_225_270(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const int32_t              nOSFactors,
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

        auto pidx_d_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tdx_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 51); 

            auto tdy_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 52); 

            auto tdy_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 53); 

            auto tdy_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 53); 

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

            auto tdx_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 25); 

            auto tdy_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 25); 

            auto tdz_yzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 25); 

            auto tdx_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 26); 

            auto tdy_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 26); 

            auto tdz_yzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 26); 

            auto tdx_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 27); 

            auto tdy_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 27); 

            auto tdz_zzz_x_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 27); 

            auto tdx_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 28); 

            auto tdy_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 28); 

            auto tdz_zzz_y_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 28); 

            auto tdx_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * idx + 29); 

            auto tdy_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 30 * bdim + 30 * idx + 29); 

            auto tdz_zzz_z_0 = primBuffer.data(pidx_d_3_1_m0 + 60 * bdim + 30 * idx + 29); 

            auto ts_yzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 51); 

            auto ts_yzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 52); 

            auto ts_yzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 53); 

            auto ts_zzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 54); 

            auto ts_zzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 55); 

            auto ts_zzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 56); 

            auto ts_zzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 57); 

            auto ts_zzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 58); 

            auto ts_zzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 59); 

            // set up pointers to integrals

            auto tdx_yyzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 75); 

            auto tdy_yyzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 75); 

            auto tdz_yyzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 75); 

            auto tdx_yyzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 76); 

            auto tdy_yyzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 76); 

            auto tdz_yyzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 76); 

            auto tdx_yyzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 77); 

            auto tdy_yyzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 77); 

            auto tdz_yyzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 77); 

            auto tdx_yzzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 78); 

            auto tdy_yzzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 78); 

            auto tdz_yzzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 78); 

            auto tdx_yzzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 79); 

            auto tdy_yzzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 79); 

            auto tdz_yzzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 79); 

            auto tdx_yzzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 80); 

            auto tdy_yzzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 80); 

            auto tdz_yzzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 80); 

            auto tdx_yzzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 81); 

            auto tdy_yzzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 81); 

            auto tdz_yzzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 81); 

            auto tdx_yzzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 82); 

            auto tdy_yzzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 82); 

            auto tdz_yzzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 82); 

            auto tdx_yzzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 83); 

            auto tdy_yzzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 83); 

            auto tdz_yzzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 83); 

            auto tdx_zzzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 84); 

            auto tdy_zzzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 84); 

            auto tdz_zzzz_xx_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 84); 

            auto tdx_zzzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 85); 

            auto tdy_zzzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 85); 

            auto tdz_zzzz_xy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 85); 

            auto tdx_zzzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 86); 

            auto tdy_zzzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 86); 

            auto tdz_zzzz_xz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 86); 

            auto tdx_zzzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 87); 

            auto tdy_zzzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 87); 

            auto tdz_zzzz_yy_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 87); 

            auto tdx_zzzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 88); 

            auto tdy_zzzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 88); 

            auto tdz_zzzz_yz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 88); 

            auto tdx_zzzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * idx + 89); 

            auto tdy_zzzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 90 * bdim + 90 * idx + 89); 

            auto tdz_zzzz_zz_0 = primBuffer.data(pidx_d_4_2_m0 + 180 * bdim + 90 * idx + 89); 

            // Batch of Integrals (225,270)

            #pragma omp simd aligned(fx, pa_y, pa_z, tdx_yyzz_yy_0, tdx_yyzz_yz_0, tdx_yyzz_zz_0, tdx_yzz_y_0, \
                                     tdx_yzz_yy_0, tdx_yzz_yz_0, tdx_yzz_z_0, tdx_yzz_zz_0, tdx_yzzz_xx_0, \
                                     tdx_yzzz_xy_0, tdx_yzzz_xz_0, tdx_yzzz_yy_0, tdx_yzzz_yz_0, tdx_yzzz_zz_0, \
                                     tdx_zz_xx_0, tdx_zz_xy_0, tdx_zz_xz_0, tdx_zz_yy_0, tdx_zz_yz_0, tdx_zz_zz_0, \
                                     tdx_zzz_x_0, tdx_zzz_xx_0, tdx_zzz_xy_0, tdx_zzz_xz_0, tdx_zzz_y_0, tdx_zzz_yy_0, \
                                     tdx_zzz_yz_0, tdx_zzz_z_0, tdx_zzz_zz_0, tdx_zzzz_xx_0, tdx_zzzz_xy_0, \
                                     tdx_zzzz_xz_0, tdx_zzzz_yy_0, tdx_zzzz_yz_0, tdx_zzzz_zz_0, tdy_yyzz_yy_0, \
                                     tdy_yyzz_yz_0, tdy_yyzz_zz_0, tdy_yzz_y_0, tdy_yzz_yy_0, tdy_yzz_yz_0, tdy_yzz_z_0, \
                                     tdy_yzz_zz_0, tdy_yzzz_xx_0, tdy_yzzz_xy_0, tdy_yzzz_xz_0, tdy_yzzz_yy_0, \
                                     tdy_yzzz_yz_0, tdy_yzzz_zz_0, tdy_zz_xx_0, tdy_zz_xy_0, tdy_zz_xz_0, tdy_zz_yy_0, \
                                     tdy_zz_yz_0, tdy_zz_zz_0, tdy_zzz_x_0, tdy_zzz_xx_0, tdy_zzz_xy_0, tdy_zzz_xz_0, \
                                     tdy_zzz_y_0, tdy_zzz_yy_0, tdy_zzz_yz_0, tdy_zzz_z_0, tdy_zzz_zz_0, tdy_zzzz_xx_0, \
                                     tdy_zzzz_xy_0, tdy_zzzz_xz_0, tdy_zzzz_yy_0, tdy_zzzz_yz_0, tdy_zzzz_zz_0, \
                                     tdz_yyzz_yy_0, tdz_yyzz_yz_0, tdz_yyzz_zz_0, tdz_yzz_y_0, tdz_yzz_yy_0, \
                                     tdz_yzz_yz_0, tdz_yzz_z_0, tdz_yzz_zz_0, tdz_yzzz_xx_0, tdz_yzzz_xy_0, \
                                     tdz_yzzz_xz_0, tdz_yzzz_yy_0, tdz_yzzz_yz_0, tdz_yzzz_zz_0, tdz_zz_xx_0, \
                                     tdz_zz_xy_0, tdz_zz_xz_0, tdz_zz_yy_0, tdz_zz_yz_0, tdz_zz_zz_0, tdz_zzz_x_0, \
                                     tdz_zzz_xx_0, tdz_zzz_xy_0, tdz_zzz_xz_0, tdz_zzz_y_0, tdz_zzz_yy_0, tdz_zzz_yz_0, \
                                     tdz_zzz_z_0, tdz_zzz_zz_0, tdz_zzzz_xx_0, tdz_zzzz_xy_0, tdz_zzzz_xz_0, \
                                     tdz_zzzz_yy_0, tdz_zzzz_yz_0, tdz_zzzz_zz_0, ts_yzz_yy_0, ts_yzz_yz_0, ts_yzz_zz_0, \
                                     ts_zzz_xx_0, ts_zzz_xy_0, ts_zzz_xz_0, ts_zzz_yy_0, ts_zzz_yz_0, ts_zzz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yyzz_yy_0[j] = pa_y[j] * tdx_yzz_yy_0[j] + 0.5 * fl1_fx * tdx_zz_yy_0[j] + fl1_fx * tdx_yzz_y_0[j];

                tdy_yyzz_yy_0[j] = pa_y[j] * tdy_yzz_yy_0[j] + 0.5 * fl1_fx * tdy_zz_yy_0[j] + fl1_fx * tdy_yzz_y_0[j] + 0.5 * fl1_fx * ts_yzz_yy_0[j];

                tdz_yyzz_yy_0[j] = pa_y[j] * tdz_yzz_yy_0[j] + 0.5 * fl1_fx * tdz_zz_yy_0[j] + fl1_fx * tdz_yzz_y_0[j];

                tdx_yyzz_yz_0[j] = pa_y[j] * tdx_yzz_yz_0[j] + 0.5 * fl1_fx * tdx_zz_yz_0[j] + 0.5 * fl1_fx * tdx_yzz_z_0[j];

                tdy_yyzz_yz_0[j] = pa_y[j] * tdy_yzz_yz_0[j] + 0.5 * fl1_fx * tdy_zz_yz_0[j] + 0.5 * fl1_fx * tdy_yzz_z_0[j] + 0.5 * fl1_fx * ts_yzz_yz_0[j];

                tdz_yyzz_yz_0[j] = pa_y[j] * tdz_yzz_yz_0[j] + 0.5 * fl1_fx * tdz_zz_yz_0[j] + 0.5 * fl1_fx * tdz_yzz_z_0[j];

                tdx_yyzz_zz_0[j] = pa_y[j] * tdx_yzz_zz_0[j] + 0.5 * fl1_fx * tdx_zz_zz_0[j];

                tdy_yyzz_zz_0[j] = pa_y[j] * tdy_yzz_zz_0[j] + 0.5 * fl1_fx * tdy_zz_zz_0[j] + 0.5 * fl1_fx * ts_yzz_zz_0[j];

                tdz_yyzz_zz_0[j] = pa_y[j] * tdz_yzz_zz_0[j] + 0.5 * fl1_fx * tdz_zz_zz_0[j];

                tdx_yzzz_xx_0[j] = pa_y[j] * tdx_zzz_xx_0[j];

                tdy_yzzz_xx_0[j] = pa_y[j] * tdy_zzz_xx_0[j] + 0.5 * fl1_fx * ts_zzz_xx_0[j];

                tdz_yzzz_xx_0[j] = pa_y[j] * tdz_zzz_xx_0[j];

                tdx_yzzz_xy_0[j] = pa_y[j] * tdx_zzz_xy_0[j] + 0.5 * fl1_fx * tdx_zzz_x_0[j];

                tdy_yzzz_xy_0[j] = pa_y[j] * tdy_zzz_xy_0[j] + 0.5 * fl1_fx * tdy_zzz_x_0[j] + 0.5 * fl1_fx * ts_zzz_xy_0[j];

                tdz_yzzz_xy_0[j] = pa_y[j] * tdz_zzz_xy_0[j] + 0.5 * fl1_fx * tdz_zzz_x_0[j];

                tdx_yzzz_xz_0[j] = pa_y[j] * tdx_zzz_xz_0[j];

                tdy_yzzz_xz_0[j] = pa_y[j] * tdy_zzz_xz_0[j] + 0.5 * fl1_fx * ts_zzz_xz_0[j];

                tdz_yzzz_xz_0[j] = pa_y[j] * tdz_zzz_xz_0[j];

                tdx_yzzz_yy_0[j] = pa_y[j] * tdx_zzz_yy_0[j] + fl1_fx * tdx_zzz_y_0[j];

                tdy_yzzz_yy_0[j] = pa_y[j] * tdy_zzz_yy_0[j] + fl1_fx * tdy_zzz_y_0[j] + 0.5 * fl1_fx * ts_zzz_yy_0[j];

                tdz_yzzz_yy_0[j] = pa_y[j] * tdz_zzz_yy_0[j] + fl1_fx * tdz_zzz_y_0[j];

                tdx_yzzz_yz_0[j] = pa_y[j] * tdx_zzz_yz_0[j] + 0.5 * fl1_fx * tdx_zzz_z_0[j];

                tdy_yzzz_yz_0[j] = pa_y[j] * tdy_zzz_yz_0[j] + 0.5 * fl1_fx * tdy_zzz_z_0[j] + 0.5 * fl1_fx * ts_zzz_yz_0[j];

                tdz_yzzz_yz_0[j] = pa_y[j] * tdz_zzz_yz_0[j] + 0.5 * fl1_fx * tdz_zzz_z_0[j];

                tdx_yzzz_zz_0[j] = pa_y[j] * tdx_zzz_zz_0[j];

                tdy_yzzz_zz_0[j] = pa_y[j] * tdy_zzz_zz_0[j] + 0.5 * fl1_fx * ts_zzz_zz_0[j];

                tdz_yzzz_zz_0[j] = pa_y[j] * tdz_zzz_zz_0[j];

                tdx_zzzz_xx_0[j] = pa_z[j] * tdx_zzz_xx_0[j] + 1.5 * fl1_fx * tdx_zz_xx_0[j];

                tdy_zzzz_xx_0[j] = pa_z[j] * tdy_zzz_xx_0[j] + 1.5 * fl1_fx * tdy_zz_xx_0[j];

                tdz_zzzz_xx_0[j] = pa_z[j] * tdz_zzz_xx_0[j] + 1.5 * fl1_fx * tdz_zz_xx_0[j] + 0.5 * fl1_fx * ts_zzz_xx_0[j];

                tdx_zzzz_xy_0[j] = pa_z[j] * tdx_zzz_xy_0[j] + 1.5 * fl1_fx * tdx_zz_xy_0[j];

                tdy_zzzz_xy_0[j] = pa_z[j] * tdy_zzz_xy_0[j] + 1.5 * fl1_fx * tdy_zz_xy_0[j];

                tdz_zzzz_xy_0[j] = pa_z[j] * tdz_zzz_xy_0[j] + 1.5 * fl1_fx * tdz_zz_xy_0[j] + 0.5 * fl1_fx * ts_zzz_xy_0[j];

                tdx_zzzz_xz_0[j] = pa_z[j] * tdx_zzz_xz_0[j] + 1.5 * fl1_fx * tdx_zz_xz_0[j] + 0.5 * fl1_fx * tdx_zzz_x_0[j];

                tdy_zzzz_xz_0[j] = pa_z[j] * tdy_zzz_xz_0[j] + 1.5 * fl1_fx * tdy_zz_xz_0[j] + 0.5 * fl1_fx * tdy_zzz_x_0[j];

                tdz_zzzz_xz_0[j] = pa_z[j] * tdz_zzz_xz_0[j] + 1.5 * fl1_fx * tdz_zz_xz_0[j] + 0.5 * fl1_fx * tdz_zzz_x_0[j] + 0.5 * fl1_fx * ts_zzz_xz_0[j];

                tdx_zzzz_yy_0[j] = pa_z[j] * tdx_zzz_yy_0[j] + 1.5 * fl1_fx * tdx_zz_yy_0[j];

                tdy_zzzz_yy_0[j] = pa_z[j] * tdy_zzz_yy_0[j] + 1.5 * fl1_fx * tdy_zz_yy_0[j];

                tdz_zzzz_yy_0[j] = pa_z[j] * tdz_zzz_yy_0[j] + 1.5 * fl1_fx * tdz_zz_yy_0[j] + 0.5 * fl1_fx * ts_zzz_yy_0[j];

                tdx_zzzz_yz_0[j] = pa_z[j] * tdx_zzz_yz_0[j] + 1.5 * fl1_fx * tdx_zz_yz_0[j] + 0.5 * fl1_fx * tdx_zzz_y_0[j];

                tdy_zzzz_yz_0[j] = pa_z[j] * tdy_zzz_yz_0[j] + 1.5 * fl1_fx * tdy_zz_yz_0[j] + 0.5 * fl1_fx * tdy_zzz_y_0[j];

                tdz_zzzz_yz_0[j] = pa_z[j] * tdz_zzz_yz_0[j] + 1.5 * fl1_fx * tdz_zz_yz_0[j] + 0.5 * fl1_fx * tdz_zzz_y_0[j] + 0.5 * fl1_fx * ts_zzz_yz_0[j];

                tdx_zzzz_zz_0[j] = pa_z[j] * tdx_zzz_zz_0[j] + 1.5 * fl1_fx * tdx_zz_zz_0[j] + fl1_fx * tdx_zzz_z_0[j];

                tdy_zzzz_zz_0[j] = pa_z[j] * tdy_zzz_zz_0[j] + 1.5 * fl1_fx * tdy_zz_zz_0[j] + fl1_fx * tdy_zzz_z_0[j];

                tdz_zzzz_zz_0[j] = pa_z[j] * tdz_zzz_zz_0[j] + 1.5 * fl1_fx * tdz_zz_zz_0[j] + fl1_fx * tdz_zzz_z_0[j] + 0.5 * fl1_fx * ts_zzz_zz_0[j];
            }

            idx++;
        }
    }


} // ediprecfunc namespace

