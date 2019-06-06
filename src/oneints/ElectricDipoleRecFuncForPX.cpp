//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricDipoleRecFuncForPX.hpp"

namespace ediprecfunc { // ediprecfunc namespace

    void
    compElectricDipoleForPP(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_1_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tdx_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx); 

            auto tdy_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx); 

            auto tdz_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx); 

            auto tdx_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 1); 

            auto tdy_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 2); 

            auto tdy_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 2); 

            auto tdx_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + idx); 

            auto tdy_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + bdim + idx); 

            auto tdz_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + 2 * bdim + idx); 

            auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx); 

            auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1); 

            auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2); 

            // set up pointers to integrals

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

            auto tdx_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 6); 

            auto tdy_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 6); 

            auto tdz_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 6); 

            auto tdx_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 7); 

            auto tdy_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 7); 

            auto tdz_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 7); 

            auto tdx_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 8); 

            auto tdy_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 8); 

            auto tdz_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 8); 

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, tdx_0_0_0, tdx_0_x_0, tdx_0_y_0, tdx_0_z_0, tdx_x_x_0, \
                                     tdx_x_y_0, tdx_x_z_0, tdx_y_x_0, tdx_y_y_0, tdx_y_z_0, tdx_z_x_0, tdx_z_y_0, \
                                     tdx_z_z_0, tdy_0_0_0, tdy_0_x_0, tdy_0_y_0, tdy_0_z_0, tdy_x_x_0, tdy_x_y_0, \
                                     tdy_x_z_0, tdy_y_x_0, tdy_y_y_0, tdy_y_z_0, tdy_z_x_0, tdy_z_y_0, tdy_z_z_0, \
                                     tdz_0_0_0, tdz_0_x_0, tdz_0_y_0, tdz_0_z_0, tdz_x_x_0, tdz_x_y_0, tdz_x_z_0, \
                                     tdz_y_x_0, tdz_y_y_0, tdz_y_z_0, tdz_z_x_0, tdz_z_y_0, tdz_z_z_0, ts_0_x_0, \
                                     ts_0_y_0, ts_0_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_x_x_0[j] = pa_x[j] * tdx_0_x_0[j] + 0.5 * fl1_fx * tdx_0_0_0[j] + 0.5 * fl1_fx * ts_0_x_0[j];

                tdy_x_x_0[j] = pa_x[j] * tdy_0_x_0[j] + 0.5 * fl1_fx * tdy_0_0_0[j];

                tdz_x_x_0[j] = pa_x[j] * tdz_0_x_0[j] + 0.5 * fl1_fx * tdz_0_0_0[j];

                tdx_x_y_0[j] = pa_x[j] * tdx_0_y_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

                tdy_x_y_0[j] = pa_x[j] * tdy_0_y_0[j];

                tdz_x_y_0[j] = pa_x[j] * tdz_0_y_0[j];

                tdx_x_z_0[j] = pa_x[j] * tdx_0_z_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

                tdy_x_z_0[j] = pa_x[j] * tdy_0_z_0[j];

                tdz_x_z_0[j] = pa_x[j] * tdz_0_z_0[j];

                tdx_y_x_0[j] = pa_y[j] * tdx_0_x_0[j];

                tdy_y_x_0[j] = pa_y[j] * tdy_0_x_0[j] + 0.5 * fl1_fx * ts_0_x_0[j];

                tdz_y_x_0[j] = pa_y[j] * tdz_0_x_0[j];

                tdx_y_y_0[j] = pa_y[j] * tdx_0_y_0[j] + 0.5 * fl1_fx * tdx_0_0_0[j];

                tdy_y_y_0[j] = pa_y[j] * tdy_0_y_0[j] + 0.5 * fl1_fx * tdy_0_0_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

                tdz_y_y_0[j] = pa_y[j] * tdz_0_y_0[j] + 0.5 * fl1_fx * tdz_0_0_0[j];

                tdx_y_z_0[j] = pa_y[j] * tdx_0_z_0[j];

                tdy_y_z_0[j] = pa_y[j] * tdy_0_z_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

                tdz_y_z_0[j] = pa_y[j] * tdz_0_z_0[j];

                tdx_z_x_0[j] = pa_z[j] * tdx_0_x_0[j];

                tdy_z_x_0[j] = pa_z[j] * tdy_0_x_0[j];

                tdz_z_x_0[j] = pa_z[j] * tdz_0_x_0[j] + 0.5 * fl1_fx * ts_0_x_0[j];

                tdx_z_y_0[j] = pa_z[j] * tdx_0_y_0[j];

                tdy_z_y_0[j] = pa_z[j] * tdy_0_y_0[j];

                tdz_z_y_0[j] = pa_z[j] * tdz_0_y_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

                tdx_z_z_0[j] = pa_z[j] * tdx_0_z_0[j] + 0.5 * fl1_fx * tdx_0_0_0[j];

                tdy_z_z_0[j] = pa_z[j] * tdy_0_z_0[j] + 0.5 * fl1_fx * tdy_0_0_0[j];

                tdz_z_z_0[j] = pa_z[j] * tdz_0_z_0[j] + 0.5 * fl1_fx * tdz_0_0_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForPD(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForPD_0_27(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForPD_27_54(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 
    }

    void
    compElectricDipoleForPD_0_27(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (0,27)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_1_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

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

            auto tdx_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx); 

            auto tdy_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx); 

            auto tdz_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx); 

            auto tdx_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 1); 

            auto tdy_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 2); 

            auto tdy_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 2); 

            auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx); 

            auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1); 

            auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2); 

            auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3); 

            auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4); 

            auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5); 

            // set up pointers to integrals

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

            // Batch of Integrals (0,27)

            #pragma omp simd aligned(fx, pa_x, pa_y, tdx_0_x_0, tdx_0_xx_0, tdx_0_xy_0, tdx_0_xz_0, tdx_0_y_0, \
                                     tdx_0_yy_0, tdx_0_yz_0, tdx_0_z_0, tdx_0_zz_0, tdx_x_xx_0, tdx_x_xy_0, tdx_x_xz_0, \
                                     tdx_x_yy_0, tdx_x_yz_0, tdx_x_zz_0, tdx_y_xx_0, tdx_y_xy_0, tdx_y_xz_0, tdy_0_x_0, \
                                     tdy_0_xx_0, tdy_0_xy_0, tdy_0_xz_0, tdy_0_y_0, tdy_0_yy_0, tdy_0_yz_0, tdy_0_z_0, \
                                     tdy_0_zz_0, tdy_x_xx_0, tdy_x_xy_0, tdy_x_xz_0, tdy_x_yy_0, tdy_x_yz_0, tdy_x_zz_0, \
                                     tdy_y_xx_0, tdy_y_xy_0, tdy_y_xz_0, tdz_0_x_0, tdz_0_xx_0, tdz_0_xy_0, tdz_0_xz_0, \
                                     tdz_0_y_0, tdz_0_yy_0, tdz_0_yz_0, tdz_0_z_0, tdz_0_zz_0, tdz_x_xx_0, tdz_x_xy_0, \
                                     tdz_x_xz_0, tdz_x_yy_0, tdz_x_yz_0, tdz_x_zz_0, tdz_y_xx_0, tdz_y_xy_0, tdz_y_xz_0, \
                                     ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, ts_0_yy_0, ts_0_yz_0, ts_0_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_x_xx_0[j] = pa_x[j] * tdx_0_xx_0[j] + fl1_fx * tdx_0_x_0[j] + 0.5 * fl1_fx * ts_0_xx_0[j];

                tdy_x_xx_0[j] = pa_x[j] * tdy_0_xx_0[j] + fl1_fx * tdy_0_x_0[j];

                tdz_x_xx_0[j] = pa_x[j] * tdz_0_xx_0[j] + fl1_fx * tdz_0_x_0[j];

                tdx_x_xy_0[j] = pa_x[j] * tdx_0_xy_0[j] + 0.5 * fl1_fx * tdx_0_y_0[j] + 0.5 * fl1_fx * ts_0_xy_0[j];

                tdy_x_xy_0[j] = pa_x[j] * tdy_0_xy_0[j] + 0.5 * fl1_fx * tdy_0_y_0[j];

                tdz_x_xy_0[j] = pa_x[j] * tdz_0_xy_0[j] + 0.5 * fl1_fx * tdz_0_y_0[j];

                tdx_x_xz_0[j] = pa_x[j] * tdx_0_xz_0[j] + 0.5 * fl1_fx * tdx_0_z_0[j] + 0.5 * fl1_fx * ts_0_xz_0[j];

                tdy_x_xz_0[j] = pa_x[j] * tdy_0_xz_0[j] + 0.5 * fl1_fx * tdy_0_z_0[j];

                tdz_x_xz_0[j] = pa_x[j] * tdz_0_xz_0[j] + 0.5 * fl1_fx * tdz_0_z_0[j];

                tdx_x_yy_0[j] = pa_x[j] * tdx_0_yy_0[j] + 0.5 * fl1_fx * ts_0_yy_0[j];

                tdy_x_yy_0[j] = pa_x[j] * tdy_0_yy_0[j];

                tdz_x_yy_0[j] = pa_x[j] * tdz_0_yy_0[j];

                tdx_x_yz_0[j] = pa_x[j] * tdx_0_yz_0[j] + 0.5 * fl1_fx * ts_0_yz_0[j];

                tdy_x_yz_0[j] = pa_x[j] * tdy_0_yz_0[j];

                tdz_x_yz_0[j] = pa_x[j] * tdz_0_yz_0[j];

                tdx_x_zz_0[j] = pa_x[j] * tdx_0_zz_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];

                tdy_x_zz_0[j] = pa_x[j] * tdy_0_zz_0[j];

                tdz_x_zz_0[j] = pa_x[j] * tdz_0_zz_0[j];

                tdx_y_xx_0[j] = pa_y[j] * tdx_0_xx_0[j];

                tdy_y_xx_0[j] = pa_y[j] * tdy_0_xx_0[j] + 0.5 * fl1_fx * ts_0_xx_0[j];

                tdz_y_xx_0[j] = pa_y[j] * tdz_0_xx_0[j];

                tdx_y_xy_0[j] = pa_y[j] * tdx_0_xy_0[j] + 0.5 * fl1_fx * tdx_0_x_0[j];

                tdy_y_xy_0[j] = pa_y[j] * tdy_0_xy_0[j] + 0.5 * fl1_fx * tdy_0_x_0[j] + 0.5 * fl1_fx * ts_0_xy_0[j];

                tdz_y_xy_0[j] = pa_y[j] * tdz_0_xy_0[j] + 0.5 * fl1_fx * tdz_0_x_0[j];

                tdx_y_xz_0[j] = pa_y[j] * tdx_0_xz_0[j];

                tdy_y_xz_0[j] = pa_y[j] * tdy_0_xz_0[j] + 0.5 * fl1_fx * ts_0_xz_0[j];

                tdz_y_xz_0[j] = pa_y[j] * tdz_0_xz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForPD_27_54(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (27,54)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_1_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto tdx_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx); 

            auto tdy_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx); 

            auto tdz_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx); 

            auto tdx_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 1); 

            auto tdy_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 2); 

            auto tdy_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 2); 

            auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx); 

            auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1); 

            auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2); 

            auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3); 

            auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4); 

            auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5); 

            // set up pointers to integrals

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

            // Batch of Integrals (27,54)

            #pragma omp simd aligned(fx, pa_y, pa_z, tdx_0_x_0, tdx_0_xx_0, tdx_0_xy_0, tdx_0_xz_0, tdx_0_y_0, \
                                     tdx_0_yy_0, tdx_0_yz_0, tdx_0_z_0, tdx_0_zz_0, tdx_y_yy_0, tdx_y_yz_0, tdx_y_zz_0, \
                                     tdx_z_xx_0, tdx_z_xy_0, tdx_z_xz_0, tdx_z_yy_0, tdx_z_yz_0, tdx_z_zz_0, tdy_0_x_0, \
                                     tdy_0_xx_0, tdy_0_xy_0, tdy_0_xz_0, tdy_0_y_0, tdy_0_yy_0, tdy_0_yz_0, tdy_0_z_0, \
                                     tdy_0_zz_0, tdy_y_yy_0, tdy_y_yz_0, tdy_y_zz_0, tdy_z_xx_0, tdy_z_xy_0, tdy_z_xz_0, \
                                     tdy_z_yy_0, tdy_z_yz_0, tdy_z_zz_0, tdz_0_x_0, tdz_0_xx_0, tdz_0_xy_0, tdz_0_xz_0, \
                                     tdz_0_y_0, tdz_0_yy_0, tdz_0_yz_0, tdz_0_z_0, tdz_0_zz_0, tdz_y_yy_0, tdz_y_yz_0, \
                                     tdz_y_zz_0, tdz_z_xx_0, tdz_z_xy_0, tdz_z_xz_0, tdz_z_yy_0, tdz_z_yz_0, tdz_z_zz_0, \
                                     ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, ts_0_yy_0, ts_0_yz_0, ts_0_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_y_yy_0[j] = pa_y[j] * tdx_0_yy_0[j] + fl1_fx * tdx_0_y_0[j];

                tdy_y_yy_0[j] = pa_y[j] * tdy_0_yy_0[j] + fl1_fx * tdy_0_y_0[j] + 0.5 * fl1_fx * ts_0_yy_0[j];

                tdz_y_yy_0[j] = pa_y[j] * tdz_0_yy_0[j] + fl1_fx * tdz_0_y_0[j];

                tdx_y_yz_0[j] = pa_y[j] * tdx_0_yz_0[j] + 0.5 * fl1_fx * tdx_0_z_0[j];

                tdy_y_yz_0[j] = pa_y[j] * tdy_0_yz_0[j] + 0.5 * fl1_fx * tdy_0_z_0[j] + 0.5 * fl1_fx * ts_0_yz_0[j];

                tdz_y_yz_0[j] = pa_y[j] * tdz_0_yz_0[j] + 0.5 * fl1_fx * tdz_0_z_0[j];

                tdx_y_zz_0[j] = pa_y[j] * tdx_0_zz_0[j];

                tdy_y_zz_0[j] = pa_y[j] * tdy_0_zz_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];

                tdz_y_zz_0[j] = pa_y[j] * tdz_0_zz_0[j];

                tdx_z_xx_0[j] = pa_z[j] * tdx_0_xx_0[j];

                tdy_z_xx_0[j] = pa_z[j] * tdy_0_xx_0[j];

                tdz_z_xx_0[j] = pa_z[j] * tdz_0_xx_0[j] + 0.5 * fl1_fx * ts_0_xx_0[j];

                tdx_z_xy_0[j] = pa_z[j] * tdx_0_xy_0[j];

                tdy_z_xy_0[j] = pa_z[j] * tdy_0_xy_0[j];

                tdz_z_xy_0[j] = pa_z[j] * tdz_0_xy_0[j] + 0.5 * fl1_fx * ts_0_xy_0[j];

                tdx_z_xz_0[j] = pa_z[j] * tdx_0_xz_0[j] + 0.5 * fl1_fx * tdx_0_x_0[j];

                tdy_z_xz_0[j] = pa_z[j] * tdy_0_xz_0[j] + 0.5 * fl1_fx * tdy_0_x_0[j];

                tdz_z_xz_0[j] = pa_z[j] * tdz_0_xz_0[j] + 0.5 * fl1_fx * tdz_0_x_0[j] + 0.5 * fl1_fx * ts_0_xz_0[j];

                tdx_z_yy_0[j] = pa_z[j] * tdx_0_yy_0[j];

                tdy_z_yy_0[j] = pa_z[j] * tdy_0_yy_0[j];

                tdz_z_yy_0[j] = pa_z[j] * tdz_0_yy_0[j] + 0.5 * fl1_fx * ts_0_yy_0[j];

                tdx_z_yz_0[j] = pa_z[j] * tdx_0_yz_0[j] + 0.5 * fl1_fx * tdx_0_y_0[j];

                tdy_z_yz_0[j] = pa_z[j] * tdy_0_yz_0[j] + 0.5 * fl1_fx * tdy_0_y_0[j];

                tdz_z_yz_0[j] = pa_z[j] * tdz_0_yz_0[j] + 0.5 * fl1_fx * tdz_0_y_0[j] + 0.5 * fl1_fx * ts_0_yz_0[j];

                tdx_z_zz_0[j] = pa_z[j] * tdx_0_zz_0[j] + fl1_fx * tdx_0_z_0[j];

                tdy_z_zz_0[j] = pa_z[j] * tdy_0_zz_0[j] + fl1_fx * tdy_0_z_0[j];

                tdz_z_zz_0[j] = pa_z[j] * tdz_0_zz_0[j] + fl1_fx * tdz_0_z_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDP(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForDP_0_27(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForDP_27_54(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 
    }

    void
    compElectricDipoleForDP_0_27(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (0,27)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

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

            auto tdx_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 6); 

            auto tdy_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 6); 

            auto tdz_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 6); 

            auto tdx_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 7); 

            auto tdy_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 7); 

            auto tdz_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 7); 

            auto tdx_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 8); 

            auto tdy_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 8); 

            auto tdz_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 8); 

            auto tdx_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx); 

            auto tdy_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx); 

            auto tdz_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx); 

            auto tdx_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 1); 

            auto tdy_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 2); 

            auto tdy_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 2); 

            auto tdx_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx); 

            auto tdy_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx); 

            auto tdz_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx); 

            auto tdx_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 1); 

            auto tdy_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 2); 

            auto tdy_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 2); 

            auto ts_x_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx); 

            auto ts_x_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 1); 

            auto ts_x_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 2); 

            auto ts_y_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 3); 

            auto ts_y_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 4); 

            auto ts_y_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 5); 

            auto ts_z_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 6); 

            auto ts_z_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 7); 

            auto ts_z_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 8); 

            // set up pointers to integrals

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

            // Batch of Integrals (0,27)

            #pragma omp simd aligned(fx, pa_x, tdx_0_x_0, tdx_0_y_0, tdx_0_z_0, tdx_x_0_0, tdx_x_x_0, tdx_x_y_0, \
                                     tdx_x_z_0, tdx_xx_x_0, tdx_xx_y_0, tdx_xx_z_0, tdx_xy_x_0, tdx_xy_y_0, tdx_xy_z_0, \
                                     tdx_xz_x_0, tdx_xz_y_0, tdx_xz_z_0, tdx_y_0_0, tdx_y_x_0, tdx_y_y_0, tdx_y_z_0, \
                                     tdx_z_0_0, tdx_z_x_0, tdx_z_y_0, tdx_z_z_0, tdy_0_x_0, tdy_0_y_0, tdy_0_z_0, \
                                     tdy_x_0_0, tdy_x_x_0, tdy_x_y_0, tdy_x_z_0, tdy_xx_x_0, tdy_xx_y_0, tdy_xx_z_0, \
                                     tdy_xy_x_0, tdy_xy_y_0, tdy_xy_z_0, tdy_xz_x_0, tdy_xz_y_0, tdy_xz_z_0, tdy_y_0_0, \
                                     tdy_y_x_0, tdy_y_y_0, tdy_y_z_0, tdy_z_0_0, tdy_z_x_0, tdy_z_y_0, tdy_z_z_0, \
                                     tdz_0_x_0, tdz_0_y_0, tdz_0_z_0, tdz_x_0_0, tdz_x_x_0, tdz_x_y_0, tdz_x_z_0, \
                                     tdz_xx_x_0, tdz_xx_y_0, tdz_xx_z_0, tdz_xy_x_0, tdz_xy_y_0, tdz_xy_z_0, tdz_xz_x_0, \
                                     tdz_xz_y_0, tdz_xz_z_0, tdz_y_0_0, tdz_y_x_0, tdz_y_y_0, tdz_y_z_0, tdz_z_0_0, \
                                     tdz_z_x_0, tdz_z_y_0, tdz_z_z_0, ts_x_x_0, ts_x_y_0, ts_x_z_0, ts_y_x_0, ts_y_y_0, \
                                     ts_y_z_0, ts_z_x_0, ts_z_y_0, ts_z_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xx_x_0[j] = pa_x[j] * tdx_x_x_0[j] + 0.5 * fl1_fx * tdx_0_x_0[j] + 0.5 * fl1_fx * tdx_x_0_0[j] + 0.5 * fl1_fx * ts_x_x_0[j];

                tdy_xx_x_0[j] = pa_x[j] * tdy_x_x_0[j] + 0.5 * fl1_fx * tdy_0_x_0[j] + 0.5 * fl1_fx * tdy_x_0_0[j];

                tdz_xx_x_0[j] = pa_x[j] * tdz_x_x_0[j] + 0.5 * fl1_fx * tdz_0_x_0[j] + 0.5 * fl1_fx * tdz_x_0_0[j];

                tdx_xx_y_0[j] = pa_x[j] * tdx_x_y_0[j] + 0.5 * fl1_fx * tdx_0_y_0[j] + 0.5 * fl1_fx * ts_x_y_0[j];

                tdy_xx_y_0[j] = pa_x[j] * tdy_x_y_0[j] + 0.5 * fl1_fx * tdy_0_y_0[j];

                tdz_xx_y_0[j] = pa_x[j] * tdz_x_y_0[j] + 0.5 * fl1_fx * tdz_0_y_0[j];

                tdx_xx_z_0[j] = pa_x[j] * tdx_x_z_0[j] + 0.5 * fl1_fx * tdx_0_z_0[j] + 0.5 * fl1_fx * ts_x_z_0[j];

                tdy_xx_z_0[j] = pa_x[j] * tdy_x_z_0[j] + 0.5 * fl1_fx * tdy_0_z_0[j];

                tdz_xx_z_0[j] = pa_x[j] * tdz_x_z_0[j] + 0.5 * fl1_fx * tdz_0_z_0[j];

                tdx_xy_x_0[j] = pa_x[j] * tdx_y_x_0[j] + 0.5 * fl1_fx * tdx_y_0_0[j] + 0.5 * fl1_fx * ts_y_x_0[j];

                tdy_xy_x_0[j] = pa_x[j] * tdy_y_x_0[j] + 0.5 * fl1_fx * tdy_y_0_0[j];

                tdz_xy_x_0[j] = pa_x[j] * tdz_y_x_0[j] + 0.5 * fl1_fx * tdz_y_0_0[j];

                tdx_xy_y_0[j] = pa_x[j] * tdx_y_y_0[j] + 0.5 * fl1_fx * ts_y_y_0[j];

                tdy_xy_y_0[j] = pa_x[j] * tdy_y_y_0[j];

                tdz_xy_y_0[j] = pa_x[j] * tdz_y_y_0[j];

                tdx_xy_z_0[j] = pa_x[j] * tdx_y_z_0[j] + 0.5 * fl1_fx * ts_y_z_0[j];

                tdy_xy_z_0[j] = pa_x[j] * tdy_y_z_0[j];

                tdz_xy_z_0[j] = pa_x[j] * tdz_y_z_0[j];

                tdx_xz_x_0[j] = pa_x[j] * tdx_z_x_0[j] + 0.5 * fl1_fx * tdx_z_0_0[j] + 0.5 * fl1_fx * ts_z_x_0[j];

                tdy_xz_x_0[j] = pa_x[j] * tdy_z_x_0[j] + 0.5 * fl1_fx * tdy_z_0_0[j];

                tdz_xz_x_0[j] = pa_x[j] * tdz_z_x_0[j] + 0.5 * fl1_fx * tdz_z_0_0[j];

                tdx_xz_y_0[j] = pa_x[j] * tdx_z_y_0[j] + 0.5 * fl1_fx * ts_z_y_0[j];

                tdy_xz_y_0[j] = pa_x[j] * tdy_z_y_0[j];

                tdz_xz_y_0[j] = pa_x[j] * tdz_z_y_0[j];

                tdx_xz_z_0[j] = pa_x[j] * tdx_z_z_0[j] + 0.5 * fl1_fx * ts_z_z_0[j];

                tdy_xz_z_0[j] = pa_x[j] * tdy_z_z_0[j];

                tdz_xz_z_0[j] = pa_x[j] * tdz_z_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDP_27_54(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (27,54)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto tdx_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx); 

            auto tdy_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx); 

            auto tdz_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx); 

            auto tdx_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 1); 

            auto tdy_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 2); 

            auto tdy_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 2); 

            auto tdx_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 1); 

            auto tdy_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 2); 

            auto tdy_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 2); 

            auto ts_y_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 3); 

            auto ts_y_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 4); 

            auto ts_y_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 5); 

            auto ts_z_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 6); 

            auto ts_z_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 7); 

            auto ts_z_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 8); 

            // set up pointers to integrals

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

            auto tdx_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 15); 

            auto tdy_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 15); 

            auto tdz_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 15); 

            auto tdx_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 16); 

            auto tdy_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 16); 

            auto tdz_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 16); 

            auto tdx_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 17); 

            auto tdy_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 17); 

            auto tdz_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 17); 

            // Batch of Integrals (27,54)

            #pragma omp simd aligned(fx, pa_y, pa_z, tdx_0_x_0, tdx_0_y_0, tdx_0_z_0, tdx_y_0_0, tdx_y_x_0, \
                                     tdx_y_y_0, tdx_y_z_0, tdx_yy_x_0, tdx_yy_y_0, tdx_yy_z_0, tdx_yz_x_0, tdx_yz_y_0, \
                                     tdx_yz_z_0, tdx_z_0_0, tdx_z_x_0, tdx_z_y_0, tdx_z_z_0, tdx_zz_x_0, tdx_zz_y_0, \
                                     tdx_zz_z_0, tdy_0_x_0, tdy_0_y_0, tdy_0_z_0, tdy_y_0_0, tdy_y_x_0, tdy_y_y_0, \
                                     tdy_y_z_0, tdy_yy_x_0, tdy_yy_y_0, tdy_yy_z_0, tdy_yz_x_0, tdy_yz_y_0, tdy_yz_z_0, \
                                     tdy_z_0_0, tdy_z_x_0, tdy_z_y_0, tdy_z_z_0, tdy_zz_x_0, tdy_zz_y_0, tdy_zz_z_0, \
                                     tdz_0_x_0, tdz_0_y_0, tdz_0_z_0, tdz_y_0_0, tdz_y_x_0, tdz_y_y_0, tdz_y_z_0, \
                                     tdz_yy_x_0, tdz_yy_y_0, tdz_yy_z_0, tdz_yz_x_0, tdz_yz_y_0, tdz_yz_z_0, tdz_z_0_0, \
                                     tdz_z_x_0, tdz_z_y_0, tdz_z_z_0, tdz_zz_x_0, tdz_zz_y_0, tdz_zz_z_0, ts_y_x_0, \
                                     ts_y_y_0, ts_y_z_0, ts_z_x_0, ts_z_y_0, ts_z_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yy_x_0[j] = pa_y[j] * tdx_y_x_0[j] + 0.5 * fl1_fx * tdx_0_x_0[j];

                tdy_yy_x_0[j] = pa_y[j] * tdy_y_x_0[j] + 0.5 * fl1_fx * tdy_0_x_0[j] + 0.5 * fl1_fx * ts_y_x_0[j];

                tdz_yy_x_0[j] = pa_y[j] * tdz_y_x_0[j] + 0.5 * fl1_fx * tdz_0_x_0[j];

                tdx_yy_y_0[j] = pa_y[j] * tdx_y_y_0[j] + 0.5 * fl1_fx * tdx_0_y_0[j] + 0.5 * fl1_fx * tdx_y_0_0[j];

                tdy_yy_y_0[j] = pa_y[j] * tdy_y_y_0[j] + 0.5 * fl1_fx * tdy_0_y_0[j] + 0.5 * fl1_fx * tdy_y_0_0[j] + 0.5 * fl1_fx * ts_y_y_0[j];

                tdz_yy_y_0[j] = pa_y[j] * tdz_y_y_0[j] + 0.5 * fl1_fx * tdz_0_y_0[j] + 0.5 * fl1_fx * tdz_y_0_0[j];

                tdx_yy_z_0[j] = pa_y[j] * tdx_y_z_0[j] + 0.5 * fl1_fx * tdx_0_z_0[j];

                tdy_yy_z_0[j] = pa_y[j] * tdy_y_z_0[j] + 0.5 * fl1_fx * tdy_0_z_0[j] + 0.5 * fl1_fx * ts_y_z_0[j];

                tdz_yy_z_0[j] = pa_y[j] * tdz_y_z_0[j] + 0.5 * fl1_fx * tdz_0_z_0[j];

                tdx_yz_x_0[j] = pa_y[j] * tdx_z_x_0[j];

                tdy_yz_x_0[j] = pa_y[j] * tdy_z_x_0[j] + 0.5 * fl1_fx * ts_z_x_0[j];

                tdz_yz_x_0[j] = pa_y[j] * tdz_z_x_0[j];

                tdx_yz_y_0[j] = pa_y[j] * tdx_z_y_0[j] + 0.5 * fl1_fx * tdx_z_0_0[j];

                tdy_yz_y_0[j] = pa_y[j] * tdy_z_y_0[j] + 0.5 * fl1_fx * tdy_z_0_0[j] + 0.5 * fl1_fx * ts_z_y_0[j];

                tdz_yz_y_0[j] = pa_y[j] * tdz_z_y_0[j] + 0.5 * fl1_fx * tdz_z_0_0[j];

                tdx_yz_z_0[j] = pa_y[j] * tdx_z_z_0[j];

                tdy_yz_z_0[j] = pa_y[j] * tdy_z_z_0[j] + 0.5 * fl1_fx * ts_z_z_0[j];

                tdz_yz_z_0[j] = pa_y[j] * tdz_z_z_0[j];

                tdx_zz_x_0[j] = pa_z[j] * tdx_z_x_0[j] + 0.5 * fl1_fx * tdx_0_x_0[j];

                tdy_zz_x_0[j] = pa_z[j] * tdy_z_x_0[j] + 0.5 * fl1_fx * tdy_0_x_0[j];

                tdz_zz_x_0[j] = pa_z[j] * tdz_z_x_0[j] + 0.5 * fl1_fx * tdz_0_x_0[j] + 0.5 * fl1_fx * ts_z_x_0[j];

                tdx_zz_y_0[j] = pa_z[j] * tdx_z_y_0[j] + 0.5 * fl1_fx * tdx_0_y_0[j];

                tdy_zz_y_0[j] = pa_z[j] * tdy_z_y_0[j] + 0.5 * fl1_fx * tdy_0_y_0[j];

                tdz_zz_y_0[j] = pa_z[j] * tdz_z_y_0[j] + 0.5 * fl1_fx * tdz_0_y_0[j] + 0.5 * fl1_fx * ts_z_y_0[j];

                tdx_zz_z_0[j] = pa_z[j] * tdx_z_z_0[j] + 0.5 * fl1_fx * tdx_0_z_0[j] + 0.5 * fl1_fx * tdx_z_0_0[j];

                tdy_zz_z_0[j] = pa_z[j] * tdy_z_z_0[j] + 0.5 * fl1_fx * tdy_0_z_0[j] + 0.5 * fl1_fx * tdy_z_0_0[j];

                tdz_zz_z_0[j] = pa_z[j] * tdz_z_z_0[j] + 0.5 * fl1_fx * tdz_0_z_0[j] + 0.5 * fl1_fx * tdz_z_0_0[j] + 0.5 * fl1_fx * ts_z_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForPF(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForPF_0_45(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForPF_45_90(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 
    }

    void
    compElectricDipoleForPF_0_45(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_1_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

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

            auto ts_0_xxx_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx); 

            auto ts_0_xxy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 1); 

            auto ts_0_xxz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 2); 

            auto ts_0_xyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 3); 

            auto ts_0_xyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 4); 

            auto ts_0_xzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 5); 

            auto ts_0_yyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 6); 

            auto ts_0_yyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 7); 

            auto ts_0_yzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 8); 

            auto ts_0_zzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 9); 

            // set up pointers to integrals

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

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, pa_y, tdx_0_xx_0, tdx_0_xxx_0, tdx_0_xxy_0, tdx_0_xxz_0, \
                                     tdx_0_xy_0, tdx_0_xyy_0, tdx_0_xyz_0, tdx_0_xz_0, tdx_0_xzz_0, tdx_0_yy_0, \
                                     tdx_0_yyy_0, tdx_0_yyz_0, tdx_0_yz_0, tdx_0_yzz_0, tdx_0_zz_0, tdx_0_zzz_0, \
                                     tdx_x_xxx_0, tdx_x_xxy_0, tdx_x_xxz_0, tdx_x_xyy_0, tdx_x_xyz_0, tdx_x_xzz_0, \
                                     tdx_x_yyy_0, tdx_x_yyz_0, tdx_x_yzz_0, tdx_x_zzz_0, tdx_y_xxx_0, tdx_y_xxy_0, \
                                     tdx_y_xxz_0, tdx_y_xyy_0, tdx_y_xyz_0, tdy_0_xx_0, tdy_0_xxx_0, tdy_0_xxy_0, \
                                     tdy_0_xxz_0, tdy_0_xy_0, tdy_0_xyy_0, tdy_0_xyz_0, tdy_0_xz_0, tdy_0_xzz_0, \
                                     tdy_0_yy_0, tdy_0_yyy_0, tdy_0_yyz_0, tdy_0_yz_0, tdy_0_yzz_0, tdy_0_zz_0, \
                                     tdy_0_zzz_0, tdy_x_xxx_0, tdy_x_xxy_0, tdy_x_xxz_0, tdy_x_xyy_0, tdy_x_xyz_0, \
                                     tdy_x_xzz_0, tdy_x_yyy_0, tdy_x_yyz_0, tdy_x_yzz_0, tdy_x_zzz_0, tdy_y_xxx_0, \
                                     tdy_y_xxy_0, tdy_y_xxz_0, tdy_y_xyy_0, tdy_y_xyz_0, tdz_0_xx_0, tdz_0_xxx_0, \
                                     tdz_0_xxy_0, tdz_0_xxz_0, tdz_0_xy_0, tdz_0_xyy_0, tdz_0_xyz_0, tdz_0_xz_0, \
                                     tdz_0_xzz_0, tdz_0_yy_0, tdz_0_yyy_0, tdz_0_yyz_0, tdz_0_yz_0, tdz_0_yzz_0, \
                                     tdz_0_zz_0, tdz_0_zzz_0, tdz_x_xxx_0, tdz_x_xxy_0, tdz_x_xxz_0, tdz_x_xyy_0, \
                                     tdz_x_xyz_0, tdz_x_xzz_0, tdz_x_yyy_0, tdz_x_yyz_0, tdz_x_yzz_0, tdz_x_zzz_0, \
                                     tdz_y_xxx_0, tdz_y_xxy_0, tdz_y_xxz_0, tdz_y_xyy_0, tdz_y_xyz_0, ts_0_xxx_0, \
                                     ts_0_xxy_0, ts_0_xxz_0, ts_0_xyy_0, ts_0_xyz_0, ts_0_xzz_0, ts_0_yyy_0, ts_0_yyz_0, \
                                     ts_0_yzz_0, ts_0_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_x_xxx_0[j] = pa_x[j] * tdx_0_xxx_0[j] + 1.5 * fl1_fx * tdx_0_xx_0[j] + 0.5 * fl1_fx * ts_0_xxx_0[j];

                tdy_x_xxx_0[j] = pa_x[j] * tdy_0_xxx_0[j] + 1.5 * fl1_fx * tdy_0_xx_0[j];

                tdz_x_xxx_0[j] = pa_x[j] * tdz_0_xxx_0[j] + 1.5 * fl1_fx * tdz_0_xx_0[j];

                tdx_x_xxy_0[j] = pa_x[j] * tdx_0_xxy_0[j] + fl1_fx * tdx_0_xy_0[j] + 0.5 * fl1_fx * ts_0_xxy_0[j];

                tdy_x_xxy_0[j] = pa_x[j] * tdy_0_xxy_0[j] + fl1_fx * tdy_0_xy_0[j];

                tdz_x_xxy_0[j] = pa_x[j] * tdz_0_xxy_0[j] + fl1_fx * tdz_0_xy_0[j];

                tdx_x_xxz_0[j] = pa_x[j] * tdx_0_xxz_0[j] + fl1_fx * tdx_0_xz_0[j] + 0.5 * fl1_fx * ts_0_xxz_0[j];

                tdy_x_xxz_0[j] = pa_x[j] * tdy_0_xxz_0[j] + fl1_fx * tdy_0_xz_0[j];

                tdz_x_xxz_0[j] = pa_x[j] * tdz_0_xxz_0[j] + fl1_fx * tdz_0_xz_0[j];

                tdx_x_xyy_0[j] = pa_x[j] * tdx_0_xyy_0[j] + 0.5 * fl1_fx * tdx_0_yy_0[j] + 0.5 * fl1_fx * ts_0_xyy_0[j];

                tdy_x_xyy_0[j] = pa_x[j] * tdy_0_xyy_0[j] + 0.5 * fl1_fx * tdy_0_yy_0[j];

                tdz_x_xyy_0[j] = pa_x[j] * tdz_0_xyy_0[j] + 0.5 * fl1_fx * tdz_0_yy_0[j];

                tdx_x_xyz_0[j] = pa_x[j] * tdx_0_xyz_0[j] + 0.5 * fl1_fx * tdx_0_yz_0[j] + 0.5 * fl1_fx * ts_0_xyz_0[j];

                tdy_x_xyz_0[j] = pa_x[j] * tdy_0_xyz_0[j] + 0.5 * fl1_fx * tdy_0_yz_0[j];

                tdz_x_xyz_0[j] = pa_x[j] * tdz_0_xyz_0[j] + 0.5 * fl1_fx * tdz_0_yz_0[j];

                tdx_x_xzz_0[j] = pa_x[j] * tdx_0_xzz_0[j] + 0.5 * fl1_fx * tdx_0_zz_0[j] + 0.5 * fl1_fx * ts_0_xzz_0[j];

                tdy_x_xzz_0[j] = pa_x[j] * tdy_0_xzz_0[j] + 0.5 * fl1_fx * tdy_0_zz_0[j];

                tdz_x_xzz_0[j] = pa_x[j] * tdz_0_xzz_0[j] + 0.5 * fl1_fx * tdz_0_zz_0[j];

                tdx_x_yyy_0[j] = pa_x[j] * tdx_0_yyy_0[j] + 0.5 * fl1_fx * ts_0_yyy_0[j];

                tdy_x_yyy_0[j] = pa_x[j] * tdy_0_yyy_0[j];

                tdz_x_yyy_0[j] = pa_x[j] * tdz_0_yyy_0[j];

                tdx_x_yyz_0[j] = pa_x[j] * tdx_0_yyz_0[j] + 0.5 * fl1_fx * ts_0_yyz_0[j];

                tdy_x_yyz_0[j] = pa_x[j] * tdy_0_yyz_0[j];

                tdz_x_yyz_0[j] = pa_x[j] * tdz_0_yyz_0[j];

                tdx_x_yzz_0[j] = pa_x[j] * tdx_0_yzz_0[j] + 0.5 * fl1_fx * ts_0_yzz_0[j];

                tdy_x_yzz_0[j] = pa_x[j] * tdy_0_yzz_0[j];

                tdz_x_yzz_0[j] = pa_x[j] * tdz_0_yzz_0[j];

                tdx_x_zzz_0[j] = pa_x[j] * tdx_0_zzz_0[j] + 0.5 * fl1_fx * ts_0_zzz_0[j];

                tdy_x_zzz_0[j] = pa_x[j] * tdy_0_zzz_0[j];

                tdz_x_zzz_0[j] = pa_x[j] * tdz_0_zzz_0[j];

                tdx_y_xxx_0[j] = pa_y[j] * tdx_0_xxx_0[j];

                tdy_y_xxx_0[j] = pa_y[j] * tdy_0_xxx_0[j] + 0.5 * fl1_fx * ts_0_xxx_0[j];

                tdz_y_xxx_0[j] = pa_y[j] * tdz_0_xxx_0[j];

                tdx_y_xxy_0[j] = pa_y[j] * tdx_0_xxy_0[j] + 0.5 * fl1_fx * tdx_0_xx_0[j];

                tdy_y_xxy_0[j] = pa_y[j] * tdy_0_xxy_0[j] + 0.5 * fl1_fx * tdy_0_xx_0[j] + 0.5 * fl1_fx * ts_0_xxy_0[j];

                tdz_y_xxy_0[j] = pa_y[j] * tdz_0_xxy_0[j] + 0.5 * fl1_fx * tdz_0_xx_0[j];

                tdx_y_xxz_0[j] = pa_y[j] * tdx_0_xxz_0[j];

                tdy_y_xxz_0[j] = pa_y[j] * tdy_0_xxz_0[j] + 0.5 * fl1_fx * ts_0_xxz_0[j];

                tdz_y_xxz_0[j] = pa_y[j] * tdz_0_xxz_0[j];

                tdx_y_xyy_0[j] = pa_y[j] * tdx_0_xyy_0[j] + fl1_fx * tdx_0_xy_0[j];

                tdy_y_xyy_0[j] = pa_y[j] * tdy_0_xyy_0[j] + fl1_fx * tdy_0_xy_0[j] + 0.5 * fl1_fx * ts_0_xyy_0[j];

                tdz_y_xyy_0[j] = pa_y[j] * tdz_0_xyy_0[j] + fl1_fx * tdz_0_xy_0[j];

                tdx_y_xyz_0[j] = pa_y[j] * tdx_0_xyz_0[j] + 0.5 * fl1_fx * tdx_0_xz_0[j];

                tdy_y_xyz_0[j] = pa_y[j] * tdy_0_xyz_0[j] + 0.5 * fl1_fx * tdy_0_xz_0[j] + 0.5 * fl1_fx * ts_0_xyz_0[j];

                tdz_y_xyz_0[j] = pa_y[j] * tdz_0_xyz_0[j] + 0.5 * fl1_fx * tdz_0_xz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForPF_45_90(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_1_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto ts_0_xxx_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx); 

            auto ts_0_xxy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 1); 

            auto ts_0_xxz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 2); 

            auto ts_0_xyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 3); 

            auto ts_0_xyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 4); 

            auto ts_0_xzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 5); 

            auto ts_0_yyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 6); 

            auto ts_0_yyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 7); 

            auto ts_0_yzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 8); 

            auto ts_0_zzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 9); 

            // set up pointers to integrals

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

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_y, pa_z, tdx_0_xx_0, tdx_0_xxx_0, tdx_0_xxy_0, tdx_0_xxz_0, \
                                     tdx_0_xy_0, tdx_0_xyy_0, tdx_0_xyz_0, tdx_0_xz_0, tdx_0_xzz_0, tdx_0_yy_0, \
                                     tdx_0_yyy_0, tdx_0_yyz_0, tdx_0_yz_0, tdx_0_yzz_0, tdx_0_zz_0, tdx_0_zzz_0, \
                                     tdx_y_xzz_0, tdx_y_yyy_0, tdx_y_yyz_0, tdx_y_yzz_0, tdx_y_zzz_0, tdx_z_xxx_0, \
                                     tdx_z_xxy_0, tdx_z_xxz_0, tdx_z_xyy_0, tdx_z_xyz_0, tdx_z_xzz_0, tdx_z_yyy_0, \
                                     tdx_z_yyz_0, tdx_z_yzz_0, tdx_z_zzz_0, tdy_0_xx_0, tdy_0_xxx_0, tdy_0_xxy_0, \
                                     tdy_0_xxz_0, tdy_0_xy_0, tdy_0_xyy_0, tdy_0_xyz_0, tdy_0_xz_0, tdy_0_xzz_0, \
                                     tdy_0_yy_0, tdy_0_yyy_0, tdy_0_yyz_0, tdy_0_yz_0, tdy_0_yzz_0, tdy_0_zz_0, \
                                     tdy_0_zzz_0, tdy_y_xzz_0, tdy_y_yyy_0, tdy_y_yyz_0, tdy_y_yzz_0, tdy_y_zzz_0, \
                                     tdy_z_xxx_0, tdy_z_xxy_0, tdy_z_xxz_0, tdy_z_xyy_0, tdy_z_xyz_0, tdy_z_xzz_0, \
                                     tdy_z_yyy_0, tdy_z_yyz_0, tdy_z_yzz_0, tdy_z_zzz_0, tdz_0_xx_0, tdz_0_xxx_0, \
                                     tdz_0_xxy_0, tdz_0_xxz_0, tdz_0_xy_0, tdz_0_xyy_0, tdz_0_xyz_0, tdz_0_xz_0, \
                                     tdz_0_xzz_0, tdz_0_yy_0, tdz_0_yyy_0, tdz_0_yyz_0, tdz_0_yz_0, tdz_0_yzz_0, \
                                     tdz_0_zz_0, tdz_0_zzz_0, tdz_y_xzz_0, tdz_y_yyy_0, tdz_y_yyz_0, tdz_y_yzz_0, \
                                     tdz_y_zzz_0, tdz_z_xxx_0, tdz_z_xxy_0, tdz_z_xxz_0, tdz_z_xyy_0, tdz_z_xyz_0, \
                                     tdz_z_xzz_0, tdz_z_yyy_0, tdz_z_yyz_0, tdz_z_yzz_0, tdz_z_zzz_0, ts_0_xxx_0, \
                                     ts_0_xxy_0, ts_0_xxz_0, ts_0_xyy_0, ts_0_xyz_0, ts_0_xzz_0, ts_0_yyy_0, ts_0_yyz_0, \
                                     ts_0_yzz_0, ts_0_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_y_xzz_0[j] = pa_y[j] * tdx_0_xzz_0[j];

                tdy_y_xzz_0[j] = pa_y[j] * tdy_0_xzz_0[j] + 0.5 * fl1_fx * ts_0_xzz_0[j];

                tdz_y_xzz_0[j] = pa_y[j] * tdz_0_xzz_0[j];

                tdx_y_yyy_0[j] = pa_y[j] * tdx_0_yyy_0[j] + 1.5 * fl1_fx * tdx_0_yy_0[j];

                tdy_y_yyy_0[j] = pa_y[j] * tdy_0_yyy_0[j] + 1.5 * fl1_fx * tdy_0_yy_0[j] + 0.5 * fl1_fx * ts_0_yyy_0[j];

                tdz_y_yyy_0[j] = pa_y[j] * tdz_0_yyy_0[j] + 1.5 * fl1_fx * tdz_0_yy_0[j];

                tdx_y_yyz_0[j] = pa_y[j] * tdx_0_yyz_0[j] + fl1_fx * tdx_0_yz_0[j];

                tdy_y_yyz_0[j] = pa_y[j] * tdy_0_yyz_0[j] + fl1_fx * tdy_0_yz_0[j] + 0.5 * fl1_fx * ts_0_yyz_0[j];

                tdz_y_yyz_0[j] = pa_y[j] * tdz_0_yyz_0[j] + fl1_fx * tdz_0_yz_0[j];

                tdx_y_yzz_0[j] = pa_y[j] * tdx_0_yzz_0[j] + 0.5 * fl1_fx * tdx_0_zz_0[j];

                tdy_y_yzz_0[j] = pa_y[j] * tdy_0_yzz_0[j] + 0.5 * fl1_fx * tdy_0_zz_0[j] + 0.5 * fl1_fx * ts_0_yzz_0[j];

                tdz_y_yzz_0[j] = pa_y[j] * tdz_0_yzz_0[j] + 0.5 * fl1_fx * tdz_0_zz_0[j];

                tdx_y_zzz_0[j] = pa_y[j] * tdx_0_zzz_0[j];

                tdy_y_zzz_0[j] = pa_y[j] * tdy_0_zzz_0[j] + 0.5 * fl1_fx * ts_0_zzz_0[j];

                tdz_y_zzz_0[j] = pa_y[j] * tdz_0_zzz_0[j];

                tdx_z_xxx_0[j] = pa_z[j] * tdx_0_xxx_0[j];

                tdy_z_xxx_0[j] = pa_z[j] * tdy_0_xxx_0[j];

                tdz_z_xxx_0[j] = pa_z[j] * tdz_0_xxx_0[j] + 0.5 * fl1_fx * ts_0_xxx_0[j];

                tdx_z_xxy_0[j] = pa_z[j] * tdx_0_xxy_0[j];

                tdy_z_xxy_0[j] = pa_z[j] * tdy_0_xxy_0[j];

                tdz_z_xxy_0[j] = pa_z[j] * tdz_0_xxy_0[j] + 0.5 * fl1_fx * ts_0_xxy_0[j];

                tdx_z_xxz_0[j] = pa_z[j] * tdx_0_xxz_0[j] + 0.5 * fl1_fx * tdx_0_xx_0[j];

                tdy_z_xxz_0[j] = pa_z[j] * tdy_0_xxz_0[j] + 0.5 * fl1_fx * tdy_0_xx_0[j];

                tdz_z_xxz_0[j] = pa_z[j] * tdz_0_xxz_0[j] + 0.5 * fl1_fx * tdz_0_xx_0[j] + 0.5 * fl1_fx * ts_0_xxz_0[j];

                tdx_z_xyy_0[j] = pa_z[j] * tdx_0_xyy_0[j];

                tdy_z_xyy_0[j] = pa_z[j] * tdy_0_xyy_0[j];

                tdz_z_xyy_0[j] = pa_z[j] * tdz_0_xyy_0[j] + 0.5 * fl1_fx * ts_0_xyy_0[j];

                tdx_z_xyz_0[j] = pa_z[j] * tdx_0_xyz_0[j] + 0.5 * fl1_fx * tdx_0_xy_0[j];

                tdy_z_xyz_0[j] = pa_z[j] * tdy_0_xyz_0[j] + 0.5 * fl1_fx * tdy_0_xy_0[j];

                tdz_z_xyz_0[j] = pa_z[j] * tdz_0_xyz_0[j] + 0.5 * fl1_fx * tdz_0_xy_0[j] + 0.5 * fl1_fx * ts_0_xyz_0[j];

                tdx_z_xzz_0[j] = pa_z[j] * tdx_0_xzz_0[j] + fl1_fx * tdx_0_xz_0[j];

                tdy_z_xzz_0[j] = pa_z[j] * tdy_0_xzz_0[j] + fl1_fx * tdy_0_xz_0[j];

                tdz_z_xzz_0[j] = pa_z[j] * tdz_0_xzz_0[j] + fl1_fx * tdz_0_xz_0[j] + 0.5 * fl1_fx * ts_0_xzz_0[j];

                tdx_z_yyy_0[j] = pa_z[j] * tdx_0_yyy_0[j];

                tdy_z_yyy_0[j] = pa_z[j] * tdy_0_yyy_0[j];

                tdz_z_yyy_0[j] = pa_z[j] * tdz_0_yyy_0[j] + 0.5 * fl1_fx * ts_0_yyy_0[j];

                tdx_z_yyz_0[j] = pa_z[j] * tdx_0_yyz_0[j] + 0.5 * fl1_fx * tdx_0_yy_0[j];

                tdy_z_yyz_0[j] = pa_z[j] * tdy_0_yyz_0[j] + 0.5 * fl1_fx * tdy_0_yy_0[j];

                tdz_z_yyz_0[j] = pa_z[j] * tdz_0_yyz_0[j] + 0.5 * fl1_fx * tdz_0_yy_0[j] + 0.5 * fl1_fx * ts_0_yyz_0[j];

                tdx_z_yzz_0[j] = pa_z[j] * tdx_0_yzz_0[j] + fl1_fx * tdx_0_yz_0[j];

                tdy_z_yzz_0[j] = pa_z[j] * tdy_0_yzz_0[j] + fl1_fx * tdy_0_yz_0[j];

                tdz_z_yzz_0[j] = pa_z[j] * tdz_0_yzz_0[j] + fl1_fx * tdz_0_yz_0[j] + 0.5 * fl1_fx * ts_0_yzz_0[j];

                tdx_z_zzz_0[j] = pa_z[j] * tdx_0_zzz_0[j] + 1.5 * fl1_fx * tdx_0_zz_0[j];

                tdy_z_zzz_0[j] = pa_z[j] * tdy_0_zzz_0[j] + 1.5 * fl1_fx * tdy_0_zz_0[j];

                tdz_z_zzz_0[j] = pa_z[j] * tdz_0_zzz_0[j] + 1.5 * fl1_fx * tdz_0_zz_0[j] + 0.5 * fl1_fx * ts_0_zzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForFP(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForFP_0_45(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForFP_45_90(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 
    }

    void
    compElectricDipoleForFP_0_45(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_3_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

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

            auto tdx_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 6); 

            auto tdy_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 6); 

            auto tdz_z_x_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 6); 

            auto tdx_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 7); 

            auto tdy_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 7); 

            auto tdz_z_y_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 7); 

            auto tdx_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * idx + 8); 

            auto tdy_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 9 * bdim + 9 * idx + 8); 

            auto tdz_z_z_0 = primBuffer.data(pidx_d_1_1_m0 + 18 * bdim + 9 * idx + 8); 

            auto tdx_xx_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx); 

            auto tdy_xx_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx); 

            auto tdz_xx_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx); 

            auto tdx_xy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 1); 

            auto tdy_xy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 1); 

            auto tdz_xy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 1); 

            auto tdx_xz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 2); 

            auto tdy_xz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 2); 

            auto tdz_xz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 2); 

            auto tdx_yy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 3); 

            auto tdy_yy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 3); 

            auto tdz_yy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 3); 

            auto tdx_yz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 4); 

            auto tdy_yz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 4); 

            auto tdz_yz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 4); 

            auto ts_xx_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx); 

            auto ts_xx_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 1); 

            auto ts_xx_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 2); 

            auto ts_xy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 3); 

            auto ts_xy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 4); 

            auto ts_xy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 5); 

            auto ts_xz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 6); 

            auto ts_xz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 7); 

            auto ts_xz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 8); 

            auto ts_yy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 9); 

            auto ts_yy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 10); 

            auto ts_yy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 11); 

            auto ts_yz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 12); 

            auto ts_yz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 13); 

            auto ts_yz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 14); 

            // set up pointers to integrals

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

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, tdx_x_x_0, tdx_x_y_0, tdx_x_z_0, tdx_xx_0_0, tdx_xx_x_0, \
                                     tdx_xx_y_0, tdx_xx_z_0, tdx_xxx_x_0, tdx_xxx_y_0, tdx_xxx_z_0, tdx_xxy_x_0, \
                                     tdx_xxy_y_0, tdx_xxy_z_0, tdx_xxz_x_0, tdx_xxz_y_0, tdx_xxz_z_0, tdx_xy_0_0, \
                                     tdx_xy_x_0, tdx_xy_y_0, tdx_xy_z_0, tdx_xyy_x_0, tdx_xyy_y_0, tdx_xyy_z_0, \
                                     tdx_xyz_x_0, tdx_xyz_y_0, tdx_xyz_z_0, tdx_xz_0_0, tdx_xz_x_0, tdx_xz_y_0, \
                                     tdx_xz_z_0, tdx_y_x_0, tdx_y_y_0, tdx_y_z_0, tdx_yy_0_0, tdx_yy_x_0, tdx_yy_y_0, \
                                     tdx_yy_z_0, tdx_yz_0_0, tdx_yz_x_0, tdx_yz_y_0, tdx_yz_z_0, tdx_z_x_0, tdx_z_y_0, \
                                     tdx_z_z_0, tdy_x_x_0, tdy_x_y_0, tdy_x_z_0, tdy_xx_0_0, tdy_xx_x_0, tdy_xx_y_0, \
                                     tdy_xx_z_0, tdy_xxx_x_0, tdy_xxx_y_0, tdy_xxx_z_0, tdy_xxy_x_0, tdy_xxy_y_0, \
                                     tdy_xxy_z_0, tdy_xxz_x_0, tdy_xxz_y_0, tdy_xxz_z_0, tdy_xy_0_0, tdy_xy_x_0, \
                                     tdy_xy_y_0, tdy_xy_z_0, tdy_xyy_x_0, tdy_xyy_y_0, tdy_xyy_z_0, tdy_xyz_x_0, \
                                     tdy_xyz_y_0, tdy_xyz_z_0, tdy_xz_0_0, tdy_xz_x_0, tdy_xz_y_0, tdy_xz_z_0, tdy_y_x_0, \
                                     tdy_y_y_0, tdy_y_z_0, tdy_yy_0_0, tdy_yy_x_0, tdy_yy_y_0, tdy_yy_z_0, tdy_yz_0_0, \
                                     tdy_yz_x_0, tdy_yz_y_0, tdy_yz_z_0, tdy_z_x_0, tdy_z_y_0, tdy_z_z_0, tdz_x_x_0, \
                                     tdz_x_y_0, tdz_x_z_0, tdz_xx_0_0, tdz_xx_x_0, tdz_xx_y_0, tdz_xx_z_0, tdz_xxx_x_0, \
                                     tdz_xxx_y_0, tdz_xxx_z_0, tdz_xxy_x_0, tdz_xxy_y_0, tdz_xxy_z_0, tdz_xxz_x_0, \
                                     tdz_xxz_y_0, tdz_xxz_z_0, tdz_xy_0_0, tdz_xy_x_0, tdz_xy_y_0, tdz_xy_z_0, \
                                     tdz_xyy_x_0, tdz_xyy_y_0, tdz_xyy_z_0, tdz_xyz_x_0, tdz_xyz_y_0, tdz_xyz_z_0, \
                                     tdz_xz_0_0, tdz_xz_x_0, tdz_xz_y_0, tdz_xz_z_0, tdz_y_x_0, tdz_y_y_0, tdz_y_z_0, \
                                     tdz_yy_0_0, tdz_yy_x_0, tdz_yy_y_0, tdz_yy_z_0, tdz_yz_0_0, tdz_yz_x_0, tdz_yz_y_0, \
                                     tdz_yz_z_0, tdz_z_x_0, tdz_z_y_0, tdz_z_z_0, ts_xx_x_0, ts_xx_y_0, ts_xx_z_0, \
                                     ts_xy_x_0, ts_xy_y_0, ts_xy_z_0, ts_xz_x_0, ts_xz_y_0, ts_xz_z_0, ts_yy_x_0, \
                                     ts_yy_y_0, ts_yy_z_0, ts_yz_x_0, ts_yz_y_0, ts_yz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxx_x_0[j] = pa_x[j] * tdx_xx_x_0[j] + fl1_fx * tdx_x_x_0[j] + 0.5 * fl1_fx * tdx_xx_0_0[j] + 0.5 * fl1_fx * ts_xx_x_0[j];

                tdy_xxx_x_0[j] = pa_x[j] * tdy_xx_x_0[j] + fl1_fx * tdy_x_x_0[j] + 0.5 * fl1_fx * tdy_xx_0_0[j];

                tdz_xxx_x_0[j] = pa_x[j] * tdz_xx_x_0[j] + fl1_fx * tdz_x_x_0[j] + 0.5 * fl1_fx * tdz_xx_0_0[j];

                tdx_xxx_y_0[j] = pa_x[j] * tdx_xx_y_0[j] + fl1_fx * tdx_x_y_0[j] + 0.5 * fl1_fx * ts_xx_y_0[j];

                tdy_xxx_y_0[j] = pa_x[j] * tdy_xx_y_0[j] + fl1_fx * tdy_x_y_0[j];

                tdz_xxx_y_0[j] = pa_x[j] * tdz_xx_y_0[j] + fl1_fx * tdz_x_y_0[j];

                tdx_xxx_z_0[j] = pa_x[j] * tdx_xx_z_0[j] + fl1_fx * tdx_x_z_0[j] + 0.5 * fl1_fx * ts_xx_z_0[j];

                tdy_xxx_z_0[j] = pa_x[j] * tdy_xx_z_0[j] + fl1_fx * tdy_x_z_0[j];

                tdz_xxx_z_0[j] = pa_x[j] * tdz_xx_z_0[j] + fl1_fx * tdz_x_z_0[j];

                tdx_xxy_x_0[j] = pa_x[j] * tdx_xy_x_0[j] + 0.5 * fl1_fx * tdx_y_x_0[j] + 0.5 * fl1_fx * tdx_xy_0_0[j] + 0.5 * fl1_fx * ts_xy_x_0[j];

                tdy_xxy_x_0[j] = pa_x[j] * tdy_xy_x_0[j] + 0.5 * fl1_fx * tdy_y_x_0[j] + 0.5 * fl1_fx * tdy_xy_0_0[j];

                tdz_xxy_x_0[j] = pa_x[j] * tdz_xy_x_0[j] + 0.5 * fl1_fx * tdz_y_x_0[j] + 0.5 * fl1_fx * tdz_xy_0_0[j];

                tdx_xxy_y_0[j] = pa_x[j] * tdx_xy_y_0[j] + 0.5 * fl1_fx * tdx_y_y_0[j] + 0.5 * fl1_fx * ts_xy_y_0[j];

                tdy_xxy_y_0[j] = pa_x[j] * tdy_xy_y_0[j] + 0.5 * fl1_fx * tdy_y_y_0[j];

                tdz_xxy_y_0[j] = pa_x[j] * tdz_xy_y_0[j] + 0.5 * fl1_fx * tdz_y_y_0[j];

                tdx_xxy_z_0[j] = pa_x[j] * tdx_xy_z_0[j] + 0.5 * fl1_fx * tdx_y_z_0[j] + 0.5 * fl1_fx * ts_xy_z_0[j];

                tdy_xxy_z_0[j] = pa_x[j] * tdy_xy_z_0[j] + 0.5 * fl1_fx * tdy_y_z_0[j];

                tdz_xxy_z_0[j] = pa_x[j] * tdz_xy_z_0[j] + 0.5 * fl1_fx * tdz_y_z_0[j];

                tdx_xxz_x_0[j] = pa_x[j] * tdx_xz_x_0[j] + 0.5 * fl1_fx * tdx_z_x_0[j] + 0.5 * fl1_fx * tdx_xz_0_0[j] + 0.5 * fl1_fx * ts_xz_x_0[j];

                tdy_xxz_x_0[j] = pa_x[j] * tdy_xz_x_0[j] + 0.5 * fl1_fx * tdy_z_x_0[j] + 0.5 * fl1_fx * tdy_xz_0_0[j];

                tdz_xxz_x_0[j] = pa_x[j] * tdz_xz_x_0[j] + 0.5 * fl1_fx * tdz_z_x_0[j] + 0.5 * fl1_fx * tdz_xz_0_0[j];

                tdx_xxz_y_0[j] = pa_x[j] * tdx_xz_y_0[j] + 0.5 * fl1_fx * tdx_z_y_0[j] + 0.5 * fl1_fx * ts_xz_y_0[j];

                tdy_xxz_y_0[j] = pa_x[j] * tdy_xz_y_0[j] + 0.5 * fl1_fx * tdy_z_y_0[j];

                tdz_xxz_y_0[j] = pa_x[j] * tdz_xz_y_0[j] + 0.5 * fl1_fx * tdz_z_y_0[j];

                tdx_xxz_z_0[j] = pa_x[j] * tdx_xz_z_0[j] + 0.5 * fl1_fx * tdx_z_z_0[j] + 0.5 * fl1_fx * ts_xz_z_0[j];

                tdy_xxz_z_0[j] = pa_x[j] * tdy_xz_z_0[j] + 0.5 * fl1_fx * tdy_z_z_0[j];

                tdz_xxz_z_0[j] = pa_x[j] * tdz_xz_z_0[j] + 0.5 * fl1_fx * tdz_z_z_0[j];

                tdx_xyy_x_0[j] = pa_x[j] * tdx_yy_x_0[j] + 0.5 * fl1_fx * tdx_yy_0_0[j] + 0.5 * fl1_fx * ts_yy_x_0[j];

                tdy_xyy_x_0[j] = pa_x[j] * tdy_yy_x_0[j] + 0.5 * fl1_fx * tdy_yy_0_0[j];

                tdz_xyy_x_0[j] = pa_x[j] * tdz_yy_x_0[j] + 0.5 * fl1_fx * tdz_yy_0_0[j];

                tdx_xyy_y_0[j] = pa_x[j] * tdx_yy_y_0[j] + 0.5 * fl1_fx * ts_yy_y_0[j];

                tdy_xyy_y_0[j] = pa_x[j] * tdy_yy_y_0[j];

                tdz_xyy_y_0[j] = pa_x[j] * tdz_yy_y_0[j];

                tdx_xyy_z_0[j] = pa_x[j] * tdx_yy_z_0[j] + 0.5 * fl1_fx * ts_yy_z_0[j];

                tdy_xyy_z_0[j] = pa_x[j] * tdy_yy_z_0[j];

                tdz_xyy_z_0[j] = pa_x[j] * tdz_yy_z_0[j];

                tdx_xyz_x_0[j] = pa_x[j] * tdx_yz_x_0[j] + 0.5 * fl1_fx * tdx_yz_0_0[j] + 0.5 * fl1_fx * ts_yz_x_0[j];

                tdy_xyz_x_0[j] = pa_x[j] * tdy_yz_x_0[j] + 0.5 * fl1_fx * tdy_yz_0_0[j];

                tdz_xyz_x_0[j] = pa_x[j] * tdz_yz_x_0[j] + 0.5 * fl1_fx * tdz_yz_0_0[j];

                tdx_xyz_y_0[j] = pa_x[j] * tdx_yz_y_0[j] + 0.5 * fl1_fx * ts_yz_y_0[j];

                tdy_xyz_y_0[j] = pa_x[j] * tdy_yz_y_0[j];

                tdz_xyz_y_0[j] = pa_x[j] * tdz_yz_y_0[j];

                tdx_xyz_z_0[j] = pa_x[j] * tdx_yz_z_0[j] + 0.5 * fl1_fx * ts_yz_z_0[j];

                tdy_xyz_z_0[j] = pa_x[j] * tdy_yz_z_0[j];

                tdz_xyz_z_0[j] = pa_x[j] * tdz_yz_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForFP_45_90(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_3_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto tdx_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 15); 

            auto tdy_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 15); 

            auto tdz_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 15); 

            auto tdx_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 16); 

            auto tdy_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 16); 

            auto tdz_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 16); 

            auto tdx_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 17); 

            auto tdy_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 17); 

            auto tdz_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 17); 

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

            auto tdx_yy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 3); 

            auto tdy_yy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 3); 

            auto tdz_yy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 3); 

            auto tdx_yz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 4); 

            auto tdy_yz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 4); 

            auto tdz_yz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 4); 

            auto tdx_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 5); 

            auto tdy_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 5); 

            auto tdz_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 5); 

            auto ts_yy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 9); 

            auto ts_yy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 10); 

            auto ts_yy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 11); 

            auto ts_yz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 12); 

            auto ts_yz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 13); 

            auto ts_yz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 14); 

            auto ts_zz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 15); 

            auto ts_zz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 16); 

            auto ts_zz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 17); 

            // set up pointers to integrals

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

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, tdx_xzz_x_0, tdx_xzz_y_0, tdx_xzz_z_0, tdx_y_x_0, \
                                     tdx_y_y_0, tdx_y_z_0, tdx_yy_0_0, tdx_yy_x_0, tdx_yy_y_0, tdx_yy_z_0, tdx_yyy_x_0, \
                                     tdx_yyy_y_0, tdx_yyy_z_0, tdx_yyz_x_0, tdx_yyz_y_0, tdx_yyz_z_0, tdx_yz_0_0, \
                                     tdx_yz_x_0, tdx_yz_y_0, tdx_yz_z_0, tdx_yzz_x_0, tdx_yzz_y_0, tdx_yzz_z_0, \
                                     tdx_z_x_0, tdx_z_y_0, tdx_z_z_0, tdx_zz_0_0, tdx_zz_x_0, tdx_zz_y_0, tdx_zz_z_0, \
                                     tdx_zzz_x_0, tdx_zzz_y_0, tdx_zzz_z_0, tdy_xzz_x_0, tdy_xzz_y_0, tdy_xzz_z_0, \
                                     tdy_y_x_0, tdy_y_y_0, tdy_y_z_0, tdy_yy_0_0, tdy_yy_x_0, tdy_yy_y_0, tdy_yy_z_0, \
                                     tdy_yyy_x_0, tdy_yyy_y_0, tdy_yyy_z_0, tdy_yyz_x_0, tdy_yyz_y_0, tdy_yyz_z_0, \
                                     tdy_yz_0_0, tdy_yz_x_0, tdy_yz_y_0, tdy_yz_z_0, tdy_yzz_x_0, tdy_yzz_y_0, \
                                     tdy_yzz_z_0, tdy_z_x_0, tdy_z_y_0, tdy_z_z_0, tdy_zz_0_0, tdy_zz_x_0, tdy_zz_y_0, \
                                     tdy_zz_z_0, tdy_zzz_x_0, tdy_zzz_y_0, tdy_zzz_z_0, tdz_xzz_x_0, tdz_xzz_y_0, \
                                     tdz_xzz_z_0, tdz_y_x_0, tdz_y_y_0, tdz_y_z_0, tdz_yy_0_0, tdz_yy_x_0, tdz_yy_y_0, \
                                     tdz_yy_z_0, tdz_yyy_x_0, tdz_yyy_y_0, tdz_yyy_z_0, tdz_yyz_x_0, tdz_yyz_y_0, \
                                     tdz_yyz_z_0, tdz_yz_0_0, tdz_yz_x_0, tdz_yz_y_0, tdz_yz_z_0, tdz_yzz_x_0, \
                                     tdz_yzz_y_0, tdz_yzz_z_0, tdz_z_x_0, tdz_z_y_0, tdz_z_z_0, tdz_zz_0_0, tdz_zz_x_0, \
                                     tdz_zz_y_0, tdz_zz_z_0, tdz_zzz_x_0, tdz_zzz_y_0, tdz_zzz_z_0, ts_yy_x_0, \
                                     ts_yy_y_0, ts_yy_z_0, ts_yz_x_0, ts_yz_y_0, ts_yz_z_0, ts_zz_x_0, ts_zz_y_0, \
                                     ts_zz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xzz_x_0[j] = pa_x[j] * tdx_zz_x_0[j] + 0.5 * fl1_fx * tdx_zz_0_0[j] + 0.5 * fl1_fx * ts_zz_x_0[j];

                tdy_xzz_x_0[j] = pa_x[j] * tdy_zz_x_0[j] + 0.5 * fl1_fx * tdy_zz_0_0[j];

                tdz_xzz_x_0[j] = pa_x[j] * tdz_zz_x_0[j] + 0.5 * fl1_fx * tdz_zz_0_0[j];

                tdx_xzz_y_0[j] = pa_x[j] * tdx_zz_y_0[j] + 0.5 * fl1_fx * ts_zz_y_0[j];

                tdy_xzz_y_0[j] = pa_x[j] * tdy_zz_y_0[j];

                tdz_xzz_y_0[j] = pa_x[j] * tdz_zz_y_0[j];

                tdx_xzz_z_0[j] = pa_x[j] * tdx_zz_z_0[j] + 0.5 * fl1_fx * ts_zz_z_0[j];

                tdy_xzz_z_0[j] = pa_x[j] * tdy_zz_z_0[j];

                tdz_xzz_z_0[j] = pa_x[j] * tdz_zz_z_0[j];

                tdx_yyy_x_0[j] = pa_y[j] * tdx_yy_x_0[j] + fl1_fx * tdx_y_x_0[j];

                tdy_yyy_x_0[j] = pa_y[j] * tdy_yy_x_0[j] + fl1_fx * tdy_y_x_0[j] + 0.5 * fl1_fx * ts_yy_x_0[j];

                tdz_yyy_x_0[j] = pa_y[j] * tdz_yy_x_0[j] + fl1_fx * tdz_y_x_0[j];

                tdx_yyy_y_0[j] = pa_y[j] * tdx_yy_y_0[j] + fl1_fx * tdx_y_y_0[j] + 0.5 * fl1_fx * tdx_yy_0_0[j];

                tdy_yyy_y_0[j] = pa_y[j] * tdy_yy_y_0[j] + fl1_fx * tdy_y_y_0[j] + 0.5 * fl1_fx * tdy_yy_0_0[j] + 0.5 * fl1_fx * ts_yy_y_0[j];

                tdz_yyy_y_0[j] = pa_y[j] * tdz_yy_y_0[j] + fl1_fx * tdz_y_y_0[j] + 0.5 * fl1_fx * tdz_yy_0_0[j];

                tdx_yyy_z_0[j] = pa_y[j] * tdx_yy_z_0[j] + fl1_fx * tdx_y_z_0[j];

                tdy_yyy_z_0[j] = pa_y[j] * tdy_yy_z_0[j] + fl1_fx * tdy_y_z_0[j] + 0.5 * fl1_fx * ts_yy_z_0[j];

                tdz_yyy_z_0[j] = pa_y[j] * tdz_yy_z_0[j] + fl1_fx * tdz_y_z_0[j];

                tdx_yyz_x_0[j] = pa_y[j] * tdx_yz_x_0[j] + 0.5 * fl1_fx * tdx_z_x_0[j];

                tdy_yyz_x_0[j] = pa_y[j] * tdy_yz_x_0[j] + 0.5 * fl1_fx * tdy_z_x_0[j] + 0.5 * fl1_fx * ts_yz_x_0[j];

                tdz_yyz_x_0[j] = pa_y[j] * tdz_yz_x_0[j] + 0.5 * fl1_fx * tdz_z_x_0[j];

                tdx_yyz_y_0[j] = pa_y[j] * tdx_yz_y_0[j] + 0.5 * fl1_fx * tdx_z_y_0[j] + 0.5 * fl1_fx * tdx_yz_0_0[j];

                tdy_yyz_y_0[j] = pa_y[j] * tdy_yz_y_0[j] + 0.5 * fl1_fx * tdy_z_y_0[j] + 0.5 * fl1_fx * tdy_yz_0_0[j] + 0.5 * fl1_fx * ts_yz_y_0[j];

                tdz_yyz_y_0[j] = pa_y[j] * tdz_yz_y_0[j] + 0.5 * fl1_fx * tdz_z_y_0[j] + 0.5 * fl1_fx * tdz_yz_0_0[j];

                tdx_yyz_z_0[j] = pa_y[j] * tdx_yz_z_0[j] + 0.5 * fl1_fx * tdx_z_z_0[j];

                tdy_yyz_z_0[j] = pa_y[j] * tdy_yz_z_0[j] + 0.5 * fl1_fx * tdy_z_z_0[j] + 0.5 * fl1_fx * ts_yz_z_0[j];

                tdz_yyz_z_0[j] = pa_y[j] * tdz_yz_z_0[j] + 0.5 * fl1_fx * tdz_z_z_0[j];

                tdx_yzz_x_0[j] = pa_y[j] * tdx_zz_x_0[j];

                tdy_yzz_x_0[j] = pa_y[j] * tdy_zz_x_0[j] + 0.5 * fl1_fx * ts_zz_x_0[j];

                tdz_yzz_x_0[j] = pa_y[j] * tdz_zz_x_0[j];

                tdx_yzz_y_0[j] = pa_y[j] * tdx_zz_y_0[j] + 0.5 * fl1_fx * tdx_zz_0_0[j];

                tdy_yzz_y_0[j] = pa_y[j] * tdy_zz_y_0[j] + 0.5 * fl1_fx * tdy_zz_0_0[j] + 0.5 * fl1_fx * ts_zz_y_0[j];

                tdz_yzz_y_0[j] = pa_y[j] * tdz_zz_y_0[j] + 0.5 * fl1_fx * tdz_zz_0_0[j];

                tdx_yzz_z_0[j] = pa_y[j] * tdx_zz_z_0[j];

                tdy_yzz_z_0[j] = pa_y[j] * tdy_zz_z_0[j] + 0.5 * fl1_fx * ts_zz_z_0[j];

                tdz_yzz_z_0[j] = pa_y[j] * tdz_zz_z_0[j];

                tdx_zzz_x_0[j] = pa_z[j] * tdx_zz_x_0[j] + fl1_fx * tdx_z_x_0[j];

                tdy_zzz_x_0[j] = pa_z[j] * tdy_zz_x_0[j] + fl1_fx * tdy_z_x_0[j];

                tdz_zzz_x_0[j] = pa_z[j] * tdz_zz_x_0[j] + fl1_fx * tdz_z_x_0[j] + 0.5 * fl1_fx * ts_zz_x_0[j];

                tdx_zzz_y_0[j] = pa_z[j] * tdx_zz_y_0[j] + fl1_fx * tdx_z_y_0[j];

                tdy_zzz_y_0[j] = pa_z[j] * tdy_zz_y_0[j] + fl1_fx * tdy_z_y_0[j];

                tdz_zzz_y_0[j] = pa_z[j] * tdz_zz_y_0[j] + fl1_fx * tdz_z_y_0[j] + 0.5 * fl1_fx * ts_zz_y_0[j];

                tdx_zzz_z_0[j] = pa_z[j] * tdx_zz_z_0[j] + fl1_fx * tdx_z_z_0[j] + 0.5 * fl1_fx * tdx_zz_0_0[j];

                tdy_zzz_z_0[j] = pa_z[j] * tdy_zz_z_0[j] + fl1_fx * tdy_z_z_0[j] + 0.5 * fl1_fx * tdy_zz_0_0[j];

                tdz_zzz_z_0[j] = pa_z[j] * tdz_zz_z_0[j] + fl1_fx * tdz_z_z_0[j] + 0.5 * fl1_fx * tdz_zz_0_0[j] + 0.5 * fl1_fx * ts_zz_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForPG(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForPG_0_45(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForPG_45_90(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        ediprecfunc::compElectricDipoleForPG_90_135(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 
    }

    void
    compElectricDipoleForPG_0_45(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_1_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

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

            auto ts_0_xxxx_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx); 

            auto ts_0_xxxy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 1); 

            auto ts_0_xxxz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 2); 

            auto ts_0_xxyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 3); 

            auto ts_0_xxyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 4); 

            auto ts_0_xxzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 5); 

            auto ts_0_xyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 6); 

            auto ts_0_xyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 7); 

            auto ts_0_xyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 8); 

            auto ts_0_xzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 9); 

            auto ts_0_yyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 10); 

            auto ts_0_yyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 11); 

            auto ts_0_yyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 12); 

            auto ts_0_yzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 13); 

            auto ts_0_zzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 14); 

            // set up pointers to integrals

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

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, tdx_0_xxx_0, tdx_0_xxxx_0, tdx_0_xxxy_0, tdx_0_xxxz_0, \
                                     tdx_0_xxy_0, tdx_0_xxyy_0, tdx_0_xxyz_0, tdx_0_xxz_0, tdx_0_xxzz_0, tdx_0_xyy_0, \
                                     tdx_0_xyyy_0, tdx_0_xyyz_0, tdx_0_xyz_0, tdx_0_xyzz_0, tdx_0_xzz_0, tdx_0_xzzz_0, \
                                     tdx_0_yyy_0, tdx_0_yyyy_0, tdx_0_yyyz_0, tdx_0_yyz_0, tdx_0_yyzz_0, tdx_0_yzz_0, \
                                     tdx_0_yzzz_0, tdx_0_zzz_0, tdx_0_zzzz_0, tdx_x_xxxx_0, tdx_x_xxxy_0, tdx_x_xxxz_0, \
                                     tdx_x_xxyy_0, tdx_x_xxyz_0, tdx_x_xxzz_0, tdx_x_xyyy_0, tdx_x_xyyz_0, tdx_x_xyzz_0, \
                                     tdx_x_xzzz_0, tdx_x_yyyy_0, tdx_x_yyyz_0, tdx_x_yyzz_0, tdx_x_yzzz_0, tdx_x_zzzz_0, \
                                     tdy_0_xxx_0, tdy_0_xxxx_0, tdy_0_xxxy_0, tdy_0_xxxz_0, tdy_0_xxy_0, tdy_0_xxyy_0, \
                                     tdy_0_xxyz_0, tdy_0_xxz_0, tdy_0_xxzz_0, tdy_0_xyy_0, tdy_0_xyyy_0, tdy_0_xyyz_0, \
                                     tdy_0_xyz_0, tdy_0_xyzz_0, tdy_0_xzz_0, tdy_0_xzzz_0, tdy_0_yyy_0, tdy_0_yyyy_0, \
                                     tdy_0_yyyz_0, tdy_0_yyz_0, tdy_0_yyzz_0, tdy_0_yzz_0, tdy_0_yzzz_0, tdy_0_zzz_0, \
                                     tdy_0_zzzz_0, tdy_x_xxxx_0, tdy_x_xxxy_0, tdy_x_xxxz_0, tdy_x_xxyy_0, tdy_x_xxyz_0, \
                                     tdy_x_xxzz_0, tdy_x_xyyy_0, tdy_x_xyyz_0, tdy_x_xyzz_0, tdy_x_xzzz_0, tdy_x_yyyy_0, \
                                     tdy_x_yyyz_0, tdy_x_yyzz_0, tdy_x_yzzz_0, tdy_x_zzzz_0, tdz_0_xxx_0, tdz_0_xxxx_0, \
                                     tdz_0_xxxy_0, tdz_0_xxxz_0, tdz_0_xxy_0, tdz_0_xxyy_0, tdz_0_xxyz_0, tdz_0_xxz_0, \
                                     tdz_0_xxzz_0, tdz_0_xyy_0, tdz_0_xyyy_0, tdz_0_xyyz_0, tdz_0_xyz_0, tdz_0_xyzz_0, \
                                     tdz_0_xzz_0, tdz_0_xzzz_0, tdz_0_yyy_0, tdz_0_yyyy_0, tdz_0_yyyz_0, tdz_0_yyz_0, \
                                     tdz_0_yyzz_0, tdz_0_yzz_0, tdz_0_yzzz_0, tdz_0_zzz_0, tdz_0_zzzz_0, tdz_x_xxxx_0, \
                                     tdz_x_xxxy_0, tdz_x_xxxz_0, tdz_x_xxyy_0, tdz_x_xxyz_0, tdz_x_xxzz_0, tdz_x_xyyy_0, \
                                     tdz_x_xyyz_0, tdz_x_xyzz_0, tdz_x_xzzz_0, tdz_x_yyyy_0, tdz_x_yyyz_0, tdz_x_yyzz_0, \
                                     tdz_x_yzzz_0, tdz_x_zzzz_0, ts_0_xxxx_0, ts_0_xxxy_0, ts_0_xxxz_0, ts_0_xxyy_0, \
                                     ts_0_xxyz_0, ts_0_xxzz_0, ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyzz_0, ts_0_xzzz_0, \
                                     ts_0_yyyy_0, ts_0_yyyz_0, ts_0_yyzz_0, ts_0_yzzz_0, ts_0_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_x_xxxx_0[j] = pa_x[j] * tdx_0_xxxx_0[j] + 2.0 * fl1_fx * tdx_0_xxx_0[j] + 0.5 * fl1_fx * ts_0_xxxx_0[j];

                tdy_x_xxxx_0[j] = pa_x[j] * tdy_0_xxxx_0[j] + 2.0 * fl1_fx * tdy_0_xxx_0[j];

                tdz_x_xxxx_0[j] = pa_x[j] * tdz_0_xxxx_0[j] + 2.0 * fl1_fx * tdz_0_xxx_0[j];

                tdx_x_xxxy_0[j] = pa_x[j] * tdx_0_xxxy_0[j] + 1.5 * fl1_fx * tdx_0_xxy_0[j] + 0.5 * fl1_fx * ts_0_xxxy_0[j];

                tdy_x_xxxy_0[j] = pa_x[j] * tdy_0_xxxy_0[j] + 1.5 * fl1_fx * tdy_0_xxy_0[j];

                tdz_x_xxxy_0[j] = pa_x[j] * tdz_0_xxxy_0[j] + 1.5 * fl1_fx * tdz_0_xxy_0[j];

                tdx_x_xxxz_0[j] = pa_x[j] * tdx_0_xxxz_0[j] + 1.5 * fl1_fx * tdx_0_xxz_0[j] + 0.5 * fl1_fx * ts_0_xxxz_0[j];

                tdy_x_xxxz_0[j] = pa_x[j] * tdy_0_xxxz_0[j] + 1.5 * fl1_fx * tdy_0_xxz_0[j];

                tdz_x_xxxz_0[j] = pa_x[j] * tdz_0_xxxz_0[j] + 1.5 * fl1_fx * tdz_0_xxz_0[j];

                tdx_x_xxyy_0[j] = pa_x[j] * tdx_0_xxyy_0[j] + fl1_fx * tdx_0_xyy_0[j] + 0.5 * fl1_fx * ts_0_xxyy_0[j];

                tdy_x_xxyy_0[j] = pa_x[j] * tdy_0_xxyy_0[j] + fl1_fx * tdy_0_xyy_0[j];

                tdz_x_xxyy_0[j] = pa_x[j] * tdz_0_xxyy_0[j] + fl1_fx * tdz_0_xyy_0[j];

                tdx_x_xxyz_0[j] = pa_x[j] * tdx_0_xxyz_0[j] + fl1_fx * tdx_0_xyz_0[j] + 0.5 * fl1_fx * ts_0_xxyz_0[j];

                tdy_x_xxyz_0[j] = pa_x[j] * tdy_0_xxyz_0[j] + fl1_fx * tdy_0_xyz_0[j];

                tdz_x_xxyz_0[j] = pa_x[j] * tdz_0_xxyz_0[j] + fl1_fx * tdz_0_xyz_0[j];

                tdx_x_xxzz_0[j] = pa_x[j] * tdx_0_xxzz_0[j] + fl1_fx * tdx_0_xzz_0[j] + 0.5 * fl1_fx * ts_0_xxzz_0[j];

                tdy_x_xxzz_0[j] = pa_x[j] * tdy_0_xxzz_0[j] + fl1_fx * tdy_0_xzz_0[j];

                tdz_x_xxzz_0[j] = pa_x[j] * tdz_0_xxzz_0[j] + fl1_fx * tdz_0_xzz_0[j];

                tdx_x_xyyy_0[j] = pa_x[j] * tdx_0_xyyy_0[j] + 0.5 * fl1_fx * tdx_0_yyy_0[j] + 0.5 * fl1_fx * ts_0_xyyy_0[j];

                tdy_x_xyyy_0[j] = pa_x[j] * tdy_0_xyyy_0[j] + 0.5 * fl1_fx * tdy_0_yyy_0[j];

                tdz_x_xyyy_0[j] = pa_x[j] * tdz_0_xyyy_0[j] + 0.5 * fl1_fx * tdz_0_yyy_0[j];

                tdx_x_xyyz_0[j] = pa_x[j] * tdx_0_xyyz_0[j] + 0.5 * fl1_fx * tdx_0_yyz_0[j] + 0.5 * fl1_fx * ts_0_xyyz_0[j];

                tdy_x_xyyz_0[j] = pa_x[j] * tdy_0_xyyz_0[j] + 0.5 * fl1_fx * tdy_0_yyz_0[j];

                tdz_x_xyyz_0[j] = pa_x[j] * tdz_0_xyyz_0[j] + 0.5 * fl1_fx * tdz_0_yyz_0[j];

                tdx_x_xyzz_0[j] = pa_x[j] * tdx_0_xyzz_0[j] + 0.5 * fl1_fx * tdx_0_yzz_0[j] + 0.5 * fl1_fx * ts_0_xyzz_0[j];

                tdy_x_xyzz_0[j] = pa_x[j] * tdy_0_xyzz_0[j] + 0.5 * fl1_fx * tdy_0_yzz_0[j];

                tdz_x_xyzz_0[j] = pa_x[j] * tdz_0_xyzz_0[j] + 0.5 * fl1_fx * tdz_0_yzz_0[j];

                tdx_x_xzzz_0[j] = pa_x[j] * tdx_0_xzzz_0[j] + 0.5 * fl1_fx * tdx_0_zzz_0[j] + 0.5 * fl1_fx * ts_0_xzzz_0[j];

                tdy_x_xzzz_0[j] = pa_x[j] * tdy_0_xzzz_0[j] + 0.5 * fl1_fx * tdy_0_zzz_0[j];

                tdz_x_xzzz_0[j] = pa_x[j] * tdz_0_xzzz_0[j] + 0.5 * fl1_fx * tdz_0_zzz_0[j];

                tdx_x_yyyy_0[j] = pa_x[j] * tdx_0_yyyy_0[j] + 0.5 * fl1_fx * ts_0_yyyy_0[j];

                tdy_x_yyyy_0[j] = pa_x[j] * tdy_0_yyyy_0[j];

                tdz_x_yyyy_0[j] = pa_x[j] * tdz_0_yyyy_0[j];

                tdx_x_yyyz_0[j] = pa_x[j] * tdx_0_yyyz_0[j] + 0.5 * fl1_fx * ts_0_yyyz_0[j];

                tdy_x_yyyz_0[j] = pa_x[j] * tdy_0_yyyz_0[j];

                tdz_x_yyyz_0[j] = pa_x[j] * tdz_0_yyyz_0[j];

                tdx_x_yyzz_0[j] = pa_x[j] * tdx_0_yyzz_0[j] + 0.5 * fl1_fx * ts_0_yyzz_0[j];

                tdy_x_yyzz_0[j] = pa_x[j] * tdy_0_yyzz_0[j];

                tdz_x_yyzz_0[j] = pa_x[j] * tdz_0_yyzz_0[j];

                tdx_x_yzzz_0[j] = pa_x[j] * tdx_0_yzzz_0[j] + 0.5 * fl1_fx * ts_0_yzzz_0[j];

                tdy_x_yzzz_0[j] = pa_x[j] * tdy_0_yzzz_0[j];

                tdz_x_yzzz_0[j] = pa_x[j] * tdz_0_yzzz_0[j];

                tdx_x_zzzz_0[j] = pa_x[j] * tdx_0_zzzz_0[j] + 0.5 * fl1_fx * ts_0_zzzz_0[j];

                tdy_x_zzzz_0[j] = pa_x[j] * tdy_0_zzzz_0[j];

                tdz_x_zzzz_0[j] = pa_x[j] * tdz_0_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForPG_45_90(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_1_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

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

            auto ts_0_xxxx_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx); 

            auto ts_0_xxxy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 1); 

            auto ts_0_xxxz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 2); 

            auto ts_0_xxyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 3); 

            auto ts_0_xxyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 4); 

            auto ts_0_xxzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 5); 

            auto ts_0_xyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 6); 

            auto ts_0_xyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 7); 

            auto ts_0_xyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 8); 

            auto ts_0_xzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 9); 

            auto ts_0_yyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 10); 

            auto ts_0_yyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 11); 

            auto ts_0_yyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 12); 

            auto ts_0_yzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 13); 

            auto ts_0_zzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 14); 

            // set up pointers to integrals

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

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_y, tdx_0_xxx_0, tdx_0_xxxx_0, tdx_0_xxxy_0, tdx_0_xxxz_0, \
                                     tdx_0_xxy_0, tdx_0_xxyy_0, tdx_0_xxyz_0, tdx_0_xxz_0, tdx_0_xxzz_0, tdx_0_xyy_0, \
                                     tdx_0_xyyy_0, tdx_0_xyyz_0, tdx_0_xyz_0, tdx_0_xyzz_0, tdx_0_xzz_0, tdx_0_xzzz_0, \
                                     tdx_0_yyy_0, tdx_0_yyyy_0, tdx_0_yyyz_0, tdx_0_yyz_0, tdx_0_yyzz_0, tdx_0_yzz_0, \
                                     tdx_0_yzzz_0, tdx_0_zzz_0, tdx_0_zzzz_0, tdx_y_xxxx_0, tdx_y_xxxy_0, tdx_y_xxxz_0, \
                                     tdx_y_xxyy_0, tdx_y_xxyz_0, tdx_y_xxzz_0, tdx_y_xyyy_0, tdx_y_xyyz_0, tdx_y_xyzz_0, \
                                     tdx_y_xzzz_0, tdx_y_yyyy_0, tdx_y_yyyz_0, tdx_y_yyzz_0, tdx_y_yzzz_0, tdx_y_zzzz_0, \
                                     tdy_0_xxx_0, tdy_0_xxxx_0, tdy_0_xxxy_0, tdy_0_xxxz_0, tdy_0_xxy_0, tdy_0_xxyy_0, \
                                     tdy_0_xxyz_0, tdy_0_xxz_0, tdy_0_xxzz_0, tdy_0_xyy_0, tdy_0_xyyy_0, tdy_0_xyyz_0, \
                                     tdy_0_xyz_0, tdy_0_xyzz_0, tdy_0_xzz_0, tdy_0_xzzz_0, tdy_0_yyy_0, tdy_0_yyyy_0, \
                                     tdy_0_yyyz_0, tdy_0_yyz_0, tdy_0_yyzz_0, tdy_0_yzz_0, tdy_0_yzzz_0, tdy_0_zzz_0, \
                                     tdy_0_zzzz_0, tdy_y_xxxx_0, tdy_y_xxxy_0, tdy_y_xxxz_0, tdy_y_xxyy_0, tdy_y_xxyz_0, \
                                     tdy_y_xxzz_0, tdy_y_xyyy_0, tdy_y_xyyz_0, tdy_y_xyzz_0, tdy_y_xzzz_0, tdy_y_yyyy_0, \
                                     tdy_y_yyyz_0, tdy_y_yyzz_0, tdy_y_yzzz_0, tdy_y_zzzz_0, tdz_0_xxx_0, tdz_0_xxxx_0, \
                                     tdz_0_xxxy_0, tdz_0_xxxz_0, tdz_0_xxy_0, tdz_0_xxyy_0, tdz_0_xxyz_0, tdz_0_xxz_0, \
                                     tdz_0_xxzz_0, tdz_0_xyy_0, tdz_0_xyyy_0, tdz_0_xyyz_0, tdz_0_xyz_0, tdz_0_xyzz_0, \
                                     tdz_0_xzz_0, tdz_0_xzzz_0, tdz_0_yyy_0, tdz_0_yyyy_0, tdz_0_yyyz_0, tdz_0_yyz_0, \
                                     tdz_0_yyzz_0, tdz_0_yzz_0, tdz_0_yzzz_0, tdz_0_zzz_0, tdz_0_zzzz_0, tdz_y_xxxx_0, \
                                     tdz_y_xxxy_0, tdz_y_xxxz_0, tdz_y_xxyy_0, tdz_y_xxyz_0, tdz_y_xxzz_0, tdz_y_xyyy_0, \
                                     tdz_y_xyyz_0, tdz_y_xyzz_0, tdz_y_xzzz_0, tdz_y_yyyy_0, tdz_y_yyyz_0, tdz_y_yyzz_0, \
                                     tdz_y_yzzz_0, tdz_y_zzzz_0, ts_0_xxxx_0, ts_0_xxxy_0, ts_0_xxxz_0, ts_0_xxyy_0, \
                                     ts_0_xxyz_0, ts_0_xxzz_0, ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyzz_0, ts_0_xzzz_0, \
                                     ts_0_yyyy_0, ts_0_yyyz_0, ts_0_yyzz_0, ts_0_yzzz_0, ts_0_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_y_xxxx_0[j] = pa_y[j] * tdx_0_xxxx_0[j];

                tdy_y_xxxx_0[j] = pa_y[j] * tdy_0_xxxx_0[j] + 0.5 * fl1_fx * ts_0_xxxx_0[j];

                tdz_y_xxxx_0[j] = pa_y[j] * tdz_0_xxxx_0[j];

                tdx_y_xxxy_0[j] = pa_y[j] * tdx_0_xxxy_0[j] + 0.5 * fl1_fx * tdx_0_xxx_0[j];

                tdy_y_xxxy_0[j] = pa_y[j] * tdy_0_xxxy_0[j] + 0.5 * fl1_fx * tdy_0_xxx_0[j] + 0.5 * fl1_fx * ts_0_xxxy_0[j];

                tdz_y_xxxy_0[j] = pa_y[j] * tdz_0_xxxy_0[j] + 0.5 * fl1_fx * tdz_0_xxx_0[j];

                tdx_y_xxxz_0[j] = pa_y[j] * tdx_0_xxxz_0[j];

                tdy_y_xxxz_0[j] = pa_y[j] * tdy_0_xxxz_0[j] + 0.5 * fl1_fx * ts_0_xxxz_0[j];

                tdz_y_xxxz_0[j] = pa_y[j] * tdz_0_xxxz_0[j];

                tdx_y_xxyy_0[j] = pa_y[j] * tdx_0_xxyy_0[j] + fl1_fx * tdx_0_xxy_0[j];

                tdy_y_xxyy_0[j] = pa_y[j] * tdy_0_xxyy_0[j] + fl1_fx * tdy_0_xxy_0[j] + 0.5 * fl1_fx * ts_0_xxyy_0[j];

                tdz_y_xxyy_0[j] = pa_y[j] * tdz_0_xxyy_0[j] + fl1_fx * tdz_0_xxy_0[j];

                tdx_y_xxyz_0[j] = pa_y[j] * tdx_0_xxyz_0[j] + 0.5 * fl1_fx * tdx_0_xxz_0[j];

                tdy_y_xxyz_0[j] = pa_y[j] * tdy_0_xxyz_0[j] + 0.5 * fl1_fx * tdy_0_xxz_0[j] + 0.5 * fl1_fx * ts_0_xxyz_0[j];

                tdz_y_xxyz_0[j] = pa_y[j] * tdz_0_xxyz_0[j] + 0.5 * fl1_fx * tdz_0_xxz_0[j];

                tdx_y_xxzz_0[j] = pa_y[j] * tdx_0_xxzz_0[j];

                tdy_y_xxzz_0[j] = pa_y[j] * tdy_0_xxzz_0[j] + 0.5 * fl1_fx * ts_0_xxzz_0[j];

                tdz_y_xxzz_0[j] = pa_y[j] * tdz_0_xxzz_0[j];

                tdx_y_xyyy_0[j] = pa_y[j] * tdx_0_xyyy_0[j] + 1.5 * fl1_fx * tdx_0_xyy_0[j];

                tdy_y_xyyy_0[j] = pa_y[j] * tdy_0_xyyy_0[j] + 1.5 * fl1_fx * tdy_0_xyy_0[j] + 0.5 * fl1_fx * ts_0_xyyy_0[j];

                tdz_y_xyyy_0[j] = pa_y[j] * tdz_0_xyyy_0[j] + 1.5 * fl1_fx * tdz_0_xyy_0[j];

                tdx_y_xyyz_0[j] = pa_y[j] * tdx_0_xyyz_0[j] + fl1_fx * tdx_0_xyz_0[j];

                tdy_y_xyyz_0[j] = pa_y[j] * tdy_0_xyyz_0[j] + fl1_fx * tdy_0_xyz_0[j] + 0.5 * fl1_fx * ts_0_xyyz_0[j];

                tdz_y_xyyz_0[j] = pa_y[j] * tdz_0_xyyz_0[j] + fl1_fx * tdz_0_xyz_0[j];

                tdx_y_xyzz_0[j] = pa_y[j] * tdx_0_xyzz_0[j] + 0.5 * fl1_fx * tdx_0_xzz_0[j];

                tdy_y_xyzz_0[j] = pa_y[j] * tdy_0_xyzz_0[j] + 0.5 * fl1_fx * tdy_0_xzz_0[j] + 0.5 * fl1_fx * ts_0_xyzz_0[j];

                tdz_y_xyzz_0[j] = pa_y[j] * tdz_0_xyzz_0[j] + 0.5 * fl1_fx * tdz_0_xzz_0[j];

                tdx_y_xzzz_0[j] = pa_y[j] * tdx_0_xzzz_0[j];

                tdy_y_xzzz_0[j] = pa_y[j] * tdy_0_xzzz_0[j] + 0.5 * fl1_fx * ts_0_xzzz_0[j];

                tdz_y_xzzz_0[j] = pa_y[j] * tdz_0_xzzz_0[j];

                tdx_y_yyyy_0[j] = pa_y[j] * tdx_0_yyyy_0[j] + 2.0 * fl1_fx * tdx_0_yyy_0[j];

                tdy_y_yyyy_0[j] = pa_y[j] * tdy_0_yyyy_0[j] + 2.0 * fl1_fx * tdy_0_yyy_0[j] + 0.5 * fl1_fx * ts_0_yyyy_0[j];

                tdz_y_yyyy_0[j] = pa_y[j] * tdz_0_yyyy_0[j] + 2.0 * fl1_fx * tdz_0_yyy_0[j];

                tdx_y_yyyz_0[j] = pa_y[j] * tdx_0_yyyz_0[j] + 1.5 * fl1_fx * tdx_0_yyz_0[j];

                tdy_y_yyyz_0[j] = pa_y[j] * tdy_0_yyyz_0[j] + 1.5 * fl1_fx * tdy_0_yyz_0[j] + 0.5 * fl1_fx * ts_0_yyyz_0[j];

                tdz_y_yyyz_0[j] = pa_y[j] * tdz_0_yyyz_0[j] + 1.5 * fl1_fx * tdz_0_yyz_0[j];

                tdx_y_yyzz_0[j] = pa_y[j] * tdx_0_yyzz_0[j] + fl1_fx * tdx_0_yzz_0[j];

                tdy_y_yyzz_0[j] = pa_y[j] * tdy_0_yyzz_0[j] + fl1_fx * tdy_0_yzz_0[j] + 0.5 * fl1_fx * ts_0_yyzz_0[j];

                tdz_y_yyzz_0[j] = pa_y[j] * tdz_0_yyzz_0[j] + fl1_fx * tdz_0_yzz_0[j];

                tdx_y_yzzz_0[j] = pa_y[j] * tdx_0_yzzz_0[j] + 0.5 * fl1_fx * tdx_0_zzz_0[j];

                tdy_y_yzzz_0[j] = pa_y[j] * tdy_0_yzzz_0[j] + 0.5 * fl1_fx * tdy_0_zzz_0[j] + 0.5 * fl1_fx * ts_0_yzzz_0[j];

                tdz_y_yzzz_0[j] = pa_y[j] * tdz_0_yzzz_0[j] + 0.5 * fl1_fx * tdz_0_zzz_0[j];

                tdx_y_zzzz_0[j] = pa_y[j] * tdx_0_zzzz_0[j];

                tdy_y_zzzz_0[j] = pa_y[j] * tdy_0_zzzz_0[j] + 0.5 * fl1_fx * ts_0_zzzz_0[j];

                tdz_y_zzzz_0[j] = pa_y[j] * tdz_0_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForPG_90_135(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_1_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto ts_0_xxxx_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx); 

            auto ts_0_xxxy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 1); 

            auto ts_0_xxxz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 2); 

            auto ts_0_xxyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 3); 

            auto ts_0_xxyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 4); 

            auto ts_0_xxzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 5); 

            auto ts_0_xyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 6); 

            auto ts_0_xyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 7); 

            auto ts_0_xyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 8); 

            auto ts_0_xzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 9); 

            auto ts_0_yyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 10); 

            auto ts_0_yyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 11); 

            auto ts_0_yyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 12); 

            auto ts_0_yzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 13); 

            auto ts_0_zzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 14); 

            // set up pointers to integrals

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

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fx, pa_z, tdx_0_xxx_0, tdx_0_xxxx_0, tdx_0_xxxy_0, tdx_0_xxxz_0, \
                                     tdx_0_xxy_0, tdx_0_xxyy_0, tdx_0_xxyz_0, tdx_0_xxz_0, tdx_0_xxzz_0, tdx_0_xyy_0, \
                                     tdx_0_xyyy_0, tdx_0_xyyz_0, tdx_0_xyz_0, tdx_0_xyzz_0, tdx_0_xzz_0, tdx_0_xzzz_0, \
                                     tdx_0_yyy_0, tdx_0_yyyy_0, tdx_0_yyyz_0, tdx_0_yyz_0, tdx_0_yyzz_0, tdx_0_yzz_0, \
                                     tdx_0_yzzz_0, tdx_0_zzz_0, tdx_0_zzzz_0, tdx_z_xxxx_0, tdx_z_xxxy_0, tdx_z_xxxz_0, \
                                     tdx_z_xxyy_0, tdx_z_xxyz_0, tdx_z_xxzz_0, tdx_z_xyyy_0, tdx_z_xyyz_0, tdx_z_xyzz_0, \
                                     tdx_z_xzzz_0, tdx_z_yyyy_0, tdx_z_yyyz_0, tdx_z_yyzz_0, tdx_z_yzzz_0, tdx_z_zzzz_0, \
                                     tdy_0_xxx_0, tdy_0_xxxx_0, tdy_0_xxxy_0, tdy_0_xxxz_0, tdy_0_xxy_0, tdy_0_xxyy_0, \
                                     tdy_0_xxyz_0, tdy_0_xxz_0, tdy_0_xxzz_0, tdy_0_xyy_0, tdy_0_xyyy_0, tdy_0_xyyz_0, \
                                     tdy_0_xyz_0, tdy_0_xyzz_0, tdy_0_xzz_0, tdy_0_xzzz_0, tdy_0_yyy_0, tdy_0_yyyy_0, \
                                     tdy_0_yyyz_0, tdy_0_yyz_0, tdy_0_yyzz_0, tdy_0_yzz_0, tdy_0_yzzz_0, tdy_0_zzz_0, \
                                     tdy_0_zzzz_0, tdy_z_xxxx_0, tdy_z_xxxy_0, tdy_z_xxxz_0, tdy_z_xxyy_0, tdy_z_xxyz_0, \
                                     tdy_z_xxzz_0, tdy_z_xyyy_0, tdy_z_xyyz_0, tdy_z_xyzz_0, tdy_z_xzzz_0, tdy_z_yyyy_0, \
                                     tdy_z_yyyz_0, tdy_z_yyzz_0, tdy_z_yzzz_0, tdy_z_zzzz_0, tdz_0_xxx_0, tdz_0_xxxx_0, \
                                     tdz_0_xxxy_0, tdz_0_xxxz_0, tdz_0_xxy_0, tdz_0_xxyy_0, tdz_0_xxyz_0, tdz_0_xxz_0, \
                                     tdz_0_xxzz_0, tdz_0_xyy_0, tdz_0_xyyy_0, tdz_0_xyyz_0, tdz_0_xyz_0, tdz_0_xyzz_0, \
                                     tdz_0_xzz_0, tdz_0_xzzz_0, tdz_0_yyy_0, tdz_0_yyyy_0, tdz_0_yyyz_0, tdz_0_yyz_0, \
                                     tdz_0_yyzz_0, tdz_0_yzz_0, tdz_0_yzzz_0, tdz_0_zzz_0, tdz_0_zzzz_0, tdz_z_xxxx_0, \
                                     tdz_z_xxxy_0, tdz_z_xxxz_0, tdz_z_xxyy_0, tdz_z_xxyz_0, tdz_z_xxzz_0, tdz_z_xyyy_0, \
                                     tdz_z_xyyz_0, tdz_z_xyzz_0, tdz_z_xzzz_0, tdz_z_yyyy_0, tdz_z_yyyz_0, tdz_z_yyzz_0, \
                                     tdz_z_yzzz_0, tdz_z_zzzz_0, ts_0_xxxx_0, ts_0_xxxy_0, ts_0_xxxz_0, ts_0_xxyy_0, \
                                     ts_0_xxyz_0, ts_0_xxzz_0, ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyzz_0, ts_0_xzzz_0, \
                                     ts_0_yyyy_0, ts_0_yyyz_0, ts_0_yyzz_0, ts_0_yzzz_0, ts_0_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_z_xxxx_0[j] = pa_z[j] * tdx_0_xxxx_0[j];

                tdy_z_xxxx_0[j] = pa_z[j] * tdy_0_xxxx_0[j];

                tdz_z_xxxx_0[j] = pa_z[j] * tdz_0_xxxx_0[j] + 0.5 * fl1_fx * ts_0_xxxx_0[j];

                tdx_z_xxxy_0[j] = pa_z[j] * tdx_0_xxxy_0[j];

                tdy_z_xxxy_0[j] = pa_z[j] * tdy_0_xxxy_0[j];

                tdz_z_xxxy_0[j] = pa_z[j] * tdz_0_xxxy_0[j] + 0.5 * fl1_fx * ts_0_xxxy_0[j];

                tdx_z_xxxz_0[j] = pa_z[j] * tdx_0_xxxz_0[j] + 0.5 * fl1_fx * tdx_0_xxx_0[j];

                tdy_z_xxxz_0[j] = pa_z[j] * tdy_0_xxxz_0[j] + 0.5 * fl1_fx * tdy_0_xxx_0[j];

                tdz_z_xxxz_0[j] = pa_z[j] * tdz_0_xxxz_0[j] + 0.5 * fl1_fx * tdz_0_xxx_0[j] + 0.5 * fl1_fx * ts_0_xxxz_0[j];

                tdx_z_xxyy_0[j] = pa_z[j] * tdx_0_xxyy_0[j];

                tdy_z_xxyy_0[j] = pa_z[j] * tdy_0_xxyy_0[j];

                tdz_z_xxyy_0[j] = pa_z[j] * tdz_0_xxyy_0[j] + 0.5 * fl1_fx * ts_0_xxyy_0[j];

                tdx_z_xxyz_0[j] = pa_z[j] * tdx_0_xxyz_0[j] + 0.5 * fl1_fx * tdx_0_xxy_0[j];

                tdy_z_xxyz_0[j] = pa_z[j] * tdy_0_xxyz_0[j] + 0.5 * fl1_fx * tdy_0_xxy_0[j];

                tdz_z_xxyz_0[j] = pa_z[j] * tdz_0_xxyz_0[j] + 0.5 * fl1_fx * tdz_0_xxy_0[j] + 0.5 * fl1_fx * ts_0_xxyz_0[j];

                tdx_z_xxzz_0[j] = pa_z[j] * tdx_0_xxzz_0[j] + fl1_fx * tdx_0_xxz_0[j];

                tdy_z_xxzz_0[j] = pa_z[j] * tdy_0_xxzz_0[j] + fl1_fx * tdy_0_xxz_0[j];

                tdz_z_xxzz_0[j] = pa_z[j] * tdz_0_xxzz_0[j] + fl1_fx * tdz_0_xxz_0[j] + 0.5 * fl1_fx * ts_0_xxzz_0[j];

                tdx_z_xyyy_0[j] = pa_z[j] * tdx_0_xyyy_0[j];

                tdy_z_xyyy_0[j] = pa_z[j] * tdy_0_xyyy_0[j];

                tdz_z_xyyy_0[j] = pa_z[j] * tdz_0_xyyy_0[j] + 0.5 * fl1_fx * ts_0_xyyy_0[j];

                tdx_z_xyyz_0[j] = pa_z[j] * tdx_0_xyyz_0[j] + 0.5 * fl1_fx * tdx_0_xyy_0[j];

                tdy_z_xyyz_0[j] = pa_z[j] * tdy_0_xyyz_0[j] + 0.5 * fl1_fx * tdy_0_xyy_0[j];

                tdz_z_xyyz_0[j] = pa_z[j] * tdz_0_xyyz_0[j] + 0.5 * fl1_fx * tdz_0_xyy_0[j] + 0.5 * fl1_fx * ts_0_xyyz_0[j];

                tdx_z_xyzz_0[j] = pa_z[j] * tdx_0_xyzz_0[j] + fl1_fx * tdx_0_xyz_0[j];

                tdy_z_xyzz_0[j] = pa_z[j] * tdy_0_xyzz_0[j] + fl1_fx * tdy_0_xyz_0[j];

                tdz_z_xyzz_0[j] = pa_z[j] * tdz_0_xyzz_0[j] + fl1_fx * tdz_0_xyz_0[j] + 0.5 * fl1_fx * ts_0_xyzz_0[j];

                tdx_z_xzzz_0[j] = pa_z[j] * tdx_0_xzzz_0[j] + 1.5 * fl1_fx * tdx_0_xzz_0[j];

                tdy_z_xzzz_0[j] = pa_z[j] * tdy_0_xzzz_0[j] + 1.5 * fl1_fx * tdy_0_xzz_0[j];

                tdz_z_xzzz_0[j] = pa_z[j] * tdz_0_xzzz_0[j] + 1.5 * fl1_fx * tdz_0_xzz_0[j] + 0.5 * fl1_fx * ts_0_xzzz_0[j];

                tdx_z_yyyy_0[j] = pa_z[j] * tdx_0_yyyy_0[j];

                tdy_z_yyyy_0[j] = pa_z[j] * tdy_0_yyyy_0[j];

                tdz_z_yyyy_0[j] = pa_z[j] * tdz_0_yyyy_0[j] + 0.5 * fl1_fx * ts_0_yyyy_0[j];

                tdx_z_yyyz_0[j] = pa_z[j] * tdx_0_yyyz_0[j] + 0.5 * fl1_fx * tdx_0_yyy_0[j];

                tdy_z_yyyz_0[j] = pa_z[j] * tdy_0_yyyz_0[j] + 0.5 * fl1_fx * tdy_0_yyy_0[j];

                tdz_z_yyyz_0[j] = pa_z[j] * tdz_0_yyyz_0[j] + 0.5 * fl1_fx * tdz_0_yyy_0[j] + 0.5 * fl1_fx * ts_0_yyyz_0[j];

                tdx_z_yyzz_0[j] = pa_z[j] * tdx_0_yyzz_0[j] + fl1_fx * tdx_0_yyz_0[j];

                tdy_z_yyzz_0[j] = pa_z[j] * tdy_0_yyzz_0[j] + fl1_fx * tdy_0_yyz_0[j];

                tdz_z_yyzz_0[j] = pa_z[j] * tdz_0_yyzz_0[j] + fl1_fx * tdz_0_yyz_0[j] + 0.5 * fl1_fx * ts_0_yyzz_0[j];

                tdx_z_yzzz_0[j] = pa_z[j] * tdx_0_yzzz_0[j] + 1.5 * fl1_fx * tdx_0_yzz_0[j];

                tdy_z_yzzz_0[j] = pa_z[j] * tdy_0_yzzz_0[j] + 1.5 * fl1_fx * tdy_0_yzz_0[j];

                tdz_z_yzzz_0[j] = pa_z[j] * tdz_0_yzzz_0[j] + 1.5 * fl1_fx * tdz_0_yzz_0[j] + 0.5 * fl1_fx * ts_0_yzzz_0[j];

                tdx_z_zzzz_0[j] = pa_z[j] * tdx_0_zzzz_0[j] + 2.0 * fl1_fx * tdx_0_zzz_0[j];

                tdy_z_zzzz_0[j] = pa_z[j] * tdy_0_zzzz_0[j] + 2.0 * fl1_fx * tdy_0_zzz_0[j];

                tdz_z_zzzz_0[j] = pa_z[j] * tdz_0_zzzz_0[j] + 2.0 * fl1_fx * tdz_0_zzz_0[j] + 0.5 * fl1_fx * ts_0_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGP(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        ediprecfunc::compElectricDipoleForGP_0_45(primBuffer,
                                                  recursionMap,
                                                  osFactors,
                                                  paDistances, 
                                                  braGtoBlock,
                                                  ketGtoBlock,
                                                  iContrGto); 

        ediprecfunc::compElectricDipoleForGP_45_90(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        ediprecfunc::compElectricDipoleForGP_90_135(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 
    }

    void
    compElectricDipoleForGP_0_45(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

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

            auto tdx_xxx_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx); 

            auto tdy_xxx_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx); 

            auto tdz_xxx_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx); 

            auto tdx_xxy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 1); 

            auto tdy_xxy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 1); 

            auto tdz_xxy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 1); 

            auto tdx_xxz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 2); 

            auto tdy_xxz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 2); 

            auto tdz_xxz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 2); 

            auto tdx_xyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 3); 

            auto tdy_xyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 3); 

            auto tdz_xyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 3); 

            auto tdx_xyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 4); 

            auto tdy_xyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 4); 

            auto tdz_xyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 4); 

            auto ts_xxx_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx); 

            auto ts_xxx_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 1); 

            auto ts_xxx_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 2); 

            auto ts_xxy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 3); 

            auto ts_xxy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 4); 

            auto ts_xxy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 5); 

            auto ts_xxz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 6); 

            auto ts_xxz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 7); 

            auto ts_xxz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 8); 

            auto ts_xyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 9); 

            auto ts_xyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 10); 

            auto ts_xyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 11); 

            auto ts_xyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 12); 

            auto ts_xyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 13); 

            auto ts_xyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 14); 

            // set up pointers to integrals

            auto tdx_xxxx_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx); 

            auto tdy_xxxx_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx); 

            auto tdz_xxxx_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx); 

            auto tdx_xxxx_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 1); 

            auto tdy_xxxx_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 1); 

            auto tdz_xxxx_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 1); 

            auto tdx_xxxx_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 2); 

            auto tdy_xxxx_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 2); 

            auto tdz_xxxx_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 2); 

            auto tdx_xxxy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 3); 

            auto tdy_xxxy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 3); 

            auto tdz_xxxy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 3); 

            auto tdx_xxxy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 4); 

            auto tdy_xxxy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 4); 

            auto tdz_xxxy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 4); 

            auto tdx_xxxy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 5); 

            auto tdy_xxxy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 5); 

            auto tdz_xxxy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 5); 

            auto tdx_xxxz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 6); 

            auto tdy_xxxz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 6); 

            auto tdz_xxxz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 6); 

            auto tdx_xxxz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 7); 

            auto tdy_xxxz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 7); 

            auto tdz_xxxz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 7); 

            auto tdx_xxxz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 8); 

            auto tdy_xxxz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 8); 

            auto tdz_xxxz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 8); 

            auto tdx_xxyy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 9); 

            auto tdy_xxyy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 9); 

            auto tdz_xxyy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 9); 

            auto tdx_xxyy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 10); 

            auto tdy_xxyy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 10); 

            auto tdz_xxyy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 10); 

            auto tdx_xxyy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 11); 

            auto tdy_xxyy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 11); 

            auto tdz_xxyy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 11); 

            auto tdx_xxyz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 12); 

            auto tdy_xxyz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 12); 

            auto tdz_xxyz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 12); 

            auto tdx_xxyz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 13); 

            auto tdy_xxyz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 13); 

            auto tdz_xxyz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 13); 

            auto tdx_xxyz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 14); 

            auto tdy_xxyz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 14); 

            auto tdz_xxyz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 14); 

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, tdx_xx_x_0, tdx_xx_y_0, tdx_xx_z_0, tdx_xxx_0_0, tdx_xxx_x_0, \
                                     tdx_xxx_y_0, tdx_xxx_z_0, tdx_xxxx_x_0, tdx_xxxx_y_0, tdx_xxxx_z_0, tdx_xxxy_x_0, \
                                     tdx_xxxy_y_0, tdx_xxxy_z_0, tdx_xxxz_x_0, tdx_xxxz_y_0, tdx_xxxz_z_0, tdx_xxy_0_0, \
                                     tdx_xxy_x_0, tdx_xxy_y_0, tdx_xxy_z_0, tdx_xxyy_x_0, tdx_xxyy_y_0, tdx_xxyy_z_0, \
                                     tdx_xxyz_x_0, tdx_xxyz_y_0, tdx_xxyz_z_0, tdx_xxz_0_0, tdx_xxz_x_0, tdx_xxz_y_0, \
                                     tdx_xxz_z_0, tdx_xy_x_0, tdx_xy_y_0, tdx_xy_z_0, tdx_xyy_0_0, tdx_xyy_x_0, \
                                     tdx_xyy_y_0, tdx_xyy_z_0, tdx_xyz_0_0, tdx_xyz_x_0, tdx_xyz_y_0, tdx_xyz_z_0, \
                                     tdx_xz_x_0, tdx_xz_y_0, tdx_xz_z_0, tdx_yy_x_0, tdx_yy_y_0, tdx_yy_z_0, tdx_yz_x_0, \
                                     tdx_yz_y_0, tdx_yz_z_0, tdy_xx_x_0, tdy_xx_y_0, tdy_xx_z_0, tdy_xxx_0_0, \
                                     tdy_xxx_x_0, tdy_xxx_y_0, tdy_xxx_z_0, tdy_xxxx_x_0, tdy_xxxx_y_0, tdy_xxxx_z_0, \
                                     tdy_xxxy_x_0, tdy_xxxy_y_0, tdy_xxxy_z_0, tdy_xxxz_x_0, tdy_xxxz_y_0, tdy_xxxz_z_0, \
                                     tdy_xxy_0_0, tdy_xxy_x_0, tdy_xxy_y_0, tdy_xxy_z_0, tdy_xxyy_x_0, tdy_xxyy_y_0, \
                                     tdy_xxyy_z_0, tdy_xxyz_x_0, tdy_xxyz_y_0, tdy_xxyz_z_0, tdy_xxz_0_0, tdy_xxz_x_0, \
                                     tdy_xxz_y_0, tdy_xxz_z_0, tdy_xy_x_0, tdy_xy_y_0, tdy_xy_z_0, tdy_xyy_0_0, \
                                     tdy_xyy_x_0, tdy_xyy_y_0, tdy_xyy_z_0, tdy_xyz_0_0, tdy_xyz_x_0, tdy_xyz_y_0, \
                                     tdy_xyz_z_0, tdy_xz_x_0, tdy_xz_y_0, tdy_xz_z_0, tdy_yy_x_0, tdy_yy_y_0, tdy_yy_z_0, \
                                     tdy_yz_x_0, tdy_yz_y_0, tdy_yz_z_0, tdz_xx_x_0, tdz_xx_y_0, tdz_xx_z_0, \
                                     tdz_xxx_0_0, tdz_xxx_x_0, tdz_xxx_y_0, tdz_xxx_z_0, tdz_xxxx_x_0, tdz_xxxx_y_0, \
                                     tdz_xxxx_z_0, tdz_xxxy_x_0, tdz_xxxy_y_0, tdz_xxxy_z_0, tdz_xxxz_x_0, tdz_xxxz_y_0, \
                                     tdz_xxxz_z_0, tdz_xxy_0_0, tdz_xxy_x_0, tdz_xxy_y_0, tdz_xxy_z_0, tdz_xxyy_x_0, \
                                     tdz_xxyy_y_0, tdz_xxyy_z_0, tdz_xxyz_x_0, tdz_xxyz_y_0, tdz_xxyz_z_0, tdz_xxz_0_0, \
                                     tdz_xxz_x_0, tdz_xxz_y_0, tdz_xxz_z_0, tdz_xy_x_0, tdz_xy_y_0, tdz_xy_z_0, \
                                     tdz_xyy_0_0, tdz_xyy_x_0, tdz_xyy_y_0, tdz_xyy_z_0, tdz_xyz_0_0, tdz_xyz_x_0, \
                                     tdz_xyz_y_0, tdz_xyz_z_0, tdz_xz_x_0, tdz_xz_y_0, tdz_xz_z_0, tdz_yy_x_0, \
                                     tdz_yy_y_0, tdz_yy_z_0, tdz_yz_x_0, tdz_yz_y_0, tdz_yz_z_0, ts_xxx_x_0, ts_xxx_y_0, \
                                     ts_xxx_z_0, ts_xxy_x_0, ts_xxy_y_0, ts_xxy_z_0, ts_xxz_x_0, ts_xxz_y_0, ts_xxz_z_0, \
                                     ts_xyy_x_0, ts_xyy_y_0, ts_xyy_z_0, ts_xyz_x_0, ts_xyz_y_0, ts_xyz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxxx_x_0[j] = pa_x[j] * tdx_xxx_x_0[j] + 1.5 * fl1_fx * tdx_xx_x_0[j] + 0.5 * fl1_fx * tdx_xxx_0_0[j] + 0.5 * fl1_fx * ts_xxx_x_0[j];

                tdy_xxxx_x_0[j] = pa_x[j] * tdy_xxx_x_0[j] + 1.5 * fl1_fx * tdy_xx_x_0[j] + 0.5 * fl1_fx * tdy_xxx_0_0[j];

                tdz_xxxx_x_0[j] = pa_x[j] * tdz_xxx_x_0[j] + 1.5 * fl1_fx * tdz_xx_x_0[j] + 0.5 * fl1_fx * tdz_xxx_0_0[j];

                tdx_xxxx_y_0[j] = pa_x[j] * tdx_xxx_y_0[j] + 1.5 * fl1_fx * tdx_xx_y_0[j] + 0.5 * fl1_fx * ts_xxx_y_0[j];

                tdy_xxxx_y_0[j] = pa_x[j] * tdy_xxx_y_0[j] + 1.5 * fl1_fx * tdy_xx_y_0[j];

                tdz_xxxx_y_0[j] = pa_x[j] * tdz_xxx_y_0[j] + 1.5 * fl1_fx * tdz_xx_y_0[j];

                tdx_xxxx_z_0[j] = pa_x[j] * tdx_xxx_z_0[j] + 1.5 * fl1_fx * tdx_xx_z_0[j] + 0.5 * fl1_fx * ts_xxx_z_0[j];

                tdy_xxxx_z_0[j] = pa_x[j] * tdy_xxx_z_0[j] + 1.5 * fl1_fx * tdy_xx_z_0[j];

                tdz_xxxx_z_0[j] = pa_x[j] * tdz_xxx_z_0[j] + 1.5 * fl1_fx * tdz_xx_z_0[j];

                tdx_xxxy_x_0[j] = pa_x[j] * tdx_xxy_x_0[j] + fl1_fx * tdx_xy_x_0[j] + 0.5 * fl1_fx * tdx_xxy_0_0[j] + 0.5 * fl1_fx * ts_xxy_x_0[j];

                tdy_xxxy_x_0[j] = pa_x[j] * tdy_xxy_x_0[j] + fl1_fx * tdy_xy_x_0[j] + 0.5 * fl1_fx * tdy_xxy_0_0[j];

                tdz_xxxy_x_0[j] = pa_x[j] * tdz_xxy_x_0[j] + fl1_fx * tdz_xy_x_0[j] + 0.5 * fl1_fx * tdz_xxy_0_0[j];

                tdx_xxxy_y_0[j] = pa_x[j] * tdx_xxy_y_0[j] + fl1_fx * tdx_xy_y_0[j] + 0.5 * fl1_fx * ts_xxy_y_0[j];

                tdy_xxxy_y_0[j] = pa_x[j] * tdy_xxy_y_0[j] + fl1_fx * tdy_xy_y_0[j];

                tdz_xxxy_y_0[j] = pa_x[j] * tdz_xxy_y_0[j] + fl1_fx * tdz_xy_y_0[j];

                tdx_xxxy_z_0[j] = pa_x[j] * tdx_xxy_z_0[j] + fl1_fx * tdx_xy_z_0[j] + 0.5 * fl1_fx * ts_xxy_z_0[j];

                tdy_xxxy_z_0[j] = pa_x[j] * tdy_xxy_z_0[j] + fl1_fx * tdy_xy_z_0[j];

                tdz_xxxy_z_0[j] = pa_x[j] * tdz_xxy_z_0[j] + fl1_fx * tdz_xy_z_0[j];

                tdx_xxxz_x_0[j] = pa_x[j] * tdx_xxz_x_0[j] + fl1_fx * tdx_xz_x_0[j] + 0.5 * fl1_fx * tdx_xxz_0_0[j] + 0.5 * fl1_fx * ts_xxz_x_0[j];

                tdy_xxxz_x_0[j] = pa_x[j] * tdy_xxz_x_0[j] + fl1_fx * tdy_xz_x_0[j] + 0.5 * fl1_fx * tdy_xxz_0_0[j];

                tdz_xxxz_x_0[j] = pa_x[j] * tdz_xxz_x_0[j] + fl1_fx * tdz_xz_x_0[j] + 0.5 * fl1_fx * tdz_xxz_0_0[j];

                tdx_xxxz_y_0[j] = pa_x[j] * tdx_xxz_y_0[j] + fl1_fx * tdx_xz_y_0[j] + 0.5 * fl1_fx * ts_xxz_y_0[j];

                tdy_xxxz_y_0[j] = pa_x[j] * tdy_xxz_y_0[j] + fl1_fx * tdy_xz_y_0[j];

                tdz_xxxz_y_0[j] = pa_x[j] * tdz_xxz_y_0[j] + fl1_fx * tdz_xz_y_0[j];

                tdx_xxxz_z_0[j] = pa_x[j] * tdx_xxz_z_0[j] + fl1_fx * tdx_xz_z_0[j] + 0.5 * fl1_fx * ts_xxz_z_0[j];

                tdy_xxxz_z_0[j] = pa_x[j] * tdy_xxz_z_0[j] + fl1_fx * tdy_xz_z_0[j];

                tdz_xxxz_z_0[j] = pa_x[j] * tdz_xxz_z_0[j] + fl1_fx * tdz_xz_z_0[j];

                tdx_xxyy_x_0[j] = pa_x[j] * tdx_xyy_x_0[j] + 0.5 * fl1_fx * tdx_yy_x_0[j] + 0.5 * fl1_fx * tdx_xyy_0_0[j] + 0.5 * fl1_fx * ts_xyy_x_0[j];

                tdy_xxyy_x_0[j] = pa_x[j] * tdy_xyy_x_0[j] + 0.5 * fl1_fx * tdy_yy_x_0[j] + 0.5 * fl1_fx * tdy_xyy_0_0[j];

                tdz_xxyy_x_0[j] = pa_x[j] * tdz_xyy_x_0[j] + 0.5 * fl1_fx * tdz_yy_x_0[j] + 0.5 * fl1_fx * tdz_xyy_0_0[j];

                tdx_xxyy_y_0[j] = pa_x[j] * tdx_xyy_y_0[j] + 0.5 * fl1_fx * tdx_yy_y_0[j] + 0.5 * fl1_fx * ts_xyy_y_0[j];

                tdy_xxyy_y_0[j] = pa_x[j] * tdy_xyy_y_0[j] + 0.5 * fl1_fx * tdy_yy_y_0[j];

                tdz_xxyy_y_0[j] = pa_x[j] * tdz_xyy_y_0[j] + 0.5 * fl1_fx * tdz_yy_y_0[j];

                tdx_xxyy_z_0[j] = pa_x[j] * tdx_xyy_z_0[j] + 0.5 * fl1_fx * tdx_yy_z_0[j] + 0.5 * fl1_fx * ts_xyy_z_0[j];

                tdy_xxyy_z_0[j] = pa_x[j] * tdy_xyy_z_0[j] + 0.5 * fl1_fx * tdy_yy_z_0[j];

                tdz_xxyy_z_0[j] = pa_x[j] * tdz_xyy_z_0[j] + 0.5 * fl1_fx * tdz_yy_z_0[j];

                tdx_xxyz_x_0[j] = pa_x[j] * tdx_xyz_x_0[j] + 0.5 * fl1_fx * tdx_yz_x_0[j] + 0.5 * fl1_fx * tdx_xyz_0_0[j] + 0.5 * fl1_fx * ts_xyz_x_0[j];

                tdy_xxyz_x_0[j] = pa_x[j] * tdy_xyz_x_0[j] + 0.5 * fl1_fx * tdy_yz_x_0[j] + 0.5 * fl1_fx * tdy_xyz_0_0[j];

                tdz_xxyz_x_0[j] = pa_x[j] * tdz_xyz_x_0[j] + 0.5 * fl1_fx * tdz_yz_x_0[j] + 0.5 * fl1_fx * tdz_xyz_0_0[j];

                tdx_xxyz_y_0[j] = pa_x[j] * tdx_xyz_y_0[j] + 0.5 * fl1_fx * tdx_yz_y_0[j] + 0.5 * fl1_fx * ts_xyz_y_0[j];

                tdy_xxyz_y_0[j] = pa_x[j] * tdy_xyz_y_0[j] + 0.5 * fl1_fx * tdy_yz_y_0[j];

                tdz_xxyz_y_0[j] = pa_x[j] * tdz_xyz_y_0[j] + 0.5 * fl1_fx * tdz_yz_y_0[j];

                tdx_xxyz_z_0[j] = pa_x[j] * tdx_xyz_z_0[j] + 0.5 * fl1_fx * tdx_yz_z_0[j] + 0.5 * fl1_fx * ts_xyz_z_0[j];

                tdy_xxyz_z_0[j] = pa_x[j] * tdy_xyz_z_0[j] + 0.5 * fl1_fx * tdy_yz_z_0[j];

                tdz_xxyz_z_0[j] = pa_x[j] * tdz_xyz_z_0[j] + 0.5 * fl1_fx * tdz_yz_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGP_45_90(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

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

            auto tdx_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 15); 

            auto tdy_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 15); 

            auto tdz_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 15); 

            auto tdx_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 16); 

            auto tdy_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 16); 

            auto tdz_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 16); 

            auto tdx_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 17); 

            auto tdy_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 17); 

            auto tdz_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 17); 

            auto tdx_xzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 5); 

            auto tdy_xzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 5); 

            auto tdz_xzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 5); 

            auto tdx_yyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 6); 

            auto tdy_yyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 6); 

            auto tdz_yyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 6); 

            auto tdx_yyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 7); 

            auto tdy_yyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 7); 

            auto tdz_yyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 7); 

            auto tdx_yzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 8); 

            auto tdy_yzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 8); 

            auto tdz_yzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 8); 

            auto tdx_zzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 9); 

            auto tdy_zzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 9); 

            auto tdz_zzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 9); 

            auto ts_xzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 15); 

            auto ts_xzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 16); 

            auto ts_xzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 17); 

            auto ts_yyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 18); 

            auto ts_yyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 19); 

            auto ts_yyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 20); 

            auto ts_yyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 21); 

            auto ts_yyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 22); 

            auto ts_yyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 23); 

            auto ts_yzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 24); 

            auto ts_yzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 25); 

            auto ts_yzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 26); 

            auto ts_zzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 27); 

            auto ts_zzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 28); 

            auto ts_zzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 29); 

            // set up pointers to integrals

            auto tdx_xxzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 15); 

            auto tdy_xxzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 15); 

            auto tdz_xxzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 15); 

            auto tdx_xxzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 16); 

            auto tdy_xxzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 16); 

            auto tdz_xxzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 16); 

            auto tdx_xxzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 17); 

            auto tdy_xxzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 17); 

            auto tdz_xxzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 17); 

            auto tdx_xyyy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 18); 

            auto tdy_xyyy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 18); 

            auto tdz_xyyy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 18); 

            auto tdx_xyyy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 19); 

            auto tdy_xyyy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 19); 

            auto tdz_xyyy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 19); 

            auto tdx_xyyy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 20); 

            auto tdy_xyyy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 20); 

            auto tdz_xyyy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 20); 

            auto tdx_xyyz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 21); 

            auto tdy_xyyz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 21); 

            auto tdz_xyyz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 21); 

            auto tdx_xyyz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 22); 

            auto tdy_xyyz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 22); 

            auto tdz_xyyz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 22); 

            auto tdx_xyyz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 23); 

            auto tdy_xyyz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 23); 

            auto tdz_xyyz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 23); 

            auto tdx_xyzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 24); 

            auto tdy_xyzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 24); 

            auto tdz_xyzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 24); 

            auto tdx_xyzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 25); 

            auto tdy_xyzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 25); 

            auto tdz_xyzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 25); 

            auto tdx_xyzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 26); 

            auto tdy_xyzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 26); 

            auto tdz_xyzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 26); 

            auto tdx_xzzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 27); 

            auto tdy_xzzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 27); 

            auto tdz_xzzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 27); 

            auto tdx_xzzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 28); 

            auto tdy_xzzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 28); 

            auto tdz_xzzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 28); 

            auto tdx_xzzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 29); 

            auto tdy_xzzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 29); 

            auto tdz_xzzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 29); 

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, tdx_xxzz_x_0, tdx_xxzz_y_0, tdx_xxzz_z_0, tdx_xyyy_x_0, \
                                     tdx_xyyy_y_0, tdx_xyyy_z_0, tdx_xyyz_x_0, tdx_xyyz_y_0, tdx_xyyz_z_0, tdx_xyzz_x_0, \
                                     tdx_xyzz_y_0, tdx_xyzz_z_0, tdx_xzz_0_0, tdx_xzz_x_0, tdx_xzz_y_0, tdx_xzz_z_0, \
                                     tdx_xzzz_x_0, tdx_xzzz_y_0, tdx_xzzz_z_0, tdx_yyy_0_0, tdx_yyy_x_0, tdx_yyy_y_0, \
                                     tdx_yyy_z_0, tdx_yyz_0_0, tdx_yyz_x_0, tdx_yyz_y_0, tdx_yyz_z_0, tdx_yzz_0_0, \
                                     tdx_yzz_x_0, tdx_yzz_y_0, tdx_yzz_z_0, tdx_zz_x_0, tdx_zz_y_0, tdx_zz_z_0, \
                                     tdx_zzz_0_0, tdx_zzz_x_0, tdx_zzz_y_0, tdx_zzz_z_0, tdy_xxzz_x_0, tdy_xxzz_y_0, \
                                     tdy_xxzz_z_0, tdy_xyyy_x_0, tdy_xyyy_y_0, tdy_xyyy_z_0, tdy_xyyz_x_0, tdy_xyyz_y_0, \
                                     tdy_xyyz_z_0, tdy_xyzz_x_0, tdy_xyzz_y_0, tdy_xyzz_z_0, tdy_xzz_0_0, tdy_xzz_x_0, \
                                     tdy_xzz_y_0, tdy_xzz_z_0, tdy_xzzz_x_0, tdy_xzzz_y_0, tdy_xzzz_z_0, tdy_yyy_0_0, \
                                     tdy_yyy_x_0, tdy_yyy_y_0, tdy_yyy_z_0, tdy_yyz_0_0, tdy_yyz_x_0, tdy_yyz_y_0, \
                                     tdy_yyz_z_0, tdy_yzz_0_0, tdy_yzz_x_0, tdy_yzz_y_0, tdy_yzz_z_0, tdy_zz_x_0, \
                                     tdy_zz_y_0, tdy_zz_z_0, tdy_zzz_0_0, tdy_zzz_x_0, tdy_zzz_y_0, tdy_zzz_z_0, \
                                     tdz_xxzz_x_0, tdz_xxzz_y_0, tdz_xxzz_z_0, tdz_xyyy_x_0, tdz_xyyy_y_0, tdz_xyyy_z_0, \
                                     tdz_xyyz_x_0, tdz_xyyz_y_0, tdz_xyyz_z_0, tdz_xyzz_x_0, tdz_xyzz_y_0, tdz_xyzz_z_0, \
                                     tdz_xzz_0_0, tdz_xzz_x_0, tdz_xzz_y_0, tdz_xzz_z_0, tdz_xzzz_x_0, tdz_xzzz_y_0, \
                                     tdz_xzzz_z_0, tdz_yyy_0_0, tdz_yyy_x_0, tdz_yyy_y_0, tdz_yyy_z_0, tdz_yyz_0_0, \
                                     tdz_yyz_x_0, tdz_yyz_y_0, tdz_yyz_z_0, tdz_yzz_0_0, tdz_yzz_x_0, tdz_yzz_y_0, \
                                     tdz_yzz_z_0, tdz_zz_x_0, tdz_zz_y_0, tdz_zz_z_0, tdz_zzz_0_0, tdz_zzz_x_0, \
                                     tdz_zzz_y_0, tdz_zzz_z_0, ts_xzz_x_0, ts_xzz_y_0, ts_xzz_z_0, ts_yyy_x_0, \
                                     ts_yyy_y_0, ts_yyy_z_0, ts_yyz_x_0, ts_yyz_y_0, ts_yyz_z_0, ts_yzz_x_0, ts_yzz_y_0, \
                                     ts_yzz_z_0, ts_zzz_x_0, ts_zzz_y_0, ts_zzz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxzz_x_0[j] = pa_x[j] * tdx_xzz_x_0[j] + 0.5 * fl1_fx * tdx_zz_x_0[j] + 0.5 * fl1_fx * tdx_xzz_0_0[j] + 0.5 * fl1_fx * ts_xzz_x_0[j];

                tdy_xxzz_x_0[j] = pa_x[j] * tdy_xzz_x_0[j] + 0.5 * fl1_fx * tdy_zz_x_0[j] + 0.5 * fl1_fx * tdy_xzz_0_0[j];

                tdz_xxzz_x_0[j] = pa_x[j] * tdz_xzz_x_0[j] + 0.5 * fl1_fx * tdz_zz_x_0[j] + 0.5 * fl1_fx * tdz_xzz_0_0[j];

                tdx_xxzz_y_0[j] = pa_x[j] * tdx_xzz_y_0[j] + 0.5 * fl1_fx * tdx_zz_y_0[j] + 0.5 * fl1_fx * ts_xzz_y_0[j];

                tdy_xxzz_y_0[j] = pa_x[j] * tdy_xzz_y_0[j] + 0.5 * fl1_fx * tdy_zz_y_0[j];

                tdz_xxzz_y_0[j] = pa_x[j] * tdz_xzz_y_0[j] + 0.5 * fl1_fx * tdz_zz_y_0[j];

                tdx_xxzz_z_0[j] = pa_x[j] * tdx_xzz_z_0[j] + 0.5 * fl1_fx * tdx_zz_z_0[j] + 0.5 * fl1_fx * ts_xzz_z_0[j];

                tdy_xxzz_z_0[j] = pa_x[j] * tdy_xzz_z_0[j] + 0.5 * fl1_fx * tdy_zz_z_0[j];

                tdz_xxzz_z_0[j] = pa_x[j] * tdz_xzz_z_0[j] + 0.5 * fl1_fx * tdz_zz_z_0[j];

                tdx_xyyy_x_0[j] = pa_x[j] * tdx_yyy_x_0[j] + 0.5 * fl1_fx * tdx_yyy_0_0[j] + 0.5 * fl1_fx * ts_yyy_x_0[j];

                tdy_xyyy_x_0[j] = pa_x[j] * tdy_yyy_x_0[j] + 0.5 * fl1_fx * tdy_yyy_0_0[j];

                tdz_xyyy_x_0[j] = pa_x[j] * tdz_yyy_x_0[j] + 0.5 * fl1_fx * tdz_yyy_0_0[j];

                tdx_xyyy_y_0[j] = pa_x[j] * tdx_yyy_y_0[j] + 0.5 * fl1_fx * ts_yyy_y_0[j];

                tdy_xyyy_y_0[j] = pa_x[j] * tdy_yyy_y_0[j];

                tdz_xyyy_y_0[j] = pa_x[j] * tdz_yyy_y_0[j];

                tdx_xyyy_z_0[j] = pa_x[j] * tdx_yyy_z_0[j] + 0.5 * fl1_fx * ts_yyy_z_0[j];

                tdy_xyyy_z_0[j] = pa_x[j] * tdy_yyy_z_0[j];

                tdz_xyyy_z_0[j] = pa_x[j] * tdz_yyy_z_0[j];

                tdx_xyyz_x_0[j] = pa_x[j] * tdx_yyz_x_0[j] + 0.5 * fl1_fx * tdx_yyz_0_0[j] + 0.5 * fl1_fx * ts_yyz_x_0[j];

                tdy_xyyz_x_0[j] = pa_x[j] * tdy_yyz_x_0[j] + 0.5 * fl1_fx * tdy_yyz_0_0[j];

                tdz_xyyz_x_0[j] = pa_x[j] * tdz_yyz_x_0[j] + 0.5 * fl1_fx * tdz_yyz_0_0[j];

                tdx_xyyz_y_0[j] = pa_x[j] * tdx_yyz_y_0[j] + 0.5 * fl1_fx * ts_yyz_y_0[j];

                tdy_xyyz_y_0[j] = pa_x[j] * tdy_yyz_y_0[j];

                tdz_xyyz_y_0[j] = pa_x[j] * tdz_yyz_y_0[j];

                tdx_xyyz_z_0[j] = pa_x[j] * tdx_yyz_z_0[j] + 0.5 * fl1_fx * ts_yyz_z_0[j];

                tdy_xyyz_z_0[j] = pa_x[j] * tdy_yyz_z_0[j];

                tdz_xyyz_z_0[j] = pa_x[j] * tdz_yyz_z_0[j];

                tdx_xyzz_x_0[j] = pa_x[j] * tdx_yzz_x_0[j] + 0.5 * fl1_fx * tdx_yzz_0_0[j] + 0.5 * fl1_fx * ts_yzz_x_0[j];

                tdy_xyzz_x_0[j] = pa_x[j] * tdy_yzz_x_0[j] + 0.5 * fl1_fx * tdy_yzz_0_0[j];

                tdz_xyzz_x_0[j] = pa_x[j] * tdz_yzz_x_0[j] + 0.5 * fl1_fx * tdz_yzz_0_0[j];

                tdx_xyzz_y_0[j] = pa_x[j] * tdx_yzz_y_0[j] + 0.5 * fl1_fx * ts_yzz_y_0[j];

                tdy_xyzz_y_0[j] = pa_x[j] * tdy_yzz_y_0[j];

                tdz_xyzz_y_0[j] = pa_x[j] * tdz_yzz_y_0[j];

                tdx_xyzz_z_0[j] = pa_x[j] * tdx_yzz_z_0[j] + 0.5 * fl1_fx * ts_yzz_z_0[j];

                tdy_xyzz_z_0[j] = pa_x[j] * tdy_yzz_z_0[j];

                tdz_xyzz_z_0[j] = pa_x[j] * tdz_yzz_z_0[j];

                tdx_xzzz_x_0[j] = pa_x[j] * tdx_zzz_x_0[j] + 0.5 * fl1_fx * tdx_zzz_0_0[j] + 0.5 * fl1_fx * ts_zzz_x_0[j];

                tdy_xzzz_x_0[j] = pa_x[j] * tdy_zzz_x_0[j] + 0.5 * fl1_fx * tdy_zzz_0_0[j];

                tdz_xzzz_x_0[j] = pa_x[j] * tdz_zzz_x_0[j] + 0.5 * fl1_fx * tdz_zzz_0_0[j];

                tdx_xzzz_y_0[j] = pa_x[j] * tdx_zzz_y_0[j] + 0.5 * fl1_fx * ts_zzz_y_0[j];

                tdy_xzzz_y_0[j] = pa_x[j] * tdy_zzz_y_0[j];

                tdz_xzzz_y_0[j] = pa_x[j] * tdz_zzz_y_0[j];

                tdx_xzzz_z_0[j] = pa_x[j] * tdx_zzz_z_0[j] + 0.5 * fl1_fx * ts_zzz_z_0[j];

                tdy_xzzz_z_0[j] = pa_x[j] * tdy_zzz_z_0[j];

                tdz_xzzz_z_0[j] = pa_x[j] * tdz_zzz_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGP_90_135(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_d_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto tdx_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 15); 

            auto tdy_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 15); 

            auto tdz_zz_x_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 15); 

            auto tdx_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 16); 

            auto tdy_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 16); 

            auto tdz_zz_y_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 16); 

            auto tdx_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * idx + 17); 

            auto tdy_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 18 * bdim + 18 * idx + 17); 

            auto tdz_zz_z_0 = primBuffer.data(pidx_d_2_1_m0 + 36 * bdim + 18 * idx + 17); 

            auto tdx_yyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 6); 

            auto tdy_yyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 6); 

            auto tdz_yyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 6); 

            auto tdx_yyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 7); 

            auto tdy_yyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 7); 

            auto tdz_yyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 7); 

            auto tdx_yzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 8); 

            auto tdy_yzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 8); 

            auto tdz_yzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 8); 

            auto tdx_zzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 9); 

            auto tdy_zzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 9); 

            auto tdz_zzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 9); 

            auto ts_yyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 18); 

            auto ts_yyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 19); 

            auto ts_yyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 20); 

            auto ts_yyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 21); 

            auto ts_yyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 22); 

            auto ts_yyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 23); 

            auto ts_yzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 24); 

            auto ts_yzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 25); 

            auto ts_yzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 26); 

            auto ts_zzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 27); 

            auto ts_zzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 28); 

            auto ts_zzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 29); 

            // set up pointers to integrals

            auto tdx_yyyy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 30); 

            auto tdy_yyyy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 30); 

            auto tdz_yyyy_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 30); 

            auto tdx_yyyy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 31); 

            auto tdy_yyyy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 31); 

            auto tdz_yyyy_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 31); 

            auto tdx_yyyy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 32); 

            auto tdy_yyyy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 32); 

            auto tdz_yyyy_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 32); 

            auto tdx_yyyz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 33); 

            auto tdy_yyyz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 33); 

            auto tdz_yyyz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 33); 

            auto tdx_yyyz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 34); 

            auto tdy_yyyz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 34); 

            auto tdz_yyyz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 34); 

            auto tdx_yyyz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 35); 

            auto tdy_yyyz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 35); 

            auto tdz_yyyz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 35); 

            auto tdx_yyzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 36); 

            auto tdy_yyzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 36); 

            auto tdz_yyzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 36); 

            auto tdx_yyzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 37); 

            auto tdy_yyzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 37); 

            auto tdz_yyzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 37); 

            auto tdx_yyzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 38); 

            auto tdy_yyzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 38); 

            auto tdz_yyzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 38); 

            auto tdx_yzzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 39); 

            auto tdy_yzzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 39); 

            auto tdz_yzzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 39); 

            auto tdx_yzzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 40); 

            auto tdy_yzzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 40); 

            auto tdz_yzzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 40); 

            auto tdx_yzzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 41); 

            auto tdy_yzzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 41); 

            auto tdz_yzzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 41); 

            auto tdx_zzzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 42); 

            auto tdy_zzzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 42); 

            auto tdz_zzzz_x_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 42); 

            auto tdx_zzzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 43); 

            auto tdy_zzzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 43); 

            auto tdz_zzzz_y_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 43); 

            auto tdx_zzzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * idx + 44); 

            auto tdy_zzzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 45 * bdim + 45 * idx + 44); 

            auto tdz_zzzz_z_0 = primBuffer.data(pidx_d_4_1_m0 + 90 * bdim + 45 * idx + 44); 

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fx, pa_y, pa_z, tdx_yy_x_0, tdx_yy_y_0, tdx_yy_z_0, tdx_yyy_0_0, \
                                     tdx_yyy_x_0, tdx_yyy_y_0, tdx_yyy_z_0, tdx_yyyy_x_0, tdx_yyyy_y_0, tdx_yyyy_z_0, \
                                     tdx_yyyz_x_0, tdx_yyyz_y_0, tdx_yyyz_z_0, tdx_yyz_0_0, tdx_yyz_x_0, tdx_yyz_y_0, \
                                     tdx_yyz_z_0, tdx_yyzz_x_0, tdx_yyzz_y_0, tdx_yyzz_z_0, tdx_yz_x_0, tdx_yz_y_0, \
                                     tdx_yz_z_0, tdx_yzz_0_0, tdx_yzz_x_0, tdx_yzz_y_0, tdx_yzz_z_0, tdx_yzzz_x_0, \
                                     tdx_yzzz_y_0, tdx_yzzz_z_0, tdx_zz_x_0, tdx_zz_y_0, tdx_zz_z_0, tdx_zzz_0_0, \
                                     tdx_zzz_x_0, tdx_zzz_y_0, tdx_zzz_z_0, tdx_zzzz_x_0, tdx_zzzz_y_0, tdx_zzzz_z_0, \
                                     tdy_yy_x_0, tdy_yy_y_0, tdy_yy_z_0, tdy_yyy_0_0, tdy_yyy_x_0, tdy_yyy_y_0, \
                                     tdy_yyy_z_0, tdy_yyyy_x_0, tdy_yyyy_y_0, tdy_yyyy_z_0, tdy_yyyz_x_0, tdy_yyyz_y_0, \
                                     tdy_yyyz_z_0, tdy_yyz_0_0, tdy_yyz_x_0, tdy_yyz_y_0, tdy_yyz_z_0, tdy_yyzz_x_0, \
                                     tdy_yyzz_y_0, tdy_yyzz_z_0, tdy_yz_x_0, tdy_yz_y_0, tdy_yz_z_0, tdy_yzz_0_0, \
                                     tdy_yzz_x_0, tdy_yzz_y_0, tdy_yzz_z_0, tdy_yzzz_x_0, tdy_yzzz_y_0, tdy_yzzz_z_0, \
                                     tdy_zz_x_0, tdy_zz_y_0, tdy_zz_z_0, tdy_zzz_0_0, tdy_zzz_x_0, tdy_zzz_y_0, \
                                     tdy_zzz_z_0, tdy_zzzz_x_0, tdy_zzzz_y_0, tdy_zzzz_z_0, tdz_yy_x_0, tdz_yy_y_0, \
                                     tdz_yy_z_0, tdz_yyy_0_0, tdz_yyy_x_0, tdz_yyy_y_0, tdz_yyy_z_0, tdz_yyyy_x_0, \
                                     tdz_yyyy_y_0, tdz_yyyy_z_0, tdz_yyyz_x_0, tdz_yyyz_y_0, tdz_yyyz_z_0, tdz_yyz_0_0, \
                                     tdz_yyz_x_0, tdz_yyz_y_0, tdz_yyz_z_0, tdz_yyzz_x_0, tdz_yyzz_y_0, tdz_yyzz_z_0, \
                                     tdz_yz_x_0, tdz_yz_y_0, tdz_yz_z_0, tdz_yzz_0_0, tdz_yzz_x_0, tdz_yzz_y_0, \
                                     tdz_yzz_z_0, tdz_yzzz_x_0, tdz_yzzz_y_0, tdz_yzzz_z_0, tdz_zz_x_0, tdz_zz_y_0, \
                                     tdz_zz_z_0, tdz_zzz_0_0, tdz_zzz_x_0, tdz_zzz_y_0, tdz_zzz_z_0, tdz_zzzz_x_0, \
                                     tdz_zzzz_y_0, tdz_zzzz_z_0, ts_yyy_x_0, ts_yyy_y_0, ts_yyy_z_0, ts_yyz_x_0, \
                                     ts_yyz_y_0, ts_yyz_z_0, ts_yzz_x_0, ts_yzz_y_0, ts_yzz_z_0, ts_zzz_x_0, ts_zzz_y_0, \
                                     ts_zzz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_yyyy_x_0[j] = pa_y[j] * tdx_yyy_x_0[j] + 1.5 * fl1_fx * tdx_yy_x_0[j];

                tdy_yyyy_x_0[j] = pa_y[j] * tdy_yyy_x_0[j] + 1.5 * fl1_fx * tdy_yy_x_0[j] + 0.5 * fl1_fx * ts_yyy_x_0[j];

                tdz_yyyy_x_0[j] = pa_y[j] * tdz_yyy_x_0[j] + 1.5 * fl1_fx * tdz_yy_x_0[j];

                tdx_yyyy_y_0[j] = pa_y[j] * tdx_yyy_y_0[j] + 1.5 * fl1_fx * tdx_yy_y_0[j] + 0.5 * fl1_fx * tdx_yyy_0_0[j];

                tdy_yyyy_y_0[j] = pa_y[j] * tdy_yyy_y_0[j] + 1.5 * fl1_fx * tdy_yy_y_0[j] + 0.5 * fl1_fx * tdy_yyy_0_0[j] + 0.5 * fl1_fx * ts_yyy_y_0[j];

                tdz_yyyy_y_0[j] = pa_y[j] * tdz_yyy_y_0[j] + 1.5 * fl1_fx * tdz_yy_y_0[j] + 0.5 * fl1_fx * tdz_yyy_0_0[j];

                tdx_yyyy_z_0[j] = pa_y[j] * tdx_yyy_z_0[j] + 1.5 * fl1_fx * tdx_yy_z_0[j];

                tdy_yyyy_z_0[j] = pa_y[j] * tdy_yyy_z_0[j] + 1.5 * fl1_fx * tdy_yy_z_0[j] + 0.5 * fl1_fx * ts_yyy_z_0[j];

                tdz_yyyy_z_0[j] = pa_y[j] * tdz_yyy_z_0[j] + 1.5 * fl1_fx * tdz_yy_z_0[j];

                tdx_yyyz_x_0[j] = pa_y[j] * tdx_yyz_x_0[j] + fl1_fx * tdx_yz_x_0[j];

                tdy_yyyz_x_0[j] = pa_y[j] * tdy_yyz_x_0[j] + fl1_fx * tdy_yz_x_0[j] + 0.5 * fl1_fx * ts_yyz_x_0[j];

                tdz_yyyz_x_0[j] = pa_y[j] * tdz_yyz_x_0[j] + fl1_fx * tdz_yz_x_0[j];

                tdx_yyyz_y_0[j] = pa_y[j] * tdx_yyz_y_0[j] + fl1_fx * tdx_yz_y_0[j] + 0.5 * fl1_fx * tdx_yyz_0_0[j];

                tdy_yyyz_y_0[j] = pa_y[j] * tdy_yyz_y_0[j] + fl1_fx * tdy_yz_y_0[j] + 0.5 * fl1_fx * tdy_yyz_0_0[j] + 0.5 * fl1_fx * ts_yyz_y_0[j];

                tdz_yyyz_y_0[j] = pa_y[j] * tdz_yyz_y_0[j] + fl1_fx * tdz_yz_y_0[j] + 0.5 * fl1_fx * tdz_yyz_0_0[j];

                tdx_yyyz_z_0[j] = pa_y[j] * tdx_yyz_z_0[j] + fl1_fx * tdx_yz_z_0[j];

                tdy_yyyz_z_0[j] = pa_y[j] * tdy_yyz_z_0[j] + fl1_fx * tdy_yz_z_0[j] + 0.5 * fl1_fx * ts_yyz_z_0[j];

                tdz_yyyz_z_0[j] = pa_y[j] * tdz_yyz_z_0[j] + fl1_fx * tdz_yz_z_0[j];

                tdx_yyzz_x_0[j] = pa_y[j] * tdx_yzz_x_0[j] + 0.5 * fl1_fx * tdx_zz_x_0[j];

                tdy_yyzz_x_0[j] = pa_y[j] * tdy_yzz_x_0[j] + 0.5 * fl1_fx * tdy_zz_x_0[j] + 0.5 * fl1_fx * ts_yzz_x_0[j];

                tdz_yyzz_x_0[j] = pa_y[j] * tdz_yzz_x_0[j] + 0.5 * fl1_fx * tdz_zz_x_0[j];

                tdx_yyzz_y_0[j] = pa_y[j] * tdx_yzz_y_0[j] + 0.5 * fl1_fx * tdx_zz_y_0[j] + 0.5 * fl1_fx * tdx_yzz_0_0[j];

                tdy_yyzz_y_0[j] = pa_y[j] * tdy_yzz_y_0[j] + 0.5 * fl1_fx * tdy_zz_y_0[j] + 0.5 * fl1_fx * tdy_yzz_0_0[j] + 0.5 * fl1_fx * ts_yzz_y_0[j];

                tdz_yyzz_y_0[j] = pa_y[j] * tdz_yzz_y_0[j] + 0.5 * fl1_fx * tdz_zz_y_0[j] + 0.5 * fl1_fx * tdz_yzz_0_0[j];

                tdx_yyzz_z_0[j] = pa_y[j] * tdx_yzz_z_0[j] + 0.5 * fl1_fx * tdx_zz_z_0[j];

                tdy_yyzz_z_0[j] = pa_y[j] * tdy_yzz_z_0[j] + 0.5 * fl1_fx * tdy_zz_z_0[j] + 0.5 * fl1_fx * ts_yzz_z_0[j];

                tdz_yyzz_z_0[j] = pa_y[j] * tdz_yzz_z_0[j] + 0.5 * fl1_fx * tdz_zz_z_0[j];

                tdx_yzzz_x_0[j] = pa_y[j] * tdx_zzz_x_0[j];

                tdy_yzzz_x_0[j] = pa_y[j] * tdy_zzz_x_0[j] + 0.5 * fl1_fx * ts_zzz_x_0[j];

                tdz_yzzz_x_0[j] = pa_y[j] * tdz_zzz_x_0[j];

                tdx_yzzz_y_0[j] = pa_y[j] * tdx_zzz_y_0[j] + 0.5 * fl1_fx * tdx_zzz_0_0[j];

                tdy_yzzz_y_0[j] = pa_y[j] * tdy_zzz_y_0[j] + 0.5 * fl1_fx * tdy_zzz_0_0[j] + 0.5 * fl1_fx * ts_zzz_y_0[j];

                tdz_yzzz_y_0[j] = pa_y[j] * tdz_zzz_y_0[j] + 0.5 * fl1_fx * tdz_zzz_0_0[j];

                tdx_yzzz_z_0[j] = pa_y[j] * tdx_zzz_z_0[j];

                tdy_yzzz_z_0[j] = pa_y[j] * tdy_zzz_z_0[j] + 0.5 * fl1_fx * ts_zzz_z_0[j];

                tdz_yzzz_z_0[j] = pa_y[j] * tdz_zzz_z_0[j];

                tdx_zzzz_x_0[j] = pa_z[j] * tdx_zzz_x_0[j] + 1.5 * fl1_fx * tdx_zz_x_0[j];

                tdy_zzzz_x_0[j] = pa_z[j] * tdy_zzz_x_0[j] + 1.5 * fl1_fx * tdy_zz_x_0[j];

                tdz_zzzz_x_0[j] = pa_z[j] * tdz_zzz_x_0[j] + 1.5 * fl1_fx * tdz_zz_x_0[j] + 0.5 * fl1_fx * ts_zzz_x_0[j];

                tdx_zzzz_y_0[j] = pa_z[j] * tdx_zzz_y_0[j] + 1.5 * fl1_fx * tdx_zz_y_0[j];

                tdy_zzzz_y_0[j] = pa_z[j] * tdy_zzz_y_0[j] + 1.5 * fl1_fx * tdy_zz_y_0[j];

                tdz_zzzz_y_0[j] = pa_z[j] * tdz_zzz_y_0[j] + 1.5 * fl1_fx * tdz_zz_y_0[j] + 0.5 * fl1_fx * ts_zzz_y_0[j];

                tdx_zzzz_z_0[j] = pa_z[j] * tdx_zzz_z_0[j] + 1.5 * fl1_fx * tdx_zz_z_0[j] + 0.5 * fl1_fx * tdx_zzz_0_0[j];

                tdy_zzzz_z_0[j] = pa_z[j] * tdy_zzz_z_0[j] + 1.5 * fl1_fx * tdy_zz_z_0[j] + 0.5 * fl1_fx * tdy_zzz_0_0[j];

                tdz_zzzz_z_0[j] = pa_z[j] * tdz_zzz_z_0[j] + 1.5 * fl1_fx * tdz_zz_z_0[j] + 0.5 * fl1_fx * tdz_zzz_0_0[j] + 0.5 * fl1_fx * ts_zzz_z_0[j];
            }

            idx++;
        }
    }


} // ediprecfunc namespace

