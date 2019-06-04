//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForPX.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForPP(      CMemBlock2D<double>& primBuffer,
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

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_t_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_1_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_0_x_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx); 

            auto tt_0_y_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 1); 

            auto tt_0_z_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 2); 

            auto tt_0_0_0 = primBuffer.data(pidx_t_0_0_m0 + idx); 

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

            auto tt_x_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx); 

            auto tt_x_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 1); 

            auto tt_x_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 2); 

            auto tt_y_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 3); 

            auto tt_y_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 4); 

            auto tt_y_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 5); 

            auto tt_z_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 6); 

            auto tt_z_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 7); 

            auto tt_z_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 8); 

            #pragma omp simd aligned(fx, fz, pa_x, pa_y, pa_z, ts_x_x_0, ts_x_y_0, ts_x_z_0, ts_y_x_0, ts_y_y_0, \
                                     ts_y_z_0, ts_z_x_0, ts_z_y_0, ts_z_z_0, tt_0_0_0, tt_0_x_0, tt_0_y_0, tt_0_z_0, \
                                     tt_x_x_0, tt_x_y_0, tt_x_z_0, tt_y_x_0, tt_y_y_0, tt_y_z_0, tt_z_x_0, tt_z_y_0, \
                                     tt_z_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_x_x_0[j] = pa_x[j] * tt_0_x_0[j] + 0.5 * fl1_fx * tt_0_0_0[j] + 2.0 * fl1_fz * ts_x_x_0[j];

                tt_x_y_0[j] = pa_x[j] * tt_0_y_0[j] + 2.0 * fl1_fz * ts_x_y_0[j];

                tt_x_z_0[j] = pa_x[j] * tt_0_z_0[j] + 2.0 * fl1_fz * ts_x_z_0[j];

                tt_y_x_0[j] = pa_y[j] * tt_0_x_0[j] + 2.0 * fl1_fz * ts_y_x_0[j];

                tt_y_y_0[j] = pa_y[j] * tt_0_y_0[j] + 0.5 * fl1_fx * tt_0_0_0[j] + 2.0 * fl1_fz * ts_y_y_0[j];

                tt_y_z_0[j] = pa_y[j] * tt_0_z_0[j] + 2.0 * fl1_fz * ts_y_z_0[j];

                tt_z_x_0[j] = pa_z[j] * tt_0_x_0[j] + 2.0 * fl1_fz * ts_z_x_0[j];

                tt_z_y_0[j] = pa_z[j] * tt_0_y_0[j] + 2.0 * fl1_fz * ts_z_y_0[j];

                tt_z_z_0[j] = pa_z[j] * tt_0_z_0[j] + 0.5 * fl1_fx * tt_0_0_0[j] + 2.0 * fl1_fz * ts_z_z_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForPD(      CMemBlock2D<double>& primBuffer,
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

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_t_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_1_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_0_xx_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx); 

            auto tt_0_xy_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 1); 

            auto tt_0_xz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 2); 

            auto tt_0_yy_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 3); 

            auto tt_0_yz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 4); 

            auto tt_0_zz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 5); 

            auto tt_0_x_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx); 

            auto tt_0_y_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 1); 

            auto tt_0_z_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 2); 

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

            auto ts_z_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 12); 

            auto ts_z_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 13); 

            auto ts_z_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 14); 

            auto ts_z_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 15); 

            auto ts_z_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 16); 

            auto ts_z_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 17); 

            // set up pointers to integrals

            auto tt_x_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx); 

            auto tt_x_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 1); 

            auto tt_x_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 2); 

            auto tt_x_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 3); 

            auto tt_x_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 4); 

            auto tt_x_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 5); 

            auto tt_y_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 6); 

            auto tt_y_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 7); 

            auto tt_y_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 8); 

            auto tt_y_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 9); 

            auto tt_y_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 10); 

            auto tt_y_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 11); 

            auto tt_z_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 12); 

            auto tt_z_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 13); 

            auto tt_z_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 14); 

            auto tt_z_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 15); 

            auto tt_z_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 16); 

            auto tt_z_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 17); 

            #pragma omp simd aligned(fx, fz, pa_x, pa_y, pa_z, ts_x_xx_0, ts_x_xy_0, ts_x_xz_0, ts_x_yy_0, \
                                     ts_x_yz_0, ts_x_zz_0, ts_y_xx_0, ts_y_xy_0, ts_y_xz_0, ts_y_yy_0, ts_y_yz_0, \
                                     ts_y_zz_0, ts_z_xx_0, ts_z_xy_0, ts_z_xz_0, ts_z_yy_0, ts_z_yz_0, ts_z_zz_0, \
                                     tt_0_x_0, tt_0_xx_0, tt_0_xy_0, tt_0_xz_0, tt_0_y_0, tt_0_yy_0, tt_0_yz_0, \
                                     tt_0_z_0, tt_0_zz_0, tt_x_xx_0, tt_x_xy_0, tt_x_xz_0, tt_x_yy_0, tt_x_yz_0, \
                                     tt_x_zz_0, tt_y_xx_0, tt_y_xy_0, tt_y_xz_0, tt_y_yy_0, tt_y_yz_0, tt_y_zz_0, \
                                     tt_z_xx_0, tt_z_xy_0, tt_z_xz_0, tt_z_yy_0, tt_z_yz_0, tt_z_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_x_xx_0[j] = pa_x[j] * tt_0_xx_0[j] + fl1_fx * tt_0_x_0[j] + 2.0 * fl1_fz * ts_x_xx_0[j];

                tt_x_xy_0[j] = pa_x[j] * tt_0_xy_0[j] + 0.5 * fl1_fx * tt_0_y_0[j] + 2.0 * fl1_fz * ts_x_xy_0[j];

                tt_x_xz_0[j] = pa_x[j] * tt_0_xz_0[j] + 0.5 * fl1_fx * tt_0_z_0[j] + 2.0 * fl1_fz * ts_x_xz_0[j];

                tt_x_yy_0[j] = pa_x[j] * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_x_yy_0[j];

                tt_x_yz_0[j] = pa_x[j] * tt_0_yz_0[j] + 2.0 * fl1_fz * ts_x_yz_0[j];

                tt_x_zz_0[j] = pa_x[j] * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_x_zz_0[j];

                tt_y_xx_0[j] = pa_y[j] * tt_0_xx_0[j] + 2.0 * fl1_fz * ts_y_xx_0[j];

                tt_y_xy_0[j] = pa_y[j] * tt_0_xy_0[j] + 0.5 * fl1_fx * tt_0_x_0[j] + 2.0 * fl1_fz * ts_y_xy_0[j];

                tt_y_xz_0[j] = pa_y[j] * tt_0_xz_0[j] + 2.0 * fl1_fz * ts_y_xz_0[j];

                tt_y_yy_0[j] = pa_y[j] * tt_0_yy_0[j] + fl1_fx * tt_0_y_0[j] + 2.0 * fl1_fz * ts_y_yy_0[j];

                tt_y_yz_0[j] = pa_y[j] * tt_0_yz_0[j] + 0.5 * fl1_fx * tt_0_z_0[j] + 2.0 * fl1_fz * ts_y_yz_0[j];

                tt_y_zz_0[j] = pa_y[j] * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_y_zz_0[j];

                tt_z_xx_0[j] = pa_z[j] * tt_0_xx_0[j] + 2.0 * fl1_fz * ts_z_xx_0[j];

                tt_z_xy_0[j] = pa_z[j] * tt_0_xy_0[j] + 2.0 * fl1_fz * ts_z_xy_0[j];

                tt_z_xz_0[j] = pa_z[j] * tt_0_xz_0[j] + 0.5 * fl1_fx * tt_0_x_0[j] + 2.0 * fl1_fz * ts_z_xz_0[j];

                tt_z_yy_0[j] = pa_z[j] * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_z_yy_0[j];

                tt_z_yz_0[j] = pa_z[j] * tt_0_yz_0[j] + 0.5 * fl1_fx * tt_0_y_0[j] + 2.0 * fl1_fz * ts_z_yz_0[j];

                tt_z_zz_0[j] = pa_z[j] * tt_0_zz_0[j] + fl1_fx * tt_0_z_0[j] + 2.0 * fl1_fz * ts_z_zz_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForDP(      CMemBlock2D<double>& primBuffer,
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

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_t_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_2_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_x_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx); 

            auto tt_x_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 1); 

            auto tt_x_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 2); 

            auto tt_y_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 3); 

            auto tt_y_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 4); 

            auto tt_y_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 5); 

            auto tt_z_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 6); 

            auto tt_z_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 7); 

            auto tt_z_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 8); 

            auto tt_0_x_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx); 

            auto tt_0_y_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 1); 

            auto tt_0_z_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 2); 

            auto tt_x_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx); 

            auto tt_y_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx + 1); 

            auto tt_z_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx + 2); 

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

            auto ts_zz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 15); 

            auto ts_zz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 16); 

            auto ts_zz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 17); 

            auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx); 

            auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1); 

            auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2); 

            // set up pointers to integrals

            auto tt_xx_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx); 

            auto tt_xx_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 1); 

            auto tt_xx_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 2); 

            auto tt_xy_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 3); 

            auto tt_xy_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 4); 

            auto tt_xy_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 5); 

            auto tt_xz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 6); 

            auto tt_xz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 7); 

            auto tt_xz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 8); 

            auto tt_yy_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 9); 

            auto tt_yy_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 10); 

            auto tt_yy_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 11); 

            auto tt_yz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 12); 

            auto tt_yz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 13); 

            auto tt_yz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 14); 

            auto tt_zz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 15); 

            auto tt_zz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 16); 

            auto tt_zz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 17); 

            #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, pa_z, ts_0_x_0, ts_0_y_0, ts_0_z_0, ts_xx_x_0, \
                                     ts_xx_y_0, ts_xx_z_0, ts_xy_x_0, ts_xy_y_0, ts_xy_z_0, ts_xz_x_0, ts_xz_y_0, \
                                     ts_xz_z_0, ts_yy_x_0, ts_yy_y_0, ts_yy_z_0, ts_yz_x_0, ts_yz_y_0, ts_yz_z_0, \
                                     ts_zz_x_0, ts_zz_y_0, ts_zz_z_0, tt_0_x_0, tt_0_y_0, tt_0_z_0, tt_x_0_0, tt_x_x_0, \
                                     tt_x_y_0, tt_x_z_0, tt_xx_x_0, tt_xx_y_0, tt_xx_z_0, tt_xy_x_0, tt_xy_y_0, \
                                     tt_xy_z_0, tt_xz_x_0, tt_xz_y_0, tt_xz_z_0, tt_y_0_0, tt_y_x_0, tt_y_y_0, tt_y_z_0, \
                                     tt_yy_x_0, tt_yy_y_0, tt_yy_z_0, tt_yz_x_0, tt_yz_y_0, tt_yz_z_0, tt_z_0_0, \
                                     tt_z_x_0, tt_z_y_0, tt_z_z_0, tt_zz_x_0, tt_zz_y_0, tt_zz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fga = fga[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_xx_x_0[j] = pa_x[j] * tt_x_x_0[j] + 0.5 * fl1_fx * tt_0_x_0[j] + 0.5 * fl1_fx * tt_x_0_0[j] + 2.0 * fl1_fz * ts_xx_x_0[j] - fl1_fz * fl1_fga * ts_0_x_0[j];

                tt_xx_y_0[j] = pa_x[j] * tt_x_y_0[j] + 0.5 * fl1_fx * tt_0_y_0[j] + 2.0 * fl1_fz * ts_xx_y_0[j] - fl1_fz * fl1_fga * ts_0_y_0[j];

                tt_xx_z_0[j] = pa_x[j] * tt_x_z_0[j] + 0.5 * fl1_fx * tt_0_z_0[j] + 2.0 * fl1_fz * ts_xx_z_0[j] - fl1_fz * fl1_fga * ts_0_z_0[j];

                tt_xy_x_0[j] = pa_x[j] * tt_y_x_0[j] + 0.5 * fl1_fx * tt_y_0_0[j] + 2.0 * fl1_fz * ts_xy_x_0[j];

                tt_xy_y_0[j] = pa_x[j] * tt_y_y_0[j] + 2.0 * fl1_fz * ts_xy_y_0[j];

                tt_xy_z_0[j] = pa_x[j] * tt_y_z_0[j] + 2.0 * fl1_fz * ts_xy_z_0[j];

                tt_xz_x_0[j] = pa_x[j] * tt_z_x_0[j] + 0.5 * fl1_fx * tt_z_0_0[j] + 2.0 * fl1_fz * ts_xz_x_0[j];

                tt_xz_y_0[j] = pa_x[j] * tt_z_y_0[j] + 2.0 * fl1_fz * ts_xz_y_0[j];

                tt_xz_z_0[j] = pa_x[j] * tt_z_z_0[j] + 2.0 * fl1_fz * ts_xz_z_0[j];

                tt_yy_x_0[j] = pa_y[j] * tt_y_x_0[j] + 0.5 * fl1_fx * tt_0_x_0[j] + 2.0 * fl1_fz * ts_yy_x_0[j] - fl1_fz * fl1_fga * ts_0_x_0[j];

                tt_yy_y_0[j] = pa_y[j] * tt_y_y_0[j] + 0.5 * fl1_fx * tt_0_y_0[j] + 0.5 * fl1_fx * tt_y_0_0[j] + 2.0 * fl1_fz * ts_yy_y_0[j] - fl1_fz * fl1_fga * ts_0_y_0[j];

                tt_yy_z_0[j] = pa_y[j] * tt_y_z_0[j] + 0.5 * fl1_fx * tt_0_z_0[j] + 2.0 * fl1_fz * ts_yy_z_0[j] - fl1_fz * fl1_fga * ts_0_z_0[j];

                tt_yz_x_0[j] = pa_y[j] * tt_z_x_0[j] + 2.0 * fl1_fz * ts_yz_x_0[j];

                tt_yz_y_0[j] = pa_y[j] * tt_z_y_0[j] + 0.5 * fl1_fx * tt_z_0_0[j] + 2.0 * fl1_fz * ts_yz_y_0[j];

                tt_yz_z_0[j] = pa_y[j] * tt_z_z_0[j] + 2.0 * fl1_fz * ts_yz_z_0[j];

                tt_zz_x_0[j] = pa_z[j] * tt_z_x_0[j] + 0.5 * fl1_fx * tt_0_x_0[j] + 2.0 * fl1_fz * ts_zz_x_0[j] - fl1_fz * fl1_fga * ts_0_x_0[j];

                tt_zz_y_0[j] = pa_z[j] * tt_z_y_0[j] + 0.5 * fl1_fx * tt_0_y_0[j] + 2.0 * fl1_fz * ts_zz_y_0[j] - fl1_fz * fl1_fga * ts_0_y_0[j];

                tt_zz_z_0[j] = pa_z[j] * tt_z_z_0[j] + 0.5 * fl1_fx * tt_0_z_0[j] + 0.5 * fl1_fx * tt_z_0_0[j] + 2.0 * fl1_fz * ts_zz_z_0[j] - fl1_fz * fl1_fga * ts_0_z_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForPF(      CMemBlock2D<double>& primBuffer,
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

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_t_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_1_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_0_xxx_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx); 

            auto tt_0_xxy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 1); 

            auto tt_0_xxz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 2); 

            auto tt_0_xyy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 3); 

            auto tt_0_xyz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 4); 

            auto tt_0_xzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 5); 

            auto tt_0_yyy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 6); 

            auto tt_0_yyz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 7); 

            auto tt_0_yzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 8); 

            auto tt_0_zzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 9); 

            auto tt_0_xx_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx); 

            auto tt_0_xy_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 1); 

            auto tt_0_xz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 2); 

            auto tt_0_yy_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 3); 

            auto tt_0_yz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 4); 

            auto tt_0_zz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 5); 

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

            auto tt_x_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx); 

            auto tt_x_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 1); 

            auto tt_x_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 2); 

            auto tt_x_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 3); 

            auto tt_x_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 4); 

            auto tt_x_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 5); 

            auto tt_x_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 6); 

            auto tt_x_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 7); 

            auto tt_x_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 8); 

            auto tt_x_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 9); 

            auto tt_y_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 10); 

            auto tt_y_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 11); 

            auto tt_y_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 12); 

            auto tt_y_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 13); 

            auto tt_y_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 14); 

            auto tt_y_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 15); 

            auto tt_y_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 16); 

            auto tt_y_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 17); 

            auto tt_y_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 18); 

            auto tt_y_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 19); 

            auto tt_z_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 20); 

            auto tt_z_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 21); 

            auto tt_z_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 22); 

            auto tt_z_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 23); 

            auto tt_z_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 24); 

            auto tt_z_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 25); 

            auto tt_z_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 26); 

            auto tt_z_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 27); 

            auto tt_z_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 28); 

            auto tt_z_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 29); 

            #pragma omp simd aligned(fx, fz, pa_x, pa_y, pa_z, ts_x_xxx_0, ts_x_xxy_0, ts_x_xxz_0, ts_x_xyy_0, \
                                     ts_x_xyz_0, ts_x_xzz_0, ts_x_yyy_0, ts_x_yyz_0, ts_x_yzz_0, ts_x_zzz_0, ts_y_xxx_0, \
                                     ts_y_xxy_0, ts_y_xxz_0, ts_y_xyy_0, ts_y_xyz_0, ts_y_xzz_0, ts_y_yyy_0, ts_y_yyz_0, \
                                     ts_y_yzz_0, ts_y_zzz_0, ts_z_xxx_0, ts_z_xxy_0, ts_z_xxz_0, ts_z_xyy_0, ts_z_xyz_0, \
                                     ts_z_xzz_0, ts_z_yyy_0, ts_z_yyz_0, ts_z_yzz_0, ts_z_zzz_0, tt_0_xx_0, tt_0_xxx_0, \
                                     tt_0_xxy_0, tt_0_xxz_0, tt_0_xy_0, tt_0_xyy_0, tt_0_xyz_0, tt_0_xz_0, tt_0_xzz_0, \
                                     tt_0_yy_0, tt_0_yyy_0, tt_0_yyz_0, tt_0_yz_0, tt_0_yzz_0, tt_0_zz_0, tt_0_zzz_0, \
                                     tt_x_xxx_0, tt_x_xxy_0, tt_x_xxz_0, tt_x_xyy_0, tt_x_xyz_0, tt_x_xzz_0, tt_x_yyy_0, \
                                     tt_x_yyz_0, tt_x_yzz_0, tt_x_zzz_0, tt_y_xxx_0, tt_y_xxy_0, tt_y_xxz_0, tt_y_xyy_0, \
                                     tt_y_xyz_0, tt_y_xzz_0, tt_y_yyy_0, tt_y_yyz_0, tt_y_yzz_0, tt_y_zzz_0, tt_z_xxx_0, \
                                     tt_z_xxy_0, tt_z_xxz_0, tt_z_xyy_0, tt_z_xyz_0, tt_z_xzz_0, tt_z_yyy_0, tt_z_yyz_0, \
                                     tt_z_yzz_0, tt_z_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_x_xxx_0[j] = pa_x[j] * tt_0_xxx_0[j] + 1.5 * fl1_fx * tt_0_xx_0[j] + 2.0 * fl1_fz * ts_x_xxx_0[j];

                tt_x_xxy_0[j] = pa_x[j] * tt_0_xxy_0[j] + fl1_fx * tt_0_xy_0[j] + 2.0 * fl1_fz * ts_x_xxy_0[j];

                tt_x_xxz_0[j] = pa_x[j] * tt_0_xxz_0[j] + fl1_fx * tt_0_xz_0[j] + 2.0 * fl1_fz * ts_x_xxz_0[j];

                tt_x_xyy_0[j] = pa_x[j] * tt_0_xyy_0[j] + 0.5 * fl1_fx * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_x_xyy_0[j];

                tt_x_xyz_0[j] = pa_x[j] * tt_0_xyz_0[j] + 0.5 * fl1_fx * tt_0_yz_0[j] + 2.0 * fl1_fz * ts_x_xyz_0[j];

                tt_x_xzz_0[j] = pa_x[j] * tt_0_xzz_0[j] + 0.5 * fl1_fx * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_x_xzz_0[j];

                tt_x_yyy_0[j] = pa_x[j] * tt_0_yyy_0[j] + 2.0 * fl1_fz * ts_x_yyy_0[j];

                tt_x_yyz_0[j] = pa_x[j] * tt_0_yyz_0[j] + 2.0 * fl1_fz * ts_x_yyz_0[j];

                tt_x_yzz_0[j] = pa_x[j] * tt_0_yzz_0[j] + 2.0 * fl1_fz * ts_x_yzz_0[j];

                tt_x_zzz_0[j] = pa_x[j] * tt_0_zzz_0[j] + 2.0 * fl1_fz * ts_x_zzz_0[j];

                tt_y_xxx_0[j] = pa_y[j] * tt_0_xxx_0[j] + 2.0 * fl1_fz * ts_y_xxx_0[j];

                tt_y_xxy_0[j] = pa_y[j] * tt_0_xxy_0[j] + 0.5 * fl1_fx * tt_0_xx_0[j] + 2.0 * fl1_fz * ts_y_xxy_0[j];

                tt_y_xxz_0[j] = pa_y[j] * tt_0_xxz_0[j] + 2.0 * fl1_fz * ts_y_xxz_0[j];

                tt_y_xyy_0[j] = pa_y[j] * tt_0_xyy_0[j] + fl1_fx * tt_0_xy_0[j] + 2.0 * fl1_fz * ts_y_xyy_0[j];

                tt_y_xyz_0[j] = pa_y[j] * tt_0_xyz_0[j] + 0.5 * fl1_fx * tt_0_xz_0[j] + 2.0 * fl1_fz * ts_y_xyz_0[j];

                tt_y_xzz_0[j] = pa_y[j] * tt_0_xzz_0[j] + 2.0 * fl1_fz * ts_y_xzz_0[j];

                tt_y_yyy_0[j] = pa_y[j] * tt_0_yyy_0[j] + 1.5 * fl1_fx * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_y_yyy_0[j];

                tt_y_yyz_0[j] = pa_y[j] * tt_0_yyz_0[j] + fl1_fx * tt_0_yz_0[j] + 2.0 * fl1_fz * ts_y_yyz_0[j];

                tt_y_yzz_0[j] = pa_y[j] * tt_0_yzz_0[j] + 0.5 * fl1_fx * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_y_yzz_0[j];

                tt_y_zzz_0[j] = pa_y[j] * tt_0_zzz_0[j] + 2.0 * fl1_fz * ts_y_zzz_0[j];

                tt_z_xxx_0[j] = pa_z[j] * tt_0_xxx_0[j] + 2.0 * fl1_fz * ts_z_xxx_0[j];

                tt_z_xxy_0[j] = pa_z[j] * tt_0_xxy_0[j] + 2.0 * fl1_fz * ts_z_xxy_0[j];

                tt_z_xxz_0[j] = pa_z[j] * tt_0_xxz_0[j] + 0.5 * fl1_fx * tt_0_xx_0[j] + 2.0 * fl1_fz * ts_z_xxz_0[j];

                tt_z_xyy_0[j] = pa_z[j] * tt_0_xyy_0[j] + 2.0 * fl1_fz * ts_z_xyy_0[j];

                tt_z_xyz_0[j] = pa_z[j] * tt_0_xyz_0[j] + 0.5 * fl1_fx * tt_0_xy_0[j] + 2.0 * fl1_fz * ts_z_xyz_0[j];

                tt_z_xzz_0[j] = pa_z[j] * tt_0_xzz_0[j] + fl1_fx * tt_0_xz_0[j] + 2.0 * fl1_fz * ts_z_xzz_0[j];

                tt_z_yyy_0[j] = pa_z[j] * tt_0_yyy_0[j] + 2.0 * fl1_fz * ts_z_yyy_0[j];

                tt_z_yyz_0[j] = pa_z[j] * tt_0_yyz_0[j] + 0.5 * fl1_fx * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_z_yyz_0[j];

                tt_z_yzz_0[j] = pa_z[j] * tt_0_yzz_0[j] + fl1_fx * tt_0_yz_0[j] + 2.0 * fl1_fz * ts_z_yzz_0[j];

                tt_z_zzz_0[j] = pa_z[j] * tt_0_zzz_0[j] + 1.5 * fl1_fx * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_z_zzz_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFP(      CMemBlock2D<double>& primBuffer,
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

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_t_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_3_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_xx_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx); 

            auto tt_xx_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 1); 

            auto tt_xx_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 2); 

            auto tt_xy_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 3); 

            auto tt_xy_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 4); 

            auto tt_xy_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 5); 

            auto tt_xz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 6); 

            auto tt_xz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 7); 

            auto tt_xz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 8); 

            auto tt_yy_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 9); 

            auto tt_yy_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 10); 

            auto tt_yy_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 11); 

            auto tt_yz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 12); 

            auto tt_yz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 13); 

            auto tt_yz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 14); 

            auto tt_zz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 15); 

            auto tt_zz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 16); 

            auto tt_zz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 17); 

            auto tt_x_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx); 

            auto tt_x_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 1); 

            auto tt_x_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 2); 

            auto tt_y_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 3); 

            auto tt_y_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 4); 

            auto tt_y_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 5); 

            auto tt_z_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 6); 

            auto tt_z_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 7); 

            auto tt_z_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 8); 

            auto tt_xx_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx); 

            auto tt_xy_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 1); 

            auto tt_xz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 2); 

            auto tt_yy_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 3); 

            auto tt_yz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 4); 

            auto tt_zz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 5); 

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

            auto tt_xxx_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx); 

            auto tt_xxx_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 1); 

            auto tt_xxx_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 2); 

            auto tt_xxy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 3); 

            auto tt_xxy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 4); 

            auto tt_xxy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 5); 

            auto tt_xxz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 6); 

            auto tt_xxz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 7); 

            auto tt_xxz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 8); 

            auto tt_xyy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 9); 

            auto tt_xyy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 10); 

            auto tt_xyy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 11); 

            auto tt_xyz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 12); 

            auto tt_xyz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 13); 

            auto tt_xyz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 14); 

            auto tt_xzz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 15); 

            auto tt_xzz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 16); 

            auto tt_xzz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 17); 

            auto tt_yyy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 18); 

            auto tt_yyy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 19); 

            auto tt_yyy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 20); 

            auto tt_yyz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 21); 

            auto tt_yyz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 22); 

            auto tt_yyz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 23); 

            auto tt_yzz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 24); 

            auto tt_yzz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 25); 

            auto tt_yzz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 26); 

            auto tt_zzz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 27); 

            auto tt_zzz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 28); 

            auto tt_zzz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 29); 

            #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, pa_z, ts_x_x_0, ts_x_y_0, ts_x_z_0, ts_xxx_x_0, \
                                     ts_xxx_y_0, ts_xxx_z_0, ts_xxy_x_0, ts_xxy_y_0, ts_xxy_z_0, ts_xxz_x_0, ts_xxz_y_0, \
                                     ts_xxz_z_0, ts_xyy_x_0, ts_xyy_y_0, ts_xyy_z_0, ts_xyz_x_0, ts_xyz_y_0, ts_xyz_z_0, \
                                     ts_xzz_x_0, ts_xzz_y_0, ts_xzz_z_0, ts_y_x_0, ts_y_y_0, ts_y_z_0, ts_yyy_x_0, \
                                     ts_yyy_y_0, ts_yyy_z_0, ts_yyz_x_0, ts_yyz_y_0, ts_yyz_z_0, ts_yzz_x_0, ts_yzz_y_0, \
                                     ts_yzz_z_0, ts_z_x_0, ts_z_y_0, ts_z_z_0, ts_zzz_x_0, ts_zzz_y_0, ts_zzz_z_0, \
                                     tt_x_x_0, tt_x_y_0, tt_x_z_0, tt_xx_0_0, tt_xx_x_0, tt_xx_y_0, tt_xx_z_0, \
                                     tt_xxx_x_0, tt_xxx_y_0, tt_xxx_z_0, tt_xxy_x_0, tt_xxy_y_0, tt_xxy_z_0, tt_xxz_x_0, \
                                     tt_xxz_y_0, tt_xxz_z_0, tt_xy_0_0, tt_xy_x_0, tt_xy_y_0, tt_xy_z_0, tt_xyy_x_0, \
                                     tt_xyy_y_0, tt_xyy_z_0, tt_xyz_x_0, tt_xyz_y_0, tt_xyz_z_0, tt_xz_0_0, tt_xz_x_0, \
                                     tt_xz_y_0, tt_xz_z_0, tt_xzz_x_0, tt_xzz_y_0, tt_xzz_z_0, tt_y_x_0, tt_y_y_0, \
                                     tt_y_z_0, tt_yy_0_0, tt_yy_x_0, tt_yy_y_0, tt_yy_z_0, tt_yyy_x_0, tt_yyy_y_0, \
                                     tt_yyy_z_0, tt_yyz_x_0, tt_yyz_y_0, tt_yyz_z_0, tt_yz_0_0, tt_yz_x_0, tt_yz_y_0, \
                                     tt_yz_z_0, tt_yzz_x_0, tt_yzz_y_0, tt_yzz_z_0, tt_z_x_0, tt_z_y_0, tt_z_z_0, \
                                     tt_zz_0_0, tt_zz_x_0, tt_zz_y_0, tt_zz_z_0, tt_zzz_x_0, tt_zzz_y_0, tt_zzz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fga = fga[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_xxx_x_0[j] = pa_x[j] * tt_xx_x_0[j] + fl1_fx * tt_x_x_0[j] + 0.5 * fl1_fx * tt_xx_0_0[j] + 2.0 * fl1_fz * ts_xxx_x_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_x_0[j];

                tt_xxx_y_0[j] = pa_x[j] * tt_xx_y_0[j] + fl1_fx * tt_x_y_0[j] + 2.0 * fl1_fz * ts_xxx_y_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_y_0[j];

                tt_xxx_z_0[j] = pa_x[j] * tt_xx_z_0[j] + fl1_fx * tt_x_z_0[j] + 2.0 * fl1_fz * ts_xxx_z_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_z_0[j];

                tt_xxy_x_0[j] = pa_x[j] * tt_xy_x_0[j] + 0.5 * fl1_fx * tt_y_x_0[j] + 0.5 * fl1_fx * tt_xy_0_0[j] + 2.0 * fl1_fz * ts_xxy_x_0[j] - fl1_fz * fl1_fga * ts_y_x_0[j];

                tt_xxy_y_0[j] = pa_x[j] * tt_xy_y_0[j] + 0.5 * fl1_fx * tt_y_y_0[j] + 2.0 * fl1_fz * ts_xxy_y_0[j] - fl1_fz * fl1_fga * ts_y_y_0[j];

                tt_xxy_z_0[j] = pa_x[j] * tt_xy_z_0[j] + 0.5 * fl1_fx * tt_y_z_0[j] + 2.0 * fl1_fz * ts_xxy_z_0[j] - fl1_fz * fl1_fga * ts_y_z_0[j];

                tt_xxz_x_0[j] = pa_x[j] * tt_xz_x_0[j] + 0.5 * fl1_fx * tt_z_x_0[j] + 0.5 * fl1_fx * tt_xz_0_0[j] + 2.0 * fl1_fz * ts_xxz_x_0[j] - fl1_fz * fl1_fga * ts_z_x_0[j];

                tt_xxz_y_0[j] = pa_x[j] * tt_xz_y_0[j] + 0.5 * fl1_fx * tt_z_y_0[j] + 2.0 * fl1_fz * ts_xxz_y_0[j] - fl1_fz * fl1_fga * ts_z_y_0[j];

                tt_xxz_z_0[j] = pa_x[j] * tt_xz_z_0[j] + 0.5 * fl1_fx * tt_z_z_0[j] + 2.0 * fl1_fz * ts_xxz_z_0[j] - fl1_fz * fl1_fga * ts_z_z_0[j];

                tt_xyy_x_0[j] = pa_x[j] * tt_yy_x_0[j] + 0.5 * fl1_fx * tt_yy_0_0[j] + 2.0 * fl1_fz * ts_xyy_x_0[j];

                tt_xyy_y_0[j] = pa_x[j] * tt_yy_y_0[j] + 2.0 * fl1_fz * ts_xyy_y_0[j];

                tt_xyy_z_0[j] = pa_x[j] * tt_yy_z_0[j] + 2.0 * fl1_fz * ts_xyy_z_0[j];

                tt_xyz_x_0[j] = pa_x[j] * tt_yz_x_0[j] + 0.5 * fl1_fx * tt_yz_0_0[j] + 2.0 * fl1_fz * ts_xyz_x_0[j];

                tt_xyz_y_0[j] = pa_x[j] * tt_yz_y_0[j] + 2.0 * fl1_fz * ts_xyz_y_0[j];

                tt_xyz_z_0[j] = pa_x[j] * tt_yz_z_0[j] + 2.0 * fl1_fz * ts_xyz_z_0[j];

                tt_xzz_x_0[j] = pa_x[j] * tt_zz_x_0[j] + 0.5 * fl1_fx * tt_zz_0_0[j] + 2.0 * fl1_fz * ts_xzz_x_0[j];

                tt_xzz_y_0[j] = pa_x[j] * tt_zz_y_0[j] + 2.0 * fl1_fz * ts_xzz_y_0[j];

                tt_xzz_z_0[j] = pa_x[j] * tt_zz_z_0[j] + 2.0 * fl1_fz * ts_xzz_z_0[j];

                tt_yyy_x_0[j] = pa_y[j] * tt_yy_x_0[j] + fl1_fx * tt_y_x_0[j] + 2.0 * fl1_fz * ts_yyy_x_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_x_0[j];

                tt_yyy_y_0[j] = pa_y[j] * tt_yy_y_0[j] + fl1_fx * tt_y_y_0[j] + 0.5 * fl1_fx * tt_yy_0_0[j] + 2.0 * fl1_fz * ts_yyy_y_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_y_0[j];

                tt_yyy_z_0[j] = pa_y[j] * tt_yy_z_0[j] + fl1_fx * tt_y_z_0[j] + 2.0 * fl1_fz * ts_yyy_z_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_z_0[j];

                tt_yyz_x_0[j] = pa_y[j] * tt_yz_x_0[j] + 0.5 * fl1_fx * tt_z_x_0[j] + 2.0 * fl1_fz * ts_yyz_x_0[j] - fl1_fz * fl1_fga * ts_z_x_0[j];

                tt_yyz_y_0[j] = pa_y[j] * tt_yz_y_0[j] + 0.5 * fl1_fx * tt_z_y_0[j] + 0.5 * fl1_fx * tt_yz_0_0[j] + 2.0 * fl1_fz * ts_yyz_y_0[j] - fl1_fz * fl1_fga * ts_z_y_0[j];

                tt_yyz_z_0[j] = pa_y[j] * tt_yz_z_0[j] + 0.5 * fl1_fx * tt_z_z_0[j] + 2.0 * fl1_fz * ts_yyz_z_0[j] - fl1_fz * fl1_fga * ts_z_z_0[j];

                tt_yzz_x_0[j] = pa_y[j] * tt_zz_x_0[j] + 2.0 * fl1_fz * ts_yzz_x_0[j];

                tt_yzz_y_0[j] = pa_y[j] * tt_zz_y_0[j] + 0.5 * fl1_fx * tt_zz_0_0[j] + 2.0 * fl1_fz * ts_yzz_y_0[j];

                tt_yzz_z_0[j] = pa_y[j] * tt_zz_z_0[j] + 2.0 * fl1_fz * ts_yzz_z_0[j];

                tt_zzz_x_0[j] = pa_z[j] * tt_zz_x_0[j] + fl1_fx * tt_z_x_0[j] + 2.0 * fl1_fz * ts_zzz_x_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_x_0[j];

                tt_zzz_y_0[j] = pa_z[j] * tt_zz_y_0[j] + fl1_fx * tt_z_y_0[j] + 2.0 * fl1_fz * ts_zzz_y_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_y_0[j];

                tt_zzz_z_0[j] = pa_z[j] * tt_zz_z_0[j] + fl1_fx * tt_z_z_0[j] + 0.5 * fl1_fx * tt_zz_0_0[j] + 2.0 * fl1_fz * ts_zzz_z_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_z_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForPG(      CMemBlock2D<double>& primBuffer,
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

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_t_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_1_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_0_xxxx_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx); 

            auto tt_0_xxxy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 1); 

            auto tt_0_xxxz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 2); 

            auto tt_0_xxyy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 3); 

            auto tt_0_xxyz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 4); 

            auto tt_0_xxzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 5); 

            auto tt_0_xyyy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 6); 

            auto tt_0_xyyz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 7); 

            auto tt_0_xyzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 8); 

            auto tt_0_xzzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 9); 

            auto tt_0_yyyy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 10); 

            auto tt_0_yyyz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 11); 

            auto tt_0_yyzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 12); 

            auto tt_0_yzzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 13); 

            auto tt_0_zzzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 14); 

            auto tt_0_xxx_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx); 

            auto tt_0_xxy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 1); 

            auto tt_0_xxz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 2); 

            auto tt_0_xyy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 3); 

            auto tt_0_xyz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 4); 

            auto tt_0_xzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 5); 

            auto tt_0_yyy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 6); 

            auto tt_0_yyz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 7); 

            auto tt_0_yzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 8); 

            auto tt_0_zzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 9); 

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

            auto tt_x_xxxx_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx); 

            auto tt_x_xxxy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 1); 

            auto tt_x_xxxz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 2); 

            auto tt_x_xxyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 3); 

            auto tt_x_xxyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 4); 

            auto tt_x_xxzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 5); 

            auto tt_x_xyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 6); 

            auto tt_x_xyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 7); 

            auto tt_x_xyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 8); 

            auto tt_x_xzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 9); 

            auto tt_x_yyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 10); 

            auto tt_x_yyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 11); 

            auto tt_x_yyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 12); 

            auto tt_x_yzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 13); 

            auto tt_x_zzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 14); 

            auto tt_y_xxxx_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 15); 

            auto tt_y_xxxy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 16); 

            auto tt_y_xxxz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 17); 

            auto tt_y_xxyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 18); 

            auto tt_y_xxyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 19); 

            auto tt_y_xxzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 20); 

            auto tt_y_xyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 21); 

            auto tt_y_xyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 22); 

            auto tt_y_xyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 23); 

            auto tt_y_xzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 24); 

            auto tt_y_yyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 25); 

            auto tt_y_yyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 26); 

            auto tt_y_yyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 27); 

            auto tt_y_yzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 28); 

            auto tt_y_zzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 29); 

            auto tt_z_xxxx_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 30); 

            auto tt_z_xxxy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 31); 

            auto tt_z_xxxz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 32); 

            auto tt_z_xxyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 33); 

            auto tt_z_xxyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 34); 

            auto tt_z_xxzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 35); 

            auto tt_z_xyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 36); 

            auto tt_z_xyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 37); 

            auto tt_z_xyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 38); 

            auto tt_z_xzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 39); 

            auto tt_z_yyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 40); 

            auto tt_z_yyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 41); 

            auto tt_z_yyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 42); 

            auto tt_z_yzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 43); 

            auto tt_z_zzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 44); 

            #pragma omp simd aligned(fx, fz, pa_x, pa_y, pa_z, ts_x_xxxx_0, ts_x_xxxy_0, ts_x_xxxz_0, ts_x_xxyy_0, \
                                     ts_x_xxyz_0, ts_x_xxzz_0, ts_x_xyyy_0, ts_x_xyyz_0, ts_x_xyzz_0, ts_x_xzzz_0, \
                                     ts_x_yyyy_0, ts_x_yyyz_0, ts_x_yyzz_0, ts_x_yzzz_0, ts_x_zzzz_0, ts_y_xxxx_0, \
                                     ts_y_xxxy_0, ts_y_xxxz_0, ts_y_xxyy_0, ts_y_xxyz_0, ts_y_xxzz_0, ts_y_xyyy_0, \
                                     ts_y_xyyz_0, ts_y_xyzz_0, ts_y_xzzz_0, ts_y_yyyy_0, ts_y_yyyz_0, ts_y_yyzz_0, \
                                     ts_y_yzzz_0, ts_y_zzzz_0, ts_z_xxxx_0, ts_z_xxxy_0, ts_z_xxxz_0, ts_z_xxyy_0, \
                                     ts_z_xxyz_0, ts_z_xxzz_0, ts_z_xyyy_0, ts_z_xyyz_0, ts_z_xyzz_0, ts_z_xzzz_0, \
                                     ts_z_yyyy_0, ts_z_yyyz_0, ts_z_yyzz_0, ts_z_yzzz_0, ts_z_zzzz_0, tt_0_xxx_0, \
                                     tt_0_xxxx_0, tt_0_xxxy_0, tt_0_xxxz_0, tt_0_xxy_0, tt_0_xxyy_0, tt_0_xxyz_0, \
                                     tt_0_xxz_0, tt_0_xxzz_0, tt_0_xyy_0, tt_0_xyyy_0, tt_0_xyyz_0, tt_0_xyz_0, \
                                     tt_0_xyzz_0, tt_0_xzz_0, tt_0_xzzz_0, tt_0_yyy_0, tt_0_yyyy_0, tt_0_yyyz_0, \
                                     tt_0_yyz_0, tt_0_yyzz_0, tt_0_yzz_0, tt_0_yzzz_0, tt_0_zzz_0, tt_0_zzzz_0, \
                                     tt_x_xxxx_0, tt_x_xxxy_0, tt_x_xxxz_0, tt_x_xxyy_0, tt_x_xxyz_0, tt_x_xxzz_0, \
                                     tt_x_xyyy_0, tt_x_xyyz_0, tt_x_xyzz_0, tt_x_xzzz_0, tt_x_yyyy_0, tt_x_yyyz_0, \
                                     tt_x_yyzz_0, tt_x_yzzz_0, tt_x_zzzz_0, tt_y_xxxx_0, tt_y_xxxy_0, tt_y_xxxz_0, \
                                     tt_y_xxyy_0, tt_y_xxyz_0, tt_y_xxzz_0, tt_y_xyyy_0, tt_y_xyyz_0, tt_y_xyzz_0, \
                                     tt_y_xzzz_0, tt_y_yyyy_0, tt_y_yyyz_0, tt_y_yyzz_0, tt_y_yzzz_0, tt_y_zzzz_0, \
                                     tt_z_xxxx_0, tt_z_xxxy_0, tt_z_xxxz_0, tt_z_xxyy_0, tt_z_xxyz_0, tt_z_xxzz_0, \
                                     tt_z_xyyy_0, tt_z_xyyz_0, tt_z_xyzz_0, tt_z_xzzz_0, tt_z_yyyy_0, tt_z_yyyz_0, \
                                     tt_z_yyzz_0, tt_z_yzzz_0, tt_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_x_xxxx_0[j] = pa_x[j] * tt_0_xxxx_0[j] + 2.0 * fl1_fx * tt_0_xxx_0[j] + 2.0 * fl1_fz * ts_x_xxxx_0[j];

                tt_x_xxxy_0[j] = pa_x[j] * tt_0_xxxy_0[j] + 1.5 * fl1_fx * tt_0_xxy_0[j] + 2.0 * fl1_fz * ts_x_xxxy_0[j];

                tt_x_xxxz_0[j] = pa_x[j] * tt_0_xxxz_0[j] + 1.5 * fl1_fx * tt_0_xxz_0[j] + 2.0 * fl1_fz * ts_x_xxxz_0[j];

                tt_x_xxyy_0[j] = pa_x[j] * tt_0_xxyy_0[j] + fl1_fx * tt_0_xyy_0[j] + 2.0 * fl1_fz * ts_x_xxyy_0[j];

                tt_x_xxyz_0[j] = pa_x[j] * tt_0_xxyz_0[j] + fl1_fx * tt_0_xyz_0[j] + 2.0 * fl1_fz * ts_x_xxyz_0[j];

                tt_x_xxzz_0[j] = pa_x[j] * tt_0_xxzz_0[j] + fl1_fx * tt_0_xzz_0[j] + 2.0 * fl1_fz * ts_x_xxzz_0[j];

                tt_x_xyyy_0[j] = pa_x[j] * tt_0_xyyy_0[j] + 0.5 * fl1_fx * tt_0_yyy_0[j] + 2.0 * fl1_fz * ts_x_xyyy_0[j];

                tt_x_xyyz_0[j] = pa_x[j] * tt_0_xyyz_0[j] + 0.5 * fl1_fx * tt_0_yyz_0[j] + 2.0 * fl1_fz * ts_x_xyyz_0[j];

                tt_x_xyzz_0[j] = pa_x[j] * tt_0_xyzz_0[j] + 0.5 * fl1_fx * tt_0_yzz_0[j] + 2.0 * fl1_fz * ts_x_xyzz_0[j];

                tt_x_xzzz_0[j] = pa_x[j] * tt_0_xzzz_0[j] + 0.5 * fl1_fx * tt_0_zzz_0[j] + 2.0 * fl1_fz * ts_x_xzzz_0[j];

                tt_x_yyyy_0[j] = pa_x[j] * tt_0_yyyy_0[j] + 2.0 * fl1_fz * ts_x_yyyy_0[j];

                tt_x_yyyz_0[j] = pa_x[j] * tt_0_yyyz_0[j] + 2.0 * fl1_fz * ts_x_yyyz_0[j];

                tt_x_yyzz_0[j] = pa_x[j] * tt_0_yyzz_0[j] + 2.0 * fl1_fz * ts_x_yyzz_0[j];

                tt_x_yzzz_0[j] = pa_x[j] * tt_0_yzzz_0[j] + 2.0 * fl1_fz * ts_x_yzzz_0[j];

                tt_x_zzzz_0[j] = pa_x[j] * tt_0_zzzz_0[j] + 2.0 * fl1_fz * ts_x_zzzz_0[j];

                tt_y_xxxx_0[j] = pa_y[j] * tt_0_xxxx_0[j] + 2.0 * fl1_fz * ts_y_xxxx_0[j];

                tt_y_xxxy_0[j] = pa_y[j] * tt_0_xxxy_0[j] + 0.5 * fl1_fx * tt_0_xxx_0[j] + 2.0 * fl1_fz * ts_y_xxxy_0[j];

                tt_y_xxxz_0[j] = pa_y[j] * tt_0_xxxz_0[j] + 2.0 * fl1_fz * ts_y_xxxz_0[j];

                tt_y_xxyy_0[j] = pa_y[j] * tt_0_xxyy_0[j] + fl1_fx * tt_0_xxy_0[j] + 2.0 * fl1_fz * ts_y_xxyy_0[j];

                tt_y_xxyz_0[j] = pa_y[j] * tt_0_xxyz_0[j] + 0.5 * fl1_fx * tt_0_xxz_0[j] + 2.0 * fl1_fz * ts_y_xxyz_0[j];

                tt_y_xxzz_0[j] = pa_y[j] * tt_0_xxzz_0[j] + 2.0 * fl1_fz * ts_y_xxzz_0[j];

                tt_y_xyyy_0[j] = pa_y[j] * tt_0_xyyy_0[j] + 1.5 * fl1_fx * tt_0_xyy_0[j] + 2.0 * fl1_fz * ts_y_xyyy_0[j];

                tt_y_xyyz_0[j] = pa_y[j] * tt_0_xyyz_0[j] + fl1_fx * tt_0_xyz_0[j] + 2.0 * fl1_fz * ts_y_xyyz_0[j];

                tt_y_xyzz_0[j] = pa_y[j] * tt_0_xyzz_0[j] + 0.5 * fl1_fx * tt_0_xzz_0[j] + 2.0 * fl1_fz * ts_y_xyzz_0[j];

                tt_y_xzzz_0[j] = pa_y[j] * tt_0_xzzz_0[j] + 2.0 * fl1_fz * ts_y_xzzz_0[j];

                tt_y_yyyy_0[j] = pa_y[j] * tt_0_yyyy_0[j] + 2.0 * fl1_fx * tt_0_yyy_0[j] + 2.0 * fl1_fz * ts_y_yyyy_0[j];

                tt_y_yyyz_0[j] = pa_y[j] * tt_0_yyyz_0[j] + 1.5 * fl1_fx * tt_0_yyz_0[j] + 2.0 * fl1_fz * ts_y_yyyz_0[j];

                tt_y_yyzz_0[j] = pa_y[j] * tt_0_yyzz_0[j] + fl1_fx * tt_0_yzz_0[j] + 2.0 * fl1_fz * ts_y_yyzz_0[j];

                tt_y_yzzz_0[j] = pa_y[j] * tt_0_yzzz_0[j] + 0.5 * fl1_fx * tt_0_zzz_0[j] + 2.0 * fl1_fz * ts_y_yzzz_0[j];

                tt_y_zzzz_0[j] = pa_y[j] * tt_0_zzzz_0[j] + 2.0 * fl1_fz * ts_y_zzzz_0[j];

                tt_z_xxxx_0[j] = pa_z[j] * tt_0_xxxx_0[j] + 2.0 * fl1_fz * ts_z_xxxx_0[j];

                tt_z_xxxy_0[j] = pa_z[j] * tt_0_xxxy_0[j] + 2.0 * fl1_fz * ts_z_xxxy_0[j];

                tt_z_xxxz_0[j] = pa_z[j] * tt_0_xxxz_0[j] + 0.5 * fl1_fx * tt_0_xxx_0[j] + 2.0 * fl1_fz * ts_z_xxxz_0[j];

                tt_z_xxyy_0[j] = pa_z[j] * tt_0_xxyy_0[j] + 2.0 * fl1_fz * ts_z_xxyy_0[j];

                tt_z_xxyz_0[j] = pa_z[j] * tt_0_xxyz_0[j] + 0.5 * fl1_fx * tt_0_xxy_0[j] + 2.0 * fl1_fz * ts_z_xxyz_0[j];

                tt_z_xxzz_0[j] = pa_z[j] * tt_0_xxzz_0[j] + fl1_fx * tt_0_xxz_0[j] + 2.0 * fl1_fz * ts_z_xxzz_0[j];

                tt_z_xyyy_0[j] = pa_z[j] * tt_0_xyyy_0[j] + 2.0 * fl1_fz * ts_z_xyyy_0[j];

                tt_z_xyyz_0[j] = pa_z[j] * tt_0_xyyz_0[j] + 0.5 * fl1_fx * tt_0_xyy_0[j] + 2.0 * fl1_fz * ts_z_xyyz_0[j];

                tt_z_xyzz_0[j] = pa_z[j] * tt_0_xyzz_0[j] + fl1_fx * tt_0_xyz_0[j] + 2.0 * fl1_fz * ts_z_xyzz_0[j];

                tt_z_xzzz_0[j] = pa_z[j] * tt_0_xzzz_0[j] + 1.5 * fl1_fx * tt_0_xzz_0[j] + 2.0 * fl1_fz * ts_z_xzzz_0[j];

                tt_z_yyyy_0[j] = pa_z[j] * tt_0_yyyy_0[j] + 2.0 * fl1_fz * ts_z_yyyy_0[j];

                tt_z_yyyz_0[j] = pa_z[j] * tt_0_yyyz_0[j] + 0.5 * fl1_fx * tt_0_yyy_0[j] + 2.0 * fl1_fz * ts_z_yyyz_0[j];

                tt_z_yyzz_0[j] = pa_z[j] * tt_0_yyzz_0[j] + fl1_fx * tt_0_yyz_0[j] + 2.0 * fl1_fz * ts_z_yyzz_0[j];

                tt_z_yzzz_0[j] = pa_z[j] * tt_0_yzzz_0[j] + 1.5 * fl1_fx * tt_0_yzz_0[j] + 2.0 * fl1_fz * ts_z_yzzz_0[j];

                tt_z_zzzz_0[j] = pa_z[j] * tt_0_zzzz_0[j] + 2.0 * fl1_fx * tt_0_zzz_0[j] + 2.0 * fl1_fz * ts_z_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGP(      CMemBlock2D<double>& primBuffer,
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

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up index of integral

        auto pidx_t_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_4_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_xxx_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx); 

            auto tt_xxx_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 1); 

            auto tt_xxx_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 2); 

            auto tt_xxy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 3); 

            auto tt_xxy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 4); 

            auto tt_xxy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 5); 

            auto tt_xxz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 6); 

            auto tt_xxz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 7); 

            auto tt_xxz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 8); 

            auto tt_xyy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 9); 

            auto tt_xyy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 10); 

            auto tt_xyy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 11); 

            auto tt_xyz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 12); 

            auto tt_xyz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 13); 

            auto tt_xyz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 14); 

            auto tt_xzz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 15); 

            auto tt_xzz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 16); 

            auto tt_xzz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 17); 

            auto tt_yyy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 18); 

            auto tt_yyy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 19); 

            auto tt_yyy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 20); 

            auto tt_yyz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 21); 

            auto tt_yyz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 22); 

            auto tt_yyz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 23); 

            auto tt_yzz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 24); 

            auto tt_yzz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 25); 

            auto tt_yzz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 26); 

            auto tt_zzz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 27); 

            auto tt_zzz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 28); 

            auto tt_zzz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 29); 

            auto tt_xx_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx); 

            auto tt_xx_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 1); 

            auto tt_xx_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 2); 

            auto tt_xy_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 3); 

            auto tt_xy_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 4); 

            auto tt_xy_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 5); 

            auto tt_xz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 6); 

            auto tt_xz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 7); 

            auto tt_xz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 8); 

            auto tt_yy_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 9); 

            auto tt_yy_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 10); 

            auto tt_yy_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 11); 

            auto tt_yz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 12); 

            auto tt_yz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 13); 

            auto tt_yz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 14); 

            auto tt_zz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 15); 

            auto tt_zz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 16); 

            auto tt_zz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 17); 

            auto tt_xxx_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx); 

            auto tt_xxy_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx + 1); 

            auto tt_xxz_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx + 2); 

            auto tt_xyy_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx + 3); 

            auto tt_xyz_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx + 4); 

            auto tt_xzz_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx + 5); 

            auto tt_yyy_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx + 6); 

            auto tt_yyz_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx + 7); 

            auto tt_yzz_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx + 8); 

            auto tt_zzz_0_0 = primBuffer.data(pidx_t_3_0_m0 + 10 * idx + 9); 

            auto ts_xxxx_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx); 

            auto ts_xxxx_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 1); 

            auto ts_xxxx_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 2); 

            auto ts_xxxy_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 3); 

            auto ts_xxxy_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 4); 

            auto ts_xxxy_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 5); 

            auto ts_xxxz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 6); 

            auto ts_xxxz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 7); 

            auto ts_xxxz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 8); 

            auto ts_xxyy_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 9); 

            auto ts_xxyy_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 10); 

            auto ts_xxyy_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 11); 

            auto ts_xxyz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 12); 

            auto ts_xxyz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 13); 

            auto ts_xxyz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 14); 

            auto ts_xxzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 15); 

            auto ts_xxzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 16); 

            auto ts_xxzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 17); 

            auto ts_xyyy_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 18); 

            auto ts_xyyy_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 19); 

            auto ts_xyyy_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 20); 

            auto ts_xyyz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 21); 

            auto ts_xyyz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 22); 

            auto ts_xyyz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 23); 

            auto ts_xyzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 24); 

            auto ts_xyzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 25); 

            auto ts_xyzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 26); 

            auto ts_xzzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 27); 

            auto ts_xzzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 28); 

            auto ts_xzzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 29); 

            auto ts_yyyy_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 30); 

            auto ts_yyyy_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 31); 

            auto ts_yyyy_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 32); 

            auto ts_yyyz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 33); 

            auto ts_yyyz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 34); 

            auto ts_yyyz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 35); 

            auto ts_yyzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 36); 

            auto ts_yyzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 37); 

            auto ts_yyzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 38); 

            auto ts_yzzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 39); 

            auto ts_yzzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 40); 

            auto ts_yzzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 41); 

            auto ts_zzzz_x_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 42); 

            auto ts_zzzz_y_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 43); 

            auto ts_zzzz_z_0 = primBuffer.data(pidx_s_4_1_m0 + 45 * idx + 44); 

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

            auto ts_zz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 15); 

            auto ts_zz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 16); 

            auto ts_zz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 17); 

            // set up pointers to integrals

            auto tt_xxxx_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx); 

            auto tt_xxxx_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 1); 

            auto tt_xxxx_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 2); 

            auto tt_xxxy_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 3); 

            auto tt_xxxy_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 4); 

            auto tt_xxxy_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 5); 

            auto tt_xxxz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 6); 

            auto tt_xxxz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 7); 

            auto tt_xxxz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 8); 

            auto tt_xxyy_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 9); 

            auto tt_xxyy_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 10); 

            auto tt_xxyy_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 11); 

            auto tt_xxyz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 12); 

            auto tt_xxyz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 13); 

            auto tt_xxyz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 14); 

            auto tt_xxzz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 15); 

            auto tt_xxzz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 16); 

            auto tt_xxzz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 17); 

            auto tt_xyyy_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 18); 

            auto tt_xyyy_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 19); 

            auto tt_xyyy_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 20); 

            auto tt_xyyz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 21); 

            auto tt_xyyz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 22); 

            auto tt_xyyz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 23); 

            auto tt_xyzz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 24); 

            auto tt_xyzz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 25); 

            auto tt_xyzz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 26); 

            auto tt_xzzz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 27); 

            auto tt_xzzz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 28); 

            auto tt_xzzz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 29); 

            auto tt_yyyy_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 30); 

            auto tt_yyyy_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 31); 

            auto tt_yyyy_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 32); 

            auto tt_yyyz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 33); 

            auto tt_yyyz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 34); 

            auto tt_yyyz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 35); 

            auto tt_yyzz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 36); 

            auto tt_yyzz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 37); 

            auto tt_yyzz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 38); 

            auto tt_yzzz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 39); 

            auto tt_yzzz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 40); 

            auto tt_yzzz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 41); 

            auto tt_zzzz_x_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 42); 

            auto tt_zzzz_y_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 43); 

            auto tt_zzzz_z_0 = primBuffer.data(pidx_t_4_1_m0 + 45 * idx + 44); 

            #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, pa_z, ts_xx_x_0, ts_xx_y_0, ts_xx_z_0, ts_xxxx_x_0, \
                                     ts_xxxx_y_0, ts_xxxx_z_0, ts_xxxy_x_0, ts_xxxy_y_0, ts_xxxy_z_0, ts_xxxz_x_0, \
                                     ts_xxxz_y_0, ts_xxxz_z_0, ts_xxyy_x_0, ts_xxyy_y_0, ts_xxyy_z_0, ts_xxyz_x_0, \
                                     ts_xxyz_y_0, ts_xxyz_z_0, ts_xxzz_x_0, ts_xxzz_y_0, ts_xxzz_z_0, ts_xy_x_0, \
                                     ts_xy_y_0, ts_xy_z_0, ts_xyyy_x_0, ts_xyyy_y_0, ts_xyyy_z_0, ts_xyyz_x_0, \
                                     ts_xyyz_y_0, ts_xyyz_z_0, ts_xyzz_x_0, ts_xyzz_y_0, ts_xyzz_z_0, ts_xz_x_0, \
                                     ts_xz_y_0, ts_xz_z_0, ts_xzzz_x_0, ts_xzzz_y_0, ts_xzzz_z_0, ts_yy_x_0, ts_yy_y_0, \
                                     ts_yy_z_0, ts_yyyy_x_0, ts_yyyy_y_0, ts_yyyy_z_0, ts_yyyz_x_0, ts_yyyz_y_0, \
                                     ts_yyyz_z_0, ts_yyzz_x_0, ts_yyzz_y_0, ts_yyzz_z_0, ts_yz_x_0, ts_yz_y_0, ts_yz_z_0, \
                                     ts_yzzz_x_0, ts_yzzz_y_0, ts_yzzz_z_0, ts_zz_x_0, ts_zz_y_0, ts_zz_z_0, ts_zzzz_x_0, \
                                     ts_zzzz_y_0, ts_zzzz_z_0, tt_xx_x_0, tt_xx_y_0, tt_xx_z_0, tt_xxx_0_0, tt_xxx_x_0, \
                                     tt_xxx_y_0, tt_xxx_z_0, tt_xxxx_x_0, tt_xxxx_y_0, tt_xxxx_z_0, tt_xxxy_x_0, \
                                     tt_xxxy_y_0, tt_xxxy_z_0, tt_xxxz_x_0, tt_xxxz_y_0, tt_xxxz_z_0, tt_xxy_0_0, \
                                     tt_xxy_x_0, tt_xxy_y_0, tt_xxy_z_0, tt_xxyy_x_0, tt_xxyy_y_0, tt_xxyy_z_0, \
                                     tt_xxyz_x_0, tt_xxyz_y_0, tt_xxyz_z_0, tt_xxz_0_0, tt_xxz_x_0, tt_xxz_y_0, \
                                     tt_xxz_z_0, tt_xxzz_x_0, tt_xxzz_y_0, tt_xxzz_z_0, tt_xy_x_0, tt_xy_y_0, tt_xy_z_0, \
                                     tt_xyy_0_0, tt_xyy_x_0, tt_xyy_y_0, tt_xyy_z_0, tt_xyyy_x_0, tt_xyyy_y_0, \
                                     tt_xyyy_z_0, tt_xyyz_x_0, tt_xyyz_y_0, tt_xyyz_z_0, tt_xyz_0_0, tt_xyz_x_0, \
                                     tt_xyz_y_0, tt_xyz_z_0, tt_xyzz_x_0, tt_xyzz_y_0, tt_xyzz_z_0, tt_xz_x_0, \
                                     tt_xz_y_0, tt_xz_z_0, tt_xzz_0_0, tt_xzz_x_0, tt_xzz_y_0, tt_xzz_z_0, tt_xzzz_x_0, \
                                     tt_xzzz_y_0, tt_xzzz_z_0, tt_yy_x_0, tt_yy_y_0, tt_yy_z_0, tt_yyy_0_0, tt_yyy_x_0, \
                                     tt_yyy_y_0, tt_yyy_z_0, tt_yyyy_x_0, tt_yyyy_y_0, tt_yyyy_z_0, tt_yyyz_x_0, \
                                     tt_yyyz_y_0, tt_yyyz_z_0, tt_yyz_0_0, tt_yyz_x_0, tt_yyz_y_0, tt_yyz_z_0, \
                                     tt_yyzz_x_0, tt_yyzz_y_0, tt_yyzz_z_0, tt_yz_x_0, tt_yz_y_0, tt_yz_z_0, tt_yzz_0_0, \
                                     tt_yzz_x_0, tt_yzz_y_0, tt_yzz_z_0, tt_yzzz_x_0, tt_yzzz_y_0, tt_yzzz_z_0, \
                                     tt_zz_x_0, tt_zz_y_0, tt_zz_z_0, tt_zzz_0_0, tt_zzz_x_0, tt_zzz_y_0, tt_zzz_z_0, \
                                     tt_zzzz_x_0, tt_zzzz_y_0, tt_zzzz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fga = fga[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_xxxx_x_0[j] = pa_x[j] * tt_xxx_x_0[j] + 1.5 * fl1_fx * tt_xx_x_0[j] + 0.5 * fl1_fx * tt_xxx_0_0[j] + 2.0 * fl1_fz * ts_xxxx_x_0[j] - 3.0 * fl1_fz * fl1_fga * ts_xx_x_0[j];

                tt_xxxx_y_0[j] = pa_x[j] * tt_xxx_y_0[j] + 1.5 * fl1_fx * tt_xx_y_0[j] + 2.0 * fl1_fz * ts_xxxx_y_0[j] - 3.0 * fl1_fz * fl1_fga * ts_xx_y_0[j];

                tt_xxxx_z_0[j] = pa_x[j] * tt_xxx_z_0[j] + 1.5 * fl1_fx * tt_xx_z_0[j] + 2.0 * fl1_fz * ts_xxxx_z_0[j] - 3.0 * fl1_fz * fl1_fga * ts_xx_z_0[j];

                tt_xxxy_x_0[j] = pa_x[j] * tt_xxy_x_0[j] + fl1_fx * tt_xy_x_0[j] + 0.5 * fl1_fx * tt_xxy_0_0[j] + 2.0 * fl1_fz * ts_xxxy_x_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xy_x_0[j];

                tt_xxxy_y_0[j] = pa_x[j] * tt_xxy_y_0[j] + fl1_fx * tt_xy_y_0[j] + 2.0 * fl1_fz * ts_xxxy_y_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xy_y_0[j];

                tt_xxxy_z_0[j] = pa_x[j] * tt_xxy_z_0[j] + fl1_fx * tt_xy_z_0[j] + 2.0 * fl1_fz * ts_xxxy_z_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xy_z_0[j];

                tt_xxxz_x_0[j] = pa_x[j] * tt_xxz_x_0[j] + fl1_fx * tt_xz_x_0[j] + 0.5 * fl1_fx * tt_xxz_0_0[j] + 2.0 * fl1_fz * ts_xxxz_x_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xz_x_0[j];

                tt_xxxz_y_0[j] = pa_x[j] * tt_xxz_y_0[j] + fl1_fx * tt_xz_y_0[j] + 2.0 * fl1_fz * ts_xxxz_y_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xz_y_0[j];

                tt_xxxz_z_0[j] = pa_x[j] * tt_xxz_z_0[j] + fl1_fx * tt_xz_z_0[j] + 2.0 * fl1_fz * ts_xxxz_z_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xz_z_0[j];

                tt_xxyy_x_0[j] = pa_x[j] * tt_xyy_x_0[j] + 0.5 * fl1_fx * tt_yy_x_0[j] + 0.5 * fl1_fx * tt_xyy_0_0[j] + 2.0 * fl1_fz * ts_xxyy_x_0[j] - fl1_fz * fl1_fga * ts_yy_x_0[j];

                tt_xxyy_y_0[j] = pa_x[j] * tt_xyy_y_0[j] + 0.5 * fl1_fx * tt_yy_y_0[j] + 2.0 * fl1_fz * ts_xxyy_y_0[j] - fl1_fz * fl1_fga * ts_yy_y_0[j];

                tt_xxyy_z_0[j] = pa_x[j] * tt_xyy_z_0[j] + 0.5 * fl1_fx * tt_yy_z_0[j] + 2.0 * fl1_fz * ts_xxyy_z_0[j] - fl1_fz * fl1_fga * ts_yy_z_0[j];

                tt_xxyz_x_0[j] = pa_x[j] * tt_xyz_x_0[j] + 0.5 * fl1_fx * tt_yz_x_0[j] + 0.5 * fl1_fx * tt_xyz_0_0[j] + 2.0 * fl1_fz * ts_xxyz_x_0[j] - fl1_fz * fl1_fga * ts_yz_x_0[j];

                tt_xxyz_y_0[j] = pa_x[j] * tt_xyz_y_0[j] + 0.5 * fl1_fx * tt_yz_y_0[j] + 2.0 * fl1_fz * ts_xxyz_y_0[j] - fl1_fz * fl1_fga * ts_yz_y_0[j];

                tt_xxyz_z_0[j] = pa_x[j] * tt_xyz_z_0[j] + 0.5 * fl1_fx * tt_yz_z_0[j] + 2.0 * fl1_fz * ts_xxyz_z_0[j] - fl1_fz * fl1_fga * ts_yz_z_0[j];

                tt_xxzz_x_0[j] = pa_x[j] * tt_xzz_x_0[j] + 0.5 * fl1_fx * tt_zz_x_0[j] + 0.5 * fl1_fx * tt_xzz_0_0[j] + 2.0 * fl1_fz * ts_xxzz_x_0[j] - fl1_fz * fl1_fga * ts_zz_x_0[j];

                tt_xxzz_y_0[j] = pa_x[j] * tt_xzz_y_0[j] + 0.5 * fl1_fx * tt_zz_y_0[j] + 2.0 * fl1_fz * ts_xxzz_y_0[j] - fl1_fz * fl1_fga * ts_zz_y_0[j];

                tt_xxzz_z_0[j] = pa_x[j] * tt_xzz_z_0[j] + 0.5 * fl1_fx * tt_zz_z_0[j] + 2.0 * fl1_fz * ts_xxzz_z_0[j] - fl1_fz * fl1_fga * ts_zz_z_0[j];

                tt_xyyy_x_0[j] = pa_x[j] * tt_yyy_x_0[j] + 0.5 * fl1_fx * tt_yyy_0_0[j] + 2.0 * fl1_fz * ts_xyyy_x_0[j];

                tt_xyyy_y_0[j] = pa_x[j] * tt_yyy_y_0[j] + 2.0 * fl1_fz * ts_xyyy_y_0[j];

                tt_xyyy_z_0[j] = pa_x[j] * tt_yyy_z_0[j] + 2.0 * fl1_fz * ts_xyyy_z_0[j];

                tt_xyyz_x_0[j] = pa_x[j] * tt_yyz_x_0[j] + 0.5 * fl1_fx * tt_yyz_0_0[j] + 2.0 * fl1_fz * ts_xyyz_x_0[j];

                tt_xyyz_y_0[j] = pa_x[j] * tt_yyz_y_0[j] + 2.0 * fl1_fz * ts_xyyz_y_0[j];

                tt_xyyz_z_0[j] = pa_x[j] * tt_yyz_z_0[j] + 2.0 * fl1_fz * ts_xyyz_z_0[j];

                tt_xyzz_x_0[j] = pa_x[j] * tt_yzz_x_0[j] + 0.5 * fl1_fx * tt_yzz_0_0[j] + 2.0 * fl1_fz * ts_xyzz_x_0[j];

                tt_xyzz_y_0[j] = pa_x[j] * tt_yzz_y_0[j] + 2.0 * fl1_fz * ts_xyzz_y_0[j];

                tt_xyzz_z_0[j] = pa_x[j] * tt_yzz_z_0[j] + 2.0 * fl1_fz * ts_xyzz_z_0[j];

                tt_xzzz_x_0[j] = pa_x[j] * tt_zzz_x_0[j] + 0.5 * fl1_fx * tt_zzz_0_0[j] + 2.0 * fl1_fz * ts_xzzz_x_0[j];

                tt_xzzz_y_0[j] = pa_x[j] * tt_zzz_y_0[j] + 2.0 * fl1_fz * ts_xzzz_y_0[j];

                tt_xzzz_z_0[j] = pa_x[j] * tt_zzz_z_0[j] + 2.0 * fl1_fz * ts_xzzz_z_0[j];

                tt_yyyy_x_0[j] = pa_y[j] * tt_yyy_x_0[j] + 1.5 * fl1_fx * tt_yy_x_0[j] + 2.0 * fl1_fz * ts_yyyy_x_0[j] - 3.0 * fl1_fz * fl1_fga * ts_yy_x_0[j];

                tt_yyyy_y_0[j] = pa_y[j] * tt_yyy_y_0[j] + 1.5 * fl1_fx * tt_yy_y_0[j] + 0.5 * fl1_fx * tt_yyy_0_0[j] + 2.0 * fl1_fz * ts_yyyy_y_0[j] - 3.0 * fl1_fz * fl1_fga * ts_yy_y_0[j];

                tt_yyyy_z_0[j] = pa_y[j] * tt_yyy_z_0[j] + 1.5 * fl1_fx * tt_yy_z_0[j] + 2.0 * fl1_fz * ts_yyyy_z_0[j] - 3.0 * fl1_fz * fl1_fga * ts_yy_z_0[j];

                tt_yyyz_x_0[j] = pa_y[j] * tt_yyz_x_0[j] + fl1_fx * tt_yz_x_0[j] + 2.0 * fl1_fz * ts_yyyz_x_0[j] - 2.0 * fl1_fz * fl1_fga * ts_yz_x_0[j];

                tt_yyyz_y_0[j] = pa_y[j] * tt_yyz_y_0[j] + fl1_fx * tt_yz_y_0[j] + 0.5 * fl1_fx * tt_yyz_0_0[j] + 2.0 * fl1_fz * ts_yyyz_y_0[j] - 2.0 * fl1_fz * fl1_fga * ts_yz_y_0[j];

                tt_yyyz_z_0[j] = pa_y[j] * tt_yyz_z_0[j] + fl1_fx * tt_yz_z_0[j] + 2.0 * fl1_fz * ts_yyyz_z_0[j] - 2.0 * fl1_fz * fl1_fga * ts_yz_z_0[j];

                tt_yyzz_x_0[j] = pa_y[j] * tt_yzz_x_0[j] + 0.5 * fl1_fx * tt_zz_x_0[j] + 2.0 * fl1_fz * ts_yyzz_x_0[j] - fl1_fz * fl1_fga * ts_zz_x_0[j];

                tt_yyzz_y_0[j] = pa_y[j] * tt_yzz_y_0[j] + 0.5 * fl1_fx * tt_zz_y_0[j] + 0.5 * fl1_fx * tt_yzz_0_0[j] + 2.0 * fl1_fz * ts_yyzz_y_0[j] - fl1_fz * fl1_fga * ts_zz_y_0[j];

                tt_yyzz_z_0[j] = pa_y[j] * tt_yzz_z_0[j] + 0.5 * fl1_fx * tt_zz_z_0[j] + 2.0 * fl1_fz * ts_yyzz_z_0[j] - fl1_fz * fl1_fga * ts_zz_z_0[j];

                tt_yzzz_x_0[j] = pa_y[j] * tt_zzz_x_0[j] + 2.0 * fl1_fz * ts_yzzz_x_0[j];

                tt_yzzz_y_0[j] = pa_y[j] * tt_zzz_y_0[j] + 0.5 * fl1_fx * tt_zzz_0_0[j] + 2.0 * fl1_fz * ts_yzzz_y_0[j];

                tt_yzzz_z_0[j] = pa_y[j] * tt_zzz_z_0[j] + 2.0 * fl1_fz * ts_yzzz_z_0[j];

                tt_zzzz_x_0[j] = pa_z[j] * tt_zzz_x_0[j] + 1.5 * fl1_fx * tt_zz_x_0[j] + 2.0 * fl1_fz * ts_zzzz_x_0[j] - 3.0 * fl1_fz * fl1_fga * ts_zz_x_0[j];

                tt_zzzz_y_0[j] = pa_z[j] * tt_zzz_y_0[j] + 1.5 * fl1_fx * tt_zz_y_0[j] + 2.0 * fl1_fz * ts_zzzz_y_0[j] - 3.0 * fl1_fz * fl1_fga * ts_zz_y_0[j];

                tt_zzzz_z_0[j] = pa_z[j] * tt_zzz_z_0[j] + 1.5 * fl1_fx * tt_zz_z_0[j] + 0.5 * fl1_fx * tt_zzz_0_0[j] + 2.0 * fl1_fz * ts_zzzz_z_0[j] - 3.0 * fl1_fz * fl1_fga * ts_zz_z_0[j];
            }

            idx++;
        }
    }


} // kinrecfunc namespace

