//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForSX.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForSS(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& abDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // set up pointers to primitives data on bra side
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // set up pointers to R(AB) distances
        
        auto abx = abDistances.data(0);
        
        auto aby = abDistances.data(1);
        
        auto abz = abDistances.data(2);
        
        // set up pointer to overlap integrals
        
        auto sidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                               {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));
        
        // set up pointer to kinetic energy integrals
        
        auto tidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true,
                                                               {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fz = osFactors.data(4 * idx + 1);
            
            // set up primitive integrals data
            
            auto s_0_0 = primBuffer.data(sidx + idx);
            
            auto t_0_0 = primBuffer.data(tidx + idx);
            
            #pragma omp simd aligned(s_0_0, t_0_0, fz, abx, aby, abz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double r2ab = abx[j] * abx[j] + aby[j] * aby[j] + abz[j] * abz[j];
                
                t_0_0[j] = fz[j] * (3.0 - 2.0 * fz[j] * r2ab) * s_0_0[j];
            }
            
            idx++;
        }

    }

    void
    compKineticEnergyForSP(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& pbDistances,
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

        auto pidx_t_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_0_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
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

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_0_0_0 = primBuffer.data(pidx_t_0_0_m0 + idx); 

            auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx); 

            auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1); 

            auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2); 

            // set up pointers to integrals

            auto tt_0_x_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx); 

            auto tt_0_y_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 1); 

            auto tt_0_z_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 2); 

            #pragma omp simd aligned(fz, pb_x, pb_y, pb_z, ts_0_x_0, ts_0_y_0, ts_0_z_0, tt_0_0_0, tt_0_x_0, \
                                     tt_0_y_0, tt_0_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fz = fz[j];

                tt_0_x_0[j] = pb_x[j] * tt_0_0_0[j] + 2.0 * fl1_fz * ts_0_x_0[j];

                tt_0_y_0[j] = pb_y[j] * tt_0_0_0[j] + 2.0 * fl1_fz * ts_0_y_0[j];

                tt_0_z_0[j] = pb_z[j] * tt_0_0_0[j] + 2.0 * fl1_fz * ts_0_z_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForPS(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_t_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_1_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fz = osFactors.data(4 * idx + 1);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_0_0_0 = primBuffer.data(pidx_t_0_0_m0 + idx); 

            auto ts_x_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx); 

            auto ts_y_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 1); 

            auto ts_z_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 2); 

            // set up pointers to integrals

            auto tt_x_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx); 

            auto tt_y_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx + 1); 

            auto tt_z_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx + 2); 

            #pragma omp simd aligned(fz, pa_x, pa_y, pa_z, ts_x_0_0, ts_y_0_0, ts_z_0_0, tt_0_0_0, tt_x_0_0, \
                                     tt_y_0_0, tt_z_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fz = fz[j];

                tt_x_0_0[j] = pa_x[j] * tt_0_0_0[j] + 2.0 * fl1_fz * ts_x_0_0[j];

                tt_y_0_0[j] = pa_y[j] * tt_0_0_0[j] + 2.0 * fl1_fz * ts_y_0_0[j];

                tt_z_0_0[j] = pa_z[j] * tt_0_0_0[j] + 2.0 * fl1_fz * ts_z_0_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForSD(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& pbDistances,
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

        auto pidx_t_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_0_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tt_0_x_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx); 

            auto tt_0_y_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 1); 

            auto tt_0_z_0 = primBuffer.data(pidx_t_0_1_m0 + 3 * idx + 2); 

            auto tt_0_0_0 = primBuffer.data(pidx_t_0_0_m0 + idx); 

            auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx); 

            auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1); 

            auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2); 

            auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3); 

            auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4); 

            auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5); 

            auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx); 

            // set up pointers to integrals

            auto tt_0_xx_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx); 

            auto tt_0_xy_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 1); 

            auto tt_0_xz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 2); 

            auto tt_0_yy_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 3); 

            auto tt_0_yz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 4); 

            auto tt_0_zz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 5); 

            #pragma omp simd aligned(fgb, fx, fz, pb_x, pb_y, pb_z, ts_0_0_0, ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, \
                                     ts_0_yy_0, ts_0_yz_0, ts_0_zz_0, tt_0_0_0, tt_0_x_0, tt_0_xx_0, tt_0_xy_0, \
                                     tt_0_xz_0, tt_0_y_0, tt_0_yy_0, tt_0_yz_0, tt_0_z_0, tt_0_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_0_xx_0[j] = pb_x[j] * tt_0_x_0[j] + 0.5 * fl1_fx * tt_0_0_0[j] + 2.0 * fl1_fz * ts_0_xx_0[j] - fl1_fz * fl1_fgb * ts_0_0_0[j];

                tt_0_xy_0[j] = pb_x[j] * tt_0_y_0[j] + 2.0 * fl1_fz * ts_0_xy_0[j];

                tt_0_xz_0[j] = pb_x[j] * tt_0_z_0[j] + 2.0 * fl1_fz * ts_0_xz_0[j];

                tt_0_yy_0[j] = pb_y[j] * tt_0_y_0[j] + 0.5 * fl1_fx * tt_0_0_0[j] + 2.0 * fl1_fz * ts_0_yy_0[j] - fl1_fz * fl1_fgb * ts_0_0_0[j];

                tt_0_yz_0[j] = pb_y[j] * tt_0_z_0[j] + 2.0 * fl1_fz * ts_0_yz_0[j];

                tt_0_zz_0[j] = pb_z[j] * tt_0_z_0[j] + 0.5 * fl1_fx * tt_0_0_0[j] + 2.0 * fl1_fz * ts_0_zz_0[j] - fl1_fz * fl1_fgb * ts_0_0_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForDS(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_t_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_2_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
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

            auto tt_x_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx); 

            auto tt_y_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx + 1); 

            auto tt_z_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx + 2); 

            auto tt_0_0_0 = primBuffer.data(pidx_t_0_0_m0 + idx); 

            auto ts_xx_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx); 

            auto ts_xy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 1); 

            auto ts_xz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 2); 

            auto ts_yy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 3); 

            auto ts_yz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 4); 

            auto ts_zz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 5); 

            auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx); 

            // set up pointers to integrals

            auto tt_xx_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx); 

            auto tt_xy_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 1); 

            auto tt_xz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 2); 

            auto tt_yy_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 3); 

            auto tt_yz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 4); 

            auto tt_zz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 5); 

            #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, pa_z, ts_0_0_0, ts_xx_0_0, ts_xy_0_0, ts_xz_0_0, \
                                     ts_yy_0_0, ts_yz_0_0, ts_zz_0_0, tt_0_0_0, tt_x_0_0, tt_xx_0_0, tt_xy_0_0, \
                                     tt_xz_0_0, tt_y_0_0, tt_yy_0_0, tt_yz_0_0, tt_z_0_0, tt_zz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fga = fga[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_xx_0_0[j] = pa_x[j] * tt_x_0_0[j] + 0.5 * fl1_fx * tt_0_0_0[j] + 2.0 * fl1_fz * ts_xx_0_0[j] - fl1_fz * fl1_fga * ts_0_0_0[j];

                tt_xy_0_0[j] = pa_x[j] * tt_y_0_0[j] + 2.0 * fl1_fz * ts_xy_0_0[j];

                tt_xz_0_0[j] = pa_x[j] * tt_z_0_0[j] + 2.0 * fl1_fz * ts_xz_0_0[j];

                tt_yy_0_0[j] = pa_y[j] * tt_y_0_0[j] + 0.5 * fl1_fx * tt_0_0_0[j] + 2.0 * fl1_fz * ts_yy_0_0[j] - fl1_fz * fl1_fga * ts_0_0_0[j];

                tt_yz_0_0[j] = pa_y[j] * tt_z_0_0[j] + 2.0 * fl1_fz * ts_yz_0_0[j];

                tt_zz_0_0[j] = pa_z[j] * tt_z_0_0[j] + 0.5 * fl1_fx * tt_0_0_0[j] + 2.0 * fl1_fz * ts_zz_0_0[j] - fl1_fz * fl1_fga * ts_0_0_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForSF(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& pbDistances,
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

        auto pidx_t_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_0_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
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

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx); 

            auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1); 

            auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(fgb, fx, fz, pb_x, pb_y, pb_z, ts_0_x_0, ts_0_xxx_0, ts_0_xxy_0, ts_0_xxz_0, \
                                     ts_0_xyy_0, ts_0_xyz_0, ts_0_xzz_0, ts_0_y_0, ts_0_yyy_0, ts_0_yyz_0, ts_0_yzz_0, \
                                     ts_0_z_0, ts_0_zzz_0, tt_0_x_0, tt_0_xx_0, tt_0_xxx_0, tt_0_xxy_0, tt_0_xxz_0, \
                                     tt_0_xy_0, tt_0_xyy_0, tt_0_xyz_0, tt_0_xz_0, tt_0_xzz_0, tt_0_y_0, tt_0_yy_0, \
                                     tt_0_yyy_0, tt_0_yyz_0, tt_0_yz_0, tt_0_yzz_0, tt_0_z_0, tt_0_zz_0, tt_0_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_0_xxx_0[j] = pb_x[j] * tt_0_xx_0[j] + fl1_fx * tt_0_x_0[j] + 2.0 * fl1_fz * ts_0_xxx_0[j] - 2.0 * fl1_fz * fl1_fgb * ts_0_x_0[j];

                tt_0_xxy_0[j] = pb_x[j] * tt_0_xy_0[j] + 0.5 * fl1_fx * tt_0_y_0[j] + 2.0 * fl1_fz * ts_0_xxy_0[j] - fl1_fz * fl1_fgb * ts_0_y_0[j];

                tt_0_xxz_0[j] = pb_x[j] * tt_0_xz_0[j] + 0.5 * fl1_fx * tt_0_z_0[j] + 2.0 * fl1_fz * ts_0_xxz_0[j] - fl1_fz * fl1_fgb * ts_0_z_0[j];

                tt_0_xyy_0[j] = pb_x[j] * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_0_xyy_0[j];

                tt_0_xyz_0[j] = pb_x[j] * tt_0_yz_0[j] + 2.0 * fl1_fz * ts_0_xyz_0[j];

                tt_0_xzz_0[j] = pb_x[j] * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_0_xzz_0[j];

                tt_0_yyy_0[j] = pb_y[j] * tt_0_yy_0[j] + fl1_fx * tt_0_y_0[j] + 2.0 * fl1_fz * ts_0_yyy_0[j] - 2.0 * fl1_fz * fl1_fgb * ts_0_y_0[j];

                tt_0_yyz_0[j] = pb_y[j] * tt_0_yz_0[j] + 0.5 * fl1_fx * tt_0_z_0[j] + 2.0 * fl1_fz * ts_0_yyz_0[j] - fl1_fz * fl1_fgb * ts_0_z_0[j];

                tt_0_yzz_0[j] = pb_y[j] * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_0_yzz_0[j];

                tt_0_zzz_0[j] = pb_z[j] * tt_0_zz_0[j] + fl1_fx * tt_0_z_0[j] + 2.0 * fl1_fz * ts_0_zzz_0[j] - 2.0 * fl1_fz * fl1_fgb * ts_0_z_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFS(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_t_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_3_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
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

            auto tt_xx_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx); 

            auto tt_xy_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 1); 

            auto tt_xz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 2); 

            auto tt_yy_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 3); 

            auto tt_yz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 4); 

            auto tt_zz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 5); 

            auto tt_x_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx); 

            auto tt_y_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx + 1); 

            auto tt_z_0_0 = primBuffer.data(pidx_t_1_0_m0 + 3 * idx + 2); 

            auto ts_xxx_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx); 

            auto ts_xxy_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 1); 

            auto ts_xxz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 2); 

            auto ts_xyy_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 3); 

            auto ts_xyz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 4); 

            auto ts_xzz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 5); 

            auto ts_yyy_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 6); 

            auto ts_yyz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 7); 

            auto ts_yzz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 8); 

            auto ts_zzz_0_0 = primBuffer.data(pidx_s_3_0_m0 + 10 * idx + 9); 

            auto ts_x_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx); 

            auto ts_y_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 1); 

            auto ts_z_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 2); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, pa_z, ts_x_0_0, ts_xxx_0_0, ts_xxy_0_0, ts_xxz_0_0, \
                                     ts_xyy_0_0, ts_xyz_0_0, ts_xzz_0_0, ts_y_0_0, ts_yyy_0_0, ts_yyz_0_0, ts_yzz_0_0, \
                                     ts_z_0_0, ts_zzz_0_0, tt_x_0_0, tt_xx_0_0, tt_xxx_0_0, tt_xxy_0_0, tt_xxz_0_0, \
                                     tt_xy_0_0, tt_xyy_0_0, tt_xyz_0_0, tt_xz_0_0, tt_xzz_0_0, tt_y_0_0, tt_yy_0_0, \
                                     tt_yyy_0_0, tt_yyz_0_0, tt_yz_0_0, tt_yzz_0_0, tt_z_0_0, tt_zz_0_0, tt_zzz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fga = fga[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_xxx_0_0[j] = pa_x[j] * tt_xx_0_0[j] + fl1_fx * tt_x_0_0[j] + 2.0 * fl1_fz * ts_xxx_0_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_0_0[j];

                tt_xxy_0_0[j] = pa_x[j] * tt_xy_0_0[j] + 0.5 * fl1_fx * tt_y_0_0[j] + 2.0 * fl1_fz * ts_xxy_0_0[j] - fl1_fz * fl1_fga * ts_y_0_0[j];

                tt_xxz_0_0[j] = pa_x[j] * tt_xz_0_0[j] + 0.5 * fl1_fx * tt_z_0_0[j] + 2.0 * fl1_fz * ts_xxz_0_0[j] - fl1_fz * fl1_fga * ts_z_0_0[j];

                tt_xyy_0_0[j] = pa_x[j] * tt_yy_0_0[j] + 2.0 * fl1_fz * ts_xyy_0_0[j];

                tt_xyz_0_0[j] = pa_x[j] * tt_yz_0_0[j] + 2.0 * fl1_fz * ts_xyz_0_0[j];

                tt_xzz_0_0[j] = pa_x[j] * tt_zz_0_0[j] + 2.0 * fl1_fz * ts_xzz_0_0[j];

                tt_yyy_0_0[j] = pa_y[j] * tt_yy_0_0[j] + fl1_fx * tt_y_0_0[j] + 2.0 * fl1_fz * ts_yyy_0_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_0_0[j];

                tt_yyz_0_0[j] = pa_y[j] * tt_yz_0_0[j] + 0.5 * fl1_fx * tt_z_0_0[j] + 2.0 * fl1_fz * ts_yyz_0_0[j] - fl1_fz * fl1_fga * ts_z_0_0[j];

                tt_yzz_0_0[j] = pa_y[j] * tt_zz_0_0[j] + 2.0 * fl1_fz * ts_yzz_0_0[j];

                tt_zzz_0_0[j] = pa_z[j] * tt_zz_0_0[j] + fl1_fx * tt_z_0_0[j] + 2.0 * fl1_fz * ts_zzz_0_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_0_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForSG(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& pbDistances,
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

        auto pidx_t_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_0_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx); 

            auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1); 

            auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2); 

            auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3); 

            auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4); 

            auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(fgb, fx, fz, pb_x, pb_y, pb_z, ts_0_xx_0, ts_0_xxxx_0, ts_0_xxxy_0, \
                                     ts_0_xxxz_0, ts_0_xxyy_0, ts_0_xxyz_0, ts_0_xxzz_0, ts_0_xy_0, ts_0_xyyy_0, \
                                     ts_0_xyyz_0, ts_0_xyzz_0, ts_0_xz_0, ts_0_xzzz_0, ts_0_yy_0, ts_0_yyyy_0, \
                                     ts_0_yyyz_0, ts_0_yyzz_0, ts_0_yz_0, ts_0_yzzz_0, ts_0_zz_0, ts_0_zzzz_0, tt_0_xx_0, \
                                     tt_0_xxx_0, tt_0_xxxx_0, tt_0_xxxy_0, tt_0_xxxz_0, tt_0_xxy_0, tt_0_xxyy_0, \
                                     tt_0_xxyz_0, tt_0_xxz_0, tt_0_xxzz_0, tt_0_xy_0, tt_0_xyy_0, tt_0_xyyy_0, \
                                     tt_0_xyyz_0, tt_0_xyz_0, tt_0_xyzz_0, tt_0_xz_0, tt_0_xzz_0, tt_0_xzzz_0, tt_0_yy_0, \
                                     tt_0_yyy_0, tt_0_yyyy_0, tt_0_yyyz_0, tt_0_yyz_0, tt_0_yyzz_0, tt_0_yz_0, \
                                     tt_0_yzz_0, tt_0_yzzz_0, tt_0_zz_0, tt_0_zzz_0, tt_0_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_0_xxxx_0[j] = pb_x[j] * tt_0_xxx_0[j] + 1.5 * fl1_fx * tt_0_xx_0[j] + 2.0 * fl1_fz * ts_0_xxxx_0[j] - 3.0 * fl1_fz * fl1_fgb * ts_0_xx_0[j];

                tt_0_xxxy_0[j] = pb_x[j] * tt_0_xxy_0[j] + fl1_fx * tt_0_xy_0[j] + 2.0 * fl1_fz * ts_0_xxxy_0[j] - 2.0 * fl1_fz * fl1_fgb * ts_0_xy_0[j];

                tt_0_xxxz_0[j] = pb_x[j] * tt_0_xxz_0[j] + fl1_fx * tt_0_xz_0[j] + 2.0 * fl1_fz * ts_0_xxxz_0[j] - 2.0 * fl1_fz * fl1_fgb * ts_0_xz_0[j];

                tt_0_xxyy_0[j] = pb_x[j] * tt_0_xyy_0[j] + 0.5 * fl1_fx * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_0_xxyy_0[j] - fl1_fz * fl1_fgb * ts_0_yy_0[j];

                tt_0_xxyz_0[j] = pb_x[j] * tt_0_xyz_0[j] + 0.5 * fl1_fx * tt_0_yz_0[j] + 2.0 * fl1_fz * ts_0_xxyz_0[j] - fl1_fz * fl1_fgb * ts_0_yz_0[j];

                tt_0_xxzz_0[j] = pb_x[j] * tt_0_xzz_0[j] + 0.5 * fl1_fx * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_0_xxzz_0[j] - fl1_fz * fl1_fgb * ts_0_zz_0[j];

                tt_0_xyyy_0[j] = pb_x[j] * tt_0_yyy_0[j] + 2.0 * fl1_fz * ts_0_xyyy_0[j];

                tt_0_xyyz_0[j] = pb_x[j] * tt_0_yyz_0[j] + 2.0 * fl1_fz * ts_0_xyyz_0[j];

                tt_0_xyzz_0[j] = pb_x[j] * tt_0_yzz_0[j] + 2.0 * fl1_fz * ts_0_xyzz_0[j];

                tt_0_xzzz_0[j] = pb_x[j] * tt_0_zzz_0[j] + 2.0 * fl1_fz * ts_0_xzzz_0[j];

                tt_0_yyyy_0[j] = pb_y[j] * tt_0_yyy_0[j] + 1.5 * fl1_fx * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_0_yyyy_0[j] - 3.0 * fl1_fz * fl1_fgb * ts_0_yy_0[j];

                tt_0_yyyz_0[j] = pb_y[j] * tt_0_yyz_0[j] + fl1_fx * tt_0_yz_0[j] + 2.0 * fl1_fz * ts_0_yyyz_0[j] - 2.0 * fl1_fz * fl1_fgb * ts_0_yz_0[j];

                tt_0_yyzz_0[j] = pb_y[j] * tt_0_yzz_0[j] + 0.5 * fl1_fx * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_0_yyzz_0[j] - fl1_fz * fl1_fgb * ts_0_zz_0[j];

                tt_0_yzzz_0[j] = pb_y[j] * tt_0_zzz_0[j] + 2.0 * fl1_fz * ts_0_yzzz_0[j];

                tt_0_zzzz_0[j] = pb_z[j] * tt_0_zzz_0[j] + 1.5 * fl1_fx * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_0_zzzz_0[j] - 3.0 * fl1_fz * fl1_fgb * ts_0_zz_0[j];
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGS(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_t_4_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {4, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_t_4_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_t_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_t_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_4_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {4, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
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

            auto tt_xx_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx); 

            auto tt_xy_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 1); 

            auto tt_xz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 2); 

            auto tt_yy_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 3); 

            auto tt_yz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 4); 

            auto tt_zz_0_0 = primBuffer.data(pidx_t_2_0_m0 + 6 * idx + 5); 

            auto ts_xxxx_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx); 

            auto ts_xxxy_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 1); 

            auto ts_xxxz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 2); 

            auto ts_xxyy_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 3); 

            auto ts_xxyz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 4); 

            auto ts_xxzz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 5); 

            auto ts_xyyy_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 6); 

            auto ts_xyyz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 7); 

            auto ts_xyzz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 8); 

            auto ts_xzzz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 9); 

            auto ts_yyyy_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 10); 

            auto ts_yyyz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 11); 

            auto ts_yyzz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 12); 

            auto ts_yzzz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 13); 

            auto ts_zzzz_0_0 = primBuffer.data(pidx_s_4_0_m0 + 15 * idx + 14); 

            auto ts_xx_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx); 

            auto ts_xy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 1); 

            auto ts_xz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 2); 

            auto ts_yy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 3); 

            auto ts_yz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 4); 

            auto ts_zz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 5); 

            // set up pointers to integrals

            auto tt_xxxx_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx); 

            auto tt_xxxy_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 1); 

            auto tt_xxxz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 2); 

            auto tt_xxyy_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 3); 

            auto tt_xxyz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 4); 

            auto tt_xxzz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 5); 

            auto tt_xyyy_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 6); 

            auto tt_xyyz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 7); 

            auto tt_xyzz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 8); 

            auto tt_xzzz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 9); 

            auto tt_yyyy_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 10); 

            auto tt_yyyz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 11); 

            auto tt_yyzz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 12); 

            auto tt_yzzz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 13); 

            auto tt_zzzz_0_0 = primBuffer.data(pidx_t_4_0_m0 + 15 * idx + 14); 

            #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, pa_z, ts_xx_0_0, ts_xxxx_0_0, ts_xxxy_0_0, \
                                     ts_xxxz_0_0, ts_xxyy_0_0, ts_xxyz_0_0, ts_xxzz_0_0, ts_xy_0_0, ts_xyyy_0_0, \
                                     ts_xyyz_0_0, ts_xyzz_0_0, ts_xz_0_0, ts_xzzz_0_0, ts_yy_0_0, ts_yyyy_0_0, \
                                     ts_yyyz_0_0, ts_yyzz_0_0, ts_yz_0_0, ts_yzzz_0_0, ts_zz_0_0, ts_zzzz_0_0, tt_xx_0_0, \
                                     tt_xxx_0_0, tt_xxxx_0_0, tt_xxxy_0_0, tt_xxxz_0_0, tt_xxy_0_0, tt_xxyy_0_0, \
                                     tt_xxyz_0_0, tt_xxz_0_0, tt_xxzz_0_0, tt_xy_0_0, tt_xyy_0_0, tt_xyyy_0_0, \
                                     tt_xyyz_0_0, tt_xyz_0_0, tt_xyzz_0_0, tt_xz_0_0, tt_xzz_0_0, tt_xzzz_0_0, tt_yy_0_0, \
                                     tt_yyy_0_0, tt_yyyy_0_0, tt_yyyz_0_0, tt_yyz_0_0, tt_yyzz_0_0, tt_yz_0_0, \
                                     tt_yzz_0_0, tt_yzzz_0_0, tt_zz_0_0, tt_zzz_0_0, tt_zzzz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fga = fga[j];

                double fl1_fx = fx[j];

                double fl1_fz = fz[j];

                tt_xxxx_0_0[j] = pa_x[j] * tt_xxx_0_0[j] + 1.5 * fl1_fx * tt_xx_0_0[j] + 2.0 * fl1_fz * ts_xxxx_0_0[j] - 3.0 * fl1_fz * fl1_fga * ts_xx_0_0[j];

                tt_xxxy_0_0[j] = pa_x[j] * tt_xxy_0_0[j] + fl1_fx * tt_xy_0_0[j] + 2.0 * fl1_fz * ts_xxxy_0_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xy_0_0[j];

                tt_xxxz_0_0[j] = pa_x[j] * tt_xxz_0_0[j] + fl1_fx * tt_xz_0_0[j] + 2.0 * fl1_fz * ts_xxxz_0_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xz_0_0[j];

                tt_xxyy_0_0[j] = pa_x[j] * tt_xyy_0_0[j] + 0.5 * fl1_fx * tt_yy_0_0[j] + 2.0 * fl1_fz * ts_xxyy_0_0[j] - fl1_fz * fl1_fga * ts_yy_0_0[j];

                tt_xxyz_0_0[j] = pa_x[j] * tt_xyz_0_0[j] + 0.5 * fl1_fx * tt_yz_0_0[j] + 2.0 * fl1_fz * ts_xxyz_0_0[j] - fl1_fz * fl1_fga * ts_yz_0_0[j];

                tt_xxzz_0_0[j] = pa_x[j] * tt_xzz_0_0[j] + 0.5 * fl1_fx * tt_zz_0_0[j] + 2.0 * fl1_fz * ts_xxzz_0_0[j] - fl1_fz * fl1_fga * ts_zz_0_0[j];

                tt_xyyy_0_0[j] = pa_x[j] * tt_yyy_0_0[j] + 2.0 * fl1_fz * ts_xyyy_0_0[j];

                tt_xyyz_0_0[j] = pa_x[j] * tt_yyz_0_0[j] + 2.0 * fl1_fz * ts_xyyz_0_0[j];

                tt_xyzz_0_0[j] = pa_x[j] * tt_yzz_0_0[j] + 2.0 * fl1_fz * ts_xyzz_0_0[j];

                tt_xzzz_0_0[j] = pa_x[j] * tt_zzz_0_0[j] + 2.0 * fl1_fz * ts_xzzz_0_0[j];

                tt_yyyy_0_0[j] = pa_y[j] * tt_yyy_0_0[j] + 1.5 * fl1_fx * tt_yy_0_0[j] + 2.0 * fl1_fz * ts_yyyy_0_0[j] - 3.0 * fl1_fz * fl1_fga * ts_yy_0_0[j];

                tt_yyyz_0_0[j] = pa_y[j] * tt_yyz_0_0[j] + fl1_fx * tt_yz_0_0[j] + 2.0 * fl1_fz * ts_yyyz_0_0[j] - 2.0 * fl1_fz * fl1_fga * ts_yz_0_0[j];

                tt_yyzz_0_0[j] = pa_y[j] * tt_yzz_0_0[j] + 0.5 * fl1_fx * tt_zz_0_0[j] + 2.0 * fl1_fz * ts_yyzz_0_0[j] - fl1_fz * fl1_fga * ts_zz_0_0[j];

                tt_yzzz_0_0[j] = pa_y[j] * tt_zzz_0_0[j] + 2.0 * fl1_fz * ts_yzzz_0_0[j];

                tt_zzzz_0_0[j] = pa_z[j] * tt_zzz_0_0[j] + 1.5 * fl1_fx * tt_zz_0_0[j] + 2.0 * fl1_fz * ts_zzzz_0_0[j] - 3.0 * fl1_fz * fl1_fga * ts_zz_0_0[j];
            }

            idx++;
        }
    }
    
} // kinrecfunc namespace

