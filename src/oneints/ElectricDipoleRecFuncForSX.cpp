//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricDipoleRecFuncForSX.hpp"

namespace ediprecfunc { // ediprecfunc namespace

    void
    compElectricDipoleForSS(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& pcDistances,
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
        
        // set up pointer to overlap integrals
        
        auto sidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                               {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));
        
        // set up pointer to electric dipole integrals
        
        auto didx = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true,
                                                               {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));
        
        if (didx == -1) return; 
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to ditances R(PC)
            
            auto pcx = pcDistances.data(3 * idx);
            
            auto pcy = pcDistances.data(3 * idx + 1);
            
            auto pcz = pcDistances.data(3 * idx + 2);
            
            // set up primitives buffer data
            
            auto fovl = primBuffer.data(sidx + idx);
            
            auto fdipx = primBuffer.data(didx + idx);
            
            auto fdipy = primBuffer.data(didx + bdim + idx);
            
            auto fdipz = primBuffer.data(didx + 2 * bdim + idx);
            
            #pragma omp simd aligned(fovl, fdipx, fdipy, fdipz, pcx,\
                                     pcy, pcz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fx = fovl[j];
                
                fdipx[j] = pcx[j] * fx;
                
                fdipy[j] = pcy[j] * fx;
                
                fdipz[j] = pcz[j] * fx;
            }
            
            idx++;
        }

    }

    void
    compElectricDipoleForSP(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& pbDistances,
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

        auto pidx_d_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_0_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tdx_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + idx); 

            auto tdy_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + bdim + idx); 

            auto tdz_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + 2 * bdim + idx); 

            auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx); 

            // set up pointers to integrals

            auto tdx_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx); 

            auto tdy_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx); 

            auto tdz_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx); 

            auto tdx_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 1); 

            auto tdy_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 2); 

            auto tdy_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 2); 

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, tdx_0_0_0, tdx_0_x_0, tdx_0_y_0, tdx_0_z_0, tdy_0_0_0, \
                                     tdy_0_x_0, tdy_0_y_0, tdy_0_z_0, tdz_0_0_0, tdz_0_x_0, tdz_0_y_0, tdz_0_z_0, \
                                     ts_0_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_0_x_0[j] = pb_x[j] * tdx_0_0_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

                tdy_0_x_0[j] = pb_x[j] * tdy_0_0_0[j];

                tdz_0_x_0[j] = pb_x[j] * tdz_0_0_0[j];

                tdx_0_y_0[j] = pb_y[j] * tdx_0_0_0[j];

                tdy_0_y_0[j] = pb_y[j] * tdy_0_0_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

                tdz_0_y_0[j] = pb_y[j] * tdz_0_0_0[j];

                tdx_0_z_0[j] = pb_z[j] * tdx_0_0_0[j];

                tdy_0_z_0[j] = pb_z[j] * tdy_0_0_0[j];

                tdz_0_z_0[j] = pb_z[j] * tdz_0_0_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForPS(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
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

        auto pidx_d_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_1_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
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

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tdx_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + idx); 

            auto tdy_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + bdim + idx); 

            auto tdz_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + 2 * bdim + idx); 

            auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx); 

            // set up pointers to integrals

            auto tdx_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx); 

            auto tdy_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx); 

            auto tdz_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx); 

            auto tdx_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 1); 

            auto tdy_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 2); 

            auto tdy_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 2); 

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, tdx_0_0_0, tdx_x_0_0, tdx_y_0_0, tdx_z_0_0, tdy_0_0_0, \
                                     tdy_x_0_0, tdy_y_0_0, tdy_z_0_0, tdz_0_0_0, tdz_x_0_0, tdz_y_0_0, tdz_z_0_0, \
                                     ts_0_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_x_0_0[j] = pa_x[j] * tdx_0_0_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

                tdy_x_0_0[j] = pa_x[j] * tdy_0_0_0[j];

                tdz_x_0_0[j] = pa_x[j] * tdz_0_0_0[j];

                tdx_y_0_0[j] = pa_y[j] * tdx_0_0_0[j];

                tdy_y_0_0[j] = pa_y[j] * tdy_0_0_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

                tdz_y_0_0[j] = pa_y[j] * tdz_0_0_0[j];

                tdx_z_0_0[j] = pa_z[j] * tdx_0_0_0[j];

                tdy_z_0_0[j] = pa_z[j] * tdy_0_0_0[j];

                tdz_z_0_0[j] = pa_z[j] * tdz_0_0_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForSD(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& pbDistances,
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

        auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_0_2_m0 == -1) return;

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

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, tdx_0_0_0, tdx_0_x_0, tdx_0_xx_0, tdx_0_xy_0, \
                                     tdx_0_xz_0, tdx_0_y_0, tdx_0_yy_0, tdx_0_yz_0, tdx_0_z_0, tdx_0_zz_0, tdy_0_0_0, \
                                     tdy_0_x_0, tdy_0_xx_0, tdy_0_xy_0, tdy_0_xz_0, tdy_0_y_0, tdy_0_yy_0, tdy_0_yz_0, \
                                     tdy_0_z_0, tdy_0_zz_0, tdz_0_0_0, tdz_0_x_0, tdz_0_xx_0, tdz_0_xy_0, tdz_0_xz_0, \
                                     tdz_0_y_0, tdz_0_yy_0, tdz_0_yz_0, tdz_0_z_0, tdz_0_zz_0, ts_0_x_0, ts_0_y_0, \
                                     ts_0_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_0_xx_0[j] = pb_x[j] * tdx_0_x_0[j] + 0.5 * fl1_fx * tdx_0_0_0[j] + 0.5 * fl1_fx * ts_0_x_0[j];

                tdy_0_xx_0[j] = pb_x[j] * tdy_0_x_0[j] + 0.5 * fl1_fx * tdy_0_0_0[j];

                tdz_0_xx_0[j] = pb_x[j] * tdz_0_x_0[j] + 0.5 * fl1_fx * tdz_0_0_0[j];

                tdx_0_xy_0[j] = pb_x[j] * tdx_0_y_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

                tdy_0_xy_0[j] = pb_x[j] * tdy_0_y_0[j];

                tdz_0_xy_0[j] = pb_x[j] * tdz_0_y_0[j];

                tdx_0_xz_0[j] = pb_x[j] * tdx_0_z_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

                tdy_0_xz_0[j] = pb_x[j] * tdy_0_z_0[j];

                tdz_0_xz_0[j] = pb_x[j] * tdz_0_z_0[j];

                tdx_0_yy_0[j] = pb_y[j] * tdx_0_y_0[j] + 0.5 * fl1_fx * tdx_0_0_0[j];

                tdy_0_yy_0[j] = pb_y[j] * tdy_0_y_0[j] + 0.5 * fl1_fx * tdy_0_0_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

                tdz_0_yy_0[j] = pb_y[j] * tdz_0_y_0[j] + 0.5 * fl1_fx * tdz_0_0_0[j];

                tdx_0_yz_0[j] = pb_y[j] * tdx_0_z_0[j];

                tdy_0_yz_0[j] = pb_y[j] * tdy_0_z_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

                tdz_0_yz_0[j] = pb_y[j] * tdz_0_z_0[j];

                tdx_0_zz_0[j] = pb_z[j] * tdx_0_z_0[j] + 0.5 * fl1_fx * tdx_0_0_0[j];

                tdy_0_zz_0[j] = pb_z[j] * tdy_0_z_0[j] + 0.5 * fl1_fx * tdy_0_0_0[j];

                tdz_0_zz_0[j] = pb_z[j] * tdz_0_z_0[j] + 0.5 * fl1_fx * tdz_0_0_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForDS(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
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

        auto pidx_d_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_2_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
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

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tdx_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx); 

            auto tdy_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx); 

            auto tdz_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx); 

            auto tdx_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 1); 

            auto tdy_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 2); 

            auto tdy_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 2); 

            auto tdx_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + idx); 

            auto tdy_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + bdim + idx); 

            auto tdz_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + 2 * bdim + idx); 

            auto ts_x_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx); 

            auto ts_y_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 1); 

            auto ts_z_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 2); 

            // set up pointers to integrals

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

            auto tdx_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 5); 

            auto tdy_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 5); 

            auto tdz_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 5); 

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, tdx_0_0_0, tdx_x_0_0, tdx_xx_0_0, tdx_xy_0_0, \
                                     tdx_xz_0_0, tdx_y_0_0, tdx_yy_0_0, tdx_yz_0_0, tdx_z_0_0, tdx_zz_0_0, tdy_0_0_0, \
                                     tdy_x_0_0, tdy_xx_0_0, tdy_xy_0_0, tdy_xz_0_0, tdy_y_0_0, tdy_yy_0_0, tdy_yz_0_0, \
                                     tdy_z_0_0, tdy_zz_0_0, tdz_0_0_0, tdz_x_0_0, tdz_xx_0_0, tdz_xy_0_0, tdz_xz_0_0, \
                                     tdz_y_0_0, tdz_yy_0_0, tdz_yz_0_0, tdz_z_0_0, tdz_zz_0_0, ts_x_0_0, ts_y_0_0, \
                                     ts_z_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xx_0_0[j] = pa_x[j] * tdx_x_0_0[j] + 0.5 * fl1_fx * tdx_0_0_0[j] + 0.5 * fl1_fx * ts_x_0_0[j];

                tdy_xx_0_0[j] = pa_x[j] * tdy_x_0_0[j] + 0.5 * fl1_fx * tdy_0_0_0[j];

                tdz_xx_0_0[j] = pa_x[j] * tdz_x_0_0[j] + 0.5 * fl1_fx * tdz_0_0_0[j];

                tdx_xy_0_0[j] = pa_x[j] * tdx_y_0_0[j] + 0.5 * fl1_fx * ts_y_0_0[j];

                tdy_xy_0_0[j] = pa_x[j] * tdy_y_0_0[j];

                tdz_xy_0_0[j] = pa_x[j] * tdz_y_0_0[j];

                tdx_xz_0_0[j] = pa_x[j] * tdx_z_0_0[j] + 0.5 * fl1_fx * ts_z_0_0[j];

                tdy_xz_0_0[j] = pa_x[j] * tdy_z_0_0[j];

                tdz_xz_0_0[j] = pa_x[j] * tdz_z_0_0[j];

                tdx_yy_0_0[j] = pa_y[j] * tdx_y_0_0[j] + 0.5 * fl1_fx * tdx_0_0_0[j];

                tdy_yy_0_0[j] = pa_y[j] * tdy_y_0_0[j] + 0.5 * fl1_fx * tdy_0_0_0[j] + 0.5 * fl1_fx * ts_y_0_0[j];

                tdz_yy_0_0[j] = pa_y[j] * tdz_y_0_0[j] + 0.5 * fl1_fx * tdz_0_0_0[j];

                tdx_yz_0_0[j] = pa_y[j] * tdx_z_0_0[j];

                tdy_yz_0_0[j] = pa_y[j] * tdy_z_0_0[j] + 0.5 * fl1_fx * ts_z_0_0[j];

                tdz_yz_0_0[j] = pa_y[j] * tdz_z_0_0[j];

                tdx_zz_0_0[j] = pa_z[j] * tdx_z_0_0[j] + 0.5 * fl1_fx * tdx_0_0_0[j];

                tdy_zz_0_0[j] = pa_z[j] * tdy_z_0_0[j] + 0.5 * fl1_fx * tdy_0_0_0[j];

                tdz_zz_0_0[j] = pa_z[j] * tdz_z_0_0[j] + 0.5 * fl1_fx * tdz_0_0_0[j] + 0.5 * fl1_fx * ts_z_0_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForSF(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& pbDistances,
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

        auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_0_3_m0 == -1) return;

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

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, tdx_0_x_0, tdx_0_xx_0, tdx_0_xxx_0, tdx_0_xxy_0, \
                                     tdx_0_xxz_0, tdx_0_xy_0, tdx_0_xyy_0, tdx_0_xyz_0, tdx_0_xz_0, tdx_0_xzz_0, \
                                     tdx_0_y_0, tdx_0_yy_0, tdx_0_yyy_0, tdx_0_yyz_0, tdx_0_yz_0, tdx_0_yzz_0, \
                                     tdx_0_z_0, tdx_0_zz_0, tdx_0_zzz_0, tdy_0_x_0, tdy_0_xx_0, tdy_0_xxx_0, \
                                     tdy_0_xxy_0, tdy_0_xxz_0, tdy_0_xy_0, tdy_0_xyy_0, tdy_0_xyz_0, tdy_0_xz_0, \
                                     tdy_0_xzz_0, tdy_0_y_0, tdy_0_yy_0, tdy_0_yyy_0, tdy_0_yyz_0, tdy_0_yz_0, \
                                     tdy_0_yzz_0, tdy_0_z_0, tdy_0_zz_0, tdy_0_zzz_0, tdz_0_x_0, tdz_0_xx_0, tdz_0_xxx_0, \
                                     tdz_0_xxy_0, tdz_0_xxz_0, tdz_0_xy_0, tdz_0_xyy_0, tdz_0_xyz_0, tdz_0_xz_0, \
                                     tdz_0_xzz_0, tdz_0_y_0, tdz_0_yy_0, tdz_0_yyy_0, tdz_0_yyz_0, tdz_0_yz_0, \
                                     tdz_0_yzz_0, tdz_0_z_0, tdz_0_zz_0, tdz_0_zzz_0, ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, \
                                     ts_0_yy_0, ts_0_yz_0, ts_0_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_0_xxx_0[j] = pb_x[j] * tdx_0_xx_0[j] + fl1_fx * tdx_0_x_0[j] + 0.5 * fl1_fx * ts_0_xx_0[j];

                tdy_0_xxx_0[j] = pb_x[j] * tdy_0_xx_0[j] + fl1_fx * tdy_0_x_0[j];

                tdz_0_xxx_0[j] = pb_x[j] * tdz_0_xx_0[j] + fl1_fx * tdz_0_x_0[j];

                tdx_0_xxy_0[j] = pb_x[j] * tdx_0_xy_0[j] + 0.5 * fl1_fx * tdx_0_y_0[j] + 0.5 * fl1_fx * ts_0_xy_0[j];

                tdy_0_xxy_0[j] = pb_x[j] * tdy_0_xy_0[j] + 0.5 * fl1_fx * tdy_0_y_0[j];

                tdz_0_xxy_0[j] = pb_x[j] * tdz_0_xy_0[j] + 0.5 * fl1_fx * tdz_0_y_0[j];

                tdx_0_xxz_0[j] = pb_x[j] * tdx_0_xz_0[j] + 0.5 * fl1_fx * tdx_0_z_0[j] + 0.5 * fl1_fx * ts_0_xz_0[j];

                tdy_0_xxz_0[j] = pb_x[j] * tdy_0_xz_0[j] + 0.5 * fl1_fx * tdy_0_z_0[j];

                tdz_0_xxz_0[j] = pb_x[j] * tdz_0_xz_0[j] + 0.5 * fl1_fx * tdz_0_z_0[j];

                tdx_0_xyy_0[j] = pb_x[j] * tdx_0_yy_0[j] + 0.5 * fl1_fx * ts_0_yy_0[j];

                tdy_0_xyy_0[j] = pb_x[j] * tdy_0_yy_0[j];

                tdz_0_xyy_0[j] = pb_x[j] * tdz_0_yy_0[j];

                tdx_0_xyz_0[j] = pb_x[j] * tdx_0_yz_0[j] + 0.5 * fl1_fx * ts_0_yz_0[j];

                tdy_0_xyz_0[j] = pb_x[j] * tdy_0_yz_0[j];

                tdz_0_xyz_0[j] = pb_x[j] * tdz_0_yz_0[j];

                tdx_0_xzz_0[j] = pb_x[j] * tdx_0_zz_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];

                tdy_0_xzz_0[j] = pb_x[j] * tdy_0_zz_0[j];

                tdz_0_xzz_0[j] = pb_x[j] * tdz_0_zz_0[j];

                tdx_0_yyy_0[j] = pb_y[j] * tdx_0_yy_0[j] + fl1_fx * tdx_0_y_0[j];

                tdy_0_yyy_0[j] = pb_y[j] * tdy_0_yy_0[j] + fl1_fx * tdy_0_y_0[j] + 0.5 * fl1_fx * ts_0_yy_0[j];

                tdz_0_yyy_0[j] = pb_y[j] * tdz_0_yy_0[j] + fl1_fx * tdz_0_y_0[j];

                tdx_0_yyz_0[j] = pb_y[j] * tdx_0_yz_0[j] + 0.5 * fl1_fx * tdx_0_z_0[j];

                tdy_0_yyz_0[j] = pb_y[j] * tdy_0_yz_0[j] + 0.5 * fl1_fx * tdy_0_z_0[j] + 0.5 * fl1_fx * ts_0_yz_0[j];

                tdz_0_yyz_0[j] = pb_y[j] * tdz_0_yz_0[j] + 0.5 * fl1_fx * tdz_0_z_0[j];

                tdx_0_yzz_0[j] = pb_y[j] * tdx_0_zz_0[j];

                tdy_0_yzz_0[j] = pb_y[j] * tdy_0_zz_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];

                tdz_0_yzz_0[j] = pb_y[j] * tdz_0_zz_0[j];

                tdx_0_zzz_0[j] = pb_z[j] * tdx_0_zz_0[j] + fl1_fx * tdx_0_z_0[j];

                tdy_0_zzz_0[j] = pb_z[j] * tdy_0_zz_0[j] + fl1_fx * tdy_0_z_0[j];

                tdz_0_zzz_0[j] = pb_z[j] * tdz_0_zz_0[j] + fl1_fx * tdz_0_z_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForFS(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
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

        auto pidx_d_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_3_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
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

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto tdx_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 5); 

            auto tdy_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 5); 

            auto tdz_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 5); 

            auto tdx_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx); 

            auto tdy_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx); 

            auto tdz_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx); 

            auto tdx_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 1); 

            auto tdy_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 1); 

            auto tdz_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 1); 

            auto tdx_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 2); 

            auto tdy_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 2); 

            auto tdz_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 2); 

            auto ts_xx_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx); 

            auto ts_xy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 1); 

            auto ts_xz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 2); 

            auto ts_yy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 3); 

            auto ts_yz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 4); 

            auto ts_zz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 5); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, tdx_x_0_0, tdx_xx_0_0, tdx_xxx_0_0, tdx_xxy_0_0, \
                                     tdx_xxz_0_0, tdx_xy_0_0, tdx_xyy_0_0, tdx_xyz_0_0, tdx_xz_0_0, tdx_xzz_0_0, \
                                     tdx_y_0_0, tdx_yy_0_0, tdx_yyy_0_0, tdx_yyz_0_0, tdx_yz_0_0, tdx_yzz_0_0, \
                                     tdx_z_0_0, tdx_zz_0_0, tdx_zzz_0_0, tdy_x_0_0, tdy_xx_0_0, tdy_xxx_0_0, \
                                     tdy_xxy_0_0, tdy_xxz_0_0, tdy_xy_0_0, tdy_xyy_0_0, tdy_xyz_0_0, tdy_xz_0_0, \
                                     tdy_xzz_0_0, tdy_y_0_0, tdy_yy_0_0, tdy_yyy_0_0, tdy_yyz_0_0, tdy_yz_0_0, \
                                     tdy_yzz_0_0, tdy_z_0_0, tdy_zz_0_0, tdy_zzz_0_0, tdz_x_0_0, tdz_xx_0_0, tdz_xxx_0_0, \
                                     tdz_xxy_0_0, tdz_xxz_0_0, tdz_xy_0_0, tdz_xyy_0_0, tdz_xyz_0_0, tdz_xz_0_0, \
                                     tdz_xzz_0_0, tdz_y_0_0, tdz_yy_0_0, tdz_yyy_0_0, tdz_yyz_0_0, tdz_yz_0_0, \
                                     tdz_yzz_0_0, tdz_z_0_0, tdz_zz_0_0, tdz_zzz_0_0, ts_xx_0_0, ts_xy_0_0, ts_xz_0_0, \
                                     ts_yy_0_0, ts_yz_0_0, ts_zz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxx_0_0[j] = pa_x[j] * tdx_xx_0_0[j] + fl1_fx * tdx_x_0_0[j] + 0.5 * fl1_fx * ts_xx_0_0[j];

                tdy_xxx_0_0[j] = pa_x[j] * tdy_xx_0_0[j] + fl1_fx * tdy_x_0_0[j];

                tdz_xxx_0_0[j] = pa_x[j] * tdz_xx_0_0[j] + fl1_fx * tdz_x_0_0[j];

                tdx_xxy_0_0[j] = pa_x[j] * tdx_xy_0_0[j] + 0.5 * fl1_fx * tdx_y_0_0[j] + 0.5 * fl1_fx * ts_xy_0_0[j];

                tdy_xxy_0_0[j] = pa_x[j] * tdy_xy_0_0[j] + 0.5 * fl1_fx * tdy_y_0_0[j];

                tdz_xxy_0_0[j] = pa_x[j] * tdz_xy_0_0[j] + 0.5 * fl1_fx * tdz_y_0_0[j];

                tdx_xxz_0_0[j] = pa_x[j] * tdx_xz_0_0[j] + 0.5 * fl1_fx * tdx_z_0_0[j] + 0.5 * fl1_fx * ts_xz_0_0[j];

                tdy_xxz_0_0[j] = pa_x[j] * tdy_xz_0_0[j] + 0.5 * fl1_fx * tdy_z_0_0[j];

                tdz_xxz_0_0[j] = pa_x[j] * tdz_xz_0_0[j] + 0.5 * fl1_fx * tdz_z_0_0[j];

                tdx_xyy_0_0[j] = pa_x[j] * tdx_yy_0_0[j] + 0.5 * fl1_fx * ts_yy_0_0[j];

                tdy_xyy_0_0[j] = pa_x[j] * tdy_yy_0_0[j];

                tdz_xyy_0_0[j] = pa_x[j] * tdz_yy_0_0[j];

                tdx_xyz_0_0[j] = pa_x[j] * tdx_yz_0_0[j] + 0.5 * fl1_fx * ts_yz_0_0[j];

                tdy_xyz_0_0[j] = pa_x[j] * tdy_yz_0_0[j];

                tdz_xyz_0_0[j] = pa_x[j] * tdz_yz_0_0[j];

                tdx_xzz_0_0[j] = pa_x[j] * tdx_zz_0_0[j] + 0.5 * fl1_fx * ts_zz_0_0[j];

                tdy_xzz_0_0[j] = pa_x[j] * tdy_zz_0_0[j];

                tdz_xzz_0_0[j] = pa_x[j] * tdz_zz_0_0[j];

                tdx_yyy_0_0[j] = pa_y[j] * tdx_yy_0_0[j] + fl1_fx * tdx_y_0_0[j];

                tdy_yyy_0_0[j] = pa_y[j] * tdy_yy_0_0[j] + fl1_fx * tdy_y_0_0[j] + 0.5 * fl1_fx * ts_yy_0_0[j];

                tdz_yyy_0_0[j] = pa_y[j] * tdz_yy_0_0[j] + fl1_fx * tdz_y_0_0[j];

                tdx_yyz_0_0[j] = pa_y[j] * tdx_yz_0_0[j] + 0.5 * fl1_fx * tdx_z_0_0[j];

                tdy_yyz_0_0[j] = pa_y[j] * tdy_yz_0_0[j] + 0.5 * fl1_fx * tdy_z_0_0[j] + 0.5 * fl1_fx * ts_yz_0_0[j];

                tdz_yyz_0_0[j] = pa_y[j] * tdz_yz_0_0[j] + 0.5 * fl1_fx * tdz_z_0_0[j];

                tdx_yzz_0_0[j] = pa_y[j] * tdx_zz_0_0[j];

                tdy_yzz_0_0[j] = pa_y[j] * tdy_zz_0_0[j] + 0.5 * fl1_fx * ts_zz_0_0[j];

                tdz_yzz_0_0[j] = pa_y[j] * tdz_zz_0_0[j];

                tdx_zzz_0_0[j] = pa_z[j] * tdx_zz_0_0[j] + fl1_fx * tdx_z_0_0[j];

                tdy_zzz_0_0[j] = pa_z[j] * tdy_zz_0_0[j] + fl1_fx * tdy_z_0_0[j];

                tdz_zzz_0_0[j] = pa_z[j] * tdz_zz_0_0[j] + fl1_fx * tdz_z_0_0[j] + 0.5 * fl1_fx * ts_zz_0_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForSG(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& pbDistances,
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

        auto pidx_d_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_0_4_m0 == -1) return;

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

            auto fx = osFactors.data(nOSFactors * idx);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, tdx_0_xx_0, tdx_0_xxx_0, tdx_0_xxxx_0, tdx_0_xxxy_0, \
                                     tdx_0_xxxz_0, tdx_0_xxy_0, tdx_0_xxyy_0, tdx_0_xxyz_0, tdx_0_xxz_0, tdx_0_xxzz_0, \
                                     tdx_0_xy_0, tdx_0_xyy_0, tdx_0_xyyy_0, tdx_0_xyyz_0, tdx_0_xyz_0, tdx_0_xyzz_0, \
                                     tdx_0_xz_0, tdx_0_xzz_0, tdx_0_xzzz_0, tdx_0_yy_0, tdx_0_yyy_0, tdx_0_yyyy_0, \
                                     tdx_0_yyyz_0, tdx_0_yyz_0, tdx_0_yyzz_0, tdx_0_yz_0, tdx_0_yzz_0, tdx_0_yzzz_0, \
                                     tdx_0_zz_0, tdx_0_zzz_0, tdx_0_zzzz_0, tdy_0_xx_0, tdy_0_xxx_0, tdy_0_xxxx_0, \
                                     tdy_0_xxxy_0, tdy_0_xxxz_0, tdy_0_xxy_0, tdy_0_xxyy_0, tdy_0_xxyz_0, tdy_0_xxz_0, \
                                     tdy_0_xxzz_0, tdy_0_xy_0, tdy_0_xyy_0, tdy_0_xyyy_0, tdy_0_xyyz_0, tdy_0_xyz_0, \
                                     tdy_0_xyzz_0, tdy_0_xz_0, tdy_0_xzz_0, tdy_0_xzzz_0, tdy_0_yy_0, tdy_0_yyy_0, \
                                     tdy_0_yyyy_0, tdy_0_yyyz_0, tdy_0_yyz_0, tdy_0_yyzz_0, tdy_0_yz_0, tdy_0_yzz_0, \
                                     tdy_0_yzzz_0, tdy_0_zz_0, tdy_0_zzz_0, tdy_0_zzzz_0, tdz_0_xx_0, tdz_0_xxx_0, \
                                     tdz_0_xxxx_0, tdz_0_xxxy_0, tdz_0_xxxz_0, tdz_0_xxy_0, tdz_0_xxyy_0, tdz_0_xxyz_0, \
                                     tdz_0_xxz_0, tdz_0_xxzz_0, tdz_0_xy_0, tdz_0_xyy_0, tdz_0_xyyy_0, tdz_0_xyyz_0, \
                                     tdz_0_xyz_0, tdz_0_xyzz_0, tdz_0_xz_0, tdz_0_xzz_0, tdz_0_xzzz_0, tdz_0_yy_0, \
                                     tdz_0_yyy_0, tdz_0_yyyy_0, tdz_0_yyyz_0, tdz_0_yyz_0, tdz_0_yyzz_0, tdz_0_yz_0, \
                                     tdz_0_yzz_0, tdz_0_yzzz_0, tdz_0_zz_0, tdz_0_zzz_0, tdz_0_zzzz_0, ts_0_xxx_0, \
                                     ts_0_xxy_0, ts_0_xxz_0, ts_0_xyy_0, ts_0_xyz_0, ts_0_xzz_0, ts_0_yyy_0, ts_0_yyz_0, \
                                     ts_0_yzz_0, ts_0_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_0_xxxx_0[j] = pb_x[j] * tdx_0_xxx_0[j] + 1.5 * fl1_fx * tdx_0_xx_0[j] + 0.5 * fl1_fx * ts_0_xxx_0[j];

                tdy_0_xxxx_0[j] = pb_x[j] * tdy_0_xxx_0[j] + 1.5 * fl1_fx * tdy_0_xx_0[j];

                tdz_0_xxxx_0[j] = pb_x[j] * tdz_0_xxx_0[j] + 1.5 * fl1_fx * tdz_0_xx_0[j];

                tdx_0_xxxy_0[j] = pb_x[j] * tdx_0_xxy_0[j] + fl1_fx * tdx_0_xy_0[j] + 0.5 * fl1_fx * ts_0_xxy_0[j];

                tdy_0_xxxy_0[j] = pb_x[j] * tdy_0_xxy_0[j] + fl1_fx * tdy_0_xy_0[j];

                tdz_0_xxxy_0[j] = pb_x[j] * tdz_0_xxy_0[j] + fl1_fx * tdz_0_xy_0[j];

                tdx_0_xxxz_0[j] = pb_x[j] * tdx_0_xxz_0[j] + fl1_fx * tdx_0_xz_0[j] + 0.5 * fl1_fx * ts_0_xxz_0[j];

                tdy_0_xxxz_0[j] = pb_x[j] * tdy_0_xxz_0[j] + fl1_fx * tdy_0_xz_0[j];

                tdz_0_xxxz_0[j] = pb_x[j] * tdz_0_xxz_0[j] + fl1_fx * tdz_0_xz_0[j];

                tdx_0_xxyy_0[j] = pb_x[j] * tdx_0_xyy_0[j] + 0.5 * fl1_fx * tdx_0_yy_0[j] + 0.5 * fl1_fx * ts_0_xyy_0[j];

                tdy_0_xxyy_0[j] = pb_x[j] * tdy_0_xyy_0[j] + 0.5 * fl1_fx * tdy_0_yy_0[j];

                tdz_0_xxyy_0[j] = pb_x[j] * tdz_0_xyy_0[j] + 0.5 * fl1_fx * tdz_0_yy_0[j];

                tdx_0_xxyz_0[j] = pb_x[j] * tdx_0_xyz_0[j] + 0.5 * fl1_fx * tdx_0_yz_0[j] + 0.5 * fl1_fx * ts_0_xyz_0[j];

                tdy_0_xxyz_0[j] = pb_x[j] * tdy_0_xyz_0[j] + 0.5 * fl1_fx * tdy_0_yz_0[j];

                tdz_0_xxyz_0[j] = pb_x[j] * tdz_0_xyz_0[j] + 0.5 * fl1_fx * tdz_0_yz_0[j];

                tdx_0_xxzz_0[j] = pb_x[j] * tdx_0_xzz_0[j] + 0.5 * fl1_fx * tdx_0_zz_0[j] + 0.5 * fl1_fx * ts_0_xzz_0[j];

                tdy_0_xxzz_0[j] = pb_x[j] * tdy_0_xzz_0[j] + 0.5 * fl1_fx * tdy_0_zz_0[j];

                tdz_0_xxzz_0[j] = pb_x[j] * tdz_0_xzz_0[j] + 0.5 * fl1_fx * tdz_0_zz_0[j];

                tdx_0_xyyy_0[j] = pb_x[j] * tdx_0_yyy_0[j] + 0.5 * fl1_fx * ts_0_yyy_0[j];

                tdy_0_xyyy_0[j] = pb_x[j] * tdy_0_yyy_0[j];

                tdz_0_xyyy_0[j] = pb_x[j] * tdz_0_yyy_0[j];

                tdx_0_xyyz_0[j] = pb_x[j] * tdx_0_yyz_0[j] + 0.5 * fl1_fx * ts_0_yyz_0[j];

                tdy_0_xyyz_0[j] = pb_x[j] * tdy_0_yyz_0[j];

                tdz_0_xyyz_0[j] = pb_x[j] * tdz_0_yyz_0[j];

                tdx_0_xyzz_0[j] = pb_x[j] * tdx_0_yzz_0[j] + 0.5 * fl1_fx * ts_0_yzz_0[j];

                tdy_0_xyzz_0[j] = pb_x[j] * tdy_0_yzz_0[j];

                tdz_0_xyzz_0[j] = pb_x[j] * tdz_0_yzz_0[j];

                tdx_0_xzzz_0[j] = pb_x[j] * tdx_0_zzz_0[j] + 0.5 * fl1_fx * ts_0_zzz_0[j];

                tdy_0_xzzz_0[j] = pb_x[j] * tdy_0_zzz_0[j];

                tdz_0_xzzz_0[j] = pb_x[j] * tdz_0_zzz_0[j];

                tdx_0_yyyy_0[j] = pb_y[j] * tdx_0_yyy_0[j] + 1.5 * fl1_fx * tdx_0_yy_0[j];

                tdy_0_yyyy_0[j] = pb_y[j] * tdy_0_yyy_0[j] + 1.5 * fl1_fx * tdy_0_yy_0[j] + 0.5 * fl1_fx * ts_0_yyy_0[j];

                tdz_0_yyyy_0[j] = pb_y[j] * tdz_0_yyy_0[j] + 1.5 * fl1_fx * tdz_0_yy_0[j];

                tdx_0_yyyz_0[j] = pb_y[j] * tdx_0_yyz_0[j] + fl1_fx * tdx_0_yz_0[j];

                tdy_0_yyyz_0[j] = pb_y[j] * tdy_0_yyz_0[j] + fl1_fx * tdy_0_yz_0[j] + 0.5 * fl1_fx * ts_0_yyz_0[j];

                tdz_0_yyyz_0[j] = pb_y[j] * tdz_0_yyz_0[j] + fl1_fx * tdz_0_yz_0[j];

                tdx_0_yyzz_0[j] = pb_y[j] * tdx_0_yzz_0[j] + 0.5 * fl1_fx * tdx_0_zz_0[j];

                tdy_0_yyzz_0[j] = pb_y[j] * tdy_0_yzz_0[j] + 0.5 * fl1_fx * tdy_0_zz_0[j] + 0.5 * fl1_fx * ts_0_yzz_0[j];

                tdz_0_yyzz_0[j] = pb_y[j] * tdz_0_yzz_0[j] + 0.5 * fl1_fx * tdz_0_zz_0[j];

                tdx_0_yzzz_0[j] = pb_y[j] * tdx_0_zzz_0[j];

                tdy_0_yzzz_0[j] = pb_y[j] * tdy_0_zzz_0[j] + 0.5 * fl1_fx * ts_0_zzz_0[j];

                tdz_0_yzzz_0[j] = pb_y[j] * tdz_0_zzz_0[j];

                tdx_0_zzzz_0[j] = pb_z[j] * tdx_0_zzz_0[j] + 1.5 * fl1_fx * tdx_0_zz_0[j];

                tdy_0_zzzz_0[j] = pb_z[j] * tdy_0_zzz_0[j] + 1.5 * fl1_fx * tdy_0_zz_0[j];

                tdz_0_zzzz_0[j] = pb_z[j] * tdz_0_zzz_0[j] + 1.5 * fl1_fx * tdz_0_zz_0[j] + 0.5 * fl1_fx * ts_0_zzz_0[j];
            }

            idx++;
        }
    }

    void
    compElectricDipoleForGS(      CMemBlock2D<double>& primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
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

        auto pidx_d_4_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {4, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_d_4_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_d_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, 
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
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

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto tdx_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 5); 

            auto tdy_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 5); 

            auto tdz_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 5); 

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

            // set up pointers to integrals

            auto tdx_xxxx_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx); 

            auto tdy_xxxx_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx); 

            auto tdz_xxxx_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx); 

            auto tdx_xxxy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 1); 

            auto tdy_xxxy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 1); 

            auto tdz_xxxy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 1); 

            auto tdx_xxxz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 2); 

            auto tdy_xxxz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 2); 

            auto tdz_xxxz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 2); 

            auto tdx_xxyy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 3); 

            auto tdy_xxyy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 3); 

            auto tdz_xxyy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 3); 

            auto tdx_xxyz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 4); 

            auto tdy_xxyz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 4); 

            auto tdz_xxyz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 4); 

            auto tdx_xxzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 5); 

            auto tdy_xxzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 5); 

            auto tdz_xxzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 5); 

            auto tdx_xyyy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 6); 

            auto tdy_xyyy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 6); 

            auto tdz_xyyy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 6); 

            auto tdx_xyyz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 7); 

            auto tdy_xyyz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 7); 

            auto tdz_xyyz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 7); 

            auto tdx_xyzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 8); 

            auto tdy_xyzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 8); 

            auto tdz_xyzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 8); 

            auto tdx_xzzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 9); 

            auto tdy_xzzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 9); 

            auto tdz_xzzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 9); 

            auto tdx_yyyy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 10); 

            auto tdy_yyyy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 10); 

            auto tdz_yyyy_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 10); 

            auto tdx_yyyz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 11); 

            auto tdy_yyyz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 11); 

            auto tdz_yyyz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 11); 

            auto tdx_yyzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 12); 

            auto tdy_yyzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 12); 

            auto tdz_yyzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 12); 

            auto tdx_yzzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 13); 

            auto tdy_yzzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 13); 

            auto tdz_yzzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 13); 

            auto tdx_zzzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * idx + 14); 

            auto tdy_zzzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 15 * bdim + 15 * idx + 14); 

            auto tdz_zzzz_0_0 = primBuffer.data(pidx_d_4_0_m0 + 30 * bdim + 15 * idx + 14); 

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, tdx_xx_0_0, tdx_xxx_0_0, tdx_xxxx_0_0, tdx_xxxy_0_0, \
                                     tdx_xxxz_0_0, tdx_xxy_0_0, tdx_xxyy_0_0, tdx_xxyz_0_0, tdx_xxz_0_0, tdx_xxzz_0_0, \
                                     tdx_xy_0_0, tdx_xyy_0_0, tdx_xyyy_0_0, tdx_xyyz_0_0, tdx_xyz_0_0, tdx_xyzz_0_0, \
                                     tdx_xz_0_0, tdx_xzz_0_0, tdx_xzzz_0_0, tdx_yy_0_0, tdx_yyy_0_0, tdx_yyyy_0_0, \
                                     tdx_yyyz_0_0, tdx_yyz_0_0, tdx_yyzz_0_0, tdx_yz_0_0, tdx_yzz_0_0, tdx_yzzz_0_0, \
                                     tdx_zz_0_0, tdx_zzz_0_0, tdx_zzzz_0_0, tdy_xx_0_0, tdy_xxx_0_0, tdy_xxxx_0_0, \
                                     tdy_xxxy_0_0, tdy_xxxz_0_0, tdy_xxy_0_0, tdy_xxyy_0_0, tdy_xxyz_0_0, tdy_xxz_0_0, \
                                     tdy_xxzz_0_0, tdy_xy_0_0, tdy_xyy_0_0, tdy_xyyy_0_0, tdy_xyyz_0_0, tdy_xyz_0_0, \
                                     tdy_xyzz_0_0, tdy_xz_0_0, tdy_xzz_0_0, tdy_xzzz_0_0, tdy_yy_0_0, tdy_yyy_0_0, \
                                     tdy_yyyy_0_0, tdy_yyyz_0_0, tdy_yyz_0_0, tdy_yyzz_0_0, tdy_yz_0_0, tdy_yzz_0_0, \
                                     tdy_yzzz_0_0, tdy_zz_0_0, tdy_zzz_0_0, tdy_zzzz_0_0, tdz_xx_0_0, tdz_xxx_0_0, \
                                     tdz_xxxx_0_0, tdz_xxxy_0_0, tdz_xxxz_0_0, tdz_xxy_0_0, tdz_xxyy_0_0, tdz_xxyz_0_0, \
                                     tdz_xxz_0_0, tdz_xxzz_0_0, tdz_xy_0_0, tdz_xyy_0_0, tdz_xyyy_0_0, tdz_xyyz_0_0, \
                                     tdz_xyz_0_0, tdz_xyzz_0_0, tdz_xz_0_0, tdz_xzz_0_0, tdz_xzzz_0_0, tdz_yy_0_0, \
                                     tdz_yyy_0_0, tdz_yyyy_0_0, tdz_yyyz_0_0, tdz_yyz_0_0, tdz_yyzz_0_0, tdz_yz_0_0, \
                                     tdz_yzz_0_0, tdz_yzzz_0_0, tdz_zz_0_0, tdz_zzz_0_0, tdz_zzzz_0_0, ts_xxx_0_0, \
                                     ts_xxy_0_0, ts_xxz_0_0, ts_xyy_0_0, ts_xyz_0_0, ts_xzz_0_0, ts_yyy_0_0, ts_yyz_0_0, \
                                     ts_yzz_0_0, ts_zzz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tdx_xxxx_0_0[j] = pa_x[j] * tdx_xxx_0_0[j] + 1.5 * fl1_fx * tdx_xx_0_0[j] + 0.5 * fl1_fx * ts_xxx_0_0[j];

                tdy_xxxx_0_0[j] = pa_x[j] * tdy_xxx_0_0[j] + 1.5 * fl1_fx * tdy_xx_0_0[j];

                tdz_xxxx_0_0[j] = pa_x[j] * tdz_xxx_0_0[j] + 1.5 * fl1_fx * tdz_xx_0_0[j];

                tdx_xxxy_0_0[j] = pa_x[j] * tdx_xxy_0_0[j] + fl1_fx * tdx_xy_0_0[j] + 0.5 * fl1_fx * ts_xxy_0_0[j];

                tdy_xxxy_0_0[j] = pa_x[j] * tdy_xxy_0_0[j] + fl1_fx * tdy_xy_0_0[j];

                tdz_xxxy_0_0[j] = pa_x[j] * tdz_xxy_0_0[j] + fl1_fx * tdz_xy_0_0[j];

                tdx_xxxz_0_0[j] = pa_x[j] * tdx_xxz_0_0[j] + fl1_fx * tdx_xz_0_0[j] + 0.5 * fl1_fx * ts_xxz_0_0[j];

                tdy_xxxz_0_0[j] = pa_x[j] * tdy_xxz_0_0[j] + fl1_fx * tdy_xz_0_0[j];

                tdz_xxxz_0_0[j] = pa_x[j] * tdz_xxz_0_0[j] + fl1_fx * tdz_xz_0_0[j];

                tdx_xxyy_0_0[j] = pa_x[j] * tdx_xyy_0_0[j] + 0.5 * fl1_fx * tdx_yy_0_0[j] + 0.5 * fl1_fx * ts_xyy_0_0[j];

                tdy_xxyy_0_0[j] = pa_x[j] * tdy_xyy_0_0[j] + 0.5 * fl1_fx * tdy_yy_0_0[j];

                tdz_xxyy_0_0[j] = pa_x[j] * tdz_xyy_0_0[j] + 0.5 * fl1_fx * tdz_yy_0_0[j];

                tdx_xxyz_0_0[j] = pa_x[j] * tdx_xyz_0_0[j] + 0.5 * fl1_fx * tdx_yz_0_0[j] + 0.5 * fl1_fx * ts_xyz_0_0[j];

                tdy_xxyz_0_0[j] = pa_x[j] * tdy_xyz_0_0[j] + 0.5 * fl1_fx * tdy_yz_0_0[j];

                tdz_xxyz_0_0[j] = pa_x[j] * tdz_xyz_0_0[j] + 0.5 * fl1_fx * tdz_yz_0_0[j];

                tdx_xxzz_0_0[j] = pa_x[j] * tdx_xzz_0_0[j] + 0.5 * fl1_fx * tdx_zz_0_0[j] + 0.5 * fl1_fx * ts_xzz_0_0[j];

                tdy_xxzz_0_0[j] = pa_x[j] * tdy_xzz_0_0[j] + 0.5 * fl1_fx * tdy_zz_0_0[j];

                tdz_xxzz_0_0[j] = pa_x[j] * tdz_xzz_0_0[j] + 0.5 * fl1_fx * tdz_zz_0_0[j];

                tdx_xyyy_0_0[j] = pa_x[j] * tdx_yyy_0_0[j] + 0.5 * fl1_fx * ts_yyy_0_0[j];

                tdy_xyyy_0_0[j] = pa_x[j] * tdy_yyy_0_0[j];

                tdz_xyyy_0_0[j] = pa_x[j] * tdz_yyy_0_0[j];

                tdx_xyyz_0_0[j] = pa_x[j] * tdx_yyz_0_0[j] + 0.5 * fl1_fx * ts_yyz_0_0[j];

                tdy_xyyz_0_0[j] = pa_x[j] * tdy_yyz_0_0[j];

                tdz_xyyz_0_0[j] = pa_x[j] * tdz_yyz_0_0[j];

                tdx_xyzz_0_0[j] = pa_x[j] * tdx_yzz_0_0[j] + 0.5 * fl1_fx * ts_yzz_0_0[j];

                tdy_xyzz_0_0[j] = pa_x[j] * tdy_yzz_0_0[j];

                tdz_xyzz_0_0[j] = pa_x[j] * tdz_yzz_0_0[j];

                tdx_xzzz_0_0[j] = pa_x[j] * tdx_zzz_0_0[j] + 0.5 * fl1_fx * ts_zzz_0_0[j];

                tdy_xzzz_0_0[j] = pa_x[j] * tdy_zzz_0_0[j];

                tdz_xzzz_0_0[j] = pa_x[j] * tdz_zzz_0_0[j];

                tdx_yyyy_0_0[j] = pa_y[j] * tdx_yyy_0_0[j] + 1.5 * fl1_fx * tdx_yy_0_0[j];

                tdy_yyyy_0_0[j] = pa_y[j] * tdy_yyy_0_0[j] + 1.5 * fl1_fx * tdy_yy_0_0[j] + 0.5 * fl1_fx * ts_yyy_0_0[j];

                tdz_yyyy_0_0[j] = pa_y[j] * tdz_yyy_0_0[j] + 1.5 * fl1_fx * tdz_yy_0_0[j];

                tdx_yyyz_0_0[j] = pa_y[j] * tdx_yyz_0_0[j] + fl1_fx * tdx_yz_0_0[j];

                tdy_yyyz_0_0[j] = pa_y[j] * tdy_yyz_0_0[j] + fl1_fx * tdy_yz_0_0[j] + 0.5 * fl1_fx * ts_yyz_0_0[j];

                tdz_yyyz_0_0[j] = pa_y[j] * tdz_yyz_0_0[j] + fl1_fx * tdz_yz_0_0[j];

                tdx_yyzz_0_0[j] = pa_y[j] * tdx_yzz_0_0[j] + 0.5 * fl1_fx * tdx_zz_0_0[j];

                tdy_yyzz_0_0[j] = pa_y[j] * tdy_yzz_0_0[j] + 0.5 * fl1_fx * tdy_zz_0_0[j] + 0.5 * fl1_fx * ts_yzz_0_0[j];

                tdz_yyzz_0_0[j] = pa_y[j] * tdz_yzz_0_0[j] + 0.5 * fl1_fx * tdz_zz_0_0[j];

                tdx_yzzz_0_0[j] = pa_y[j] * tdx_zzz_0_0[j];

                tdy_yzzz_0_0[j] = pa_y[j] * tdy_zzz_0_0[j] + 0.5 * fl1_fx * ts_zzz_0_0[j];

                tdz_yzzz_0_0[j] = pa_y[j] * tdz_zzz_0_0[j];

                tdx_zzzz_0_0[j] = pa_z[j] * tdx_zzz_0_0[j] + 1.5 * fl1_fx * tdx_zz_0_0[j];

                tdy_zzzz_0_0[j] = pa_z[j] * tdy_zzz_0_0[j] + 1.5 * fl1_fx * tdy_zz_0_0[j];

                tdz_zzzz_0_0[j] = pa_z[j] * tdz_zzz_0_0[j] + 1.5 * fl1_fx * tdz_zz_0_0[j] + 0.5 * fl1_fx * ts_zzz_0_0[j];
            }

            idx++;
        }
    }


} // ediprecfunc namespace

