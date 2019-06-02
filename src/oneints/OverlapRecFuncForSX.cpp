//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapRecFuncForSX.hpp"

#include "MathConst.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace

    void
    compOverlapForSS(      CMemBlock2D<double>& primBuffer,
                     const CRecursionMap&       recursionMap,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& abDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        // set up index of auxilary integrals
        
        auto pidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));
        
        if (pidx == -1) return;
        
        // set up pointers to primitives data on bra side
        
        auto bnorm = braGtoBlock.getNormFactors();
        
        auto spos = braGtoBlock.getStartPositions();
        
        auto epos = braGtoBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto knorm = ketGtoBlock.getNormFactors();
        
        auto nprim = ketGtoBlock.getNumberOfPrimGtos();
        
        // set up pointers to R(AB) distances
        
        auto abx = abDistances.data(0);
        
        auto aby = abDistances.data(1);
        
        auto abz = abDistances.data(2);
        
        // fetch up pi value
        
        auto fpi = mathconst::getPiValue();
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(2 * idx);
            
            auto fz = osFactors.data(2 * idx + 1);
            
            auto fb = bnorm[i];
            
            // set up primitives buffer data
            
            auto fovl = primBuffer.data(pidx + idx);
            
            #pragma omp simd aligned(fovl, fx, fz, knorm, abx, aby, abz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                fovl[j] = fb * knorm[j] * std::pow(fpi * fx[j], 1.5)
                
                        * std::exp(-fz[j] * (abx[j] * abx[j] + aby[j] * aby[j] +
                                     
                                   abz[j] * abz[j]));
            }
            
            idx++;
        }
    }

    void
    compOverlapForSP(      CMemBlock2D<double>& primBuffer,
                     const CRecursionMap&       recursionMap,
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

        auto pidx_s_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_s_0_1_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx); 

            // set up pointers to integrals

            auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx); 

            auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1); 

            auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2); 

            #pragma omp simd aligned(pb_x, pb_y, pb_z, ts_0_0_0, ts_0_x_0, ts_0_y_0, ts_0_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                ts_0_x_0[j] = pb_x[j] * ts_0_0_0[j];

                ts_0_y_0[j] = pb_y[j] * ts_0_0_0[j];

                ts_0_z_0[j] = pb_z[j] * ts_0_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForPS(      CMemBlock2D<double>& primBuffer,
                     const CRecursionMap&       recursionMap,
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

        auto pidx_s_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_s_1_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx); 

            // set up pointers to integrals

            auto ts_x_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx); 

            auto ts_y_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 1); 

            auto ts_z_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 2); 

            #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_0_0_0, ts_x_0_0, ts_y_0_0, ts_z_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                ts_x_0_0[j] = pa_x[j] * ts_0_0_0[j];

                ts_y_0_0[j] = pa_y[j] * ts_0_0_0[j];

                ts_z_0_0[j] = pa_z[j] * ts_0_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForSD(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_s_0_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_s_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx); 

            auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1); 

            auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2); 

            auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx); 

            // set up pointers to integrals

            auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx); 

            auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1); 

            auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2); 

            auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3); 

            auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4); 

            auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5); 

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, ts_0_0_0, ts_0_x_0, ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, \
                                     ts_0_y_0, ts_0_yy_0, ts_0_yz_0, ts_0_z_0, ts_0_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ts_0_xx_0[j] = pb_x[j] * ts_0_x_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

                ts_0_xy_0[j] = pb_x[j] * ts_0_y_0[j];

                ts_0_xz_0[j] = pb_x[j] * ts_0_z_0[j];

                ts_0_yy_0[j] = pb_y[j] * ts_0_y_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

                ts_0_yz_0[j] = pb_y[j] * ts_0_z_0[j];

                ts_0_zz_0[j] = pb_z[j] * ts_0_z_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForDS(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_s_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_s_2_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_s_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {0, -1, -1, -1}, 
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

            auto ts_x_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx); 

            auto ts_y_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 1); 

            auto ts_z_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 2); 

            auto ts_0_0_0 = primBuffer.data(pidx_s_0_0_m0 + idx); 

            // set up pointers to integrals

            auto ts_xx_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx); 

            auto ts_xy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 1); 

            auto ts_xz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 2); 

            auto ts_yy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 3); 

            auto ts_yz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 4); 

            auto ts_zz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 5); 

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_0_0_0, ts_x_0_0, ts_xx_0_0, ts_xy_0_0, ts_xz_0_0, \
                                     ts_y_0_0, ts_yy_0_0, ts_yz_0_0, ts_z_0_0, ts_zz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ts_xx_0_0[j] = pa_x[j] * ts_x_0_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

                ts_xy_0_0[j] = pa_x[j] * ts_y_0_0[j];

                ts_xz_0_0[j] = pa_x[j] * ts_z_0_0[j];

                ts_yy_0_0[j] = pa_y[j] * ts_y_0_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];

                ts_yz_0_0[j] = pa_y[j] * ts_z_0_0[j];

                ts_zz_0_0[j] = pa_z[j] * ts_z_0_0[j] + 0.5 * fl1_fx * ts_0_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForSF(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_s_0_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx); 

            auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1); 

            auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2); 

            auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3); 

            auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4); 

            auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5); 

            auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx); 

            auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1); 

            auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, ts_0_x_0, ts_0_xx_0, ts_0_xxx_0, ts_0_xxy_0, \
                                     ts_0_xxz_0, ts_0_xy_0, ts_0_xyy_0, ts_0_xyz_0, ts_0_xz_0, ts_0_xzz_0, ts_0_y_0, \
                                     ts_0_yy_0, ts_0_yyy_0, ts_0_yyz_0, ts_0_yz_0, ts_0_yzz_0, ts_0_z_0, ts_0_zz_0, \
                                     ts_0_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ts_0_xxx_0[j] = pb_x[j] * ts_0_xx_0[j] + fl1_fx * ts_0_x_0[j];

                ts_0_xxy_0[j] = pb_x[j] * ts_0_xy_0[j] + 0.5 * fl1_fx * ts_0_y_0[j];

                ts_0_xxz_0[j] = pb_x[j] * ts_0_xz_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

                ts_0_xyy_0[j] = pb_x[j] * ts_0_yy_0[j];

                ts_0_xyz_0[j] = pb_x[j] * ts_0_yz_0[j];

                ts_0_xzz_0[j] = pb_x[j] * ts_0_zz_0[j];

                ts_0_yyy_0[j] = pb_y[j] * ts_0_yy_0[j] + fl1_fx * ts_0_y_0[j];

                ts_0_yyz_0[j] = pb_y[j] * ts_0_yz_0[j] + 0.5 * fl1_fx * ts_0_z_0[j];

                ts_0_yzz_0[j] = pb_y[j] * ts_0_zz_0[j];

                ts_0_zzz_0[j] = pb_z[j] * ts_0_zz_0[j] + fl1_fx * ts_0_z_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForFS(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_s_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_s_3_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_s_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {1, -1, -1, -1}, {0, -1, -1, -1}, 
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

            auto ts_xx_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx); 

            auto ts_xy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 1); 

            auto ts_xz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 2); 

            auto ts_yy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 3); 

            auto ts_yz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 4); 

            auto ts_zz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 5); 

            auto ts_x_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx); 

            auto ts_y_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 1); 

            auto ts_z_0_0 = primBuffer.data(pidx_s_1_0_m0 + 3 * idx + 2); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_x_0_0, ts_xx_0_0, ts_xxx_0_0, ts_xxy_0_0, \
                                     ts_xxz_0_0, ts_xy_0_0, ts_xyy_0_0, ts_xyz_0_0, ts_xz_0_0, ts_xzz_0_0, ts_y_0_0, \
                                     ts_yy_0_0, ts_yyy_0_0, ts_yyz_0_0, ts_yz_0_0, ts_yzz_0_0, ts_z_0_0, ts_zz_0_0, \
                                     ts_zzz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ts_xxx_0_0[j] = pa_x[j] * ts_xx_0_0[j] + fl1_fx * ts_x_0_0[j];

                ts_xxy_0_0[j] = pa_x[j] * ts_xy_0_0[j] + 0.5 * fl1_fx * ts_y_0_0[j];

                ts_xxz_0_0[j] = pa_x[j] * ts_xz_0_0[j] + 0.5 * fl1_fx * ts_z_0_0[j];

                ts_xyy_0_0[j] = pa_x[j] * ts_yy_0_0[j];

                ts_xyz_0_0[j] = pa_x[j] * ts_yz_0_0[j];

                ts_xzz_0_0[j] = pa_x[j] * ts_zz_0_0[j];

                ts_yyy_0_0[j] = pa_y[j] * ts_yy_0_0[j] + fl1_fx * ts_y_0_0[j];

                ts_yyz_0_0[j] = pa_y[j] * ts_yz_0_0[j] + 0.5 * fl1_fx * ts_z_0_0[j];

                ts_yzz_0_0[j] = pa_y[j] * ts_zz_0_0[j];

                ts_zzz_0_0[j] = pa_z[j] * ts_zz_0_0[j] + fl1_fx * ts_z_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForSG(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_s_0_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

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

            auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx); 

            auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1); 

            auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2); 

            auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3); 

            auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4); 

            auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, ts_0_xx_0, ts_0_xxx_0, ts_0_xxxx_0, ts_0_xxxy_0, \
                                     ts_0_xxxz_0, ts_0_xxy_0, ts_0_xxyy_0, ts_0_xxyz_0, ts_0_xxz_0, ts_0_xxzz_0, \
                                     ts_0_xy_0, ts_0_xyy_0, ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyz_0, ts_0_xyzz_0, \
                                     ts_0_xz_0, ts_0_xzz_0, ts_0_xzzz_0, ts_0_yy_0, ts_0_yyy_0, ts_0_yyyy_0, \
                                     ts_0_yyyz_0, ts_0_yyz_0, ts_0_yyzz_0, ts_0_yz_0, ts_0_yzz_0, ts_0_yzzz_0, ts_0_zz_0, \
                                     ts_0_zzz_0, ts_0_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ts_0_xxxx_0[j] = pb_x[j] * ts_0_xxx_0[j] + 1.5 * fl1_fx * ts_0_xx_0[j];

                ts_0_xxxy_0[j] = pb_x[j] * ts_0_xxy_0[j] + fl1_fx * ts_0_xy_0[j];

                ts_0_xxxz_0[j] = pb_x[j] * ts_0_xxz_0[j] + fl1_fx * ts_0_xz_0[j];

                ts_0_xxyy_0[j] = pb_x[j] * ts_0_xyy_0[j] + 0.5 * fl1_fx * ts_0_yy_0[j];

                ts_0_xxyz_0[j] = pb_x[j] * ts_0_xyz_0[j] + 0.5 * fl1_fx * ts_0_yz_0[j];

                ts_0_xxzz_0[j] = pb_x[j] * ts_0_xzz_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];

                ts_0_xyyy_0[j] = pb_x[j] * ts_0_yyy_0[j];

                ts_0_xyyz_0[j] = pb_x[j] * ts_0_yyz_0[j];

                ts_0_xyzz_0[j] = pb_x[j] * ts_0_yzz_0[j];

                ts_0_xzzz_0[j] = pb_x[j] * ts_0_zzz_0[j];

                ts_0_yyyy_0[j] = pb_y[j] * ts_0_yyy_0[j] + 1.5 * fl1_fx * ts_0_yy_0[j];

                ts_0_yyyz_0[j] = pb_y[j] * ts_0_yyz_0[j] + fl1_fx * ts_0_yz_0[j];

                ts_0_yyzz_0[j] = pb_y[j] * ts_0_yzz_0[j] + 0.5 * fl1_fx * ts_0_zz_0[j];

                ts_0_yzzz_0[j] = pb_y[j] * ts_0_zzz_0[j];

                ts_0_zzzz_0[j] = pb_z[j] * ts_0_zzz_0[j] + 1.5 * fl1_fx * ts_0_zz_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForGS(      CMemBlock2D<double>& primBuffer,
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

        auto pidx_s_4_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {4, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_s_4_0_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_s_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {3, -1, -1, -1}, {0, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_s_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                         {2, -1, -1, -1}, {0, -1, -1, -1}, 
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

            auto ts_xx_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx); 

            auto ts_xy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 1); 

            auto ts_xz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 2); 

            auto ts_yy_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 3); 

            auto ts_yz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 4); 

            auto ts_zz_0_0 = primBuffer.data(pidx_s_2_0_m0 + 6 * idx + 5); 

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_xx_0_0, ts_xxx_0_0, ts_xxxx_0_0, ts_xxxy_0_0, \
                                     ts_xxxz_0_0, ts_xxy_0_0, ts_xxyy_0_0, ts_xxyz_0_0, ts_xxz_0_0, ts_xxzz_0_0, \
                                     ts_xy_0_0, ts_xyy_0_0, ts_xyyy_0_0, ts_xyyz_0_0, ts_xyz_0_0, ts_xyzz_0_0, \
                                     ts_xz_0_0, ts_xzz_0_0, ts_xzzz_0_0, ts_yy_0_0, ts_yyy_0_0, ts_yyyy_0_0, \
                                     ts_yyyz_0_0, ts_yyz_0_0, ts_yyzz_0_0, ts_yz_0_0, ts_yzz_0_0, ts_yzzz_0_0, ts_zz_0_0, \
                                     ts_zzz_0_0, ts_zzzz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ts_xxxx_0_0[j] = pa_x[j] * ts_xxx_0_0[j] + 1.5 * fl1_fx * ts_xx_0_0[j];

                ts_xxxy_0_0[j] = pa_x[j] * ts_xxy_0_0[j] + fl1_fx * ts_xy_0_0[j];

                ts_xxxz_0_0[j] = pa_x[j] * ts_xxz_0_0[j] + fl1_fx * ts_xz_0_0[j];

                ts_xxyy_0_0[j] = pa_x[j] * ts_xyy_0_0[j] + 0.5 * fl1_fx * ts_yy_0_0[j];

                ts_xxyz_0_0[j] = pa_x[j] * ts_xyz_0_0[j] + 0.5 * fl1_fx * ts_yz_0_0[j];

                ts_xxzz_0_0[j] = pa_x[j] * ts_xzz_0_0[j] + 0.5 * fl1_fx * ts_zz_0_0[j];

                ts_xyyy_0_0[j] = pa_x[j] * ts_yyy_0_0[j];

                ts_xyyz_0_0[j] = pa_x[j] * ts_yyz_0_0[j];

                ts_xyzz_0_0[j] = pa_x[j] * ts_yzz_0_0[j];

                ts_xzzz_0_0[j] = pa_x[j] * ts_zzz_0_0[j];

                ts_yyyy_0_0[j] = pa_y[j] * ts_yyy_0_0[j] + 1.5 * fl1_fx * ts_yy_0_0[j];

                ts_yyyz_0_0[j] = pa_y[j] * ts_yyz_0_0[j] + fl1_fx * ts_yz_0_0[j];

                ts_yyzz_0_0[j] = pa_y[j] * ts_yzz_0_0[j] + 0.5 * fl1_fx * ts_zz_0_0[j];

                ts_yzzz_0_0[j] = pa_y[j] * ts_zzz_0_0[j];

                ts_zzzz_0_0[j] = pa_z[j] * ts_zzz_0_0[j] + 1.5 * fl1_fx * ts_zz_0_0[j];
            }

            idx++;
        }
    }


} // ovlrecfunc namespace

