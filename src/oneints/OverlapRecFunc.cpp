//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapRecFunc.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace

    void
    compOverlapForSS(      CMemBlock2D<double>& primBuffer,
                           CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& abDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
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
        
        // flag for (s||s) integrals generation
        
        bool doints = ((braGtoBlock.getAngularMomentum() == 0) &&
                       (ketGtoBlock.getAngularMomentum() == 0));
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors
            
            auto fx = osFactors.data(2 * idx);
            
            auto fz = osFactors.data(2 * idx + 1);
            
            auto fb = bnorm[i];
            
            // set up primitives buffer data
            
            auto fovl = (doints) ? primBuffer.data(idx) : auxBuffer.data(idx);
            
            #pragma omp simd aligned(fovl, fx, fz, knorm, abx, aby,\
                                     abz: VLX_ALIGN)
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
                     const CMemBlock2D<double>& auxBuffer,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_0_x = primBuffer.data(3 * idx);

            auto t_0_y = primBuffer.data(3 * idx + 1);

            auto t_0_z = primBuffer.data(3 * idx + 2);

            // Batch of Integrals (0) = (0,3)

            #pragma omp simd aligned(pb_x, pb_y, pb_z, s_0_0, t_0_x, t_0_y, t_0_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_0_x[j] = pb_x[j] * s_0_0[j];

                t_0_y[j] = pb_y[j] * s_0_0[j];

                t_0_z[j] = pb_z[j] * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForSD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_0_xx = primBuffer.data(6 * idx);

            auto t_0_xy = primBuffer.data(6 * idx + 1);

            auto t_0_xz = primBuffer.data(6 * idx + 2);

            auto t_0_yy = primBuffer.data(6 * idx + 3);

            auto t_0_yz = primBuffer.data(6 * idx + 4);

            auto t_0_zz = primBuffer.data(6 * idx + 5);

            // Batch of Integrals (0) = (0,6)

            #pragma omp simd aligned(fx, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, t_0_xx, t_0_xy, t_0_xz, \
                                     t_0_yy, t_0_yz, t_0_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_0_xx[j] = (0.5 * fx[j] + pb_xx[j]) * s_0_0[j];

                t_0_xy[j] = pb_xy[j] * s_0_0[j];

                t_0_xz[j] = pb_xz[j] * s_0_0[j];

                t_0_yy[j] = (0.5 * fx[j] + pb_yy[j]) * s_0_0[j];

                t_0_yz[j] = pb_yz[j] * s_0_0[j];

                t_0_zz[j] = (0.5 * fx[j] + pb_zz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForSF(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_0_xxx = primBuffer.data(10 * idx);

            auto t_0_xxy = primBuffer.data(10 * idx + 1);

            auto t_0_xxz = primBuffer.data(10 * idx + 2);

            auto t_0_xyy = primBuffer.data(10 * idx + 3);

            auto t_0_xyz = primBuffer.data(10 * idx + 4);

            auto t_0_xzz = primBuffer.data(10 * idx + 5);

            auto t_0_yyy = primBuffer.data(10 * idx + 6);

            auto t_0_yyz = primBuffer.data(10 * idx + 7);

            auto t_0_yzz = primBuffer.data(10 * idx + 8);

            auto t_0_zzz = primBuffer.data(10 * idx + 9);

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fx, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_xzz, pb_y, pb_yyy, pb_yyz, \
                                     pb_yzz, pb_z, pb_zzz, s_0_0, t_0_xxx, t_0_xxy, t_0_xxz, t_0_xyy, t_0_xyz, t_0_xzz, \
                                     t_0_yyy, t_0_yyz, t_0_yzz, t_0_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_0_xxx[j] = (1.5 * pb_x[j] * fx[j] + pb_xxx[j]) * s_0_0[j];

                t_0_xxy[j] = (0.5 * fx[j] * pb_y[j] + pb_xxy[j]) * s_0_0[j];

                t_0_xxz[j] = (0.5 * fx[j] * pb_z[j] + pb_xxz[j]) * s_0_0[j];

                t_0_xyy[j] = (0.5 * pb_x[j] * fx[j] + pb_xyy[j]) * s_0_0[j];

                t_0_xyz[j] = pb_xyz[j] * s_0_0[j];

                t_0_xzz[j] = (0.5 * pb_x[j] * fx[j] + pb_xzz[j]) * s_0_0[j];

                t_0_yyy[j] = (1.5 * pb_y[j] * fx[j] + pb_yyy[j]) * s_0_0[j];

                t_0_yyz[j] = (0.5 * fx[j] * pb_z[j] + pb_yyz[j]) * s_0_0[j];

                t_0_yzz[j] = (0.5 * pb_y[j] * fx[j] + pb_yzz[j]) * s_0_0[j];

                t_0_zzz[j] = (1.5 * pb_z[j] * fx[j] + pb_zzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForSG(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_0_xxxx = primBuffer.data(15 * idx);

            auto t_0_xxxy = primBuffer.data(15 * idx + 1);

            auto t_0_xxxz = primBuffer.data(15 * idx + 2);

            auto t_0_xxyy = primBuffer.data(15 * idx + 3);

            auto t_0_xxyz = primBuffer.data(15 * idx + 4);

            auto t_0_xxzz = primBuffer.data(15 * idx + 5);

            auto t_0_xyyy = primBuffer.data(15 * idx + 6);

            auto t_0_xyyz = primBuffer.data(15 * idx + 7);

            auto t_0_xyzz = primBuffer.data(15 * idx + 8);

            auto t_0_xzzz = primBuffer.data(15 * idx + 9);

            auto t_0_yyyy = primBuffer.data(15 * idx + 10);

            auto t_0_yyyz = primBuffer.data(15 * idx + 11);

            auto t_0_yyzz = primBuffer.data(15 * idx + 12);

            auto t_0_yzzz = primBuffer.data(15 * idx + 13);

            auto t_0_zzzz = primBuffer.data(15 * idx + 14);

            // Batch of Integrals (0) = (0,8)

            #pragma omp simd aligned(fx, pb_xx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxyy, pb_xxyz, pb_xxzz, pb_xy, \
                                     pb_xyyy, pb_xyyz, pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, t_0_xxxx, t_0_xxxy, t_0_xxxz, \
                                     t_0_xxyy, t_0_xxyz, t_0_xxzz, t_0_xyyy, t_0_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_0_xxxx[j] = (0.75 * fx[j] * fx[j] + 3.0 * pb_xx[j] * fx[j] + pb_xxxx[j]) * s_0_0[j];

                t_0_xxxy[j] = (1.5 * pb_xy[j] * fx[j] + pb_xxxy[j]) * s_0_0[j];

                t_0_xxxz[j] = (1.5 * pb_xz[j] * fx[j] + pb_xxxz[j]) * s_0_0[j];

                t_0_xxyy[j] = (0.25 * fx[j] * fx[j] + 0.5 * pb_xx[j] * fx[j] + 0.5 * fx[j] * pb_yy[j] + pb_xxyy[j]) * s_0_0[j];

                t_0_xxyz[j] = (0.5 * fx[j] * pb_yz[j] + pb_xxyz[j]) * s_0_0[j];

                t_0_xxzz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pb_xx[j] * fx[j] + 0.5 * fx[j] * pb_zz[j] + pb_xxzz[j]) * s_0_0[j];

                t_0_xyyy[j] = (1.5 * pb_xy[j] * fx[j] + pb_xyyy[j]) * s_0_0[j];

                t_0_xyyz[j] = (0.5 * pb_xz[j] * fx[j] + pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (8,15)

            #pragma omp simd aligned(fx, pb_xy, pb_xyzz, pb_xz, pb_xzzz, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, \
                                     pb_yzzz, pb_zz, pb_zzzz, s_0_0, t_0_xyzz, t_0_xzzz, t_0_yyyy, t_0_yyyz, t_0_yyzz, \
                                     t_0_yzzz, t_0_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_0_xyzz[j] = (0.5 * pb_xy[j] * fx[j] + pb_xyzz[j]) * s_0_0[j];

                t_0_xzzz[j] = (1.5 * pb_xz[j] * fx[j] + pb_xzzz[j]) * s_0_0[j];

                t_0_yyyy[j] = (0.75 * fx[j] * fx[j] + 3.0 * pb_yy[j] * fx[j] + pb_yyyy[j]) * s_0_0[j];

                t_0_yyyz[j] = (1.5 * pb_yz[j] * fx[j] + pb_yyyz[j]) * s_0_0[j];

                t_0_yyzz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pb_yy[j] * fx[j] + 0.5 * fx[j] * pb_zz[j] + pb_yyzz[j]) * s_0_0[j];

                t_0_yzzz[j] = (1.5 * pb_yz[j] * fx[j] + pb_yzzz[j]) * s_0_0[j];

                t_0_zzzz[j] = (0.75 * fx[j] * fx[j] + 3.0 * pb_zz[j] * fx[j] + pb_zzzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForPS(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_0 = primBuffer.data(3 * idx);

            auto t_y_0 = primBuffer.data(3 * idx + 1);

            auto t_z_0 = primBuffer.data(3 * idx + 2);

            // Batch of Integrals (0) = (0,3)

            #pragma omp simd aligned(pa_x, pa_y, pa_z, s_0_0, t_x_0, t_y_0, t_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_x_0[j] = pa_x[j] * s_0_0[j];

                t_y_0[j] = pa_y[j] * s_0_0[j];

                t_z_0[j] = pa_z[j] * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForPP(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_x = primBuffer.data(9 * idx);

            auto t_x_y = primBuffer.data(9 * idx + 1);

            auto t_x_z = primBuffer.data(9 * idx + 2);

            auto t_y_x = primBuffer.data(9 * idx + 3);

            auto t_y_y = primBuffer.data(9 * idx + 4);

            auto t_y_z = primBuffer.data(9 * idx + 5);

            auto t_z_x = primBuffer.data(9 * idx + 6);

            auto t_z_y = primBuffer.data(9 * idx + 7);

            auto t_z_z = primBuffer.data(9 * idx + 8);

            // Batch of Integrals (0) = (0,9)

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pb_x, pb_y, pb_z, s_0_0, t_x_x, t_x_y, t_x_z, t_y_x, t_y_y, \
                                     t_y_z, t_z_x, t_z_y, t_z_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_x_x[j] = (0.5 * fx[j] + pa_x[j] * pb_x[j]) * s_0_0[j];

                t_x_y[j] = pa_x[j] * pb_y[j] * s_0_0[j];

                t_x_z[j] = pa_x[j] * pb_z[j] * s_0_0[j];

                t_y_x[j] = pa_y[j] * pb_x[j] * s_0_0[j];

                t_y_y[j] = (0.5 * fx[j] + pa_y[j] * pb_y[j]) * s_0_0[j];

                t_y_z[j] = pa_y[j] * pb_z[j] * s_0_0[j];

                t_z_x[j] = pa_z[j] * pb_x[j] * s_0_0[j];

                t_z_y[j] = pa_z[j] * pb_y[j] * s_0_0[j];

                t_z_z[j] = (0.5 * fx[j] + pa_z[j] * pb_z[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForPD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_xx = primBuffer.data(18 * idx);

            auto t_x_xy = primBuffer.data(18 * idx + 1);

            auto t_x_xz = primBuffer.data(18 * idx + 2);

            auto t_x_yy = primBuffer.data(18 * idx + 3);

            auto t_x_yz = primBuffer.data(18 * idx + 4);

            auto t_x_zz = primBuffer.data(18 * idx + 5);

            auto t_y_xx = primBuffer.data(18 * idx + 6);

            auto t_y_xy = primBuffer.data(18 * idx + 7);

            auto t_y_xz = primBuffer.data(18 * idx + 8);

            auto t_y_yy = primBuffer.data(18 * idx + 9);

            auto t_y_yz = primBuffer.data(18 * idx + 10);

            auto t_y_zz = primBuffer.data(18 * idx + 11);

            auto t_z_xx = primBuffer.data(18 * idx + 12);

            auto t_z_xy = primBuffer.data(18 * idx + 13);

            auto t_z_xz = primBuffer.data(18 * idx + 14);

            auto t_z_yy = primBuffer.data(18 * idx + 15);

            auto t_z_yz = primBuffer.data(18 * idx + 16);

            auto t_z_zz = primBuffer.data(18 * idx + 17);

            // Batch of Integrals (0) = (0,9)

            #pragma omp simd aligned(fx, pa_x, pa_y, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, \
                                     t_x_xx, t_x_xy, t_x_xz, t_x_yy, t_x_yz, t_x_zz, t_y_xx, t_y_xy, t_y_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_x_xx[j] = (0.5 * pa_x[j] * fx[j] + fx[j] * pb_x[j] + pa_x[j] * pb_xx[j]) * s_0_0[j];

                t_x_xy[j] = (0.5 * fx[j] * pb_y[j] + pa_x[j] * pb_xy[j]) * s_0_0[j];

                t_x_xz[j] = (0.5 * fx[j] * pb_z[j] + pa_x[j] * pb_xz[j]) * s_0_0[j];

                t_x_yy[j] = (0.5 * pa_x[j] * fx[j] + pa_x[j] * pb_yy[j]) * s_0_0[j];

                t_x_yz[j] = pa_x[j] * pb_yz[j] * s_0_0[j];

                t_x_zz[j] = (0.5 * pa_x[j] * fx[j] + pa_x[j] * pb_zz[j]) * s_0_0[j];

                t_y_xx[j] = (0.5 * pa_y[j] * fx[j] + pa_y[j] * pb_xx[j]) * s_0_0[j];

                t_y_xy[j] = (0.5 * fx[j] * pb_x[j] + pa_y[j] * pb_xy[j]) * s_0_0[j];

                t_y_xz[j] = pa_y[j] * pb_xz[j] * s_0_0[j];
            }

            // Batch of Integrals (1) = (9,18)

            #pragma omp simd aligned(fx, pa_y, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, \
                                     t_y_yy, t_y_yz, t_y_zz, t_z_xx, t_z_xy, t_z_xz, t_z_yy, t_z_yz, t_z_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_y_yy[j] = (0.5 * pa_y[j] * fx[j] + fx[j] * pb_y[j] + pa_y[j] * pb_yy[j]) * s_0_0[j];

                t_y_yz[j] = (0.5 * fx[j] * pb_z[j] + pa_y[j] * pb_yz[j]) * s_0_0[j];

                t_y_zz[j] = (0.5 * pa_y[j] * fx[j] + pa_y[j] * pb_zz[j]) * s_0_0[j];

                t_z_xx[j] = (0.5 * pa_z[j] * fx[j] + pa_z[j] * pb_xx[j]) * s_0_0[j];

                t_z_xy[j] = pa_z[j] * pb_xy[j] * s_0_0[j];

                t_z_xz[j] = (0.5 * fx[j] * pb_x[j] + pa_z[j] * pb_xz[j]) * s_0_0[j];

                t_z_yy[j] = (0.5 * pa_z[j] * fx[j] + pa_z[j] * pb_yy[j]) * s_0_0[j];

                t_z_yz[j] = (0.5 * fx[j] * pb_y[j] + pa_z[j] * pb_yz[j]) * s_0_0[j];

                t_z_zz[j] = (0.5 * pa_z[j] * fx[j] + fx[j] * pb_z[j] + pa_z[j] * pb_zz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForPF(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_xxx = primBuffer.data(30 * idx);

            auto t_x_xxy = primBuffer.data(30 * idx + 1);

            auto t_x_xxz = primBuffer.data(30 * idx + 2);

            auto t_x_xyy = primBuffer.data(30 * idx + 3);

            auto t_x_xyz = primBuffer.data(30 * idx + 4);

            auto t_x_xzz = primBuffer.data(30 * idx + 5);

            auto t_x_yyy = primBuffer.data(30 * idx + 6);

            auto t_x_yyz = primBuffer.data(30 * idx + 7);

            auto t_x_yzz = primBuffer.data(30 * idx + 8);

            auto t_x_zzz = primBuffer.data(30 * idx + 9);

            auto t_y_xxx = primBuffer.data(30 * idx + 10);

            auto t_y_xxy = primBuffer.data(30 * idx + 11);

            auto t_y_xxz = primBuffer.data(30 * idx + 12);

            auto t_y_xyy = primBuffer.data(30 * idx + 13);

            auto t_y_xyz = primBuffer.data(30 * idx + 14);

            auto t_y_xzz = primBuffer.data(30 * idx + 15);

            auto t_y_yyy = primBuffer.data(30 * idx + 16);

            auto t_y_yyz = primBuffer.data(30 * idx + 17);

            auto t_y_yzz = primBuffer.data(30 * idx + 18);

            auto t_y_zzz = primBuffer.data(30 * idx + 19);

            auto t_z_xxx = primBuffer.data(30 * idx + 20);

            auto t_z_xxy = primBuffer.data(30 * idx + 21);

            auto t_z_xxz = primBuffer.data(30 * idx + 22);

            auto t_z_xyy = primBuffer.data(30 * idx + 23);

            auto t_z_xyz = primBuffer.data(30 * idx + 24);

            auto t_z_xzz = primBuffer.data(30 * idx + 25);

            auto t_z_yyy = primBuffer.data(30 * idx + 26);

            auto t_z_yyz = primBuffer.data(30 * idx + 27);

            auto t_z_yzz = primBuffer.data(30 * idx + 28);

            auto t_z_zzz = primBuffer.data(30 * idx + 29);

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fx, pa_x, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_x_xxx, \
                                     t_x_xxy, t_x_xxz, t_x_xyy, t_x_xyz, t_x_xzz, t_x_yyy, t_x_yyz, t_x_yzz, t_x_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_x_xxx[j] = (0.75 * fx[j] * fx[j] + 1.5 * pa_x[j] * pb_x[j] * fx[j] + 1.5 * fx[j] * pb_xx[j] +

                             pa_x[j] * pb_xxx[j]) * s_0_0[j];

                t_x_xxy[j] = (0.5 * pa_x[j] * fx[j] * pb_y[j] + fx[j] * pb_xy[j] + pa_x[j] * pb_xxy[j]) * s_0_0[j];

                t_x_xxz[j] = (0.5 * pa_x[j] * fx[j] * pb_z[j] + fx[j] * pb_xz[j] + pa_x[j] * pb_xxz[j]) * s_0_0[j];

                t_x_xyy[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_x[j] * pb_x[j] * fx[j] + 0.5 * fx[j] * pb_yy[j] +

                             pa_x[j] * pb_xyy[j]) * s_0_0[j];

                t_x_xyz[j] = (0.5 * fx[j] * pb_yz[j] + pa_x[j] * pb_xyz[j]) * s_0_0[j];

                t_x_xzz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_x[j] * pb_x[j] * fx[j] + 0.5 * fx[j] * pb_zz[j] +

                             pa_x[j] * pb_xzz[j]) * s_0_0[j];

                t_x_yyy[j] = (1.5 * pa_x[j] * pb_y[j] * fx[j] + pa_x[j] * pb_yyy[j]) * s_0_0[j];

                t_x_yyz[j] = (0.5 * pa_x[j] * fx[j] * pb_z[j] + pa_x[j] * pb_yyz[j]) * s_0_0[j];

                t_x_yzz[j] = (0.5 * pa_x[j] * pb_y[j] * fx[j] + pa_x[j] * pb_yzz[j]) * s_0_0[j];

                t_x_zzz[j] = (1.5 * pa_x[j] * pb_z[j] * fx[j] + pa_x[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fx, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_y_xxx, \
                                     t_y_xxy, t_y_xxz, t_y_xyy, t_y_xyz, t_y_xzz, t_y_yyy, t_y_yyz, t_y_yzz, t_y_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_y_xxx[j] = (1.5 * pa_y[j] * pb_x[j] * fx[j] + pa_y[j] * pb_xxx[j]) * s_0_0[j];

                t_y_xxy[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pb_xx[j] +

                             pa_y[j] * pb_xxy[j]) * s_0_0[j];

                t_y_xxz[j] = (0.5 * pa_y[j] * fx[j] * pb_z[j] + pa_y[j] * pb_xxz[j]) * s_0_0[j];

                t_y_xyy[j] = (0.5 * pa_y[j] * pb_x[j] * fx[j] + fx[j] * pb_xy[j] + pa_y[j] * pb_xyy[j]) * s_0_0[j];

                t_y_xyz[j] = (0.5 * fx[j] * pb_xz[j] + pa_y[j] * pb_xyz[j]) * s_0_0[j];

                t_y_xzz[j] = (0.5 * pa_y[j] * pb_x[j] * fx[j] + pa_y[j] * pb_xzz[j]) * s_0_0[j];

                t_y_yyy[j] = (0.75 * fx[j] * fx[j] + 1.5 * pa_y[j] * pb_y[j] * fx[j] + 1.5 * fx[j] * pb_yy[j] +

                             pa_y[j] * pb_yyy[j]) * s_0_0[j];

                t_y_yyz[j] = (0.5 * pa_y[j] * fx[j] * pb_z[j] + fx[j] * pb_yz[j] + pa_y[j] * pb_yyz[j]) * s_0_0[j];

                t_y_yzz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_y[j] * pb_y[j] * fx[j] + 0.5 * fx[j] * pb_zz[j] +

                             pa_y[j] * pb_yzz[j]) * s_0_0[j];

                t_y_zzz[j] = (1.5 * pa_y[j] * pb_z[j] * fx[j] + pa_y[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fx, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_z_xxx, \
                                     t_z_xxy, t_z_xxz, t_z_xyy, t_z_xyz, t_z_xzz, t_z_yyy, t_z_yyz, t_z_yzz, t_z_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_z_xxx[j] = (1.5 * pa_z[j] * pb_x[j] * fx[j] + pa_z[j] * pb_xxx[j]) * s_0_0[j];

                t_z_xxy[j] = (0.5 * pa_z[j] * fx[j] * pb_y[j] + pa_z[j] * pb_xxy[j]) * s_0_0[j];

                t_z_xxz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_z[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pb_xx[j] +

                             pa_z[j] * pb_xxz[j]) * s_0_0[j];

                t_z_xyy[j] = (0.5 * pa_z[j] * pb_x[j] * fx[j] + pa_z[j] * pb_xyy[j]) * s_0_0[j];

                t_z_xyz[j] = (0.5 * fx[j] * pb_xy[j] + pa_z[j] * pb_xyz[j]) * s_0_0[j];

                t_z_xzz[j] = (0.5 * pa_z[j] * pb_x[j] * fx[j] + fx[j] * pb_xz[j] + pa_z[j] * pb_xzz[j]) * s_0_0[j];

                t_z_yyy[j] = (1.5 * pa_z[j] * pb_y[j] * fx[j] + pa_z[j] * pb_yyy[j]) * s_0_0[j];

                t_z_yyz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_z[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pb_yy[j] +

                             pa_z[j] * pb_yyz[j]) * s_0_0[j];

                t_z_yzz[j] = (0.5 * pa_z[j] * pb_y[j] * fx[j] + fx[j] * pb_yz[j] + pa_z[j] * pb_yzz[j]) * s_0_0[j];

                t_z_zzz[j] = (0.75 * fx[j] * fx[j] + 1.5 * pa_z[j] * pb_z[j] * fx[j] + 1.5 * fx[j] * pb_zz[j] +

                             pa_z[j] * pb_zzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForPG(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_x_xxxx = primBuffer.data(45 * idx);

            auto t_x_xxxy = primBuffer.data(45 * idx + 1);

            auto t_x_xxxz = primBuffer.data(45 * idx + 2);

            auto t_x_xxyy = primBuffer.data(45 * idx + 3);

            auto t_x_xxyz = primBuffer.data(45 * idx + 4);

            auto t_x_xxzz = primBuffer.data(45 * idx + 5);

            auto t_x_xyyy = primBuffer.data(45 * idx + 6);

            auto t_x_xyyz = primBuffer.data(45 * idx + 7);

            auto t_x_xyzz = primBuffer.data(45 * idx + 8);

            auto t_x_xzzz = primBuffer.data(45 * idx + 9);

            auto t_x_yyyy = primBuffer.data(45 * idx + 10);

            auto t_x_yyyz = primBuffer.data(45 * idx + 11);

            auto t_x_yyzz = primBuffer.data(45 * idx + 12);

            auto t_x_yzzz = primBuffer.data(45 * idx + 13);

            auto t_x_zzzz = primBuffer.data(45 * idx + 14);

            auto t_y_xxxx = primBuffer.data(45 * idx + 15);

            auto t_y_xxxy = primBuffer.data(45 * idx + 16);

            auto t_y_xxxz = primBuffer.data(45 * idx + 17);

            auto t_y_xxyy = primBuffer.data(45 * idx + 18);

            auto t_y_xxyz = primBuffer.data(45 * idx + 19);

            auto t_y_xxzz = primBuffer.data(45 * idx + 20);

            auto t_y_xyyy = primBuffer.data(45 * idx + 21);

            auto t_y_xyyz = primBuffer.data(45 * idx + 22);

            auto t_y_xyzz = primBuffer.data(45 * idx + 23);

            auto t_y_xzzz = primBuffer.data(45 * idx + 24);

            auto t_y_yyyy = primBuffer.data(45 * idx + 25);

            auto t_y_yyyz = primBuffer.data(45 * idx + 26);

            auto t_y_yyzz = primBuffer.data(45 * idx + 27);

            auto t_y_yzzz = primBuffer.data(45 * idx + 28);

            auto t_y_zzzz = primBuffer.data(45 * idx + 29);

            auto t_z_xxxx = primBuffer.data(45 * idx + 30);

            auto t_z_xxxy = primBuffer.data(45 * idx + 31);

            auto t_z_xxxz = primBuffer.data(45 * idx + 32);

            auto t_z_xxyy = primBuffer.data(45 * idx + 33);

            auto t_z_xxyz = primBuffer.data(45 * idx + 34);

            auto t_z_xxzz = primBuffer.data(45 * idx + 35);

            auto t_z_xyyy = primBuffer.data(45 * idx + 36);

            auto t_z_xyyz = primBuffer.data(45 * idx + 37);

            auto t_z_xyzz = primBuffer.data(45 * idx + 38);

            auto t_z_xzzz = primBuffer.data(45 * idx + 39);

            auto t_z_yyyy = primBuffer.data(45 * idx + 40);

            auto t_z_yyyz = primBuffer.data(45 * idx + 41);

            auto t_z_yyzz = primBuffer.data(45 * idx + 42);

            auto t_z_yzzz = primBuffer.data(45 * idx + 43);

            auto t_z_zzzz = primBuffer.data(45 * idx + 44);

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fx, pa_x, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_x_xxxx, \
                                     t_x_xxxy, t_x_xxxz, t_x_xxyy, t_x_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_x_xxxx[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pb_x[j] +

                              3.0 * pa_x[j] * pb_xx[j] * fx[j] + 2.0 * fx[j] * pb_xxx[j] + pa_x[j] * pb_xxxx[j]) * s_0_0[j];

                t_x_xxxy[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_x[j] * pb_xy[j] * fx[j] +

                              1.5 * fx[j] * pb_xxy[j] + pa_x[j] * pb_xxxy[j]) * s_0_0[j];

                t_x_xxxz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_x[j] * pb_xz[j] * fx[j] +

                              1.5 * fx[j] * pb_xxz[j] + pa_x[j] * pb_xxxz[j]) * s_0_0[j];

                t_x_xxyy[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_x[j] * pb_xx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_yy[j] + fx[j] * pb_xyy[j] + pa_x[j] * pb_xxyy[j]) * s_0_0[j];

                t_x_xxyz[j] = (0.5 * pa_x[j] * fx[j] * pb_yz[j] + fx[j] * pb_xyz[j] + pa_x[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fx, pa_x, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, pb_xyzz, pb_xz, pb_xzz, \
                                     pb_xzzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_x_xxzz, t_x_xyyy, \
                                     t_x_xyyz, t_x_xyzz, t_x_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_x_xxzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_x[j] * pb_xx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_zz[j] + fx[j] * pb_xzz[j] + pa_x[j] * pb_xxzz[j]) * s_0_0[j];

                t_x_xyyy[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_x[j] * pb_xy[j] * fx[j] +

                              0.5 * fx[j] * pb_yyy[j] + pa_x[j] * pb_xyyy[j]) * s_0_0[j];

                t_x_xyyz[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * pb_xz[j] * fx[j] +

                              0.5 * fx[j] * pb_yyz[j] + pa_x[j] * pb_xyyz[j]) * s_0_0[j];

                t_x_xyzz[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * pb_xy[j] * fx[j] +

                              0.5 * fx[j] * pb_yzz[j] + pa_x[j] * pb_xyzz[j]) * s_0_0[j];

                t_x_xzzz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_x[j] * pb_xz[j] * fx[j] +

                              0.5 * fx[j] * pb_zzz[j] + pa_x[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fx, pa_x, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_zz, pb_zzzz, \
                                     s_0_0, t_x_yyyy, t_x_yyyz, t_x_yyzz, t_x_yzzz, t_x_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_x_yyyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 3.0 * pa_x[j] * pb_yy[j] * fx[j] +

                              pa_x[j] * pb_yyyy[j]) * s_0_0[j];

                t_x_yyyz[j] = (1.5 * pa_x[j] * pb_yz[j] * fx[j] + pa_x[j] * pb_yyyz[j]) * s_0_0[j];

                t_x_yyzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * pb_yy[j] * fx[j] +

                              0.5 * pa_x[j] * fx[j] * pb_zz[j] + pa_x[j] * pb_yyzz[j]) * s_0_0[j];

                t_x_yzzz[j] = (1.5 * pa_x[j] * pb_yz[j] * fx[j] + pa_x[j] * pb_yzzz[j]) * s_0_0[j];

                t_x_zzzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 3.0 * pa_x[j] * pb_zz[j] * fx[j] +

                              pa_x[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fx, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_y_xxxx, t_y_xxxy, \
                                     t_y_xxxz, t_y_xxyy, t_y_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_y_xxxx[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 3.0 * pa_y[j] * pb_xx[j] * fx[j] +

                              pa_y[j] * pb_xxxx[j]) * s_0_0[j];

                t_y_xxxy[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_y[j] * pb_xy[j] * fx[j] +

                              0.5 * fx[j] * pb_xxx[j] + pa_y[j] * pb_xxxy[j]) * s_0_0[j];

                t_y_xxxz[j] = (1.5 * pa_y[j] * pb_xz[j] * fx[j] + pa_y[j] * pb_xxxz[j]) * s_0_0[j];

                t_y_xxyy[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_y[j] * pb_xx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_yy[j] + fx[j] * pb_xxy[j] + pa_y[j] * pb_xxyy[j]) * s_0_0[j];

                t_y_xxyz[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_y[j] * fx[j] * pb_yz[j] +

                              0.5 * fx[j] * pb_xxz[j] + pa_y[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fx, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, \
                                     pb_xz, pb_xzz, pb_xzzz, pb_zz, s_0_0, t_y_xxzz, t_y_xyyy, t_y_xyyz, t_y_xyzz, \
                                     t_y_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_y_xxzz[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * pa_y[j] * pb_xx[j] * fx[j] +

                              0.5 * pa_y[j] * fx[j] * pb_zz[j] + pa_y[j] * pb_xxzz[j]) * s_0_0[j];

                t_y_xyyy[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_y[j] * pb_xy[j] * fx[j] +

                              1.5 * fx[j] * pb_xyy[j] + pa_y[j] * pb_xyyy[j]) * s_0_0[j];

                t_y_xyyz[j] = (0.5 * pa_y[j] * pb_xz[j] * fx[j] + fx[j] * pb_xyz[j] + pa_y[j] * pb_xyyz[j]) * s_0_0[j];

                t_y_xyzz[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_y[j] * pb_xy[j] * fx[j] +

                              0.5 * fx[j] * pb_xzz[j] + pa_y[j] * pb_xyzz[j]) * s_0_0[j];

                t_y_xzzz[j] = (1.5 * pa_y[j] * pb_xz[j] * fx[j] + pa_y[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fx, pa_y, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_y_yyyy, t_y_yyyz, t_y_yyzz, t_y_yzzz, \
                                     t_y_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_y_yyyy[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pb_y[j] +

                              3.0 * pa_y[j] * pb_yy[j] * fx[j] + 2.0 * fx[j] * pb_yyy[j] + pa_y[j] * pb_yyyy[j]) * s_0_0[j];

                t_y_yyyz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_y[j] * pb_yz[j] * fx[j] +

                              1.5 * fx[j] * pb_yyz[j] + pa_y[j] * pb_yyyz[j]) * s_0_0[j];

                t_y_yyzz[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_y[j] * pb_yy[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_zz[j] + fx[j] * pb_yzz[j] + pa_y[j] * pb_yyzz[j]) * s_0_0[j];

                t_y_yzzz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_y[j] * pb_yz[j] * fx[j] +

                              0.5 * fx[j] * pb_zzz[j] + pa_y[j] * pb_yzzz[j]) * s_0_0[j];

                t_y_zzzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 3.0 * pa_y[j] * pb_zz[j] * fx[j] +

                              pa_y[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fx, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, s_0_0, t_z_xxxx, t_z_xxxy, t_z_xxxz, \
                                     t_z_xxyy, t_z_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_z_xxxx[j] = (0.75 * pa_z[j] * fx[j] * fx[j] + 3.0 * pa_z[j] * pb_xx[j] * fx[j] +

                              pa_z[j] * pb_xxxx[j]) * s_0_0[j];

                t_z_xxxy[j] = (1.5 * pa_z[j] * pb_xy[j] * fx[j] + pa_z[j] * pb_xxxy[j]) * s_0_0[j];

                t_z_xxxz[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_z[j] * pb_xz[j] * fx[j] +

                              0.5 * fx[j] * pb_xxx[j] + pa_z[j] * pb_xxxz[j]) * s_0_0[j];

                t_z_xxyy[j] = (0.25 * pa_z[j] * fx[j] * fx[j] + 0.5 * pa_z[j] * pb_xx[j] * fx[j] +

                              0.5 * pa_z[j] * fx[j] * pb_yy[j] + pa_z[j] * pb_xxyy[j]) * s_0_0[j];

                t_z_xxyz[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_z[j] * fx[j] * pb_yz[j] +

                              0.5 * fx[j] * pb_xxy[j] + pa_z[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fx, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, s_0_0, t_z_xxzz, t_z_xyyy, t_z_xyyz, \
                                     t_z_xyzz, t_z_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_z_xxzz[j] = (0.25 * pa_z[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_z[j] * pb_xx[j] * fx[j] + 0.5 * pa_z[j] * fx[j] * pb_zz[j] + fx[j] * pb_xxz[j] + pa_z[j] * pb_xxzz[j]) * s_0_0[j];

                t_z_xyyy[j] = (1.5 * pa_z[j] * pb_xy[j] * fx[j] + pa_z[j] * pb_xyyy[j]) * s_0_0[j];

                t_z_xyyz[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_z[j] * pb_xz[j] * fx[j] +

                              0.5 * fx[j] * pb_xyy[j] + pa_z[j] * pb_xyyz[j]) * s_0_0[j];

                t_z_xyzz[j] = (0.5 * pa_z[j] * pb_xy[j] * fx[j] + fx[j] * pb_xyz[j] + pa_z[j] * pb_xyzz[j]) * s_0_0[j];

                t_z_xzzz[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_z[j] * pb_xz[j] * fx[j] +

                              1.5 * fx[j] * pb_xzz[j] + pa_z[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fx, pa_z, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_z_yyyy, t_z_yyyz, t_z_yyzz, t_z_yzzz, \
                                     t_z_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_z_yyyy[j] = (0.75 * pa_z[j] * fx[j] * fx[j] + 3.0 * pa_z[j] * pb_yy[j] * fx[j] +

                              pa_z[j] * pb_yyyy[j]) * s_0_0[j];

                t_z_yyyz[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_z[j] * pb_yz[j] * fx[j] +

                              0.5 * fx[j] * pb_yyy[j] + pa_z[j] * pb_yyyz[j]) * s_0_0[j];

                t_z_yyzz[j] = (0.25 * pa_z[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_z[j] * pb_yy[j] * fx[j] + 0.5 * pa_z[j] * fx[j] * pb_zz[j] + fx[j] * pb_yyz[j] + pa_z[j] * pb_yyzz[j]) * s_0_0[j];

                t_z_yzzz[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_z[j] * pb_yz[j] * fx[j] +

                              1.5 * fx[j] * pb_yzz[j] + pa_z[j] * pb_yzzz[j]) * s_0_0[j];

                t_z_zzzz[j] = (0.75 * pa_z[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pb_z[j] +

                              3.0 * pa_z[j] * pb_zz[j] * fx[j] + 2.0 * fx[j] * pb_zzz[j] + pa_z[j] * pb_zzzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForDS(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            auto pa_xz = paDistances.data(9 * idx + 5);

            auto pa_yy = paDistances.data(9 * idx + 6);

            auto pa_yz = paDistances.data(9 * idx + 7);

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_0 = primBuffer.data(6 * idx);

            auto t_xy_0 = primBuffer.data(6 * idx + 1);

            auto t_xz_0 = primBuffer.data(6 * idx + 2);

            auto t_yy_0 = primBuffer.data(6 * idx + 3);

            auto t_yz_0 = primBuffer.data(6 * idx + 4);

            auto t_zz_0 = primBuffer.data(6 * idx + 5);

            // Batch of Integrals (0) = (0,6)

            #pragma omp simd aligned(fx, pa_xx, pa_xy, pa_xz, pa_yy, pa_yz, pa_zz, s_0_0, t_xx_0, t_xy_0, t_xz_0, \
                                     t_yy_0, t_yz_0, t_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_0[j] = (0.5 * fx[j] + pa_xx[j]) * s_0_0[j];

                t_xy_0[j] = pa_xy[j] * s_0_0[j];

                t_xz_0[j] = pa_xz[j] * s_0_0[j];

                t_yy_0[j] = (0.5 * fx[j] + pa_yy[j]) * s_0_0[j];

                t_yz_0[j] = pa_yz[j] * s_0_0[j];

                t_zz_0[j] = (0.5 * fx[j] + pa_zz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForDP(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(9 * idx);

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            auto pa_xz = paDistances.data(9 * idx + 5);

            auto pa_yy = paDistances.data(9 * idx + 6);

            auto pa_yz = paDistances.data(9 * idx + 7);

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_x = primBuffer.data(18 * idx);

            auto t_xx_y = primBuffer.data(18 * idx + 1);

            auto t_xx_z = primBuffer.data(18 * idx + 2);

            auto t_xy_x = primBuffer.data(18 * idx + 3);

            auto t_xy_y = primBuffer.data(18 * idx + 4);

            auto t_xy_z = primBuffer.data(18 * idx + 5);

            auto t_xz_x = primBuffer.data(18 * idx + 6);

            auto t_xz_y = primBuffer.data(18 * idx + 7);

            auto t_xz_z = primBuffer.data(18 * idx + 8);

            auto t_yy_x = primBuffer.data(18 * idx + 9);

            auto t_yy_y = primBuffer.data(18 * idx + 10);

            auto t_yy_z = primBuffer.data(18 * idx + 11);

            auto t_yz_x = primBuffer.data(18 * idx + 12);

            auto t_yz_y = primBuffer.data(18 * idx + 13);

            auto t_yz_z = primBuffer.data(18 * idx + 14);

            auto t_zz_x = primBuffer.data(18 * idx + 15);

            auto t_zz_y = primBuffer.data(18 * idx + 16);

            auto t_zz_z = primBuffer.data(18 * idx + 17);

            // Batch of Integrals (0) = (0,9)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xy, pa_xz, pa_y, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xx_x, t_xx_y, \
                                     t_xx_z, t_xy_x, t_xy_y, t_xy_z, t_xz_x, t_xz_y, t_xz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_x[j] = (pa_x[j] * fx[j] + 0.5 * fx[j] * pb_x[j] + pa_xx[j] * pb_x[j]) * s_0_0[j];

                t_xx_y[j] = (0.5 * fx[j] * pb_y[j] + pa_xx[j] * pb_y[j]) * s_0_0[j];

                t_xx_z[j] = (0.5 * fx[j] * pb_z[j] + pa_xx[j] * pb_z[j]) * s_0_0[j];

                t_xy_x[j] = (0.5 * fx[j] * pa_y[j] + pa_xy[j] * pb_x[j]) * s_0_0[j];

                t_xy_y[j] = (0.5 * pa_x[j] * fx[j] + pa_xy[j] * pb_y[j]) * s_0_0[j];

                t_xy_z[j] = pa_xy[j] * pb_z[j] * s_0_0[j];

                t_xz_x[j] = (0.5 * fx[j] * pa_z[j] + pa_xz[j] * pb_x[j]) * s_0_0[j];

                t_xz_y[j] = pa_xz[j] * pb_y[j] * s_0_0[j];

                t_xz_z[j] = (0.5 * pa_x[j] * fx[j] + pa_xz[j] * pb_z[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (9,18)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yz, pa_z, pa_zz, pb_x, pb_y, pb_z, s_0_0, t_yy_x, t_yy_y, t_yy_z, \
                                     t_yz_x, t_yz_y, t_yz_z, t_zz_x, t_zz_y, t_zz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yy_x[j] = (0.5 * fx[j] * pb_x[j] + pa_yy[j] * pb_x[j]) * s_0_0[j];

                t_yy_y[j] = (pa_y[j] * fx[j] + 0.5 * fx[j] * pb_y[j] + pa_yy[j] * pb_y[j]) * s_0_0[j];

                t_yy_z[j] = (0.5 * fx[j] * pb_z[j] + pa_yy[j] * pb_z[j]) * s_0_0[j];

                t_yz_x[j] = pa_yz[j] * pb_x[j] * s_0_0[j];

                t_yz_y[j] = (0.5 * fx[j] * pa_z[j] + pa_yz[j] * pb_y[j]) * s_0_0[j];

                t_yz_z[j] = (0.5 * pa_y[j] * fx[j] + pa_yz[j] * pb_z[j]) * s_0_0[j];

                t_zz_x[j] = (0.5 * fx[j] * pb_x[j] + pa_zz[j] * pb_x[j]) * s_0_0[j];

                t_zz_y[j] = (0.5 * fx[j] * pb_y[j] + pa_zz[j] * pb_y[j]) * s_0_0[j];

                t_zz_z[j] = (pa_z[j] * fx[j] + 0.5 * fx[j] * pb_z[j] + pa_zz[j] * pb_z[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForDD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(9 * idx);

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            auto pa_xz = paDistances.data(9 * idx + 5);

            auto pa_yy = paDistances.data(9 * idx + 6);

            auto pa_yz = paDistances.data(9 * idx + 7);

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_xx = primBuffer.data(36 * idx);

            auto t_xx_xy = primBuffer.data(36 * idx + 1);

            auto t_xx_xz = primBuffer.data(36 * idx + 2);

            auto t_xx_yy = primBuffer.data(36 * idx + 3);

            auto t_xx_yz = primBuffer.data(36 * idx + 4);

            auto t_xx_zz = primBuffer.data(36 * idx + 5);

            auto t_xy_xx = primBuffer.data(36 * idx + 6);

            auto t_xy_xy = primBuffer.data(36 * idx + 7);

            auto t_xy_xz = primBuffer.data(36 * idx + 8);

            auto t_xy_yy = primBuffer.data(36 * idx + 9);

            auto t_xy_yz = primBuffer.data(36 * idx + 10);

            auto t_xy_zz = primBuffer.data(36 * idx + 11);

            auto t_xz_xx = primBuffer.data(36 * idx + 12);

            auto t_xz_xy = primBuffer.data(36 * idx + 13);

            auto t_xz_xz = primBuffer.data(36 * idx + 14);

            auto t_xz_yy = primBuffer.data(36 * idx + 15);

            auto t_xz_yz = primBuffer.data(36 * idx + 16);

            auto t_xz_zz = primBuffer.data(36 * idx + 17);

            auto t_yy_xx = primBuffer.data(36 * idx + 18);

            auto t_yy_xy = primBuffer.data(36 * idx + 19);

            auto t_yy_xz = primBuffer.data(36 * idx + 20);

            auto t_yy_yy = primBuffer.data(36 * idx + 21);

            auto t_yy_yz = primBuffer.data(36 * idx + 22);

            auto t_yy_zz = primBuffer.data(36 * idx + 23);

            auto t_yz_xx = primBuffer.data(36 * idx + 24);

            auto t_yz_xy = primBuffer.data(36 * idx + 25);

            auto t_yz_xz = primBuffer.data(36 * idx + 26);

            auto t_yz_yy = primBuffer.data(36 * idx + 27);

            auto t_yz_yz = primBuffer.data(36 * idx + 28);

            auto t_yz_zz = primBuffer.data(36 * idx + 29);

            auto t_zz_xx = primBuffer.data(36 * idx + 30);

            auto t_zz_xy = primBuffer.data(36 * idx + 31);

            auto t_zz_xz = primBuffer.data(36 * idx + 32);

            auto t_zz_yy = primBuffer.data(36 * idx + 33);

            auto t_zz_yz = primBuffer.data(36 * idx + 34);

            auto t_zz_zz = primBuffer.data(36 * idx + 35);

            // Batch of Integrals (0) = (0,9)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xy, pa_y, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_xx_xx, t_xx_xy, t_xx_xz, t_xx_yy, t_xx_yz, t_xx_zz, t_xy_xx, \
                                     t_xy_xy, t_xy_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_xx[j] = (0.75 * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] + 2.0 * pa_x[j] * fx[j] * pb_x[j] +

                             0.5 * fx[j] * pb_xx[j] + pa_xx[j] * pb_xx[j]) * s_0_0[j];

                t_xx_xy[j] = (pa_x[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pb_xy[j] + pa_xx[j] * pb_xy[j]) * s_0_0[j];

                t_xx_xz[j] = (pa_x[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pb_xz[j] + pa_xx[j] * pb_xz[j]) * s_0_0[j];

                t_xx_yy[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] + 0.5 * fx[j] * pb_yy[j] +

                             pa_xx[j] * pb_yy[j]) * s_0_0[j];

                t_xx_yz[j] = (0.5 * fx[j] * pb_yz[j] + pa_xx[j] * pb_yz[j]) * s_0_0[j];

                t_xx_zz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] + 0.5 * fx[j] * pb_zz[j] +

                             pa_xx[j] * pb_zz[j]) * s_0_0[j];

                t_xy_xx[j] = (0.5 * pa_xy[j] * fx[j] + fx[j] * pa_y[j] * pb_x[j] + pa_xy[j] * pb_xx[j]) * s_0_0[j];

                t_xy_xy[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_x[j] +

                             0.5 * fx[j] * pa_y[j] * pb_y[j] + pa_xy[j] * pb_xy[j]) * s_0_0[j];

                t_xy_xz[j] = (0.5 * fx[j] * pa_y[j] * pb_z[j] + pa_xy[j] * pb_xz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (9,18)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xz, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_xy_yy, t_xy_yz, t_xy_zz, t_xz_xx, t_xz_xy, t_xz_xz, t_xz_yy, \
                                     t_xz_yz, t_xz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xy_yy[j] = (0.5 * pa_xy[j] * fx[j] + pa_x[j] * fx[j] * pb_y[j] + pa_xy[j] * pb_yy[j]) * s_0_0[j];

                t_xy_yz[j] = (0.5 * pa_x[j] * fx[j] * pb_z[j] + pa_xy[j] * pb_yz[j]) * s_0_0[j];

                t_xy_zz[j] = (0.5 * pa_xy[j] * fx[j] + pa_xy[j] * pb_zz[j]) * s_0_0[j];

                t_xz_xx[j] = (0.5 * pa_xz[j] * fx[j] + fx[j] * pa_z[j] * pb_x[j] + pa_xz[j] * pb_xx[j]) * s_0_0[j];

                t_xz_xy[j] = (0.5 * fx[j] * pa_z[j] * pb_y[j] + pa_xz[j] * pb_xy[j]) * s_0_0[j];

                t_xz_xz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_x[j] +

                             0.5 * fx[j] * pa_z[j] * pb_z[j] + pa_xz[j] * pb_xz[j]) * s_0_0[j];

                t_xz_yy[j] = (0.5 * pa_xz[j] * fx[j] + pa_xz[j] * pb_yy[j]) * s_0_0[j];

                t_xz_yz[j] = (0.5 * pa_x[j] * fx[j] * pb_y[j] + pa_xz[j] * pb_yz[j]) * s_0_0[j];

                t_xz_zz[j] = (0.5 * pa_xz[j] * fx[j] + pa_x[j] * fx[j] * pb_z[j] + pa_xz[j] * pb_zz[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (18,27)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yz, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_yy_xx, t_yy_xy, t_yy_xz, t_yy_yy, t_yy_yz, t_yy_zz, t_yz_xx, \
                                     t_yz_xy, t_yz_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yy_xx[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] + 0.5 * fx[j] * pb_xx[j] +

                             pa_yy[j] * pb_xx[j]) * s_0_0[j];

                t_yy_xy[j] = (pa_y[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pb_xy[j] + pa_yy[j] * pb_xy[j]) * s_0_0[j];

                t_yy_xz[j] = (0.5 * fx[j] * pb_xz[j] + pa_yy[j] * pb_xz[j]) * s_0_0[j];

                t_yy_yy[j] = (0.75 * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] + 2.0 * pa_y[j] * fx[j] * pb_y[j] +

                             0.5 * fx[j] * pb_yy[j] + pa_yy[j] * pb_yy[j]) * s_0_0[j];

                t_yy_yz[j] = (pa_y[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pb_yz[j] + pa_yy[j] * pb_yz[j]) * s_0_0[j];

                t_yy_zz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] + 0.5 * fx[j] * pb_zz[j] +

                             pa_yy[j] * pb_zz[j]) * s_0_0[j];

                t_yz_xx[j] = (0.5 * pa_yz[j] * fx[j] + pa_yz[j] * pb_xx[j]) * s_0_0[j];

                t_yz_xy[j] = (0.5 * fx[j] * pa_z[j] * pb_x[j] + pa_yz[j] * pb_xy[j]) * s_0_0[j];

                t_yz_xz[j] = (0.5 * pa_y[j] * fx[j] * pb_x[j] + pa_yz[j] * pb_xz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (27,36)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_yz_yy, t_yz_yz, t_yz_zz, t_zz_xx, t_zz_xy, t_zz_xz, t_zz_yy, \
                                     t_zz_yz, t_zz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_yy[j] = (0.5 * pa_yz[j] * fx[j] + fx[j] * pa_z[j] * pb_y[j] + pa_yz[j] * pb_yy[j]) * s_0_0[j];

                t_yz_yz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_y[j] +

                             0.5 * fx[j] * pa_z[j] * pb_z[j] + pa_yz[j] * pb_yz[j]) * s_0_0[j];

                t_yz_zz[j] = (0.5 * pa_yz[j] * fx[j] + pa_y[j] * fx[j] * pb_z[j] + pa_yz[j] * pb_zz[j]) * s_0_0[j];

                t_zz_xx[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_zz[j] * fx[j] + 0.5 * fx[j] * pb_xx[j] +

                             pa_zz[j] * pb_xx[j]) * s_0_0[j];

                t_zz_xy[j] = (0.5 * fx[j] * pb_xy[j] + pa_zz[j] * pb_xy[j]) * s_0_0[j];

                t_zz_xz[j] = (pa_z[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pb_xz[j] + pa_zz[j] * pb_xz[j]) * s_0_0[j];

                t_zz_yy[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_zz[j] * fx[j] + 0.5 * fx[j] * pb_yy[j] +

                             pa_zz[j] * pb_yy[j]) * s_0_0[j];

                t_zz_yz[j] = (pa_z[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pb_yz[j] + pa_zz[j] * pb_yz[j]) * s_0_0[j];

                t_zz_zz[j] = (0.75 * fx[j] * fx[j] + 0.5 * pa_zz[j] * fx[j] + 2.0 * pa_z[j] * fx[j] * pb_z[j] +

                             0.5 * fx[j] * pb_zz[j] + pa_zz[j] * pb_zz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForDF(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(9 * idx);

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            auto pa_xz = paDistances.data(9 * idx + 5);

            auto pa_yy = paDistances.data(9 * idx + 6);

            auto pa_yz = paDistances.data(9 * idx + 7);

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_xxx = primBuffer.data(60 * idx);

            auto t_xx_xxy = primBuffer.data(60 * idx + 1);

            auto t_xx_xxz = primBuffer.data(60 * idx + 2);

            auto t_xx_xyy = primBuffer.data(60 * idx + 3);

            auto t_xx_xyz = primBuffer.data(60 * idx + 4);

            auto t_xx_xzz = primBuffer.data(60 * idx + 5);

            auto t_xx_yyy = primBuffer.data(60 * idx + 6);

            auto t_xx_yyz = primBuffer.data(60 * idx + 7);

            auto t_xx_yzz = primBuffer.data(60 * idx + 8);

            auto t_xx_zzz = primBuffer.data(60 * idx + 9);

            auto t_xy_xxx = primBuffer.data(60 * idx + 10);

            auto t_xy_xxy = primBuffer.data(60 * idx + 11);

            auto t_xy_xxz = primBuffer.data(60 * idx + 12);

            auto t_xy_xyy = primBuffer.data(60 * idx + 13);

            auto t_xy_xyz = primBuffer.data(60 * idx + 14);

            auto t_xy_xzz = primBuffer.data(60 * idx + 15);

            auto t_xy_yyy = primBuffer.data(60 * idx + 16);

            auto t_xy_yyz = primBuffer.data(60 * idx + 17);

            auto t_xy_yzz = primBuffer.data(60 * idx + 18);

            auto t_xy_zzz = primBuffer.data(60 * idx + 19);

            auto t_xz_xxx = primBuffer.data(60 * idx + 20);

            auto t_xz_xxy = primBuffer.data(60 * idx + 21);

            auto t_xz_xxz = primBuffer.data(60 * idx + 22);

            auto t_xz_xyy = primBuffer.data(60 * idx + 23);

            auto t_xz_xyz = primBuffer.data(60 * idx + 24);

            auto t_xz_xzz = primBuffer.data(60 * idx + 25);

            auto t_xz_yyy = primBuffer.data(60 * idx + 26);

            auto t_xz_yyz = primBuffer.data(60 * idx + 27);

            auto t_xz_yzz = primBuffer.data(60 * idx + 28);

            auto t_xz_zzz = primBuffer.data(60 * idx + 29);

            auto t_yy_xxx = primBuffer.data(60 * idx + 30);

            auto t_yy_xxy = primBuffer.data(60 * idx + 31);

            auto t_yy_xxz = primBuffer.data(60 * idx + 32);

            auto t_yy_xyy = primBuffer.data(60 * idx + 33);

            auto t_yy_xyz = primBuffer.data(60 * idx + 34);

            auto t_yy_xzz = primBuffer.data(60 * idx + 35);

            auto t_yy_yyy = primBuffer.data(60 * idx + 36);

            auto t_yy_yyz = primBuffer.data(60 * idx + 37);

            auto t_yy_yzz = primBuffer.data(60 * idx + 38);

            auto t_yy_zzz = primBuffer.data(60 * idx + 39);

            auto t_yz_xxx = primBuffer.data(60 * idx + 40);

            auto t_yz_xxy = primBuffer.data(60 * idx + 41);

            auto t_yz_xxz = primBuffer.data(60 * idx + 42);

            auto t_yz_xyy = primBuffer.data(60 * idx + 43);

            auto t_yz_xyz = primBuffer.data(60 * idx + 44);

            auto t_yz_xzz = primBuffer.data(60 * idx + 45);

            auto t_yz_yyy = primBuffer.data(60 * idx + 46);

            auto t_yz_yyz = primBuffer.data(60 * idx + 47);

            auto t_yz_yzz = primBuffer.data(60 * idx + 48);

            auto t_yz_zzz = primBuffer.data(60 * idx + 49);

            auto t_zz_xxx = primBuffer.data(60 * idx + 50);

            auto t_zz_xxy = primBuffer.data(60 * idx + 51);

            auto t_zz_xxz = primBuffer.data(60 * idx + 52);

            auto t_zz_xyy = primBuffer.data(60 * idx + 53);

            auto t_zz_xyz = primBuffer.data(60 * idx + 54);

            auto t_zz_xzz = primBuffer.data(60 * idx + 55);

            auto t_zz_yyy = primBuffer.data(60 * idx + 56);

            auto t_zz_yyz = primBuffer.data(60 * idx + 57);

            auto t_zz_yzz = primBuffer.data(60 * idx + 58);

            auto t_zz_zzz = primBuffer.data(60 * idx + 59);

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xx_xxx, t_xx_xxy, t_xx_xxz, t_xx_xyy, t_xx_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_xxx[j] = (1.5 * pa_x[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pb_x[j] +

                              1.5 * pa_xx[j] * pb_x[j] * fx[j] + 3.0 * pa_x[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pb_xxx[j] + pa_xx[j] * pb_xxx[j]) * s_0_0[j];

                t_xx_xxy[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xx[j] * fx[j] * pb_y[j] +

                              2.0 * pa_x[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pb_xxy[j] + pa_xx[j] * pb_xxy[j]) * s_0_0[j];

                t_xx_xxz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xx[j] * fx[j] * pb_z[j] +

                              2.0 * pa_x[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pb_xxz[j] + pa_xx[j] * pb_xxz[j]) * s_0_0[j];

                t_xx_xyy[j] = (0.5 * pa_x[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xx[j] * pb_x[j] * fx[j] + pa_x[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pb_xyy[j] + pa_xx[j] * pb_xyy[j]) * s_0_0[j];

                t_xx_xyz[j] = (pa_x[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pb_xyz[j] + pa_xx[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     s_0_0, t_xx_xzz, t_xx_yyy, t_xx_yyz, t_xx_yzz, t_xx_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_xzz[j] = (0.5 * pa_x[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xx[j] * pb_x[j] * fx[j] + pa_x[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pb_xzz[j] + pa_xx[j] * pb_xzz[j]) * s_0_0[j];

                t_xx_yyy[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xx[j] * pb_y[j] * fx[j] +

                              0.5 * fx[j] * pb_yyy[j] + pa_xx[j] * pb_yyy[j]) * s_0_0[j];

                t_xx_yyz[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xx[j] * fx[j] * pb_z[j] +

                              0.5 * fx[j] * pb_yyz[j] + pa_xx[j] * pb_yyz[j]) * s_0_0[j];

                t_xx_yzz[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xx[j] * pb_y[j] * fx[j] +

                              0.5 * fx[j] * pb_yzz[j] + pa_xx[j] * pb_yzz[j]) * s_0_0[j];

                t_xx_zzz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xx[j] * pb_z[j] * fx[j] +

                              0.5 * fx[j] * pb_zzz[j] + pa_xx[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xy_xxx, t_xy_xxy, t_xy_xxz, t_xy_xyy, \
                                     t_xy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xy_xxx[j] = (0.75 * fx[j] * fx[j] * pa_y[j] + 1.5 * pa_xy[j] * pb_x[j] * fx[j] +

                              1.5 * fx[j] * pa_y[j] * pb_xx[j] + pa_xy[j] * pb_xxx[j]) * s_0_0[j];

                t_xy_xxy[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xy[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * pb_xx[j] + fx[j] * pa_y[j] * pb_xy[j] + pa_xy[j] * pb_xxy[j]) * s_0_0[j];

                t_xy_xxz[j] = (0.5 * pa_xy[j] * fx[j] * pb_z[j] + fx[j] * pa_y[j] * pb_xz[j] + pa_xy[j] * pb_xxz[j]) * s_0_0[j];

                t_xy_xyy[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_xy[j] * pb_x[j] * fx[j] + pa_x[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_y[j] * pb_yy[j] + pa_xy[j] * pb_xyy[j]) * s_0_0[j];

                t_xy_xyz[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * pb_xz[j] +

                              0.5 * fx[j] * pa_y[j] * pb_yz[j] + pa_xy[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_y, pb_x, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, s_0_0, t_xy_xzz, t_xy_yyy, t_xy_yyz, t_xy_yzz, t_xy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xy_xzz[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xy[j] * pb_x[j] * fx[j] +

                              0.5 * fx[j] * pa_y[j] * pb_zz[j] + pa_xy[j] * pb_xzz[j]) * s_0_0[j];

                t_xy_yyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 1.5 * pa_xy[j] * pb_y[j] * fx[j] +

                              1.5 * pa_x[j] * fx[j] * pb_yy[j] + pa_xy[j] * pb_yyy[j]) * s_0_0[j];

                t_xy_yyz[j] = (0.5 * pa_xy[j] * fx[j] * pb_z[j] + pa_x[j] * fx[j] * pb_yz[j] + pa_xy[j] * pb_yyz[j]) * s_0_0[j];

                t_xy_yzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xy[j] * pb_y[j] * fx[j] +

                              0.5 * pa_x[j] * fx[j] * pb_zz[j] + pa_xy[j] * pb_yzz[j]) * s_0_0[j];

                t_xy_zzz[j] = (1.5 * pa_xy[j] * pb_z[j] * fx[j] + pa_xy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xz_xxx, t_xz_xxy, t_xz_xxz, t_xz_xyy, \
                                     t_xz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_xxx[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 1.5 * pa_xz[j] * pb_x[j] * fx[j] +

                              1.5 * fx[j] * pa_z[j] * pb_xx[j] + pa_xz[j] * pb_xxx[j]) * s_0_0[j];

                t_xz_xxy[j] = (0.5 * pa_xz[j] * fx[j] * pb_y[j] + fx[j] * pa_z[j] * pb_xy[j] + pa_xz[j] * pb_xxy[j]) * s_0_0[j];

                t_xz_xxz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xz[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * pb_xx[j] + fx[j] * pa_z[j] * pb_xz[j] + pa_xz[j] * pb_xxz[j]) * s_0_0[j];

                t_xz_xyy[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_xz[j] * pb_x[j] * fx[j] +

                              0.5 * fx[j] * pa_z[j] * pb_yy[j] + pa_xz[j] * pb_xyy[j]) * s_0_0[j];

                t_xz_xyz[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * pb_xy[j] +

                              0.5 * fx[j] * pa_z[j] * pb_yz[j] + pa_xz[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_z, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xz_xzz, t_xz_yyy, t_xz_yyz, t_xz_yzz, t_xz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_xzz[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_xz[j] * pb_x[j] * fx[j] + pa_x[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_z[j] * pb_zz[j] + pa_xz[j] * pb_xzz[j]) * s_0_0[j];

                t_xz_yyy[j] = (1.5 * pa_xz[j] * pb_y[j] * fx[j] + pa_xz[j] * pb_yyy[j]) * s_0_0[j];

                t_xz_yyz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xz[j] * fx[j] * pb_z[j] +

                              0.5 * pa_x[j] * fx[j] * pb_yy[j] + pa_xz[j] * pb_yyz[j]) * s_0_0[j];

                t_xz_yzz[j] = (0.5 * pa_xz[j] * pb_y[j] * fx[j] + pa_x[j] * fx[j] * pb_yz[j] + pa_xz[j] * pb_yzz[j]) * s_0_0[j];

                t_xz_zzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 1.5 * pa_xz[j] * pb_z[j] * fx[j] +

                              1.5 * pa_x[j] * fx[j] * pb_zz[j] + pa_xz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_y, pb_z, s_0_0, t_yy_xxx, t_yy_xxy, t_yy_xxz, t_yy_xyy, t_yy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yy_xxx[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yy[j] * pb_x[j] * fx[j] +

                              0.5 * fx[j] * pb_xxx[j] + pa_yy[j] * pb_xxx[j]) * s_0_0[j];

                t_yy_xxy[j] = (0.5 * pa_y[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_yy[j] * fx[j] * pb_y[j] + pa_y[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pb_xxy[j] + pa_yy[j] * pb_xxy[j]) * s_0_0[j];

                t_yy_xxz[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_yy[j] * fx[j] * pb_z[j] +

                              0.5 * fx[j] * pb_xxz[j] + pa_yy[j] * pb_xxz[j]) * s_0_0[j];

                t_yy_xyy[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_yy[j] * pb_x[j] * fx[j] +

                              2.0 * pa_y[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pb_xyy[j] + pa_yy[j] * pb_xyy[j]) * s_0_0[j];

                t_yy_xyz[j] = (pa_y[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pb_xyz[j] + pa_yy[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pb_x, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, \
                                     pb_zz, pb_zzz, s_0_0, t_yy_xzz, t_yy_yyy, t_yy_yyz, t_yy_yzz, t_yy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yy_xzz[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_yy[j] * pb_x[j] * fx[j] +

                              0.5 * fx[j] * pb_xzz[j] + pa_yy[j] * pb_xzz[j]) * s_0_0[j];

                t_yy_yyy[j] = (1.5 * pa_y[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pb_y[j] +

                              1.5 * pa_yy[j] * pb_y[j] * fx[j] + 3.0 * pa_y[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pb_yyy[j] + pa_yy[j] * pb_yyy[j]) * s_0_0[j];

                t_yy_yyz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_yy[j] * fx[j] * pb_z[j] +

                              2.0 * pa_y[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pb_yyz[j] + pa_yy[j] * pb_yyz[j]) * s_0_0[j];

                t_yy_yzz[j] = (0.5 * pa_y[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_yy[j] * pb_y[j] * fx[j] + pa_y[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pb_yzz[j] + pa_yy[j] * pb_yzz[j]) * s_0_0[j];

                t_yy_zzz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_yy[j] * pb_z[j] * fx[j] +

                              0.5 * fx[j] * pb_zzz[j] + pa_yy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_z, s_0_0, t_yz_xxx, t_yz_xxy, t_yz_xxz, t_yz_xyy, t_yz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_xxx[j] = (1.5 * pa_yz[j] * pb_x[j] * fx[j] + pa_yz[j] * pb_xxx[j]) * s_0_0[j];

                t_yz_xxy[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_yz[j] * fx[j] * pb_y[j] +

                              0.5 * fx[j] * pa_z[j] * pb_xx[j] + pa_yz[j] * pb_xxy[j]) * s_0_0[j];

                t_yz_xxz[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * pa_yz[j] * fx[j] * pb_z[j] +

                              0.5 * pa_y[j] * fx[j] * pb_xx[j] + pa_yz[j] * pb_xxz[j]) * s_0_0[j];

                t_yz_xyy[j] = (0.5 * pa_yz[j] * pb_x[j] * fx[j] + fx[j] * pa_z[j] * pb_xy[j] + pa_yz[j] * pb_xyy[j]) * s_0_0[j];

                t_yz_xyz[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_y[j] * fx[j] * pb_xy[j] +

                              0.5 * fx[j] * pa_z[j] * pb_xz[j] + pa_yz[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_yz_xzz, t_yz_yyy, t_yz_yyz, t_yz_yzz, t_yz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_xzz[j] = (0.5 * pa_yz[j] * pb_x[j] * fx[j] + pa_y[j] * fx[j] * pb_xz[j] + pa_yz[j] * pb_xzz[j]) * s_0_0[j];

                t_yz_yyy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 1.5 * pa_yz[j] * pb_y[j] * fx[j] +

                              1.5 * fx[j] * pa_z[j] * pb_yy[j] + pa_yz[j] * pb_yyy[j]) * s_0_0[j];

                t_yz_yyz[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_yz[j] * fx[j] * pb_z[j] + 0.5 * pa_y[j] * fx[j] * pb_yy[j] + fx[j] * pa_z[j] * pb_yz[j] + pa_yz[j] * pb_yyz[j]) * s_0_0[j];

                t_yz_yzz[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_yz[j] * pb_y[j] * fx[j] + pa_y[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_z[j] * pb_zz[j] + pa_yz[j] * pb_yzz[j]) * s_0_0[j];

                t_yz_zzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 1.5 * pa_yz[j] * pb_z[j] * fx[j] +

                              1.5 * pa_y[j] * fx[j] * pb_zz[j] + pa_yz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_y, \
                                     pb_z, s_0_0, t_zz_xxx, t_zz_xxy, t_zz_xxz, t_zz_xyy, t_zz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zz_xxx[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_zz[j] * pb_x[j] * fx[j] +

                              0.5 * fx[j] * pb_xxx[j] + pa_zz[j] * pb_xxx[j]) * s_0_0[j];

                t_zz_xxy[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_zz[j] * fx[j] * pb_y[j] +

                              0.5 * fx[j] * pb_xxy[j] + pa_zz[j] * pb_xxy[j]) * s_0_0[j];

                t_zz_xxz[j] = (0.5 * pa_z[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_zz[j] * fx[j] * pb_z[j] + pa_z[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pb_xxz[j] + pa_zz[j] * pb_xxz[j]) * s_0_0[j];

                t_zz_xyy[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_zz[j] * pb_x[j] * fx[j] +

                              0.5 * fx[j] * pb_xyy[j] + pa_zz[j] * pb_xyy[j]) * s_0_0[j];

                t_zz_xyz[j] = (pa_z[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pb_xyz[j] + pa_zz[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, s_0_0, t_zz_xzz, t_zz_yyy, t_zz_yyz, t_zz_yzz, t_zz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zz_xzz[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_zz[j] * pb_x[j] * fx[j] +

                              2.0 * pa_z[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pb_xzz[j] + pa_zz[j] * pb_xzz[j]) * s_0_0[j];

                t_zz_yyy[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_zz[j] * pb_y[j] * fx[j] +

                              0.5 * fx[j] * pb_yyy[j] + pa_zz[j] * pb_yyy[j]) * s_0_0[j];

                t_zz_yyz[j] = (0.5 * pa_z[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_zz[j] * fx[j] * pb_z[j] + pa_z[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pb_yyz[j] + pa_zz[j] * pb_yyz[j]) * s_0_0[j];

                t_zz_yzz[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_zz[j] * pb_y[j] * fx[j] +

                              2.0 * pa_z[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pb_yzz[j] + pa_zz[j] * pb_yzz[j]) * s_0_0[j];

                t_zz_zzz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pb_z[j] +

                              1.5 * pa_zz[j] * pb_z[j] * fx[j] + 3.0 * pa_z[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pb_zzz[j] + pa_zz[j] * pb_zzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForDG(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(9 * idx);

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            auto pa_xz = paDistances.data(9 * idx + 5);

            auto pa_yy = paDistances.data(9 * idx + 6);

            auto pa_yz = paDistances.data(9 * idx + 7);

            auto pa_zz = paDistances.data(9 * idx + 8);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xx_xxxx = primBuffer.data(90 * idx);

            auto t_xx_xxxy = primBuffer.data(90 * idx + 1);

            auto t_xx_xxxz = primBuffer.data(90 * idx + 2);

            auto t_xx_xxyy = primBuffer.data(90 * idx + 3);

            auto t_xx_xxyz = primBuffer.data(90 * idx + 4);

            auto t_xx_xxzz = primBuffer.data(90 * idx + 5);

            auto t_xx_xyyy = primBuffer.data(90 * idx + 6);

            auto t_xx_xyyz = primBuffer.data(90 * idx + 7);

            auto t_xx_xyzz = primBuffer.data(90 * idx + 8);

            auto t_xx_xzzz = primBuffer.data(90 * idx + 9);

            auto t_xx_yyyy = primBuffer.data(90 * idx + 10);

            auto t_xx_yyyz = primBuffer.data(90 * idx + 11);

            auto t_xx_yyzz = primBuffer.data(90 * idx + 12);

            auto t_xx_yzzz = primBuffer.data(90 * idx + 13);

            auto t_xx_zzzz = primBuffer.data(90 * idx + 14);

            auto t_xy_xxxx = primBuffer.data(90 * idx + 15);

            auto t_xy_xxxy = primBuffer.data(90 * idx + 16);

            auto t_xy_xxxz = primBuffer.data(90 * idx + 17);

            auto t_xy_xxyy = primBuffer.data(90 * idx + 18);

            auto t_xy_xxyz = primBuffer.data(90 * idx + 19);

            auto t_xy_xxzz = primBuffer.data(90 * idx + 20);

            auto t_xy_xyyy = primBuffer.data(90 * idx + 21);

            auto t_xy_xyyz = primBuffer.data(90 * idx + 22);

            auto t_xy_xyzz = primBuffer.data(90 * idx + 23);

            auto t_xy_xzzz = primBuffer.data(90 * idx + 24);

            auto t_xy_yyyy = primBuffer.data(90 * idx + 25);

            auto t_xy_yyyz = primBuffer.data(90 * idx + 26);

            auto t_xy_yyzz = primBuffer.data(90 * idx + 27);

            auto t_xy_yzzz = primBuffer.data(90 * idx + 28);

            auto t_xy_zzzz = primBuffer.data(90 * idx + 29);

            auto t_xz_xxxx = primBuffer.data(90 * idx + 30);

            auto t_xz_xxxy = primBuffer.data(90 * idx + 31);

            auto t_xz_xxxz = primBuffer.data(90 * idx + 32);

            auto t_xz_xxyy = primBuffer.data(90 * idx + 33);

            auto t_xz_xxyz = primBuffer.data(90 * idx + 34);

            auto t_xz_xxzz = primBuffer.data(90 * idx + 35);

            auto t_xz_xyyy = primBuffer.data(90 * idx + 36);

            auto t_xz_xyyz = primBuffer.data(90 * idx + 37);

            auto t_xz_xyzz = primBuffer.data(90 * idx + 38);

            auto t_xz_xzzz = primBuffer.data(90 * idx + 39);

            auto t_xz_yyyy = primBuffer.data(90 * idx + 40);

            auto t_xz_yyyz = primBuffer.data(90 * idx + 41);

            auto t_xz_yyzz = primBuffer.data(90 * idx + 42);

            auto t_xz_yzzz = primBuffer.data(90 * idx + 43);

            auto t_xz_zzzz = primBuffer.data(90 * idx + 44);

            auto t_yy_xxxx = primBuffer.data(90 * idx + 45);

            auto t_yy_xxxy = primBuffer.data(90 * idx + 46);

            auto t_yy_xxxz = primBuffer.data(90 * idx + 47);

            auto t_yy_xxyy = primBuffer.data(90 * idx + 48);

            auto t_yy_xxyz = primBuffer.data(90 * idx + 49);

            auto t_yy_xxzz = primBuffer.data(90 * idx + 50);

            auto t_yy_xyyy = primBuffer.data(90 * idx + 51);

            auto t_yy_xyyz = primBuffer.data(90 * idx + 52);

            auto t_yy_xyzz = primBuffer.data(90 * idx + 53);

            auto t_yy_xzzz = primBuffer.data(90 * idx + 54);

            auto t_yy_yyyy = primBuffer.data(90 * idx + 55);

            auto t_yy_yyyz = primBuffer.data(90 * idx + 56);

            auto t_yy_yyzz = primBuffer.data(90 * idx + 57);

            auto t_yy_yzzz = primBuffer.data(90 * idx + 58);

            auto t_yy_zzzz = primBuffer.data(90 * idx + 59);

            auto t_yz_xxxx = primBuffer.data(90 * idx + 60);

            auto t_yz_xxxy = primBuffer.data(90 * idx + 61);

            auto t_yz_xxxz = primBuffer.data(90 * idx + 62);

            auto t_yz_xxyy = primBuffer.data(90 * idx + 63);

            auto t_yz_xxyz = primBuffer.data(90 * idx + 64);

            auto t_yz_xxzz = primBuffer.data(90 * idx + 65);

            auto t_yz_xyyy = primBuffer.data(90 * idx + 66);

            auto t_yz_xyyz = primBuffer.data(90 * idx + 67);

            auto t_yz_xyzz = primBuffer.data(90 * idx + 68);

            auto t_yz_xzzz = primBuffer.data(90 * idx + 69);

            auto t_yz_yyyy = primBuffer.data(90 * idx + 70);

            auto t_yz_yyyz = primBuffer.data(90 * idx + 71);

            auto t_yz_yyzz = primBuffer.data(90 * idx + 72);

            auto t_yz_yzzz = primBuffer.data(90 * idx + 73);

            auto t_yz_zzzz = primBuffer.data(90 * idx + 74);

            auto t_zz_xxxx = primBuffer.data(90 * idx + 75);

            auto t_zz_xxxy = primBuffer.data(90 * idx + 76);

            auto t_zz_xxxz = primBuffer.data(90 * idx + 77);

            auto t_zz_xxyy = primBuffer.data(90 * idx + 78);

            auto t_zz_xxyz = primBuffer.data(90 * idx + 79);

            auto t_zz_xxzz = primBuffer.data(90 * idx + 80);

            auto t_zz_xyyy = primBuffer.data(90 * idx + 81);

            auto t_zz_xyyz = primBuffer.data(90 * idx + 82);

            auto t_zz_xyzz = primBuffer.data(90 * idx + 83);

            auto t_zz_xzzz = primBuffer.data(90 * idx + 84);

            auto t_zz_yyyy = primBuffer.data(90 * idx + 85);

            auto t_zz_yyyz = primBuffer.data(90 * idx + 86);

            auto t_zz_yyzz = primBuffer.data(90 * idx + 87);

            auto t_zz_yzzz = primBuffer.data(90 * idx + 88);

            auto t_zz_zzzz = primBuffer.data(90 * idx + 89);

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xx_xxxx, \
                                     t_xx_xxxy, t_xx_xxxz, t_xx_xxyy, t_xx_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_xxxx[j] = (1.875 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               6.0 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_xx[j] * pb_xx[j] * fx[j] +

                               4.0 * pa_x[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pb_xxxx[j] + pa_xx[j] * pb_xxxx[j]) * s_0_0[j];

                t_xx_xxxy[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * fx[j] * fx[j] * pb_xy[j] +

                               1.5 * pa_xx[j] * pb_xy[j] * fx[j] + 3.0 * pa_x[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pb_xxxy[j] + pa_xx[j] * pb_xxxy[j]) * s_0_0[j];

                t_xx_xxxz[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * fx[j] * fx[j] * pb_xz[j] +

                               1.5 * pa_xx[j] * pb_xz[j] * fx[j] + 3.0 * pa_x[j] * fx[j] * pb_xxz[j] + 0.5 * fx[j] * pb_xxxz[j] + pa_xx[j] * pb_xxxz[j]) * s_0_0[j];

                t_xx_xxyy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] +

                               0.5 * pa_xx[j] * pb_xx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * pb_yy[j] + 2.0 * pa_x[j] * fx[j] * pb_xyy[j] +

                               0.5 * fx[j] * pb_xxyy[j] + pa_xx[j] * pb_xxyy[j]) * s_0_0[j];

                t_xx_xxyz[j] = (0.75 * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yz[j] +

                               2.0 * pa_x[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pb_xxyz[j] + pa_xx[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, pb_xyzz, pb_xz, \
                                     pb_xzz, pb_xzzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xx_xxzz, \
                                     t_xx_xyyy, t_xx_xyyz, t_xx_xyzz, t_xx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_xxzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] +

                               0.5 * pa_xx[j] * pb_xx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * pb_zz[j] + 2.0 * pa_x[j] * fx[j] * pb_xzz[j] +

                               0.5 * fx[j] * pb_xxzz[j] + pa_xx[j] * pb_xxzz[j]) * s_0_0[j];

                t_xx_xyyy[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_xy[j] +

                               1.5 * pa_xx[j] * pb_xy[j] * fx[j] + pa_x[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pb_xyyy[j] + pa_xx[j] * pb_xyyy[j]) * s_0_0[j];

                t_xx_xyyz[j] = (0.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_xz[j] +

                               0.5 * pa_xx[j] * pb_xz[j] * fx[j] + pa_x[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pb_xyyz[j] + pa_xx[j] * pb_xyyz[j]) * s_0_0[j];

                t_xx_xyzz[j] = (0.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_xx[j] * pb_xy[j] * fx[j] + pa_x[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pb_xyzz[j] + pa_xx[j] * pb_xyzz[j]) * s_0_0[j];

                t_xx_xzzz[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_xz[j] +

                               1.5 * pa_xx[j] * pb_xz[j] * fx[j] + pa_x[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pb_xzzz[j] + pa_xx[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fx, pa_xx, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_zz, pb_zzzz, \
                                     s_0_0, t_xx_yyyy, t_xx_yyyz, t_xx_yyzz, t_xx_yzzz, t_xx_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_yyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               1.5 * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xx[j] * pb_yy[j] * fx[j] + 0.5 * fx[j] * pb_yyyy[j] + pa_xx[j] * pb_yyyy[j]) * s_0_0[j];

                t_xx_yyyz[j] = (0.75 * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xx[j] * pb_yz[j] * fx[j] +

                               0.5 * fx[j] * pb_yyyz[j] + pa_xx[j] * pb_yyyz[j]) * s_0_0[j];

                t_xx_yyzz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xx[j] * pb_yy[j] * fx[j] +

                               0.5 * pa_xx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pb_yyzz[j] + pa_xx[j] * pb_yyzz[j]) * s_0_0[j];

                t_xx_yzzz[j] = (0.75 * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xx[j] * pb_yz[j] * fx[j] +

                               0.5 * fx[j] * pb_yzzz[j] + pa_xx[j] * pb_yzzz[j]) * s_0_0[j];

                t_xx_zzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               1.5 * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xx[j] * pb_zz[j] * fx[j] + 0.5 * fx[j] * pb_zzzz[j] + pa_xx[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, \
                                     t_xy_xxxx, t_xy_xxxy, t_xy_xxxz, t_xy_xxyy, t_xy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xy_xxxx[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               3.0 * pa_xy[j] * pb_xx[j] * fx[j] + 2.0 * fx[j] * pa_y[j] * pb_xxx[j] + pa_xy[j] * pb_xxxx[j]) * s_0_0[j];

                t_xy_xxxy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] +

                               0.75 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xy[j] * pb_xy[j] * fx[j] +

                               0.5 * pa_x[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_y[j] * pb_xxy[j] + pa_xy[j] * pb_xxxy[j]) * s_0_0[j];

                t_xy_xxxz[j] = (0.75 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 1.5 * pa_xy[j] * pb_xz[j] * fx[j] +

                               1.5 * fx[j] * pa_y[j] * pb_xxz[j] + pa_xy[j] * pb_xxxz[j]) * s_0_0[j];

                t_xy_xxyy[j] = (0.25 * pa_xy[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.5 * fx[j] * fx[j] * pa_y[j] * pb_x[j] + fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xy[j] * pb_xx[j] * fx[j] +

                               0.5 * pa_xy[j] * fx[j] * pb_yy[j] + pa_x[j] * fx[j] * pb_xxy[j] + fx[j] * pa_y[j] * pb_xyy[j] + pa_xy[j] * pb_xxyy[j]) * s_0_0[j];

                t_xy_xxyz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pb_xz[j] +

                               0.5 * pa_xy[j] * fx[j] * pb_yz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxz[j] + fx[j] * pa_y[j] * pb_xyz[j] +

                               pa_xy[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, \
                                     pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, s_0_0, t_xy_xxzz, t_xy_xyyy, t_xy_xyyz, t_xy_xyzz, \
                                     t_xy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xy_xxzz[j] = (0.25 * pa_xy[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               0.5 * pa_xy[j] * pb_xx[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * pb_zz[j] + fx[j] * pa_y[j] * pb_xzz[j] +

                               pa_xy[j] * pb_xxzz[j]) * s_0_0[j];

                t_xy_xyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] +

                               0.75 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xy[j] * pb_xy[j] * fx[j] +

                               1.5 * pa_x[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_y[j] * pb_yyy[j] + pa_xy[j] * pb_xyyy[j]) * s_0_0[j];

                t_xy_xyyz[j] = (0.25 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pb_yz[j] +

                               0.5 * pa_xy[j] * pb_xz[j] * fx[j] + pa_x[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_y[j] * pb_yyz[j] +

                               pa_xy[j] * pb_xyyz[j]) * s_0_0[j];

                t_xy_xyzz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_x[j] +

                               0.25 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xy[j] * pb_xy[j] * fx[j] +

                               0.5 * pa_x[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_y[j] * pb_yzz[j] + pa_xy[j] * pb_xyzz[j]) * s_0_0[j];

                t_xy_xzzz[j] = (0.75 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 1.5 * pa_xy[j] * pb_xz[j] * fx[j] +

                               0.5 * fx[j] * pa_y[j] * pb_zzz[j] + pa_xy[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, \
                                     pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xy_yyyy, t_xy_yyyz, \
                                     t_xy_yyzz, t_xy_yzzz, t_xy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xy_yyyy[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               3.0 * pa_xy[j] * pb_yy[j] * fx[j] + 2.0 * pa_x[j] * fx[j] * pb_yyy[j] + pa_xy[j] * pb_yyyy[j]) * s_0_0[j];

                t_xy_yyyz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xy[j] * pb_yz[j] * fx[j] +

                               1.5 * pa_x[j] * fx[j] * pb_yyz[j] + pa_xy[j] * pb_yyyz[j]) * s_0_0[j];

                t_xy_yyzz[j] = (0.25 * pa_xy[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.5 * pa_xy[j] * pb_yy[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * pb_zz[j] + pa_x[j] * fx[j] * pb_yzz[j] +

                               pa_xy[j] * pb_yyzz[j]) * s_0_0[j];

                t_xy_yzzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xy[j] * pb_yz[j] * fx[j] +

                               0.5 * pa_x[j] * fx[j] * pb_zzz[j] + pa_xy[j] * pb_yzzz[j]) * s_0_0[j];

                t_xy_zzzz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] + 3.0 * pa_xy[j] * pb_zz[j] * fx[j] +

                               pa_xy[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, \
                                     t_xz_xxxx, t_xz_xxxy, t_xz_xxxz, t_xz_xxyy, t_xz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_xxxx[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               3.0 * pa_xz[j] * pb_xx[j] * fx[j] + 2.0 * fx[j] * pa_z[j] * pb_xxx[j] + pa_xz[j] * pb_xxxx[j]) * s_0_0[j];

                t_xz_xxxy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 1.5 * pa_xz[j] * pb_xy[j] * fx[j] +

                               1.5 * fx[j] * pa_z[j] * pb_xxy[j] + pa_xz[j] * pb_xxxy[j]) * s_0_0[j];

                t_xz_xxxz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] +

                               0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xz[j] * pb_xz[j] * fx[j] +

                               0.5 * pa_x[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_z[j] * pb_xxz[j] + pa_xz[j] * pb_xxxz[j]) * s_0_0[j];

                t_xz_xxyy[j] = (0.25 * pa_xz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               0.5 * pa_xz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xz[j] * fx[j] * pb_yy[j] + fx[j] * pa_z[j] * pb_xyy[j] +

                               pa_xz[j] * pb_xxyy[j]) * s_0_0[j];

                t_xz_xxyz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_xz[j] * fx[j] * pb_yz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxy[j] + fx[j] * pa_z[j] * pb_xyz[j] +

                               pa_xz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xz_xxzz, t_xz_xyyy, t_xz_xyyz, t_xz_xyzz, \
                                     t_xz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_xxzz[j] = (0.25 * pa_xz[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xz[j] * pb_xx[j] * fx[j] +

                               0.5 * pa_xz[j] * fx[j] * pb_zz[j] + pa_x[j] * fx[j] * pb_xxz[j] + fx[j] * pa_z[j] * pb_xzz[j] + pa_xz[j] * pb_xxzz[j]) * s_0_0[j];

                t_xz_xyyy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 1.5 * pa_xz[j] * pb_xy[j] * fx[j] +

                               0.5 * fx[j] * pa_z[j] * pb_yyy[j] + pa_xz[j] * pb_xyyy[j]) * s_0_0[j];

                t_xz_xyyz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_x[j] +

                               0.25 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xz[j] * pb_xz[j] * fx[j] +

                               0.5 * pa_x[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_z[j] * pb_yyz[j] + pa_xz[j] * pb_xyyz[j]) * s_0_0[j];

                t_xz_xyzz[j] = (0.25 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pb_yz[j] +

                               0.5 * pa_xz[j] * pb_xy[j] * fx[j] + pa_x[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_z[j] * pb_yzz[j] +

                               pa_xz[j] * pb_xyzz[j]) * s_0_0[j];

                t_xz_xzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] +

                               0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xz[j] * pb_xz[j] * fx[j] +

                               1.5 * pa_x[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_z[j] * pb_zzz[j] + pa_xz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, \
                                     pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xz_yyyy, t_xz_yyyz, \
                                     t_xz_yyzz, t_xz_yzzz, t_xz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_yyyy[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] + 3.0 * pa_xz[j] * pb_yy[j] * fx[j] +

                               pa_xz[j] * pb_yyyy[j]) * s_0_0[j];

                t_xz_yyyz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xz[j] * pb_yz[j] * fx[j] +

                               0.5 * pa_x[j] * fx[j] * pb_yyy[j] + pa_xz[j] * pb_yyyz[j]) * s_0_0[j];

                t_xz_yyzz[j] = (0.25 * pa_xz[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_xz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xz[j] * fx[j] * pb_zz[j] + pa_x[j] * fx[j] * pb_yyz[j] +

                               pa_xz[j] * pb_yyzz[j]) * s_0_0[j];

                t_xz_yzzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xz[j] * pb_yz[j] * fx[j] +

                               1.5 * pa_x[j] * fx[j] * pb_yzz[j] + pa_xz[j] * pb_yzzz[j]) * s_0_0[j];

                t_xz_zzzz[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               3.0 * pa_xz[j] * pb_zz[j] * fx[j] + 2.0 * pa_x[j] * fx[j] * pb_zzz[j] + pa_xz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_yy_xxxx, t_yy_xxxy, \
                                     t_yy_xxxz, t_yy_xxyy, t_yy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yy_xxxx[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] +

                               1.5 * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_yy[j] * pb_xx[j] * fx[j] + 0.5 * fx[j] * pb_xxxx[j] + pa_yy[j] * pb_xxxx[j]) * s_0_0[j];

                t_yy_xxxy[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pb_xy[j] +

                               1.5 * pa_yy[j] * pb_xy[j] * fx[j] + pa_y[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pb_xxxy[j] + pa_yy[j] * pb_xxxy[j]) * s_0_0[j];

                t_yy_xxxz[j] = (0.75 * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yy[j] * pb_xz[j] * fx[j] +

                               0.5 * fx[j] * pb_xxxz[j] + pa_yy[j] * pb_xxxz[j]) * s_0_0[j];

                t_yy_xxyy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] +

                               pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] +

                               0.5 * pa_yy[j] * pb_xx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] * pb_yy[j] + 2.0 * pa_y[j] * fx[j] * pb_xxy[j] +

                               0.5 * fx[j] * pb_xxyy[j] + pa_yy[j] * pb_xxyy[j]) * s_0_0[j];

                t_yy_xxyz[j] = (0.5 * pa_y[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_yz[j] +

                               0.5 * pa_yy[j] * fx[j] * pb_yz[j] + pa_y[j] * fx[j] * pb_xxz[j] + 0.5 * fx[j] * pb_xxyz[j] + pa_yy[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_zz, s_0_0, t_yy_xxzz, t_yy_xyyy, t_yy_xyyz, \
                                     t_yy_xyzz, t_yy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yy_xxzz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_yy[j] * pb_xx[j] * fx[j] +

                               0.5 * pa_yy[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pb_xxzz[j] + pa_yy[j] * pb_xxzz[j]) * s_0_0[j];

                t_yy_xyyy[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pb_xy[j] +

                               1.5 * pa_yy[j] * pb_xy[j] * fx[j] + 3.0 * pa_y[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pb_xyyy[j] + pa_yy[j] * pb_xyyy[j]) * s_0_0[j];

                t_yy_xyyz[j] = (0.75 * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_yy[j] * pb_xz[j] * fx[j] +

                               2.0 * pa_y[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pb_xyyz[j] + pa_yy[j] * pb_xyyz[j]) * s_0_0[j];

                t_yy_xyzz[j] = (0.5 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_yy[j] * pb_xy[j] * fx[j] + pa_y[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pb_xyzz[j] + pa_yy[j] * pb_xyzz[j]) * s_0_0[j];

                t_yy_xzzz[j] = (0.75 * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yy[j] * pb_xz[j] * fx[j] +

                               0.5 * fx[j] * pb_xzzz[j] + pa_yy[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, \
                                     pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_yy_yyyy, t_yy_yyyz, \
                                     t_yy_yyzz, t_yy_yzzz, t_yy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yy_yyyy[j] = (1.875 * fx[j] * fx[j] * fx[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] +

                               6.0 * pa_y[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_yy[j] * pb_yy[j] * fx[j] +

                               4.0 * pa_y[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pb_yyyy[j] + pa_yy[j] * pb_yyyy[j]) * s_0_0[j];

                t_yy_yyyz[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * fx[j] * fx[j] * pb_yz[j] +

                               1.5 * pa_yy[j] * pb_yz[j] * fx[j] + 3.0 * pa_y[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pb_yyyz[j] + pa_yy[j] * pb_yyyz[j]) * s_0_0[j];

                t_yy_yyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] +

                               pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] +

                               0.5 * pa_yy[j] * pb_yy[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] * pb_zz[j] + 2.0 * pa_y[j] * fx[j] * pb_yzz[j] +

                               0.5 * fx[j] * pb_yyzz[j] + pa_yy[j] * pb_yyzz[j]) * s_0_0[j];

                t_yy_yzzz[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_yz[j] +

                               1.5 * pa_yy[j] * pb_yz[j] * fx[j] + pa_y[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pb_yzzz[j] + pa_yy[j] * pb_yzzz[j]) * s_0_0[j];

                t_yy_zzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] +

                               1.5 * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_yy[j] * pb_zz[j] * fx[j] + 0.5 * fx[j] * pb_zzzz[j] + pa_yy[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (12) = (60,65)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_yz_xxxx, \
                                     t_yz_xxxy, t_yz_xxxz, t_yz_xxyy, t_yz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_xxxx[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] + 3.0 * pa_yz[j] * pb_xx[j] * fx[j] +

                               pa_yz[j] * pb_xxxx[j]) * s_0_0[j];

                t_yz_xxxy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 1.5 * pa_yz[j] * pb_xy[j] * fx[j] +

                               0.5 * fx[j] * pa_z[j] * pb_xxx[j] + pa_yz[j] * pb_xxxy[j]) * s_0_0[j];

                t_yz_xxxz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yz[j] * pb_xz[j] * fx[j] +

                               0.5 * pa_y[j] * fx[j] * pb_xxx[j] + pa_yz[j] * pb_xxxz[j]) * s_0_0[j];

                t_yz_xxyy[j] = (0.25 * pa_yz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               0.5 * pa_yz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yz[j] * fx[j] * pb_yy[j] + fx[j] * pa_z[j] * pb_xxy[j] +

                               pa_yz[j] * pb_xxyy[j]) * s_0_0[j];

                t_yz_xxyz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_y[j] * fx[j] * fx[j] * pb_y[j] +

                               0.25 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_yz[j] * fx[j] * pb_yz[j] +

                               0.5 * pa_y[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pa_z[j] * pb_xxz[j] + pa_yz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (13) = (65,70)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, s_0_0, t_yz_xxzz, \
                                     t_yz_xyyy, t_yz_xyyz, t_yz_xyzz, t_yz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_xxzz[j] = (0.25 * pa_yz[j] * fx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_yz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yz[j] * fx[j] * pb_zz[j] + pa_y[j] * fx[j] * pb_xxz[j] +

                               pa_yz[j] * pb_xxzz[j]) * s_0_0[j];

                t_yz_xyyy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 1.5 * pa_yz[j] * pb_xy[j] * fx[j] +

                               1.5 * fx[j] * pa_z[j] * pb_xyy[j] + pa_yz[j] * pb_xyyy[j]) * s_0_0[j];

                t_yz_xyyz[j] = (0.25 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_yz[j] * pb_xz[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_xyy[j] + fx[j] * pa_z[j] * pb_xyz[j] +

                               pa_yz[j] * pb_xyyz[j]) * s_0_0[j];

                t_yz_xyzz[j] = (0.25 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pb_xz[j] +

                               0.5 * pa_yz[j] * pb_xy[j] * fx[j] + pa_y[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_z[j] * pb_xzz[j] +

                               pa_yz[j] * pb_xyzz[j]) * s_0_0[j];

                t_yz_xzzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yz[j] * pb_xz[j] * fx[j] +

                               1.5 * pa_y[j] * fx[j] * pb_xzz[j] + pa_yz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (14) = (70,75)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, \
                                     pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_yz_yyyy, t_yz_yyyz, \
                                     t_yz_yyzz, t_yz_yzzz, t_yz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_yyyy[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               3.0 * pa_yz[j] * pb_yy[j] * fx[j] + 2.0 * fx[j] * pa_z[j] * pb_yyy[j] + pa_yz[j] * pb_yyyy[j]) * s_0_0[j];

                t_yz_yyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_y[j] +

                               0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_yz[j] * pb_yz[j] * fx[j] +

                               0.5 * pa_y[j] * fx[j] * pb_yyy[j] + 1.5 * fx[j] * pa_z[j] * pb_yyz[j] + pa_yz[j] * pb_yyyz[j]) * s_0_0[j];

                t_yz_yyzz[j] = (0.25 * pa_yz[j] * fx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_yz[j] * pb_yy[j] * fx[j] +

                               0.5 * pa_yz[j] * fx[j] * pb_zz[j] + pa_y[j] * fx[j] * pb_yyz[j] + fx[j] * pa_z[j] * pb_yzz[j] + pa_yz[j] * pb_yyzz[j]) * s_0_0[j];

                t_yz_yzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_y[j] +

                               0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_yz[j] * pb_yz[j] * fx[j] +

                               1.5 * pa_y[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_z[j] * pb_zzz[j] + pa_yz[j] * pb_yzzz[j]) * s_0_0[j];

                t_yz_zzzz[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_z[j] +

                               3.0 * pa_yz[j] * pb_zz[j] * fx[j] + 2.0 * pa_y[j] * fx[j] * pb_zzz[j] + pa_yz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (15) = (75,80)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, s_0_0, t_zz_xxxx, t_zz_xxxy, t_zz_xxxz, \
                                     t_zz_xxyy, t_zz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zz_xxxx[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] +

                               1.5 * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_zz[j] * pb_xx[j] * fx[j] + 0.5 * fx[j] * pb_xxxx[j] + pa_zz[j] * pb_xxxx[j]) * s_0_0[j];

                t_zz_xxxy[j] = (0.75 * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_zz[j] * pb_xy[j] * fx[j] +

                               0.5 * fx[j] * pb_xxxy[j] + pa_zz[j] * pb_xxxy[j]) * s_0_0[j];

                t_zz_xxxz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pb_xz[j] +

                               1.5 * pa_zz[j] * pb_xz[j] * fx[j] + pa_z[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pb_xxxz[j] + pa_zz[j] * pb_xxxz[j]) * s_0_0[j];

                t_zz_xxyy[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_zz[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_zz[j] * pb_xx[j] * fx[j] +

                               0.5 * pa_zz[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pb_xxyy[j] + pa_zz[j] * pb_xxyy[j]) * s_0_0[j];

                t_zz_xxyz[j] = (0.5 * pa_z[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_yz[j] +

                               0.5 * pa_zz[j] * fx[j] * pb_yz[j] + pa_z[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pb_xxyz[j] + pa_zz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (16) = (80,85)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, \
                                     pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, s_0_0, t_zz_xxzz, t_zz_xyyy, \
                                     t_zz_xyyz, t_zz_xyzz, t_zz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zz_xxzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.25 * pa_zz[j] * fx[j] * fx[j] +

                               pa_z[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] +

                               0.5 * pa_zz[j] * pb_xx[j] * fx[j] + 0.5 * pa_zz[j] * fx[j] * pb_zz[j] + 2.0 * pa_z[j] * fx[j] * pb_xxz[j] +

                               0.5 * fx[j] * pb_xxzz[j] + pa_zz[j] * pb_xxzz[j]) * s_0_0[j];

                t_zz_xyyy[j] = (0.75 * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_zz[j] * pb_xy[j] * fx[j] +

                               0.5 * fx[j] * pb_xyyy[j] + pa_zz[j] * pb_xyyy[j]) * s_0_0[j];

                t_zz_xyyz[j] = (0.5 * pa_z[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pb_xz[j] +

                               0.5 * pa_zz[j] * pb_xz[j] * fx[j] + pa_z[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pb_xyyz[j] + pa_zz[j] * pb_xyyz[j]) * s_0_0[j];

                t_zz_xyzz[j] = (0.75 * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_zz[j] * pb_xy[j] * fx[j] +

                               2.0 * pa_z[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pb_xyzz[j] + pa_zz[j] * pb_xyzz[j]) * s_0_0[j];

                t_zz_xzzz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pb_xz[j] +

                               1.5 * pa_zz[j] * pb_xz[j] * fx[j] + 3.0 * pa_z[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pb_xzzz[j] + pa_zz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (17) = (85,90)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, \
                                     pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_zz_yyyy, t_zz_yyyz, \
                                     t_zz_yyzz, t_zz_yzzz, t_zz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zz_yyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] +

                               1.5 * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_zz[j] * pb_yy[j] * fx[j] + 0.5 * fx[j] * pb_yyyy[j] + pa_zz[j] * pb_yyyy[j]) * s_0_0[j];

                t_zz_yyyz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_yz[j] +

                               1.5 * pa_zz[j] * pb_yz[j] * fx[j] + pa_z[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pb_yyyz[j] + pa_zz[j] * pb_yyyz[j]) * s_0_0[j];

                t_zz_yyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.25 * pa_zz[j] * fx[j] * fx[j] +

                               pa_z[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] +

                               0.5 * pa_zz[j] * pb_yy[j] * fx[j] + 0.5 * pa_zz[j] * fx[j] * pb_zz[j] + 2.0 * pa_z[j] * fx[j] * pb_yyz[j] +

                               0.5 * fx[j] * pb_yyzz[j] + pa_zz[j] * pb_yyzz[j]) * s_0_0[j];

                t_zz_yzzz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * fx[j] * fx[j] * pb_yz[j] +

                               1.5 * pa_zz[j] * pb_yz[j] * fx[j] + 3.0 * pa_z[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pb_yzzz[j] + pa_zz[j] * pb_yzzz[j]) * s_0_0[j];

                t_zz_zzzz[j] = (1.875 * fx[j] * fx[j] * fx[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] +

                               6.0 * pa_z[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_zz[j] * pb_zz[j] * fx[j] +

                               4.0 * pa_z[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pb_zzzz[j] + pa_zz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForFS(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxx_0 = primBuffer.data(10 * idx);

            auto t_xxy_0 = primBuffer.data(10 * idx + 1);

            auto t_xxz_0 = primBuffer.data(10 * idx + 2);

            auto t_xyy_0 = primBuffer.data(10 * idx + 3);

            auto t_xyz_0 = primBuffer.data(10 * idx + 4);

            auto t_xzz_0 = primBuffer.data(10 * idx + 5);

            auto t_yyy_0 = primBuffer.data(10 * idx + 6);

            auto t_yyz_0 = primBuffer.data(10 * idx + 7);

            auto t_yzz_0 = primBuffer.data(10 * idx + 8);

            auto t_zzz_0 = primBuffer.data(10 * idx + 9);

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxy, pa_xxz, pa_xyy, pa_xyz, pa_xzz, pa_y, pa_yyy, pa_yyz, \
                                     pa_yzz, pa_z, pa_zzz, s_0_0, t_xxx_0, t_xxy_0, t_xxz_0, t_xyy_0, t_xyz_0, t_xzz_0, \
                                     t_yyy_0, t_yyz_0, t_yzz_0, t_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_0[j] = (1.5 * pa_x[j] * fx[j] + pa_xxx[j]) * s_0_0[j];

                t_xxy_0[j] = (0.5 * fx[j] * pa_y[j] + pa_xxy[j]) * s_0_0[j];

                t_xxz_0[j] = (0.5 * fx[j] * pa_z[j] + pa_xxz[j]) * s_0_0[j];

                t_xyy_0[j] = (0.5 * pa_x[j] * fx[j] + pa_xyy[j]) * s_0_0[j];

                t_xyz_0[j] = pa_xyz[j] * s_0_0[j];

                t_xzz_0[j] = (0.5 * pa_x[j] * fx[j] + pa_xzz[j]) * s_0_0[j];

                t_yyy_0[j] = (1.5 * pa_y[j] * fx[j] + pa_yyy[j]) * s_0_0[j];

                t_yyz_0[j] = (0.5 * fx[j] * pa_z[j] + pa_yyz[j]) * s_0_0[j];

                t_yzz_0[j] = (0.5 * pa_y[j] * fx[j] + pa_yzz[j]) * s_0_0[j];

                t_zzz_0[j] = (1.5 * pa_z[j] * fx[j] + pa_zzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForFP(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxx_x = primBuffer.data(30 * idx);

            auto t_xxx_y = primBuffer.data(30 * idx + 1);

            auto t_xxx_z = primBuffer.data(30 * idx + 2);

            auto t_xxy_x = primBuffer.data(30 * idx + 3);

            auto t_xxy_y = primBuffer.data(30 * idx + 4);

            auto t_xxy_z = primBuffer.data(30 * idx + 5);

            auto t_xxz_x = primBuffer.data(30 * idx + 6);

            auto t_xxz_y = primBuffer.data(30 * idx + 7);

            auto t_xxz_z = primBuffer.data(30 * idx + 8);

            auto t_xyy_x = primBuffer.data(30 * idx + 9);

            auto t_xyy_y = primBuffer.data(30 * idx + 10);

            auto t_xyy_z = primBuffer.data(30 * idx + 11);

            auto t_xyz_x = primBuffer.data(30 * idx + 12);

            auto t_xyz_y = primBuffer.data(30 * idx + 13);

            auto t_xyz_z = primBuffer.data(30 * idx + 14);

            auto t_xzz_x = primBuffer.data(30 * idx + 15);

            auto t_xzz_y = primBuffer.data(30 * idx + 16);

            auto t_xzz_z = primBuffer.data(30 * idx + 17);

            auto t_yyy_x = primBuffer.data(30 * idx + 18);

            auto t_yyy_y = primBuffer.data(30 * idx + 19);

            auto t_yyy_z = primBuffer.data(30 * idx + 20);

            auto t_yyz_x = primBuffer.data(30 * idx + 21);

            auto t_yyz_y = primBuffer.data(30 * idx + 22);

            auto t_yyz_z = primBuffer.data(30 * idx + 23);

            auto t_yzz_x = primBuffer.data(30 * idx + 24);

            auto t_yzz_y = primBuffer.data(30 * idx + 25);

            auto t_yzz_z = primBuffer.data(30 * idx + 26);

            auto t_zzz_x = primBuffer.data(30 * idx + 27);

            auto t_zzz_y = primBuffer.data(30 * idx + 28);

            auto t_zzz_z = primBuffer.data(30 * idx + 29);

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxy, pa_xxz, pa_xy, pa_xyy, pa_xz, pa_y, pa_yy, pa_z, \
                                     pb_x, pb_y, pb_z, s_0_0, t_xxx_x, t_xxx_y, t_xxx_z, t_xxy_x, t_xxy_y, t_xxy_z, \
                                     t_xxz_x, t_xxz_y, t_xxz_z, t_xyy_x: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_x[j] = (0.75 * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] + 1.5 * pa_x[j] * fx[j] * pb_x[j] +

                             pa_xxx[j] * pb_x[j]) * s_0_0[j];

                t_xxx_y[j] = (1.5 * pa_x[j] * fx[j] * pb_y[j] + pa_xxx[j] * pb_y[j]) * s_0_0[j];

                t_xxx_z[j] = (1.5 * pa_x[j] * fx[j] * pb_z[j] + pa_xxx[j] * pb_z[j]) * s_0_0[j];

                t_xxy_x[j] = (pa_xy[j] * fx[j] + 0.5 * fx[j] * pa_y[j] * pb_x[j] + pa_xxy[j] * pb_x[j]) * s_0_0[j];

                t_xxy_y[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] + 0.5 * fx[j] * pa_y[j] * pb_y[j] +

                             pa_xxy[j] * pb_y[j]) * s_0_0[j];

                t_xxy_z[j] = (0.5 * fx[j] * pa_y[j] * pb_z[j] + pa_xxy[j] * pb_z[j]) * s_0_0[j];

                t_xxz_x[j] = (pa_xz[j] * fx[j] + 0.5 * fx[j] * pa_z[j] * pb_x[j] + pa_xxz[j] * pb_x[j]) * s_0_0[j];

                t_xxz_y[j] = (0.5 * fx[j] * pa_z[j] * pb_y[j] + pa_xxz[j] * pb_y[j]) * s_0_0[j];

                t_xxz_z[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] + 0.5 * fx[j] * pa_z[j] * pb_z[j] +

                             pa_xxz[j] * pb_z[j]) * s_0_0[j];

                t_xyy_x[j] = (0.25 * fx[j] * fx[j] + 0.5 * fx[j] * pa_yy[j] + 0.5 * pa_x[j] * fx[j] * pb_x[j] +

                             pa_xyy[j] * pb_x[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyz, pa_xz, pa_xzz, pa_y, pa_yy, pa_yyy, pa_yz, pa_zz, \
                                     pb_x, pb_y, pb_z, s_0_0, t_xyy_y, t_xyy_z, t_xyz_x, t_xyz_y, t_xyz_z, t_xzz_x, \
                                     t_xzz_y, t_xzz_z, t_yyy_x, t_yyy_y: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_y[j] = (pa_xy[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_y[j] + pa_xyy[j] * pb_y[j]) * s_0_0[j];

                t_xyy_z[j] = (0.5 * pa_x[j] * fx[j] * pb_z[j] + pa_xyy[j] * pb_z[j]) * s_0_0[j];

                t_xyz_x[j] = (0.5 * fx[j] * pa_yz[j] + pa_xyz[j] * pb_x[j]) * s_0_0[j];

                t_xyz_y[j] = (0.5 * pa_xz[j] * fx[j] + pa_xyz[j] * pb_y[j]) * s_0_0[j];

                t_xyz_z[j] = (0.5 * pa_xy[j] * fx[j] + pa_xyz[j] * pb_z[j]) * s_0_0[j];

                t_xzz_x[j] = (0.25 * fx[j] * fx[j] + 0.5 * fx[j] * pa_zz[j] + 0.5 * pa_x[j] * fx[j] * pb_x[j] +

                             pa_xzz[j] * pb_x[j]) * s_0_0[j];

                t_xzz_y[j] = (0.5 * pa_x[j] * fx[j] * pb_y[j] + pa_xzz[j] * pb_y[j]) * s_0_0[j];

                t_xzz_z[j] = (pa_xz[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_z[j] + pa_xzz[j] * pb_z[j]) * s_0_0[j];

                t_yyy_x[j] = (1.5 * pa_y[j] * fx[j] * pb_x[j] + pa_yyy[j] * pb_x[j]) * s_0_0[j];

                t_yyy_y[j] = (0.75 * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] + 1.5 * pa_y[j] * fx[j] * pb_y[j] +

                             pa_yyy[j] * pb_y[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, pa_zzz, pb_x, pb_y, \
                                     pb_z, s_0_0, t_yyy_z, t_yyz_x, t_yyz_y, t_yyz_z, t_yzz_x, t_yzz_y, t_yzz_z, \
                                     t_zzz_x, t_zzz_y, t_zzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_z[j] = (1.5 * pa_y[j] * fx[j] * pb_z[j] + pa_yyy[j] * pb_z[j]) * s_0_0[j];

                t_yyz_x[j] = (0.5 * fx[j] * pa_z[j] * pb_x[j] + pa_yyz[j] * pb_x[j]) * s_0_0[j];

                t_yyz_y[j] = (pa_yz[j] * fx[j] + 0.5 * fx[j] * pa_z[j] * pb_y[j] + pa_yyz[j] * pb_y[j]) * s_0_0[j];

                t_yyz_z[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] + 0.5 * fx[j] * pa_z[j] * pb_z[j] +

                             pa_yyz[j] * pb_z[j]) * s_0_0[j];

                t_yzz_x[j] = (0.5 * pa_y[j] * fx[j] * pb_x[j] + pa_yzz[j] * pb_x[j]) * s_0_0[j];

                t_yzz_y[j] = (0.25 * fx[j] * fx[j] + 0.5 * fx[j] * pa_zz[j] + 0.5 * pa_y[j] * fx[j] * pb_y[j] +

                             pa_yzz[j] * pb_y[j]) * s_0_0[j];

                t_yzz_z[j] = (pa_yz[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_z[j] + pa_yzz[j] * pb_z[j]) * s_0_0[j];

                t_zzz_x[j] = (1.5 * pa_z[j] * fx[j] * pb_x[j] + pa_zzz[j] * pb_x[j]) * s_0_0[j];

                t_zzz_y[j] = (1.5 * pa_z[j] * fx[j] * pb_y[j] + pa_zzz[j] * pb_y[j]) * s_0_0[j];

                t_zzz_z[j] = (0.75 * fx[j] * fx[j] + 1.5 * pa_zz[j] * fx[j] + 1.5 * pa_z[j] * fx[j] * pb_z[j] +

                             pa_zzz[j] * pb_z[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForFD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxx_xx = primBuffer.data(60 * idx);

            auto t_xxx_xy = primBuffer.data(60 * idx + 1);

            auto t_xxx_xz = primBuffer.data(60 * idx + 2);

            auto t_xxx_yy = primBuffer.data(60 * idx + 3);

            auto t_xxx_yz = primBuffer.data(60 * idx + 4);

            auto t_xxx_zz = primBuffer.data(60 * idx + 5);

            auto t_xxy_xx = primBuffer.data(60 * idx + 6);

            auto t_xxy_xy = primBuffer.data(60 * idx + 7);

            auto t_xxy_xz = primBuffer.data(60 * idx + 8);

            auto t_xxy_yy = primBuffer.data(60 * idx + 9);

            auto t_xxy_yz = primBuffer.data(60 * idx + 10);

            auto t_xxy_zz = primBuffer.data(60 * idx + 11);

            auto t_xxz_xx = primBuffer.data(60 * idx + 12);

            auto t_xxz_xy = primBuffer.data(60 * idx + 13);

            auto t_xxz_xz = primBuffer.data(60 * idx + 14);

            auto t_xxz_yy = primBuffer.data(60 * idx + 15);

            auto t_xxz_yz = primBuffer.data(60 * idx + 16);

            auto t_xxz_zz = primBuffer.data(60 * idx + 17);

            auto t_xyy_xx = primBuffer.data(60 * idx + 18);

            auto t_xyy_xy = primBuffer.data(60 * idx + 19);

            auto t_xyy_xz = primBuffer.data(60 * idx + 20);

            auto t_xyy_yy = primBuffer.data(60 * idx + 21);

            auto t_xyy_yz = primBuffer.data(60 * idx + 22);

            auto t_xyy_zz = primBuffer.data(60 * idx + 23);

            auto t_xyz_xx = primBuffer.data(60 * idx + 24);

            auto t_xyz_xy = primBuffer.data(60 * idx + 25);

            auto t_xyz_xz = primBuffer.data(60 * idx + 26);

            auto t_xyz_yy = primBuffer.data(60 * idx + 27);

            auto t_xyz_yz = primBuffer.data(60 * idx + 28);

            auto t_xyz_zz = primBuffer.data(60 * idx + 29);

            auto t_xzz_xx = primBuffer.data(60 * idx + 30);

            auto t_xzz_xy = primBuffer.data(60 * idx + 31);

            auto t_xzz_xz = primBuffer.data(60 * idx + 32);

            auto t_xzz_yy = primBuffer.data(60 * idx + 33);

            auto t_xzz_yz = primBuffer.data(60 * idx + 34);

            auto t_xzz_zz = primBuffer.data(60 * idx + 35);

            auto t_yyy_xx = primBuffer.data(60 * idx + 36);

            auto t_yyy_xy = primBuffer.data(60 * idx + 37);

            auto t_yyy_xz = primBuffer.data(60 * idx + 38);

            auto t_yyy_yy = primBuffer.data(60 * idx + 39);

            auto t_yyy_yz = primBuffer.data(60 * idx + 40);

            auto t_yyy_zz = primBuffer.data(60 * idx + 41);

            auto t_yyz_xx = primBuffer.data(60 * idx + 42);

            auto t_yyz_xy = primBuffer.data(60 * idx + 43);

            auto t_yyz_xz = primBuffer.data(60 * idx + 44);

            auto t_yyz_yy = primBuffer.data(60 * idx + 45);

            auto t_yyz_yz = primBuffer.data(60 * idx + 46);

            auto t_yyz_zz = primBuffer.data(60 * idx + 47);

            auto t_yzz_xx = primBuffer.data(60 * idx + 48);

            auto t_yzz_xy = primBuffer.data(60 * idx + 49);

            auto t_yzz_xz = primBuffer.data(60 * idx + 50);

            auto t_yzz_yy = primBuffer.data(60 * idx + 51);

            auto t_yzz_yz = primBuffer.data(60 * idx + 52);

            auto t_yzz_zz = primBuffer.data(60 * idx + 53);

            auto t_zzz_xx = primBuffer.data(60 * idx + 54);

            auto t_zzz_xy = primBuffer.data(60 * idx + 55);

            auto t_zzz_xz = primBuffer.data(60 * idx + 56);

            auto t_zzz_yy = primBuffer.data(60 * idx + 57);

            auto t_zzz_yz = primBuffer.data(60 * idx + 58);

            auto t_zzz_zz = primBuffer.data(60 * idx + 59);

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, \
                                     t_xxx_xx, t_xxx_xy, t_xxx_xz, t_xxx_yy, t_xxx_yz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xx[j] = (2.25 * pa_x[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xxx[j] * fx[j] + 3.0 * pa_xx[j] * fx[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * pb_xx[j] +

                              pa_xxx[j] * pb_xx[j]) * s_0_0[j];

                t_xxx_xy[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xx[j] * fx[j] * pb_y[j] +

                              1.5 * pa_x[j] * fx[j] * pb_xy[j] + pa_xxx[j] * pb_xy[j]) * s_0_0[j];

                t_xxx_xz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xx[j] * fx[j] * pb_z[j] +

                              1.5 * pa_x[j] * fx[j] * pb_xz[j] + pa_xxx[j] * pb_xz[j]) * s_0_0[j];

                t_xxx_yy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] +

                              1.5 * pa_x[j] * fx[j] * pb_yy[j] + pa_xxx[j] * pb_yy[j]) * s_0_0[j];

                t_xxx_yz[j] = (1.5 * pa_x[j] * fx[j] * pb_yz[j] + pa_xxx[j] * pb_yz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_z, pb_zz, s_0_0, t_xxx_zz, t_xxy_xx, t_xxy_xy, t_xxy_xz, t_xxy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_zz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] +

                              1.5 * pa_x[j] * fx[j] * pb_zz[j] + pa_xxx[j] * pb_zz[j]) * s_0_0[j];

                t_xxy_xx[j] = (0.75 * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xxy[j] * fx[j] +

                              2.0 * pa_xy[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_y[j] * pb_xx[j] + pa_xxy[j] * pb_xx[j]) * s_0_0[j];

                t_xxy_xy[j] = (0.5 * pa_x[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xx[j] * fx[j] * pb_x[j] + pa_xy[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_y[j] * pb_xy[j] + pa_xxy[j] * pb_xy[j]) * s_0_0[j];

                t_xxy_xz[j] = (pa_xy[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_y[j] * pb_xz[j] + pa_xxy[j] * pb_xz[j]) * s_0_0[j];

                t_xxy_yy[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_xxy[j] * fx[j] + pa_xx[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_y[j] * pb_yy[j] + pa_xxy[j] * pb_yy[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxz, pa_xz, pa_y, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_xxy_yz, t_xxy_zz, t_xxz_xx, t_xxz_xy, t_xxz_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_yz[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xx[j] * fx[j] * pb_z[j] +

                              0.5 * fx[j] * pa_y[j] * pb_yz[j] + pa_xxy[j] * pb_yz[j]) * s_0_0[j];

                t_xxy_zz[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xxy[j] * fx[j] +

                              0.5 * fx[j] * pa_y[j] * pb_zz[j] + pa_xxy[j] * pb_zz[j]) * s_0_0[j];

                t_xxz_xx[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_xxz[j] * fx[j] +

                              2.0 * pa_xz[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_z[j] * pb_xx[j] + pa_xxz[j] * pb_xx[j]) * s_0_0[j];

                t_xxz_xy[j] = (pa_xz[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_z[j] * pb_xy[j] + pa_xxz[j] * pb_xy[j]) * s_0_0[j];

                t_xxz_xz[j] = (0.5 * pa_x[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xx[j] * fx[j] * pb_x[j] + pa_xz[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_z[j] * pb_xz[j] + pa_xxz[j] * pb_xz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xy, pa_xyy, pa_y, pa_yy, pa_z, pb_x, pb_xx, pb_xy, pb_y, \
                                     pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xxz_yy, t_xxz_yz, t_xxz_zz, t_xyy_xx, t_xyy_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_yy[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_xxz[j] * fx[j] +

                              0.5 * fx[j] * pa_z[j] * pb_yy[j] + pa_xxz[j] * pb_yy[j]) * s_0_0[j];

                t_xxz_yz[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xx[j] * fx[j] * pb_y[j] +

                              0.5 * fx[j] * pa_z[j] * pb_yz[j] + pa_xxz[j] * pb_yz[j]) * s_0_0[j];

                t_xxz_zz[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_xxz[j] * fx[j] + pa_xx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_z[j] * pb_zz[j] + pa_xxz[j] * pb_zz[j]) * s_0_0[j];

                t_xyy_xx[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xyy[j] * fx[j] + fx[j] * pa_yy[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * pb_xx[j] + pa_xyy[j] * pb_xx[j]) * s_0_0[j];

                t_xyy_xy[j] = (0.5 * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * pb_y[j] +

                              pa_xy[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_yy[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * pb_xy[j] +

                              pa_xyy[j] * pb_xy[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyz, pa_yy, pa_yz, pb_x, pb_xx, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_xyy_xz, t_xyy_yy, t_xyy_yz, t_xyy_zz, t_xyz_xx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xz[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_yy[j] * pb_z[j] +

                              0.5 * pa_x[j] * fx[j] * pb_xz[j] + pa_xyy[j] * pb_xz[j]) * s_0_0[j];

                t_xyy_yy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] +

                              2.0 * pa_xy[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * pb_yy[j] + pa_xyy[j] * pb_yy[j]) * s_0_0[j];

                t_xyy_yz[j] = (pa_xy[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * pb_yz[j] + pa_xyy[j] * pb_yz[j]) * s_0_0[j];

                t_xyy_zz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] +

                              0.5 * pa_x[j] * fx[j] * pb_zz[j] + pa_xyy[j] * pb_zz[j]) * s_0_0[j];

                t_xyz_xx[j] = (0.5 * pa_xyz[j] * fx[j] + fx[j] * pa_yz[j] * pb_x[j] + pa_xyz[j] * pb_xx[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_xyz_xy, t_xyz_xz, t_xyz_yy, t_xyz_yz, t_xyz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xy[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_xz[j] * fx[j] * pb_x[j] +

                              0.5 * fx[j] * pa_yz[j] * pb_y[j] + pa_xyz[j] * pb_xy[j]) * s_0_0[j];

                t_xyz_xz[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xy[j] * fx[j] * pb_x[j] +

                              0.5 * fx[j] * pa_yz[j] * pb_z[j] + pa_xyz[j] * pb_xz[j]) * s_0_0[j];

                t_xyz_yy[j] = (0.5 * pa_xyz[j] * fx[j] + pa_xz[j] * fx[j] * pb_y[j] + pa_xyz[j] * pb_yy[j]) * s_0_0[j];

                t_xyz_yz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * pb_y[j] +

                              0.5 * pa_xz[j] * fx[j] * pb_z[j] + pa_xyz[j] * pb_yz[j]) * s_0_0[j];

                t_xyz_zz[j] = (0.5 * pa_xyz[j] * fx[j] + pa_xy[j] * fx[j] * pb_z[j] + pa_xyz[j] * pb_zz[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, s_0_0, t_xzz_xx, t_xzz_xy, t_xzz_xz, t_xzz_yy, t_xzz_yz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xx[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xzz[j] * fx[j] + fx[j] * pa_zz[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * pb_xx[j] + pa_xzz[j] * pb_xx[j]) * s_0_0[j];

                t_xzz_xy[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_zz[j] * pb_y[j] +

                              0.5 * pa_x[j] * fx[j] * pb_xy[j] + pa_xzz[j] * pb_xy[j]) * s_0_0[j];

                t_xzz_xz[j] = (0.5 * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * pb_z[j] +

                              pa_xz[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_zz[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * pb_xz[j] +

                              pa_xzz[j] * pb_xz[j]) * s_0_0[j];

                t_xzz_yy[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] +

                              0.5 * pa_x[j] * fx[j] * pb_yy[j] + pa_xzz[j] * pb_yy[j]) * s_0_0[j];

                t_xzz_yz[j] = (pa_xz[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * pb_yz[j] + pa_xzz[j] * pb_yz[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_z, pb_zz, s_0_0, t_xzz_zz, t_yyy_xx, t_yyy_xy, t_yyy_xz, t_yyy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_zz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] +

                              2.0 * pa_xz[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * pb_zz[j] + pa_xzz[j] * pb_zz[j]) * s_0_0[j];

                t_yyy_xx[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] +

                              1.5 * pa_y[j] * fx[j] * pb_xx[j] + pa_yyy[j] * pb_xx[j]) * s_0_0[j];

                t_yyy_xy[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yy[j] * fx[j] * pb_x[j] +

                              1.5 * pa_y[j] * fx[j] * pb_xy[j] + pa_yyy[j] * pb_xy[j]) * s_0_0[j];

                t_yyy_xz[j] = (1.5 * pa_y[j] * fx[j] * pb_xz[j] + pa_yyy[j] * pb_xz[j]) * s_0_0[j];

                t_yyy_yy[j] = (2.25 * pa_y[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_yyy[j] * fx[j] + 3.0 * pa_yy[j] * fx[j] * pb_y[j] + 1.5 * pa_y[j] * fx[j] * pb_yy[j] +

                              pa_yyy[j] * pb_yy[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_yyy_yz, t_yyy_zz, t_yyz_xx, t_yyz_xy, t_yyz_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_yz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_yy[j] * fx[j] * pb_z[j] +

                              1.5 * pa_y[j] * fx[j] * pb_yz[j] + pa_yyy[j] * pb_yz[j]) * s_0_0[j];

                t_yyy_zz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] +

                              1.5 * pa_y[j] * fx[j] * pb_zz[j] + pa_yyy[j] * pb_zz[j]) * s_0_0[j];

                t_yyz_xx[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_yyz[j] * fx[j] +

                              0.5 * fx[j] * pa_z[j] * pb_xx[j] + pa_yyz[j] * pb_xx[j]) * s_0_0[j];

                t_yyz_xy[j] = (pa_yz[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_z[j] * pb_xy[j] + pa_yyz[j] * pb_xy[j]) * s_0_0[j];

                t_yyz_xz[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_yy[j] * fx[j] * pb_x[j] +

                              0.5 * fx[j] * pa_z[j] * pb_xz[j] + pa_yyz[j] * pb_xz[j]) * s_0_0[j];
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_y, pb_yy, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_yyz_yy, t_yyz_yz, t_yyz_zz, t_yzz_xx, t_yzz_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_yy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_yyz[j] * fx[j] +

                              2.0 * pa_yz[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_z[j] * pb_yy[j] + pa_yyz[j] * pb_yy[j]) * s_0_0[j];

                t_yyz_yz[j] = (0.5 * pa_y[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_yy[j] * fx[j] * pb_y[j] + pa_yz[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_z[j] * pb_yz[j] + pa_yyz[j] * pb_yz[j]) * s_0_0[j];

                t_yyz_zz[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_yyz[j] * fx[j] + pa_yy[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_z[j] * pb_zz[j] + pa_yyz[j] * pb_zz[j]) * s_0_0[j];

                t_yzz_xx[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * pa_yzz[j] * fx[j] +

                              0.5 * pa_y[j] * fx[j] * pb_xx[j] + pa_yzz[j] * pb_xx[j]) * s_0_0[j];

                t_yzz_xy[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_zz[j] * pb_x[j] +

                              0.5 * pa_y[j] * fx[j] * pb_xy[j] + pa_yzz[j] * pb_xy[j]) * s_0_0[j];
            }

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, pb_zz, s_0_0, t_yzz_xz, t_yzz_yy, t_yzz_yz, t_yzz_zz, t_zzz_xx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xz[j] = (pa_yz[j] * fx[j] * pb_x[j] + 0.5 * pa_y[j] * fx[j] * pb_xz[j] + pa_yzz[j] * pb_xz[j]) * s_0_0[j];

                t_yzz_yy[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_yzz[j] * fx[j] + fx[j] * pa_zz[j] * pb_y[j] + 0.5 * pa_y[j] * fx[j] * pb_yy[j] + pa_yzz[j] * pb_yy[j]) * s_0_0[j];

                t_yzz_yz[j] = (0.5 * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * pb_z[j] +

                              pa_yz[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_zz[j] * pb_z[j] + 0.5 * pa_y[j] * fx[j] * pb_yz[j] +

                              pa_yzz[j] * pb_yz[j]) * s_0_0[j];

                t_yzz_zz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 0.5 * pa_yzz[j] * fx[j] +

                              2.0 * pa_yz[j] * fx[j] * pb_z[j] + 0.5 * pa_y[j] * fx[j] * pb_zz[j] + pa_yzz[j] * pb_zz[j]) * s_0_0[j];

                t_zzz_xx[j] = (0.75 * pa_z[j] * fx[j] * fx[j] + 0.5 * pa_zzz[j] * fx[j] +

                              1.5 * pa_z[j] * fx[j] * pb_xx[j] + pa_zzz[j] * pb_xx[j]) * s_0_0[j];
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_x, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, \
                                     t_zzz_xy, t_zzz_xz, t_zzz_yy, t_zzz_yz, t_zzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xy[j] = (1.5 * pa_z[j] * fx[j] * pb_xy[j] + pa_zzz[j] * pb_xy[j]) * s_0_0[j];

                t_zzz_xz[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_zz[j] * fx[j] * pb_x[j] +

                              1.5 * pa_z[j] * fx[j] * pb_xz[j] + pa_zzz[j] * pb_xz[j]) * s_0_0[j];

                t_zzz_yy[j] = (0.75 * pa_z[j] * fx[j] * fx[j] + 0.5 * pa_zzz[j] * fx[j] +

                              1.5 * pa_z[j] * fx[j] * pb_yy[j] + pa_zzz[j] * pb_yy[j]) * s_0_0[j];

                t_zzz_yz[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_zz[j] * fx[j] * pb_y[j] +

                              1.5 * pa_z[j] * fx[j] * pb_yz[j] + pa_zzz[j] * pb_yz[j]) * s_0_0[j];

                t_zzz_zz[j] = (2.25 * pa_z[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_zzz[j] * fx[j] + 3.0 * pa_zz[j] * fx[j] * pb_z[j] + 1.5 * pa_z[j] * fx[j] * pb_zz[j] +

                              pa_zzz[j] * pb_zz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForFF(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxx_xxx = primBuffer.data(100 * idx);

            auto t_xxx_xxy = primBuffer.data(100 * idx + 1);

            auto t_xxx_xxz = primBuffer.data(100 * idx + 2);

            auto t_xxx_xyy = primBuffer.data(100 * idx + 3);

            auto t_xxx_xyz = primBuffer.data(100 * idx + 4);

            auto t_xxx_xzz = primBuffer.data(100 * idx + 5);

            auto t_xxx_yyy = primBuffer.data(100 * idx + 6);

            auto t_xxx_yyz = primBuffer.data(100 * idx + 7);

            auto t_xxx_yzz = primBuffer.data(100 * idx + 8);

            auto t_xxx_zzz = primBuffer.data(100 * idx + 9);

            auto t_xxy_xxx = primBuffer.data(100 * idx + 10);

            auto t_xxy_xxy = primBuffer.data(100 * idx + 11);

            auto t_xxy_xxz = primBuffer.data(100 * idx + 12);

            auto t_xxy_xyy = primBuffer.data(100 * idx + 13);

            auto t_xxy_xyz = primBuffer.data(100 * idx + 14);

            auto t_xxy_xzz = primBuffer.data(100 * idx + 15);

            auto t_xxy_yyy = primBuffer.data(100 * idx + 16);

            auto t_xxy_yyz = primBuffer.data(100 * idx + 17);

            auto t_xxy_yzz = primBuffer.data(100 * idx + 18);

            auto t_xxy_zzz = primBuffer.data(100 * idx + 19);

            auto t_xxz_xxx = primBuffer.data(100 * idx + 20);

            auto t_xxz_xxy = primBuffer.data(100 * idx + 21);

            auto t_xxz_xxz = primBuffer.data(100 * idx + 22);

            auto t_xxz_xyy = primBuffer.data(100 * idx + 23);

            auto t_xxz_xyz = primBuffer.data(100 * idx + 24);

            auto t_xxz_xzz = primBuffer.data(100 * idx + 25);

            auto t_xxz_yyy = primBuffer.data(100 * idx + 26);

            auto t_xxz_yyz = primBuffer.data(100 * idx + 27);

            auto t_xxz_yzz = primBuffer.data(100 * idx + 28);

            auto t_xxz_zzz = primBuffer.data(100 * idx + 29);

            auto t_xyy_xxx = primBuffer.data(100 * idx + 30);

            auto t_xyy_xxy = primBuffer.data(100 * idx + 31);

            auto t_xyy_xxz = primBuffer.data(100 * idx + 32);

            auto t_xyy_xyy = primBuffer.data(100 * idx + 33);

            auto t_xyy_xyz = primBuffer.data(100 * idx + 34);

            auto t_xyy_xzz = primBuffer.data(100 * idx + 35);

            auto t_xyy_yyy = primBuffer.data(100 * idx + 36);

            auto t_xyy_yyz = primBuffer.data(100 * idx + 37);

            auto t_xyy_yzz = primBuffer.data(100 * idx + 38);

            auto t_xyy_zzz = primBuffer.data(100 * idx + 39);

            auto t_xyz_xxx = primBuffer.data(100 * idx + 40);

            auto t_xyz_xxy = primBuffer.data(100 * idx + 41);

            auto t_xyz_xxz = primBuffer.data(100 * idx + 42);

            auto t_xyz_xyy = primBuffer.data(100 * idx + 43);

            auto t_xyz_xyz = primBuffer.data(100 * idx + 44);

            auto t_xyz_xzz = primBuffer.data(100 * idx + 45);

            auto t_xyz_yyy = primBuffer.data(100 * idx + 46);

            auto t_xyz_yyz = primBuffer.data(100 * idx + 47);

            auto t_xyz_yzz = primBuffer.data(100 * idx + 48);

            auto t_xyz_zzz = primBuffer.data(100 * idx + 49);

            auto t_xzz_xxx = primBuffer.data(100 * idx + 50);

            auto t_xzz_xxy = primBuffer.data(100 * idx + 51);

            auto t_xzz_xxz = primBuffer.data(100 * idx + 52);

            auto t_xzz_xyy = primBuffer.data(100 * idx + 53);

            auto t_xzz_xyz = primBuffer.data(100 * idx + 54);

            auto t_xzz_xzz = primBuffer.data(100 * idx + 55);

            auto t_xzz_yyy = primBuffer.data(100 * idx + 56);

            auto t_xzz_yyz = primBuffer.data(100 * idx + 57);

            auto t_xzz_yzz = primBuffer.data(100 * idx + 58);

            auto t_xzz_zzz = primBuffer.data(100 * idx + 59);

            auto t_yyy_xxx = primBuffer.data(100 * idx + 60);

            auto t_yyy_xxy = primBuffer.data(100 * idx + 61);

            auto t_yyy_xxz = primBuffer.data(100 * idx + 62);

            auto t_yyy_xyy = primBuffer.data(100 * idx + 63);

            auto t_yyy_xyz = primBuffer.data(100 * idx + 64);

            auto t_yyy_xzz = primBuffer.data(100 * idx + 65);

            auto t_yyy_yyy = primBuffer.data(100 * idx + 66);

            auto t_yyy_yyz = primBuffer.data(100 * idx + 67);

            auto t_yyy_yzz = primBuffer.data(100 * idx + 68);

            auto t_yyy_zzz = primBuffer.data(100 * idx + 69);

            auto t_yyz_xxx = primBuffer.data(100 * idx + 70);

            auto t_yyz_xxy = primBuffer.data(100 * idx + 71);

            auto t_yyz_xxz = primBuffer.data(100 * idx + 72);

            auto t_yyz_xyy = primBuffer.data(100 * idx + 73);

            auto t_yyz_xyz = primBuffer.data(100 * idx + 74);

            auto t_yyz_xzz = primBuffer.data(100 * idx + 75);

            auto t_yyz_yyy = primBuffer.data(100 * idx + 76);

            auto t_yyz_yyz = primBuffer.data(100 * idx + 77);

            auto t_yyz_yzz = primBuffer.data(100 * idx + 78);

            auto t_yyz_zzz = primBuffer.data(100 * idx + 79);

            auto t_yzz_xxx = primBuffer.data(100 * idx + 80);

            auto t_yzz_xxy = primBuffer.data(100 * idx + 81);

            auto t_yzz_xxz = primBuffer.data(100 * idx + 82);

            auto t_yzz_xyy = primBuffer.data(100 * idx + 83);

            auto t_yzz_xyz = primBuffer.data(100 * idx + 84);

            auto t_yzz_xzz = primBuffer.data(100 * idx + 85);

            auto t_yzz_yyy = primBuffer.data(100 * idx + 86);

            auto t_yzz_yyz = primBuffer.data(100 * idx + 87);

            auto t_yzz_yzz = primBuffer.data(100 * idx + 88);

            auto t_yzz_zzz = primBuffer.data(100 * idx + 89);

            auto t_zzz_xxx = primBuffer.data(100 * idx + 90);

            auto t_zzz_xxy = primBuffer.data(100 * idx + 91);

            auto t_zzz_xxz = primBuffer.data(100 * idx + 92);

            auto t_zzz_xyy = primBuffer.data(100 * idx + 93);

            auto t_zzz_xyz = primBuffer.data(100 * idx + 94);

            auto t_zzz_xzz = primBuffer.data(100 * idx + 95);

            auto t_zzz_yyy = primBuffer.data(100 * idx + 96);

            auto t_zzz_yyz = primBuffer.data(100 * idx + 97);

            auto t_zzz_yzz = primBuffer.data(100 * idx + 98);

            auto t_zzz_zzz = primBuffer.data(100 * idx + 99);

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxx_xxx, t_xxx_xxy, t_xxx_xxz, \
                                     t_xxx_xyy, t_xxx_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xxx[j] = (1.875 * fx[j] * fx[j] * fx[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] +

                               6.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xxx[j] * pb_x[j] * fx[j] +

                               4.5 * pa_xx[j] * fx[j] * pb_xx[j] + 1.5 * pa_x[j] * fx[j] * pb_xxx[j] + pa_xxx[j] * pb_xxx[j]) * s_0_0[j];

                t_xxx_xxy[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_xxx[j] * fx[j] * pb_y[j] + 3.0 * pa_xx[j] * fx[j] * pb_xy[j] + 1.5 * pa_x[j] * fx[j] * pb_xxy[j] +

                               pa_xxx[j] * pb_xxy[j]) * s_0_0[j];

                t_xxx_xxz[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pb_xz[j] +

                               0.5 * pa_xxx[j] * fx[j] * pb_z[j] + 3.0 * pa_xx[j] * fx[j] * pb_xz[j] + 1.5 * pa_x[j] * fx[j] * pb_xxz[j] +

                               pa_xxx[j] * pb_xxz[j]) * s_0_0[j];

                t_xxx_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xxx[j] * pb_x[j] * fx[j] +

                               1.5 * pa_xx[j] * fx[j] * pb_yy[j] + 1.5 * pa_x[j] * fx[j] * pb_xyy[j] + pa_xxx[j] * pb_xyy[j]) * s_0_0[j];

                t_xxx_xyz[j] = (0.75 * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xx[j] * fx[j] * pb_yz[j] +

                               1.5 * pa_x[j] * fx[j] * pb_xyz[j] + pa_xxx[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, s_0_0, t_xxx_xzz, t_xxx_yyy, t_xxx_yyz, t_xxx_yzz, t_xxx_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xxx[j] * pb_x[j] * fx[j] +

                               1.5 * pa_xx[j] * fx[j] * pb_zz[j] + 1.5 * pa_x[j] * fx[j] * pb_xzz[j] + pa_xxx[j] * pb_xzz[j]) * s_0_0[j];

                t_xxx_yyy[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xxx[j] * pb_y[j] * fx[j] +

                               1.5 * pa_x[j] * fx[j] * pb_yyy[j] + pa_xxx[j] * pb_yyy[j]) * s_0_0[j];

                t_xxx_yyz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xxx[j] * fx[j] * pb_z[j] +

                               1.5 * pa_x[j] * fx[j] * pb_yyz[j] + pa_xxx[j] * pb_yyz[j]) * s_0_0[j];

                t_xxx_yzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xxx[j] * pb_y[j] * fx[j] +

                               1.5 * pa_x[j] * fx[j] * pb_yzz[j] + pa_xxx[j] * pb_yzz[j]) * s_0_0[j];

                t_xxx_zzz[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xxx[j] * pb_z[j] * fx[j] +

                               1.5 * pa_x[j] * fx[j] * pb_zzz[j] + pa_xxx[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxy_xxx, t_xxy_xxy, \
                                     t_xxy_xxz, t_xxy_xyy, t_xxy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xxx[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               1.5 * pa_xxy[j] * pb_x[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pa_y[j] * pb_xxx[j] +

                               pa_xxy[j] * pb_xxx[j]) * s_0_0[j];

                t_xxy_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] +

                               0.5 * pa_xxy[j] * fx[j] * pb_y[j] + 0.5 * pa_xx[j] * fx[j] * pb_xx[j] + 2.0 * pa_xy[j] * fx[j] * pb_xy[j] +

                               0.5 * fx[j] * pa_y[j] * pb_xxy[j] + pa_xxy[j] * pb_xxy[j]) * s_0_0[j];

                t_xxy_xxz[j] = (0.75 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 0.5 * pa_xxy[j] * fx[j] * pb_z[j] +

                               2.0 * pa_xy[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_y[j] * pb_xxz[j] + pa_xxy[j] * pb_xxz[j]) * s_0_0[j];

                t_xxy_xyy[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] + pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.25 * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xxy[j] * pb_x[j] * fx[j] +

                               pa_xx[j] * fx[j] * pb_xy[j] + pa_xy[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_y[j] * pb_xyy[j] +

                               pa_xxy[j] * pb_xyy[j]) * s_0_0[j];

                t_xxy_xyz[j] = (0.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_xz[j] +

                               0.5 * pa_xx[j] * fx[j] * pb_xz[j] + pa_xy[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_y[j] * pb_xyz[j] +

                               pa_xxy[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxy_xzz, t_xxy_yyy, t_xxy_yyz, t_xxy_yzz, \
                                     t_xxy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xzz[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               0.5 * pa_xxy[j] * pb_x[j] * fx[j] + pa_xy[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_y[j] * pb_xzz[j] +

                               pa_xxy[j] * pb_xzz[j]) * s_0_0[j];

                t_xxy_yyy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               0.75 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xxy[j] * pb_y[j] * fx[j] +

                               1.5 * pa_xx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_y[j] * pb_yyy[j] + pa_xxy[j] * pb_yyy[j]) * s_0_0[j];

                t_xxy_yyz[j] = (0.25 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pb_yz[j] +

                               0.5 * pa_xxy[j] * fx[j] * pb_z[j] + pa_xx[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_y[j] * pb_yyz[j] +

                               pa_xxy[j] * pb_yyz[j]) * s_0_0[j];

                t_xxy_yzz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xxy[j] * pb_y[j] * fx[j] +

                               0.5 * pa_xx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_y[j] * pb_yzz[j] + pa_xxy[j] * pb_yzz[j]) * s_0_0[j];

                t_xxy_zzz[j] = (0.75 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 1.5 * pa_xxy[j] * pb_z[j] * fx[j] +

                               0.5 * fx[j] * pa_y[j] * pb_zzz[j] + pa_xxy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxz_xxx, t_xxz_xxy, \
                                     t_xxz_xxz, t_xxz_xyy, t_xxz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xxx[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               1.5 * pa_xxz[j] * pb_x[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pa_z[j] * pb_xxx[j] +

                               pa_xxz[j] * pb_xxx[j]) * s_0_0[j];

                t_xxz_xxy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.5 * pa_xxz[j] * fx[j] * pb_y[j] +

                               2.0 * pa_xz[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_z[j] * pb_xxy[j] + pa_xxz[j] * pb_xxy[j]) * s_0_0[j];

                t_xxz_xxz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] +

                               0.5 * pa_xxz[j] * fx[j] * pb_z[j] + 0.5 * pa_xx[j] * fx[j] * pb_xx[j] + 2.0 * pa_xz[j] * fx[j] * pb_xz[j] +

                               0.5 * fx[j] * pa_z[j] * pb_xxz[j] + pa_xxz[j] * pb_xxz[j]) * s_0_0[j];

                t_xxz_xyy[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               0.5 * pa_xxz[j] * pb_x[j] * fx[j] + pa_xz[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_z[j] * pb_xyy[j] +

                               pa_xxz[j] * pb_xyy[j]) * s_0_0[j];

                t_xxz_xyz[j] = (0.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_xx[j] * fx[j] * pb_xy[j] + pa_xz[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_z[j] * pb_xyz[j] +

                               pa_xxz[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxz_xzz, t_xxz_yyy, t_xxz_yyz, \
                                     t_xxz_yzz, t_xxz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xzz[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] + pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               0.25 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xxz[j] * pb_x[j] * fx[j] +

                               pa_xx[j] * fx[j] * pb_xz[j] + pa_xz[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_z[j] * pb_xzz[j] +

                               pa_xxz[j] * pb_xzz[j]) * s_0_0[j];

                t_xxz_yyy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 1.5 * pa_xxz[j] * pb_y[j] * fx[j] +

                               0.5 * fx[j] * pa_z[j] * pb_yyy[j] + pa_xxz[j] * pb_yyy[j]) * s_0_0[j];

                t_xxz_yyz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xxz[j] * fx[j] * pb_z[j] +

                               0.5 * pa_xx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_z[j] * pb_yyz[j] + pa_xxz[j] * pb_yyz[j]) * s_0_0[j];

                t_xxz_yzz[j] = (0.25 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pb_yz[j] +

                               0.5 * pa_xxz[j] * pb_y[j] * fx[j] + pa_xx[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_z[j] * pb_yzz[j] +

                               pa_xxz[j] * pb_yzz[j]) * s_0_0[j];

                t_xxz_zzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xxz[j] * pb_z[j] * fx[j] +

                               1.5 * pa_xx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_z[j] * pb_zzz[j] + pa_xxz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xyy_xxx, t_xyy_xxy, \
                                     t_xyy_xxz, t_xyy_xyy, t_xyy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xxx[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] +

                               0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xyy[j] * pb_x[j] * fx[j] +

                               1.5 * fx[j] * pa_yy[j] * pb_xx[j] + 0.5 * pa_x[j] * fx[j] * pb_xxx[j] + pa_xyy[j] * pb_xxx[j]) * s_0_0[j];

                t_xyy_xxy[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] + fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               0.25 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xyy[j] * fx[j] * pb_y[j] +

                               pa_xy[j] * fx[j] * pb_xx[j] + fx[j] * pa_yy[j] * pb_xy[j] + 0.5 * pa_x[j] * fx[j] * pb_xxy[j] +

                               pa_xyy[j] * pb_xxy[j]) * s_0_0[j];

                t_xyy_xxz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pb_xz[j] +

                               0.5 * pa_xyy[j] * fx[j] * pb_z[j] + fx[j] * pa_yy[j] * pb_xz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxz[j] +

                               pa_xyy[j] * pb_xxz[j]) * s_0_0[j];

                t_xyy_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] +

                               0.25 * fx[j] * fx[j] * pa_yy[j] + fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] +

                               0.5 * pa_xyy[j] * pb_x[j] * fx[j] + 2.0 * pa_xy[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_yy[j] * pb_yy[j] +

                               0.5 * pa_x[j] * fx[j] * pb_xyy[j] + pa_xyy[j] * pb_xyy[j]) * s_0_0[j];

                t_xyy_xyz[j] = (0.5 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_yz[j] +

                               pa_xy[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_yy[j] * pb_yz[j] + 0.5 * pa_x[j] * fx[j] * pb_xyz[j] +

                               pa_xyy[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_yy, pb_x, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xyy_xzz, t_xyy_yyy, t_xyy_yyz, t_xyy_yzz, \
                                     t_xyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xzz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] +

                               0.25 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xyy[j] * pb_x[j] * fx[j] +

                               0.5 * fx[j] * pa_yy[j] * pb_zz[j] + 0.5 * pa_x[j] * fx[j] * pb_xzz[j] + pa_xyy[j] * pb_xzz[j]) * s_0_0[j];

                t_xyy_yyy[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               1.5 * pa_xyy[j] * pb_y[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * pb_yy[j] + 0.5 * pa_x[j] * fx[j] * pb_yyy[j] +

                               pa_xyy[j] * pb_yyy[j]) * s_0_0[j];

                t_xyy_yyz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xyy[j] * fx[j] * pb_z[j] +

                               2.0 * pa_xy[j] * fx[j] * pb_yz[j] + 0.5 * pa_x[j] * fx[j] * pb_yyz[j] + pa_xyy[j] * pb_yyz[j]) * s_0_0[j];

                t_xyy_yzz[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.5 * pa_xyy[j] * pb_y[j] * fx[j] + pa_xy[j] * fx[j] * pb_zz[j] + 0.5 * pa_x[j] * fx[j] * pb_yzz[j] +

                               pa_xyy[j] * pb_yzz[j]) * s_0_0[j];

                t_xyy_zzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xyy[j] * pb_z[j] * fx[j] +

                               0.5 * pa_x[j] * fx[j] * pb_zzz[j] + pa_xyy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xyz_xxx, \
                                     t_xyz_xxy, t_xyz_xxz, t_xyz_xyy, t_xyz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xxx[j] = (0.75 * fx[j] * fx[j] * pa_yz[j] + 1.5 * pa_xyz[j] * pb_x[j] * fx[j] +

                               1.5 * fx[j] * pa_yz[j] * pb_xx[j] + pa_xyz[j] * pb_xxx[j]) * s_0_0[j];

                t_xyz_xxy[j] = (0.25 * pa_xz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               0.5 * pa_xyz[j] * fx[j] * pb_y[j] + 0.5 * pa_xz[j] * fx[j] * pb_xx[j] + fx[j] * pa_yz[j] * pb_xy[j] +

                               pa_xyz[j] * pb_xxy[j]) * s_0_0[j];

                t_xyz_xxz[j] = (0.25 * pa_xy[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               0.5 * pa_xyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xy[j] * fx[j] * pb_xx[j] + fx[j] * pa_yz[j] * pb_xz[j] +

                               pa_xyz[j] * pb_xxz[j]) * s_0_0[j];

                t_xyz_xyy[j] = (0.25 * fx[j] * fx[j] * pa_yz[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               0.5 * pa_xyz[j] * pb_x[j] * fx[j] + pa_xz[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_yz[j] * pb_yy[j] +

                               pa_xyz[j] * pb_xyy[j]) * s_0_0[j];

                t_xyz_xyz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_x[j] +

                               0.25 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.5 * pa_xy[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_xz[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yz[j] + pa_xyz[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xyz_xzz, t_xyz_yyy, \
                                     t_xyz_yyz, t_xyz_yzz, t_xyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xzz[j] = (0.25 * fx[j] * fx[j] * pa_yz[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                               0.5 * pa_xyz[j] * pb_x[j] * fx[j] + pa_xy[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_yz[j] * pb_zz[j] +

                               pa_xyz[j] * pb_xzz[j]) * s_0_0[j];

                t_xyz_yyy[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] + 1.5 * pa_xyz[j] * pb_y[j] * fx[j] +

                               1.5 * pa_xz[j] * fx[j] * pb_yy[j] + pa_xyz[j] * pb_yyy[j]) * s_0_0[j];

                t_xyz_yyz[j] = (0.25 * pa_xy[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.5 * pa_xyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xy[j] * fx[j] * pb_yy[j] + pa_xz[j] * fx[j] * pb_yz[j] +

                               pa_xyz[j] * pb_yyz[j]) * s_0_0[j];

                t_xyz_yzz[j] = (0.25 * pa_xz[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_xyz[j] * pb_y[j] * fx[j] + pa_xy[j] * fx[j] * pb_yz[j] + 0.5 * pa_xz[j] * fx[j] * pb_zz[j] +

                               pa_xyz[j] * pb_yzz[j]) * s_0_0[j];

                t_xyz_zzz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] + 1.5 * pa_xyz[j] * pb_z[j] * fx[j] +

                               1.5 * pa_xy[j] * fx[j] * pb_zz[j] + pa_xyz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xzz_xxx, t_xzz_xxy, \
                                     t_xzz_xxz, t_xzz_xyy, t_xzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xxx[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] +

                               0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xzz[j] * pb_x[j] * fx[j] +

                               1.5 * fx[j] * pa_zz[j] * pb_xx[j] + 0.5 * pa_x[j] * fx[j] * pb_xxx[j] + pa_xzz[j] * pb_xxx[j]) * s_0_0[j];

                t_xzz_xxy[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_xzz[j] * fx[j] * pb_y[j] + fx[j] * pa_zz[j] * pb_xy[j] + 0.5 * pa_x[j] * fx[j] * pb_xxy[j] +

                               pa_xzz[j] * pb_xxy[j]) * s_0_0[j];

                t_xzz_xxz[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] + fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               0.25 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xzz[j] * fx[j] * pb_z[j] +

                               pa_xz[j] * fx[j] * pb_xx[j] + fx[j] * pa_zz[j] * pb_xz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxz[j] +

                               pa_xzz[j] * pb_xxz[j]) * s_0_0[j];

                t_xzz_xyy[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] +

                               0.25 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xzz[j] * pb_x[j] * fx[j] +

                               0.5 * fx[j] * pa_zz[j] * pb_yy[j] + 0.5 * pa_x[j] * fx[j] * pb_xyy[j] + pa_xzz[j] * pb_xyy[j]) * s_0_0[j];

                t_xzz_xyz[j] = (0.5 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_yz[j] +

                               pa_xz[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_zz[j] * pb_yz[j] + 0.5 * pa_x[j] * fx[j] * pb_xyz[j] +

                               pa_xzz[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xzz_xzz, t_xzz_yyy, t_xzz_yyz, \
                                     t_xzz_yzz, t_xzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] +

                               0.25 * fx[j] * fx[j] * pa_zz[j] + fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] +

                               0.5 * pa_xzz[j] * pb_x[j] * fx[j] + 2.0 * pa_xz[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_zz[j] * pb_zz[j] +

                               0.5 * pa_x[j] * fx[j] * pb_xzz[j] + pa_xzz[j] * pb_xzz[j]) * s_0_0[j];

                t_xzz_yyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xzz[j] * pb_y[j] * fx[j] +

                               0.5 * pa_x[j] * fx[j] * pb_yyy[j] + pa_xzz[j] * pb_yyy[j]) * s_0_0[j];

                t_xzz_yyz[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_xzz[j] * fx[j] * pb_z[j] + pa_xz[j] * fx[j] * pb_yy[j] + 0.5 * pa_x[j] * fx[j] * pb_yyz[j] +

                               pa_xzz[j] * pb_yyz[j]) * s_0_0[j];

                t_xzz_yzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xzz[j] * pb_y[j] * fx[j] +

                               2.0 * pa_xz[j] * fx[j] * pb_yz[j] + 0.5 * pa_x[j] * fx[j] * pb_yzz[j] + pa_xzz[j] * pb_yzz[j]) * s_0_0[j];

                t_xzz_zzz[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               1.5 * pa_xzz[j] * pb_z[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * pb_zz[j] + 0.5 * pa_x[j] * fx[j] * pb_zzz[j] +

                               pa_xzz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (12) = (60,65)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_y, pb_z, s_0_0, t_yyy_xxx, t_yyy_xxy, t_yyy_xxz, t_yyy_xyy, \
                                     t_yyy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xxx[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yyy[j] * pb_x[j] * fx[j] +

                               1.5 * pa_y[j] * fx[j] * pb_xxx[j] + pa_yyy[j] * pb_xxx[j]) * s_0_0[j];

                t_yyy_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] +

                               0.75 * pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_yyy[j] * fx[j] * pb_y[j] +

                               1.5 * pa_yy[j] * fx[j] * pb_xx[j] + 1.5 * pa_y[j] * fx[j] * pb_xxy[j] + pa_yyy[j] * pb_xxy[j]) * s_0_0[j];

                t_yyy_xxz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_yyy[j] * fx[j] * pb_z[j] +

                               1.5 * pa_y[j] * fx[j] * pb_xxz[j] + pa_yyy[j] * pb_xxz[j]) * s_0_0[j];

                t_yyy_xyy[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_yyy[j] * pb_x[j] * fx[j] + 3.0 * pa_yy[j] * fx[j] * pb_xy[j] + 1.5 * pa_y[j] * fx[j] * pb_xyy[j] +

                               pa_yyy[j] * pb_xyy[j]) * s_0_0[j];

                t_yyy_xyz[j] = (0.75 * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yy[j] * fx[j] * pb_xz[j] +

                               1.5 * pa_y[j] * fx[j] * pb_xyz[j] + pa_yyy[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (13) = (65,70)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pb_x, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, s_0_0, t_yyy_xzz, t_yyy_yyy, t_yyy_yyz, t_yyy_yzz, \
                                     t_yyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_yyy[j] * pb_x[j] * fx[j] +

                               1.5 * pa_y[j] * fx[j] * pb_xzz[j] + pa_yyy[j] * pb_xzz[j]) * s_0_0[j];

                t_yyy_yyy[j] = (1.875 * fx[j] * fx[j] * fx[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] +

                               6.75 * pa_y[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_yyy[j] * pb_y[j] * fx[j] +

                               4.5 * pa_yy[j] * fx[j] * pb_yy[j] + 1.5 * pa_y[j] * fx[j] * pb_yyy[j] + pa_yyy[j] * pb_yyy[j]) * s_0_0[j];

                t_yyy_yyz[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pb_yz[j] +

                               0.5 * pa_yyy[j] * fx[j] * pb_z[j] + 3.0 * pa_yy[j] * fx[j] * pb_yz[j] + 1.5 * pa_y[j] * fx[j] * pb_yyz[j] +

                               pa_yyy[j] * pb_yyz[j]) * s_0_0[j];

                t_yyy_yzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] +

                               0.75 * pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_yyy[j] * pb_y[j] * fx[j] +

                               1.5 * pa_yy[j] * fx[j] * pb_zz[j] + 1.5 * pa_y[j] * fx[j] * pb_yzz[j] + pa_yyy[j] * pb_yzz[j]) * s_0_0[j];

                t_yyy_zzz[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_yyy[j] * pb_z[j] * fx[j] +

                               1.5 * pa_y[j] * fx[j] * pb_zzz[j] + pa_yyy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (14) = (70,75)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, s_0_0, t_yyz_xxx, t_yyz_xxy, t_yyz_xxz, t_yyz_xyy, \
                                     t_yyz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xxx[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 1.5 * pa_yyz[j] * pb_x[j] * fx[j] +

                               0.5 * fx[j] * pa_z[j] * pb_xxx[j] + pa_yyz[j] * pb_xxx[j]) * s_0_0[j];

                t_yyz_xxy[j] = (0.5 * pa_yz[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               0.5 * pa_yyz[j] * fx[j] * pb_y[j] + pa_yz[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pa_z[j] * pb_xxy[j] +

                               pa_yyz[j] * pb_xxy[j]) * s_0_0[j];

                t_yyz_xxz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_yyz[j] * fx[j] * pb_z[j] +

                               0.5 * pa_yy[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pa_z[j] * pb_xxz[j] + pa_yyz[j] * pb_xxz[j]) * s_0_0[j];

                t_yyz_xyy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.5 * pa_yyz[j] * pb_x[j] * fx[j] +

                               2.0 * pa_yz[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_z[j] * pb_xyy[j] + pa_yyz[j] * pb_xyy[j]) * s_0_0[j];

                t_yyz_xyz[j] = (0.5 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_yy[j] * fx[j] * pb_xy[j] + pa_yz[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_z[j] * pb_xyz[j] +

                               pa_yyz[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (15) = (75,80)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_yyz_xzz, t_yyz_yyy, t_yyz_yyz, \
                                     t_yyz_yzz, t_yyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xzz[j] = (0.25 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pb_xz[j] +

                               0.5 * pa_yyz[j] * pb_x[j] * fx[j] + pa_yy[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_z[j] * pb_xzz[j] +

                               pa_yyz[j] * pb_xzz[j]) * s_0_0[j];

                t_yyz_yyy[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               1.5 * pa_yyz[j] * pb_y[j] * fx[j] + 3.0 * pa_yz[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_z[j] * pb_yyy[j] +

                               pa_yyz[j] * pb_yyy[j]) * s_0_0[j];

                t_yyz_yyz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] +

                               pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] +

                               0.5 * pa_yyz[j] * fx[j] * pb_z[j] + 0.5 * pa_yy[j] * fx[j] * pb_yy[j] + 2.0 * pa_yz[j] * fx[j] * pb_yz[j] +

                               0.5 * fx[j] * pa_z[j] * pb_yyz[j] + pa_yyz[j] * pb_yyz[j]) * s_0_0[j];

                t_yyz_yzz[j] = (0.5 * pa_yz[j] * fx[j] * fx[j] + pa_y[j] * fx[j] * fx[j] * pb_z[j] +

                               0.25 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_yyz[j] * pb_y[j] * fx[j] +

                               pa_yy[j] * fx[j] * pb_yz[j] + pa_yz[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_z[j] * pb_yzz[j] +

                               pa_yyz[j] * pb_yzz[j]) * s_0_0[j];

                t_yyz_zzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] +

                               0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_yyz[j] * pb_z[j] * fx[j] +

                               1.5 * pa_yy[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_z[j] * pb_zzz[j] + pa_yyz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (16) = (80,85)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, s_0_0, t_yzz_xxx, t_yzz_xxy, t_yzz_xxz, t_yzz_xyy, \
                                     t_yzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xxx[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yzz[j] * pb_x[j] * fx[j] +

                               0.5 * pa_y[j] * fx[j] * pb_xxx[j] + pa_yzz[j] * pb_xxx[j]) * s_0_0[j];

                t_yzz_xxy[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] +

                               0.25 * pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_yzz[j] * fx[j] * pb_y[j] +

                               0.5 * fx[j] * pa_zz[j] * pb_xx[j] + 0.5 * pa_y[j] * fx[j] * pb_xxy[j] + pa_yzz[j] * pb_xxy[j]) * s_0_0[j];

                t_yzz_xxz[j] = (0.5 * pa_yz[j] * fx[j] * fx[j] + 0.25 * pa_y[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_yzz[j] * fx[j] * pb_z[j] + pa_yz[j] * fx[j] * pb_xx[j] + 0.5 * pa_y[j] * fx[j] * pb_xxz[j] +

                               pa_yzz[j] * pb_xxz[j]) * s_0_0[j];

                t_yzz_xyy[j] = (0.25 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pb_xy[j] +

                               0.5 * pa_yzz[j] * pb_x[j] * fx[j] + fx[j] * pa_zz[j] * pb_xy[j] + 0.5 * pa_y[j] * fx[j] * pb_xyy[j] +

                               pa_yzz[j] * pb_xyy[j]) * s_0_0[j];

                t_yzz_xyz[j] = (0.5 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pb_xz[j] +

                               pa_yz[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_zz[j] * pb_xz[j] + 0.5 * pa_y[j] * fx[j] * pb_xyz[j] +

                               pa_yzz[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (17) = (85,90)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_yzz_xzz, t_yzz_yyy, t_yzz_yyz, \
                                     t_yzz_yzz, t_yzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_yzz[j] * pb_x[j] * fx[j] +

                               2.0 * pa_yz[j] * fx[j] * pb_xz[j] + 0.5 * pa_y[j] * fx[j] * pb_xzz[j] + pa_yzz[j] * pb_xzz[j]) * s_0_0[j];

                t_yzz_yyy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] +

                               0.75 * pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_yzz[j] * pb_y[j] * fx[j] +

                               1.5 * fx[j] * pa_zz[j] * pb_yy[j] + 0.5 * pa_y[j] * fx[j] * pb_yyy[j] + pa_yzz[j] * pb_yyy[j]) * s_0_0[j];

                t_yzz_yyz[j] = (0.5 * pa_yz[j] * fx[j] * fx[j] + fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               0.25 * pa_y[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_yzz[j] * fx[j] * pb_z[j] +

                               pa_yz[j] * fx[j] * pb_yy[j] + fx[j] * pa_zz[j] * pb_yz[j] + 0.5 * pa_y[j] * fx[j] * pb_yyz[j] +

                               pa_yzz[j] * pb_yyz[j]) * s_0_0[j];

                t_yzz_yzz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_y[j] +

                               0.25 * fx[j] * fx[j] * pa_zz[j] + fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] +

                               0.5 * pa_yzz[j] * pb_y[j] * fx[j] + 2.0 * pa_yz[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_zz[j] * pb_zz[j] +

                               0.5 * pa_y[j] * fx[j] * pb_yzz[j] + pa_yzz[j] * pb_yzz[j]) * s_0_0[j];

                t_yzz_zzz[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] + 2.25 * pa_y[j] * fx[j] * fx[j] * pb_z[j] +

                               1.5 * pa_yzz[j] * pb_z[j] * fx[j] + 3.0 * pa_yz[j] * fx[j] * pb_zz[j] + 0.5 * pa_y[j] * fx[j] * pb_zzz[j] +

                               pa_yzz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (18) = (90,95)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_y, pb_z, s_0_0, t_zzz_xxx, t_zzz_xxy, t_zzz_xxz, t_zzz_xyy, t_zzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xxx[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_zzz[j] * pb_x[j] * fx[j] +

                               1.5 * pa_z[j] * fx[j] * pb_xxx[j] + pa_zzz[j] * pb_xxx[j]) * s_0_0[j];

                t_zzz_xxy[j] = (0.75 * pa_z[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_zzz[j] * fx[j] * pb_y[j] +

                               1.5 * pa_z[j] * fx[j] * pb_xxy[j] + pa_zzz[j] * pb_xxy[j]) * s_0_0[j];

                t_zzz_xxz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] +

                               0.75 * pa_z[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_zzz[j] * fx[j] * pb_z[j] +

                               1.5 * pa_zz[j] * fx[j] * pb_xx[j] + 1.5 * pa_z[j] * fx[j] * pb_xxz[j] + pa_zzz[j] * pb_xxz[j]) * s_0_0[j];

                t_zzz_xyy[j] = (0.75 * pa_z[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_zzz[j] * pb_x[j] * fx[j] +

                               1.5 * pa_z[j] * fx[j] * pb_xyy[j] + pa_zzz[j] * pb_xyy[j]) * s_0_0[j];

                t_zzz_xyz[j] = (0.75 * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_zz[j] * fx[j] * pb_xy[j] +

                               1.5 * pa_z[j] * fx[j] * pb_xyz[j] + pa_zzz[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (19) = (95,100)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_zzz_xzz, t_zzz_yyy, t_zzz_yyz, t_zzz_yzz, \
                                     t_zzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xzz[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * fx[j] * fx[j] * pb_xz[j] +

                               0.5 * pa_zzz[j] * pb_x[j] * fx[j] + 3.0 * pa_zz[j] * fx[j] * pb_xz[j] + 1.5 * pa_z[j] * fx[j] * pb_xzz[j] +

                               pa_zzz[j] * pb_xzz[j]) * s_0_0[j];

                t_zzz_yyy[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_zzz[j] * pb_y[j] * fx[j] +

                               1.5 * pa_z[j] * fx[j] * pb_yyy[j] + pa_zzz[j] * pb_yyy[j]) * s_0_0[j];

                t_zzz_yyz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] +

                               0.75 * pa_z[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_zzz[j] * fx[j] * pb_z[j] +

                               1.5 * pa_zz[j] * fx[j] * pb_yy[j] + 1.5 * pa_z[j] * fx[j] * pb_yyz[j] + pa_zzz[j] * pb_yyz[j]) * s_0_0[j];

                t_zzz_yzz[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pb_yz[j] +

                               0.5 * pa_zzz[j] * pb_y[j] * fx[j] + 3.0 * pa_zz[j] * fx[j] * pb_yz[j] + 1.5 * pa_z[j] * fx[j] * pb_yzz[j] +

                               pa_zzz[j] * pb_yzz[j]) * s_0_0[j];

                t_zzz_zzz[j] = (1.875 * fx[j] * fx[j] * fx[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] +

                               6.75 * pa_z[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_zzz[j] * pb_z[j] * fx[j] +

                               4.5 * pa_zz[j] * fx[j] * pb_zz[j] + 1.5 * pa_z[j] * fx[j] * pb_zzz[j] + pa_zzz[j] * pb_zzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForFG(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(19 * idx);

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(19 * idx + 3);

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            auto pa_yyy = paDistances.data(19 * idx + 15);

            auto pa_yyz = paDistances.data(19 * idx + 16);

            auto pa_yzz = paDistances.data(19 * idx + 17);

            auto pa_zzz = paDistances.data(19 * idx + 18);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxx_xxxx = primBuffer.data(150 * idx);

            auto t_xxx_xxxy = primBuffer.data(150 * idx + 1);

            auto t_xxx_xxxz = primBuffer.data(150 * idx + 2);

            auto t_xxx_xxyy = primBuffer.data(150 * idx + 3);

            auto t_xxx_xxyz = primBuffer.data(150 * idx + 4);

            auto t_xxx_xxzz = primBuffer.data(150 * idx + 5);

            auto t_xxx_xyyy = primBuffer.data(150 * idx + 6);

            auto t_xxx_xyyz = primBuffer.data(150 * idx + 7);

            auto t_xxx_xyzz = primBuffer.data(150 * idx + 8);

            auto t_xxx_xzzz = primBuffer.data(150 * idx + 9);

            auto t_xxx_yyyy = primBuffer.data(150 * idx + 10);

            auto t_xxx_yyyz = primBuffer.data(150 * idx + 11);

            auto t_xxx_yyzz = primBuffer.data(150 * idx + 12);

            auto t_xxx_yzzz = primBuffer.data(150 * idx + 13);

            auto t_xxx_zzzz = primBuffer.data(150 * idx + 14);

            auto t_xxy_xxxx = primBuffer.data(150 * idx + 15);

            auto t_xxy_xxxy = primBuffer.data(150 * idx + 16);

            auto t_xxy_xxxz = primBuffer.data(150 * idx + 17);

            auto t_xxy_xxyy = primBuffer.data(150 * idx + 18);

            auto t_xxy_xxyz = primBuffer.data(150 * idx + 19);

            auto t_xxy_xxzz = primBuffer.data(150 * idx + 20);

            auto t_xxy_xyyy = primBuffer.data(150 * idx + 21);

            auto t_xxy_xyyz = primBuffer.data(150 * idx + 22);

            auto t_xxy_xyzz = primBuffer.data(150 * idx + 23);

            auto t_xxy_xzzz = primBuffer.data(150 * idx + 24);

            auto t_xxy_yyyy = primBuffer.data(150 * idx + 25);

            auto t_xxy_yyyz = primBuffer.data(150 * idx + 26);

            auto t_xxy_yyzz = primBuffer.data(150 * idx + 27);

            auto t_xxy_yzzz = primBuffer.data(150 * idx + 28);

            auto t_xxy_zzzz = primBuffer.data(150 * idx + 29);

            auto t_xxz_xxxx = primBuffer.data(150 * idx + 30);

            auto t_xxz_xxxy = primBuffer.data(150 * idx + 31);

            auto t_xxz_xxxz = primBuffer.data(150 * idx + 32);

            auto t_xxz_xxyy = primBuffer.data(150 * idx + 33);

            auto t_xxz_xxyz = primBuffer.data(150 * idx + 34);

            auto t_xxz_xxzz = primBuffer.data(150 * idx + 35);

            auto t_xxz_xyyy = primBuffer.data(150 * idx + 36);

            auto t_xxz_xyyz = primBuffer.data(150 * idx + 37);

            auto t_xxz_xyzz = primBuffer.data(150 * idx + 38);

            auto t_xxz_xzzz = primBuffer.data(150 * idx + 39);

            auto t_xxz_yyyy = primBuffer.data(150 * idx + 40);

            auto t_xxz_yyyz = primBuffer.data(150 * idx + 41);

            auto t_xxz_yyzz = primBuffer.data(150 * idx + 42);

            auto t_xxz_yzzz = primBuffer.data(150 * idx + 43);

            auto t_xxz_zzzz = primBuffer.data(150 * idx + 44);

            auto t_xyy_xxxx = primBuffer.data(150 * idx + 45);

            auto t_xyy_xxxy = primBuffer.data(150 * idx + 46);

            auto t_xyy_xxxz = primBuffer.data(150 * idx + 47);

            auto t_xyy_xxyy = primBuffer.data(150 * idx + 48);

            auto t_xyy_xxyz = primBuffer.data(150 * idx + 49);

            auto t_xyy_xxzz = primBuffer.data(150 * idx + 50);

            auto t_xyy_xyyy = primBuffer.data(150 * idx + 51);

            auto t_xyy_xyyz = primBuffer.data(150 * idx + 52);

            auto t_xyy_xyzz = primBuffer.data(150 * idx + 53);

            auto t_xyy_xzzz = primBuffer.data(150 * idx + 54);

            auto t_xyy_yyyy = primBuffer.data(150 * idx + 55);

            auto t_xyy_yyyz = primBuffer.data(150 * idx + 56);

            auto t_xyy_yyzz = primBuffer.data(150 * idx + 57);

            auto t_xyy_yzzz = primBuffer.data(150 * idx + 58);

            auto t_xyy_zzzz = primBuffer.data(150 * idx + 59);

            auto t_xyz_xxxx = primBuffer.data(150 * idx + 60);

            auto t_xyz_xxxy = primBuffer.data(150 * idx + 61);

            auto t_xyz_xxxz = primBuffer.data(150 * idx + 62);

            auto t_xyz_xxyy = primBuffer.data(150 * idx + 63);

            auto t_xyz_xxyz = primBuffer.data(150 * idx + 64);

            auto t_xyz_xxzz = primBuffer.data(150 * idx + 65);

            auto t_xyz_xyyy = primBuffer.data(150 * idx + 66);

            auto t_xyz_xyyz = primBuffer.data(150 * idx + 67);

            auto t_xyz_xyzz = primBuffer.data(150 * idx + 68);

            auto t_xyz_xzzz = primBuffer.data(150 * idx + 69);

            auto t_xyz_yyyy = primBuffer.data(150 * idx + 70);

            auto t_xyz_yyyz = primBuffer.data(150 * idx + 71);

            auto t_xyz_yyzz = primBuffer.data(150 * idx + 72);

            auto t_xyz_yzzz = primBuffer.data(150 * idx + 73);

            auto t_xyz_zzzz = primBuffer.data(150 * idx + 74);

            auto t_xzz_xxxx = primBuffer.data(150 * idx + 75);

            auto t_xzz_xxxy = primBuffer.data(150 * idx + 76);

            auto t_xzz_xxxz = primBuffer.data(150 * idx + 77);

            auto t_xzz_xxyy = primBuffer.data(150 * idx + 78);

            auto t_xzz_xxyz = primBuffer.data(150 * idx + 79);

            auto t_xzz_xxzz = primBuffer.data(150 * idx + 80);

            auto t_xzz_xyyy = primBuffer.data(150 * idx + 81);

            auto t_xzz_xyyz = primBuffer.data(150 * idx + 82);

            auto t_xzz_xyzz = primBuffer.data(150 * idx + 83);

            auto t_xzz_xzzz = primBuffer.data(150 * idx + 84);

            auto t_xzz_yyyy = primBuffer.data(150 * idx + 85);

            auto t_xzz_yyyz = primBuffer.data(150 * idx + 86);

            auto t_xzz_yyzz = primBuffer.data(150 * idx + 87);

            auto t_xzz_yzzz = primBuffer.data(150 * idx + 88);

            auto t_xzz_zzzz = primBuffer.data(150 * idx + 89);

            auto t_yyy_xxxx = primBuffer.data(150 * idx + 90);

            auto t_yyy_xxxy = primBuffer.data(150 * idx + 91);

            auto t_yyy_xxxz = primBuffer.data(150 * idx + 92);

            auto t_yyy_xxyy = primBuffer.data(150 * idx + 93);

            auto t_yyy_xxyz = primBuffer.data(150 * idx + 94);

            auto t_yyy_xxzz = primBuffer.data(150 * idx + 95);

            auto t_yyy_xyyy = primBuffer.data(150 * idx + 96);

            auto t_yyy_xyyz = primBuffer.data(150 * idx + 97);

            auto t_yyy_xyzz = primBuffer.data(150 * idx + 98);

            auto t_yyy_xzzz = primBuffer.data(150 * idx + 99);

            auto t_yyy_yyyy = primBuffer.data(150 * idx + 100);

            auto t_yyy_yyyz = primBuffer.data(150 * idx + 101);

            auto t_yyy_yyzz = primBuffer.data(150 * idx + 102);

            auto t_yyy_yzzz = primBuffer.data(150 * idx + 103);

            auto t_yyy_zzzz = primBuffer.data(150 * idx + 104);

            auto t_yyz_xxxx = primBuffer.data(150 * idx + 105);

            auto t_yyz_xxxy = primBuffer.data(150 * idx + 106);

            auto t_yyz_xxxz = primBuffer.data(150 * idx + 107);

            auto t_yyz_xxyy = primBuffer.data(150 * idx + 108);

            auto t_yyz_xxyz = primBuffer.data(150 * idx + 109);

            auto t_yyz_xxzz = primBuffer.data(150 * idx + 110);

            auto t_yyz_xyyy = primBuffer.data(150 * idx + 111);

            auto t_yyz_xyyz = primBuffer.data(150 * idx + 112);

            auto t_yyz_xyzz = primBuffer.data(150 * idx + 113);

            auto t_yyz_xzzz = primBuffer.data(150 * idx + 114);

            auto t_yyz_yyyy = primBuffer.data(150 * idx + 115);

            auto t_yyz_yyyz = primBuffer.data(150 * idx + 116);

            auto t_yyz_yyzz = primBuffer.data(150 * idx + 117);

            auto t_yyz_yzzz = primBuffer.data(150 * idx + 118);

            auto t_yyz_zzzz = primBuffer.data(150 * idx + 119);

            auto t_yzz_xxxx = primBuffer.data(150 * idx + 120);

            auto t_yzz_xxxy = primBuffer.data(150 * idx + 121);

            auto t_yzz_xxxz = primBuffer.data(150 * idx + 122);

            auto t_yzz_xxyy = primBuffer.data(150 * idx + 123);

            auto t_yzz_xxyz = primBuffer.data(150 * idx + 124);

            auto t_yzz_xxzz = primBuffer.data(150 * idx + 125);

            auto t_yzz_xyyy = primBuffer.data(150 * idx + 126);

            auto t_yzz_xyyz = primBuffer.data(150 * idx + 127);

            auto t_yzz_xyzz = primBuffer.data(150 * idx + 128);

            auto t_yzz_xzzz = primBuffer.data(150 * idx + 129);

            auto t_yzz_yyyy = primBuffer.data(150 * idx + 130);

            auto t_yzz_yyyz = primBuffer.data(150 * idx + 131);

            auto t_yzz_yyzz = primBuffer.data(150 * idx + 132);

            auto t_yzz_yzzz = primBuffer.data(150 * idx + 133);

            auto t_yzz_zzzz = primBuffer.data(150 * idx + 134);

            auto t_zzz_xxxx = primBuffer.data(150 * idx + 135);

            auto t_zzz_xxxy = primBuffer.data(150 * idx + 136);

            auto t_zzz_xxxz = primBuffer.data(150 * idx + 137);

            auto t_zzz_xxyy = primBuffer.data(150 * idx + 138);

            auto t_zzz_xxyz = primBuffer.data(150 * idx + 139);

            auto t_zzz_xxzz = primBuffer.data(150 * idx + 140);

            auto t_zzz_xyyy = primBuffer.data(150 * idx + 141);

            auto t_zzz_xyyz = primBuffer.data(150 * idx + 142);

            auto t_zzz_xyzz = primBuffer.data(150 * idx + 143);

            auto t_zzz_xzzz = primBuffer.data(150 * idx + 144);

            auto t_zzz_yyyy = primBuffer.data(150 * idx + 145);

            auto t_zzz_yyyz = primBuffer.data(150 * idx + 146);

            auto t_zzz_yyzz = primBuffer.data(150 * idx + 147);

            auto t_zzz_yzzz = primBuffer.data(150 * idx + 148);

            auto t_zzz_zzzz = primBuffer.data(150 * idx + 149);

            // Batch of Integrals (0) = (0,2)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxy, pb_xy, \
                                     pb_y, s_0_0, t_xxx_xxxx, t_xxx_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xxxx[j] = (5.625 * pa_x[j] * fx[j] * fx[j] * fx[j] + 7.5 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xxx[j] * fx[j] * fx[j] + 9.0 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] +

                                13.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xxx[j] * pb_xx[j] * fx[j] +

                                6.0 * pa_xx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_x[j] * fx[j] * pb_xxxx[j] + pa_xxx[j] * pb_xxxx[j]) * s_0_0[j];

                t_xxx_xxxy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                6.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 2.25 * fx[j] * fx[j] * pb_xxy[j] + 1.5 * pa_xxx[j] * pb_xy[j] * fx[j] +

                                4.5 * pa_xx[j] * fx[j] * pb_xxy[j] + 1.5 * pa_x[j] * fx[j] * pb_xxxy[j] + pa_xxx[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (2,4)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxxz, pb_xxyy, pb_xxz, pb_xyy, pb_xz, \
                                     pb_yy, pb_z, s_0_0, t_xxx_xxxz, t_xxx_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xxxz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                6.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 2.25 * fx[j] * fx[j] * pb_xxz[j] + 1.5 * pa_xxx[j] * pb_xz[j] * fx[j] +

                                4.5 * pa_xx[j] * fx[j] * pb_xxz[j] + 1.5 * pa_x[j] * fx[j] * pb_xxxz[j] + pa_xxx[j] * pb_xxxz[j]) * s_0_0[j];

                t_xxx_xxyy[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xxx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] +

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * fx[j] * fx[j] * pb_xyy[j] +

                                0.5 * pa_xxx[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xx[j] * fx[j] * pb_xyy[j] +

                                1.5 * pa_x[j] * fx[j] * pb_xxyy[j] + pa_xxx[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (4,6)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxyz, pb_xxzz, pb_xyz, pb_xzz, pb_yz, \
                                     pb_zz, s_0_0, t_xxx_xxyz, t_xxx_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xxyz[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * fx[j] * fx[j] * pb_xyz[j] +

                                0.5 * pa_xxx[j] * fx[j] * pb_yz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xyz[j] + 1.5 * pa_x[j] * fx[j] * pb_xxyz[j] +

                                pa_xxx[j] * pb_xxyz[j]) * s_0_0[j];

                t_xxx_xxzz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xxx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] +

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * fx[j] * fx[j] * pb_xzz[j] +

                                0.5 * pa_xxx[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xzz[j] +

                                1.5 * pa_x[j] * fx[j] * pb_xxzz[j] + pa_xxx[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (6,8)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pb_xy, pb_xyyy, pb_xyyz, pb_xz, pb_y, pb_yyy, pb_yyz, \
                                     pb_z, s_0_0, t_xxx_xyyy, t_xxx_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xyyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xxx[j] * pb_xy[j] * fx[j] +

                                1.5 * pa_xx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_x[j] * fx[j] * pb_xyyy[j] + pa_xxx[j] * pb_xyyy[j]) * s_0_0[j];

                t_xxx_xyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xxx[j] * pb_xz[j] * fx[j] +

                                1.5 * pa_xx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_x[j] * fx[j] * pb_xyyz[j] + pa_xxx[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (8,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pb_xy, pb_xyzz, pb_xz, pb_xzzz, pb_y, pb_yzz, pb_z, \
                                     pb_zzz, s_0_0, t_xxx_xyzz, t_xxx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xxx[j] * pb_xy[j] * fx[j] +

                                1.5 * pa_xx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_x[j] * fx[j] * pb_xyzz[j] + pa_xxx[j] * pb_xyzz[j]) * s_0_0[j];

                t_xxx_xzzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xxx[j] * pb_xz[j] * fx[j] +

                                1.5 * pa_xx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_x[j] * fx[j] * pb_xzzz[j] + pa_xxx[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (10,12)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pb_yy, pb_yyyy, pb_yyyz, pb_yz, s_0_0, t_xxx_yyyy, \
                                     t_xxx_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_yyyy[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] +

                                4.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xxx[j] * pb_yy[j] * fx[j] + 1.5 * pa_x[j] * fx[j] * pb_yyyy[j] +

                                pa_xxx[j] * pb_yyyy[j]) * s_0_0[j];

                t_xxx_yyyz[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xxx[j] * pb_yz[j] * fx[j] +

                                1.5 * pa_x[j] * fx[j] * pb_yyyz[j] + pa_xxx[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (12,14)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pb_yy, pb_yyzz, pb_yz, pb_yzzz, pb_zz, s_0_0, t_xxx_yyzz, \
                                     t_xxx_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_yyzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xxx[j] * pb_yy[j] * fx[j] +

                                0.5 * pa_xxx[j] * fx[j] * pb_zz[j] + 1.5 * pa_x[j] * fx[j] * pb_yyzz[j] + pa_xxx[j] * pb_yyzz[j]) * s_0_0[j];

                t_xxx_yzzz[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xxx[j] * pb_yz[j] * fx[j] +

                                1.5 * pa_x[j] * fx[j] * pb_yzzz[j] + pa_xxx[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (14,16)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_zz, \
                                     pb_zzzz, s_0_0, t_xxx_zzzz, t_xxy_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_zzzz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] +

                                4.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xxx[j] * pb_zz[j] * fx[j] + 1.5 * pa_x[j] * fx[j] * pb_zzzz[j] +

                                pa_xxx[j] * pb_zzzz[j]) * s_0_0[j];

                t_xxy_xxxx[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] +

                                6.0 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 3.0 * pa_xxy[j] * pb_xx[j] * fx[j] +

                                4.0 * pa_xy[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pa_y[j] * pb_xxxx[j] + pa_xxy[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (16,18)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxz, pb_xy, pb_xz, pb_y, pb_z, s_0_0, t_xxy_xxxy, t_xxy_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xxxy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                2.25 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_xxy[j] * pb_xy[j] * fx[j] +

                                0.5 * pa_xx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xy[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pa_y[j] * pb_xxxy[j] +

                                pa_xxy[j] * pb_xxxy[j]) * s_0_0[j];

                t_xxy_xxxz[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                2.25 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 1.5 * pa_xxy[j] * pb_xz[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * pb_xxz[j] +

                                0.5 * fx[j] * pa_y[j] * pb_xxxz[j] + pa_xxy[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (9) = (18,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxy, pb_xxyy, pb_xxyz, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxy_xxyy, \
                                     t_xxy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xxyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_xxy[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + pa_xy[j] * fx[j] * fx[j] * pb_x[j] +

                                2.0 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] +

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_xxy[j] * pb_xx[j] * fx[j] +

                                0.5 * pa_xxy[j] * fx[j] * pb_yy[j] + pa_xx[j] * fx[j] * pb_xxy[j] + 2.0 * pa_xy[j] * fx[j] * pb_xyy[j] +

                                0.5 * fx[j] * pa_y[j] * pb_xxyy[j] + pa_xxy[j] * pb_xxyy[j]) * s_0_0[j];

                t_xxy_xxyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] +

                                0.5 * pa_xxy[j] * fx[j] * pb_yz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxz[j] + 2.0 * pa_xy[j] * fx[j] * pb_xyz[j] +

                                0.5 * fx[j] * pa_y[j] * pb_xxyz[j] + pa_xxy[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (10) = (20,22)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_zz, s_0_0, t_xxy_xxzz, t_xxy_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xxzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * pa_xxy[j] * fx[j] * fx[j] +

                                pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] +

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.5 * pa_xxy[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxy[j] * fx[j] * pb_zz[j] +

                                2.0 * pa_xy[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_y[j] * pb_xxzz[j] + pa_xxy[j] * pb_xxzz[j]) * s_0_0[j];

                t_xxy_xyyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] +

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] + 1.5 * pa_xxy[j] * pb_xy[j] * fx[j] +

                                1.5 * pa_xx[j] * fx[j] * pb_xyy[j] + pa_xy[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_y[j] * pb_xyyy[j] +

                                pa_xxy[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (11) = (22,24)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, \
                                     pb_xzz, pb_y, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, t_xxy_xyyz, t_xxy_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xyyz[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + pa_x[j] * fx[j] * fx[j] * pb_yz[j] +

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xxy[j] * pb_xz[j] * fx[j] +

                                pa_xx[j] * fx[j] * pb_xyz[j] + pa_xy[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pa_y[j] * pb_xyyz[j] +

                                pa_xxy[j] * pb_xyyz[j]) * s_0_0[j];

                t_xxy_xyzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] +

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_xxy[j] * pb_xy[j] * fx[j] +

                                0.5 * pa_xx[j] * fx[j] * pb_xzz[j] + pa_xy[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_y[j] * pb_xyzz[j] +

                                pa_xxy[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (12) = (24,26)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xy, pa_y, pb_xz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_z, \
                                     pb_zzz, s_0_0, t_xxy_xzzz, t_xxy_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xzzz[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 1.5 * pa_xxy[j] * pb_xz[j] * fx[j] + pa_xy[j] * fx[j] * pb_zzz[j] +

                                0.5 * fx[j] * pa_y[j] * pb_xzzz[j] + pa_xxy[j] * pb_xzzz[j]) * s_0_0[j];

                t_xxy_yyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_xxy[j] * fx[j] * fx[j] + 3.0 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] +

                                fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_xxy[j] * pb_yy[j] * fx[j] + 2.0 * pa_xx[j] * fx[j] * pb_yyy[j] +

                                0.5 * fx[j] * pa_y[j] * pb_yyyy[j] + pa_xxy[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (13) = (26,28)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_y, pb_y, pb_yy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, s_0_0, t_xxy_yyyz, t_xxy_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_yyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_xxy[j] * pb_yz[j] * fx[j] +

                                1.5 * pa_xx[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pa_y[j] * pb_yyyz[j] + pa_xxy[j] * pb_yyyz[j]) * s_0_0[j];

                t_xxy_yyzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_xxy[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_yzz[j] +

                                0.5 * pa_xxy[j] * pb_yy[j] * fx[j] + 0.5 * pa_xxy[j] * fx[j] * pb_zz[j] + pa_xx[j] * fx[j] * pb_yzz[j] +

                                0.5 * fx[j] * pa_y[j] * pb_yyzz[j] + pa_xxy[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (14) = (28,30)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_y, pb_yz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, \
                                     t_xxy_yzzz, t_xxy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_yzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xxy[j] * pb_yz[j] * fx[j] +

                                0.5 * pa_xx[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_y[j] * pb_yzzz[j] + pa_xxy[j] * pb_yzzz[j]) * s_0_0[j];

                t_xxy_zzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] +

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 3.0 * pa_xxy[j] * pb_zz[j] * fx[j] + 0.5 * fx[j] * pa_y[j] * pb_zzzz[j] +

                                pa_xxy[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (15) = (30,32)

            #pragma omp simd aligned(fx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxy, pb_xy, \
                                     pb_y, s_0_0, t_xxz_xxxx, t_xxz_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xxxx[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] +

                                6.0 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 3.0 * pa_xxz[j] * pb_xx[j] * fx[j] +

                                4.0 * pa_xz[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pa_z[j] * pb_xxxx[j] + pa_xxz[j] * pb_xxxx[j]) * s_0_0[j];

                t_xxz_xxxy[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                2.25 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 1.5 * pa_xxz[j] * pb_xy[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * pb_xxy[j] +

                                0.5 * fx[j] * pa_z[j] * pb_xxxy[j] + pa_xxz[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (16) = (32,34)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxz, pb_xxyy, \
                                     pb_xxz, pb_xyy, pb_xz, pb_yy, pb_z, s_0_0, t_xxz_xxxz, t_xxz_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xxxz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                2.25 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_xxz[j] * pb_xz[j] * fx[j] +

                                0.5 * pa_xx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xz[j] * fx[j] * pb_xxz[j] + 0.5 * fx[j] * pa_z[j] * pb_xxxz[j] +

                                pa_xxz[j] * pb_xxxz[j]) * s_0_0[j];

                t_xxz_xxyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * pa_xxz[j] * fx[j] * fx[j] +

                                pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_xxz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxz[j] * fx[j] * pb_yy[j] +

                                2.0 * pa_xz[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_z[j] * pb_xxyy[j] + pa_xxz[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (17) = (34,36)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxy, pb_xxyz, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_xxz_xxyz, \
                                     t_xxz_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xxyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] +

                                0.5 * pa_xxz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxy[j] + 2.0 * pa_xz[j] * fx[j] * pb_xyz[j] +

                                0.5 * fx[j] * pa_z[j] * pb_xxyz[j] + pa_xxz[j] * pb_xxyz[j]) * s_0_0[j];

                t_xxz_xxzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_xxz[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + pa_xz[j] * fx[j] * fx[j] * pb_x[j] +

                                2.0 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_xxz[j] * pb_xx[j] * fx[j] +

                                0.5 * pa_xxz[j] * fx[j] * pb_zz[j] + pa_xx[j] * fx[j] * pb_xxz[j] + 2.0 * pa_xz[j] * fx[j] * pb_xzz[j] +

                                0.5 * fx[j] * pa_z[j] * pb_xxzz[j] + pa_xxz[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (18) = (36,38)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xz, \
                                     pb_y, pb_yy, pb_yyy, pb_yyz, pb_z, s_0_0, t_xxz_xyyy, t_xxz_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xyyy[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 1.5 * pa_xxz[j] * pb_xy[j] * fx[j] + pa_xz[j] * fx[j] * pb_yyy[j] +

                                0.5 * fx[j] * pa_z[j] * pb_xyyy[j] + pa_xxz[j] * pb_xyyy[j]) * s_0_0[j];

                t_xxz_xyyz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_xxz[j] * pb_xz[j] * fx[j] +

                                0.5 * pa_xx[j] * fx[j] * pb_xyy[j] + pa_xz[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pa_z[j] * pb_xyyz[j] +

                                pa_xxz[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (19) = (38,40)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xy, pb_xyz, pb_xyzz, pb_xz, pb_xzz, \
                                     pb_xzzz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxz_xyzz, t_xxz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xyzz[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + pa_x[j] * fx[j] * fx[j] * pb_yz[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xxz[j] * pb_xy[j] * fx[j] +

                                pa_xx[j] * fx[j] * pb_xyz[j] + pa_xz[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_z[j] * pb_xyzz[j] +

                                pa_xxz[j] * pb_xyzz[j]) * s_0_0[j];

                t_xxz_xzzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] + 1.5 * pa_xxz[j] * pb_xz[j] * fx[j] +

                                1.5 * pa_xx[j] * fx[j] * pb_xzz[j] + pa_xz[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_z[j] * pb_xzzz[j] +

                                pa_xxz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (20) = (40,42)

            #pragma omp simd aligned(fx, pa_xx, pa_xxz, pa_z, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yz, s_0_0, \
                                     t_xxz_yyyy, t_xxz_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_yyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 3.0 * pa_xxz[j] * pb_yy[j] * fx[j] + 0.5 * fx[j] * pa_z[j] * pb_yyyy[j] +

                                pa_xxz[j] * pb_yyyy[j]) * s_0_0[j];

                t_xxz_yyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xxz[j] * pb_yz[j] * fx[j] +

                                0.5 * pa_xx[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_z[j] * pb_yyyz[j] + pa_xxz[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (21) = (42,44)

            #pragma omp simd aligned(fx, pa_xx, pa_xxz, pa_z, pb_y, pb_yy, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, \
                                     pb_z, pb_zz, s_0_0, t_xxz_yyzz, t_xxz_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_yyzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_xxz[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_yyz[j] +

                                0.5 * pa_xxz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xxz[j] * fx[j] * pb_zz[j] + pa_xx[j] * fx[j] * pb_yyz[j] +

                                0.5 * fx[j] * pa_z[j] * pb_yyzz[j] + pa_xxz[j] * pb_yyzz[j]) * s_0_0[j];

                t_xxz_yzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_xxz[j] * pb_yz[j] * fx[j] +

                                1.5 * pa_xx[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_z[j] * pb_yzzz[j] + pa_xxz[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (22) = (44,46)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xyy, pa_yy, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_z, \
                                     pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xxz_zzzz, t_xyy_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_zzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_xxz[j] * fx[j] * fx[j] + 3.0 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] +

                                fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_xxz[j] * pb_zz[j] * fx[j] + 2.0 * pa_xx[j] * fx[j] * pb_zzz[j] +

                                0.5 * fx[j] * pa_z[j] * pb_zzzz[j] + pa_xxz[j] * pb_zzzz[j]) * s_0_0[j];

                t_xyy_xxxx[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xyy[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xyy[j] * pb_xx[j] * fx[j] + 2.0 * fx[j] * pa_yy[j] * pb_xxx[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xxxx[j] + pa_xyy[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (23) = (46,48)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxz, pb_xy, pb_xz, pb_y, pb_z, s_0_0, t_xyy_xxxy, t_xyy_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xxxy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] +

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] +

                                1.5 * pa_xyy[j] * pb_xy[j] * fx[j] + pa_xy[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_yy[j] * pb_xxy[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xxxy[j] + pa_xyy[j] * pb_xxxy[j]) * s_0_0[j];

                t_xyy_xxxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] + 1.5 * pa_xyy[j] * pb_xz[j] * fx[j] +

                                1.5 * fx[j] * pa_yy[j] * pb_xxz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxxz[j] + pa_xyy[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (24) = (48,50)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxy, pb_xxyy, pb_xxyz, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xyy_xxyy, \
                                     t_xyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xxyy[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xyy[j] * fx[j] * fx[j] + pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                0.5 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 2.0 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_xyy[j] * pb_xx[j] * fx[j] +

                                0.5 * pa_xyy[j] * fx[j] * pb_yy[j] + 2.0 * pa_xy[j] * fx[j] * pb_xxy[j] + fx[j] * pa_yy[j] * pb_xyy[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xxyy[j] + pa_xyy[j] * pb_xxyy[j]) * s_0_0[j];

                t_xyy_xxyz[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + fx[j] * fx[j] * pa_y[j] * pb_xz[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xyy[j] * fx[j] * pb_yz[j] +

                                pa_xy[j] * fx[j] * pb_xxz[j] + fx[j] * pa_yy[j] * pb_xyz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxyz[j] +

                                pa_xyy[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (25) = (50,52)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_zz, s_0_0, t_xyy_xxzz, t_xyy_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xxzz[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xyy[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_xzz[j] +

                                0.5 * pa_xyy[j] * pb_xx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * pb_zz[j] + fx[j] * pa_yy[j] * pb_xzz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xxzz[j] + pa_xyy[j] * pb_xxzz[j]) * s_0_0[j];

                t_xyy_xyyy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] +

                                0.75 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] +

                                1.5 * pa_xyy[j] * pb_xy[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_yy[j] * pb_yyy[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xyyy[j] + pa_xyy[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (26) = (52,54)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, \
                                     pb_xzz, pb_y, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, t_xyy_xyyz, t_xyy_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] +

                                0.25 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] +

                                0.5 * pa_xyy[j] * pb_xz[j] * fx[j] + 2.0 * pa_xy[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yy[j] * pb_yyz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xyyz[j] + pa_xyy[j] * pb_xyyz[j]) * s_0_0[j];

                t_xyy_xyzz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.5 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] +

                                0.5 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_yzz[j] +

                                0.5 * pa_xyy[j] * pb_xy[j] * fx[j] + pa_xy[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yy[j] * pb_yzz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xyzz[j] + pa_xyy[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (27) = (54,56)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_yy, pb_xz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_z, \
                                     pb_zzz, s_0_0, t_xyy_xzzz, t_xyy_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xyy[j] * pb_xz[j] * fx[j] +

                                0.5 * fx[j] * pa_yy[j] * pb_zzz[j] + 0.5 * pa_x[j] * fx[j] * pb_xzzz[j] + pa_xyy[j] * pb_xzzz[j]) * s_0_0[j];

                t_xyy_yyyy[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] +

                                6.0 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xyy[j] * pb_yy[j] * fx[j] +

                                4.0 * pa_xy[j] * fx[j] * pb_yyy[j] + 0.5 * pa_x[j] * fx[j] * pb_yyyy[j] + pa_xyy[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (28) = (56,58)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pb_y, pb_yy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, s_0_0, t_xyy_yyyz, t_xyy_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_yyyz[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xyy[j] * pb_yz[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * pb_yyz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_yyyz[j] + pa_xyy[j] * pb_yyyz[j]) * s_0_0[j];

                t_xyy_yyzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xyy[j] * fx[j] * fx[j] +

                                pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xyy[j] * pb_yy[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * pb_zz[j] +

                                2.0 * pa_xy[j] * fx[j] * pb_yzz[j] + 0.5 * pa_x[j] * fx[j] * pb_yyzz[j] + pa_xyy[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (29) = (58,60)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pb_yz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, \
                                     t_xyy_yzzz, t_xyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_yzzz[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xyy[j] * pb_yz[j] * fx[j] + pa_xy[j] * fx[j] * pb_zzz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_yzzz[j] + pa_xyy[j] * pb_yzzz[j]) * s_0_0[j];

                t_xyy_zzzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] +

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xyy[j] * pb_zz[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_zzzz[j] +

                                pa_xyy[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (30) = (60,62)

            #pragma omp simd aligned(fx, pa_xyz, pa_xz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxy, \
                                     pb_xy, pb_y, s_0_0, t_xyz_xxxx, t_xyz_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xxxx[j] = (0.75 * pa_xyz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] +

                                3.0 * pa_xyz[j] * pb_xx[j] * fx[j] + 2.0 * fx[j] * pa_yz[j] * pb_xxx[j] + pa_xyz[j] * pb_xxxx[j]) * s_0_0[j];

                t_xyz_xxxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 1.5 * pa_xyz[j] * pb_xy[j] * fx[j] +

                                0.5 * pa_xz[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_yz[j] * pb_xxy[j] + pa_xyz[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (31) = (62,64)

            #pragma omp simd aligned(fx, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxz, pb_xxy, \
                                     pb_xxyy, pb_xxz, pb_xy, pb_xyy, pb_xz, pb_y, pb_yy, pb_z, s_0_0, t_xyz_xxxz, t_xyz_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xxxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 1.5 * pa_xyz[j] * pb_xz[j] * fx[j] +

                                0.5 * pa_xy[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_yz[j] * pb_xxz[j] + pa_xyz[j] * pb_xxxz[j]) * s_0_0[j];

                t_xyz_xxyy[j] = (0.25 * pa_xyz[j] * fx[j] * fx[j] + 0.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                0.5 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xyz[j] * pb_xx[j] * fx[j] +

                                0.5 * pa_xyz[j] * fx[j] * pb_yy[j] + pa_xz[j] * fx[j] * pb_xxy[j] + fx[j] * pa_yz[j] * pb_xyy[j] + pa_xyz[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (32) = (64,66)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxy, pb_xxyz, \
                                     pb_xxz, pb_xxzz, pb_xy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, \
                                     t_xyz_xxyz, t_xyz_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xxyz[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] +

                                0.5 * pa_xyz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxy[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxz[j] +

                                fx[j] * pa_yz[j] * pb_xyz[j] + pa_xyz[j] * pb_xxyz[j]) * s_0_0[j];

                t_xyz_xxzz[j] = (0.25 * pa_xyz[j] * fx[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xyz[j] * pb_xx[j] * fx[j] +

                                0.5 * pa_xyz[j] * fx[j] * pb_zz[j] + pa_xy[j] * fx[j] * pb_xxz[j] + fx[j] * pa_yz[j] * pb_xzz[j] + pa_xyz[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (33) = (66,68)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_z, s_0_0, t_xyz_xyyy, \
                                     t_xyz_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 1.5 * pa_xyz[j] * pb_xy[j] * fx[j] +

                                1.5 * pa_xz[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyy[j] + pa_xyz[j] * pb_xyyy[j]) * s_0_0[j];

                t_xyz_xyyz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] +

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] +

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xyz[j] * pb_xz[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * pb_xyy[j] +

                                pa_xz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyz[j] + pa_xyz[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (34) = (68,70)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xy, pb_xyz, pb_xyzz, \
                                     pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xyz_xyzz, \
                                     t_xyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xyzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] +

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xyz[j] * pb_xy[j] * fx[j] + pa_xy[j] * fx[j] * pb_xyz[j] +

                                0.5 * pa_xz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yzz[j] + pa_xyz[j] * pb_xyzz[j]) * s_0_0[j];

                t_xyz_xzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 1.5 * pa_xyz[j] * pb_xz[j] * fx[j] +

                                1.5 * pa_xy[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yz[j] * pb_zzz[j] + pa_xyz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (35) = (70,72)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_xyz_yyyy, t_xyz_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_yyyy[j] = (0.75 * pa_xyz[j] * fx[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                3.0 * pa_xyz[j] * pb_yy[j] * fx[j] + 2.0 * pa_xz[j] * fx[j] * pb_yyy[j] + pa_xyz[j] * pb_yyyy[j]) * s_0_0[j];

                t_xyz_yyyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xyz[j] * pb_yz[j] * fx[j] +

                                0.5 * pa_xy[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyz[j] + pa_xyz[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (36) = (72,74)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xz, pb_y, pb_yy, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xyz_yyzz, t_xyz_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_yyzz[j] = (0.25 * pa_xyz[j] * fx[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyz[j] * pb_yy[j] * fx[j] +

                                0.5 * pa_xyz[j] * fx[j] * pb_zz[j] + pa_xy[j] * fx[j] * pb_yyz[j] + pa_xz[j] * fx[j] * pb_yzz[j] + pa_xyz[j] * pb_yyzz[j]) * s_0_0[j];

                t_xyz_yzzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xyz[j] * pb_yz[j] * fx[j] +

                                1.5 * pa_xy[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xz[j] * fx[j] * pb_zzz[j] + pa_xyz[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (37) = (74,76)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xzz, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_z, pb_zz, \
                                     pb_zzz, pb_zzzz, s_0_0, t_xyz_zzzz, t_xzz_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_zzzz[j] = (0.75 * pa_xyz[j] * fx[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                3.0 * pa_xyz[j] * pb_zz[j] * fx[j] + 2.0 * pa_xy[j] * fx[j] * pb_zzz[j] + pa_xyz[j] * pb_zzzz[j]) * s_0_0[j];

                t_xzz_xxxx[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xzz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xzz[j] * pb_xx[j] * fx[j] + 2.0 * fx[j] * pa_zz[j] * pb_xxx[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xxxx[j] + pa_xzz[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (38) = (76,78)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxz, pb_xy, pb_xz, pb_y, pb_z, s_0_0, t_xzz_xxxy, t_xzz_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xxxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] + 1.5 * pa_xzz[j] * pb_xy[j] * fx[j] +

                                1.5 * fx[j] * pa_zz[j] * pb_xxy[j] + 0.5 * pa_x[j] * fx[j] * pb_xxxy[j] + pa_xzz[j] * pb_xxxy[j]) * s_0_0[j];

                t_xzz_xxxz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] +

                                1.5 * pa_xzz[j] * pb_xz[j] * fx[j] + pa_xz[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_zz[j] * pb_xxz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xxxz[j] + pa_xzz[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (39) = (78,80)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxy, pb_xxyy, pb_xxyz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_y, pb_yy, pb_yz, s_0_0, t_xzz_xxyy, t_xzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xxyy[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xzz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * fx[j] * pb_xyy[j] +

                                0.5 * pa_xzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] * pb_yy[j] + fx[j] * pa_zz[j] * pb_xyy[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xxyy[j] + pa_xzz[j] * pb_xxyy[j]) * s_0_0[j];

                t_xzz_xxyz[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + fx[j] * fx[j] * pa_z[j] * pb_xy[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xzz[j] * fx[j] * pb_yz[j] +

                                pa_xz[j] * fx[j] * pb_xxy[j] + fx[j] * pa_zz[j] * pb_xyz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxyz[j] +

                                pa_xzz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (40) = (80,82)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyyy, \
                                     pb_xz, pb_xzz, pb_y, pb_yyy, pb_z, pb_zz, s_0_0, t_xzz_xxzz, t_xzz_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xxzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xzz[j] * fx[j] * fx[j] + pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                0.5 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 2.0 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_xzz[j] * pb_xx[j] * fx[j] +

                                0.5 * pa_xzz[j] * fx[j] * pb_zz[j] + 2.0 * pa_xz[j] * fx[j] * pb_xxz[j] + fx[j] * pa_zz[j] * pb_xzz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xxzz[j] + pa_xzz[j] * pb_xxzz[j]) * s_0_0[j];

                t_xzz_xyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xzz[j] * pb_xy[j] * fx[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_yyy[j] + 0.5 * pa_x[j] * fx[j] * pb_xyyy[j] + pa_xzz[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (41) = (82,84)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xy, pb_xyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_y, pb_yy, pb_yyz, pb_yz, pb_yzz, pb_z, s_0_0, t_xzz_xyyz, t_xzz_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xyyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] +

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] +

                                0.5 * pa_xzz[j] * pb_xz[j] * fx[j] + pa_xz[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xyyz[j] + pa_xzz[j] * pb_xyyz[j]) * s_0_0[j];

                t_xzz_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] +

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yzz[j] +

                                0.5 * pa_xzz[j] * pb_xy[j] * fx[j] + 2.0 * pa_xz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yzz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xyzz[j] + pa_xzz[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (42) = (84,86)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xz, pb_xzz, pb_xzzz, pb_yy, pb_yyyy, \
                                     pb_z, pb_zz, pb_zzz, s_0_0, t_xzz_xzzz, t_xzz_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xzzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] +

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] +

                                1.5 * pa_xzz[j] * pb_xz[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_zzz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_xzzz[j] + pa_xzz[j] * pb_xzzz[j]) * s_0_0[j];

                t_xzz_yyyy[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] +

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xzz[j] * pb_yy[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_yyyy[j] +

                                pa_xzz[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (43) = (86,88)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pb_y, pb_yy, pb_yyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, \
                                     pb_z, pb_zz, s_0_0, t_xzz_yyyz, t_xzz_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_yyyz[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xzz[j] * pb_yz[j] * fx[j] + pa_xz[j] * fx[j] * pb_yyy[j] +

                                0.5 * pa_x[j] * fx[j] * pb_yyyz[j] + pa_xzz[j] * pb_yyyz[j]) * s_0_0[j];

                t_xzz_yyzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xzz[j] * fx[j] * fx[j] +

                                pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xzz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] * pb_zz[j] +

                                2.0 * pa_xz[j] * fx[j] * pb_yyz[j] + 0.5 * pa_x[j] * fx[j] * pb_yyzz[j] + pa_xzz[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (44) = (88,90)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pb_y, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, \
                                     pb_zzzz, s_0_0, t_xzz_yzzz, t_xzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_yzzz[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xzz[j] * pb_yz[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * pb_yzz[j] +

                                0.5 * pa_x[j] * fx[j] * pb_yzzz[j] + pa_xzz[j] * pb_yzzz[j]) * s_0_0[j];

                t_xzz_zzzz[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] +

                                6.0 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xzz[j] * pb_zz[j] * fx[j] +

                                4.0 * pa_xz[j] * fx[j] * pb_zzz[j] + 0.5 * pa_x[j] * fx[j] * pb_zzzz[j] + pa_xzz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (45) = (90,92)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xy, s_0_0, \
                                     t_yyy_xxxx, t_yyy_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xxxx[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yyy[j] * fx[j] * fx[j] +

                                4.5 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_yyy[j] * pb_xx[j] * fx[j] + 1.5 * pa_y[j] * fx[j] * pb_xxxx[j] +

                                pa_yyy[j] * pb_xxxx[j]) * s_0_0[j];

                t_yyy_xxxy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_yyy[j] * pb_xy[j] * fx[j] +

                                1.5 * pa_yy[j] * fx[j] * pb_xxx[j] + 1.5 * pa_y[j] * fx[j] * pb_xxxy[j] + pa_yyy[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (46) = (92,94)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pb_xx, pb_xxxz, pb_xxy, pb_xxyy, pb_xz, pb_y, pb_yy, \
                                     s_0_0, t_yyy_xxxz, t_yyy_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xxxz[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yyy[j] * pb_xz[j] * fx[j] +

                                1.5 * pa_y[j] * fx[j] * pb_xxxz[j] + pa_yyy[j] * pb_xxxz[j]) * s_0_0[j];

                t_yyy_xxyy[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_yyy[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] +

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * fx[j] * fx[j] * pb_xxy[j] +

                                0.5 * pa_yyy[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] * pb_yy[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxy[j] +

                                1.5 * pa_y[j] * fx[j] * pb_xxyy[j] + pa_yyy[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (47) = (94,96)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pb_xx, pb_xxyz, pb_xxz, pb_xxzz, pb_yz, pb_z, pb_zz, \
                                     s_0_0, t_yyy_xxyz, t_yyy_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xxyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_yyy[j] * fx[j] * pb_yz[j] +

                                1.5 * pa_yy[j] * fx[j] * pb_xxz[j] + 1.5 * pa_y[j] * fx[j] * pb_xxyz[j] + pa_yyy[j] * pb_xxyz[j]) * s_0_0[j];

                t_yyy_xxzz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_yyy[j] * fx[j] * fx[j] +

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_yyy[j] * pb_xx[j] * fx[j] +

                                0.5 * pa_yyy[j] * fx[j] * pb_zz[j] + 1.5 * pa_y[j] * fx[j] * pb_xxzz[j] + pa_yyy[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (48) = (96,98)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pb_x, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xz, \
                                     s_0_0, t_yyy_xyyy, t_yyy_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xyyy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                6.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 2.25 * fx[j] * fx[j] * pb_xyy[j] + 1.5 * pa_yyy[j] * pb_xy[j] * fx[j] +

                                4.5 * pa_yy[j] * fx[j] * pb_xyy[j] + 1.5 * pa_y[j] * fx[j] * pb_xyyy[j] + pa_yyy[j] * pb_xyyy[j]) * s_0_0[j];

                t_yyy_xyyz[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * fx[j] * fx[j] * pb_xyz[j] +

                                0.5 * pa_yyy[j] * pb_xz[j] * fx[j] + 3.0 * pa_yy[j] * fx[j] * pb_xyz[j] + 1.5 * pa_y[j] * fx[j] * pb_xyyz[j] +

                                pa_yyy[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (49) = (98,100)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pb_x, pb_xy, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, s_0_0, \
                                     t_yyy_xyzz, t_yyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_yyy[j] * pb_xy[j] * fx[j] +

                                1.5 * pa_yy[j] * fx[j] * pb_xzz[j] + 1.5 * pa_y[j] * fx[j] * pb_xyzz[j] + pa_yyy[j] * pb_xyzz[j]) * s_0_0[j];

                t_yyy_xzzz[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yyy[j] * pb_xz[j] * fx[j] +

                                1.5 * pa_y[j] * fx[j] * pb_xzzz[j] + pa_yyy[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (50) = (100,102)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yz, \
                                     pb_z, s_0_0, t_yyy_yyyy, t_yyy_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_yyyy[j] = (5.625 * pa_y[j] * fx[j] * fx[j] * fx[j] + 7.5 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_yyy[j] * fx[j] * fx[j] + 9.0 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] +

                                13.5 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yyy[j] * pb_yy[j] * fx[j] +

                                6.0 * pa_yy[j] * fx[j] * pb_yyy[j] + 1.5 * pa_y[j] * fx[j] * pb_yyyy[j] + pa_yyy[j] * pb_yyyy[j]) * s_0_0[j];

                t_yyy_yyyz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                6.75 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 2.25 * fx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_yyy[j] * pb_yz[j] * fx[j] +

                                4.5 * pa_yy[j] * fx[j] * pb_yyz[j] + 1.5 * pa_y[j] * fx[j] * pb_yyyz[j] + pa_yyy[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (51) = (102,104)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pb_y, pb_yy, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, \
                                     pb_zzz, s_0_0, t_yyy_yyzz, t_yyy_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_yyzz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_yyy[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] +

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * fx[j] * fx[j] * pb_yzz[j] +

                                0.5 * pa_yyy[j] * pb_yy[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] * pb_zz[j] + 3.0 * pa_yy[j] * fx[j] * pb_yzz[j] +

                                1.5 * pa_y[j] * fx[j] * pb_yyzz[j] + pa_yyy[j] * pb_yyzz[j]) * s_0_0[j];

                t_yyy_yzzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_yyy[j] * pb_yz[j] * fx[j] +

                                1.5 * pa_yy[j] * fx[j] * pb_zzz[j] + 1.5 * pa_y[j] * fx[j] * pb_yzzz[j] + pa_yyy[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (52) = (104,106)

            #pragma omp simd aligned(fx, pa_y, pa_yyy, pa_yyz, pa_z, pb_xx, pb_xxxx, pb_zz, pb_zzzz, s_0_0, \
                                     t_yyy_zzzz, t_yyz_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_zzzz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yyy[j] * fx[j] * fx[j] +

                                4.5 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_yyy[j] * pb_zz[j] * fx[j] + 1.5 * pa_y[j] * fx[j] * pb_zzzz[j] +

                                pa_yyy[j] * pb_zzzz[j]) * s_0_0[j];

                t_yyz_xxxx[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_yyz[j] * fx[j] * fx[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 3.0 * pa_yyz[j] * pb_xx[j] * fx[j] + 0.5 * fx[j] * pa_z[j] * pb_xxxx[j] +

                                pa_yyz[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (53) = (106,108)

            #pragma omp simd aligned(fx, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xxx, pb_xxxy, pb_xxxz, pb_xy, pb_xz, \
                                     s_0_0, t_yyz_xxxy, t_yyz_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xxxy[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 1.5 * pa_yyz[j] * pb_xy[j] * fx[j] + pa_yz[j] * fx[j] * pb_xxx[j] +

                                0.5 * fx[j] * pa_z[j] * pb_xxxy[j] + pa_yyz[j] * pb_xxxy[j]) * s_0_0[j];

                t_yyz_xxxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_yyz[j] * pb_xz[j] * fx[j] +

                                0.5 * pa_yy[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pa_z[j] * pb_xxxz[j] + pa_yyz[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (54) = (108,110)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_xx, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, \
                                     pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_yyz_xxyy, t_yyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xxyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * pa_yyz[j] * fx[j] * fx[j] +

                                pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_yyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyz[j] * fx[j] * pb_yy[j] +

                                2.0 * pa_yz[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pa_z[j] * pb_xxyy[j] + pa_yyz[j] * pb_xxyy[j]) * s_0_0[j];

                t_yyz_xxyz[j] = (0.25 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_yyz[j] * fx[j] * pb_yz[j] +

                                0.5 * pa_yy[j] * fx[j] * pb_xxy[j] + pa_yz[j] * fx[j] * pb_xxz[j] + 0.5 * fx[j] * pa_z[j] * pb_xxyz[j] +

                                pa_yyz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (55) = (110,112)

            #pragma omp simd aligned(fx, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_z, pb_zz, s_0_0, t_yyz_xxzz, t_yyz_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xxzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_yyz[j] * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_xxz[j] +

                                0.5 * pa_yyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyz[j] * fx[j] * pb_zz[j] + pa_yy[j] * fx[j] * pb_xxz[j] +

                                0.5 * fx[j] * pa_z[j] * pb_xxzz[j] + pa_yyz[j] * pb_xxzz[j]) * s_0_0[j];

                t_yyz_xyyy[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] +

                                2.25 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 1.5 * pa_yyz[j] * pb_xy[j] * fx[j] + 3.0 * pa_yz[j] * fx[j] * pb_xyy[j] +

                                0.5 * fx[j] * pa_z[j] * pb_xyyy[j] + pa_yyz[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (56) = (112,114)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xy, pb_xyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, s_0_0, t_yyz_xyyz, t_yyz_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xyy[j] +

                                0.5 * pa_yyz[j] * pb_xz[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] * pb_xyy[j] + 2.0 * pa_yz[j] * fx[j] * pb_xyz[j] +

                                0.5 * fx[j] * pa_z[j] * pb_xyyz[j] + pa_yyz[j] * pb_xyyz[j]) * s_0_0[j];

                t_yyz_xyzz[j] = (0.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + pa_y[j] * fx[j] * fx[j] * pb_xz[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_yyz[j] * pb_xy[j] * fx[j] +

                                pa_yy[j] * fx[j] * pb_xyz[j] + pa_yz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_z[j] * pb_xyzz[j] +

                                pa_yyz[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (57) = (114,116)

            #pragma omp simd aligned(fx, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, s_0_0, t_yyz_xzzz, t_yyz_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] + 1.5 * pa_yyz[j] * pb_xz[j] * fx[j] +

                                1.5 * pa_yy[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_z[j] * pb_xzzz[j] + pa_yyz[j] * pb_xzzz[j]) * s_0_0[j];

                t_yyz_yyyy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_yyz[j] * fx[j] * fx[j] +

                                6.0 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 3.0 * pa_yyz[j] * pb_yy[j] * fx[j] +

                                4.0 * pa_yz[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_z[j] * pb_yyyy[j] + pa_yyz[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (58) = (116,118)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, t_yyz_yyyz, t_yyz_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_yyyz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] +

                                2.25 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_yyz[j] * pb_yz[j] * fx[j] +

                                0.5 * pa_yy[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yz[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pa_z[j] * pb_yyyz[j] +

                                pa_yyz[j] * pb_yyyz[j]) * s_0_0[j];

                t_yyz_yyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_yyz[j] * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + pa_yz[j] * fx[j] * fx[j] * pb_y[j] +

                                2.0 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_yyz[j] * pb_yy[j] * fx[j] +

                                0.5 * pa_yyz[j] * fx[j] * pb_zz[j] + pa_yy[j] * fx[j] * pb_yyz[j] + 2.0 * pa_yz[j] * fx[j] * pb_yzz[j] +

                                0.5 * fx[j] * pa_z[j] * pb_yyzz[j] + pa_yyz[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (59) = (118,120)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_y, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, \
                                     pb_zzz, pb_zzzz, s_0_0, t_yyz_yzzz, t_yyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_yzzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_yyz[j] * pb_yz[j] * fx[j] +

                                1.5 * pa_yy[j] * fx[j] * pb_yzz[j] + pa_yz[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_z[j] * pb_yzzz[j] +

                                pa_yyz[j] * pb_yzzz[j]) * s_0_0[j];

                t_yyz_zzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_yyz[j] * fx[j] * fx[j] + 3.0 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] +

                                fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_yyz[j] * pb_zz[j] * fx[j] + 2.0 * pa_yy[j] * fx[j] * pb_zzz[j] +

                                0.5 * fx[j] * pa_z[j] * pb_zzzz[j] + pa_yyz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (60) = (120,122)

            #pragma omp simd aligned(fx, pa_y, pa_yzz, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xy, s_0_0, \
                                     t_yzz_xxxx, t_yzz_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xxxx[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yzz[j] * fx[j] * fx[j] +

                                1.5 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_yzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_xxxx[j] +

                                pa_yzz[j] * pb_xxxx[j]) * s_0_0[j];

                t_yzz_xxxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] +

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_yzz[j] * pb_xy[j] * fx[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_xxx[j] + 0.5 * pa_y[j] * fx[j] * pb_xxxy[j] + pa_yzz[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (61) = (122,124)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xz, pb_y, pb_yy, s_0_0, t_yzz_xxxz, t_yzz_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xxxz[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yzz[j] * pb_xz[j] * fx[j] + pa_yz[j] * fx[j] * pb_xxx[j] +

                                0.5 * pa_y[j] * fx[j] * pb_xxxz[j] + pa_yzz[j] * pb_xxxz[j]) * s_0_0[j];

                t_yzz_xxyy[j] = (0.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_yzz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] +

                                0.25 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * fx[j] * pb_xxy[j] +

                                0.5 * pa_yzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yzz[j] * fx[j] * pb_yy[j] + fx[j] * pa_zz[j] * pb_xxy[j] +

                                0.5 * pa_y[j] * fx[j] * pb_xxyy[j] + pa_yzz[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (62) = (124,126)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_xx, pb_xxy, pb_xxyz, pb_xxz, pb_xxzz, \
                                     pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_yzz_xxyz, t_yzz_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xxyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] +

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.25 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] +

                                0.5 * pa_yzz[j] * fx[j] * pb_yz[j] + pa_yz[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxz[j] +

                                0.5 * pa_y[j] * fx[j] * pb_xxyz[j] + pa_yzz[j] * pb_xxyz[j]) * s_0_0[j];

                t_yzz_xxzz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_yzz[j] * fx[j] * fx[j] +

                                pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] +

                                0.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_yzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yzz[j] * fx[j] * pb_zz[j] +

                                2.0 * pa_yz[j] * fx[j] * pb_xxz[j] + 0.5 * pa_y[j] * fx[j] * pb_xxzz[j] + pa_yzz[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (63) = (126,128)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, \
                                     pb_xyz, pb_xz, s_0_0, t_yzz_xyyy, t_yzz_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] +

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] + 1.5 * pa_yzz[j] * pb_xy[j] * fx[j] +

                                1.5 * fx[j] * pa_zz[j] * pb_xyy[j] + 0.5 * pa_y[j] * fx[j] * pb_xyyy[j] + pa_yzz[j] * pb_xyyy[j]) * s_0_0[j];

                t_yzz_xyyz[j] = (0.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xy[j] +

                                0.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_yzz[j] * pb_xz[j] * fx[j] +

                                pa_yz[j] * fx[j] * pb_xyy[j] + fx[j] * pa_zz[j] * pb_xyz[j] + 0.5 * pa_y[j] * fx[j] * pb_xyyz[j] +

                                pa_yzz[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (64) = (128,130)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xy, pb_xyz, pb_xyzz, pb_xz, pb_xzz, \
                                     pb_xzzz, s_0_0, t_yzz_xyzz, t_yzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] +

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xzz[j] +

                                0.5 * pa_yzz[j] * pb_xy[j] * fx[j] + 2.0 * pa_yz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xzz[j] +

                                0.5 * pa_y[j] * fx[j] * pb_xyzz[j] + pa_yzz[j] * pb_xyzz[j]) * s_0_0[j];

                t_yzz_xzzz[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] +

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yzz[j] * pb_xz[j] * fx[j] + 3.0 * pa_yz[j] * fx[j] * pb_xzz[j] +

                                0.5 * pa_y[j] * fx[j] * pb_xzzz[j] + pa_yzz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (65) = (130,132)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yz, pb_z, s_0_0, t_yzz_yyyy, t_yzz_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_yyyy[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_yzz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] +

                                fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yzz[j] * pb_yy[j] * fx[j] + 2.0 * fx[j] * pa_zz[j] * pb_yyy[j] +

                                0.5 * pa_y[j] * fx[j] * pb_yyyy[j] + pa_yzz[j] * pb_yyyy[j]) * s_0_0[j];

                t_yzz_yyyz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] +

                                1.5 * pa_yzz[j] * pb_yz[j] * fx[j] + pa_yz[j] * fx[j] * pb_yyy[j] + 1.5 * fx[j] * pa_zz[j] * pb_yyz[j] +

                                0.5 * pa_y[j] * fx[j] * pb_yyyz[j] + pa_yzz[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (66) = (132,134)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, s_0_0, t_yzz_yyzz, t_yzz_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_yyzz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_yzz[j] * fx[j] * fx[j] + pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] +

                                0.5 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 2.0 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] +

                                0.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_yzz[j] * pb_yy[j] * fx[j] +

                                0.5 * pa_yzz[j] * fx[j] * pb_zz[j] + 2.0 * pa_yz[j] * fx[j] * pb_yyz[j] + fx[j] * pa_zz[j] * pb_yzz[j] +

                                0.5 * pa_y[j] * fx[j] * pb_yyzz[j] + pa_yzz[j] * pb_yyzz[j]) * s_0_0[j];

                t_yzz_yzzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] +

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] +

                                1.5 * pa_yzz[j] * pb_yz[j] * fx[j] + 3.0 * pa_yz[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_zzz[j] +

                                0.5 * pa_y[j] * fx[j] * pb_yzzz[j] + pa_yzz[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (67) = (134,136)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zzz, pb_xx, pb_xxxx, pb_z, pb_zz, pb_zzz, \
                                     pb_zzzz, s_0_0, t_yzz_zzzz, t_zzz_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_zzzz[j] = (1.875 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yzz[j] * fx[j] * fx[j] +

                                6.0 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_yzz[j] * pb_zz[j] * fx[j] +

                                4.0 * pa_yz[j] * fx[j] * pb_zzz[j] + 0.5 * pa_y[j] * fx[j] * pb_zzzz[j] + pa_yzz[j] * pb_zzzz[j]) * s_0_0[j];

                t_zzz_xxxx[j] = (1.125 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_zzz[j] * fx[j] * fx[j] +

                                4.5 * pa_z[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_zzz[j] * pb_xx[j] * fx[j] + 1.5 * pa_z[j] * fx[j] * pb_xxxx[j] +

                                pa_zzz[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (68) = (136,138)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_x, pb_xxx, pb_xxxy, pb_xxxz, pb_xy, pb_xz, s_0_0, \
                                     t_zzz_xxxy, t_zzz_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xxxy[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_zzz[j] * pb_xy[j] * fx[j] +

                                1.5 * pa_z[j] * fx[j] * pb_xxxy[j] + pa_zzz[j] * pb_xxxy[j]) * s_0_0[j];

                t_zzz_xxxz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] +

                                2.25 * pa_z[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_zzz[j] * pb_xz[j] * fx[j] +

                                1.5 * pa_zz[j] * fx[j] * pb_xxx[j] + 1.5 * pa_z[j] * fx[j] * pb_xxxz[j] + pa_zzz[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (69) = (138,140)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_xx, pb_xxy, pb_xxyy, pb_xxyz, pb_y, pb_yy, pb_yz, \
                                     s_0_0, t_zzz_xxyy, t_zzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xxyy[j] = (0.375 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_zzz[j] * fx[j] * fx[j] +

                                0.75 * pa_z[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_z[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_zzz[j] * pb_xx[j] * fx[j] +

                                0.5 * pa_zzz[j] * fx[j] * pb_yy[j] + 1.5 * pa_z[j] * fx[j] * pb_xxyy[j] + pa_zzz[j] * pb_xxyy[j]) * s_0_0[j];

                t_zzz_xxyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_z[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_zzz[j] * fx[j] * pb_yz[j] +

                                1.5 * pa_zz[j] * fx[j] * pb_xxy[j] + 1.5 * pa_z[j] * fx[j] * pb_xxyz[j] + pa_zzz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (70) = (140,142)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyyy, pb_z, pb_zz, \
                                     s_0_0, t_zzz_xxzz, t_zzz_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xxzz[j] = (1.125 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_zzz[j] * fx[j] * fx[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] +

                                2.25 * pa_z[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_z[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * fx[j] * fx[j] * pb_xxz[j] +

                                0.5 * pa_zzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_zzz[j] * fx[j] * pb_zz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxz[j] +

                                1.5 * pa_z[j] * fx[j] * pb_xxzz[j] + pa_zzz[j] * pb_xxzz[j]) * s_0_0[j];

                t_zzz_xyyy[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_zzz[j] * pb_xy[j] * fx[j] +

                                1.5 * pa_z[j] * fx[j] * pb_xyyy[j] + pa_zzz[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (71) = (142,144)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_x, pb_xy, pb_xyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, \
                                     s_0_0, t_zzz_xyyz, t_zzz_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_z[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_zzz[j] * pb_xz[j] * fx[j] +

                                1.5 * pa_zz[j] * fx[j] * pb_xyy[j] + 1.5 * pa_z[j] * fx[j] * pb_xyyz[j] + pa_zzz[j] * pb_xyyz[j]) * s_0_0[j];

                t_zzz_xyzz[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * fx[j] * fx[j] * pb_xyz[j] +

                                0.5 * pa_zzz[j] * pb_xy[j] * fx[j] + 3.0 * pa_zz[j] * fx[j] * pb_xyz[j] + 1.5 * pa_z[j] * fx[j] * pb_xyzz[j] +

                                pa_zzz[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (72) = (144,146)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_x, pb_xz, pb_xzz, pb_xzzz, pb_yy, pb_yyyy, s_0_0, \
                                     t_zzz_xzzz, t_zzz_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xzzz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] +

                                6.75 * pa_z[j] * fx[j] * fx[j] * pb_xz[j] + 2.25 * fx[j] * fx[j] * pb_xzz[j] + 1.5 * pa_zzz[j] * pb_xz[j] * fx[j] +

                                4.5 * pa_zz[j] * fx[j] * pb_xzz[j] + 1.5 * pa_z[j] * fx[j] * pb_xzzz[j] + pa_zzz[j] * pb_xzzz[j]) * s_0_0[j];

                t_zzz_yyyy[j] = (1.125 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_zzz[j] * fx[j] * fx[j] +

                                4.5 * pa_z[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_zzz[j] * pb_yy[j] * fx[j] + 1.5 * pa_z[j] * fx[j] * pb_yyyy[j] +

                                pa_zzz[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (73) = (146,148)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_y, pb_yy, pb_yyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, \
                                     pb_z, pb_zz, s_0_0, t_zzz_yyyz, t_zzz_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_yyyz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] +

                                2.25 * pa_z[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_zzz[j] * pb_yz[j] * fx[j] +

                                1.5 * pa_zz[j] * fx[j] * pb_yyy[j] + 1.5 * pa_z[j] * fx[j] * pb_yyyz[j] + pa_zzz[j] * pb_yyyz[j]) * s_0_0[j];

                t_zzz_yyzz[j] = (1.125 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_zzz[j] * fx[j] * fx[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] +

                                2.25 * pa_z[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_z[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * fx[j] * fx[j] * pb_yyz[j] +

                                0.5 * pa_zzz[j] * pb_yy[j] * fx[j] + 0.5 * pa_zzz[j] * fx[j] * pb_zz[j] + 3.0 * pa_zz[j] * fx[j] * pb_yyz[j] +

                                1.5 * pa_z[j] * fx[j] * pb_yyzz[j] + pa_zzz[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (74) = (148,150)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pb_y, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, \
                                     pb_zzzz, s_0_0, t_zzz_yzzz, t_zzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_yzzz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] +

                                6.75 * pa_z[j] * fx[j] * fx[j] * pb_yz[j] + 2.25 * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_zzz[j] * pb_yz[j] * fx[j] +

                                4.5 * pa_zz[j] * fx[j] * pb_yzz[j] + 1.5 * pa_z[j] * fx[j] * pb_yzzz[j] + pa_zzz[j] * pb_yzzz[j]) * s_0_0[j];

                t_zzz_zzzz[j] = (5.625 * pa_z[j] * fx[j] * fx[j] * fx[j] + 7.5 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_zzz[j] * fx[j] * fx[j] + 9.0 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] +

                                13.5 * pa_z[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_zzz[j] * pb_zz[j] * fx[j] +

                                6.0 * pa_zz[j] * fx[j] * pb_zzz[j] + 1.5 * pa_z[j] * fx[j] * pb_zzzz[j] + pa_zzz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForGS(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxxx_0 = primBuffer.data(15 * idx);

            auto t_xxxy_0 = primBuffer.data(15 * idx + 1);

            auto t_xxxz_0 = primBuffer.data(15 * idx + 2);

            auto t_xxyy_0 = primBuffer.data(15 * idx + 3);

            auto t_xxyz_0 = primBuffer.data(15 * idx + 4);

            auto t_xxzz_0 = primBuffer.data(15 * idx + 5);

            auto t_xyyy_0 = primBuffer.data(15 * idx + 6);

            auto t_xyyz_0 = primBuffer.data(15 * idx + 7);

            auto t_xyzz_0 = primBuffer.data(15 * idx + 8);

            auto t_xzzz_0 = primBuffer.data(15 * idx + 9);

            auto t_yyyy_0 = primBuffer.data(15 * idx + 10);

            auto t_yyyz_0 = primBuffer.data(15 * idx + 11);

            auto t_yyzz_0 = primBuffer.data(15 * idx + 12);

            auto t_yzzz_0 = primBuffer.data(15 * idx + 13);

            auto t_zzzz_0 = primBuffer.data(15 * idx + 14);

            // Batch of Integrals (0) = (0,8)

            #pragma omp simd aligned(fx, pa_xx, pa_xxxx, pa_xxxy, pa_xxxz, pa_xxyy, pa_xxyz, pa_xxzz, pa_xy, \
                                     pa_xyyy, pa_xyyz, pa_xz, pa_yy, pa_yz, pa_zz, s_0_0, t_xxxx_0, t_xxxy_0, t_xxxz_0, \
                                     t_xxyy_0, t_xxyz_0, t_xxzz_0, t_xyyy_0, t_xyyz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_0[j] = (0.75 * fx[j] * fx[j] + 3.0 * pa_xx[j] * fx[j] + pa_xxxx[j]) * s_0_0[j];

                t_xxxy_0[j] = (1.5 * pa_xy[j] * fx[j] + pa_xxxy[j]) * s_0_0[j];

                t_xxxz_0[j] = (1.5 * pa_xz[j] * fx[j] + pa_xxxz[j]) * s_0_0[j];

                t_xxyy_0[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] + 0.5 * fx[j] * pa_yy[j] + pa_xxyy[j]) * s_0_0[j];

                t_xxyz_0[j] = (0.5 * fx[j] * pa_yz[j] + pa_xxyz[j]) * s_0_0[j];

                t_xxzz_0[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] + 0.5 * fx[j] * pa_zz[j] + pa_xxzz[j]) * s_0_0[j];

                t_xyyy_0[j] = (1.5 * pa_xy[j] * fx[j] + pa_xyyy[j]) * s_0_0[j];

                t_xyyz_0[j] = (0.5 * pa_xz[j] * fx[j] + pa_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (8,15)

            #pragma omp simd aligned(fx, pa_xy, pa_xyzz, pa_xz, pa_xzzz, pa_yy, pa_yyyy, pa_yyyz, pa_yyzz, pa_yz, \
                                     pa_yzzz, pa_zz, pa_zzzz, s_0_0, t_xyzz_0, t_xzzz_0, t_yyyy_0, t_yyyz_0, t_yyzz_0, \
                                     t_yzzz_0, t_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_0[j] = (0.5 * pa_xy[j] * fx[j] + pa_xyzz[j]) * s_0_0[j];

                t_xzzz_0[j] = (1.5 * pa_xz[j] * fx[j] + pa_xzzz[j]) * s_0_0[j];

                t_yyyy_0[j] = (0.75 * fx[j] * fx[j] + 3.0 * pa_yy[j] * fx[j] + pa_yyyy[j]) * s_0_0[j];

                t_yyyz_0[j] = (1.5 * pa_yz[j] * fx[j] + pa_yyyz[j]) * s_0_0[j];

                t_yyzz_0[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] + 0.5 * fx[j] * pa_zz[j] + pa_yyzz[j]) * s_0_0[j];

                t_yzzz_0[j] = (1.5 * pa_yz[j] * fx[j] + pa_yzzz[j]) * s_0_0[j];

                t_zzzz_0[j] = (0.75 * fx[j] * fx[j] + 3.0 * pa_zz[j] * fx[j] + pa_zzzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForGP(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxxx_x = primBuffer.data(45 * idx);

            auto t_xxxx_y = primBuffer.data(45 * idx + 1);

            auto t_xxxx_z = primBuffer.data(45 * idx + 2);

            auto t_xxxy_x = primBuffer.data(45 * idx + 3);

            auto t_xxxy_y = primBuffer.data(45 * idx + 4);

            auto t_xxxy_z = primBuffer.data(45 * idx + 5);

            auto t_xxxz_x = primBuffer.data(45 * idx + 6);

            auto t_xxxz_y = primBuffer.data(45 * idx + 7);

            auto t_xxxz_z = primBuffer.data(45 * idx + 8);

            auto t_xxyy_x = primBuffer.data(45 * idx + 9);

            auto t_xxyy_y = primBuffer.data(45 * idx + 10);

            auto t_xxyy_z = primBuffer.data(45 * idx + 11);

            auto t_xxyz_x = primBuffer.data(45 * idx + 12);

            auto t_xxyz_y = primBuffer.data(45 * idx + 13);

            auto t_xxyz_z = primBuffer.data(45 * idx + 14);

            auto t_xxzz_x = primBuffer.data(45 * idx + 15);

            auto t_xxzz_y = primBuffer.data(45 * idx + 16);

            auto t_xxzz_z = primBuffer.data(45 * idx + 17);

            auto t_xyyy_x = primBuffer.data(45 * idx + 18);

            auto t_xyyy_y = primBuffer.data(45 * idx + 19);

            auto t_xyyy_z = primBuffer.data(45 * idx + 20);

            auto t_xyyz_x = primBuffer.data(45 * idx + 21);

            auto t_xyyz_y = primBuffer.data(45 * idx + 22);

            auto t_xyyz_z = primBuffer.data(45 * idx + 23);

            auto t_xyzz_x = primBuffer.data(45 * idx + 24);

            auto t_xyzz_y = primBuffer.data(45 * idx + 25);

            auto t_xyzz_z = primBuffer.data(45 * idx + 26);

            auto t_xzzz_x = primBuffer.data(45 * idx + 27);

            auto t_xzzz_y = primBuffer.data(45 * idx + 28);

            auto t_xzzz_z = primBuffer.data(45 * idx + 29);

            auto t_yyyy_x = primBuffer.data(45 * idx + 30);

            auto t_yyyy_y = primBuffer.data(45 * idx + 31);

            auto t_yyyy_z = primBuffer.data(45 * idx + 32);

            auto t_yyyz_x = primBuffer.data(45 * idx + 33);

            auto t_yyyz_y = primBuffer.data(45 * idx + 34);

            auto t_yyyz_z = primBuffer.data(45 * idx + 35);

            auto t_yyzz_x = primBuffer.data(45 * idx + 36);

            auto t_yyzz_y = primBuffer.data(45 * idx + 37);

            auto t_yyzz_z = primBuffer.data(45 * idx + 38);

            auto t_yzzz_x = primBuffer.data(45 * idx + 39);

            auto t_yzzz_y = primBuffer.data(45 * idx + 40);

            auto t_yzzz_z = primBuffer.data(45 * idx + 41);

            auto t_zzzz_x = primBuffer.data(45 * idx + 42);

            auto t_zzzz_y = primBuffer.data(45 * idx + 43);

            auto t_zzzz_z = primBuffer.data(45 * idx + 44);

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_y, pb_z, \
                                     s_0_0, t_xxxx_x, t_xxxx_y, t_xxxx_z, t_xxxy_x, t_xxxy_y: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_x[j] = (3.0 * pa_x[j] * fx[j] * fx[j] + 2.0 * pa_xxx[j] * fx[j] +

                              0.75 * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_xx[j] * fx[j] * pb_x[j] + pa_xxxx[j] * pb_x[j]) * s_0_0[j];

                t_xxxx_y[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_xx[j] * fx[j] * pb_y[j] +

                              pa_xxxx[j] * pb_y[j]) * s_0_0[j];

                t_xxxx_z[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_xx[j] * fx[j] * pb_z[j] +

                              pa_xxxx[j] * pb_z[j]) * s_0_0[j];

                t_xxxy_x[j] = (0.75 * fx[j] * fx[j] * pa_y[j] + 1.5 * pa_xxy[j] * fx[j] +

                              1.5 * pa_xy[j] * fx[j] * pb_x[j] + pa_xxxy[j] * pb_x[j]) * s_0_0[j];

                t_xxxy_y[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] +

                              1.5 * pa_xy[j] * fx[j] * pb_y[j] + pa_xxxy[j] * pb_y[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxxz, pa_xxyy, pa_xxz, pa_xy, pa_xyy, pa_xz, \
                                     pa_yy, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xxxy_z, t_xxxz_x, t_xxxz_y, t_xxxz_z, \
                                     t_xxyy_x: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_z[j] = (1.5 * pa_xy[j] * fx[j] * pb_z[j] + pa_xxxy[j] * pb_z[j]) * s_0_0[j];

                t_xxxz_x[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 1.5 * pa_xxz[j] * fx[j] +

                              1.5 * pa_xz[j] * fx[j] * pb_x[j] + pa_xxxz[j] * pb_x[j]) * s_0_0[j];

                t_xxxz_y[j] = (1.5 * pa_xz[j] * fx[j] * pb_y[j] + pa_xxxz[j] * pb_y[j]) * s_0_0[j];

                t_xxxz_z[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] +

                              1.5 * pa_xz[j] * fx[j] * pb_z[j] + pa_xxxz[j] * pb_z[j]) * s_0_0[j];

                t_xxyy_x[j] = (0.5 * pa_x[j] * fx[j] * fx[j] + pa_xyy[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xx[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_yy[j] * pb_x[j] + pa_xxyy[j] * pb_x[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyy, pa_xxyz, pa_xxz, pa_xyz, pa_y, pa_yy, pa_yz, pa_z, \
                                     pb_x, pb_y, pb_z, s_0_0, t_xxyy_y, t_xxyy_z, t_xxyz_x, t_xxyz_y, t_xxyz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_y[j] = (0.5 * fx[j] * fx[j] * pa_y[j] + pa_xxy[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_xx[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_yy[j] * pb_y[j] + pa_xxyy[j] * pb_y[j]) * s_0_0[j];

                t_xxyy_z[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xx[j] * fx[j] * pb_z[j] +

                              0.5 * fx[j] * pa_yy[j] * pb_z[j] + pa_xxyy[j] * pb_z[j]) * s_0_0[j];

                t_xxyz_x[j] = (pa_xyz[j] * fx[j] + 0.5 * fx[j] * pa_yz[j] * pb_x[j] + pa_xxyz[j] * pb_x[j]) * s_0_0[j];

                t_xxyz_y[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_xxz[j] * fx[j] +

                              0.5 * fx[j] * pa_yz[j] * pb_y[j] + pa_xxyz[j] * pb_y[j]) * s_0_0[j];

                t_xxyz_z[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xxy[j] * fx[j] +

                              0.5 * fx[j] * pa_yz[j] * pb_z[j] + pa_xxyz[j] * pb_z[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xy, pa_xyy, pa_xyyy, pa_xzz, pa_y, pa_yyy, \
                                     pa_z, pa_zz, pb_x, pb_y, pb_z, s_0_0, t_xxzz_x, t_xxzz_y, t_xxzz_z, t_xyyy_x, \
                                     t_xyyy_y: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_x[j] = (0.5 * pa_x[j] * fx[j] * fx[j] + pa_xzz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_x[j] +

                              0.5 * pa_xx[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_zz[j] * pb_x[j] + pa_xxzz[j] * pb_x[j]) * s_0_0[j];

                t_xxzz_y[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xx[j] * fx[j] * pb_y[j] +

                              0.5 * fx[j] * pa_zz[j] * pb_y[j] + pa_xxzz[j] * pb_y[j]) * s_0_0[j];

                t_xxzz_z[j] = (0.5 * fx[j] * fx[j] * pa_z[j] + pa_xxz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_xx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_zz[j] * pb_z[j] + pa_xxzz[j] * pb_z[j]) * s_0_0[j];

                t_xyyy_x[j] = (0.75 * fx[j] * fx[j] * pa_y[j] + 0.5 * fx[j] * pa_yyy[j] +

                              1.5 * pa_xy[j] * fx[j] * pb_x[j] + pa_xyyy[j] * pb_x[j]) * s_0_0[j];

                t_xyyy_y[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 1.5 * pa_xyy[j] * fx[j] +

                              1.5 * pa_xy[j] * fx[j] * pb_y[j] + pa_xyyy[j] * pb_y[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_xyyz, pa_xyz, pa_xyzz, pa_xz, pa_y, pa_yyz, \
                                     pa_yzz, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xyyy_z, t_xyyz_x, t_xyyz_y, t_xyyz_z, \
                                     t_xyzz_x: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_z[j] = (1.5 * pa_xy[j] * fx[j] * pb_z[j] + pa_xyyy[j] * pb_z[j]) * s_0_0[j];

                t_xyyz_x[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * pa_yyz[j] +

                              0.5 * pa_xz[j] * fx[j] * pb_x[j] + pa_xyyz[j] * pb_x[j]) * s_0_0[j];

                t_xyyz_y[j] = (pa_xyz[j] * fx[j] + 0.5 * pa_xz[j] * fx[j] * pb_y[j] + pa_xyyz[j] * pb_y[j]) * s_0_0[j];

                t_xyyz_z[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] +

                              0.5 * pa_xz[j] * fx[j] * pb_z[j] + pa_xyyz[j] * pb_z[j]) * s_0_0[j];

                t_xyzz_x[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * fx[j] * pa_yzz[j] +

                              0.5 * pa_xy[j] * fx[j] * pb_x[j] + pa_xyzz[j] * pb_x[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zzz, pb_x, pb_y, \
                                     pb_z, s_0_0, t_xyzz_y, t_xyzz_z, t_xzzz_x, t_xzzz_y, t_xzzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_y[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] +

                              0.5 * pa_xy[j] * fx[j] * pb_y[j] + pa_xyzz[j] * pb_y[j]) * s_0_0[j];

                t_xyzz_z[j] = (pa_xyz[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * pb_z[j] + pa_xyzz[j] * pb_z[j]) * s_0_0[j];

                t_xzzz_x[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * pa_zzz[j] +

                              1.5 * pa_xz[j] * fx[j] * pb_x[j] + pa_xzzz[j] * pb_x[j]) * s_0_0[j];

                t_xzzz_y[j] = (1.5 * pa_xz[j] * fx[j] * pb_y[j] + pa_xzzz[j] * pb_y[j]) * s_0_0[j];

                t_xzzz_z[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 1.5 * pa_xzz[j] * fx[j] +

                              1.5 * pa_xz[j] * fx[j] * pb_z[j] + pa_xzzz[j] * pb_z[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_y, pb_z, \
                                     s_0_0, t_yyyy_x, t_yyyy_y, t_yyyy_z, t_yyyz_x, t_yyyz_y: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_x[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_yy[j] * fx[j] * pb_x[j] +

                              pa_yyyy[j] * pb_x[j]) * s_0_0[j];

                t_yyyy_y[j] = (3.0 * pa_y[j] * fx[j] * fx[j] + 2.0 * pa_yyy[j] * fx[j] +

                              0.75 * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_yy[j] * fx[j] * pb_y[j] + pa_yyyy[j] * pb_y[j]) * s_0_0[j];

                t_yyyy_z[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_yy[j] * fx[j] * pb_z[j] +

                              pa_yyyy[j] * pb_z[j]) * s_0_0[j];

                t_yyyz_x[j] = (1.5 * pa_yz[j] * fx[j] * pb_x[j] + pa_yyyz[j] * pb_x[j]) * s_0_0[j];

                t_yyyz_y[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 1.5 * pa_yyz[j] * fx[j] +

                              1.5 * pa_yz[j] * fx[j] * pb_y[j] + pa_yyyz[j] * pb_y[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_yzzz, pa_z, \
                                     pa_zz, pb_x, pb_y, pb_z, s_0_0, t_yyyz_z, t_yyzz_x, t_yyzz_y, t_yyzz_z, t_yzzz_x: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_z[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] +

                              1.5 * pa_yz[j] * fx[j] * pb_z[j] + pa_yyyz[j] * pb_z[j]) * s_0_0[j];

                t_yyzz_x[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_yy[j] * fx[j] * pb_x[j] +

                              0.5 * fx[j] * pa_zz[j] * pb_x[j] + pa_yyzz[j] * pb_x[j]) * s_0_0[j];

                t_yyzz_y[j] = (0.5 * pa_y[j] * fx[j] * fx[j] + pa_yzz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_y[j] +

                              0.5 * pa_yy[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_zz[j] * pb_y[j] + pa_yyzz[j] * pb_y[j]) * s_0_0[j];

                t_yyzz_z[j] = (0.5 * fx[j] * fx[j] * pa_z[j] + pa_yyz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_z[j] +

                              0.5 * pa_yy[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_zz[j] * pb_z[j] + pa_yyzz[j] * pb_z[j]) * s_0_0[j];

                t_yzzz_x[j] = (1.5 * pa_yz[j] * fx[j] * pb_x[j] + pa_yzzz[j] * pb_x[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_y, pb_z, \
                                     s_0_0, t_yzzz_y, t_yzzz_z, t_zzzz_x, t_zzzz_y, t_zzzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_y[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * pa_zzz[j] +

                              1.5 * pa_yz[j] * fx[j] * pb_y[j] + pa_yzzz[j] * pb_y[j]) * s_0_0[j];

                t_yzzz_z[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 1.5 * pa_yzz[j] * fx[j] +

                              1.5 * pa_yz[j] * fx[j] * pb_z[j] + pa_yzzz[j] * pb_z[j]) * s_0_0[j];

                t_zzzz_x[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_zz[j] * fx[j] * pb_x[j] +

                              pa_zzzz[j] * pb_x[j]) * s_0_0[j];

                t_zzzz_y[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_zz[j] * fx[j] * pb_y[j] +

                              pa_zzzz[j] * pb_y[j]) * s_0_0[j];

                t_zzzz_z[j] = (3.0 * pa_z[j] * fx[j] * fx[j] + 2.0 * pa_zzz[j] * fx[j] +

                              0.75 * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_zz[j] * fx[j] * pb_z[j] + pa_zzzz[j] * pb_z[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForGD(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(9 * idx);

            auto pb_y = pbDistances.data(9 * idx + 1);

            auto pb_z = pbDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(9 * idx + 3);

            auto pb_xy = pbDistances.data(9 * idx + 4);

            auto pb_xz = pbDistances.data(9 * idx + 5);

            auto pb_yy = pbDistances.data(9 * idx + 6);

            auto pb_yz = pbDistances.data(9 * idx + 7);

            auto pb_zz = pbDistances.data(9 * idx + 8);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxxx_xx = primBuffer.data(90 * idx);

            auto t_xxxx_xy = primBuffer.data(90 * idx + 1);

            auto t_xxxx_xz = primBuffer.data(90 * idx + 2);

            auto t_xxxx_yy = primBuffer.data(90 * idx + 3);

            auto t_xxxx_yz = primBuffer.data(90 * idx + 4);

            auto t_xxxx_zz = primBuffer.data(90 * idx + 5);

            auto t_xxxy_xx = primBuffer.data(90 * idx + 6);

            auto t_xxxy_xy = primBuffer.data(90 * idx + 7);

            auto t_xxxy_xz = primBuffer.data(90 * idx + 8);

            auto t_xxxy_yy = primBuffer.data(90 * idx + 9);

            auto t_xxxy_yz = primBuffer.data(90 * idx + 10);

            auto t_xxxy_zz = primBuffer.data(90 * idx + 11);

            auto t_xxxz_xx = primBuffer.data(90 * idx + 12);

            auto t_xxxz_xy = primBuffer.data(90 * idx + 13);

            auto t_xxxz_xz = primBuffer.data(90 * idx + 14);

            auto t_xxxz_yy = primBuffer.data(90 * idx + 15);

            auto t_xxxz_yz = primBuffer.data(90 * idx + 16);

            auto t_xxxz_zz = primBuffer.data(90 * idx + 17);

            auto t_xxyy_xx = primBuffer.data(90 * idx + 18);

            auto t_xxyy_xy = primBuffer.data(90 * idx + 19);

            auto t_xxyy_xz = primBuffer.data(90 * idx + 20);

            auto t_xxyy_yy = primBuffer.data(90 * idx + 21);

            auto t_xxyy_yz = primBuffer.data(90 * idx + 22);

            auto t_xxyy_zz = primBuffer.data(90 * idx + 23);

            auto t_xxyz_xx = primBuffer.data(90 * idx + 24);

            auto t_xxyz_xy = primBuffer.data(90 * idx + 25);

            auto t_xxyz_xz = primBuffer.data(90 * idx + 26);

            auto t_xxyz_yy = primBuffer.data(90 * idx + 27);

            auto t_xxyz_yz = primBuffer.data(90 * idx + 28);

            auto t_xxyz_zz = primBuffer.data(90 * idx + 29);

            auto t_xxzz_xx = primBuffer.data(90 * idx + 30);

            auto t_xxzz_xy = primBuffer.data(90 * idx + 31);

            auto t_xxzz_xz = primBuffer.data(90 * idx + 32);

            auto t_xxzz_yy = primBuffer.data(90 * idx + 33);

            auto t_xxzz_yz = primBuffer.data(90 * idx + 34);

            auto t_xxzz_zz = primBuffer.data(90 * idx + 35);

            auto t_xyyy_xx = primBuffer.data(90 * idx + 36);

            auto t_xyyy_xy = primBuffer.data(90 * idx + 37);

            auto t_xyyy_xz = primBuffer.data(90 * idx + 38);

            auto t_xyyy_yy = primBuffer.data(90 * idx + 39);

            auto t_xyyy_yz = primBuffer.data(90 * idx + 40);

            auto t_xyyy_zz = primBuffer.data(90 * idx + 41);

            auto t_xyyz_xx = primBuffer.data(90 * idx + 42);

            auto t_xyyz_xy = primBuffer.data(90 * idx + 43);

            auto t_xyyz_xz = primBuffer.data(90 * idx + 44);

            auto t_xyyz_yy = primBuffer.data(90 * idx + 45);

            auto t_xyyz_yz = primBuffer.data(90 * idx + 46);

            auto t_xyyz_zz = primBuffer.data(90 * idx + 47);

            auto t_xyzz_xx = primBuffer.data(90 * idx + 48);

            auto t_xyzz_xy = primBuffer.data(90 * idx + 49);

            auto t_xyzz_xz = primBuffer.data(90 * idx + 50);

            auto t_xyzz_yy = primBuffer.data(90 * idx + 51);

            auto t_xyzz_yz = primBuffer.data(90 * idx + 52);

            auto t_xyzz_zz = primBuffer.data(90 * idx + 53);

            auto t_xzzz_xx = primBuffer.data(90 * idx + 54);

            auto t_xzzz_xy = primBuffer.data(90 * idx + 55);

            auto t_xzzz_xz = primBuffer.data(90 * idx + 56);

            auto t_xzzz_yy = primBuffer.data(90 * idx + 57);

            auto t_xzzz_yz = primBuffer.data(90 * idx + 58);

            auto t_xzzz_zz = primBuffer.data(90 * idx + 59);

            auto t_yyyy_xx = primBuffer.data(90 * idx + 60);

            auto t_yyyy_xy = primBuffer.data(90 * idx + 61);

            auto t_yyyy_xz = primBuffer.data(90 * idx + 62);

            auto t_yyyy_yy = primBuffer.data(90 * idx + 63);

            auto t_yyyy_yz = primBuffer.data(90 * idx + 64);

            auto t_yyyy_zz = primBuffer.data(90 * idx + 65);

            auto t_yyyz_xx = primBuffer.data(90 * idx + 66);

            auto t_yyyz_xy = primBuffer.data(90 * idx + 67);

            auto t_yyyz_xz = primBuffer.data(90 * idx + 68);

            auto t_yyyz_yy = primBuffer.data(90 * idx + 69);

            auto t_yyyz_yz = primBuffer.data(90 * idx + 70);

            auto t_yyyz_zz = primBuffer.data(90 * idx + 71);

            auto t_yyzz_xx = primBuffer.data(90 * idx + 72);

            auto t_yyzz_xy = primBuffer.data(90 * idx + 73);

            auto t_yyzz_xz = primBuffer.data(90 * idx + 74);

            auto t_yyzz_yy = primBuffer.data(90 * idx + 75);

            auto t_yyzz_yz = primBuffer.data(90 * idx + 76);

            auto t_yyzz_zz = primBuffer.data(90 * idx + 77);

            auto t_yzzz_xx = primBuffer.data(90 * idx + 78);

            auto t_yzzz_xy = primBuffer.data(90 * idx + 79);

            auto t_yzzz_xz = primBuffer.data(90 * idx + 80);

            auto t_yzzz_yy = primBuffer.data(90 * idx + 81);

            auto t_yzzz_yz = primBuffer.data(90 * idx + 82);

            auto t_yzzz_zz = primBuffer.data(90 * idx + 83);

            auto t_zzzz_xx = primBuffer.data(90 * idx + 84);

            auto t_zzzz_xy = primBuffer.data(90 * idx + 85);

            auto t_zzzz_xz = primBuffer.data(90 * idx + 86);

            auto t_zzzz_yy = primBuffer.data(90 * idx + 87);

            auto t_zzzz_yz = primBuffer.data(90 * idx + 88);

            auto t_zzzz_zz = primBuffer.data(90 * idx + 89);

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, s_0_0, t_xxxx_xx, t_xxxx_xy, t_xxxx_xz, t_xxxx_yy, t_xxxx_yz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xx[j] = (1.875 * fx[j] * fx[j] * fx[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] +

                               6.0 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xxxx[j] * fx[j] + 4.0 * pa_xxx[j] * fx[j] * pb_x[j] +

                               0.75 * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_xx[j] * fx[j] * pb_xx[j] + pa_xxxx[j] * pb_xx[j]) * s_0_0[j];

                t_xxxx_xy[j] = (3.0 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 2.0 * pa_xxx[j] * fx[j] * pb_y[j] +

                               0.75 * fx[j] * fx[j] * pb_xy[j] + 3.0 * pa_xx[j] * fx[j] * pb_xy[j] + pa_xxxx[j] * pb_xy[j]) * s_0_0[j];

                t_xxxx_xz[j] = (3.0 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 2.0 * pa_xxx[j] * fx[j] * pb_z[j] +

                               0.75 * fx[j] * fx[j] * pb_xz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xz[j] + pa_xxxx[j] * pb_xz[j]) * s_0_0[j];

                t_xxxx_yy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] +

                               0.5 * pa_xxxx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xx[j] * fx[j] * pb_yy[j] +

                               pa_xxxx[j] * pb_yy[j]) * s_0_0[j];

                t_xxxx_yz[j] = (0.75 * fx[j] * fx[j] * pb_yz[j] + 3.0 * pa_xx[j] * fx[j] * pb_yz[j] +

                               pa_xxxx[j] * pb_yz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xy, \
                                     pb_xz, pb_y, pb_yy, pb_z, pb_zz, s_0_0, t_xxxx_zz, t_xxxy_xx, t_xxxy_xy, t_xxxy_xz, \
                                     t_xxxy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_zz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] +

                               0.5 * pa_xxxx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xx[j] * fx[j] * pb_zz[j] +

                               pa_xxxx[j] * pb_zz[j]) * s_0_0[j];

                t_xxxy_xx[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               0.5 * pa_xxxy[j] * fx[j] + 3.0 * pa_xxy[j] * fx[j] * pb_x[j] + 1.5 * pa_xy[j] * fx[j] * pb_xx[j] +

                               pa_xxxy[j] * pb_xx[j]) * s_0_0[j];

                t_xxxy_xy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.5 * pa_xxx[j] * fx[j] * pb_x[j] +

                               1.5 * pa_xxy[j] * fx[j] * pb_y[j] + 1.5 * pa_xy[j] * fx[j] * pb_xy[j] + pa_xxxy[j] * pb_xy[j]) * s_0_0[j];

                t_xxxy_xz[j] = (0.75 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 1.5 * pa_xxy[j] * fx[j] * pb_z[j] +

                               1.5 * pa_xy[j] * fx[j] * pb_xz[j] + pa_xxxy[j] * pb_xz[j]) * s_0_0[j];

                t_xxxy_yy[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.5 * pa_xxxy[j] * fx[j] + pa_xxx[j] * fx[j] * pb_y[j] + 1.5 * pa_xy[j] * fx[j] * pb_yy[j] +

                               pa_xxxy[j] * pb_yy[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxxz, pa_xxz, pa_xy, pa_xz, pa_z, pb_x, pb_xx, \
                                     pb_xy, pb_xz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_xxxy_yz, t_xxxy_zz, t_xxxz_xx, \
                                     t_xxxz_xy, t_xxxz_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_yz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xxx[j] * fx[j] * pb_z[j] +

                               1.5 * pa_xy[j] * fx[j] * pb_yz[j] + pa_xxxy[j] * pb_yz[j]) * s_0_0[j];

                t_xxxy_zz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] + 0.5 * pa_xxxy[j] * fx[j] +

                               1.5 * pa_xy[j] * fx[j] * pb_zz[j] + pa_xxxy[j] * pb_zz[j]) * s_0_0[j];

                t_xxxz_xx[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               0.5 * pa_xxxz[j] * fx[j] + 3.0 * pa_xxz[j] * fx[j] * pb_x[j] + 1.5 * pa_xz[j] * fx[j] * pb_xx[j] +

                               pa_xxxz[j] * pb_xx[j]) * s_0_0[j];

                t_xxxz_xy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 1.5 * pa_xxz[j] * fx[j] * pb_y[j] +

                               1.5 * pa_xz[j] * fx[j] * pb_xy[j] + pa_xxxz[j] * pb_xy[j]) * s_0_0[j];

                t_xxxz_xz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.5 * pa_xxx[j] * fx[j] * pb_x[j] +

                               1.5 * pa_xxz[j] * fx[j] * pb_z[j] + 1.5 * pa_xz[j] * fx[j] * pb_xz[j] + pa_xxxz[j] * pb_xz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_xz, pa_y, \
                                     pa_yy, pb_x, pb_xx, pb_xy, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xxxz_yy, t_xxxz_yz, \
                                     t_xxxz_zz, t_xxyy_xx, t_xxyy_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_yy[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] + 0.5 * pa_xxxz[j] * fx[j] +

                               1.5 * pa_xz[j] * fx[j] * pb_yy[j] + pa_xxxz[j] * pb_yy[j]) * s_0_0[j];

                t_xxxz_yz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xxx[j] * fx[j] * pb_y[j] +

                               1.5 * pa_xz[j] * fx[j] * pb_yz[j] + pa_xxxz[j] * pb_yz[j]) * s_0_0[j];

                t_xxxz_zz[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_xxxz[j] * fx[j] + pa_xxx[j] * fx[j] * pb_z[j] + 1.5 * pa_xz[j] * fx[j] * pb_zz[j] +

                               pa_xxxz[j] * pb_zz[j]) * s_0_0[j];

                t_xxyy_xx[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] +

                               0.25 * pa_xx[j] * fx[j] * fx[j] + pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xxyy[j] * fx[j] +

                               2.0 * pa_xyy[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_xx[j] * fx[j] * pb_xx[j] +

                               0.5 * fx[j] * pa_yy[j] * pb_xx[j] + pa_xxyy[j] * pb_xx[j]) * s_0_0[j];

                t_xxyy_xy[j] = (pa_xy[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.5 * fx[j] * fx[j] * pa_y[j] * pb_x[j] + pa_xxy[j] * fx[j] * pb_x[j] + pa_xyy[j] * fx[j] * pb_y[j] +

                               0.25 * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_yy[j] * pb_xy[j] +

                               pa_xxyy[j] * pb_xy[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xxyz, pa_xyy, pa_xyz, pa_y, pa_yy, pa_yz, \
                                     pb_x, pb_xx, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xxyy_xz, t_xxyy_yy, \
                                     t_xxyy_yz, t_xxyy_zz, t_xxyz_xx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xz[j] = (0.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + pa_xyy[j] * fx[j] * pb_z[j] +

                               0.25 * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_yy[j] * pb_xz[j] +

                               pa_xxyy[j] * pb_xz[j]) * s_0_0[j];

                t_xxyy_yy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_yy[j] + fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.5 * pa_xxyy[j] * fx[j] +

                               2.0 * pa_xxy[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xx[j] * fx[j] * pb_yy[j] +

                               0.5 * fx[j] * pa_yy[j] * pb_yy[j] + pa_xxyy[j] * pb_yy[j]) * s_0_0[j];

                t_xxyy_yz[j] = (0.5 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + pa_xxy[j] * fx[j] * pb_z[j] +

                               0.25 * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_yy[j] * pb_yz[j] +

                               pa_xxyy[j] * pb_yz[j]) * s_0_0[j];

                t_xxyy_zz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_yy[j] + 0.5 * pa_xxyy[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] +

                               0.5 * pa_xx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_yy[j] * pb_zz[j] + pa_xxyy[j] * pb_zz[j]) * s_0_0[j];

                t_xxyz_xx[j] = (0.75 * fx[j] * fx[j] * pa_yz[j] + 0.5 * pa_xxyz[j] * fx[j] +

                               2.0 * pa_xyz[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_yz[j] * pb_xx[j] + pa_xxyz[j] * pb_xx[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, \
                                     pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xxyz_xy, t_xxyz_xz, t_xxyz_yy, \
                                     t_xxyz_yz, t_xxyz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xy[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               0.5 * pa_xxz[j] * fx[j] * pb_x[j] + pa_xyz[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_yz[j] * pb_xy[j] +

                               pa_xxyz[j] * pb_xy[j]) * s_0_0[j];

                t_xxyz_xz[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               0.5 * pa_xxy[j] * fx[j] * pb_x[j] + pa_xyz[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_yz[j] * pb_xz[j] +

                               pa_xxyz[j] * pb_xz[j]) * s_0_0[j];

                t_xxyz_yy[j] = (0.25 * fx[j] * fx[j] * pa_yz[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               0.5 * pa_xxyz[j] * fx[j] + pa_xxz[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_yz[j] * pb_yy[j] +

                               pa_xxyz[j] * pb_yy[j]) * s_0_0[j];

                t_xxyz_yz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.5 * pa_xxy[j] * fx[j] * pb_y[j] +

                               0.5 * pa_xxz[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_yz[j] * pb_yz[j] + pa_xxyz[j] * pb_yz[j]) * s_0_0[j];

                t_xxyz_zz[j] = (0.25 * fx[j] * fx[j] * pa_yz[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                               0.5 * pa_xxyz[j] * fx[j] + pa_xxy[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_yz[j] * pb_zz[j] +

                               pa_xxyz[j] * pb_zz[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xy, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxzz_xx, t_xxzz_xy, t_xxzz_xz, t_xxzz_yy, \
                                     t_xxzz_yz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xx[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] +

                               0.25 * pa_xx[j] * fx[j] * fx[j] + pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xxzz[j] * fx[j] +

                               2.0 * pa_xzz[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_xx[j] * fx[j] * pb_xx[j] +

                               0.5 * fx[j] * pa_zz[j] * pb_xx[j] + pa_xxzz[j] * pb_xx[j]) * s_0_0[j];

                t_xxzz_xy[j] = (0.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + pa_xzz[j] * fx[j] * pb_y[j] +

                               0.25 * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_zz[j] * pb_xy[j] +

                               pa_xxzz[j] * pb_xy[j]) * s_0_0[j];

                t_xxzz_xz[j] = (pa_xz[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + pa_xxz[j] * fx[j] * pb_x[j] + pa_xzz[j] * fx[j] * pb_z[j] +

                               0.25 * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xz[j] +

                               pa_xxzz[j] * pb_xz[j]) * s_0_0[j];

                t_xxzz_yy[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_zz[j] + 0.5 * pa_xxzz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] +

                               0.5 * pa_xx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_zz[j] * pb_yy[j] + pa_xxzz[j] * pb_yy[j]) * s_0_0[j];

                t_xxzz_yz[j] = (0.5 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + pa_xxz[j] * fx[j] * pb_y[j] +

                               0.25 * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yz[j] +

                               pa_xxzz[j] * pb_yz[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, \
                                     pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_z, pb_zz, s_0_0, t_xxzz_zz, \
                                     t_xyyy_xx, t_xyyy_xy, t_xyyy_xz, t_xyyy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_zz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_zz[j] + fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.5 * pa_xxzz[j] * fx[j] +

                               2.0 * pa_xxz[j] * fx[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xx[j] * fx[j] * pb_zz[j] +

                               0.5 * fx[j] * pa_zz[j] * pb_zz[j] + pa_xxzz[j] * pb_zz[j]) * s_0_0[j];

                t_xyyy_xx[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               0.5 * pa_xyyy[j] * fx[j] + fx[j] * pa_yyy[j] * pb_x[j] + 1.5 * pa_xy[j] * fx[j] * pb_xx[j] +

                               pa_xyyy[j] * pb_xx[j]) * s_0_0[j];

                t_xyyy_xy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] +

                               0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 1.5 * pa_xyy[j] * fx[j] * pb_x[j] +

                               0.5 * fx[j] * pa_yyy[j] * pb_y[j] + 1.5 * pa_xy[j] * fx[j] * pb_xy[j] + pa_xyyy[j] * pb_xy[j]) * s_0_0[j];

                t_xyyy_xz[j] = (0.75 * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 0.5 * fx[j] * pa_yyy[j] * pb_z[j] +

                               1.5 * pa_xy[j] * fx[j] * pb_xz[j] + pa_xyyy[j] * pb_xz[j]) * s_0_0[j];

                t_xyyy_yy[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.5 * pa_xyyy[j] * fx[j] + 3.0 * pa_xyy[j] * fx[j] * pb_y[j] + 1.5 * pa_xy[j] * fx[j] * pb_yy[j] +

                               pa_xyyy[j] * pb_yy[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_xyyz, pa_xyz, pa_xz, pa_yy, pa_yyz, pa_yz, \
                                     pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_xyyy_yz, t_xyyy_zz, \
                                     t_xyyz_xx, t_xyyz_xy, t_xyyz_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_yz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xyy[j] * fx[j] * pb_z[j] +

                               1.5 * pa_xy[j] * fx[j] * pb_yz[j] + pa_xyyy[j] * pb_yz[j]) * s_0_0[j];

                t_xyyy_zz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] + 0.5 * pa_xyyy[j] * fx[j] +

                               1.5 * pa_xy[j] * fx[j] * pb_zz[j] + pa_xyyy[j] * pb_zz[j]) * s_0_0[j];

                t_xyyz_xx[j] = (0.25 * pa_xz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               0.5 * pa_xyyz[j] * fx[j] + fx[j] * pa_yyz[j] * pb_x[j] + 0.5 * pa_xz[j] * fx[j] * pb_xx[j] +

                               pa_xyyz[j] * pb_xx[j]) * s_0_0[j];

                t_xyyz_xy[j] = (0.5 * fx[j] * fx[j] * pa_yz[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               pa_xyz[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_yyz[j] * pb_y[j] + 0.5 * pa_xz[j] * fx[j] * pb_xy[j] +

                               pa_xyyz[j] * pb_xy[j]) * s_0_0[j];

                t_xyyz_xz[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] +

                               0.25 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.5 * pa_xyy[j] * fx[j] * pb_x[j] +

                               0.5 * fx[j] * pa_yyz[j] * pb_z[j] + 0.5 * pa_xz[j] * fx[j] * pb_xz[j] + pa_xyyz[j] * pb_xz[j]) * s_0_0[j];
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yzz, \
                                     pa_zz, pb_x, pb_xx, pb_xy, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xyyz_yy, t_xyyz_yz, \
                                     t_xyyz_zz, t_xyzz_xx, t_xyzz_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_yy[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] + 0.5 * pa_xyyz[j] * fx[j] +

                               2.0 * pa_xyz[j] * fx[j] * pb_y[j] + 0.5 * pa_xz[j] * fx[j] * pb_yy[j] + pa_xyyz[j] * pb_yy[j]) * s_0_0[j];

                t_xyyz_yz[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.5 * pa_xyy[j] * fx[j] * pb_y[j] + pa_xyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xz[j] * fx[j] * pb_yz[j] +

                               pa_xyyz[j] * pb_yz[j]) * s_0_0[j];

                t_xyyz_zz[j] = (0.25 * pa_xz[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_xyyz[j] * fx[j] + pa_xyy[j] * fx[j] * pb_z[j] + 0.5 * pa_xz[j] * fx[j] * pb_zz[j] +

                               pa_xyyz[j] * pb_zz[j]) * s_0_0[j];

                t_xyzz_xx[j] = (0.25 * pa_xy[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                               0.5 * pa_xyzz[j] * fx[j] + fx[j] * pa_yzz[j] * pb_x[j] + 0.5 * pa_xy[j] * fx[j] * pb_xx[j] +

                               pa_xyzz[j] * pb_xx[j]) * s_0_0[j];

                t_xyzz_xy[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] +

                               0.25 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.5 * pa_xzz[j] * fx[j] * pb_x[j] +

                               0.5 * fx[j] * pa_yzz[j] * pb_y[j] + 0.5 * pa_xy[j] * fx[j] * pb_xy[j] + pa_xyzz[j] * pb_xy[j]) * s_0_0[j];
            }

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_xzzz, pa_y, pa_yz, pa_yzz, \
                                     pa_z, pa_zzz, pb_x, pb_xx, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xyzz_xz, \
                                     t_xyzz_yy, t_xyzz_yz, t_xyzz_zz, t_xzzz_xx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xz[j] = (0.5 * fx[j] * fx[j] * pa_yz[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                               pa_xyz[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_yzz[j] * pb_z[j] + 0.5 * pa_xy[j] * fx[j] * pb_xz[j] +

                               pa_xyzz[j] * pb_xz[j]) * s_0_0[j];

                t_xyzz_yy[j] = (0.25 * pa_xy[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_y[j] +

                               0.5 * pa_xyzz[j] * fx[j] + pa_xzz[j] * fx[j] * pb_y[j] + 0.5 * pa_xy[j] * fx[j] * pb_yy[j] +

                               pa_xyzz[j] * pb_yy[j]) * s_0_0[j];

                t_xyzz_yz[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               pa_xyz[j] * fx[j] * pb_y[j] + 0.5 * pa_xzz[j] * fx[j] * pb_z[j] + 0.5 * pa_xy[j] * fx[j] * pb_yz[j] +

                               pa_xyzz[j] * pb_yz[j]) * s_0_0[j];

                t_xyzz_zz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] + 0.5 * pa_xyzz[j] * fx[j] +

                               2.0 * pa_xyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xy[j] * fx[j] * pb_zz[j] + pa_xyzz[j] * pb_zz[j]) * s_0_0[j];

                t_xzzz_xx[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                               0.5 * pa_xzzz[j] * fx[j] + fx[j] * pa_zzz[j] * pb_x[j] + 1.5 * pa_xz[j] * fx[j] * pb_xx[j] +

                               pa_xzzz[j] * pb_xx[j]) * s_0_0[j];
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xy, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xzzz_xy, t_xzzz_xz, t_xzzz_yy, t_xzzz_yz, \
                                     t_xzzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.5 * fx[j] * pa_zzz[j] * pb_y[j] +

                               1.5 * pa_xz[j] * fx[j] * pb_xy[j] + pa_xzzz[j] * pb_xy[j]) * s_0_0[j];

                t_xzzz_xz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] +

                               0.75 * pa_x[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 1.5 * pa_xzz[j] * fx[j] * pb_x[j] +

                               0.5 * fx[j] * pa_zzz[j] * pb_z[j] + 1.5 * pa_xz[j] * fx[j] * pb_xz[j] + pa_xzzz[j] * pb_xz[j]) * s_0_0[j];

                t_xzzz_yy[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] + 0.5 * pa_xzzz[j] * fx[j] +

                               1.5 * pa_xz[j] * fx[j] * pb_yy[j] + pa_xzzz[j] * pb_yy[j]) * s_0_0[j];

                t_xzzz_yz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xzz[j] * fx[j] * pb_y[j] +

                               1.5 * pa_xz[j] * fx[j] * pb_yz[j] + pa_xzzz[j] * pb_yz[j]) * s_0_0[j];

                t_xzzz_zz[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_xzzz[j] * fx[j] + 3.0 * pa_xzz[j] * fx[j] * pb_z[j] + 1.5 * pa_xz[j] * fx[j] * pb_zz[j] +

                               pa_xzzz[j] * pb_zz[j]) * s_0_0[j];
            }

            // Batch of Integrals (12) = (60,65)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, s_0_0, t_yyyy_xx, t_yyyy_xy, t_yyyy_xz, t_yyyy_yy, t_yyyy_yz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xx[j] = (0.375 * fx[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] +

                               0.5 * pa_yyyy[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_yy[j] * fx[j] * pb_xx[j] +

                               pa_yyyy[j] * pb_xx[j]) * s_0_0[j];

                t_yyyy_xy[j] = (3.0 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 2.0 * pa_yyy[j] * fx[j] * pb_x[j] +

                               0.75 * fx[j] * fx[j] * pb_xy[j] + 3.0 * pa_yy[j] * fx[j] * pb_xy[j] + pa_yyyy[j] * pb_xy[j]) * s_0_0[j];

                t_yyyy_xz[j] = (0.75 * fx[j] * fx[j] * pb_xz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xz[j] +

                               pa_yyyy[j] * pb_xz[j]) * s_0_0[j];

                t_yyyy_yy[j] = (1.875 * fx[j] * fx[j] * fx[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] +

                               6.0 * pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_yyyy[j] * fx[j] + 4.0 * pa_yyy[j] * fx[j] * pb_y[j] +

                               0.75 * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_yy[j] * fx[j] * pb_yy[j] + pa_yyyy[j] * pb_yy[j]) * s_0_0[j];

                t_yyyy_yz[j] = (3.0 * pa_y[j] * fx[j] * fx[j] * pb_z[j] + 2.0 * pa_yyy[j] * fx[j] * pb_z[j] +

                               0.75 * fx[j] * fx[j] * pb_yz[j] + 3.0 * pa_yy[j] * fx[j] * pb_yz[j] + pa_yyyy[j] * pb_yz[j]) * s_0_0[j];
            }

            // Batch of Integrals (13) = (65,70)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xy, \
                                     pb_xz, pb_y, pb_yy, pb_zz, s_0_0, t_yyyy_zz, t_yyyz_xx, t_yyyz_xy, t_yyyz_xz, \
                                     t_yyyz_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_zz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] +

                               0.5 * pa_yyyy[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_yy[j] * fx[j] * pb_zz[j] +

                               pa_yyyy[j] * pb_zz[j]) * s_0_0[j];

                t_yyyz_xx[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] + 0.5 * pa_yyyz[j] * fx[j] +

                               1.5 * pa_yz[j] * fx[j] * pb_xx[j] + pa_yyyz[j] * pb_xx[j]) * s_0_0[j];

                t_yyyz_xy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 1.5 * pa_yyz[j] * fx[j] * pb_x[j] +

                               1.5 * pa_yz[j] * fx[j] * pb_xy[j] + pa_yyyz[j] * pb_xy[j]) * s_0_0[j];

                t_yyyz_xz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_yyy[j] * fx[j] * pb_x[j] +

                               1.5 * pa_yz[j] * fx[j] * pb_xz[j] + pa_yyyz[j] * pb_xz[j]) * s_0_0[j];

                t_yyyz_yy[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               0.5 * pa_yyyz[j] * fx[j] + 3.0 * pa_yyz[j] * fx[j] * pb_y[j] + 1.5 * pa_yz[j] * fx[j] * pb_yy[j] +

                               pa_yyyz[j] * pb_yy[j]) * s_0_0[j];
            }

            // Batch of Integrals (14) = (70,75)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_yyyz_yz, t_yyyz_zz, \
                                     t_yyzz_xx, t_yyzz_xy, t_yyzz_xz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_yz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] +

                               0.75 * pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.5 * pa_yyy[j] * fx[j] * pb_y[j] +

                               1.5 * pa_yyz[j] * fx[j] * pb_z[j] + 1.5 * pa_yz[j] * fx[j] * pb_yz[j] + pa_yyyz[j] * pb_yz[j]) * s_0_0[j];

                t_yyyz_zz[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_yyyz[j] * fx[j] + pa_yyy[j] * fx[j] * pb_z[j] + 1.5 * pa_yz[j] * fx[j] * pb_zz[j] +

                               pa_yyyz[j] * pb_zz[j]) * s_0_0[j];

                t_yyzz_xx[j] = (0.125 * fx[j] * fx[j] * fx[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_zz[j] + 0.5 * pa_yyzz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_xx[j] +

                               0.5 * pa_yy[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pa_zz[j] * pb_xx[j] + pa_yyzz[j] * pb_xx[j]) * s_0_0[j];

                t_yyzz_xy[j] = (0.5 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + pa_yzz[j] * fx[j] * pb_x[j] +

                               0.25 * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_yy[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_zz[j] * pb_xy[j] +

                               pa_yyzz[j] * pb_xy[j]) * s_0_0[j];

                t_yyzz_xz[j] = (0.5 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + pa_yyz[j] * fx[j] * pb_x[j] +

                               0.25 * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xz[j] +

                               pa_yyzz[j] * pb_xz[j]) * s_0_0[j];
            }

            // Batch of Integrals (15) = (75,80)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, \
                                     pb_x, pb_xx, pb_xy, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_yyzz_yy, t_yyzz_yz, \
                                     t_yyzz_zz, t_yzzz_xx, t_yzzz_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_yy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] +

                               0.25 * pa_yy[j] * fx[j] * fx[j] + pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_yyzz[j] * fx[j] +

                               2.0 * pa_yzz[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_yy[j] * fx[j] * pb_yy[j] +

                               0.5 * fx[j] * pa_zz[j] * pb_yy[j] + pa_yyzz[j] * pb_yy[j]) * s_0_0[j];

                t_yyzz_yz[j] = (pa_yz[j] * fx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * fx[j] * fx[j] * pa_z[j] * pb_y[j] + pa_yyz[j] * fx[j] * pb_y[j] + pa_yzz[j] * fx[j] * pb_z[j] +

                               0.25 * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_yy[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yz[j] +

                               pa_yyzz[j] * pb_yz[j]) * s_0_0[j];

                t_yyzz_zz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] +

                               0.25 * fx[j] * fx[j] * pa_zz[j] + fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.5 * pa_yyzz[j] * fx[j] +

                               2.0 * pa_yyz[j] * fx[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_yy[j] * fx[j] * pb_zz[j] +

                               0.5 * fx[j] * pa_zz[j] * pb_zz[j] + pa_yyzz[j] * pb_zz[j]) * s_0_0[j];

                t_yzzz_xx[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] + 0.5 * pa_yzzz[j] * fx[j] +

                               1.5 * pa_yz[j] * fx[j] * pb_xx[j] + pa_yzzz[j] * pb_xx[j]) * s_0_0[j];

                t_yzzz_xy[j] = (0.75 * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.5 * fx[j] * pa_zzz[j] * pb_x[j] +

                               1.5 * pa_yz[j] * fx[j] * pb_xy[j] + pa_yzzz[j] * pb_xy[j]) * s_0_0[j];
            }

            // Batch of Integrals (16) = (80,85)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xx, pb_xz, \
                                     pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_yzzz_xz, t_yzzz_yy, t_yzzz_yz, t_yzzz_zz, \
                                     t_zzzz_xx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yzz[j] * fx[j] * pb_x[j] +

                               1.5 * pa_yz[j] * fx[j] * pb_xz[j] + pa_yzzz[j] * pb_xz[j]) * s_0_0[j];

                t_yzzz_yy[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                               0.5 * pa_yzzz[j] * fx[j] + fx[j] * pa_zzz[j] * pb_y[j] + 1.5 * pa_yz[j] * fx[j] * pb_yy[j] +

                               pa_yzzz[j] * pb_yy[j]) * s_0_0[j];

                t_yzzz_yz[j] = (0.375 * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] +

                               0.75 * pa_y[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 1.5 * pa_yzz[j] * fx[j] * pb_y[j] +

                               0.5 * fx[j] * pa_zzz[j] * pb_z[j] + 1.5 * pa_yz[j] * fx[j] * pb_yz[j] + pa_yzzz[j] * pb_yz[j]) * s_0_0[j];

                t_yzzz_zz[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_z[j] +

                               0.5 * pa_yzzz[j] * fx[j] + 3.0 * pa_yzz[j] * fx[j] * pb_z[j] + 1.5 * pa_yz[j] * fx[j] * pb_zz[j] +

                               pa_yzzz[j] * pb_zz[j]) * s_0_0[j];

                t_zzzz_xx[j] = (0.375 * fx[j] * fx[j] * fx[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] +

                               0.5 * pa_zzzz[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_zz[j] * fx[j] * pb_xx[j] +

                               pa_zzzz[j] * pb_xx[j]) * s_0_0[j];
            }

            // Batch of Integrals (17) = (85,90)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_zzzz_xy, t_zzzz_xz, t_zzzz_yy, t_zzzz_yz, t_zzzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xy[j] = (0.75 * fx[j] * fx[j] * pb_xy[j] + 3.0 * pa_zz[j] * fx[j] * pb_xy[j] +

                               pa_zzzz[j] * pb_xy[j]) * s_0_0[j];

                t_zzzz_xz[j] = (3.0 * pa_z[j] * fx[j] * fx[j] * pb_x[j] + 2.0 * pa_zzz[j] * fx[j] * pb_x[j] +

                               0.75 * fx[j] * fx[j] * pb_xz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xz[j] + pa_zzzz[j] * pb_xz[j]) * s_0_0[j];

                t_zzzz_yy[j] = (0.375 * fx[j] * fx[j] * fx[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] +

                               0.5 * pa_zzzz[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_zz[j] * fx[j] * pb_yy[j] +

                               pa_zzzz[j] * pb_yy[j]) * s_0_0[j];

                t_zzzz_yz[j] = (3.0 * pa_z[j] * fx[j] * fx[j] * pb_y[j] + 2.0 * pa_zzz[j] * fx[j] * pb_y[j] +

                               0.75 * fx[j] * fx[j] * pb_yz[j] + 3.0 * pa_zz[j] * fx[j] * pb_yz[j] + pa_zzzz[j] * pb_yz[j]) * s_0_0[j];

                t_zzzz_zz[j] = (1.875 * fx[j] * fx[j] * fx[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] +

                               6.0 * pa_z[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_zzzz[j] * fx[j] + 4.0 * pa_zzz[j] * fx[j] * pb_z[j] +

                               0.75 * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_zz[j] * fx[j] * pb_zz[j] + pa_zzzz[j] * pb_zz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForGF(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(19 * idx);

            auto pb_y = pbDistances.data(19 * idx + 1);

            auto pb_z = pbDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(19 * idx + 3);

            auto pb_xy = pbDistances.data(19 * idx + 4);

            auto pb_xz = pbDistances.data(19 * idx + 5);

            auto pb_yy = pbDistances.data(19 * idx + 6);

            auto pb_yz = pbDistances.data(19 * idx + 7);

            auto pb_zz = pbDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(19 * idx + 9);

            auto pb_xxy = pbDistances.data(19 * idx + 10);

            auto pb_xxz = pbDistances.data(19 * idx + 11);

            auto pb_xyy = pbDistances.data(19 * idx + 12);

            auto pb_xyz = pbDistances.data(19 * idx + 13);

            auto pb_xzz = pbDistances.data(19 * idx + 14);

            auto pb_yyy = pbDistances.data(19 * idx + 15);

            auto pb_yyz = pbDistances.data(19 * idx + 16);

            auto pb_yzz = pbDistances.data(19 * idx + 17);

            auto pb_zzz = pbDistances.data(19 * idx + 18);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxxx_xxx = primBuffer.data(150 * idx);

            auto t_xxxx_xxy = primBuffer.data(150 * idx + 1);

            auto t_xxxx_xxz = primBuffer.data(150 * idx + 2);

            auto t_xxxx_xyy = primBuffer.data(150 * idx + 3);

            auto t_xxxx_xyz = primBuffer.data(150 * idx + 4);

            auto t_xxxx_xzz = primBuffer.data(150 * idx + 5);

            auto t_xxxx_yyy = primBuffer.data(150 * idx + 6);

            auto t_xxxx_yyz = primBuffer.data(150 * idx + 7);

            auto t_xxxx_yzz = primBuffer.data(150 * idx + 8);

            auto t_xxxx_zzz = primBuffer.data(150 * idx + 9);

            auto t_xxxy_xxx = primBuffer.data(150 * idx + 10);

            auto t_xxxy_xxy = primBuffer.data(150 * idx + 11);

            auto t_xxxy_xxz = primBuffer.data(150 * idx + 12);

            auto t_xxxy_xyy = primBuffer.data(150 * idx + 13);

            auto t_xxxy_xyz = primBuffer.data(150 * idx + 14);

            auto t_xxxy_xzz = primBuffer.data(150 * idx + 15);

            auto t_xxxy_yyy = primBuffer.data(150 * idx + 16);

            auto t_xxxy_yyz = primBuffer.data(150 * idx + 17);

            auto t_xxxy_yzz = primBuffer.data(150 * idx + 18);

            auto t_xxxy_zzz = primBuffer.data(150 * idx + 19);

            auto t_xxxz_xxx = primBuffer.data(150 * idx + 20);

            auto t_xxxz_xxy = primBuffer.data(150 * idx + 21);

            auto t_xxxz_xxz = primBuffer.data(150 * idx + 22);

            auto t_xxxz_xyy = primBuffer.data(150 * idx + 23);

            auto t_xxxz_xyz = primBuffer.data(150 * idx + 24);

            auto t_xxxz_xzz = primBuffer.data(150 * idx + 25);

            auto t_xxxz_yyy = primBuffer.data(150 * idx + 26);

            auto t_xxxz_yyz = primBuffer.data(150 * idx + 27);

            auto t_xxxz_yzz = primBuffer.data(150 * idx + 28);

            auto t_xxxz_zzz = primBuffer.data(150 * idx + 29);

            auto t_xxyy_xxx = primBuffer.data(150 * idx + 30);

            auto t_xxyy_xxy = primBuffer.data(150 * idx + 31);

            auto t_xxyy_xxz = primBuffer.data(150 * idx + 32);

            auto t_xxyy_xyy = primBuffer.data(150 * idx + 33);

            auto t_xxyy_xyz = primBuffer.data(150 * idx + 34);

            auto t_xxyy_xzz = primBuffer.data(150 * idx + 35);

            auto t_xxyy_yyy = primBuffer.data(150 * idx + 36);

            auto t_xxyy_yyz = primBuffer.data(150 * idx + 37);

            auto t_xxyy_yzz = primBuffer.data(150 * idx + 38);

            auto t_xxyy_zzz = primBuffer.data(150 * idx + 39);

            auto t_xxyz_xxx = primBuffer.data(150 * idx + 40);

            auto t_xxyz_xxy = primBuffer.data(150 * idx + 41);

            auto t_xxyz_xxz = primBuffer.data(150 * idx + 42);

            auto t_xxyz_xyy = primBuffer.data(150 * idx + 43);

            auto t_xxyz_xyz = primBuffer.data(150 * idx + 44);

            auto t_xxyz_xzz = primBuffer.data(150 * idx + 45);

            auto t_xxyz_yyy = primBuffer.data(150 * idx + 46);

            auto t_xxyz_yyz = primBuffer.data(150 * idx + 47);

            auto t_xxyz_yzz = primBuffer.data(150 * idx + 48);

            auto t_xxyz_zzz = primBuffer.data(150 * idx + 49);

            auto t_xxzz_xxx = primBuffer.data(150 * idx + 50);

            auto t_xxzz_xxy = primBuffer.data(150 * idx + 51);

            auto t_xxzz_xxz = primBuffer.data(150 * idx + 52);

            auto t_xxzz_xyy = primBuffer.data(150 * idx + 53);

            auto t_xxzz_xyz = primBuffer.data(150 * idx + 54);

            auto t_xxzz_xzz = primBuffer.data(150 * idx + 55);

            auto t_xxzz_yyy = primBuffer.data(150 * idx + 56);

            auto t_xxzz_yyz = primBuffer.data(150 * idx + 57);

            auto t_xxzz_yzz = primBuffer.data(150 * idx + 58);

            auto t_xxzz_zzz = primBuffer.data(150 * idx + 59);

            auto t_xyyy_xxx = primBuffer.data(150 * idx + 60);

            auto t_xyyy_xxy = primBuffer.data(150 * idx + 61);

            auto t_xyyy_xxz = primBuffer.data(150 * idx + 62);

            auto t_xyyy_xyy = primBuffer.data(150 * idx + 63);

            auto t_xyyy_xyz = primBuffer.data(150 * idx + 64);

            auto t_xyyy_xzz = primBuffer.data(150 * idx + 65);

            auto t_xyyy_yyy = primBuffer.data(150 * idx + 66);

            auto t_xyyy_yyz = primBuffer.data(150 * idx + 67);

            auto t_xyyy_yzz = primBuffer.data(150 * idx + 68);

            auto t_xyyy_zzz = primBuffer.data(150 * idx + 69);

            auto t_xyyz_xxx = primBuffer.data(150 * idx + 70);

            auto t_xyyz_xxy = primBuffer.data(150 * idx + 71);

            auto t_xyyz_xxz = primBuffer.data(150 * idx + 72);

            auto t_xyyz_xyy = primBuffer.data(150 * idx + 73);

            auto t_xyyz_xyz = primBuffer.data(150 * idx + 74);

            auto t_xyyz_xzz = primBuffer.data(150 * idx + 75);

            auto t_xyyz_yyy = primBuffer.data(150 * idx + 76);

            auto t_xyyz_yyz = primBuffer.data(150 * idx + 77);

            auto t_xyyz_yzz = primBuffer.data(150 * idx + 78);

            auto t_xyyz_zzz = primBuffer.data(150 * idx + 79);

            auto t_xyzz_xxx = primBuffer.data(150 * idx + 80);

            auto t_xyzz_xxy = primBuffer.data(150 * idx + 81);

            auto t_xyzz_xxz = primBuffer.data(150 * idx + 82);

            auto t_xyzz_xyy = primBuffer.data(150 * idx + 83);

            auto t_xyzz_xyz = primBuffer.data(150 * idx + 84);

            auto t_xyzz_xzz = primBuffer.data(150 * idx + 85);

            auto t_xyzz_yyy = primBuffer.data(150 * idx + 86);

            auto t_xyzz_yyz = primBuffer.data(150 * idx + 87);

            auto t_xyzz_yzz = primBuffer.data(150 * idx + 88);

            auto t_xyzz_zzz = primBuffer.data(150 * idx + 89);

            auto t_xzzz_xxx = primBuffer.data(150 * idx + 90);

            auto t_xzzz_xxy = primBuffer.data(150 * idx + 91);

            auto t_xzzz_xxz = primBuffer.data(150 * idx + 92);

            auto t_xzzz_xyy = primBuffer.data(150 * idx + 93);

            auto t_xzzz_xyz = primBuffer.data(150 * idx + 94);

            auto t_xzzz_xzz = primBuffer.data(150 * idx + 95);

            auto t_xzzz_yyy = primBuffer.data(150 * idx + 96);

            auto t_xzzz_yyz = primBuffer.data(150 * idx + 97);

            auto t_xzzz_yzz = primBuffer.data(150 * idx + 98);

            auto t_xzzz_zzz = primBuffer.data(150 * idx + 99);

            auto t_yyyy_xxx = primBuffer.data(150 * idx + 100);

            auto t_yyyy_xxy = primBuffer.data(150 * idx + 101);

            auto t_yyyy_xxz = primBuffer.data(150 * idx + 102);

            auto t_yyyy_xyy = primBuffer.data(150 * idx + 103);

            auto t_yyyy_xyz = primBuffer.data(150 * idx + 104);

            auto t_yyyy_xzz = primBuffer.data(150 * idx + 105);

            auto t_yyyy_yyy = primBuffer.data(150 * idx + 106);

            auto t_yyyy_yyz = primBuffer.data(150 * idx + 107);

            auto t_yyyy_yzz = primBuffer.data(150 * idx + 108);

            auto t_yyyy_zzz = primBuffer.data(150 * idx + 109);

            auto t_yyyz_xxx = primBuffer.data(150 * idx + 110);

            auto t_yyyz_xxy = primBuffer.data(150 * idx + 111);

            auto t_yyyz_xxz = primBuffer.data(150 * idx + 112);

            auto t_yyyz_xyy = primBuffer.data(150 * idx + 113);

            auto t_yyyz_xyz = primBuffer.data(150 * idx + 114);

            auto t_yyyz_xzz = primBuffer.data(150 * idx + 115);

            auto t_yyyz_yyy = primBuffer.data(150 * idx + 116);

            auto t_yyyz_yyz = primBuffer.data(150 * idx + 117);

            auto t_yyyz_yzz = primBuffer.data(150 * idx + 118);

            auto t_yyyz_zzz = primBuffer.data(150 * idx + 119);

            auto t_yyzz_xxx = primBuffer.data(150 * idx + 120);

            auto t_yyzz_xxy = primBuffer.data(150 * idx + 121);

            auto t_yyzz_xxz = primBuffer.data(150 * idx + 122);

            auto t_yyzz_xyy = primBuffer.data(150 * idx + 123);

            auto t_yyzz_xyz = primBuffer.data(150 * idx + 124);

            auto t_yyzz_xzz = primBuffer.data(150 * idx + 125);

            auto t_yyzz_yyy = primBuffer.data(150 * idx + 126);

            auto t_yyzz_yyz = primBuffer.data(150 * idx + 127);

            auto t_yyzz_yzz = primBuffer.data(150 * idx + 128);

            auto t_yyzz_zzz = primBuffer.data(150 * idx + 129);

            auto t_yzzz_xxx = primBuffer.data(150 * idx + 130);

            auto t_yzzz_xxy = primBuffer.data(150 * idx + 131);

            auto t_yzzz_xxz = primBuffer.data(150 * idx + 132);

            auto t_yzzz_xyy = primBuffer.data(150 * idx + 133);

            auto t_yzzz_xyz = primBuffer.data(150 * idx + 134);

            auto t_yzzz_xzz = primBuffer.data(150 * idx + 135);

            auto t_yzzz_yyy = primBuffer.data(150 * idx + 136);

            auto t_yzzz_yyz = primBuffer.data(150 * idx + 137);

            auto t_yzzz_yzz = primBuffer.data(150 * idx + 138);

            auto t_yzzz_zzz = primBuffer.data(150 * idx + 139);

            auto t_zzzz_xxx = primBuffer.data(150 * idx + 140);

            auto t_zzzz_xxy = primBuffer.data(150 * idx + 141);

            auto t_zzzz_xxz = primBuffer.data(150 * idx + 142);

            auto t_zzzz_xyy = primBuffer.data(150 * idx + 143);

            auto t_zzzz_xyz = primBuffer.data(150 * idx + 144);

            auto t_zzzz_xzz = primBuffer.data(150 * idx + 145);

            auto t_zzzz_yyy = primBuffer.data(150 * idx + 146);

            auto t_zzzz_yyz = primBuffer.data(150 * idx + 147);

            auto t_zzzz_yzz = primBuffer.data(150 * idx + 148);

            auto t_zzzz_zzz = primBuffer.data(150 * idx + 149);

            // Batch of Integrals (0) = (0,2)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xy, pb_y, s_0_0, \
                                     t_xxxx_xxx, t_xxxx_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xxx[j] = (7.5 * pa_x[j] * fx[j] * fx[j] * fx[j] + 5.625 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                3.0 * pa_xxx[j] * fx[j] * fx[j] + 13.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] +

                                9.0 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xxxx[j] * pb_x[j] * fx[j] + 6.0 * pa_xxx[j] * fx[j] * pb_xx[j] +

                                0.75 * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxx[j] + pa_xxxx[j] * pb_xxx[j]) * s_0_0[j];

                t_xxxx_xxy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                6.0 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xxxx[j] * fx[j] * pb_y[j] + 4.0 * pa_xxx[j] * fx[j] * pb_xy[j] +

                                0.75 * fx[j] * fx[j] * pb_xxy[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxy[j] + pa_xxxx[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (2,4)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xxz, pb_xyy, pb_xz, pb_yy, pb_z, s_0_0, \
                                     t_xxxx_xxz, t_xxxx_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xxz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                6.0 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xxxx[j] * fx[j] * pb_z[j] + 4.0 * pa_xxx[j] * fx[j] * pb_xz[j] +

                                0.75 * fx[j] * fx[j] * pb_xxz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxz[j] + pa_xxxx[j] * pb_xxz[j]) * s_0_0[j];

                t_xxxx_xyy[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] + pa_xxx[j] * fx[j] * fx[j] +

                                0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] +

                                0.5 * pa_xxxx[j] * pb_x[j] * fx[j] + 2.0 * pa_xxx[j] * fx[j] * pb_yy[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] +

                                3.0 * pa_xx[j] * fx[j] * pb_xyy[j] + pa_xxxx[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (4,6)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xyz, pb_xzz, pb_yz, pb_zz, s_0_0, \
                                     t_xxxx_xyz, t_xxxx_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xyz[j] = (3.0 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 2.0 * pa_xxx[j] * fx[j] * pb_yz[j] +

                                0.75 * fx[j] * fx[j] * pb_xyz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xyz[j] + pa_xxxx[j] * pb_xyz[j]) * s_0_0[j];

                t_xxxx_xzz[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] + pa_xxx[j] * fx[j] * fx[j] +

                                0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] +

                                0.5 * pa_xxxx[j] * pb_x[j] * fx[j] + 2.0 * pa_xxx[j] * fx[j] * pb_zz[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] +

                                3.0 * pa_xx[j] * fx[j] * pb_xzz[j] + pa_xxxx[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (6,8)

            #pragma omp simd aligned(fx, pa_xx, pa_xxxx, pb_y, pb_yyy, pb_yyz, pb_z, s_0_0, t_xxxx_yyy, t_xxxx_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_yyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * pa_xxxx[j] * pb_y[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_xx[j] * fx[j] * pb_yyy[j] +

                                pa_xxxx[j] * pb_yyy[j]) * s_0_0[j];

                t_xxxx_yyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_xxxx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] + 3.0 * pa_xx[j] * fx[j] * pb_yyz[j] +

                                pa_xxxx[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (8,10)

            #pragma omp simd aligned(fx, pa_xx, pa_xxxx, pb_y, pb_yzz, pb_z, pb_zzz, s_0_0, t_xxxx_yzz, t_xxxx_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_yzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.5 * pa_xxxx[j] * pb_y[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] + 3.0 * pa_xx[j] * fx[j] * pb_yzz[j] +

                                pa_xxxx[j] * pb_yzz[j]) * s_0_0[j];

                t_xxxx_zzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_xxxx[j] * pb_z[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_xx[j] * fx[j] * pb_zzz[j] +

                                pa_xxxx[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (10,12)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xy, pb_y, s_0_0, t_xxxy_xxx, t_xxxy_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xxx[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_y[j] + 2.25 * pa_xxy[j] * fx[j] * fx[j] +

                                6.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 1.5 * pa_xxxy[j] * pb_x[j] * fx[j] +

                                4.5 * pa_xxy[j] * fx[j] * pb_xx[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxx[j] + pa_xxxy[j] * pb_xxx[j]) * s_0_0[j];

                t_xxxy_xxy[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xxx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] +

                                2.25 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.5 * pa_xxxy[j] * fx[j] * pb_y[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xx[j] +

                                3.0 * pa_xxy[j] * fx[j] * pb_xy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxy[j] + pa_xxxy[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (12,14)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xz, pb_y, pb_yy, pb_z, s_0_0, t_xxxy_xxz, t_xxxy_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xxz[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xxxy[j] * fx[j] * pb_z[j] + 3.0 * pa_xxy[j] * fx[j] * pb_xz[j] +

                                1.5 * pa_xy[j] * fx[j] * pb_xxz[j] + pa_xxxy[j] * pb_xxz[j]) * s_0_0[j];

                t_xxxy_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_xxy[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] +

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.5 * pa_xxxy[j] * pb_x[j] * fx[j] + pa_xxx[j] * fx[j] * pb_xy[j] +

                                1.5 * pa_xxy[j] * fx[j] * pb_yy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyy[j] + pa_xxxy[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (14,16)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xyz, pb_xz, pb_xzz, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_xxxy_xyz, t_xxxy_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xz[j] +

                                1.5 * pa_xxy[j] * fx[j] * pb_yz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyz[j] + pa_xxxy[j] * pb_xyz[j]) * s_0_0[j];

                t_xxxy_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] +

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.5 * pa_xxxy[j] * pb_x[j] * fx[j] +

                                1.5 * pa_xxy[j] * fx[j] * pb_zz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xzz[j] + pa_xxxy[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (16,18)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxy, pa_xy, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_z, s_0_0, \
                                     t_xxxy_yyy, t_xxxy_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_yyy[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] +

                                2.25 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xxxy[j] * pb_y[j] * fx[j] +

                                1.5 * pa_xxx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xy[j] * fx[j] * pb_yyy[j] + pa_xxxy[j] * pb_yyy[j]) * s_0_0[j];

                t_xxxy_yyz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xxxy[j] * fx[j] * pb_z[j] + pa_xxx[j] * fx[j] * pb_yz[j] +

                                1.5 * pa_xy[j] * fx[j] * pb_yyz[j] + pa_xxxy[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (9) = (18,20)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxy, pa_xy, pb_y, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xxxy_yzz, t_xxxy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_yzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] +

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xxxy[j] * pb_y[j] * fx[j] +

                                0.5 * pa_xxx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xy[j] * fx[j] * pb_yzz[j] + pa_xxxy[j] * pb_yzz[j]) * s_0_0[j];

                t_xxxy_zzz[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xxxy[j] * pb_z[j] * fx[j] +

                                1.5 * pa_xy[j] * fx[j] * pb_zzz[j] + pa_xxxy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (10) = (20,22)

            #pragma omp simd aligned(fx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xy, pb_y, s_0_0, \
                                     t_xxxz_xxx, t_xxxz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xxx[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] + 2.25 * pa_xxz[j] * fx[j] * fx[j] +

                                6.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 1.5 * pa_xxxz[j] * pb_x[j] * fx[j] +

                                4.5 * pa_xxz[j] * fx[j] * pb_xx[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxx[j] + pa_xxxz[j] * pb_xxx[j]) * s_0_0[j];

                t_xxxz_xxy[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xxxz[j] * fx[j] * pb_y[j] + 3.0 * pa_xxz[j] * fx[j] * pb_xy[j] +

                                1.5 * pa_xz[j] * fx[j] * pb_xxy[j] + pa_xxxz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (11) = (22,24)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxz, pb_xyy, \
                                     pb_xz, pb_yy, pb_z, s_0_0, t_xxxz_xxz, t_xxxz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xxz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xxx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] +

                                2.25 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_xxxz[j] * fx[j] * pb_z[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xx[j] +

                                3.0 * pa_xxz[j] * fx[j] * pb_xz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxz[j] + pa_xxxz[j] * pb_xxz[j]) * s_0_0[j];

                t_xxxz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_xxxz[j] * pb_x[j] * fx[j] +

                                1.5 * pa_xxz[j] * fx[j] * pb_yy[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyy[j] + pa_xxxz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (12) = (24,26)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_xxxz_xyz, t_xxxz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xy[j] +

                                1.5 * pa_xxz[j] * fx[j] * pb_yz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyz[j] + pa_xxxz[j] * pb_xyz[j]) * s_0_0[j];

                t_xxxz_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_xxz[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xxxz[j] * pb_x[j] * fx[j] + pa_xxx[j] * fx[j] * pb_xz[j] +

                                1.5 * pa_xxz[j] * fx[j] * pb_zz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xzz[j] + pa_xxxz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (13) = (26,28)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxz, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_z, s_0_0, \
                                     t_xxxz_yyy, t_xxxz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_yyy[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xxxz[j] * pb_y[j] * fx[j] +

                                1.5 * pa_xz[j] * fx[j] * pb_yyy[j] + pa_xxxz[j] * pb_yyy[j]) * s_0_0[j];

                t_xxxz_yyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xxxz[j] * fx[j] * pb_z[j] +

                                0.5 * pa_xxx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyz[j] + pa_xxxz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (14) = (28,30)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxz, pa_xz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xxxz_yzz, t_xxxz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_yzz[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xxxz[j] * pb_y[j] * fx[j] + pa_xxx[j] * fx[j] * pb_yz[j] +

                                1.5 * pa_xz[j] * fx[j] * pb_yzz[j] + pa_xxxz[j] * pb_yzz[j]) * s_0_0[j];

                t_xxxz_zzz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] +

                                2.25 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xxxz[j] * pb_z[j] * fx[j] +

                                1.5 * pa_xxx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xz[j] * fx[j] * pb_zzz[j] + pa_xxxz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (15) = (30,32)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, \
                                     pb_xxy, pb_xy, pb_y, s_0_0, t_xxyy_xxx, t_xxyy_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xxx[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                1.5 * pa_xyy[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] +

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xxyy[j] * pb_x[j] * fx[j] +

                                3.0 * pa_xyy[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxx[j] +

                                0.5 * fx[j] * pa_yy[j] * pb_xxx[j] + pa_xxyy[j] * pb_xxx[j]) * s_0_0[j];

                t_xxyy_xxy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.5 * pa_xxy[j] * fx[j] * fx[j] + 2.0 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + pa_x[j] * fx[j] * fx[j] * pb_xy[j] +

                                0.5 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.5 * pa_xxyy[j] * fx[j] * pb_y[j] + pa_xxy[j] * fx[j] * pb_xx[j] +

                                2.0 * pa_xyy[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxy[j] +

                                0.5 * fx[j] * pa_yy[j] * pb_xxy[j] + pa_xxyy[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (16) = (32,34)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xz, pb_y, pb_yy, pb_z, s_0_0, t_xxyy_xxz, t_xxyy_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] +

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xxyy[j] * fx[j] * pb_z[j] +

                                2.0 * pa_xyy[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxz[j] +

                                0.5 * fx[j] * pa_yy[j] * pb_xxz[j] + pa_xxyy[j] * pb_xxz[j]) * s_0_0[j];

                t_xxyy_xyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xyy[j] * fx[j] * fx[j] + 2.0 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] +

                                0.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + fx[j] * fx[j] * pa_y[j] * pb_xy[j] +

                                0.5 * pa_xxyy[j] * pb_x[j] * fx[j] + 2.0 * pa_xxy[j] * fx[j] * pb_xy[j] + pa_xyy[j] * fx[j] * pb_yy[j] +

                                0.25 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_yy[j] * pb_xyy[j] +

                                pa_xxyy[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (17) = (34,36)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xyz, pb_xz, \
                                     pb_xzz, pb_yz, pb_z, pb_zz, s_0_0, t_xxyy_xyz, t_xxyy_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xyz[j] = (pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] +

                                0.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + pa_xxy[j] * fx[j] * pb_xz[j] + pa_xyy[j] * fx[j] * pb_yz[j] +

                                0.25 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yy[j] * pb_xyz[j] +

                                pa_xxyy[j] * pb_xyz[j]) * s_0_0[j];

                t_xxyy_xzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * fx[j] +

                                0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 0.5 * pa_xxyy[j] * pb_x[j] * fx[j] +

                                pa_xyy[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xzz[j] +

                                0.5 * fx[j] * pa_yy[j] * pb_xzz[j] + pa_xxyy[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (18) = (36,38)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyy, pa_y, pa_yy, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_z, \
                                     s_0_0, t_xxyy_yyy, t_xxyy_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_yyy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * pa_xxy[j] * fx[j] * fx[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 1.5 * pa_xxyy[j] * pb_y[j] * fx[j] +

                                3.0 * pa_xxy[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyy[j] +

                                0.5 * fx[j] * pa_yy[j] * pb_yyy[j] + pa_xxyy[j] * pb_yyy[j]) * s_0_0[j];

                t_xxyy_yyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.5 * pa_xxyy[j] * fx[j] * pb_z[j] +

                                2.0 * pa_xxy[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyz[j] +

                                0.5 * fx[j] * pa_yy[j] * pb_yyz[j] + pa_xxyy[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (19) = (38,40)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyy, pa_y, pa_yy, pb_y, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xxyy_yzz, t_xxyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_yzz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xxy[j] * fx[j] * fx[j] +

                                0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.5 * pa_xxyy[j] * pb_y[j] * fx[j] +

                                pa_xxy[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yzz[j] +

                                0.5 * fx[j] * pa_yy[j] * pb_yzz[j] + pa_xxyy[j] * pb_yzz[j]) * s_0_0[j];

                t_xxyy_zzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + 1.5 * pa_xxyy[j] * pb_z[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] +

                                0.5 * pa_xx[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_yy[j] * pb_zzz[j] + pa_xxyy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (20) = (40,42)

            #pragma omp simd aligned(fx, pa_xxyz, pa_xxz, pa_xyz, pa_xz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xy, pb_y, s_0_0, t_xxyz_xxx, t_xxyz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xxx[j] = (1.5 * pa_xyz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] +

                                1.5 * pa_xxyz[j] * pb_x[j] * fx[j] + 3.0 * pa_xyz[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxx[j] +

                                pa_xxyz[j] * pb_xxx[j]) * s_0_0[j];

                t_xxyz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * pa_xxz[j] * fx[j] * fx[j] +

                                pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_y[j] + 0.5 * pa_xxz[j] * fx[j] * pb_xx[j] +

                                2.0 * pa_xyz[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxy[j] + pa_xxyz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (21) = (42,44)

            #pragma omp simd aligned(fx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xz, pb_y, pb_yy, pb_z, s_0_0, t_xxyz_xxz, t_xxyz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * pa_xxy[j] * fx[j] * fx[j] +

                                pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xxy[j] * fx[j] * pb_xx[j] +

                                2.0 * pa_xyz[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxz[j] + pa_xxyz[j] * pb_xxz[j]) * s_0_0[j];

                t_xxyz_xyy[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xxyz[j] * pb_x[j] * fx[j] +

                                pa_xxz[j] * fx[j] * pb_xy[j] + pa_xyz[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_yz[j] * pb_xyy[j] +

                                pa_xxyz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (22) = (44,46)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, \
                                     pb_x, pb_xy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_xxyz_xyz, \
                                     t_xxyz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xyz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_xxy[j] * fx[j] * pb_xy[j] +

                                0.5 * pa_xxz[j] * fx[j] * pb_xz[j] + pa_xyz[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xyz[j] +

                                pa_xxyz[j] * pb_xyz[j]) * s_0_0[j];

                t_xxyz_xzz[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xxyz[j] * pb_x[j] * fx[j] +

                                pa_xxy[j] * fx[j] * pb_xz[j] + pa_xyz[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xzz[j] +

                                pa_xxyz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (23) = (46,48)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_y, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_xxyz_yyy, t_xxyz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_yyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] +

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 1.5 * pa_xxyz[j] * pb_y[j] * fx[j] +

                                1.5 * pa_xxz[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyy[j] + pa_xxyz[j] * pb_yyy[j]) * s_0_0[j];

                t_xxyz_yyz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_xxy[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] +

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xxy[j] * fx[j] * pb_yy[j] +

                                pa_xxz[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyz[j] + pa_xxyz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (24) = (48,50)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_y, pa_yz, pa_z, pb_y, pb_yz, pb_yzz, pb_z, \
                                     pb_zz, pb_zzz, s_0_0, t_xxyz_yzz, t_xxyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_yzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_xxz[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xxyz[j] * pb_y[j] * fx[j] + pa_xxy[j] * fx[j] * pb_yz[j] +

                                0.5 * pa_xxz[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yzz[j] + pa_xxyz[j] * pb_yzz[j]) * s_0_0[j];

                t_xxyz_zzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] +

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 1.5 * pa_xxyz[j] * pb_z[j] * fx[j] +

                                1.5 * pa_xxy[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_yz[j] * pb_zzz[j] + pa_xxyz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (25) = (50,52)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxzz, pa_xzz, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xy, pb_y, \
                                     s_0_0, t_xxzz_xxx, t_xxzz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xxx[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                1.5 * pa_xzz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] +

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xxzz[j] * pb_x[j] * fx[j] +

                                3.0 * pa_xzz[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxx[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_xxx[j] + pa_xxzz[j] * pb_xxx[j]) * s_0_0[j];

                t_xxzz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] +

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xxzz[j] * fx[j] * pb_y[j] +

                                2.0 * pa_xzz[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxy[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_xxy[j] + pa_xxzz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (26) = (52,54)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, \
                                     pb_xyy, pb_xz, pb_yy, pb_z, s_0_0, t_xxzz_xxz, t_xxzz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xxz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_xxz[j] * fx[j] * fx[j] + 2.0 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + pa_x[j] * fx[j] * fx[j] * pb_xz[j] +

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_xxzz[j] * fx[j] * pb_z[j] + pa_xxz[j] * fx[j] * pb_xx[j] +

                                2.0 * pa_xzz[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxz[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_xxz[j] + pa_xxzz[j] * pb_xxz[j]) * s_0_0[j];

                t_xxzz_xyy[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] * fx[j] +

                                0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 0.5 * pa_xxzz[j] * pb_x[j] * fx[j] +

                                pa_xzz[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyy[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_xyy[j] + pa_xxzz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (27) = (54,56)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xy, pb_xyz, \
                                     pb_xz, pb_xzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_xxzz_xyz, t_xxzz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xyz[j] = (pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] +

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + pa_xxz[j] * fx[j] * pb_xy[j] + pa_xzz[j] * fx[j] * pb_yz[j] +

                                0.25 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xyz[j] +

                                pa_xxzz[j] * pb_xyz[j]) * s_0_0[j];

                t_xxzz_xzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xzz[j] * fx[j] * fx[j] + 2.0 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xz[j] +

                                0.5 * pa_xxzz[j] * pb_x[j] * fx[j] + 2.0 * pa_xxz[j] * fx[j] * pb_xz[j] + pa_xzz[j] * fx[j] * pb_zz[j] +

                                0.25 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xzz[j] +

                                pa_xxzz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (28) = (56,58)

            #pragma omp simd aligned(fx, pa_xx, pa_xxz, pa_xxzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_z, s_0_0, \
                                     t_xxzz_yyy, t_xxzz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_yyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 1.5 * pa_xxzz[j] * pb_y[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] +

                                0.5 * pa_xx[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyy[j] + pa_xxzz[j] * pb_yyy[j]) * s_0_0[j];

                t_xxzz_yyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_xxz[j] * fx[j] * fx[j] +

                                0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_xxzz[j] * fx[j] * pb_z[j] +

                                pa_xxz[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyz[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_yyz[j] + pa_xxzz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (29) = (58,60)

            #pragma omp simd aligned(fx, pa_xx, pa_xxz, pa_xxzz, pa_z, pa_zz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     s_0_0, t_xxzz_yzz, t_xxzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_yzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xxzz[j] * pb_y[j] * fx[j] +

                                2.0 * pa_xxz[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yzz[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_yzz[j] + pa_xxzz[j] * pb_yzz[j]) * s_0_0[j];

                t_xxzz_zzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_xxz[j] * fx[j] * fx[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 1.5 * pa_xxzz[j] * pb_z[j] * fx[j] +

                                3.0 * pa_xxz[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_zzz[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_zzz[j] + pa_xxzz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (30) = (60,62)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xy, pb_y, s_0_0, t_xyyy_xxx, t_xyyy_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xxx[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * pa_yyy[j] +

                                2.25 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 1.5 * pa_xyyy[j] * pb_x[j] * fx[j] +

                                1.5 * fx[j] * pa_yyy[j] * pb_xx[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxx[j] + pa_xyyy[j] * pb_xxx[j]) * s_0_0[j];

                t_xyyy_xxy[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xyy[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] +

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.5 * pa_xyyy[j] * fx[j] * pb_y[j] + 1.5 * pa_xyy[j] * fx[j] * pb_xx[j] +

                                fx[j] * pa_yyy[j] * pb_xy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxy[j] + pa_xyyy[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (31) = (62,64)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xz, pb_y, pb_yy, pb_z, s_0_0, t_xyyy_xxz, t_xyyy_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xxz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xyyy[j] * fx[j] * pb_z[j] + fx[j] * pa_yyy[j] * pb_xz[j] +

                                1.5 * pa_xy[j] * fx[j] * pb_xxz[j] + pa_xyyy[j] * pb_xxz[j]) * s_0_0[j];

                t_xyyy_xyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                2.25 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_yyy[j] + 1.5 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] +

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.5 * pa_xyyy[j] * pb_x[j] * fx[j] +

                                3.0 * pa_xyy[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_yyy[j] * pb_yy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyy[j] +

                                pa_xyyy[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (32) = (64,66)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xyz, pb_xz, pb_xzz, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_xyyy_xyz, t_xyyy_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 1.5 * pa_xyy[j] * fx[j] * pb_xz[j] +

                                0.5 * fx[j] * pa_yyy[j] * pb_yz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyz[j] + pa_xyyy[j] * pb_xyz[j]) * s_0_0[j];

                t_xyyy_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * pa_yyy[j] +

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.5 * pa_xyyy[j] * pb_x[j] * fx[j] +

                                0.5 * fx[j] * pa_yyy[j] * pb_zz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xzz[j] + pa_xyyy[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (33) = (66,68)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_z, s_0_0, \
                                     t_xyyy_yyy, t_xyyy_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_yyy[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_xyy[j] * fx[j] * fx[j] +

                                6.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xyyy[j] * pb_y[j] * fx[j] +

                                4.5 * pa_xyy[j] * fx[j] * pb_yy[j] + 1.5 * pa_xy[j] * fx[j] * pb_yyy[j] + pa_xyyy[j] * pb_yyy[j]) * s_0_0[j];

                t_xyyy_yyz[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyyy[j] * fx[j] * pb_z[j] + 3.0 * pa_xyy[j] * fx[j] * pb_yz[j] +

                                1.5 * pa_xy[j] * fx[j] * pb_yyz[j] + pa_xyyy[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (34) = (68,70)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pb_y, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xyyy_yzz, t_xyyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_yzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] +

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xyyy[j] * pb_y[j] * fx[j] +

                                1.5 * pa_xyy[j] * fx[j] * pb_zz[j] + 1.5 * pa_xy[j] * fx[j] * pb_yzz[j] + pa_xyyy[j] * pb_yzz[j]) * s_0_0[j];

                t_xyyy_zzz[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xyyy[j] * pb_z[j] * fx[j] +

                                1.5 * pa_xy[j] * fx[j] * pb_zzz[j] + pa_xyyy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (35) = (70,72)

            #pragma omp simd aligned(fx, pa_xyyz, pa_xyz, pa_xz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xy, pb_y, s_0_0, t_xyyz_xxx, t_xyyz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xxx[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * pa_yyz[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 1.5 * pa_xyyz[j] * pb_x[j] * fx[j] +

                                1.5 * fx[j] * pa_yyz[j] * pb_xx[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxx[j] + pa_xyyz[j] * pb_xxx[j]) * s_0_0[j];

                t_xyyz_xxy[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + fx[j] * fx[j] * pa_yz[j] * pb_x[j] +

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_y[j] +

                                pa_xyz[j] * fx[j] * pb_xx[j] + fx[j] * pa_yyz[j] * pb_xy[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxy[j] +

                                pa_xyyz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (36) = (72,74)

            #pragma omp simd aligned(fx, pa_x, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xz, pb_y, pb_yy, pb_z, s_0_0, t_xyyz_xxz, t_xyyz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xxz[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xyy[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] +

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xyy[j] * fx[j] * pb_xx[j] +

                                fx[j] * pa_yyz[j] * pb_xz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxz[j] + pa_xyyz[j] * pb_xxz[j]) * s_0_0[j];

                t_xyyz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * fx[j] * fx[j] * pa_yyz[j] + fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] +

                                0.5 * pa_xyyz[j] * pb_x[j] * fx[j] + 2.0 * pa_xyz[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_yyz[j] * pb_yy[j] +

                                0.5 * pa_xz[j] * fx[j] * pb_xyy[j] + pa_xyyz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (37) = (74,76)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, \
                                     pb_x, pb_xy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_xyyz_xyz, \
                                     t_xyyz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.5 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] +

                                0.5 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xyy[j] * fx[j] * pb_xy[j] + pa_xyz[j] * fx[j] * pb_xz[j] +

                                0.5 * fx[j] * pa_yyz[j] * pb_yz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xyz[j] + pa_xyyz[j] * pb_xyz[j]) * s_0_0[j];

                t_xyyz_xzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_yyz[j] + 0.5 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] +

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] +

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xyyz[j] * pb_x[j] * fx[j] + pa_xyy[j] * fx[j] * pb_xz[j] +

                                0.5 * fx[j] * pa_yyz[j] * pb_zz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xzz[j] + pa_xyyz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (38) = (76,78)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_xyyz_yyy, t_xyyz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_yyy[j] = (1.5 * pa_xyz[j] * fx[j] * fx[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * pa_xyyz[j] * pb_y[j] * fx[j] + 3.0 * pa_xyz[j] * fx[j] * pb_yy[j] + 0.5 * pa_xz[j] * fx[j] * pb_yyy[j] +

                                pa_xyyz[j] * pb_yyy[j]) * s_0_0[j];

                t_xyyz_yyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xyy[j] * fx[j] * fx[j] +

                                pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xyy[j] * fx[j] * pb_yy[j] +

                                2.0 * pa_xyz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xz[j] * fx[j] * pb_yyz[j] + pa_xyyz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (39) = (78,80)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, s_0_0, t_xyyz_yzz, t_xyyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_yzz[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyyz[j] * pb_y[j] * fx[j] +

                                pa_xyy[j] * fx[j] * pb_yz[j] + pa_xyz[j] * fx[j] * pb_zz[j] + 0.5 * pa_xz[j] * fx[j] * pb_yzz[j] +

                                pa_xyyz[j] * pb_yzz[j]) * s_0_0[j];

                t_xyyz_zzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xyyz[j] * pb_z[j] * fx[j] +

                                1.5 * pa_xyy[j] * fx[j] * pb_zz[j] + 0.5 * pa_xz[j] * fx[j] * pb_zzz[j] + pa_xyyz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (40) = (80,82)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyzz, pa_xzz, pa_y, pa_yzz, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xy, pb_y, s_0_0, t_xyzz_xxx, t_xyzz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xxx[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * pa_yzz[j] +

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 1.5 * pa_xyzz[j] * pb_x[j] * fx[j] +

                                1.5 * fx[j] * pa_yzz[j] * pb_xx[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxx[j] + pa_xyzz[j] * pb_xxx[j]) * s_0_0[j];

                t_xyzz_xxy[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * pa_xzz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] +

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                0.5 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_y[j] + 0.5 * pa_xzz[j] * fx[j] * pb_xx[j] +

                                fx[j] * pa_yzz[j] * pb_xy[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxy[j] + pa_xyzz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (41) = (82,84)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xzz, pa_y, pa_yz, pa_yzz, pa_zz, pb_x, pb_xx, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xz, pb_y, pb_yy, pb_z, s_0_0, t_xyzz_xxz, t_xyzz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xxz[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + fx[j] * fx[j] * pa_yz[j] * pb_x[j] +

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_z[j] +

                                pa_xyz[j] * fx[j] * pb_xx[j] + fx[j] * pa_yzz[j] * pb_xz[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxz[j] +

                                pa_xyzz[j] * pb_xxz[j]) * s_0_0[j];

                t_xyzz_xyy[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * fx[j] * fx[j] * pa_yzz[j] + 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] +

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] +

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.5 * pa_xyzz[j] * pb_x[j] * fx[j] + pa_xzz[j] * fx[j] * pb_xy[j] +

                                0.5 * fx[j] * pa_yzz[j] * pb_yy[j] + 0.5 * pa_xy[j] * fx[j] * pb_xyy[j] + pa_xyzz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (42) = (84,86)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pb_x, pb_xy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_xyzz_xyz, \
                                     t_xyzz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] +

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] +

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + pa_xyz[j] * fx[j] * pb_xy[j] + 0.5 * pa_xzz[j] * fx[j] * pb_xz[j] +

                                0.5 * fx[j] * pa_yzz[j] * pb_yz[j] + 0.5 * pa_xy[j] * fx[j] * pb_xyz[j] + pa_xyzz[j] * pb_xyz[j]) * s_0_0[j];

                t_xyzz_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * fx[j] * fx[j] * pa_yzz[j] + fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] +

                                0.5 * pa_xyzz[j] * pb_x[j] * fx[j] + 2.0 * pa_xyz[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_yzz[j] * pb_zz[j] +

                                0.5 * pa_xy[j] * fx[j] * pb_xzz[j] + pa_xyzz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (43) = (86,88)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_xyzz_yyy, t_xyzz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_yyy[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] +

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xyzz[j] * pb_y[j] * fx[j] +

                                1.5 * pa_xzz[j] * fx[j] * pb_yy[j] + 0.5 * pa_xy[j] * fx[j] * pb_yyy[j] + pa_xyzz[j] * pb_yyy[j]) * s_0_0[j];

                t_xyzz_yyz[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_z[j] +

                                pa_xyz[j] * fx[j] * pb_yy[j] + pa_xzz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xy[j] * fx[j] * pb_yyz[j] +

                                pa_xyzz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (44) = (88,90)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, s_0_0, t_xyzz_yzz, t_xyzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_yzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_xzz[j] * fx[j] * fx[j] + pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] +

                                0.5 * pa_xyzz[j] * pb_y[j] * fx[j] + 2.0 * pa_xyz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xzz[j] * fx[j] * pb_zz[j] +

                                0.5 * pa_xy[j] * fx[j] * pb_yzz[j] + pa_xyzz[j] * pb_yzz[j]) * s_0_0[j];

                t_xyzz_zzz[j] = (1.5 * pa_xyz[j] * fx[j] * fx[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_xyzz[j] * pb_z[j] * fx[j] + 3.0 * pa_xyz[j] * fx[j] * pb_zz[j] + 0.5 * pa_xy[j] * fx[j] * pb_zzz[j] +

                                pa_xyzz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (45) = (90,92)

            #pragma omp simd aligned(fx, pa_xz, pa_xzzz, pa_z, pa_zzz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xy, pb_y, s_0_0, \
                                     t_xzzz_xxx, t_xzzz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xxx[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * pa_zzz[j] +

                                2.25 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 1.5 * pa_xzzz[j] * pb_x[j] * fx[j] +

                                1.5 * fx[j] * pa_zzz[j] * pb_xx[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxx[j] + pa_xzzz[j] * pb_xxx[j]) * s_0_0[j];

                t_xzzz_xxy[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xzzz[j] * fx[j] * pb_y[j] + fx[j] * pa_zzz[j] * pb_xy[j] +

                                1.5 * pa_xz[j] * fx[j] * pb_xxy[j] + pa_xzzz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (46) = (92,94)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxz, pb_xyy, \
                                     pb_xz, pb_yy, pb_z, s_0_0, t_xzzz_xxz, t_xzzz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xxz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_xzz[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_xzzz[j] * fx[j] * pb_z[j] + 1.5 * pa_xzz[j] * fx[j] * pb_xx[j] +

                                fx[j] * pa_zzz[j] * pb_xz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxz[j] + pa_xzzz[j] * pb_xxz[j]) * s_0_0[j];

                t_xzzz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * pa_zzz[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_xzzz[j] * pb_x[j] * fx[j] +

                                0.5 * fx[j] * pa_zzz[j] * pb_yy[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyy[j] + pa_xzzz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (47) = (94,96)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_xzzz_xyz, t_xzzz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] +

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 1.5 * pa_xzz[j] * fx[j] * pb_xy[j] +

                                0.5 * fx[j] * pa_zzz[j] * pb_yz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyz[j] + pa_xzzz[j] * pb_xyz[j]) * s_0_0[j];

                t_xzzz_xzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                2.25 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_zzz[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] +

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xzzz[j] * pb_x[j] * fx[j] +

                                3.0 * pa_xzz[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_zzz[j] * pb_zz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xzz[j] +

                                pa_xzzz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (48) = (96,98)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_z, s_0_0, \
                                     t_xzzz_yyy, t_xzzz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_yyy[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xzzz[j] * pb_y[j] * fx[j] +

                                1.5 * pa_xz[j] * fx[j] * pb_yyy[j] + pa_xzzz[j] * pb_yyy[j]) * s_0_0[j];

                t_xzzz_yyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] +

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xzzz[j] * fx[j] * pb_z[j] +

                                1.5 * pa_xzz[j] * fx[j] * pb_yy[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyz[j] + pa_xzzz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (49) = (98,100)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xzzz_yzz, t_xzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_yzz[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xzzz[j] * pb_y[j] * fx[j] + 3.0 * pa_xzz[j] * fx[j] * pb_yz[j] +

                                1.5 * pa_xz[j] * fx[j] * pb_yzz[j] + pa_xzzz[j] * pb_yzz[j]) * s_0_0[j];

                t_xzzz_zzz[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_xzz[j] * fx[j] * fx[j] +

                                6.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xzzz[j] * pb_z[j] * fx[j] +

                                4.5 * pa_xzz[j] * fx[j] * pb_zz[j] + 1.5 * pa_xz[j] * fx[j] * pb_zzz[j] + pa_xzzz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (50) = (100,102)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_y, s_0_0, \
                                     t_yyyy_xxx, t_yyyy_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xxx[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                1.5 * pa_yyyy[j] * pb_x[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxx[j] +

                                pa_yyyy[j] * pb_xxx[j]) * s_0_0[j];

                t_yyyy_xxy[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * fx[j] + pa_yyy[j] * fx[j] * fx[j] +

                                0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] +

                                0.5 * pa_yyyy[j] * fx[j] * pb_y[j] + 2.0 * pa_yyy[j] * fx[j] * pb_xx[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] +

                                3.0 * pa_yy[j] * fx[j] * pb_xxy[j] + pa_yyyy[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (51) = (102,104)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xxz, pb_xy, pb_xyy, pb_z, s_0_0, \
                                     t_yyyy_xxz, t_yyyy_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_yyyy[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxz[j] +

                                pa_yyyy[j] * pb_xxz[j]) * s_0_0[j];

                t_yyyy_xyy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                6.0 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_yyyy[j] * pb_x[j] * fx[j] + 4.0 * pa_yyy[j] * fx[j] * pb_xy[j] +

                                0.75 * fx[j] * fx[j] * pb_xyy[j] + 3.0 * pa_yy[j] * fx[j] * pb_xyy[j] + pa_yyyy[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (52) = (104,106)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xyz, pb_xz, pb_xzz, s_0_0, t_yyyy_xyz, \
                                     t_yyyy_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xyz[j] = (3.0 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 2.0 * pa_yyy[j] * fx[j] * pb_xz[j] +

                                0.75 * fx[j] * fx[j] * pb_xyz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xyz[j] + pa_yyyy[j] * pb_xyz[j]) * s_0_0[j];

                t_yyyy_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.5 * pa_yyyy[j] * pb_x[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xzz[j] +

                                pa_yyyy[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (53) = (106,108)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_z, s_0_0, \
                                     t_yyyy_yyy, t_yyyy_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_yyy[j] = (7.5 * pa_y[j] * fx[j] * fx[j] * fx[j] + 5.625 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                3.0 * pa_yyy[j] * fx[j] * fx[j] + 13.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] +

                                9.0 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_yyyy[j] * pb_y[j] * fx[j] + 6.0 * pa_yyy[j] * fx[j] * pb_yy[j] +

                                0.75 * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yy[j] * fx[j] * pb_yyy[j] + pa_yyyy[j] * pb_yyy[j]) * s_0_0[j];

                t_yyyy_yyz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                6.0 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_yyyy[j] * fx[j] * pb_z[j] + 4.0 * pa_yyy[j] * fx[j] * pb_yz[j] +

                                0.75 * fx[j] * fx[j] * pb_yyz[j] + 3.0 * pa_yy[j] * fx[j] * pb_yyz[j] + pa_yyyy[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (54) = (108,110)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_y, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_yyyy_yzz, t_yyyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_yzz[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * fx[j] + pa_yyy[j] * fx[j] * fx[j] +

                                0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] +

                                0.5 * pa_yyyy[j] * pb_y[j] * fx[j] + 2.0 * pa_yyy[j] * fx[j] * pb_zz[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] +

                                3.0 * pa_yy[j] * fx[j] * pb_yzz[j] + pa_yyyy[j] * pb_yzz[j]) * s_0_0[j];

                t_yyyy_zzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_yyyy[j] * pb_z[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_yy[j] * fx[j] * pb_zzz[j] +

                                pa_yyyy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (55) = (110,112)

            #pragma omp simd aligned(fx, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_y, s_0_0, \
                                     t_yyyz_xxx, t_yyyz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xxx[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yyyz[j] * pb_x[j] * fx[j] +

                                1.5 * pa_yz[j] * fx[j] * pb_xxx[j] + pa_yyyz[j] * pb_xxx[j]) * s_0_0[j];

                t_yyyz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_yyz[j] * fx[j] * fx[j] +

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_yyyz[j] * fx[j] * pb_y[j] +

                                1.5 * pa_yyz[j] * fx[j] * pb_xx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxy[j] + pa_yyyz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (56) = (112,114)

            #pragma omp simd aligned(fx, pa_y, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, pb_xy, pb_xyy, \
                                     pb_z, s_0_0, t_yyyz_xxz, t_yyyz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xxz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_yyy[j] * fx[j] * fx[j] +

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_yyyz[j] * fx[j] * pb_z[j] +

                                0.5 * pa_yyy[j] * fx[j] * pb_xx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxz[j] + pa_yyyz[j] * pb_xxz[j]) * s_0_0[j];

                t_yyyz_xyy[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_yyyz[j] * pb_x[j] * fx[j] + 3.0 * pa_yyz[j] * fx[j] * pb_xy[j] +

                                1.5 * pa_yz[j] * fx[j] * pb_xyy[j] + pa_yyyz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (57) = (114,116)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xy, pb_xyz, pb_xz, \
                                     pb_xzz, s_0_0, t_yyyz_xyz, t_yyyz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_yyy[j] * fx[j] * pb_xy[j] +

                                1.5 * pa_yyz[j] * fx[j] * pb_xz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xyz[j] + pa_yyyz[j] * pb_xyz[j]) * s_0_0[j];

                t_yyyz_xzz[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] +

                                1.5 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_yyyz[j] * pb_x[j] * fx[j] + pa_yyy[j] * fx[j] * pb_xz[j] +

                                1.5 * pa_yz[j] * fx[j] * pb_xzz[j] + pa_yyyz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (58) = (116,118)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_yyyz_yyy, t_yyyz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_yyy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] + 2.25 * pa_yyz[j] * fx[j] * fx[j] +

                                6.75 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 1.5 * pa_yyyz[j] * pb_y[j] * fx[j] +

                                4.5 * pa_yyz[j] * fx[j] * pb_yy[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyy[j] + pa_yyyz[j] * pb_yyy[j]) * s_0_0[j];

                t_yyyz_yyz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.25 * pa_yyy[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] +

                                2.25 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_yyyz[j] * fx[j] * pb_z[j] + 0.5 * pa_yyy[j] * fx[j] * pb_yy[j] +

                                3.0 * pa_yyz[j] * fx[j] * pb_yz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyz[j] + pa_yyyz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (59) = (118,120)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_y, pb_yz, pb_yzz, pb_z, \
                                     pb_zz, pb_zzz, s_0_0, t_yyyz_yzz, t_yyyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_yzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_yyz[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] +

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_yyyz[j] * pb_y[j] * fx[j] + pa_yyy[j] * fx[j] * pb_yz[j] +

                                1.5 * pa_yyz[j] * fx[j] * pb_zz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yzz[j] + pa_yyyz[j] * pb_yzz[j]) * s_0_0[j];

                t_yyyz_zzz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yyy[j] * fx[j] * fx[j] +

                                2.25 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_yyyz[j] * pb_z[j] * fx[j] +

                                1.5 * pa_yyy[j] * fx[j] * pb_zz[j] + 1.5 * pa_yz[j] * fx[j] * pb_zzz[j] + pa_yyyz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (60) = (120,122)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyzz, pa_yzz, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_y, s_0_0, \
                                     t_yyzz_xxx, t_yyzz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xxx[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 1.5 * pa_yyzz[j] * pb_x[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] +

                                0.5 * pa_yy[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxx[j] + pa_yyzz[j] * pb_xxx[j]) * s_0_0[j];

                t_yyzz_xxy[j] = (0.25 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.5 * pa_yzz[j] * fx[j] * fx[j] +

                                0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] +

                                0.5 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 0.5 * pa_yyzz[j] * fx[j] * pb_y[j] +

                                pa_yzz[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_yy[j] * fx[j] * pb_xxy[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_xxy[j] + pa_yyzz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (61) = (122,124)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xy, \
                                     pb_xyy, pb_z, s_0_0, t_yyzz_xxz, t_yyzz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xxz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_yyz[j] * fx[j] * fx[j] +

                                0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_yyzz[j] * fx[j] * pb_z[j] +

                                pa_yyz[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xxz[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_xxz[j] + pa_yyzz[j] * pb_xxz[j]) * s_0_0[j];

                t_yyzz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] +

                                0.25 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_yyzz[j] * pb_x[j] * fx[j] +

                                2.0 * pa_yzz[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_yy[j] * fx[j] * pb_xyy[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_xyy[j] + pa_yyzz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (62) = (124,126)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xy, pb_xyz, \
                                     pb_xz, pb_xzz, s_0_0, t_yyzz_xyz, t_yyzz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xyz[j] = (pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] +

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + pa_yyz[j] * fx[j] * pb_xy[j] + pa_yzz[j] * fx[j] * pb_xz[j] +

                                0.25 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xyz[j] +

                                pa_yyzz[j] * pb_xyz[j]) * s_0_0[j];

                t_yyzz_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] +

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_yyzz[j] * pb_x[j] * fx[j] +

                                2.0 * pa_yyz[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xzz[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_xzz[j] + pa_yyzz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (63) = (126,128)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_z, s_0_0, t_yyzz_yyy, t_yyzz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_yyy[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * pa_yzz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] +

                                0.75 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_yyzz[j] * pb_y[j] * fx[j] +

                                3.0 * pa_yzz[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 0.5 * pa_yy[j] * fx[j] * pb_yyy[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_yyy[j] + pa_yyzz[j] * pb_yyy[j]) * s_0_0[j];

                t_yyzz_yyz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_yyz[j] * fx[j] * fx[j] + 2.0 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + pa_y[j] * fx[j] * fx[j] * pb_yz[j] +

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_yyzz[j] * fx[j] * pb_z[j] + pa_yyz[j] * fx[j] * pb_yy[j] +

                                2.0 * pa_yzz[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_yy[j] * fx[j] * pb_yyz[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_yyz[j] + pa_yyzz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (64) = (128,130)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_y, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, s_0_0, t_yyzz_yzz, t_yyzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_yzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_yzz[j] * fx[j] * fx[j] + 2.0 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] +

                                0.5 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + fx[j] * fx[j] * pa_z[j] * pb_yz[j] +

                                0.5 * pa_yyzz[j] * pb_y[j] * fx[j] + 2.0 * pa_yyz[j] * fx[j] * pb_yz[j] + pa_yzz[j] * fx[j] * pb_zz[j] +

                                0.25 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yzz[j] +

                                pa_yyzz[j] * pb_yzz[j]) * s_0_0[j];

                t_yyzz_zzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                1.5 * pa_yyz[j] * fx[j] * fx[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] +

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 1.5 * pa_yyzz[j] * pb_z[j] * fx[j] +

                                3.0 * pa_yyz[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_zzz[j] +

                                0.5 * fx[j] * pa_zz[j] * pb_zzz[j] + pa_yyzz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (65) = (130,132)

            #pragma omp simd aligned(fx, pa_yz, pa_yzzz, pa_z, pa_zzz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_y, s_0_0, \
                                     t_yzzz_xxx, t_yzzz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xxx[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yzzz[j] * pb_x[j] * fx[j] +

                                1.5 * pa_yz[j] * fx[j] * pb_xxx[j] + pa_yzzz[j] * pb_xxx[j]) * s_0_0[j];

                t_yzzz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * pa_zzz[j] +

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_yzzz[j] * fx[j] * pb_y[j] +

                                0.5 * fx[j] * pa_zzz[j] * pb_xx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxy[j] + pa_yzzz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (66) = (132,134)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zzz, pb_x, pb_xx, pb_xxz, pb_xy, pb_xyy, \
                                     pb_z, s_0_0, t_yzzz_xxz, t_yzzz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xxz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yzz[j] * fx[j] * fx[j] +

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_yzzz[j] * fx[j] * pb_z[j] +

                                1.5 * pa_yzz[j] * fx[j] * pb_xx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxz[j] + pa_yzzz[j] * pb_xxz[j]) * s_0_0[j];

                t_yzzz_xyy[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_yzzz[j] * pb_x[j] * fx[j] + fx[j] * pa_zzz[j] * pb_xy[j] +

                                1.5 * pa_yz[j] * fx[j] * pb_xyy[j] + pa_yzzz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (67) = (134,136)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xy, pb_xyz, pb_xz, \
                                     pb_xzz, s_0_0, t_yzzz_xyz, t_yzzz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] +

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 1.5 * pa_yzz[j] * fx[j] * pb_xy[j] +

                                0.5 * fx[j] * pa_zzz[j] * pb_xz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xyz[j] + pa_yzzz[j] * pb_xyz[j]) * s_0_0[j];

                t_yzzz_xzz[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] +

                                1.5 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_yzzz[j] * pb_x[j] * fx[j] + 3.0 * pa_yzz[j] * fx[j] * pb_xz[j] +

                                1.5 * pa_yz[j] * fx[j] * pb_xzz[j] + pa_yzzz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (68) = (136,138)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_yzzz_yyy, t_yzzz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_yyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * pa_zzz[j] +

                                2.25 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 1.5 * pa_yzzz[j] * pb_y[j] * fx[j] +

                                1.5 * fx[j] * pa_zzz[j] * pb_yy[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyy[j] + pa_yzzz[j] * pb_yyy[j]) * s_0_0[j];

                t_yzzz_yyz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                0.75 * pa_yzz[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] +

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] +

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_yzzz[j] * fx[j] * pb_z[j] + 1.5 * pa_yzz[j] * fx[j] * pb_yy[j] +

                                fx[j] * pa_zzz[j] * pb_yz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyz[j] + pa_yzzz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (69) = (138,140)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_y, pb_yz, pb_yzz, pb_z, \
                                     pb_zz, pb_zzz, s_0_0, t_yzzz_yzz, t_yzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_yzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                2.25 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pa_zzz[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] +

                                1.5 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_yzzz[j] * pb_y[j] * fx[j] +

                                3.0 * pa_yzz[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_zzz[j] * pb_zz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yzz[j] +

                                pa_yzzz[j] * pb_yzz[j]) * s_0_0[j];

                t_yzzz_zzz[j] = (1.875 * pa_y[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_yzz[j] * fx[j] * fx[j] +

                                6.75 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_yzzz[j] * pb_z[j] * fx[j] +

                                4.5 * pa_yzz[j] * fx[j] * pb_zz[j] + 1.5 * pa_yz[j] * fx[j] * pb_zzz[j] + pa_yzzz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (70) = (140,142)

            #pragma omp simd aligned(fx, pa_zz, pa_zzzz, pb_x, pb_xxx, pb_xxy, pb_y, s_0_0, t_zzzz_xxx, t_zzzz_xxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xxx[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] +

                                1.5 * pa_zzzz[j] * pb_x[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxx[j] +

                                pa_zzzz[j] * pb_xxx[j]) * s_0_0[j];

                t_zzzz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] +

                                0.5 * pa_zzzz[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxy[j] +

                                pa_zzzz[j] * pb_xxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (71) = (142,144)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xx, pb_xxz, pb_xyy, pb_z, s_0_0, \
                                     t_zzzz_xxz, t_zzzz_xyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xxz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * fx[j] + pa_zzz[j] * fx[j] * fx[j] +

                                0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_z[j] * fx[j] * fx[j] * pb_xx[j] +

                                0.5 * pa_zzzz[j] * fx[j] * pb_z[j] + 2.0 * pa_zzz[j] * fx[j] * pb_xx[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] +

                                3.0 * pa_zz[j] * fx[j] * pb_xxz[j] + pa_zzzz[j] * pb_xxz[j]) * s_0_0[j];

                t_zzzz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] +

                                0.5 * pa_zzzz[j] * pb_x[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] + 3.0 * pa_zz[j] * fx[j] * pb_xyy[j] +

                                pa_zzzz[j] * pb_xyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (72) = (144,146)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xy, pb_xyz, pb_xz, pb_xzz, s_0_0, \
                                     t_zzzz_xyz, t_zzzz_xzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xyz[j] = (3.0 * pa_z[j] * fx[j] * fx[j] * pb_xy[j] + 2.0 * pa_zzz[j] * fx[j] * pb_xy[j] +

                                0.75 * fx[j] * fx[j] * pb_xyz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xyz[j] + pa_zzzz[j] * pb_xyz[j]) * s_0_0[j];

                t_zzzz_xzz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] +

                                6.0 * pa_z[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_zzzz[j] * pb_x[j] * fx[j] + 4.0 * pa_zzz[j] * fx[j] * pb_xz[j] +

                                0.75 * fx[j] * fx[j] * pb_xzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xzz[j] + pa_zzzz[j] * pb_xzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (73) = (146,148)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_z, s_0_0, \
                                     t_zzzz_yyy, t_zzzz_yyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_yyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] +

                                1.5 * pa_zzzz[j] * pb_y[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_zz[j] * fx[j] * pb_yyy[j] +

                                pa_zzzz[j] * pb_yyy[j]) * s_0_0[j];

                t_zzzz_yyz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * fx[j] + pa_zzz[j] * fx[j] * fx[j] +

                                0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_z[j] * fx[j] * fx[j] * pb_yy[j] +

                                0.5 * pa_zzzz[j] * fx[j] * pb_z[j] + 2.0 * pa_zzz[j] * fx[j] * pb_yy[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] +

                                3.0 * pa_zz[j] * fx[j] * pb_yyz[j] + pa_zzzz[j] * pb_yyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (74) = (148,150)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_zzzz_yzz, t_zzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_yzz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] +

                                6.0 * pa_z[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_zzzz[j] * pb_y[j] * fx[j] + 4.0 * pa_zzz[j] * fx[j] * pb_yz[j] +

                                0.75 * fx[j] * fx[j] * pb_yzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_yzz[j] + pa_zzzz[j] * pb_yzz[j]) * s_0_0[j];

                t_zzzz_zzz[j] = (7.5 * pa_z[j] * fx[j] * fx[j] * fx[j] + 5.625 * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                3.0 * pa_zzz[j] * fx[j] * fx[j] + 13.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] +

                                9.0 * pa_z[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_zzzz[j] * pb_z[j] * fx[j] + 6.0 * pa_zzz[j] * fx[j] * pb_zz[j] +

                                0.75 * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_zzz[j] + pa_zzzz[j] * pb_zzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    void
    compOverlapForGG(      CMemBlock2D<double>& primBuffer,
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& osFactors,
                     const CMemBlock2D<double>& paDistances,
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

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(2 * idx);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xxx = pbDistances.data(34 * idx + 9);

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_xxxx_xxxx = primBuffer.data(225 * idx);

            auto t_xxxx_xxxy = primBuffer.data(225 * idx + 1);

            auto t_xxxx_xxxz = primBuffer.data(225 * idx + 2);

            auto t_xxxx_xxyy = primBuffer.data(225 * idx + 3);

            auto t_xxxx_xxyz = primBuffer.data(225 * idx + 4);

            auto t_xxxx_xxzz = primBuffer.data(225 * idx + 5);

            auto t_xxxx_xyyy = primBuffer.data(225 * idx + 6);

            auto t_xxxx_xyyz = primBuffer.data(225 * idx + 7);

            auto t_xxxx_xyzz = primBuffer.data(225 * idx + 8);

            auto t_xxxx_xzzz = primBuffer.data(225 * idx + 9);

            auto t_xxxx_yyyy = primBuffer.data(225 * idx + 10);

            auto t_xxxx_yyyz = primBuffer.data(225 * idx + 11);

            auto t_xxxx_yyzz = primBuffer.data(225 * idx + 12);

            auto t_xxxx_yzzz = primBuffer.data(225 * idx + 13);

            auto t_xxxx_zzzz = primBuffer.data(225 * idx + 14);

            auto t_xxxy_xxxx = primBuffer.data(225 * idx + 15);

            auto t_xxxy_xxxy = primBuffer.data(225 * idx + 16);

            auto t_xxxy_xxxz = primBuffer.data(225 * idx + 17);

            auto t_xxxy_xxyy = primBuffer.data(225 * idx + 18);

            auto t_xxxy_xxyz = primBuffer.data(225 * idx + 19);

            auto t_xxxy_xxzz = primBuffer.data(225 * idx + 20);

            auto t_xxxy_xyyy = primBuffer.data(225 * idx + 21);

            auto t_xxxy_xyyz = primBuffer.data(225 * idx + 22);

            auto t_xxxy_xyzz = primBuffer.data(225 * idx + 23);

            auto t_xxxy_xzzz = primBuffer.data(225 * idx + 24);

            auto t_xxxy_yyyy = primBuffer.data(225 * idx + 25);

            auto t_xxxy_yyyz = primBuffer.data(225 * idx + 26);

            auto t_xxxy_yyzz = primBuffer.data(225 * idx + 27);

            auto t_xxxy_yzzz = primBuffer.data(225 * idx + 28);

            auto t_xxxy_zzzz = primBuffer.data(225 * idx + 29);

            auto t_xxxz_xxxx = primBuffer.data(225 * idx + 30);

            auto t_xxxz_xxxy = primBuffer.data(225 * idx + 31);

            auto t_xxxz_xxxz = primBuffer.data(225 * idx + 32);

            auto t_xxxz_xxyy = primBuffer.data(225 * idx + 33);

            auto t_xxxz_xxyz = primBuffer.data(225 * idx + 34);

            auto t_xxxz_xxzz = primBuffer.data(225 * idx + 35);

            auto t_xxxz_xyyy = primBuffer.data(225 * idx + 36);

            auto t_xxxz_xyyz = primBuffer.data(225 * idx + 37);

            auto t_xxxz_xyzz = primBuffer.data(225 * idx + 38);

            auto t_xxxz_xzzz = primBuffer.data(225 * idx + 39);

            auto t_xxxz_yyyy = primBuffer.data(225 * idx + 40);

            auto t_xxxz_yyyz = primBuffer.data(225 * idx + 41);

            auto t_xxxz_yyzz = primBuffer.data(225 * idx + 42);

            auto t_xxxz_yzzz = primBuffer.data(225 * idx + 43);

            auto t_xxxz_zzzz = primBuffer.data(225 * idx + 44);

            auto t_xxyy_xxxx = primBuffer.data(225 * idx + 45);

            auto t_xxyy_xxxy = primBuffer.data(225 * idx + 46);

            auto t_xxyy_xxxz = primBuffer.data(225 * idx + 47);

            auto t_xxyy_xxyy = primBuffer.data(225 * idx + 48);

            auto t_xxyy_xxyz = primBuffer.data(225 * idx + 49);

            auto t_xxyy_xxzz = primBuffer.data(225 * idx + 50);

            auto t_xxyy_xyyy = primBuffer.data(225 * idx + 51);

            auto t_xxyy_xyyz = primBuffer.data(225 * idx + 52);

            auto t_xxyy_xyzz = primBuffer.data(225 * idx + 53);

            auto t_xxyy_xzzz = primBuffer.data(225 * idx + 54);

            auto t_xxyy_yyyy = primBuffer.data(225 * idx + 55);

            auto t_xxyy_yyyz = primBuffer.data(225 * idx + 56);

            auto t_xxyy_yyzz = primBuffer.data(225 * idx + 57);

            auto t_xxyy_yzzz = primBuffer.data(225 * idx + 58);

            auto t_xxyy_zzzz = primBuffer.data(225 * idx + 59);

            auto t_xxyz_xxxx = primBuffer.data(225 * idx + 60);

            auto t_xxyz_xxxy = primBuffer.data(225 * idx + 61);

            auto t_xxyz_xxxz = primBuffer.data(225 * idx + 62);

            auto t_xxyz_xxyy = primBuffer.data(225 * idx + 63);

            auto t_xxyz_xxyz = primBuffer.data(225 * idx + 64);

            auto t_xxyz_xxzz = primBuffer.data(225 * idx + 65);

            auto t_xxyz_xyyy = primBuffer.data(225 * idx + 66);

            auto t_xxyz_xyyz = primBuffer.data(225 * idx + 67);

            auto t_xxyz_xyzz = primBuffer.data(225 * idx + 68);

            auto t_xxyz_xzzz = primBuffer.data(225 * idx + 69);

            auto t_xxyz_yyyy = primBuffer.data(225 * idx + 70);

            auto t_xxyz_yyyz = primBuffer.data(225 * idx + 71);

            auto t_xxyz_yyzz = primBuffer.data(225 * idx + 72);

            auto t_xxyz_yzzz = primBuffer.data(225 * idx + 73);

            auto t_xxyz_zzzz = primBuffer.data(225 * idx + 74);

            auto t_xxzz_xxxx = primBuffer.data(225 * idx + 75);

            auto t_xxzz_xxxy = primBuffer.data(225 * idx + 76);

            auto t_xxzz_xxxz = primBuffer.data(225 * idx + 77);

            auto t_xxzz_xxyy = primBuffer.data(225 * idx + 78);

            auto t_xxzz_xxyz = primBuffer.data(225 * idx + 79);

            auto t_xxzz_xxzz = primBuffer.data(225 * idx + 80);

            auto t_xxzz_xyyy = primBuffer.data(225 * idx + 81);

            auto t_xxzz_xyyz = primBuffer.data(225 * idx + 82);

            auto t_xxzz_xyzz = primBuffer.data(225 * idx + 83);

            auto t_xxzz_xzzz = primBuffer.data(225 * idx + 84);

            auto t_xxzz_yyyy = primBuffer.data(225 * idx + 85);

            auto t_xxzz_yyyz = primBuffer.data(225 * idx + 86);

            auto t_xxzz_yyzz = primBuffer.data(225 * idx + 87);

            auto t_xxzz_yzzz = primBuffer.data(225 * idx + 88);

            auto t_xxzz_zzzz = primBuffer.data(225 * idx + 89);

            auto t_xyyy_xxxx = primBuffer.data(225 * idx + 90);

            auto t_xyyy_xxxy = primBuffer.data(225 * idx + 91);

            auto t_xyyy_xxxz = primBuffer.data(225 * idx + 92);

            auto t_xyyy_xxyy = primBuffer.data(225 * idx + 93);

            auto t_xyyy_xxyz = primBuffer.data(225 * idx + 94);

            auto t_xyyy_xxzz = primBuffer.data(225 * idx + 95);

            auto t_xyyy_xyyy = primBuffer.data(225 * idx + 96);

            auto t_xyyy_xyyz = primBuffer.data(225 * idx + 97);

            auto t_xyyy_xyzz = primBuffer.data(225 * idx + 98);

            auto t_xyyy_xzzz = primBuffer.data(225 * idx + 99);

            auto t_xyyy_yyyy = primBuffer.data(225 * idx + 100);

            auto t_xyyy_yyyz = primBuffer.data(225 * idx + 101);

            auto t_xyyy_yyzz = primBuffer.data(225 * idx + 102);

            auto t_xyyy_yzzz = primBuffer.data(225 * idx + 103);

            auto t_xyyy_zzzz = primBuffer.data(225 * idx + 104);

            auto t_xyyz_xxxx = primBuffer.data(225 * idx + 105);

            auto t_xyyz_xxxy = primBuffer.data(225 * idx + 106);

            auto t_xyyz_xxxz = primBuffer.data(225 * idx + 107);

            auto t_xyyz_xxyy = primBuffer.data(225 * idx + 108);

            auto t_xyyz_xxyz = primBuffer.data(225 * idx + 109);

            auto t_xyyz_xxzz = primBuffer.data(225 * idx + 110);

            auto t_xyyz_xyyy = primBuffer.data(225 * idx + 111);

            auto t_xyyz_xyyz = primBuffer.data(225 * idx + 112);

            auto t_xyyz_xyzz = primBuffer.data(225 * idx + 113);

            auto t_xyyz_xzzz = primBuffer.data(225 * idx + 114);

            auto t_xyyz_yyyy = primBuffer.data(225 * idx + 115);

            auto t_xyyz_yyyz = primBuffer.data(225 * idx + 116);

            auto t_xyyz_yyzz = primBuffer.data(225 * idx + 117);

            auto t_xyyz_yzzz = primBuffer.data(225 * idx + 118);

            auto t_xyyz_zzzz = primBuffer.data(225 * idx + 119);

            auto t_xyzz_xxxx = primBuffer.data(225 * idx + 120);

            auto t_xyzz_xxxy = primBuffer.data(225 * idx + 121);

            auto t_xyzz_xxxz = primBuffer.data(225 * idx + 122);

            auto t_xyzz_xxyy = primBuffer.data(225 * idx + 123);

            auto t_xyzz_xxyz = primBuffer.data(225 * idx + 124);

            auto t_xyzz_xxzz = primBuffer.data(225 * idx + 125);

            auto t_xyzz_xyyy = primBuffer.data(225 * idx + 126);

            auto t_xyzz_xyyz = primBuffer.data(225 * idx + 127);

            auto t_xyzz_xyzz = primBuffer.data(225 * idx + 128);

            auto t_xyzz_xzzz = primBuffer.data(225 * idx + 129);

            auto t_xyzz_yyyy = primBuffer.data(225 * idx + 130);

            auto t_xyzz_yyyz = primBuffer.data(225 * idx + 131);

            auto t_xyzz_yyzz = primBuffer.data(225 * idx + 132);

            auto t_xyzz_yzzz = primBuffer.data(225 * idx + 133);

            auto t_xyzz_zzzz = primBuffer.data(225 * idx + 134);

            auto t_xzzz_xxxx = primBuffer.data(225 * idx + 135);

            auto t_xzzz_xxxy = primBuffer.data(225 * idx + 136);

            auto t_xzzz_xxxz = primBuffer.data(225 * idx + 137);

            auto t_xzzz_xxyy = primBuffer.data(225 * idx + 138);

            auto t_xzzz_xxyz = primBuffer.data(225 * idx + 139);

            auto t_xzzz_xxzz = primBuffer.data(225 * idx + 140);

            auto t_xzzz_xyyy = primBuffer.data(225 * idx + 141);

            auto t_xzzz_xyyz = primBuffer.data(225 * idx + 142);

            auto t_xzzz_xyzz = primBuffer.data(225 * idx + 143);

            auto t_xzzz_xzzz = primBuffer.data(225 * idx + 144);

            auto t_xzzz_yyyy = primBuffer.data(225 * idx + 145);

            auto t_xzzz_yyyz = primBuffer.data(225 * idx + 146);

            auto t_xzzz_yyzz = primBuffer.data(225 * idx + 147);

            auto t_xzzz_yzzz = primBuffer.data(225 * idx + 148);

            auto t_xzzz_zzzz = primBuffer.data(225 * idx + 149);

            auto t_yyyy_xxxx = primBuffer.data(225 * idx + 150);

            auto t_yyyy_xxxy = primBuffer.data(225 * idx + 151);

            auto t_yyyy_xxxz = primBuffer.data(225 * idx + 152);

            auto t_yyyy_xxyy = primBuffer.data(225 * idx + 153);

            auto t_yyyy_xxyz = primBuffer.data(225 * idx + 154);

            auto t_yyyy_xxzz = primBuffer.data(225 * idx + 155);

            auto t_yyyy_xyyy = primBuffer.data(225 * idx + 156);

            auto t_yyyy_xyyz = primBuffer.data(225 * idx + 157);

            auto t_yyyy_xyzz = primBuffer.data(225 * idx + 158);

            auto t_yyyy_xzzz = primBuffer.data(225 * idx + 159);

            auto t_yyyy_yyyy = primBuffer.data(225 * idx + 160);

            auto t_yyyy_yyyz = primBuffer.data(225 * idx + 161);

            auto t_yyyy_yyzz = primBuffer.data(225 * idx + 162);

            auto t_yyyy_yzzz = primBuffer.data(225 * idx + 163);

            auto t_yyyy_zzzz = primBuffer.data(225 * idx + 164);

            auto t_yyyz_xxxx = primBuffer.data(225 * idx + 165);

            auto t_yyyz_xxxy = primBuffer.data(225 * idx + 166);

            auto t_yyyz_xxxz = primBuffer.data(225 * idx + 167);

            auto t_yyyz_xxyy = primBuffer.data(225 * idx + 168);

            auto t_yyyz_xxyz = primBuffer.data(225 * idx + 169);

            auto t_yyyz_xxzz = primBuffer.data(225 * idx + 170);

            auto t_yyyz_xyyy = primBuffer.data(225 * idx + 171);

            auto t_yyyz_xyyz = primBuffer.data(225 * idx + 172);

            auto t_yyyz_xyzz = primBuffer.data(225 * idx + 173);

            auto t_yyyz_xzzz = primBuffer.data(225 * idx + 174);

            auto t_yyyz_yyyy = primBuffer.data(225 * idx + 175);

            auto t_yyyz_yyyz = primBuffer.data(225 * idx + 176);

            auto t_yyyz_yyzz = primBuffer.data(225 * idx + 177);

            auto t_yyyz_yzzz = primBuffer.data(225 * idx + 178);

            auto t_yyyz_zzzz = primBuffer.data(225 * idx + 179);

            auto t_yyzz_xxxx = primBuffer.data(225 * idx + 180);

            auto t_yyzz_xxxy = primBuffer.data(225 * idx + 181);

            auto t_yyzz_xxxz = primBuffer.data(225 * idx + 182);

            auto t_yyzz_xxyy = primBuffer.data(225 * idx + 183);

            auto t_yyzz_xxyz = primBuffer.data(225 * idx + 184);

            auto t_yyzz_xxzz = primBuffer.data(225 * idx + 185);

            auto t_yyzz_xyyy = primBuffer.data(225 * idx + 186);

            auto t_yyzz_xyyz = primBuffer.data(225 * idx + 187);

            auto t_yyzz_xyzz = primBuffer.data(225 * idx + 188);

            auto t_yyzz_xzzz = primBuffer.data(225 * idx + 189);

            auto t_yyzz_yyyy = primBuffer.data(225 * idx + 190);

            auto t_yyzz_yyyz = primBuffer.data(225 * idx + 191);

            auto t_yyzz_yyzz = primBuffer.data(225 * idx + 192);

            auto t_yyzz_yzzz = primBuffer.data(225 * idx + 193);

            auto t_yyzz_zzzz = primBuffer.data(225 * idx + 194);

            auto t_yzzz_xxxx = primBuffer.data(225 * idx + 195);

            auto t_yzzz_xxxy = primBuffer.data(225 * idx + 196);

            auto t_yzzz_xxxz = primBuffer.data(225 * idx + 197);

            auto t_yzzz_xxyy = primBuffer.data(225 * idx + 198);

            auto t_yzzz_xxyz = primBuffer.data(225 * idx + 199);

            auto t_yzzz_xxzz = primBuffer.data(225 * idx + 200);

            auto t_yzzz_xyyy = primBuffer.data(225 * idx + 201);

            auto t_yzzz_xyyz = primBuffer.data(225 * idx + 202);

            auto t_yzzz_xyzz = primBuffer.data(225 * idx + 203);

            auto t_yzzz_xzzz = primBuffer.data(225 * idx + 204);

            auto t_yzzz_yyyy = primBuffer.data(225 * idx + 205);

            auto t_yzzz_yyyz = primBuffer.data(225 * idx + 206);

            auto t_yzzz_yyzz = primBuffer.data(225 * idx + 207);

            auto t_yzzz_yzzz = primBuffer.data(225 * idx + 208);

            auto t_yzzz_zzzz = primBuffer.data(225 * idx + 209);

            auto t_zzzz_xxxx = primBuffer.data(225 * idx + 210);

            auto t_zzzz_xxxy = primBuffer.data(225 * idx + 211);

            auto t_zzzz_xxxz = primBuffer.data(225 * idx + 212);

            auto t_zzzz_xxyy = primBuffer.data(225 * idx + 213);

            auto t_zzzz_xxyz = primBuffer.data(225 * idx + 214);

            auto t_zzzz_xxzz = primBuffer.data(225 * idx + 215);

            auto t_zzzz_xyyy = primBuffer.data(225 * idx + 216);

            auto t_zzzz_xyyz = primBuffer.data(225 * idx + 217);

            auto t_zzzz_xyzz = primBuffer.data(225 * idx + 218);

            auto t_zzzz_xzzz = primBuffer.data(225 * idx + 219);

            auto t_zzzz_yyyy = primBuffer.data(225 * idx + 220);

            auto t_zzzz_yyyz = primBuffer.data(225 * idx + 221);

            auto t_zzzz_yyzz = primBuffer.data(225 * idx + 222);

            auto t_zzzz_yzzz = primBuffer.data(225 * idx + 223);

            auto t_zzzz_zzzz = primBuffer.data(225 * idx + 224);

            // Batch of Integrals (0) = (0,2)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxy, \
                                     pb_xy, pb_y, s_0_0, t_xxxx_xxxx, t_xxxx_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xxxx[j] = (6.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 11.25 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 30.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 11.25 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_xxxx[j] * fx[j] * fx[j] +

                                 12.0 * pa_xxx[j] * fx[j] * fx[j] * pb_x[j] + 27.0 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 12.0 * pa_x[j] * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xxxx[j] * pb_xx[j] * fx[j] + 8.0 * pa_xxx[j] * fx[j] * pb_xxx[j] +

                                 0.75 * fx[j] * fx[j] * pb_xxxx[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxxx[j] + pa_xxxx[j] * pb_xxxx[j]) * s_0_0[j];

                t_xxxx_xxxy[j] = (7.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 5.625 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 3.0 * pa_xxx[j] * fx[j] * fx[j] * pb_y[j] +

                                 13.5 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] + 9.0 * pa_x[j] * fx[j] * fx[j] * pb_xxy[j] + 1.5 * pa_xxxx[j] * pb_xy[j] * fx[j] +

                                 6.0 * pa_xxx[j] * fx[j] * pb_xxy[j] + 0.75 * fx[j] * fx[j] * pb_xxxy[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxxy[j] +

                                 pa_xxxx[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (2,4)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xx, pb_xxxz, pb_xxyy, pb_xxz, pb_xyy, \
                                     pb_xz, pb_yy, pb_z, s_0_0, t_xxxx_xxxz, t_xxxx_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xxxz[j] = (7.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 5.625 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 3.0 * pa_xxx[j] * fx[j] * fx[j] * pb_z[j] +

                                 13.5 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] + 9.0 * pa_x[j] * fx[j] * fx[j] * pb_xxz[j] + 1.5 * pa_xxxx[j] * pb_xz[j] * fx[j] +

                                 6.0 * pa_xxx[j] * fx[j] * pb_xxz[j] + 0.75 * fx[j] * fx[j] * pb_xxxz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxxz[j] +

                                 pa_xxxx[j] * pb_xxxz[j]) * s_0_0[j];

                t_xxxx_xxyy[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.875 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * pa_xxxx[j] * fx[j] * fx[j] +

                                 2.0 * pa_xxx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 6.0 * pa_x[j] * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_xxxx[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxxx[j] * fx[j] * pb_yy[j] +

                                 4.0 * pa_xxx[j] * fx[j] * pb_xyy[j] + 0.75 * fx[j] * fx[j] * pb_xxyy[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxyy[j] +

                                 pa_xxxx[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (4,6)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xx, pb_xxyz, pb_xxzz, pb_xyz, pb_xzz, \
                                     pb_yz, pb_zz, s_0_0, t_xxxx_xxyz, t_xxxx_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xxyz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_yz[j] +

                                 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 6.0 * pa_x[j] * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xxxx[j] * fx[j] * pb_yz[j] +

                                 4.0 * pa_xxx[j] * fx[j] * pb_xyz[j] + 0.75 * fx[j] * fx[j] * pb_xxyz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxyz[j] +

                                 pa_xxxx[j] * pb_xxyz[j]) * s_0_0[j];

                t_xxxx_xxzz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.875 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * pa_xxxx[j] * fx[j] * fx[j] +

                                 2.0 * pa_xxx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 6.0 * pa_x[j] * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_xxxx[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxxx[j] * fx[j] * pb_zz[j] +

                                 4.0 * pa_xxx[j] * fx[j] * pb_xzz[j] + 0.75 * fx[j] * fx[j] * pb_xxzz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxzz[j] +

                                 pa_xxxx[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (6,8)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_xy, pb_xyyy, pb_xyyz, pb_xz, pb_y, pb_yyy, \
                                     pb_yyz, pb_z, s_0_0, t_xxxx_xyyy, t_xxxx_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xyyy[j] = (4.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 3.0 * pa_xxx[j] * fx[j] * fx[j] * pb_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xxxx[j] * pb_xy[j] * fx[j] +

                                 2.0 * pa_xxx[j] * fx[j] * pb_yyy[j] + 0.75 * fx[j] * fx[j] * pb_xyyy[j] + 3.0 * pa_xx[j] * fx[j] * pb_xyyy[j] +

                                 pa_xxxx[j] * pb_xyyy[j]) * s_0_0[j];

                t_xxxx_xyyz[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 pa_xxx[j] * fx[j] * fx[j] * pb_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xxxx[j] * pb_xz[j] * fx[j] +

                                 2.0 * pa_xxx[j] * fx[j] * pb_yyz[j] + 0.75 * fx[j] * fx[j] * pb_xyyz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xyyz[j] +

                                 pa_xxxx[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (8,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_xy, pb_xyzz, pb_xz, pb_xzzz, pb_y, pb_yzz, \
                                     pb_z, pb_zzz, s_0_0, t_xxxx_xyzz, t_xxxx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xyzz[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 pa_xxx[j] * fx[j] * fx[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xxxx[j] * pb_xy[j] * fx[j] +

                                 2.0 * pa_xxx[j] * fx[j] * pb_yzz[j] + 0.75 * fx[j] * fx[j] * pb_xyzz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xyzz[j] +

                                 pa_xxxx[j] * pb_xyzz[j]) * s_0_0[j];

                t_xxxx_xzzz[j] = (4.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 3.0 * pa_xxx[j] * fx[j] * fx[j] * pb_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xxxx[j] * pb_xz[j] * fx[j] +

                                 2.0 * pa_xxx[j] * fx[j] * pb_zzz[j] + 0.75 * fx[j] * fx[j] * pb_xzzz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xzzz[j] +

                                 pa_xxxx[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (5) = (10,12)

            #pragma omp simd aligned(fx, pa_xx, pa_xxxx, pb_yy, pb_yyyy, pb_yyyz, pb_yz, s_0_0, t_xxxx_yyyy, \
                                     t_xxxx_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_yyyy[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_xxxx[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 9.0 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 3.0 * pa_xxxx[j] * pb_yy[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yyyy[j] + 3.0 * pa_xx[j] * fx[j] * pb_yyyy[j] +

                                 pa_xxxx[j] * pb_yyyy[j]) * s_0_0[j];

                t_xxxx_yyyz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_yz[j] +

                                 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xxxx[j] * pb_yz[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yyyz[j] +

                                 3.0 * pa_xx[j] * fx[j] * pb_yyyz[j] + pa_xxxx[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (6) = (12,14)

            #pragma omp simd aligned(fx, pa_xx, pa_xxxx, pb_yy, pb_yyzz, pb_yz, pb_yzzz, pb_zz, s_0_0, t_xxxx_yyzz, \
                                     t_xxxx_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_yyzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.25 * pa_xxxx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xxxx[j] * pb_yy[j] * fx[j] +

                                 0.5 * pa_xxxx[j] * fx[j] * pb_zz[j] + 0.75 * fx[j] * fx[j] * pb_yyzz[j] + 3.0 * pa_xx[j] * fx[j] * pb_yyzz[j] +

                                 pa_xxxx[j] * pb_yyzz[j]) * s_0_0[j];

                t_xxxx_yzzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_yz[j] +

                                 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xxxx[j] * pb_yz[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yzzz[j] +

                                 3.0 * pa_xx[j] * fx[j] * pb_yzzz[j] + pa_xxxx[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (7) = (14,16)

            #pragma omp simd aligned(fx, pa_xx, pa_xxxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_zz, pb_zzzz, s_0_0, t_xxxx_zzzz, t_xxxy_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_zzzz[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_xxxx[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 9.0 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 3.0 * pa_xxxx[j] * pb_zz[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_zzzz[j] + 3.0 * pa_xx[j] * fx[j] * pb_zzzz[j] +

                                 pa_xxxx[j] * pb_zzzz[j]) * s_0_0[j];

                t_xxxy_xxxx[j] = (5.625 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 7.5 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.75 * pa_xxxy[j] * fx[j] * fx[j] + 9.0 * pa_xxy[j] * fx[j] * fx[j] * pb_x[j] +

                                 13.5 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * fx[j] * fx[j] * pa_y[j] * pb_xxx[j] + 3.0 * pa_xxxy[j] * pb_xx[j] * fx[j] +

                                 6.0 * pa_xxy[j] * fx[j] * pb_xxx[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxxx[j] + pa_xxxy[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (8) = (16,18)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxz, pb_xy, pb_xz, pb_y, pb_z, s_0_0, t_xxxy_xxxy, \
                                     t_xxxy_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xxxy[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 1.125 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 3.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.875 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xxy[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 6.75 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xxx[j] +

                                 2.25 * fx[j] * fx[j] * pa_y[j] * pb_xxy[j] + 1.5 * pa_xxxy[j] * pb_xy[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xxx[j] +

                                 4.5 * pa_xxy[j] * fx[j] * pb_xxy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxxy[j] + pa_xxxy[j] * pb_xxxy[j]) * s_0_0[j];

                t_xxxy_xxxz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 2.25 * pa_xxy[j] * fx[j] * fx[j] * pb_z[j] + 6.75 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] +

                                 2.25 * fx[j] * fx[j] * pa_y[j] * pb_xxz[j] + 1.5 * pa_xxxy[j] * pb_xz[j] * fx[j] + 4.5 * pa_xxy[j] * fx[j] * pb_xxz[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_xxxz[j] + pa_xxxy[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (9) = (18,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, \
                                     t_xxxy_xxyy, t_xxxy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xxyy[j] = (1.125 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * pa_xxxy[j] * fx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * pa_xxy[j] * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xxy[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_xyy[j] + 0.5 * pa_xxxy[j] * pb_xx[j] * fx[j] +

                                 0.5 * pa_xxxy[j] * fx[j] * pb_yy[j] + pa_xxx[j] * fx[j] * pb_xxy[j] + 3.0 * pa_xxy[j] * fx[j] * pb_xyy[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_xxyy[j] + pa_xxxy[j] * pb_xxyy[j]) * s_0_0[j];

                t_xxxy_xxyz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xxz[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_xyz[j] + 0.5 * pa_xxxy[j] * fx[j] * pb_yz[j] +

                                 0.5 * pa_xxx[j] * fx[j] * pb_xxz[j] + 3.0 * pa_xxy[j] * fx[j] * pb_xyz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxyz[j] +

                                 pa_xxxy[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (10) = (20,22)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xzz, pb_y, pb_yy, pb_yyy, pb_zz, s_0_0, t_xxxy_xxzz, t_xxxy_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xxzz[j] = (1.125 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.25 * pa_xxxy[j] * fx[j] * fx[j] + 1.5 * pa_xxy[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] +

                                 1.5 * fx[j] * fx[j] * pa_y[j] * pb_xzz[j] + 0.5 * pa_xxxy[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxxy[j] * fx[j] * pb_zz[j] +

                                 3.0 * pa_xxy[j] * fx[j] * pb_xzz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxzz[j] + pa_xxxy[j] * pb_xxzz[j]) * s_0_0[j];

                t_xxxy_xyyy[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 1.125 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xxy[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_xyy[j] +

                                 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yyy[j] + 1.5 * pa_xxxy[j] * pb_xy[j] * fx[j] + 1.5 * pa_xxx[j] * fx[j] * pb_xyy[j] +

                                 1.5 * pa_xxy[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyyy[j] + pa_xxxy[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (11) = (22,24)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xy, pb_xyyz, \
                                     pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_y, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, \
                                     t_xxxy_xyyz, t_xxxy_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xyz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yyz[j] + 0.5 * pa_xxxy[j] * pb_xz[j] * fx[j] +

                                 pa_xxx[j] * fx[j] * pb_xyz[j] + 1.5 * pa_xxy[j] * fx[j] * pb_yyz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyyz[j] +

                                 pa_xxxy[j] * pb_xyyz[j]) * s_0_0[j];

                t_xxxy_xyzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * pa_xxy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xzz[j] +

                                 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yzz[j] + 0.5 * pa_xxxy[j] * pb_xy[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xzz[j] +

                                 1.5 * pa_xxy[j] * fx[j] * pb_yzz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyzz[j] + pa_xxxy[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (12) = (24,26)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_xz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyyy, pb_z, pb_zzz, s_0_0, t_xxxy_xzzz, t_xxxy_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xzzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 2.25 * pa_xxy[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zzz[j] + 1.5 * pa_xxxy[j] * pb_xz[j] * fx[j] + 1.5 * pa_xxy[j] * fx[j] * pb_zzz[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_xzzz[j] + pa_xxxy[j] * pb_xzzz[j]) * s_0_0[j];

                t_xxxy_yyyy[j] = (1.125 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 4.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xxxy[j] * fx[j] * fx[j] + 3.0 * pa_xxx[j] * fx[j] * fx[j] * pb_y[j] +

                                 4.5 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_xxxy[j] * pb_yy[j] * fx[j] +

                                 2.0 * pa_xxx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xy[j] * fx[j] * pb_yyyy[j] + pa_xxxy[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (13) = (26,28)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxy, pa_xy, pb_y, pb_yy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, s_0_0, t_xxxy_yyyz, t_xxxy_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_yyyz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * pa_xxx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_xxxy[j] * pb_yz[j] * fx[j] + 1.5 * pa_xxx[j] * fx[j] * pb_yyz[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_yyyz[j] + pa_xxxy[j] * pb_yyyz[j]) * s_0_0[j];

                t_xxxy_yyzz[j] = (0.375 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xxxy[j] * fx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xxxy[j] * pb_yy[j] * fx[j] + 0.5 * pa_xxxy[j] * fx[j] * pb_zz[j] +

                                 pa_xxx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_xy[j] * fx[j] * pb_yyzz[j] + pa_xxxy[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (14) = (28,30)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxy, pa_xy, pb_yz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     s_0_0, t_xxxy_yzzz, t_xxxy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_yzzz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * pa_xxx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xxxy[j] * pb_yz[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_zzz[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_yzzz[j] + pa_xxxy[j] * pb_yzzz[j]) * s_0_0[j];

                t_xxxy_zzzz[j] = (1.125 * pa_xy[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxxy[j] * fx[j] * fx[j] +

                                 4.5 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xxxy[j] * pb_zz[j] * fx[j] + 1.5 * pa_xy[j] * fx[j] * pb_zzzz[j] +

                                 pa_xxxy[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (15) = (30,32)

            #pragma omp simd aligned(fx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxy, \
                                     pb_xy, pb_y, s_0_0, t_xxxz_xxxx, t_xxxz_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xxxx[j] = (5.625 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 7.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.75 * pa_xxxz[j] * fx[j] * fx[j] + 9.0 * pa_xxz[j] * fx[j] * fx[j] * pb_x[j] +

                                 13.5 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * fx[j] * fx[j] * pa_z[j] * pb_xxx[j] + 3.0 * pa_xxxz[j] * pb_xx[j] * fx[j] +

                                 6.0 * pa_xxz[j] * fx[j] * pb_xxx[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxxx[j] + pa_xxxz[j] * pb_xxxx[j]) * s_0_0[j];

                t_xxxz_xxxy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 2.25 * pa_xxz[j] * fx[j] * fx[j] * pb_y[j] + 6.75 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xxy[j] + 1.5 * pa_xxxz[j] * pb_xy[j] * fx[j] + 4.5 * pa_xxz[j] * fx[j] * pb_xxy[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_xxxy[j] + pa_xxxz[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (16) = (32,34)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxz, pb_xxyy, pb_xxz, pb_xyy, pb_xz, pb_yy, pb_z, s_0_0, t_xxxz_xxxz, t_xxxz_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xxxz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 1.125 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 3.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xxz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 6.75 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xxx[j] +

                                 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xxz[j] + 1.5 * pa_xxxz[j] * pb_xz[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xxx[j] +

                                 4.5 * pa_xxz[j] * fx[j] * pb_xxz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxxz[j] + pa_xxxz[j] * pb_xxxz[j]) * s_0_0[j];

                t_xxxz_xxyy[j] = (1.125 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.25 * pa_xxxz[j] * fx[j] * fx[j] + 1.5 * pa_xxz[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] +

                                 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xyy[j] + 0.5 * pa_xxxz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxxz[j] * fx[j] * pb_yy[j] +

                                 3.0 * pa_xxz[j] * fx[j] * pb_xyy[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxyy[j] + pa_xxxz[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (17) = (34,36)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxy, \
                                     pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, \
                                     t_xxxz_xxyz, t_xxxz_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xxyz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xxy[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xyz[j] + 0.5 * pa_xxxz[j] * fx[j] * pb_yz[j] +

                                 0.5 * pa_xxx[j] * fx[j] * pb_xxy[j] + 3.0 * pa_xxz[j] * fx[j] * pb_xyz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxyz[j] +

                                 pa_xxxz[j] * pb_xxyz[j]) * s_0_0[j];

                t_xxxz_xxzz[j] = (1.125 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * pa_xxxz[j] * fx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * pa_xxz[j] * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xxz[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xzz[j] + 0.5 * pa_xxxz[j] * pb_xx[j] * fx[j] +

                                 0.5 * pa_xxxz[j] * fx[j] * pb_zz[j] + pa_xxx[j] * fx[j] * pb_xxz[j] + 3.0 * pa_xxz[j] * fx[j] * pb_xzz[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_xxzz[j] + pa_xxxz[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (18) = (36,38)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_z, s_0_0, t_xxxz_xyyy, \
                                     t_xxxz_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xyyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 2.25 * pa_xxz[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yyy[j] + 1.5 * pa_xxxz[j] * pb_xy[j] * fx[j] + 1.5 * pa_xxz[j] * fx[j] * pb_yyy[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_xyyy[j] + pa_xxxz[j] * pb_xyyy[j]) * s_0_0[j];

                t_xxxz_xyyz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * pa_xxz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xyy[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yyz[j] + 0.5 * pa_xxxz[j] * pb_xz[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xyy[j] +

                                 1.5 * pa_xxz[j] * fx[j] * pb_yyz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyyz[j] + pa_xxxz[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (19) = (38,40)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xy, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xxxz_xyzz, t_xxxz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xyz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yzz[j] + 0.5 * pa_xxxz[j] * pb_xy[j] * fx[j] +

                                 pa_xxx[j] * fx[j] * pb_xyz[j] + 1.5 * pa_xxz[j] * fx[j] * pb_yzz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyzz[j] +

                                 pa_xxxz[j] * pb_xyzz[j]) * s_0_0[j];

                t_xxxz_xzzz[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 1.125 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xxz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_xzz[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zzz[j] + 1.5 * pa_xxxz[j] * pb_xz[j] * fx[j] + 1.5 * pa_xxx[j] * fx[j] * pb_xzz[j] +

                                 1.5 * pa_xxz[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xzzz[j] + pa_xxxz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (20) = (40,42)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxz, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yz, \
                                     s_0_0, t_xxxz_yyyy, t_xxxz_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_yyyy[j] = (1.125 * pa_xz[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxxz[j] * fx[j] * fx[j] +

                                 4.5 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xxxz[j] * pb_yy[j] * fx[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyyy[j] +

                                 pa_xxxz[j] * pb_yyyy[j]) * s_0_0[j];

                t_xxxz_yyyz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * pa_xxx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xxxz[j] * pb_yz[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_yyy[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_yyyz[j] + pa_xxxz[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (21) = (42,44)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxz, pa_xz, pb_y, pb_yy, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, s_0_0, t_xxxz_yyzz, t_xxxz_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_yyzz[j] = (0.375 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_xxxz[j] * fx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xxxz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xxxz[j] * fx[j] * pb_zz[j] +

                                 pa_xxx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyzz[j] + pa_xxxz[j] * pb_yyzz[j]) * s_0_0[j];

                t_xxxz_yzzz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * pa_xxx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_xxxz[j] * pb_yz[j] * fx[j] + 1.5 * pa_xxx[j] * fx[j] * pb_yzz[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_yzzz[j] + pa_xxxz[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (22) = (44,46)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxyy, pa_xyy, pa_xz, pa_yy, pb_x, pb_xx, \
                                     pb_xxx, pb_xxxx, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xxxz_zzzz, t_xxyy_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_zzzz[j] = (1.125 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 4.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xxxz[j] * fx[j] * fx[j] + 3.0 * pa_xxx[j] * fx[j] * fx[j] * pb_z[j] +

                                 4.5 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_xxxz[j] * pb_zz[j] * fx[j] +

                                 2.0 * pa_xxx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xz[j] * fx[j] * pb_zzzz[j] + pa_xxxz[j] * pb_zzzz[j]) * s_0_0[j];

                t_xxyy_xxxx[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 1.875 * fx[j] * fx[j] * fx[j] * pa_yy[j] +

                                 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_xxyy[j] * fx[j] * fx[j] + 6.0 * pa_xyy[j] * fx[j] * fx[j] * pb_x[j] +

                                 4.5 * fx[j] * fx[j] * pa_yy[j] * pb_xx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 2.0 * pa_x[j] * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xxyy[j] * pb_xx[j] * fx[j] + 4.0 * pa_xyy[j] * fx[j] * pb_xxx[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxxx[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxxx[j] + 0.5 * fx[j] * pa_yy[j] * pb_xxxx[j] +

                                 pa_xxyy[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (23) = (46,48)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxz, pb_xy, pb_xz, pb_y, pb_z, s_0_0, t_xxyy_xxxy, \
                                     t_xxyy_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xxxy[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_xxy[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * pa_xyy[j] * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] +

                                 2.25 * fx[j] * fx[j] * pa_yy[j] * pb_xy[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xxx[j] + 1.5 * pa_xxyy[j] * pb_xy[j] * fx[j] +

                                 pa_xxy[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xyy[j] * fx[j] * pb_xxy[j] + 0.25 * fx[j] * fx[j] * pb_xxxy[j] +

                                 0.5 * pa_xx[j] * fx[j] * pb_xxxy[j] + 0.5 * fx[j] * pa_yy[j] * pb_xxxy[j] + pa_xxyy[j] * pb_xxxy[j]) * s_0_0[j];

                t_xxyy_xxxz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_xyy[j] * fx[j] * fx[j] * pb_z[j] +

                                 2.25 * fx[j] * fx[j] * pa_yy[j] * pb_xz[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xxz[j] + 1.5 * pa_xxyy[j] * pb_xz[j] * fx[j] + 3.0 * pa_xyy[j] * fx[j] * pb_xxz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxxz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxxz[j] + 0.5 * fx[j] * pa_yy[j] * pb_xxxz[j] +

                                 pa_xxyy[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (24) = (48,50)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, \
                                     t_xxyy_xxyy, t_xxyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xxyy[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_yy[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.25 * pa_xxyy[j] * fx[j] * fx[j] + pa_xxy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 pa_xyy[j] * fx[j] * fx[j] * pb_x[j] + 4.0 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_yy[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] + pa_x[j] * fx[j] * fx[j] * pb_xyy[j] +

                                 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_xx[j] + fx[j] * fx[j] * pa_y[j] * pb_xxy[j] + 0.5 * pa_xxyy[j] * pb_xx[j] * fx[j] +

                                 0.5 * pa_xxyy[j] * fx[j] * pb_yy[j] + 2.0 * pa_xxy[j] * fx[j] * pb_xxy[j] + 2.0 * pa_xyy[j] * fx[j] * pb_xyy[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxyy[j] + 0.5 * fx[j] * pa_yy[j] * pb_xxyy[j] +

                                 pa_xxyy[j] * pb_xxyy[j]) * s_0_0[j];

                t_xxyy_xxyz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xxy[j] * fx[j] * fx[j] * pb_z[j] +

                                 2.0 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_yz[j] +

                                 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + pa_x[j] * fx[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xxz[j] +

                                 0.5 * pa_xxyy[j] * fx[j] * pb_yz[j] + pa_xxy[j] * fx[j] * pb_xxz[j] + 2.0 * pa_xyy[j] * fx[j] * pb_xyz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxyz[j] + 0.5 * fx[j] * pa_yy[j] * pb_xxyz[j] +

                                 pa_xxyy[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (25) = (50,52)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xzz, pb_y, pb_yy, pb_yyy, pb_zz, s_0_0, t_xxyy_xxzz, \
                                     t_xxyy_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xxzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_yy[j] +

                                 0.125 * pa_xx[j] * fx[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * pa_xxyy[j] * fx[j] * fx[j] + pa_xyy[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_zz[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] + pa_x[j] * fx[j] * fx[j] * pb_xzz[j] +

                                 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_xx[j] + 0.5 * pa_xxyy[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxyy[j] * fx[j] * pb_zz[j] +

                                 2.0 * pa_xyy[j] * fx[j] * pb_xzz[j] + 0.25 * fx[j] * fx[j] * pb_xxzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxzz[j] +

                                 0.5 * fx[j] * pa_yy[j] * pb_xxzz[j] + pa_xxyy[j] * pb_xxzz[j]) * s_0_0[j];

                t_xxyy_xyyy[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_xxy[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_xyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 3.0 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yyy[j] +

                                 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_xy[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_xyy[j] + 1.5 * pa_xxyy[j] * pb_xy[j] * fx[j] +

                                 3.0 * pa_xxy[j] * fx[j] * pb_xyy[j] + pa_xyy[j] * fx[j] * pb_yyy[j] + 0.25 * fx[j] * fx[j] * pb_xyyy[j] +

                                 0.5 * pa_xx[j] * fx[j] * pb_xyyy[j] + 0.5 * fx[j] * pa_yy[j] * pb_xyyy[j] + pa_xxyy[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (26) = (52,54)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xy, pb_xyyz, \
                                     pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_y, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, \
                                     t_xxyy_xyyz, t_xxyy_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xyyz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.5 * pa_xyy[j] * fx[j] * fx[j] * pb_z[j] + 2.0 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yyz[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_xz[j] + fx[j] * fx[j] * pa_y[j] * pb_xyz[j] +

                                 0.5 * pa_xxyy[j] * pb_xz[j] * fx[j] + 2.0 * pa_xxy[j] * fx[j] * pb_xyz[j] + pa_xyy[j] * fx[j] * pb_yyz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xyyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyyz[j] + 0.5 * fx[j] * pa_yy[j] * pb_xyyz[j] +

                                 pa_xxyy[j] * pb_xyyz[j]) * s_0_0[j];

                t_xxyy_xyzz[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                                 0.5 * pa_xxy[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xyy[j] * fx[j] * fx[j] * pb_y[j] + pa_xy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.125 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yzz[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_xy[j] +

                                 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xzz[j] + 0.5 * pa_xxyy[j] * pb_xy[j] * fx[j] + pa_xxy[j] * fx[j] * pb_xzz[j] +

                                 pa_xyy[j] * fx[j] * pb_yzz[j] + 0.25 * fx[j] * fx[j] * pb_xyzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyzz[j] +

                                 0.5 * fx[j] * pa_yy[j] * pb_xyzz[j] + pa_xxyy[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (27) = (54,56)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xyy, pa_y, pa_yy, pb_xz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyyy, pb_z, pb_zzz, s_0_0, t_xxyy_xzzz, t_xxyy_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xzzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * pa_xyy[j] * fx[j] * fx[j] * pb_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_zzz[j] +

                                 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_xz[j] + 1.5 * pa_xxyy[j] * pb_xz[j] * fx[j] + pa_xyy[j] * fx[j] * pb_zzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xzzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xzzz[j] + 0.5 * fx[j] * pa_yy[j] * pb_xzzz[j] +

                                 pa_xxyy[j] * pb_xzzz[j]) * s_0_0[j];

                t_xxyy_yyyy[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 1.875 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_yy[j] + 3.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xxyy[j] * fx[j] * fx[j] + 6.0 * pa_xxy[j] * fx[j] * fx[j] * pb_y[j] +

                                 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * fx[j] * fx[j] * pa_yy[j] * pb_yy[j] +

                                 2.0 * fx[j] * fx[j] * pa_y[j] * pb_yyy[j] + 3.0 * pa_xxyy[j] * pb_yy[j] * fx[j] + 4.0 * pa_xxy[j] * fx[j] * pb_yyy[j] +

                                 0.25 * fx[j] * fx[j] * pb_yyyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyyy[j] + 0.5 * fx[j] * pa_yy[j] * pb_yyyy[j] +

                                 pa_xxyy[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (28) = (56,58)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyy, pa_y, pa_yy, pb_y, pb_yy, pb_yyyz, pb_yyz, pb_yyzz, \
                                     pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, t_xxyy_yyyz, t_xxyy_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_yyyz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xxy[j] * fx[j] * fx[j] * pb_z[j] +

                                 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_yz[j] +

                                 1.5 * fx[j] * fx[j] * pa_y[j] * pb_yyz[j] + 1.5 * pa_xxyy[j] * pb_yz[j] * fx[j] + 3.0 * pa_xxy[j] * fx[j] * pb_yyz[j] +

                                 0.25 * fx[j] * fx[j] * pb_yyyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyyz[j] + 0.5 * fx[j] * pa_yy[j] * pb_yyyz[j] +

                                 pa_xxyy[j] * pb_yyyz[j]) * s_0_0[j];

                t_xxyy_yyzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.125 * fx[j] * fx[j] * fx[j] * pa_yy[j] + 0.5 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * pa_xxyy[j] * fx[j] * fx[j] + pa_xxy[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_yy[j] +

                                 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_zz[j] + fx[j] * fx[j] * pa_y[j] * pb_yzz[j] + 0.5 * pa_xxyy[j] * pb_yy[j] * fx[j] +

                                 0.5 * pa_xxyy[j] * fx[j] * pb_zz[j] + 2.0 * pa_xxy[j] * fx[j] * pb_yzz[j] + 0.25 * fx[j] * fx[j] * pb_yyzz[j] +

                                 0.5 * pa_xx[j] * fx[j] * pb_yyzz[j] + 0.5 * fx[j] * pa_yy[j] * pb_yyzz[j] + pa_xxyy[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (29) = (58,60)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyy, pa_y, pa_yy, pb_yz, pb_yzzz, pb_z, pb_zz, pb_zzz, \
                                     pb_zzzz, s_0_0, t_xxyy_yzzz, t_xxyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_yzzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 1.5 * pa_xxy[j] * fx[j] * fx[j] * pb_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_yz[j] +

                                 0.5 * fx[j] * fx[j] * pa_y[j] * pb_zzz[j] + 1.5 * pa_xxyy[j] * pb_yz[j] * fx[j] + pa_xxy[j] * fx[j] * pb_zzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_yzzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yzzz[j] + 0.5 * fx[j] * pa_yy[j] * pb_yzzz[j] +

                                 pa_xxyy[j] * pb_yzzz[j]) * s_0_0[j];

                t_xxyy_zzzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_yy[j] + 0.75 * pa_xxyy[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * fx[j] * fx[j] * pa_yy[j] * pb_zz[j] + 3.0 * pa_xxyy[j] * pb_zz[j] * fx[j] +

                                 0.25 * fx[j] * fx[j] * pb_zzzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_zzzz[j] + 0.5 * fx[j] * pa_yy[j] * pb_zzzz[j] +

                                 pa_xxyy[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (30) = (60,62)

            #pragma omp simd aligned(fx, pa_xxyz, pa_xxz, pa_xyz, pa_xz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxy, pb_xy, pb_y, s_0_0, t_xxyz_xxxx, t_xxyz_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xxxx[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_yz[j] + 0.75 * pa_xxyz[j] * fx[j] * fx[j] +

                                 6.0 * pa_xyz[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * fx[j] * fx[j] * pa_yz[j] * pb_xx[j] + 3.0 * pa_xxyz[j] * pb_xx[j] * fx[j] +

                                 4.0 * pa_xyz[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxxx[j] + pa_xxyz[j] * pb_xxxx[j]) * s_0_0[j];

                t_xxyz_xxxy[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] +

                                 2.25 * fx[j] * fx[j] * pa_yz[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_xxx[j] + 1.5 * pa_xxyz[j] * pb_xy[j] * fx[j] +

                                 0.5 * pa_xxz[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xyz[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxxy[j] +

                                 pa_xxyz[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (31) = (62,64)

            #pragma omp simd aligned(fx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, \
                                     pb_xxx, pb_xxxz, pb_xxy, pb_xxyy, pb_xxz, pb_xy, pb_xyy, pb_xz, pb_y, pb_yy, pb_z, s_0_0, \
                                     t_xxyz_xxxz, t_xxyz_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xxxz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] +

                                 2.25 * fx[j] * fx[j] * pa_yz[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_xxx[j] + 1.5 * pa_xxyz[j] * pb_xz[j] * fx[j] +

                                 0.5 * pa_xxy[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xyz[j] * fx[j] * pb_xxz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxxz[j] +

                                 pa_xxyz[j] * pb_xxxz[j]) * s_0_0[j];

                t_xxyz_xxyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.25 * pa_xxyz[j] * fx[j] * fx[j] + 0.5 * pa_xxz[j] * fx[j] * fx[j] * pb_y[j] +

                                 pa_xyz[j] * fx[j] * fx[j] * pb_x[j] + 2.0 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_yz[j] * pb_xx[j] +

                                 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xxy[j] + 0.5 * pa_xxyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_yy[j] +

                                 pa_xxz[j] * fx[j] * pb_xxy[j] + 2.0 * pa_xyz[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxyy[j] +

                                 pa_xxyz[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (32) = (64,66)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, \
                                     pb_x, pb_xx, pb_xxy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyz, pb_xz, pb_xzz, pb_y, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_xxyz_xxyz, t_xxyz_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xxyz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.125 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.25 * pa_xxy[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xxz[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] + pa_xy[j] * fx[j] * fx[j] * pb_xy[j] + pa_xz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_xxy[j] +

                                 0.25 * fx[j] * fx[j] * pa_z[j] * pb_xxz[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xxy[j] * fx[j] * pb_xxy[j] +

                                 0.5 * pa_xxz[j] * fx[j] * pb_xxz[j] + 2.0 * pa_xyz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxyz[j] +

                                 pa_xxyz[j] * pb_xxyz[j]) * s_0_0[j];

                t_xxyz_xxzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 0.25 * pa_xxyz[j] * fx[j] * fx[j] + 0.5 * pa_xxy[j] * fx[j] * fx[j] * pb_z[j] +

                                 pa_xyz[j] * fx[j] * fx[j] * pb_x[j] + 2.0 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_yz[j] * pb_xx[j] +

                                 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xxz[j] + 0.5 * pa_xxyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_zz[j] +

                                 pa_xxy[j] * fx[j] * pb_xxz[j] + 2.0 * pa_xyz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxzz[j] +

                                 pa_xxyz[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (33) = (66,68)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, \
                                     pb_x, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_xxyz_xyyy, t_xxyz_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xyyy[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xyy[j] + 1.5 * pa_xxyz[j] * pb_xy[j] * fx[j] +

                                 1.5 * pa_xxz[j] * fx[j] * pb_xyy[j] + pa_xyz[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_yz[j] * pb_xyyy[j] +

                                 pa_xxyz[j] * pb_xyyy[j]) * s_0_0[j];

                t_xxyz_xyyz[j] = (0.25 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                                 0.25 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * pa_xxy[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xyz[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] + pa_xz[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pa_yz[j] * pb_xz[j] +

                                 0.25 * fx[j] * fx[j] * pa_y[j] * pb_xyy[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xyz[j] + 0.5 * pa_xxyz[j] * pb_xz[j] * fx[j] +

                                 0.5 * pa_xxy[j] * fx[j] * pb_xyy[j] + pa_xxz[j] * fx[j] * pb_xyz[j] + pa_xyz[j] * fx[j] * pb_yyz[j] +

                                 0.5 * fx[j] * pa_yz[j] * pb_xyyz[j] + pa_xxyz[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (34) = (68,70)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, \
                                     pb_x, pb_xy, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, s_0_0, t_xxyz_xyzz, t_xxyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xyzz[j] = (0.25 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 0.25 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * pa_xxz[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xyz[j] * fx[j] * fx[j] * pb_y[j] + pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.5 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_yz[j] * pb_xy[j] +

                                 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xyz[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_xzz[j] + 0.5 * pa_xxyz[j] * pb_xy[j] * fx[j] +

                                 pa_xxy[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xxz[j] * fx[j] * pb_xzz[j] + pa_xyz[j] * fx[j] * pb_yzz[j] +

                                 0.5 * fx[j] * pa_yz[j] * pb_xyzz[j] + pa_xxyz[j] * pb_xyzz[j]) * s_0_0[j];

                t_xxyz_xzzz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_xzz[j] + 1.5 * pa_xxyz[j] * pb_xz[j] * fx[j] +

                                 1.5 * pa_xxy[j] * fx[j] * pb_xzz[j] + pa_xyz[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xzzz[j] +

                                 pa_xxyz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (35) = (70,72)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_y, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yz, pb_z, s_0_0, t_xxyz_yyyy, t_xxyz_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_yyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.75 * pa_xxyz[j] * fx[j] * fx[j] + 3.0 * pa_xxz[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * fx[j] * fx[j] * pa_yz[j] * pb_yy[j] + fx[j] * fx[j] * pa_z[j] * pb_yyy[j] + 3.0 * pa_xxyz[j] * pb_yy[j] * fx[j] +

                                 2.0 * pa_xxz[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyyy[j] + pa_xxyz[j] * pb_yyyy[j]) * s_0_0[j];

                t_xxyz_yyyz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * pa_xxz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_yyy[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yyz[j] + 1.5 * pa_xxyz[j] * pb_yz[j] * fx[j] + 0.5 * pa_xxy[j] * fx[j] * pb_yyy[j] +

                                 1.5 * pa_xxz[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyyz[j] + pa_xxyz[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (36) = (72,74)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_y, pa_yz, pa_z, pb_y, pb_yy, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxyz_yyzz, t_xxyz_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_yyzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 0.25 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 0.5 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * pa_xxyz[j] * fx[j] * fx[j] + 0.5 * pa_xxy[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.5 * pa_xxz[j] * fx[j] * fx[j] * pb_y[j] + pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pa_yz[j] * pb_yy[j] +

                                 0.25 * fx[j] * fx[j] * pa_yz[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_yyz[j] +

                                 0.5 * fx[j] * fx[j] * pa_z[j] * pb_yzz[j] + 0.5 * pa_xxyz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_zz[j] +

                                 pa_xxy[j] * fx[j] * pb_yyz[j] + pa_xxz[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyzz[j] +

                                 pa_xxyz[j] * pb_yyzz[j]) * s_0_0[j];

                t_xxyz_yzzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * pa_xxz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yzz[j] +

                                 0.25 * fx[j] * fx[j] * pa_z[j] * pb_zzz[j] + 1.5 * pa_xxyz[j] * pb_yz[j] * fx[j] + 1.5 * pa_xxy[j] * fx[j] * pb_yzz[j] +

                                 0.5 * pa_xxz[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yzzz[j] + pa_xxyz[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (37) = (74,76)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxzz, pa_xzz, pa_y, pa_yz, pa_zz, pb_x, pb_xx, \
                                     pb_xxx, pb_xxxx, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xxyz_zzzz, t_xxzz_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_zzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 0.75 * pa_xxyz[j] * fx[j] * fx[j] + 3.0 * pa_xxy[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * fx[j] * fx[j] * pa_yz[j] * pb_zz[j] + fx[j] * fx[j] * pa_y[j] * pb_zzz[j] + 3.0 * pa_xxyz[j] * pb_zz[j] * fx[j] +

                                 2.0 * pa_xxy[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_yz[j] * pb_zzzz[j] + pa_xxyz[j] * pb_zzzz[j]) * s_0_0[j];

                t_xxzz_xxxx[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 1.875 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_xxzz[j] * fx[j] * fx[j] + 6.0 * pa_xzz[j] * fx[j] * fx[j] * pb_x[j] +

                                 4.5 * fx[j] * fx[j] * pa_zz[j] * pb_xx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 2.0 * pa_x[j] * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xxzz[j] * pb_xx[j] * fx[j] + 4.0 * pa_xzz[j] * fx[j] * pb_xxx[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxxx[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxxx[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxxx[j] +

                                 pa_xxzz[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (38) = (76,78)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxz, pb_xy, pb_xz, pb_y, pb_z, s_0_0, t_xxzz_xxxy, \
                                     t_xxzz_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xxxy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_xzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_xy[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xxy[j] + 1.5 * pa_xxzz[j] * pb_xy[j] * fx[j] + 3.0 * pa_xzz[j] * fx[j] * pb_xxy[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxxy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxxy[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxxy[j] +

                                 pa_xxzz[j] * pb_xxxy[j]) * s_0_0[j];

                t_xxzz_xxxz[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_xxz[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * pa_xzz[j] * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] +

                                 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_xz[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xxz[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xxx[j] + 1.5 * pa_xxzz[j] * pb_xz[j] * fx[j] +

                                 pa_xxz[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xzz[j] * fx[j] * pb_xxz[j] + 0.25 * fx[j] * fx[j] * pb_xxxz[j] +

                                 0.5 * pa_xx[j] * fx[j] * pb_xxxz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxxz[j] + pa_xxzz[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (39) = (78,80)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xy, pb_xyy, pb_xyz, pb_y, pb_yy, pb_yz, s_0_0, t_xxzz_xxyy, \
                                     t_xxzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xxyy[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 0.125 * pa_xx[j] * fx[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * pa_xxzz[j] * fx[j] * fx[j] + pa_xzz[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_yy[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] + pa_x[j] * fx[j] * fx[j] * pb_xyy[j] +

                                 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_xx[j] + 0.5 * pa_xxzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxzz[j] * fx[j] * pb_yy[j] +

                                 2.0 * pa_xzz[j] * fx[j] * pb_xyy[j] + 0.25 * fx[j] * fx[j] * pb_xxyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxyy[j] +

                                 0.5 * fx[j] * pa_zz[j] * pb_xxyy[j] + pa_xxzz[j] * pb_xxyy[j]) * s_0_0[j];

                t_xxzz_xxyz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xxz[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.0 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_yz[j] +

                                 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + pa_x[j] * fx[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xxy[j] +

                                 0.5 * pa_xxzz[j] * fx[j] * pb_yz[j] + pa_xxz[j] * fx[j] * pb_xxy[j] + 2.0 * pa_xzz[j] * fx[j] * pb_xyz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxyz[j] +

                                 pa_xxzz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (40) = (80,82)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyyy, pb_xz, pb_xzz, pb_y, pb_yyy, pb_z, pb_zz, s_0_0, t_xxzz_xxzz, \
                                     t_xxzz_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xxzz[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.25 * pa_xxzz[j] * fx[j] * fx[j] + pa_xxz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 pa_xzz[j] * fx[j] * fx[j] * pb_x[j] + 4.0 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_zz[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] + pa_x[j] * fx[j] * fx[j] * pb_xzz[j] +

                                 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_xx[j] + fx[j] * fx[j] * pa_z[j] * pb_xxz[j] + 0.5 * pa_xxzz[j] * pb_xx[j] * fx[j] +

                                 0.5 * pa_xxzz[j] * fx[j] * pb_zz[j] + 2.0 * pa_xxz[j] * fx[j] * pb_xxz[j] + 2.0 * pa_xzz[j] * fx[j] * pb_xzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxzz[j] +

                                 pa_xxzz[j] * pb_xxzz[j]) * s_0_0[j];

                t_xxzz_xyyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * pa_xzz[j] * fx[j] * fx[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yyy[j] +

                                 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_xy[j] + 1.5 * pa_xxzz[j] * pb_xy[j] * fx[j] + pa_xzz[j] * fx[j] * pb_yyy[j] +

                                 0.25 * fx[j] * fx[j] * pb_xyyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyyy[j] + 0.5 * fx[j] * pa_zz[j] * pb_xyyy[j] +

                                 pa_xxzz[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (41) = (82,84)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xy, pb_xyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_y, pb_yy, pb_yyz, pb_yz, pb_yzz, pb_z, s_0_0, \
                                     t_xxzz_xyyz, t_xxzz_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xyyz[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 0.5 * pa_xxz[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xzz[j] * fx[j] * fx[j] * pb_z[j] + pa_xz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.125 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yyz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_xz[j] +

                                 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xyy[j] + 0.5 * pa_xxzz[j] * pb_xz[j] * fx[j] + pa_xxz[j] * fx[j] * pb_xyy[j] +

                                 pa_xzz[j] * fx[j] * pb_yyz[j] + 0.25 * fx[j] * fx[j] * pb_xyyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyyz[j] +

                                 0.5 * fx[j] * pa_zz[j] * pb_xyyz[j] + pa_xxzz[j] * pb_xyyz[j]) * s_0_0[j];

                t_xxzz_xyzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.5 * pa_xzz[j] * fx[j] * fx[j] * pb_y[j] + 2.0 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yzz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_xy[j] + fx[j] * fx[j] * pa_z[j] * pb_xyz[j] +

                                 0.5 * pa_xxzz[j] * pb_xy[j] * fx[j] + 2.0 * pa_xxz[j] * fx[j] * pb_xyz[j] + pa_xzz[j] * fx[j] * pb_yzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xyzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xyzz[j] +

                                 pa_xxzz[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (42) = (84,86)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xz, pb_xzz, \
                                     pb_xzzz, pb_yy, pb_yyyy, pb_z, pb_zz, pb_zzz, s_0_0, t_xxzz_xzzz, t_xxzz_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xzzz[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_xxz[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_xzz[j] * fx[j] * fx[j] * pb_z[j] +

                                 3.0 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_zzz[j] +

                                 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_xz[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xzz[j] + 1.5 * pa_xxzz[j] * pb_xz[j] * fx[j] +

                                 3.0 * pa_xxz[j] * fx[j] * pb_xzz[j] + pa_xzz[j] * fx[j] * pb_zzz[j] + 0.25 * fx[j] * fx[j] * pb_xzzz[j] +

                                 0.5 * pa_xx[j] * fx[j] * pb_xzzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xzzz[j] + pa_xxzz[j] * pb_xzzz[j]) * s_0_0[j];

                t_xxzz_yyyy[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] + 0.75 * pa_xxzz[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_yy[j] + 3.0 * pa_xxzz[j] * pb_yy[j] * fx[j] +

                                 0.25 * fx[j] * fx[j] * pb_yyyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyyy[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyyy[j] +

                                 pa_xxzz[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (43) = (86,88)

            #pragma omp simd aligned(fx, pa_xx, pa_xxz, pa_xxzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_z, pb_zz, s_0_0, t_xxzz_yyyz, t_xxzz_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_yyyz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 1.5 * pa_xxz[j] * fx[j] * fx[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_yz[j] +

                                 0.5 * fx[j] * fx[j] * pa_z[j] * pb_yyy[j] + 1.5 * pa_xxzz[j] * pb_yz[j] * fx[j] + pa_xxz[j] * fx[j] * pb_yyy[j] +

                                 0.25 * fx[j] * fx[j] * pb_yyyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyyz[j] +

                                 pa_xxzz[j] * pb_yyyz[j]) * s_0_0[j];

                t_xxzz_yyzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.125 * fx[j] * fx[j] * fx[j] * pa_zz[j] + 0.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * pa_xxzz[j] * fx[j] * fx[j] + pa_xxz[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_yy[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_yy[j] +

                                 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_zz[j] + fx[j] * fx[j] * pa_z[j] * pb_yyz[j] + 0.5 * pa_xxzz[j] * pb_yy[j] * fx[j] +

                                 0.5 * pa_xxzz[j] * fx[j] * pb_zz[j] + 2.0 * pa_xxz[j] * fx[j] * pb_yyz[j] + 0.25 * fx[j] * fx[j] * pb_yyzz[j] +

                                 0.5 * pa_xx[j] * fx[j] * pb_yyzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyzz[j] + pa_xxzz[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (44) = (88,90)

            #pragma omp simd aligned(fx, pa_xx, pa_xxz, pa_xxzz, pa_z, pa_zz, pb_y, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, \
                                     pb_zzz, pb_zzzz, s_0_0, t_xxzz_yzzz, t_xxzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_yzzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xxz[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_yz[j] +

                                 1.5 * fx[j] * fx[j] * pa_z[j] * pb_yzz[j] + 1.5 * pa_xxzz[j] * pb_yz[j] * fx[j] + 3.0 * pa_xxz[j] * fx[j] * pb_yzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_yzzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yzzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yzzz[j] +

                                 pa_xxzz[j] * pb_yzzz[j]) * s_0_0[j];

                t_xxzz_zzzz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 1.875 * pa_xx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] + 3.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_xxzz[j] * fx[j] * fx[j] + 6.0 * pa_xxz[j] * fx[j] * fx[j] * pb_z[j] +

                                 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_zz[j] +

                                 2.0 * fx[j] * fx[j] * pa_z[j] * pb_zzz[j] + 3.0 * pa_xxzz[j] * pb_zz[j] * fx[j] + 4.0 * pa_xxz[j] * fx[j] * pb_zzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_zzzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_zzzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_zzzz[j] +

                                 pa_xxzz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (45) = (90,92)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxy, pb_xy, pb_y, s_0_0, t_xyyy_xxxx, t_xyyy_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xxxx[j] = (1.125 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 4.5 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.75 * pa_xyyy[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_yyy[j] * pb_x[j] +

                                 4.5 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * fx[j] * fx[j] * pa_y[j] * pb_xxx[j] + 3.0 * pa_xyyy[j] * pb_xx[j] * fx[j] +

                                 2.0 * fx[j] * pa_yyy[j] * pb_xxx[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxxx[j] + pa_xyyy[j] * pb_xxxx[j]) * s_0_0[j];

                t_xyyy_xxxy[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_yy[j] +

                                 1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 2.25 * pa_xyy[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_yyy[j] * pb_y[j] + 2.25 * fx[j] * fx[j] * pa_yy[j] * pb_xx[j] +

                                 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xxx[j] +

                                 2.25 * fx[j] * fx[j] * pa_y[j] * pb_xxy[j] + 1.5 * pa_xyyy[j] * pb_xy[j] * fx[j] + 1.5 * pa_xyy[j] * fx[j] * pb_xxx[j] +

                                 1.5 * fx[j] * pa_yyy[j] * pb_xxy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxxy[j] + pa_xyyy[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (46) = (92,94)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxz, pb_xy, pb_xyy, pb_xz, pb_y, pb_yy, pb_z, s_0_0, t_xyyy_xxxz, \
                                     t_xyyy_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xxxz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 0.75 * fx[j] * fx[j] * pa_yyy[j] * pb_z[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] +

                                 2.25 * fx[j] * fx[j] * pa_y[j] * pb_xxz[j] + 1.5 * pa_xyyy[j] * pb_xz[j] * fx[j] + 1.5 * fx[j] * pa_yyy[j] * pb_xxz[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_xxxz[j] + pa_xyyy[j] * pb_xxxz[j]) * s_0_0[j];

                t_xyyy_xxyy[j] = (1.125 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * pa_xyyy[j] * fx[j] * fx[j] + 1.5 * pa_xyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pa_yyy[j] * pb_x[j] +

                                 3.0 * fx[j] * fx[j] * pa_yy[j] * pb_xy[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xxy[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_xyy[j] + 0.5 * pa_xyyy[j] * pb_xx[j] * fx[j] +

                                 0.5 * pa_xyyy[j] * fx[j] * pb_yy[j] + 3.0 * pa_xyy[j] * fx[j] * pb_xxy[j] + fx[j] * pa_yyy[j] * pb_xyy[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_xxyy[j] + pa_xyyy[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (47) = (94,96)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxyz, \
                                     pb_xxz, pb_xxzz, pb_xyz, pb_xz, pb_xzz, pb_yz, pb_z, pb_zz, s_0_0, t_xyyy_xxyz, \
                                     t_xyyy_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xxyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * fx[j] * fx[j] * pa_yy[j] * pb_xz[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xxz[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_xyz[j] + 0.5 * pa_xyyy[j] * fx[j] * pb_yz[j] +

                                 1.5 * pa_xyy[j] * fx[j] * pb_xxz[j] + fx[j] * pa_yyy[j] * pb_xyz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxyz[j] +

                                 pa_xyyy[j] * pb_xxyz[j]) * s_0_0[j];

                t_xyyy_xxzz[j] = (0.375 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.25 * pa_xyyy[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_yyy[j] * pb_x[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * fx[j] * fx[j] * pa_y[j] * pb_xzz[j] + 0.5 * pa_xyyy[j] * pb_xx[j] * fx[j] + 0.5 * pa_xyyy[j] * fx[j] * pb_zz[j] +

                                 fx[j] * pa_yyy[j] * pb_xzz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxzz[j] + pa_xyyy[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (48) = (96,98)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_z, s_0_0, \
                                     t_xyyy_xyyy, t_xyyy_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xyyy[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] +

                                 1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_yy[j] +

                                 3.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 2.25 * pa_xyy[j] * fx[j] * fx[j] * pb_x[j] + 6.75 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * fx[j] * fx[j] * pa_yyy[j] * pb_y[j] + 2.25 * fx[j] * fx[j] * pa_yy[j] * pb_yy[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * pb_xyy[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yyy[j] + 1.5 * pa_xyyy[j] * pb_xy[j] * fx[j] +

                                 4.5 * pa_xyy[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_yyy[j] * pb_yyy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyyy[j] +

                                 pa_xyyy[j] * pb_xyyy[j]) * s_0_0[j];

                t_xyyy_xyyz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.25 * fx[j] * fx[j] * pa_yyy[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_yy[j] * pb_yz[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xyz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yyz[j] + 0.5 * pa_xyyy[j] * pb_xz[j] * fx[j] +

                                 3.0 * pa_xyy[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yyy[j] * pb_yyz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyyz[j] +

                                 pa_xyyy[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (49) = (98,100)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xy, pb_xyzz, pb_xz, \
                                     pb_xzz, pb_xzzz, pb_y, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xyyy_xyzz, t_xyyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xyzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_yy[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.25 * fx[j] * fx[j] * pa_yyy[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_zz[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xzz[j] +

                                 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yzz[j] + 0.5 * pa_xyyy[j] * pb_xy[j] * fx[j] + 1.5 * pa_xyy[j] * fx[j] * pb_xzz[j] +

                                 0.5 * fx[j] * pa_yyy[j] * pb_yzz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyzz[j] + pa_xyyy[j] * pb_xyzz[j]) * s_0_0[j];

                t_xyyy_xzzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 0.75 * fx[j] * fx[j] * pa_yyy[j] * pb_z[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zzz[j] + 1.5 * pa_xyyy[j] * pb_xz[j] * fx[j] + 0.5 * fx[j] * pa_yyy[j] * pb_zzz[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_xzzz[j] + pa_xyyy[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (50) = (100,102)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_xyyy_yyyy, t_xyyy_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_yyyy[j] = (5.625 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 7.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xyyy[j] * fx[j] * fx[j] + 9.0 * pa_xyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 13.5 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_xyyy[j] * pb_yy[j] * fx[j] +

                                 6.0 * pa_xyy[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xy[j] * fx[j] * pb_yyyy[j] + pa_xyyy[j] * pb_yyyy[j]) * s_0_0[j];

                t_xyyy_yyyz[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 2.25 * pa_xyy[j] * fx[j] * fx[j] * pb_z[j] + 6.75 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_xyyy[j] * pb_yz[j] * fx[j] + 4.5 * pa_xyy[j] * fx[j] * pb_yyz[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_yyyz[j] + pa_xyyy[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (51) = (102,104)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pb_y, pb_yy, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, \
                                     pb_z, pb_zz, pb_zzz, s_0_0, t_xyyy_yyzz, t_xyyy_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_yyzz[j] = (1.125 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xyyy[j] * fx[j] * fx[j] + 1.5 * pa_xyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xyyy[j] * pb_yy[j] * fx[j] + 0.5 * pa_xyyy[j] * fx[j] * pb_zz[j] +

                                 3.0 * pa_xyy[j] * fx[j] * pb_yzz[j] + 1.5 * pa_xy[j] * fx[j] * pb_yyzz[j] + pa_xyyy[j] * pb_yyzz[j]) * s_0_0[j];

                t_xyyy_yzzz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 2.25 * pa_xyy[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xyyy[j] * pb_yz[j] * fx[j] + 1.5 * pa_xyy[j] * fx[j] * pb_zzz[j] +

                                 1.5 * pa_xy[j] * fx[j] * pb_yzzz[j] + pa_xyyy[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (52) = (104,106)

            #pragma omp simd aligned(fx, pa_xy, pa_xyyy, pa_xyyz, pa_xz, pa_yyz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_zz, pb_zzzz, s_0_0, t_xyyy_zzzz, t_xyyz_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_zzzz[j] = (1.125 * pa_xy[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyyy[j] * fx[j] * fx[j] +

                                 4.5 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xyyy[j] * pb_zz[j] * fx[j] + 1.5 * pa_xy[j] * fx[j] * pb_zzzz[j] +

                                 pa_xyyy[j] * pb_zzzz[j]) * s_0_0[j];

                t_xyyz_xxxx[j] = (0.375 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.75 * pa_xyyz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_yyz[j] * pb_x[j] +

                                 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] + fx[j] * fx[j] * pa_z[j] * pb_xxx[j] + 3.0 * pa_xyyz[j] * pb_xx[j] * fx[j] +

                                 2.0 * fx[j] * pa_yyz[j] * pb_xxx[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxxx[j] + pa_xyyz[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (53) = (106,108)

            #pragma omp simd aligned(fx, pa_x, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, \
                                     pb_xxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxz, pb_xy, pb_xz, pb_y, pb_z, s_0_0, t_xyyz_xxxy, \
                                     t_xyyz_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xxxy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_yyz[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pa_yz[j] * pb_xx[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xxy[j] + 1.5 * pa_xyyz[j] * pb_xy[j] * fx[j] +

                                 pa_xyz[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_yyz[j] * pb_xxy[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxxy[j] +

                                 pa_xyyz[j] * pb_xxxy[j]) * s_0_0[j];

                t_xyyz_xxxz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_yy[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_yyz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_xx[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xxx[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xxz[j] + 1.5 * pa_xyyz[j] * pb_xz[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * pb_xxx[j] +

                                 1.5 * fx[j] * pa_yyz[j] * pb_xxz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxxz[j] + pa_xyyz[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (54) = (108,110)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, \
                                     pb_x, pb_xx, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, s_0_0, t_xyyz_xxyy, t_xyyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xxyy[j] = (0.375 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.25 * pa_xyyz[j] * fx[j] * fx[j] + pa_xyz[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pa_yyz[j] * pb_x[j] +

                                 2.0 * fx[j] * fx[j] * pa_yz[j] * pb_xy[j] + 0.25 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xyy[j] + 0.5 * pa_xyyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_yy[j] +

                                 2.0 * pa_xyz[j] * fx[j] * pb_xxy[j] + fx[j] * pa_yyz[j] * pb_xyy[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxyy[j] +

                                 pa_xyyz[j] * pb_xxyy[j]) * s_0_0[j];

                t_xyyz_xxyz[j] = (0.25 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.5 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.25 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * pa_xyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.5 * pa_xyz[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.5 * fx[j] * fx[j] * pa_yy[j] * pb_xy[j] + fx[j] * fx[j] * pa_yz[j] * pb_xz[j] + 0.25 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xyz[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_yz[j] +

                                 0.5 * pa_xyy[j] * fx[j] * pb_xxy[j] + pa_xyz[j] * fx[j] * pb_xxz[j] + fx[j] * pa_yyz[j] * pb_xyz[j] +

                                 0.5 * pa_xz[j] * fx[j] * pb_xxyz[j] + pa_xyyz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (55) = (110,112)

            #pragma omp simd aligned(fx, pa_x, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, \
                                     pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_z, pb_zz, \
                                     s_0_0, t_xyyz_xxzz, t_xyyz_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xxzz[j] = (0.125 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 0.5 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * pa_xyyz[j] * fx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.5 * fx[j] * fx[j] * pa_yyz[j] * pb_x[j] + fx[j] * fx[j] * pa_yy[j] * pb_xz[j] + 0.25 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.25 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xxz[j] +

                                 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xzz[j] + 0.5 * pa_xyyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_zz[j] +

                                 pa_xyy[j] * fx[j] * pb_xxz[j] + fx[j] * pa_yyz[j] * pb_xzz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxzz[j] +

                                 pa_xyyz[j] * pb_xxzz[j]) * s_0_0[j];

                t_xyyz_xyyy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_yyz[j] * pb_y[j] +

                                 1.5 * fx[j] * fx[j] * pa_yz[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_yyy[j] + 1.5 * pa_xyyz[j] * pb_xy[j] * fx[j] +

                                 3.0 * pa_xyz[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_yyz[j] * pb_yyy[j] + 0.5 * pa_xz[j] * fx[j] * pb_xyyy[j] +

                                 pa_xyyz[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (56) = (112,114)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, \
                                     pb_x, pb_xy, pb_xyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyz, \
                                     pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, t_xyyz_xyyz, t_xyyz_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xyyz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.125 * fx[j] * fx[j] * fx[j] * pa_yy[j] +

                                 0.5 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.125 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * pa_xyy[j] * fx[j] * fx[j] * pb_x[j] + pa_xy[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pa_yyz[j] * pb_z[j] +

                                 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_yy[j] + fx[j] * fx[j] * pa_yz[j] * pb_yz[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xyy[j] +

                                 0.25 * fx[j] * fx[j] * pa_z[j] * pb_yyz[j] + 0.5 * pa_xyyz[j] * pb_xz[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * pb_xyy[j] +

                                 2.0 * pa_xyz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yyz[j] * pb_yyz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xyyz[j] +

                                 pa_xyyz[j] * pb_xyyz[j]) * s_0_0[j];

                t_xyyz_xyzz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 0.5 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 0.25 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyz[j] * fx[j] * fx[j] * pb_x[j] + pa_xy[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.25 * fx[j] * fx[j] * pa_yyz[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pa_yy[j] * pb_yz[j] +

                                 0.5 * fx[j] * fx[j] * pa_yz[j] * pb_zz[j] + 0.25 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xyz[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_yzz[j] + 0.5 * pa_xyyz[j] * pb_xy[j] * fx[j] +

                                 pa_xyy[j] * fx[j] * pb_xyz[j] + pa_xyz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yyz[j] * pb_yzz[j] +

                                 0.5 * pa_xz[j] * fx[j] * pb_xyzz[j] + pa_xyyz[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (57) = (114,116)

            #pragma omp simd aligned(fx, pa_x, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_yy, pa_yyz, pa_z, pb_x, pb_xz, pb_xzz, \
                                     pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_z, pb_zz, pb_zzz, s_0_0, t_xyyz_xzzz, \
                                     t_xyyz_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xzzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_yy[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_yyz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_zz[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xzz[j] +

                                 0.25 * fx[j] * fx[j] * pa_z[j] * pb_zzz[j] + 1.5 * pa_xyyz[j] * pb_xz[j] * fx[j] + 1.5 * pa_xyy[j] * fx[j] * pb_xzz[j] +

                                 0.5 * fx[j] * pa_yyz[j] * pb_zzz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xzzz[j] + pa_xyyz[j] * pb_xzzz[j]) * s_0_0[j];

                t_xyyz_yyyy[j] = (1.875 * pa_xz[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyyz[j] * fx[j] * fx[j] +

                                 6.0 * pa_xyz[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xyyz[j] * pb_yy[j] * fx[j] +

                                 4.0 * pa_xyz[j] * fx[j] * pb_yyy[j] + 0.5 * pa_xz[j] * fx[j] * pb_yyyy[j] + pa_xyyz[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (58) = (116,118)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, t_xyyz_yyyz, t_xyyz_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_yyyz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] +

                                 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xyyz[j] * pb_yz[j] * fx[j] +

                                 0.5 * pa_xyy[j] * fx[j] * pb_yyy[j] + 3.0 * pa_xyz[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xz[j] * fx[j] * pb_yyyz[j] +

                                 pa_xyyz[j] * pb_yyyz[j]) * s_0_0[j];

                t_xyyz_yyzz[j] = (0.375 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_xyyz[j] * fx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * fx[j] * pb_z[j] +

                                 pa_xyz[j] * fx[j] * fx[j] * pb_y[j] + 2.0 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xyyz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_zz[j] +

                                 pa_xyy[j] * fx[j] * pb_yyz[j] + 2.0 * pa_xyz[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xz[j] * fx[j] * pb_yyzz[j] +

                                 pa_xyyz[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (59) = (118,120)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pb_y, pb_yz, pb_yzz, pb_yzzz, \
                                     pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xyyz_yzzz, t_xyyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_yzzz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_xyyz[j] * pb_yz[j] * fx[j] +

                                 1.5 * pa_xyy[j] * fx[j] * pb_yzz[j] + pa_xyz[j] * fx[j] * pb_zzz[j] + 0.5 * pa_xz[j] * fx[j] * pb_yzzz[j] +

                                 pa_xyyz[j] * pb_yzzz[j]) * s_0_0[j];

                t_xyyz_zzzz[j] = (0.375 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xyyz[j] * fx[j] * fx[j] + 3.0 * pa_xyy[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] + pa_x[j] * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_xyyz[j] * pb_zz[j] * fx[j] +

                                 2.0 * pa_xyy[j] * fx[j] * pb_zzz[j] + 0.5 * pa_xz[j] * fx[j] * pb_zzzz[j] + pa_xyyz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (60) = (120,122)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyzz, pa_xzz, pa_y, pa_yzz, pa_zz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxy, pb_xy, pb_y, s_0_0, t_xyzz_xxxx, t_xyzz_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xxxx[j] = (0.375 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.75 * pa_xyzz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_yzz[j] * pb_x[j] +

                                 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] + fx[j] * fx[j] * pa_y[j] * pb_xxx[j] + 3.0 * pa_xyzz[j] * pb_xx[j] * fx[j] +

                                 2.0 * fx[j] * pa_yzz[j] * pb_xxx[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxxx[j] + pa_xyzz[j] * pb_xxxx[j]) * s_0_0[j];

                t_xyzz_xxxy[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_yzz[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_xx[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xxx[j] +

                                 0.75 * fx[j] * fx[j] * pa_y[j] * pb_xxy[j] + 1.5 * pa_xyzz[j] * pb_xy[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] * pb_xxx[j] +

                                 1.5 * fx[j] * pa_yzz[j] * pb_xxy[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxxy[j] + pa_xyzz[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (61) = (122,124)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xzz, pa_y, pa_yz, pa_yzz, pa_zz, pb_x, pb_xx, \
                                     pb_xxx, pb_xxxz, pb_xxy, pb_xxyy, pb_xxz, pb_xy, pb_xyy, pb_xz, pb_y, pb_yy, pb_z, s_0_0, \
                                     t_xyzz_xxxz, t_xyzz_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xxxz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_yzz[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_yz[j] * pb_xx[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_xxz[j] + 1.5 * pa_xyzz[j] * pb_xz[j] * fx[j] +

                                 pa_xyz[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_yzz[j] * pb_xxz[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxxz[j] +

                                 pa_xyzz[j] * pb_xxxz[j]) * s_0_0[j];

                t_xyzz_xxyy[j] = (0.125 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] +

                                 0.5 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * pa_xyzz[j] * fx[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.5 * fx[j] * fx[j] * pa_yzz[j] * pb_x[j] + fx[j] * fx[j] * pa_zz[j] * pb_xy[j] + 0.25 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.25 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xxy[j] +

                                 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xyy[j] + 0.5 * pa_xyzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_yy[j] +

                                 pa_xzz[j] * fx[j] * pb_xxy[j] + fx[j] * pa_yzz[j] * pb_xyy[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxyy[j] +

                                 pa_xyzz[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (62) = (124,126)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pb_x, pb_xx, pb_xxy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyz, pb_xz, pb_xzz, pb_y, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_xyzz_xxyz, t_xyzz_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xxyz[j] = (0.25 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.25 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xyz[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.25 * pa_xzz[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] + fx[j] * fx[j] * pa_yz[j] * pb_xy[j] +

                                 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_xz[j] + 0.25 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xxz[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xyz[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_yz[j] +

                                 pa_xyz[j] * fx[j] * pb_xxy[j] + 0.5 * pa_xzz[j] * fx[j] * pb_xxz[j] + fx[j] * pa_yzz[j] * pb_xyz[j] +

                                 0.5 * pa_xy[j] * fx[j] * pb_xxyz[j] + pa_xyzz[j] * pb_xxyz[j]) * s_0_0[j];

                t_xyzz_xxzz[j] = (0.375 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_x[j] + 0.25 * pa_xyzz[j] * fx[j] * fx[j] + pa_xyz[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pa_yzz[j] * pb_x[j] +

                                 2.0 * fx[j] * fx[j] * pa_yz[j] * pb_xz[j] + 0.25 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xzz[j] + 0.5 * pa_xyzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_zz[j] +

                                 2.0 * pa_xyz[j] * fx[j] * pb_xxz[j] + fx[j] * pa_yzz[j] * pb_xzz[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxzz[j] +

                                 pa_xyzz[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (63) = (126,128)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pb_x, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_xyzz_xyyy, t_xyzz_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xyyy[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_yzz[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_yy[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xyy[j] +

                                 0.25 * fx[j] * fx[j] * pa_y[j] * pb_yyy[j] + 1.5 * pa_xyzz[j] * pb_xy[j] * fx[j] + 1.5 * pa_xzz[j] * fx[j] * pb_xyy[j] +

                                 0.5 * fx[j] * pa_yzz[j] * pb_yyy[j] + 0.5 * pa_xy[j] * fx[j] * pb_xyyy[j] + pa_xyzz[j] * pb_xyyy[j]) * s_0_0[j];

                t_xyzz_xyyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 0.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] +

                                 0.25 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyz[j] * fx[j] * fx[j] * pb_x[j] + pa_xz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.25 * fx[j] * fx[j] * pa_yzz[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pa_yz[j] * pb_yy[j] +

                                 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_yz[j] + 0.25 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xyz[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_yyz[j] + 0.5 * pa_xyzz[j] * pb_xz[j] * fx[j] +

                                 pa_xyz[j] * fx[j] * pb_xyy[j] + pa_xzz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yzz[j] * pb_yyz[j] +

                                 0.5 * pa_xy[j] * fx[j] * pb_xyyz[j] + pa_xyzz[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (64) = (128,130)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pb_x, pb_xy, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, s_0_0, t_xyzz_xyzz, t_xyzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xyzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_y[j] +

                                 0.125 * fx[j] * fx[j] * fx[j] * pa_zz[j] + 0.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.125 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.25 * pa_xzz[j] * fx[j] * fx[j] * pb_x[j] + pa_xz[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pa_yzz[j] * pb_y[j] +

                                 fx[j] * fx[j] * pa_yz[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_zz[j] +

                                 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xzz[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_yzz[j] + 0.5 * pa_xyzz[j] * pb_xy[j] * fx[j] +

                                 2.0 * pa_xyz[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xzz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yzz[j] * pb_yzz[j] +

                                 0.5 * pa_xy[j] * fx[j] * pb_xyzz[j] + pa_xyzz[j] * pb_xyzz[j]) * s_0_0[j];

                t_xyzz_xzzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_yz[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] * pb_z[j] + 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_yzz[j] * pb_z[j] +

                                 1.5 * fx[j] * fx[j] * pa_yz[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_zzz[j] + 1.5 * pa_xyzz[j] * pb_xz[j] * fx[j] +

                                 3.0 * pa_xyz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yzz[j] * pb_zzz[j] + 0.5 * pa_xy[j] * fx[j] * pb_xzzz[j] +

                                 pa_xyzz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (65) = (130,132)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yz, pb_z, s_0_0, t_xyzz_yyyy, t_xyzz_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_yyyy[j] = (0.375 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xyzz[j] * fx[j] * fx[j] + 3.0 * pa_xzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] + pa_x[j] * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_xyzz[j] * pb_yy[j] * fx[j] +

                                 2.0 * pa_xzz[j] * fx[j] * pb_yyy[j] + 0.5 * pa_xy[j] * fx[j] * pb_yyyy[j] + pa_xyzz[j] * pb_yyyy[j]) * s_0_0[j];

                t_xyzz_yyyz[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * pa_xzz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_xyzz[j] * pb_yz[j] * fx[j] +

                                 pa_xyz[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xzz[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xy[j] * fx[j] * pb_yyyz[j] +

                                 pa_xyzz[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (66) = (132,134)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pb_y, pb_yy, pb_yyz, pb_yyzz, \
                                     pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xyzz_yyzz, t_xyzz_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_yyzz[j] = (0.375 * pa_xy[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xyzz[j] * fx[j] * fx[j] + pa_xyz[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.0 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xyzz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_zz[j] +

                                 2.0 * pa_xyz[j] * fx[j] * pb_yyz[j] + pa_xzz[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xy[j] * fx[j] * pb_yyzz[j] +

                                 pa_xyzz[j] * pb_yyzz[j]) * s_0_0[j];

                t_xyzz_yzzz[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xyz[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xyzz[j] * pb_yz[j] * fx[j] +

                                 3.0 * pa_xyz[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xzz[j] * fx[j] * pb_zzz[j] + 0.5 * pa_xy[j] * fx[j] * pb_yzzz[j] +

                                 pa_xyzz[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (67) = (134,136)

            #pragma omp simd aligned(fx, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzzz, pa_z, pa_zzz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xyzz_zzzz, t_xzzz_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_zzzz[j] = (1.875 * pa_xy[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyzz[j] * fx[j] * fx[j] +

                                 6.0 * pa_xyz[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_xy[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xyzz[j] * pb_zz[j] * fx[j] +

                                 4.0 * pa_xyz[j] * fx[j] * pb_zzz[j] + 0.5 * pa_xy[j] * fx[j] * pb_zzzz[j] + pa_xyzz[j] * pb_zzzz[j]) * s_0_0[j];

                t_xzzz_xxxx[j] = (1.125 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 4.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.75 * pa_xzzz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_zzz[j] * pb_x[j] +

                                 4.5 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * fx[j] * fx[j] * pa_z[j] * pb_xxx[j] + 3.0 * pa_xzzz[j] * pb_xx[j] * fx[j] +

                                 2.0 * fx[j] * pa_zzz[j] * pb_xxx[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxxx[j] + pa_xzzz[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (68) = (136,138)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxz, pb_xy, pb_xz, pb_y, pb_z, s_0_0, t_xzzz_xxxy, \
                                     t_xzzz_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xxxy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 0.75 * fx[j] * fx[j] * pa_zzz[j] * pb_y[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xxy[j] + 1.5 * pa_xzzz[j] * pb_xy[j] * fx[j] + 1.5 * fx[j] * pa_zzz[j] * pb_xxy[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_xxxy[j] + pa_xzzz[j] * pb_xxxy[j]) * s_0_0[j];

                t_xzzz_xxxz[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 2.25 * pa_xzz[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_zzz[j] * pb_z[j] + 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_xx[j] +

                                 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xxx[j] +

                                 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xxz[j] + 1.5 * pa_xzzz[j] * pb_xz[j] * fx[j] + 1.5 * pa_xzz[j] * fx[j] * pb_xxx[j] +

                                 1.5 * fx[j] * pa_zzz[j] * pb_xxz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxxz[j] + pa_xzzz[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (69) = (138,140)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xy, pb_xyy, pb_xyz, pb_y, pb_yy, pb_yz, s_0_0, t_xzzz_xxyy, \
                                     t_xzzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xxyy[j] = (0.375 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.25 * pa_xzzz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_zzz[j] * pb_x[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xyy[j] + 0.5 * pa_xzzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xzzz[j] * fx[j] * pb_yy[j] +

                                 fx[j] * pa_zzz[j] * pb_xyy[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxyy[j] + pa_xzzz[j] * pb_xxyy[j]) * s_0_0[j];

                t_xzzz_xxyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_xy[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xxy[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xyz[j] + 0.5 * pa_xzzz[j] * fx[j] * pb_yz[j] +

                                 1.5 * pa_xzz[j] * fx[j] * pb_xxy[j] + fx[j] * pa_zzz[j] * pb_xyz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxyz[j] +

                                 pa_xzzz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (70) = (140,142)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyyy, pb_xz, pb_xzz, pb_y, pb_yyy, pb_z, pb_zz, s_0_0, t_xzzz_xxzz, \
                                     t_xzzz_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xxzz[j] = (1.125 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * pa_xzzz[j] * fx[j] * fx[j] + 1.5 * pa_xzz[j] * fx[j] * fx[j] * pb_z[j] +

                                 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pa_zzz[j] * pb_x[j] +

                                 3.0 * fx[j] * fx[j] * pa_zz[j] * pb_xz[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xxz[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xzz[j] + 0.5 * pa_xzzz[j] * pb_xx[j] * fx[j] +

                                 0.5 * pa_xzzz[j] * fx[j] * pb_zz[j] + 3.0 * pa_xzz[j] * fx[j] * pb_xxz[j] + fx[j] * pa_zzz[j] * pb_xzz[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_xxzz[j] + pa_xzzz[j] * pb_xxzz[j]) * s_0_0[j];

                t_xzzz_xyyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 0.75 * fx[j] * fx[j] * pa_zzz[j] * pb_y[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yyy[j] + 1.5 * pa_xzzz[j] * pb_xy[j] * fx[j] + 0.5 * fx[j] * pa_zzz[j] * pb_yyy[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_xyyy[j] + pa_xzzz[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (71) = (142,144)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xy, pb_xyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_y, pb_yy, pb_yyz, pb_yz, pb_yzz, pb_z, s_0_0, \
                                     t_xzzz_xyyz, t_xzzz_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xyyz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.25 * fx[j] * fx[j] * pa_zzz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_yy[j] +

                                 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xyy[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yyz[j] + 0.5 * pa_xzzz[j] * pb_xz[j] * fx[j] + 1.5 * pa_xzz[j] * fx[j] * pb_xyy[j] +

                                 0.5 * fx[j] * pa_zzz[j] * pb_yyz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyyz[j] + pa_xzzz[j] * pb_xyyz[j]) * s_0_0[j];

                t_xzzz_xyzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.25 * fx[j] * fx[j] * pa_zzz[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_yz[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xyz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yzz[j] + 0.5 * pa_xzzz[j] * pb_xy[j] * fx[j] +

                                 3.0 * pa_xzz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zzz[j] * pb_yzz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyzz[j] +

                                 pa_xzzz[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (72) = (144,146)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xz, pb_xzz, \
                                     pb_xzzz, pb_yy, pb_yyyy, pb_z, pb_zz, pb_zzz, s_0_0, t_xzzz_xzzz, t_xzzz_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xzzz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] +

                                 1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 3.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 2.25 * pa_xzz[j] * fx[j] * fx[j] * pb_x[j] + 6.75 * pa_xz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * fx[j] * fx[j] * pa_zzz[j] * pb_z[j] + 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_zz[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * pb_xzz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zzz[j] + 1.5 * pa_xzzz[j] * pb_xz[j] * fx[j] +

                                 4.5 * pa_xzz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_zzz[j] * pb_zzz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xzzz[j] +

                                 pa_xzzz[j] * pb_xzzz[j]) * s_0_0[j];

                t_xzzz_yyyy[j] = (1.125 * pa_xz[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xzzz[j] * fx[j] * fx[j] +

                                 4.5 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xzzz[j] * pb_yy[j] * fx[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyyy[j] +

                                 pa_xzzz[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (73) = (146,148)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyz, pb_yyz, pb_yyzz, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_xzzz_yyyz, t_xzzz_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_yyyz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_xzz[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xzzz[j] * pb_yz[j] * fx[j] + 1.5 * pa_xzz[j] * fx[j] * pb_yyy[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_yyyz[j] + pa_xzzz[j] * pb_yyyz[j]) * s_0_0[j];

                t_xzzz_yyzz[j] = (1.125 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_xzzz[j] * fx[j] * fx[j] + 1.5 * pa_xzz[j] * fx[j] * fx[j] * pb_z[j] +

                                 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_x[j] * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xzzz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xzzz[j] * fx[j] * pb_zz[j] +

                                 3.0 * pa_xzz[j] * fx[j] * pb_yyz[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyzz[j] + pa_xzzz[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (74) = (148,150)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pb_y, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, \
                                     pb_zzz, pb_zzzz, s_0_0, t_xzzz_yzzz, t_xzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_yzzz[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_xzz[j] * fx[j] * fx[j] * pb_y[j] + 6.75 * pa_xz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 2.25 * pa_x[j] * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_xzzz[j] * pb_yz[j] * fx[j] + 4.5 * pa_xzz[j] * fx[j] * pb_yzz[j] +

                                 1.5 * pa_xz[j] * fx[j] * pb_yzzz[j] + pa_xzzz[j] * pb_yzzz[j]) * s_0_0[j];

                t_xzzz_zzzz[j] = (5.625 * pa_xz[j] * fx[j] * fx[j] * fx[j] +

                                 7.5 * pa_x[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xzzz[j] * fx[j] * fx[j] + 9.0 * pa_xzz[j] * fx[j] * fx[j] * pb_z[j] +

                                 13.5 * pa_xz[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_xzzz[j] * pb_zz[j] * fx[j] +

                                 6.0 * pa_xzz[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xz[j] * fx[j] * pb_zzzz[j] + pa_xzzz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (75) = (150,152)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xy, \
                                     s_0_0, t_yyyy_xxxx, t_yyyy_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xxxx[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_yyyy[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 9.0 * pa_yy[j] * fx[j] * fx[j] * pb_xx[j] +

                                 3.0 * pa_yyyy[j] * pb_xx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xxxx[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxxx[j] +

                                 pa_yyyy[j] * pb_xxxx[j]) * s_0_0[j];

                t_yyyy_xxxy[j] = (4.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 3.0 * pa_yyy[j] * fx[j] * fx[j] * pb_x[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_xy[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_yyyy[j] * pb_xy[j] * fx[j] +

                                 2.0 * pa_yyy[j] * fx[j] * pb_xxx[j] + 0.75 * fx[j] * fx[j] * pb_xxxy[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxxy[j] +

                                 pa_yyyy[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (76) = (152,154)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_xx, pb_xxxz, pb_xxy, pb_xxyy, pb_xz, pb_y, \
                                     pb_yy, s_0_0, t_yyyy_xxxz, t_yyyy_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xxxz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yyyy[j] * pb_xz[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xxxz[j] +

                                 3.0 * pa_yy[j] * fx[j] * pb_xxxz[j] + pa_yyyy[j] * pb_xxxz[j]) * s_0_0[j];

                t_yyyy_xxyy[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 3.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.875 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_yyyy[j] * fx[j] * fx[j] +

                                 2.0 * pa_yyy[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_yy[j] +

                                 6.0 * pa_y[j] * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_yyyy[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyyy[j] * fx[j] * pb_yy[j] +

                                 4.0 * pa_yyy[j] * fx[j] * pb_xxy[j] + 0.75 * fx[j] * fx[j] * pb_xxyy[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxyy[j] +

                                 pa_yyyy[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (77) = (154,156)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_xx, pb_xxyz, pb_xxz, pb_xxzz, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_yyyy_xxyz, t_yyyy_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xxyz[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 pa_yyy[j] * fx[j] * fx[j] * pb_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_yz[j] +

                                 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_yz[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_yyyy[j] * fx[j] * pb_yz[j] +

                                 2.0 * pa_yyy[j] * fx[j] * pb_xxz[j] + 0.75 * fx[j] * fx[j] * pb_xxyz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxyz[j] +

                                 pa_yyyy[j] * pb_xxyz[j]) * s_0_0[j];

                t_yyyy_xxzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 0.25 * pa_yyyy[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_yyyy[j] * pb_xx[j] * fx[j] +

                                 0.5 * pa_yyyy[j] * fx[j] * pb_zz[j] + 0.75 * fx[j] * fx[j] * pb_xxzz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxzz[j] +

                                 pa_yyyy[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (78) = (156,158)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xz, s_0_0, t_yyyy_xyyy, t_yyyy_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xyyy[j] = (7.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 5.625 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 3.0 * pa_yyy[j] * fx[j] * fx[j] * pb_x[j] +

                                 13.5 * pa_yy[j] * fx[j] * fx[j] * pb_xy[j] + 9.0 * pa_y[j] * fx[j] * fx[j] * pb_xyy[j] + 1.5 * pa_yyyy[j] * pb_xy[j] * fx[j] +

                                 6.0 * pa_yyy[j] * fx[j] * pb_xyy[j] + 0.75 * fx[j] * fx[j] * pb_xyyy[j] + 3.0 * pa_yy[j] * fx[j] * pb_xyyy[j] +

                                 pa_yyyy[j] * pb_xyyy[j]) * s_0_0[j];

                t_yyyy_xyyz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_xz[j] + 6.0 * pa_y[j] * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_yyyy[j] * pb_xz[j] * fx[j] +

                                 4.0 * pa_yyy[j] * fx[j] * pb_xyz[j] + 0.75 * fx[j] * fx[j] * pb_xyyz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xyyz[j] +

                                 pa_yyyy[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (79) = (158,160)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xy, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, \
                                     s_0_0, t_yyyy_xyzz, t_yyyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xyzz[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 pa_yyy[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_xy[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_yyyy[j] * pb_xy[j] * fx[j] +

                                 2.0 * pa_yyy[j] * fx[j] * pb_xzz[j] + 0.75 * fx[j] * fx[j] * pb_xyzz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xyzz[j] +

                                 pa_yyyy[j] * pb_xyzz[j]) * s_0_0[j];

                t_yyyy_xzzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yyyy[j] * pb_xz[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xzzz[j] +

                                 3.0 * pa_yy[j] * fx[j] * pb_xzzz[j] + pa_yyyy[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (80) = (160,162)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yz, pb_z, s_0_0, t_yyyy_yyyy, t_yyyy_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_yyyy[j] = (6.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 11.25 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 30.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 11.25 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_yyyy[j] * fx[j] * fx[j] +

                                 12.0 * pa_yyy[j] * fx[j] * fx[j] * pb_y[j] + 27.0 * pa_yy[j] * fx[j] * fx[j] * pb_yy[j] +

                                 12.0 * pa_y[j] * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yyyy[j] * pb_yy[j] * fx[j] + 8.0 * pa_yyy[j] * fx[j] * pb_yyy[j] +

                                 0.75 * fx[j] * fx[j] * pb_yyyy[j] + 3.0 * pa_yy[j] * fx[j] * pb_yyyy[j] + pa_yyyy[j] * pb_yyyy[j]) * s_0_0[j];

                t_yyyy_yyyz[j] = (7.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 5.625 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 3.0 * pa_yyy[j] * fx[j] * fx[j] * pb_z[j] +

                                 13.5 * pa_yy[j] * fx[j] * fx[j] * pb_yz[j] + 9.0 * pa_y[j] * fx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_yyyy[j] * pb_yz[j] * fx[j] +

                                 6.0 * pa_yyy[j] * fx[j] * pb_yyz[j] + 0.75 * fx[j] * fx[j] * pb_yyyz[j] + 3.0 * pa_yy[j] * fx[j] * pb_yyyz[j] +

                                 pa_yyyy[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (81) = (162,164)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_y, pb_yy, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, \
                                     pb_z, pb_zz, pb_zzz, s_0_0, t_yyyy_yyzz, t_yyyy_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_yyzz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 3.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.875 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * pa_yyyy[j] * fx[j] * fx[j] +

                                 2.0 * pa_yyy[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_yy[j] +

                                 6.0 * pa_y[j] * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_yyyy[j] * pb_yy[j] * fx[j] + 0.5 * pa_yyyy[j] * fx[j] * pb_zz[j] +

                                 4.0 * pa_yyy[j] * fx[j] * pb_yzz[j] + 0.75 * fx[j] * fx[j] * pb_yyzz[j] + 3.0 * pa_yy[j] * fx[j] * pb_yyzz[j] +

                                 pa_yyyy[j] * pb_yyzz[j]) * s_0_0[j];

                t_yyyy_yzzz[j] = (4.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 3.0 * pa_yyy[j] * fx[j] * fx[j] * pb_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_yz[j] +

                                 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_yz[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_yyyy[j] * pb_yz[j] * fx[j] +

                                 2.0 * pa_yyy[j] * fx[j] * pb_zzz[j] + 0.75 * fx[j] * fx[j] * pb_yzzz[j] + 3.0 * pa_yy[j] * fx[j] * pb_yzzz[j] +

                                 pa_yyyy[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (82) = (164,166)

            #pragma omp simd aligned(fx, pa_yy, pa_yyyy, pa_yyyz, pa_yz, pb_xx, pb_xxxx, pb_zz, pb_zzzz, s_0_0, \
                                     t_yyyy_zzzz, t_yyyz_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_zzzz[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_yyyy[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 9.0 * pa_yy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 3.0 * pa_yyyy[j] * pb_zz[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_zzzz[j] + 3.0 * pa_yy[j] * fx[j] * pb_zzzz[j] +

                                 pa_yyyy[j] * pb_zzzz[j]) * s_0_0[j];

                t_yyyz_xxxx[j] = (1.125 * pa_yz[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yyyz[j] * fx[j] * fx[j] +

                                 4.5 * pa_yz[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_yyyz[j] * pb_xx[j] * fx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxxx[j] +

                                 pa_yyyz[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (83) = (166,168)

            #pragma omp simd aligned(fx, pa_y, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xxx, pb_xxxy, pb_xxxz, \
                                     pb_xy, pb_xz, s_0_0, t_yyyz_xxxy, t_yyyz_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xxxy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 2.25 * pa_yyz[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xxx[j] + 1.5 * pa_yyyz[j] * pb_xy[j] * fx[j] + 1.5 * pa_yyz[j] * fx[j] * pb_xxx[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_xxxy[j] + pa_yyyz[j] * pb_xxxy[j]) * s_0_0[j];

                t_yyyz_xxxz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * pa_yyy[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_yyyz[j] * pb_xz[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] * pb_xxx[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_xxxz[j] + pa_yyyz[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (84) = (168,170)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_xx, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_yyyz_xxyy, t_yyyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xxyy[j] = (1.125 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.25 * pa_yyyz[j] * fx[j] * fx[j] + 1.5 * pa_yyz[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xxy[j] + 0.5 * pa_yyyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyyz[j] * fx[j] * pb_yy[j] +

                                 3.0 * pa_yyz[j] * fx[j] * pb_xxy[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxyy[j] + pa_yyyz[j] * pb_xxyy[j]) * s_0_0[j];

                t_yyyz_xxyz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_yyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * pa_yyz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xxy[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xxz[j] + 0.5 * pa_yyyz[j] * fx[j] * pb_yz[j] + 0.5 * pa_yyy[j] * fx[j] * pb_xxy[j] +

                                 1.5 * pa_yyz[j] * fx[j] * pb_xxz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxyz[j] + pa_yyyz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (85) = (170,172)

            #pragma omp simd aligned(fx, pa_y, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_z, pb_zz, s_0_0, t_yyyz_xxzz, t_yyyz_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xxzz[j] = (0.375 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_yyyz[j] * fx[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_y[j] * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_yyyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyyz[j] * fx[j] * pb_zz[j] +

                                 pa_yyy[j] * fx[j] * pb_xxz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxzz[j] + pa_yyyz[j] * pb_xxzz[j]) * s_0_0[j];

                t_yyyz_xyyy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 2.25 * pa_yyz[j] * fx[j] * fx[j] * pb_x[j] + 6.75 * pa_yz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xyy[j] + 1.5 * pa_yyyz[j] * pb_xy[j] * fx[j] + 4.5 * pa_yyz[j] * fx[j] * pb_xyy[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_xyyy[j] + pa_yyyz[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (86) = (172,174)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xy, pb_xyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, s_0_0, t_yyyz_xyyz, t_yyyz_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xyyz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * pa_yyy[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_xy[j] + 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xyy[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xyz[j] + 0.5 * pa_yyyz[j] * pb_xz[j] * fx[j] +

                                 0.5 * pa_yyy[j] * fx[j] * pb_xyy[j] + 3.0 * pa_yyz[j] * fx[j] * pb_xyz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xyyz[j] +

                                 pa_yyyz[j] * pb_xyyz[j]) * s_0_0[j];

                t_yyyz_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * pa_yyz[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 1.5 * pa_y[j] * fx[j] * fx[j] * pb_xyz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xzz[j] + 0.5 * pa_yyyz[j] * pb_xy[j] * fx[j] +

                                 pa_yyy[j] * fx[j] * pb_xyz[j] + 1.5 * pa_yyz[j] * fx[j] * pb_xzz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xyzz[j] +

                                 pa_yyyz[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (87) = (174,176)

            #pragma omp simd aligned(fx, pa_y, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyyy, s_0_0, t_yyyz_xzzz, t_yyyz_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xzzz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * pa_yyy[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 2.25 * pa_y[j] * fx[j] * fx[j] * pb_xzz[j] + 1.5 * pa_yyyz[j] * pb_xz[j] * fx[j] + 1.5 * pa_yyy[j] * fx[j] * pb_xzz[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_xzzz[j] + pa_yyyz[j] * pb_xzzz[j]) * s_0_0[j];

                t_yyyz_yyyy[j] = (5.625 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 7.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.75 * pa_yyyz[j] * fx[j] * fx[j] + 9.0 * pa_yyz[j] * fx[j] * fx[j] * pb_y[j] +

                                 13.5 * pa_yz[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * fx[j] * fx[j] * pa_z[j] * pb_yyy[j] + 3.0 * pa_yyyz[j] * pb_yy[j] * fx[j] +

                                 6.0 * pa_yyz[j] * fx[j] * pb_yyy[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyyy[j] + pa_yyyz[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (88) = (176,178)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, t_yyyz_yyyz, t_yyyz_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_yyyz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 1.125 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 3.375 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_yyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_yyz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_yy[j] +

                                 6.75 * pa_yz[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yyy[j] +

                                 2.25 * fx[j] * fx[j] * pa_z[j] * pb_yyz[j] + 1.5 * pa_yyyz[j] * pb_yz[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] * pb_yyy[j] +

                                 4.5 * pa_yyz[j] * fx[j] * pb_yyz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyyz[j] + pa_yyyz[j] * pb_yyyz[j]) * s_0_0[j];

                t_yyyz_yyzz[j] = (1.125 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * pa_yyyz[j] * fx[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * pa_yyz[j] * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_yy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yyz[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_yzz[j] + 0.5 * pa_yyyz[j] * pb_yy[j] * fx[j] +

                                 0.5 * pa_yyyz[j] * fx[j] * pb_zz[j] + pa_yyy[j] * fx[j] * pb_yyz[j] + 3.0 * pa_yyz[j] * fx[j] * pb_yzz[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_yyzz[j] + pa_yyyz[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (89) = (178,180)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_y, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_yyyz_yzzz, t_yyyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_yzzz[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 1.125 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_yyy[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_yyz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_zz[j] +

                                 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_yz[j] + 2.25 * pa_y[j] * fx[j] * fx[j] * pb_yzz[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zzz[j] + 1.5 * pa_yyyz[j] * pb_yz[j] * fx[j] + 1.5 * pa_yyy[j] * fx[j] * pb_yzz[j] +

                                 1.5 * pa_yyz[j] * fx[j] * pb_zzz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yzzz[j] + pa_yyyz[j] * pb_yzzz[j]) * s_0_0[j];

                t_yyyz_zzzz[j] = (1.125 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 4.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_yyyz[j] * fx[j] * fx[j] + 3.0 * pa_yyy[j] * fx[j] * fx[j] * pb_z[j] +

                                 4.5 * pa_yz[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_yyyz[j] * pb_zz[j] * fx[j] +

                                 2.0 * pa_yyy[j] * fx[j] * pb_zzz[j] + 1.5 * pa_yz[j] * fx[j] * pb_zzzz[j] + pa_yyyz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (90) = (180,182)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyzz, pa_yzz, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xy, s_0_0, t_yyzz_xxxx, t_yyzz_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xxxx[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] + 0.75 * pa_yyzz[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_xx[j] +

                                 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_xx[j] + 3.0 * pa_yyzz[j] * pb_xx[j] * fx[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxxx[j] + 0.5 * pa_yy[j] * fx[j] * pb_xxxx[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxxx[j] +

                                 pa_yyzz[j] * pb_xxxx[j]) * s_0_0[j];

                t_yyzz_xxxy[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * pa_yzz[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * pb_xxx[j] +

                                 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_xy[j] + 1.5 * pa_yyzz[j] * pb_xy[j] * fx[j] + pa_yzz[j] * fx[j] * pb_xxx[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxxy[j] + 0.5 * pa_yy[j] * fx[j] * pb_xxxy[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxxy[j] +

                                 pa_yyzz[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (91) = (182,184)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xz, pb_y, pb_yy, s_0_0, t_yyzz_xxxz, t_yyzz_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xxxz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 1.5 * pa_yyz[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_xz[j] +

                                 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xxx[j] + 1.5 * pa_yyzz[j] * pb_xz[j] * fx[j] + pa_yyz[j] * fx[j] * pb_xxx[j] +

                                 0.25 * fx[j] * fx[j] * pb_xxxz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xxxz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxxz[j] +

                                 pa_yyzz[j] * pb_xxxz[j]) * s_0_0[j];

                t_yyzz_xxyy[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 0.125 * pa_yy[j] * fx[j] * fx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_yyzz[j] * fx[j] * fx[j] + pa_yzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_xx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_yy[j] + pa_y[j] * fx[j] * fx[j] * pb_xxy[j] +

                                 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_yy[j] + 0.5 * pa_yyzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyzz[j] * fx[j] * pb_yy[j] +

                                 2.0 * pa_yzz[j] * fx[j] * pb_xxy[j] + 0.25 * fx[j] * fx[j] * pb_xxyy[j] + 0.5 * pa_yy[j] * fx[j] * pb_xxyy[j] +

                                 0.5 * fx[j] * pa_zz[j] * pb_xxyy[j] + pa_yyzz[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (92) = (184,186)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_xx, pb_xxy, \
                                     pb_xxyz, pb_xxz, pb_xxzz, pb_y, pb_yz, pb_z, pb_zz, s_0_0, t_yyzz_xxyz, t_yyzz_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xxyz[j] = (0.5 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 0.25 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 0.5 * pa_yyz[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_yzz[j] * fx[j] * fx[j] * pb_z[j] + pa_yz[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.125 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.5 * pa_y[j] * fx[j] * fx[j] * pb_xxz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_yz[j] +

                                 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xxy[j] + 0.5 * pa_yyzz[j] * fx[j] * pb_yz[j] + pa_yyz[j] * fx[j] * pb_xxy[j] +

                                 pa_yzz[j] * fx[j] * pb_xxz[j] + 0.25 * fx[j] * fx[j] * pb_xxyz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xxyz[j] +

                                 0.5 * fx[j] * pa_zz[j] * pb_xxyz[j] + pa_yyzz[j] * pb_xxyz[j]) * s_0_0[j];

                t_yyzz_xxzz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 0.125 * fx[j] * fx[j] * fx[j] * pa_zz[j] + 0.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_yyzz[j] * fx[j] * fx[j] + pa_yyz[j] * fx[j] * fx[j] * pb_z[j] +

                                 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_xx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_xx[j] +

                                 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_zz[j] + fx[j] * fx[j] * pa_z[j] * pb_xxz[j] + 0.5 * pa_yyzz[j] * pb_xx[j] * fx[j] +

                                 0.5 * pa_yyzz[j] * fx[j] * pb_zz[j] + 2.0 * pa_yyz[j] * fx[j] * pb_xxz[j] + 0.25 * fx[j] * fx[j] * pb_xxzz[j] +

                                 0.5 * pa_yy[j] * fx[j] * pb_xxzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxzz[j] + pa_yyzz[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (93) = (186,188)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xz, s_0_0, t_yyzz_xyyy, t_yyzz_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xyyy[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_yzz[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_xy[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_xy[j] +

                                 1.5 * pa_y[j] * fx[j] * fx[j] * pb_xyy[j] + 1.5 * pa_yyzz[j] * pb_xy[j] * fx[j] + 3.0 * pa_yzz[j] * fx[j] * pb_xyy[j] +

                                 0.25 * fx[j] * fx[j] * pb_xyyy[j] + 0.5 * pa_yy[j] * fx[j] * pb_xyyy[j] + 0.5 * fx[j] * pa_zz[j] * pb_xyyy[j] +

                                 pa_yyzz[j] * pb_xyyy[j]) * s_0_0[j];

                t_yyzz_xyyz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_yyz[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.0 * pa_yz[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_xz[j] +

                                 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_xz[j] + pa_y[j] * fx[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xyy[j] +

                                 0.5 * pa_yyzz[j] * pb_xz[j] * fx[j] + pa_yyz[j] * fx[j] * pb_xyy[j] + 2.0 * pa_yzz[j] * fx[j] * pb_xyz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xyyz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xyyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xyyz[j] +

                                 pa_yyzz[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (94) = (188,190)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xy, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, s_0_0, t_yyzz_xyzz, t_yyzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xyzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.5 * pa_yzz[j] * fx[j] * fx[j] * pb_x[j] + 2.0 * pa_yz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.5 * pa_y[j] * fx[j] * fx[j] * pb_xzz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_xy[j] + fx[j] * fx[j] * pa_z[j] * pb_xyz[j] +

                                 0.5 * pa_yyzz[j] * pb_xy[j] * fx[j] + 2.0 * pa_yyz[j] * fx[j] * pb_xyz[j] + pa_yzz[j] * fx[j] * pb_xzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xyzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xyzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xyzz[j] +

                                 pa_yyzz[j] * pb_xyzz[j]) * s_0_0[j];

                t_yyzz_xzzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yyz[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_xz[j] +

                                 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xzz[j] + 1.5 * pa_yyzz[j] * pb_xz[j] * fx[j] + 3.0 * pa_yyz[j] * fx[j] * pb_xzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_xzzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xzzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xzzz[j] +

                                 pa_yyzz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (95) = (190,192)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yz, pb_z, s_0_0, t_yyzz_yyyy, t_yyzz_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_yyyy[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 1.875 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 0.375 * pa_yy[j] * fx[j] * fx[j] * fx[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_yyzz[j] * fx[j] * fx[j] + 6.0 * pa_yzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 4.5 * fx[j] * fx[j] * pa_zz[j] * pb_yy[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_yy[j] +

                                 2.0 * pa_y[j] * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yyzz[j] * pb_yy[j] * fx[j] + 4.0 * pa_yzz[j] * fx[j] * pb_yyy[j] +

                                 0.25 * fx[j] * fx[j] * pb_yyyy[j] + 0.5 * pa_yy[j] * fx[j] * pb_yyyy[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyyy[j] +

                                 pa_yyzz[j] * pb_yyyy[j]) * s_0_0[j];

                t_yyzz_yyyz[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_yyz[j] * fx[j] * fx[j] * pb_y[j] +

                                 1.5 * pa_yzz[j] * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_yz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_yz[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_yz[j] +

                                 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_yyy[j] + 1.5 * pa_yyzz[j] * pb_yz[j] * fx[j] +

                                 pa_yyz[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yzz[j] * fx[j] * pb_yyz[j] + 0.25 * fx[j] * fx[j] * pb_yyyz[j] +

                                 0.5 * pa_yy[j] * fx[j] * pb_yyyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyyz[j] + pa_yyzz[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (96) = (192,194)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, s_0_0, t_yyzz_yyzz, t_yyzz_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_yyzz[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 1.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 0.25 * pa_yyzz[j] * fx[j] * fx[j] + pa_yyz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_yy[j] +

                                 pa_yzz[j] * fx[j] * fx[j] * pb_y[j] + 4.0 * pa_yz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_zz[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_zz[j] + pa_y[j] * fx[j] * fx[j] * pb_yzz[j] +

                                 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_yy[j] + fx[j] * fx[j] * pa_z[j] * pb_yyz[j] + 0.5 * pa_yyzz[j] * pb_yy[j] * fx[j] +

                                 0.5 * pa_yyzz[j] * fx[j] * pb_zz[j] + 2.0 * pa_yyz[j] * fx[j] * pb_yyz[j] + 2.0 * pa_yzz[j] * fx[j] * pb_yzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_yyzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_yyzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyzz[j] +

                                 pa_yyzz[j] * pb_yyzz[j]) * s_0_0[j];

                t_yyzz_yzzz[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_yyz[j] * fx[j] * fx[j] * pb_y[j] +

                                 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_yzz[j] * fx[j] * fx[j] * pb_z[j] +

                                 3.0 * pa_yz[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * pb_zzz[j] +

                                 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_yz[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_yzz[j] + 1.5 * pa_yyzz[j] * pb_yz[j] * fx[j] +

                                 3.0 * pa_yyz[j] * fx[j] * pb_yzz[j] + pa_yzz[j] * fx[j] * pb_zzz[j] + 0.25 * fx[j] * fx[j] * pb_yzzz[j] +

                                 0.5 * pa_yy[j] * fx[j] * pb_yzzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yzzz[j] + pa_yyzz[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (97) = (194,196)

            #pragma omp simd aligned(fx, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzzz, pa_z, pa_zz, pb_xx, pb_xxxx, pb_z, \
                                     pb_zz, pb_zzz, pb_zzzz, s_0_0, t_yyzz_zzzz, t_yzzz_xxxx: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_zzzz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 1.875 * pa_yy[j] * fx[j] * fx[j] * fx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] + 3.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_yyzz[j] * fx[j] * fx[j] + 6.0 * pa_yyz[j] * fx[j] * fx[j] * pb_z[j] +

                                 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_zz[j] +

                                 2.0 * fx[j] * fx[j] * pa_z[j] * pb_zzz[j] + 3.0 * pa_yyzz[j] * pb_zz[j] * fx[j] + 4.0 * pa_yyz[j] * fx[j] * pb_zzz[j] +

                                 0.25 * fx[j] * fx[j] * pb_zzzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_zzzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_zzzz[j] +

                                 pa_yyzz[j] * pb_zzzz[j]) * s_0_0[j];

                t_yzzz_xxxx[j] = (1.125 * pa_yz[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yzzz[j] * fx[j] * fx[j] +

                                 4.5 * pa_yz[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_yzzz[j] * pb_xx[j] * fx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxxx[j] +

                                 pa_yzzz[j] * pb_xxxx[j]) * s_0_0[j];
            }

            // Batch of Integrals (98) = (196,198)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zzz, pb_x, pb_xxx, pb_xxxy, pb_xxxz, \
                                     pb_xy, pb_xz, s_0_0, t_yzzz_xxxy, t_yzzz_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xxxy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_zzz[j] * pb_x[j] + 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xxx[j] + 1.5 * pa_yzzz[j] * pb_xy[j] * fx[j] + 0.5 * fx[j] * pa_zzz[j] * pb_xxx[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_xxxy[j] + pa_yzzz[j] * pb_xxxy[j]) * s_0_0[j];

                t_yzzz_xxxz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_yzz[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_yzzz[j] * pb_xz[j] * fx[j] + 1.5 * pa_yzz[j] * fx[j] * pb_xxx[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_xxxz[j] + pa_yzzz[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (99) = (198,200)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_xx, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_yzzz_xxyy, t_yzzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xxyy[j] = (0.375 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.25 * pa_yzzz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_zzz[j] * pb_y[j] +

                                 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xxy[j] + 0.5 * pa_yzzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yzzz[j] * fx[j] * pb_yy[j] +

                                 fx[j] * pa_zzz[j] * pb_xxy[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxyy[j] + pa_yzzz[j] * pb_xxyy[j]) * s_0_0[j];

                t_yzzz_xxyz[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_yzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.25 * fx[j] * fx[j] * pa_zzz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_xx[j] +

                                 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xxy[j] +

                                 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xxz[j] + 0.5 * pa_yzzz[j] * fx[j] * pb_yz[j] + 1.5 * pa_yzz[j] * fx[j] * pb_xxy[j] +

                                 0.5 * fx[j] * pa_zzz[j] * pb_xxz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxyz[j] + pa_yzzz[j] * pb_xxyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (100) = (200,202)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zzz, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_z, pb_zz, s_0_0, t_yzzz_xxzz, t_yzzz_xyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xxzz[j] = (1.125 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_yzzz[j] * fx[j] * fx[j] + 1.5 * pa_yzz[j] * fx[j] * fx[j] * pb_z[j] +

                                 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_y[j] * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_yzzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yzzz[j] * fx[j] * pb_zz[j] +

                                 3.0 * pa_yzz[j] * fx[j] * pb_xxz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxzz[j] + pa_yzzz[j] * pb_xxzz[j]) * s_0_0[j];

                t_yzzz_xyyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * pa_zzz[j] * pb_x[j] + 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xyy[j] + 1.5 * pa_yzzz[j] * pb_xy[j] * fx[j] + 1.5 * fx[j] * pa_zzz[j] * pb_xyy[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_xyyy[j] + pa_yzzz[j] * pb_xyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (101) = (202,204)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xy, pb_xyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, s_0_0, t_yzzz_xyyz, t_yzzz_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xyyz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * pa_yzz[j] * fx[j] * fx[j] * pb_x[j] +

                                 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_xy[j] + 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xyy[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_xyz[j] + 0.5 * pa_yzzz[j] * pb_xz[j] * fx[j] +

                                 1.5 * pa_yzz[j] * fx[j] * pb_xyy[j] + fx[j] * pa_zzz[j] * pb_xyz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xyyz[j] +

                                 pa_yzzz[j] * pb_xyyz[j]) * s_0_0[j];

                t_yzzz_xyzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_x[j] +

                                 0.75 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_xy[j] +

                                 0.25 * fx[j] * fx[j] * pa_zzz[j] * pb_x[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_xz[j] +

                                 1.5 * pa_y[j] * fx[j] * fx[j] * pb_xyz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xzz[j] + 0.5 * pa_yzzz[j] * pb_xy[j] * fx[j] +

                                 3.0 * pa_yzz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zzz[j] * pb_xzz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xyzz[j] +

                                 pa_yzzz[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (102) = (204,206)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zzz, pb_x, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyyy, s_0_0, t_yzzz_xzzz, t_yzzz_yyyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xzzz[j] = (1.875 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 2.25 * pa_yzz[j] * fx[j] * fx[j] * pb_x[j] + 6.75 * pa_yz[j] * fx[j] * fx[j] * pb_xz[j] +

                                 2.25 * pa_y[j] * fx[j] * fx[j] * pb_xzz[j] + 1.5 * pa_yzzz[j] * pb_xz[j] * fx[j] + 4.5 * pa_yzz[j] * fx[j] * pb_xzz[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_xzzz[j] + pa_yzzz[j] * pb_xzzz[j]) * s_0_0[j];

                t_yzzz_yyyy[j] = (1.125 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 4.5 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.75 * pa_yzzz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_zzz[j] * pb_y[j] +

                                 4.5 * pa_yz[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * fx[j] * fx[j] * pa_z[j] * pb_yyy[j] + 3.0 * pa_yzzz[j] * pb_yy[j] * fx[j] +

                                 2.0 * fx[j] * pa_zzz[j] * pb_yyy[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyyy[j] + pa_yzzz[j] * pb_yyyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (103) = (206,208)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, t_yzzz_yyyz, t_yzzz_yyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_yyyz[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] +

                                 1.125 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 2.25 * pa_yzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 0.75 * fx[j] * fx[j] * pa_zzz[j] * pb_z[j] + 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_yy[j] +

                                 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yyy[j] +

                                 2.25 * fx[j] * fx[j] * pa_z[j] * pb_yyz[j] + 1.5 * pa_yzzz[j] * pb_yz[j] * fx[j] + 1.5 * pa_yzz[j] * fx[j] * pb_yyy[j] +

                                 1.5 * fx[j] * pa_zzz[j] * pb_yyz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyyz[j] + pa_yzzz[j] * pb_yyyz[j]) * s_0_0[j];

                t_yzzz_yyzz[j] = (1.125 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 2.25 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_y[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] +

                                 1.5 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * pa_yzzz[j] * fx[j] * fx[j] + 1.5 * pa_yzz[j] * fx[j] * fx[j] * pb_z[j] +

                                 2.25 * pa_yz[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * fx[j] * pa_zzz[j] * pb_y[j] +

                                 3.0 * fx[j] * fx[j] * pa_zz[j] * pb_yz[j] + 0.75 * pa_yz[j] * fx[j] * fx[j] * pb_zz[j] +

                                 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yyz[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_yzz[j] + 0.5 * pa_yzzz[j] * pb_yy[j] * fx[j] +

                                 0.5 * pa_yzzz[j] * fx[j] * pb_zz[j] + 3.0 * pa_yzz[j] * fx[j] * pb_yyz[j] + fx[j] * pa_zzz[j] * pb_yzz[j] +

                                 1.5 * pa_yz[j] * fx[j] * pb_yyzz[j] + pa_yzzz[j] * pb_yyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (104) = (208,210)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_y, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_yzzz_yzzz, t_yzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_yzzz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] +

                                 1.875 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pa_zz[j] +

                                 3.375 * fx[j] * fx[j] * fx[j] * pa_z[j] * pb_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_zz[j] +

                                 2.25 * pa_yzz[j] * fx[j] * fx[j] * pb_y[j] + 6.75 * pa_yz[j] * fx[j] * fx[j] * pb_yz[j] +

                                 0.75 * fx[j] * fx[j] * pa_zzz[j] * pb_z[j] + 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_zz[j] +

                                 2.25 * pa_y[j] * fx[j] * fx[j] * pb_yzz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zzz[j] + 1.5 * pa_yzzz[j] * pb_yz[j] * fx[j] +

                                 4.5 * pa_yzz[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_zzz[j] * pb_zzz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yzzz[j] +

                                 pa_yzzz[j] * pb_yzzz[j]) * s_0_0[j];

                t_yzzz_zzzz[j] = (5.625 * pa_yz[j] * fx[j] * fx[j] * fx[j] +

                                 7.5 * pa_y[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_yzzz[j] * fx[j] * fx[j] + 9.0 * pa_yzz[j] * fx[j] * fx[j] * pb_z[j] +

                                 13.5 * pa_yz[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_yzzz[j] * pb_zz[j] * fx[j] +

                                 6.0 * pa_yzz[j] * fx[j] * pb_zzz[j] + 1.5 * pa_yz[j] * fx[j] * pb_zzzz[j] + pa_yzzz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (105) = (210,212)

            #pragma omp simd aligned(fx, pa_zz, pa_zzzz, pb_xx, pb_xxxx, pb_xxxy, pb_xy, s_0_0, t_zzzz_xxxx, \
                                     t_zzzz_xxxy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xxxx[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_zzzz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 9.0 * pa_zz[j] * fx[j] * fx[j] * pb_xx[j] +

                                 3.0 * pa_zzzz[j] * pb_xx[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xxxx[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxxx[j] +

                                 pa_zzzz[j] * pb_xxxx[j]) * s_0_0[j];

                t_zzzz_xxxy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_zzzz[j] * pb_xy[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xxxy[j] +

                                 3.0 * pa_zz[j] * fx[j] * pb_xxxy[j] + pa_zzzz[j] * pb_xxxy[j]) * s_0_0[j];
            }

            // Batch of Integrals (106) = (212,214)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xx, pb_xxx, pb_xxxz, pb_xxyy, pb_xz, \
                                     pb_yy, s_0_0, t_zzzz_xxxz, t_zzzz_xxyy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xxxz[j] = (4.5 * pa_z[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 3.0 * pa_zzz[j] * fx[j] * fx[j] * pb_x[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_xz[j] + 3.0 * pa_z[j] * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_zzzz[j] * pb_xz[j] * fx[j] +

                                 2.0 * pa_zzz[j] * fx[j] * pb_xxx[j] + 0.75 * fx[j] * fx[j] * pb_xxxz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxxz[j] +

                                 pa_zzzz[j] * pb_xxxz[j]) * s_0_0[j];

                t_zzzz_xxyy[j] = (0.1875 * fx[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] * fx[j] +

                                 0.25 * pa_zzzz[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_yy[j] +

                                 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_zzzz[j] * pb_xx[j] * fx[j] +

                                 0.5 * pa_zzzz[j] * fx[j] * pb_yy[j] + 0.75 * fx[j] * fx[j] * pb_xxyy[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxyy[j] +

                                 pa_zzzz[j] * pb_xxyy[j]) * s_0_0[j];
            }

            // Batch of Integrals (107) = (214,216)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_xx, pb_xxy, pb_xxyz, pb_xxz, pb_xxzz, pb_y, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_zzzz_xxyz, t_zzzz_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xxyz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 pa_zzz[j] * fx[j] * fx[j] * pb_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_yz[j] +

                                 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_yz[j] + 3.0 * pa_z[j] * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_zzzz[j] * fx[j] * pb_yz[j] +

                                 2.0 * pa_zzz[j] * fx[j] * pb_xxy[j] + 0.75 * fx[j] * fx[j] * pb_xxyz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxyz[j] +

                                 pa_zzzz[j] * pb_xxyz[j]) * s_0_0[j];

                t_zzzz_xxzz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * fx[j] +

                                 3.0 * pa_z[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.875 * fx[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_zzzz[j] * fx[j] * fx[j] +

                                 2.0 * pa_zzz[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_xx[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_zz[j] +

                                 6.0 * pa_z[j] * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_zzzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_zzzz[j] * fx[j] * pb_zz[j] +

                                 4.0 * pa_zzz[j] * fx[j] * pb_xxz[j] + 0.75 * fx[j] * fx[j] * pb_xxzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxzz[j] +

                                 pa_zzzz[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (108) = (216,218)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xz, \
                                     s_0_0, t_zzzz_xyyy, t_zzzz_xyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xyyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_zzzz[j] * pb_xy[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xyyy[j] +

                                 3.0 * pa_zz[j] * fx[j] * pb_xyyy[j] + pa_zzzz[j] * pb_xyyy[j]) * s_0_0[j];

                t_zzzz_xyyz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 pa_zzz[j] * fx[j] * fx[j] * pb_x[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_xz[j] +

                                 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_xz[j] + 3.0 * pa_z[j] * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_zzzz[j] * pb_xz[j] * fx[j] +

                                 2.0 * pa_zzz[j] * fx[j] * pb_xyy[j] + 0.75 * fx[j] * fx[j] * pb_xyyz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xyyz[j] +

                                 pa_zzzz[j] * pb_xyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (109) = (218,220)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xy, pb_xyz, pb_xyzz, pb_xz, pb_xzz, \
                                     pb_xzzz, s_0_0, t_zzzz_xyzz, t_zzzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xyzz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_xy[j] +

                                 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_xy[j] + 6.0 * pa_z[j] * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_zzzz[j] * pb_xy[j] * fx[j] +

                                 4.0 * pa_zzz[j] * fx[j] * pb_xyz[j] + 0.75 * fx[j] * fx[j] * pb_xyzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xyzz[j] +

                                 pa_zzzz[j] * pb_xyzz[j]) * s_0_0[j];

                t_zzzz_xzzz[j] = (7.5 * pa_z[j] * fx[j] * fx[j] * fx[j] * pb_x[j] +

                                 5.625 * fx[j] * fx[j] * fx[j] * pb_xz[j] + 3.0 * pa_zzz[j] * fx[j] * fx[j] * pb_x[j] +

                                 13.5 * pa_zz[j] * fx[j] * fx[j] * pb_xz[j] + 9.0 * pa_z[j] * fx[j] * fx[j] * pb_xzz[j] + 1.5 * pa_zzzz[j] * pb_xz[j] * fx[j] +

                                 6.0 * pa_zzz[j] * fx[j] * pb_xzz[j] + 0.75 * fx[j] * fx[j] * pb_xzzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xzzz[j] +

                                 pa_zzzz[j] * pb_xzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (110) = (220,222)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yz, \
                                     s_0_0, t_zzzz_yyyy, t_zzzz_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_yyyy[j] = (0.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * fx[j] +

                                 0.75 * pa_zzzz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 9.0 * pa_zz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 3.0 * pa_zzzz[j] * pb_yy[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yyyy[j] + 3.0 * pa_zz[j] * fx[j] * pb_yyyy[j] +

                                 pa_zzzz[j] * pb_yyyy[j]) * s_0_0[j];

                t_zzzz_yyyz[j] = (4.5 * pa_z[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 3.0 * pa_zzz[j] * fx[j] * fx[j] * pb_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_yz[j] +

                                 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_yz[j] + 3.0 * pa_z[j] * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_zzzz[j] * pb_yz[j] * fx[j] +

                                 2.0 * pa_zzz[j] * fx[j] * pb_yyy[j] + 0.75 * fx[j] * fx[j] * pb_yyyz[j] + 3.0 * pa_zz[j] * fx[j] * pb_yyyz[j] +

                                 pa_zzzz[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (111) = (222,224)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_y, pb_yy, pb_yyz, pb_yyzz, pb_yz, pb_yzz, \
                                     pb_yzzz, pb_z, pb_zz, s_0_0, t_zzzz_yyzz, t_zzzz_yzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_yyzz[j] = (0.9375 * fx[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * fx[j] +

                                 3.0 * pa_z[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.875 * fx[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * pa_zzzz[j] * fx[j] * fx[j] +

                                 2.0 * pa_zzz[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_yy[j] +

                                 0.375 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_zz[j] +

                                 6.0 * pa_z[j] * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_zzzz[j] * pb_yy[j] * fx[j] + 0.5 * pa_zzzz[j] * fx[j] * pb_zz[j] +

                                 4.0 * pa_zzz[j] * fx[j] * pb_yyz[j] + 0.75 * fx[j] * fx[j] * pb_yyzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_yyzz[j] +

                                 pa_zzzz[j] * pb_yyzz[j]) * s_0_0[j];

                t_zzzz_yzzz[j] = (7.5 * pa_z[j] * fx[j] * fx[j] * fx[j] * pb_y[j] +

                                 5.625 * fx[j] * fx[j] * fx[j] * pb_yz[j] + 3.0 * pa_zzz[j] * fx[j] * fx[j] * pb_y[j] +

                                 13.5 * pa_zz[j] * fx[j] * fx[j] * pb_yz[j] + 9.0 * pa_z[j] * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_zzzz[j] * pb_yz[j] * fx[j] +

                                 6.0 * pa_zzz[j] * fx[j] * pb_yzz[j] + 0.75 * fx[j] * fx[j] * pb_yzzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_yzzz[j] +

                                 pa_zzzz[j] * pb_yzzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (112) = (224,225)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, \
                                     t_zzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_zzzz[j] = (6.5625 * fx[j] * fx[j] * fx[j] * fx[j] + 11.25 * pa_zz[j] * fx[j] * fx[j] * fx[j] +

                                 30.0 * pa_z[j] * fx[j] * fx[j] * fx[j] * pb_z[j] + 11.25 * fx[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_zzzz[j] * fx[j] * fx[j] +

                                 12.0 * pa_zzz[j] * fx[j] * fx[j] * pb_z[j] + 27.0 * pa_zz[j] * fx[j] * fx[j] * pb_zz[j] +

                                 12.0 * pa_z[j] * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_zzzz[j] * pb_zz[j] * fx[j] + 8.0 * pa_zzz[j] * fx[j] * pb_zzz[j] +

                                 0.75 * fx[j] * fx[j] * pb_zzzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_zzzz[j] + pa_zzzz[j] * pb_zzzz[j]) * s_0_0[j];
            }

            idx++;
        }
    }

    
    
} // ovlrecfunc namespace

