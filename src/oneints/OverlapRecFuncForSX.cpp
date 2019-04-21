//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapRecFuncForSX.hpp"

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
                     const CMemBlock2D<double>& auxBuffer,
                     const CMemBlock2D<double>& pbDistances,
                     const CGtoBlock&           braGtoBlock,
                     const CGtoBlock&           ketGtoBlock,
                     const int32_t              iContrGto)
    {
        // Batch of Integrals (0,3)

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

            // Batch of Integrals (0,3)

            #pragma omp simd aligned(pb_x, pb_y, pb_z, s_0_0, t_0_x, t_0_y, t_0_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                t_0_x[j] = fl_s_0_0 * pb_x[j];

                t_0_y[j] = fl_s_0_0 * pb_y[j];

                t_0_z[j] = fl_s_0_0 * pb_z[j];
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
        // Batch of Integrals (0,3)

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

            // Batch of Integrals (0,3)

            #pragma omp simd aligned(pa_x, pa_y, pa_z, s_0_0, t_x_0, t_y_0, t_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                t_x_0[j] = fl_s_0_0 * pa_x[j];

                t_y_0[j] = fl_s_0_0 * pa_y[j];

                t_z_0[j] = fl_s_0_0 * pa_z[j];
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
        // Batch of Integrals (0,6)

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

            // Batch of Integrals (0,6)

            #pragma omp simd aligned(fx, pb_xx, pb_xy, pb_xz, pb_yy, pb_yz, pb_zz, s_0_0, t_0_xx, t_0_xy, t_0_xz, \
                                     t_0_yy, t_0_yz, t_0_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                t_0_xx[j] = fl_s_0_0 * (0.5 * fl1_fx + pb_xx[j]);

                t_0_xy[j] = fl_s_0_0 * pb_xy[j];

                t_0_xz[j] = fl_s_0_0 * pb_xz[j];

                t_0_yy[j] = fl_s_0_0 * (0.5 * fl1_fx + pb_yy[j]);

                t_0_yz[j] = fl_s_0_0 * pb_yz[j];

                t_0_zz[j] = fl_s_0_0 * (0.5 * fl1_fx + pb_zz[j]);
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
        // Batch of Integrals (0,6)

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

            // Batch of Integrals (0,6)

            #pragma omp simd aligned(fx, pa_xx, pa_xy, pa_xz, pa_yy, pa_yz, pa_zz, s_0_0, t_xx_0, t_xy_0, t_xz_0, \
                                     t_yy_0, t_yz_0, t_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                t_xx_0[j] = fl_s_0_0 * (0.5 * fl1_fx + pa_xx[j]);

                t_xy_0[j] = fl_s_0_0 * pa_xy[j];

                t_xz_0[j] = fl_s_0_0 * pa_xz[j];

                t_yy_0[j] = fl_s_0_0 * (0.5 * fl1_fx + pa_yy[j]);

                t_yz_0[j] = fl_s_0_0 * pa_yz[j];

                t_zz_0[j] = fl_s_0_0 * (0.5 * fl1_fx + pa_zz[j]);
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
        // Batch of Integrals (0,10)

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

            // set up pointers to 2-th order tensor of distance R(PB)

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

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pb_x, pb_xxx, pb_xxy, pb_xxz, pb_xyy, pb_xyz, pb_xzz, pb_y, pb_yyy, pb_yyz, \
                                     pb_yzz, pb_z, pb_zzz, s_0_0, t_0_xxx, t_0_xxy, t_0_xxz, t_0_xyy, t_0_xyz, t_0_xzz, \
                                     t_0_yyy, t_0_yyz, t_0_yzz, t_0_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                t_0_xxx[j] = fl_s_0_0 * (1.5 * pb_x[j] * fl1_fx + pb_xxx[j]);

                t_0_xxy[j] = fl_s_0_0 * (0.5 * fl1_fx * pb_y[j] + pb_xxy[j]);

                t_0_xxz[j] = fl_s_0_0 * (0.5 * fl1_fx * pb_z[j] + pb_xxz[j]);

                t_0_xyy[j] = fl_s_0_0 * (0.5 * pb_x[j] * fl1_fx + pb_xyy[j]);

                t_0_xyz[j] = fl_s_0_0 * pb_xyz[j];

                t_0_xzz[j] = fl_s_0_0 * (0.5 * pb_x[j] * fl1_fx + pb_xzz[j]);

                t_0_yyy[j] = fl_s_0_0 * (1.5 * pb_y[j] * fl1_fx + pb_yyy[j]);

                t_0_yyz[j] = fl_s_0_0 * (0.5 * fl1_fx * pb_z[j] + pb_yyz[j]);

                t_0_yzz[j] = fl_s_0_0 * (0.5 * pb_y[j] * fl1_fx + pb_yzz[j]);

                t_0_zzz[j] = fl_s_0_0 * (1.5 * pb_z[j] * fl1_fx + pb_zzz[j]);
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
        // Batch of Integrals (0,10)

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

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxy, pa_xxz, pa_xyy, pa_xyz, pa_xzz, pa_y, pa_yyy, pa_yyz, \
                                     pa_yzz, pa_z, pa_zzz, s_0_0, t_xxx_0, t_xxy_0, t_xxz_0, t_xyy_0, t_xyz_0, t_xzz_0, \
                                     t_yyy_0, t_yyz_0, t_yzz_0, t_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                t_xxx_0[j] = fl_s_0_0 * (1.5 * pa_x[j] * fl1_fx + pa_xxx[j]);

                t_xxy_0[j] = fl_s_0_0 * (0.5 * fl1_fx * pa_y[j] + pa_xxy[j]);

                t_xxz_0[j] = fl_s_0_0 * (0.5 * fl1_fx * pa_z[j] + pa_xxz[j]);

                t_xyy_0[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl1_fx + pa_xyy[j]);

                t_xyz_0[j] = fl_s_0_0 * pa_xyz[j];

                t_xzz_0[j] = fl_s_0_0 * (0.5 * pa_x[j] * fl1_fx + pa_xzz[j]);

                t_yyy_0[j] = fl_s_0_0 * (1.5 * pa_y[j] * fl1_fx + pa_yyy[j]);

                t_yyz_0[j] = fl_s_0_0 * (0.5 * fl1_fx * pa_z[j] + pa_yyz[j]);

                t_yzz_0[j] = fl_s_0_0 * (0.5 * pa_y[j] * fl1_fx + pa_yzz[j]);

                t_zzz_0[j] = fl_s_0_0 * (1.5 * pa_z[j] * fl1_fx + pa_zzz[j]);
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
        // Batch of Integrals (0,15)

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

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xx = pbDistances.data(34 * idx + 3);

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

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

            // Batch of Integrals (0,15)

            #pragma omp simd aligned(fx, pb_xx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxyy, pb_xxyz, pb_xxzz, pb_xy, \
                                     pb_xyyy, pb_xyyz, pb_xyzz, pb_xz, pb_xzzz, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, \
                                     pb_yzzz, pb_zz, pb_zzzz, s_0_0, t_0_xxxx, t_0_xxxy, t_0_xxxz, t_0_xxyy, t_0_xxyz, \
                                     t_0_xxzz, t_0_xyyy, t_0_xyyz, t_0_xyzz, t_0_xzzz, t_0_yyyy, t_0_yyyz, t_0_yyzz, \
                                     t_0_yzzz, t_0_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_0_xxxx[j] = fl_s_0_0 * (0.75 * fl2_fx + 3.0 * pb_xx[j] * fl1_fx + pb_xxxx[j]);

                t_0_xxxy[j] = fl_s_0_0 * (1.5 * pb_xy[j] * fl1_fx + pb_xxxy[j]);

                t_0_xxxz[j] = fl_s_0_0 * (1.5 * pb_xz[j] * fl1_fx + pb_xxxz[j]);

                t_0_xxyy[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pb_xx[j] * fl1_fx + 0.5 * fl1_fx * pb_yy[j] + pb_xxyy[j]);

                t_0_xxyz[j] = fl_s_0_0 * (0.5 * fl1_fx * pb_yz[j] + pb_xxyz[j]);

                t_0_xxzz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pb_xx[j] * fl1_fx + 0.5 * fl1_fx * pb_zz[j] + pb_xxzz[j]);

                t_0_xyyy[j] = fl_s_0_0 * (1.5 * pb_xy[j] * fl1_fx + pb_xyyy[j]);

                t_0_xyyz[j] = fl_s_0_0 * (0.5 * pb_xz[j] * fl1_fx + pb_xyyz[j]);

                t_0_xyzz[j] = fl_s_0_0 * (0.5 * pb_xy[j] * fl1_fx + pb_xyzz[j]);

                t_0_xzzz[j] = fl_s_0_0 * (1.5 * pb_xz[j] * fl1_fx + pb_xzzz[j]);

                t_0_yyyy[j] = fl_s_0_0 * (0.75 * fl2_fx + 3.0 * pb_yy[j] * fl1_fx + pb_yyyy[j]);

                t_0_yyyz[j] = fl_s_0_0 * (1.5 * pb_yz[j] * fl1_fx + pb_yyyz[j]);

                t_0_yyzz[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pb_yy[j] * fl1_fx + 0.5 * fl1_fx * pb_zz[j] + pb_yyzz[j]);

                t_0_yzzz[j] = fl_s_0_0 * (1.5 * pb_yz[j] * fl1_fx + pb_yzzz[j]);

                t_0_zzzz[j] = fl_s_0_0 * (0.75 * fl2_fx + 3.0 * pb_zz[j] * fl1_fx + pb_zzzz[j]);
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
        // Batch of Integrals (0,15)

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

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

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

            // Batch of Integrals (0,15)

            #pragma omp simd aligned(fx, pa_xx, pa_xxxx, pa_xxxy, pa_xxxz, pa_xxyy, pa_xxyz, pa_xxzz, pa_xy, \
                                     pa_xyyy, pa_xyyz, pa_xyzz, pa_xz, pa_xzzz, pa_yy, pa_yyyy, pa_yyyz, pa_yyzz, pa_yz, \
                                     pa_yzzz, pa_zz, pa_zzzz, s_0_0, t_xxxx_0, t_xxxy_0, t_xxxz_0, t_xxyy_0, t_xxyz_0, \
                                     t_xxzz_0, t_xyyy_0, t_xyyz_0, t_xyzz_0, t_xzzz_0, t_yyyy_0, t_yyyz_0, t_yyzz_0, \
                                     t_yzzz_0, t_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                t_xxxx_0[j] = fl_s_0_0 * (0.75 * fl2_fx + 3.0 * pa_xx[j] * fl1_fx + pa_xxxx[j]);

                t_xxxy_0[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl1_fx + pa_xxxy[j]);

                t_xxxz_0[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl1_fx + pa_xxxz[j]);

                t_xxyy_0[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * fl1_fx * pa_yy[j] + pa_xxyy[j]);

                t_xxyz_0[j] = fl_s_0_0 * (0.5 * fl1_fx * pa_yz[j] + pa_xxyz[j]);

                t_xxzz_0[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_xx[j] * fl1_fx + 0.5 * fl1_fx * pa_zz[j] + pa_xxzz[j]);

                t_xyyy_0[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl1_fx + pa_xyyy[j]);

                t_xyyz_0[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl1_fx + pa_xyyz[j]);

                t_xyzz_0[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl1_fx + pa_xyzz[j]);

                t_xzzz_0[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl1_fx + pa_xzzz[j]);

                t_yyyy_0[j] = fl_s_0_0 * (0.75 * fl2_fx + 3.0 * pa_yy[j] * fl1_fx + pa_yyyy[j]);

                t_yyyz_0[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl1_fx + pa_yyyz[j]);

                t_yyzz_0[j] = fl_s_0_0 * (0.25 * fl2_fx + 0.5 * pa_yy[j] * fl1_fx + 0.5 * fl1_fx * pa_zz[j] + pa_yyzz[j]);

                t_yzzz_0[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl1_fx + pa_yzzz[j]);

                t_zzzz_0[j] = fl_s_0_0 * (0.75 * fl2_fx + 3.0 * pa_zz[j] * fl1_fx + pa_zzzz[j]);
            }

            idx++;
        }
    }


} // ovlrecfunc namespace

