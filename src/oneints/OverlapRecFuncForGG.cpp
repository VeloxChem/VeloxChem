//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapRecFuncForGG.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace

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
        ovlrecfunc::compOverlapForGG_0_10(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                          braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_10_20(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_20_30(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_30_40(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_40_50(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_50_60(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_60_70(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_70_80(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_80_90(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                           braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_90_100(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                            braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_100_110(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_110_120(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_120_130(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_130_140(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_140_150(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_150_160(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_160_170(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_170_180(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_180_189(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_189_198(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_198_207(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_207_216(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 

        ovlrecfunc::compOverlapForGG_216_225(primBuffer, auxBuffer, osFactors, paDistances, pbDistances, 
                                             braGtoBlock, ketGtoBlock, iContrGto); 
    }

    void
    compOverlapForGG_0_10(      CMemBlock2D<double>& primBuffer,
                          const CMemBlock2D<double>& auxBuffer,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
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

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxx = paDistances.data(34 * idx + 19);

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

            // Batch of Integrals (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, s_0_0, t_xxxx_xxxx, t_xxxx_xxxy, t_xxxx_xxxz, t_xxxx_xxyy, t_xxxx_xxyz, \
                                     t_xxxx_xxzz, t_xxxx_xyyy, t_xxxx_xyyz, t_xxxx_xyzz, t_xxxx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxx_xxxx[j] = fl_s_0_0 * (6.5625 * fl4_fx + 11.25 * pa_xx[j] * fl3_fx + 30.0 * pa_x[j] * fl3_fx * pb_x[j] + 11.25 * fl3_fx * pb_xx[j] + 0.75 * pa_xxxx[j] * fl2_fx + 12.0 * pa_xxx[j] * fl2_fx * pb_x[j] + 27.0 * pa_xx[j] * fl2_fx * pb_xx[j] + 12.0 * pa_x[j] * fl2_fx * pb_xxx[j] + 3.0 * pa_xxxx[j] * pb_xx[j] * fl1_fx + 8.0 * pa_xxx[j] * fl1_fx * pb_xxx[j] + 0.75 * fl2_fx * pb_xxxx[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xxxx[j] + pa_xxxx[j] * pb_xxxx[j]);

                t_xxxx_xxxy[j] = fl_s_0_0 * (7.5 * pa_x[j] * fl3_fx * pb_y[j] + 5.625 * fl3_fx * pb_xy[j] + 3.0 * pa_xxx[j] * fl2_fx * pb_y[j] + 13.5 * pa_xx[j] * fl2_fx * pb_xy[j] + 9.0 * pa_x[j] * fl2_fx * pb_xxy[j] + 1.5 * pa_xxxx[j] * pb_xy[j] * fl1_fx + 6.0 * pa_xxx[j] * fl1_fx * pb_xxy[j] + 0.75 * fl2_fx * pb_xxxy[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xxxy[j] + pa_xxxx[j] * pb_xxxy[j]);

                t_xxxx_xxxz[j] = fl_s_0_0 * (7.5 * pa_x[j] * fl3_fx * pb_z[j] + 5.625 * fl3_fx * pb_xz[j] + 3.0 * pa_xxx[j] * fl2_fx * pb_z[j] + 13.5 * pa_xx[j] * fl2_fx * pb_xz[j] + 9.0 * pa_x[j] * fl2_fx * pb_xxz[j] + 1.5 * pa_xxxx[j] * pb_xz[j] * fl1_fx + 6.0 * pa_xxx[j] * fl1_fx * pb_xxz[j] + 0.75 * fl2_fx * pb_xxxz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xxxz[j] + pa_xxxx[j] * pb_xxxz[j]);

                t_xxxx_xxyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_xx[j] * fl3_fx + 3.0 * pa_x[j] * fl3_fx * pb_x[j] + 1.875 * fl3_fx * pb_yy[j] + 0.25 * pa_xxxx[j] * fl2_fx + 2.0 * pa_xxx[j] * fl2_fx * pb_x[j] + 4.5 * pa_xx[j] * fl2_fx * pb_yy[j] + 0.375 * fl3_fx * pb_xx[j] + 1.5 * pa_xx[j] * fl2_fx * pb_xx[j] + 6.0 * pa_x[j] * fl2_fx * pb_xyy[j] + 0.5 * pa_xxxx[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxxx[j] * fl1_fx * pb_yy[j] + 4.0 * pa_xxx[j] * fl1_fx * pb_xyy[j] + 0.75 * fl2_fx * pb_xxyy[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xxyy[j] + pa_xxxx[j] * pb_xxyy[j]);

                t_xxxx_xxyz[j] = fl_s_0_0 * (1.875 * fl3_fx * pb_yz[j] + 4.5 * pa_xx[j] * fl2_fx * pb_yz[j] + 6.0 * pa_x[j] * fl2_fx * pb_xyz[j] + 0.5 * pa_xxxx[j] * fl1_fx * pb_yz[j] + 4.0 * pa_xxx[j] * fl1_fx * pb_xyz[j] + 0.75 * fl2_fx * pb_xxyz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xxyz[j] + pa_xxxx[j] * pb_xxyz[j]);

                t_xxxx_xxzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_xx[j] * fl3_fx + 3.0 * pa_x[j] * fl3_fx * pb_x[j] + 1.875 * fl3_fx * pb_zz[j] + 0.25 * pa_xxxx[j] * fl2_fx + 2.0 * pa_xxx[j] * fl2_fx * pb_x[j] + 4.5 * pa_xx[j] * fl2_fx * pb_zz[j] + 0.375 * fl3_fx * pb_xx[j] + 1.5 * pa_xx[j] * fl2_fx * pb_xx[j] + 6.0 * pa_x[j] * fl2_fx * pb_xzz[j] + 0.5 * pa_xxxx[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxxx[j] * fl1_fx * pb_zz[j] + 4.0 * pa_xxx[j] * fl1_fx * pb_xzz[j] + 0.75 * fl2_fx * pb_xxzz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xxzz[j] + pa_xxxx[j] * pb_xxzz[j]);

                t_xxxx_xyyy[j] = fl_s_0_0 * (4.5 * pa_x[j] * fl3_fx * pb_y[j] + 3.0 * pa_xxx[j] * fl2_fx * pb_y[j] + 1.125 * fl3_fx * pb_xy[j] + 4.5 * pa_xx[j] * fl2_fx * pb_xy[j] + 3.0 * pa_x[j] * fl2_fx * pb_yyy[j] + 1.5 * pa_xxxx[j] * pb_xy[j] * fl1_fx + 2.0 * pa_xxx[j] * fl1_fx * pb_yyy[j] + 0.75 * fl2_fx * pb_xyyy[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xyyy[j] + pa_xxxx[j] * pb_xyyy[j]);

                t_xxxx_xyyz[j] = fl_s_0_0 * (1.5 * pa_x[j] * fl3_fx * pb_z[j] + pa_xxx[j] * fl2_fx * pb_z[j] + 0.375 * fl3_fx * pb_xz[j] + 1.5 * pa_xx[j] * fl2_fx * pb_xz[j] + 3.0 * pa_x[j] * fl2_fx * pb_yyz[j] + 0.5 * pa_xxxx[j] * pb_xz[j] * fl1_fx + 2.0 * pa_xxx[j] * fl1_fx * pb_yyz[j] + 0.75 * fl2_fx * pb_xyyz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xyyz[j] + pa_xxxx[j] * pb_xyyz[j]);

                t_xxxx_xyzz[j] = fl_s_0_0 * (1.5 * pa_x[j] * fl3_fx * pb_y[j] + pa_xxx[j] * fl2_fx * pb_y[j] + 0.375 * fl3_fx * pb_xy[j] + 1.5 * pa_xx[j] * fl2_fx * pb_xy[j] + 3.0 * pa_x[j] * fl2_fx * pb_yzz[j] + 0.5 * pa_xxxx[j] * pb_xy[j] * fl1_fx + 2.0 * pa_xxx[j] * fl1_fx * pb_yzz[j] + 0.75 * fl2_fx * pb_xyzz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xyzz[j] + pa_xxxx[j] * pb_xyzz[j]);

                t_xxxx_xzzz[j] = fl_s_0_0 * (4.5 * pa_x[j] * fl3_fx * pb_z[j] + 3.0 * pa_xxx[j] * fl2_fx * pb_z[j] + 1.125 * fl3_fx * pb_xz[j] + 4.5 * pa_xx[j] * fl2_fx * pb_xz[j] + 3.0 * pa_x[j] * fl2_fx * pb_zzz[j] + 1.5 * pa_xxxx[j] * pb_xz[j] * fl1_fx + 2.0 * pa_xxx[j] * fl1_fx * pb_zzz[j] + 0.75 * fl2_fx * pb_xzzz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_xzzz[j] + pa_xxxx[j] * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_10_20(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (10,20)

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

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxx = paDistances.data(34 * idx + 19);

            auto pa_xxxy = paDistances.data(34 * idx + 20);

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

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (10,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, \
                                     pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_y, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_z, pb_zz, \
                                     pb_zzzz, s_0_0, t_xxxx_yyyy, t_xxxx_yyyz, t_xxxx_yyzz, t_xxxx_yzzz, t_xxxx_zzzz, \
                                     t_xxxy_xxxx, t_xxxy_xxxy, t_xxxy_xxxz, t_xxxy_xxyy, t_xxxy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxx_yyyy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_xx[j] * fl3_fx + 0.75 * pa_xxxx[j] * fl2_fx + 2.25 * fl3_fx * pb_yy[j] + 9.0 * pa_xx[j] * fl2_fx * pb_yy[j] + 3.0 * pa_xxxx[j] * pb_yy[j] * fl1_fx + 0.75 * fl2_fx * pb_yyyy[j] + 3.0 * pa_xx[j] * fl1_fx * pb_yyyy[j] + pa_xxxx[j] * pb_yyyy[j]);

                t_xxxx_yyyz[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_yz[j] + 4.5 * pa_xx[j] * fl2_fx * pb_yz[j] + 1.5 * pa_xxxx[j] * pb_yz[j] * fl1_fx + 0.75 * fl2_fx * pb_yyyz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_yyyz[j] + pa_xxxx[j] * pb_yyyz[j]);

                t_xxxx_yyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.75 * pa_xx[j] * fl3_fx + 0.25 * pa_xxxx[j] * fl2_fx + 0.375 * fl3_fx * pb_yy[j] + 0.375 * fl3_fx * pb_zz[j] + 1.5 * pa_xx[j] * fl2_fx * pb_yy[j] + 1.5 * pa_xx[j] * fl2_fx * pb_zz[j] + 0.5 * pa_xxxx[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xxxx[j] * fl1_fx * pb_zz[j] + 0.75 * fl2_fx * pb_yyzz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_yyzz[j] + pa_xxxx[j] * pb_yyzz[j]);

                t_xxxx_yzzz[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_yz[j] + 4.5 * pa_xx[j] * fl2_fx * pb_yz[j] + 1.5 * pa_xxxx[j] * pb_yz[j] * fl1_fx + 0.75 * fl2_fx * pb_yzzz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_yzzz[j] + pa_xxxx[j] * pb_yzzz[j]);

                t_xxxx_zzzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_xx[j] * fl3_fx + 0.75 * pa_xxxx[j] * fl2_fx + 2.25 * fl3_fx * pb_zz[j] + 9.0 * pa_xx[j] * fl2_fx * pb_zz[j] + 3.0 * pa_xxxx[j] * pb_zz[j] * fl1_fx + 0.75 * fl2_fx * pb_zzzz[j] + 3.0 * pa_xx[j] * fl1_fx * pb_zzzz[j] + pa_xxxx[j] * pb_zzzz[j]);

                t_xxxy_xxxx[j] = fl_s_0_0 * (5.625 * pa_xy[j] * fl3_fx + 7.5 * fl3_fx * pa_y[j] * pb_x[j] + 0.75 * pa_xxxy[j] * fl2_fx + 9.0 * pa_xxy[j] * fl2_fx * pb_x[j] + 13.5 * pa_xy[j] * fl2_fx * pb_xx[j] + 3.0 * fl2_fx * pa_y[j] * pb_xxx[j] + 3.0 * pa_xxxy[j] * pb_xx[j] * fl1_fx + 6.0 * pa_xxy[j] * fl1_fx * pb_xxx[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxxx[j] + pa_xxxy[j] * pb_xxxx[j]);

                t_xxxy_xxxy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.125 * pa_xx[j] * fl3_fx + 3.375 * pa_x[j] * fl3_fx * pb_x[j] + 1.875 * fl3_fx * pa_y[j] * pb_y[j] + 1.125 * fl3_fx * pb_xx[j] + 0.75 * pa_xxx[j] * fl2_fx * pb_x[j] + 2.25 * pa_xxy[j] * fl2_fx * pb_y[j] + 2.25 * pa_xx[j] * fl2_fx * pb_xx[j] + 6.75 * pa_xy[j] * fl2_fx * pb_xy[j] + 0.75 * pa_x[j] * fl2_fx * pb_xxx[j] + 2.25 * fl2_fx * pa_y[j] * pb_xxy[j] + 1.5 * pa_xxxy[j] * pb_xy[j] * fl1_fx + 0.5 * pa_xxx[j] * fl1_fx * pb_xxx[j] + 4.5 * pa_xxy[j] * fl1_fx * pb_xxy[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxxy[j] + pa_xxxy[j] * pb_xxxy[j]);

                t_xxxy_xxxz[j] = fl_s_0_0 * (1.875 * fl3_fx * pa_y[j] * pb_z[j] + 2.25 * pa_xxy[j] * fl2_fx * pb_z[j] + 6.75 * pa_xy[j] * fl2_fx * pb_xz[j] + 2.25 * fl2_fx * pa_y[j] * pb_xxz[j] + 1.5 * pa_xxxy[j] * pb_xz[j] * fl1_fx + 4.5 * pa_xxy[j] * fl1_fx * pb_xxz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxxz[j] + pa_xxxy[j] * pb_xxxz[j]);

                t_xxxy_xxyy[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 2.25 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * fl3_fx * pa_y[j] * pb_x[j] + 1.5 * fl3_fx * pb_xy[j] + 0.25 * pa_xxxy[j] * fl2_fx + 0.5 * pa_xxx[j] * fl2_fx * pb_y[j] + 1.5 * pa_xxy[j] * fl2_fx * pb_x[j] + 3.0 * pa_xx[j] * fl2_fx * pb_xy[j] + 2.25 * pa_xy[j] * fl2_fx * pb_yy[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xx[j] + 1.5 * pa_x[j] * fl2_fx * pb_xxy[j] + 1.5 * fl2_fx * pa_y[j] * pb_xyy[j] + 0.5 * pa_xxxy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxxy[j] * fl1_fx * pb_yy[j] + pa_xxx[j] * fl1_fx * pb_xxy[j] + 3.0 * pa_xxy[j] * fl1_fx * pb_xyy[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxyy[j] + pa_xxxy[j] * pb_xxyy[j]);

                t_xxxy_xxyz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx * pb_z[j] + 0.75 * fl3_fx * pb_xz[j] + 0.25 * pa_xxx[j] * fl2_fx * pb_z[j] + 1.5 * pa_xx[j] * fl2_fx * pb_xz[j] + 2.25 * pa_xy[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xxz[j] + 1.5 * fl2_fx * pa_y[j] * pb_xyz[j] + 0.5 * pa_xxxy[j] * fl1_fx * pb_yz[j] + 0.5 * pa_xxx[j] * fl1_fx * pb_xxz[j] + 3.0 * pa_xxy[j] * fl1_fx * pb_xyz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxyz[j] + pa_xxxy[j] * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_20_30(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (20,30)

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

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxy = paDistances.data(34 * idx + 20);

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

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // Batch of Integrals (20,30)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     s_0_0, t_xxxy_xxzz, t_xxxy_xyyy, t_xxxy_xyyz, t_xxxy_xyzz, t_xxxy_xzzz, \
                                     t_xxxy_yyyy, t_xxxy_yyyz, t_xxxy_yyzz, t_xxxy_yzzz, t_xxxy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxy_xxzz[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 0.75 * fl3_fx * pa_y[j] * pb_x[j] + 0.25 * pa_xxxy[j] * fl2_fx + 1.5 * pa_xxy[j] * fl2_fx * pb_x[j] + 2.25 * pa_xy[j] * fl2_fx * pb_zz[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xx[j] + 1.5 * fl2_fx * pa_y[j] * pb_xzz[j] + 0.5 * pa_xxxy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxxy[j] * fl1_fx * pb_zz[j] + 3.0 * pa_xxy[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxzz[j] + pa_xxxy[j] * pb_xxzz[j]);

                t_xxxy_xyyy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * pa_xx[j] * fl3_fx + 1.125 * pa_x[j] * fl3_fx * pb_x[j] + 1.125 * fl3_fx * pa_y[j] * pb_y[j] + 1.125 * fl3_fx * pb_yy[j] + 0.75 * pa_xxx[j] * fl2_fx * pb_x[j] + 2.25 * pa_xxy[j] * fl2_fx * pb_y[j] + 2.25 * pa_xx[j] * fl2_fx * pb_yy[j] + 2.25 * pa_xy[j] * fl2_fx * pb_xy[j] + 2.25 * pa_x[j] * fl2_fx * pb_xyy[j] + 0.75 * fl2_fx * pa_y[j] * pb_yyy[j] + 1.5 * pa_xxxy[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xxx[j] * fl1_fx * pb_xyy[j] + 1.5 * pa_xxy[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xyyy[j] + pa_xxxy[j] * pb_xyyy[j]);

                t_xxxy_xyyz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_y[j] * pb_z[j] + 0.75 * fl3_fx * pb_yz[j] + 0.75 * pa_xxy[j] * fl2_fx * pb_z[j] + 1.5 * pa_xx[j] * fl2_fx * pb_yz[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xz[j] + 1.5 * pa_x[j] * fl2_fx * pb_xyz[j] + 0.75 * fl2_fx * pa_y[j] * pb_yyz[j] + 0.5 * pa_xxxy[j] * pb_xz[j] * fl1_fx + pa_xxx[j] * fl1_fx * pb_xyz[j] + 1.5 * pa_xxy[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xyyz[j] + pa_xxxy[j] * pb_xyyz[j]);

                t_xxxy_xyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pb_zz[j] + 0.25 * pa_xxx[j] * fl2_fx * pb_x[j] + 0.75 * pa_xxy[j] * fl2_fx * pb_y[j] + 0.75 * pa_xx[j] * fl2_fx * pb_zz[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xy[j] + 0.75 * pa_x[j] * fl2_fx * pb_xzz[j] + 0.75 * fl2_fx * pa_y[j] * pb_yzz[j] + 0.5 * pa_xxxy[j] * pb_xy[j] * fl1_fx + 0.5 * pa_xxx[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_xxy[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xyzz[j] + pa_xxxy[j] * pb_xyzz[j]);

                t_xxxy_xzzz[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_y[j] * pb_z[j] + 2.25 * pa_xxy[j] * fl2_fx * pb_z[j] + 2.25 * pa_xy[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_y[j] * pb_zzz[j] + 1.5 * pa_xxxy[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xxy[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xzzz[j] + pa_xxxy[j] * pb_xzzz[j]);

                t_xxxy_yyyy[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 4.5 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * pa_xxxy[j] * fl2_fx + 3.0 * pa_xxx[j] * fl2_fx * pb_y[j] + 4.5 * pa_xy[j] * fl2_fx * pb_yy[j] + 3.0 * pa_x[j] * fl2_fx * pb_yyy[j] + 3.0 * pa_xxxy[j] * pb_yy[j] * fl1_fx + 2.0 * pa_xxx[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_xy[j] * fl1_fx * pb_yyyy[j] + pa_xxxy[j] * pb_yyyy[j]);

                t_xxxy_yyyz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx * pb_z[j] + 0.75 * pa_xxx[j] * fl2_fx * pb_z[j] + 2.25 * pa_xy[j] * fl2_fx * pb_yz[j] + 2.25 * pa_x[j] * fl2_fx * pb_yyz[j] + 1.5 * pa_xxxy[j] * pb_yz[j] * fl1_fx + 1.5 * pa_xxx[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_yyyz[j] + pa_xxxy[j] * pb_yyyz[j]);

                t_xxxy_yyzz[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 0.75 * pa_x[j] * fl3_fx * pb_y[j] + 0.25 * pa_xxxy[j] * fl2_fx + 0.5 * pa_xxx[j] * fl2_fx * pb_y[j] + 0.75 * pa_xy[j] * fl2_fx * pb_yy[j] + 0.75 * pa_xy[j] * fl2_fx * pb_zz[j] + 1.5 * pa_x[j] * fl2_fx * pb_yzz[j] + 0.5 * pa_xxxy[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xxxy[j] * fl1_fx * pb_zz[j] + pa_xxx[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_yyzz[j] + pa_xxxy[j] * pb_yyzz[j]);

                t_xxxy_yzzz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx * pb_z[j] + 0.75 * pa_xxx[j] * fl2_fx * pb_z[j] + 2.25 * pa_xy[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_zzz[j] + 1.5 * pa_xxxy[j] * pb_yz[j] * fl1_fx + 0.5 * pa_xxx[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_yzzz[j] + pa_xxxy[j] * pb_yzzz[j]);

                t_xxxy_zzzz[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 0.75 * pa_xxxy[j] * fl2_fx + 4.5 * pa_xy[j] * fl2_fx * pb_zz[j] + 3.0 * pa_xxxy[j] * pb_zz[j] * fl1_fx + 1.5 * pa_xy[j] * fl1_fx * pb_zzzz[j] + pa_xxxy[j] * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_30_40(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (30,40)

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

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xz = paDistances.data(34 * idx + 5);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxz = paDistances.data(34 * idx + 21);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (30,40)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxxz_xxxx, t_xxxz_xxxy, t_xxxz_xxxz, \
                                     t_xxxz_xxyy, t_xxxz_xxyz, t_xxxz_xxzz, t_xxxz_xyyy, t_xxxz_xyyz, t_xxxz_xyzz, \
                                     t_xxxz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxz_xxxx[j] = fl_s_0_0 * (5.625 * pa_xz[j] * fl3_fx + 7.5 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * pa_xxxz[j] * fl2_fx + 9.0 * pa_xxz[j] * fl2_fx * pb_x[j] + 13.5 * pa_xz[j] * fl2_fx * pb_xx[j] + 3.0 * fl2_fx * pa_z[j] * pb_xxx[j] + 3.0 * pa_xxxz[j] * pb_xx[j] * fl1_fx + 6.0 * pa_xxz[j] * fl1_fx * pb_xxx[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxxx[j] + pa_xxxz[j] * pb_xxxx[j]);

                t_xxxz_xxxy[j] = fl_s_0_0 * (1.875 * fl3_fx * pa_z[j] * pb_y[j] + 2.25 * pa_xxz[j] * fl2_fx * pb_y[j] + 6.75 * pa_xz[j] * fl2_fx * pb_xy[j] + 2.25 * fl2_fx * pa_z[j] * pb_xxy[j] + 1.5 * pa_xxxz[j] * pb_xy[j] * fl1_fx + 4.5 * pa_xxz[j] * fl1_fx * pb_xxy[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxxy[j] + pa_xxxz[j] * pb_xxxy[j]);

                t_xxxz_xxxz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.125 * pa_xx[j] * fl3_fx + 3.375 * pa_x[j] * fl3_fx * pb_x[j] + 1.875 * fl3_fx * pa_z[j] * pb_z[j] + 1.125 * fl3_fx * pb_xx[j] + 0.75 * pa_xxx[j] * fl2_fx * pb_x[j] + 2.25 * pa_xxz[j] * fl2_fx * pb_z[j] + 2.25 * pa_xx[j] * fl2_fx * pb_xx[j] + 6.75 * pa_xz[j] * fl2_fx * pb_xz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xxx[j] + 2.25 * fl2_fx * pa_z[j] * pb_xxz[j] + 1.5 * pa_xxxz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xxx[j] * fl1_fx * pb_xxx[j] + 4.5 * pa_xxz[j] * fl1_fx * pb_xxz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxxz[j] + pa_xxxz[j] * pb_xxxz[j]);

                t_xxxz_xxyy[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 0.75 * fl3_fx * pa_z[j] * pb_x[j] + 0.25 * pa_xxxz[j] * fl2_fx + 1.5 * pa_xxz[j] * fl2_fx * pb_x[j] + 2.25 * pa_xz[j] * fl2_fx * pb_yy[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xx[j] + 1.5 * fl2_fx * pa_z[j] * pb_xyy[j] + 0.5 * pa_xxxz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxxz[j] * fl1_fx * pb_yy[j] + 3.0 * pa_xxz[j] * fl1_fx * pb_xyy[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxyy[j] + pa_xxxz[j] * pb_xxyy[j]);

                t_xxxz_xxyz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * fl3_fx * pb_xy[j] + 0.25 * pa_xxx[j] * fl2_fx * pb_y[j] + 1.5 * pa_xx[j] * fl2_fx * pb_xy[j] + 2.25 * pa_xz[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xxy[j] + 1.5 * fl2_fx * pa_z[j] * pb_xyz[j] + 0.5 * pa_xxxz[j] * fl1_fx * pb_yz[j] + 0.5 * pa_xxx[j] * fl1_fx * pb_xxy[j] + 3.0 * pa_xxz[j] * fl1_fx * pb_xyz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxyz[j] + pa_xxxz[j] * pb_xxyz[j]);

                t_xxxz_xxzz[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 2.25 * pa_x[j] * fl3_fx * pb_z[j] + 0.75 * fl3_fx * pa_z[j] * pb_x[j] + 1.5 * fl3_fx * pb_xz[j] + 0.25 * pa_xxxz[j] * fl2_fx + 0.5 * pa_xxx[j] * fl2_fx * pb_z[j] + 1.5 * pa_xxz[j] * fl2_fx * pb_x[j] + 3.0 * pa_xx[j] * fl2_fx * pb_xz[j] + 2.25 * pa_xz[j] * fl2_fx * pb_zz[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xx[j] + 1.5 * pa_x[j] * fl2_fx * pb_xxz[j] + 1.5 * fl2_fx * pa_z[j] * pb_xzz[j] + 0.5 * pa_xxxz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxxz[j] * fl1_fx * pb_zz[j] + pa_xxx[j] * fl1_fx * pb_xxz[j] + 3.0 * pa_xxz[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxzz[j] + pa_xxxz[j] * pb_xxzz[j]);

                t_xxxz_xyyy[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_z[j] * pb_y[j] + 2.25 * pa_xxz[j] * fl2_fx * pb_y[j] + 2.25 * pa_xz[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_z[j] * pb_yyy[j] + 1.5 * pa_xxxz[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xxz[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xyyy[j] + pa_xxxz[j] * pb_xyyy[j]);

                t_xxxz_xyyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_yy[j] + 0.25 * pa_xxx[j] * fl2_fx * pb_x[j] + 0.75 * pa_xxz[j] * fl2_fx * pb_z[j] + 0.75 * pa_xx[j] * fl2_fx * pb_yy[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xyy[j] + 0.75 * fl2_fx * pa_z[j] * pb_yyz[j] + 0.5 * pa_xxxz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xxx[j] * fl1_fx * pb_xyy[j] + 1.5 * pa_xxz[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xyyz[j] + pa_xxxz[j] * pb_xyyz[j]);

                t_xxxz_xyzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] * pb_y[j] + 0.75 * fl3_fx * pb_yz[j] + 0.75 * pa_xxz[j] * fl2_fx * pb_y[j] + 1.5 * pa_xx[j] * fl2_fx * pb_yz[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xy[j] + 1.5 * pa_x[j] * fl2_fx * pb_xyz[j] + 0.75 * fl2_fx * pa_z[j] * pb_yzz[j] + 0.5 * pa_xxxz[j] * pb_xy[j] * fl1_fx + pa_xxx[j] * fl1_fx * pb_xyz[j] + 1.5 * pa_xxz[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xyzz[j] + pa_xxxz[j] * pb_xyzz[j]);

                t_xxxz_xzzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * pa_xx[j] * fl3_fx + 1.125 * pa_x[j] * fl3_fx * pb_x[j] + 1.125 * fl3_fx * pa_z[j] * pb_z[j] + 1.125 * fl3_fx * pb_zz[j] + 0.75 * pa_xxx[j] * fl2_fx * pb_x[j] + 2.25 * pa_xxz[j] * fl2_fx * pb_z[j] + 2.25 * pa_xx[j] * fl2_fx * pb_zz[j] + 2.25 * pa_xz[j] * fl2_fx * pb_xz[j] + 2.25 * pa_x[j] * fl2_fx * pb_xzz[j] + 0.75 * fl2_fx * pa_z[j] * pb_zzz[j] + 1.5 * pa_xxxz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xxx[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_xxz[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xzzz[j] + pa_xxxz[j] * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_40_50(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (40,50)

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

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(34 * idx + 9);

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxxz = paDistances.data(34 * idx + 21);

            auto pa_xxyy = paDistances.data(34 * idx + 22);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (40,50)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_xz, pa_y, \
                                     pa_yy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xxxz_yyyy, \
                                     t_xxxz_yyyz, t_xxxz_yyzz, t_xxxz_yzzz, t_xxxz_zzzz, t_xxyy_xxxx, t_xxyy_xxxy, \
                                     t_xxyy_xxxz, t_xxyy_xxyy, t_xxyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxxz_yyyy[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 0.75 * pa_xxxz[j] * fl2_fx + 4.5 * pa_xz[j] * fl2_fx * pb_yy[j] + 3.0 * pa_xxxz[j] * pb_yy[j] * fl1_fx + 1.5 * pa_xz[j] * fl1_fx * pb_yyyy[j] + pa_xxxz[j] * pb_yyyy[j]);

                t_xxxz_yyyz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * pa_xxx[j] * fl2_fx * pb_y[j] + 2.25 * pa_xz[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_yyy[j] + 1.5 * pa_xxxz[j] * pb_yz[j] * fl1_fx + 0.5 * pa_xxx[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_xz[j] * fl1_fx * pb_yyyz[j] + pa_xxxz[j] * pb_yyyz[j]);

                t_xxxz_yyzz[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 0.75 * pa_x[j] * fl3_fx * pb_z[j] + 0.25 * pa_xxxz[j] * fl2_fx + 0.5 * pa_xxx[j] * fl2_fx * pb_z[j] + 0.75 * pa_xz[j] * fl2_fx * pb_yy[j] + 0.75 * pa_xz[j] * fl2_fx * pb_zz[j] + 1.5 * pa_x[j] * fl2_fx * pb_yyz[j] + 0.5 * pa_xxxz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xxxz[j] * fl1_fx * pb_zz[j] + pa_xxx[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_yyzz[j] + pa_xxxz[j] * pb_yyzz[j]);

                t_xxxz_yzzz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * pa_xxx[j] * fl2_fx * pb_y[j] + 2.25 * pa_xz[j] * fl2_fx * pb_yz[j] + 2.25 * pa_x[j] * fl2_fx * pb_yzz[j] + 1.5 * pa_xxxz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_xxx[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_yzzz[j] + pa_xxxz[j] * pb_yzzz[j]);

                t_xxxz_zzzz[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 4.5 * pa_x[j] * fl3_fx * pb_z[j] + 0.75 * pa_xxxz[j] * fl2_fx + 3.0 * pa_xxx[j] * fl2_fx * pb_z[j] + 4.5 * pa_xz[j] * fl2_fx * pb_zz[j] + 3.0 * pa_x[j] * fl2_fx * pb_zzz[j] + 3.0 * pa_xxxz[j] * pb_zz[j] * fl1_fx + 2.0 * pa_xxx[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_zzzz[j] + pa_xxxz[j] * pb_zzzz[j]);

                t_xxyy_xxxx[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * fl3_fx * pa_yy[j] + 0.375 * pa_xx[j] * fl3_fx + 3.0 * pa_x[j] * fl3_fx * pb_x[j] + 2.25 * fl3_fx * pb_xx[j] + 0.75 * pa_xxyy[j] * fl2_fx + 6.0 * pa_xyy[j] * fl2_fx * pb_x[j] + 4.5 * fl2_fx * pa_yy[j] * pb_xx[j] + 1.5 * pa_xx[j] * fl2_fx * pb_xx[j] + 2.0 * pa_x[j] * fl2_fx * pb_xxx[j] + 3.0 * pa_xxyy[j] * pb_xx[j] * fl1_fx + 4.0 * pa_xyy[j] * fl1_fx * pb_xxx[j] + 0.25 * fl2_fx * pb_xxxx[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxxx[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xxxx[j] + pa_xxyy[j] * pb_xxxx[j]);

                t_xxyy_xxxy[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl3_fx + 2.25 * fl3_fx * pa_y[j] * pb_x[j] + 0.75 * pa_x[j] * fl3_fx * pb_y[j] + 1.125 * fl3_fx * pb_xy[j] + 1.5 * pa_xxy[j] * fl2_fx * pb_x[j] + 1.5 * pa_xyy[j] * fl2_fx * pb_y[j] + 3.0 * pa_xy[j] * fl2_fx * pb_xx[j] + 2.25 * fl2_fx * pa_yy[j] * pb_xy[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xy[j] + 1.5 * pa_x[j] * fl2_fx * pb_xxy[j] + 0.5 * fl2_fx * pa_y[j] * pb_xxx[j] + 1.5 * pa_xxyy[j] * pb_xy[j] * fl1_fx + pa_xxy[j] * fl1_fx * pb_xxx[j] + 3.0 * pa_xyy[j] * fl1_fx * pb_xxy[j] + 0.25 * fl2_fx * pb_xxxy[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxxy[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xxxy[j] + pa_xxyy[j] * pb_xxxy[j]);

                t_xxyy_xxxz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx * pb_z[j] + 1.125 * fl3_fx * pb_xz[j] + 1.5 * pa_xyy[j] * fl2_fx * pb_z[j] + 2.25 * fl2_fx * pa_yy[j] * pb_xz[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xz[j] + 1.5 * pa_x[j] * fl2_fx * pb_xxz[j] + 1.5 * pa_xxyy[j] * pb_xz[j] * fl1_fx + 3.0 * pa_xyy[j] * fl1_fx * pb_xxz[j] + 0.25 * fl2_fx * pb_xxxz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxxz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xxxz[j] + pa_xxyy[j] * pb_xxxz[j]);

                t_xxyy_xxyy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 1.5 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_yy[j] + 1.5 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pb_xx[j] + 0.375 * fl3_fx * pb_yy[j] + 0.25 * pa_xxyy[j] * fl2_fx + pa_xxy[j] * fl2_fx * pb_y[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xx[j] + pa_xyy[j] * fl2_fx * pb_x[j] + 4.0 * pa_xy[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_yy[j] * pb_yy[j] + 0.25 * pa_xx[j] * fl2_fx * pb_yy[j] + pa_x[j] * fl2_fx * pb_xyy[j] + 0.25 * fl2_fx * pa_yy[j] * pb_xx[j] + fl2_fx * pa_y[j] * pb_xxy[j] + 0.5 * pa_xxyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxyy[j] * fl1_fx * pb_yy[j] + 2.0 * pa_xxy[j] * fl1_fx * pb_xxy[j] + 2.0 * pa_xyy[j] * fl1_fx * pb_xyy[j] + 0.25 * fl2_fx * pb_xxyy[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxyy[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xxyy[j] + pa_xxyy[j] * pb_xxyy[j]);

                t_xxyy_xxyz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_y[j] * pb_z[j] + 0.375 * fl3_fx * pb_yz[j] + 0.5 * pa_xxy[j] * fl2_fx * pb_z[j] + 2.0 * pa_xy[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_yy[j] * pb_yz[j] + 0.25 * pa_xx[j] * fl2_fx * pb_yz[j] + pa_x[j] * fl2_fx * pb_xyz[j] + 0.5 * fl2_fx * pa_y[j] * pb_xxz[j] + 0.5 * pa_xxyy[j] * fl1_fx * pb_yz[j] + pa_xxy[j] * fl1_fx * pb_xxz[j] + 2.0 * pa_xyy[j] * fl1_fx * pb_xyz[j] + 0.25 * fl2_fx * pb_xxyz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxyz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xxyz[j] + pa_xxyy[j] * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_50_60(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (50,60)

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

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_yy = paDistances.data(34 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xyy = paDistances.data(34 * idx + 12);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxyy = paDistances.data(34 * idx + 22);

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

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // Batch of Integrals (50,60)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, \
                                     pb_zzzz, s_0_0, t_xxyy_xxzz, t_xxyy_xyyy, t_xxyy_xyyz, t_xxyy_xyzz, t_xxyy_xzzz, \
                                     t_xxyy_yyyy, t_xxyy_yyyz, t_xxyy_yyzz, t_xxyy_yzzz, t_xxyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxyy_xxzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_yy[j] + 0.125 * pa_xx[j] * fl3_fx + 0.5 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pb_zz[j] + 0.25 * pa_xxyy[j] * fl2_fx + pa_xyy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yy[j] * pb_zz[j] + 0.125 * fl3_fx * pb_xx[j] + 0.25 * pa_xx[j] * fl2_fx * pb_xx[j] + 0.25 * pa_xx[j] * fl2_fx * pb_zz[j] + pa_x[j] * fl2_fx * pb_xzz[j] + 0.25 * fl2_fx * pa_yy[j] * pb_xx[j] + 0.5 * pa_xxyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxyy[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xyy[j] * fl1_fx * pb_xzz[j] + 0.25 * fl2_fx * pb_xxzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxzz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xxzz[j] + pa_xxyy[j] * pb_xxzz[j]);

                t_xxyy_xyyy[j] = fl_s_0_0 * (1.5 * pa_xy[j] * fl3_fx + 2.25 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * fl3_fx * pa_y[j] * pb_x[j] + 1.125 * fl3_fx * pb_xy[j] + 1.5 * pa_xxy[j] * fl2_fx * pb_x[j] + 2.25 * pa_xx[j] * fl2_fx * pb_xy[j] + 1.5 * pa_xyy[j] * fl2_fx * pb_y[j] + 3.0 * pa_xy[j] * fl2_fx * pb_yy[j] + 0.5 * pa_x[j] * fl2_fx * pb_yyy[j] + 0.75 * fl2_fx * pa_yy[j] * pb_xy[j] + 1.5 * fl2_fx * pa_y[j] * pb_xyy[j] + 1.5 * pa_xxyy[j] * pb_xy[j] * fl1_fx + 3.0 * pa_xxy[j] * fl1_fx * pb_xyy[j] + pa_xyy[j] * fl1_fx * pb_yyy[j] + 0.25 * fl2_fx * pb_xyyy[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xyyy[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xyyy[j] + pa_xxyy[j] * pb_xyyy[j]);

                t_xxyy_xyyz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx * pb_z[j] + 0.375 * fl3_fx * pb_xz[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xz[j] + 0.5 * pa_xyy[j] * fl2_fx * pb_z[j] + 2.0 * pa_xy[j] * fl2_fx * pb_yz[j] + 0.5 * pa_x[j] * fl2_fx * pb_yyz[j] + 0.25 * fl2_fx * pa_yy[j] * pb_xz[j] + fl2_fx * pa_y[j] * pb_xyz[j] + 0.5 * pa_xxyy[j] * pb_xz[j] * fl1_fx + 2.0 * pa_xxy[j] * fl1_fx * pb_xyz[j] + pa_xyy[j] * fl1_fx * pb_yyz[j] + 0.25 * fl2_fx * pb_xyyz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xyyz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xyyz[j] + pa_xxyy[j] * pb_xyyz[j]);

                t_xxyy_xyzz[j] = fl_s_0_0 * (0.5 * pa_xy[j] * fl3_fx + 0.25 * pa_x[j] * fl3_fx * pb_y[j] + 0.25 * fl3_fx * pa_y[j] * pb_x[j] + 0.5 * pa_xxy[j] * fl2_fx * pb_x[j] + 0.5 * pa_xyy[j] * fl2_fx * pb_y[j] + pa_xy[j] * fl2_fx * pb_zz[j] + 0.125 * fl3_fx * pb_xy[j] + 0.25 * pa_xx[j] * fl2_fx * pb_xy[j] + 0.5 * pa_x[j] * fl2_fx * pb_yzz[j] + 0.25 * fl2_fx * pa_yy[j] * pb_xy[j] + 0.5 * fl2_fx * pa_y[j] * pb_xzz[j] + 0.5 * pa_xxyy[j] * pb_xy[j] * fl1_fx + pa_xxy[j] * fl1_fx * pb_xzz[j] + pa_xyy[j] * fl1_fx * pb_yzz[j] + 0.25 * fl2_fx * pb_xyzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xyzz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xyzz[j] + pa_xxyy[j] * pb_xyzz[j]);

                t_xxyy_xzzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx * pb_z[j] + 1.5 * pa_xyy[j] * fl2_fx * pb_z[j] + 0.375 * fl3_fx * pb_xz[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xz[j] + 0.5 * pa_x[j] * fl2_fx * pb_zzz[j] + 0.75 * fl2_fx * pa_yy[j] * pb_xz[j] + 1.5 * pa_xxyy[j] * pb_xz[j] * fl1_fx + pa_xyy[j] * fl1_fx * pb_zzz[j] + 0.25 * fl2_fx * pb_xzzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xzzz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_xzzz[j] + pa_xxyy[j] * pb_xzzz[j]);

                t_xxyy_yyyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_xx[j] * fl3_fx + 0.375 * fl3_fx * pa_yy[j] + 3.0 * fl3_fx * pa_y[j] * pb_y[j] + 2.25 * fl3_fx * pb_yy[j] + 0.75 * pa_xxyy[j] * fl2_fx + 6.0 * pa_xxy[j] * fl2_fx * pb_y[j] + 4.5 * pa_xx[j] * fl2_fx * pb_yy[j] + 1.5 * fl2_fx * pa_yy[j] * pb_yy[j] + 2.0 * fl2_fx * pa_y[j] * pb_yyy[j] + 3.0 * pa_xxyy[j] * pb_yy[j] * fl1_fx + 4.0 * pa_xxy[j] * fl1_fx * pb_yyy[j] + 0.25 * fl2_fx * pb_yyyy[j] + 0.5 * pa_xx[j] * fl1_fx * pb_yyyy[j] + 0.5 * fl1_fx * pa_yy[j] * pb_yyyy[j] + pa_xxyy[j] * pb_yyyy[j]);

                t_xxyy_yyyz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_y[j] * pb_z[j] + 1.125 * fl3_fx * pb_yz[j] + 1.5 * pa_xxy[j] * fl2_fx * pb_z[j] + 2.25 * pa_xx[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pa_yy[j] * pb_yz[j] + 1.5 * fl2_fx * pa_y[j] * pb_yyz[j] + 1.5 * pa_xxyy[j] * pb_yz[j] * fl1_fx + 3.0 * pa_xxy[j] * fl1_fx * pb_yyz[j] + 0.25 * fl2_fx * pb_yyyz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_yyyz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_yyyz[j] + pa_xxyy[j] * pb_yyyz[j]);

                t_xxyy_yyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.125 * fl3_fx * pa_yy[j] + 0.5 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pb_zz[j] + 0.25 * pa_xxyy[j] * fl2_fx + pa_xxy[j] * fl2_fx * pb_y[j] + 0.75 * pa_xx[j] * fl2_fx * pb_zz[j] + 0.125 * fl3_fx * pb_yy[j] + 0.25 * pa_xx[j] * fl2_fx * pb_yy[j] + 0.25 * fl2_fx * pa_yy[j] * pb_yy[j] + 0.25 * fl2_fx * pa_yy[j] * pb_zz[j] + fl2_fx * pa_y[j] * pb_yzz[j] + 0.5 * pa_xxyy[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xxyy[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xxy[j] * fl1_fx * pb_yzz[j] + 0.25 * fl2_fx * pb_yyzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_yyzz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_yyzz[j] + pa_xxyy[j] * pb_yyzz[j]);

                t_xxyy_yzzz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_y[j] * pb_z[j] + 1.5 * pa_xxy[j] * fl2_fx * pb_z[j] + 0.375 * fl3_fx * pb_yz[j] + 0.75 * pa_xx[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pa_yy[j] * pb_yz[j] + 0.5 * fl2_fx * pa_y[j] * pb_zzz[j] + 1.5 * pa_xxyy[j] * pb_yz[j] * fl1_fx + pa_xxy[j] * fl1_fx * pb_zzz[j] + 0.25 * fl2_fx * pb_yzzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_yzzz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_yzzz[j] + pa_xxyy[j] * pb_yzzz[j]);

                t_xxyy_zzzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * fl3_fx * pa_yy[j] + 0.75 * pa_xxyy[j] * fl2_fx + 0.75 * fl3_fx * pb_zz[j] + 1.5 * pa_xx[j] * fl2_fx * pb_zz[j] + 1.5 * fl2_fx * pa_yy[j] * pb_zz[j] + 3.0 * pa_xxyy[j] * pb_zz[j] * fl1_fx + 0.25 * fl2_fx * pb_zzzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_zzzz[j] + 0.5 * fl1_fx * pa_yy[j] * pb_zzzz[j] + pa_xxyy[j] * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_60_70(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (60,70)

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

            auto pa_yz = paDistances.data(34 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxyz = paDistances.data(34 * idx + 23);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (60,70)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, \
                                     pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxyz_xxxx, \
                                     t_xxyz_xxxy, t_xxyz_xxxz, t_xxyz_xxyy, t_xxyz_xxyz, t_xxyz_xxzz, t_xxyz_xyyy, \
                                     t_xxyz_xyyz, t_xxyz_xyzz, t_xxyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxyz_xxxx[j] = fl_s_0_0 * (1.875 * fl3_fx * pa_yz[j] + 0.75 * pa_xxyz[j] * fl2_fx + 6.0 * pa_xyz[j] * fl2_fx * pb_x[j] + 4.5 * fl2_fx * pa_yz[j] * pb_xx[j] + 3.0 * pa_xxyz[j] * pb_xx[j] * fl1_fx + 4.0 * pa_xyz[j] * fl1_fx * pb_xxx[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xxxx[j] + pa_xxyz[j] * pb_xxxx[j]);

                t_xxyz_xxxy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl3_fx + 1.125 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * pa_xxz[j] * fl2_fx * pb_x[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_y[j] + 1.5 * pa_xz[j] * fl2_fx * pb_xx[j] + 2.25 * fl2_fx * pa_yz[j] * pb_xy[j] + 0.25 * fl2_fx * pa_z[j] * pb_xxx[j] + 1.5 * pa_xxyz[j] * pb_xy[j] * fl1_fx + 0.5 * pa_xxz[j] * fl1_fx * pb_xxx[j] + 3.0 * pa_xyz[j] * fl1_fx * pb_xxy[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xxxy[j] + pa_xxyz[j] * pb_xxxy[j]);

                t_xxyz_xxxz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl3_fx + 1.125 * fl3_fx * pa_y[j] * pb_x[j] + 0.75 * pa_xxy[j] * fl2_fx * pb_x[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_z[j] + 1.5 * pa_xy[j] * fl2_fx * pb_xx[j] + 2.25 * fl2_fx * pa_yz[j] * pb_xz[j] + 0.25 * fl2_fx * pa_y[j] * pb_xxx[j] + 1.5 * pa_xxyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xxy[j] * fl1_fx * pb_xxx[j] + 3.0 * pa_xyz[j] * fl1_fx * pb_xxz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xxxz[j] + pa_xxyz[j] * pb_xxxz[j]);

                t_xxyz_xxyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_yz[j] + 0.75 * fl3_fx * pa_z[j] * pb_y[j] + 0.25 * pa_xxyz[j] * fl2_fx + 0.5 * pa_xxz[j] * fl2_fx * pb_y[j] + pa_xyz[j] * fl2_fx * pb_x[j] + 2.0 * pa_xz[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_yz[j] * pb_yy[j] + 0.25 * fl2_fx * pa_yz[j] * pb_xx[j] + 0.5 * fl2_fx * pa_z[j] * pb_xxy[j] + 0.5 * pa_xxyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxyz[j] * fl1_fx * pb_yy[j] + pa_xxz[j] * fl1_fx * pb_xxy[j] + 2.0 * pa_xyz[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xxyy[j] + pa_xxyz[j] * pb_xxyy[j]);

                t_xxyz_xxyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.125 * pa_xx[j] * fl3_fx + 0.5 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.125 * fl3_fx * pb_xx[j] + 0.25 * pa_xxy[j] * fl2_fx * pb_y[j] + 0.25 * pa_xxz[j] * fl2_fx * pb_z[j] + 0.25 * pa_xx[j] * fl2_fx * pb_xx[j] + pa_xy[j] * fl2_fx * pb_xy[j] + pa_xz[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_yz[j] * pb_yz[j] + 0.25 * fl2_fx * pa_y[j] * pb_xxy[j] + 0.25 * fl2_fx * pa_z[j] * pb_xxz[j] + 0.5 * pa_xxyz[j] * fl1_fx * pb_yz[j] + 0.5 * pa_xxy[j] * fl1_fx * pb_xxy[j] + 0.5 * pa_xxz[j] * fl1_fx * pb_xxz[j] + 2.0 * pa_xyz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xxyz[j] + pa_xxyz[j] * pb_xxyz[j]);

                t_xxyz_xxzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_yz[j] + 0.75 * fl3_fx * pa_y[j] * pb_z[j] + 0.25 * pa_xxyz[j] * fl2_fx + 0.5 * pa_xxy[j] * fl2_fx * pb_z[j] + pa_xyz[j] * fl2_fx * pb_x[j] + 2.0 * pa_xy[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_yz[j] * pb_zz[j] + 0.25 * fl2_fx * pa_yz[j] * pb_xx[j] + 0.5 * fl2_fx * pa_y[j] * pb_xxz[j] + 0.5 * pa_xxyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxyz[j] * fl1_fx * pb_zz[j] + pa_xxy[j] * fl1_fx * pb_xxz[j] + 2.0 * pa_xyz[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xxzz[j] + pa_xxyz[j] * pb_xxzz[j]);

                t_xxyz_xyyy[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl3_fx + 0.375 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * pa_xxz[j] * fl2_fx * pb_x[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_y[j] + 1.5 * pa_xz[j] * fl2_fx * pb_yy[j] + 0.75 * fl2_fx * pa_yz[j] * pb_xy[j] + 0.75 * fl2_fx * pa_z[j] * pb_xyy[j] + 1.5 * pa_xxyz[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xxz[j] * fl1_fx * pb_xyy[j] + pa_xyz[j] * fl1_fx * pb_yyy[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xyyy[j] + pa_xxyz[j] * pb_xyyy[j]);

                t_xxyz_xyyz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl3_fx + 0.5 * pa_x[j] * fl3_fx * pb_y[j] + 0.125 * fl3_fx * pa_y[j] * pb_x[j] + 0.25 * fl3_fx * pb_xy[j] + 0.25 * pa_xxy[j] * fl2_fx * pb_x[j] + 0.5 * pa_xx[j] * fl2_fx * pb_xy[j] + 0.5 * pa_xyz[j] * fl2_fx * pb_z[j] + 0.5 * pa_xy[j] * fl2_fx * pb_yy[j] + pa_xz[j] * fl2_fx * pb_yz[j] + 0.25 * fl2_fx * pa_yz[j] * pb_xz[j] + 0.25 * fl2_fx * pa_y[j] * pb_xyy[j] + 0.5 * fl2_fx * pa_z[j] * pb_xyz[j] + 0.5 * pa_xxyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xxy[j] * fl1_fx * pb_xyy[j] + pa_xxz[j] * fl1_fx * pb_xyz[j] + pa_xyz[j] * fl1_fx * pb_yyz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xyyz[j] + pa_xxyz[j] * pb_xyyz[j]);

                t_xxyz_xyzz[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl3_fx + 0.5 * pa_x[j] * fl3_fx * pb_z[j] + 0.125 * fl3_fx * pa_z[j] * pb_x[j] + 0.25 * fl3_fx * pb_xz[j] + 0.25 * pa_xxz[j] * fl2_fx * pb_x[j] + 0.5 * pa_xx[j] * fl2_fx * pb_xz[j] + 0.5 * pa_xyz[j] * fl2_fx * pb_y[j] + pa_xy[j] * fl2_fx * pb_yz[j] + 0.5 * pa_xz[j] * fl2_fx * pb_zz[j] + 0.25 * fl2_fx * pa_yz[j] * pb_xy[j] + 0.5 * fl2_fx * pa_y[j] * pb_xyz[j] + 0.25 * fl2_fx * pa_z[j] * pb_xzz[j] + 0.5 * pa_xxyz[j] * pb_xy[j] * fl1_fx + pa_xxy[j] * fl1_fx * pb_xyz[j] + 0.5 * pa_xxz[j] * fl1_fx * pb_xzz[j] + pa_xyz[j] * fl1_fx * pb_yzz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xyzz[j] + pa_xxyz[j] * pb_xyzz[j]);

                t_xxyz_xzzz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl3_fx + 0.375 * fl3_fx * pa_y[j] * pb_x[j] + 0.75 * pa_xxy[j] * fl2_fx * pb_x[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_z[j] + 1.5 * pa_xy[j] * fl2_fx * pb_zz[j] + 0.75 * fl2_fx * pa_yz[j] * pb_xz[j] + 0.75 * fl2_fx * pa_y[j] * pb_xzz[j] + 1.5 * pa_xxyz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xxy[j] * fl1_fx * pb_xzz[j] + pa_xyz[j] * fl1_fx * pb_zzz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_xzzz[j] + pa_xxyz[j] * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_70_80(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (70,80)

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

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxy = paDistances.data(34 * idx + 10);

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxyz = paDistances.data(34 * idx + 23);

            auto pa_xxzz = paDistances.data(34 * idx + 24);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (70,80)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_y, pa_yz, \
                                     pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, \
                                     t_xxyz_yyyy, t_xxyz_yyyz, t_xxyz_yyzz, t_xxyz_yzzz, t_xxyz_zzzz, t_xxzz_xxxx, \
                                     t_xxzz_xxxy, t_xxzz_xxxz, t_xxzz_xxyy, t_xxzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxyz_yyyy[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_yz[j] + 1.5 * fl3_fx * pa_z[j] * pb_y[j] + 0.75 * pa_xxyz[j] * fl2_fx + 3.0 * pa_xxz[j] * fl2_fx * pb_y[j] + 1.5 * fl2_fx * pa_yz[j] * pb_yy[j] + fl2_fx * pa_z[j] * pb_yyy[j] + 3.0 * pa_xxyz[j] * pb_yy[j] * fl1_fx + 2.0 * pa_xxz[j] * fl1_fx * pb_yyy[j] + 0.5 * fl1_fx * pa_yz[j] * pb_yyyy[j] + pa_xxyz[j] * pb_yyyy[j]);

                t_xxyz_yyyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_yy[j] + 0.75 * pa_xxy[j] * fl2_fx * pb_y[j] + 0.75 * pa_xxz[j] * fl2_fx * pb_z[j] + 0.75 * pa_xx[j] * fl2_fx * pb_yy[j] + 0.75 * fl2_fx * pa_yz[j] * pb_yz[j] + 0.25 * fl2_fx * pa_y[j] * pb_yyy[j] + 0.75 * fl2_fx * pa_z[j] * pb_yyz[j] + 1.5 * pa_xxyz[j] * pb_yz[j] * fl1_fx + 0.5 * pa_xxy[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_xxz[j] * fl1_fx * pb_yyz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_yyyz[j] + pa_xxyz[j] * pb_yyyz[j]);

                t_xxyz_yyzz[j] = fl_s_0_0 * (0.125 * fl3_fx * pa_yz[j] + 0.25 * fl3_fx * pa_y[j] * pb_z[j] + 0.25 * fl3_fx * pa_z[j] * pb_y[j] + 0.5 * fl3_fx * pb_yz[j] + 0.25 * pa_xxyz[j] * fl2_fx + 0.5 * pa_xxy[j] * fl2_fx * pb_z[j] + 0.5 * pa_xxz[j] * fl2_fx * pb_y[j] + pa_xx[j] * fl2_fx * pb_yz[j] + 0.25 * fl2_fx * pa_yz[j] * pb_yy[j] + 0.25 * fl2_fx * pa_yz[j] * pb_zz[j] + 0.5 * fl2_fx * pa_y[j] * pb_yyz[j] + 0.5 * fl2_fx * pa_z[j] * pb_yzz[j] + 0.5 * pa_xxyz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xxyz[j] * fl1_fx * pb_zz[j] + pa_xxy[j] * fl1_fx * pb_yyz[j] + pa_xxz[j] * fl1_fx * pb_yzz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_yyzz[j] + pa_xxyz[j] * pb_yyzz[j]);

                t_xxyz_yzzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_zz[j] + 0.75 * pa_xxy[j] * fl2_fx * pb_y[j] + 0.75 * pa_xxz[j] * fl2_fx * pb_z[j] + 0.75 * pa_xx[j] * fl2_fx * pb_zz[j] + 0.75 * fl2_fx * pa_yz[j] * pb_yz[j] + 0.75 * fl2_fx * pa_y[j] * pb_yzz[j] + 0.25 * fl2_fx * pa_z[j] * pb_zzz[j] + 1.5 * pa_xxyz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_xxy[j] * fl1_fx * pb_yzz[j] + 0.5 * pa_xxz[j] * fl1_fx * pb_zzz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_yzzz[j] + pa_xxyz[j] * pb_yzzz[j]);

                t_xxyz_zzzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_yz[j] + 1.5 * fl3_fx * pa_y[j] * pb_z[j] + 0.75 * pa_xxyz[j] * fl2_fx + 3.0 * pa_xxy[j] * fl2_fx * pb_z[j] + 1.5 * fl2_fx * pa_yz[j] * pb_zz[j] + fl2_fx * pa_y[j] * pb_zzz[j] + 3.0 * pa_xxyz[j] * pb_zz[j] * fl1_fx + 2.0 * pa_xxy[j] * fl1_fx * pb_zzz[j] + 0.5 * fl1_fx * pa_yz[j] * pb_zzzz[j] + pa_xxyz[j] * pb_zzzz[j]);

                t_xxzz_xxxx[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * fl3_fx * pa_zz[j] + 0.375 * pa_xx[j] * fl3_fx + 3.0 * pa_x[j] * fl3_fx * pb_x[j] + 2.25 * fl3_fx * pb_xx[j] + 0.75 * pa_xxzz[j] * fl2_fx + 6.0 * pa_xzz[j] * fl2_fx * pb_x[j] + 4.5 * fl2_fx * pa_zz[j] * pb_xx[j] + 1.5 * pa_xx[j] * fl2_fx * pb_xx[j] + 2.0 * pa_x[j] * fl2_fx * pb_xxx[j] + 3.0 * pa_xxzz[j] * pb_xx[j] * fl1_fx + 4.0 * pa_xzz[j] * fl1_fx * pb_xxx[j] + 0.25 * fl2_fx * pb_xxxx[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxxx[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxxx[j] + pa_xxzz[j] * pb_xxxx[j]);

                t_xxzz_xxxy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx * pb_y[j] + 1.125 * fl3_fx * pb_xy[j] + 1.5 * pa_xzz[j] * fl2_fx * pb_y[j] + 2.25 * fl2_fx * pa_zz[j] * pb_xy[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xy[j] + 1.5 * pa_x[j] * fl2_fx * pb_xxy[j] + 1.5 * pa_xxzz[j] * pb_xy[j] * fl1_fx + 3.0 * pa_xzz[j] * fl1_fx * pb_xxy[j] + 0.25 * fl2_fx * pb_xxxy[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxxy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxxy[j] + pa_xxzz[j] * pb_xxxy[j]);

                t_xxzz_xxxz[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl3_fx + 2.25 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * pa_x[j] * fl3_fx * pb_z[j] + 1.125 * fl3_fx * pb_xz[j] + 1.5 * pa_xxz[j] * fl2_fx * pb_x[j] + 1.5 * pa_xzz[j] * fl2_fx * pb_z[j] + 3.0 * pa_xz[j] * fl2_fx * pb_xx[j] + 2.25 * fl2_fx * pa_zz[j] * pb_xz[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xz[j] + 1.5 * pa_x[j] * fl2_fx * pb_xxz[j] + 0.5 * fl2_fx * pa_z[j] * pb_xxx[j] + 1.5 * pa_xxzz[j] * pb_xz[j] * fl1_fx + pa_xxz[j] * fl1_fx * pb_xxx[j] + 3.0 * pa_xzz[j] * fl1_fx * pb_xxz[j] + 0.25 * fl2_fx * pb_xxxz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxxz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxxz[j] + pa_xxzz[j] * pb_xxxz[j]);

                t_xxzz_xxyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_zz[j] + 0.125 * pa_xx[j] * fl3_fx + 0.5 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pb_yy[j] + 0.25 * pa_xxzz[j] * fl2_fx + pa_xzz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_zz[j] * pb_yy[j] + 0.125 * fl3_fx * pb_xx[j] + 0.25 * pa_xx[j] * fl2_fx * pb_xx[j] + 0.25 * pa_xx[j] * fl2_fx * pb_yy[j] + pa_x[j] * fl2_fx * pb_xyy[j] + 0.25 * fl2_fx * pa_zz[j] * pb_xx[j] + 0.5 * pa_xxzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxzz[j] * fl1_fx * pb_yy[j] + 2.0 * pa_xzz[j] * fl1_fx * pb_xyy[j] + 0.25 * fl2_fx * pb_xxyy[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxyy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxyy[j] + pa_xxzz[j] * pb_xxyy[j]);

                t_xxzz_xxyz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] * pb_y[j] + 0.375 * fl3_fx * pb_yz[j] + 0.5 * pa_xxz[j] * fl2_fx * pb_y[j] + 2.0 * pa_xz[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_zz[j] * pb_yz[j] + 0.25 * pa_xx[j] * fl2_fx * pb_yz[j] + pa_x[j] * fl2_fx * pb_xyz[j] + 0.5 * fl2_fx * pa_z[j] * pb_xxy[j] + 0.5 * pa_xxzz[j] * fl1_fx * pb_yz[j] + pa_xxz[j] * fl1_fx * pb_xxy[j] + 2.0 * pa_xzz[j] * fl1_fx * pb_xyz[j] + 0.25 * fl2_fx * pb_xxyz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxyz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxyz[j] + pa_xxzz[j] * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_80_90(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        // Batch of Integrals (80,90)

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

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(34 * idx + 3);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxz = paDistances.data(34 * idx + 11);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xxzz = paDistances.data(34 * idx + 24);

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

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // Batch of Integrals (80,90)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, \
                                     pb_zzz, pb_zzzz, s_0_0, t_xxzz_xxzz, t_xxzz_xyyy, t_xxzz_xyyz, t_xxzz_xyzz, \
                                     t_xxzz_xzzz, t_xxzz_yyyy, t_xxzz_yyyz, t_xxzz_yyzz, t_xxzz_yzzz, t_xxzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xxzz_xxzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 1.5 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_zz[j] + 1.5 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_xx[j] + 0.375 * fl3_fx * pb_zz[j] + 0.25 * pa_xxzz[j] * fl2_fx + pa_xxz[j] * fl2_fx * pb_z[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xx[j] + pa_xzz[j] * fl2_fx * pb_x[j] + 4.0 * pa_xz[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_zz[j] + 0.25 * pa_xx[j] * fl2_fx * pb_zz[j] + pa_x[j] * fl2_fx * pb_xzz[j] + 0.25 * fl2_fx * pa_zz[j] * pb_xx[j] + fl2_fx * pa_z[j] * pb_xxz[j] + 0.5 * pa_xxzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xxzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xxz[j] * fl1_fx * pb_xxz[j] + 2.0 * pa_xzz[j] * fl1_fx * pb_xzz[j] + 0.25 * fl2_fx * pb_xxzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xxzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxzz[j] + pa_xxzz[j] * pb_xxzz[j]);

                t_xxzz_xyyy[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx * pb_y[j] + 1.5 * pa_xzz[j] * fl2_fx * pb_y[j] + 0.375 * fl3_fx * pb_xy[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xy[j] + 0.5 * pa_x[j] * fl2_fx * pb_yyy[j] + 0.75 * fl2_fx * pa_zz[j] * pb_xy[j] + 1.5 * pa_xxzz[j] * pb_xy[j] * fl1_fx + pa_xzz[j] * fl1_fx * pb_yyy[j] + 0.25 * fl2_fx * pb_xyyy[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xyyy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xyyy[j] + pa_xxzz[j] * pb_xyyy[j]);

                t_xxzz_xyyz[j] = fl_s_0_0 * (0.5 * pa_xz[j] * fl3_fx + 0.25 * pa_x[j] * fl3_fx * pb_z[j] + 0.25 * fl3_fx * pa_z[j] * pb_x[j] + 0.5 * pa_xxz[j] * fl2_fx * pb_x[j] + 0.5 * pa_xzz[j] * fl2_fx * pb_z[j] + pa_xz[j] * fl2_fx * pb_yy[j] + 0.125 * fl3_fx * pb_xz[j] + 0.25 * pa_xx[j] * fl2_fx * pb_xz[j] + 0.5 * pa_x[j] * fl2_fx * pb_yyz[j] + 0.25 * fl2_fx * pa_zz[j] * pb_xz[j] + 0.5 * fl2_fx * pa_z[j] * pb_xyy[j] + 0.5 * pa_xxzz[j] * pb_xz[j] * fl1_fx + pa_xxz[j] * fl1_fx * pb_xyy[j] + pa_xzz[j] * fl1_fx * pb_yyz[j] + 0.25 * fl2_fx * pb_xyyz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xyyz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xyyz[j] + pa_xxzz[j] * pb_xyyz[j]);

                t_xxzz_xyzz[j] = fl_s_0_0 * (0.75 * pa_x[j] * fl3_fx * pb_y[j] + 0.375 * fl3_fx * pb_xy[j] + 0.75 * pa_xx[j] * fl2_fx * pb_xy[j] + 0.5 * pa_xzz[j] * fl2_fx * pb_y[j] + 2.0 * pa_xz[j] * fl2_fx * pb_yz[j] + 0.5 * pa_x[j] * fl2_fx * pb_yzz[j] + 0.25 * fl2_fx * pa_zz[j] * pb_xy[j] + fl2_fx * pa_z[j] * pb_xyz[j] + 0.5 * pa_xxzz[j] * pb_xy[j] * fl1_fx + 2.0 * pa_xxz[j] * fl1_fx * pb_xyz[j] + pa_xzz[j] * fl1_fx * pb_yzz[j] + 0.25 * fl2_fx * pb_xyzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xyzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xyzz[j] + pa_xxzz[j] * pb_xyzz[j]);

                t_xxzz_xzzz[j] = fl_s_0_0 * (1.5 * pa_xz[j] * fl3_fx + 2.25 * pa_x[j] * fl3_fx * pb_z[j] + 0.75 * fl3_fx * pa_z[j] * pb_x[j] + 1.125 * fl3_fx * pb_xz[j] + 1.5 * pa_xxz[j] * fl2_fx * pb_x[j] + 2.25 * pa_xx[j] * fl2_fx * pb_xz[j] + 1.5 * pa_xzz[j] * fl2_fx * pb_z[j] + 3.0 * pa_xz[j] * fl2_fx * pb_zz[j] + 0.5 * pa_x[j] * fl2_fx * pb_zzz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_xz[j] + 1.5 * fl2_fx * pa_z[j] * pb_xzz[j] + 1.5 * pa_xxzz[j] * pb_xz[j] * fl1_fx + 3.0 * pa_xxz[j] * fl1_fx * pb_xzz[j] + pa_xzz[j] * fl1_fx * pb_zzz[j] + 0.25 * fl2_fx * pb_xzzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_xzzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xzzz[j] + pa_xxzz[j] * pb_xzzz[j]);

                t_xxzz_yyyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.375 * fl3_fx * pa_zz[j] + 0.75 * pa_xxzz[j] * fl2_fx + 0.75 * fl3_fx * pb_yy[j] + 1.5 * pa_xx[j] * fl2_fx * pb_yy[j] + 1.5 * fl2_fx * pa_zz[j] * pb_yy[j] + 3.0 * pa_xxzz[j] * pb_yy[j] * fl1_fx + 0.25 * fl2_fx * pb_yyyy[j] + 0.5 * pa_xx[j] * fl1_fx * pb_yyyy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yyyy[j] + pa_xxzz[j] * pb_yyyy[j]);

                t_xxzz_yyyz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] * pb_y[j] + 1.5 * pa_xxz[j] * fl2_fx * pb_y[j] + 0.375 * fl3_fx * pb_yz[j] + 0.75 * pa_xx[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_yz[j] + 0.5 * fl2_fx * pa_z[j] * pb_yyy[j] + 1.5 * pa_xxzz[j] * pb_yz[j] * fl1_fx + pa_xxz[j] * fl1_fx * pb_yyy[j] + 0.25 * fl2_fx * pb_yyyz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_yyyz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yyyz[j] + pa_xxzz[j] * pb_yyyz[j]);

                t_xxzz_yyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_xx[j] * fl3_fx + 0.125 * fl3_fx * pa_zz[j] + 0.5 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_yy[j] + 0.25 * pa_xxzz[j] * fl2_fx + pa_xxz[j] * fl2_fx * pb_z[j] + 0.75 * pa_xx[j] * fl2_fx * pb_yy[j] + 0.125 * fl3_fx * pb_zz[j] + 0.25 * pa_xx[j] * fl2_fx * pb_zz[j] + 0.25 * fl2_fx * pa_zz[j] * pb_yy[j] + 0.25 * fl2_fx * pa_zz[j] * pb_zz[j] + fl2_fx * pa_z[j] * pb_yyz[j] + 0.5 * pa_xxzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xxzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xxz[j] * fl1_fx * pb_yyz[j] + 0.25 * fl2_fx * pb_yyzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_yyzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yyzz[j] + pa_xxzz[j] * pb_yyzz[j]);

                t_xxzz_yzzz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] * pb_y[j] + 1.125 * fl3_fx * pb_yz[j] + 1.5 * pa_xxz[j] * fl2_fx * pb_y[j] + 2.25 * pa_xx[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_yz[j] + 1.5 * fl2_fx * pa_z[j] * pb_yzz[j] + 1.5 * pa_xxzz[j] * pb_yz[j] * fl1_fx + 3.0 * pa_xxz[j] * fl1_fx * pb_yzz[j] + 0.25 * fl2_fx * pb_yzzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_yzzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yzzz[j] + pa_xxzz[j] * pb_yzzz[j]);

                t_xxzz_zzzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_xx[j] * fl3_fx + 0.375 * fl3_fx * pa_zz[j] + 3.0 * fl3_fx * pa_z[j] * pb_z[j] + 2.25 * fl3_fx * pb_zz[j] + 0.75 * pa_xxzz[j] * fl2_fx + 6.0 * pa_xxz[j] * fl2_fx * pb_z[j] + 4.5 * pa_xx[j] * fl2_fx * pb_zz[j] + 1.5 * fl2_fx * pa_zz[j] * pb_zz[j] + 2.0 * fl2_fx * pa_z[j] * pb_zzz[j] + 3.0 * pa_xxzz[j] * pb_zz[j] * fl1_fx + 4.0 * pa_xxz[j] * fl1_fx * pb_zzz[j] + 0.25 * fl2_fx * pb_zzzz[j] + 0.5 * pa_xx[j] * fl1_fx * pb_zzzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_zzzz[j] + pa_xxzz[j] * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_90_100(      CMemBlock2D<double>& primBuffer,
                            const CMemBlock2D<double>& auxBuffer,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CMemBlock2D<double>& pbDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
    {
        // Batch of Integrals (90,100)

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

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_yy = paDistances.data(34 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xyyy = paDistances.data(34 * idx + 25);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (90,100)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, \
                                     pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xyyy_xxxx, t_xyyy_xxxy, t_xyyy_xxxz, \
                                     t_xyyy_xxyy, t_xyyy_xxyz, t_xyyy_xxzz, t_xyyy_xyyy, t_xyyy_xyyz, t_xyyy_xyzz, \
                                     t_xyyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyyy_xxxx[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 4.5 * fl3_fx * pa_y[j] * pb_x[j] + 0.75 * pa_xyyy[j] * fl2_fx + 3.0 * fl2_fx * pa_yyy[j] * pb_x[j] + 4.5 * pa_xy[j] * fl2_fx * pb_xx[j] + 3.0 * fl2_fx * pa_y[j] * pb_xxx[j] + 3.0 * pa_xyyy[j] * pb_xx[j] * fl1_fx + 2.0 * fl1_fx * pa_yyy[j] * pb_xxx[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxxx[j] + pa_xyyy[j] * pb_xxxx[j]);

                t_xyyy_xxxy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * fl3_fx * pa_yy[j] + 1.125 * pa_x[j] * fl3_fx * pb_x[j] + 1.125 * fl3_fx * pa_y[j] * pb_y[j] + 1.125 * fl3_fx * pb_xx[j] + 2.25 * pa_xyy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yyy[j] * pb_y[j] + 2.25 * fl2_fx * pa_yy[j] * pb_xx[j] + 2.25 * pa_xy[j] * fl2_fx * pb_xy[j] + 0.75 * pa_x[j] * fl2_fx * pb_xxx[j] + 2.25 * fl2_fx * pa_y[j] * pb_xxy[j] + 1.5 * pa_xyyy[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xyy[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_yyy[j] * pb_xxy[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxxy[j] + pa_xyyy[j] * pb_xxxy[j]);

                t_xyyy_xxxz[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_y[j] * pb_z[j] + 0.75 * fl2_fx * pa_yyy[j] * pb_z[j] + 2.25 * pa_xy[j] * fl2_fx * pb_xz[j] + 2.25 * fl2_fx * pa_y[j] * pb_xxz[j] + 1.5 * pa_xyyy[j] * pb_xz[j] * fl1_fx + 1.5 * fl1_fx * pa_yyy[j] * pb_xxz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxxz[j] + pa_xyyy[j] * pb_xxxz[j]);

                t_xyyy_xxyy[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 2.25 * fl3_fx * pa_y[j] * pb_x[j] + 0.75 * pa_x[j] * fl3_fx * pb_y[j] + 1.5 * fl3_fx * pb_xy[j] + 0.25 * pa_xyyy[j] * fl2_fx + 1.5 * pa_xyy[j] * fl2_fx * pb_y[j] + 2.25 * pa_xy[j] * fl2_fx * pb_xx[j] + 0.5 * fl2_fx * pa_yyy[j] * pb_x[j] + 3.0 * fl2_fx * pa_yy[j] * pb_xy[j] + 0.75 * pa_xy[j] * fl2_fx * pb_yy[j] + 1.5 * pa_x[j] * fl2_fx * pb_xxy[j] + 1.5 * fl2_fx * pa_y[j] * pb_xyy[j] + 0.5 * pa_xyyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyyy[j] * fl1_fx * pb_yy[j] + 3.0 * pa_xyy[j] * fl1_fx * pb_xxy[j] + fl1_fx * pa_yyy[j] * pb_xyy[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxyy[j] + pa_xyyy[j] * pb_xxyy[j]);

                t_xyyy_xxyz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx * pb_z[j] + 0.75 * fl3_fx * pb_xz[j] + 0.75 * pa_xyy[j] * fl2_fx * pb_z[j] + 1.5 * fl2_fx * pa_yy[j] * pb_xz[j] + 0.75 * pa_xy[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xxz[j] + 1.5 * fl2_fx * pa_y[j] * pb_xyz[j] + 0.5 * pa_xyyy[j] * fl1_fx * pb_yz[j] + 1.5 * pa_xyy[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_yyy[j] * pb_xyz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxyz[j] + pa_xyyy[j] * pb_xxyz[j]);

                t_xyyy_xxzz[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 0.75 * fl3_fx * pa_y[j] * pb_x[j] + 0.25 * pa_xyyy[j] * fl2_fx + 0.5 * fl2_fx * pa_yyy[j] * pb_x[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xx[j] + 0.75 * pa_xy[j] * fl2_fx * pb_zz[j] + 1.5 * fl2_fx * pa_y[j] * pb_xzz[j] + 0.5 * pa_xyyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyyy[j] * fl1_fx * pb_zz[j] + fl1_fx * pa_yyy[j] * pb_xzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xxzz[j] + pa_xyyy[j] * pb_xxzz[j]);

                t_xyyy_xyyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_x[j] * fl3_fx * pb_x[j] + 1.125 * fl3_fx * pa_yy[j] + 3.375 * fl3_fx * pa_y[j] * pb_y[j] + 1.125 * fl3_fx * pb_yy[j] + 2.25 * pa_xyy[j] * fl2_fx * pb_x[j] + 6.75 * pa_xy[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_yyy[j] * pb_y[j] + 2.25 * fl2_fx * pa_yy[j] * pb_yy[j] + 2.25 * pa_x[j] * fl2_fx * pb_xyy[j] + 0.75 * fl2_fx * pa_y[j] * pb_yyy[j] + 1.5 * pa_xyyy[j] * pb_xy[j] * fl1_fx + 4.5 * pa_xyy[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_yyy[j] * pb_yyy[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xyyy[j] + pa_xyyy[j] * pb_xyyy[j]);

                t_xyyy_xyyz[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_y[j] * pb_z[j] + 0.75 * fl3_fx * pb_yz[j] + 2.25 * pa_xy[j] * fl2_fx * pb_xz[j] + 0.25 * fl2_fx * pa_yyy[j] * pb_z[j] + 1.5 * fl2_fx * pa_yy[j] * pb_yz[j] + 1.5 * pa_x[j] * fl2_fx * pb_xyz[j] + 0.75 * fl2_fx * pa_y[j] * pb_yyz[j] + 0.5 * pa_xyyy[j] * pb_xz[j] * fl1_fx + 3.0 * pa_xyy[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_yyy[j] * pb_yyz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xyyz[j] + pa_xyyy[j] * pb_xyyz[j]);

                t_xyyy_xyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_yy[j] + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pb_zz[j] + 0.75 * pa_xyy[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pa_yyy[j] * pb_y[j] + 0.75 * fl2_fx * pa_yy[j] * pb_zz[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xy[j] + 0.75 * pa_x[j] * fl2_fx * pb_xzz[j] + 0.75 * fl2_fx * pa_y[j] * pb_yzz[j] + 0.5 * pa_xyyy[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xyy[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_yyy[j] * pb_yzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xyzz[j] + pa_xyyy[j] * pb_xyzz[j]);

                t_xyyy_xzzz[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_y[j] * pb_z[j] + 0.75 * fl2_fx * pa_yyy[j] * pb_z[j] + 2.25 * pa_xy[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_y[j] * pb_zzz[j] + 1.5 * pa_xyyy[j] * pb_xz[j] * fl1_fx + 0.5 * fl1_fx * pa_yyy[j] * pb_zzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_xzzz[j] + pa_xyyy[j] * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_100_110(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (100,110)

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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xyyy = paDistances.data(34 * idx + 25);

            auto pa_xyyz = paDistances.data(34 * idx + 26);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (100,110)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_xyyz, pa_xyz, pa_xz, pa_y, pa_yy, pa_yyz, \
                                     pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xyyy_yyyy, \
                                     t_xyyy_yyyz, t_xyyy_yyzz, t_xyyy_yzzz, t_xyyy_zzzz, t_xyyz_xxxx, t_xyyz_xxxy, \
                                     t_xyyz_xxxz, t_xyyz_xxyy, t_xyyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyyy_yyyy[j] = fl_s_0_0 * (5.625 * pa_xy[j] * fl3_fx + 7.5 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * pa_xyyy[j] * fl2_fx + 9.0 * pa_xyy[j] * fl2_fx * pb_y[j] + 13.5 * pa_xy[j] * fl2_fx * pb_yy[j] + 3.0 * pa_x[j] * fl2_fx * pb_yyy[j] + 3.0 * pa_xyyy[j] * pb_yy[j] * fl1_fx + 6.0 * pa_xyy[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_xy[j] * fl1_fx * pb_yyyy[j] + pa_xyyy[j] * pb_yyyy[j]);

                t_xyyy_yyyz[j] = fl_s_0_0 * (1.875 * pa_x[j] * fl3_fx * pb_z[j] + 2.25 * pa_xyy[j] * fl2_fx * pb_z[j] + 6.75 * pa_xy[j] * fl2_fx * pb_yz[j] + 2.25 * pa_x[j] * fl2_fx * pb_yyz[j] + 1.5 * pa_xyyy[j] * pb_yz[j] * fl1_fx + 4.5 * pa_xyy[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_yyyz[j] + pa_xyyy[j] * pb_yyyz[j]);

                t_xyyy_yyzz[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 0.75 * pa_x[j] * fl3_fx * pb_y[j] + 0.25 * pa_xyyy[j] * fl2_fx + 1.5 * pa_xyy[j] * fl2_fx * pb_y[j] + 2.25 * pa_xy[j] * fl2_fx * pb_zz[j] + 0.75 * pa_xy[j] * fl2_fx * pb_yy[j] + 1.5 * pa_x[j] * fl2_fx * pb_yzz[j] + 0.5 * pa_xyyy[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xyyy[j] * fl1_fx * pb_zz[j] + 3.0 * pa_xyy[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_yyzz[j] + pa_xyyy[j] * pb_yyzz[j]);

                t_xyyy_yzzz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx * pb_z[j] + 2.25 * pa_xyy[j] * fl2_fx * pb_z[j] + 2.25 * pa_xy[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_zzz[j] + 1.5 * pa_xyyy[j] * pb_yz[j] * fl1_fx + 1.5 * pa_xyy[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_xy[j] * fl1_fx * pb_yzzz[j] + pa_xyyy[j] * pb_yzzz[j]);

                t_xyyy_zzzz[j] = fl_s_0_0 * (1.125 * pa_xy[j] * fl3_fx + 0.75 * pa_xyyy[j] * fl2_fx + 4.5 * pa_xy[j] * fl2_fx * pb_zz[j] + 3.0 * pa_xyyy[j] * pb_zz[j] * fl1_fx + 1.5 * pa_xy[j] * fl1_fx * pb_zzzz[j] + pa_xyyy[j] * pb_zzzz[j]);

                t_xyyz_xxxx[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 1.5 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * pa_xyyz[j] * fl2_fx + 3.0 * fl2_fx * pa_yyz[j] * pb_x[j] + 1.5 * pa_xz[j] * fl2_fx * pb_xx[j] + fl2_fx * pa_z[j] * pb_xxx[j] + 3.0 * pa_xyyz[j] * pb_xx[j] * fl1_fx + 2.0 * fl1_fx * pa_yyz[j] * pb_xxx[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xxxx[j] + pa_xyyz[j] * pb_xxxx[j]);

                t_xyyz_xxxy[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_yz[j] + 0.375 * fl3_fx * pa_z[j] * pb_y[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yyz[j] * pb_y[j] + 1.5 * fl2_fx * pa_yz[j] * pb_xx[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_z[j] * pb_xxy[j] + 1.5 * pa_xyyz[j] * pb_xy[j] * fl1_fx + pa_xyz[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_yyz[j] * pb_xxy[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xxxy[j] + pa_xyyz[j] * pb_xxxy[j]);

                t_xyyz_xxxz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_yy[j] + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_xx[j] + 0.75 * pa_xyy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yyz[j] * pb_z[j] + 0.75 * fl2_fx * pa_yy[j] * pb_xx[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xz[j] + 0.25 * pa_x[j] * fl2_fx * pb_xxx[j] + 0.75 * fl2_fx * pa_z[j] * pb_xxz[j] + 1.5 * pa_xyyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xyy[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_yyz[j] * pb_xxz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xxxz[j] + pa_xyyz[j] * pb_xxxz[j]);

                t_xyyz_xxyy[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 0.75 * fl3_fx * pa_z[j] * pb_x[j] + 0.25 * pa_xyyz[j] * fl2_fx + pa_xyz[j] * fl2_fx * pb_y[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xx[j] + 0.5 * fl2_fx * pa_yyz[j] * pb_x[j] + 2.0 * fl2_fx * pa_yz[j] * pb_xy[j] + 0.25 * pa_xz[j] * fl2_fx * pb_yy[j] + 0.5 * fl2_fx * pa_z[j] * pb_xyy[j] + 0.5 * pa_xyyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyyz[j] * fl1_fx * pb_yy[j] + 2.0 * pa_xyz[j] * fl1_fx * pb_xxy[j] + fl1_fx * pa_yyz[j] * pb_xyy[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xxyy[j] + pa_xyyz[j] * pb_xxyy[j]);

                t_xyyz_xxyz[j] = fl_s_0_0 * (0.25 * pa_xy[j] * fl3_fx + 0.5 * fl3_fx * pa_y[j] * pb_x[j] + 0.125 * pa_x[j] * fl3_fx * pb_y[j] + 0.25 * fl3_fx * pb_xy[j] + 0.25 * pa_xyy[j] * fl2_fx * pb_y[j] + 0.5 * pa_xyz[j] * fl2_fx * pb_z[j] + 0.5 * pa_xy[j] * fl2_fx * pb_xx[j] + 0.5 * fl2_fx * pa_yy[j] * pb_xy[j] + fl2_fx * pa_yz[j] * pb_xz[j] + 0.25 * pa_xz[j] * fl2_fx * pb_yz[j] + 0.25 * pa_x[j] * fl2_fx * pb_xxy[j] + 0.5 * fl2_fx * pa_z[j] * pb_xyz[j] + 0.5 * pa_xyyz[j] * fl1_fx * pb_yz[j] + 0.5 * pa_xyy[j] * fl1_fx * pb_xxy[j] + pa_xyz[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_yyz[j] * pb_xyz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xxyz[j] + pa_xyyz[j] * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_110_120(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (110,120)

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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyy = paDistances.data(34 * idx + 12);

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xyyz = paDistances.data(34 * idx + 26);

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

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // Batch of Integrals (110,120)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, \
                                     pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, \
                                     pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, \
                                     pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xyyz_xxzz, t_xyyz_xyyy, \
                                     t_xyyz_xyyz, t_xyyz_xyzz, t_xyyz_xzzz, t_xyyz_yyyy, t_xyyz_yyyz, t_xyyz_yyzz, \
                                     t_xyyz_yzzz, t_xyyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyyz_xxzz[j] = fl_s_0_0 * (0.125 * pa_xz[j] * fl3_fx + 0.25 * pa_x[j] * fl3_fx * pb_z[j] + 0.25 * fl3_fx * pa_z[j] * pb_x[j] + 0.5 * fl3_fx * pb_xz[j] + 0.25 * pa_xyyz[j] * fl2_fx + 0.5 * pa_xyy[j] * fl2_fx * pb_z[j] + 0.5 * fl2_fx * pa_yyz[j] * pb_x[j] + fl2_fx * pa_yy[j] * pb_xz[j] + 0.25 * pa_xz[j] * fl2_fx * pb_xx[j] + 0.25 * pa_xz[j] * fl2_fx * pb_zz[j] + 0.5 * pa_x[j] * fl2_fx * pb_xxz[j] + 0.5 * fl2_fx * pa_z[j] * pb_xzz[j] + 0.5 * pa_xyyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyyz[j] * fl1_fx * pb_zz[j] + pa_xyy[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_yyz[j] * pb_xzz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xxzz[j] + pa_xyyz[j] * pb_xxzz[j]);

                t_xyyz_xyyy[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_yz[j] + 1.125 * fl3_fx * pa_z[j] * pb_y[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_x[j] + 2.25 * pa_xz[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_yyz[j] * pb_y[j] + 1.5 * fl2_fx * pa_yz[j] * pb_yy[j] + 0.25 * fl2_fx * pa_z[j] * pb_yyy[j] + 1.5 * pa_xyyz[j] * pb_xy[j] * fl1_fx + 3.0 * pa_xyz[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_yyz[j] * pb_yyy[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xyyy[j] + pa_xyyz[j] * pb_xyyy[j]);

                t_xyyz_xyyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.125 * fl3_fx * pa_yy[j] + 0.5 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.125 * fl3_fx * pb_yy[j] + 0.25 * pa_xyy[j] * fl2_fx * pb_x[j] + pa_xy[j] * fl2_fx * pb_xy[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xz[j] + 0.25 * fl2_fx * pa_yyz[j] * pb_z[j] + 0.25 * fl2_fx * pa_yy[j] * pb_yy[j] + fl2_fx * pa_yz[j] * pb_yz[j] + 0.25 * pa_x[j] * fl2_fx * pb_xyy[j] + 0.25 * fl2_fx * pa_z[j] * pb_yyz[j] + 0.5 * pa_xyyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_xyy[j] * fl1_fx * pb_xyy[j] + 2.0 * pa_xyz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_yyz[j] * pb_yyz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xyyz[j] + pa_xyyz[j] * pb_xyyz[j]);

                t_xyyz_xyzz[j] = fl_s_0_0 * (0.25 * fl3_fx * pa_yz[j] + 0.5 * fl3_fx * pa_y[j] * pb_z[j] + 0.125 * fl3_fx * pa_z[j] * pb_y[j] + 0.25 * fl3_fx * pb_yz[j] + 0.5 * pa_xyz[j] * fl2_fx * pb_x[j] + pa_xy[j] * fl2_fx * pb_xz[j] + 0.25 * fl2_fx * pa_yyz[j] * pb_y[j] + 0.5 * fl2_fx * pa_yy[j] * pb_yz[j] + 0.5 * fl2_fx * pa_yz[j] * pb_zz[j] + 0.25 * pa_xz[j] * fl2_fx * pb_xy[j] + 0.5 * pa_x[j] * fl2_fx * pb_xyz[j] + 0.25 * fl2_fx * pa_z[j] * pb_yzz[j] + 0.5 * pa_xyyz[j] * pb_xy[j] * fl1_fx + pa_xyy[j] * fl1_fx * pb_xyz[j] + pa_xyz[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_yyz[j] * pb_yzz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xyzz[j] + pa_xyyz[j] * pb_xyzz[j]);

                t_xyyz_xzzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_yy[j] + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_zz[j] + 0.75 * pa_xyy[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yyz[j] * pb_z[j] + 0.75 * fl2_fx * pa_yy[j] * pb_zz[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xzz[j] + 0.25 * fl2_fx * pa_z[j] * pb_zzz[j] + 1.5 * pa_xyyz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xyy[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_yyz[j] * pb_zzz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_xzzz[j] + pa_xyyz[j] * pb_xzzz[j]);

                t_xyyz_yyyy[j] = fl_s_0_0 * (1.875 * pa_xz[j] * fl3_fx + 0.75 * pa_xyyz[j] * fl2_fx + 6.0 * pa_xyz[j] * fl2_fx * pb_y[j] + 4.5 * pa_xz[j] * fl2_fx * pb_yy[j] + 3.0 * pa_xyyz[j] * pb_yy[j] * fl1_fx + 4.0 * pa_xyz[j] * fl1_fx * pb_yyy[j] + 0.5 * pa_xz[j] * fl1_fx * pb_yyyy[j] + pa_xyyz[j] * pb_yyyy[j]);

                t_xyyz_yyyz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl3_fx + 1.125 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * pa_xyy[j] * fl2_fx * pb_y[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_z[j] + 1.5 * pa_xy[j] * fl2_fx * pb_yy[j] + 2.25 * pa_xz[j] * fl2_fx * pb_yz[j] + 0.25 * pa_x[j] * fl2_fx * pb_yyy[j] + 1.5 * pa_xyyz[j] * pb_yz[j] * fl1_fx + 0.5 * pa_xyy[j] * fl1_fx * pb_yyy[j] + 3.0 * pa_xyz[j] * fl1_fx * pb_yyz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_yyyz[j] + pa_xyyz[j] * pb_yyyz[j]);

                t_xyyz_yyzz[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 0.75 * pa_x[j] * fl3_fx * pb_z[j] + 0.25 * pa_xyyz[j] * fl2_fx + 0.5 * pa_xyy[j] * fl2_fx * pb_z[j] + pa_xyz[j] * fl2_fx * pb_y[j] + 2.0 * pa_xy[j] * fl2_fx * pb_yz[j] + 0.75 * pa_xz[j] * fl2_fx * pb_zz[j] + 0.25 * pa_xz[j] * fl2_fx * pb_yy[j] + 0.5 * pa_x[j] * fl2_fx * pb_yyz[j] + 0.5 * pa_xyyz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xyyz[j] * fl1_fx * pb_zz[j] + pa_xyy[j] * fl1_fx * pb_yyz[j] + 2.0 * pa_xyz[j] * fl1_fx * pb_yzz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_yyzz[j] + pa_xyyz[j] * pb_yyzz[j]);

                t_xyyz_yzzz[j] = fl_s_0_0 * (0.75 * pa_xy[j] * fl3_fx + 0.375 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * pa_xyy[j] * fl2_fx * pb_y[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_z[j] + 1.5 * pa_xy[j] * fl2_fx * pb_zz[j] + 0.75 * pa_xz[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_yzz[j] + 1.5 * pa_xyyz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_xyy[j] * fl1_fx * pb_yzz[j] + pa_xyz[j] * fl1_fx * pb_zzz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_yzzz[j] + pa_xyyz[j] * pb_yzzz[j]);

                t_xyyz_zzzz[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 1.5 * pa_x[j] * fl3_fx * pb_z[j] + 0.75 * pa_xyyz[j] * fl2_fx + 3.0 * pa_xyy[j] * fl2_fx * pb_z[j] + 1.5 * pa_xz[j] * fl2_fx * pb_zz[j] + pa_x[j] * fl2_fx * pb_zzz[j] + 3.0 * pa_xyyz[j] * pb_zz[j] * fl1_fx + 2.0 * pa_xyy[j] * fl1_fx * pb_zzz[j] + 0.5 * pa_xz[j] * fl1_fx * pb_zzzz[j] + pa_xyyz[j] * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_120_130(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (120,130)

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

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xyzz = paDistances.data(34 * idx + 27);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (120,130)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xyzz_xxxx, \
                                     t_xyzz_xxxy, t_xyzz_xxxz, t_xyzz_xxyy, t_xyzz_xxyz, t_xyzz_xxzz, t_xyzz_xyyy, \
                                     t_xyzz_xyyz, t_xyzz_xyzz, t_xyzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyzz_xxxx[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 1.5 * fl3_fx * pa_y[j] * pb_x[j] + 0.75 * pa_xyzz[j] * fl2_fx + 3.0 * fl2_fx * pa_yzz[j] * pb_x[j] + 1.5 * pa_xy[j] * fl2_fx * pb_xx[j] + fl2_fx * pa_y[j] * pb_xxx[j] + 3.0 * pa_xyzz[j] * pb_xx[j] * fl1_fx + 2.0 * fl1_fx * pa_yzz[j] * pb_xxx[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xxxx[j] + pa_xyzz[j] * pb_xxxx[j]);

                t_xyzz_xxxy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_zz[j] + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pb_xx[j] + 0.75 * pa_xzz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yzz[j] * pb_y[j] + 0.75 * fl2_fx * pa_zz[j] * pb_xx[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xy[j] + 0.25 * pa_x[j] * fl2_fx * pb_xxx[j] + 0.75 * fl2_fx * pa_y[j] * pb_xxy[j] + 1.5 * pa_xyzz[j] * pb_xy[j] * fl1_fx + 0.5 * pa_xzz[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_yzz[j] * pb_xxy[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xxxy[j] + pa_xyzz[j] * pb_xxxy[j]);

                t_xyzz_xxxz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_yz[j] + 0.375 * fl3_fx * pa_y[j] * pb_z[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yzz[j] * pb_z[j] + 1.5 * fl2_fx * pa_yz[j] * pb_xx[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_y[j] * pb_xxz[j] + 1.5 * pa_xyzz[j] * pb_xz[j] * fl1_fx + pa_xyz[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_yzz[j] * pb_xxz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xxxz[j] + pa_xyzz[j] * pb_xxxz[j]);

                t_xyzz_xxyy[j] = fl_s_0_0 * (0.125 * pa_xy[j] * fl3_fx + 0.25 * pa_x[j] * fl3_fx * pb_y[j] + 0.25 * fl3_fx * pa_y[j] * pb_x[j] + 0.5 * fl3_fx * pb_xy[j] + 0.25 * pa_xyzz[j] * fl2_fx + 0.5 * pa_xzz[j] * fl2_fx * pb_y[j] + 0.5 * fl2_fx * pa_yzz[j] * pb_x[j] + fl2_fx * pa_zz[j] * pb_xy[j] + 0.25 * pa_xy[j] * fl2_fx * pb_xx[j] + 0.25 * pa_xy[j] * fl2_fx * pb_yy[j] + 0.5 * pa_x[j] * fl2_fx * pb_xxy[j] + 0.5 * fl2_fx * pa_y[j] * pb_xyy[j] + 0.5 * pa_xyzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyzz[j] * fl1_fx * pb_yy[j] + pa_xzz[j] * fl1_fx * pb_xxy[j] + fl1_fx * pa_yzz[j] * pb_xyy[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xxyy[j] + pa_xyzz[j] * pb_xxyy[j]);

                t_xyzz_xxyz[j] = fl_s_0_0 * (0.25 * pa_xz[j] * fl3_fx + 0.5 * fl3_fx * pa_z[j] * pb_x[j] + 0.125 * pa_x[j] * fl3_fx * pb_z[j] + 0.25 * fl3_fx * pb_xz[j] + 0.5 * pa_xyz[j] * fl2_fx * pb_y[j] + 0.25 * pa_xzz[j] * fl2_fx * pb_z[j] + 0.5 * pa_xz[j] * fl2_fx * pb_xx[j] + fl2_fx * pa_yz[j] * pb_xy[j] + 0.5 * fl2_fx * pa_zz[j] * pb_xz[j] + 0.25 * pa_xy[j] * fl2_fx * pb_yz[j] + 0.25 * pa_x[j] * fl2_fx * pb_xxz[j] + 0.5 * fl2_fx * pa_y[j] * pb_xyz[j] + 0.5 * pa_xyzz[j] * fl1_fx * pb_yz[j] + pa_xyz[j] * fl1_fx * pb_xxy[j] + 0.5 * pa_xzz[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_yzz[j] * pb_xyz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xxyz[j] + pa_xyzz[j] * pb_xxyz[j]);

                t_xyzz_xxzz[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 0.75 * fl3_fx * pa_y[j] * pb_x[j] + 0.25 * pa_xyzz[j] * fl2_fx + pa_xyz[j] * fl2_fx * pb_z[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xx[j] + 0.5 * fl2_fx * pa_yzz[j] * pb_x[j] + 2.0 * fl2_fx * pa_yz[j] * pb_xz[j] + 0.25 * pa_xy[j] * fl2_fx * pb_zz[j] + 0.5 * fl2_fx * pa_y[j] * pb_xzz[j] + 0.5 * pa_xyzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xyzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xyz[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_yzz[j] * pb_xzz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xxzz[j] + pa_xyzz[j] * pb_xxzz[j]);

                t_xyzz_xyyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_zz[j] + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_y[j] * pb_y[j] + 0.375 * fl3_fx * pb_yy[j] + 0.75 * pa_xzz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_yzz[j] * pb_y[j] + 0.75 * fl2_fx * pa_zz[j] * pb_yy[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xy[j] + 0.75 * pa_x[j] * fl2_fx * pb_xyy[j] + 0.25 * fl2_fx * pa_y[j] * pb_yyy[j] + 1.5 * pa_xyzz[j] * pb_xy[j] * fl1_fx + 1.5 * pa_xzz[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_yzz[j] * pb_yyy[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xyyy[j] + pa_xyzz[j] * pb_xyyy[j]);

                t_xyzz_xyyz[j] = fl_s_0_0 * (0.25 * fl3_fx * pa_yz[j] + 0.5 * fl3_fx * pa_z[j] * pb_y[j] + 0.125 * fl3_fx * pa_y[j] * pb_z[j] + 0.25 * fl3_fx * pb_yz[j] + 0.5 * pa_xyz[j] * fl2_fx * pb_x[j] + pa_xz[j] * fl2_fx * pb_xy[j] + 0.25 * fl2_fx * pa_yzz[j] * pb_z[j] + 0.5 * fl2_fx * pa_yz[j] * pb_yy[j] + 0.5 * fl2_fx * pa_zz[j] * pb_yz[j] + 0.25 * pa_xy[j] * fl2_fx * pb_xz[j] + 0.5 * pa_x[j] * fl2_fx * pb_xyz[j] + 0.25 * fl2_fx * pa_y[j] * pb_yyz[j] + 0.5 * pa_xyzz[j] * pb_xz[j] * fl1_fx + pa_xyz[j] * fl1_fx * pb_xyy[j] + pa_xzz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_yzz[j] * pb_yyz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xyyz[j] + pa_xyzz[j] * pb_xyyz[j]);

                t_xyzz_xyzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_y[j] * pb_y[j] + 0.125 * fl3_fx * pa_zz[j] + 0.5 * fl3_fx * pa_z[j] * pb_z[j] + 0.125 * fl3_fx * pb_zz[j] + 0.75 * pa_xy[j] * fl2_fx * pb_xy[j] + 0.25 * pa_xzz[j] * fl2_fx * pb_x[j] + pa_xz[j] * fl2_fx * pb_xz[j] + 0.25 * fl2_fx * pa_yzz[j] * pb_y[j] + fl2_fx * pa_yz[j] * pb_yz[j] + 0.25 * fl2_fx * pa_zz[j] * pb_zz[j] + 0.25 * pa_x[j] * fl2_fx * pb_xzz[j] + 0.25 * fl2_fx * pa_y[j] * pb_yzz[j] + 0.5 * pa_xyzz[j] * pb_xy[j] * fl1_fx + 2.0 * pa_xyz[j] * fl1_fx * pb_xyz[j] + 0.5 * pa_xzz[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_yzz[j] * pb_yzz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xyzz[j] + pa_xyzz[j] * pb_xyzz[j]);

                t_xyzz_xzzz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_yz[j] + 1.125 * fl3_fx * pa_y[j] * pb_z[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_x[j] + 2.25 * pa_xy[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_yzz[j] * pb_z[j] + 1.5 * fl2_fx * pa_yz[j] * pb_zz[j] + 0.25 * fl2_fx * pa_y[j] * pb_zzz[j] + 1.5 * pa_xyzz[j] * pb_xz[j] * fl1_fx + 3.0 * pa_xyz[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_yzz[j] * pb_zzz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_xzzz[j] + pa_xyzz[j] * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_130_140(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (130,140)

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

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xyzz = paDistances.data(34 * idx + 27);

            auto pa_xzzz = paDistances.data(34 * idx + 28);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (130,140)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, \
                                     pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, \
                                     pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, \
                                     pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xyzz_yyyy, \
                                     t_xyzz_yyyz, t_xyzz_yyzz, t_xyzz_yzzz, t_xyzz_zzzz, t_xzzz_xxxx, t_xzzz_xxxy, \
                                     t_xzzz_xxxz, t_xzzz_xxyy, t_xzzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xyzz_yyyy[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 1.5 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * pa_xyzz[j] * fl2_fx + 3.0 * pa_xzz[j] * fl2_fx * pb_y[j] + 1.5 * pa_xy[j] * fl2_fx * pb_yy[j] + pa_x[j] * fl2_fx * pb_yyy[j] + 3.0 * pa_xyzz[j] * pb_yy[j] * fl1_fx + 2.0 * pa_xzz[j] * fl1_fx * pb_yyy[j] + 0.5 * pa_xy[j] * fl1_fx * pb_yyyy[j] + pa_xyzz[j] * pb_yyyy[j]);

                t_xyzz_yyyz[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl3_fx + 0.375 * pa_x[j] * fl3_fx * pb_z[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_y[j] + 0.75 * pa_xzz[j] * fl2_fx * pb_z[j] + 1.5 * pa_xz[j] * fl2_fx * pb_yy[j] + 0.75 * pa_xy[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_yyz[j] + 1.5 * pa_xyzz[j] * pb_yz[j] * fl1_fx + pa_xyz[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_xzz[j] * fl1_fx * pb_yyz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_yyyz[j] + pa_xyzz[j] * pb_yyyz[j]);

                t_xyzz_yyzz[j] = fl_s_0_0 * (0.375 * pa_xy[j] * fl3_fx + 0.75 * pa_x[j] * fl3_fx * pb_y[j] + 0.25 * pa_xyzz[j] * fl2_fx + pa_xyz[j] * fl2_fx * pb_z[j] + 0.75 * pa_xy[j] * fl2_fx * pb_yy[j] + 0.5 * pa_xzz[j] * fl2_fx * pb_y[j] + 2.0 * pa_xz[j] * fl2_fx * pb_yz[j] + 0.25 * pa_xy[j] * fl2_fx * pb_zz[j] + 0.5 * pa_x[j] * fl2_fx * pb_yzz[j] + 0.5 * pa_xyzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xyzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_xyz[j] * fl1_fx * pb_yyz[j] + pa_xzz[j] * fl1_fx * pb_yzz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_yyzz[j] + pa_xyzz[j] * pb_yyzz[j]);

                t_xyzz_yzzz[j] = fl_s_0_0 * (0.75 * pa_xz[j] * fl3_fx + 1.125 * pa_x[j] * fl3_fx * pb_z[j] + 1.5 * pa_xyz[j] * fl2_fx * pb_y[j] + 2.25 * pa_xy[j] * fl2_fx * pb_yz[j] + 0.75 * pa_xzz[j] * fl2_fx * pb_z[j] + 1.5 * pa_xz[j] * fl2_fx * pb_zz[j] + 0.25 * pa_x[j] * fl2_fx * pb_zzz[j] + 1.5 * pa_xyzz[j] * pb_yz[j] * fl1_fx + 3.0 * pa_xyz[j] * fl1_fx * pb_yzz[j] + 0.5 * pa_xzz[j] * fl1_fx * pb_zzz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_yzzz[j] + pa_xyzz[j] * pb_yzzz[j]);

                t_xyzz_zzzz[j] = fl_s_0_0 * (1.875 * pa_xy[j] * fl3_fx + 0.75 * pa_xyzz[j] * fl2_fx + 6.0 * pa_xyz[j] * fl2_fx * pb_z[j] + 4.5 * pa_xy[j] * fl2_fx * pb_zz[j] + 3.0 * pa_xyzz[j] * pb_zz[j] * fl1_fx + 4.0 * pa_xyz[j] * fl1_fx * pb_zzz[j] + 0.5 * pa_xy[j] * fl1_fx * pb_zzzz[j] + pa_xyzz[j] * pb_zzzz[j]);

                t_xzzz_xxxx[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 4.5 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * pa_xzzz[j] * fl2_fx + 3.0 * fl2_fx * pa_zzz[j] * pb_x[j] + 4.5 * pa_xz[j] * fl2_fx * pb_xx[j] + 3.0 * fl2_fx * pa_z[j] * pb_xxx[j] + 3.0 * pa_xzzz[j] * pb_xx[j] * fl1_fx + 2.0 * fl1_fx * pa_zzz[j] * pb_xxx[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxxx[j] + pa_xzzz[j] * pb_xxxx[j]);

                t_xzzz_xxxy[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_z[j] * pb_y[j] + 0.75 * fl2_fx * pa_zzz[j] * pb_y[j] + 2.25 * pa_xz[j] * fl2_fx * pb_xy[j] + 2.25 * fl2_fx * pa_z[j] * pb_xxy[j] + 1.5 * pa_xzzz[j] * pb_xy[j] * fl1_fx + 1.5 * fl1_fx * pa_zzz[j] * pb_xxy[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxxy[j] + pa_xzzz[j] * pb_xxxy[j]);

                t_xzzz_xxxz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * fl3_fx * pa_zz[j] + 1.125 * pa_x[j] * fl3_fx * pb_x[j] + 1.125 * fl3_fx * pa_z[j] * pb_z[j] + 1.125 * fl3_fx * pb_xx[j] + 2.25 * pa_xzz[j] * fl2_fx * pb_x[j] + 0.75 * fl2_fx * pa_zzz[j] * pb_z[j] + 2.25 * fl2_fx * pa_zz[j] * pb_xx[j] + 2.25 * pa_xz[j] * fl2_fx * pb_xz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xxx[j] + 2.25 * fl2_fx * pa_z[j] * pb_xxz[j] + 1.5 * pa_xzzz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xzz[j] * fl1_fx * pb_xxx[j] + 1.5 * fl1_fx * pa_zzz[j] * pb_xxz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxxz[j] + pa_xzzz[j] * pb_xxxz[j]);

                t_xzzz_xxyy[j] = fl_s_0_0 * (0.375 * pa_xz[j] * fl3_fx + 0.75 * fl3_fx * pa_z[j] * pb_x[j] + 0.25 * pa_xzzz[j] * fl2_fx + 0.5 * fl2_fx * pa_zzz[j] * pb_x[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xx[j] + 0.75 * pa_xz[j] * fl2_fx * pb_yy[j] + 1.5 * fl2_fx * pa_z[j] * pb_xyy[j] + 0.5 * pa_xzzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xzzz[j] * fl1_fx * pb_yy[j] + fl1_fx * pa_zzz[j] * pb_xyy[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxyy[j] + pa_xzzz[j] * pb_xxyy[j]);

                t_xzzz_xxyz[j] = fl_s_0_0 * (0.375 * pa_x[j] * fl3_fx * pb_y[j] + 0.75 * fl3_fx * pb_xy[j] + 0.75 * pa_xzz[j] * fl2_fx * pb_y[j] + 1.5 * fl2_fx * pa_zz[j] * pb_xy[j] + 0.75 * pa_xz[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xxy[j] + 1.5 * fl2_fx * pa_z[j] * pb_xyz[j] + 0.5 * pa_xzzz[j] * fl1_fx * pb_yz[j] + 1.5 * pa_xzz[j] * fl1_fx * pb_xxy[j] + fl1_fx * pa_zzz[j] * pb_xyz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxyz[j] + pa_xzzz[j] * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_140_150(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (140,150)

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

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_xzzz = paDistances.data(34 * idx + 28);

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

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // Batch of Integrals (140,150)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, \
                                     pb_zzz, pb_zzzz, s_0_0, t_xzzz_xxzz, t_xzzz_xyyy, t_xzzz_xyyz, t_xzzz_xyzz, \
                                     t_xzzz_xzzz, t_xzzz_yyyy, t_xzzz_yyyz, t_xzzz_yyzz, t_xzzz_yzzz, t_xzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_xzzz_xxzz[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 2.25 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * pa_x[j] * fl3_fx * pb_z[j] + 1.5 * fl3_fx * pb_xz[j] + 0.25 * pa_xzzz[j] * fl2_fx + 1.5 * pa_xzz[j] * fl2_fx * pb_z[j] + 2.25 * pa_xz[j] * fl2_fx * pb_xx[j] + 0.5 * fl2_fx * pa_zzz[j] * pb_x[j] + 3.0 * fl2_fx * pa_zz[j] * pb_xz[j] + 0.75 * pa_xz[j] * fl2_fx * pb_zz[j] + 1.5 * pa_x[j] * fl2_fx * pb_xxz[j] + 1.5 * fl2_fx * pa_z[j] * pb_xzz[j] + 0.5 * pa_xzzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_xzzz[j] * fl1_fx * pb_zz[j] + 3.0 * pa_xzz[j] * fl1_fx * pb_xxz[j] + fl1_fx * pa_zzz[j] * pb_xzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xxzz[j] + pa_xzzz[j] * pb_xxzz[j]);

                t_xzzz_xyyy[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_z[j] * pb_y[j] + 0.75 * fl2_fx * pa_zzz[j] * pb_y[j] + 2.25 * pa_xz[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_z[j] * pb_yyy[j] + 1.5 * pa_xzzz[j] * pb_xy[j] * fl1_fx + 0.5 * fl1_fx * pa_zzz[j] * pb_yyy[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xyyy[j] + pa_xzzz[j] * pb_xyyy[j]);

                t_xzzz_xyyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_zz[j] + 0.375 * pa_x[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_yy[j] + 0.75 * pa_xzz[j] * fl2_fx * pb_x[j] + 0.25 * fl2_fx * pa_zzz[j] * pb_z[j] + 0.75 * fl2_fx * pa_zz[j] * pb_yy[j] + 0.75 * pa_xz[j] * fl2_fx * pb_xz[j] + 0.75 * pa_x[j] * fl2_fx * pb_xyy[j] + 0.75 * fl2_fx * pa_z[j] * pb_yyz[j] + 0.5 * pa_xzzz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_xzz[j] * fl1_fx * pb_xyy[j] + 0.5 * fl1_fx * pa_zzz[j] * pb_yyz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xyyz[j] + pa_xzzz[j] * pb_xyyz[j]);

                t_xzzz_xyzz[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_z[j] * pb_y[j] + 0.75 * fl3_fx * pb_yz[j] + 2.25 * pa_xz[j] * fl2_fx * pb_xy[j] + 0.25 * fl2_fx * pa_zzz[j] * pb_y[j] + 1.5 * fl2_fx * pa_zz[j] * pb_yz[j] + 1.5 * pa_x[j] * fl2_fx * pb_xyz[j] + 0.75 * fl2_fx * pa_z[j] * pb_yzz[j] + 0.5 * pa_xzzz[j] * pb_xy[j] * fl1_fx + 3.0 * pa_xzz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_zzz[j] * pb_yzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xyzz[j] + pa_xzzz[j] * pb_xyzz[j]);

                t_xzzz_xzzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_x[j] * fl3_fx * pb_x[j] + 1.125 * fl3_fx * pa_zz[j] + 3.375 * fl3_fx * pa_z[j] * pb_z[j] + 1.125 * fl3_fx * pb_zz[j] + 2.25 * pa_xzz[j] * fl2_fx * pb_x[j] + 6.75 * pa_xz[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_zzz[j] * pb_z[j] + 2.25 * fl2_fx * pa_zz[j] * pb_zz[j] + 2.25 * pa_x[j] * fl2_fx * pb_xzz[j] + 0.75 * fl2_fx * pa_z[j] * pb_zzz[j] + 1.5 * pa_xzzz[j] * pb_xz[j] * fl1_fx + 4.5 * pa_xzz[j] * fl1_fx * pb_xzz[j] + 0.5 * fl1_fx * pa_zzz[j] * pb_zzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_xzzz[j] + pa_xzzz[j] * pb_xzzz[j]);

                t_xzzz_yyyy[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 0.75 * pa_xzzz[j] * fl2_fx + 4.5 * pa_xz[j] * fl2_fx * pb_yy[j] + 3.0 * pa_xzzz[j] * pb_yy[j] * fl1_fx + 1.5 * pa_xz[j] * fl1_fx * pb_yyyy[j] + pa_xzzz[j] * pb_yyyy[j]);

                t_xzzz_yyyz[j] = fl_s_0_0 * (1.125 * pa_x[j] * fl3_fx * pb_y[j] + 2.25 * pa_xzz[j] * fl2_fx * pb_y[j] + 2.25 * pa_xz[j] * fl2_fx * pb_yz[j] + 0.75 * pa_x[j] * fl2_fx * pb_yyy[j] + 1.5 * pa_xzzz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_xzz[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_xz[j] * fl1_fx * pb_yyyz[j] + pa_xzzz[j] * pb_yyyz[j]);

                t_xzzz_yyzz[j] = fl_s_0_0 * (1.125 * pa_xz[j] * fl3_fx + 0.75 * pa_x[j] * fl3_fx * pb_z[j] + 0.25 * pa_xzzz[j] * fl2_fx + 1.5 * pa_xzz[j] * fl2_fx * pb_z[j] + 2.25 * pa_xz[j] * fl2_fx * pb_yy[j] + 0.75 * pa_xz[j] * fl2_fx * pb_zz[j] + 1.5 * pa_x[j] * fl2_fx * pb_yyz[j] + 0.5 * pa_xzzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_xzzz[j] * fl1_fx * pb_zz[j] + 3.0 * pa_xzz[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_yyzz[j] + pa_xzzz[j] * pb_yyzz[j]);

                t_xzzz_yzzz[j] = fl_s_0_0 * (1.875 * pa_x[j] * fl3_fx * pb_y[j] + 2.25 * pa_xzz[j] * fl2_fx * pb_y[j] + 6.75 * pa_xz[j] * fl2_fx * pb_yz[j] + 2.25 * pa_x[j] * fl2_fx * pb_yzz[j] + 1.5 * pa_xzzz[j] * pb_yz[j] * fl1_fx + 4.5 * pa_xzz[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_yzzz[j] + pa_xzzz[j] * pb_yzzz[j]);

                t_xzzz_zzzz[j] = fl_s_0_0 * (5.625 * pa_xz[j] * fl3_fx + 7.5 * pa_x[j] * fl3_fx * pb_z[j] + 0.75 * pa_xzzz[j] * fl2_fx + 9.0 * pa_xzz[j] * fl2_fx * pb_z[j] + 13.5 * pa_xz[j] * fl2_fx * pb_zz[j] + 3.0 * pa_x[j] * fl2_fx * pb_zzz[j] + 3.0 * pa_xzzz[j] * pb_zz[j] * fl1_fx + 6.0 * pa_xzz[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_xz[j] * fl1_fx * pb_zzzz[j] + pa_xzzz[j] * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_150_160(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (150,160)

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

            auto pa_y = paDistances.data(34 * idx + 1);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(34 * idx + 6);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyy = paDistances.data(34 * idx + 15);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_yyyy = paDistances.data(34 * idx + 29);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (150,160)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_yyyy_xxxx, \
                                     t_yyyy_xxxy, t_yyyy_xxxz, t_yyyy_xxyy, t_yyyy_xxyz, t_yyyy_xxzz, t_yyyy_xyyy, \
                                     t_yyyy_xyyz, t_yyyy_xyzz, t_yyyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyyy_xxxx[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_yy[j] * fl3_fx + 0.75 * pa_yyyy[j] * fl2_fx + 2.25 * fl3_fx * pb_xx[j] + 9.0 * pa_yy[j] * fl2_fx * pb_xx[j] + 3.0 * pa_yyyy[j] * pb_xx[j] * fl1_fx + 0.75 * fl2_fx * pb_xxxx[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xxxx[j] + pa_yyyy[j] * pb_xxxx[j]);

                t_yyyy_xxxy[j] = fl_s_0_0 * (4.5 * pa_y[j] * fl3_fx * pb_x[j] + 3.0 * pa_yyy[j] * fl2_fx * pb_x[j] + 1.125 * fl3_fx * pb_xy[j] + 4.5 * pa_yy[j] * fl2_fx * pb_xy[j] + 3.0 * pa_y[j] * fl2_fx * pb_xxx[j] + 1.5 * pa_yyyy[j] * pb_xy[j] * fl1_fx + 2.0 * pa_yyy[j] * fl1_fx * pb_xxx[j] + 0.75 * fl2_fx * pb_xxxy[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xxxy[j] + pa_yyyy[j] * pb_xxxy[j]);

                t_yyyy_xxxz[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_xz[j] + 4.5 * pa_yy[j] * fl2_fx * pb_xz[j] + 1.5 * pa_yyyy[j] * pb_xz[j] * fl1_fx + 0.75 * fl2_fx * pb_xxxz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xxxz[j] + pa_yyyy[j] * pb_xxxz[j]);

                t_yyyy_xxyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_yy[j] * fl3_fx + 3.0 * pa_y[j] * fl3_fx * pb_y[j] + 1.875 * fl3_fx * pb_xx[j] + 0.25 * pa_yyyy[j] * fl2_fx + 2.0 * pa_yyy[j] * fl2_fx * pb_y[j] + 4.5 * pa_yy[j] * fl2_fx * pb_xx[j] + 0.375 * fl3_fx * pb_yy[j] + 1.5 * pa_yy[j] * fl2_fx * pb_yy[j] + 6.0 * pa_y[j] * fl2_fx * pb_xxy[j] + 0.5 * pa_yyyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyyy[j] * fl1_fx * pb_yy[j] + 4.0 * pa_yyy[j] * fl1_fx * pb_xxy[j] + 0.75 * fl2_fx * pb_xxyy[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xxyy[j] + pa_yyyy[j] * pb_xxyy[j]);

                t_yyyy_xxyz[j] = fl_s_0_0 * (1.5 * pa_y[j] * fl3_fx * pb_z[j] + pa_yyy[j] * fl2_fx * pb_z[j] + 0.375 * fl3_fx * pb_yz[j] + 1.5 * pa_yy[j] * fl2_fx * pb_yz[j] + 3.0 * pa_y[j] * fl2_fx * pb_xxz[j] + 0.5 * pa_yyyy[j] * fl1_fx * pb_yz[j] + 2.0 * pa_yyy[j] * fl1_fx * pb_xxz[j] + 0.75 * fl2_fx * pb_xxyz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xxyz[j] + pa_yyyy[j] * pb_xxyz[j]);

                t_yyyy_xxzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.75 * pa_yy[j] * fl3_fx + 0.25 * pa_yyyy[j] * fl2_fx + 0.375 * fl3_fx * pb_xx[j] + 0.375 * fl3_fx * pb_zz[j] + 1.5 * pa_yy[j] * fl2_fx * pb_xx[j] + 1.5 * pa_yy[j] * fl2_fx * pb_zz[j] + 0.5 * pa_yyyy[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyyy[j] * fl1_fx * pb_zz[j] + 0.75 * fl2_fx * pb_xxzz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xxzz[j] + pa_yyyy[j] * pb_xxzz[j]);

                t_yyyy_xyyy[j] = fl_s_0_0 * (7.5 * pa_y[j] * fl3_fx * pb_x[j] + 5.625 * fl3_fx * pb_xy[j] + 3.0 * pa_yyy[j] * fl2_fx * pb_x[j] + 13.5 * pa_yy[j] * fl2_fx * pb_xy[j] + 9.0 * pa_y[j] * fl2_fx * pb_xyy[j] + 1.5 * pa_yyyy[j] * pb_xy[j] * fl1_fx + 6.0 * pa_yyy[j] * fl1_fx * pb_xyy[j] + 0.75 * fl2_fx * pb_xyyy[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xyyy[j] + pa_yyyy[j] * pb_xyyy[j]);

                t_yyyy_xyyz[j] = fl_s_0_0 * (1.875 * fl3_fx * pb_xz[j] + 4.5 * pa_yy[j] * fl2_fx * pb_xz[j] + 6.0 * pa_y[j] * fl2_fx * pb_xyz[j] + 0.5 * pa_yyyy[j] * pb_xz[j] * fl1_fx + 4.0 * pa_yyy[j] * fl1_fx * pb_xyz[j] + 0.75 * fl2_fx * pb_xyyz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xyyz[j] + pa_yyyy[j] * pb_xyyz[j]);

                t_yyyy_xyzz[j] = fl_s_0_0 * (1.5 * pa_y[j] * fl3_fx * pb_x[j] + pa_yyy[j] * fl2_fx * pb_x[j] + 0.375 * fl3_fx * pb_xy[j] + 1.5 * pa_yy[j] * fl2_fx * pb_xy[j] + 3.0 * pa_y[j] * fl2_fx * pb_xzz[j] + 0.5 * pa_yyyy[j] * pb_xy[j] * fl1_fx + 2.0 * pa_yyy[j] * fl1_fx * pb_xzz[j] + 0.75 * fl2_fx * pb_xyzz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xyzz[j] + pa_yyyy[j] * pb_xyzz[j]);

                t_yyyy_xzzz[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_xz[j] + 4.5 * pa_yy[j] * fl2_fx * pb_xz[j] + 1.5 * pa_yyyy[j] * pb_xz[j] * fl1_fx + 0.75 * fl2_fx * pb_xzzz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_xzzz[j] + pa_yyyy[j] * pb_xzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_160_170(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (160,170)

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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_yyyy = paDistances.data(34 * idx + 29);

            auto pa_yyyz = paDistances.data(34 * idx + 30);

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

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

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

            // Batch of Integrals (160,170)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, \
                                     pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, \
                                     pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, \
                                     pb_zz, pb_zzz, pb_zzzz, s_0_0, t_yyyy_yyyy, t_yyyy_yyyz, t_yyyy_yyzz, \
                                     t_yyyy_yzzz, t_yyyy_zzzz, t_yyyz_xxxx, t_yyyz_xxxy, t_yyyz_xxxz, t_yyyz_xxyy, \
                                     t_yyyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyyy_yyyy[j] = fl_s_0_0 * (6.5625 * fl4_fx + 11.25 * pa_yy[j] * fl3_fx + 30.0 * pa_y[j] * fl3_fx * pb_y[j] + 11.25 * fl3_fx * pb_yy[j] + 0.75 * pa_yyyy[j] * fl2_fx + 12.0 * pa_yyy[j] * fl2_fx * pb_y[j] + 27.0 * pa_yy[j] * fl2_fx * pb_yy[j] + 12.0 * pa_y[j] * fl2_fx * pb_yyy[j] + 3.0 * pa_yyyy[j] * pb_yy[j] * fl1_fx + 8.0 * pa_yyy[j] * fl1_fx * pb_yyy[j] + 0.75 * fl2_fx * pb_yyyy[j] + 3.0 * pa_yy[j] * fl1_fx * pb_yyyy[j] + pa_yyyy[j] * pb_yyyy[j]);

                t_yyyy_yyyz[j] = fl_s_0_0 * (7.5 * pa_y[j] * fl3_fx * pb_z[j] + 5.625 * fl3_fx * pb_yz[j] + 3.0 * pa_yyy[j] * fl2_fx * pb_z[j] + 13.5 * pa_yy[j] * fl2_fx * pb_yz[j] + 9.0 * pa_y[j] * fl2_fx * pb_yyz[j] + 1.5 * pa_yyyy[j] * pb_yz[j] * fl1_fx + 6.0 * pa_yyy[j] * fl1_fx * pb_yyz[j] + 0.75 * fl2_fx * pb_yyyz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_yyyz[j] + pa_yyyy[j] * pb_yyyz[j]);

                t_yyyy_yyzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_yy[j] * fl3_fx + 3.0 * pa_y[j] * fl3_fx * pb_y[j] + 1.875 * fl3_fx * pb_zz[j] + 0.25 * pa_yyyy[j] * fl2_fx + 2.0 * pa_yyy[j] * fl2_fx * pb_y[j] + 4.5 * pa_yy[j] * fl2_fx * pb_zz[j] + 0.375 * fl3_fx * pb_yy[j] + 1.5 * pa_yy[j] * fl2_fx * pb_yy[j] + 6.0 * pa_y[j] * fl2_fx * pb_yzz[j] + 0.5 * pa_yyyy[j] * pb_yy[j] * fl1_fx + 0.5 * pa_yyyy[j] * fl1_fx * pb_zz[j] + 4.0 * pa_yyy[j] * fl1_fx * pb_yzz[j] + 0.75 * fl2_fx * pb_yyzz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_yyzz[j] + pa_yyyy[j] * pb_yyzz[j]);

                t_yyyy_yzzz[j] = fl_s_0_0 * (4.5 * pa_y[j] * fl3_fx * pb_z[j] + 3.0 * pa_yyy[j] * fl2_fx * pb_z[j] + 1.125 * fl3_fx * pb_yz[j] + 4.5 * pa_yy[j] * fl2_fx * pb_yz[j] + 3.0 * pa_y[j] * fl2_fx * pb_zzz[j] + 1.5 * pa_yyyy[j] * pb_yz[j] * fl1_fx + 2.0 * pa_yyy[j] * fl1_fx * pb_zzz[j] + 0.75 * fl2_fx * pb_yzzz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_yzzz[j] + pa_yyyy[j] * pb_yzzz[j]);

                t_yyyy_zzzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_yy[j] * fl3_fx + 0.75 * pa_yyyy[j] * fl2_fx + 2.25 * fl3_fx * pb_zz[j] + 9.0 * pa_yy[j] * fl2_fx * pb_zz[j] + 3.0 * pa_yyyy[j] * pb_zz[j] * fl1_fx + 0.75 * fl2_fx * pb_zzzz[j] + 3.0 * pa_yy[j] * fl1_fx * pb_zzzz[j] + pa_yyyy[j] * pb_zzzz[j]);

                t_yyyz_xxxx[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 0.75 * pa_yyyz[j] * fl2_fx + 4.5 * pa_yz[j] * fl2_fx * pb_xx[j] + 3.0 * pa_yyyz[j] * pb_xx[j] * fl1_fx + 1.5 * pa_yz[j] * fl1_fx * pb_xxxx[j] + pa_yyyz[j] * pb_xxxx[j]);

                t_yyyz_xxxy[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_z[j] * pb_x[j] + 2.25 * pa_yyz[j] * fl2_fx * pb_x[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_z[j] * pb_xxx[j] + 1.5 * pa_yyyz[j] * pb_xy[j] * fl1_fx + 1.5 * pa_yyz[j] * fl1_fx * pb_xxx[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxxy[j] + pa_yyyz[j] * pb_xxxy[j]);

                t_yyyz_xxxz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx * pb_x[j] + 0.75 * pa_yyy[j] * fl2_fx * pb_x[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xz[j] + 0.75 * pa_y[j] * fl2_fx * pb_xxx[j] + 1.5 * pa_yyyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_yyy[j] * fl1_fx * pb_xxx[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxxz[j] + pa_yyyz[j] * pb_xxxz[j]);

                t_yyyz_xxyy[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 0.75 * fl3_fx * pa_z[j] * pb_y[j] + 0.25 * pa_yyyz[j] * fl2_fx + 1.5 * pa_yyz[j] * fl2_fx * pb_y[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xx[j] + 0.75 * pa_yz[j] * fl2_fx * pb_yy[j] + 1.5 * fl2_fx * pa_z[j] * pb_xxy[j] + 0.5 * pa_yyyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyyz[j] * fl1_fx * pb_yy[j] + 3.0 * pa_yyz[j] * fl1_fx * pb_xxy[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxyy[j] + pa_yyyz[j] * pb_xxyy[j]);

                t_yyyz_xxyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.375 * pa_y[j] * fl3_fx * pb_y[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_xx[j] + 0.25 * pa_yyy[j] * fl2_fx * pb_y[j] + 0.75 * pa_yyz[j] * fl2_fx * pb_z[j] + 0.75 * pa_yy[j] * fl2_fx * pb_xx[j] + 0.75 * pa_yz[j] * fl2_fx * pb_yz[j] + 0.75 * pa_y[j] * fl2_fx * pb_xxy[j] + 0.75 * fl2_fx * pa_z[j] * pb_xxz[j] + 0.5 * pa_yyyz[j] * fl1_fx * pb_yz[j] + 0.5 * pa_yyy[j] * fl1_fx * pb_xxy[j] + 1.5 * pa_yyz[j] * fl1_fx * pb_xxz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxyz[j] + pa_yyyz[j] * pb_xxyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_170_180(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (170,180)

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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_yyyz = paDistances.data(34 * idx + 30);

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

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            // Batch of Integrals (170,180)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, \
                                     pb_zzz, pb_zzzz, s_0_0, t_yyyz_xxzz, t_yyyz_xyyy, t_yyyz_xyyz, t_yyyz_xyzz, \
                                     t_yyyz_xzzz, t_yyyz_yyyy, t_yyyz_yyyz, t_yyyz_yyzz, t_yyyz_yzzz, t_yyyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyyz_xxzz[j] = fl_s_0_0 * (0.375 * pa_yz[j] * fl3_fx + 0.75 * pa_y[j] * fl3_fx * pb_z[j] + 0.25 * pa_yyyz[j] * fl2_fx + 0.5 * pa_yyy[j] * fl2_fx * pb_z[j] + 0.75 * pa_yz[j] * fl2_fx * pb_xx[j] + 0.75 * pa_yz[j] * fl2_fx * pb_zz[j] + 1.5 * pa_y[j] * fl2_fx * pb_xxz[j] + 0.5 * pa_yyyz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyyz[j] * fl1_fx * pb_zz[j] + pa_yyy[j] * fl1_fx * pb_xxz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxzz[j] + pa_yyyz[j] * pb_xxzz[j]);

                t_yyyz_xyyy[j] = fl_s_0_0 * (1.875 * fl3_fx * pa_z[j] * pb_x[j] + 2.25 * pa_yyz[j] * fl2_fx * pb_x[j] + 6.75 * pa_yz[j] * fl2_fx * pb_xy[j] + 2.25 * fl2_fx * pa_z[j] * pb_xyy[j] + 1.5 * pa_yyyz[j] * pb_xy[j] * fl1_fx + 4.5 * pa_yyz[j] * fl1_fx * pb_xyy[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xyyy[j] + pa_yyyz[j] * pb_xyyy[j]);

                t_yyyz_xyyz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx * pb_x[j] + 0.75 * fl3_fx * pb_xy[j] + 0.25 * pa_yyy[j] * fl2_fx * pb_x[j] + 1.5 * pa_yy[j] * fl2_fx * pb_xy[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xz[j] + 0.75 * pa_y[j] * fl2_fx * pb_xyy[j] + 1.5 * fl2_fx * pa_z[j] * pb_xyz[j] + 0.5 * pa_yyyz[j] * pb_xz[j] * fl1_fx + 0.5 * pa_yyy[j] * fl1_fx * pb_xyy[j] + 3.0 * pa_yyz[j] * fl1_fx * pb_xyz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xyyz[j] + pa_yyyz[j] * pb_xyyz[j]);

                t_yyyz_xyzz[j] = fl_s_0_0 * (0.375 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * fl3_fx * pb_xz[j] + 0.75 * pa_yyz[j] * fl2_fx * pb_x[j] + 1.5 * pa_yy[j] * fl2_fx * pb_xz[j] + 0.75 * pa_yz[j] * fl2_fx * pb_xy[j] + 1.5 * pa_y[j] * fl2_fx * pb_xyz[j] + 0.75 * fl2_fx * pa_z[j] * pb_xzz[j] + 0.5 * pa_yyyz[j] * pb_xy[j] * fl1_fx + pa_yyy[j] * fl1_fx * pb_xyz[j] + 1.5 * pa_yyz[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xyzz[j] + pa_yyyz[j] * pb_xyzz[j]);

                t_yyyz_xzzz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx * pb_x[j] + 0.75 * pa_yyy[j] * fl2_fx * pb_x[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xz[j] + 2.25 * pa_y[j] * fl2_fx * pb_xzz[j] + 1.5 * pa_yyyz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_yyy[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xzzz[j] + pa_yyyz[j] * pb_xzzz[j]);

                t_yyyz_yyyy[j] = fl_s_0_0 * (5.625 * pa_yz[j] * fl3_fx + 7.5 * fl3_fx * pa_z[j] * pb_y[j] + 0.75 * pa_yyyz[j] * fl2_fx + 9.0 * pa_yyz[j] * fl2_fx * pb_y[j] + 13.5 * pa_yz[j] * fl2_fx * pb_yy[j] + 3.0 * fl2_fx * pa_z[j] * pb_yyy[j] + 3.0 * pa_yyyz[j] * pb_yy[j] * fl1_fx + 6.0 * pa_yyz[j] * fl1_fx * pb_yyy[j] + 1.5 * pa_yz[j] * fl1_fx * pb_yyyy[j] + pa_yyyz[j] * pb_yyyy[j]);

                t_yyyz_yyyz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.125 * pa_yy[j] * fl3_fx + 3.375 * pa_y[j] * fl3_fx * pb_y[j] + 1.875 * fl3_fx * pa_z[j] * pb_z[j] + 1.125 * fl3_fx * pb_yy[j] + 0.75 * pa_yyy[j] * fl2_fx * pb_y[j] + 2.25 * pa_yyz[j] * fl2_fx * pb_z[j] + 2.25 * pa_yy[j] * fl2_fx * pb_yy[j] + 6.75 * pa_yz[j] * fl2_fx * pb_yz[j] + 0.75 * pa_y[j] * fl2_fx * pb_yyy[j] + 2.25 * fl2_fx * pa_z[j] * pb_yyz[j] + 1.5 * pa_yyyz[j] * pb_yz[j] * fl1_fx + 0.5 * pa_yyy[j] * fl1_fx * pb_yyy[j] + 4.5 * pa_yyz[j] * fl1_fx * pb_yyz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_yyyz[j] + pa_yyyz[j] * pb_yyyz[j]);

                t_yyyz_yyzz[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 2.25 * pa_y[j] * fl3_fx * pb_z[j] + 0.75 * fl3_fx * pa_z[j] * pb_y[j] + 1.5 * fl3_fx * pb_yz[j] + 0.25 * pa_yyyz[j] * fl2_fx + 0.5 * pa_yyy[j] * fl2_fx * pb_z[j] + 1.5 * pa_yyz[j] * fl2_fx * pb_y[j] + 3.0 * pa_yy[j] * fl2_fx * pb_yz[j] + 2.25 * pa_yz[j] * fl2_fx * pb_zz[j] + 0.75 * pa_yz[j] * fl2_fx * pb_yy[j] + 1.5 * pa_y[j] * fl2_fx * pb_yyz[j] + 1.5 * fl2_fx * pa_z[j] * pb_yzz[j] + 0.5 * pa_yyyz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_yyyz[j] * fl1_fx * pb_zz[j] + pa_yyy[j] * fl1_fx * pb_yyz[j] + 3.0 * pa_yyz[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_yyzz[j] + pa_yyyz[j] * pb_yyzz[j]);

                t_yyyz_yzzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * pa_yy[j] * fl3_fx + 1.125 * pa_y[j] * fl3_fx * pb_y[j] + 1.125 * fl3_fx * pa_z[j] * pb_z[j] + 1.125 * fl3_fx * pb_zz[j] + 0.75 * pa_yyy[j] * fl2_fx * pb_y[j] + 2.25 * pa_yyz[j] * fl2_fx * pb_z[j] + 2.25 * pa_yy[j] * fl2_fx * pb_zz[j] + 2.25 * pa_yz[j] * fl2_fx * pb_yz[j] + 2.25 * pa_y[j] * fl2_fx * pb_yzz[j] + 0.75 * fl2_fx * pa_z[j] * pb_zzz[j] + 1.5 * pa_yyyz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_yyy[j] * fl1_fx * pb_yzz[j] + 1.5 * pa_yyz[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_yzzz[j] + pa_yyyz[j] * pb_yzzz[j]);

                t_yyyz_zzzz[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 4.5 * pa_y[j] * fl3_fx * pb_z[j] + 0.75 * pa_yyyz[j] * fl2_fx + 3.0 * pa_yyy[j] * fl2_fx * pb_z[j] + 4.5 * pa_yz[j] * fl2_fx * pb_zz[j] + 3.0 * pa_y[j] * fl2_fx * pb_zzz[j] + 3.0 * pa_yyyz[j] * pb_zz[j] * fl1_fx + 2.0 * pa_yyy[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_zzzz[j] + pa_yyyz[j] * pb_zzzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_180_189(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (180,189)

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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_yyzz = paDistances.data(34 * idx + 31);

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

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yyzz_xxxx = primBuffer.data(225 * idx + 180);

            auto t_yyzz_xxxy = primBuffer.data(225 * idx + 181);

            auto t_yyzz_xxxz = primBuffer.data(225 * idx + 182);

            auto t_yyzz_xxyy = primBuffer.data(225 * idx + 183);

            auto t_yyzz_xxyz = primBuffer.data(225 * idx + 184);

            auto t_yyzz_xxzz = primBuffer.data(225 * idx + 185);

            auto t_yyzz_xyyy = primBuffer.data(225 * idx + 186);

            auto t_yyzz_xyyz = primBuffer.data(225 * idx + 187);

            auto t_yyzz_xyzz = primBuffer.data(225 * idx + 188);

            // Batch of Integrals (180,189)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, \
                                     t_yyzz_xxxx, t_yyzz_xxxy, t_yyzz_xxxz, t_yyzz_xxyy, t_yyzz_xxyz, t_yyzz_xxzz, \
                                     t_yyzz_xyyy, t_yyzz_xyyz, t_yyzz_xyzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyzz_xxxx[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.375 * fl3_fx * pa_zz[j] + 0.75 * pa_yyzz[j] * fl2_fx + 0.75 * fl3_fx * pb_xx[j] + 1.5 * pa_yy[j] * fl2_fx * pb_xx[j] + 1.5 * fl2_fx * pa_zz[j] * pb_xx[j] + 3.0 * pa_yyzz[j] * pb_xx[j] * fl1_fx + 0.25 * fl2_fx * pb_xxxx[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xxxx[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxxx[j] + pa_yyzz[j] * pb_xxxx[j]);

                t_yyzz_xxxy[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx * pb_x[j] + 1.5 * pa_yzz[j] * fl2_fx * pb_x[j] + 0.375 * fl3_fx * pb_xy[j] + 0.75 * pa_yy[j] * fl2_fx * pb_xy[j] + 0.5 * pa_y[j] * fl2_fx * pb_xxx[j] + 0.75 * fl2_fx * pa_zz[j] * pb_xy[j] + 1.5 * pa_yyzz[j] * pb_xy[j] * fl1_fx + pa_yzz[j] * fl1_fx * pb_xxx[j] + 0.25 * fl2_fx * pb_xxxy[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xxxy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxxy[j] + pa_yyzz[j] * pb_xxxy[j]);

                t_yyzz_xxxz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] * pb_x[j] + 1.5 * pa_yyz[j] * fl2_fx * pb_x[j] + 0.375 * fl3_fx * pb_xz[j] + 0.75 * pa_yy[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_xz[j] + 0.5 * fl2_fx * pa_z[j] * pb_xxx[j] + 1.5 * pa_yyzz[j] * pb_xz[j] * fl1_fx + pa_yyz[j] * fl1_fx * pb_xxx[j] + 0.25 * fl2_fx * pb_xxxz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xxxz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxxz[j] + pa_yyzz[j] * pb_xxxz[j]);

                t_yyzz_xxyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_zz[j] + 0.125 * pa_yy[j] * fl3_fx + 0.5 * pa_y[j] * fl3_fx * pb_y[j] + 0.375 * fl3_fx * pb_xx[j] + 0.25 * pa_yyzz[j] * fl2_fx + pa_yzz[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pa_zz[j] * pb_xx[j] + 0.125 * fl3_fx * pb_yy[j] + 0.25 * pa_yy[j] * fl2_fx * pb_xx[j] + 0.25 * pa_yy[j] * fl2_fx * pb_yy[j] + pa_y[j] * fl2_fx * pb_xxy[j] + 0.25 * fl2_fx * pa_zz[j] * pb_yy[j] + 0.5 * pa_yyzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyzz[j] * fl1_fx * pb_yy[j] + 2.0 * pa_yzz[j] * fl1_fx * pb_xxy[j] + 0.25 * fl2_fx * pb_xxyy[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xxyy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxyy[j] + pa_yyzz[j] * pb_xxyy[j]);

                t_yyzz_xxyz[j] = fl_s_0_0 * (0.5 * pa_yz[j] * fl3_fx + 0.25 * pa_y[j] * fl3_fx * pb_z[j] + 0.25 * fl3_fx * pa_z[j] * pb_y[j] + 0.5 * pa_yyz[j] * fl2_fx * pb_y[j] + 0.5 * pa_yzz[j] * fl2_fx * pb_z[j] + pa_yz[j] * fl2_fx * pb_xx[j] + 0.125 * fl3_fx * pb_yz[j] + 0.25 * pa_yy[j] * fl2_fx * pb_yz[j] + 0.5 * pa_y[j] * fl2_fx * pb_xxz[j] + 0.25 * fl2_fx * pa_zz[j] * pb_yz[j] + 0.5 * fl2_fx * pa_z[j] * pb_xxy[j] + 0.5 * pa_yyzz[j] * fl1_fx * pb_yz[j] + pa_yyz[j] * fl1_fx * pb_xxy[j] + pa_yzz[j] * fl1_fx * pb_xxz[j] + 0.25 * fl2_fx * pb_xxyz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xxyz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxyz[j] + pa_yyzz[j] * pb_xxyz[j]);

                t_yyzz_xxzz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 0.125 * fl3_fx * pa_zz[j] + 0.5 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_xx[j] + 0.25 * pa_yyzz[j] * fl2_fx + pa_yyz[j] * fl2_fx * pb_z[j] + 0.75 * pa_yy[j] * fl2_fx * pb_xx[j] + 0.125 * fl3_fx * pb_zz[j] + 0.25 * pa_yy[j] * fl2_fx * pb_zz[j] + 0.25 * fl2_fx * pa_zz[j] * pb_xx[j] + 0.25 * fl2_fx * pa_zz[j] * pb_zz[j] + fl2_fx * pa_z[j] * pb_xxz[j] + 0.5 * pa_yyzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yyzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_yyz[j] * fl1_fx * pb_xxz[j] + 0.25 * fl2_fx * pb_xxzz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xxzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xxzz[j] + pa_yyzz[j] * pb_xxzz[j]);

                t_yyzz_xyyy[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx * pb_x[j] + 1.125 * fl3_fx * pb_xy[j] + 1.5 * pa_yzz[j] * fl2_fx * pb_x[j] + 2.25 * fl2_fx * pa_zz[j] * pb_xy[j] + 0.75 * pa_yy[j] * fl2_fx * pb_xy[j] + 1.5 * pa_y[j] * fl2_fx * pb_xyy[j] + 1.5 * pa_yyzz[j] * pb_xy[j] * fl1_fx + 3.0 * pa_yzz[j] * fl1_fx * pb_xyy[j] + 0.25 * fl2_fx * pb_xyyy[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xyyy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xyyy[j] + pa_yyzz[j] * pb_xyyy[j]);

                t_yyzz_xyyz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] * pb_x[j] + 0.375 * fl3_fx * pb_xz[j] + 0.5 * pa_yyz[j] * fl2_fx * pb_x[j] + 2.0 * pa_yz[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_zz[j] * pb_xz[j] + 0.25 * pa_yy[j] * fl2_fx * pb_xz[j] + pa_y[j] * fl2_fx * pb_xyz[j] + 0.5 * fl2_fx * pa_z[j] * pb_xyy[j] + 0.5 * pa_yyzz[j] * pb_xz[j] * fl1_fx + pa_yyz[j] * fl1_fx * pb_xyy[j] + 2.0 * pa_yzz[j] * fl1_fx * pb_xyz[j] + 0.25 * fl2_fx * pb_xyyz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xyyz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xyyz[j] + pa_yyzz[j] * pb_xyyz[j]);

                t_yyzz_xyzz[j] = fl_s_0_0 * (0.75 * pa_y[j] * fl3_fx * pb_x[j] + 0.375 * fl3_fx * pb_xy[j] + 0.75 * pa_yy[j] * fl2_fx * pb_xy[j] + 0.5 * pa_yzz[j] * fl2_fx * pb_x[j] + 2.0 * pa_yz[j] * fl2_fx * pb_xz[j] + 0.5 * pa_y[j] * fl2_fx * pb_xzz[j] + 0.25 * fl2_fx * pa_zz[j] * pb_xy[j] + fl2_fx * pa_z[j] * pb_xyz[j] + 0.5 * pa_yyzz[j] * pb_xy[j] * fl1_fx + 2.0 * pa_yyz[j] * fl1_fx * pb_xyz[j] + pa_yzz[j] * fl1_fx * pb_xzz[j] + 0.25 * fl2_fx * pb_xyzz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xyzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xyzz[j] + pa_yyzz[j] * pb_xyzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_189_198(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (189,198)

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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_yyzz = paDistances.data(34 * idx + 31);

            auto pa_yzzz = paDistances.data(34 * idx + 32);

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

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxxx = pbDistances.data(34 * idx + 19);

            auto pb_xxxy = pbDistances.data(34 * idx + 20);

            auto pb_xxxz = pbDistances.data(34 * idx + 21);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yyzz_xzzz = primBuffer.data(225 * idx + 189);

            auto t_yyzz_yyyy = primBuffer.data(225 * idx + 190);

            auto t_yyzz_yyyz = primBuffer.data(225 * idx + 191);

            auto t_yyzz_yyzz = primBuffer.data(225 * idx + 192);

            auto t_yyzz_yzzz = primBuffer.data(225 * idx + 193);

            auto t_yyzz_zzzz = primBuffer.data(225 * idx + 194);

            auto t_yzzz_xxxx = primBuffer.data(225 * idx + 195);

            auto t_yzzz_xxxy = primBuffer.data(225 * idx + 196);

            auto t_yzzz_xxxz = primBuffer.data(225 * idx + 197);

            // Batch of Integrals (189,198)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, \
                                     pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xy, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, \
                                     pb_zzz, pb_zzzz, s_0_0, t_yyzz_xzzz, t_yyzz_yyyy, t_yyzz_yyyz, t_yyzz_yyzz, \
                                     t_yyzz_yzzz, t_yyzz_zzzz, t_yzzz_xxxx, t_yzzz_xxxy, t_yzzz_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yyzz_xzzz[j] = fl_s_0_0 * (0.75 * fl3_fx * pa_z[j] * pb_x[j] + 1.125 * fl3_fx * pb_xz[j] + 1.5 * pa_yyz[j] * fl2_fx * pb_x[j] + 2.25 * pa_yy[j] * fl2_fx * pb_xz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_xz[j] + 1.5 * fl2_fx * pa_z[j] * pb_xzz[j] + 1.5 * pa_yyzz[j] * pb_xz[j] * fl1_fx + 3.0 * pa_yyz[j] * fl1_fx * pb_xzz[j] + 0.25 * fl2_fx * pb_xzzz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_xzzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_xzzz[j] + pa_yyzz[j] * pb_xzzz[j]);

                t_yyzz_yyyy[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * fl3_fx * pa_zz[j] + 0.375 * pa_yy[j] * fl3_fx + 3.0 * pa_y[j] * fl3_fx * pb_y[j] + 2.25 * fl3_fx * pb_yy[j] + 0.75 * pa_yyzz[j] * fl2_fx + 6.0 * pa_yzz[j] * fl2_fx * pb_y[j] + 4.5 * fl2_fx * pa_zz[j] * pb_yy[j] + 1.5 * pa_yy[j] * fl2_fx * pb_yy[j] + 2.0 * pa_y[j] * fl2_fx * pb_yyy[j] + 3.0 * pa_yyzz[j] * pb_yy[j] * fl1_fx + 4.0 * pa_yzz[j] * fl1_fx * pb_yyy[j] + 0.25 * fl2_fx * pb_yyyy[j] + 0.5 * pa_yy[j] * fl1_fx * pb_yyyy[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yyyy[j] + pa_yyzz[j] * pb_yyyy[j]);

                t_yyzz_yyyz[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl3_fx + 2.25 * fl3_fx * pa_z[j] * pb_y[j] + 0.75 * pa_y[j] * fl3_fx * pb_z[j] + 1.125 * fl3_fx * pb_yz[j] + 1.5 * pa_yyz[j] * fl2_fx * pb_y[j] + 1.5 * pa_yzz[j] * fl2_fx * pb_z[j] + 3.0 * pa_yz[j] * fl2_fx * pb_yy[j] + 2.25 * fl2_fx * pa_zz[j] * pb_yz[j] + 0.75 * pa_yy[j] * fl2_fx * pb_yz[j] + 1.5 * pa_y[j] * fl2_fx * pb_yyz[j] + 0.5 * fl2_fx * pa_z[j] * pb_yyy[j] + 1.5 * pa_yyzz[j] * pb_yz[j] * fl1_fx + pa_yyz[j] * fl1_fx * pb_yyy[j] + 3.0 * pa_yzz[j] * fl1_fx * pb_yyz[j] + 0.25 * fl2_fx * pb_yyyz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_yyyz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yyyz[j] + pa_yyzz[j] * pb_yyyz[j]);

                t_yyzz_yyzz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 0.375 * pa_yy[j] * fl3_fx + 1.5 * pa_y[j] * fl3_fx * pb_y[j] + 0.375 * fl3_fx * pa_zz[j] + 1.5 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_yy[j] + 0.375 * fl3_fx * pb_zz[j] + 0.25 * pa_yyzz[j] * fl2_fx + pa_yyz[j] * fl2_fx * pb_z[j] + 0.75 * pa_yy[j] * fl2_fx * pb_yy[j] + pa_yzz[j] * fl2_fx * pb_y[j] + 4.0 * pa_yz[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_zz[j] + 0.25 * pa_yy[j] * fl2_fx * pb_zz[j] + pa_y[j] * fl2_fx * pb_yzz[j] + 0.25 * fl2_fx * pa_zz[j] * pb_yy[j] + fl2_fx * pa_z[j] * pb_yyz[j] + 0.5 * pa_yyzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_yyzz[j] * fl1_fx * pb_zz[j] + 2.0 * pa_yyz[j] * fl1_fx * pb_yyz[j] + 2.0 * pa_yzz[j] * fl1_fx * pb_yzz[j] + 0.25 * fl2_fx * pb_yyzz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_yyzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yyzz[j] + pa_yyzz[j] * pb_yyzz[j]);

                t_yyzz_yzzz[j] = fl_s_0_0 * (1.5 * pa_yz[j] * fl3_fx + 2.25 * pa_y[j] * fl3_fx * pb_z[j] + 0.75 * fl3_fx * pa_z[j] * pb_y[j] + 1.125 * fl3_fx * pb_yz[j] + 1.5 * pa_yyz[j] * fl2_fx * pb_y[j] + 2.25 * pa_yy[j] * fl2_fx * pb_yz[j] + 1.5 * pa_yzz[j] * fl2_fx * pb_z[j] + 3.0 * pa_yz[j] * fl2_fx * pb_zz[j] + 0.5 * pa_y[j] * fl2_fx * pb_zzz[j] + 0.75 * fl2_fx * pa_zz[j] * pb_yz[j] + 1.5 * fl2_fx * pa_z[j] * pb_yzz[j] + 1.5 * pa_yyzz[j] * pb_yz[j] * fl1_fx + 3.0 * pa_yyz[j] * fl1_fx * pb_yzz[j] + pa_yzz[j] * fl1_fx * pb_zzz[j] + 0.25 * fl2_fx * pb_yzzz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_yzzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_yzzz[j] + pa_yyzz[j] * pb_yzzz[j]);

                t_yyzz_zzzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_yy[j] * fl3_fx + 0.375 * fl3_fx * pa_zz[j] + 3.0 * fl3_fx * pa_z[j] * pb_z[j] + 2.25 * fl3_fx * pb_zz[j] + 0.75 * pa_yyzz[j] * fl2_fx + 6.0 * pa_yyz[j] * fl2_fx * pb_z[j] + 4.5 * pa_yy[j] * fl2_fx * pb_zz[j] + 1.5 * fl2_fx * pa_zz[j] * pb_zz[j] + 2.0 * fl2_fx * pa_z[j] * pb_zzz[j] + 3.0 * pa_yyzz[j] * pb_zz[j] * fl1_fx + 4.0 * pa_yyz[j] * fl1_fx * pb_zzz[j] + 0.25 * fl2_fx * pb_zzzz[j] + 0.5 * pa_yy[j] * fl1_fx * pb_zzzz[j] + 0.5 * fl1_fx * pa_zz[j] * pb_zzzz[j] + pa_yyzz[j] * pb_zzzz[j]);

                t_yzzz_xxxx[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 0.75 * pa_yzzz[j] * fl2_fx + 4.5 * pa_yz[j] * fl2_fx * pb_xx[j] + 3.0 * pa_yzzz[j] * pb_xx[j] * fl1_fx + 1.5 * pa_yz[j] * fl1_fx * pb_xxxx[j] + pa_yzzz[j] * pb_xxxx[j]);

                t_yzzz_xxxy[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * fl2_fx * pa_zzz[j] * pb_x[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xy[j] + 0.75 * fl2_fx * pa_z[j] * pb_xxx[j] + 1.5 * pa_yzzz[j] * pb_xy[j] * fl1_fx + 0.5 * fl1_fx * pa_zzz[j] * pb_xxx[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxxy[j] + pa_yzzz[j] * pb_xxxy[j]);

                t_yzzz_xxxz[j] = fl_s_0_0 * (1.125 * pa_y[j] * fl3_fx * pb_x[j] + 2.25 * pa_yzz[j] * fl2_fx * pb_x[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xz[j] + 0.75 * pa_y[j] * fl2_fx * pb_xxx[j] + 1.5 * pa_yzzz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_yzz[j] * fl1_fx * pb_xxx[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxxz[j] + pa_yzzz[j] * pb_xxxz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_198_207(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (198,207)

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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_yzzz = paDistances.data(34 * idx + 32);

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

            auto pb_xxy = pbDistances.data(34 * idx + 10);

            auto pb_xxz = pbDistances.data(34 * idx + 11);

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            // set up pointers to 4-th order tensor of distance R(PB)

            auto pb_xxyy = pbDistances.data(34 * idx + 22);

            auto pb_xxyz = pbDistances.data(34 * idx + 23);

            auto pb_xxzz = pbDistances.data(34 * idx + 24);

            auto pb_xyyy = pbDistances.data(34 * idx + 25);

            auto pb_xyyz = pbDistances.data(34 * idx + 26);

            auto pb_xyzz = pbDistances.data(34 * idx + 27);

            auto pb_xzzz = pbDistances.data(34 * idx + 28);

            auto pb_yyyy = pbDistances.data(34 * idx + 29);

            auto pb_yyyz = pbDistances.data(34 * idx + 30);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yzzz_xxyy = primBuffer.data(225 * idx + 198);

            auto t_yzzz_xxyz = primBuffer.data(225 * idx + 199);

            auto t_yzzz_xxzz = primBuffer.data(225 * idx + 200);

            auto t_yzzz_xyyy = primBuffer.data(225 * idx + 201);

            auto t_yzzz_xyyz = primBuffer.data(225 * idx + 202);

            auto t_yzzz_xyzz = primBuffer.data(225 * idx + 203);

            auto t_yzzz_xzzz = primBuffer.data(225 * idx + 204);

            auto t_yzzz_yyyy = primBuffer.data(225 * idx + 205);

            auto t_yzzz_yyyz = primBuffer.data(225 * idx + 206);

            // Batch of Integrals (198,207)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, \
                                     pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_yzzz_xxyy, t_yzzz_xxyz, t_yzzz_xxzz, t_yzzz_xyyy, t_yzzz_xyyz, \
                                     t_yzzz_xyzz, t_yzzz_xzzz, t_yzzz_yyyy, t_yzzz_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yzzz_xxyy[j] = fl_s_0_0 * (0.375 * pa_yz[j] * fl3_fx + 0.75 * fl3_fx * pa_z[j] * pb_y[j] + 0.25 * pa_yzzz[j] * fl2_fx + 0.5 * fl2_fx * pa_zzz[j] * pb_y[j] + 0.75 * pa_yz[j] * fl2_fx * pb_xx[j] + 0.75 * pa_yz[j] * fl2_fx * pb_yy[j] + 1.5 * fl2_fx * pa_z[j] * pb_xxy[j] + 0.5 * pa_yzzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yzzz[j] * fl1_fx * pb_yy[j] + fl1_fx * pa_zzz[j] * pb_xxy[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxyy[j] + pa_yzzz[j] * pb_xxyy[j]);

                t_yzzz_xxyz[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.375 * fl3_fx * pa_zz[j] + 0.375 * pa_y[j] * fl3_fx * pb_y[j] + 0.375 * fl3_fx * pa_z[j] * pb_z[j] + 0.375 * fl3_fx * pb_xx[j] + 0.75 * pa_yzz[j] * fl2_fx * pb_y[j] + 0.25 * fl2_fx * pa_zzz[j] * pb_z[j] + 0.75 * fl2_fx * pa_zz[j] * pb_xx[j] + 0.75 * pa_yz[j] * fl2_fx * pb_yz[j] + 0.75 * pa_y[j] * fl2_fx * pb_xxy[j] + 0.75 * fl2_fx * pa_z[j] * pb_xxz[j] + 0.5 * pa_yzzz[j] * fl1_fx * pb_yz[j] + 1.5 * pa_yzz[j] * fl1_fx * pb_xxy[j] + 0.5 * fl1_fx * pa_zzz[j] * pb_xxz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxyz[j] + pa_yzzz[j] * pb_xxyz[j]);

                t_yzzz_xxzz[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 0.75 * pa_y[j] * fl3_fx * pb_z[j] + 0.25 * pa_yzzz[j] * fl2_fx + 1.5 * pa_yzz[j] * fl2_fx * pb_z[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xx[j] + 0.75 * pa_yz[j] * fl2_fx * pb_zz[j] + 1.5 * pa_y[j] * fl2_fx * pb_xxz[j] + 0.5 * pa_yzzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_yzzz[j] * fl1_fx * pb_zz[j] + 3.0 * pa_yzz[j] * fl1_fx * pb_xxz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xxzz[j] + pa_yzzz[j] * pb_xxzz[j]);

                t_yzzz_xyyy[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * fl2_fx * pa_zzz[j] * pb_x[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xy[j] + 2.25 * fl2_fx * pa_z[j] * pb_xyy[j] + 1.5 * pa_yzzz[j] * pb_xy[j] * fl1_fx + 1.5 * fl1_fx * pa_zzz[j] * pb_xyy[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xyyy[j] + pa_yzzz[j] * pb_xyyy[j]);

                t_yzzz_xyyz[j] = fl_s_0_0 * (0.375 * pa_y[j] * fl3_fx * pb_x[j] + 0.75 * fl3_fx * pb_xy[j] + 0.75 * pa_yzz[j] * fl2_fx * pb_x[j] + 1.5 * fl2_fx * pa_zz[j] * pb_xy[j] + 0.75 * pa_yz[j] * fl2_fx * pb_xz[j] + 0.75 * pa_y[j] * fl2_fx * pb_xyy[j] + 1.5 * fl2_fx * pa_z[j] * pb_xyz[j] + 0.5 * pa_yzzz[j] * pb_xz[j] * fl1_fx + 1.5 * pa_yzz[j] * fl1_fx * pb_xyy[j] + fl1_fx * pa_zzz[j] * pb_xyz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xyyz[j] + pa_yzzz[j] * pb_xyyz[j]);

                t_yzzz_xyzz[j] = fl_s_0_0 * (1.125 * fl3_fx * pa_z[j] * pb_x[j] + 0.75 * fl3_fx * pb_xz[j] + 2.25 * pa_yz[j] * fl2_fx * pb_xy[j] + 0.25 * fl2_fx * pa_zzz[j] * pb_x[j] + 1.5 * fl2_fx * pa_zz[j] * pb_xz[j] + 1.5 * pa_y[j] * fl2_fx * pb_xyz[j] + 0.75 * fl2_fx * pa_z[j] * pb_xzz[j] + 0.5 * pa_yzzz[j] * pb_xy[j] * fl1_fx + 3.0 * pa_yzz[j] * fl1_fx * pb_xyz[j] + 0.5 * fl1_fx * pa_zzz[j] * pb_xzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xyzz[j] + pa_yzzz[j] * pb_xyzz[j]);

                t_yzzz_xzzz[j] = fl_s_0_0 * (1.875 * pa_y[j] * fl3_fx * pb_x[j] + 2.25 * pa_yzz[j] * fl2_fx * pb_x[j] + 6.75 * pa_yz[j] * fl2_fx * pb_xz[j] + 2.25 * pa_y[j] * fl2_fx * pb_xzz[j] + 1.5 * pa_yzzz[j] * pb_xz[j] * fl1_fx + 4.5 * pa_yzz[j] * fl1_fx * pb_xzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_xzzz[j] + pa_yzzz[j] * pb_xzzz[j]);

                t_yzzz_yyyy[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 4.5 * fl3_fx * pa_z[j] * pb_y[j] + 0.75 * pa_yzzz[j] * fl2_fx + 3.0 * fl2_fx * pa_zzz[j] * pb_y[j] + 4.5 * pa_yz[j] * fl2_fx * pb_yy[j] + 3.0 * fl2_fx * pa_z[j] * pb_yyy[j] + 3.0 * pa_yzzz[j] * pb_yy[j] * fl1_fx + 2.0 * fl1_fx * pa_zzz[j] * pb_yyy[j] + 1.5 * pa_yz[j] * fl1_fx * pb_yyyy[j] + pa_yzzz[j] * pb_yyyy[j]);

                t_yzzz_yyyz[j] = fl_s_0_0 * (0.5625 * fl4_fx + 1.125 * fl3_fx * pa_zz[j] + 1.125 * pa_y[j] * fl3_fx * pb_y[j] + 1.125 * fl3_fx * pa_z[j] * pb_z[j] + 1.125 * fl3_fx * pb_yy[j] + 2.25 * pa_yzz[j] * fl2_fx * pb_y[j] + 0.75 * fl2_fx * pa_zzz[j] * pb_z[j] + 2.25 * fl2_fx * pa_zz[j] * pb_yy[j] + 2.25 * pa_yz[j] * fl2_fx * pb_yz[j] + 0.75 * pa_y[j] * fl2_fx * pb_yyy[j] + 2.25 * fl2_fx * pa_z[j] * pb_yyz[j] + 1.5 * pa_yzzz[j] * pb_yz[j] * fl1_fx + 1.5 * pa_yzz[j] * fl1_fx * pb_yyy[j] + 1.5 * fl1_fx * pa_zzz[j] * pb_yyz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_yyyz[j] + pa_yzzz[j] * pb_yyyz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_207_216(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (207,216)

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

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

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

            auto pb_yyzz = pbDistances.data(34 * idx + 31);

            auto pb_yzzz = pbDistances.data(34 * idx + 32);

            auto pb_zzzz = pbDistances.data(34 * idx + 33);

            // set up pointers to auxilary integrals

            auto s_0_0 = auxBuffer.data(idx);

            // set up pointers to integrals

            auto t_yzzz_yyzz = primBuffer.data(225 * idx + 207);

            auto t_yzzz_yzzz = primBuffer.data(225 * idx + 208);

            auto t_yzzz_zzzz = primBuffer.data(225 * idx + 209);

            auto t_zzzz_xxxx = primBuffer.data(225 * idx + 210);

            auto t_zzzz_xxxy = primBuffer.data(225 * idx + 211);

            auto t_zzzz_xxxz = primBuffer.data(225 * idx + 212);

            auto t_zzzz_xxyy = primBuffer.data(225 * idx + 213);

            auto t_zzzz_xxyz = primBuffer.data(225 * idx + 214);

            auto t_zzzz_xxzz = primBuffer.data(225 * idx + 215);

            // Batch of Integrals (207,216)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xx, \
                                     pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, \
                                     pb_xz, pb_y, pb_yy, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, \
                                     pb_zzzz, s_0_0, t_yzzz_yyzz, t_yzzz_yzzz, t_yzzz_zzzz, t_zzzz_xxxx, t_zzzz_xxxy, \
                                     t_zzzz_xxxz, t_zzzz_xxyy, t_zzzz_xxyz, t_zzzz_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_yzzz_yyzz[j] = fl_s_0_0 * (1.125 * pa_yz[j] * fl3_fx + 2.25 * fl3_fx * pa_z[j] * pb_y[j] + 0.75 * pa_y[j] * fl3_fx * pb_z[j] + 1.5 * fl3_fx * pb_yz[j] + 0.25 * pa_yzzz[j] * fl2_fx + 1.5 * pa_yzz[j] * fl2_fx * pb_z[j] + 2.25 * pa_yz[j] * fl2_fx * pb_yy[j] + 0.5 * fl2_fx * pa_zzz[j] * pb_y[j] + 3.0 * fl2_fx * pa_zz[j] * pb_yz[j] + 0.75 * pa_yz[j] * fl2_fx * pb_zz[j] + 1.5 * pa_y[j] * fl2_fx * pb_yyz[j] + 1.5 * fl2_fx * pa_z[j] * pb_yzz[j] + 0.5 * pa_yzzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_yzzz[j] * fl1_fx * pb_zz[j] + 3.0 * pa_yzz[j] * fl1_fx * pb_yyz[j] + fl1_fx * pa_zzz[j] * pb_yzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_yyzz[j] + pa_yzzz[j] * pb_yyzz[j]);

                t_yzzz_yzzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 1.875 * pa_y[j] * fl3_fx * pb_y[j] + 1.125 * fl3_fx * pa_zz[j] + 3.375 * fl3_fx * pa_z[j] * pb_z[j] + 1.125 * fl3_fx * pb_zz[j] + 2.25 * pa_yzz[j] * fl2_fx * pb_y[j] + 6.75 * pa_yz[j] * fl2_fx * pb_yz[j] + 0.75 * fl2_fx * pa_zzz[j] * pb_z[j] + 2.25 * fl2_fx * pa_zz[j] * pb_zz[j] + 2.25 * pa_y[j] * fl2_fx * pb_yzz[j] + 0.75 * fl2_fx * pa_z[j] * pb_zzz[j] + 1.5 * pa_yzzz[j] * pb_yz[j] * fl1_fx + 4.5 * pa_yzz[j] * fl1_fx * pb_yzz[j] + 0.5 * fl1_fx * pa_zzz[j] * pb_zzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_yzzz[j] + pa_yzzz[j] * pb_yzzz[j]);

                t_yzzz_zzzz[j] = fl_s_0_0 * (5.625 * pa_yz[j] * fl3_fx + 7.5 * pa_y[j] * fl3_fx * pb_z[j] + 0.75 * pa_yzzz[j] * fl2_fx + 9.0 * pa_yzz[j] * fl2_fx * pb_z[j] + 13.5 * pa_yz[j] * fl2_fx * pb_zz[j] + 3.0 * pa_y[j] * fl2_fx * pb_zzz[j] + 3.0 * pa_yzzz[j] * pb_zz[j] * fl1_fx + 6.0 * pa_yzz[j] * fl1_fx * pb_zzz[j] + 1.5 * pa_yz[j] * fl1_fx * pb_zzzz[j] + pa_yzzz[j] * pb_zzzz[j]);

                t_zzzz_xxxx[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_zz[j] * fl3_fx + 0.75 * pa_zzzz[j] * fl2_fx + 2.25 * fl3_fx * pb_xx[j] + 9.0 * pa_zz[j] * fl2_fx * pb_xx[j] + 3.0 * pa_zzzz[j] * pb_xx[j] * fl1_fx + 0.75 * fl2_fx * pb_xxxx[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xxxx[j] + pa_zzzz[j] * pb_xxxx[j]);

                t_zzzz_xxxy[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_xy[j] + 4.5 * pa_zz[j] * fl2_fx * pb_xy[j] + 1.5 * pa_zzzz[j] * pb_xy[j] * fl1_fx + 0.75 * fl2_fx * pb_xxxy[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xxxy[j] + pa_zzzz[j] * pb_xxxy[j]);

                t_zzzz_xxxz[j] = fl_s_0_0 * (4.5 * pa_z[j] * fl3_fx * pb_x[j] + 3.0 * pa_zzz[j] * fl2_fx * pb_x[j] + 1.125 * fl3_fx * pb_xz[j] + 4.5 * pa_zz[j] * fl2_fx * pb_xz[j] + 3.0 * pa_z[j] * fl2_fx * pb_xxx[j] + 1.5 * pa_zzzz[j] * pb_xz[j] * fl1_fx + 2.0 * pa_zzz[j] * fl1_fx * pb_xxx[j] + 0.75 * fl2_fx * pb_xxxz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xxxz[j] + pa_zzzz[j] * pb_xxxz[j]);

                t_zzzz_xxyy[j] = fl_s_0_0 * (0.1875 * fl4_fx + 0.75 * pa_zz[j] * fl3_fx + 0.25 * pa_zzzz[j] * fl2_fx + 0.375 * fl3_fx * pb_xx[j] + 0.375 * fl3_fx * pb_yy[j] + 1.5 * pa_zz[j] * fl2_fx * pb_xx[j] + 1.5 * pa_zz[j] * fl2_fx * pb_yy[j] + 0.5 * pa_zzzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_zzzz[j] * fl1_fx * pb_yy[j] + 0.75 * fl2_fx * pb_xxyy[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xxyy[j] + pa_zzzz[j] * pb_xxyy[j]);

                t_zzzz_xxyz[j] = fl_s_0_0 * (1.5 * pa_z[j] * fl3_fx * pb_y[j] + pa_zzz[j] * fl2_fx * pb_y[j] + 0.375 * fl3_fx * pb_yz[j] + 1.5 * pa_zz[j] * fl2_fx * pb_yz[j] + 3.0 * pa_z[j] * fl2_fx * pb_xxy[j] + 0.5 * pa_zzzz[j] * fl1_fx * pb_yz[j] + 2.0 * pa_zzz[j] * fl1_fx * pb_xxy[j] + 0.75 * fl2_fx * pb_xxyz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xxyz[j] + pa_zzzz[j] * pb_xxyz[j]);

                t_zzzz_xxzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_zz[j] * fl3_fx + 3.0 * pa_z[j] * fl3_fx * pb_z[j] + 1.875 * fl3_fx * pb_xx[j] + 0.25 * pa_zzzz[j] * fl2_fx + 2.0 * pa_zzz[j] * fl2_fx * pb_z[j] + 4.5 * pa_zz[j] * fl2_fx * pb_xx[j] + 0.375 * fl3_fx * pb_zz[j] + 1.5 * pa_zz[j] * fl2_fx * pb_zz[j] + 6.0 * pa_z[j] * fl2_fx * pb_xxz[j] + 0.5 * pa_zzzz[j] * pb_xx[j] * fl1_fx + 0.5 * pa_zzzz[j] * fl1_fx * pb_zz[j] + 4.0 * pa_zzz[j] * fl1_fx * pb_xxz[j] + 0.75 * fl2_fx * pb_xxzz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xxzz[j] + pa_zzzz[j] * pb_xxzz[j]);
            }

            idx++;
        }
    }

    void
    compOverlapForGG_216_225(      CMemBlock2D<double>& primBuffer,
                             const CMemBlock2D<double>& auxBuffer,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        // Batch of Integrals (216,225)

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

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

            auto pa_zzzz = paDistances.data(34 * idx + 33);

            // set up pointers to 1-th order tensor of distance R(PB)

            auto pb_x = pbDistances.data(34 * idx);

            auto pb_y = pbDistances.data(34 * idx + 1);

            auto pb_z = pbDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PB)

            auto pb_xy = pbDistances.data(34 * idx + 4);

            auto pb_xz = pbDistances.data(34 * idx + 5);

            auto pb_yy = pbDistances.data(34 * idx + 6);

            auto pb_yz = pbDistances.data(34 * idx + 7);

            auto pb_zz = pbDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PB)

            auto pb_xyy = pbDistances.data(34 * idx + 12);

            auto pb_xyz = pbDistances.data(34 * idx + 13);

            auto pb_xzz = pbDistances.data(34 * idx + 14);

            auto pb_yyy = pbDistances.data(34 * idx + 15);

            auto pb_yyz = pbDistances.data(34 * idx + 16);

            auto pb_yzz = pbDistances.data(34 * idx + 17);

            auto pb_zzz = pbDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PB)

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

            auto t_zzzz_xyyy = primBuffer.data(225 * idx + 216);

            auto t_zzzz_xyyz = primBuffer.data(225 * idx + 217);

            auto t_zzzz_xyzz = primBuffer.data(225 * idx + 218);

            auto t_zzzz_xzzz = primBuffer.data(225 * idx + 219);

            auto t_zzzz_yyyy = primBuffer.data(225 * idx + 220);

            auto t_zzzz_yyyz = primBuffer.data(225 * idx + 221);

            auto t_zzzz_yyzz = primBuffer.data(225 * idx + 222);

            auto t_zzzz_yzzz = primBuffer.data(225 * idx + 223);

            auto t_zzzz_zzzz = primBuffer.data(225 * idx + 224);

            // Batch of Integrals (216,225)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, \
                                     pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_zzzz_xyyy, \
                                     t_zzzz_xyyz, t_zzzz_xyzz, t_zzzz_xzzz, t_zzzz_yyyy, t_zzzz_yyyz, t_zzzz_yyzz, \
                                     t_zzzz_yzzz, t_zzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl_s_0_0 = s_0_0[j];

                double fl1_fx = fx[j];

                double fl2_fx = fx[j] * fx[j];

                double fl3_fx = fx[j] * fx[j] * fx[j];

                double fl4_fx = fx[j] * fx[j] * fx[j] * fx[j];

                t_zzzz_xyyy[j] = fl_s_0_0 * (1.125 * fl3_fx * pb_xy[j] + 4.5 * pa_zz[j] * fl2_fx * pb_xy[j] + 1.5 * pa_zzzz[j] * pb_xy[j] * fl1_fx + 0.75 * fl2_fx * pb_xyyy[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xyyy[j] + pa_zzzz[j] * pb_xyyy[j]);

                t_zzzz_xyyz[j] = fl_s_0_0 * (1.5 * pa_z[j] * fl3_fx * pb_x[j] + pa_zzz[j] * fl2_fx * pb_x[j] + 0.375 * fl3_fx * pb_xz[j] + 1.5 * pa_zz[j] * fl2_fx * pb_xz[j] + 3.0 * pa_z[j] * fl2_fx * pb_xyy[j] + 0.5 * pa_zzzz[j] * pb_xz[j] * fl1_fx + 2.0 * pa_zzz[j] * fl1_fx * pb_xyy[j] + 0.75 * fl2_fx * pb_xyyz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xyyz[j] + pa_zzzz[j] * pb_xyyz[j]);

                t_zzzz_xyzz[j] = fl_s_0_0 * (1.875 * fl3_fx * pb_xy[j] + 4.5 * pa_zz[j] * fl2_fx * pb_xy[j] + 6.0 * pa_z[j] * fl2_fx * pb_xyz[j] + 0.5 * pa_zzzz[j] * pb_xy[j] * fl1_fx + 4.0 * pa_zzz[j] * fl1_fx * pb_xyz[j] + 0.75 * fl2_fx * pb_xyzz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xyzz[j] + pa_zzzz[j] * pb_xyzz[j]);

                t_zzzz_xzzz[j] = fl_s_0_0 * (7.5 * pa_z[j] * fl3_fx * pb_x[j] + 5.625 * fl3_fx * pb_xz[j] + 3.0 * pa_zzz[j] * fl2_fx * pb_x[j] + 13.5 * pa_zz[j] * fl2_fx * pb_xz[j] + 9.0 * pa_z[j] * fl2_fx * pb_xzz[j] + 1.5 * pa_zzzz[j] * pb_xz[j] * fl1_fx + 6.0 * pa_zzz[j] * fl1_fx * pb_xzz[j] + 0.75 * fl2_fx * pb_xzzz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_xzzz[j] + pa_zzzz[j] * pb_xzzz[j]);

                t_zzzz_yyyy[j] = fl_s_0_0 * (0.5625 * fl4_fx + 2.25 * pa_zz[j] * fl3_fx + 0.75 * pa_zzzz[j] * fl2_fx + 2.25 * fl3_fx * pb_yy[j] + 9.0 * pa_zz[j] * fl2_fx * pb_yy[j] + 3.0 * pa_zzzz[j] * pb_yy[j] * fl1_fx + 0.75 * fl2_fx * pb_yyyy[j] + 3.0 * pa_zz[j] * fl1_fx * pb_yyyy[j] + pa_zzzz[j] * pb_yyyy[j]);

                t_zzzz_yyyz[j] = fl_s_0_0 * (4.5 * pa_z[j] * fl3_fx * pb_y[j] + 3.0 * pa_zzz[j] * fl2_fx * pb_y[j] + 1.125 * fl3_fx * pb_yz[j] + 4.5 * pa_zz[j] * fl2_fx * pb_yz[j] + 3.0 * pa_z[j] * fl2_fx * pb_yyy[j] + 1.5 * pa_zzzz[j] * pb_yz[j] * fl1_fx + 2.0 * pa_zzz[j] * fl1_fx * pb_yyy[j] + 0.75 * fl2_fx * pb_yyyz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_yyyz[j] + pa_zzzz[j] * pb_yyyz[j]);

                t_zzzz_yyzz[j] = fl_s_0_0 * (0.9375 * fl4_fx + 2.25 * pa_zz[j] * fl3_fx + 3.0 * pa_z[j] * fl3_fx * pb_z[j] + 1.875 * fl3_fx * pb_yy[j] + 0.25 * pa_zzzz[j] * fl2_fx + 2.0 * pa_zzz[j] * fl2_fx * pb_z[j] + 4.5 * pa_zz[j] * fl2_fx * pb_yy[j] + 0.375 * fl3_fx * pb_zz[j] + 1.5 * pa_zz[j] * fl2_fx * pb_zz[j] + 6.0 * pa_z[j] * fl2_fx * pb_yyz[j] + 0.5 * pa_zzzz[j] * pb_yy[j] * fl1_fx + 0.5 * pa_zzzz[j] * fl1_fx * pb_zz[j] + 4.0 * pa_zzz[j] * fl1_fx * pb_yyz[j] + 0.75 * fl2_fx * pb_yyzz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_yyzz[j] + pa_zzzz[j] * pb_yyzz[j]);

                t_zzzz_yzzz[j] = fl_s_0_0 * (7.5 * pa_z[j] * fl3_fx * pb_y[j] + 5.625 * fl3_fx * pb_yz[j] + 3.0 * pa_zzz[j] * fl2_fx * pb_y[j] + 13.5 * pa_zz[j] * fl2_fx * pb_yz[j] + 9.0 * pa_z[j] * fl2_fx * pb_yzz[j] + 1.5 * pa_zzzz[j] * pb_yz[j] * fl1_fx + 6.0 * pa_zzz[j] * fl1_fx * pb_yzz[j] + 0.75 * fl2_fx * pb_yzzz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_yzzz[j] + pa_zzzz[j] * pb_yzzz[j]);

                t_zzzz_zzzz[j] = fl_s_0_0 * (6.5625 * fl4_fx + 11.25 * pa_zz[j] * fl3_fx + 30.0 * pa_z[j] * fl3_fx * pb_z[j] + 11.25 * fl3_fx * pb_zz[j] + 0.75 * pa_zzzz[j] * fl2_fx + 12.0 * pa_zzz[j] * fl2_fx * pb_z[j] + 27.0 * pa_zz[j] * fl2_fx * pb_zz[j] + 12.0 * pa_z[j] * fl2_fx * pb_zzz[j] + 3.0 * pa_zzzz[j] * pb_zz[j] * fl1_fx + 8.0 * pa_zzz[j] * fl1_fx * pb_zzz[j] + 0.75 * fl2_fx * pb_zzzz[j] + 3.0 * pa_zz[j] * fl1_fx * pb_zzzz[j] + pa_zzzz[j] * pb_zzzz[j]);
            }

            idx++;
        }
    }


} // ovlrecfunc namespace

