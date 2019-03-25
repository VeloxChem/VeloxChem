//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapRecFuncForDX.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace

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

            // Batch of Integrals (0) = (0,12)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xy, pa_y, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_xx_xx, t_xx_xy, t_xx_xz, t_xx_yy, t_xx_yz, t_xx_zz, t_xy_xx, \
                                     t_xy_xy, t_xy_xz, t_xy_yy, t_xy_yz, t_xy_zz: VLX_ALIGN)
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

                t_xy_yy[j] = (0.5 * pa_xy[j] * fx[j] + pa_x[j] * fx[j] * pb_y[j] + pa_xy[j] * pb_yy[j]) * s_0_0[j];

                t_xy_yz[j] = (0.5 * pa_x[j] * fx[j] * pb_z[j] + pa_xy[j] * pb_yz[j]) * s_0_0[j];

                t_xy_zz[j] = (0.5 * pa_xy[j] * fx[j] + pa_xy[j] * pb_zz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (12,24)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_y, pa_yy, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, pb_zz, s_0_0, t_xz_xx, t_xz_xy, t_xz_xz, t_xz_yy, t_xz_yz, t_xz_zz, t_yy_xx, \
                                     t_yy_xy, t_yy_xz, t_yy_yy, t_yy_yz, t_yy_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_xx[j] = (0.5 * pa_xz[j] * fx[j] + fx[j] * pa_z[j] * pb_x[j] + pa_xz[j] * pb_xx[j]) * s_0_0[j];

                t_xz_xy[j] = (0.5 * fx[j] * pa_z[j] * pb_y[j] + pa_xz[j] * pb_xy[j]) * s_0_0[j];

                t_xz_xz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_x[j] + 

                             0.5 * fx[j] * pa_z[j] * pb_z[j] + pa_xz[j] * pb_xz[j]) * s_0_0[j];

                t_xz_yy[j] = (0.5 * pa_xz[j] * fx[j] + pa_xz[j] * pb_yy[j]) * s_0_0[j];

                t_xz_yz[j] = (0.5 * pa_x[j] * fx[j] * pb_y[j] + pa_xz[j] * pb_yz[j]) * s_0_0[j];

                t_xz_zz[j] = (0.5 * pa_xz[j] * fx[j] + pa_x[j] * fx[j] * pb_z[j] + pa_xz[j] * pb_zz[j]) * s_0_0[j];

                t_yy_xx[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] + 0.5 * fx[j] * pb_xx[j] + 

                             pa_yy[j] * pb_xx[j]) * s_0_0[j];

                t_yy_xy[j] = (pa_y[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pb_xy[j] + pa_yy[j] * pb_xy[j]) * s_0_0[j];

                t_yy_xz[j] = (0.5 * fx[j] * pb_xz[j] + pa_yy[j] * pb_xz[j]) * s_0_0[j];

                t_yy_yy[j] = (0.75 * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] + 2.0 * pa_y[j] * fx[j] * pb_y[j] + 

                             0.5 * fx[j] * pb_yy[j] + pa_yy[j] * pb_yy[j]) * s_0_0[j];

                t_yy_yz[j] = (pa_y[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pb_yz[j] + pa_yy[j] * pb_yz[j]) * s_0_0[j];

                t_yy_zz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] + 0.5 * fx[j] * pb_zz[j] + 

                             pa_yy[j] * pb_zz[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (24,36)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     pb_zz, s_0_0, t_yz_xx, t_yz_xy, t_yz_xz, t_yz_yy, t_yz_yz, t_yz_zz, t_zz_xx, \
                                     t_zz_xy, t_zz_xz, t_zz_yy, t_zz_yz, t_zz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_xx[j] = (0.5 * pa_yz[j] * fx[j] + pa_yz[j] * pb_xx[j]) * s_0_0[j];

                t_yz_xy[j] = (0.5 * fx[j] * pa_z[j] * pb_x[j] + pa_yz[j] * pb_xy[j]) * s_0_0[j];

                t_yz_xz[j] = (0.5 * pa_y[j] * fx[j] * pb_x[j] + pa_yz[j] * pb_xz[j]) * s_0_0[j];

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

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xx_xxx, \
                                     t_xx_xxy, t_xx_xxz, t_xx_xyy, t_xx_xyz, t_xx_xzz, t_xx_yyy, t_xx_yyz, t_xx_yzz, \
                                     t_xx_zzz: VLX_ALIGN)
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

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xy_xxx, t_xy_xxy, t_xy_xxz, t_xy_xyy, t_xy_xyz, t_xy_xzz, t_xy_yyy, t_xy_yyz, \
                                     t_xy_yzz, t_xy_zzz: VLX_ALIGN)
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

                t_xy_xzz[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xy[j] * pb_x[j] * fx[j] + 

                              0.5 * fx[j] * pa_y[j] * pb_zz[j] + pa_xy[j] * pb_xzz[j]) * s_0_0[j];

                t_xy_yyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 1.5 * pa_xy[j] * pb_y[j] * fx[j] + 

                              1.5 * pa_x[j] * fx[j] * pb_yy[j] + pa_xy[j] * pb_yyy[j]) * s_0_0[j];

                t_xy_yyz[j] = (0.5 * pa_xy[j] * fx[j] * pb_z[j] + pa_x[j] * fx[j] * pb_yz[j] + pa_xy[j] * pb_yyz[j]) * s_0_0[j];

                t_xy_yzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xy[j] * pb_y[j] * fx[j] + 

                              0.5 * pa_x[j] * fx[j] * pb_zz[j] + pa_xy[j] * pb_yzz[j]) * s_0_0[j];

                t_xy_zzz[j] = (1.5 * pa_xy[j] * pb_z[j] * fx[j] + pa_xy[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xz_xxx, t_xz_xxy, t_xz_xxz, t_xz_xyy, t_xz_xyz, t_xz_xzz, t_xz_yyy, t_xz_yyz, \
                                     t_xz_yzz, t_xz_zzz: VLX_ALIGN)
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

                t_xz_xzz[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * fx[j] * pb_z[j] + 

                              0.5 * pa_xz[j] * pb_x[j] * fx[j] + pa_x[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_z[j] * pb_zz[j] + pa_xz[j] * pb_xzz[j]) * s_0_0[j];

                t_xz_yyy[j] = (1.5 * pa_xz[j] * pb_y[j] * fx[j] + pa_xz[j] * pb_yyy[j]) * s_0_0[j];

                t_xz_yyz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xz[j] * fx[j] * pb_z[j] + 

                              0.5 * pa_x[j] * fx[j] * pb_yy[j] + pa_xz[j] * pb_yyz[j]) * s_0_0[j];

                t_xz_yzz[j] = (0.5 * pa_xz[j] * pb_y[j] * fx[j] + pa_x[j] * fx[j] * pb_yz[j] + pa_xz[j] * pb_yzz[j]) * s_0_0[j];

                t_xz_zzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 1.5 * pa_xz[j] * pb_z[j] * fx[j] + 

                              1.5 * pa_x[j] * fx[j] * pb_zz[j] + pa_xz[j] * pb_zzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (30,40)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_yy_xxx, \
                                     t_yy_xxy, t_yy_xxz, t_yy_xyy, t_yy_xyz, t_yy_xzz, t_yy_yyy, t_yy_yyz, t_yy_yzz, \
                                     t_yy_zzz: VLX_ALIGN)
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

            // Batch of Integrals (4) = (40,50)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_yz_xxx, t_yz_xxy, t_yz_xxz, t_yz_xyy, t_yz_xyz, t_yz_xzz, t_yz_yyy, t_yz_yyz, \
                                     t_yz_yzz, t_yz_zzz: VLX_ALIGN)
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

            // Batch of Integrals (5) = (50,60)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_zz_xxx, \
                                     t_zz_xxy, t_zz_xxz, t_zz_xyy, t_zz_xyz, t_zz_xzz, t_zz_yyy, t_zz_yyz, t_zz_yzz, \
                                     t_zz_zzz: VLX_ALIGN)
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

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_xxx_xx, t_xxx_xy, t_xxx_xz, t_xxx_yy, t_xxx_yz, \
                                     t_xxx_zz, t_xxy_xx, t_xxy_xy, t_xxy_xz, t_xxy_yy: VLX_ALIGN)
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

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxz, pa_xy, pa_xyy, pa_xz, pa_y, pa_yy, pa_z, pb_x, \
                                     pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xxy_yz, t_xxy_zz, \
                                     t_xxz_xx, t_xxz_xy, t_xxz_xz, t_xxz_yy, t_xxz_yz, t_xxz_zz, t_xyy_xx, t_xyy_xy: VLX_ALIGN)
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

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyz, pa_xz, pa_y, pa_yy, pa_yz, pa_z, pb_x, pb_xx, pb_xy, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xyy_xz, t_xyy_yy, t_xyy_yz, t_xyy_zz, \
                                     t_xyz_xx, t_xyz_xy, t_xyz_xz, t_xyz_yy, t_xyz_yz, t_xyz_zz: VLX_ALIGN)
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

                t_xyz_xy[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_xz[j] * fx[j] * pb_x[j] + 

                              0.5 * fx[j] * pa_yz[j] * pb_y[j] + pa_xyz[j] * pb_xy[j]) * s_0_0[j];

                t_xyz_xz[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xy[j] * fx[j] * pb_x[j] + 

                              0.5 * fx[j] * pa_yz[j] * pb_z[j] + pa_xyz[j] * pb_xz[j]) * s_0_0[j];

                t_xyz_yy[j] = (0.5 * pa_xyz[j] * fx[j] + pa_xz[j] * fx[j] * pb_y[j] + pa_xyz[j] * pb_yy[j]) * s_0_0[j];

                t_xyz_yz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * pb_y[j] + 

                              0.5 * pa_xz[j] * fx[j] * pb_z[j] + pa_xyz[j] * pb_yz[j]) * s_0_0[j];

                t_xyz_zz[j] = (0.5 * pa_xyz[j] * fx[j] + pa_xy[j] * fx[j] * pb_z[j] + pa_xyz[j] * pb_zz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (30,40)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_y, pa_yy, pa_yyy, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, \
                                     pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xzz_xx, t_xzz_xy, t_xzz_xz, t_xzz_yy, \
                                     t_xzz_yz, t_xzz_zz, t_yyy_xx, t_yyy_xy, t_yyy_xz, t_yyy_yy: VLX_ALIGN)
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

            // Batch of Integrals (4) = (40,50)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xy, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_yyy_yz, t_yyy_zz, t_yyz_xx, t_yyz_xy, \
                                     t_yyz_xz, t_yyz_yy, t_yyz_yz, t_yyz_zz, t_yzz_xx, t_yzz_xy: VLX_ALIGN)
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

            // Batch of Integrals (5) = (50,60)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, pb_zz, s_0_0, t_yzz_xz, t_yzz_yy, t_yzz_yz, t_yzz_zz, t_zzz_xx, \
                                     t_zzz_xy, t_zzz_xz, t_zzz_yy, t_zzz_yz, t_zzz_zz: VLX_ALIGN)
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

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, \
                                     pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xx_xxxx, t_xx_xxxy, t_xx_xxxz, t_xx_xxyy, t_xx_xxyz, t_xx_xxzz, t_xx_xyyy, \
                                     t_xx_xyyz, t_xx_xyzz, t_xx_xzzz: VLX_ALIGN)
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

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yyyy, \
                                     pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_z, pb_zz, pb_zzzz, s_0_0, t_xx_yyyy, t_xx_yyyz, \
                                     t_xx_yyzz, t_xx_yzzz, t_xx_zzzz, t_xy_xxxx, t_xy_xxxy, t_xy_xxxz, t_xy_xxyy, \
                                     t_xy_xxyz: VLX_ALIGN)
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

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, \
                                     pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xy_xxzz, \
                                     t_xy_xyyy, t_xy_xyyz, t_xy_xyzz, t_xy_xzzz, t_xy_yyyy, t_xy_yyyz, t_xy_yyzz, \
                                     t_xy_yzzz, t_xy_zzzz: VLX_ALIGN)
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

            // Batch of Integrals (3) = (30,40)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, \
                                     pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     s_0_0, t_xz_xxxx, t_xz_xxxy, t_xz_xxxz, t_xz_xxyy, t_xz_xxyz, t_xz_xxzz, \
                                     t_xz_xyyy, t_xz_xyyz, t_xz_xyzz, t_xz_xzzz: VLX_ALIGN)
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

            // Batch of Integrals (4) = (40,50)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, \
                                     t_xz_yyyy, t_xz_yyyz, t_xz_yyzz, t_xz_yzzz, t_xz_zzzz, t_yy_xxxx, t_yy_xxxy, \
                                     t_yy_xxxz, t_yy_xxyy, t_yy_xxyz: VLX_ALIGN)
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

            // Batch of Integrals (5) = (50,60)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, \
                                     pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_yy_xxzz, t_yy_xyyy, \
                                     t_yy_xyyz, t_yy_xyzz, t_yy_xzzz, t_yy_yyyy, t_yy_yyyz, t_yy_yyzz, t_yy_yzzz, \
                                     t_yy_zzzz: VLX_ALIGN)
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

            // Batch of Integrals (6) = (60,70)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, \
                                     pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, \
                                     pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_yz_xxxx, t_yz_xxxy, \
                                     t_yz_xxxz, t_yz_xxyy, t_yz_xxyz, t_yz_xxzz, t_yz_xyyy, t_yz_xyyz, t_yz_xyzz, \
                                     t_yz_xzzz: VLX_ALIGN)
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

            // Batch of Integrals (7) = (70,80)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_yz_yyyy, \
                                     t_yz_yyyz, t_yz_yyzz, t_yz_yzzz, t_yz_zzzz, t_zz_xxxx, t_zz_xxxy, t_zz_xxxz, \
                                     t_zz_xxyy, t_zz_xxyz: VLX_ALIGN)
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

            // Batch of Integrals (8) = (80,90)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, \
                                     pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_zz_xxzz, \
                                     t_zz_xyyy, t_zz_xyyz, t_zz_xyzz, t_zz_xzzz, t_zz_yyyy, t_zz_yyyz, t_zz_yyzz, \
                                     t_zz_yzzz, t_zz_zzzz: VLX_ALIGN)
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

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xy, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_xxxx_xx, t_xxxx_xy, t_xxxx_xz, \
                                     t_xxxx_yy, t_xxxx_yz, t_xxxx_zz, t_xxxy_xx, t_xxxy_xy, t_xxxy_xz, t_xxxy_yy: VLX_ALIGN)
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

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxxz, pa_xxy, pa_xxyy, pa_xxz, pa_xy, \
                                     pa_xyy, pa_xz, pa_y, pa_yy, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, \
                                     s_0_0, t_xxxy_yz, t_xxxy_zz, t_xxxz_xx, t_xxxz_xy, t_xxxz_xz, t_xxxz_yy, \
                                     t_xxxz_yz, t_xxxz_zz, t_xxyy_xx, t_xxyy_xy: VLX_ALIGN)
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

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xxyz, pa_xxz, pa_xy, pa_xyy, pa_xyz, pa_xz, \
                                     pa_y, pa_yy, pa_yz, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, \
                                     s_0_0, t_xxyy_xz, t_xxyy_yy, t_xxyy_yz, t_xxyy_zz, t_xxyz_xx, t_xxyz_xy, \
                                     t_xxyz_xz, t_xxyz_yy, t_xxyz_yz, t_xxyz_zz: VLX_ALIGN)
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

            // Batch of Integrals (3) = (30,40)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xy, pa_xyy, pa_xyyy, pa_xz, pa_xzz, pa_y, \
                                     pa_yy, pa_yyy, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, \
                                     s_0_0, t_xxzz_xx, t_xxzz_xy, t_xxzz_xz, t_xxzz_yy, t_xxzz_yz, t_xxzz_zz, \
                                     t_xyyy_xx, t_xyyy_xy, t_xyyy_xz, t_xyyy_yy: VLX_ALIGN)
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

            // Batch of Integrals (4) = (40,50)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_xyyz, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, \
                                     pa_yy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, pb_zz, s_0_0, t_xyyy_yz, t_xyyy_zz, t_xyyz_xx, t_xyyz_xy, t_xyyz_xz, \
                                     t_xyyz_yy, t_xyyz_yz, t_xyyz_zz, t_xyzz_xx, t_xyzz_xy: VLX_ALIGN)
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

            // Batch of Integrals (5) = (50,60)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_xzzz, pa_y, pa_yz, pa_yzz, \
                                     pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, \
                                     t_xyzz_xz, t_xyzz_yy, t_xyzz_yz, t_xyzz_zz, t_xzzz_xx, t_xzzz_xy, t_xzzz_xz, \
                                     t_xzzz_yy, t_xzzz_yz, t_xzzz_zz: VLX_ALIGN)
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

            // Batch of Integrals (6) = (60,70)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xy, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_yyyy_xx, t_yyyy_xy, t_yyyy_xz, \
                                     t_yyyy_yy, t_yyyy_yz, t_yyyy_zz, t_yyyz_xx, t_yyyz_xy, t_yyyz_xz, t_yyyz_yy: VLX_ALIGN)
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

            // Batch of Integrals (7) = (70,80)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_yzzz, pa_z, \
                                     pa_zz, pa_zzz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, \
                                     t_yyyz_yz, t_yyyz_zz, t_yyzz_xx, t_yyzz_xy, t_yyzz_xz, t_yyzz_yy, t_yyzz_yz, \
                                     t_yyzz_zz, t_yzzz_xx, t_yzzz_xy: VLX_ALIGN)
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

            // Batch of Integrals (8) = (80,90)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xx, pb_xy, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, s_0_0, t_yzzz_xz, t_yzzz_yy, t_yzzz_yz, \
                                     t_yzzz_zz, t_zzzz_xx, t_zzzz_xy, t_zzzz_xz, t_zzzz_yy, t_zzzz_yz, t_zzzz_zz: VLX_ALIGN)
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


} // ovlrecfunc namespace

