//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapRecFuncForPX.hpp"

namespace ovlrecfunc { // ovlrecfunc namespace

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

            // Batch of Integrals (0) = (0,15)

            #pragma omp simd aligned(fx, pa_x, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_x_xxx, \
                                     t_x_xxy, t_x_xxz, t_x_xyy, t_x_xyz, t_x_xzz, t_x_yyy, t_x_yyz, t_x_yzz, t_x_zzz, \
                                     t_y_xxx, t_y_xxy, t_y_xxz, t_y_xyy, t_y_xyz: VLX_ALIGN)
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

                t_y_xxx[j] = (1.5 * pa_y[j] * pb_x[j] * fx[j] + pa_y[j] * pb_xxx[j]) * s_0_0[j];

                t_y_xxy[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pb_xx[j] + 

                             pa_y[j] * pb_xxy[j]) * s_0_0[j];

                t_y_xxz[j] = (0.5 * pa_y[j] * fx[j] * pb_z[j] + pa_y[j] * pb_xxz[j]) * s_0_0[j];

                t_y_xyy[j] = (0.5 * pa_y[j] * pb_x[j] * fx[j] + fx[j] * pb_xy[j] + pa_y[j] * pb_xyy[j]) * s_0_0[j];

                t_y_xyz[j] = (0.5 * fx[j] * pb_xz[j] + pa_y[j] * pb_xyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (15,30)

            #pragma omp simd aligned(fx, pa_y, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_y_xzz, \
                                     t_y_yyy, t_y_yyz, t_y_yzz, t_y_zzz, t_z_xxx, t_z_xxy, t_z_xxz, t_z_xyy, t_z_xyz, \
                                     t_z_xzz, t_z_yyy, t_z_yyz, t_z_yzz, t_z_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_y_xzz[j] = (0.5 * pa_y[j] * pb_x[j] * fx[j] + pa_y[j] * pb_xzz[j]) * s_0_0[j];

                t_y_yyy[j] = (0.75 * fx[j] * fx[j] + 1.5 * pa_y[j] * pb_y[j] * fx[j] + 1.5 * fx[j] * pb_yy[j] + 

                             pa_y[j] * pb_yyy[j]) * s_0_0[j];

                t_y_yyz[j] = (0.5 * pa_y[j] * fx[j] * pb_z[j] + fx[j] * pb_yz[j] + pa_y[j] * pb_yyz[j]) * s_0_0[j];

                t_y_yzz[j] = (0.25 * fx[j] * fx[j] + 0.5 * pa_y[j] * pb_y[j] * fx[j] + 0.5 * fx[j] * pb_zz[j] + 

                             pa_y[j] * pb_yzz[j]) * s_0_0[j];

                t_y_zzz[j] = (1.5 * pa_y[j] * pb_z[j] * fx[j] + pa_y[j] * pb_zzz[j]) * s_0_0[j];

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

            // Batch of Integrals (0) = (0,15)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxy, pa_xxz, pa_xy, pa_xyy, pa_xyz, pa_xz, pa_y, pa_yy, \
                                     pa_yz, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xxx_x, t_xxx_y, t_xxx_z, t_xxy_x, t_xxy_y, \
                                     t_xxy_z, t_xxz_x, t_xxz_y, t_xxz_z, t_xyy_x, t_xyy_y, t_xyy_z, t_xyz_x, t_xyz_y, \
                                     t_xyz_z: VLX_ALIGN)
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

                t_xyy_y[j] = (pa_xy[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_y[j] + pa_xyy[j] * pb_y[j]) * s_0_0[j];

                t_xyy_z[j] = (0.5 * pa_x[j] * fx[j] * pb_z[j] + pa_xyy[j] * pb_z[j]) * s_0_0[j];

                t_xyz_x[j] = (0.5 * fx[j] * pa_yz[j] + pa_xyz[j] * pb_x[j]) * s_0_0[j];

                t_xyz_y[j] = (0.5 * pa_xz[j] * fx[j] + pa_xyz[j] * pb_y[j]) * s_0_0[j];

                t_xyz_z[j] = (0.5 * pa_xy[j] * fx[j] + pa_xyz[j] * pb_z[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (15,30)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_y, pa_yy, pa_yyy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pa_zzz, pb_x, pb_y, pb_z, s_0_0, t_xzz_x, t_xzz_y, t_xzz_z, t_yyy_x, t_yyy_y, t_yyy_z, \
                                     t_yyz_x, t_yyz_y, t_yyz_z, t_yzz_x, t_yzz_y, t_yzz_z, t_zzz_x, t_zzz_y, t_zzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_x[j] = (0.25 * fx[j] * fx[j] + 0.5 * fx[j] * pa_zz[j] + 0.5 * pa_x[j] * fx[j] * pb_x[j] + 

                             pa_xzz[j] * pb_x[j]) * s_0_0[j];

                t_xzz_y[j] = (0.5 * pa_x[j] * fx[j] * pb_y[j] + pa_xzz[j] * pb_y[j]) * s_0_0[j];

                t_xzz_z[j] = (pa_xz[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_z[j] + pa_xzz[j] * pb_z[j]) * s_0_0[j];

                t_yyy_x[j] = (1.5 * pa_y[j] * fx[j] * pb_x[j] + pa_yyy[j] * pb_x[j]) * s_0_0[j];

                t_yyy_y[j] = (0.75 * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] + 1.5 * pa_y[j] * fx[j] * pb_y[j] + 

                             pa_yyy[j] * pb_y[j]) * s_0_0[j];

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

            // Batch of Integrals (0) = (0,9)

            #pragma omp simd aligned(fx, pa_x, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, s_0_0, t_x_xxxx, \
                                     t_x_xxxy, t_x_xxxz, t_x_xxyy, t_x_xxyz, t_x_xxzz, t_x_xyyy, t_x_xyyz, t_x_xyzz: VLX_ALIGN)
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

                t_x_xxzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_x[j] + 

                              0.5 * pa_x[j] * pb_xx[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_zz[j] + fx[j] * pb_xzz[j] + pa_x[j] * pb_xxzz[j]) * s_0_0[j];

                t_x_xyyy[j] = (0.75 * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_x[j] * pb_xy[j] * fx[j] + 

                              0.5 * fx[j] * pb_yyy[j] + pa_x[j] * pb_xyyy[j]) * s_0_0[j];

                t_x_xyyz[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * pb_xz[j] * fx[j] + 

                              0.5 * fx[j] * pb_yyz[j] + pa_x[j] * pb_xyyz[j]) * s_0_0[j];

                t_x_xyzz[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * pb_xy[j] * fx[j] + 

                              0.5 * fx[j] * pb_yzz[j] + pa_x[j] * pb_xyzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (9,18)

            #pragma omp simd aligned(fx, pa_x, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xy, pb_xz, \
                                     pb_xzzz, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     s_0_0, t_x_xzzz, t_x_yyyy, t_x_yyyz, t_x_yyzz, t_x_yzzz, t_x_zzzz, t_y_xxxx, \
                                     t_y_xxxy, t_y_xxxz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_x_xzzz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_x[j] * pb_xz[j] * fx[j] + 

                              0.5 * fx[j] * pb_zzz[j] + pa_x[j] * pb_xzzz[j]) * s_0_0[j];

                t_x_yyyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 3.0 * pa_x[j] * pb_yy[j] * fx[j] + 

                              pa_x[j] * pb_yyyy[j]) * s_0_0[j];

                t_x_yyyz[j] = (1.5 * pa_x[j] * pb_yz[j] * fx[j] + pa_x[j] * pb_yyyz[j]) * s_0_0[j];

                t_x_yyzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_x[j] * pb_yy[j] * fx[j] + 

                              0.5 * pa_x[j] * fx[j] * pb_zz[j] + pa_x[j] * pb_yyzz[j]) * s_0_0[j];

                t_x_yzzz[j] = (1.5 * pa_x[j] * pb_yz[j] * fx[j] + pa_x[j] * pb_yzzz[j]) * s_0_0[j];

                t_x_zzzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 3.0 * pa_x[j] * pb_zz[j] * fx[j] + 

                              pa_x[j] * pb_zzzz[j]) * s_0_0[j];

                t_y_xxxx[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 3.0 * pa_y[j] * pb_xx[j] * fx[j] + 

                              pa_y[j] * pb_xxxx[j]) * s_0_0[j];

                t_y_xxxy[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_y[j] * pb_xy[j] * fx[j] + 

                              0.5 * fx[j] * pb_xxx[j] + pa_y[j] * pb_xxxy[j]) * s_0_0[j];

                t_y_xxxz[j] = (1.5 * pa_y[j] * pb_xz[j] * fx[j] + pa_y[j] * pb_xxxz[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (18,27)

            #pragma omp simd aligned(fx, pa_y, pb_x, pb_xx, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yz, pb_z, pb_zz, s_0_0, t_y_xxyy, t_y_xxyz, t_y_xxzz, t_y_xyyy, \
                                     t_y_xyyz, t_y_xyzz, t_y_xzzz, t_y_yyyy, t_y_yyyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_y_xxyy[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_y[j] + 

                              0.5 * pa_y[j] * pb_xx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_yy[j] + fx[j] * pb_xxy[j] + pa_y[j] * pb_xxyy[j]) * s_0_0[j];

                t_y_xxyz[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_y[j] * fx[j] * pb_yz[j] + 

                              0.5 * fx[j] * pb_xxz[j] + pa_y[j] * pb_xxyz[j]) * s_0_0[j];

                t_y_xxzz[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * pa_y[j] * pb_xx[j] * fx[j] + 

                              0.5 * pa_y[j] * fx[j] * pb_zz[j] + pa_y[j] * pb_xxzz[j]) * s_0_0[j];

                t_y_xyyy[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_y[j] * pb_xy[j] * fx[j] + 

                              1.5 * fx[j] * pb_xyy[j] + pa_y[j] * pb_xyyy[j]) * s_0_0[j];

                t_y_xyyz[j] = (0.5 * pa_y[j] * pb_xz[j] * fx[j] + fx[j] * pb_xyz[j] + pa_y[j] * pb_xyyz[j]) * s_0_0[j];

                t_y_xyzz[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_y[j] * pb_xy[j] * fx[j] + 

                              0.5 * fx[j] * pb_xzz[j] + pa_y[j] * pb_xyzz[j]) * s_0_0[j];

                t_y_xzzz[j] = (1.5 * pa_y[j] * pb_xz[j] * fx[j] + pa_y[j] * pb_xzzz[j]) * s_0_0[j];

                t_y_yyyy[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pb_y[j] + 

                              3.0 * pa_y[j] * pb_yy[j] * fx[j] + 2.0 * fx[j] * pb_yyy[j] + pa_y[j] * pb_yyyy[j]) * s_0_0[j];

                t_y_yyyz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_y[j] * pb_yz[j] * fx[j] + 

                              1.5 * fx[j] * pb_yyz[j] + pa_y[j] * pb_yyyz[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (27,36)

            #pragma omp simd aligned(fx, pa_y, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, \
                                     pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xz, pb_y, pb_yy, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, \
                                     pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_y_yyzz, t_y_yzzz, t_y_zzzz, t_z_xxxx, \
                                     t_z_xxxy, t_z_xxxz, t_z_xxyy, t_z_xxyz, t_z_xxzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_y_yyzz[j] = (0.25 * pa_y[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_y[j] + 

                              0.5 * pa_y[j] * pb_yy[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_zz[j] + fx[j] * pb_yzz[j] + pa_y[j] * pb_yyzz[j]) * s_0_0[j];

                t_y_yzzz[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_y[j] * pb_yz[j] * fx[j] + 

                              0.5 * fx[j] * pb_zzz[j] + pa_y[j] * pb_yzzz[j]) * s_0_0[j];

                t_y_zzzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 3.0 * pa_y[j] * pb_zz[j] * fx[j] + 

                              pa_y[j] * pb_zzzz[j]) * s_0_0[j];

                t_z_xxxx[j] = (0.75 * pa_z[j] * fx[j] * fx[j] + 3.0 * pa_z[j] * pb_xx[j] * fx[j] + 

                              pa_z[j] * pb_xxxx[j]) * s_0_0[j];

                t_z_xxxy[j] = (1.5 * pa_z[j] * pb_xy[j] * fx[j] + pa_z[j] * pb_xxxy[j]) * s_0_0[j];

                t_z_xxxz[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_z[j] * pb_xz[j] * fx[j] + 

                              0.5 * fx[j] * pb_xxx[j] + pa_z[j] * pb_xxxz[j]) * s_0_0[j];

                t_z_xxyy[j] = (0.25 * pa_z[j] * fx[j] * fx[j] + 0.5 * pa_z[j] * pb_xx[j] * fx[j] + 

                              0.5 * pa_z[j] * fx[j] * pb_yy[j] + pa_z[j] * pb_xxyy[j]) * s_0_0[j];

                t_z_xxyz[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_z[j] * fx[j] * pb_yz[j] + 

                              0.5 * fx[j] * pb_xxy[j] + pa_z[j] * pb_xxyz[j]) * s_0_0[j];

                t_z_xxzz[j] = (0.25 * pa_z[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pb_z[j] + 

                              0.5 * pa_z[j] * pb_xx[j] * fx[j] + 0.5 * pa_z[j] * fx[j] * pb_zz[j] + fx[j] * pb_xxz[j] + pa_z[j] * pb_xxzz[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (36,45)

            #pragma omp simd aligned(fx, pa_z, pb_x, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, \
                                     pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, \
                                     pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_z_xyyy, t_z_xyyz, t_z_xyzz, t_z_xzzz, \
                                     t_z_yyyy, t_z_yyyz, t_z_yyzz, t_z_yzzz, t_z_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_z_xyyy[j] = (1.5 * pa_z[j] * pb_xy[j] * fx[j] + pa_z[j] * pb_xyyy[j]) * s_0_0[j];

                t_z_xyyz[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_z[j] * pb_xz[j] * fx[j] + 

                              0.5 * fx[j] * pb_xyy[j] + pa_z[j] * pb_xyyz[j]) * s_0_0[j];

                t_z_xyzz[j] = (0.5 * pa_z[j] * pb_xy[j] * fx[j] + fx[j] * pb_xyz[j] + pa_z[j] * pb_xyzz[j]) * s_0_0[j];

                t_z_xzzz[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_z[j] * pb_xz[j] * fx[j] + 

                              1.5 * fx[j] * pb_xzz[j] + pa_z[j] * pb_xzzz[j]) * s_0_0[j];

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

            // Batch of Integrals (0) = (0,9)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pa_xxxy, pa_xxxz, pa_xxy, pa_xxz, pa_xy, pa_xz, \
                                     pa_y, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xxxx_x, t_xxxx_y, t_xxxx_z, t_xxxy_x, \
                                     t_xxxy_y, t_xxxy_z, t_xxxz_x, t_xxxz_y, t_xxxz_z: VLX_ALIGN)
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

                t_xxxy_z[j] = (1.5 * pa_xy[j] * fx[j] * pb_z[j] + pa_xxxy[j] * pb_z[j]) * s_0_0[j];

                t_xxxz_x[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 1.5 * pa_xxz[j] * fx[j] + 

                              1.5 * pa_xz[j] * fx[j] * pb_x[j] + pa_xxxz[j] * pb_x[j]) * s_0_0[j];

                t_xxxz_y[j] = (1.5 * pa_xz[j] * fx[j] * pb_y[j] + pa_xxxz[j] * pb_y[j]) * s_0_0[j];

                t_xxxz_z[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] + 

                              1.5 * pa_xz[j] * fx[j] * pb_z[j] + pa_xxxz[j] * pb_z[j]) * s_0_0[j];
            }

            // Batch of Integrals (1) = (9,18)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xxyz, pa_xxz, pa_xxzz, pa_xyy, pa_xyz, \
                                     pa_xzz, pa_y, pa_yy, pa_yz, pa_z, pa_zz, pb_x, pb_y, pb_z, s_0_0, t_xxyy_x, t_xxyy_y, \
                                     t_xxyy_z, t_xxyz_x, t_xxyz_y, t_xxyz_z, t_xxzz_x, t_xxzz_y, t_xxzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_x[j] = (0.5 * pa_x[j] * fx[j] * fx[j] + pa_xyy[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_x[j] + 

                              0.5 * pa_xx[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_yy[j] * pb_x[j] + pa_xxyy[j] * pb_x[j]) * s_0_0[j];

                t_xxyy_y[j] = (0.5 * fx[j] * fx[j] * pa_y[j] + pa_xxy[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_y[j] + 

                              0.5 * pa_xx[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_yy[j] * pb_y[j] + pa_xxyy[j] * pb_y[j]) * s_0_0[j];

                t_xxyy_z[j] = (0.25 * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_xx[j] * fx[j] * pb_z[j] + 

                              0.5 * fx[j] * pa_yy[j] * pb_z[j] + pa_xxyy[j] * pb_z[j]) * s_0_0[j];

                t_xxyz_x[j] = (pa_xyz[j] * fx[j] + 0.5 * fx[j] * pa_yz[j] * pb_x[j] + pa_xxyz[j] * pb_x[j]) * s_0_0[j];

                t_xxyz_y[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_xxz[j] * fx[j] + 

                              0.5 * fx[j] * pa_yz[j] * pb_y[j] + pa_xxyz[j] * pb_y[j]) * s_0_0[j];

                t_xxyz_z[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xxy[j] * fx[j] + 

                              0.5 * fx[j] * pa_yz[j] * pb_z[j] + pa_xxyz[j] * pb_z[j]) * s_0_0[j];

                t_xxzz_x[j] = (0.5 * pa_x[j] * fx[j] * fx[j] + pa_xzz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_x[j] + 

                              0.5 * pa_xx[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * pa_zz[j] * pb_x[j] + pa_xxzz[j] * pb_x[j]) * s_0_0[j];

                t_xxzz_y[j] = (0.25 * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xx[j] * fx[j] * pb_y[j] + 

                              0.5 * fx[j] * pa_zz[j] * pb_y[j] + pa_xxzz[j] * pb_y[j]) * s_0_0[j];

                t_xxzz_z[j] = (0.5 * fx[j] * fx[j] * pa_z[j] + pa_xxz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_z[j] + 

                              0.5 * pa_xx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_zz[j] * pb_z[j] + pa_xxzz[j] * pb_z[j]) * s_0_0[j];
            }

            // Batch of Integrals (2) = (18,27)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_xyyz, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, \
                                     pa_yyy, pa_yyz, pa_yzz, pa_z, pb_x, pb_y, pb_z, s_0_0, t_xyyy_x, t_xyyy_y, t_xyyy_z, \
                                     t_xyyz_x, t_xyyz_y, t_xyyz_z, t_xyzz_x, t_xyzz_y, t_xyzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_x[j] = (0.75 * fx[j] * fx[j] * pa_y[j] + 0.5 * fx[j] * pa_yyy[j] + 

                              1.5 * pa_xy[j] * fx[j] * pb_x[j] + pa_xyyy[j] * pb_x[j]) * s_0_0[j];

                t_xyyy_y[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 1.5 * pa_xyy[j] * fx[j] + 

                              1.5 * pa_xy[j] * fx[j] * pb_y[j] + pa_xyyy[j] * pb_y[j]) * s_0_0[j];

                t_xyyy_z[j] = (1.5 * pa_xy[j] * fx[j] * pb_z[j] + pa_xyyy[j] * pb_z[j]) * s_0_0[j];

                t_xyyz_x[j] = (0.25 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * pa_yyz[j] + 

                              0.5 * pa_xz[j] * fx[j] * pb_x[j] + pa_xyyz[j] * pb_x[j]) * s_0_0[j];

                t_xyyz_y[j] = (pa_xyz[j] * fx[j] + 0.5 * pa_xz[j] * fx[j] * pb_y[j] + pa_xyyz[j] * pb_y[j]) * s_0_0[j];

                t_xyyz_z[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] + 

                              0.5 * pa_xz[j] * fx[j] * pb_z[j] + pa_xyyz[j] * pb_z[j]) * s_0_0[j];

                t_xyzz_x[j] = (0.25 * fx[j] * fx[j] * pa_y[j] + 0.5 * fx[j] * pa_yzz[j] + 

                              0.5 * pa_xy[j] * fx[j] * pb_x[j] + pa_xyzz[j] * pb_x[j]) * s_0_0[j];

                t_xyzz_y[j] = (0.25 * pa_x[j] * fx[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] + 

                              0.5 * pa_xy[j] * fx[j] * pb_y[j] + pa_xyzz[j] * pb_y[j]) * s_0_0[j];

                t_xyzz_z[j] = (pa_xyz[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * pb_z[j] + pa_xyzz[j] * pb_z[j]) * s_0_0[j];
            }

            // Batch of Integrals (3) = (27,36)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_y, pa_yy, pa_yyy, pa_yyyy, pa_yyyz, pa_yyz, \
                                     pa_yz, pa_z, pa_zzz, pb_x, pb_y, pb_z, s_0_0, t_xzzz_x, t_xzzz_y, t_xzzz_z, t_yyyy_x, \
                                     t_yyyy_y, t_yyyy_z, t_yyyz_x, t_yyyz_y, t_yyyz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_x[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 0.5 * fx[j] * pa_zzz[j] + 

                              1.5 * pa_xz[j] * fx[j] * pb_x[j] + pa_xzzz[j] * pb_x[j]) * s_0_0[j];

                t_xzzz_y[j] = (1.5 * pa_xz[j] * fx[j] * pb_y[j] + pa_xzzz[j] * pb_y[j]) * s_0_0[j];

                t_xzzz_z[j] = (0.75 * pa_x[j] * fx[j] * fx[j] + 1.5 * pa_xzz[j] * fx[j] + 

                              1.5 * pa_xz[j] * fx[j] * pb_z[j] + pa_xzzz[j] * pb_z[j]) * s_0_0[j];

                t_yyyy_x[j] = (0.75 * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_yy[j] * fx[j] * pb_x[j] + 

                              pa_yyyy[j] * pb_x[j]) * s_0_0[j];

                t_yyyy_y[j] = (3.0 * pa_y[j] * fx[j] * fx[j] + 2.0 * pa_yyy[j] * fx[j] + 

                              0.75 * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_yy[j] * fx[j] * pb_y[j] + pa_yyyy[j] * pb_y[j]) * s_0_0[j];

                t_yyyy_z[j] = (0.75 * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_yy[j] * fx[j] * pb_z[j] + 

                              pa_yyyy[j] * pb_z[j]) * s_0_0[j];

                t_yyyz_x[j] = (1.5 * pa_yz[j] * fx[j] * pb_x[j] + pa_yyyz[j] * pb_x[j]) * s_0_0[j];

                t_yyyz_y[j] = (0.75 * fx[j] * fx[j] * pa_z[j] + 1.5 * pa_yyz[j] * fx[j] + 

                              1.5 * pa_yz[j] * fx[j] * pb_y[j] + pa_yyyz[j] * pb_y[j]) * s_0_0[j];

                t_yyyz_z[j] = (0.75 * pa_y[j] * fx[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] + 

                              1.5 * pa_yz[j] * fx[j] * pb_z[j] + pa_yyyz[j] * pb_z[j]) * s_0_0[j];
            }

            // Batch of Integrals (4) = (36,45)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, \
                                     pa_zzzz, pb_x, pb_y, pb_z, s_0_0, t_yyzz_x, t_yyzz_y, t_yyzz_z, t_yzzz_x, t_yzzz_y, \
                                     t_yzzz_z, t_zzzz_x, t_zzzz_y, t_zzzz_z: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_x[j] = (0.25 * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_yy[j] * fx[j] * pb_x[j] + 

                              0.5 * fx[j] * pa_zz[j] * pb_x[j] + pa_yyzz[j] * pb_x[j]) * s_0_0[j];

                t_yyzz_y[j] = (0.5 * pa_y[j] * fx[j] * fx[j] + pa_yzz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_y[j] + 

                              0.5 * pa_yy[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * pa_zz[j] * pb_y[j] + pa_yyzz[j] * pb_y[j]) * s_0_0[j];

                t_yyzz_z[j] = (0.5 * fx[j] * fx[j] * pa_z[j] + pa_yyz[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_z[j] + 

                              0.5 * pa_yy[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * pa_zz[j] * pb_z[j] + pa_yyzz[j] * pb_z[j]) * s_0_0[j];

                t_yzzz_x[j] = (1.5 * pa_yz[j] * fx[j] * pb_x[j] + pa_yzzz[j] * pb_x[j]) * s_0_0[j];

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


} // ovlrecfunc namespace

