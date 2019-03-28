//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForDX.hpp"

#include "KineticEnergyVecFuncForDX.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForDD(      CMemBlock2D<double>& primBuffer,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

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

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xy, pa_y, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_xx_xx, t_xx_xy, t_xx_xz, t_xx_yy, t_xx_yz, \
                                     t_xx_zz, t_xy_xx, t_xy_xy, t_xy_xz, t_xy_yy, t_xy_yz, t_xy_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_xx[j] = kinvecfunc::fvec_xx_xx_s_0(fx[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xx_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xx_xy[j] = kinvecfunc::fvec_xx_xy_s_0(fx[j], pa_x[j], pa_xx[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xx_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xx_xz[j] = kinvecfunc::fvec_xx_xz_s_0(fx[j], pa_x[j], pa_xx[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xx_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xx_yy[j] = kinvecfunc::fvec_xx_yy_s_0(fx[j], pa_xx[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xx_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_yy[j], r_0_0[j]);

                t_xx_yz[j] = kinvecfunc::fvec_xx_yz_s_0(fx[j], pa_xx[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xx_yz_r_0(fga[j], fx[j], fz[j], pa_xx[j], pb_yz[j], r_0_0[j]);

                t_xx_zz[j] = kinvecfunc::fvec_xx_zz_s_0(fx[j], pa_xx[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xx_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_zz[j], r_0_0[j]);

                t_xy_xx[j] = kinvecfunc::fvec_xy_xx_s_0(fx[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xy_xx_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xy_xy[j] = kinvecfunc::fvec_xy_xy_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xy_xy_r_0(fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xy_xz[j] = kinvecfunc::fvec_xy_xz_s_0(fx[j], pa_xy[j], pa_y[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xy_xz_r_0(fx[j], fz[j], pa_xy[j], pa_y[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xy_yy[j] = kinvecfunc::fvec_xy_yy_s_0(fx[j], pa_x[j], pa_xy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xy_yy_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xy_yz[j] = kinvecfunc::fvec_xy_yz_s_0(fx[j], pa_x[j], pa_xy[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xy_yz_r_0(fx[j], fz[j], pa_x[j], pa_xy[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xy_zz[j] = kinvecfunc::fvec_xy_zz_s_0(fx[j], pa_xy[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xy_zz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pb_zz[j], r_0_0[j]);
            }

            // Batch of Integrals (1) = (12,24)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_y, pa_yy, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_xz_xx, t_xz_xy, t_xz_xz, t_xz_yy, t_xz_yz, \
                                     t_xz_zz, t_yy_xx, t_yy_xy, t_yy_xz, t_yy_yy, t_yy_yz, t_yy_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_xx[j] = kinvecfunc::fvec_xz_xx_s_0(fx[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xz_xx_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xz_xy[j] = kinvecfunc::fvec_xz_xy_s_0(fx[j], pa_xz[j], pa_z[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xz_xy_r_0(fx[j], fz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xz_xz[j] = kinvecfunc::fvec_xz_xz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xz_xz_r_0(fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xz_yy[j] = kinvecfunc::fvec_xz_yy_s_0(fx[j], pa_xz[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xz_yy_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pb_yy[j], r_0_0[j]);

                t_xz_yz[j] = kinvecfunc::fvec_xz_yz_s_0(fx[j], pa_x[j], pa_xz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xz_yz_r_0(fx[j], fz[j], pa_x[j], pa_xz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xz_zz[j] = kinvecfunc::fvec_xz_zz_s_0(fx[j], pa_x[j], pa_xz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xz_zz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yy_xx[j] = kinvecfunc::fvec_yy_xx_s_0(fx[j], pa_yy[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_yy_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_xx[j], r_0_0[j]);

                t_yy_xy[j] = kinvecfunc::fvec_yy_xy_s_0(fx[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yy_xy_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], r_0_0[j]);

                t_yy_xz[j] = kinvecfunc::fvec_yy_xz_s_0(fx[j], pa_yy[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yy_xz_r_0(fga[j], fx[j], fz[j], pa_yy[j], pb_xz[j], r_0_0[j]);

                t_yy_yy[j] = kinvecfunc::fvec_yy_yy_s_0(fx[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yy_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yy_yz[j] = kinvecfunc::fvec_yy_yz_s_0(fx[j], pa_y[j], pa_yy[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yy_yz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yy_zz[j] = kinvecfunc::fvec_yy_zz_s_0(fx[j], pa_yy[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yy_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_zz[j], r_0_0[j]);
            }

            // Batch of Integrals (2) = (24,36)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_yz_xx, t_yz_xy, t_yz_xz, t_yz_yy, t_yz_yz, \
                                     t_yz_zz, t_zz_xx, t_zz_xy, t_zz_xz, t_zz_yy, t_zz_yz, t_zz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_xx[j] = kinvecfunc::fvec_yz_xx_s_0(fx[j], pa_yz[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_yz_xx_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pb_xx[j], r_0_0[j]);

                t_yz_xy[j] = kinvecfunc::fvec_yz_xy_s_0(fx[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yz_xy_r_0(fx[j], fz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], r_0_0[j]);

                t_yz_xz[j] = kinvecfunc::fvec_yz_xz_s_0(fx[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yz_xz_r_0(fx[j], fz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], r_0_0[j]);

                t_yz_yy[j] = kinvecfunc::fvec_yz_yy_s_0(fx[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yz_yy_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yz_yz[j] = kinvecfunc::fvec_yz_yz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yz_yz_r_0(fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yz_zz[j] = kinvecfunc::fvec_yz_zz_s_0(fx[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yz_zz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_zz_xx[j] = kinvecfunc::fvec_zz_xx_s_0(fx[j], pa_zz[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_zz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_xx[j], r_0_0[j]);

                t_zz_xy[j] = kinvecfunc::fvec_zz_xy_s_0(fx[j], pa_zz[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_zz_xy_r_0(fga[j], fx[j], fz[j], pa_zz[j], pb_xy[j], r_0_0[j]);

                t_zz_xz[j] = kinvecfunc::fvec_zz_xz_s_0(fx[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_zz_xz_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], r_0_0[j]);

                t_zz_yy[j] = kinvecfunc::fvec_zz_yy_s_0(fx[j], pa_zz[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_zz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_yy[j], r_0_0[j]);

                t_zz_yz[j] = kinvecfunc::fvec_zz_yz_s_0(fx[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_zz_yz_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_zz_zz[j] = kinvecfunc::fvec_zz_zz_s_0(fx[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_zz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], r_0_0[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForDF(      CMemBlock2D<double>& primBuffer,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

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

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_xx_xxx, t_xx_xxy, t_xx_xxz, t_xx_xyy, t_xx_xyz, t_xx_xzz, \
                                     t_xx_yyy, t_xx_yyz, t_xx_yzz, t_xx_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_xxx[j] = kinvecfunc::fvec_xx_xxx_s_0(fx[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xx_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xx_xxy[j] = kinvecfunc::fvec_xx_xxy_s_0(fx[j], pa_x[j], pa_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xx_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xx_xxz[j] = kinvecfunc::fvec_xx_xxz_s_0(fx[j], pa_x[j], pa_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xx_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xx_xyy[j] = kinvecfunc::fvec_xx_xyy_s_0(fx[j], pa_x[j], pa_xx[j], pb_x[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xx_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_x[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xx_xyz[j] = kinvecfunc::fvec_xx_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pb_xyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xx_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xyz[j], pb_yz[j], r_0_0[j]);

                t_xx_xzz[j] = kinvecfunc::fvec_xx_xzz_s_0(fx[j], pa_x[j], pa_xx[j], pb_x[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xx_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_x[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xx_yyy[j] = kinvecfunc::fvec_xx_yyy_s_0(fx[j], pa_xx[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xx_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xx_yyz[j] = kinvecfunc::fvec_xx_yyz_s_0(fx[j], pa_xx[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xx_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xx_yzz[j] = kinvecfunc::fvec_xx_yzz_s_0(fx[j], pa_xx[j], pb_y[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xx_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_y[j], pb_yzz[j], r_0_0[j]);

                t_xx_zzz[j] = kinvecfunc::fvec_xx_zzz_s_0(fx[j], pa_xx[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xx_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_xy_xxx, t_xy_xxy, t_xy_xxz, t_xy_xyy, t_xy_xyz, t_xy_xzz, \
                                     t_xy_yyy, t_xy_yyz, t_xy_yzz, t_xy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xy_xxx[j] = kinvecfunc::fvec_xy_xxx_s_0(fx[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xy_xxx_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xy_xxy[j] = kinvecfunc::fvec_xy_xxy_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xy_xxy_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xy_xxz[j] = kinvecfunc::fvec_xy_xxz_s_0(fx[j], pa_xy[j], pa_y[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xy_xxz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_y[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xy_xyy[j] = kinvecfunc::fvec_xy_xyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xy_xyy_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xy_xyz[j] = kinvecfunc::fvec_xy_xyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xy_xyz_r_0(fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xy_xzz[j] = kinvecfunc::fvec_xy_xzz_s_0(fx[j], pa_xy[j], pa_y[j], pb_x[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xy_xzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_y[j], pb_x[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xy_yyy[j] = kinvecfunc::fvec_xy_yyy_s_0(fx[j], pa_x[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xy_yyy_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xy_yyz[j] = kinvecfunc::fvec_xy_yyz_s_0(fx[j], pa_x[j], pa_xy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xy_yyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xy_yzz[j] = kinvecfunc::fvec_xy_yzz_s_0(fx[j], pa_x[j], pa_xy[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xy_yzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xy_zzz[j] = kinvecfunc::fvec_xy_zzz_s_0(fx[j], pa_xy[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xy_zzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_xz_xxx, t_xz_xxy, t_xz_xxz, t_xz_xyy, t_xz_xyz, t_xz_xzz, \
                                     t_xz_yyy, t_xz_yyz, t_xz_yzz, t_xz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_xxx[j] = kinvecfunc::fvec_xz_xxx_s_0(fx[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xz_xxx_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xz_xxy[j] = kinvecfunc::fvec_xz_xxy_s_0(fx[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xz_xxy_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xz_xxz[j] = kinvecfunc::fvec_xz_xxz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xz_xxz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xz_xyy[j] = kinvecfunc::fvec_xz_xyy_s_0(fx[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xz_xyy_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xz_xyz[j] = kinvecfunc::fvec_xz_xyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xz_xyz_r_0(fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xz_xzz[j] = kinvecfunc::fvec_xz_xzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xz_xzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xz_yyy[j] = kinvecfunc::fvec_xz_yyy_s_0(fx[j], pa_xz[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xz_yyy_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xz_yyz[j] = kinvecfunc::fvec_xz_yyz_s_0(fx[j], pa_x[j], pa_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xz_yyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xz_yzz[j] = kinvecfunc::fvec_xz_yzz_s_0(fx[j], pa_x[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xz_yzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_xz_zzz[j] = kinvecfunc::fvec_xz_zzz_s_0(fx[j], pa_x[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xz_zzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (3) = (30,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_yy_xxx, t_yy_xxy, t_yy_xxz, t_yy_xyy, t_yy_xyz, t_yy_xzz, \
                                     t_yy_yyy, t_yy_yyz, t_yy_yzz, t_yy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yy_xxx[j] = kinvecfunc::fvec_yy_xxx_s_0(fx[j], pa_yy[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_yy_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_yy_xxy[j] = kinvecfunc::fvec_yy_xxy_s_0(fx[j], pa_y[j], pa_yy[j], pb_xx[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_yy_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_xx[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_yy_xxz[j] = kinvecfunc::fvec_yy_xxz_s_0(fx[j], pa_yy[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yy_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_yy_xyy[j] = kinvecfunc::fvec_yy_xyy_s_0(fx[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_yy_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], r_0_0[j]);

                t_yy_xyz[j] = kinvecfunc::fvec_yy_xyz_s_0(fx[j], pa_y[j], pa_yy[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yy_xyz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yy_xzz[j] = kinvecfunc::fvec_yy_xzz_s_0(fx[j], pa_yy[j], pb_x[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yy_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_x[j], pb_xzz[j], r_0_0[j]);

                t_yy_yyy[j] = kinvecfunc::fvec_yy_yyy_s_0(fx[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_yy_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_yy_yyz[j] = kinvecfunc::fvec_yy_yyz_s_0(fx[j], pa_y[j], pa_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yy_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yy_yzz[j] = kinvecfunc::fvec_yy_yzz_s_0(fx[j], pa_y[j], pa_yy[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yy_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_yy_zzz[j] = kinvecfunc::fvec_yy_zzz_s_0(fx[j], pa_yy[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yy_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (4) = (40,50)

            #pragma omp simd aligned(fgb, fx, fz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_yz_xxx, t_yz_xxy, t_yz_xxz, t_yz_xyy, t_yz_xyz, t_yz_xzz, \
                                     t_yz_yyy, t_yz_yyz, t_yz_yzz, t_yz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_xxx[j] = kinvecfunc::fvec_yz_xxx_s_0(fx[j], pa_yz[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_yz_xxx_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_yz_xxy[j] = kinvecfunc::fvec_yz_xxy_s_0(fx[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_yz_xxy_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_yz_xxz[j] = kinvecfunc::fvec_yz_xxz_s_0(fx[j], pa_y[j], pa_yz[j], pb_xx[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yz_xxz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pb_xx[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_yz_xyy[j] = kinvecfunc::fvec_yz_xyy_s_0(fx[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_yz_xyy_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], r_0_0[j]);

                t_yz_xyz[j] = kinvecfunc::fvec_yz_xyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yz_xyz_r_0(fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yz_xzz[j] = kinvecfunc::fvec_yz_xzz_s_0(fx[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yz_xzz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_yz_yyy[j] = kinvecfunc::fvec_yz_yyy_s_0(fx[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_yz_yyy_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_yz_yyz[j] = kinvecfunc::fvec_yz_yyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yz_yyz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yz_yzz[j] = kinvecfunc::fvec_yz_yzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yz_yzz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yz_zzz[j] = kinvecfunc::fvec_yz_zzz_s_0(fx[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yz_zzz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (5) = (50,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, \
                                     pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_zz_xxx, t_zz_xxy, t_zz_xxz, t_zz_xyy, t_zz_xyz, t_zz_xzz, \
                                     t_zz_yyy, t_zz_yyz, t_zz_yzz, t_zz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zz_xxx[j] = kinvecfunc::fvec_zz_xxx_s_0(fx[j], pa_zz[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_zz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_zz_xxy[j] = kinvecfunc::fvec_zz_xxy_s_0(fx[j], pa_zz[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_zz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_zz_xxz[j] = kinvecfunc::fvec_zz_xxz_s_0(fx[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_zz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_zz_xyy[j] = kinvecfunc::fvec_zz_xyy_s_0(fx[j], pa_zz[j], pb_x[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_zz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_x[j], pb_xyy[j], r_0_0[j]);

                t_zz_xyz[j] = kinvecfunc::fvec_zz_xyz_s_0(fx[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], s_0_0[j]) + kinvecfunc::fvec_zz_xyz_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], r_0_0[j]);

                t_zz_xzz[j] = kinvecfunc::fvec_zz_xzz_s_0(fx[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_zz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_zz_yyy[j] = kinvecfunc::fvec_zz_yyy_s_0(fx[j], pa_zz[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_zz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_zz_yyz[j] = kinvecfunc::fvec_zz_yyz_s_0(fx[j], pa_z[j], pa_zz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_zz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_zz_yzz[j] = kinvecfunc::fvec_zz_yzz_s_0(fx[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_zz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_zz_zzz[j] = kinvecfunc::fvec_zz_zzz_s_0(fx[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_zz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForFD(      CMemBlock2D<double>& primBuffer,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

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

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

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

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xy, pb_xz, \
                                     pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_xxx_xx, t_xxx_xy, t_xxx_xz, \
                                     t_xxx_yy, t_xxx_yz, t_xxx_zz, t_xxy_xx, t_xxy_xy, t_xxy_xz, t_xxy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xx[j] = kinvecfunc::fvec_xxx_xx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xxx_xy[j] = kinvecfunc::fvec_xxx_xy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxx_xz[j] = kinvecfunc::fvec_xxx_xz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxx_yy[j] = kinvecfunc::fvec_xxx_yy_s_0(fx[j], pa_x[j], pa_xxx[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxx_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_yy[j], r_0_0[j]);

                t_xxx_yz[j] = kinvecfunc::fvec_xxx_yz_s_0(fx[j], pa_x[j], pa_xxx[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_yz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_yz[j], r_0_0[j]);

                t_xxx_zz[j] = kinvecfunc::fvec_xxx_zz_s_0(fx[j], pa_x[j], pa_xxx[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_zz[j], r_0_0[j]);

                t_xxy_xx[j] = kinvecfunc::fvec_xxy_xx_s_0(fx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xxy_xy[j] = kinvecfunc::fvec_xxy_xy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxy_xz[j] = kinvecfunc::fvec_xxy_xz_s_0(fx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xz_r_0(fga[j], fx[j], fz[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxy_yy[j] = kinvecfunc::fvec_xxy_yy_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxy_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yy[j], r_0_0[j]);
            }

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xxz, pa_xy, pa_xyy, pa_xz, pa_y, pa_yy, pa_z, \
                                     pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_xxy_yz, \
                                     t_xxy_zz, t_xxz_xx, t_xxz_xy, t_xxz_xz, t_xxz_yy, t_xxz_yz, t_xxz_zz, t_xyy_xx, \
                                     t_xyy_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_yz[j] = kinvecfunc::fvec_xxy_yz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxy_yz_r_0(fga[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxy_zz[j] = kinvecfunc::fvec_xxy_zz_s_0(fx[j], pa_xxy[j], pa_y[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_y[j], pb_zz[j], r_0_0[j]);

                t_xxz_xx[j] = kinvecfunc::fvec_xxz_xx_s_0(fx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xxz_xy[j] = kinvecfunc::fvec_xxz_xy_s_0(fx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xy_r_0(fga[j], fx[j], fz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxz_xz[j] = kinvecfunc::fvec_xxz_xz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxz_yy[j] = kinvecfunc::fvec_xxz_yy_s_0(fx[j], pa_xxz[j], pa_z[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_z[j], pb_yy[j], r_0_0[j]);

                t_xxz_yz[j] = kinvecfunc::fvec_xxz_yz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_yz_r_0(fga[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xxz_zz[j] = kinvecfunc::fvec_xxz_zz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyy_xx[j] = kinvecfunc::fvec_xyy_xx_s_0(fx[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xyy_xy[j] = kinvecfunc::fvec_xyy_xy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);
            }

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyz, pa_xz, pa_y, pa_yy, pa_yz, pa_z, pb_x, \
                                     pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_xyy_xz, t_xyy_yy, \
                                     t_xyy_yz, t_xyy_zz, t_xyz_xx, t_xyz_xy, t_xyz_xz, t_xyz_yy, t_xyz_yz, t_xyz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xz[j] = kinvecfunc::fvec_xyy_xz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyy_yy[j] = kinvecfunc::fvec_xyy_yy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyy_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyy_yz[j] = kinvecfunc::fvec_xyy_yz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyy_yz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyy_zz[j] = kinvecfunc::fvec_xyy_zz_s_0(fx[j], pa_x[j], pa_xyy[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pb_zz[j], r_0_0[j]);

                t_xyz_xx[j] = kinvecfunc::fvec_xyz_xx_s_0(fx[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xx_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xyz_xy[j] = kinvecfunc::fvec_xyz_xy_s_0(fx[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xy_r_0(fx[j], fz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyz_xz[j] = kinvecfunc::fvec_xyz_xz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xz_r_0(fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyz_yy[j] = kinvecfunc::fvec_xyz_yy_s_0(fx[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyz_yy_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyz_yz[j] = kinvecfunc::fvec_xyz_yz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyz_yz_r_0(fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyz_zz[j] = kinvecfunc::fvec_xyz_zz_s_0(fx[j], pa_xy[j], pa_xyz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_zz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pb_z[j], pb_zz[j], r_0_0[j]);
            }

            // Batch of Integrals (3) = (30,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_y, pa_yy, pa_yyy, pa_z, pa_zz, pb_x, pb_xx, \
                                     pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_xzz_xx, t_xzz_xy, \
                                     t_xzz_xz, t_xzz_yy, t_xzz_yz, t_xzz_zz, t_yyy_xx, t_yyy_xy, t_yyy_xz, t_yyy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xx[j] = kinvecfunc::fvec_xzz_xx_s_0(fx[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xzz_xy[j] = kinvecfunc::fvec_xzz_xy_s_0(fx[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xzz_xz[j] = kinvecfunc::fvec_xzz_xz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xzz_yy[j] = kinvecfunc::fvec_xzz_yy_s_0(fx[j], pa_x[j], pa_xzz[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xzz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pb_yy[j], r_0_0[j]);

                t_xzz_yz[j] = kinvecfunc::fvec_xzz_yz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_yz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xzz_zz[j] = kinvecfunc::fvec_xzz_zz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yyy_xx[j] = kinvecfunc::fvec_yyy_xx_s_0(fx[j], pa_y[j], pa_yyy[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_xx[j], r_0_0[j]);

                t_yyy_xy[j] = kinvecfunc::fvec_yyy_xy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xy_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], r_0_0[j]);

                t_yyy_xz[j] = kinvecfunc::fvec_yyy_xz_s_0(fx[j], pa_y[j], pa_yyy[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_xz[j], r_0_0[j]);

                t_yyy_yy[j] = kinvecfunc::fvec_yyy_yy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yyy_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yy[j], r_0_0[j]);
            }

            // Batch of Integrals (4) = (40,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, \
                                     pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_yyy_yz, t_yyy_zz, \
                                     t_yyz_xx, t_yyz_xy, t_yyz_xz, t_yyz_yy, t_yyz_yz, t_yyz_zz, t_yzz_xx, t_yzz_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_yz[j] = kinvecfunc::fvec_yyy_yz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyy_yz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyy_zz[j] = kinvecfunc::fvec_yyy_zz_s_0(fx[j], pa_y[j], pa_yyy[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_zz[j], r_0_0[j]);

                t_yyz_xx[j] = kinvecfunc::fvec_yyz_xx_s_0(fx[j], pa_yyz[j], pa_z[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_z[j], pb_xx[j], r_0_0[j]);

                t_yyz_xy[j] = kinvecfunc::fvec_yyz_xy_s_0(fx[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xy_r_0(fga[j], fx[j], fz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], r_0_0[j]);

                t_yyz_xz[j] = kinvecfunc::fvec_yyz_xz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xz_r_0(fga[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], r_0_0[j]);

                t_yyz_yy[j] = kinvecfunc::fvec_yyz_yy_s_0(fx[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yyz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yyz_yz[j] = kinvecfunc::fvec_yyz_yz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyz_yz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyz_zz[j] = kinvecfunc::fvec_yyz_zz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yzz_xx[j] = kinvecfunc::fvec_yzz_xx_s_0(fx[j], pa_y[j], pa_yzz[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pb_xx[j], r_0_0[j]);

                t_yzz_xy[j] = kinvecfunc::fvec_yzz_xy_s_0(fx[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xy_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], r_0_0[j]);
            }

            // Batch of Integrals (5) = (50,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xy, pb_xz, \
                                     pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_yzz_xz, t_yzz_yy, t_yzz_yz, \
                                     t_yzz_zz, t_zzz_xx, t_zzz_xy, t_zzz_xz, t_zzz_yy, t_zzz_yz, t_zzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xz[j] = kinvecfunc::fvec_yzz_xz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], r_0_0[j]);

                t_yzz_yy[j] = kinvecfunc::fvec_yzz_yy_s_0(fx[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yzz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yzz_yz[j] = kinvecfunc::fvec_yzz_yz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yzz_yz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yzz_zz[j] = kinvecfunc::fvec_yzz_zz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_zzz_xx[j] = kinvecfunc::fvec_zzz_xx_s_0(fx[j], pa_z[j], pa_zzz[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_xx[j], r_0_0[j]);

                t_zzz_xy[j] = kinvecfunc::fvec_zzz_xy_s_0(fx[j], pa_z[j], pa_zzz[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xy_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_xy[j], r_0_0[j]);

                t_zzz_xz[j] = kinvecfunc::fvec_zzz_xz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xz_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], r_0_0[j]);

                t_zzz_yy[j] = kinvecfunc::fvec_zzz_yy_s_0(fx[j], pa_z[j], pa_zzz[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_zzz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_yy[j], r_0_0[j]);

                t_zzz_yz[j] = kinvecfunc::fvec_zzz_yz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_yz_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_zzz_zz[j] = kinvecfunc::fvec_zzz_zz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_z[j], pb_zz[j], r_0_0[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForDG(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        kinrecfunc::compKineticEnergyForDG_0_4(primBuffer, auxBuffer, osFactors, paDistances, pbDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        kinrecfunc::compKineticEnergyForDG_5_8(primBuffer, auxBuffer, osFactors, paDistances, pbDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
    }

    void
    compKineticEnergyForDG_0_4(      CMemBlock2D<double>& primBuffer,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(9 * idx);

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xx = paDistances.data(9 * idx + 3);

            auto pa_xy = paDistances.data(9 * idx + 4);

            auto pa_xz = paDistances.data(9 * idx + 5);

            auto pa_yy = paDistances.data(9 * idx + 6);

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

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

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

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, r_0_0, s_0_0, t_xx_xxxx, t_xx_xxxy, t_xx_xxxz, t_xx_xxyy, t_xx_xxyz, \
                                     t_xx_xxzz, t_xx_xyyy, t_xx_xyyz, t_xx_xyzz, t_xx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_xxxx[j] = kinvecfunc::fvec_xx_xxxx_s_0(fx[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_xx_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], r_0_0[j]);

                t_xx_xxxy[j] = kinvecfunc::fvec_xx_xxxy_s_0(fx[j], pa_x[j], pa_xx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xx_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xx_xxxz[j] = kinvecfunc::fvec_xx_xxxz_s_0(fx[j], pa_x[j], pa_xx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xx_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xx_xxyy[j] = kinvecfunc::fvec_xx_xxyy_s_0(fx[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xx_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xx_xxyz[j] = kinvecfunc::fvec_xx_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pb_xxyz[j], pb_xyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xx_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xxyz[j], pb_xyz[j], pb_yz[j], r_0_0[j]);

                t_xx_xxzz[j] = kinvecfunc::fvec_xx_xxzz_s_0(fx[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xx_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xx_xyyy[j] = kinvecfunc::fvec_xx_xyyy_s_0(fx[j], pa_x[j], pa_xx[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xx_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xx_xyyz[j] = kinvecfunc::fvec_xx_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pb_xyyz[j], pb_xz[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xx_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xyyz[j], pb_xz[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xx_xyzz[j] = kinvecfunc::fvec_xx_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pb_xy[j], pb_xyzz[j], pb_y[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xx_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xy[j], pb_xyzz[j], pb_y[j], pb_yzz[j], r_0_0[j]);

                t_xx_xzzz[j] = kinvecfunc::fvec_xx_xzzz_s_0(fx[j], pa_x[j], pa_xx[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xx_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, \
                                     pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_z, pb_zz, pb_zzzz, r_0_0, s_0_0, \
                                     t_xx_yyyy, t_xx_yyyz, t_xx_yyzz, t_xx_yzzz, t_xx_zzzz, t_xy_xxxx, t_xy_xxxy, \
                                     t_xy_xxxz, t_xy_xxyy, t_xy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xx_yyyy[j] = kinvecfunc::fvec_xx_yyyy_s_0(fx[j], pa_xx[j], pb_yy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_xx_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_yy[j], pb_yyyy[j], r_0_0[j]);

                t_xx_yyyz[j] = kinvecfunc::fvec_xx_yyyz_s_0(fx[j], pa_xx[j], pb_yyyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xx_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_yyyz[j], pb_yz[j], r_0_0[j]);

                t_xx_yyzz[j] = kinvecfunc::fvec_xx_yyzz_s_0(fx[j], pa_xx[j], pb_yy[j], pb_yyzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xx_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_yy[j], pb_yyzz[j], pb_zz[j], r_0_0[j]);

                t_xx_yzzz[j] = kinvecfunc::fvec_xx_yzzz_s_0(fx[j], pa_xx[j], pb_yz[j], pb_yzzz[j], s_0_0[j]) + kinvecfunc::fvec_xx_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_yz[j], pb_yzzz[j], r_0_0[j]);

                t_xx_zzzz[j] = kinvecfunc::fvec_xx_zzzz_s_0(fx[j], pa_xx[j], pb_zz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_xx_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pb_zz[j], pb_zzzz[j], r_0_0[j]);

                t_xy_xxxx[j] = kinvecfunc::fvec_xy_xxxx_s_0(fx[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_xy_xxxx_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], r_0_0[j]);

                t_xy_xxxy[j] = kinvecfunc::fvec_xy_xxxy_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xy_xxxy_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xy_xxxz[j] = kinvecfunc::fvec_xy_xxxz_s_0(fx[j], pa_xy[j], pa_y[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xy_xxxz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_y[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xy_xxyy[j] = kinvecfunc::fvec_xy_xxyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xy_xxyy_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xy_xxyz[j] = kinvecfunc::fvec_xy_xxyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_xxyz[j], pb_xxz[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xy_xxyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_xxyz[j], pb_xxz[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, \
                                     t_xy_xxzz, t_xy_xyyy, t_xy_xyyz, t_xy_xyzz, t_xy_xzzz, t_xy_yyyy, t_xy_yyyz, \
                                     t_xy_yyzz, t_xy_yzzz, t_xy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xy_xxzz[j] = kinvecfunc::fvec_xy_xxzz_s_0(fx[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xy_xxzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xy_xyyy[j] = kinvecfunc::fvec_xy_xyyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xy_xyyy_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xy_xyyz[j] = kinvecfunc::fvec_xy_xyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xy_xyyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xy_xyzz[j] = kinvecfunc::fvec_xy_xyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xy_xyzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xy_xzzz[j] = kinvecfunc::fvec_xy_xzzz_s_0(fx[j], pa_xy[j], pa_y[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xy_xzzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_y[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);

                t_xy_yyyy[j] = kinvecfunc::fvec_xy_yyyy_s_0(fx[j], pa_x[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_xy_yyyy_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], r_0_0[j]);

                t_xy_yyyz[j] = kinvecfunc::fvec_xy_yyyz_s_0(fx[j], pa_x[j], pa_xy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xy_yyyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xy_yyzz[j] = kinvecfunc::fvec_xy_yyzz_s_0(fx[j], pa_x[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xy_yyzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xy_yzzz[j] = kinvecfunc::fvec_xy_yzzz_s_0(fx[j], pa_x[j], pa_xy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xy_yzzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);

                t_xy_zzzz[j] = kinvecfunc::fvec_xy_zzzz_s_0(fx[j], pa_xy[j], pb_zz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_xy_zzzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pb_zz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (3) = (30,40)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, r_0_0, s_0_0, t_xz_xxxx, t_xz_xxxy, t_xz_xxxz, t_xz_xxyy, t_xz_xxyz, \
                                     t_xz_xxzz, t_xz_xyyy, t_xz_xyyz, t_xz_xyzz, t_xz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_xxxx[j] = kinvecfunc::fvec_xz_xxxx_s_0(fx[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_xz_xxxx_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], r_0_0[j]);

                t_xz_xxxy[j] = kinvecfunc::fvec_xz_xxxy_s_0(fx[j], pa_xz[j], pa_z[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xz_xxxy_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pa_z[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xz_xxxz[j] = kinvecfunc::fvec_xz_xxxz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xz_xxxz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xz_xxyy[j] = kinvecfunc::fvec_xz_xxyy_s_0(fx[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xz_xxyy_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xz_xxyz[j] = kinvecfunc::fvec_xz_xxyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xxyz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xz_xxyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xxyz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xz_xxzz[j] = kinvecfunc::fvec_xz_xxzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xz_xxzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xz_xyyy[j] = kinvecfunc::fvec_xz_xyyy_s_0(fx[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xz_xyyy_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xz_xyyz[j] = kinvecfunc::fvec_xz_xyyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xz_xyyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xz_xyzz[j] = kinvecfunc::fvec_xz_xyzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xz_xyzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_xz_xzzz[j] = kinvecfunc::fvec_xz_xzzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xz_xzzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (4) = (40,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, \
                                     s_0_0, t_xz_yyyy, t_xz_yyyz, t_xz_yyzz, t_xz_yzzz, t_xz_zzzz, t_yy_xxxx, \
                                     t_yy_xxxy, t_yy_xxxz, t_yy_xxyy, t_yy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xz_yyyy[j] = kinvecfunc::fvec_xz_yyyy_s_0(fx[j], pa_xz[j], pb_yy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_xz_yyyy_r_0(fgb[j], fx[j], fz[j], pa_xz[j], pb_yy[j], pb_yyyy[j], r_0_0[j]);

                t_xz_yyyz[j] = kinvecfunc::fvec_xz_yyyz_s_0(fx[j], pa_x[j], pa_xz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xz_yyyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], r_0_0[j]);

                t_xz_yyzz[j] = kinvecfunc::fvec_xz_yyzz_s_0(fx[j], pa_x[j], pa_xz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xz_yyzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xz_yzzz[j] = kinvecfunc::fvec_xz_yzzz_s_0(fx[j], pa_x[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], s_0_0[j]) + kinvecfunc::fvec_xz_yzzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], r_0_0[j]);

                t_xz_zzzz[j] = kinvecfunc::fvec_xz_zzzz_s_0(fx[j], pa_x[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_xz_zzzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], r_0_0[j]);

                t_yy_xxxx[j] = kinvecfunc::fvec_yy_xxxx_s_0(fx[j], pa_yy[j], pb_xx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_yy_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_xx[j], pb_xxxx[j], r_0_0[j]);

                t_yy_xxxy[j] = kinvecfunc::fvec_yy_xxxy_s_0(fx[j], pa_y[j], pa_yy[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yy_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], r_0_0[j]);

                t_yy_xxxz[j] = kinvecfunc::fvec_yy_xxxz_s_0(fx[j], pa_yy[j], pb_xxxz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yy_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_xxxz[j], pb_xz[j], r_0_0[j]);

                t_yy_xxyy[j] = kinvecfunc::fvec_yy_xxyy_s_0(fx[j], pa_y[j], pa_yy[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yy_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yy_xxyz[j] = kinvecfunc::fvec_yy_xxyz_s_0(fx[j], pa_y[j], pa_yy[j], pb_xxyz[j], pb_xxz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yy_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_xxyz[j], pb_xxz[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForDG_5_8(      CMemBlock2D<double>& primBuffer,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_y = paDistances.data(9 * idx + 1);

            auto pa_z = paDistances.data(9 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

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

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

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

            // Batch of Integrals (5) = (50,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, \
                                     pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_yy_xxzz, \
                                     t_yy_xyyy, t_yy_xyyz, t_yy_xyzz, t_yy_xzzz, t_yy_yyyy, t_yy_yyyz, t_yy_yyzz, \
                                     t_yy_yzzz, t_yy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yy_xxzz[j] = kinvecfunc::fvec_yy_xxzz_s_0(fx[j], pa_yy[j], pb_xx[j], pb_xxzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yy_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_xx[j], pb_xxzz[j], pb_zz[j], r_0_0[j]);

                t_yy_xyyy[j] = kinvecfunc::fvec_yy_xyyy_s_0(fx[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], s_0_0[j]) + kinvecfunc::fvec_yy_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], r_0_0[j]);

                t_yy_xyyz[j] = kinvecfunc::fvec_yy_xyyz_s_0(fx[j], pa_y[j], pa_yy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yy_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yy_xyzz[j] = kinvecfunc::fvec_yy_xyzz_s_0(fx[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yy_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], r_0_0[j]);

                t_yy_xzzz[j] = kinvecfunc::fvec_yy_xzzz_s_0(fx[j], pa_yy[j], pb_xz[j], pb_xzzz[j], s_0_0[j]) + kinvecfunc::fvec_yy_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_xz[j], pb_xzzz[j], r_0_0[j]);

                t_yy_yyyy[j] = kinvecfunc::fvec_yy_yyyy_s_0(fx[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_yy_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], r_0_0[j]);

                t_yy_yyyz[j] = kinvecfunc::fvec_yy_yyyz_s_0(fx[j], pa_y[j], pa_yy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yy_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yy_yyzz[j] = kinvecfunc::fvec_yy_yyzz_s_0(fx[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yy_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_yy_yzzz[j] = kinvecfunc::fvec_yy_yzzz_s_0(fx[j], pa_y[j], pa_yy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yy_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);

                t_yy_zzzz[j] = kinvecfunc::fvec_yy_zzzz_s_0(fx[j], pa_yy[j], pb_zz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_yy_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pb_zz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (6) = (60,70)

            #pragma omp simd aligned(fgb, fx, fz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_yz_xxxx, \
                                     t_yz_xxxy, t_yz_xxxz, t_yz_xxyy, t_yz_xxyz, t_yz_xxzz, t_yz_xyyy, t_yz_xyyz, \
                                     t_yz_xyzz, t_yz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_xxxx[j] = kinvecfunc::fvec_yz_xxxx_s_0(fx[j], pa_yz[j], pb_xx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_yz_xxxx_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pb_xx[j], pb_xxxx[j], r_0_0[j]);

                t_yz_xxxy[j] = kinvecfunc::fvec_yz_xxxy_s_0(fx[j], pa_yz[j], pa_z[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yz_xxxy_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], r_0_0[j]);

                t_yz_xxxz[j] = kinvecfunc::fvec_yz_xxxz_s_0(fx[j], pa_y[j], pa_yz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yz_xxxz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], r_0_0[j]);

                t_yz_xxyy[j] = kinvecfunc::fvec_yz_xxyy_s_0(fx[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yz_xxyy_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yz_xxyz[j] = kinvecfunc::fvec_yz_xxyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yz_xxyz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yz_xxzz[j] = kinvecfunc::fvec_yz_xxzz_s_0(fx[j], pa_y[j], pa_yz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yz_xxzz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yz_xyyy[j] = kinvecfunc::fvec_yz_xyyy_s_0(fx[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], s_0_0[j]) + kinvecfunc::fvec_yz_xyyy_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], r_0_0[j]);

                t_yz_xyyz[j] = kinvecfunc::fvec_yz_xyyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yz_xyyz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yz_xyzz[j] = kinvecfunc::fvec_yz_xyzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yz_xyzz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_yz_xzzz[j] = kinvecfunc::fvec_yz_xzzz_s_0(fx[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], s_0_0[j]) + kinvecfunc::fvec_yz_xzzz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (7) = (70,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, \
                                     t_yz_yyyy, t_yz_yyyz, t_yz_yyzz, t_yz_yzzz, t_yz_zzzz, t_zz_xxxx, t_zz_xxxy, \
                                     t_zz_xxxz, t_zz_xxyy, t_zz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yz_yyyy[j] = kinvecfunc::fvec_yz_yyyy_s_0(fx[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_yz_yyyy_r_0(fgb[j], fx[j], fz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], r_0_0[j]);

                t_yz_yyyz[j] = kinvecfunc::fvec_yz_yyyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yz_yyyz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yz_yyzz[j] = kinvecfunc::fvec_yz_yyzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yz_yyzz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yz_yzzz[j] = kinvecfunc::fvec_yz_yzzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yz_yzzz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);

                t_yz_zzzz[j] = kinvecfunc::fvec_yz_zzzz_s_0(fx[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_yz_zzzz_r_0(fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], r_0_0[j]);

                t_zz_xxxx[j] = kinvecfunc::fvec_zz_xxxx_s_0(fx[j], pa_zz[j], pb_xx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_zz_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_xx[j], pb_xxxx[j], r_0_0[j]);

                t_zz_xxxy[j] = kinvecfunc::fvec_zz_xxxy_s_0(fx[j], pa_zz[j], pb_xxxy[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_zz_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_xxxy[j], pb_xy[j], r_0_0[j]);

                t_zz_xxxz[j] = kinvecfunc::fvec_zz_xxxz_s_0(fx[j], pa_z[j], pa_zz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_zz_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], r_0_0[j]);

                t_zz_xxyy[j] = kinvecfunc::fvec_zz_xxyy_s_0(fx[j], pa_zz[j], pb_xx[j], pb_xxyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_zz_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_xx[j], pb_xxyy[j], pb_yy[j], r_0_0[j]);

                t_zz_xxyz[j] = kinvecfunc::fvec_zz_xxyz_s_0(fx[j], pa_z[j], pa_zz[j], pb_xxy[j], pb_xxyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_zz_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_xxy[j], pb_xxyz[j], pb_y[j], pb_yz[j], r_0_0[j]);
            }

            // Batch of Integrals (8) = (80,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, \
                                     t_zz_xxzz, t_zz_xyyy, t_zz_xyyz, t_zz_xyzz, t_zz_xzzz, t_zz_yyyy, t_zz_yyyz, \
                                     t_zz_yyzz, t_zz_yzzz, t_zz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zz_xxzz[j] = kinvecfunc::fvec_zz_xxzz_s_0(fx[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_zz_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_zz_xyyy[j] = kinvecfunc::fvec_zz_xyyy_s_0(fx[j], pa_zz[j], pb_xy[j], pb_xyyy[j], s_0_0[j]) + kinvecfunc::fvec_zz_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_xy[j], pb_xyyy[j], r_0_0[j]);

                t_zz_xyyz[j] = kinvecfunc::fvec_zz_xyyz_s_0(fx[j], pa_z[j], pa_zz[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_zz_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], r_0_0[j]);

                t_zz_xyzz[j] = kinvecfunc::fvec_zz_xyzz_s_0(fx[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], s_0_0[j]) + kinvecfunc::fvec_zz_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], r_0_0[j]);

                t_zz_xzzz[j] = kinvecfunc::fvec_zz_xzzz_s_0(fx[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], s_0_0[j]) + kinvecfunc::fvec_zz_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], r_0_0[j]);

                t_zz_yyyy[j] = kinvecfunc::fvec_zz_yyyy_s_0(fx[j], pa_zz[j], pb_yy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_zz_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pb_yy[j], pb_yyyy[j], r_0_0[j]);

                t_zz_yyyz[j] = kinvecfunc::fvec_zz_yyyz_s_0(fx[j], pa_z[j], pa_zz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_zz_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], r_0_0[j]);

                t_zz_yyzz[j] = kinvecfunc::fvec_zz_yyzz_s_0(fx[j], pa_z[j], pa_zz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_zz_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_zz_yzzz[j] = kinvecfunc::fvec_zz_yzzz_s_0(fx[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], s_0_0[j]) + kinvecfunc::fvec_zz_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], r_0_0[j]);

                t_zz_zzzz[j] = kinvecfunc::fvec_zz_zzzz_s_0(fx[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_zz_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], r_0_0[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForGD(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        kinrecfunc::compKineticEnergyForGD_0_4(primBuffer, auxBuffer, osFactors, paDistances, pbDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
       
        kinrecfunc::compKineticEnergyForGD_5_8(primBuffer, auxBuffer, osFactors, paDistances, pbDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
    }
    
    void
    compKineticEnergyForGD_0_4(      CMemBlock2D<double>& primBuffer,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

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

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

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

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, \
                                     pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_xxxx_xx, \
                                     t_xxxx_xy, t_xxxx_xz, t_xxxx_yy, t_xxxx_yz, t_xxxx_zz, t_xxxy_xx, t_xxxy_xy, \
                                     t_xxxy_xz, t_xxxy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xx[j] = kinvecfunc::fvec_xxxx_xx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xxxx_xy[j] = kinvecfunc::fvec_xxxx_xy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxxx_xz[j] = kinvecfunc::fvec_xxxx_xz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxxx_yy[j] = kinvecfunc::fvec_xxxx_yy_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxxx[j], pb_yy[j], r_0_0[j]);

                t_xxxx_yz[j] = kinvecfunc::fvec_xxxx_yz_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_yz_r_0(fga[j], fx[j], fz[j], pa_xx[j], pa_xxxx[j], pb_yz[j], r_0_0[j]);

                t_xxxx_zz[j] = kinvecfunc::fvec_xxxx_zz_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxxx[j], pb_zz[j], r_0_0[j]);

                t_xxxy_xx[j] = kinvecfunc::fvec_xxxy_xx_s_0(fx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xxxy_xy[j] = kinvecfunc::fvec_xxxy_xy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxxy_xz[j] = kinvecfunc::fvec_xxxy_xz_s_0(fx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_xz_r_0(fga[j], fx[j], fz[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxxy_yy[j] = kinvecfunc::fvec_xxxy_yy_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_y[j], pb_yy[j], r_0_0[j]);
            }

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxxz, pa_xxy, pa_xxyy, pa_xxz, \
                                     pa_xy, pa_xyy, pa_xz, pa_y, pa_yy, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, pb_zz, r_0_0, s_0_0, t_xxxy_yz, t_xxxy_zz, t_xxxz_xx, t_xxxz_xy, t_xxxz_xz, \
                                     t_xxxz_yy, t_xxxz_yz, t_xxxz_zz, t_xxyy_xx, t_xxyy_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_yz[j] = kinvecfunc::fvec_xxxy_yz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_yz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxxy_zz[j] = kinvecfunc::fvec_xxxy_zz_s_0(fx[j], pa_xxxy[j], pa_xy[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxy[j], pa_xy[j], pb_zz[j], r_0_0[j]);

                t_xxxz_xx[j] = kinvecfunc::fvec_xxxz_xx_s_0(fx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xxxz_xy[j] = kinvecfunc::fvec_xxxz_xy_s_0(fx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_xy_r_0(fga[j], fx[j], fz[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxxz_xz[j] = kinvecfunc::fvec_xxxz_xz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxxz_yy[j] = kinvecfunc::fvec_xxxz_yy_s_0(fx[j], pa_xxxz[j], pa_xz[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxz[j], pa_xz[j], pb_yy[j], r_0_0[j]);

                t_xxxz_yz[j] = kinvecfunc::fvec_xxxz_yz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_yz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xxxz_zz[j] = kinvecfunc::fvec_xxxz_zz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xxyy_xx[j] = kinvecfunc::fvec_xxyy_xx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xxyy_xy[j] = kinvecfunc::fvec_xxyy_xy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);
            }

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xxyz, pa_xxz, pa_xy, pa_xyy, \
                                     pa_xyz, pa_xz, pa_y, pa_yy, pa_yz, pa_z, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     pb_zz, r_0_0, s_0_0, t_xxyy_xz, t_xxyy_yy, t_xxyy_yz, t_xxyy_zz, t_xxyz_xx, \
                                     t_xxyz_xy, t_xxyz_xz, t_xxyz_yy, t_xxyz_yz, t_xxyz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xz[j] = kinvecfunc::fvec_xxyy_xz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxyy_yy[j] = kinvecfunc::fvec_xxyy_yy_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xxyy_yz[j] = kinvecfunc::fvec_xxyy_yz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_yz_r_0(fga[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxyy_zz[j] = kinvecfunc::fvec_xxyy_zz_s_0(fx[j], pa_xx[j], pa_xxyy[j], pa_yy[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxyy[j], pa_yy[j], pb_zz[j], r_0_0[j]);

                t_xxyz_xx[j] = kinvecfunc::fvec_xxyz_xx_s_0(fx[j], pa_xxyz[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxyz[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xxyz_xy[j] = kinvecfunc::fvec_xxyz_xy_s_0(fx[j], pa_xxyz[j], pa_xxz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_xy_r_0(fga[j], fx[j], fz[j], pa_xxyz[j], pa_xxz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxyz_xz[j] = kinvecfunc::fvec_xxyz_xz_s_0(fx[j], pa_xxy[j], pa_xxyz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_xz_r_0(fga[j], fx[j], fz[j], pa_xxy[j], pa_xxyz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxyz_yy[j] = kinvecfunc::fvec_xxyz_yy_s_0(fx[j], pa_xxyz[j], pa_xxz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxyz[j], pa_xxz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xxyz_yz[j] = kinvecfunc::fvec_xxyz_yz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_yz_r_0(fga[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxyz_zz[j] = kinvecfunc::fvec_xxyz_zz_s_0(fx[j], pa_xxy[j], pa_xxyz[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xxyz[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], r_0_0[j]);
            }

            // Batch of Integrals (3) = (30,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xy, pa_xyy, pa_xyyy, pa_xz, \
                                     pa_xzz, pa_y, pa_yy, pa_yyy, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, pb_zz, r_0_0, s_0_0, t_xxzz_xx, t_xxzz_xy, t_xxzz_xz, t_xxzz_yy, t_xxzz_yz, \
                                     t_xxzz_zz, t_xyyy_xx, t_xyyy_xy, t_xyyy_xz, t_xyyy_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xx[j] = kinvecfunc::fvec_xxzz_xx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xxzz_xy[j] = kinvecfunc::fvec_xxzz_xy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxzz_xz[j] = kinvecfunc::fvec_xxzz_xz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxzz_yy[j] = kinvecfunc::fvec_xxzz_yy_s_0(fx[j], pa_xx[j], pa_xxzz[j], pa_zz[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxzz[j], pa_zz[j], pb_yy[j], r_0_0[j]);

                t_xxzz_yz[j] = kinvecfunc::fvec_xxzz_yz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_yz_r_0(fga[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xxzz_zz[j] = kinvecfunc::fvec_xxzz_zz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyyy_xx[j] = kinvecfunc::fvec_xyyy_xx_s_0(fx[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xyyy_xy[j] = kinvecfunc::fvec_xyyy_xy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyyy_xz[j] = kinvecfunc::fvec_xyyy_xz_s_0(fx[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_xz_r_0(fga[j], fx[j], fz[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyyy_yy[j] = kinvecfunc::fvec_xyyy_yy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_y[j], pb_yy[j], r_0_0[j]);
            }

            // Batch of Integrals (4) = (40,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_xyyz, pa_xyz, pa_xyzz, pa_xz, \
                                     pa_xzz, pa_y, pa_yy, pa_yyz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_xyyy_yz, t_xyyy_zz, t_xyyz_xx, t_xyyz_xy, \
                                     t_xyyz_xz, t_xyyz_yy, t_xyyz_yz, t_xyyz_zz, t_xyzz_xx, t_xyzz_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_yz[j] = kinvecfunc::fvec_xyyy_yz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_yz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyyy_zz[j] = kinvecfunc::fvec_xyyy_zz_s_0(fx[j], pa_xy[j], pa_xyyy[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyyy[j], pb_zz[j], r_0_0[j]);

                t_xyyz_xx[j] = kinvecfunc::fvec_xyyz_xx_s_0(fx[j], pa_xyyz[j], pa_xz[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xyyz[j], pa_xz[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xyyz_xy[j] = kinvecfunc::fvec_xyyz_xy_s_0(fx[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_xy_r_0(fga[j], fx[j], fz[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyyz_xz[j] = kinvecfunc::fvec_xyyz_xz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyyz_yy[j] = kinvecfunc::fvec_xyyz_yy_s_0(fx[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyyz_yz[j] = kinvecfunc::fvec_xyyz_yz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_yz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyyz_zz[j] = kinvecfunc::fvec_xyyz_zz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyzz_xx[j] = kinvecfunc::fvec_xyzz_xx_s_0(fx[j], pa_xy[j], pa_xyzz[j], pa_y[j], pa_yzz[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyzz[j], pa_y[j], pa_yzz[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xyzz_xy[j] = kinvecfunc::fvec_xyzz_xy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_xy_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_y[j], r_0_0[j]);
            }

            idx++;
        }
    }

    void
    compKineticEnergyForGD_5_8(      CMemBlock2D<double>& primBuffer,
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

            auto fx = osFactors.data(4 * idx);

            auto fz = osFactors.data(4 * idx + 1);

            auto fga = osFactors.data(4 * idx + 2);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to 1-th order tensor of distance R(PA)

            auto pa_x = paDistances.data(34 * idx);

            auto pa_y = paDistances.data(34 * idx + 1);

            auto pa_z = paDistances.data(34 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_xy = paDistances.data(34 * idx + 4);

            auto pa_xz = paDistances.data(34 * idx + 5);

            auto pa_yy = paDistances.data(34 * idx + 6);

            auto pa_yz = paDistances.data(34 * idx + 7);

            auto pa_zz = paDistances.data(34 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyz = paDistances.data(34 * idx + 13);

            auto pa_xzz = paDistances.data(34 * idx + 14);

            auto pa_yyy = paDistances.data(34 * idx + 15);

            auto pa_yyz = paDistances.data(34 * idx + 16);

            auto pa_yzz = paDistances.data(34 * idx + 17);

            auto pa_zzz = paDistances.data(34 * idx + 18);

            // set up pointers to 4-th order tensor of distance R(PA)

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

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

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

            // Batch of Integrals (5) = (50,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_xzzz, pa_y, pa_yz, \
                                     pa_yzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, \
                                     r_0_0, s_0_0, t_xyzz_xz, t_xyzz_yy, t_xyzz_yz, t_xyzz_zz, t_xzzz_xx, t_xzzz_xy, \
                                     t_xzzz_xz, t_xzzz_yy, t_xzzz_yz, t_xzzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xz[j] = kinvecfunc::fvec_xyzz_xz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_xz_r_0(fga[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyzz_yy[j] = kinvecfunc::fvec_xyzz_yy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyzz_yz[j] = kinvecfunc::fvec_xyzz_yz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_yz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyzz_zz[j] = kinvecfunc::fvec_xyzz_zz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xzzz_xx[j] = kinvecfunc::fvec_xzzz_xx_s_0(fx[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xx[j], r_0_0[j]);

                t_xzzz_xy[j] = kinvecfunc::fvec_xzzz_xy_s_0(fx[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_xy_r_0(fga[j], fx[j], fz[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xzzz_xz[j] = kinvecfunc::fvec_xzzz_xz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_xz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xzzz_yy[j] = kinvecfunc::fvec_xzzz_yy_s_0(fx[j], pa_xz[j], pa_xzzz[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xz[j], pa_xzzz[j], pb_yy[j], r_0_0[j]);

                t_xzzz_yz[j] = kinvecfunc::fvec_xzzz_yz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_yz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xzzz_zz[j] = kinvecfunc::fvec_xzzz_zz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_z[j], pb_zz[j], r_0_0[j]);
            }

            // Batch of Integrals (6) = (60,70)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, \
                                     pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_yyyy_xx, \
                                     t_yyyy_xy, t_yyyy_xz, t_yyyy_yy, t_yyyy_yz, t_yyyy_zz, t_yyyz_xx, t_yyyz_xy, \
                                     t_yyyz_xz, t_yyyz_yy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xx[j] = kinvecfunc::fvec_yyyy_xx_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyyy[j], pb_xx[j], r_0_0[j]);

                t_yyyy_xy[j] = kinvecfunc::fvec_yyyy_xy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_x[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_xy_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_x[j], pb_xy[j], r_0_0[j]);

                t_yyyy_xz[j] = kinvecfunc::fvec_yyyy_xz_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_xz_r_0(fga[j], fx[j], fz[j], pa_yy[j], pa_yyyy[j], pb_xz[j], r_0_0[j]);

                t_yyyy_yy[j] = kinvecfunc::fvec_yyyy_yy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yyyy_yz[j] = kinvecfunc::fvec_yyyy_yz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_yz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyyy_zz[j] = kinvecfunc::fvec_yyyy_zz_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyyy[j], pb_zz[j], r_0_0[j]);

                t_yyyz_xx[j] = kinvecfunc::fvec_yyyz_xx_s_0(fx[j], pa_yyyz[j], pa_yz[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyyz[j], pa_yz[j], pb_xx[j], r_0_0[j]);

                t_yyyz_xy[j] = kinvecfunc::fvec_yyyz_xy_s_0(fx[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_xy_r_0(fga[j], fx[j], fz[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], r_0_0[j]);

                t_yyyz_xz[j] = kinvecfunc::fvec_yyyz_xz_s_0(fx[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_x[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_xz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_x[j], pb_xz[j], r_0_0[j]);

                t_yyyz_yy[j] = kinvecfunc::fvec_yyyz_yy_s_0(fx[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], r_0_0[j]);
            }

            // Batch of Integrals (7) = (70,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yyzz, pa_yz, pa_yzz, \
                                     pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, \
                                     r_0_0, s_0_0, t_yyyz_yz, t_yyyz_zz, t_yyzz_xx, t_yyzz_xy, t_yyzz_xz, t_yyzz_yy, \
                                     t_yyzz_yz, t_yyzz_zz, t_yzzz_xx, t_yzzz_xy: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_yz[j] = kinvecfunc::fvec_yyyz_yz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_yz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyyz_zz[j] = kinvecfunc::fvec_yyyz_zz_s_0(fx[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yyzz_xx[j] = kinvecfunc::fvec_yyzz_xx_s_0(fx[j], pa_yy[j], pa_yyzz[j], pa_zz[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyzz[j], pa_zz[j], pb_xx[j], r_0_0[j]);

                t_yyzz_xy[j] = kinvecfunc::fvec_yyzz_xy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_xy_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], r_0_0[j]);

                t_yyzz_xz[j] = kinvecfunc::fvec_yyzz_xz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_xz_r_0(fga[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], r_0_0[j]);

                t_yyzz_yy[j] = kinvecfunc::fvec_yyzz_yy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yyzz_yz[j] = kinvecfunc::fvec_yyzz_yz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_yz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyzz_zz[j] = kinvecfunc::fvec_yyzz_zz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yzzz_xx[j] = kinvecfunc::fvec_yzzz_xx_s_0(fx[j], pa_yz[j], pa_yzzz[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yz[j], pa_yzzz[j], pb_xx[j], r_0_0[j]);

                t_yzzz_xy[j] = kinvecfunc::fvec_yzzz_xy_s_0(fx[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_xy_r_0(fga[j], fx[j], fz[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xy[j], r_0_0[j]);
            }

            // Batch of Integrals (8) = (80,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, \
                                     pb_xx, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, pb_zz, r_0_0, s_0_0, t_yzzz_xz, \
                                     t_yzzz_yy, t_yzzz_yz, t_yzzz_zz, t_zzzz_xx, t_zzzz_xy, t_zzzz_xz, t_zzzz_yy, \
                                     t_zzzz_yz, t_zzzz_zz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xz[j] = kinvecfunc::fvec_yzzz_xz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_x[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_xz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_x[j], pb_xz[j], r_0_0[j]);

                t_yzzz_yy[j] = kinvecfunc::fvec_yzzz_yy_s_0(fx[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yzzz_yz[j] = kinvecfunc::fvec_yzzz_yz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_yz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yzzz_zz[j] = kinvecfunc::fvec_yzzz_zz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_zzzz_xx[j] = kinvecfunc::fvec_zzzz_xx_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_xx[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_xx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pa_zzzz[j], pb_xx[j], r_0_0[j]);

                t_zzzz_xy[j] = kinvecfunc::fvec_zzzz_xy_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_xy_r_0(fga[j], fx[j], fz[j], pa_zz[j], pa_zzzz[j], pb_xy[j], r_0_0[j]);

                t_zzzz_xz[j] = kinvecfunc::fvec_zzzz_xz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_x[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_xz_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_x[j], pb_xz[j], r_0_0[j]);

                t_zzzz_yy[j] = kinvecfunc::fvec_zzzz_yy_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_yy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pa_zzzz[j], pb_yy[j], r_0_0[j]);

                t_zzzz_yz[j] = kinvecfunc::fvec_zzzz_yz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_yz_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_zzzz_zz[j] = kinvecfunc::fvec_zzzz_zz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_zz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_z[j], pb_zz[j], r_0_0[j]);
            }

            idx++;
        }
    }

} // kinrecfunc namespace

