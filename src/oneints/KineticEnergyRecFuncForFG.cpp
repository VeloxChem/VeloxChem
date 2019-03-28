//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForFG.hpp"

#include "KineticEnergyVecFuncForFG.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForFG(      CMemBlock2D<double>& primBuffer,
                           const CMemBlock2D<double>& auxBuffer,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pbDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        kinrecfunc::compKineticEnergyForFG_0_9(primBuffer, auxBuffer, osFactors, paDistances, pbDistances,
                                               braGtoBlock, ketGtoBlock, iContrGto);
        
        kinrecfunc::compKineticEnergyForFG_10_18(primBuffer, auxBuffer, osFactors, paDistances, pbDistances,
                                                 braGtoBlock, ketGtoBlock, iContrGto);
        
        kinrecfunc::compKineticEnergyForFG_19_27(primBuffer, auxBuffer, osFactors, paDistances, pbDistances,
                                                 braGtoBlock, ketGtoBlock, iContrGto);
    }

    void
    compKineticEnergyForFG_0_9(      CMemBlock2D<double>& primBuffer,
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

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xxx = paDistances.data(19 * idx + 9);

            auto pa_xxy = paDistances.data(19 * idx + 10);

            auto pa_xxz = paDistances.data(19 * idx + 11);

            auto pa_xyy = paDistances.data(19 * idx + 12);

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

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, r_0_0, s_0_0, t_xxx_xxxx, t_xxx_xxxy, t_xxx_xxxz, t_xxx_xxyy, \
                                     t_xxx_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xxxx[j] = kinvecfunc::fvec_xxx_xxxx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], r_0_0[j]);

                t_xxx_xxxy[j] = kinvecfunc::fvec_xxx_xxxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxx_xxxz[j] = kinvecfunc::fvec_xxx_xxxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxx_xxyy[j] = kinvecfunc::fvec_xxx_xxyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xxx_xxyz[j] = kinvecfunc::fvec_xxx_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxyz[j], pb_xyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxyz[j], pb_xyz[j], pb_yz[j], r_0_0[j]);
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, \
                                     s_0_0, t_xxx_xxzz, t_xxx_xyyy, t_xxx_xyyz, t_xxx_xyzz, t_xxx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xxzz[j] = kinvecfunc::fvec_xxx_xxzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xxx_xyyy[j] = kinvecfunc::fvec_xxx_xyyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xxx_xyyz[j] = kinvecfunc::fvec_xxx_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xyyz[j], pb_xz[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xyyz[j], pb_xz[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xxx_xyzz[j] = kinvecfunc::fvec_xxx_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xy[j], pb_xyzz[j], pb_y[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xy[j], pb_xyzz[j], pb_y[j], pb_yzz[j], r_0_0[j]);

                t_xxx_xzzz[j] = kinvecfunc::fvec_xxx_xzzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xxx, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, \
                                     pb_zz, pb_zzzz, r_0_0, s_0_0, t_xxx_yyyy, t_xxx_yyyz, t_xxx_yyzz, t_xxx_yzzz, \
                                     t_xxx_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_yyyy[j] = kinvecfunc::fvec_xxx_yyyy_s_0(fx[j], pa_x[j], pa_xxx[j], pb_yy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_xxx_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_yy[j], pb_yyyy[j], r_0_0[j]);

                t_xxx_yyyz[j] = kinvecfunc::fvec_xxx_yyyz_s_0(fx[j], pa_x[j], pa_xxx[j], pb_yyyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_yyyz[j], pb_yz[j], r_0_0[j]);

                t_xxx_yyzz[j] = kinvecfunc::fvec_xxx_yyzz_s_0(fx[j], pa_x[j], pa_xxx[j], pb_yy[j], pb_yyzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_yy[j], pb_yyzz[j], pb_zz[j], r_0_0[j]);

                t_xxx_yzzz[j] = kinvecfunc::fvec_xxx_yzzz_s_0(fx[j], pa_x[j], pa_xxx[j], pb_yz[j], pb_yzzz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_yz[j], pb_yzzz[j], r_0_0[j]);

                t_xxx_zzzz[j] = kinvecfunc::fvec_xxx_zzzz_s_0(fx[j], pa_x[j], pa_xxx[j], pb_zz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_zz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xxy_xxxx, t_xxy_xxxy, t_xxy_xxxz, t_xxy_xxyy, \
                                     t_xxy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xxxx[j] = kinvecfunc::fvec_xxy_xxxx_s_0(fx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], r_0_0[j]);

                t_xxy_xxxy[j] = kinvecfunc::fvec_xxy_xxxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxy_xxxz[j] = kinvecfunc::fvec_xxy_xxxz_s_0(fx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxy_xxyy[j] = kinvecfunc::fvec_xxy_xxyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xxy_xxyz[j] = kinvecfunc::fvec_xxy_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxyz[j], pb_xxz[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxyz[j], pb_xxz[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxy_xxzz, t_xxy_xyyy, \
                                     t_xxy_xyyz, t_xxy_xyzz, t_xxy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xxzz[j] = kinvecfunc::fvec_xxy_xxzz_s_0(fx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xxy_xyyy[j] = kinvecfunc::fvec_xxy_xyyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xxy_xyyz[j] = kinvecfunc::fvec_xxy_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxy_xyzz[j] = kinvecfunc::fvec_xxy_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xxy_xzzz[j] = kinvecfunc::fvec_xxy_xzzz_s_0(fx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_xx, pa_xxy, pa_y, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xxy_yyyy, \
                                     t_xxy_yyyz, t_xxy_yyzz, t_xxy_yzzz, t_xxy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_yyyy[j] = kinvecfunc::fvec_xxy_yyyy_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_xxy_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], r_0_0[j]);

                t_xxy_yyyz[j] = kinvecfunc::fvec_xxy_yyyz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxy_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxy_yyzz[j] = kinvecfunc::fvec_xxy_yyzz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xxy_yzzz[j] = kinvecfunc::fvec_xxy_yzzz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);

                t_xxy_zzzz[j] = kinvecfunc::fvec_xxy_zzzz_s_0(fx[j], pa_xxy[j], pa_y[j], pb_zz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_y[j], pb_zz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xxz_xxxx, t_xxz_xxxy, t_xxz_xxxz, t_xxz_xxyy, \
                                     t_xxz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xxxx[j] = kinvecfunc::fvec_xxz_xxxx_s_0(fx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], r_0_0[j]);

                t_xxz_xxxy[j] = kinvecfunc::fvec_xxz_xxxy_s_0(fx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxz_xxxz[j] = kinvecfunc::fvec_xxz_xxxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxz_xxyy[j] = kinvecfunc::fvec_xxz_xxyy_s_0(fx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xxz_xxyz[j] = kinvecfunc::fvec_xxz_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xxyz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xxyz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], r_0_0[j]);
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxz_xxzz, \
                                     t_xxz_xyyy, t_xxz_xyyz, t_xxz_xyzz, t_xxz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xxzz[j] = kinvecfunc::fvec_xxz_xxzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xxz_xyyy[j] = kinvecfunc::fvec_xxz_xyyy_s_0(fx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xxz_xyyz[j] = kinvecfunc::fvec_xxz_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xxz_xyzz[j] = kinvecfunc::fvec_xxz_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_xxz_xzzz[j] = kinvecfunc::fvec_xxz_xzzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_xx, pa_xxz, pa_z, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xxz_yyyy, \
                                     t_xxz_yyyz, t_xxz_yyzz, t_xxz_yzzz, t_xxz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_yyyy[j] = kinvecfunc::fvec_xxz_yyyy_s_0(fx[j], pa_xxz[j], pa_z[j], pb_yy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_xxz_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_z[j], pb_yy[j], pb_yyyy[j], r_0_0[j]);

                t_xxz_yyyz[j] = kinvecfunc::fvec_xxz_yyyz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], r_0_0[j]);

                t_xxz_yyzz[j] = kinvecfunc::fvec_xxz_yyzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xxz_yzzz[j] = kinvecfunc::fvec_xxz_yzzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], r_0_0[j]);

                t_xxz_zzzz[j] = kinvecfunc::fvec_xxz_zzzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xyy_xxxx, t_xyy_xxxy, t_xyy_xxxz, t_xyy_xxyy, \
                                     t_xyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xxxx[j] = kinvecfunc::fvec_xyy_xxxx_s_0(fx[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], r_0_0[j]);

                t_xyy_xxxy[j] = kinvecfunc::fvec_xyy_xxxy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyy_xxxz[j] = kinvecfunc::fvec_xyy_xxxz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyy_xxyy[j] = kinvecfunc::fvec_xyy_xxyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyy_xxyz[j] = kinvecfunc::fvec_xyy_xxyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xxyz[j], pb_xxz[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xxyz[j], pb_xxz[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForFG_10_18(      CMemBlock2D<double>& primBuffer,
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

            auto pa_xy = paDistances.data(19 * idx + 4);

            auto pa_xz = paDistances.data(19 * idx + 5);

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

            auto pa_xyy = paDistances.data(19 * idx + 12);

            auto pa_xyz = paDistances.data(19 * idx + 13);

            auto pa_xzz = paDistances.data(19 * idx + 14);

            auto pa_yyy = paDistances.data(19 * idx + 15);

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

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyy_xxzz, t_xyy_xyyy, \
                                     t_xyy_xyyz, t_xyy_xyzz, t_xyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xxzz[j] = kinvecfunc::fvec_xyy_xxzz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xyy_xyyy[j] = kinvecfunc::fvec_xyy_xyyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xyy_xyyz[j] = kinvecfunc::fvec_xyy_xyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyy_xyzz[j] = kinvecfunc::fvec_xyy_xyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xyy_xzzz[j] = kinvecfunc::fvec_xyy_xzzz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xyy_yyyy, \
                                     t_xyy_yyyz, t_xyy_yyzz, t_xyy_yzzz, t_xyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_yyyy[j] = kinvecfunc::fvec_xyy_yyyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_xyy_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], r_0_0[j]);

                t_xyy_yyyz[j] = kinvecfunc::fvec_xyy_yyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyy_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyy_yyzz[j] = kinvecfunc::fvec_xyy_yyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xyy_yzzz[j] = kinvecfunc::fvec_xyy_yzzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);

                t_xyy_zzzz[j] = kinvecfunc::fvec_xyy_zzzz_s_0(fx[j], pa_x[j], pa_xyy[j], pb_zz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pb_zz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (12) = (60,65)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xyz_xxxx, t_xyz_xxxy, t_xyz_xxxz, \
                                     t_xyz_xxyy, t_xyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xxxx[j] = kinvecfunc::fvec_xyz_xxxx_s_0(fx[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xxxx_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], r_0_0[j]);

                t_xyz_xxxy[j] = kinvecfunc::fvec_xyz_xxxy_s_0(fx[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xxxy_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyz_xxxz[j] = kinvecfunc::fvec_xyz_xxxz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xxxz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyz_xxyy[j] = kinvecfunc::fvec_xyz_xxyy_s_0(fx[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xxyy_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyz_xxyz[j] = kinvecfunc::fvec_xyz_xxyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xxyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (13) = (65,70)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyz_xxzz, \
                                     t_xyz_xyyy, t_xyz_xyyz, t_xyz_xyzz, t_xyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xxzz[j] = kinvecfunc::fvec_xyz_xxzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xxzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyz_xyyy[j] = kinvecfunc::fvec_xyz_xyyy_s_0(fx[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xyyy_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xyz_xyyz[j] = kinvecfunc::fvec_xyz_xyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xyyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyz_xyzz[j] = kinvecfunc::fvec_xyz_xyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xyzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyz_xzzz[j] = kinvecfunc::fvec_xyz_xzzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xzzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (14) = (70,75)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, \
                                     t_xyz_yyyy, t_xyz_yyyz, t_xyz_yyzz, t_xyz_yzzz, t_xyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_yyyy[j] = kinvecfunc::fvec_xyz_yyyy_s_0(fx[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_xyz_yyyy_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], r_0_0[j]);

                t_xyz_yyyz[j] = kinvecfunc::fvec_xyz_yyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyz_yyyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyz_yyzz[j] = kinvecfunc::fvec_xyz_yyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_yyzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyz_yzzz[j] = kinvecfunc::fvec_xyz_yzzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_yzzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);

                t_xyz_zzzz[j] = kinvecfunc::fvec_xyz_zzzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_zzzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (15) = (75,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xzz_xxxx, t_xzz_xxxy, t_xzz_xxxz, t_xzz_xxyy, \
                                     t_xzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xxxx[j] = kinvecfunc::fvec_xzz_xxxx_s_0(fx[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], r_0_0[j]);

                t_xzz_xxxy[j] = kinvecfunc::fvec_xzz_xxxy_s_0(fx[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xzz_xxxz[j] = kinvecfunc::fvec_xzz_xxxz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xzz_xxyy[j] = kinvecfunc::fvec_xzz_xxyy_s_0(fx[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xzz_xxyz[j] = kinvecfunc::fvec_xzz_xxyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xxy[j], pb_xxyz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xxy[j], pb_xxyz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], r_0_0[j]);
            }

            // Batch of Integrals (16) = (80,85)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xzz_xxzz, \
                                     t_xzz_xyyy, t_xzz_xyyz, t_xzz_xyzz, t_xzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xxzz[j] = kinvecfunc::fvec_xzz_xxzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xzz_xyyy[j] = kinvecfunc::fvec_xzz_xyyy_s_0(fx[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xzz_xyyz[j] = kinvecfunc::fvec_xzz_xyyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xzz_xyzz[j] = kinvecfunc::fvec_xzz_xyzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_xzz_xzzz[j] = kinvecfunc::fvec_xzz_xzzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (17) = (85,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xzz_yyyy, \
                                     t_xzz_yyyz, t_xzz_yyzz, t_xzz_yzzz, t_xzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_yyyy[j] = kinvecfunc::fvec_xzz_yyyy_s_0(fx[j], pa_x[j], pa_xzz[j], pb_yy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_xzz_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pb_yy[j], pb_yyyy[j], r_0_0[j]);

                t_xzz_yyyz[j] = kinvecfunc::fvec_xzz_yyyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], r_0_0[j]);

                t_xzz_yyzz[j] = kinvecfunc::fvec_xzz_yyzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xzz_yzzz[j] = kinvecfunc::fvec_xzz_yzzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], r_0_0[j]);

                t_xzz_zzzz[j] = kinvecfunc::fvec_xzz_zzzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (18) = (90,95)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_yyy_xxxx, t_yyy_xxxy, t_yyy_xxxz, t_yyy_xxyy, t_yyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xxxx[j] = kinvecfunc::fvec_yyy_xxxx_s_0(fx[j], pa_y[j], pa_yyy[j], pb_xx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_xx[j], pb_xxxx[j], r_0_0[j]);

                t_yyy_xxxy[j] = kinvecfunc::fvec_yyy_xxxy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], r_0_0[j]);

                t_yyy_xxxz[j] = kinvecfunc::fvec_yyy_xxxz_s_0(fx[j], pa_y[j], pa_yyy[j], pb_xxxz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_xxxz[j], pb_xz[j], r_0_0[j]);

                t_yyy_xxyy[j] = kinvecfunc::fvec_yyy_xxyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yyy_xxyz[j] = kinvecfunc::fvec_yyy_xxyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xxyz[j], pb_xxz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xxyz[j], pb_xxz[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            idx++;
        }
    }
    
    void
    compKineticEnergyForFG_19_27(      CMemBlock2D<double>& primBuffer,
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

            auto pa_y = paDistances.data(19 * idx + 1);

            auto pa_z = paDistances.data(19 * idx + 2);

            // set up pointers to 2-th order tensor of distance R(PA)

            auto pa_yy = paDistances.data(19 * idx + 6);

            auto pa_yz = paDistances.data(19 * idx + 7);

            auto pa_zz = paDistances.data(19 * idx + 8);

            // set up pointers to 3-th order tensor of distance R(PA)

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

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

            // set up pointers to integrals

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

            // Batch of Integrals (19) = (95,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_zz, r_0_0, s_0_0, t_yyy_xxzz, \
                                     t_yyy_xyyy, t_yyy_xyyz, t_yyy_xyzz, t_yyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xxzz[j] = kinvecfunc::fvec_yyy_xxzz_s_0(fx[j], pa_y[j], pa_yyy[j], pb_xx[j], pb_xxzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_xx[j], pb_xxzz[j], pb_zz[j], r_0_0[j]);

                t_yyy_xyyy[j] = kinvecfunc::fvec_yyy_xyyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], r_0_0[j]);

                t_yyy_xyyz[j] = kinvecfunc::fvec_yyy_xyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yyy_xyzz[j] = kinvecfunc::fvec_yyy_xyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], r_0_0[j]);

                t_yyy_xzzz[j] = kinvecfunc::fvec_yyy_xzzz_s_0(fx[j], pa_y[j], pa_yyy[j], pb_xz[j], pb_xzzz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_xz[j], pb_xzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (20) = (100,105)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_yyy_yyyy, \
                                     t_yyy_yyyz, t_yyy_yyzz, t_yyy_yzzz, t_yyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_yyyy[j] = kinvecfunc::fvec_yyy_yyyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_yyy_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], r_0_0[j]);

                t_yyy_yyyz[j] = kinvecfunc::fvec_yyy_yyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyy_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyy_yyzz[j] = kinvecfunc::fvec_yyy_yyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_yyy_yzzz[j] = kinvecfunc::fvec_yyy_yzzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], r_0_0[j]);

                t_yyy_zzzz[j] = kinvecfunc::fvec_yyy_zzzz_s_0(fx[j], pa_y[j], pa_yyy[j], pb_zz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_zz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (21) = (105,110)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, r_0_0, s_0_0, t_yyz_xxxx, t_yyz_xxxy, t_yyz_xxxz, t_yyz_xxyy, t_yyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xxxx[j] = kinvecfunc::fvec_yyz_xxxx_s_0(fx[j], pa_yyz[j], pa_z[j], pb_xx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_z[j], pb_xx[j], pb_xxxx[j], r_0_0[j]);

                t_yyz_xxxy[j] = kinvecfunc::fvec_yyz_xxxy_s_0(fx[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], r_0_0[j]);

                t_yyz_xxxz[j] = kinvecfunc::fvec_yyz_xxxz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], r_0_0[j]);

                t_yyz_xxyy[j] = kinvecfunc::fvec_yyz_xxyy_s_0(fx[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yyz_xxyz[j] = kinvecfunc::fvec_yyz_xxyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (22) = (110,115)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, \
                                     r_0_0, s_0_0, t_yyz_xxzz, t_yyz_xyyy, t_yyz_xyyz, t_yyz_xyzz, t_yyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xxzz[j] = kinvecfunc::fvec_yyz_xxzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yyz_xyyy[j] = kinvecfunc::fvec_yyz_xyyy_s_0(fx[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], r_0_0[j]);

                t_yyz_xyyz[j] = kinvecfunc::fvec_yyz_xyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yyz_xyzz[j] = kinvecfunc::fvec_yyz_xyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_yyz_xzzz[j] = kinvecfunc::fvec_yyz_xzzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (23) = (115,120)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, \
                                     s_0_0, t_yyz_yyyy, t_yyz_yyyz, t_yyz_yyzz, t_yyz_yzzz, t_yyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_yyyy[j] = kinvecfunc::fvec_yyz_yyyy_s_0(fx[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_yyz_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], r_0_0[j]);

                t_yyz_yyyz[j] = kinvecfunc::fvec_yyz_yyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyz_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyz_yyzz[j] = kinvecfunc::fvec_yyz_yyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yyz_yzzz[j] = kinvecfunc::fvec_yyz_yzzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);

                t_yyz_zzzz[j] = kinvecfunc::fvec_yyz_zzzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (24) = (120,125)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, r_0_0, s_0_0, t_yzz_xxxx, t_yzz_xxxy, t_yzz_xxxz, t_yzz_xxyy, t_yzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xxxx[j] = kinvecfunc::fvec_yzz_xxxx_s_0(fx[j], pa_y[j], pa_yzz[j], pb_xx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pb_xx[j], pb_xxxx[j], r_0_0[j]);

                t_yzz_xxxy[j] = kinvecfunc::fvec_yzz_xxxy_s_0(fx[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], r_0_0[j]);

                t_yzz_xxxz[j] = kinvecfunc::fvec_yzz_xxxz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], r_0_0[j]);

                t_yzz_xxyy[j] = kinvecfunc::fvec_yzz_xxyy_s_0(fx[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_yzz_xxyz[j] = kinvecfunc::fvec_yzz_xxyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (25) = (125,130)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, \
                                     r_0_0, s_0_0, t_yzz_xxzz, t_yzz_xyyy, t_yzz_xyyz, t_yzz_xyzz, t_yzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xxzz[j] = kinvecfunc::fvec_yzz_xxzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yzz_xyyy[j] = kinvecfunc::fvec_yzz_xyyy_s_0(fx[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], r_0_0[j]);

                t_yzz_xyyz[j] = kinvecfunc::fvec_yzz_xyyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yzz_xyzz[j] = kinvecfunc::fvec_yzz_xyzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_yzz_xzzz[j] = kinvecfunc::fvec_yzz_xzzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (26) = (130,135)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, \
                                     s_0_0, t_yzz_yyyy, t_yzz_yyyz, t_yzz_yyzz, t_yzz_yzzz, t_yzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_yyyy[j] = kinvecfunc::fvec_yzz_yyyy_s_0(fx[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_yzz_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], r_0_0[j]);

                t_yzz_yyyz[j] = kinvecfunc::fvec_yzz_yyyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yzz_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yzz_yyzz[j] = kinvecfunc::fvec_yzz_yyzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yzz_yzzz[j] = kinvecfunc::fvec_yzz_yzzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);

                t_yzz_zzzz[j] = kinvecfunc::fvec_yzz_zzzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (27) = (135,140)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, r_0_0, s_0_0, \
                                     t_zzz_xxxx, t_zzz_xxxy, t_zzz_xxxz, t_zzz_xxyy, t_zzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xxxx[j] = kinvecfunc::fvec_zzz_xxxx_s_0(fx[j], pa_z[j], pa_zzz[j], pb_xx[j], pb_xxxx[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xxxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_xx[j], pb_xxxx[j], r_0_0[j]);

                t_zzz_xxxy[j] = kinvecfunc::fvec_zzz_xxxy_s_0(fx[j], pa_z[j], pa_zzz[j], pb_xxxy[j], pb_xy[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xxxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_xxxy[j], pb_xy[j], r_0_0[j]);

                t_zzz_xxxz[j] = kinvecfunc::fvec_zzz_xxxz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xxxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], r_0_0[j]);

                t_zzz_xxyy[j] = kinvecfunc::fvec_zzz_xxyy_s_0(fx[j], pa_z[j], pa_zzz[j], pb_xx[j], pb_xxyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xxyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_xx[j], pb_xxyy[j], pb_yy[j], r_0_0[j]);

                t_zzz_xxyz[j] = kinvecfunc::fvec_zzz_xxyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xxy[j], pb_xxyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xxyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xxy[j], pb_xxyz[j], pb_y[j], pb_yz[j], r_0_0[j]);
            }

            // Batch of Integrals (28) = (140,145)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, r_0_0, s_0_0, \
                                     t_zzz_xxzz, t_zzz_xyyy, t_zzz_xyyz, t_zzz_xyzz, t_zzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xxzz[j] = kinvecfunc::fvec_zzz_xxzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xxzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_zzz_xyyy[j] = kinvecfunc::fvec_zzz_xyyy_s_0(fx[j], pa_z[j], pa_zzz[j], pb_xy[j], pb_xyyy[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_xy[j], pb_xyyy[j], r_0_0[j]);

                t_zzz_xyyz[j] = kinvecfunc::fvec_zzz_xyyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], r_0_0[j]);

                t_zzz_xyzz[j] = kinvecfunc::fvec_zzz_xyzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], r_0_0[j]);

                t_zzz_xzzz[j] = kinvecfunc::fvec_zzz_xzzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], r_0_0[j]);
            }

            // Batch of Integrals (29) = (145,150)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_zzz_yyyy, \
                                     t_zzz_yyyz, t_zzz_yyzz, t_zzz_yzzz, t_zzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_yyyy[j] = kinvecfunc::fvec_zzz_yyyy_s_0(fx[j], pa_z[j], pa_zzz[j], pb_yy[j], pb_yyyy[j], s_0_0[j]) + kinvecfunc::fvec_zzz_yyyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_yy[j], pb_yyyy[j], r_0_0[j]);

                t_zzz_yyyz[j] = kinvecfunc::fvec_zzz_yyyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_yyyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], r_0_0[j]);

                t_zzz_yyzz[j] = kinvecfunc::fvec_zzz_yyzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_yyzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_zzz_yzzz[j] = kinvecfunc::fvec_zzz_yzzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_yzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], r_0_0[j]);

                t_zzz_zzzz[j] = kinvecfunc::fvec_zzz_zzzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_zzzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], r_0_0[j]);
            }

            idx++;
        }
    }

} // kinrecfunc namespace

