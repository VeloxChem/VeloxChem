//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForFF.hpp"

#include "KineticEnergyVecFuncForFF.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForFF(      CMemBlock2D<double>& primBuffer,
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

            // Batch of Integrals (0) = (0,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, r_0_0, s_0_0, t_xxx_xxx, t_xxx_xxy, t_xxx_xxz, t_xxx_xyy, t_xxx_xyz, \
                                     t_xxx_xzz, t_xxx_yyy, t_xxx_yyz, t_xxx_yzz, t_xxx_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xxx[j] = kinvecfunc::fvec_xxx_xxx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xxx_xxy[j] = kinvecfunc::fvec_xxx_xxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxx_xxz[j] = kinvecfunc::fvec_xxx_xxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxx_xyy[j] = kinvecfunc::fvec_xxx_xyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xxx_xyz[j] = kinvecfunc::fvec_xxx_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_xyz[j], pb_yz[j], r_0_0[j]);

                t_xxx_xzz[j] = kinvecfunc::fvec_xxx_xzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pb_x[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xxx_yyy[j] = kinvecfunc::fvec_xxx_yyy_s_0(fx[j], pa_x[j], pa_xxx[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxx_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xxx_yyz[j] = kinvecfunc::fvec_xxx_yyz_s_0(fx[j], pa_x[j], pa_xxx[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxx_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xxx_yzz[j] = kinvecfunc::fvec_xxx_yzz_s_0(fx[j], pa_x[j], pa_xxx[j], pb_y[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_y[j], pb_yzz[j], r_0_0[j]);

                t_xxx_zzz[j] = kinvecfunc::fvec_xxx_zzz_s_0(fx[j], pa_x[j], pa_xxx[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxx_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (1) = (10,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxy_xxx, t_xxy_xxy, t_xxy_xxz, t_xxy_xyy, \
                                     t_xxy_xyz, t_xxy_xzz, t_xxy_yyy, t_xxy_yyz, t_xxy_yzz, t_xxy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xxx[j] = kinvecfunc::fvec_xxy_xxx_s_0(fx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xxy_xxy[j] = kinvecfunc::fvec_xxy_xxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxy_xxz[j] = kinvecfunc::fvec_xxy_xxz_s_0(fx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxy_xyy[j] = kinvecfunc::fvec_xxy_xyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xxy_xyz[j] = kinvecfunc::fvec_xxy_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxy_xzz[j] = kinvecfunc::fvec_xxy_xzz_s_0(fx[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xxy_yyy[j] = kinvecfunc::fvec_xxy_yyy_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxy_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xxy_yyz[j] = kinvecfunc::fvec_xxy_yyz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxy_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxy_yzz[j] = kinvecfunc::fvec_xxy_yzz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_y[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xxy_zzz[j] = kinvecfunc::fvec_xxy_zzz_s_0(fx[j], pa_xxy[j], pa_y[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxy_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_y[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (2) = (20,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxz_xxx, t_xxz_xxy, t_xxz_xxz, t_xxz_xyy, \
                                     t_xxz_xyz, t_xxz_xzz, t_xxz_yyy, t_xxz_yyz, t_xxz_yzz, t_xxz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xxx[j] = kinvecfunc::fvec_xxz_xxx_s_0(fx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xxz_xxy[j] = kinvecfunc::fvec_xxz_xxy_s_0(fx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxz_xxz[j] = kinvecfunc::fvec_xxz_xxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxz_xyy[j] = kinvecfunc::fvec_xxz_xyy_s_0(fx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xxz_xyz[j] = kinvecfunc::fvec_xxz_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xxz_xzz[j] = kinvecfunc::fvec_xxz_xzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xxz_yyy[j] = kinvecfunc::fvec_xxz_yyy_s_0(fx[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xxz_yyz[j] = kinvecfunc::fvec_xxz_yyz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xxz_yzz[j] = kinvecfunc::fvec_xxz_yzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_xxz_zzz[j] = kinvecfunc::fvec_xxz_zzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_z[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (3) = (30,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyy_xxx, t_xyy_xxy, t_xyy_xxz, t_xyy_xyy, \
                                     t_xyy_xyz, t_xyy_xzz, t_xyy_yyy, t_xyy_yyz, t_xyy_yzz, t_xyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xxx[j] = kinvecfunc::fvec_xyy_xxx_s_0(fx[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xyy_xxy[j] = kinvecfunc::fvec_xyy_xxy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyy_xxz[j] = kinvecfunc::fvec_xyy_xxz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyy_xyy[j] = kinvecfunc::fvec_xyy_xyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyy_xyz[j] = kinvecfunc::fvec_xyy_xyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyy_xzz[j] = kinvecfunc::fvec_xyy_xzz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xyy_yyy[j] = kinvecfunc::fvec_xyy_yyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xyy_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xyy_yyz[j] = kinvecfunc::fvec_xyy_yyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyy_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyy_yzz[j] = kinvecfunc::fvec_xyy_yzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xyy_zzz[j] = kinvecfunc::fvec_xyy_zzz_s_0(fx[j], pa_x[j], pa_xyy[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xyy_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (4) = (40,50)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyz_xxx, t_xyz_xxy, t_xyz_xxz, \
                                     t_xyz_xyy, t_xyz_xyz, t_xyz_xzz, t_xyz_yyy, t_xyz_yyz, t_xyz_yzz, t_xyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xxx[j] = kinvecfunc::fvec_xyz_xxx_s_0(fx[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xxx_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xyz_xxy[j] = kinvecfunc::fvec_xyz_xxy_s_0(fx[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xxy_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyz_xxz[j] = kinvecfunc::fvec_xyz_xxz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xxz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyz_xyy[j] = kinvecfunc::fvec_xyz_xyy_s_0(fx[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xyy_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyz_xyz[j] = kinvecfunc::fvec_xyz_xyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xyz_r_0(fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyz_xzz[j] = kinvecfunc::fvec_xyz_xzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_xzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyz_yyy[j] = kinvecfunc::fvec_xyz_yyy_s_0(fx[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xyz_yyy_r_0(fgb[j], fx[j], fz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xyz_yyz[j] = kinvecfunc::fvec_xyz_yyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyz_yyz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyz_yzz[j] = kinvecfunc::fvec_xyz_yzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_yzz_r_0(fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyz_zzz[j] = kinvecfunc::fvec_xyz_zzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xyz_zzz_r_0(fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (5) = (50,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xzz_xxx, t_xzz_xxy, t_xzz_xxz, t_xzz_xyy, \
                                     t_xzz_xyz, t_xzz_xzz, t_xzz_yyy, t_xzz_yyz, t_xzz_yzz, t_xzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xxx[j] = kinvecfunc::fvec_xzz_xxx_s_0(fx[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xzz_xxy[j] = kinvecfunc::fvec_xzz_xxy_s_0(fx[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xzz_xxz[j] = kinvecfunc::fvec_xzz_xxz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xzz_xyy[j] = kinvecfunc::fvec_xzz_xyy_s_0(fx[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xzz_xyz[j] = kinvecfunc::fvec_xzz_xyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], r_0_0[j]);

                t_xzz_xzz[j] = kinvecfunc::fvec_xzz_xzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xzz_yyy[j] = kinvecfunc::fvec_xzz_yyy_s_0(fx[j], pa_x[j], pa_xzz[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xzz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xzz[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xzz_yyz[j] = kinvecfunc::fvec_xzz_yyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xzz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xzz_yzz[j] = kinvecfunc::fvec_xzz_yzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_xzz_zzz[j] = kinvecfunc::fvec_xzz_zzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xzz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (6) = (60,70)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, r_0_0, s_0_0, t_yyy_xxx, t_yyy_xxy, t_yyy_xxz, t_yyy_xyy, t_yyy_xyz, \
                                     t_yyy_xzz, t_yyy_yyy, t_yyy_yyz, t_yyy_yzz, t_yyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xxx[j] = kinvecfunc::fvec_yyy_xxx_s_0(fx[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_yyy_xxy[j] = kinvecfunc::fvec_yyy_xxy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xx[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xx[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_yyy_xxz[j] = kinvecfunc::fvec_yyy_xxz_s_0(fx[j], pa_y[j], pa_yyy[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_yyy_xyy[j] = kinvecfunc::fvec_yyy_xyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyy[j], r_0_0[j]);

                t_yyy_xyz[j] = kinvecfunc::fvec_yyy_xyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xyz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yyy_xzz[j] = kinvecfunc::fvec_yyy_xzz_s_0(fx[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xzz[j], r_0_0[j]);

                t_yyy_yyy[j] = kinvecfunc::fvec_yyy_yyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_yyy_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_yyy_yyz[j] = kinvecfunc::fvec_yyy_yyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyy_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyy_yzz[j] = kinvecfunc::fvec_yyy_yzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_yyy_zzz[j] = kinvecfunc::fvec_yyy_zzz_s_0(fx[j], pa_y[j], pa_yyy[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yyy_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (7) = (70,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_yyz_xxx, t_yyz_xxy, t_yyz_xxz, t_yyz_xyy, \
                                     t_yyz_xyz, t_yyz_xzz, t_yyz_yyy, t_yyz_yyz, t_yyz_yzz, t_yyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xxx[j] = kinvecfunc::fvec_yyz_xxx_s_0(fx[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_yyz_xxy[j] = kinvecfunc::fvec_yyz_xxy_s_0(fx[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_yyz_xxz[j] = kinvecfunc::fvec_yyz_xxz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_xx[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_xx[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_yyz_xyy[j] = kinvecfunc::fvec_yyz_xyy_s_0(fx[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], r_0_0[j]);

                t_yyz_xyz[j] = kinvecfunc::fvec_yyz_xyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xyz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yyz_xzz[j] = kinvecfunc::fvec_yyz_xzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_yyz_yyy[j] = kinvecfunc::fvec_yyz_yyy_s_0(fx[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_yyz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_yyz_yyz[j] = kinvecfunc::fvec_yyz_yyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyz_yzz[j] = kinvecfunc::fvec_yyz_yzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yyz_zzz[j] = kinvecfunc::fvec_yyz_zzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yyz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (8) = (80,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, \
                                     pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, \
                                     pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_yzz_xxx, t_yzz_xxy, t_yzz_xxz, t_yzz_xyy, \
                                     t_yzz_xyz, t_yzz_xzz, t_yzz_yyy, t_yzz_yyz, t_yzz_yzz, t_yzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xxx[j] = kinvecfunc::fvec_yzz_xxx_s_0(fx[j], pa_y[j], pa_yzz[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_yzz_xxy[j] = kinvecfunc::fvec_yzz_xxy_s_0(fx[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_yzz_xxz[j] = kinvecfunc::fvec_yzz_xxz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_xx[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_xx[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_yzz_xyy[j] = kinvecfunc::fvec_yzz_xyy_s_0(fx[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], r_0_0[j]);

                t_yzz_xyz[j] = kinvecfunc::fvec_yzz_xyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xyz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], r_0_0[j]);

                t_yzz_xzz[j] = kinvecfunc::fvec_yzz_xzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_yzz_yyy[j] = kinvecfunc::fvec_yzz_yyy_s_0(fx[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_yzz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_yzz_yyz[j] = kinvecfunc::fvec_yzz_yyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yzz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yzz_yzz[j] = kinvecfunc::fvec_yzz_yzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yzz_zzz[j] = kinvecfunc::fvec_yzz_zzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yzz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (9) = (90,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, \
                                     pb_xyy, pb_xyz, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, \
                                     pb_zzz, r_0_0, s_0_0, t_zzz_xxx, t_zzz_xxy, t_zzz_xxz, t_zzz_xyy, t_zzz_xyz, \
                                     t_zzz_xzz, t_zzz_yyy, t_zzz_yyz, t_zzz_yzz, t_zzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xxx[j] = kinvecfunc::fvec_zzz_xxx_s_0(fx[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_zzz_xxy[j] = kinvecfunc::fvec_zzz_xxy_s_0(fx[j], pa_z[j], pa_zzz[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_zzz_xxz[j] = kinvecfunc::fvec_zzz_xxz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xx[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xx[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_zzz_xyy[j] = kinvecfunc::fvec_zzz_xyy_s_0(fx[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xyy[j], r_0_0[j]);

                t_zzz_xyz[j] = kinvecfunc::fvec_zzz_xyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xy[j], pb_xyz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xyz_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xy[j], pb_xyz[j], r_0_0[j]);

                t_zzz_xzz[j] = kinvecfunc::fvec_zzz_xzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_zzz_yyy[j] = kinvecfunc::fvec_zzz_yyy_s_0(fx[j], pa_z[j], pa_zzz[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_zzz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zzz[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_zzz_yyz[j] = kinvecfunc::fvec_zzz_yyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_zzz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_zzz_yzz[j] = kinvecfunc::fvec_zzz_yzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_zzz_zzz[j] = kinvecfunc::fvec_zzz_zzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_zzz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            idx++;
        }
    }


} // kinrecfunc namespace

