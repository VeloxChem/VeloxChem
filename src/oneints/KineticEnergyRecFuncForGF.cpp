//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForGF.hpp"

#include "KineticEnergyVecFuncForGF.hpp"

namespace kinrecfunc { // kinrecfunc namespace

    void
    compKineticEnergyForGF(      CMemBlock2D<double>& primBuffer,
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

            auto s_0_0 = auxBuffer.data(2 * idx);

            auto r_0_0 = auxBuffer.data(2 * idx + 1);

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

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, \
                                     pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xxxx_xxx, \
                                     t_xxxx_xxy, t_xxxx_xxz, t_xxxx_xyy, t_xxxx_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xxx[j] = kinvecfunc::fvec_xxxx_xxx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xxxx_xxy[j] = kinvecfunc::fvec_xxxx_xxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxxx_xxz[j] = kinvecfunc::fvec_xxxx_xxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxxx_xyy[j] = kinvecfunc::fvec_xxxx_xyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xxxx_xyz[j] = kinvecfunc::fvec_xxxx_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xyz[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xyz[j], pb_yz[j], r_0_0[j]);
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxxx_xzz, t_xxxx_yyy, t_xxxx_yyz, \
                                     t_xxxx_yzz, t_xxxx_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xzz[j] = kinvecfunc::fvec_xxxx_xzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xxxx_yyy[j] = kinvecfunc::fvec_xxxx_yyy_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxxx[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xxxx_yyz[j] = kinvecfunc::fvec_xxxx_yyz_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxxx[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xxxx_yzz[j] = kinvecfunc::fvec_xxxx_yzz_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_y[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxxx[j], pb_y[j], pb_yzz[j], r_0_0[j]);

                t_xxxx_zzz[j] = kinvecfunc::fvec_xxxx_zzz_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxxx_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxxx[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_xxxy_xxx, t_xxxy_xxy, t_xxxy_xxz, t_xxxy_xyy, t_xxxy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xxx[j] = kinvecfunc::fvec_xxxy_xxx_s_0(fx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xxxy_xxy[j] = kinvecfunc::fvec_xxxy_xxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxxy_xxz[j] = kinvecfunc::fvec_xxxy_xxz_s_0(fx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxxy_xyy[j] = kinvecfunc::fvec_xxxy_xyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xxxy_xyz[j] = kinvecfunc::fvec_xxxy_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxxy_xzz, \
                                     t_xxxy_yyy, t_xxxy_yyz, t_xxxy_yzz, t_xxxy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xzz[j] = kinvecfunc::fvec_xxxy_xzz_s_0(fx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xxxy_yyy[j] = kinvecfunc::fvec_xxxy_yyy_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xxxy_yyz[j] = kinvecfunc::fvec_xxxy_yyz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxxy_yzz[j] = kinvecfunc::fvec_xxxy_yzz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xxxy_zzz[j] = kinvecfunc::fvec_xxxy_zzz_s_0(fx[j], pa_xxxy[j], pa_xy[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxxy_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxy[j], pa_xy[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_xxxz_xxx, t_xxxz_xxy, t_xxxz_xxz, t_xxxz_xyy, t_xxxz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xxx[j] = kinvecfunc::fvec_xxxz_xxx_s_0(fx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xxxz_xxy[j] = kinvecfunc::fvec_xxxz_xxy_s_0(fx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxxz_xxz[j] = kinvecfunc::fvec_xxxz_xxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxxz_xyy[j] = kinvecfunc::fvec_xxxz_xyy_s_0(fx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xxxz_xyz[j] = kinvecfunc::fvec_xxxz_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], r_0_0[j]);
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_xxxz_xzz, t_xxxz_yyy, t_xxxz_yyz, t_xxxz_yzz, t_xxxz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xzz[j] = kinvecfunc::fvec_xxxz_xzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xxxz_yyy[j] = kinvecfunc::fvec_xxxz_yyy_s_0(fx[j], pa_xxxz[j], pa_xz[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxxz[j], pa_xz[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xxxz_yyz[j] = kinvecfunc::fvec_xxxz_yyz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xxxz_yzz[j] = kinvecfunc::fvec_xxxz_yzz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_xxxz_zzz[j] = kinvecfunc::fvec_xxxz_zzz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxxz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, \
                                     pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     r_0_0, s_0_0, t_xxyy_xxx, t_xxyy_xxy, t_xxyy_xxz, t_xxyy_xyy, t_xxyy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xxx[j] = kinvecfunc::fvec_xxyy_xxx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xxyy_xxy[j] = kinvecfunc::fvec_xxyy_xxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxyy_xxz[j] = kinvecfunc::fvec_xxyy_xxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxyy_xyy[j] = kinvecfunc::fvec_xxyy_xyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xxyy_xyz[j] = kinvecfunc::fvec_xxyy_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xyy, pa_y, pa_yy, pb_x, pb_xzz, \
                                     pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_xxyy_xzz, t_xxyy_yyy, t_xxyy_yyz, t_xxyy_yzz, t_xxyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xzz[j] = kinvecfunc::fvec_xxyy_xzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xxyy_yyy[j] = kinvecfunc::fvec_xxyy_yyy_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xxyy_yyz[j] = kinvecfunc::fvec_xxyy_yyz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxyy_yzz[j] = kinvecfunc::fvec_xxyy_yzz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xxyy_zzz[j] = kinvecfunc::fvec_xxyy_zzz_s_0(fx[j], pa_xx[j], pa_xxyy[j], pa_yy[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxyy_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxyy[j], pa_yy[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, \
                                     pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xxyz_xxx, t_xxyz_xxy, t_xxyz_xxz, t_xxyz_xyy, \
                                     t_xxyz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xxx[j] = kinvecfunc::fvec_xxyz_xxx_s_0(fx[j], pa_xxyz[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxyz[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xxyz_xxy[j] = kinvecfunc::fvec_xxyz_xxy_s_0(fx[j], pa_xxyz[j], pa_xxz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxyz[j], pa_xxz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxyz_xxz[j] = kinvecfunc::fvec_xxyz_xxz_s_0(fx[j], pa_xxy[j], pa_xxyz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xxyz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxyz_xyy[j] = kinvecfunc::fvec_xxyz_xyy_s_0(fx[j], pa_xxyz[j], pa_xxz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxyz[j], pa_xxz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xxyz_xyz[j] = kinvecfunc::fvec_xxyz_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_y, pa_yz, pa_z, \
                                     pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_xxyz_xzz, t_xxyz_yyy, t_xxyz_yyz, t_xxyz_yzz, t_xxyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xzz[j] = kinvecfunc::fvec_xxyz_xzz_s_0(fx[j], pa_xxy[j], pa_xxyz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xxyz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xxyz_yyy[j] = kinvecfunc::fvec_xxyz_yyy_s_0(fx[j], pa_xxyz[j], pa_xxz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxyz[j], pa_xxz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xxyz_yyz[j] = kinvecfunc::fvec_xxyz_yyz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xxyz_yzz[j] = kinvecfunc::fvec_xxyz_yzz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xxyz_zzz[j] = kinvecfunc::fvec_xxyz_zzz_s_0(fx[j], pa_xxy[j], pa_xxyz[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxyz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xxy[j], pa_xxyz[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, \
                                     pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     r_0_0, s_0_0, t_xxzz_xxx, t_xxzz_xxy, t_xxzz_xxz, t_xxzz_xyy, t_xxzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xxx[j] = kinvecfunc::fvec_xxzz_xxx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xxzz_xxy[j] = kinvecfunc::fvec_xxzz_xxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xxzz_xxz[j] = kinvecfunc::fvec_xxzz_xxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xxzz_xyy[j] = kinvecfunc::fvec_xxzz_xyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xxzz_xyz[j] = kinvecfunc::fvec_xxzz_xyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], r_0_0[j]);
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, \
                                     pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, \
                                     s_0_0, t_xxzz_xzz, t_xxzz_yyy, t_xxzz_yyz, t_xxzz_yzz, t_xxzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xzz[j] = kinvecfunc::fvec_xxzz_xzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xxzz_yyy[j] = kinvecfunc::fvec_xxzz_yyy_s_0(fx[j], pa_xx[j], pa_xxzz[j], pa_zz[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxzz[j], pa_zz[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xxzz_yyz[j] = kinvecfunc::fvec_xxzz_yyz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xxzz_yzz[j] = kinvecfunc::fvec_xxzz_yzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_xxzz_zzz[j] = kinvecfunc::fvec_xxzz_zzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xxzz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (12) = (60,65)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_xyyy_xxx, t_xyyy_xxy, t_xyyy_xxz, t_xyyy_xyy, t_xyyy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xxx[j] = kinvecfunc::fvec_xyyy_xxx_s_0(fx[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xyyy_xxy[j] = kinvecfunc::fvec_xyyy_xxy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyyy_xxz[j] = kinvecfunc::fvec_xyyy_xxz_s_0(fx[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyyy_xyy[j] = kinvecfunc::fvec_xyyy_xyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyyy_xyz[j] = kinvecfunc::fvec_xyyy_xyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (13) = (65,70)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yyy, pb_x, pb_xzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyyy_xzz, \
                                     t_xyyy_yyy, t_xyyy_yyz, t_xyyy_yzz, t_xyyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xzz[j] = kinvecfunc::fvec_xyyy_xzz_s_0(fx[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xzz[j], pb_zz[j], r_0_0[j]);

                t_xyyy_yyy[j] = kinvecfunc::fvec_xyyy_yyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xyyy_yyz[j] = kinvecfunc::fvec_xyyy_yyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyyy_yzz[j] = kinvecfunc::fvec_xyyy_yzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_xyyy_zzz[j] = kinvecfunc::fvec_xyyy_zzz_s_0(fx[j], pa_xy[j], pa_xyyy[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xyyy_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyyy[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (14) = (70,75)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_y, pa_yy, pa_yyz, \
                                     pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xyyz_xxx, t_xyyz_xxy, t_xyyz_xxz, t_xyyz_xyy, \
                                     t_xyyz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xxx[j] = kinvecfunc::fvec_xyyz_xxx_s_0(fx[j], pa_xyyz[j], pa_xz[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xyyz[j], pa_xz[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xyyz_xxy[j] = kinvecfunc::fvec_xyyz_xxy_s_0(fx[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyyz_xxz[j] = kinvecfunc::fvec_xyyz_xxz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyyz_xyy[j] = kinvecfunc::fvec_xyyz_xyy_s_0(fx[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyyz_xyz[j] = kinvecfunc::fvec_xyyz_xyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (15) = (75,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_yy, pa_yyz, pa_z, \
                                     pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_xyyz_xzz, t_xyyz_yyy, t_xyyz_yyz, t_xyyz_yzz, t_xyyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xzz[j] = kinvecfunc::fvec_xyyz_xzz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyyz_yyy[j] = kinvecfunc::fvec_xyyz_yyy_s_0(fx[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xyyz_yyz[j] = kinvecfunc::fvec_xyyz_yyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyyz_yzz[j] = kinvecfunc::fvec_xyyz_yzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyyz_zzz[j] = kinvecfunc::fvec_xyyz_zzz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xyyz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (16) = (80,85)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, \
                                     pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xyzz_xxx, t_xyzz_xxy, t_xyzz_xxz, t_xyzz_xyy, \
                                     t_xyzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xxx[j] = kinvecfunc::fvec_xyzz_xxx_s_0(fx[j], pa_xy[j], pa_xyzz[j], pa_y[j], pa_yzz[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyzz[j], pa_y[j], pa_yzz[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xyzz_xxy[j] = kinvecfunc::fvec_xyzz_xxy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xyzz_xxz[j] = kinvecfunc::fvec_xyzz_xxz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xyzz_xyy[j] = kinvecfunc::fvec_xyzz_xyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], r_0_0[j]);

                t_xyzz_xyz[j] = kinvecfunc::fvec_xyzz_xyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], r_0_0[j]);
            }

            // Batch of Integrals (17) = (85,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, \
                                     pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_xyzz_xzz, t_xyzz_yyy, t_xyzz_yyz, t_xyzz_yzz, t_xyzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xzz[j] = kinvecfunc::fvec_xyzz_xzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyzz_yyy[j] = kinvecfunc::fvec_xyzz_yyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_xyzz_yyz[j] = kinvecfunc::fvec_xyzz_yyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_xyzz_yzz[j] = kinvecfunc::fvec_xyzz_yzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xyzz_zzz[j] = kinvecfunc::fvec_xyzz_zzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xyzz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (18) = (90,95)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_xzzz_xxx, t_xzzz_xxy, t_xzzz_xxz, t_xzzz_xyy, t_xzzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xxx[j] = kinvecfunc::fvec_xzzz_xxx_s_0(fx[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xx[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xx[j], pb_xxx[j], r_0_0[j]);

                t_xzzz_xxy[j] = kinvecfunc::fvec_xzzz_xxy_s_0(fx[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_xxy[j], pb_xy[j], pb_y[j], r_0_0[j]);

                t_xzzz_xxz[j] = kinvecfunc::fvec_xzzz_xxz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xz[j], pb_z[j], r_0_0[j]);

                t_xzzz_xyy[j] = kinvecfunc::fvec_xzzz_xyy_s_0(fx[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xyy[j], pb_yy[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xyy[j], pb_yy[j], r_0_0[j]);

                t_xzzz_xyz[j] = kinvecfunc::fvec_xzzz_xyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_xyz_r_0(fga[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], r_0_0[j]);
            }

            // Batch of Integrals (19) = (95,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_xzzz_xzz, t_xzzz_yyy, t_xzzz_yyz, t_xzzz_yzz, t_xzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xzz[j] = kinvecfunc::fvec_xzzz_xzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_xzzz_yyy[j] = kinvecfunc::fvec_xzzz_yyy_s_0(fx[j], pa_xz[j], pa_xzzz[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_xz[j], pa_xzzz[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_xzzz_yyz[j] = kinvecfunc::fvec_xzzz_yyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_xzzz_yzz[j] = kinvecfunc::fvec_xzzz_yzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_xzzz_zzz[j] = kinvecfunc::fvec_xzzz_zzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_xzzz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (20) = (100,105)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, \
                                     pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, r_0_0, s_0_0, t_yyyy_xxx, t_yyyy_xxy, \
                                     t_yyyy_xxz, t_yyyy_xyy, t_yyyy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xxx[j] = kinvecfunc::fvec_yyyy_xxx_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyyy[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_yyyy_xxy[j] = kinvecfunc::fvec_yyyy_xxy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_xx[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_xx[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_yyyy_xxz[j] = kinvecfunc::fvec_yyyy_xxz_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyyy[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_yyyy_xyy[j] = kinvecfunc::fvec_yyyy_xyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_x[j], pb_xy[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_x[j], pb_xy[j], pb_xyy[j], r_0_0[j]);

                t_yyyy_xyz[j] = kinvecfunc::fvec_yyyy_xyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_xyz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_xyz[j], pb_xz[j], r_0_0[j]);
            }

            // Batch of Integrals (21) = (105,110)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_yyyy_xzz, t_yyyy_yyy, \
                                     t_yyyy_yyz, t_yyyy_yzz, t_yyyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xzz[j] = kinvecfunc::fvec_yyyy_xzz_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_x[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyyy[j], pb_x[j], pb_xzz[j], r_0_0[j]);

                t_yyyy_yyy[j] = kinvecfunc::fvec_yyyy_yyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_yyyy_yyz[j] = kinvecfunc::fvec_yyyy_yyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyyy_yzz[j] = kinvecfunc::fvec_yyyy_yzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_y[j], pb_yzz[j], pb_zz[j], r_0_0[j]);

                t_yyyy_zzz[j] = kinvecfunc::fvec_yyyy_zzz_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_z[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yyyy_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyyy[j], pb_z[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (22) = (110,115)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, r_0_0, s_0_0, \
                                     t_yyyz_xxx, t_yyyz_xxy, t_yyyz_xxz, t_yyyz_xyy, t_yyyz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xxx[j] = kinvecfunc::fvec_yyyz_xxx_s_0(fx[j], pa_yyyz[j], pa_yz[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyyz[j], pa_yz[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_yyyz_xxy[j] = kinvecfunc::fvec_yyyz_xxy_s_0(fx[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_yyyz_xxz[j] = kinvecfunc::fvec_yyyz_xxz_s_0(fx[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_xx[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_xx[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_yyyz_xyy[j] = kinvecfunc::fvec_yyyz_xyy_s_0(fx[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], r_0_0[j]);

                t_yyyz_xyz[j] = kinvecfunc::fvec_yyyz_xyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_xyz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], r_0_0[j]);
            }

            // Batch of Integrals (23) = (115,120)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_yyyz_xzz, t_yyyz_yyy, t_yyyz_yyz, t_yyyz_yzz, t_yyyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xzz[j] = kinvecfunc::fvec_yyyz_xzz_s_0(fx[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_yyyz_yyy[j] = kinvecfunc::fvec_yyyz_yyy_s_0(fx[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_yyyz_yyz[j] = kinvecfunc::fvec_yyyz_yyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyyz_yzz[j] = kinvecfunc::fvec_yyyz_yzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yyyz_zzz[j] = kinvecfunc::fvec_yyyz_zzz_s_0(fx[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yyyz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (24) = (120,125)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, \
                                     pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, r_0_0, s_0_0, \
                                     t_yyzz_xxx, t_yyzz_xxy, t_yyzz_xxz, t_yyzz_xyy, t_yyzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xxx[j] = kinvecfunc::fvec_yyzz_xxx_s_0(fx[j], pa_yy[j], pa_yyzz[j], pa_zz[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyzz[j], pa_zz[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_yyzz_xxy[j] = kinvecfunc::fvec_yyzz_xxy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_yyzz_xxz[j] = kinvecfunc::fvec_yyzz_xxz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_yyzz_xyy[j] = kinvecfunc::fvec_yyzz_xyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], r_0_0[j]);

                t_yyzz_xyz[j] = kinvecfunc::fvec_yyzz_xyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_xyz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], r_0_0[j]);
            }

            // Batch of Integrals (25) = (125,130)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, \
                                     pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, \
                                     s_0_0, t_yyzz_xzz, t_yyzz_yyy, t_yyzz_yyz, t_yyzz_yzz, t_yyzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xzz[j] = kinvecfunc::fvec_yyzz_xzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_yyzz_yyy[j] = kinvecfunc::fvec_yyzz_yyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_yyzz_yyz[j] = kinvecfunc::fvec_yyzz_yyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yyzz_yzz[j] = kinvecfunc::fvec_yyzz_yzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yyzz_zzz[j] = kinvecfunc::fvec_yyzz_zzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yyzz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (26) = (130,135)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, r_0_0, s_0_0, \
                                     t_yzzz_xxx, t_yzzz_xxy, t_yzzz_xxz, t_yzzz_xyy, t_yzzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xxx[j] = kinvecfunc::fvec_yzzz_xxx_s_0(fx[j], pa_yz[j], pa_yzzz[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yz[j], pa_yzzz[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_yzzz_xxy[j] = kinvecfunc::fvec_yzzz_xxy_s_0(fx[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_xx[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_xx[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_yzzz_xxz[j] = kinvecfunc::fvec_yzzz_xxz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_xx[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_xx[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_yzzz_xyy[j] = kinvecfunc::fvec_yzzz_xyy_s_0(fx[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xy[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xy[j], pb_xyy[j], r_0_0[j]);

                t_yzzz_xyz[j] = kinvecfunc::fvec_yzzz_xyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_xyz_r_0(fga[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xz[j], r_0_0[j]);
            }

            // Batch of Integrals (27) = (135,140)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_yzzz_xzz, t_yzzz_yyy, t_yzzz_yyz, t_yzzz_yzz, t_yzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xzz[j] = kinvecfunc::fvec_yzzz_xzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_x[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_x[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_yzzz_yyy[j] = kinvecfunc::fvec_yzzz_yyy_s_0(fx[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_y[j], pb_yy[j], pb_yyy[j], r_0_0[j]);

                t_yzzz_yyz[j] = kinvecfunc::fvec_yzzz_yyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], r_0_0[j]);

                t_yzzz_yzz[j] = kinvecfunc::fvec_yzzz_yzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], r_0_0[j]);

                t_yzzz_zzz[j] = kinvecfunc::fvec_yzzz_zzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_yzzz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            // Batch of Integrals (28) = (140,145)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, \
                                     pb_xy, pb_xyy, pb_xyz, pb_y, pb_z, r_0_0, s_0_0, t_zzzz_xxx, t_zzzz_xxy, t_zzzz_xxz, \
                                     t_zzzz_xyy, t_zzzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xxx[j] = kinvecfunc::fvec_zzzz_xxx_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_x[j], pb_xxx[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_xxx_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pa_zzzz[j], pb_x[j], pb_xxx[j], r_0_0[j]);

                t_zzzz_xxy[j] = kinvecfunc::fvec_zzzz_xxy_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_xxy[j], pb_y[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_xxy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pa_zzzz[j], pb_xxy[j], pb_y[j], r_0_0[j]);

                t_zzzz_xxz[j] = kinvecfunc::fvec_zzzz_xxz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_xx[j], pb_xxz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_xxz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_xx[j], pb_xxz[j], pb_z[j], r_0_0[j]);

                t_zzzz_xyy[j] = kinvecfunc::fvec_zzzz_xyy_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_x[j], pb_xyy[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_xyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pa_zzzz[j], pb_x[j], pb_xyy[j], r_0_0[j]);

                t_zzzz_xyz[j] = kinvecfunc::fvec_zzzz_xyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_xy[j], pb_xyz[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_xyz_r_0(fga[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_xy[j], pb_xyz[j], r_0_0[j]);
            }

            // Batch of Integrals (29) = (145,150)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_zzzz_xzz, \
                                     t_zzzz_yyy, t_zzzz_yyz, t_zzzz_yzz, t_zzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xzz[j] = kinvecfunc::fvec_zzzz_xzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_x[j], pb_xz[j], pb_xzz[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_xzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_x[j], pb_xz[j], pb_xzz[j], r_0_0[j]);

                t_zzzz_yyy[j] = kinvecfunc::fvec_zzzz_yyy_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_y[j], pb_yyy[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_yyy_r_0(fga[j], fgb[j], fx[j], fz[j], pa_zz[j], pa_zzzz[j], pb_y[j], pb_yyy[j], r_0_0[j]);

                t_zzzz_yyz[j] = kinvecfunc::fvec_zzzz_yyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_yyz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_yy[j], pb_yyz[j], pb_z[j], r_0_0[j]);

                t_zzzz_yzz[j] = kinvecfunc::fvec_zzzz_yzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_yzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_y[j], pb_yz[j], pb_yzz[j], r_0_0[j]);

                t_zzzz_zzz[j] = kinvecfunc::fvec_zzzz_zzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]) + kinvecfunc::fvec_zzzz_zzz_r_0(fga[j], fgb[j], fx[j], fz[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], r_0_0[j]);
            }

            idx++;
        }
    }


} // kinrecfunc namespace

