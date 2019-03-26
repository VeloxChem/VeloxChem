//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForGF.hpp"

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
                t_xxxx_xxx[j] = (7.5 * pa_x[j] * fx[j] * fx[j] * fx[j] + 5.625 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                3.0 * pa_xxx[j] * fx[j] * fx[j] + 13.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 

                                9.0 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xxxx[j] * pb_x[j] * fx[j] + 6.0 * pa_xxx[j] * fx[j] * pb_xx[j] + 

                                0.75 * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxx[j] + pa_xxxx[j] * pb_xxx[j]) * s_0_0[j] + (60.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                9.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 9.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                13.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 6.0 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                45.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 30.0 * pa_xxx[j] * fz[j] * fx[j] * fx[j] + 

                                135.0 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 2.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                9.0 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 9.0 * pa_xx[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                18.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] - 3.0 * pa_xxxx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                90.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 18.0 * pa_xxxx[j] * fz[j] * pb_x[j] * fx[j] + 

                                72.0 * pa_xxx[j] * fz[j] * fx[j] * pb_xx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xxx[j] - 

                                6.0 * pa_xx[j] * fz[j] * fga[j] * pb_xxx[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 

                                36.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxx[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xxxx_xxy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                6.0 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xxxx[j] * fx[j] * pb_y[j] + 4.0 * pa_xxx[j] * fx[j] * pb_xy[j] + 

                                0.75 * fx[j] * fx[j] * pb_xxy[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxy[j] + pa_xxxx[j] * pb_xxy[j]) * s_0_0[j] + (-4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] + 

                                15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 45.0 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 3.0 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                3.0 * pa_xx[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 12.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] - 

                                pa_xxxx[j] * fz[j] * fgb[j] * pb_y[j] + 60.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                6.0 * pa_xxxx[j] * fz[j] * fx[j] * pb_y[j] + 48.0 * pa_xxx[j] * fz[j] * fx[j] * pb_xy[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xxy[j] - 6.0 * pa_xx[j] * fz[j] * fga[j] * pb_xxy[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 36.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxy[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xxxx_xxz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xxxx[j] * fx[j] * pb_z[j] + 4.0 * pa_xxx[j] * fx[j] * pb_xz[j] + 

                                0.75 * fx[j] * fx[j] * pb_xxz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xxz[j] + pa_xxxx[j] * pb_xxz[j]) * s_0_0[j] + (-4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] + 

                                15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 45.0 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 3.0 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                3.0 * pa_xx[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 12.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] - 

                                pa_xxxx[j] * fz[j] * fgb[j] * pb_z[j] + 60.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                6.0 * pa_xxxx[j] * fz[j] * fx[j] * pb_z[j] + 48.0 * pa_xxx[j] * fz[j] * fx[j] * pb_xz[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xxz[j] - 6.0 * pa_xx[j] * fz[j] * fga[j] * pb_xxz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 36.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxz[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xxxx_xyy[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] + pa_xxx[j] * fx[j] * fx[j] + 

                                0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 

                                0.5 * pa_xxxx[j] * pb_x[j] * fx[j] + 2.0 * pa_xxx[j] * fx[j] * pb_yy[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] + 

                                3.0 * pa_xx[j] * fx[j] * pb_xyy[j] + pa_xxxx[j] * pb_xyy[j]) * s_0_0[j] + (-3.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 2.0 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                12.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 10.0 * pa_xxx[j] * fz[j] * fx[j] * fx[j] - 

                                0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                3.0 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 3.0 * pa_xx[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                6.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                pa_xxxx[j] * pb_x[j] * fz[j] * fgb[j] + 15.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                30.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 6.0 * pa_xxxx[j] * fz[j] * pb_x[j] * fx[j] + 

                                24.0 * pa_xxx[j] * fz[j] * fx[j] * pb_yy[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xyy[j] - 

                                6.0 * pa_xx[j] * fz[j] * fga[j] * pb_xyy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 

                                36.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyy[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xxxx_xyz[j] = (3.0 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 2.0 * pa_xxx[j] * fx[j] * pb_yz[j] + 

                                0.75 * fx[j] * fx[j] * pb_xyz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xyz[j] + pa_xxxx[j] * pb_xyz[j]) * s_0_0[j] + (-6.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] + 

                                30.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 24.0 * pa_xxx[j] * fz[j] * fx[j] * pb_yz[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xyz[j] - 6.0 * pa_xx[j] * fz[j] * fga[j] * pb_xyz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 36.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyz[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xzz, pb_y, pb_yyy, pb_yyz, \
                                     pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxxx_xzz, t_xxxx_yyy, t_xxxx_yyz, \
                                     t_xxxx_yzz, t_xxxx_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xzz[j] = (1.5 * pa_x[j] * fx[j] * fx[j] * fx[j] + pa_xxx[j] * fx[j] * fx[j] + 

                                0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 

                                0.5 * pa_xxxx[j] * pb_x[j] * fx[j] + 2.0 * pa_xxx[j] * fx[j] * pb_zz[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] + 

                                3.0 * pa_xx[j] * fx[j] * pb_xzz[j] + pa_xxxx[j] * pb_xzz[j]) * s_0_0[j] + (-3.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 2.0 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                12.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 10.0 * pa_xxx[j] * fz[j] * fx[j] * fx[j] - 

                                0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                3.0 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 3.0 * pa_xx[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                6.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                pa_xxxx[j] * pb_x[j] * fz[j] * fgb[j] + 15.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                30.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 6.0 * pa_xxxx[j] * fz[j] * pb_x[j] * fx[j] + 

                                24.0 * pa_xxx[j] * fz[j] * fx[j] * pb_zz[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xzz[j] - 

                                6.0 * pa_xx[j] * fz[j] * fga[j] * pb_xzz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 

                                36.0 * pa_xx[j] * fz[j] * fx[j] * pb_xzz[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xxxx_yyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_xxxx[j] * pb_y[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_xx[j] * fx[j] * pb_yyy[j] + 

                                pa_xxxx[j] * pb_yyy[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 9.0 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                9.0 * pa_xx[j] * fz[j] * fga[j] * pb_y[j] * fx[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                3.0 * pa_xxxx[j] * pb_y[j] * fz[j] * fgb[j] + 45.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                18.0 * pa_xxxx[j] * fz[j] * pb_y[j] * fx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_yyy[j] - 

                                6.0 * pa_xx[j] * fz[j] * fga[j] * pb_yyy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 

                                36.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyy[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xxxx_yyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_xxxx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] + 3.0 * pa_xx[j] * fx[j] * pb_yyz[j] + 

                                pa_xxxx[j] * pb_yyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 3.0 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                3.0 * pa_xx[j] * fz[j] * fga[j] * fx[j] * pb_z[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                pa_xxxx[j] * fz[j] * fgb[j] * pb_z[j] + 15.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_xxxx[j] * fz[j] * fx[j] * pb_z[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_yyz[j] - 

                                6.0 * pa_xx[j] * fz[j] * fga[j] * pb_yyz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 

                                36.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyz[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xxxx_yzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.5 * pa_xxxx[j] * pb_y[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] + 3.0 * pa_xx[j] * fx[j] * pb_yzz[j] + 

                                pa_xxxx[j] * pb_yzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 3.0 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                3.0 * pa_xx[j] * fz[j] * fga[j] * pb_y[j] * fx[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                pa_xxxx[j] * pb_y[j] * fz[j] * fgb[j] + 15.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                6.0 * pa_xxxx[j] * fz[j] * pb_y[j] * fx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_yzz[j] - 

                                6.0 * pa_xx[j] * fz[j] * fga[j] * pb_yzz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 

                                36.0 * pa_xx[j] * fz[j] * fx[j] * pb_yzz[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xxxx_zzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_xxxx[j] * pb_z[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_xx[j] * fx[j] * pb_zzz[j] + 

                                pa_xxxx[j] * pb_zzz[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 9.0 * pa_xx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                9.0 * pa_xx[j] * fz[j] * fga[j] * pb_z[j] * fx[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                3.0 * pa_xxxx[j] * pb_z[j] * fz[j] * fgb[j] + 45.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                18.0 * pa_xxxx[j] * fz[j] * pb_z[j] * fx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_zzz[j] - 

                                6.0 * pa_xx[j] * fz[j] * fga[j] * pb_zzz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 

                                36.0 * pa_xx[j] * fz[j] * fx[j] * pb_zzz[j] + 14.0 * pa_xxxx[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_xxxy_xxx, t_xxxy_xxy, t_xxxy_xxz, t_xxxy_xyy, t_xxxy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xxx[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_y[j] + 2.25 * pa_xxy[j] * fx[j] * fx[j] + 

                                6.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 1.5 * pa_xxxy[j] * pb_x[j] * fx[j] + 

                                4.5 * pa_xxy[j] * fx[j] * pb_xx[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxx[j] + pa_xxxy[j] * pb_xxx[j]) * s_0_0[j] + (15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_y[j] - 

                                2.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pa_y[j] - 

                                4.5 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] + 22.5 * pa_xxy[j] * fx[j] * fx[j] * fz[j] + 

                                67.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 4.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                4.5 * pa_xy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 4.5 * fx[j] * fz[j] * fga[j] * pa_y[j] * pb_xx[j] - 

                                3.0 * pa_xxxy[j] * pb_x[j] * fz[j] * fgb[j] + 22.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xx[j] + 

                                18.0 * pa_xxxy[j] * fz[j] * pb_x[j] * fx[j] + 54.0 * pa_xxy[j] * fx[j] * fz[j] * pb_xx[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xxx[j] + 18.0 * pa_xy[j] * fx[j] * fz[j] * pb_xxx[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xxxy_xxy[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xxx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.25 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.5 * pa_xxxy[j] * fx[j] * pb_y[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xx[j] + 

                                3.0 * pa_xxy[j] * fx[j] * pb_xy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxy[j] + pa_xxxy[j] * pb_xxy[j]) * s_0_0[j] + (9.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 0.5 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 2.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] + 

                                15.0 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 22.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                1.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_xx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pa_y[j] * pb_xy[j] - 

                                pa_xxxy[j] * fz[j] * fgb[j] * pb_y[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xy[j] + 6.0 * pa_xxxy[j] * fz[j] * fx[j] * pb_y[j] + 

                                6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_xx[j] + 36.0 * pa_xxy[j] * fx[j] * fz[j] * pb_xy[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xxy[j] + 18.0 * pa_xy[j] * fx[j] * fz[j] * pb_xxy[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xxxy_xxz[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xxxy[j] * fx[j] * pb_z[j] + 3.0 * pa_xxy[j] * fx[j] * pb_xz[j] + 

                                1.5 * pa_xy[j] * fx[j] * pb_xxz[j] + pa_xxxy[j] * pb_xxz[j]) * s_0_0[j] + (22.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pa_y[j] * pb_xz[j] - pa_xxxy[j] * fz[j] * fgb[j] * pb_z[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xz[j] + 6.0 * pa_xxxy[j] * fz[j] * fx[j] * pb_z[j] + 

                                36.0 * pa_xxy[j] * fx[j] * fz[j] * pb_xz[j] - 3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xxz[j] + 

                                18.0 * pa_xy[j] * fx[j] * fz[j] * pb_xxz[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xxxy_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_xxy[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.5 * pa_xxxy[j] * pb_x[j] * fx[j] + pa_xxx[j] * fx[j] * pb_xy[j] + 

                                1.5 * pa_xxy[j] * fx[j] * pb_yy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyy[j] + pa_xxxy[j] * pb_xyy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pa_y[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                1.5 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_y[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_xxy[j] * fx[j] * fx[j] * fz[j] + 

                                15.0 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 1.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 3.0 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_xy[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pa_y[j] * pb_yy[j] - pa_xxxy[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_yy[j] + 6.0 * pa_xxxy[j] * fz[j] * pb_x[j] * fx[j] + 

                                12.0 * pa_xxx[j] * fz[j] * fx[j] * pb_xy[j] + 18.0 * pa_xxy[j] * fx[j] * fz[j] * pb_yy[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xyy[j] + 18.0 * pa_xy[j] * fx[j] * fz[j] * pb_xyy[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xxxy_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xz[j] + 

                                1.5 * pa_xxy[j] * fx[j] * pb_yz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyz[j] + pa_xxxy[j] * pb_xyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_xz[j] - 1.5 * fx[j] * fz[j] * fga[j] * pa_y[j] * pb_yz[j] + 

                                7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_yz[j] + 

                                6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_xz[j] + 18.0 * pa_xxy[j] * fx[j] * fz[j] * pb_yz[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xyz[j] + 18.0 * pa_xy[j] * fx[j] * fz[j] * pb_xyz[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxxy_xzz, \
                                     t_xxxy_yyy, t_xxxy_yyz, t_xxxy_yzz, t_xxxy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] + 

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.5 * pa_xxxy[j] * pb_x[j] * fx[j] + 

                                1.5 * pa_xxy[j] * fx[j] * pb_zz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xzz[j] + pa_xxxy[j] * pb_xzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pa_y[j] - 1.5 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_y[j] + 7.5 * pa_xxy[j] * fx[j] * fx[j] * fz[j] - 

                                1.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pa_y[j] * pb_zz[j] - pa_xxxy[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_zz[j] + 

                                6.0 * pa_xxxy[j] * fz[j] * pb_x[j] * fx[j] + 18.0 * pa_xxy[j] * fx[j] * fz[j] * pb_zz[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xzz[j] + 18.0 * pa_xy[j] * fx[j] * fz[j] * pb_xzz[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xxxy_yyy[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] + 

                                2.25 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xxxy[j] * pb_y[j] * fx[j] + 

                                1.5 * pa_xxx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xy[j] * fx[j] * pb_yyy[j] + pa_xxxy[j] * pb_yyy[j]) * s_0_0[j] + (-2.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                2.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 1.5 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 7.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] - 

                                4.5 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 4.5 * pa_xy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                4.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - 3.0 * pa_xxxy[j] * pb_y[j] * fz[j] * fgb[j] + 

                                22.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                18.0 * pa_xxxy[j] * fz[j] * pb_y[j] * fx[j] + 18.0 * pa_xxx[j] * fz[j] * fx[j] * pb_yy[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_yyy[j] + 18.0 * pa_xy[j] * fx[j] * fz[j] * pb_yyy[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xxxy_yyz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xxxy[j] * fx[j] * pb_z[j] + pa_xxx[j] * fx[j] * pb_yz[j] + 

                                1.5 * pa_xy[j] * fx[j] * pb_yyz[j] + pa_xxxy[j] * pb_yyz[j]) * s_0_0[j] + (-1.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                1.5 * pa_xy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 3.0 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - 

                                pa_xxxy[j] * fz[j] * fgb[j] * pb_z[j] + 7.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 6.0 * pa_xxxy[j] * fz[j] * fx[j] * pb_z[j] + 

                                12.0 * pa_xxx[j] * fz[j] * fx[j] * pb_yz[j] - 3.0 * pa_xy[j] * fz[j] * fga[j] * pb_yyz[j] + 

                                18.0 * pa_xy[j] * fx[j] * fz[j] * pb_yyz[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xxxy_yzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] + 

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xxxy[j] * pb_y[j] * fx[j] + 

                                0.5 * pa_xxx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xy[j] * fx[j] * pb_yzz[j] + pa_xxxy[j] * pb_yzz[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 2.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] - 

                                1.5 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - pa_xxxy[j] * pb_y[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                6.0 * pa_xxxy[j] * fz[j] * pb_y[j] * fx[j] + 6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_zz[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_yzz[j] + 18.0 * pa_xy[j] * fx[j] * fz[j] * pb_yzz[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xxxy_zzz[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xxxy[j] * pb_z[j] * fx[j] + 

                                1.5 * pa_xy[j] * fx[j] * pb_zzz[j] + pa_xxxy[j] * pb_zzz[j]) * s_0_0[j] + (-4.5 * pa_xy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                4.5 * pa_xy[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 3.0 * pa_xxxy[j] * pb_z[j] * fz[j] * fgb[j] + 

                                22.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 18.0 * pa_xxxy[j] * fz[j] * pb_z[j] * fx[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_zzz[j] + 18.0 * pa_xy[j] * fx[j] * fz[j] * pb_zzz[j] + 14.0 * pa_xxxy[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_xxxz_xxx, t_xxxz_xxy, t_xxxz_xxz, t_xxxz_xyy, t_xxxz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xxx[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] + 2.25 * pa_xxz[j] * fx[j] * fx[j] + 

                                6.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 1.5 * pa_xxxz[j] * pb_x[j] * fx[j] + 

                                4.5 * pa_xxz[j] * fx[j] * pb_xx[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxx[j] + pa_xxxz[j] * pb_xxx[j]) * s_0_0[j] + (15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 

                                2.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pa_z[j] - 

                                4.5 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] + 22.5 * pa_xxz[j] * fx[j] * fx[j] * fz[j] + 

                                67.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 4.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                4.5 * pa_xz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 4.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_xx[j] - 

                                3.0 * pa_xxxz[j] * pb_x[j] * fz[j] * fgb[j] + 22.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] + 

                                18.0 * pa_xxxz[j] * fz[j] * pb_x[j] * fx[j] + 54.0 * pa_xxz[j] * fx[j] * fz[j] * pb_xx[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xxx[j] + 18.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxx[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xxxz_xxy[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xxxz[j] * fx[j] * pb_y[j] + 3.0 * pa_xxz[j] * fx[j] * pb_xy[j] + 

                                1.5 * pa_xz[j] * fx[j] * pb_xxy[j] + pa_xxxz[j] * pb_xxy[j]) * s_0_0[j] + (22.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                1.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_xy[j] - pa_xxxz[j] * fz[j] * fgb[j] * pb_y[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] + 6.0 * pa_xxxz[j] * fz[j] * fx[j] * pb_y[j] + 

                                36.0 * pa_xxz[j] * fx[j] * fz[j] * pb_xy[j] - 3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xxy[j] + 

                                18.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxy[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xxxz_xxz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xxx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.25 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_xxxz[j] * fx[j] * pb_z[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xx[j] + 

                                3.0 * pa_xxz[j] * fx[j] * pb_xz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxz[j] + pa_xxxz[j] * pb_xxz[j]) * s_0_0[j] + (9.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 0.5 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 2.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] + 

                                15.0 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 22.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_xx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_xz[j] - 

                                pa_xxxz[j] * fz[j] * fgb[j] * pb_z[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] + 6.0 * pa_xxxz[j] * fz[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_xx[j] + 36.0 * pa_xxz[j] * fx[j] * fz[j] * pb_xz[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xxz[j] + 18.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxz[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xxxz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_xxxz[j] * pb_x[j] * fx[j] + 

                                1.5 * pa_xxz[j] * fx[j] * pb_yy[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyy[j] + pa_xxxz[j] * pb_xyy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pa_z[j] - 1.5 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 7.5 * pa_xxz[j] * fx[j] * fx[j] * fz[j] - 

                                1.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_yy[j] - pa_xxxz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] + 

                                6.0 * pa_xxxz[j] * fz[j] * pb_x[j] * fx[j] + 18.0 * pa_xxz[j] * fx[j] * fz[j] * pb_yy[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xyy[j] + 18.0 * pa_xz[j] * fx[j] * fz[j] * pb_xyy[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xxxz_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xxx[j] * fx[j] * pb_xy[j] + 

                                1.5 * pa_xxz[j] * fx[j] * pb_yz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyz[j] + pa_xxxz[j] * pb_xyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_xy[j] - 1.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_yz[j] + 

                                7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] + 

                                6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_xy[j] + 18.0 * pa_xxz[j] * fx[j] * fz[j] * pb_yz[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xyz[j] + 18.0 * pa_xz[j] * fx[j] * fz[j] * pb_xyz[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_xxxz_xzz, t_xxxz_yyy, t_xxxz_yyz, t_xxxz_yzz, t_xxxz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_xxz[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xxxz[j] * pb_x[j] * fx[j] + pa_xxx[j] * fx[j] * pb_xz[j] + 

                                1.5 * pa_xxz[j] * fx[j] * pb_zz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xzz[j] + pa_xxxz[j] * pb_xzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pa_z[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                1.5 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_xxz[j] * fx[j] * fx[j] * fz[j] + 

                                15.0 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 1.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 3.0 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_xz[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_zz[j] - pa_xxxz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] + 6.0 * pa_xxxz[j] * fz[j] * pb_x[j] * fx[j] + 

                                12.0 * pa_xxx[j] * fz[j] * fx[j] * pb_xz[j] + 18.0 * pa_xxz[j] * fx[j] * fz[j] * pb_zz[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xzz[j] + 18.0 * pa_xz[j] * fx[j] * fz[j] * pb_xzz[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xxxz_yyy[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xxxz[j] * pb_y[j] * fx[j] + 

                                1.5 * pa_xz[j] * fx[j] * pb_yyy[j] + pa_xxxz[j] * pb_yyy[j]) * s_0_0[j] + (-4.5 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                4.5 * pa_xz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 3.0 * pa_xxxz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                22.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 18.0 * pa_xxxz[j] * fz[j] * pb_y[j] * fx[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_yyy[j] + 18.0 * pa_xz[j] * fx[j] * fz[j] * pb_yyy[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xxxz_yyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xxxz[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_xxx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyz[j] + pa_xxxz[j] * pb_yyz[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 2.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] - 

                                1.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - pa_xxxz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                6.0 * pa_xxxz[j] * fz[j] * fx[j] * pb_z[j] + 6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_yy[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_yyz[j] + 18.0 * pa_xz[j] * fx[j] * fz[j] * pb_yyz[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xxxz_yzz[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xxxz[j] * pb_y[j] * fx[j] + pa_xxx[j] * fx[j] * pb_yz[j] + 

                                1.5 * pa_xz[j] * fx[j] * pb_yzz[j] + pa_xxxz[j] * pb_yzz[j]) * s_0_0[j] + (-1.5 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 3.0 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - 

                                pa_xxxz[j] * pb_y[j] * fz[j] * fgb[j] + 7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 6.0 * pa_xxxz[j] * fz[j] * pb_y[j] * fx[j] + 

                                12.0 * pa_xxx[j] * fz[j] * fx[j] * pb_yz[j] - 3.0 * pa_xz[j] * fz[j] * fga[j] * pb_yzz[j] + 

                                18.0 * pa_xz[j] * fx[j] * fz[j] * pb_yzz[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xxxz_zzz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] + 

                                2.25 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xxxz[j] * pb_z[j] * fx[j] + 

                                1.5 * pa_xxx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xz[j] * fx[j] * pb_zzz[j] + pa_xxxz[j] * pb_zzz[j]) * s_0_0[j] + (-2.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                2.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 1.5 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 7.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] - 

                                4.5 * pa_xz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 4.5 * pa_xz[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 

                                4.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - 3.0 * pa_xxxz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                22.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                18.0 * pa_xxxz[j] * fz[j] * pb_z[j] * fx[j] + 18.0 * pa_xxx[j] * fz[j] * fx[j] * pb_zz[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_zzz[j] + 18.0 * pa_xz[j] * fx[j] * fz[j] * pb_zzz[j] + 14.0 * pa_xxxz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, \
                                     pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     r_0_0, s_0_0, t_xxyy_xxx, t_xxyy_xxy, t_xxyy_xxz, t_xxyy_xyy, t_xxyy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xxx[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                1.5 * pa_xyy[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xxyy[j] * pb_x[j] * fx[j] + 

                                3.0 * pa_xyy[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxx[j] + 

                                0.5 * fx[j] * pa_yy[j] * pb_xxx[j] + pa_xxyy[j] * pb_xxx[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 3.0 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                3.0 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 6.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 15.0 * pa_xyy[j] * fx[j] * fx[j] * fz[j] + 

                                22.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_x[j] - 0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xx[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                3.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] - 1.5 * fx[j] * pa_yy[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_yy[j] * pb_x[j] * fx[j] - 3.0 * pa_xxyy[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                18.0 * pa_xxyy[j] * fz[j] * pb_x[j] * fx[j] + 36.0 * pa_xyy[j] * fx[j] * fz[j] * pb_xx[j] - fx[j] * fz[j] * fga[j] * pb_xxx[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_xxx[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] - fz[j] * fga[j] * pa_yy[j] * pb_xxx[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxx[j] + 6.0 * fx[j] * pa_yy[j] * fz[j] * pb_xxx[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xxyy_xxy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.5 * pa_xxy[j] * fx[j] * fx[j] + 2.0 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 

                                0.5 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.5 * pa_xxyy[j] * fx[j] * pb_y[j] + pa_xxy[j] * fx[j] * pb_xx[j] + 

                                2.0 * pa_xyy[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxy[j] + 

                                0.5 * fx[j] * pa_yy[j] * pb_xxy[j] + pa_xxyy[j] * pb_xxy[j]) * s_0_0[j] + (6.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] - 

                                0.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                0.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - pa_xxy[j] * fx[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                5.0 * pa_xxy[j] * fz[j] * fx[j] * fx[j] + 20.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_y[j] - 0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 0.5 * pa_xx[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 

                                2.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] - 0.5 * fx[j] * pa_yy[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.5 * fz[j] * fga[j] * pa_yy[j] * fx[j] * pb_y[j] - fz[j] * fga[j] * pa_y[j] * fx[j] * pb_xx[j] - pa_xxyy[j] * fz[j] * fgb[j] * pb_y[j] + 

                                2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 10.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xx[j] + 6.0 * pa_xxyy[j] * fz[j] * fx[j] * pb_y[j] + 

                                12.0 * pa_xxy[j] * fz[j] * fx[j] * pb_xx[j] + 24.0 * pa_xyy[j] * fx[j] * fz[j] * pb_xy[j] - fx[j] * fz[j] * fga[j] * pb_xxy[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_xxy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] - fz[j] * fga[j] * pa_yy[j] * pb_xxy[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxy[j] + 6.0 * fx[j] * pa_yy[j] * fz[j] * pb_xxy[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xxyy_xxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + 

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_xxyy[j] * fx[j] * pb_z[j] + 

                                2.0 * pa_xyy[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxz[j] + 

                                0.5 * fx[j] * pa_yy[j] * pb_xxz[j] + pa_xxyy[j] * pb_xxz[j]) * s_0_0[j] + (-fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_z[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * pa_xx[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 2.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] - 

                                0.5 * fx[j] * pa_yy[j] * fz[j] * fgb[j] * pb_z[j] - 0.5 * fz[j] * fga[j] * pa_yy[j] * fx[j] * pb_z[j] - 

                                pa_xxyy[j] * fz[j] * fgb[j] * pb_z[j] + 2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                10.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 6.0 * pa_xxyy[j] * fz[j] * fx[j] * pb_z[j] + 

                                24.0 * pa_xyy[j] * fx[j] * fz[j] * pb_xz[j] - fx[j] * fz[j] * fga[j] * pb_xxz[j] - pa_xx[j] * fz[j] * fga[j] * pb_xxz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] - fz[j] * fga[j] * pa_yy[j] * pb_xxz[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxz[j] + 

                                6.0 * fx[j] * pa_yy[j] * fz[j] * pb_xxz[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xxyy_xyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xyy[j] * fx[j] * fx[j] + 2.0 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 

                                0.5 * pa_xxyy[j] * pb_x[j] * fx[j] + 2.0 * pa_xxy[j] * fx[j] * pb_xy[j] + pa_xyy[j] * fx[j] * pb_yy[j] + 

                                0.25 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_yy[j] * pb_xyy[j] + 

                                pa_xxyy[j] * pb_xyy[j]) * s_0_0[j] + (6.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 5.0 * pa_xyy[j] * fx[j] * fx[j] * fz[j] + 

                                20.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.5 * pa_xx[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] - 0.5 * fx[j] * pa_yy[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_yy[j] * pb_x[j] * fx[j] - 2.0 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_xy[j] - 

                                pa_xxyy[j] * pb_x[j] * fz[j] * fgb[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                2.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_x[j] + 10.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xy[j] + 

                                6.0 * pa_xxyy[j] * fz[j] * pb_x[j] * fx[j] + 24.0 * pa_xxy[j] * fz[j] * fx[j] * pb_xy[j] + 

                                12.0 * pa_xyy[j] * fx[j] * fz[j] * pb_yy[j] - fx[j] * fz[j] * fga[j] * pb_xyy[j] - pa_xx[j] * fz[j] * fga[j] * pb_xyy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] - fz[j] * fga[j] * pa_yy[j] * pb_xyy[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyy[j] + 

                                6.0 * fx[j] * pa_yy[j] * fz[j] * pb_xyy[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xxyy_xyz[j] = (pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 

                                0.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + pa_xxy[j] * fx[j] * pb_xz[j] + pa_xyy[j] * fx[j] * pb_yz[j] + 

                                0.25 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yy[j] * pb_xyz[j] + 

                                pa_xxyy[j] * pb_xyz[j]) * s_0_0[j] + (10.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] - fz[j] * fga[j] * pa_y[j] * fx[j] * pb_xz[j] + 

                                5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xz[j] + 

                                12.0 * pa_xxy[j] * fz[j] * fx[j] * pb_xz[j] + 12.0 * pa_xyy[j] * fx[j] * fz[j] * pb_yz[j] - fx[j] * fz[j] * fga[j] * pb_xyz[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_xyz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xyz[j] - fz[j] * fga[j] * pa_yy[j] * pb_xyz[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyz[j] + 6.0 * fx[j] * pa_yy[j] * fz[j] * pb_xyz[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xyy, pa_y, pa_yy, pb_x, pb_xzz, \
                                     pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_xxyy_xzz, t_xxyy_yyy, t_xxyy_yyz, t_xxyy_yzz, t_xxyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * fx[j] + 

                                0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 0.5 * pa_xxyy[j] * pb_x[j] * fx[j] + 

                                pa_xyy[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xzz[j] + 

                                0.5 * fx[j] * pa_yy[j] * pb_xzz[j] + pa_xxyy[j] * pb_xzz[j]) * s_0_0[j] + (-0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 2.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                5.0 * pa_xyy[j] * fx[j] * fx[j] * fz[j] - 0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 0.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * pa_xx[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - pa_x[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] - 

                                0.5 * fx[j] * pa_yy[j] * pb_x[j] * fz[j] * fgb[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.5 * fz[j] * fga[j] * pa_yy[j] * pb_x[j] * fx[j] - pa_xxyy[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                2.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_x[j] + 6.0 * pa_xxyy[j] * fz[j] * pb_x[j] * fx[j] + 

                                12.0 * pa_xyy[j] * fx[j] * fz[j] * pb_zz[j] - fx[j] * fz[j] * fga[j] * pb_xzz[j] - pa_xx[j] * fz[j] * fga[j] * pb_xzz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] - fz[j] * fga[j] * pa_yy[j] * pb_xzz[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xzz[j] + 

                                6.0 * fx[j] * pa_yy[j] * fz[j] * pb_xzz[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xxyy_yyy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_xxy[j] * fx[j] * fx[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 1.5 * pa_xxyy[j] * pb_y[j] * fx[j] + 

                                3.0 * pa_xxy[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyy[j] + 

                                0.5 * fx[j] * pa_yy[j] * pb_yyy[j] + pa_xxyy[j] * pb_yyy[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 3.0 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - 

                                3.0 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] + 6.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 15.0 * pa_xxy[j] * fz[j] * fx[j] * fx[j] + 

                                22.5 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * pa_xx[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                1.5 * fx[j] * pa_yy[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_yy[j] * pb_y[j] * fx[j] - 

                                3.0 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_yy[j] - 3.0 * pa_xxyy[j] * pb_y[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_y[j] + 15.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yy[j] + 

                                18.0 * pa_xxyy[j] * fz[j] * pb_y[j] * fx[j] + 36.0 * pa_xxy[j] * fz[j] * fx[j] * pb_yy[j] - fx[j] * fz[j] * fga[j] * pb_yyy[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_yyy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] - fz[j] * fga[j] * pa_yy[j] * pb_yyy[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyy[j] + 6.0 * fx[j] * pa_yy[j] * fz[j] * pb_yyy[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xxyy_yyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.5 * pa_xxyy[j] * fx[j] * pb_z[j] + 

                                2.0 * pa_xxy[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyz[j] + 

                                0.5 * fx[j] * pa_yy[j] * pb_yyz[j] + pa_xxyy[j] * pb_yyz[j]) * s_0_0[j] + (-fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * pa_xx[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 0.5 * fx[j] * pa_yy[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * fz[j] * fga[j] * pa_yy[j] * fx[j] * pb_z[j] - 2.0 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_yz[j] - 

                                pa_xxyy[j] * fz[j] * fgb[j] * pb_z[j] + 2.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_z[j] + 

                                10.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yz[j] + 6.0 * pa_xxyy[j] * fz[j] * fx[j] * pb_z[j] + 

                                24.0 * pa_xxy[j] * fz[j] * fx[j] * pb_yz[j] - fx[j] * fz[j] * fga[j] * pb_yyz[j] - pa_xx[j] * fz[j] * fga[j] * pb_yyz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] - fz[j] * fga[j] * pa_yy[j] * pb_yyz[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyz[j] + 

                                6.0 * fx[j] * pa_yy[j] * fz[j] * pb_yyz[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xxyy_yzz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.5 * pa_xxy[j] * fx[j] * fx[j] + 

                                0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.5 * pa_xxyy[j] * pb_y[j] * fx[j] + 

                                pa_xxy[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yzz[j] + 

                                0.5 * fx[j] * pa_yy[j] * pb_yzz[j] + pa_xxyy[j] * pb_yzz[j]) * s_0_0[j] + (-0.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - pa_xxy[j] * fx[j] * fz[j] * fgb[j] + 2.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 

                                5.0 * pa_xxy[j] * fz[j] * fx[j] * fx[j] - 0.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 0.5 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.5 * pa_xx[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 0.5 * fx[j] * pa_yy[j] * pb_y[j] * fz[j] * fgb[j] + 

                                fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 0.5 * fz[j] * fga[j] * pa_yy[j] * pb_y[j] * fx[j] - 

                                fz[j] * fga[j] * pa_y[j] * fx[j] * pb_zz[j] - pa_xxyy[j] * pb_y[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 2.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_y[j] + 

                                5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_zz[j] + 6.0 * pa_xxyy[j] * fz[j] * pb_y[j] * fx[j] + 

                                12.0 * pa_xxy[j] * fz[j] * fx[j] * pb_zz[j] - fx[j] * fz[j] * fga[j] * pb_yzz[j] - pa_xx[j] * fz[j] * fga[j] * pb_yzz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] - fz[j] * fga[j] * pa_yy[j] * pb_yzz[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_yzz[j] + 

                                6.0 * fx[j] * pa_yy[j] * fz[j] * pb_yzz[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xxyy_zzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + 1.5 * pa_xxyy[j] * pb_z[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 

                                0.5 * pa_xx[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_yy[j] * pb_zzz[j] + pa_xxyy[j] * pb_zzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 1.5 * pa_xx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xx[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 1.5 * fx[j] * pa_yy[j] * pb_z[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 1.5 * fz[j] * fga[j] * pa_yy[j] * pb_z[j] * fx[j] - 

                                3.0 * pa_xxyy[j] * pb_z[j] * fz[j] * fgb[j] + 7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                7.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_z[j] + 18.0 * pa_xxyy[j] * fz[j] * pb_z[j] * fx[j] - fx[j] * fz[j] * fga[j] * pb_zzz[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_zzz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] - fz[j] * fga[j] * pa_yy[j] * pb_zzz[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_zzz[j] + 6.0 * fx[j] * pa_yy[j] * fz[j] * pb_zzz[j] + 14.0 * pa_xxyy[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, \
                                     pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xxyz_xxx, t_xxyz_xxy, t_xxyz_xxz, t_xxyz_xyy, \
                                     t_xxyz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xxx[j] = (1.5 * pa_xyz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + 

                                1.5 * pa_xxyz[j] * pb_x[j] * fx[j] + 3.0 * pa_xyz[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxx[j] + 

                                pa_xxyz[j] * pb_xxx[j]) * s_0_0[j] + (-3.0 * pa_xyz[j] * fx[j] * fz[j] * fgb[j] + 15.0 * pa_xyz[j] * fx[j] * fx[j] * fz[j] + 

                                22.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_x[j] - 1.5 * fx[j] * pa_yz[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_yz[j] * pb_x[j] * fx[j] - 3.0 * pa_xxyz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                18.0 * pa_xxyz[j] * fz[j] * pb_x[j] * fx[j] + 36.0 * pa_xyz[j] * fx[j] * fz[j] * pb_xx[j] - fz[j] * fga[j] * pa_yz[j] * pb_xxx[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_xxx[j] + 14.0 * pa_xxyz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xxyz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * pa_xxz[j] * fx[j] * fx[j] + 

                                pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_y[j] + 0.5 * pa_xxz[j] * fx[j] * pb_xx[j] + 

                                2.0 * pa_xyz[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxy[j] + pa_xxyz[j] * pb_xxy[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 

                                0.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 0.25 * fz[j] * fga[j] * fx[j] * fx[j] * pa_z[j] - 

                                0.5 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] + 2.5 * pa_xxz[j] * fx[j] * fx[j] * fz[j] + 

                                10.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_y[j] - 

                                0.5 * fx[j] * pa_yz[j] * fz[j] * fgb[j] * pb_y[j] - 0.5 * fz[j] * fga[j] * pa_yz[j] * fx[j] * pb_y[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pa_z[j] * pb_xx[j] - pa_xxyz[j] * fz[j] * fgb[j] * pb_y[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] + 6.0 * pa_xxyz[j] * fz[j] * fx[j] * pb_y[j] + 

                                6.0 * pa_xxz[j] * fx[j] * fz[j] * pb_xx[j] + 24.0 * pa_xyz[j] * fx[j] * fz[j] * pb_xy[j] - fz[j] * fga[j] * pa_yz[j] * pb_xxy[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_xxy[j] + 14.0 * pa_xxyz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xxyz_xxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * pa_xxy[j] * fx[j] * fx[j] + 

                                pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xxy[j] * fx[j] * pb_xx[j] + 

                                2.0 * pa_xyz[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xxz[j] + pa_xxyz[j] * pb_xxz[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] - 

                                0.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 0.25 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 

                                0.5 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] + 2.5 * pa_xxy[j] * fz[j] * fx[j] * fx[j] + 

                                10.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_z[j] - 

                                0.5 * fx[j] * pa_yz[j] * fz[j] * fgb[j] * pb_z[j] - 0.5 * fz[j] * fga[j] * pa_yz[j] * fx[j] * pb_z[j] - 

                                0.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_xx[j] - pa_xxyz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xx[j] + 6.0 * pa_xxyz[j] * fz[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_xxy[j] * fz[j] * fx[j] * pb_xx[j] + 24.0 * pa_xyz[j] * fx[j] * fz[j] * pb_xz[j] - fz[j] * fga[j] * pa_yz[j] * pb_xxz[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_xxz[j] + 14.0 * pa_xxyz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xxyz_xyy[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xxyz[j] * pb_x[j] * fx[j] + 

                                pa_xxz[j] * fx[j] * pb_xy[j] + pa_xyz[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_yz[j] * pb_xyy[j] + 

                                pa_xxyz[j] * pb_xyy[j]) * s_0_0[j] + (-pa_xyz[j] * fx[j] * fz[j] * fgb[j] + 5.0 * pa_xyz[j] * fx[j] * fx[j] * fz[j] + 

                                10.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 0.5 * fx[j] * pa_yz[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_yz[j] * pb_x[j] * fx[j] - fz[j] * fga[j] * fx[j] * pa_z[j] * pb_xy[j] - pa_xxyz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_x[j] + 5.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] + 

                                6.0 * pa_xxyz[j] * fz[j] * pb_x[j] * fx[j] + 12.0 * pa_xxz[j] * fx[j] * fz[j] * pb_xy[j] + 

                                12.0 * pa_xyz[j] * fx[j] * fz[j] * pb_yy[j] - fz[j] * fga[j] * pa_yz[j] * pb_xyy[j] + 6.0 * fx[j] * pa_yz[j] * fz[j] * pb_xyy[j] + 

                                14.0 * pa_xxyz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xxyz_xyz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_xxy[j] * fx[j] * pb_xy[j] + 

                                0.5 * pa_xxz[j] * fx[j] * pb_xz[j] + pa_xyz[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xyz[j] + 

                                pa_xxyz[j] * pb_xyz[j]) * s_0_0[j] + (2.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                0.25 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 5.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                5.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_xy[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pa_z[j] * pb_xz[j] + 2.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] + 6.0 * pa_xxy[j] * fz[j] * fx[j] * pb_xy[j] + 

                                6.0 * pa_xxz[j] * fx[j] * fz[j] * pb_xz[j] + 12.0 * pa_xyz[j] * fx[j] * fz[j] * pb_yz[j] - fz[j] * fga[j] * pa_yz[j] * pb_xyz[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_xyz[j] + 14.0 * pa_xxyz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_y, pa_yz, pa_z, \
                                     pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_xxyz_xzz, t_xxyz_yyy, t_xxyz_yyz, t_xxyz_yzz, t_xxyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xzz[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xxyz[j] * pb_x[j] * fx[j] + 

                                pa_xxy[j] * fx[j] * pb_xz[j] + pa_xyz[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_yz[j] * pb_xzz[j] + 

                                pa_xxyz[j] * pb_xzz[j]) * s_0_0[j] + (-pa_xyz[j] * fx[j] * fz[j] * fgb[j] + 5.0 * pa_xyz[j] * fx[j] * fx[j] * fz[j] + 

                                10.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.5 * fx[j] * pa_yz[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_yz[j] * pb_x[j] * fx[j] - fz[j] * fga[j] * pa_y[j] * fx[j] * pb_xz[j] - pa_xxyz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_x[j] + 5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xz[j] + 

                                6.0 * pa_xxyz[j] * fz[j] * pb_x[j] * fx[j] + 12.0 * pa_xxy[j] * fz[j] * fx[j] * pb_xz[j] + 

                                12.0 * pa_xyz[j] * fx[j] * fz[j] * pb_zz[j] - fz[j] * fga[j] * pa_yz[j] * pb_xzz[j] + 6.0 * fx[j] * pa_yz[j] * fz[j] * pb_xzz[j] + 

                                14.0 * pa_xxyz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xxyz_yyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] + 

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 1.5 * pa_xxyz[j] * pb_y[j] * fx[j] + 

                                1.5 * pa_xxz[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyy[j] + pa_xxyz[j] * pb_yyy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pa_z[j] - 1.5 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 7.5 * pa_xxz[j] * fx[j] * fx[j] * fz[j] - 

                                1.5 * fx[j] * pa_yz[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_yz[j] * pb_y[j] * fx[j] - 

                                1.5 * fz[j] * fga[j] * fx[j] * pa_z[j] * pb_yy[j] - 3.0 * pa_xxyz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] + 

                                18.0 * pa_xxyz[j] * fz[j] * pb_y[j] * fx[j] + 18.0 * pa_xxz[j] * fx[j] * fz[j] * pb_yy[j] - fz[j] * fga[j] * pa_yz[j] * pb_yyy[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_yyy[j] + 14.0 * pa_xxyz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xxyz_yyz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_xxy[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xxyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xxy[j] * fx[j] * pb_yy[j] + 

                                pa_xxz[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyz[j] + pa_xxyz[j] * pb_yyz[j]) * s_0_0[j] + (-0.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.25 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 0.5 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - 

                                0.5 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] + fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 

                                2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 2.5 * pa_xxy[j] * fz[j] * fx[j] * fx[j] + 

                                5.0 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 0.5 * fx[j] * pa_yz[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * fz[j] * fga[j] * pa_yz[j] * fx[j] * pb_z[j] - 0.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_yy[j] - 

                                fz[j] * fga[j] * fx[j] * pa_z[j] * pb_yz[j] - pa_xxyz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_z[j] + 2.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yy[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] + 6.0 * pa_xxyz[j] * fz[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_xxy[j] * fz[j] * fx[j] * pb_yy[j] + 12.0 * pa_xxz[j] * fx[j] * fz[j] * pb_yz[j] - fz[j] * fga[j] * pa_yz[j] * pb_yyz[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_yyz[j] + 14.0 * pa_xxyz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xxyz_yzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_xxz[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xxyz[j] * pb_y[j] * fx[j] + pa_xxy[j] * fx[j] * pb_yz[j] + 

                                0.5 * pa_xxz[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yzz[j] + pa_xxyz[j] * pb_yzz[j]) * s_0_0[j] + (-0.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.25 * fz[j] * fga[j] * fx[j] * fx[j] * pa_z[j] - 0.5 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - 

                                0.5 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] + fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 

                                2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 2.5 * pa_xxz[j] * fx[j] * fx[j] * fz[j] + 

                                5.0 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.5 * fx[j] * pa_yz[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_yz[j] * pb_y[j] * fx[j] - fz[j] * fga[j] * pa_y[j] * fx[j] * pb_yz[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pa_z[j] * pb_zz[j] - pa_xxyz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_y[j] + 5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] + 6.0 * pa_xxyz[j] * fz[j] * pb_y[j] * fx[j] + 

                                12.0 * pa_xxy[j] * fz[j] * fx[j] * pb_yz[j] + 6.0 * pa_xxz[j] * fx[j] * fz[j] * pb_zz[j] - fz[j] * fga[j] * pa_yz[j] * pb_yzz[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_yzz[j] + 14.0 * pa_xxyz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xxyz_zzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] + 

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 1.5 * pa_xxyz[j] * pb_z[j] * fx[j] + 

                                1.5 * pa_xxy[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * pa_yz[j] * pb_zzz[j] + pa_xxyz[j] * pb_zzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 1.5 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 7.5 * pa_xxy[j] * fz[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * pa_yz[j] * pb_z[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_yz[j] * pb_z[j] * fx[j] - 

                                1.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_zz[j] - 3.0 * pa_xxyz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_z[j] + 7.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_zz[j] + 

                                18.0 * pa_xxyz[j] * fz[j] * pb_z[j] * fx[j] + 18.0 * pa_xxy[j] * fz[j] * fx[j] * pb_zz[j] - fz[j] * fga[j] * pa_yz[j] * pb_zzz[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_zzz[j] + 14.0 * pa_xxyz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, \
                                     pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     r_0_0, s_0_0, t_xxzz_xxx, t_xxzz_xxy, t_xxzz_xxz, t_xxzz_xyy, t_xxzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xxx[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                1.5 * pa_xzz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * pa_xxzz[j] * pb_x[j] * fx[j] + 

                                3.0 * pa_xzz[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxx[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_xxx[j] + pa_xxzz[j] * pb_xxx[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 3.0 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                3.0 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 6.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 15.0 * pa_xzz[j] * fx[j] * fx[j] * fz[j] + 

                                22.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] - 0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xx[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                3.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] - 1.5 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_zz[j] * pb_x[j] * fx[j] - 3.0 * pa_xxzz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                18.0 * pa_xxzz[j] * fz[j] * pb_x[j] * fx[j] + 36.0 * pa_xzz[j] * fx[j] * fz[j] * pb_xx[j] - fx[j] * fz[j] * fga[j] * pb_xxx[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_xxx[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] - fz[j] * fga[j] * pa_zz[j] * pb_xxx[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxx[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxx[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xxzz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_xxzz[j] * fx[j] * pb_y[j] + 

                                2.0 * pa_xzz[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxy[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_xxy[j] + pa_xxzz[j] * pb_xxy[j]) * s_0_0[j] + (-fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 0.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.5 * pa_xx[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 2.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] - 

                                0.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_y[j] - 0.5 * fz[j] * fga[j] * pa_zz[j] * fx[j] * pb_y[j] - 

                                pa_xxzz[j] * fz[j] * fgb[j] * pb_y[j] + 2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                10.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 6.0 * pa_xxzz[j] * fz[j] * fx[j] * pb_y[j] + 

                                24.0 * pa_xzz[j] * fx[j] * fz[j] * pb_xy[j] - fx[j] * fz[j] * fga[j] * pb_xxy[j] - pa_xx[j] * fz[j] * fga[j] * pb_xxy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] - fz[j] * fga[j] * pa_zz[j] * pb_xxy[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxy[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxy[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xxzz_xxz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_xxz[j] * fx[j] * fx[j] + 2.0 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_xxzz[j] * fx[j] * pb_z[j] + pa_xxz[j] * fx[j] * pb_xx[j] + 

                                2.0 * pa_xzz[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxz[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_xxz[j] + pa_xxzz[j] * pb_xxz[j]) * s_0_0[j] + (6.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] - 

                                0.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - pa_xxz[j] * fx[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                5.0 * pa_xxz[j] * fz[j] * fx[j] * fx[j] + 20.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] - 0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.5 * pa_xx[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                2.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] - 0.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * fz[j] * fga[j] * pa_zz[j] * fx[j] * pb_z[j] - fz[j] * fga[j] * pa_z[j] * fx[j] * pb_xx[j] - pa_xxzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 10.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                5.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xx[j] + 6.0 * pa_xxzz[j] * fz[j] * fx[j] * pb_z[j] + 

                                12.0 * pa_xxz[j] * fz[j] * fx[j] * pb_xx[j] + 24.0 * pa_xzz[j] * fx[j] * fz[j] * pb_xz[j] - fx[j] * fz[j] * fga[j] * pb_xxz[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_xxz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] - fz[j] * fga[j] * pa_zz[j] * pb_xxz[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxz[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxz[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xxzz_xyy[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] * fx[j] + 

                                0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 0.5 * pa_xxzz[j] * pb_x[j] * fx[j] + 

                                pa_xzz[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyy[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_xyy[j] + pa_xxzz[j] * pb_xyy[j]) * s_0_0[j] + (-0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 2.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                5.0 * pa_xzz[j] * fx[j] * fx[j] * fz[j] - 0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 0.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * pa_xx[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] - 

                                0.5 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.5 * fz[j] * fga[j] * pa_zz[j] * pb_x[j] * fx[j] - pa_xxzz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] + 6.0 * pa_xxzz[j] * fz[j] * pb_x[j] * fx[j] + 

                                12.0 * pa_xzz[j] * fx[j] * fz[j] * pb_yy[j] - fx[j] * fz[j] * fga[j] * pb_xyy[j] - pa_xx[j] * fz[j] * fga[j] * pb_xyy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] - fz[j] * fga[j] * pa_zz[j] * pb_xyy[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyy[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xyy[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xxzz_xyz[j] = (pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + pa_xxz[j] * fx[j] * pb_xy[j] + pa_xzz[j] * fx[j] * pb_yz[j] + 

                                0.25 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xyz[j] + 

                                pa_xxzz[j] * pb_xyz[j]) * s_0_0[j] + (10.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] - fz[j] * fga[j] * pa_z[j] * fx[j] * pb_xy[j] + 

                                5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 5.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xy[j] + 

                                12.0 * pa_xxz[j] * fz[j] * fx[j] * pb_xy[j] + 12.0 * pa_xzz[j] * fx[j] * fz[j] * pb_yz[j] - fx[j] * fz[j] * fga[j] * pb_xyz[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_xyz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xyz[j] - fz[j] * fga[j] * pa_zz[j] * pb_xyz[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyz[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xyz[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, \
                                     pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, \
                                     s_0_0, t_xxzz_xzz, t_xxzz_yyy, t_xxzz_yyz, t_xxzz_yzz, t_xxzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xzz[j] * fx[j] * fx[j] + 2.0 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 

                                0.5 * pa_xxzz[j] * pb_x[j] * fx[j] + 2.0 * pa_xxz[j] * fx[j] * pb_xz[j] + pa_xzz[j] * fx[j] * pb_zz[j] + 

                                0.25 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xzz[j] + 

                                pa_xxzz[j] * pb_xzz[j]) * s_0_0[j] + (6.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 5.0 * pa_xzz[j] * fx[j] * fx[j] * fz[j] + 

                                20.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.5 * pa_xx[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                pa_x[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] - 0.5 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_zz[j] * pb_x[j] * fx[j] - 2.0 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_xz[j] - 

                                pa_xxzz[j] * pb_x[j] * fz[j] * fgb[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] + 10.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xz[j] + 

                                6.0 * pa_xxzz[j] * fz[j] * pb_x[j] * fx[j] + 24.0 * pa_xxz[j] * fz[j] * fx[j] * pb_xz[j] + 

                                12.0 * pa_xzz[j] * fx[j] * fz[j] * pb_zz[j] - fx[j] * fz[j] * fga[j] * pb_xzz[j] - pa_xx[j] * fz[j] * fga[j] * pb_xzz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] - fz[j] * fga[j] * pa_zz[j] * pb_xzz[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xzz[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xzz[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xxzz_yyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 1.5 * pa_xxzz[j] * pb_y[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 

                                0.5 * pa_xx[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyy[j] + pa_xxzz[j] * pb_yyy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 1.5 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xx[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 1.5 * fx[j] * pa_zz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 1.5 * fz[j] * fga[j] * pa_zz[j] * pb_y[j] * fx[j] - 

                                3.0 * pa_xxzz[j] * pb_y[j] * fz[j] * fgb[j] + 7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] + 18.0 * pa_xxzz[j] * fz[j] * pb_y[j] * fx[j] - fx[j] * fz[j] * fga[j] * pb_yyy[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_yyy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] - fz[j] * fga[j] * pa_zz[j] * pb_yyy[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyy[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_yyy[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xxzz_yyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_xxz[j] * fx[j] * fx[j] + 

                                0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_xxzz[j] * fx[j] * pb_z[j] + 

                                pa_xxz[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yyz[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_yyz[j] + pa_xxzz[j] * pb_yyz[j]) * s_0_0[j] + (-0.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - pa_xxz[j] * fx[j] * fz[j] * fgb[j] + 2.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 

                                5.0 * pa_xxz[j] * fz[j] * fx[j] * fx[j] - 0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 0.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * pa_xx[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 0.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.5 * fz[j] * fga[j] * pa_zz[j] * fx[j] * pb_z[j] - 

                                fz[j] * fga[j] * pa_z[j] * fx[j] * pb_yy[j] - pa_xxzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 

                                5.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yy[j] + 6.0 * pa_xxzz[j] * fz[j] * fx[j] * pb_z[j] + 

                                12.0 * pa_xxz[j] * fz[j] * fx[j] * pb_yy[j] - fx[j] * fz[j] * fga[j] * pb_yyz[j] - pa_xx[j] * fz[j] * fga[j] * pb_yyz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] - fz[j] * fga[j] * pa_zz[j] * pb_yyz[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyz[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_yyz[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xxzz_yzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xxzz[j] * pb_y[j] * fx[j] + 

                                2.0 * pa_xxz[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_yzz[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_yzz[j] + pa_xxzz[j] * pb_yzz[j]) * s_0_0[j] + (-fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 0.5 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.5 * pa_xx[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 0.5 * fx[j] * pa_zz[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_zz[j] * pb_y[j] * fx[j] - 2.0 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_yz[j] - 

                                pa_xxzz[j] * pb_y[j] * fz[j] * fgb[j] + 2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] + 

                                10.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yz[j] + 6.0 * pa_xxzz[j] * fz[j] * pb_y[j] * fx[j] + 

                                24.0 * pa_xxz[j] * fz[j] * fx[j] * pb_yz[j] - fx[j] * fz[j] * fga[j] * pb_yzz[j] - pa_xx[j] * fz[j] * fga[j] * pb_yzz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] - fz[j] * fga[j] * pa_zz[j] * pb_yzz[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_yzz[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_yzz[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xxzz_zzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_xxz[j] * fx[j] * fx[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 1.5 * pa_xxzz[j] * pb_z[j] * fx[j] + 

                                3.0 * pa_xxz[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 0.5 * pa_xx[j] * fx[j] * pb_zzz[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_zzz[j] + pa_xxzz[j] * pb_zzz[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 3.0 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - 

                                3.0 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] + 6.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 15.0 * pa_xxz[j] * fz[j] * fx[j] * fx[j] + 

                                22.5 * pa_xx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.75 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 1.5 * pa_xx[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 

                                1.5 * fx[j] * pa_zz[j] * pb_z[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_zz[j] * pb_z[j] * fx[j] - 

                                3.0 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_zz[j] - 3.0 * pa_xxzz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_zz[j] + 

                                18.0 * pa_xxzz[j] * fz[j] * pb_z[j] * fx[j] + 36.0 * pa_xxz[j] * fz[j] * fx[j] * pb_zz[j] - fx[j] * fz[j] * fga[j] * pb_zzz[j] - 

                                pa_xx[j] * fz[j] * fga[j] * pb_zzz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] - fz[j] * fga[j] * pa_zz[j] * pb_zzz[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_zzz[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_zzz[j] + 14.0 * pa_xxzz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (12) = (60,65)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_xyyy_xxx, t_xyyy_xxy, t_xyyy_xxz, t_xyyy_xyy, t_xyyy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xxx[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * pa_yyy[j] + 

                                2.25 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 1.5 * pa_xyyy[j] * pb_x[j] * fx[j] + 

                                1.5 * fx[j] * pa_yyy[j] * pb_xx[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxx[j] + pa_xyyy[j] * pb_xxx[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                2.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fga[j] - 1.5 * fx[j] * pa_yyy[j] * fz[j] * fgb[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 7.5 * fx[j] * fx[j] * pa_yyy[j] * fz[j] - 

                                4.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 4.5 * pa_xy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                4.5 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_xx[j] - 3.0 * pa_xyyy[j] * pb_x[j] * fz[j] * fgb[j] + 

                                22.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 22.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xx[j] + 

                                18.0 * pa_xyyy[j] * fz[j] * pb_x[j] * fx[j] + 18.0 * fx[j] * pa_yyy[j] * fz[j] * pb_xx[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xxx[j] + 18.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxx[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xyyy_xxy[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xyy[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.5 * pa_xyyy[j] * fx[j] * pb_y[j] + 1.5 * pa_xyy[j] * fx[j] * pb_xx[j] + 

                                fx[j] * pa_yyy[j] * pb_xy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xxy[j] + pa_xyyy[j] * pb_xxy[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                1.5 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] + 

                                15.0 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_x[j] - 1.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                1.5 * pa_xy[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 1.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] - 

                                3.0 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_xy[j] - pa_xyyy[j] * fz[j] * fgb[j] * pb_y[j] + 

                                7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                15.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xy[j] + 6.0 * pa_xyyy[j] * fz[j] * fx[j] * pb_y[j] + 

                                18.0 * pa_xyy[j] * fz[j] * fx[j] * pb_xx[j] + 12.0 * fx[j] * pa_yyy[j] * fz[j] * pb_xy[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xxy[j] + 18.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxy[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xyyy_xxz[j] = (0.75 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xyyy[j] * fx[j] * pb_z[j] + fx[j] * pa_yyy[j] * pb_xz[j] + 

                                1.5 * pa_xy[j] * fx[j] * pb_xxz[j] + pa_xyyy[j] * pb_xxz[j]) * s_0_0[j] + (-1.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                1.5 * pa_xy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 3.0 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_xz[j] - 

                                pa_xyyy[j] * fz[j] * fgb[j] * pb_z[j] + 7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                15.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xz[j] + 6.0 * pa_xyyy[j] * fz[j] * fx[j] * pb_z[j] + 

                                12.0 * fx[j] * pa_yyy[j] * fz[j] * pb_xz[j] - 3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xxz[j] + 

                                18.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxz[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xyyy_xyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                2.25 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_yyy[j] + 1.5 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.5 * pa_xyyy[j] * pb_x[j] * fx[j] + 

                                3.0 * pa_xyy[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_yyy[j] * pb_yy[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyy[j] + 

                                pa_xyyy[j] * pb_xyy[j]) * s_0_0[j] + (9.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] - 

                                0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fga[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 0.5 * fx[j] * pa_yyy[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 22.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * fx[j] * fx[j] * pa_yyy[j] * fz[j] + 15.0 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_y[j] - 

                                1.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                3.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] - 1.5 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_yy[j] - 

                                pa_xyyy[j] * pb_x[j] * fz[j] * fgb[j] + 15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                7.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yy[j] + 6.0 * pa_xyyy[j] * fz[j] * pb_x[j] * fx[j] + 

                                36.0 * pa_xyy[j] * fz[j] * fx[j] * pb_xy[j] + 6.0 * fx[j] * pa_yyy[j] * fz[j] * pb_yy[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xyy[j] + 18.0 * pa_xy[j] * fz[j] * fx[j] * pb_xyy[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xyyy_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 1.5 * pa_xyy[j] * fx[j] * pb_xz[j] + 

                                0.5 * fx[j] * pa_yyy[j] * pb_yz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xyz[j] + pa_xyyy[j] * pb_xyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] - 1.5 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_yz[j] + 

                                7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yz[j] + 

                                18.0 * pa_xyy[j] * fz[j] * fx[j] * pb_xz[j] + 6.0 * fx[j] * pa_yyy[j] * fz[j] * pb_yz[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xyz[j] + 18.0 * pa_xy[j] * fz[j] * fx[j] * pb_xyz[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (13) = (65,70)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yyy, pb_x, pb_xzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyyy_xzz, \
                                     t_xyyy_yyy, t_xyyy_yyz, t_xyyy_yzz, t_xyyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * pa_yyy[j] + 

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.5 * pa_xyyy[j] * pb_x[j] * fx[j] + 

                                0.5 * fx[j] * pa_yyy[j] * pb_zz[j] + 1.5 * pa_xy[j] * fx[j] * pb_xzz[j] + pa_xyyy[j] * pb_xzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fga[j] - 0.5 * fx[j] * pa_yyy[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 2.5 * fx[j] * fx[j] * pa_yyy[j] * fz[j] - 

                                1.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                1.5 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_zz[j] - pa_xyyy[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_zz[j] + 

                                6.0 * pa_xyyy[j] * fz[j] * pb_x[j] * fx[j] + 6.0 * fx[j] * pa_yyy[j] * fz[j] * pb_zz[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_xzz[j] + 18.0 * pa_xy[j] * fz[j] * fx[j] * pb_xzz[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xyyy_yyy[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_xyy[j] * fx[j] * fx[j] + 

                                6.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xyyy[j] * pb_y[j] * fx[j] + 

                                4.5 * pa_xyy[j] * fx[j] * pb_yy[j] + 1.5 * pa_xy[j] * fx[j] * pb_yyy[j] + pa_xyyy[j] * pb_yyy[j]) * s_0_0[j] + (15.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                2.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 2.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                4.5 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 22.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] + 

                                67.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 4.5 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                4.5 * pa_xy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 4.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] - 

                                3.0 * pa_xyyy[j] * pb_y[j] * fz[j] * fgb[j] + 22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                18.0 * pa_xyyy[j] * fz[j] * pb_y[j] * fx[j] + 54.0 * pa_xyy[j] * fz[j] * fx[j] * pb_yy[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_yyy[j] + 18.0 * pa_xy[j] * fz[j] * fx[j] * pb_yyy[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xyyy_yyz[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyyy[j] * fx[j] * pb_z[j] + 3.0 * pa_xyy[j] * fx[j] * pb_yz[j] + 

                                1.5 * pa_xy[j] * fx[j] * pb_yyz[j] + pa_xyyy[j] * pb_yyz[j]) * s_0_0[j] + (22.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                3.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] - pa_xyyy[j] * fz[j] * fgb[j] * pb_z[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 6.0 * pa_xyyy[j] * fz[j] * fx[j] * pb_z[j] + 

                                36.0 * pa_xyy[j] * fz[j] * fx[j] * pb_yz[j] - 3.0 * pa_xy[j] * fz[j] * fga[j] * pb_yyz[j] + 

                                18.0 * pa_xy[j] * fz[j] * fx[j] * pb_yyz[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xyyy_yzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] + 

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xyyy[j] * pb_y[j] * fx[j] + 

                                1.5 * pa_xyy[j] * fx[j] * pb_zz[j] + 1.5 * pa_xy[j] * fx[j] * pb_yzz[j] + pa_xyyy[j] * pb_yzz[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 1.5 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 7.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] - 

                                1.5 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                1.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] - pa_xyyy[j] * pb_y[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                6.0 * pa_xyyy[j] * fz[j] * pb_y[j] * fx[j] + 18.0 * pa_xyy[j] * fz[j] * fx[j] * pb_zz[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_yzz[j] + 18.0 * pa_xy[j] * fz[j] * fx[j] * pb_yzz[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xyyy_zzz[j] = (2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_xyyy[j] * pb_z[j] * fx[j] + 

                                1.5 * pa_xy[j] * fx[j] * pb_zzz[j] + pa_xyyy[j] * pb_zzz[j]) * s_0_0[j] + (-4.5 * pa_xy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                4.5 * pa_xy[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 3.0 * pa_xyyy[j] * pb_z[j] * fz[j] * fgb[j] + 

                                22.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 18.0 * pa_xyyy[j] * fz[j] * pb_z[j] * fx[j] - 

                                3.0 * pa_xy[j] * fz[j] * fga[j] * pb_zzz[j] + 18.0 * pa_xy[j] * fz[j] * fx[j] * pb_zzz[j] + 14.0 * pa_xyyy[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (14) = (70,75)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_y, pa_yy, pa_yyz, \
                                     pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xyyz_xxx, t_xyyz_xxy, t_xyyz_xxz, t_xyyz_xyy, \
                                     t_xyyz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xxx[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * pa_yyz[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 1.5 * pa_xyyz[j] * pb_x[j] * fx[j] + 

                                1.5 * fx[j] * pa_yyz[j] * pb_xx[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxx[j] + pa_xyyz[j] * pb_xxx[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pa_z[j] - 1.5 * fx[j] * pa_yyz[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 7.5 * fx[j] * fx[j] * pa_yyz[j] * fz[j] - 

                                1.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_xx[j] - 3.0 * pa_xyyz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] + 

                                18.0 * pa_xyyz[j] * fz[j] * pb_x[j] * fx[j] + 18.0 * fx[j] * pa_yyz[j] * fz[j] * pb_xx[j] - pa_xz[j] * fz[j] * fga[j] * pb_xxx[j] + 

                                6.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxx[j] + 14.0 * pa_xyyz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xyyz_xxy[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + fx[j] * fx[j] * pa_yz[j] * pb_x[j] + 

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_y[j] + 

                                pa_xyz[j] * fx[j] * pb_xx[j] + fx[j] * pa_yyz[j] * pb_xy[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxy[j] + 

                                pa_xyyz[j] * pb_xxy[j]) * s_0_0[j] + (-pa_xyz[j] * fx[j] * fz[j] * fgb[j] + 5.0 * pa_xyz[j] * fx[j] * fx[j] * fz[j] + 

                                10.0 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_x[j] - 0.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.5 * pa_xz[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - fx[j] * fz[j] * fga[j] * pa_z[j] * pb_xy[j] - pa_xyyz[j] * fz[j] * fgb[j] * pb_y[j] + 

                                2.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 5.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] + 

                                6.0 * pa_xyyz[j] * fz[j] * fx[j] * pb_y[j] + 12.0 * pa_xyz[j] * fx[j] * fz[j] * pb_xx[j] + 

                                12.0 * fx[j] * pa_yyz[j] * fz[j] * pb_xy[j] - pa_xz[j] * fz[j] * fga[j] * pb_xxy[j] + 6.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxy[j] + 

                                14.0 * pa_xyyz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xyyz_xxz[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xyy[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xyy[j] * fx[j] * pb_xx[j] + 

                                fx[j] * pa_yyz[j] * pb_xz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxz[j] + pa_xyyz[j] * pb_xxz[j]) * s_0_0[j] + (-0.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                0.5 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] + pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 2.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] + 

                                5.0 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_x[j] - 0.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * pa_xz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_xx[j] - 

                                fx[j] * fz[j] * fga[j] * pa_z[j] * pb_xz[j] - pa_xyyz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 2.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] + 6.0 * pa_xyyz[j] * fz[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_xyy[j] * fz[j] * fx[j] * pb_xx[j] + 12.0 * fx[j] * pa_yyz[j] * fz[j] * pb_xz[j] - pa_xz[j] * fz[j] * fga[j] * pb_xxz[j] + 

                                6.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxz[j] + 14.0 * pa_xyyz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xyyz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * fx[j] * fx[j] * pa_yyz[j] + fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 

                                0.5 * pa_xyyz[j] * pb_x[j] * fx[j] + 2.0 * pa_xyz[j] * fx[j] * pb_xy[j] + 0.5 * fx[j] * pa_yyz[j] * pb_yy[j] + 

                                0.5 * pa_xz[j] * fx[j] * pb_xyy[j] + pa_xyyz[j] * pb_xyy[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 

                                0.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pa_z[j] - 

                                0.5 * fx[j] * pa_yyz[j] * fz[j] * fgb[j] + 7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * fx[j] * fx[j] * pa_yyz[j] * fz[j] + 10.0 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_y[j] - 

                                0.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.5 * pa_xz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_yy[j] - pa_xyyz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] + 6.0 * pa_xyyz[j] * fz[j] * pb_x[j] * fx[j] + 

                                24.0 * pa_xyz[j] * fx[j] * fz[j] * pb_xy[j] + 6.0 * fx[j] * pa_yyz[j] * fz[j] * pb_yy[j] - pa_xz[j] * fz[j] * fga[j] * pb_xyy[j] + 

                                6.0 * pa_xz[j] * fx[j] * fz[j] * pb_xyy[j] + 14.0 * pa_xyyz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xyyz_xyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.5 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 

                                0.5 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xyy[j] * fx[j] * pb_xy[j] + pa_xyz[j] * fx[j] * pb_xz[j] + 

                                0.5 * fx[j] * pa_yyz[j] * pb_yz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xyz[j] + pa_xyyz[j] * pb_xyz[j]) * s_0_0[j] + (2.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                5.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 2.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_y[j] + 

                                5.0 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_z[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_xy[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_yz[j] + 2.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] + 6.0 * pa_xyy[j] * fz[j] * fx[j] * pb_xy[j] + 

                                12.0 * pa_xyz[j] * fx[j] * fz[j] * pb_xz[j] + 6.0 * fx[j] * pa_yyz[j] * fz[j] * pb_yz[j] - pa_xz[j] * fz[j] * fga[j] * pb_xyz[j] + 

                                6.0 * pa_xz[j] * fx[j] * fz[j] * pb_xyz[j] + 14.0 * pa_xyyz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (15) = (75,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_yy, pa_yyz, pa_z, \
                                     pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_xyyz_xzz, t_xyyz_yyy, t_xyyz_yyz, t_xyyz_yzz, t_xyyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_yyz[j] + 0.5 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + 

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xyyz[j] * pb_x[j] * fx[j] + pa_xyy[j] * fx[j] * pb_xz[j] + 

                                0.5 * fx[j] * pa_yyz[j] * pb_zz[j] + 0.5 * pa_xz[j] * fx[j] * pb_xzz[j] + pa_xyyz[j] * pb_xzz[j]) * s_0_0[j] + (-0.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pa_z[j] - 0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                0.5 * fx[j] * pa_yyz[j] * fz[j] * fgb[j] + fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 

                                2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 2.5 * fx[j] * fx[j] * pa_yyz[j] * fz[j] + 

                                5.0 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_z[j] - 0.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * pa_xz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - pa_x[j] * fz[j] * fga[j] * fx[j] * pb_xz[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_zz[j] - pa_xyyz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] + 6.0 * pa_xyyz[j] * fz[j] * pb_x[j] * fx[j] + 

                                12.0 * pa_xyy[j] * fz[j] * fx[j] * pb_xz[j] + 6.0 * fx[j] * pa_yyz[j] * fz[j] * pb_zz[j] - pa_xz[j] * fz[j] * fga[j] * pb_xzz[j] + 

                                6.0 * pa_xz[j] * fx[j] * fz[j] * pb_xzz[j] + 14.0 * pa_xyyz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xyyz_yyy[j] = (1.5 * pa_xyz[j] * fx[j] * fx[j] + 2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_xyyz[j] * pb_y[j] * fx[j] + 3.0 * pa_xyz[j] * fx[j] * pb_yy[j] + 0.5 * pa_xz[j] * fx[j] * pb_yyy[j] + 

                                pa_xyyz[j] * pb_yyy[j]) * s_0_0[j] + (-3.0 * pa_xyz[j] * fx[j] * fz[j] * fgb[j] + 15.0 * pa_xyz[j] * fx[j] * fx[j] * fz[j] + 

                                22.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 1.5 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 3.0 * pa_xyyz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                18.0 * pa_xyyz[j] * fz[j] * pb_y[j] * fx[j] + 36.0 * pa_xyz[j] * fx[j] * fz[j] * pb_yy[j] - pa_xz[j] * fz[j] * fga[j] * pb_yyy[j] + 

                                6.0 * pa_xz[j] * fx[j] * fz[j] * pb_yyy[j] + 14.0 * pa_xyyz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xyyz_yyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xyy[j] * fx[j] * fx[j] + 

                                pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xyyz[j] * fx[j] * pb_z[j] + 0.5 * pa_xyy[j] * fx[j] * pb_yy[j] + 

                                2.0 * pa_xyz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xz[j] * fx[j] * pb_yyz[j] + pa_xyyz[j] * pb_yyz[j]) * s_0_0[j] + (3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                0.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 0.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                0.5 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 2.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] + 

                                10.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.5 * pa_xz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - pa_xyyz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 6.0 * pa_xyyz[j] * fz[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_xyy[j] * fz[j] * fx[j] * pb_yy[j] + 24.0 * pa_xyz[j] * fx[j] * fz[j] * pb_yz[j] - pa_xz[j] * fz[j] * fga[j] * pb_yyz[j] + 

                                6.0 * pa_xz[j] * fx[j] * fz[j] * pb_yyz[j] + 14.0 * pa_xyyz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xyyz_yzz[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyyz[j] * pb_y[j] * fx[j] + 

                                pa_xyy[j] * fx[j] * pb_yz[j] + pa_xyz[j] * fx[j] * pb_zz[j] + 0.5 * pa_xz[j] * fx[j] * pb_yzz[j] + 

                                pa_xyyz[j] * pb_yzz[j]) * s_0_0[j] + (-pa_xyz[j] * fx[j] * fz[j] * fgb[j] + 5.0 * pa_xyz[j] * fx[j] * fx[j] * fz[j] + 

                                10.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.5 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.5 * pa_xz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - pa_xyyz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 

                                6.0 * pa_xyyz[j] * fz[j] * pb_y[j] * fx[j] + 12.0 * pa_xyy[j] * fz[j] * fx[j] * pb_yz[j] + 

                                12.0 * pa_xyz[j] * fx[j] * fz[j] * pb_zz[j] - pa_xz[j] * fz[j] * fga[j] * pb_yzz[j] + 6.0 * pa_xz[j] * fx[j] * fz[j] * pb_yzz[j] + 

                                14.0 * pa_xyyz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xyyz_zzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xyyz[j] * pb_z[j] * fx[j] + 

                                1.5 * pa_xyy[j] * fx[j] * pb_zz[j] + 0.5 * pa_xz[j] * fx[j] * pb_zzz[j] + pa_xyyz[j] * pb_zzz[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 1.5 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 7.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] - 

                                1.5 * pa_xz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - 3.0 * pa_xyyz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                18.0 * pa_xyyz[j] * fz[j] * pb_z[j] * fx[j] + 18.0 * pa_xyy[j] * fz[j] * fx[j] * pb_zz[j] - pa_xz[j] * fz[j] * fga[j] * pb_zzz[j] + 

                                6.0 * pa_xz[j] * fx[j] * fz[j] * pb_zzz[j] + 14.0 * pa_xyyz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (16) = (80,85)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, \
                                     pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xyzz_xxx, t_xyzz_xxy, t_xyzz_xxz, t_xyzz_xyy, \
                                     t_xyzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xxx[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * pa_yzz[j] + 

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 1.5 * pa_xyzz[j] * pb_x[j] * fx[j] + 

                                1.5 * fx[j] * pa_yzz[j] * pb_xx[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxx[j] + pa_xyzz[j] * pb_xxx[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fga[j] - 1.5 * fx[j] * pa_yzz[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 7.5 * fx[j] * fx[j] * pa_yzz[j] * fz[j] - 

                                1.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                1.5 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_xx[j] - 3.0 * pa_xyzz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xx[j] + 

                                18.0 * pa_xyzz[j] * fz[j] * pb_x[j] * fx[j] + 18.0 * fx[j] * pa_yzz[j] * fz[j] * pb_xx[j] - pa_xy[j] * fz[j] * fga[j] * pb_xxx[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxx[j] + 14.0 * pa_xyzz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xyzz_xxy[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xzz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                0.5 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_y[j] + 0.5 * pa_xzz[j] * fx[j] * pb_xx[j] + 

                                fx[j] * pa_yzz[j] * pb_xy[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxy[j] + pa_xyzz[j] * pb_xxy[j]) * s_0_0[j] + (-0.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                0.5 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] + pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 2.5 * pa_xzz[j] * fx[j] * fx[j] * fz[j] + 

                                5.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] - 0.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.5 * pa_xy[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 0.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] - 

                                fx[j] * pa_y[j] * fz[j] * fga[j] * pb_xy[j] - pa_xyzz[j] * fz[j] * fgb[j] * pb_y[j] + 

                                2.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 2.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xy[j] + 6.0 * pa_xyzz[j] * fz[j] * fx[j] * pb_y[j] + 

                                6.0 * pa_xzz[j] * fx[j] * fz[j] * pb_xx[j] + 12.0 * fx[j] * pa_yzz[j] * fz[j] * pb_xy[j] - pa_xy[j] * fz[j] * fga[j] * pb_xxy[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxy[j] + 14.0 * pa_xyzz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xyzz_xxz[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + fx[j] * fx[j] * pa_yz[j] * pb_x[j] + 

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_z[j] + 

                                pa_xyz[j] * fx[j] * pb_xx[j] + fx[j] * pa_yzz[j] * pb_xz[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxz[j] + 

                                pa_xyzz[j] * pb_xxz[j]) * s_0_0[j] + (-pa_xyz[j] * fx[j] * fz[j] * fgb[j] + 5.0 * pa_xyz[j] * fz[j] * fx[j] * fx[j] + 

                                10.0 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_x[j] - 0.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * pa_xy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - fx[j] * pa_y[j] * fz[j] * fga[j] * pb_xz[j] - pa_xyzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xz[j] + 

                                6.0 * pa_xyzz[j] * fz[j] * fx[j] * pb_z[j] + 12.0 * pa_xyz[j] * fz[j] * fx[j] * pb_xx[j] + 

                                12.0 * fx[j] * pa_yzz[j] * fz[j] * pb_xz[j] - pa_xy[j] * fz[j] * fga[j] * pb_xxz[j] + 6.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxz[j] + 

                                14.0 * pa_xyzz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xyzz_xyy[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * fx[j] * fx[j] * pa_yzz[j] + 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.5 * pa_xyzz[j] * pb_x[j] * fx[j] + pa_xzz[j] * fx[j] * pb_xy[j] + 

                                0.5 * fx[j] * pa_yzz[j] * pb_yy[j] + 0.5 * pa_xy[j] * fx[j] * pb_xyy[j] + pa_xyzz[j] * pb_xyy[j]) * s_0_0[j] + (-0.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fga[j] - 0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                0.5 * fx[j] * pa_yzz[j] * fz[j] * fgb[j] + fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 

                                2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 2.5 * fx[j] * fx[j] * pa_yzz[j] * fz[j] + 

                                5.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] - 0.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * pa_xy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] - 

                                0.5 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_yy[j] - pa_xyzz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                2.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yy[j] + 6.0 * pa_xyzz[j] * fz[j] * pb_x[j] * fx[j] + 

                                12.0 * pa_xzz[j] * fx[j] * fz[j] * pb_xy[j] + 6.0 * fx[j] * pa_yzz[j] * fz[j] * pb_yy[j] - pa_xy[j] * fz[j] * fga[j] * pb_xyy[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_xyy[j] + 14.0 * pa_xyzz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xyzz_xyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + pa_xyz[j] * fx[j] * pb_xy[j] + 0.5 * pa_xzz[j] * fx[j] * pb_xz[j] + 

                                0.5 * fx[j] * pa_yzz[j] * pb_yz[j] + 0.5 * pa_xy[j] * fx[j] * pb_xyz[j] + pa_xyzz[j] * pb_xyz[j]) * s_0_0[j] + (2.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                5.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 5.0 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_y[j] + 

                                2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] - 0.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] - 

                                0.5 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_yz[j] + 2.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                2.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yz[j] + 12.0 * pa_xyz[j] * fz[j] * fx[j] * pb_xy[j] + 

                                6.0 * pa_xzz[j] * fx[j] * fz[j] * pb_xz[j] + 6.0 * fx[j] * pa_yzz[j] * fz[j] * pb_yz[j] - pa_xy[j] * fz[j] * fga[j] * pb_xyz[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_xyz[j] + 14.0 * pa_xyzz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (17) = (85,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, \
                                     pb_x, pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     r_0_0, s_0_0, t_xyzz_xzz, t_xyzz_yyy, t_xyzz_yyz, t_xyzz_yzz, t_xyzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * fx[j] * fx[j] * pa_yzz[j] + fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 

                                0.5 * pa_xyzz[j] * pb_x[j] * fx[j] + 2.0 * pa_xyz[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_yzz[j] * pb_zz[j] + 

                                0.5 * pa_xy[j] * fx[j] * pb_xzz[j] + pa_xyzz[j] * pb_xzz[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] - 

                                0.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 0.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fga[j] - 

                                0.5 * fx[j] * pa_yzz[j] * fz[j] * fgb[j] + 7.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * fx[j] * fx[j] * pa_yzz[j] * fz[j] + 10.0 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_z[j] - 

                                0.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.5 * pa_xy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                0.5 * fx[j] * pa_y[j] * fz[j] * fga[j] * pb_zz[j] - pa_xyzz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_zz[j] + 6.0 * pa_xyzz[j] * fz[j] * pb_x[j] * fx[j] + 

                                24.0 * pa_xyz[j] * fz[j] * fx[j] * pb_xz[j] + 6.0 * fx[j] * pa_yzz[j] * fz[j] * pb_zz[j] - pa_xy[j] * fz[j] * fga[j] * pb_xzz[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_xzz[j] + 14.0 * pa_xyzz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xyzz_yyy[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] + 

                                0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xyzz[j] * pb_y[j] * fx[j] + 

                                1.5 * pa_xzz[j] * fx[j] * pb_yy[j] + 0.5 * pa_xy[j] * fx[j] * pb_yyy[j] + pa_xyzz[j] * pb_yyy[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 1.5 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 7.5 * pa_xzz[j] * fx[j] * fx[j] * fz[j] - 

                                1.5 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * pa_xy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                1.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] - 3.0 * pa_xyzz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                18.0 * pa_xyzz[j] * fz[j] * pb_y[j] * fx[j] + 18.0 * pa_xzz[j] * fx[j] * fz[j] * pb_yy[j] - pa_xy[j] * fz[j] * fga[j] * pb_yyy[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_yyy[j] + 14.0 * pa_xyzz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xyzz_yyz[j] = (0.5 * pa_xyz[j] * fx[j] * fx[j] + pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyzz[j] * fx[j] * pb_z[j] + 

                                pa_xyz[j] * fx[j] * pb_yy[j] + pa_xzz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xy[j] * fx[j] * pb_yyz[j] + 

                                pa_xyzz[j] * pb_yyz[j]) * s_0_0[j] + (-pa_xyz[j] * fx[j] * fz[j] * fgb[j] + 5.0 * pa_xyz[j] * fz[j] * fx[j] * fx[j] + 

                                10.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 0.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * pa_xy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] - pa_xyzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 

                                6.0 * pa_xyzz[j] * fz[j] * fx[j] * pb_z[j] + 12.0 * pa_xyz[j] * fz[j] * fx[j] * pb_yy[j] + 

                                12.0 * pa_xzz[j] * fx[j] * fz[j] * pb_yz[j] - pa_xy[j] * fz[j] * fga[j] * pb_yyz[j] + 6.0 * pa_xy[j] * fz[j] * fx[j] * pb_yyz[j] + 

                                14.0 * pa_xyzz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xyzz_yzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_xzz[j] * fx[j] * fx[j] + pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 

                                0.5 * pa_xyzz[j] * pb_y[j] * fx[j] + 2.0 * pa_xyz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xzz[j] * fx[j] * pb_zz[j] + 

                                0.5 * pa_xy[j] * fx[j] * pb_yzz[j] + pa_xyzz[j] * pb_yzz[j]) * s_0_0[j] + (3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                0.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 0.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                0.5 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 7.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                2.5 * pa_xzz[j] * fx[j] * fx[j] * fz[j] + 10.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.5 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 0.5 * pa_xy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                0.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] - pa_xyzz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                2.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 6.0 * pa_xyzz[j] * fz[j] * pb_y[j] * fx[j] + 

                                24.0 * pa_xyz[j] * fz[j] * fx[j] * pb_yz[j] + 6.0 * pa_xzz[j] * fx[j] * fz[j] * pb_zz[j] - pa_xy[j] * fz[j] * fga[j] * pb_yzz[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_yzz[j] + 14.0 * pa_xyzz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xyzz_zzz[j] = (1.5 * pa_xyz[j] * fx[j] * fx[j] + 2.25 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_xyzz[j] * pb_z[j] * fx[j] + 3.0 * pa_xyz[j] * fx[j] * pb_zz[j] + 0.5 * pa_xy[j] * fx[j] * pb_zzz[j] + 

                                pa_xyzz[j] * pb_zzz[j]) * s_0_0[j] + (-3.0 * pa_xyz[j] * fx[j] * fz[j] * fgb[j] + 15.0 * pa_xyz[j] * fz[j] * fx[j] * fx[j] + 

                                22.5 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 1.5 * pa_xy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xy[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 3.0 * pa_xyzz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                18.0 * pa_xyzz[j] * fz[j] * pb_z[j] * fx[j] + 36.0 * pa_xyz[j] * fz[j] * fx[j] * pb_zz[j] - pa_xy[j] * fz[j] * fga[j] * pb_zzz[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_zzz[j] + 14.0 * pa_xyzz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (18) = (90,95)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_xzzz_xxx, t_xzzz_xxy, t_xzzz_xxz, t_xzzz_xyy, t_xzzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xxx[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * pa_zzz[j] + 

                                2.25 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 1.5 * pa_xzzz[j] * pb_x[j] * fx[j] + 

                                1.5 * fx[j] * pa_zzz[j] * pb_xx[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxx[j] + pa_xzzz[j] * pb_xxx[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                2.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fga[j] - 1.5 * fx[j] * pa_zzz[j] * fz[j] * fgb[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 7.5 * fx[j] * fx[j] * pa_zzz[j] * fz[j] - 

                                4.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 4.5 * pa_xz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                4.5 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_xx[j] - 3.0 * pa_xzzz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                22.5 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 22.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xx[j] + 

                                18.0 * pa_xzzz[j] * fz[j] * pb_x[j] * fx[j] + 18.0 * fx[j] * pa_zzz[j] * fz[j] * pb_xx[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xxx[j] + 18.0 * pa_xz[j] * fz[j] * fx[j] * pb_xxx[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_xzzz_xxy[j] = (0.75 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xzzz[j] * fx[j] * pb_y[j] + fx[j] * pa_zzz[j] * pb_xy[j] + 

                                1.5 * pa_xz[j] * fx[j] * pb_xxy[j] + pa_xzzz[j] * pb_xxy[j]) * s_0_0[j] + (-1.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                1.5 * pa_xz[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 3.0 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_xy[j] - 

                                pa_xzzz[j] * fz[j] * fgb[j] * pb_y[j] + 7.5 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xy[j] + 6.0 * pa_xzzz[j] * fz[j] * fx[j] * pb_y[j] + 

                                12.0 * fx[j] * pa_zzz[j] * fz[j] * pb_xy[j] - 3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xxy[j] + 

                                18.0 * pa_xz[j] * fz[j] * fx[j] * pb_xxy[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_xzzz_xxz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xzz[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_xzzz[j] * fx[j] * pb_z[j] + 1.5 * pa_xzz[j] * fx[j] * pb_xx[j] + 

                                fx[j] * pa_zzz[j] * pb_xz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xxz[j] + pa_xzzz[j] * pb_xxz[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                1.5 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_xzz[j] * fz[j] * fx[j] * fx[j] + 

                                15.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] - 1.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                1.5 * pa_xz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 1.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] - 

                                3.0 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_xz[j] - pa_xzzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                7.5 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xz[j] + 6.0 * pa_xzzz[j] * fz[j] * fx[j] * pb_z[j] + 

                                18.0 * pa_xzz[j] * fz[j] * fx[j] * pb_xx[j] + 12.0 * fx[j] * pa_zzz[j] * fz[j] * pb_xz[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xxz[j] + 18.0 * pa_xz[j] * fz[j] * fx[j] * pb_xxz[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_xzzz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * pa_zzz[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_xzzz[j] * pb_x[j] * fx[j] + 

                                0.5 * fx[j] * pa_zzz[j] * pb_yy[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyy[j] + pa_xzzz[j] * pb_xyy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fga[j] - 0.5 * fx[j] * pa_zzz[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 2.5 * fx[j] * fx[j] * pa_zzz[j] * fz[j] - 

                                1.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                1.5 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_yy[j] - pa_xzzz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yy[j] + 

                                6.0 * pa_xzzz[j] * fz[j] * pb_x[j] * fx[j] + 6.0 * fx[j] * pa_zzz[j] * fz[j] * pb_yy[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xyy[j] + 18.0 * pa_xz[j] * fz[j] * fx[j] * pb_xyy[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_xzzz_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 1.5 * pa_xzz[j] * fx[j] * pb_xy[j] + 

                                0.5 * fx[j] * pa_zzz[j] * pb_yz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xyz[j] + pa_xzzz[j] * pb_xyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] - 

                                1.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] - 1.5 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_yz[j] + 

                                7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yz[j] + 

                                18.0 * pa_xzz[j] * fz[j] * fx[j] * pb_xy[j] + 6.0 * fx[j] * pa_zzz[j] * fz[j] * pb_yz[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xyz[j] + 18.0 * pa_xz[j] * fz[j] * fx[j] * pb_xyz[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (19) = (95,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_xzzz_xzz, t_xzzz_yyy, t_xzzz_yyz, t_xzzz_yzz, t_xzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                2.25 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_zzz[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xzzz[j] * pb_x[j] * fx[j] + 

                                3.0 * pa_xzz[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * pa_zzz[j] * pb_zz[j] + 1.5 * pa_xz[j] * fx[j] * pb_xzz[j] + 

                                pa_xzzz[j] * pb_xzz[j]) * s_0_0[j] + (9.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] - 

                                0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fga[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 0.5 * fx[j] * pa_zzz[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 22.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * fx[j] * fx[j] * pa_zzz[j] * fz[j] + 15.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                3.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] - 1.5 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_zz[j] - 

                                pa_xzzz[j] * pb_x[j] * fz[j] * fgb[j] + 15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                7.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_zz[j] + 6.0 * pa_xzzz[j] * fz[j] * pb_x[j] * fx[j] + 

                                36.0 * pa_xzz[j] * fz[j] * fx[j] * pb_xz[j] + 6.0 * fx[j] * pa_zzz[j] * fz[j] * pb_zz[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_xzz[j] + 18.0 * pa_xz[j] * fz[j] * fx[j] * pb_xzz[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_xzzz_yyy[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_xzzz[j] * pb_y[j] * fx[j] + 

                                1.5 * pa_xz[j] * fx[j] * pb_yyy[j] + pa_xzzz[j] * pb_yyy[j]) * s_0_0[j] + (-4.5 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                4.5 * pa_xz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 3.0 * pa_xzzz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                22.5 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 18.0 * pa_xzzz[j] * fz[j] * pb_y[j] * fx[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_yyy[j] + 18.0 * pa_xz[j] * fz[j] * fx[j] * pb_yyy[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_xzzz_yyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xzzz[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_xzz[j] * fx[j] * pb_yy[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyz[j] + pa_xzzz[j] * pb_yyz[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 1.5 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 7.5 * pa_xzz[j] * fz[j] * fx[j] * fx[j] - 

                                1.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                1.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] - pa_xzzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                7.5 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                6.0 * pa_xzzz[j] * fz[j] * fx[j] * pb_z[j] + 18.0 * pa_xzz[j] * fz[j] * fx[j] * pb_yy[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_yyz[j] + 18.0 * pa_xz[j] * fz[j] * fx[j] * pb_yyz[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_xzzz_yzz[j] = (2.25 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xzzz[j] * pb_y[j] * fx[j] + 3.0 * pa_xzz[j] * fx[j] * pb_yz[j] + 

                                1.5 * pa_xz[j] * fx[j] * pb_yzz[j] + pa_xzzz[j] * pb_yzz[j]) * s_0_0[j] + (22.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                1.5 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * pa_xz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                3.0 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] - pa_xzzz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 6.0 * pa_xzzz[j] * fz[j] * pb_y[j] * fx[j] + 

                                36.0 * pa_xzz[j] * fz[j] * fx[j] * pb_yz[j] - 3.0 * pa_xz[j] * fz[j] * fga[j] * pb_yzz[j] + 

                                18.0 * pa_xz[j] * fz[j] * fx[j] * pb_yzz[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_xzzz_zzz[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_xzz[j] * fx[j] * fx[j] + 

                                6.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xzzz[j] * pb_z[j] * fx[j] + 

                                4.5 * pa_xzz[j] * fx[j] * pb_zz[j] + 1.5 * pa_xz[j] * fx[j] * pb_zzz[j] + pa_xzzz[j] * pb_zzz[j]) * s_0_0[j] + (15.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                2.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 2.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                4.5 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 22.5 * pa_xzz[j] * fz[j] * fx[j] * fx[j] + 

                                67.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 4.5 * pa_xz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                4.5 * pa_xz[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 4.5 * pa_x[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] - 

                                3.0 * pa_xzzz[j] * pb_z[j] * fz[j] * fgb[j] + 22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                18.0 * pa_xzzz[j] * fz[j] * pb_z[j] * fx[j] + 54.0 * pa_xzz[j] * fz[j] * fx[j] * pb_zz[j] - 

                                3.0 * pa_xz[j] * fz[j] * fga[j] * pb_zzz[j] + 18.0 * pa_xz[j] * fz[j] * fx[j] * pb_zzz[j] + 14.0 * pa_xzzz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (20) = (100,105)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, \
                                     pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, r_0_0, s_0_0, t_yyyy_xxx, t_yyyy_xxy, \
                                     t_yyyy_xxz, t_yyyy_xyy, t_yyyy_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xxx[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                1.5 * pa_yyyy[j] * pb_x[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxx[j] + 

                                pa_yyyy[j] * pb_xxx[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 9.0 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                9.0 * pa_yy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                3.0 * pa_yyyy[j] * pb_x[j] * fz[j] * fgb[j] + 45.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                18.0 * pa_yyyy[j] * fz[j] * pb_x[j] * fx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xxx[j] - 

                                6.0 * pa_yy[j] * fz[j] * fga[j] * pb_xxx[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 

                                36.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxx[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_yyyy_xxy[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * fx[j] + pa_yyy[j] * fx[j] * fx[j] + 

                                0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 

                                0.5 * pa_yyyy[j] * fx[j] * pb_y[j] + 2.0 * pa_yyy[j] * fx[j] * pb_xx[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] + 

                                3.0 * pa_yy[j] * fx[j] * pb_xxy[j] + pa_yyyy[j] * pb_xxy[j]) * s_0_0[j] + (-3.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * fga[j] - 2.0 * pa_yyy[j] * fx[j] * fz[j] * fgb[j] + 

                                12.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 10.0 * pa_yyy[j] * fz[j] * fx[j] * fx[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                3.0 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 3.0 * pa_yy[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 

                                6.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                pa_yyyy[j] * fz[j] * fgb[j] * pb_y[j] + 15.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                30.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 6.0 * pa_yyyy[j] * fz[j] * fx[j] * pb_y[j] + 

                                24.0 * pa_yyy[j] * fz[j] * fx[j] * pb_xx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xxy[j] - 

                                6.0 * pa_yy[j] * fz[j] * fga[j] * pb_xxy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 

                                36.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxy[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_yyyy_xxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_yyyy[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxz[j] + 

                                pa_yyyy[j] * pb_xxz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 3.0 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                3.0 * pa_yy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                pa_yyyy[j] * fz[j] * fgb[j] * pb_z[j] + 15.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_yyyy[j] * fz[j] * fx[j] * pb_z[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xxz[j] - 

                                6.0 * pa_yy[j] * fz[j] * fga[j] * pb_xxz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 

                                36.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxz[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_yyyy_xyy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                6.0 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_yyyy[j] * pb_x[j] * fx[j] + 4.0 * pa_yyy[j] * fx[j] * pb_xy[j] + 

                                0.75 * fx[j] * fx[j] * pb_xyy[j] + 3.0 * pa_yy[j] * fx[j] * pb_xyy[j] + pa_yyyy[j] * pb_xyy[j]) * s_0_0[j] + (-4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] + 

                                15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 45.0 * pa_yy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 3.0 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                3.0 * pa_yy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 12.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] - 

                                pa_yyyy[j] * pb_x[j] * fz[j] * fgb[j] + 60.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                6.0 * pa_yyyy[j] * fz[j] * pb_x[j] * fx[j] + 48.0 * pa_yyy[j] * fz[j] * fx[j] * pb_xy[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xyy[j] - 6.0 * pa_yy[j] * fz[j] * fga[j] * pb_xyy[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 36.0 * pa_yy[j] * fz[j] * fx[j] * pb_xyy[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_yyyy_xyz[j] = (3.0 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 2.0 * pa_yyy[j] * fx[j] * pb_xz[j] + 

                                0.75 * fx[j] * fx[j] * pb_xyz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xyz[j] + pa_yyyy[j] * pb_xyz[j]) * s_0_0[j] + (-6.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] + 

                                30.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 24.0 * pa_yyy[j] * fz[j] * fx[j] * pb_xz[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xyz[j] - 6.0 * pa_yy[j] * fz[j] * fga[j] * pb_xyz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 36.0 * pa_yy[j] * fz[j] * fx[j] * pb_xyz[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (21) = (105,110)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_yyyy_xzz, t_yyyy_yyy, \
                                     t_yyyy_yyz, t_yyyy_yzz, t_yyyy_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.5 * pa_yyyy[j] * pb_x[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] + 3.0 * pa_yy[j] * fx[j] * pb_xzz[j] + 

                                pa_yyyy[j] * pb_xzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 3.0 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                3.0 * pa_yy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                pa_yyyy[j] * pb_x[j] * fz[j] * fgb[j] + 15.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                6.0 * pa_yyyy[j] * fz[j] * pb_x[j] * fx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xzz[j] - 

                                6.0 * pa_yy[j] * fz[j] * fga[j] * pb_xzz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 

                                36.0 * pa_yy[j] * fz[j] * fx[j] * pb_xzz[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_yyyy_yyy[j] = (7.5 * pa_y[j] * fx[j] * fx[j] * fx[j] + 5.625 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                3.0 * pa_yyy[j] * fx[j] * fx[j] + 13.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 

                                9.0 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_yyyy[j] * pb_y[j] * fx[j] + 6.0 * pa_yyy[j] * fx[j] * pb_yy[j] + 

                                0.75 * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yy[j] * fx[j] * pb_yyy[j] + pa_yyyy[j] * pb_yyy[j]) * s_0_0[j] + (60.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                9.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 9.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                13.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 6.0 * pa_yyy[j] * fx[j] * fz[j] * fgb[j] + 

                                45.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 30.0 * pa_yyy[j] * fz[j] * fx[j] * fx[j] + 

                                135.0 * pa_yy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 2.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                9.0 * pa_yy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 9.0 * pa_yy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                18.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] - 3.0 * pa_yyyy[j] * pb_y[j] * fz[j] * fgb[j] + 

                                90.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 18.0 * pa_yyyy[j] * fz[j] * pb_y[j] * fx[j] + 

                                72.0 * pa_yyy[j] * fz[j] * fx[j] * pb_yy[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_yyy[j] - 

                                6.0 * pa_yy[j] * fz[j] * fga[j] * pb_yyy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 

                                36.0 * pa_yy[j] * fz[j] * fx[j] * pb_yyy[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_yyyy_yyz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_yyyy[j] * fx[j] * pb_z[j] + 4.0 * pa_yyy[j] * fx[j] * pb_yz[j] + 

                                0.75 * fx[j] * fx[j] * pb_yyz[j] + 3.0 * pa_yy[j] * fx[j] * pb_yyz[j] + pa_yyyy[j] * pb_yyz[j]) * s_0_0[j] + (-4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] + 

                                15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 45.0 * pa_yy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 3.0 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                3.0 * pa_yy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 12.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] - 

                                pa_yyyy[j] * fz[j] * fgb[j] * pb_z[j] + 60.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 

                                6.0 * pa_yyyy[j] * fz[j] * fx[j] * pb_z[j] + 48.0 * pa_yyy[j] * fz[j] * fx[j] * pb_yz[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_yyz[j] - 6.0 * pa_yy[j] * fz[j] * fga[j] * pb_yyz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 36.0 * pa_yy[j] * fz[j] * fx[j] * pb_yyz[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_yyyy_yzz[j] = (1.5 * pa_y[j] * fx[j] * fx[j] * fx[j] + pa_yyy[j] * fx[j] * fx[j] + 

                                0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 

                                0.5 * pa_yyyy[j] * pb_y[j] * fx[j] + 2.0 * pa_yyy[j] * fx[j] * pb_zz[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] + 

                                3.0 * pa_yy[j] * fx[j] * pb_yzz[j] + pa_yyyy[j] * pb_yzz[j]) * s_0_0[j] + (-3.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * fga[j] - 2.0 * pa_yyy[j] * fx[j] * fz[j] * fgb[j] + 

                                12.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 10.0 * pa_yyy[j] * fz[j] * fx[j] * fx[j] - 

                                0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                3.0 * pa_yy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 3.0 * pa_yy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                6.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                pa_yyyy[j] * pb_y[j] * fz[j] * fgb[j] + 15.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                30.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 6.0 * pa_yyyy[j] * fz[j] * pb_y[j] * fx[j] + 

                                24.0 * pa_yyy[j] * fz[j] * fx[j] * pb_zz[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_yzz[j] - 

                                6.0 * pa_yy[j] * fz[j] * fga[j] * pb_yzz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 

                                36.0 * pa_yy[j] * fz[j] * fx[j] * pb_yzz[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_yyyy_zzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_yyyy[j] * pb_z[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_yy[j] * fx[j] * pb_zzz[j] + 

                                pa_yyyy[j] * pb_zzz[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 9.0 * pa_yy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                9.0 * pa_yy[j] * fz[j] * fga[j] * pb_z[j] * fx[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                3.0 * pa_yyyy[j] * pb_z[j] * fz[j] * fgb[j] + 45.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                18.0 * pa_yyyy[j] * fz[j] * pb_z[j] * fx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_zzz[j] - 

                                6.0 * pa_yy[j] * fz[j] * fga[j] * pb_zzz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 

                                36.0 * pa_yy[j] * fz[j] * fx[j] * pb_zzz[j] + 14.0 * pa_yyyy[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (22) = (110,115)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, r_0_0, s_0_0, \
                                     t_yyyz_xxx, t_yyyz_xxy, t_yyyz_xxz, t_yyyz_xyy, t_yyyz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xxx[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yyyz[j] * pb_x[j] * fx[j] + 

                                1.5 * pa_yz[j] * fx[j] * pb_xxx[j] + pa_yyyz[j] * pb_xxx[j]) * s_0_0[j] + (-4.5 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                4.5 * pa_yz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 3.0 * pa_yyyz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                22.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 18.0 * pa_yyyz[j] * fz[j] * pb_x[j] * fx[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xxx[j] + 18.0 * pa_yz[j] * fx[j] * fz[j] * pb_xxx[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_yyyz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_yyz[j] * fx[j] * fx[j] + 

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_yyyz[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_yyz[j] * fx[j] * pb_xx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxy[j] + pa_yyyz[j] * pb_xxy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pa_z[j] - 1.5 * pa_yyz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 7.5 * pa_yyz[j] * fx[j] * fx[j] * fz[j] - 

                                1.5 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 1.5 * pa_yz[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_xx[j] - pa_yyyz[j] * fz[j] * fgb[j] * pb_y[j] + 

                                7.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] + 

                                6.0 * pa_yyyz[j] * fz[j] * fx[j] * pb_y[j] + 18.0 * pa_yyz[j] * fx[j] * fz[j] * pb_xx[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xxy[j] + 18.0 * pa_yz[j] * fx[j] * fz[j] * pb_xxy[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_yyyz_xxz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_yyy[j] * fx[j] * fx[j] + 

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_yyyz[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_yyy[j] * fx[j] * pb_xx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxz[j] + pa_yyyz[j] * pb_xxz[j]) * s_0_0[j] + (-0.75 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * pa_yyy[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 2.5 * pa_yyy[j] * fz[j] * fx[j] * fx[j] - 

                                1.5 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * pa_yz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                1.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_xx[j] - pa_yyyz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                7.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                6.0 * pa_yyyz[j] * fz[j] * fx[j] * pb_z[j] + 6.0 * pa_yyy[j] * fz[j] * fx[j] * pb_xx[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xxz[j] + 18.0 * pa_yz[j] * fx[j] * fz[j] * pb_xxz[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_yyyz_xyy[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_yyyz[j] * pb_x[j] * fx[j] + 3.0 * pa_yyz[j] * fx[j] * pb_xy[j] + 

                                1.5 * pa_yz[j] * fx[j] * pb_xyy[j] + pa_yyyz[j] * pb_xyy[j]) * s_0_0[j] + (22.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                1.5 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_yz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_xy[j] - pa_yyyz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] + 6.0 * pa_yyyz[j] * fz[j] * pb_x[j] * fx[j] + 

                                36.0 * pa_yyz[j] * fx[j] * fz[j] * pb_xy[j] - 3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xyy[j] + 

                                18.0 * pa_yz[j] * fx[j] * fz[j] * pb_xyy[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_yyyz_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_yyy[j] * fx[j] * pb_xy[j] + 

                                1.5 * pa_yyz[j] * fx[j] * pb_xz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xyz[j] + pa_yyyz[j] * pb_xyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_yy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                1.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_xy[j] - 1.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_xz[j] + 

                                7.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] + 

                                6.0 * pa_yyy[j] * fz[j] * fx[j] * pb_xy[j] + 18.0 * pa_yyz[j] * fx[j] * fz[j] * pb_xz[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xyz[j] + 18.0 * pa_yz[j] * fx[j] * fz[j] * pb_xyz[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (23) = (115,120)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_yyyz_xzz, t_yyyz_yyy, t_yyyz_yyz, t_yyyz_yzz, t_yyyz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xzz[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 

                                1.5 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_yyyz[j] * pb_x[j] * fx[j] + pa_yyy[j] * fx[j] * pb_xz[j] + 

                                1.5 * pa_yz[j] * fx[j] * pb_xzz[j] + pa_yyyz[j] * pb_xzz[j]) * s_0_0[j] + (-1.5 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * pa_yz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 3.0 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_xz[j] - 

                                pa_yyyz[j] * pb_x[j] * fz[j] * fgb[j] + 7.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                15.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 6.0 * pa_yyyz[j] * fz[j] * pb_x[j] * fx[j] + 

                                12.0 * pa_yyy[j] * fz[j] * fx[j] * pb_xz[j] - 3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xzz[j] + 

                                18.0 * pa_yz[j] * fx[j] * fz[j] * pb_xzz[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_yyyz_yyy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] + 2.25 * pa_yyz[j] * fx[j] * fx[j] + 

                                6.75 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 1.5 * pa_yyyz[j] * pb_y[j] * fx[j] + 

                                4.5 * pa_yyz[j] * fx[j] * pb_yy[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyy[j] + pa_yyyz[j] * pb_yyy[j]) * s_0_0[j] + (15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 

                                2.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pa_z[j] - 

                                4.5 * pa_yyz[j] * fx[j] * fz[j] * fgb[j] + 22.5 * pa_yyz[j] * fx[j] * fx[j] * fz[j] + 

                                67.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 4.5 * pa_yz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                4.5 * pa_yz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 4.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_yy[j] - 

                                3.0 * pa_yyyz[j] * pb_y[j] * fz[j] * fgb[j] + 22.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] + 

                                18.0 * pa_yyyz[j] * fz[j] * pb_y[j] * fx[j] + 54.0 * pa_yyz[j] * fx[j] * fz[j] * pb_yy[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_yyy[j] + 18.0 * pa_yz[j] * fx[j] * fz[j] * pb_yyy[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_yyyz_yyz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_yyy[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 

                                2.25 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_yyyz[j] * fx[j] * pb_z[j] + 0.5 * pa_yyy[j] * fx[j] * pb_yy[j] + 

                                3.0 * pa_yyz[j] * fx[j] * pb_yz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyz[j] + pa_yyyz[j] * pb_yyz[j]) * s_0_0[j] + (9.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                0.75 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 0.75 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 0.5 * pa_yyy[j] * fx[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 2.5 * pa_yyy[j] * fz[j] * fx[j] * fx[j] + 

                                15.0 * pa_yy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 22.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * pa_yz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                1.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - 3.0 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_yz[j] - 

                                pa_yyyz[j] * fz[j] * fgb[j] * pb_z[j] + 7.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] + 6.0 * pa_yyyz[j] * fz[j] * fx[j] * pb_z[j] + 

                                6.0 * pa_yyy[j] * fz[j] * fx[j] * pb_yy[j] + 36.0 * pa_yyz[j] * fx[j] * fz[j] * pb_yz[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_yyz[j] + 18.0 * pa_yz[j] * fx[j] * fz[j] * pb_yyz[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_yyyz_yzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_yyz[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_yyyz[j] * pb_y[j] * fx[j] + pa_yyy[j] * fx[j] * pb_yz[j] + 

                                1.5 * pa_yyz[j] * fx[j] * pb_zz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yzz[j] + pa_yyyz[j] * pb_yzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pa_z[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                1.5 * pa_yyz[j] * fx[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_yyz[j] * fx[j] * fx[j] * fz[j] + 

                                15.0 * pa_yy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 1.5 * pa_yz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * pa_yz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 3.0 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pa_z[j] * pb_zz[j] - pa_yyyz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                7.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 15.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] + 6.0 * pa_yyyz[j] * fz[j] * pb_y[j] * fx[j] + 

                                12.0 * pa_yyy[j] * fz[j] * fx[j] * pb_yz[j] + 18.0 * pa_yyz[j] * fx[j] * fz[j] * pb_zz[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_yzz[j] + 18.0 * pa_yz[j] * fx[j] * fz[j] * pb_yzz[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_yyyz_zzz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yyy[j] * fx[j] * fx[j] + 

                                2.25 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_yyyz[j] * pb_z[j] * fx[j] + 

                                1.5 * pa_yyy[j] * fx[j] * pb_zz[j] + 1.5 * pa_yz[j] * fx[j] * pb_zzz[j] + pa_yyyz[j] * pb_zzz[j]) * s_0_0[j] + (-2.25 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                2.25 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 1.5 * pa_yyy[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 7.5 * pa_yyy[j] * fz[j] * fx[j] * fx[j] - 

                                4.5 * pa_yz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 4.5 * pa_yz[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 

                                4.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - 3.0 * pa_yyyz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                22.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 22.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                18.0 * pa_yyyz[j] * fz[j] * pb_z[j] * fx[j] + 18.0 * pa_yyy[j] * fz[j] * fx[j] * pb_zz[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_zzz[j] + 18.0 * pa_yz[j] * fx[j] * fz[j] * pb_zzz[j] + 14.0 * pa_yyyz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (24) = (120,125)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, \
                                     pb_xx, pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, r_0_0, s_0_0, \
                                     t_yyzz_xxx, t_yyzz_xxy, t_yyzz_xxz, t_yyzz_xyy, t_yyzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xxx[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 1.5 * pa_yyzz[j] * pb_x[j] * fx[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 

                                0.5 * pa_yy[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxx[j] + pa_yyzz[j] * pb_xxx[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 1.5 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * pa_yy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 1.5 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 1.5 * fz[j] * fga[j] * pa_zz[j] * pb_x[j] * fx[j] - 

                                3.0 * pa_yyzz[j] * pb_x[j] * fz[j] * fgb[j] + 7.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] + 18.0 * pa_yyzz[j] * fz[j] * pb_x[j] * fx[j] - fx[j] * fz[j] * fga[j] * pb_xxx[j] - 

                                pa_yy[j] * fz[j] * fga[j] * pb_xxx[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] - fz[j] * fga[j] * pa_zz[j] * pb_xxx[j] + 

                                6.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxx[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxx[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_yyzz_xxy[j] = (0.25 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.5 * pa_yzz[j] * fx[j] * fx[j] + 

                                0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.5 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 0.5 * pa_yyzz[j] * fx[j] * pb_y[j] + 

                                pa_yzz[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_yy[j] * fx[j] * pb_xxy[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_xxy[j] + pa_yyzz[j] * pb_xxy[j]) * s_0_0[j] + (-0.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fga[j] - pa_yzz[j] * fx[j] * fz[j] * fgb[j] + 2.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                5.0 * pa_yzz[j] * fx[j] * fx[j] * fz[j] - 0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 0.5 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.5 * pa_yy[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - pa_y[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] - 

                                0.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_y[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.5 * fz[j] * fga[j] * pa_zz[j] * fx[j] * pb_y[j] - pa_yyzz[j] * fz[j] * fgb[j] * pb_y[j] + 

                                2.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 5.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] + 6.0 * pa_yyzz[j] * fz[j] * fx[j] * pb_y[j] + 

                                12.0 * pa_yzz[j] * fx[j] * fz[j] * pb_xx[j] - fx[j] * fz[j] * fga[j] * pb_xxy[j] - pa_yy[j] * fz[j] * fga[j] * pb_xxy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] - fz[j] * fga[j] * pa_zz[j] * pb_xxy[j] + 6.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxy[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxy[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_yyzz_xxz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.5 * pa_yyz[j] * fx[j] * fx[j] + 

                                0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_yyzz[j] * fx[j] * pb_z[j] + 

                                pa_yyz[j] * fx[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xxz[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_xxz[j] + pa_yyzz[j] * pb_xxz[j]) * s_0_0[j] + (-0.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - pa_yyz[j] * fx[j] * fz[j] * fgb[j] + 2.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 

                                5.0 * pa_yyz[j] * fz[j] * fx[j] * fx[j] - 0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 0.5 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * pa_yy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 0.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.5 * fz[j] * fga[j] * pa_zz[j] * fx[j] * pb_z[j] - 

                                fz[j] * fga[j] * pa_z[j] * fx[j] * pb_xx[j] - pa_yyzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 

                                5.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xx[j] + 6.0 * pa_yyzz[j] * fz[j] * fx[j] * pb_z[j] + 

                                12.0 * pa_yyz[j] * fz[j] * fx[j] * pb_xx[j] - fx[j] * fz[j] * fga[j] * pb_xxz[j] - pa_yy[j] * fz[j] * fga[j] * pb_xxz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] - fz[j] * fga[j] * pa_zz[j] * pb_xxz[j] + 6.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxz[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxz[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_yyzz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 

                                0.25 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.5 * pa_yyzz[j] * pb_x[j] * fx[j] + 

                                2.0 * pa_yzz[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_yy[j] * fx[j] * pb_xyy[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_xyy[j] + pa_yyzz[j] * pb_xyy[j]) * s_0_0[j] + (-fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] - 

                                0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.5 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * pa_yy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 2.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] - 

                                0.5 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] - 0.5 * fz[j] * fga[j] * pa_zz[j] * pb_x[j] * fx[j] - 

                                pa_yyzz[j] * pb_x[j] * fz[j] * fgb[j] + 2.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                10.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 6.0 * pa_yyzz[j] * fz[j] * pb_x[j] * fx[j] + 

                                24.0 * pa_yzz[j] * fx[j] * fz[j] * pb_xy[j] - fx[j] * fz[j] * fga[j] * pb_xyy[j] - pa_yy[j] * fz[j] * fga[j] * pb_xyy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] - fz[j] * fga[j] * pa_zz[j] * pb_xyy[j] + 6.0 * pa_yy[j] * fz[j] * fx[j] * pb_xyy[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xyy[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_yyzz_xyz[j] = (pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + pa_yyz[j] * fx[j] * pb_xy[j] + pa_yzz[j] * fx[j] * pb_xz[j] + 

                                0.25 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xyz[j] + 

                                pa_yyzz[j] * pb_xyz[j]) * s_0_0[j] + (10.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                pa_y[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] - fz[j] * fga[j] * pa_z[j] * fx[j] * pb_xy[j] + 

                                5.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 5.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xy[j] + 

                                12.0 * pa_yyz[j] * fz[j] * fx[j] * pb_xy[j] + 12.0 * pa_yzz[j] * fx[j] * fz[j] * pb_xz[j] - fx[j] * fz[j] * fga[j] * pb_xyz[j] - 

                                pa_yy[j] * fz[j] * fga[j] * pb_xyz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xyz[j] - fz[j] * fga[j] * pa_zz[j] * pb_xyz[j] + 

                                6.0 * pa_yy[j] * fz[j] * fx[j] * pb_xyz[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xyz[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (25) = (125,130)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, \
                                     pb_xz, pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, \
                                     s_0_0, t_yyzz_xzz, t_yyzz_yyy, t_yyzz_yyz, t_yyzz_yzz, t_yyzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.5 * pa_yyzz[j] * pb_x[j] * fx[j] + 

                                2.0 * pa_yyz[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_xzz[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_xzz[j] + pa_yyzz[j] * pb_xzz[j]) * s_0_0[j] + (-fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_yy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.5 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * pa_yy[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 0.5 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_zz[j] * pb_x[j] * fx[j] - 2.0 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_xz[j] - 

                                pa_yyzz[j] * pb_x[j] * fz[j] * fgb[j] + 2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] + 

                                10.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xz[j] + 6.0 * pa_yyzz[j] * fz[j] * pb_x[j] * fx[j] + 

                                24.0 * pa_yyz[j] * fz[j] * fx[j] * pb_xz[j] - fx[j] * fz[j] * fga[j] * pb_xzz[j] - pa_yy[j] * fz[j] * fga[j] * pb_xzz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] - fz[j] * fga[j] * pa_zz[j] * pb_xzz[j] + 6.0 * pa_yy[j] * fz[j] * fx[j] * pb_xzz[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xzz[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_yyzz_yyy[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_yzz[j] * fx[j] * fx[j] + 2.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 

                                0.75 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_yyzz[j] * pb_y[j] * fx[j] + 

                                3.0 * pa_yzz[j] * fx[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 0.5 * pa_yy[j] * fx[j] * pb_yyy[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_yyy[j] + pa_yyzz[j] * pb_yyy[j]) * s_0_0[j] + (-1.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                1.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fga[j] - 3.0 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                3.0 * pa_yzz[j] * fx[j] * fz[j] * fgb[j] + 6.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 15.0 * pa_yzz[j] * fx[j] * fx[j] * fz[j] + 

                                22.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] - 0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * pa_yy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * pa_yy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                3.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] - 1.5 * fx[j] * pa_zz[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_zz[j] * pb_y[j] * fx[j] - 3.0 * pa_yyzz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                7.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 15.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                18.0 * pa_yyzz[j] * fz[j] * pb_y[j] * fx[j] + 36.0 * pa_yzz[j] * fx[j] * fz[j] * pb_yy[j] - fx[j] * fz[j] * fga[j] * pb_yyy[j] - 

                                pa_yy[j] * fz[j] * fga[j] * pb_yyy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] - fz[j] * fga[j] * pa_zz[j] * pb_yyy[j] + 

                                6.0 * pa_yy[j] * fz[j] * fx[j] * pb_yyy[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_yyy[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_yyzz_yyz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_yyz[j] * fx[j] * fx[j] + 2.0 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_yyzz[j] * fx[j] * pb_z[j] + pa_yyz[j] * fx[j] * pb_yy[j] + 

                                2.0 * pa_yzz[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_yy[j] * fx[j] * pb_yyz[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_yyz[j] + pa_yyzz[j] * pb_yyz[j]) * s_0_0[j] + (6.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] - 

                                0.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - pa_yyz[j] * fx[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                5.0 * pa_yyz[j] * fz[j] * fx[j] * fx[j] + 20.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] - 0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.5 * pa_yy[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                2.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] - 0.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.5 * fz[j] * fga[j] * pa_zz[j] * fx[j] * pb_z[j] - fz[j] * fga[j] * pa_z[j] * fx[j] * pb_yy[j] - pa_yyzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 10.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 

                                5.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yy[j] + 6.0 * pa_yyzz[j] * fz[j] * fx[j] * pb_z[j] + 

                                12.0 * pa_yyz[j] * fz[j] * fx[j] * pb_yy[j] + 24.0 * pa_yzz[j] * fx[j] * fz[j] * pb_yz[j] - fx[j] * fz[j] * fga[j] * pb_yyz[j] - 

                                pa_yy[j] * fz[j] * fga[j] * pb_yyz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] - fz[j] * fga[j] * pa_zz[j] * pb_yyz[j] + 

                                6.0 * pa_yy[j] * fz[j] * fx[j] * pb_yyz[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_yyz[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_yyzz_yzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_yzz[j] * fx[j] * fx[j] + 2.0 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 

                                0.5 * pa_yyzz[j] * pb_y[j] * fx[j] + 2.0 * pa_yyz[j] * fx[j] * pb_yz[j] + pa_yzz[j] * fx[j] * pb_zz[j] + 

                                0.25 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yzz[j] + 

                                pa_yyzz[j] * pb_yzz[j]) * s_0_0[j] + (6.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                0.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 0.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - pa_yzz[j] * fx[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                7.5 * pa_yy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 5.0 * pa_yzz[j] * fx[j] * fx[j] * fz[j] + 

                                20.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.5 * pa_yy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 0.5 * pa_yy[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                pa_y[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] - 0.5 * fx[j] * pa_zz[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_zz[j] * pb_y[j] * fx[j] - 2.0 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_yz[j] - 

                                pa_yyzz[j] * pb_y[j] * fz[j] * fgb[j] + 5.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] + 10.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yz[j] + 

                                6.0 * pa_yyzz[j] * fz[j] * pb_y[j] * fx[j] + 24.0 * pa_yyz[j] * fz[j] * fx[j] * pb_yz[j] + 

                                12.0 * pa_yzz[j] * fx[j] * fz[j] * pb_zz[j] - fx[j] * fz[j] * fga[j] * pb_yzz[j] - pa_yy[j] * fz[j] * fga[j] * pb_yzz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] - fz[j] * fga[j] * pa_zz[j] * pb_yzz[j] + 6.0 * pa_yy[j] * fz[j] * fx[j] * pb_yzz[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_yzz[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_yyzz_zzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_yyz[j] * fx[j] * fx[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 1.5 * pa_yyzz[j] * pb_z[j] * fx[j] + 

                                3.0 * pa_yyz[j] * fx[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 0.5 * pa_yy[j] * fx[j] * pb_zzz[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_zzz[j] + pa_yyzz[j] * pb_zzz[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 3.0 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - 

                                3.0 * pa_yyz[j] * fx[j] * fz[j] * fgb[j] + 6.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 15.0 * pa_yyz[j] * fz[j] * fx[j] * fx[j] + 

                                22.5 * pa_yy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 0.75 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                1.5 * pa_yy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 1.5 * pa_yy[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 

                                1.5 * fx[j] * pa_zz[j] * pb_z[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_zz[j] * pb_z[j] * fx[j] - 

                                3.0 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_zz[j] - 3.0 * pa_yyzz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_zz[j] + 

                                18.0 * pa_yyzz[j] * fz[j] * pb_z[j] * fx[j] + 36.0 * pa_yyz[j] * fz[j] * fx[j] * pb_zz[j] - fx[j] * fz[j] * fga[j] * pb_zzz[j] - 

                                pa_yy[j] * fz[j] * fga[j] * pb_zzz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] - fz[j] * fga[j] * pa_zz[j] * pb_zzz[j] + 

                                6.0 * pa_yy[j] * fz[j] * fx[j] * pb_zzz[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_zzz[j] + 14.0 * pa_yyzz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (26) = (130,135)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, \
                                     pb_xxx, pb_xxy, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_z, r_0_0, s_0_0, \
                                     t_yzzz_xxx, t_yzzz_xxy, t_yzzz_xxz, t_yzzz_xyy, t_yzzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xxx[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_yzzz[j] * pb_x[j] * fx[j] + 

                                1.5 * pa_yz[j] * fx[j] * pb_xxx[j] + pa_yzzz[j] * pb_xxx[j]) * s_0_0[j] + (-4.5 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                4.5 * pa_yz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 3.0 * pa_yzzz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                22.5 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 18.0 * pa_yzzz[j] * fz[j] * pb_x[j] * fx[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xxx[j] + 18.0 * pa_yz[j] * fz[j] * fx[j] * pb_xxx[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_yzzz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * pa_zzz[j] + 

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_yzzz[j] * fx[j] * pb_y[j] + 

                                0.5 * fx[j] * pa_zzz[j] * pb_xx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxy[j] + pa_yzzz[j] * pb_xxy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fga[j] - 0.5 * fx[j] * pa_zzz[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 2.5 * fx[j] * fx[j] * pa_zzz[j] * fz[j] - 

                                1.5 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 1.5 * pa_yz[j] * fz[j] * fga[j] * fx[j] * pb_y[j] - 

                                1.5 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_xx[j] - pa_yzzz[j] * fz[j] * fgb[j] * pb_y[j] + 

                                7.5 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xx[j] + 

                                6.0 * pa_yzzz[j] * fz[j] * fx[j] * pb_y[j] + 6.0 * fx[j] * pa_zzz[j] * fz[j] * pb_xx[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xxy[j] + 18.0 * pa_yz[j] * fz[j] * fx[j] * pb_xxy[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_yzzz_xxz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yzz[j] * fx[j] * fx[j] + 

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * pa_yzzz[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_yzz[j] * fx[j] * pb_xx[j] + 1.5 * pa_yz[j] * fx[j] * pb_xxz[j] + pa_yzzz[j] * pb_xxz[j]) * s_0_0[j] + (-0.75 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_y[j] * fx[j] * fx[j] * fz[j] * fga[j] - 1.5 * pa_yzz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 7.5 * pa_yzz[j] * fz[j] * fx[j] * fx[j] - 

                                1.5 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * pa_yz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                1.5 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] - pa_yzzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                7.5 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 7.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 

                                6.0 * pa_yzzz[j] * fz[j] * fx[j] * pb_z[j] + 18.0 * pa_yzz[j] * fz[j] * fx[j] * pb_xx[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xxz[j] + 18.0 * pa_yz[j] * fz[j] * fx[j] * pb_xxz[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_yzzz_xyy[j] = (0.75 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_yzzz[j] * pb_x[j] * fx[j] + fx[j] * pa_zzz[j] * pb_xy[j] + 

                                1.5 * pa_yz[j] * fx[j] * pb_xyy[j] + pa_yzzz[j] * pb_xyy[j]) * s_0_0[j] + (-1.5 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * pa_yz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 3.0 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_xy[j] - 

                                pa_yzzz[j] * pb_x[j] * fz[j] * fgb[j] + 7.5 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xy[j] + 6.0 * pa_yzzz[j] * fz[j] * pb_x[j] * fx[j] + 

                                12.0 * fx[j] * pa_zzz[j] * fz[j] * pb_xy[j] - 3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xyy[j] + 

                                18.0 * pa_yz[j] * fz[j] * fx[j] * pb_xyy[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_yzzz_xyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 1.5 * pa_yzz[j] * fx[j] * pb_xy[j] + 

                                0.5 * fx[j] * pa_zzz[j] * pb_xz[j] + 1.5 * pa_yz[j] * fx[j] * pb_xyz[j] + pa_yzzz[j] * pb_xyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] - 

                                1.5 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] - 1.5 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_xz[j] + 

                                7.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xz[j] + 

                                18.0 * pa_yzz[j] * fz[j] * fx[j] * pb_xy[j] + 6.0 * fx[j] * pa_zzz[j] * fz[j] * pb_xz[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xyz[j] + 18.0 * pa_yz[j] * fz[j] * fx[j] * pb_xyz[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (27) = (135,140)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xz, \
                                     pb_xzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, \
                                     t_yzzz_xzz, t_yzzz_yyy, t_yzzz_yyz, t_yzzz_yzz, t_yzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xzz[j] = (2.25 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 

                                1.5 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_yzzz[j] * pb_x[j] * fx[j] + 3.0 * pa_yzz[j] * fx[j] * pb_xz[j] + 

                                1.5 * pa_yz[j] * fx[j] * pb_xzz[j] + pa_yzzz[j] * pb_xzz[j]) * s_0_0[j] + (22.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                1.5 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * pa_yz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 

                                3.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] - pa_yzzz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                15.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 6.0 * pa_yzzz[j] * fz[j] * pb_x[j] * fx[j] + 

                                36.0 * pa_yzz[j] * fz[j] * fx[j] * pb_xz[j] - 3.0 * pa_yz[j] * fz[j] * fga[j] * pb_xzz[j] + 

                                18.0 * pa_yz[j] * fz[j] * fx[j] * pb_xzz[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_yzzz_yyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * pa_zzz[j] + 

                                2.25 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 1.5 * pa_yzzz[j] * pb_y[j] * fx[j] + 

                                1.5 * fx[j] * pa_zzz[j] * pb_yy[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyy[j] + pa_yzzz[j] * pb_yyy[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                2.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fga[j] - 1.5 * fx[j] * pa_zzz[j] * fz[j] * fgb[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 7.5 * fx[j] * fx[j] * pa_zzz[j] * fz[j] - 

                                4.5 * pa_yz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 4.5 * pa_yz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                4.5 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_yy[j] - 3.0 * pa_yzzz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                22.5 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 22.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yy[j] + 

                                18.0 * pa_yzzz[j] * fz[j] * pb_y[j] * fx[j] + 18.0 * fx[j] * pa_zzz[j] * fz[j] * pb_yy[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_yyy[j] + 18.0 * pa_yz[j] * fz[j] * fx[j] * pb_yyy[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_yzzz_yyz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_yzz[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 

                                0.75 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_yzzz[j] * fx[j] * pb_z[j] + 1.5 * pa_yzz[j] * fx[j] * pb_yy[j] + 

                                fx[j] * pa_zzz[j] * pb_yz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yyz[j] + pa_yzzz[j] * pb_yyz[j]) * s_0_0[j] + (-0.75 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_y[j] * fx[j] * fx[j] * fz[j] * fga[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                1.5 * pa_yzz[j] * fx[j] * fz[j] * fgb[j] + 3.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_yzz[j] * fz[j] * fx[j] * fx[j] + 

                                15.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] - 1.5 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                1.5 * pa_yz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 1.5 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] - 

                                3.0 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_yz[j] - pa_yzzz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                7.5 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 7.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 

                                15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yz[j] + 6.0 * pa_yzzz[j] * fz[j] * fx[j] * pb_z[j] + 

                                18.0 * pa_yzz[j] * fz[j] * fx[j] * pb_yy[j] + 12.0 * fx[j] * pa_zzz[j] * fz[j] * pb_yz[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_yyz[j] + 18.0 * pa_yz[j] * fz[j] * fx[j] * pb_yyz[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_yzzz_yzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                2.25 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pa_zzz[j] + 1.5 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 

                                1.5 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_yzzz[j] * pb_y[j] * fx[j] + 

                                3.0 * pa_yzz[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * pa_zzz[j] * pb_zz[j] + 1.5 * pa_yz[j] * fx[j] * pb_yzz[j] + 

                                pa_yzzz[j] * pb_yzz[j]) * s_0_0[j] + (9.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] - 

                                0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fga[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 0.5 * fx[j] * pa_zzz[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 22.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                2.5 * fx[j] * fx[j] * pa_zzz[j] * fz[j] + 15.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_yz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * pa_yz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 

                                3.0 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] - 1.5 * fx[j] * pa_z[j] * fz[j] * fga[j] * pb_zz[j] - 

                                pa_yzzz[j] * pb_y[j] * fz[j] * fgb[j] + 15.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 

                                7.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_zz[j] + 6.0 * pa_yzzz[j] * fz[j] * pb_y[j] * fx[j] + 

                                36.0 * pa_yzz[j] * fz[j] * fx[j] * pb_yz[j] + 6.0 * fx[j] * pa_zzz[j] * fz[j] * pb_zz[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_yzz[j] + 18.0 * pa_yz[j] * fz[j] * fx[j] * pb_yzz[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_yzzz_zzz[j] = (1.875 * pa_y[j] * fx[j] * fx[j] * fx[j] + 2.25 * pa_yzz[j] * fx[j] * fx[j] + 

                                6.75 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_yzzz[j] * pb_z[j] * fx[j] + 

                                4.5 * pa_yzz[j] * fx[j] * pb_zz[j] + 1.5 * pa_yz[j] * fx[j] * pb_zzz[j] + pa_yzzz[j] * pb_zzz[j]) * s_0_0[j] + (15.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                2.25 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 2.25 * pa_y[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                4.5 * pa_yzz[j] * fx[j] * fz[j] * fgb[j] + 22.5 * pa_yzz[j] * fz[j] * fx[j] * fx[j] + 

                                67.5 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 4.5 * pa_yz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                4.5 * pa_yz[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 4.5 * pa_y[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] - 

                                3.0 * pa_yzzz[j] * pb_z[j] * fz[j] * fgb[j] + 22.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 

                                18.0 * pa_yzzz[j] * fz[j] * pb_z[j] * fx[j] + 54.0 * pa_yzz[j] * fz[j] * fx[j] * pb_zz[j] - 

                                3.0 * pa_yz[j] * fz[j] * fga[j] * pb_zzz[j] + 18.0 * pa_yz[j] * fz[j] * fx[j] * pb_zzz[j] + 14.0 * pa_yzzz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (28) = (140,145)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xx, pb_xxx, pb_xxy, pb_xxz, \
                                     pb_xy, pb_xyy, pb_xyz, pb_y, pb_z, r_0_0, s_0_0, t_zzzz_xxx, t_zzzz_xxy, t_zzzz_xxz, \
                                     t_zzzz_xyy, t_zzzz_xyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xxx[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] + 

                                1.5 * pa_zzzz[j] * pb_x[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxx[j] + 

                                pa_zzzz[j] * pb_xxx[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 9.0 * pa_zz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                9.0 * pa_zz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                3.0 * pa_zzzz[j] * pb_x[j] * fz[j] * fgb[j] + 45.0 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                18.0 * pa_zzzz[j] * fz[j] * pb_x[j] * fx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xxx[j] - 

                                6.0 * pa_zz[j] * fz[j] * fga[j] * pb_xxx[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 

                                36.0 * pa_zz[j] * fz[j] * fx[j] * pb_xxx[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_xxx[j]) * r_0_0[j];

                t_zzzz_xxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.5 * pa_zzzz[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxy[j] + 

                                pa_zzzz[j] * pb_xxy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 3.0 * pa_zz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                3.0 * pa_zz[j] * fz[j] * fga[j] * fx[j] * pb_y[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                pa_zzzz[j] * fz[j] * fgb[j] * pb_y[j] + 15.0 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                6.0 * pa_zzzz[j] * fz[j] * fx[j] * pb_y[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xxy[j] - 

                                6.0 * pa_zz[j] * fz[j] * fga[j] * pb_xxy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 

                                36.0 * pa_zz[j] * fz[j] * fx[j] * pb_xxy[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_xxy[j]) * r_0_0[j];

                t_zzzz_xxz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * fx[j] + pa_zzz[j] * fx[j] * fx[j] + 

                                0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_z[j] * fx[j] * fx[j] * pb_xx[j] + 

                                0.5 * pa_zzzz[j] * fx[j] * pb_z[j] + 2.0 * pa_zzz[j] * fx[j] * pb_xx[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] + 

                                3.0 * pa_zz[j] * fx[j] * pb_xxz[j] + pa_zzzz[j] * pb_xxz[j]) * s_0_0[j] + (-3.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * fga[j] - 2.0 * pa_zzz[j] * fx[j] * fz[j] * fgb[j] + 

                                12.0 * pa_z[j] * fx[j] * fx[j] * fx[j] * fz[j] + 10.0 * pa_zzz[j] * fz[j] * fx[j] * fx[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                3.0 * pa_zz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 3.0 * pa_zz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                6.0 * pa_z[j] * fx[j] * fz[j] * fga[j] * pb_xx[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                pa_zzzz[j] * fz[j] * fgb[j] * pb_z[j] + 15.0 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                30.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 6.0 * pa_zzzz[j] * fz[j] * fx[j] * pb_z[j] + 

                                24.0 * pa_zzz[j] * fz[j] * fx[j] * pb_xx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xxz[j] - 

                                6.0 * pa_zz[j] * fz[j] * fga[j] * pb_xxz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 

                                36.0 * pa_zz[j] * fz[j] * fx[j] * pb_xxz[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_xxz[j]) * r_0_0[j];

                t_zzzz_xyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.5 * pa_zzzz[j] * pb_x[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] + 3.0 * pa_zz[j] * fx[j] * pb_xyy[j] + 

                                pa_zzzz[j] * pb_xyy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 3.0 * pa_zz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                3.0 * pa_zz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                pa_zzzz[j] * pb_x[j] * fz[j] * fgb[j] + 15.0 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                6.0 * pa_zzzz[j] * fz[j] * pb_x[j] * fx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xyy[j] - 

                                6.0 * pa_zz[j] * fz[j] * fga[j] * pb_xyy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 

                                36.0 * pa_zz[j] * fz[j] * fx[j] * pb_xyy[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_xyy[j]) * r_0_0[j];

                t_zzzz_xyz[j] = (3.0 * pa_z[j] * fx[j] * fx[j] * pb_xy[j] + 2.0 * pa_zzz[j] * fx[j] * pb_xy[j] + 

                                0.75 * fx[j] * fx[j] * pb_xyz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xyz[j] + pa_zzzz[j] * pb_xyz[j]) * s_0_0[j] + (-6.0 * pa_z[j] * fx[j] * fz[j] * fga[j] * pb_xy[j] + 

                                30.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 24.0 * pa_zzz[j] * fz[j] * fx[j] * pb_xy[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xyz[j] - 6.0 * pa_zz[j] * fz[j] * fga[j] * pb_xyz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 36.0 * pa_zz[j] * fz[j] * fx[j] * pb_xyz[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_xyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (29) = (145,150)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xz, pb_xzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_zzzz_xzz, \
                                     t_zzzz_yyy, t_zzzz_yyz, t_zzzz_yzz, t_zzzz_zzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xzz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] + 

                                6.0 * pa_z[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * pa_zzzz[j] * pb_x[j] * fx[j] + 4.0 * pa_zzz[j] * fx[j] * pb_xz[j] + 

                                0.75 * fx[j] * fx[j] * pb_xzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xzz[j] + pa_zzzz[j] * pb_xzz[j]) * s_0_0[j] + (-4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] + 

                                15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 45.0 * pa_zz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 3.0 * pa_zz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                3.0 * pa_zz[j] * fz[j] * fga[j] * pb_x[j] * fx[j] - 12.0 * pa_z[j] * fx[j] * fz[j] * fga[j] * pb_xz[j] - 

                                pa_zzzz[j] * pb_x[j] * fz[j] * fgb[j] + 60.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                6.0 * pa_zzzz[j] * fz[j] * pb_x[j] * fx[j] + 48.0 * pa_zzz[j] * fz[j] * fx[j] * pb_xz[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xzz[j] - 6.0 * pa_zz[j] * fz[j] * fga[j] * pb_xzz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 36.0 * pa_zz[j] * fz[j] * fx[j] * pb_xzz[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_xzz[j]) * r_0_0[j];

                t_zzzz_yyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_zzzz[j] * pb_y[j] * fx[j] + 0.75 * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_zz[j] * fx[j] * pb_yyy[j] + 

                                pa_zzzz[j] * pb_yyy[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 9.0 * pa_zz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                9.0 * pa_zz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                3.0 * pa_zzzz[j] * pb_y[j] * fz[j] * fgb[j] + 45.0 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                18.0 * pa_zzzz[j] * fz[j] * pb_y[j] * fx[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_yyy[j] - 

                                6.0 * pa_zz[j] * fz[j] * fga[j] * pb_yyy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 

                                36.0 * pa_zz[j] * fz[j] * fx[j] * pb_yyy[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_yyy[j]) * r_0_0[j];

                t_zzzz_yyz[j] = (1.5 * pa_z[j] * fx[j] * fx[j] * fx[j] + pa_zzz[j] * fx[j] * fx[j] + 

                                0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] + 3.0 * pa_z[j] * fx[j] * fx[j] * pb_yy[j] + 

                                0.5 * pa_zzzz[j] * fx[j] * pb_z[j] + 2.0 * pa_zzz[j] * fx[j] * pb_yy[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] + 

                                3.0 * pa_zz[j] * fx[j] * pb_yyz[j] + pa_zzzz[j] * pb_yyz[j]) * s_0_0[j] + (-3.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * fga[j] - 2.0 * pa_zzz[j] * fx[j] * fz[j] * fgb[j] + 

                                12.0 * pa_z[j] * fx[j] * fx[j] * fx[j] * fz[j] + 10.0 * pa_zzz[j] * fz[j] * fx[j] * fx[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                3.0 * pa_zz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 3.0 * pa_zz[j] * fz[j] * fga[j] * fx[j] * pb_z[j] - 

                                6.0 * pa_z[j] * fx[j] * fz[j] * fga[j] * pb_yy[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                pa_zzzz[j] * fz[j] * fgb[j] * pb_z[j] + 15.0 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                30.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 6.0 * pa_zzzz[j] * fz[j] * fx[j] * pb_z[j] + 

                                24.0 * pa_zzz[j] * fz[j] * fx[j] * pb_yy[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_yyz[j] - 

                                6.0 * pa_zz[j] * fz[j] * fga[j] * pb_yyz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 

                                36.0 * pa_zz[j] * fz[j] * fx[j] * pb_yyz[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_yyz[j]) * r_0_0[j];

                t_zzzz_yzz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] + 

                                6.0 * pa_z[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_zzzz[j] * pb_y[j] * fx[j] + 4.0 * pa_zzz[j] * fx[j] * pb_yz[j] + 

                                0.75 * fx[j] * fx[j] * pb_yzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_yzz[j] + pa_zzzz[j] * pb_yzz[j]) * s_0_0[j] + (-4.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] + 

                                15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 45.0 * pa_zz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 3.0 * pa_zz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                3.0 * pa_zz[j] * fz[j] * fga[j] * pb_y[j] * fx[j] - 12.0 * pa_z[j] * fx[j] * fz[j] * fga[j] * pb_yz[j] - 

                                pa_zzzz[j] * pb_y[j] * fz[j] * fgb[j] + 60.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 

                                6.0 * pa_zzzz[j] * fz[j] * pb_y[j] * fx[j] + 48.0 * pa_zzz[j] * fz[j] * fx[j] * pb_yz[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_yzz[j] - 6.0 * pa_zz[j] * fz[j] * fga[j] * pb_yzz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 36.0 * pa_zz[j] * fz[j] * fx[j] * pb_yzz[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_yzz[j]) * r_0_0[j];

                t_zzzz_zzz[j] = (7.5 * pa_z[j] * fx[j] * fx[j] * fx[j] + 5.625 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                3.0 * pa_zzz[j] * fx[j] * fx[j] + 13.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] + 

                                9.0 * pa_z[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_zzzz[j] * pb_z[j] * fx[j] + 6.0 * pa_zzz[j] * fx[j] * pb_zz[j] + 

                                0.75 * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_zz[j] * fx[j] * pb_zzz[j] + pa_zzzz[j] * pb_zzz[j]) * s_0_0[j] + (60.0 * pa_z[j] * fx[j] * fx[j] * fx[j] * fz[j] - 

                                9.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 9.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * fga[j] - 

                                13.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 6.0 * pa_zzz[j] * fx[j] * fz[j] * fgb[j] + 

                                45.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 30.0 * pa_zzz[j] * fz[j] * fx[j] * fx[j] + 

                                135.0 * pa_zz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 2.25 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                9.0 * pa_zz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 9.0 * pa_zz[j] * fz[j] * fga[j] * pb_z[j] * fx[j] - 

                                18.0 * pa_z[j] * fx[j] * fz[j] * fga[j] * pb_zz[j] - 3.0 * pa_zzzz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                90.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] + 18.0 * pa_zzzz[j] * fz[j] * pb_z[j] * fx[j] + 

                                72.0 * pa_zzz[j] * fz[j] * fx[j] * pb_zz[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_zzz[j] - 

                                6.0 * pa_zz[j] * fz[j] * fga[j] * pb_zzz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 

                                36.0 * pa_zz[j] * fz[j] * fx[j] * pb_zzz[j] + 14.0 * pa_zzzz[j] * fz[j] * pb_zzz[j]) * r_0_0[j];
            }

            idx++;
        }
    }


} // kinrecfunc namespace

