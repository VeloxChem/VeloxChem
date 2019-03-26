//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "KineticEnergyRecFuncForFG.hpp"

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

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, r_0_0, s_0_0, t_xxx_xxxx, t_xxx_xxxy, t_xxx_xxxz, t_xxx_xxyy, \
                                     t_xxx_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xxxx[j] = (5.625 * pa_x[j] * fx[j] * fx[j] * fx[j] + 7.5 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xxx[j] * fx[j] * fx[j] + 9.0 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 

                                13.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xxx[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_xx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_x[j] * fx[j] * pb_xxxx[j] + pa_xxx[j] * pb_xxxx[j]) * s_0_0[j] + (-13.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                45.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 60.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                2.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 9.0 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                9.0 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 3.0 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] - 

                                18.0 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 7.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] + 

                                90.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 135.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] - 

                                9.0 * pa_x[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 9.0 * pa_x[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                6.0 * fx[j] * fz[j] * fga[j] * pb_xxx[j] - 6.0 * pa_xxx[j] * pb_xx[j] * fz[j] * fgb[j] + 

                                30.0 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 36.0 * pa_xxx[j] * fz[j] * pb_xx[j] * fx[j] + 

                                72.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxx[j] - 3.0 * pa_x[j] * fz[j] * fga[j] * pb_xxxx[j] + 

                                18.0 * pa_x[j] * fz[j] * fx[j] * pb_xxxx[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_xxx_xxxy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                6.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 2.25 * fx[j] * fx[j] * pb_xxy[j] + 1.5 * pa_xxx[j] * pb_xy[j] * fx[j] + 

                                4.5 * pa_xx[j] * fx[j] * pb_xxy[j] + 1.5 * pa_x[j] * fx[j] * pb_xxxy[j] + pa_xxx[j] * pb_xxxy[j]) * s_0_0[j] + (15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                2.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                4.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] + 22.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                67.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] - 4.5 * pa_x[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                4.5 * pa_x[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 4.5 * fx[j] * fz[j] * fga[j] * pb_xxy[j] - 

                                3.0 * pa_xxx[j] * pb_xy[j] * fz[j] * fgb[j] + 22.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 

                                18.0 * pa_xxx[j] * fz[j] * pb_xy[j] * fx[j] + 54.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxy[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xxxy[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_xxxy[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_xxx_xxxz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                6.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 2.25 * fx[j] * fx[j] * pb_xxz[j] + 1.5 * pa_xxx[j] * pb_xz[j] * fx[j] + 

                                4.5 * pa_xx[j] * fx[j] * pb_xxz[j] + 1.5 * pa_x[j] * fx[j] * pb_xxxz[j] + pa_xxx[j] * pb_xxxz[j]) * s_0_0[j] + (15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                2.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                4.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 22.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                67.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] - 4.5 * pa_x[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                4.5 * pa_x[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 4.5 * fx[j] * fz[j] * fga[j] * pb_xxz[j] - 

                                3.0 * pa_xxx[j] * pb_xz[j] * fz[j] * fgb[j] + 22.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 

                                18.0 * pa_xxx[j] * fz[j] * pb_xz[j] * fx[j] + 54.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxz[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xxxz[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_xxxz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_xxx_xxyy[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xxx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * fx[j] * fx[j] * pb_xyy[j] + 

                                0.5 * pa_xxx[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xx[j] * fx[j] * pb_xyy[j] + 

                                1.5 * pa_x[j] * fx[j] * pb_xxyy[j] + pa_xxx[j] * pb_xxyy[j]) * s_0_0[j] + (-3.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - pa_xxx[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] + 15.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] - 1.5 * pa_x[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                1.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_yy[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xyy[j] - pa_xxx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                pa_xxx[j] * fz[j] * fgb[j] * pb_yy[j] + 7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 6.0 * pa_xxx[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_yy[j] + 36.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyy[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xxyy[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_xxyy[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_xxx_xxyz[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * fx[j] * fx[j] * pb_xyz[j] + 

                                0.5 * pa_xxx[j] * fx[j] * pb_yz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xyz[j] + 1.5 * pa_x[j] * fx[j] * pb_xxyz[j] + 

                                pa_xxx[j] * pb_xxyz[j]) * s_0_0[j] + (22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] - 

                                1.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_yz[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xyz[j] - pa_xxx[j] * fz[j] * fgb[j] * pb_yz[j] + 15.0 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_yz[j] + 36.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyz[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xxyz[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_xxyz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxx, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, \
                                     s_0_0, t_xxx_xxzz, t_xxx_xyyy, t_xxx_xyyz, t_xxx_xyzz, t_xxx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_xxzz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xxx[j] * fx[j] * fx[j] + 1.5 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 1.5 * fx[j] * fx[j] * pb_xzz[j] + 

                                0.5 * pa_xxx[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xx[j] * fx[j] * pb_xzz[j] + 

                                1.5 * pa_x[j] * fx[j] * pb_xxzz[j] + pa_xxx[j] * pb_xxzz[j]) * s_0_0[j] + (-3.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - pa_xxx[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] + 15.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 1.5 * pa_x[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                1.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xzz[j] - pa_xxx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                pa_xxx[j] * fz[j] * fgb[j] * pb_zz[j] + 7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 6.0 * pa_xxx[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_zz[j] + 36.0 * pa_xx[j] * fz[j] * fx[j] * pb_xzz[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xxzz[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_xxzz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_xxx_xyyy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xxx[j] * pb_xy[j] * fx[j] + 

                                1.5 * pa_xx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_x[j] * fx[j] * pb_xyyy[j] + pa_xxx[j] * pb_xyyy[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 4.5 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 22.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] - 

                                4.5 * pa_x[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 4.5 * pa_x[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_yyy[j] - 3.0 * pa_xxx[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                22.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 

                                18.0 * pa_xxx[j] * fz[j] * pb_xy[j] * fx[j] + 18.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyy[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xyyy[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_xyyy[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_xxx_xyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_xxx[j] * pb_xz[j] * fx[j] + 

                                1.5 * pa_xx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_x[j] * fx[j] * pb_xyyz[j] + pa_xxx[j] * pb_xyyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 1.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_yyz[j] - pa_xxx[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 

                                6.0 * pa_xxx[j] * fz[j] * pb_xz[j] * fx[j] + 18.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyz[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xyyz[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_xyyz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_xxx_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xxx[j] * pb_xy[j] * fx[j] + 

                                1.5 * pa_xx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_x[j] * fx[j] * pb_xyzz[j] + pa_xxx[j] * pb_xyzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 1.5 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_yzz[j] - pa_xxx[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 

                                6.0 * pa_xxx[j] * fz[j] * pb_xy[j] * fx[j] + 18.0 * pa_xx[j] * fz[j] * fx[j] * pb_yzz[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xyzz[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_xyzz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_xxx_xzzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xxx[j] * pb_xz[j] * fx[j] + 

                                1.5 * pa_xx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_x[j] * fx[j] * pb_xzzz[j] + pa_xxx[j] * pb_xzzz[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 4.5 * pa_xx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 22.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 

                                4.5 * pa_x[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 4.5 * pa_x[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_zzz[j] - 3.0 * pa_xxx[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                22.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 

                                18.0 * pa_xxx[j] * fz[j] * pb_xz[j] * fx[j] + 18.0 * pa_xx[j] * fz[j] * fx[j] * pb_zzz[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xzzz[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_xzzz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xxx, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, \
                                     pb_zz, pb_zzzz, r_0_0, s_0_0, t_xxx_yyyy, t_xxx_yyyz, t_xxx_yyzz, t_xxx_yzzz, \
                                     t_xxx_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxx_yyyy[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] + 

                                4.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xxx[j] * pb_yy[j] * fx[j] + 1.5 * pa_x[j] * fx[j] * pb_yyyy[j] + 

                                pa_xxx[j] * pb_yyyy[j]) * s_0_0[j] + (-4.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                2.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * fx[j] + 7.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] - 

                                9.0 * pa_x[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 9.0 * pa_x[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                6.0 * pa_xxx[j] * pb_yy[j] * fz[j] * fgb[j] + 45.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                36.0 * pa_xxx[j] * fz[j] * pb_yy[j] * fx[j] - 3.0 * pa_x[j] * fz[j] * fga[j] * pb_yyyy[j] + 

                                18.0 * pa_x[j] * fz[j] * fx[j] * pb_yyyy[j] + 14.0 * pa_xxx[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_xxx_yyyz[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xxx[j] * pb_yz[j] * fx[j] + 

                                1.5 * pa_x[j] * fx[j] * pb_yyyz[j] + pa_xxx[j] * pb_yyyz[j]) * s_0_0[j] + (-4.5 * pa_x[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 

                                4.5 * pa_x[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 3.0 * pa_xxx[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                22.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 18.0 * pa_xxx[j] * fz[j] * pb_yz[j] * fx[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_yyyz[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_yyyz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_xxx_yyzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xxx[j] * fx[j] * fx[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xxx[j] * pb_yy[j] * fx[j] + 

                                0.5 * pa_xxx[j] * fx[j] * pb_zz[j] + 1.5 * pa_x[j] * fx[j] * pb_yyzz[j] + pa_xxx[j] * pb_yyzz[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 3.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * fx[j] + 

                                2.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] - 1.5 * pa_x[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                1.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - pa_xxx[j] * pb_yy[j] * fz[j] * fgb[j] - pa_xxx[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                6.0 * pa_xxx[j] * fz[j] * pb_yy[j] * fx[j] + 6.0 * pa_xxx[j] * fz[j] * fx[j] * pb_zz[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_yyzz[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_yyzz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_xxx_yzzz[j] = (2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xxx[j] * pb_yz[j] * fx[j] + 

                                1.5 * pa_x[j] * fx[j] * pb_yzzz[j] + pa_xxx[j] * pb_yzzz[j]) * s_0_0[j] + (-4.5 * pa_x[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 

                                4.5 * pa_x[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 3.0 * pa_xxx[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                22.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 18.0 * pa_xxx[j] * fz[j] * pb_yz[j] * fx[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_yzzz[j] + 18.0 * pa_x[j] * fz[j] * fx[j] * pb_yzzz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_xxx_zzzz[j] = (1.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xxx[j] * fx[j] * fx[j] + 

                                4.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xxx[j] * pb_zz[j] * fx[j] + 1.5 * pa_x[j] * fx[j] * pb_zzzz[j] + 

                                pa_xxx[j] * pb_zzzz[j]) * s_0_0[j] + (-4.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                2.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * pa_xxx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * fx[j] + 7.5 * pa_xxx[j] * fz[j] * fx[j] * fx[j] - 

                                9.0 * pa_x[j] * fx[j] * pb_zz[j] * fz[j] * fgb[j] - 9.0 * pa_x[j] * fz[j] * fga[j] * pb_zz[j] * fx[j] - 

                                6.0 * pa_xxx[j] * pb_zz[j] * fz[j] * fgb[j] + 45.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                36.0 * pa_xxx[j] * fz[j] * pb_zz[j] * fx[j] - 3.0 * pa_x[j] * fz[j] * fga[j] * pb_zzzz[j] + 

                                18.0 * pa_x[j] * fz[j] * fx[j] * pb_zzzz[j] + 14.0 * pa_xxx[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xxy_xxxx, t_xxy_xxxy, t_xxy_xxxz, t_xxy_xxyy, \
                                     t_xxy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xxxx[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] + 

                                6.0 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 3.0 * pa_xxy[j] * pb_xx[j] * fx[j] + 

                                4.0 * pa_xy[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pa_y[j] * pb_xxxx[j] + pa_xxy[j] * pb_xxxx[j]) * s_0_0[j] + (-4.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] + 

                                15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_y[j] - 0.75 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 

                                3.0 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] - 12.0 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xxy[j] * fz[j] * fx[j] * fx[j] + 60.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                45.0 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xx[j] - 3.0 * fx[j] * pa_y[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                3.0 * fz[j] * fga[j] * pa_y[j] * pb_xx[j] * fx[j] - 6.0 * pa_xxy[j] * pb_xx[j] * fz[j] * fgb[j] + 

                                36.0 * pa_xxy[j] * fz[j] * pb_xx[j] * fx[j] + 48.0 * pa_xy[j] * fx[j] * fz[j] * pb_xxx[j] - fz[j] * fga[j] * pa_y[j] * pb_xxxx[j] + 

                                6.0 * fx[j] * fz[j] * pa_y[j] * pb_xxxx[j] + 14.0 * pa_xxy[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_xxy_xxxy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                2.25 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_xxy[j] * pb_xy[j] * fx[j] + 

                                0.5 * pa_xx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xy[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pa_y[j] * pb_xxxy[j] + 

                                pa_xxy[j] * pb_xxxy[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                6.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - 

                                1.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 3.0 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] + 

                                7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 15.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 22.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xy[j] - 

                                1.5 * fx[j] * pa_y[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_y[j] * pb_xy[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pb_xxx[j] - 3.0 * pa_xxy[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 18.0 * pa_xxy[j] * fz[j] * pb_xy[j] * fx[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxx[j] + 36.0 * pa_xy[j] * fx[j] * fz[j] * pb_xxy[j] - fz[j] * fga[j] * pa_y[j] * pb_xxxy[j] + 

                                6.0 * fx[j] * fz[j] * pa_y[j] * pb_xxxy[j] + 14.0 * pa_xxy[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_xxy_xxxz[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                2.25 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 1.5 * pa_xxy[j] * pb_xz[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * pb_xxz[j] + 

                                0.5 * fx[j] * pa_y[j] * pb_xxxz[j] + pa_xxy[j] * pb_xxxz[j]) * s_0_0[j] + (-3.0 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                15.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 22.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xz[j] - 

                                1.5 * fx[j] * pa_y[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_y[j] * pb_xz[j] * fx[j] - 

                                3.0 * pa_xxy[j] * pb_xz[j] * fz[j] * fgb[j] + 18.0 * pa_xxy[j] * fz[j] * pb_xz[j] * fx[j] + 

                                36.0 * pa_xy[j] * fx[j] * fz[j] * pb_xxz[j] - fz[j] * fga[j] * pa_y[j] * pb_xxxz[j] + 6.0 * fx[j] * fz[j] * pa_y[j] * pb_xxxz[j] + 

                                14.0 * pa_xxy[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_xxy_xxyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_xxy[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.0 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_xxy[j] * pb_xx[j] * fx[j] + 

                                0.5 * pa_xxy[j] * fx[j] * pb_yy[j] + pa_xx[j] * fx[j] * pb_xxy[j] + 2.0 * pa_xy[j] * fx[j] * pb_xyy[j] + 

                                0.5 * fx[j] * pa_y[j] * pb_xxyy[j] + pa_xxy[j] * pb_xxyy[j]) * s_0_0[j] + (-fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_y[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 0.25 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - pa_xxy[j] * fx[j] * fz[j] * fgb[j] - pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                2.0 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 2.5 * pa_xxy[j] * fz[j] * fx[j] * fx[j] + 

                                5.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 10.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                20.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_yy[j] - 

                                0.5 * fx[j] * pa_y[j] * pb_xx[j] * fz[j] * fgb[j] - 0.5 * fx[j] * pa_y[j] * fz[j] * fgb[j] * pb_yy[j] - 

                                0.5 * fz[j] * fga[j] * pa_y[j] * pb_xx[j] * fx[j] - 0.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_yy[j] - 

                                fz[j] * fga[j] * fx[j] * pb_xxy[j] - pa_xxy[j] * pb_xx[j] * fz[j] * fgb[j] - pa_xxy[j] * fz[j] * fgb[j] * pb_yy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xx[j] + 5.0 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 

                                6.0 * pa_xxy[j] * fz[j] * pb_xx[j] * fx[j] + 6.0 * pa_xxy[j] * fz[j] * fx[j] * pb_yy[j] + 

                                12.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxy[j] + 24.0 * pa_xy[j] * fx[j] * fz[j] * pb_xyy[j] - fz[j] * fga[j] * pa_y[j] * pb_xxyy[j] + 

                                6.0 * fx[j] * fz[j] * pa_y[j] * pb_xxyy[j] + 14.0 * pa_xxy[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_xxy_xxyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] + 

                                0.5 * pa_xxy[j] * fx[j] * pb_yz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxz[j] + 2.0 * pa_xy[j] * fx[j] * pb_xyz[j] + 

                                0.5 * fx[j] * pa_y[j] * pb_xxyz[j] + pa_xxy[j] * pb_xxyz[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.25 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - 

                                0.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                10.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_yz[j] - 

                                0.5 * fx[j] * pa_y[j] * fz[j] * fgb[j] * pb_yz[j] - 0.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_yz[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pb_xxz[j] - pa_xxy[j] * fz[j] * fgb[j] * pb_yz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 

                                6.0 * pa_xxy[j] * fz[j] * fx[j] * pb_yz[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxz[j] + 

                                24.0 * pa_xy[j] * fx[j] * fz[j] * pb_xyz[j] - fz[j] * fga[j] * pa_y[j] * pb_xxyz[j] + 6.0 * fx[j] * fz[j] * pa_y[j] * pb_xxyz[j] + 

                                14.0 * pa_xxy[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxy_xxzz, t_xxy_xyyy, \
                                     t_xxy_xyyz, t_xxy_xyzz, t_xxy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_xxzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * pa_xxy[j] * fx[j] * fx[j] + 

                                pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.5 * pa_xxy[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxy[j] * fx[j] * pb_zz[j] + 

                                2.0 * pa_xy[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_y[j] * pb_xxzz[j] + pa_xxy[j] * pb_xxzz[j]) * s_0_0[j] + (-fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_y[j] - 0.25 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 

                                pa_xxy[j] * fx[j] * fz[j] * fgb[j] - 2.0 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xxy[j] * fz[j] * fx[j] * fx[j] + 10.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_zz[j] - 0.5 * fx[j] * pa_y[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * pa_y[j] * fz[j] * fgb[j] * pb_zz[j] - 0.5 * fz[j] * fga[j] * pa_y[j] * pb_xx[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_zz[j] - pa_xxy[j] * pb_xx[j] * fz[j] * fgb[j] - pa_xxy[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xx[j] + 6.0 * pa_xxy[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_xxy[j] * fz[j] * fx[j] * pb_zz[j] + 24.0 * pa_xy[j] * fx[j] * fz[j] * pb_xzz[j] - fz[j] * fga[j] * pa_y[j] * pb_xxzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_y[j] * pb_xxzz[j] + 14.0 * pa_xxy[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_xxy_xyyy[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] + 1.5 * pa_xxy[j] * pb_xy[j] * fx[j] + 

                                1.5 * pa_xx[j] * fx[j] * pb_xyy[j] + pa_xy[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_y[j] * pb_xyyy[j] + 

                                pa_xxy[j] * pb_xyyy[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                6.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - 1.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                3.0 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 15.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] - 1.5 * fx[j] * pa_y[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_y[j] * pb_xy[j] * fx[j] - 1.5 * fz[j] * fga[j] * fx[j] * pb_xyy[j] - 

                                3.0 * pa_xxy[j] * pb_xy[j] * fz[j] * fgb[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xy[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 18.0 * pa_xxy[j] * fz[j] * pb_xy[j] * fx[j] + 

                                18.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyy[j] + 12.0 * pa_xy[j] * fx[j] * fz[j] * pb_yyy[j] - fz[j] * fga[j] * pa_y[j] * pb_xyyy[j] + 

                                6.0 * fx[j] * fz[j] * pa_y[j] * pb_xyyy[j] + 14.0 * pa_xxy[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_xxy_xyyz[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xxy[j] * pb_xz[j] * fx[j] + 

                                pa_xx[j] * fx[j] * pb_xyz[j] + pa_xy[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pa_y[j] * pb_xyyz[j] + 

                                pa_xxy[j] * pb_xyyz[j]) * s_0_0[j] + (-pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                5.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 10.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] - 

                                0.5 * fx[j] * pa_y[j] * pb_xz[j] * fz[j] * fgb[j] - 0.5 * fz[j] * fga[j] * pa_y[j] * pb_xz[j] * fx[j] - 

                                fz[j] * fga[j] * fx[j] * pb_xyz[j] - pa_xxy[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xz[j] + 5.0 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * pa_xxy[j] * fz[j] * pb_xz[j] * fx[j] + 12.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyz[j] + 

                                12.0 * pa_xy[j] * fx[j] * fz[j] * pb_yyz[j] - fz[j] * fga[j] * pa_y[j] * pb_xyyz[j] + 6.0 * fx[j] * fz[j] * pa_y[j] * pb_xyyz[j] + 

                                14.0 * pa_xxy[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_xxy_xyzz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_xxy[j] * pb_xy[j] * fx[j] + 

                                0.5 * pa_xx[j] * fx[j] * pb_xzz[j] + pa_xy[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_y[j] * pb_xyzz[j] + 

                                pa_xxy[j] * pb_xyzz[j]) * s_0_0[j] + (-0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                2.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.25 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - 0.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 5.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 0.5 * fx[j] * pa_y[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_y[j] * pb_xy[j] * fx[j] - 0.5 * fz[j] * fga[j] * fx[j] * pb_xzz[j] - pa_xxy[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 

                                6.0 * pa_xxy[j] * fz[j] * pb_xy[j] * fx[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xzz[j] + 

                                12.0 * pa_xy[j] * fx[j] * fz[j] * pb_yzz[j] - fz[j] * fga[j] * pa_y[j] * pb_xyzz[j] + 6.0 * fx[j] * fz[j] * pa_y[j] * pb_xyzz[j] + 

                                14.0 * pa_xxy[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_xxy_xzzz[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 1.5 * pa_xxy[j] * pb_xz[j] * fx[j] + pa_xy[j] * fx[j] * pb_zzz[j] + 

                                0.5 * fx[j] * pa_y[j] * pb_xzzz[j] + pa_xxy[j] * pb_xzzz[j]) * s_0_0[j] + (-3.0 * pa_xy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 

                                15.0 * pa_xy[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 1.5 * fx[j] * pa_y[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_y[j] * pb_xz[j] * fx[j] - 3.0 * pa_xxy[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_xz[j] + 18.0 * pa_xxy[j] * fz[j] * pb_xz[j] * fx[j] + 

                                12.0 * pa_xy[j] * fx[j] * fz[j] * pb_zzz[j] - fz[j] * fga[j] * pa_y[j] * pb_xzzz[j] + 6.0 * fx[j] * fz[j] * pa_y[j] * pb_xzzz[j] + 

                                14.0 * pa_xxy[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_xx, pa_xxy, pa_y, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xxy_yyyy, \
                                     t_xxy_yyyz, t_xxy_yyzz, t_xxy_yzzz, t_xxy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxy_yyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_xxy[j] * fx[j] * fx[j] + 3.0 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 

                                fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_xxy[j] * pb_yy[j] * fx[j] + 2.0 * pa_xx[j] * fx[j] * pb_yyy[j] + 

                                0.5 * fx[j] * pa_y[j] * pb_yyyy[j] + pa_xxy[j] * pb_yyyy[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                3.0 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 0.75 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 

                                3.0 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - 3.0 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] - 

                                6.0 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_y[j] + 

                                12.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_xxy[j] * fz[j] * fx[j] * fx[j] + 

                                30.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] - 3.0 * fx[j] * pa_y[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                3.0 * fz[j] * fga[j] * pa_y[j] * pb_yy[j] * fx[j] - 2.0 * fz[j] * fga[j] * fx[j] * pb_yyy[j] - 

                                6.0 * pa_xxy[j] * pb_yy[j] * fz[j] * fgb[j] + 15.0 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_yy[j] + 

                                10.0 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 36.0 * pa_xxy[j] * fz[j] * pb_yy[j] * fx[j] + 

                                24.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyy[j] - fz[j] * fga[j] * pa_y[j] * pb_yyyy[j] + 6.0 * fx[j] * fz[j] * pa_y[j] * pb_yyyy[j] + 

                                14.0 * pa_xxy[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_xxy_yyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_xxy[j] * pb_yz[j] * fx[j] + 

                                1.5 * pa_xx[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pa_y[j] * pb_yyyz[j] + pa_xxy[j] * pb_yyyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - 1.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 

                                1.5 * fx[j] * pa_y[j] * pb_yz[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_y[j] * pb_yz[j] * fx[j] - 

                                1.5 * fz[j] * fga[j] * fx[j] * pb_yyz[j] - 3.0 * pa_xxy[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_yz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 

                                18.0 * pa_xxy[j] * fz[j] * pb_yz[j] * fx[j] + 18.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyz[j] - fz[j] * fga[j] * pa_y[j] * pb_yyyz[j] + 

                                6.0 * fx[j] * fz[j] * pa_y[j] * pb_yyyz[j] + 14.0 * pa_xxy[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_xxy_yyzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_xxy[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_yzz[j] + 

                                0.5 * pa_xxy[j] * pb_yy[j] * fx[j] + 0.5 * pa_xxy[j] * fx[j] * pb_zz[j] + pa_xx[j] * fx[j] * pb_yzz[j] + 

                                0.5 * fx[j] * pa_y[j] * pb_yyzz[j] + pa_xxy[j] * pb_yyzz[j]) * s_0_0[j] + (-0.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 0.25 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - pa_xxy[j] * fx[j] * fz[j] * fgb[j] - pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                fx[j] * fx[j] * fx[j] * fz[j] * pa_y[j] + 2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                2.5 * pa_xxy[j] * fz[j] * fx[j] * fx[j] + 5.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] - 

                                0.5 * fx[j] * pa_y[j] * pb_yy[j] * fz[j] * fgb[j] - 0.5 * fx[j] * pa_y[j] * fz[j] * fgb[j] * pb_zz[j] - 

                                0.5 * fz[j] * fga[j] * pa_y[j] * pb_yy[j] * fx[j] - 0.5 * fz[j] * fga[j] * pa_y[j] * fx[j] * pb_zz[j] - 

                                fz[j] * fga[j] * fx[j] * pb_yzz[j] - pa_xxy[j] * pb_yy[j] * fz[j] * fgb[j] - pa_xxy[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_yy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_zz[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 6.0 * pa_xxy[j] * fz[j] * pb_yy[j] * fx[j] + 

                                6.0 * pa_xxy[j] * fz[j] * fx[j] * pb_zz[j] + 12.0 * pa_xx[j] * fz[j] * fx[j] * pb_yzz[j] - fz[j] * fga[j] * pa_y[j] * pb_yyzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_y[j] * pb_yyzz[j] + 14.0 * pa_xxy[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_xxy_yzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xxy[j] * pb_yz[j] * fx[j] + 

                                0.5 * pa_xx[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_y[j] * pb_yzzz[j] + pa_xxy[j] * pb_yzzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - 1.5 * pa_xx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 

                                1.5 * fx[j] * pa_y[j] * pb_yz[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_y[j] * pb_yz[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pb_zzz[j] - 3.0 * pa_xxy[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_yz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 

                                18.0 * pa_xxy[j] * fz[j] * pb_yz[j] * fx[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_zzz[j] - fz[j] * fga[j] * pa_y[j] * pb_yzzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_y[j] * pb_yzzz[j] + 14.0 * pa_xxy[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_xxy_zzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xxy[j] * fx[j] * fx[j] + 

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 3.0 * pa_xxy[j] * pb_zz[j] * fx[j] + 0.5 * fx[j] * pa_y[j] * pb_zzzz[j] + 

                                pa_xxy[j] * pb_zzzz[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * pa_y[j] * fx[j] * fx[j] - 3.0 * pa_xxy[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_y[j] + 7.5 * pa_xxy[j] * fz[j] * fx[j] * fx[j] - 

                                3.0 * fx[j] * pa_y[j] * pb_zz[j] * fz[j] * fgb[j] - 3.0 * fz[j] * fga[j] * pa_y[j] * pb_zz[j] * fx[j] - 

                                6.0 * pa_xxy[j] * pb_zz[j] * fz[j] * fgb[j] + 15.0 * fx[j] * fx[j] * fz[j] * pa_y[j] * pb_zz[j] + 

                                36.0 * pa_xxy[j] * fz[j] * pb_zz[j] * fx[j] - fz[j] * fga[j] * pa_y[j] * pb_zzzz[j] + 6.0 * fx[j] * fz[j] * pa_y[j] * pb_zzzz[j] + 

                                14.0 * pa_xxy[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xxz_xxxx, t_xxz_xxxy, t_xxz_xxxz, t_xxz_xxyy, \
                                     t_xxz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xxxx[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] + 

                                6.0 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 4.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 3.0 * pa_xxz[j] * pb_xx[j] * fx[j] + 

                                4.0 * pa_xz[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pa_z[j] * pb_xxxx[j] + pa_xxz[j] * pb_xxxx[j]) * s_0_0[j] + (-4.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 0.75 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                3.0 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] - 12.0 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xxz[j] * fz[j] * fx[j] * fx[j] + 60.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                45.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] - 3.0 * fx[j] * pa_z[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                3.0 * fz[j] * fga[j] * pa_z[j] * pb_xx[j] * fx[j] - 6.0 * pa_xxz[j] * pb_xx[j] * fz[j] * fgb[j] + 

                                36.0 * pa_xxz[j] * fz[j] * pb_xx[j] * fx[j] + 48.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxx[j] - fz[j] * fga[j] * pa_z[j] * pb_xxxx[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxxx[j] + 14.0 * pa_xxz[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_xxz_xxxy[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                2.25 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 1.5 * pa_xxz[j] * pb_xy[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * pb_xxy[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_xxxy[j] + pa_xxz[j] * pb_xxxy[j]) * s_0_0[j] + (-3.0 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] + 

                                15.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 22.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] - 

                                1.5 * fx[j] * pa_z[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_z[j] * pb_xy[j] * fx[j] - 

                                3.0 * pa_xxz[j] * pb_xy[j] * fz[j] * fgb[j] + 18.0 * pa_xxz[j] * fz[j] * pb_xy[j] * fx[j] + 

                                36.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxy[j] - fz[j] * fga[j] * pa_z[j] * pb_xxxy[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxxy[j] + 

                                14.0 * pa_xxz[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_xxz_xxxz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                2.25 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_xxz[j] * pb_xz[j] * fx[j] + 

                                0.5 * pa_xx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xz[j] * fx[j] * pb_xxz[j] + 0.5 * fx[j] * pa_z[j] * pb_xxxz[j] + 

                                pa_xxz[j] * pb_xxxz[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                6.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - 

                                1.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 3.0 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 15.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 22.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] - 

                                1.5 * fx[j] * pa_z[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_z[j] * pb_xz[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pb_xxx[j] - 3.0 * pa_xxz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 18.0 * pa_xxz[j] * fz[j] * pb_xz[j] * fx[j] + 

                                6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxx[j] + 36.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxz[j] - fz[j] * fga[j] * pa_z[j] * pb_xxxz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxxz[j] + 14.0 * pa_xxz[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_xxz_xxyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * pa_xxz[j] * fx[j] * fx[j] + 

                                pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * pa_xxz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xxz[j] * fx[j] * pb_yy[j] + 

                                2.0 * pa_xz[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_z[j] * pb_xxyy[j] + pa_xxz[j] * pb_xxyy[j]) * s_0_0[j] + (-fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 0.25 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                pa_xxz[j] * fx[j] * fz[j] * fgb[j] - 2.0 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xxz[j] * fz[j] * fx[j] * fx[j] + 10.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] - 0.5 * fx[j] * pa_z[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * pa_z[j] * fz[j] * fgb[j] * pb_yy[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * pb_xx[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_yy[j] - pa_xxz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_xxz[j] * fz[j] * fgb[j] * pb_yy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] + 6.0 * pa_xxz[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_xxz[j] * fz[j] * fx[j] * pb_yy[j] + 24.0 * pa_xz[j] * fx[j] * fz[j] * pb_xyy[j] - fz[j] * fga[j] * pa_z[j] * pb_xxyy[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxyy[j] + 14.0 * pa_xxz[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_xxz_xxyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] + 

                                0.5 * pa_xxz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xx[j] * fx[j] * pb_xxy[j] + 2.0 * pa_xz[j] * fx[j] * pb_xyz[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_xxyz[j] + pa_xxz[j] * pb_xxyz[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 0.25 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - 

                                0.5 * pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] + 2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                10.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] - 

                                0.5 * fx[j] * pa_z[j] * fz[j] * fgb[j] * pb_yz[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_yz[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pb_xxy[j] - pa_xxz[j] * fz[j] * fgb[j] * pb_yz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 

                                6.0 * pa_xxz[j] * fz[j] * fx[j] * pb_yz[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxy[j] + 

                                24.0 * pa_xz[j] * fx[j] * fz[j] * pb_xyz[j] - fz[j] * fga[j] * pa_z[j] * pb_xxyz[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxyz[j] + 

                                14.0 * pa_xxz[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xx, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xxz_xxzz, \
                                     t_xxz_xyyy, t_xxz_xyyz, t_xxz_xyzz, t_xxz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_xxzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_xxz[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.0 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_xxz[j] * pb_xx[j] * fx[j] + 

                                0.5 * pa_xxz[j] * fx[j] * pb_zz[j] + pa_xx[j] * fx[j] * pb_xxz[j] + 2.0 * pa_xz[j] * fx[j] * pb_xzz[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_xxzz[j] + pa_xxz[j] * pb_xxzz[j]) * s_0_0[j] + (-fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.25 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - pa_xxz[j] * fx[j] * fz[j] * fgb[j] - pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                2.0 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 2.5 * pa_xxz[j] * fz[j] * fx[j] * fx[j] + 

                                5.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 10.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                20.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] - 

                                0.5 * fx[j] * pa_z[j] * pb_xx[j] * fz[j] * fgb[j] - 0.5 * fx[j] * pa_z[j] * fz[j] * fgb[j] * pb_zz[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * pb_xx[j] * fx[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_zz[j] - 

                                fz[j] * fga[j] * fx[j] * pb_xxz[j] - pa_xxz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_xxz[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] + 5.0 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 

                                6.0 * pa_xxz[j] * fz[j] * pb_xx[j] * fx[j] + 6.0 * pa_xxz[j] * fz[j] * fx[j] * pb_zz[j] + 

                                12.0 * pa_xx[j] * fz[j] * fx[j] * pb_xxz[j] + 24.0 * pa_xz[j] * fx[j] * fz[j] * pb_xzz[j] - fz[j] * fga[j] * pa_z[j] * pb_xxzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxzz[j] + 14.0 * pa_xxz[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_xxz_xyyy[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 1.5 * pa_xxz[j] * pb_xy[j] * fx[j] + pa_xz[j] * fx[j] * pb_yyy[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_xyyy[j] + pa_xxz[j] * pb_xyyy[j]) * s_0_0[j] + (-3.0 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                15.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 1.5 * fx[j] * pa_z[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_z[j] * pb_xy[j] * fx[j] - 3.0 * pa_xxz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] + 18.0 * pa_xxz[j] * fz[j] * pb_xy[j] * fx[j] + 

                                12.0 * pa_xz[j] * fx[j] * fz[j] * pb_yyy[j] - fz[j] * fga[j] * pa_z[j] * pb_xyyy[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xyyy[j] + 

                                14.0 * pa_xxz[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_xxz_xyyz[j] = (0.25 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_xxz[j] * pb_xz[j] * fx[j] + 

                                0.5 * pa_xx[j] * fx[j] * pb_xyy[j] + pa_xz[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pa_z[j] * pb_xyyz[j] + 

                                pa_xxz[j] * pb_xyyz[j]) * s_0_0[j] + (-0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                2.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.25 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - 0.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 5.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] - 0.5 * fx[j] * pa_z[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * pb_xz[j] * fx[j] - 0.5 * fz[j] * fga[j] * fx[j] * pb_xyy[j] - pa_xxz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 

                                6.0 * pa_xxz[j] * fz[j] * pb_xz[j] * fx[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyy[j] + 

                                12.0 * pa_xz[j] * fx[j] * fz[j] * pb_yyz[j] - fz[j] * fga[j] * pa_z[j] * pb_xyyz[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xyyz[j] + 

                                14.0 * pa_xxz[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_xxz_xyzz[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xxz[j] * pb_xy[j] * fx[j] + 

                                pa_xx[j] * fx[j] * pb_xyz[j] + pa_xz[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_z[j] * pb_xyzz[j] + 

                                pa_xxz[j] * pb_xyzz[j]) * s_0_0[j] + (-pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                5.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 10.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] - 

                                0.5 * fx[j] * pa_z[j] * pb_xy[j] * fz[j] * fgb[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * pb_xy[j] * fx[j] - 

                                fz[j] * fga[j] * fx[j] * pb_xyz[j] - pa_xxz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] + 5.0 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * pa_xxz[j] * fz[j] * pb_xy[j] * fx[j] + 12.0 * pa_xx[j] * fz[j] * fx[j] * pb_xyz[j] + 

                                12.0 * pa_xz[j] * fx[j] * fz[j] * pb_yzz[j] - fz[j] * fga[j] * pa_z[j] * pb_xyzz[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xyzz[j] + 

                                14.0 * pa_xxz[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_xxz_xzzz[j] = (0.75 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xx[j] * fx[j] * fx[j] * pb_x[j] + 1.5 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] + 1.5 * pa_xxz[j] * pb_xz[j] * fx[j] + 

                                1.5 * pa_xx[j] * fx[j] * pb_xzz[j] + pa_xz[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_z[j] * pb_xzzz[j] + 

                                pa_xxz[j] * pb_xzzz[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                6.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - 1.5 * pa_xx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                3.0 * pa_xz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 15.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 1.5 * fx[j] * pa_z[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_z[j] * pb_xz[j] * fx[j] - 1.5 * fz[j] * fga[j] * fx[j] * pb_xzz[j] - 

                                3.0 * pa_xxz[j] * pb_xz[j] * fz[j] * fgb[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 18.0 * pa_xxz[j] * fz[j] * pb_xz[j] * fx[j] + 

                                18.0 * pa_xx[j] * fz[j] * fx[j] * pb_xzz[j] + 12.0 * pa_xz[j] * fx[j] * fz[j] * pb_zzz[j] - fz[j] * fga[j] * pa_z[j] * pb_xzzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_xzzz[j] + 14.0 * pa_xxz[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_xx, pa_xxz, pa_z, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xxz_yyyy, \
                                     t_xxz_yyyz, t_xxz_yyzz, t_xxz_yzzz, t_xxz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxz_yyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xxz[j] * fx[j] * fx[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 3.0 * pa_xxz[j] * pb_yy[j] * fx[j] + 0.5 * fx[j] * pa_z[j] * pb_yyyy[j] + 

                                pa_xxz[j] * pb_yyyy[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 3.0 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 7.5 * pa_xxz[j] * fz[j] * fx[j] * fx[j] - 

                                3.0 * fx[j] * pa_z[j] * pb_yy[j] * fz[j] * fgb[j] - 3.0 * fz[j] * fga[j] * pa_z[j] * pb_yy[j] * fx[j] - 

                                6.0 * pa_xxz[j] * pb_yy[j] * fz[j] * fgb[j] + 15.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] + 

                                36.0 * pa_xxz[j] * fz[j] * pb_yy[j] * fx[j] - fz[j] * fga[j] * pa_z[j] * pb_yyyy[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_yyyy[j] + 

                                14.0 * pa_xxz[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_xxz_yyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xxz[j] * pb_yz[j] * fx[j] + 

                                0.5 * pa_xx[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_z[j] * pb_yyyz[j] + pa_xxz[j] * pb_yyyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - 1.5 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] - 

                                1.5 * fx[j] * pa_z[j] * pb_yz[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_z[j] * pb_yz[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pb_yyy[j] - 3.0 * pa_xxz[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 

                                18.0 * pa_xxz[j] * fz[j] * pb_yz[j] * fx[j] + 6.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyy[j] - fz[j] * fga[j] * pa_z[j] * pb_yyyz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_yyyz[j] + 14.0 * pa_xxz[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_xxz_yyzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_xxz[j] * fx[j] * fx[j] + 0.5 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_yyz[j] + 

                                0.5 * pa_xxz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xxz[j] * fx[j] * pb_zz[j] + pa_xx[j] * fx[j] * pb_yyz[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_yyzz[j] + pa_xxz[j] * pb_yyzz[j]) * s_0_0[j] + (-0.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.25 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - pa_xxz[j] * fx[j] * fz[j] * fgb[j] - pa_xx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                2.5 * pa_xxz[j] * fz[j] * fx[j] * fx[j] + 5.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 

                                0.5 * fx[j] * pa_z[j] * pb_yy[j] * fz[j] * fgb[j] - 0.5 * fx[j] * pa_z[j] * fz[j] * fgb[j] * pb_zz[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * pb_yy[j] * fx[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_zz[j] - 

                                fz[j] * fga[j] * fx[j] * pb_yyz[j] - pa_xxz[j] * pb_yy[j] * fz[j] * fgb[j] - pa_xxz[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 6.0 * pa_xxz[j] * fz[j] * pb_yy[j] * fx[j] + 

                                6.0 * pa_xxz[j] * fz[j] * fx[j] * pb_zz[j] + 12.0 * pa_xx[j] * fz[j] * fx[j] * pb_yyz[j] - fz[j] * fga[j] * pa_z[j] * pb_yyzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_yyzz[j] + 14.0 * pa_xxz[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_xxz_yzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_xx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_xxz[j] * pb_yz[j] * fx[j] + 

                                1.5 * pa_xx[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_z[j] * pb_yzzz[j] + pa_xxz[j] * pb_yzzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - 1.5 * pa_xx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_y[j] - 

                                1.5 * fx[j] * pa_z[j] * pb_yz[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_z[j] * pb_yz[j] * fx[j] - 

                                1.5 * fz[j] * fga[j] * fx[j] * pb_yzz[j] - 3.0 * pa_xxz[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 

                                18.0 * pa_xxz[j] * fz[j] * pb_yz[j] * fx[j] + 18.0 * pa_xx[j] * fz[j] * fx[j] * pb_yzz[j] - fz[j] * fga[j] * pa_z[j] * pb_yzzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_yzzz[j] + 14.0 * pa_xxz[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_xxz_zzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_xxz[j] * fx[j] * fx[j] + 3.0 * pa_xx[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 

                                fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_xxz[j] * pb_zz[j] * fx[j] + 2.0 * pa_xx[j] * fx[j] * pb_zzz[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_zzzz[j] + pa_xxz[j] * pb_zzzz[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                3.0 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 0.75 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                3.0 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - 3.0 * pa_xxz[j] * fx[j] * fz[j] * fgb[j] - 

                                6.0 * pa_xx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 

                                12.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_xxz[j] * fz[j] * fx[j] * fx[j] + 

                                30.0 * pa_xx[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 3.0 * fx[j] * pa_z[j] * pb_zz[j] * fz[j] * fgb[j] - 

                                3.0 * fz[j] * fga[j] * pa_z[j] * pb_zz[j] * fx[j] - 2.0 * fz[j] * fga[j] * fx[j] * pb_zzz[j] - 

                                6.0 * pa_xxz[j] * pb_zz[j] * fz[j] * fgb[j] + 15.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] + 

                                10.0 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 36.0 * pa_xxz[j] * fz[j] * pb_zz[j] * fx[j] + 

                                24.0 * pa_xx[j] * fz[j] * fx[j] * pb_zzz[j] - fz[j] * fga[j] * pa_z[j] * pb_zzzz[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_zzzz[j] + 

                                14.0 * pa_xxz[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xyy_xxxx, t_xyy_xxxy, t_xyy_xxxz, t_xyy_xxyy, \
                                     t_xyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xxxx[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xyy[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xyy[j] * pb_xx[j] * fx[j] + 2.0 * fx[j] * pa_yy[j] * pb_xxx[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xxxx[j] + pa_xyy[j] * pb_xxxx[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                3.0 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 3.0 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * fx[j] - 6.0 * fx[j] * pa_yy[j] * pb_x[j] * fz[j] * fgb[j] + 

                                12.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] + 

                                30.0 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_x[j] - 3.0 * pa_x[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 2.0 * fx[j] * fz[j] * fga[j] * pb_xxx[j] - 

                                6.0 * pa_xyy[j] * pb_xx[j] * fz[j] * fgb[j] + 15.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 

                                10.0 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 36.0 * pa_xyy[j] * fz[j] * pb_xx[j] * fx[j] + 

                                24.0 * fx[j] * pa_yy[j] * fz[j] * pb_xxx[j] - pa_x[j] * fz[j] * fga[j] * pb_xxxx[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxxx[j] + 

                                14.0 * pa_xyy[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_xyy_xxxy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 

                                1.5 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] + 

                                1.5 * pa_xyy[j] * pb_xy[j] * fx[j] + pa_xy[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_yy[j] * pb_xxy[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xxxy[j] + pa_xyy[j] * pb_xxxy[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] - 0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 3.0 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * pa_yy[j] * fz[j] * fgb[j] * pb_y[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                15.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_y[j] + 

                                15.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xx[j] - 1.5 * pa_x[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 1.5 * fx[j] * fz[j] * fga[j] * pb_xxy[j] - 

                                3.0 * pa_xyy[j] * pb_xy[j] * fz[j] * fgb[j] + 7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 18.0 * pa_xyy[j] * fz[j] * pb_xy[j] * fx[j] + 

                                12.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxx[j] + 18.0 * fx[j] * pa_yy[j] * fz[j] * pb_xxy[j] - pa_x[j] * fz[j] * fga[j] * pb_xxxy[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxxy[j] + 14.0 * pa_xyy[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_xyy_xxxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] + 1.5 * pa_xyy[j] * pb_xz[j] * fx[j] + 

                                1.5 * fx[j] * pa_yy[j] * pb_xxz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxxz[j] + pa_xyy[j] * pb_xxxz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 1.5 * fx[j] * pa_yy[j] * fz[j] * fgb[j] * pb_z[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_xxz[j] - 3.0 * pa_xyy[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 

                                18.0 * pa_xyy[j] * fz[j] * pb_xz[j] * fx[j] + 18.0 * fx[j] * pa_yy[j] * fz[j] * pb_xxz[j] - pa_x[j] * fz[j] * fga[j] * pb_xxxz[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxxz[j] + 14.0 * pa_xyy[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_xyy_xxyy[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xyy[j] * fx[j] * fx[j] + pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                0.5 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 2.0 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_xyy[j] * pb_xx[j] * fx[j] + 

                                0.5 * pa_xyy[j] * fx[j] * pb_yy[j] + 2.0 * pa_xy[j] * fx[j] * pb_xxy[j] + fx[j] * pa_yy[j] * pb_xyy[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xxyy[j] + pa_xyy[j] * pb_xxyy[j]) * s_0_0[j] + (-pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - pa_xyy[j] * fx[j] * fz[j] * fgb[j] - 

                                2.0 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - fx[j] * pa_yy[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] + 10.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 5.0 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_x[j] + 

                                20.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xy[j] - 0.5 * pa_x[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                0.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_yy[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - fx[j] * fz[j] * fga[j] * pb_xyy[j] - pa_xyy[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                pa_xyy[j] * fz[j] * fgb[j] * pb_yy[j] + 2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 6.0 * pa_xyy[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_xyy[j] * fz[j] * fx[j] * pb_yy[j] + 24.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxy[j] + 

                                12.0 * fx[j] * pa_yy[j] * fz[j] * pb_xyy[j] - pa_x[j] * fz[j] * fga[j] * pb_xxyy[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxyy[j] + 

                                14.0 * pa_xyy[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_xyy_xxyz[j] = (0.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xyy[j] * fx[j] * pb_yz[j] + 

                                pa_xy[j] * fx[j] * pb_xxz[j] + fx[j] * pa_yy[j] * pb_xyz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxyz[j] + 

                                pa_xyy[j] * pb_xxyz[j]) * s_0_0[j] + (-pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                5.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 10.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xz[j] - 

                                0.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_yz[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - 

                                fx[j] * fz[j] * fga[j] * pb_xyz[j] - pa_xyy[j] * fz[j] * fgb[j] * pb_yz[j] + 

                                2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 5.0 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * pa_xyy[j] * fz[j] * fx[j] * pb_yz[j] + 12.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxz[j] + 

                                12.0 * fx[j] * pa_yy[j] * fz[j] * pb_xyz[j] - pa_x[j] * fz[j] * fga[j] * pb_xxyz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxyz[j] + 

                                14.0 * pa_xyy[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyy_xxzz, t_xyy_xyyy, \
                                     t_xyy_xyyz, t_xyy_xyzz, t_xyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_xxzz[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xyy[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_yy[j] * pb_x[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_xzz[j] + 

                                0.5 * pa_xyy[j] * pb_xx[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * pb_zz[j] + fx[j] * pa_yy[j] * pb_xzz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xxzz[j] + pa_xyy[j] * pb_xxzz[j]) * s_0_0[j] + (-0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - pa_xyy[j] * fx[j] * fz[j] * fgb[j] + pa_x[j] * fz[j] * fx[j] * fx[j] * fx[j] - 

                                fx[j] * pa_yy[j] * pb_x[j] * fz[j] * fgb[j] + 2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] + 5.0 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_x[j] - 

                                0.5 * pa_x[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 0.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 

                                0.5 * pa_x[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - 

                                fx[j] * fz[j] * fga[j] * pb_xzz[j] - pa_xyy[j] * pb_xx[j] * fz[j] * fgb[j] - pa_xyy[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 6.0 * pa_xyy[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_xyy[j] * fz[j] * fx[j] * pb_zz[j] + 12.0 * fx[j] * pa_yy[j] * fz[j] * pb_xzz[j] - pa_x[j] * fz[j] * fga[j] * pb_xxzz[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxzz[j] + 14.0 * pa_xyy[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_xyy_xyyy[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_y[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                1.5 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 

                                0.75 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 1.5 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 

                                1.5 * pa_xyy[j] * pb_xy[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_yy[j] * pb_yyy[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xyyy[j] + pa_xyy[j] * pb_xyyy[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                3.0 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * fx[j] * pa_yy[j] * pb_y[j] * fz[j] * fgb[j] + 

                                15.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                7.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_y[j] + 15.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yy[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pb_yyy[j] - 3.0 * pa_xyy[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 18.0 * pa_xyy[j] * fz[j] * pb_xy[j] * fx[j] + 

                                36.0 * pa_xy[j] * fz[j] * fx[j] * pb_xyy[j] + 6.0 * fx[j] * pa_yy[j] * fz[j] * pb_yyy[j] - pa_x[j] * fz[j] * fga[j] * pb_xyyy[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xyyy[j] + 14.0 * pa_xyy[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_xyy_xyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 

                                0.25 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] + 

                                0.5 * pa_xyy[j] * pb_xz[j] * fx[j] + 2.0 * pa_xy[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yy[j] * pb_yyz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xyyz[j] + pa_xyy[j] * pb_xyyz[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                0.5 * fx[j] * pa_yy[j] * fz[j] * fgb[j] * pb_z[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                2.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_z[j] + 10.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yz[j] - 

                                0.5 * pa_x[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pb_yyz[j] - pa_xyy[j] * pb_xz[j] * fz[j] * fgb[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 

                                6.0 * pa_xyy[j] * fz[j] * pb_xz[j] * fx[j] + 24.0 * pa_xy[j] * fz[j] * fx[j] * pb_xyz[j] + 

                                6.0 * fx[j] * pa_yy[j] * fz[j] * pb_yyz[j] - pa_x[j] * fz[j] * fga[j] * pb_xyyz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xyyz[j] + 

                                14.0 * pa_xyy[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_xyy_xyzz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.5 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_yy[j] * pb_y[j] + 

                                0.5 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_yzz[j] + 

                                0.5 * pa_xyy[j] * pb_xy[j] * fx[j] + pa_xy[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yy[j] * pb_yzz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xyzz[j] + pa_xyy[j] * pb_xyzz[j]) * s_0_0[j] + (-0.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] + 

                                2.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] - 0.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * pa_yy[j] * pb_y[j] * fz[j] * fgb[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                5.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 2.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_y[j] + 

                                5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_zz[j] - 0.5 * pa_x[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                0.5 * pa_x[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 0.5 * fx[j] * fz[j] * fga[j] * pb_yzz[j] - pa_xyy[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 

                                6.0 * pa_xyy[j] * fz[j] * pb_xy[j] * fx[j] + 12.0 * pa_xy[j] * fz[j] * fx[j] * pb_xzz[j] + 

                                6.0 * fx[j] * pa_yy[j] * fz[j] * pb_yzz[j] - pa_x[j] * fz[j] * fga[j] * pb_xyzz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xyzz[j] + 

                                14.0 * pa_xyy[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_xyy_xzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_yy[j] * pb_z[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_xyy[j] * pb_xz[j] * fx[j] + 

                                0.5 * fx[j] * pa_yy[j] * pb_zzz[j] + 0.5 * pa_x[j] * fx[j] * pb_xzzz[j] + pa_xyy[j] * pb_xzzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 1.5 * fx[j] * pa_yy[j] * pb_z[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * fx[j] * fx[j] * pa_yy[j] * fz[j] * pb_z[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pb_zzz[j] - 3.0 * pa_xyy[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 

                                18.0 * pa_xyy[j] * fz[j] * pb_xz[j] * fx[j] + 6.0 * fx[j] * pa_yy[j] * fz[j] * pb_zzz[j] - pa_x[j] * fz[j] * fga[j] * pb_xzzz[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xzzz[j] + 14.0 * pa_xyy[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xy, pa_xyy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xyy_yyyy, \
                                     t_xyy_yyyz, t_xyy_yyzz, t_xyy_yzzz, t_xyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyy_yyyy[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] + 

                                6.0 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xyy[j] * pb_yy[j] * fx[j] + 

                                4.0 * pa_xy[j] * fx[j] * pb_yyy[j] + 0.5 * pa_x[j] * fx[j] * pb_yyyy[j] + pa_xyy[j] * pb_yyyy[j]) * s_0_0[j] + (-4.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                3.0 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] - 12.0 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] + 60.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                45.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] - 3.0 * pa_x[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 6.0 * pa_xyy[j] * pb_yy[j] * fz[j] * fgb[j] + 

                                36.0 * pa_xyy[j] * fz[j] * pb_yy[j] * fx[j] + 48.0 * pa_xy[j] * fz[j] * fx[j] * pb_yyy[j] - pa_x[j] * fz[j] * fga[j] * pb_yyyy[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_yyyy[j] + 14.0 * pa_xyy[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_xyy_yyyz[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xyy[j] * pb_yz[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * pb_yyz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_yyyz[j] + pa_xyy[j] * pb_yyyz[j]) * s_0_0[j] + (-3.0 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                15.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 

                                3.0 * pa_xyy[j] * pb_yz[j] * fz[j] * fgb[j] + 18.0 * pa_xyy[j] * fz[j] * pb_yz[j] * fx[j] + 

                                36.0 * pa_xy[j] * fz[j] * fx[j] * pb_yyz[j] - pa_x[j] * fz[j] * fga[j] * pb_yyyz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_yyyz[j] + 

                                14.0 * pa_xyy[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_xyy_yyzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xyy[j] * fx[j] * fx[j] + 

                                pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_xyy[j] * pb_yy[j] * fx[j] + 0.5 * pa_xyy[j] * fx[j] * pb_zz[j] + 

                                2.0 * pa_xy[j] * fx[j] * pb_yzz[j] + 0.5 * pa_x[j] * fx[j] * pb_yyzz[j] + pa_xyy[j] * pb_yyzz[j]) * s_0_0[j] + (-pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                pa_xyy[j] * fx[j] * fz[j] * fgb[j] - 2.0 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] + 10.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 0.5 * pa_x[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                0.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - pa_xyy[j] * pb_yy[j] * fz[j] * fgb[j] - pa_xyy[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 6.0 * pa_xyy[j] * fz[j] * pb_yy[j] * fx[j] + 

                                6.0 * pa_xyy[j] * fz[j] * fx[j] * pb_zz[j] + 24.0 * pa_xy[j] * fz[j] * fx[j] * pb_yzz[j] - pa_x[j] * fz[j] * fga[j] * pb_yyzz[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_yyzz[j] + 14.0 * pa_xyy[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_xyy_yzzz[j] = (1.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xyy[j] * pb_yz[j] * fx[j] + pa_xy[j] * fx[j] * pb_zzz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_yzzz[j] + pa_xyy[j] * pb_yzzz[j]) * s_0_0[j] + (-3.0 * pa_xy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 

                                15.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 1.5 * pa_x[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 3.0 * pa_xyy[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 18.0 * pa_xyy[j] * fz[j] * pb_yz[j] * fx[j] + 

                                12.0 * pa_xy[j] * fz[j] * fx[j] * pb_zzz[j] - pa_x[j] * fz[j] * fga[j] * pb_yzzz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_yzzz[j] + 

                                14.0 * pa_xyy[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_xyy_zzzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xyy[j] * fx[j] * fx[j] + 

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xyy[j] * pb_zz[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_zzzz[j] + 

                                pa_xyy[j] * pb_zzzz[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * pa_xyy[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * fx[j] + 7.5 * pa_xyy[j] * fz[j] * fx[j] * fx[j] - 

                                3.0 * pa_x[j] * fx[j] * pb_zz[j] * fz[j] * fgb[j] - 3.0 * pa_x[j] * fz[j] * fga[j] * pb_zz[j] * fx[j] - 

                                6.0 * pa_xyy[j] * pb_zz[j] * fz[j] * fgb[j] + 15.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                36.0 * pa_xyy[j] * fz[j] * pb_zz[j] * fx[j] - pa_x[j] * fz[j] * fga[j] * pb_zzzz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_zzzz[j] + 

                                14.0 * pa_xyy[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (12) = (60,65)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xyz_xxxx, t_xyz_xxxy, t_xyz_xxxz, \
                                     t_xyz_xxyy, t_xyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xxxx[j] = (0.75 * pa_xyz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + 

                                3.0 * pa_xyz[j] * pb_xx[j] * fx[j] + 2.0 * fx[j] * pa_yz[j] * pb_xxx[j] + pa_xyz[j] * pb_xxxx[j]) * s_0_0[j] + (-3.0 * pa_xyz[j] * fx[j] * fz[j] * fgb[j] - 

                                6.0 * fx[j] * pa_yz[j] * pb_x[j] * fz[j] * fgb[j] + 7.5 * pa_xyz[j] * fz[j] * fx[j] * fx[j] + 

                                30.0 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_x[j] - 6.0 * pa_xyz[j] * pb_xx[j] * fz[j] * fgb[j] + 

                                36.0 * pa_xyz[j] * fz[j] * pb_xx[j] * fx[j] + 24.0 * fx[j] * pa_yz[j] * fz[j] * pb_xxx[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_xyz_xxxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 1.5 * pa_xyz[j] * pb_xy[j] * fx[j] + 

                                0.5 * pa_xz[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_yz[j] * pb_xxy[j] + pa_xyz[j] * pb_xxxy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 1.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * pa_yz[j] * fz[j] * fgb[j] * pb_y[j] + 7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] - 

                                3.0 * pa_xyz[j] * pb_xy[j] * fz[j] * fgb[j] + 18.0 * pa_xyz[j] * fz[j] * pb_xy[j] * fx[j] + 

                                6.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxx[j] + 18.0 * fx[j] * pa_yz[j] * fz[j] * pb_xxy[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_xyz_xxxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_xx[j] + 1.5 * pa_xyz[j] * pb_xz[j] * fx[j] + 

                                0.5 * pa_xy[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_yz[j] * pb_xxz[j] + pa_xyz[j] * pb_xxxz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] - 1.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * pa_yz[j] * fz[j] * fgb[j] * pb_z[j] + 7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                7.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_z[j] + 7.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xx[j] - 

                                3.0 * pa_xyz[j] * pb_xz[j] * fz[j] * fgb[j] + 18.0 * pa_xyz[j] * fz[j] * pb_xz[j] * fx[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxx[j] + 18.0 * fx[j] * pa_yz[j] * fz[j] * pb_xxz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_xyz_xxyy[j] = (0.25 * pa_xyz[j] * fx[j] * fx[j] + 0.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.5 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * pa_xyz[j] * pb_xx[j] * fx[j] + 

                                0.5 * pa_xyz[j] * fx[j] * pb_yy[j] + pa_xz[j] * fx[j] * pb_xxy[j] + fx[j] * pa_yz[j] * pb_xyy[j] + pa_xyz[j] * pb_xxyy[j]) * s_0_0[j] + (-pa_xyz[j] * fx[j] * fz[j] * fgb[j] - pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                fx[j] * pa_yz[j] * pb_x[j] * fz[j] * fgb[j] + 2.5 * pa_xyz[j] * fz[j] * fx[j] * fx[j] + 

                                5.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 5.0 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_x[j] + 

                                10.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] - pa_xyz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_xyz[j] * fz[j] * fgb[j] * pb_yy[j] + 

                                6.0 * pa_xyz[j] * fz[j] * pb_xx[j] * fx[j] + 6.0 * pa_xyz[j] * fz[j] * fx[j] * pb_yy[j] + 

                                12.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxy[j] + 12.0 * fx[j] * pa_yz[j] * fz[j] * pb_xyy[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_xyz_xxyz[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_xy[j] + 0.5 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 

                                0.5 * pa_xyz[j] * fx[j] * pb_yz[j] + 0.5 * pa_xy[j] * fx[j] * pb_xxy[j] + 0.5 * pa_xz[j] * fx[j] * pb_xxz[j] + 

                                fx[j] * pa_yz[j] * pb_xyz[j] + pa_xyz[j] * pb_xxyz[j]) * s_0_0[j] + (-0.25 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.5 * pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 0.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 2.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                2.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xy[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] - pa_xyz[j] * fz[j] * fgb[j] * pb_yz[j] + 6.0 * pa_xyz[j] * fz[j] * fx[j] * pb_yz[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxy[j] + 6.0 * pa_xz[j] * fx[j] * fz[j] * pb_xxz[j] + 

                                12.0 * fx[j] * pa_yz[j] * fz[j] * pb_xyz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (13) = (65,70)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xyz_xxzz, \
                                     t_xyz_xyyy, t_xyz_xyyz, t_xyz_xyzz, t_xyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_xxzz[j] = (0.25 * pa_xyz[j] * fx[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * fx[j] * fx[j] * pa_yz[j] * pb_x[j] + fx[j] * fx[j] * pa_y[j] * pb_xz[j] + 0.5 * pa_xyz[j] * pb_xx[j] * fx[j] + 

                                0.5 * pa_xyz[j] * fx[j] * pb_zz[j] + pa_xy[j] * fx[j] * pb_xxz[j] + fx[j] * pa_yz[j] * pb_xzz[j] + pa_xyz[j] * pb_xxzz[j]) * s_0_0[j] + (-pa_xyz[j] * fx[j] * fz[j] * fgb[j] - pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                fx[j] * pa_yz[j] * pb_x[j] * fz[j] * fgb[j] + 2.5 * pa_xyz[j] * fz[j] * fx[j] * fx[j] + 

                                5.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 5.0 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_x[j] + 

                                10.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_xz[j] - pa_xyz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_xyz[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                6.0 * pa_xyz[j] * fz[j] * pb_xx[j] * fx[j] + 6.0 * pa_xyz[j] * fz[j] * fx[j] * pb_zz[j] + 

                                12.0 * pa_xy[j] * fz[j] * fx[j] * pb_xxz[j] + 12.0 * fx[j] * pa_yz[j] * fz[j] * pb_xzz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_xyz_xyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 1.5 * pa_xyz[j] * pb_xy[j] * fx[j] + 

                                1.5 * pa_xz[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyy[j] + pa_xyz[j] * pb_xyyy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 1.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * pa_yz[j] * pb_y[j] * fz[j] * fgb[j] + 7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                7.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] - 

                                3.0 * pa_xyz[j] * pb_xy[j] * fz[j] * fgb[j] + 18.0 * pa_xyz[j] * fz[j] * pb_xy[j] * fx[j] + 

                                18.0 * pa_xz[j] * fx[j] * fz[j] * pb_xyy[j] + 6.0 * fx[j] * pa_yz[j] * fz[j] * pb_yyy[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_xyz_xyyz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.25 * fx[j] * fx[j] * pa_y[j] * pb_yy[j] + 

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.5 * pa_xyz[j] * pb_xz[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * pb_xyy[j] + 

                                pa_xz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yyz[j] + pa_xyz[j] * pb_xyyz[j]) * s_0_0[j] + (-0.25 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] + 

                                fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] + 2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.5 * fx[j] * pa_yz[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                2.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_z[j] + 2.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yy[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] - pa_xyz[j] * pb_xz[j] * fz[j] * fgb[j] + 6.0 * pa_xyz[j] * fz[j] * pb_xz[j] * fx[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_xyy[j] + 12.0 * pa_xz[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_yyz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_xyz_xyzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.5 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 

                                0.25 * fx[j] * fx[j] * pa_yz[j] * pb_y[j] + 0.5 * fx[j] * fx[j] * pa_y[j] * pb_yz[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * pa_xyz[j] * pb_xy[j] * fx[j] + pa_xy[j] * fx[j] * pb_xyz[j] + 

                                0.5 * pa_xz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yz[j] * pb_yzz[j] + pa_xyz[j] * pb_xyzz[j]) * s_0_0[j] + (-0.25 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.5 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.5 * fx[j] * pa_yz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 5.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                2.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_y[j] + 5.0 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_yz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] - pa_xyz[j] * pb_xy[j] * fz[j] * fgb[j] + 6.0 * pa_xyz[j] * fz[j] * pb_xy[j] * fx[j] + 

                                12.0 * pa_xy[j] * fz[j] * fx[j] * pb_xyz[j] + 6.0 * pa_xz[j] * fx[j] * fz[j] * pb_xzz[j] + 

                                6.0 * fx[j] * pa_yz[j] * fz[j] * pb_yzz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_xyz_xzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_y[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_yz[j] * pb_z[j] + 0.75 * fx[j] * fx[j] * pa_y[j] * pb_zz[j] + 1.5 * pa_xyz[j] * pb_xz[j] * fx[j] + 

                                1.5 * pa_xy[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_yz[j] * pb_zzz[j] + pa_xyz[j] * pb_xzzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pa_y[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * pa_y[j] * fz[j] - 1.5 * pa_xy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * pa_yz[j] * pb_z[j] * fz[j] * fgb[j] + 7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                7.5 * fx[j] * fx[j] * pa_yz[j] * fz[j] * pb_z[j] + 7.5 * fx[j] * fx[j] * pa_y[j] * fz[j] * pb_zz[j] - 

                                3.0 * pa_xyz[j] * pb_xz[j] * fz[j] * fgb[j] + 18.0 * pa_xyz[j] * fz[j] * pb_xz[j] * fx[j] + 

                                18.0 * pa_xy[j] * fz[j] * fx[j] * pb_xzz[j] + 6.0 * fx[j] * pa_yz[j] * fz[j] * pb_zzz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (14) = (70,75)

            #pragma omp simd aligned(fgb, fx, fz, pa_x, pa_xy, pa_xyz, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, \
                                     t_xyz_yyyy, t_xyz_yyyz, t_xyz_yyzz, t_xyz_yzzz, t_xyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyz_yyyy[j] = (0.75 * pa_xyz[j] * fx[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                3.0 * pa_xyz[j] * pb_yy[j] * fx[j] + 2.0 * pa_xz[j] * fx[j] * pb_yyy[j] + pa_xyz[j] * pb_yyyy[j]) * s_0_0[j] + (-3.0 * pa_xyz[j] * fx[j] * fz[j] * fgb[j] - 

                                6.0 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 7.5 * pa_xyz[j] * fz[j] * fx[j] * fx[j] + 

                                30.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 6.0 * pa_xyz[j] * pb_yy[j] * fz[j] * fgb[j] + 

                                36.0 * pa_xyz[j] * fz[j] * pb_yy[j] * fx[j] + 24.0 * pa_xz[j] * fx[j] * fz[j] * pb_yyy[j] + 14.0 * pa_xyz[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_xyz_yyyz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * pa_xyz[j] * pb_yz[j] * fx[j] + 

                                0.5 * pa_xy[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xz[j] * fx[j] * pb_yyz[j] + pa_xyz[j] * pb_yyyz[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 1.5 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] - 

                                3.0 * pa_xyz[j] * pb_yz[j] * fz[j] * fgb[j] + 18.0 * pa_xyz[j] * fz[j] * pb_yz[j] * fx[j] + 

                                6.0 * pa_xy[j] * fz[j] * fx[j] * pb_yyy[j] + 18.0 * pa_xz[j] * fx[j] * fz[j] * pb_yyz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_xyz_yyzz[j] = (0.25 * pa_xyz[j] * fx[j] * fx[j] + 0.5 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * pa_xyz[j] * pb_yy[j] * fx[j] + 

                                0.5 * pa_xyz[j] * fx[j] * pb_zz[j] + pa_xy[j] * fx[j] * pb_yyz[j] + pa_xz[j] * fx[j] * pb_yzz[j] + pa_xyz[j] * pb_yyzz[j]) * s_0_0[j] + (-pa_xyz[j] * fx[j] * fz[j] * fgb[j] - pa_xy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 2.5 * pa_xyz[j] * fz[j] * fx[j] * fx[j] + 

                                5.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 5.0 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                10.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] - pa_xyz[j] * pb_yy[j] * fz[j] * fgb[j] - pa_xyz[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                6.0 * pa_xyz[j] * fz[j] * pb_yy[j] * fx[j] + 6.0 * pa_xyz[j] * fz[j] * fx[j] * pb_zz[j] + 

                                12.0 * pa_xy[j] * fz[j] * fx[j] * pb_yyz[j] + 12.0 * pa_xz[j] * fx[j] * fz[j] * pb_yzz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_xyz_yzzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xy[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * pa_xyz[j] * pb_yz[j] * fx[j] + 

                                1.5 * pa_xy[j] * fx[j] * pb_yzz[j] + 0.5 * pa_xz[j] * fx[j] * pb_zzz[j] + pa_xyz[j] * pb_yzzz[j]) * s_0_0[j] + (-0.75 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 1.5 * pa_xy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * pa_xz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 7.5 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                7.5 * pa_xz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 

                                3.0 * pa_xyz[j] * pb_yz[j] * fz[j] * fgb[j] + 18.0 * pa_xyz[j] * fz[j] * pb_yz[j] * fx[j] + 

                                18.0 * pa_xy[j] * fz[j] * fx[j] * pb_yzz[j] + 6.0 * pa_xz[j] * fx[j] * fz[j] * pb_zzz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_xyz_zzzz[j] = (0.75 * pa_xyz[j] * fx[j] * fx[j] + 3.0 * pa_xy[j] * fx[j] * fx[j] * pb_z[j] + 

                                3.0 * pa_xyz[j] * pb_zz[j] * fx[j] + 2.0 * pa_xy[j] * fx[j] * pb_zzz[j] + pa_xyz[j] * pb_zzzz[j]) * s_0_0[j] + (-3.0 * pa_xyz[j] * fx[j] * fz[j] * fgb[j] - 

                                6.0 * pa_xy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 7.5 * pa_xyz[j] * fz[j] * fx[j] * fx[j] + 

                                30.0 * pa_xy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 6.0 * pa_xyz[j] * pb_zz[j] * fz[j] * fgb[j] + 

                                36.0 * pa_xyz[j] * fz[j] * pb_zz[j] * fx[j] + 24.0 * pa_xy[j] * fz[j] * fx[j] * pb_zzz[j] + 14.0 * pa_xyz[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (15) = (75,80)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, \
                                     pb_yy, pb_yz, pb_z, r_0_0, s_0_0, t_xzz_xxxx, t_xzz_xxxy, t_xzz_xxxz, t_xzz_xxyy, \
                                     t_xzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xxxx[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_xzz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 1.5 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                fx[j] * fx[j] * pb_xxx[j] + 3.0 * pa_xzz[j] * pb_xx[j] * fx[j] + 2.0 * fx[j] * pa_zz[j] * pb_xxx[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xxxx[j] + pa_xzz[j] * pb_xxxx[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                3.0 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 3.0 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * fx[j] - 6.0 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                12.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_xzz[j] * fz[j] * fx[j] * fx[j] + 

                                30.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] - 3.0 * pa_x[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 2.0 * fx[j] * fz[j] * fga[j] * pb_xxx[j] - 

                                6.0 * pa_xzz[j] * pb_xx[j] * fz[j] * fgb[j] + 15.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 

                                10.0 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 36.0 * pa_xzz[j] * fz[j] * pb_xx[j] * fx[j] + 

                                24.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxx[j] - pa_x[j] * fz[j] * fga[j] * pb_xxxx[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxxx[j] + 

                                14.0 * pa_xzz[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_xzz_xxxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] + 1.5 * pa_xzz[j] * pb_xy[j] * fx[j] + 

                                1.5 * fx[j] * pa_zz[j] * pb_xxy[j] + 0.5 * pa_x[j] * fx[j] * pb_xxxy[j] + pa_xzz[j] * pb_xxxy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 1.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_y[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_xxy[j] - 3.0 * pa_xzz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 

                                18.0 * pa_xzz[j] * fz[j] * pb_xy[j] * fx[j] + 18.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxy[j] - pa_x[j] * fz[j] * fga[j] * pb_xxxy[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxxy[j] + 14.0 * pa_xzz[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_xzz_xxxz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] + 

                                1.5 * pa_xzz[j] * pb_xz[j] * fx[j] + pa_xz[j] * fx[j] * pb_xxx[j] + 1.5 * fx[j] * pa_zz[j] * pb_xxz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xxxz[j] + pa_xzz[j] * pb_xxxz[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] - 0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 3.0 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_z[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                15.0 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 

                                15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xx[j] - 1.5 * pa_x[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 1.5 * fx[j] * fz[j] * fga[j] * pb_xxz[j] - 

                                3.0 * pa_xzz[j] * pb_xz[j] * fz[j] * fgb[j] + 7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 18.0 * pa_xzz[j] * fz[j] * pb_xz[j] * fx[j] + 

                                12.0 * pa_xz[j] * fz[j] * fx[j] * pb_xxx[j] + 18.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxz[j] - pa_x[j] * fz[j] * fga[j] * pb_xxxz[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxxz[j] + 14.0 * pa_xzz[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_xzz_xxyy[j] = (0.125 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xzz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * fx[j] * pb_xyy[j] + 

                                0.5 * pa_xzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] * pb_yy[j] + fx[j] * pa_zz[j] * pb_xyy[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xxyy[j] + pa_xzz[j] * pb_xxyy[j]) * s_0_0[j] + (-0.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - pa_xzz[j] * fx[j] * fz[j] * fgb[j] + pa_x[j] * fz[j] * fx[j] * fx[j] * fx[j] - 

                                fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] + 2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 

                                2.5 * pa_xzz[j] * fz[j] * fx[j] * fx[j] + 5.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] - 

                                0.5 * pa_x[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 0.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_yy[j] - 

                                0.5 * pa_x[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - 

                                fx[j] * fz[j] * fga[j] * pb_xyy[j] - pa_xzz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_xzz[j] * fz[j] * fgb[j] * pb_yy[j] + 

                                2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 6.0 * pa_xzz[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_xzz[j] * fz[j] * fx[j] * pb_yy[j] + 12.0 * fx[j] * pa_zz[j] * fz[j] * pb_xyy[j] - pa_x[j] * fz[j] * fga[j] * pb_xxyy[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxyy[j] + 14.0 * pa_xzz[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_xzz_xxyz[j] = (0.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_xzz[j] * fx[j] * pb_yz[j] + 

                                pa_xz[j] * fx[j] * pb_xxy[j] + fx[j] * pa_zz[j] * pb_xyz[j] + 0.5 * pa_x[j] * fx[j] * pb_xxyz[j] + 

                                pa_xzz[j] * pb_xxyz[j]) * s_0_0[j] + (-pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] + 

                                5.0 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 10.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xy[j] - 

                                0.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_yz[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - 

                                fx[j] * fz[j] * fga[j] * pb_xyz[j] - pa_xzz[j] * fz[j] * fgb[j] * pb_yz[j] + 

                                2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 5.0 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * pa_xzz[j] * fz[j] * fx[j] * pb_yz[j] + 12.0 * pa_xz[j] * fz[j] * fx[j] * pb_xxy[j] + 

                                12.0 * fx[j] * pa_zz[j] * fz[j] * pb_xyz[j] - pa_x[j] * fz[j] * fga[j] * pb_xxyz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxyz[j] + 

                                14.0 * pa_xzz[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (16) = (80,85)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, r_0_0, s_0_0, t_xzz_xxzz, \
                                     t_xzz_xyyy, t_xzz_xyyz, t_xzz_xyzz, t_xzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_xxzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.25 * pa_xzz[j] * fx[j] * fx[j] + pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xx[j] + 

                                0.5 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 2.0 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_xzz[j] * pb_xx[j] * fx[j] + 

                                0.5 * pa_xzz[j] * fx[j] * pb_zz[j] + 2.0 * pa_xz[j] * fx[j] * pb_xxz[j] + fx[j] * pa_zz[j] * pb_xzz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xxzz[j] + pa_xzz[j] * pb_xxzz[j]) * s_0_0[j] + (-pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - pa_xzz[j] * fx[j] * fz[j] * fgb[j] - 

                                2.0 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                2.5 * pa_xzz[j] * fz[j] * fx[j] * fx[j] + 10.0 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] + 5.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] + 

                                20.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xz[j] - 0.5 * pa_x[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                0.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - fx[j] * fz[j] * fga[j] * pb_xzz[j] - pa_xzz[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                pa_xzz[j] * fz[j] * fgb[j] * pb_zz[j] + 2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 6.0 * pa_xzz[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_xzz[j] * fz[j] * fx[j] * pb_zz[j] + 24.0 * pa_xz[j] * fz[j] * fx[j] * pb_xxz[j] + 

                                12.0 * fx[j] * pa_zz[j] * fz[j] * pb_xzz[j] - pa_x[j] * fz[j] * fga[j] * pb_xxzz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xxzz[j] + 

                                14.0 * pa_xzz[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_xzz_xyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_xzz[j] * pb_xy[j] * fx[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_yyy[j] + 0.5 * pa_x[j] * fx[j] * pb_xyyy[j] + pa_xzz[j] * pb_xyyy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 1.5 * fx[j] * pa_zz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pb_yyy[j] - 3.0 * pa_xzz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 

                                18.0 * pa_xzz[j] * fz[j] * pb_xy[j] * fx[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_yyy[j] - pa_x[j] * fz[j] * fga[j] * pb_xyyy[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xyyy[j] + 14.0 * pa_xzz[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_xzz_xyyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.25 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_yyz[j] + 

                                0.5 * pa_xzz[j] * pb_xz[j] * fx[j] + pa_xz[j] * fx[j] * pb_xyy[j] + 0.5 * fx[j] * pa_zz[j] * pb_yyz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xyyz[j] + pa_xzz[j] * pb_xyyz[j]) * s_0_0[j] + (-0.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                2.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] - 0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_z[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                5.0 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 

                                5.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yy[j] - 0.5 * pa_x[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                0.5 * pa_x[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 0.5 * fx[j] * fz[j] * fga[j] * pb_yyz[j] - pa_xzz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 

                                6.0 * pa_xzz[j] * fz[j] * pb_xz[j] * fx[j] + 12.0 * pa_xz[j] * fz[j] * fx[j] * pb_xyy[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_yyz[j] - pa_x[j] * fz[j] * fga[j] * pb_xyyz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xyyz[j] + 

                                14.0 * pa_xzz[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_xzz_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_xy[j] + 

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yzz[j] + 

                                0.5 * pa_xzz[j] * pb_xy[j] * fx[j] + 2.0 * pa_xz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_yzz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xyzz[j] + pa_xzz[j] * pb_xyzz[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                0.5 * fx[j] * pa_zz[j] * pb_y[j] * fz[j] * fgb[j] + 7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] + 10.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yz[j] - 

                                0.5 * pa_x[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pb_yzz[j] - pa_xzz[j] * pb_xy[j] * fz[j] * fgb[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 

                                6.0 * pa_xzz[j] * fz[j] * pb_xy[j] * fx[j] + 24.0 * pa_xz[j] * fz[j] * fx[j] * pb_xyz[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_yzz[j] - pa_x[j] * fz[j] * fga[j] * pb_xyzz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_xyzz[j] + 

                                14.0 * pa_xzz[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_xzz_xzzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_xz[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_x[j] * fx[j] * fx[j] * pb_xz[j] + 

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 

                                1.5 * pa_xzz[j] * pb_xz[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_zzz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_xzzz[j] + pa_xzz[j] * pb_xzzz[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.75 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                3.0 * pa_xz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 1.5 * fx[j] * pa_zz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                15.0 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] + 

                                7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_zz[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pb_zzz[j] - 3.0 * pa_xzz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 18.0 * pa_xzz[j] * fz[j] * pb_xz[j] * fx[j] + 

                                36.0 * pa_xz[j] * fz[j] * fx[j] * pb_xzz[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_zzz[j] - pa_x[j] * fz[j] * fga[j] * pb_xzzz[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_xzzz[j] + 14.0 * pa_xzz[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (17) = (85,90)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_x, pa_xz, pa_xzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_xzz_yyyy, \
                                     t_xzz_yyyz, t_xzz_yyzz, t_xzz_yzzz, t_xzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzz_yyyy[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] + 

                                1.5 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_xzz[j] * pb_yy[j] * fx[j] + 0.5 * pa_x[j] * fx[j] * pb_yyyy[j] + 

                                pa_xzz[j] * pb_yyyy[j]) * s_0_0[j] + (-1.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * fx[j] + 7.5 * pa_xzz[j] * fz[j] * fx[j] * fx[j] - 

                                3.0 * pa_x[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 3.0 * pa_x[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                6.0 * pa_xzz[j] * pb_yy[j] * fz[j] * fgb[j] + 15.0 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                36.0 * pa_xzz[j] * fz[j] * pb_yy[j] * fx[j] - pa_x[j] * fz[j] * fga[j] * pb_yyyy[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_yyyy[j] + 

                                14.0 * pa_xzz[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_xzz_yyyz[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xzz[j] * pb_yz[j] * fx[j] + pa_xz[j] * fx[j] * pb_yyy[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_yyyz[j] + pa_xzz[j] * pb_yyyz[j]) * s_0_0[j] + (-3.0 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                15.0 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] - 1.5 * pa_x[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 

                                1.5 * pa_x[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 3.0 * pa_xzz[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                7.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 18.0 * pa_xzz[j] * fz[j] * pb_yz[j] * fx[j] + 

                                12.0 * pa_xz[j] * fz[j] * fx[j] * pb_yyy[j] - pa_x[j] * fz[j] * fga[j] * pb_yyyz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_yyyz[j] + 

                                14.0 * pa_xzz[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_xzz_yyzz[j] = (0.375 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_xzz[j] * fx[j] * fx[j] + 

                                pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_x[j] * fx[j] * fx[j] * pb_yy[j] + 

                                0.25 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_xzz[j] * pb_yy[j] * fx[j] + 0.5 * pa_xzz[j] * fx[j] * pb_zz[j] + 

                                2.0 * pa_xz[j] * fx[j] * pb_yyz[j] + 0.5 * pa_x[j] * fx[j] * pb_yyzz[j] + pa_xzz[j] * pb_yyzz[j]) * s_0_0[j] + (-pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.25 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                pa_xzz[j] * fx[j] * fz[j] * fgb[j] - 2.0 * pa_xz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_xzz[j] * fz[j] * fx[j] * fx[j] + 10.0 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                7.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] - 0.5 * pa_x[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                0.5 * pa_x[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 0.5 * pa_x[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                0.5 * pa_x[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - pa_xzz[j] * pb_yy[j] * fz[j] * fgb[j] - pa_xzz[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * pa_x[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 6.0 * pa_xzz[j] * fz[j] * pb_yy[j] * fx[j] + 

                                6.0 * pa_xzz[j] * fz[j] * fx[j] * pb_zz[j] + 24.0 * pa_xz[j] * fz[j] * fx[j] * pb_yyz[j] - pa_x[j] * fz[j] * fga[j] * pb_yyzz[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_yyzz[j] + 14.0 * pa_xzz[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_xzz_yzzz[j] = (1.5 * pa_xz[j] * fx[j] * fx[j] * pb_y[j] + 

                                2.25 * pa_x[j] * fx[j] * fx[j] * pb_yz[j] + 1.5 * pa_xzz[j] * pb_yz[j] * fx[j] + 3.0 * pa_xz[j] * fx[j] * pb_yzz[j] + 

                                0.5 * pa_x[j] * fx[j] * pb_yzzz[j] + pa_xzz[j] * pb_yzzz[j]) * s_0_0[j] + (-3.0 * pa_xz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                15.0 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 22.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] - 

                                1.5 * pa_x[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 1.5 * pa_x[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 

                                3.0 * pa_xzz[j] * pb_yz[j] * fz[j] * fgb[j] + 18.0 * pa_xzz[j] * fz[j] * pb_yz[j] * fx[j] + 

                                36.0 * pa_xz[j] * fz[j] * fx[j] * pb_yzz[j] - pa_x[j] * fz[j] * fga[j] * pb_yzzz[j] + 6.0 * pa_x[j] * fz[j] * fx[j] * pb_yzzz[j] + 

                                14.0 * pa_xzz[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_xzz_zzzz[j] = (1.875 * pa_x[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_xzz[j] * fx[j] * fx[j] + 

                                6.0 * pa_xz[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_x[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_xzz[j] * pb_zz[j] * fx[j] + 

                                4.0 * pa_xz[j] * fx[j] * pb_zzz[j] + 0.5 * pa_x[j] * fx[j] * pb_zzzz[j] + pa_xzz[j] * pb_zzzz[j]) * s_0_0[j] + (-4.5 * pa_x[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                15.0 * pa_x[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * pa_x[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                3.0 * pa_xzz[j] * fx[j] * fz[j] * fgb[j] - 12.0 * pa_xz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 

                                7.5 * pa_xzz[j] * fz[j] * fx[j] * fx[j] + 60.0 * pa_xz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                45.0 * pa_x[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 3.0 * pa_x[j] * fx[j] * pb_zz[j] * fz[j] * fgb[j] - 

                                3.0 * pa_x[j] * fz[j] * fga[j] * pb_zz[j] * fx[j] - 6.0 * pa_xzz[j] * pb_zz[j] * fz[j] * fgb[j] + 

                                36.0 * pa_xzz[j] * fz[j] * pb_zz[j] * fx[j] + 48.0 * pa_xz[j] * fz[j] * fx[j] * pb_zzz[j] - pa_x[j] * fz[j] * fga[j] * pb_zzzz[j] + 

                                6.0 * pa_x[j] * fz[j] * fx[j] * pb_zzzz[j] + 14.0 * pa_xzz[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (18) = (90,95)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, r_0_0, \
                                     s_0_0, t_yyy_xxxx, t_yyy_xxxy, t_yyy_xxxz, t_yyy_xxyy, t_yyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xxxx[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yyy[j] * fx[j] * fx[j] + 

                                4.5 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_yyy[j] * pb_xx[j] * fx[j] + 1.5 * pa_y[j] * fx[j] * pb_xxxx[j] + 

                                pa_yyy[j] * pb_xxxx[j]) * s_0_0[j] + (-4.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                2.25 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * pa_yyy[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_y[j] * fz[j] * fx[j] * fx[j] * fx[j] + 7.5 * pa_yyy[j] * fz[j] * fx[j] * fx[j] - 

                                9.0 * pa_y[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 9.0 * pa_y[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                6.0 * pa_yyy[j] * pb_xx[j] * fz[j] * fgb[j] + 45.0 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 

                                36.0 * pa_yyy[j] * fz[j] * pb_xx[j] * fx[j] - 3.0 * pa_y[j] * fz[j] * fga[j] * pb_xxxx[j] + 

                                18.0 * pa_y[j] * fz[j] * fx[j] * pb_xxxx[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_yyy_xxxy[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_yyy[j] * pb_xy[j] * fx[j] + 

                                1.5 * pa_yy[j] * fx[j] * pb_xxx[j] + 1.5 * pa_y[j] * fx[j] * pb_xxxy[j] + pa_yyy[j] * pb_xxxy[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 4.5 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 22.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] - 

                                4.5 * pa_y[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 4.5 * pa_y[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_xxx[j] - 3.0 * pa_yyy[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                22.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 

                                18.0 * pa_yyy[j] * fz[j] * pb_xy[j] * fx[j] + 18.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxx[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_xxxy[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_xxxy[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_yyy_xxxz[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yyy[j] * pb_xz[j] * fx[j] + 

                                1.5 * pa_y[j] * fx[j] * pb_xxxz[j] + pa_yyy[j] * pb_xxxz[j]) * s_0_0[j] + (-4.5 * pa_y[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                4.5 * pa_y[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 3.0 * pa_yyy[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                22.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 18.0 * pa_yyy[j] * fz[j] * pb_xz[j] * fx[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_xxxz[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_xxxz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_yyy_xxyy[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_yyy[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * fx[j] * fx[j] * pb_xxy[j] + 

                                0.5 * pa_yyy[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] * pb_yy[j] + 3.0 * pa_yy[j] * fx[j] * pb_xxy[j] + 

                                1.5 * pa_y[j] * fx[j] * pb_xxyy[j] + pa_yyy[j] * pb_xxyy[j]) * s_0_0[j] + (-3.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - pa_yyy[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                2.5 * pa_yyy[j] * fz[j] * fx[j] * fx[j] + 15.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                22.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] - 1.5 * pa_y[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                1.5 * pa_y[j] * fx[j] * fz[j] * fgb[j] * pb_yy[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                1.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xxy[j] - pa_yyy[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                pa_yyy[j] * fz[j] * fgb[j] * pb_yy[j] + 7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 6.0 * pa_yyy[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_yyy[j] * fz[j] * fx[j] * pb_yy[j] + 36.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxy[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_xxyy[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_xxyy[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_yyy_xxyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_xxz[j] + 0.5 * pa_yyy[j] * fx[j] * pb_yz[j] + 

                                1.5 * pa_yy[j] * fx[j] * pb_xxz[j] + 1.5 * pa_y[j] * fx[j] * pb_xxyz[j] + pa_yyy[j] * pb_xxyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 1.5 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 

                                1.5 * pa_y[j] * fx[j] * fz[j] * fgb[j] * pb_yz[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_xxz[j] - pa_yyy[j] * fz[j] * fgb[j] * pb_yz[j] + 

                                7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 

                                6.0 * pa_yyy[j] * fz[j] * fx[j] * pb_yz[j] + 18.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxz[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_xxyz[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_xxyz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (19) = (95,100)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_zz, r_0_0, s_0_0, t_yyy_xxzz, \
                                     t_yyy_xyyy, t_yyy_xyyz, t_yyy_xyzz, t_yyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_xxzz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_yyy[j] * fx[j] * fx[j] + 

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_yyy[j] * pb_xx[j] * fx[j] + 

                                0.5 * pa_yyy[j] * fx[j] * pb_zz[j] + 1.5 * pa_y[j] * fx[j] * pb_xxzz[j] + pa_yyy[j] * pb_xxzz[j]) * s_0_0[j] + (-1.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - pa_yyy[j] * fx[j] * fz[j] * fgb[j] + 3.0 * pa_y[j] * fz[j] * fx[j] * fx[j] * fx[j] + 

                                2.5 * pa_yyy[j] * fz[j] * fx[j] * fx[j] - 1.5 * pa_y[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                1.5 * pa_y[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                1.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - pa_yyy[j] * pb_xx[j] * fz[j] * fgb[j] - pa_yyy[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                6.0 * pa_yyy[j] * fz[j] * pb_xx[j] * fx[j] + 6.0 * pa_yyy[j] * fz[j] * fx[j] * pb_zz[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_xxzz[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_xxzz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_yyy_xyyy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                6.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 2.25 * fx[j] * fx[j] * pb_xyy[j] + 1.5 * pa_yyy[j] * pb_xy[j] * fx[j] + 

                                4.5 * pa_yy[j] * fx[j] * pb_xyy[j] + 1.5 * pa_y[j] * fx[j] * pb_xyyy[j] + pa_yyy[j] * pb_xyyy[j]) * s_0_0[j] + (15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                2.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                4.5 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 22.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                67.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] - 4.5 * pa_y[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                4.5 * pa_y[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 4.5 * fx[j] * fz[j] * fga[j] * pb_xyy[j] - 

                                3.0 * pa_yyy[j] * pb_xy[j] * fz[j] * fgb[j] + 22.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 

                                18.0 * pa_yyy[j] * fz[j] * pb_xy[j] * fx[j] + 54.0 * pa_yy[j] * fz[j] * fx[j] * pb_xyy[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_xyyy[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_xyyy[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_yyy_xyyz[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * fx[j] * fx[j] * pb_xyz[j] + 

                                0.5 * pa_yyy[j] * pb_xz[j] * fx[j] + 3.0 * pa_yy[j] * fx[j] * pb_xyz[j] + 1.5 * pa_y[j] * fx[j] * pb_xyyz[j] + 

                                pa_yyy[j] * pb_xyyz[j]) * s_0_0[j] + (22.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] - 

                                1.5 * pa_y[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xyz[j] - pa_yyy[j] * pb_xz[j] * fz[j] * fgb[j] + 15.0 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * pa_yyy[j] * fz[j] * pb_xz[j] * fx[j] + 36.0 * pa_yy[j] * fz[j] * fx[j] * pb_xyz[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_xyyz[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_xyyz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_yyy_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] + 0.5 * pa_yyy[j] * pb_xy[j] * fx[j] + 

                                1.5 * pa_yy[j] * fx[j] * pb_xzz[j] + 1.5 * pa_y[j] * fx[j] * pb_xyzz[j] + pa_yyy[j] * pb_xyzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 1.5 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] - 

                                1.5 * pa_y[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_xzz[j] - pa_yyy[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 

                                6.0 * pa_yyy[j] * fz[j] * pb_xy[j] * fx[j] + 18.0 * pa_yy[j] * fz[j] * fx[j] * pb_xzz[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_xyzz[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_xyzz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_yyy_xzzz[j] = (2.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yyy[j] * pb_xz[j] * fx[j] + 

                                1.5 * pa_y[j] * fx[j] * pb_xzzz[j] + pa_yyy[j] * pb_xzzz[j]) * s_0_0[j] + (-4.5 * pa_y[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                4.5 * pa_y[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 3.0 * pa_yyy[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                22.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 18.0 * pa_yyy[j] * fz[j] * pb_xz[j] * fx[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_xzzz[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_xzzz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (20) = (100,105)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_yyy_yyyy, \
                                     t_yyy_yyyz, t_yyy_yyzz, t_yyy_yzzz, t_yyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyy_yyyy[j] = (5.625 * pa_y[j] * fx[j] * fx[j] * fx[j] + 7.5 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_yyy[j] * fx[j] * fx[j] + 9.0 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 

                                13.5 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yyy[j] * pb_yy[j] * fx[j] + 

                                6.0 * pa_yy[j] * fx[j] * pb_yyy[j] + 1.5 * pa_y[j] * fx[j] * pb_yyyy[j] + pa_yyy[j] * pb_yyyy[j]) * s_0_0[j] + (-13.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                45.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 60.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                2.25 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 9.0 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                9.0 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 3.0 * pa_yyy[j] * fx[j] * fz[j] * fgb[j] - 

                                18.0 * pa_yy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 7.5 * pa_yyy[j] * fz[j] * fx[j] * fx[j] + 

                                90.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 135.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] - 

                                9.0 * pa_y[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 9.0 * pa_y[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                6.0 * fx[j] * fz[j] * fga[j] * pb_yyy[j] - 6.0 * pa_yyy[j] * pb_yy[j] * fz[j] * fgb[j] + 

                                30.0 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 36.0 * pa_yyy[j] * fz[j] * pb_yy[j] * fx[j] + 

                                72.0 * pa_yy[j] * fz[j] * fx[j] * pb_yyy[j] - 3.0 * pa_y[j] * fz[j] * fga[j] * pb_yyyy[j] + 

                                18.0 * pa_y[j] * fz[j] * fx[j] * pb_yyyy[j] + 14.0 * pa_yyy[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_yyy_yyyz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                6.75 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 2.25 * fx[j] * fx[j] * pb_yyz[j] + 1.5 * pa_yyy[j] * pb_yz[j] * fx[j] + 

                                4.5 * pa_yy[j] * fx[j] * pb_yyz[j] + 1.5 * pa_y[j] * fx[j] * pb_yyyz[j] + pa_yyy[j] * pb_yyyz[j]) * s_0_0[j] + (15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                2.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                4.5 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 22.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                67.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] - 4.5 * pa_y[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 

                                4.5 * pa_y[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 4.5 * fx[j] * fz[j] * fga[j] * pb_yyz[j] - 

                                3.0 * pa_yyy[j] * pb_yz[j] * fz[j] * fgb[j] + 22.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 

                                18.0 * pa_yyy[j] * fz[j] * pb_yz[j] * fx[j] + 54.0 * pa_yy[j] * fz[j] * fx[j] * pb_yyz[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_yyyz[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_yyyz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_yyy_yyzz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_yyy[j] * fx[j] * fx[j] + 1.5 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 1.5 * fx[j] * fx[j] * pb_yzz[j] + 

                                0.5 * pa_yyy[j] * pb_yy[j] * fx[j] + 0.5 * pa_yyy[j] * fx[j] * pb_zz[j] + 3.0 * pa_yy[j] * fx[j] * pb_yzz[j] + 

                                1.5 * pa_y[j] * fx[j] * pb_yyzz[j] + pa_yyy[j] * pb_yyzz[j]) * s_0_0[j] + (-3.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - pa_yyy[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_yy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                2.5 * pa_yyy[j] * fz[j] * fx[j] * fx[j] + 15.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                22.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 1.5 * pa_y[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                1.5 * pa_y[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                1.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_yzz[j] - pa_yyy[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                pa_yyy[j] * fz[j] * fgb[j] * pb_zz[j] + 7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 6.0 * pa_yyy[j] * fz[j] * pb_yy[j] * fx[j] + 

                                6.0 * pa_yyy[j] * fz[j] * fx[j] * pb_zz[j] + 36.0 * pa_yy[j] * fz[j] * fx[j] * pb_yzz[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_yyzz[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_yyzz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_yyy_yzzz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 2.25 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_zzz[j] + 1.5 * pa_yyy[j] * pb_yz[j] * fx[j] + 

                                1.5 * pa_yy[j] * fx[j] * pb_zzz[j] + 1.5 * pa_y[j] * fx[j] * pb_yzzz[j] + pa_yyy[j] * pb_yzzz[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 4.5 * pa_yy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 22.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 

                                4.5 * pa_y[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 4.5 * pa_y[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_zzz[j] - 3.0 * pa_yyy[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                22.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 

                                18.0 * pa_yyy[j] * fz[j] * pb_yz[j] * fx[j] + 18.0 * pa_yy[j] * fz[j] * fx[j] * pb_zzz[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_yzzz[j] + 18.0 * pa_y[j] * fz[j] * fx[j] * pb_yzzz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_yyy_zzzz[j] = (1.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yyy[j] * fx[j] * fx[j] + 

                                4.5 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_yyy[j] * pb_zz[j] * fx[j] + 1.5 * pa_y[j] * fx[j] * pb_zzzz[j] + 

                                pa_yyy[j] * pb_zzzz[j]) * s_0_0[j] + (-4.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                2.25 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * pa_yyy[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_y[j] * fz[j] * fx[j] * fx[j] * fx[j] + 7.5 * pa_yyy[j] * fz[j] * fx[j] * fx[j] - 

                                9.0 * pa_y[j] * fx[j] * pb_zz[j] * fz[j] * fgb[j] - 9.0 * pa_y[j] * fz[j] * fga[j] * pb_zz[j] * fx[j] - 

                                6.0 * pa_yyy[j] * pb_zz[j] * fz[j] * fgb[j] + 45.0 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                36.0 * pa_yyy[j] * fz[j] * pb_zz[j] * fx[j] - 3.0 * pa_y[j] * fz[j] * fga[j] * pb_zzzz[j] + 

                                18.0 * pa_y[j] * fz[j] * fx[j] * pb_zzzz[j] + 14.0 * pa_yyy[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (21) = (105,110)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, r_0_0, s_0_0, t_yyz_xxxx, t_yyz_xxxy, t_yyz_xxxz, t_yyz_xxyy, t_yyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xxxx[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_yyz[j] * fx[j] * fx[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 3.0 * pa_yyz[j] * pb_xx[j] * fx[j] + 0.5 * fx[j] * pa_z[j] * pb_xxxx[j] + 

                                pa_yyz[j] * pb_xxxx[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 3.0 * pa_yyz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 7.5 * pa_yyz[j] * fz[j] * fx[j] * fx[j] - 

                                3.0 * fx[j] * pa_z[j] * pb_xx[j] * fz[j] * fgb[j] - 3.0 * fz[j] * fga[j] * pa_z[j] * pb_xx[j] * fx[j] - 

                                6.0 * pa_yyz[j] * pb_xx[j] * fz[j] * fgb[j] + 15.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] + 

                                36.0 * pa_yyz[j] * fz[j] * pb_xx[j] * fx[j] - fz[j] * fga[j] * pa_z[j] * pb_xxxx[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxxx[j] + 

                                14.0 * pa_yyz[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_yyz_xxxy[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 1.5 * pa_yyz[j] * pb_xy[j] * fx[j] + pa_yz[j] * fx[j] * pb_xxx[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_xxxy[j] + pa_yyz[j] * pb_xxxy[j]) * s_0_0[j] + (-3.0 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                15.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 1.5 * fx[j] * pa_z[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_z[j] * pb_xy[j] * fx[j] - 3.0 * pa_yyz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] + 18.0 * pa_yyz[j] * fz[j] * pb_xy[j] * fx[j] + 

                                12.0 * pa_yz[j] * fx[j] * fz[j] * pb_xxx[j] - fz[j] * fga[j] * pa_z[j] * pb_xxxy[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxxy[j] + 

                                14.0 * pa_yyz[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_yyz_xxxz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_yyz[j] * pb_xz[j] * fx[j] + 

                                0.5 * pa_yy[j] * fx[j] * pb_xxx[j] + 0.5 * fx[j] * pa_z[j] * pb_xxxz[j] + pa_yyz[j] * pb_xxxz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - 1.5 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] - 

                                1.5 * fx[j] * pa_z[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_z[j] * pb_xz[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pb_xxx[j] - 3.0 * pa_yyz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 

                                18.0 * pa_yyz[j] * fz[j] * pb_xz[j] * fx[j] + 6.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxx[j] - fz[j] * fga[j] * pa_z[j] * pb_xxxz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxxz[j] + 14.0 * pa_yyz[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_yyz_xxyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * pa_yyz[j] * fx[j] * fx[j] + 

                                pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * pa_yyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyz[j] * fx[j] * pb_yy[j] + 

                                2.0 * pa_yz[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pa_z[j] * pb_xxyy[j] + pa_yyz[j] * pb_xxyy[j]) * s_0_0[j] + (-fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 0.25 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                pa_yyz[j] * fx[j] * fz[j] * fgb[j] - 2.0 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] + 

                                2.5 * pa_yyz[j] * fz[j] * fx[j] * fx[j] + 10.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] - 0.5 * fx[j] * pa_z[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * pa_z[j] * fz[j] * fgb[j] * pb_yy[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * pb_xx[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_yy[j] - pa_yyz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_yyz[j] * fz[j] * fgb[j] * pb_yy[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] + 6.0 * pa_yyz[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_yyz[j] * fz[j] * fx[j] * pb_yy[j] + 24.0 * pa_yz[j] * fx[j] * fz[j] * pb_xxy[j] - fz[j] * fga[j] * pa_z[j] * pb_xxyy[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxyy[j] + 14.0 * pa_yyz[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_yyz_xxyz[j] = (0.25 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 0.5 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.5 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_yyz[j] * fx[j] * pb_yz[j] + 

                                0.5 * pa_yy[j] * fx[j] * pb_xxy[j] + pa_yz[j] * fx[j] * pb_xxz[j] + 0.5 * fx[j] * pa_z[j] * pb_xxyz[j] + 

                                pa_yyz[j] * pb_xxyz[j]) * s_0_0[j] + (-0.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                2.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.25 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - 0.5 * pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                2.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 5.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                5.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] - 0.5 * fx[j] * pa_z[j] * fz[j] * fgb[j] * pb_yz[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_yz[j] - 0.5 * fz[j] * fga[j] * fx[j] * pb_xxy[j] - pa_yyz[j] * fz[j] * fgb[j] * pb_yz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 

                                6.0 * pa_yyz[j] * fz[j] * fx[j] * pb_yz[j] + 6.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxy[j] + 

                                12.0 * pa_yz[j] * fx[j] * fz[j] * pb_xxz[j] - fz[j] * fga[j] * pa_z[j] * pb_xxyz[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxyz[j] + 

                                14.0 * pa_yyz[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (22) = (110,115)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, \
                                     r_0_0, s_0_0, t_yyz_xxzz, t_yyz_xyyy, t_yyz_xyyz, t_yyz_xyzz, t_yyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_xxzz[j] = (0.125 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_yyz[j] * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.25 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_xxz[j] + 

                                0.5 * pa_yyz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yyz[j] * fx[j] * pb_zz[j] + pa_yy[j] * fx[j] * pb_xxz[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_xxzz[j] + pa_yyz[j] * pb_xxzz[j]) * s_0_0[j] + (-0.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.25 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - pa_yyz[j] * fx[j] * fz[j] * fgb[j] - pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                2.5 * pa_yyz[j] * fz[j] * fx[j] * fx[j] + 5.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 

                                0.5 * fx[j] * pa_z[j] * pb_xx[j] * fz[j] * fgb[j] - 0.5 * fx[j] * pa_z[j] * fz[j] * fgb[j] * pb_zz[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * pb_xx[j] * fx[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_zz[j] - 

                                fz[j] * fga[j] * fx[j] * pb_xxz[j] - pa_yyz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_yyz[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xx[j] + 2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 6.0 * pa_yyz[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_yyz[j] * fz[j] * fx[j] * pb_zz[j] + 12.0 * pa_yy[j] * fz[j] * fx[j] * pb_xxz[j] - fz[j] * fga[j] * pa_z[j] * pb_xxzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_xxzz[j] + 14.0 * pa_yyz[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_yyz_xyyy[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.25 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 1.5 * pa_yyz[j] * pb_xy[j] * fx[j] + 3.0 * pa_yz[j] * fx[j] * pb_xyy[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_xyyy[j] + pa_yyz[j] * pb_xyyy[j]) * s_0_0[j] + (-3.0 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                15.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 22.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] - 

                                1.5 * fx[j] * pa_z[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_z[j] * pb_xy[j] * fx[j] - 

                                3.0 * pa_yyz[j] * pb_xy[j] * fz[j] * fgb[j] + 18.0 * pa_yyz[j] * fz[j] * pb_xy[j] * fx[j] + 

                                36.0 * pa_yz[j] * fx[j] * fz[j] * pb_xyy[j] - fz[j] * fga[j] * pa_z[j] * pb_xyyy[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xyyy[j] + 

                                14.0 * pa_yyz[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_yyz_xyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.25 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xyy[j] + 

                                0.5 * pa_yyz[j] * pb_xz[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] * pb_xyy[j] + 2.0 * pa_yz[j] * fx[j] * pb_xyz[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_xyyz[j] + pa_yyz[j] * pb_xyyz[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.25 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - 

                                0.5 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 2.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                10.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] - 

                                0.5 * fx[j] * pa_z[j] * pb_xz[j] * fz[j] * fgb[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * pb_xz[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pb_xyy[j] - pa_yyz[j] * pb_xz[j] * fz[j] * fgb[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 

                                6.0 * pa_yyz[j] * fz[j] * pb_xz[j] * fx[j] + 6.0 * pa_yy[j] * fz[j] * fx[j] * pb_xyy[j] + 

                                24.0 * pa_yz[j] * fx[j] * fz[j] * pb_xyz[j] - fz[j] * fga[j] * pa_z[j] * pb_xyyz[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xyyz[j] + 

                                14.0 * pa_yyz[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_yyz_xyzz[j] = (0.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_yyz[j] * pb_xy[j] * fx[j] + 

                                pa_yy[j] * fx[j] * pb_xyz[j] + pa_yz[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_z[j] * pb_xyzz[j] + 

                                pa_yyz[j] * pb_xyzz[j]) * s_0_0[j] + (-pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                5.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 10.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] - 

                                0.5 * fx[j] * pa_z[j] * pb_xy[j] * fz[j] * fgb[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * pb_xy[j] * fx[j] - 

                                fz[j] * fga[j] * fx[j] * pb_xyz[j] - pa_yyz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xy[j] + 5.0 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * pa_yyz[j] * fz[j] * pb_xy[j] * fx[j] + 12.0 * pa_yy[j] * fz[j] * fx[j] * pb_xyz[j] + 

                                12.0 * pa_yz[j] * fx[j] * fz[j] * pb_xzz[j] - fz[j] * fga[j] * pa_z[j] * pb_xyzz[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_xyzz[j] + 

                                14.0 * pa_yyz[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_yyz_xzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_yy[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xzz[j] + 1.5 * pa_yyz[j] * pb_xz[j] * fx[j] + 

                                1.5 * pa_yy[j] * fx[j] * pb_xzz[j] + 0.5 * fx[j] * pa_z[j] * pb_xzzz[j] + pa_yyz[j] * pb_xzzz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_x[j] - 1.5 * pa_yy[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_x[j] - 

                                1.5 * fx[j] * pa_z[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_z[j] * pb_xz[j] * fx[j] - 

                                1.5 * fz[j] * fga[j] * fx[j] * pb_xzz[j] - 3.0 * pa_yyz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 

                                18.0 * pa_yyz[j] * fz[j] * pb_xz[j] * fx[j] + 18.0 * pa_yy[j] * fz[j] * fx[j] * pb_xzz[j] - fz[j] * fga[j] * pa_z[j] * pb_xzzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_xzzz[j] + 14.0 * pa_yyz[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (23) = (115,120)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, \
                                     s_0_0, t_yyz_yyyy, t_yyz_yyyz, t_yyz_yyzz, t_yyz_yzzz, t_yyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyz_yyyy[j] = (1.875 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * pa_yyz[j] * fx[j] * fx[j] + 

                                6.0 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 4.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 3.0 * pa_yyz[j] * pb_yy[j] * fx[j] + 

                                4.0 * pa_yz[j] * fx[j] * pb_yyy[j] + 0.5 * fx[j] * pa_z[j] * pb_yyyy[j] + pa_yyz[j] * pb_yyyy[j]) * s_0_0[j] + (-4.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] - 0.75 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                3.0 * pa_yyz[j] * fx[j] * fz[j] * fgb[j] - 12.0 * pa_yz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                7.5 * pa_yyz[j] * fz[j] * fx[j] * fx[j] + 60.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                45.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] - 3.0 * fx[j] * pa_z[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                3.0 * fz[j] * fga[j] * pa_z[j] * pb_yy[j] * fx[j] - 6.0 * pa_yyz[j] * pb_yy[j] * fz[j] * fgb[j] + 

                                36.0 * pa_yyz[j] * fz[j] * pb_yy[j] * fx[j] + 48.0 * pa_yz[j] * fx[j] * fz[j] * pb_yyy[j] - fz[j] * fga[j] * pa_z[j] * pb_yyyy[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_yyyy[j] + 14.0 * pa_yyz[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_yyz_yyyz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 

                                2.25 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_yyz[j] * pb_yz[j] * fx[j] + 

                                0.5 * pa_yy[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yz[j] * fx[j] * pb_yyz[j] + 0.5 * fx[j] * pa_z[j] * pb_yyyz[j] + 

                                pa_yyz[j] * pb_yyyz[j]) * s_0_0[j] + (-1.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                6.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - 

                                1.5 * pa_yy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 3.0 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                7.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 15.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                15.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 22.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] - 

                                1.5 * fx[j] * pa_z[j] * pb_yz[j] * fz[j] * fgb[j] - 1.5 * fz[j] * fga[j] * pa_z[j] * pb_yz[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * pb_yyy[j] - 3.0 * pa_yyz[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 18.0 * pa_yyz[j] * fz[j] * pb_yz[j] * fx[j] + 

                                6.0 * pa_yy[j] * fz[j] * fx[j] * pb_yyy[j] + 36.0 * pa_yz[j] * fx[j] * fz[j] * pb_yyz[j] - fz[j] * fga[j] * pa_z[j] * pb_yyyz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_yyyz[j] + 14.0 * pa_yyz[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_yyz_yyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_yyz[j] * fx[j] * fx[j] + 0.5 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 

                                2.0 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 

                                0.25 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.5 * fx[j] * fx[j] * pb_yyz[j] + 0.5 * pa_yyz[j] * pb_yy[j] * fx[j] + 

                                0.5 * pa_yyz[j] * fx[j] * pb_zz[j] + pa_yy[j] * fx[j] * pb_yyz[j] + 2.0 * pa_yz[j] * fx[j] * pb_yzz[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_yyzz[j] + pa_yyz[j] * pb_yyzz[j]) * s_0_0[j] + (-fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 0.25 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                0.5 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - pa_yyz[j] * fx[j] * fz[j] * fgb[j] - pa_yy[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                2.0 * pa_yz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 2.5 * pa_yyz[j] * fz[j] * fx[j] * fx[j] + 

                                5.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 10.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                20.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] - 

                                0.5 * fx[j] * pa_z[j] * pb_yy[j] * fz[j] * fgb[j] - 0.5 * fx[j] * pa_z[j] * fz[j] * fgb[j] * pb_zz[j] - 

                                0.5 * fz[j] * fga[j] * pa_z[j] * pb_yy[j] * fx[j] - 0.5 * fz[j] * fga[j] * pa_z[j] * fx[j] * pb_zz[j] - 

                                fz[j] * fga[j] * fx[j] * pb_yyz[j] - pa_yyz[j] * pb_yy[j] * fz[j] * fgb[j] - pa_yyz[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yy[j] + 5.0 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 

                                6.0 * pa_yyz[j] * fz[j] * pb_yy[j] * fx[j] + 6.0 * pa_yyz[j] * fz[j] * fx[j] * pb_zz[j] + 

                                12.0 * pa_yy[j] * fz[j] * fx[j] * pb_yyz[j] + 24.0 * pa_yz[j] * fx[j] * fz[j] * pb_yzz[j] - fz[j] * fga[j] * pa_z[j] * pb_yyzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_yyzz[j] + 14.0 * pa_yyz[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_yyz_yzzz[j] = (0.75 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_yy[j] * fx[j] * fx[j] * pb_y[j] + 1.5 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 

                                0.75 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_yyz[j] * pb_yz[j] * fx[j] + 

                                1.5 * pa_yy[j] * fx[j] * pb_yzz[j] + pa_yz[j] * fx[j] * pb_zzz[j] + 0.5 * fx[j] * pa_z[j] * pb_yzzz[j] + 

                                pa_yyz[j] * pb_yzzz[j]) * s_0_0[j] + (-1.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                6.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.75 * fz[j] * fga[j] * fx[j] * fx[j] * pb_y[j] - 1.5 * pa_yy[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                3.0 * pa_yz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                7.5 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 15.0 * pa_yz[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                15.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 1.5 * fx[j] * pa_z[j] * pb_yz[j] * fz[j] * fgb[j] - 

                                1.5 * fz[j] * fga[j] * pa_z[j] * pb_yz[j] * fx[j] - 1.5 * fz[j] * fga[j] * fx[j] * pb_yzz[j] - 

                                3.0 * pa_yyz[j] * pb_yz[j] * fz[j] * fgb[j] + 7.5 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_yz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 18.0 * pa_yyz[j] * fz[j] * pb_yz[j] * fx[j] + 

                                18.0 * pa_yy[j] * fz[j] * fx[j] * pb_yzz[j] + 12.0 * pa_yz[j] * fx[j] * fz[j] * pb_zzz[j] - fz[j] * fga[j] * pa_z[j] * pb_yzzz[j] + 

                                6.0 * fx[j] * fz[j] * pa_z[j] * pb_yzzz[j] + 14.0 * pa_yyz[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_yyz_zzzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_yyz[j] * fx[j] * fx[j] + 3.0 * pa_yy[j] * fx[j] * fx[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 

                                fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_yyz[j] * pb_zz[j] * fx[j] + 2.0 * pa_yy[j] * fx[j] * pb_zzz[j] + 

                                0.5 * fx[j] * pa_z[j] * pb_zzzz[j] + pa_yyz[j] * pb_zzzz[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] - 

                                3.0 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 0.75 * fz[j] * fga[j] * pa_z[j] * fx[j] * fx[j] - 

                                3.0 * fz[j] * fga[j] * fx[j] * fx[j] * pb_z[j] - 3.0 * pa_yyz[j] * fx[j] * fz[j] * fgb[j] - 

                                6.0 * pa_yy[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pa_z[j] + 

                                12.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 7.5 * pa_yyz[j] * fz[j] * fx[j] * fx[j] + 

                                30.0 * pa_yy[j] * fz[j] * fx[j] * fx[j] * pb_z[j] - 3.0 * fx[j] * pa_z[j] * pb_zz[j] * fz[j] * fgb[j] - 

                                3.0 * fz[j] * fga[j] * pa_z[j] * pb_zz[j] * fx[j] - 2.0 * fz[j] * fga[j] * fx[j] * pb_zzz[j] - 

                                6.0 * pa_yyz[j] * pb_zz[j] * fz[j] * fgb[j] + 15.0 * fx[j] * fx[j] * fz[j] * pa_z[j] * pb_zz[j] + 

                                10.0 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 36.0 * pa_yyz[j] * fz[j] * pb_zz[j] * fx[j] + 

                                24.0 * pa_yy[j] * fz[j] * fx[j] * pb_zzz[j] - fz[j] * fga[j] * pa_z[j] * pb_zzzz[j] + 6.0 * fx[j] * fz[j] * pa_z[j] * pb_zzzz[j] + 

                                14.0 * pa_yyz[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (24) = (120,125)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, pb_xxxx, \
                                     pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, \
                                     pb_z, r_0_0, s_0_0, t_yzz_xxxx, t_yzz_xxxy, t_yzz_xxxz, t_yzz_xxyy, t_yzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xxxx[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yzz[j] * fx[j] * fx[j] + 

                                1.5 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_yzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_y[j] * fx[j] * pb_xxxx[j] + 

                                pa_yzz[j] * pb_xxxx[j]) * s_0_0[j] + (-1.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * pa_yzz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_y[j] * fz[j] * fx[j] * fx[j] * fx[j] + 7.5 * pa_yzz[j] * fz[j] * fx[j] * fx[j] - 

                                3.0 * pa_y[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 3.0 * pa_y[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                6.0 * pa_yzz[j] * pb_xx[j] * fz[j] * fgb[j] + 15.0 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 

                                36.0 * pa_yzz[j] * fz[j] * pb_xx[j] * fx[j] - pa_y[j] * fz[j] * fga[j] * pb_xxxx[j] + 6.0 * pa_y[j] * fz[j] * fx[j] * pb_xxxx[j] + 

                                14.0 * pa_yzz[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_yzz_xxxy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.25 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_yzz[j] * pb_xy[j] * fx[j] + 

                                0.5 * fx[j] * pa_zz[j] * pb_xxx[j] + 0.5 * pa_y[j] * fx[j] * pb_xxxy[j] + pa_yzz[j] * pb_xxxy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 1.5 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] - 

                                1.5 * pa_y[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pb_xxx[j] - 3.0 * pa_yzz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 

                                18.0 * pa_yzz[j] * fz[j] * pb_xy[j] * fx[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxx[j] - pa_y[j] * fz[j] * fga[j] * pb_xxxy[j] + 

                                6.0 * pa_y[j] * fz[j] * fx[j] * pb_xxxy[j] + 14.0 * pa_yzz[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_yzz_xxxz[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yzz[j] * pb_xz[j] * fx[j] + pa_yz[j] * fx[j] * pb_xxx[j] + 

                                0.5 * pa_y[j] * fx[j] * pb_xxxz[j] + pa_yzz[j] * pb_xxxz[j]) * s_0_0[j] + (-3.0 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                15.0 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] - 1.5 * pa_y[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                1.5 * pa_y[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 3.0 * pa_yzz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 18.0 * pa_yzz[j] * fz[j] * pb_xz[j] * fx[j] + 

                                12.0 * pa_yz[j] * fz[j] * fx[j] * pb_xxx[j] - pa_y[j] * fz[j] * fga[j] * pb_xxxz[j] + 6.0 * pa_y[j] * fz[j] * fx[j] * pb_xxxz[j] + 

                                14.0 * pa_yzz[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_yzz_xxyy[j] = (0.125 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.25 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_yzz[j] * fx[j] * fx[j] + 0.5 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 

                                0.25 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 0.25 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * fx[j] * fx[j] * pb_xxy[j] + 

                                0.5 * pa_yzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yzz[j] * fx[j] * pb_yy[j] + fx[j] * pa_zz[j] * pb_xxy[j] + 

                                0.5 * pa_y[j] * fx[j] * pb_xxyy[j] + pa_yzz[j] * pb_xxyy[j]) * s_0_0[j] + (-0.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.25 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - pa_yzz[j] * fx[j] * fz[j] * fgb[j] + pa_y[j] * fz[j] * fx[j] * fx[j] * fx[j] - 

                                fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_y[j] + 2.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 

                                2.5 * pa_yzz[j] * fz[j] * fx[j] * fx[j] + 5.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] - 

                                0.5 * pa_y[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 0.5 * pa_y[j] * fx[j] * fz[j] * fgb[j] * pb_yy[j] - 

                                0.5 * pa_y[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 0.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - 

                                fx[j] * fz[j] * fga[j] * pb_xxy[j] - pa_yzz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_yzz[j] * fz[j] * fgb[j] * pb_yy[j] + 

                                2.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 2.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 6.0 * pa_yzz[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_yzz[j] * fz[j] * fx[j] * pb_yy[j] + 12.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxy[j] - pa_y[j] * fz[j] * fga[j] * pb_xxyy[j] + 

                                6.0 * pa_y[j] * fz[j] * fx[j] * pb_xxyy[j] + 14.0 * pa_yzz[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_yzz_xxyz[j] = (0.25 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.5 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.25 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 

                                0.5 * fx[j] * fx[j] * pa_z[j] * pb_xx[j] + 0.25 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.25 * fx[j] * fx[j] * pb_xxz[j] + 

                                0.5 * pa_yzz[j] * fx[j] * pb_yz[j] + pa_yz[j] * fx[j] * pb_xxy[j] + 0.5 * fx[j] * pa_zz[j] * pb_xxz[j] + 

                                0.5 * pa_y[j] * fx[j] * pb_xxyz[j] + pa_yzz[j] * pb_xxyz[j]) * s_0_0[j] + (-0.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                2.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] - 0.25 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_z[j] + fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                5.0 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 

                                5.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xx[j] - 0.5 * pa_y[j] * fx[j] * fz[j] * fgb[j] * pb_yz[j] - 

                                0.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - 0.5 * fx[j] * fz[j] * fga[j] * pb_xxz[j] - pa_yzz[j] * fz[j] * fgb[j] * pb_yz[j] + 

                                2.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 

                                6.0 * pa_yzz[j] * fz[j] * fx[j] * pb_yz[j] + 12.0 * pa_yz[j] * fz[j] * fx[j] * pb_xxy[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xxz[j] - pa_y[j] * fz[j] * fga[j] * pb_xxyz[j] + 6.0 * pa_y[j] * fz[j] * fx[j] * pb_xxyz[j] + 

                                14.0 * pa_yzz[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (25) = (125,130)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, \
                                     r_0_0, s_0_0, t_yzz_xxzz, t_yzz_xyyy, t_yzz_xyyz, t_yzz_xyzz, t_yzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_xxzz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_yzz[j] * fx[j] * fx[j] + 

                                pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xx[j] + 

                                0.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * pa_yzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_yzz[j] * fx[j] * pb_zz[j] + 

                                2.0 * pa_yz[j] * fx[j] * pb_xxz[j] + 0.5 * pa_y[j] * fx[j] * pb_xxzz[j] + pa_yzz[j] * pb_xxzz[j]) * s_0_0[j] + (-pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.25 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                pa_yzz[j] * fx[j] * fz[j] * fgb[j] - 2.0 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 

                                2.5 * pa_yzz[j] * fz[j] * fx[j] * fx[j] + 10.0 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                7.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] - 0.5 * pa_y[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                0.5 * pa_y[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 0.5 * pa_y[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                0.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - pa_yzz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_yzz[j] * fz[j] * fgb[j] * pb_zz[j] + 

                                2.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 6.0 * pa_yzz[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_yzz[j] * fz[j] * fx[j] * pb_zz[j] + 24.0 * pa_yz[j] * fz[j] * fx[j] * pb_xxz[j] - pa_y[j] * fz[j] * fga[j] * pb_xxzz[j] + 

                                6.0 * pa_y[j] * fz[j] * fx[j] * pb_xxzz[j] + 14.0 * pa_yzz[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_yzz_xyyy[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + 

                                0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] + 1.5 * pa_yzz[j] * pb_xy[j] * fx[j] + 

                                1.5 * fx[j] * pa_zz[j] * pb_xyy[j] + 0.5 * pa_y[j] * fx[j] * pb_xyyy[j] + pa_yzz[j] * pb_xyyy[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 1.5 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] - 

                                1.5 * pa_y[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_xyy[j] - 3.0 * pa_yzz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 

                                18.0 * pa_yzz[j] * fz[j] * pb_xy[j] * fx[j] + 18.0 * fx[j] * pa_zz[j] * fz[j] * pb_xyy[j] - pa_y[j] * fz[j] * fga[j] * pb_xyyy[j] + 

                                6.0 * pa_y[j] * fz[j] * fx[j] * pb_xyyy[j] + 14.0 * pa_yzz[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_yzz_xyyz[j] = (0.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xy[j] + 

                                0.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 0.5 * fx[j] * fx[j] * pb_xyz[j] + 0.5 * pa_yzz[j] * pb_xz[j] * fx[j] + 

                                pa_yz[j] * fx[j] * pb_xyy[j] + fx[j] * pa_zz[j] * pb_xyz[j] + 0.5 * pa_y[j] * fx[j] * pb_xyyz[j] + 

                                pa_yzz[j] * pb_xyyz[j]) * s_0_0[j] + (-pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                5.0 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 10.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xy[j] - 

                                0.5 * pa_y[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 0.5 * pa_y[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                fx[j] * fz[j] * fga[j] * pb_xyz[j] - pa_yzz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                2.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 5.0 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * pa_yzz[j] * fz[j] * pb_xz[j] * fx[j] + 12.0 * pa_yz[j] * fz[j] * fx[j] * pb_xyy[j] + 

                                12.0 * fx[j] * pa_zz[j] * fz[j] * pb_xyz[j] - pa_y[j] * fz[j] * fga[j] * pb_xyyz[j] + 6.0 * pa_y[j] * fz[j] * fx[j] * pb_xyyz[j] + 

                                14.0 * pa_yzz[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_yzz_xyzz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_xy[j] + 

                                0.25 * fx[j] * fx[j] * pa_zz[j] * pb_x[j] + fx[j] * fx[j] * pa_z[j] * pb_xz[j] + 0.25 * fx[j] * fx[j] * pb_xzz[j] + 

                                0.5 * pa_yzz[j] * pb_xy[j] * fx[j] + 2.0 * pa_yz[j] * fx[j] * pb_xyz[j] + 0.5 * fx[j] * pa_zz[j] * pb_xzz[j] + 

                                0.5 * pa_y[j] * fx[j] * pb_xyzz[j] + pa_yzz[j] * pb_xyzz[j]) * s_0_0[j] + (3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                0.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 0.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                0.5 * fx[j] * pa_zz[j] * pb_x[j] * fz[j] * fgb[j] + 7.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] + 

                                2.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_x[j] + 10.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_xz[j] - 

                                0.5 * pa_y[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 0.5 * pa_y[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pb_xzz[j] - pa_yzz[j] * pb_xy[j] * fz[j] * fgb[j] + 2.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 

                                6.0 * pa_yzz[j] * fz[j] * pb_xy[j] * fx[j] + 24.0 * pa_yz[j] * fz[j] * fx[j] * pb_xyz[j] + 

                                6.0 * fx[j] * pa_zz[j] * fz[j] * pb_xzz[j] - pa_y[j] * fz[j] * fga[j] * pb_xyzz[j] + 6.0 * pa_y[j] * fz[j] * fx[j] * pb_xyzz[j] + 

                                14.0 * pa_yzz[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_yzz_xzzz[j] = (1.5 * pa_yz[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.25 * pa_y[j] * fx[j] * fx[j] * pb_xz[j] + 1.5 * pa_yzz[j] * pb_xz[j] * fx[j] + 3.0 * pa_yz[j] * fx[j] * pb_xzz[j] + 

                                0.5 * pa_y[j] * fx[j] * pb_xzzz[j] + pa_yzz[j] * pb_xzzz[j]) * s_0_0[j] + (-3.0 * pa_yz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                15.0 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 22.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] - 

                                1.5 * pa_y[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                3.0 * pa_yzz[j] * pb_xz[j] * fz[j] * fgb[j] + 18.0 * pa_yzz[j] * fz[j] * pb_xz[j] * fx[j] + 

                                36.0 * pa_yz[j] * fz[j] * fx[j] * pb_xzz[j] - pa_y[j] * fz[j] * fga[j] * pb_xzzz[j] + 6.0 * pa_y[j] * fz[j] * fx[j] * pb_xzzz[j] + 

                                14.0 * pa_yzz[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (26) = (130,135)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, \
                                     s_0_0, t_yzz_yyyy, t_yzz_yyyz, t_yzz_yyzz, t_yzz_yzzz, t_yzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzz_yyyy[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 1.5 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_yzz[j] * fx[j] * fx[j] + 3.0 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 1.5 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 

                                fx[j] * fx[j] * pb_yyy[j] + 3.0 * pa_yzz[j] * pb_yy[j] * fx[j] + 2.0 * fx[j] * pa_zz[j] * pb_yyy[j] + 

                                0.5 * pa_y[j] * fx[j] * pb_yyyy[j] + pa_yzz[j] * pb_yyyy[j]) * s_0_0[j] + (-1.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                3.0 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 3.0 * pa_yzz[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_y[j] * fz[j] * fx[j] * fx[j] * fx[j] - 6.0 * fx[j] * pa_zz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                12.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_yzz[j] * fz[j] * fx[j] * fx[j] + 

                                30.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] - 3.0 * pa_y[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 2.0 * fx[j] * fz[j] * fga[j] * pb_yyy[j] - 

                                6.0 * pa_yzz[j] * pb_yy[j] * fz[j] * fgb[j] + 15.0 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                10.0 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 36.0 * pa_yzz[j] * fz[j] * pb_yy[j] * fx[j] + 

                                24.0 * fx[j] * pa_zz[j] * fz[j] * pb_yyy[j] - pa_y[j] * fz[j] * fga[j] * pb_yyyy[j] + 6.0 * pa_y[j] * fz[j] * fx[j] * pb_yyyy[j] + 

                                14.0 * pa_yzz[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_yzz_yyyz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 0.375 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 

                                1.5 * fx[j] * fx[j] * pa_z[j] * pb_yy[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yyz[j] + 

                                1.5 * pa_yzz[j] * pb_yz[j] * fx[j] + pa_yz[j] * fx[j] * pb_yyy[j] + 1.5 * fx[j] * pa_zz[j] * pb_yyz[j] + 

                                0.5 * pa_y[j] * fx[j] * pb_yyyz[j] + pa_yzz[j] * pb_yyyz[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] - 0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 3.0 * pa_yz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                1.5 * fx[j] * pa_zz[j] * fz[j] * fgb[j] * pb_z[j] + 3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                15.0 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 

                                15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yy[j] - 1.5 * pa_y[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 

                                1.5 * pa_y[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 1.5 * fx[j] * fz[j] * fga[j] * pb_yyz[j] - 

                                3.0 * pa_yzz[j] * pb_yz[j] * fz[j] * fgb[j] + 7.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 

                                7.5 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 18.0 * pa_yzz[j] * fz[j] * pb_yz[j] * fx[j] + 

                                12.0 * pa_yz[j] * fz[j] * fx[j] * pb_yyy[j] + 18.0 * fx[j] * pa_zz[j] * fz[j] * pb_yyz[j] - pa_y[j] * fz[j] * fga[j] * pb_yyyz[j] + 

                                6.0 * pa_y[j] * fz[j] * fx[j] * pb_yyyz[j] + 14.0 * pa_yzz[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_yzz_yyzz[j] = (0.375 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.25 * pa_yzz[j] * fx[j] * fx[j] + pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 0.75 * pa_y[j] * fx[j] * fx[j] * pb_yy[j] + 

                                0.5 * fx[j] * fx[j] * pa_zz[j] * pb_y[j] + 2.0 * fx[j] * fx[j] * pa_z[j] * pb_yz[j] + 

                                0.25 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 0.5 * fx[j] * fx[j] * pb_yzz[j] + 0.5 * pa_yzz[j] * pb_yy[j] * fx[j] + 

                                0.5 * pa_yzz[j] * fx[j] * pb_zz[j] + 2.0 * pa_yz[j] * fx[j] * pb_yyz[j] + fx[j] * pa_zz[j] * pb_yzz[j] + 

                                0.5 * pa_y[j] * fx[j] * pb_yyzz[j] + pa_yzz[j] * pb_yyzz[j]) * s_0_0[j] + (-pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                3.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                0.25 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 0.5 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                0.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - pa_yzz[j] * fx[j] * fz[j] * fgb[j] - 

                                2.0 * pa_yz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - fx[j] * pa_zz[j] * pb_y[j] * fz[j] * fgb[j] + 

                                2.5 * pa_yzz[j] * fz[j] * fx[j] * fx[j] + 10.0 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                7.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] + 5.0 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_y[j] + 

                                20.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_yz[j] - 0.5 * pa_y[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                0.5 * pa_y[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 0.5 * pa_y[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                0.5 * pa_y[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - fx[j] * fz[j] * fga[j] * pb_yzz[j] - pa_yzz[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                pa_yzz[j] * fz[j] * fgb[j] * pb_zz[j] + 2.5 * pa_y[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                5.0 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 6.0 * pa_yzz[j] * fz[j] * pb_yy[j] * fx[j] + 

                                6.0 * pa_yzz[j] * fz[j] * fx[j] * pb_zz[j] + 24.0 * pa_yz[j] * fz[j] * fx[j] * pb_yyz[j] + 

                                12.0 * fx[j] * pa_zz[j] * fz[j] * pb_yzz[j] - pa_y[j] * fz[j] * fga[j] * pb_yyzz[j] + 6.0 * pa_y[j] * fz[j] * fx[j] * pb_yyzz[j] + 

                                14.0 * pa_yzz[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_yzz_yzzz[j] = (0.75 * fx[j] * fx[j] * fx[j] * pa_z[j] + 1.125 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                1.5 * pa_yz[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_y[j] * fx[j] * fx[j] * pb_yz[j] + 

                                0.75 * fx[j] * fx[j] * pa_zz[j] * pb_z[j] + 1.5 * fx[j] * fx[j] * pa_z[j] * pb_zz[j] + 0.25 * fx[j] * fx[j] * pb_zzz[j] + 

                                1.5 * pa_yzz[j] * pb_yz[j] * fx[j] + 3.0 * pa_yz[j] * fx[j] * pb_yzz[j] + 0.5 * fx[j] * pa_zz[j] * pb_zzz[j] + 

                                0.5 * pa_y[j] * fx[j] * pb_yzzz[j] + pa_yzz[j] * pb_yzzz[j]) * s_0_0[j] + (-1.5 * fx[j] * fx[j] * pa_z[j] * fz[j] * fgb[j] + 

                                6.0 * fx[j] * fx[j] * fx[j] * pa_z[j] * fz[j] + 9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                0.75 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 

                                3.0 * pa_yz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 1.5 * fx[j] * pa_zz[j] * pb_z[j] * fz[j] * fgb[j] + 

                                15.0 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 22.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] + 

                                7.5 * fx[j] * fx[j] * pa_zz[j] * fz[j] * pb_z[j] + 15.0 * fx[j] * fx[j] * pa_z[j] * fz[j] * pb_zz[j] - 

                                1.5 * pa_y[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 1.5 * pa_y[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 

                                0.5 * fx[j] * fz[j] * fga[j] * pb_zzz[j] - 3.0 * pa_yzz[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                2.5 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 18.0 * pa_yzz[j] * fz[j] * pb_yz[j] * fx[j] + 

                                36.0 * pa_yz[j] * fz[j] * fx[j] * pb_yzz[j] + 6.0 * fx[j] * pa_zz[j] * fz[j] * pb_zzz[j] - pa_y[j] * fz[j] * fga[j] * pb_yzzz[j] + 

                                6.0 * pa_y[j] * fz[j] * fx[j] * pb_yzzz[j] + 14.0 * pa_yzz[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_yzz_zzzz[j] = (1.875 * pa_y[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_yzz[j] * fx[j] * fx[j] + 

                                6.0 * pa_yz[j] * fx[j] * fx[j] * pb_z[j] + 4.5 * pa_y[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * pa_yzz[j] * pb_zz[j] * fx[j] + 

                                4.0 * pa_yz[j] * fx[j] * pb_zzz[j] + 0.5 * pa_y[j] * fx[j] * pb_zzzz[j] + pa_yzz[j] * pb_zzzz[j]) * s_0_0[j] + (-4.5 * pa_y[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                15.0 * pa_y[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * pa_y[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                3.0 * pa_yzz[j] * fx[j] * fz[j] * fgb[j] - 12.0 * pa_yz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 

                                7.5 * pa_yzz[j] * fz[j] * fx[j] * fx[j] + 60.0 * pa_yz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                45.0 * pa_y[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 3.0 * pa_y[j] * fx[j] * pb_zz[j] * fz[j] * fgb[j] - 

                                3.0 * pa_y[j] * fz[j] * fga[j] * pb_zz[j] * fx[j] - 6.0 * pa_yzz[j] * pb_zz[j] * fz[j] * fgb[j] + 

                                36.0 * pa_yzz[j] * fz[j] * pb_zz[j] * fx[j] + 48.0 * pa_yz[j] * fz[j] * fx[j] * pb_zzz[j] - pa_y[j] * fz[j] * fga[j] * pb_zzzz[j] + 

                                6.0 * pa_y[j] * fz[j] * fx[j] * pb_zzzz[j] + 14.0 * pa_yzz[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (27) = (135,140)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, \
                                     pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, r_0_0, s_0_0, \
                                     t_zzz_xxxx, t_zzz_xxxy, t_zzz_xxxz, t_zzz_xxyy, t_zzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xxxx[j] = (1.125 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_zzz[j] * fx[j] * fx[j] + 

                                4.5 * pa_z[j] * fx[j] * fx[j] * pb_xx[j] + 3.0 * pa_zzz[j] * pb_xx[j] * fx[j] + 1.5 * pa_z[j] * fx[j] * pb_xxxx[j] + 

                                pa_zzz[j] * pb_xxxx[j]) * s_0_0[j] + (-4.5 * pa_z[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                2.25 * pa_z[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * pa_zzz[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_z[j] * fz[j] * fx[j] * fx[j] * fx[j] + 7.5 * pa_zzz[j] * fz[j] * fx[j] * fx[j] - 

                                9.0 * pa_z[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 9.0 * pa_z[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                6.0 * pa_zzz[j] * pb_xx[j] * fz[j] * fgb[j] + 45.0 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 

                                36.0 * pa_zzz[j] * fz[j] * pb_xx[j] * fx[j] - 3.0 * pa_z[j] * fz[j] * fga[j] * pb_xxxx[j] + 

                                18.0 * pa_z[j] * fz[j] * fx[j] * pb_xxxx[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xxxx[j]) * r_0_0[j];

                t_zzz_xxxy[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_zzz[j] * pb_xy[j] * fx[j] + 

                                1.5 * pa_z[j] * fx[j] * pb_xxxy[j] + pa_zzz[j] * pb_xxxy[j]) * s_0_0[j] + (-4.5 * pa_z[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                4.5 * pa_z[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 3.0 * pa_zzz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                22.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 18.0 * pa_zzz[j] * fz[j] * pb_xy[j] * fx[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_xxxy[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_xxxy[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xxxy[j]) * r_0_0[j];

                t_zzz_xxxz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] + 

                                2.25 * pa_z[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xxx[j] + 1.5 * pa_zzz[j] * pb_xz[j] * fx[j] + 

                                1.5 * pa_zz[j] * fx[j] * pb_xxx[j] + 1.5 * pa_z[j] * fx[j] * pb_xxxz[j] + pa_zzz[j] * pb_xxxz[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 4.5 * pa_zz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 22.5 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] - 

                                4.5 * pa_z[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 4.5 * pa_z[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_xxx[j] - 3.0 * pa_zzz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                22.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxx[j] + 

                                18.0 * pa_zzz[j] * fz[j] * pb_xz[j] * fx[j] + 18.0 * pa_zz[j] * fz[j] * fx[j] * pb_xxx[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_xxxz[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_xxxz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xxxz[j]) * r_0_0[j];

                t_zzz_xxyy[j] = (0.375 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.25 * pa_zzz[j] * fx[j] * fx[j] + 

                                0.75 * pa_z[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_z[j] * fx[j] * fx[j] * pb_yy[j] + 0.5 * pa_zzz[j] * pb_xx[j] * fx[j] + 

                                0.5 * pa_zzz[j] * fx[j] * pb_yy[j] + 1.5 * pa_z[j] * fx[j] * pb_xxyy[j] + pa_zzz[j] * pb_xxyy[j]) * s_0_0[j] + (-1.5 * pa_z[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                0.75 * pa_z[j] * fz[j] * fga[j] * fx[j] * fx[j] - pa_zzz[j] * fx[j] * fz[j] * fgb[j] + 3.0 * pa_z[j] * fz[j] * fx[j] * fx[j] * fx[j] + 

                                2.5 * pa_zzz[j] * fz[j] * fx[j] * fx[j] - 1.5 * pa_z[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                1.5 * pa_z[j] * fx[j] * fz[j] * fgb[j] * pb_yy[j] - 1.5 * pa_z[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                1.5 * pa_z[j] * fz[j] * fga[j] * fx[j] * pb_yy[j] - pa_zzz[j] * pb_xx[j] * fz[j] * fgb[j] - pa_zzz[j] * fz[j] * fgb[j] * pb_yy[j] + 

                                7.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_xx[j] + 7.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                6.0 * pa_zzz[j] * fz[j] * pb_xx[j] * fx[j] + 6.0 * pa_zzz[j] * fz[j] * fx[j] * pb_yy[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_xxyy[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_xxyy[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xxyy[j]) * r_0_0[j];

                t_zzz_xxyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_y[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] + 

                                0.75 * pa_z[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_xxy[j] + 0.5 * pa_zzz[j] * fx[j] * pb_yz[j] + 

                                1.5 * pa_zz[j] * fx[j] * pb_xxy[j] + 1.5 * pa_z[j] * fx[j] * pb_xxyz[j] + pa_zzz[j] * pb_xxyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 1.5 * pa_zz[j] * fx[j] * fz[j] * fgb[j] * pb_y[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 7.5 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] - 

                                1.5 * pa_z[j] * fx[j] * fz[j] * fgb[j] * pb_yz[j] - 1.5 * pa_z[j] * fz[j] * fga[j] * fx[j] * pb_yz[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_xxy[j] - pa_zzz[j] * fz[j] * fgb[j] * pb_yz[j] + 

                                7.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xxy[j] + 

                                6.0 * pa_zzz[j] * fz[j] * fx[j] * pb_yz[j] + 18.0 * pa_zz[j] * fz[j] * fx[j] * pb_xxy[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_xxyz[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_xxyz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xxyz[j]) * r_0_0[j];
            }

            // Batch of Integrals (28) = (140,145)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, r_0_0, s_0_0, \
                                     t_zzz_xxzz, t_zzz_xyyy, t_zzz_xyyz, t_zzz_xyzz, t_zzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_xxzz[j] = (1.125 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_zzz[j] * fx[j] * fx[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] + 

                                2.25 * pa_z[j] * fx[j] * fx[j] * pb_xx[j] + 0.75 * pa_z[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * fx[j] * fx[j] * pb_xxz[j] + 

                                0.5 * pa_zzz[j] * pb_xx[j] * fx[j] + 0.5 * pa_zzz[j] * fx[j] * pb_zz[j] + 3.0 * pa_zz[j] * fx[j] * pb_xxz[j] + 

                                1.5 * pa_z[j] * fx[j] * pb_xxzz[j] + pa_zzz[j] * pb_xxzz[j]) * s_0_0[j] + (-3.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_z[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * pa_z[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - pa_zzz[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_zz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                2.5 * pa_zzz[j] * fz[j] * fx[j] * fx[j] + 15.0 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                22.5 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_xx[j] - 1.5 * pa_z[j] * fx[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                1.5 * pa_z[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 1.5 * pa_z[j] * fz[j] * fga[j] * pb_xx[j] * fx[j] - 

                                1.5 * pa_z[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_xxz[j] - pa_zzz[j] * pb_xx[j] * fz[j] * fgb[j] - 

                                pa_zzz[j] * fz[j] * fgb[j] * pb_zz[j] + 7.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pb_xxz[j] + 6.0 * pa_zzz[j] * fz[j] * pb_xx[j] * fx[j] + 

                                6.0 * pa_zzz[j] * fz[j] * fx[j] * pb_zz[j] + 36.0 * pa_zz[j] * fz[j] * fx[j] * pb_xxz[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_xxzz[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_xxzz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xxzz[j]) * r_0_0[j];

                t_zzz_xyyy[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * pa_zzz[j] * pb_xy[j] * fx[j] + 

                                1.5 * pa_z[j] * fx[j] * pb_xyyy[j] + pa_zzz[j] * pb_xyyy[j]) * s_0_0[j] + (-4.5 * pa_z[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 

                                4.5 * pa_z[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 3.0 * pa_zzz[j] * pb_xy[j] * fz[j] * fgb[j] + 

                                22.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_xy[j] + 18.0 * pa_zzz[j] * fz[j] * pb_xy[j] * fx[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_xyyy[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_xyyy[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xyyy[j]) * r_0_0[j];

                t_zzz_xyyz[j] = (0.375 * fx[j] * fx[j] * fx[j] * pb_x[j] + 0.75 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] + 

                                0.75 * pa_z[j] * fx[j] * fx[j] * pb_xz[j] + 0.75 * fx[j] * fx[j] * pb_xyy[j] + 0.5 * pa_zzz[j] * pb_xz[j] * fx[j] + 

                                1.5 * pa_zz[j] * fx[j] * pb_xyy[j] + 1.5 * pa_z[j] * fx[j] * pb_xyyz[j] + pa_zzz[j] * pb_xyyz[j]) * s_0_0[j] + (-0.75 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 

                                0.75 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 1.5 * pa_zz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 

                                3.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] + 7.5 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] - 

                                1.5 * pa_z[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 1.5 * pa_z[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_xyy[j] - pa_zzz[j] * pb_xz[j] * fz[j] * fgb[j] + 

                                7.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_xz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_xyy[j] + 

                                6.0 * pa_zzz[j] * fz[j] * pb_xz[j] * fx[j] + 18.0 * pa_zz[j] * fz[j] * fx[j] * pb_xyy[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_xyyz[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_xyyz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xyyz[j]) * r_0_0[j];

                t_zzz_xyzz[j] = (2.25 * pa_z[j] * fx[j] * fx[j] * pb_xy[j] + 1.5 * fx[j] * fx[j] * pb_xyz[j] + 

                                0.5 * pa_zzz[j] * pb_xy[j] * fx[j] + 3.0 * pa_zz[j] * fx[j] * pb_xyz[j] + 1.5 * pa_z[j] * fx[j] * pb_xyzz[j] + 

                                pa_zzz[j] * pb_xyzz[j]) * s_0_0[j] + (22.5 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_xy[j] - 

                                1.5 * pa_z[j] * fx[j] * pb_xy[j] * fz[j] * fgb[j] - 1.5 * pa_z[j] * fz[j] * fga[j] * pb_xy[j] * fx[j] - 

                                3.0 * fx[j] * fz[j] * fga[j] * pb_xyz[j] - pa_zzz[j] * pb_xy[j] * fz[j] * fgb[j] + 15.0 * fx[j] * fx[j] * fz[j] * pb_xyz[j] + 

                                6.0 * pa_zzz[j] * fz[j] * pb_xy[j] * fx[j] + 36.0 * pa_zz[j] * fz[j] * fx[j] * pb_xyz[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_xyzz[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_xyzz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xyzz[j]) * r_0_0[j];

                t_zzz_xzzz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_x[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * pb_x[j] + 

                                6.75 * pa_z[j] * fx[j] * fx[j] * pb_xz[j] + 2.25 * fx[j] * fx[j] * pb_xzz[j] + 1.5 * pa_zzz[j] * pb_xz[j] * fx[j] + 

                                4.5 * pa_zz[j] * fx[j] * pb_xzz[j] + 1.5 * pa_z[j] * fx[j] * pb_xzzz[j] + pa_zzz[j] * pb_xzzz[j]) * s_0_0[j] + (15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_x[j] - 

                                2.25 * fx[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] - 2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_x[j] - 

                                4.5 * pa_zz[j] * fx[j] * pb_x[j] * fz[j] * fgb[j] + 22.5 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_x[j] + 

                                67.5 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_xz[j] - 4.5 * pa_z[j] * fx[j] * pb_xz[j] * fz[j] * fgb[j] - 

                                4.5 * pa_z[j] * fz[j] * fga[j] * pb_xz[j] * fx[j] - 4.5 * fx[j] * fz[j] * fga[j] * pb_xzz[j] - 

                                3.0 * pa_zzz[j] * pb_xz[j] * fz[j] * fgb[j] + 22.5 * fx[j] * fx[j] * fz[j] * pb_xzz[j] + 

                                18.0 * pa_zzz[j] * fz[j] * pb_xz[j] * fx[j] + 54.0 * pa_zz[j] * fz[j] * fx[j] * pb_xzz[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_xzzz[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_xzzz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_xzzz[j]) * r_0_0[j];
            }

            // Batch of Integrals (29) = (145,150)

            #pragma omp simd aligned(fga, fgb, fx, fz, pa_z, pa_zz, pa_zzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, r_0_0, s_0_0, t_zzz_yyyy, \
                                     t_zzz_yyyz, t_zzz_yyzz, t_zzz_yzzz, t_zzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzz_yyyy[j] = (1.125 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.75 * pa_zzz[j] * fx[j] * fx[j] + 

                                4.5 * pa_z[j] * fx[j] * fx[j] * pb_yy[j] + 3.0 * pa_zzz[j] * pb_yy[j] * fx[j] + 1.5 * pa_z[j] * fx[j] * pb_yyyy[j] + 

                                pa_zzz[j] * pb_yyyy[j]) * s_0_0[j] + (-4.5 * pa_z[j] * fx[j] * fx[j] * fz[j] * fgb[j] - 

                                2.25 * pa_z[j] * fz[j] * fga[j] * fx[j] * fx[j] - 3.0 * pa_zzz[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_z[j] * fz[j] * fx[j] * fx[j] * fx[j] + 7.5 * pa_zzz[j] * fz[j] * fx[j] * fx[j] - 

                                9.0 * pa_z[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 9.0 * pa_z[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                6.0 * pa_zzz[j] * pb_yy[j] * fz[j] * fgb[j] + 45.0 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_yy[j] + 

                                36.0 * pa_zzz[j] * fz[j] * pb_yy[j] * fx[j] - 3.0 * pa_z[j] * fz[j] * fga[j] * pb_yyyy[j] + 

                                18.0 * pa_z[j] * fz[j] * fx[j] * pb_yyyy[j] + 14.0 * pa_zzz[j] * fz[j] * pb_yyyy[j]) * r_0_0[j];

                t_zzz_yyyz[j] = (1.125 * fx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] + 

                                2.25 * pa_z[j] * fx[j] * fx[j] * pb_yz[j] + 0.75 * fx[j] * fx[j] * pb_yyy[j] + 1.5 * pa_zzz[j] * pb_yz[j] * fx[j] + 

                                1.5 * pa_zz[j] * fx[j] * pb_yyy[j] + 1.5 * pa_z[j] * fx[j] * pb_yyyz[j] + pa_zzz[j] * pb_yyyz[j]) * s_0_0[j] + (-2.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 

                                2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 4.5 * pa_zz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 

                                9.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] + 22.5 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] - 

                                4.5 * pa_z[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 4.5 * pa_z[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 

                                1.5 * fx[j] * fz[j] * fga[j] * pb_yyy[j] - 3.0 * pa_zzz[j] * pb_yz[j] * fz[j] * fgb[j] + 

                                22.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_yz[j] + 7.5 * fx[j] * fx[j] * fz[j] * pb_yyy[j] + 

                                18.0 * pa_zzz[j] * fz[j] * pb_yz[j] * fx[j] + 18.0 * pa_zz[j] * fz[j] * fx[j] * pb_yyy[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_yyyz[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_yyyz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_yyyz[j]) * r_0_0[j];

                t_zzz_yyzz[j] = (1.125 * pa_z[j] * fx[j] * fx[j] * fx[j] + 0.75 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.25 * pa_zzz[j] * fx[j] * fx[j] + 1.5 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] + 

                                2.25 * pa_z[j] * fx[j] * fx[j] * pb_yy[j] + 0.75 * pa_z[j] * fx[j] * fx[j] * pb_zz[j] + 1.5 * fx[j] * fx[j] * pb_yyz[j] + 

                                0.5 * pa_zzz[j] * pb_yy[j] * fx[j] + 0.5 * pa_zzz[j] * fx[j] * pb_zz[j] + 3.0 * pa_zz[j] * fx[j] * pb_yyz[j] + 

                                1.5 * pa_z[j] * fx[j] * pb_yyzz[j] + pa_zzz[j] * pb_yyzz[j]) * s_0_0[j] + (-3.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                9.0 * pa_z[j] * fx[j] * fx[j] * fx[j] * fz[j] - 0.75 * pa_z[j] * fz[j] * fga[j] * fx[j] * fx[j] - 

                                1.5 * fx[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] - 1.5 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - pa_zzz[j] * fx[j] * fz[j] * fgb[j] - 

                                3.0 * pa_zz[j] * fx[j] * fz[j] * fgb[j] * pb_z[j] + 6.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] + 

                                2.5 * pa_zzz[j] * fz[j] * fx[j] * fx[j] + 15.0 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 

                                22.5 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_yy[j] - 1.5 * pa_z[j] * fx[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                1.5 * pa_z[j] * fx[j] * fz[j] * fgb[j] * pb_zz[j] - 1.5 * pa_z[j] * fz[j] * fga[j] * pb_yy[j] * fx[j] - 

                                1.5 * pa_z[j] * fz[j] * fga[j] * fx[j] * pb_zz[j] - 3.0 * fx[j] * fz[j] * fga[j] * pb_yyz[j] - pa_zzz[j] * pb_yy[j] * fz[j] * fgb[j] - 

                                pa_zzz[j] * fz[j] * fgb[j] * pb_zz[j] + 7.5 * pa_z[j] * fz[j] * fx[j] * fx[j] * pb_zz[j] + 

                                15.0 * fx[j] * fx[j] * fz[j] * pb_yyz[j] + 6.0 * pa_zzz[j] * fz[j] * pb_yy[j] * fx[j] + 

                                6.0 * pa_zzz[j] * fz[j] * fx[j] * pb_zz[j] + 36.0 * pa_zz[j] * fz[j] * fx[j] * pb_yyz[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_yyzz[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_yyzz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_yyzz[j]) * r_0_0[j];

                t_zzz_yzzz[j] = (1.875 * fx[j] * fx[j] * fx[j] * pb_y[j] + 2.25 * pa_zz[j] * fx[j] * fx[j] * pb_y[j] + 

                                6.75 * pa_z[j] * fx[j] * fx[j] * pb_yz[j] + 2.25 * fx[j] * fx[j] * pb_yzz[j] + 1.5 * pa_zzz[j] * pb_yz[j] * fx[j] + 

                                4.5 * pa_zz[j] * fx[j] * pb_yzz[j] + 1.5 * pa_z[j] * fx[j] * pb_yzzz[j] + pa_zzz[j] * pb_yzzz[j]) * s_0_0[j] + (15.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_y[j] - 

                                2.25 * fx[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] - 2.25 * fx[j] * fx[j] * fz[j] * fga[j] * pb_y[j] - 

                                4.5 * pa_zz[j] * fx[j] * pb_y[j] * fz[j] * fgb[j] + 22.5 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_y[j] + 

                                67.5 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_yz[j] - 4.5 * pa_z[j] * fx[j] * pb_yz[j] * fz[j] * fgb[j] - 

                                4.5 * pa_z[j] * fz[j] * fga[j] * pb_yz[j] * fx[j] - 4.5 * fx[j] * fz[j] * fga[j] * pb_yzz[j] - 

                                3.0 * pa_zzz[j] * pb_yz[j] * fz[j] * fgb[j] + 22.5 * fx[j] * fx[j] * fz[j] * pb_yzz[j] + 

                                18.0 * pa_zzz[j] * fz[j] * pb_yz[j] * fx[j] + 54.0 * pa_zz[j] * fz[j] * fx[j] * pb_yzz[j] - 

                                3.0 * pa_z[j] * fz[j] * fga[j] * pb_yzzz[j] + 18.0 * pa_z[j] * fz[j] * fx[j] * pb_yzzz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_yzzz[j]) * r_0_0[j];

                t_zzz_zzzz[j] = (5.625 * pa_z[j] * fx[j] * fx[j] * fx[j] + 7.5 * fx[j] * fx[j] * fx[j] * pb_z[j] + 

                                0.75 * pa_zzz[j] * fx[j] * fx[j] + 9.0 * pa_zz[j] * fx[j] * fx[j] * pb_z[j] + 

                                13.5 * pa_z[j] * fx[j] * fx[j] * pb_zz[j] + 3.0 * fx[j] * fx[j] * pb_zzz[j] + 3.0 * pa_zzz[j] * pb_zz[j] * fx[j] + 

                                6.0 * pa_zz[j] * fx[j] * pb_zzz[j] + 1.5 * pa_z[j] * fx[j] * pb_zzzz[j] + pa_zzz[j] * pb_zzzz[j]) * s_0_0[j] + (-13.5 * pa_z[j] * fx[j] * fx[j] * fz[j] * fgb[j] + 

                                45.0 * pa_z[j] * fx[j] * fx[j] * fx[j] * fz[j] + 60.0 * fx[j] * fx[j] * fx[j] * fz[j] * pb_z[j] - 

                                2.25 * pa_z[j] * fz[j] * fga[j] * fx[j] * fx[j] - 9.0 * fx[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] - 

                                9.0 * fx[j] * fx[j] * fz[j] * fga[j] * pb_z[j] - 3.0 * pa_zzz[j] * fx[j] * fz[j] * fgb[j] - 

                                18.0 * pa_zz[j] * fx[j] * pb_z[j] * fz[j] * fgb[j] + 7.5 * pa_zzz[j] * fz[j] * fx[j] * fx[j] + 

                                90.0 * pa_zz[j] * fz[j] * fx[j] * fx[j] * pb_z[j] + 135.0 * pa_z[j] * fx[j] * fx[j] * fz[j] * pb_zz[j] - 

                                9.0 * pa_z[j] * fx[j] * pb_zz[j] * fz[j] * fgb[j] - 9.0 * pa_z[j] * fz[j] * fga[j] * pb_zz[j] * fx[j] - 

                                6.0 * fx[j] * fz[j] * fga[j] * pb_zzz[j] - 6.0 * pa_zzz[j] * pb_zz[j] * fz[j] * fgb[j] + 

                                30.0 * fx[j] * fx[j] * fz[j] * pb_zzz[j] + 36.0 * pa_zzz[j] * fz[j] * pb_zz[j] * fx[j] + 

                                72.0 * pa_zz[j] * fz[j] * fx[j] * pb_zzz[j] - 3.0 * pa_z[j] * fz[j] * fga[j] * pb_zzzz[j] + 

                                18.0 * pa_z[j] * fz[j] * fx[j] * pb_zzzz[j] + 14.0 * pa_zzz[j] * fz[j] * pb_zzzz[j]) * r_0_0[j];
            }

            idx++;
        }
    }


} // kinrecfunc namespace

