//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapRecFuncForGG.hpp"

#include "OverlapVecFuncForGG.hpp"

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

            auto t_yyzz_xxxx = primBuffer.data(225 * idx + 180);

            auto t_yyzz_xxxy = primBuffer.data(225 * idx + 181);

            auto t_yyzz_xxxz = primBuffer.data(225 * idx + 182);

            auto t_yyzz_xxyy = primBuffer.data(225 * idx + 183);

            auto t_yyzz_xxyz = primBuffer.data(225 * idx + 184);

            auto t_yyzz_xxzz = primBuffer.data(225 * idx + 185);

            auto t_yyzz_xyyy = primBuffer.data(225 * idx + 186);

            auto t_yyzz_xyyz = primBuffer.data(225 * idx + 187);

            auto t_yyzz_xyzz = primBuffer.data(225 * idx + 188);

            auto t_yyzz_xzzz = primBuffer.data(225 * idx + 189);

            auto t_yyzz_yyyy = primBuffer.data(225 * idx + 190);

            auto t_yyzz_yyyz = primBuffer.data(225 * idx + 191);

            auto t_yyzz_yyzz = primBuffer.data(225 * idx + 192);

            auto t_yyzz_yzzz = primBuffer.data(225 * idx + 193);

            auto t_yyzz_zzzz = primBuffer.data(225 * idx + 194);

            auto t_yzzz_xxxx = primBuffer.data(225 * idx + 195);

            auto t_yzzz_xxxy = primBuffer.data(225 * idx + 196);

            auto t_yzzz_xxxz = primBuffer.data(225 * idx + 197);

            auto t_yzzz_xxyy = primBuffer.data(225 * idx + 198);

            auto t_yzzz_xxyz = primBuffer.data(225 * idx + 199);

            auto t_yzzz_xxzz = primBuffer.data(225 * idx + 200);

            auto t_yzzz_xyyy = primBuffer.data(225 * idx + 201);

            auto t_yzzz_xyyz = primBuffer.data(225 * idx + 202);

            auto t_yzzz_xyzz = primBuffer.data(225 * idx + 203);

            auto t_yzzz_xzzz = primBuffer.data(225 * idx + 204);

            auto t_yzzz_yyyy = primBuffer.data(225 * idx + 205);

            auto t_yzzz_yyyz = primBuffer.data(225 * idx + 206);

            auto t_yzzz_yyzz = primBuffer.data(225 * idx + 207);

            auto t_yzzz_yzzz = primBuffer.data(225 * idx + 208);

            auto t_yzzz_zzzz = primBuffer.data(225 * idx + 209);

            auto t_zzzz_xxxx = primBuffer.data(225 * idx + 210);

            auto t_zzzz_xxxy = primBuffer.data(225 * idx + 211);

            auto t_zzzz_xxxz = primBuffer.data(225 * idx + 212);

            auto t_zzzz_xxyy = primBuffer.data(225 * idx + 213);

            auto t_zzzz_xxyz = primBuffer.data(225 * idx + 214);

            auto t_zzzz_xxzz = primBuffer.data(225 * idx + 215);

            auto t_zzzz_xyyy = primBuffer.data(225 * idx + 216);

            auto t_zzzz_xyyz = primBuffer.data(225 * idx + 217);

            auto t_zzzz_xyzz = primBuffer.data(225 * idx + 218);

            auto t_zzzz_xzzz = primBuffer.data(225 * idx + 219);

            auto t_zzzz_yyyy = primBuffer.data(225 * idx + 220);

            auto t_zzzz_yyyz = primBuffer.data(225 * idx + 221);

            auto t_zzzz_yyzz = primBuffer.data(225 * idx + 222);

            auto t_zzzz_yzzz = primBuffer.data(225 * idx + 223);

            auto t_zzzz_zzzz = primBuffer.data(225 * idx + 224);

            // Batch of Integrals (0) = (0,5)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, \
                                     s_0_0, t_xxxx_xxxx, t_xxxx_xxxy, t_xxxx_xxxz, t_xxxx_xxyy, t_xxxx_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xxxx[j] = ovlvecfunc::fvec_xxxx_xxxx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xxxx_xxxy[j] = ovlvecfunc::fvec_xxxx_xxxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xxxx_xxxz[j] = ovlvecfunc::fvec_xxxx_xxxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xxxx_xxyy[j] = ovlvecfunc::fvec_xxxx_xxyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], s_0_0[j]);

                t_xxxx_xxyz[j] = ovlvecfunc::fvec_xxxx_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xxyz[j], pb_xyz[j], pb_yz[j], s_0_0[j]);
            }

            // Batch of Integrals (1) = (5,10)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxx, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyyy, pb_xyyz, \
                                     pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yyy, pb_yyz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, \
                                     t_xxxx_xxzz, t_xxxx_xyyy, t_xxxx_xyyz, t_xxxx_xyzz, t_xxxx_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_xxzz[j] = ovlvecfunc::fvec_xxxx_xxzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], s_0_0[j]);

                t_xxxx_xyyy[j] = ovlvecfunc::fvec_xxxx_xyyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], s_0_0[j]);

                t_xxxx_xyyz[j] = ovlvecfunc::fvec_xxxx_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xyyz[j], pb_xz[j], pb_yyz[j], pb_z[j], s_0_0[j]);

                t_xxxx_xyzz[j] = ovlvecfunc::fvec_xxxx_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xy[j], pb_xyzz[j], pb_y[j], pb_yzz[j], s_0_0[j]);

                t_xxxx_xzzz[j] = ovlvecfunc::fvec_xxxx_xzzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxx[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (2) = (10,15)

            #pragma omp simd aligned(fx, pa_xx, pa_xxxx, pb_yy, pb_yyyy, pb_yyyz, pb_yyzz, pb_yz, pb_yzzz, pb_zz, \
                                     pb_zzzz, s_0_0, t_xxxx_yyyy, t_xxxx_yyyz, t_xxxx_yyzz, t_xxxx_yzzz, t_xxxx_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxx_yyyy[j] = ovlvecfunc::fvec_xxxx_yyyy_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_yy[j], pb_yyyy[j], s_0_0[j]);

                t_xxxx_yyyz[j] = ovlvecfunc::fvec_xxxx_yyyz_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_yyyz[j], pb_yz[j], s_0_0[j]);

                t_xxxx_yyzz[j] = ovlvecfunc::fvec_xxxx_yyzz_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_yy[j], pb_yyzz[j], pb_zz[j], s_0_0[j]);

                t_xxxx_yzzz[j] = ovlvecfunc::fvec_xxxx_yzzz_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_yz[j], pb_yzzz[j], s_0_0[j]);

                t_xxxx_zzzz[j] = ovlvecfunc::fvec_xxxx_zzzz_s_0(fx[j], pa_xx[j], pa_xxxx[j], pb_zz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (3) = (15,20)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxxy_xxxx, t_xxxy_xxxy, t_xxxy_xxxz, \
                                     t_xxxy_xxyy, t_xxxy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xxxx[j] = ovlvecfunc::fvec_xxxy_xxxx_s_0(fx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xxxy_xxxy[j] = ovlvecfunc::fvec_xxxy_xxxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xxxy_xxxz[j] = ovlvecfunc::fvec_xxxy_xxxz_s_0(fx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xxxy_xxyy[j] = ovlvecfunc::fvec_xxxy_xxyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_xxxy_xxyz[j] = ovlvecfunc::fvec_xxxy_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xxyz[j], pb_xxz[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (4) = (20,25)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxy, pa_xxy, pa_xy, pa_y, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxxy_xxzz, t_xxxy_xyyy, \
                                     t_xxxy_xyyz, t_xxxy_xyzz, t_xxxy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_xxzz[j] = ovlvecfunc::fvec_xxxy_xxzz_s_0(fx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], s_0_0[j]);

                t_xxxy_xyyy[j] = ovlvecfunc::fvec_xxxy_xyyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]);

                t_xxxy_xyyz[j] = ovlvecfunc::fvec_xxxy_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xxxy_xyzz[j] = ovlvecfunc::fvec_xxxy_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]);

                t_xxxy_xzzz[j] = ovlvecfunc::fvec_xxxy_xzzz_s_0(fx[j], pa_xxxy[j], pa_xxy[j], pa_xy[j], pa_y[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (5) = (25,30)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxy, pa_xy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xxxy_yyyy, \
                                     t_xxxy_yyyz, t_xxxy_yyzz, t_xxxy_yzzz, t_xxxy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxy_yyyy[j] = ovlvecfunc::fvec_xxxy_yyyy_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_xxxy_yyyz[j] = ovlvecfunc::fvec_xxxy_yyyz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xxxy_yyzz[j] = ovlvecfunc::fvec_xxxy_yyzz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], s_0_0[j]);

                t_xxxy_yzzz[j] = ovlvecfunc::fvec_xxxy_yzzz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxy[j], pa_xy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]);

                t_xxxy_zzzz[j] = ovlvecfunc::fvec_xxxy_zzzz_s_0(fx[j], pa_xxxy[j], pa_xy[j], pb_zz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (6) = (30,35)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxxz_xxxx, t_xxxz_xxxy, t_xxxz_xxxz, \
                                     t_xxxz_xxyy, t_xxxz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xxxx[j] = ovlvecfunc::fvec_xxxz_xxxx_s_0(fx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xxxz_xxxy[j] = ovlvecfunc::fvec_xxxz_xxxy_s_0(fx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xxxz_xxxz[j] = ovlvecfunc::fvec_xxxz_xxxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xxxz_xxyy[j] = ovlvecfunc::fvec_xxxz_xxyy_s_0(fx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], s_0_0[j]);

                t_xxxz_xxyz[j] = ovlvecfunc::fvec_xxxz_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xxy[j], pb_xxyz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]);
            }

            // Batch of Integrals (7) = (35,40)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxx, pa_xxxz, pa_xxz, pa_xz, pa_z, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxxz_xxzz, \
                                     t_xxxz_xyyy, t_xxxz_xyyz, t_xxxz_xyzz, t_xxxz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_xxzz[j] = ovlvecfunc::fvec_xxxz_xxzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xxxz_xyyy[j] = ovlvecfunc::fvec_xxxz_xyyy_s_0(fx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], s_0_0[j]);

                t_xxxz_xyyz[j] = ovlvecfunc::fvec_xxxz_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]);

                t_xxxz_xyzz[j] = ovlvecfunc::fvec_xxxz_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]);

                t_xxxz_xzzz[j] = ovlvecfunc::fvec_xxxz_xzzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxx[j], pa_xxxz[j], pa_xxz[j], pa_xz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (8) = (40,45)

            #pragma omp simd aligned(fx, pa_x, pa_xxx, pa_xxxz, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xxxz_yyyy, \
                                     t_xxxz_yyyz, t_xxxz_yyzz, t_xxxz_yzzz, t_xxxz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxxz_yyyy[j] = ovlvecfunc::fvec_xxxz_yyyy_s_0(fx[j], pa_xxxz[j], pa_xz[j], pb_yy[j], pb_yyyy[j], s_0_0[j]);

                t_xxxz_yyyz[j] = ovlvecfunc::fvec_xxxz_yyyz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], s_0_0[j]);

                t_xxxz_yyzz[j] = ovlvecfunc::fvec_xxxz_yyzz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xxxz_yzzz[j] = ovlvecfunc::fvec_xxxz_yzzz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], s_0_0[j]);

                t_xxxz_zzzz[j] = ovlvecfunc::fvec_xxxz_zzzz_s_0(fx[j], pa_x[j], pa_xxx[j], pa_xxxz[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (9) = (45,50)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxyy_xxxx, t_xxyy_xxxy, t_xxyy_xxxz, \
                                     t_xxyy_xxyy, t_xxyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xxxx[j] = ovlvecfunc::fvec_xxyy_xxxx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xxyy_xxxy[j] = ovlvecfunc::fvec_xxyy_xxxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xxyy_xxxz[j] = ovlvecfunc::fvec_xxyy_xxxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xxyy_xxyy[j] = ovlvecfunc::fvec_xxyy_xxyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_xxyy_xxyz[j] = ovlvecfunc::fvec_xxyy_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xxyz[j], pb_xxz[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (10) = (50,55)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyy, pa_xy, pa_xyy, pa_y, pa_yy, pb_x, pb_xx, pb_xxzz, \
                                     pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, \
                                     pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxyy_xxzz, t_xxyy_xyyy, \
                                     t_xxyy_xyyz, t_xxyy_xyzz, t_xxyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_xxzz[j] = ovlvecfunc::fvec_xxyy_xxzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], s_0_0[j]);

                t_xxyy_xyyy[j] = ovlvecfunc::fvec_xxyy_xyyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]);

                t_xxyy_xyyz[j] = ovlvecfunc::fvec_xxyy_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xxyy_xyzz[j] = ovlvecfunc::fvec_xxyy_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_xy[j], pa_xyy[j], pa_y[j], pa_yy[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]);

                t_xxyy_xzzz[j] = ovlvecfunc::fvec_xxyy_xzzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxyy[j], pa_xyy[j], pa_yy[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (11) = (55,60)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyy, pa_y, pa_yy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, \
                                     t_xxyy_yyyy, t_xxyy_yyyz, t_xxyy_yyzz, t_xxyy_yzzz, t_xxyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyy_yyyy[j] = ovlvecfunc::fvec_xxyy_yyyy_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_xxyy_yyyz[j] = ovlvecfunc::fvec_xxyy_yyyz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xxyy_yyzz[j] = ovlvecfunc::fvec_xxyy_yyzz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], s_0_0[j]);

                t_xxyy_yzzz[j] = ovlvecfunc::fvec_xxyy_yzzz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyy[j], pa_y[j], pa_yy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]);

                t_xxyy_zzzz[j] = ovlvecfunc::fvec_xxyy_zzzz_s_0(fx[j], pa_xx[j], pa_xxyy[j], pa_yy[j], pb_zz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (12) = (60,65)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, \
                                     pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, \
                                     pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxyz_xxxx, \
                                     t_xxyz_xxxy, t_xxyz_xxxz, t_xxyz_xxyy, t_xxyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xxxx[j] = ovlvecfunc::fvec_xxyz_xxxx_s_0(fx[j], pa_xxyz[j], pa_xyz[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xxyz_xxxy[j] = ovlvecfunc::fvec_xxyz_xxxy_s_0(fx[j], pa_xxyz[j], pa_xxz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xxyz_xxxz[j] = ovlvecfunc::fvec_xxyz_xxxz_s_0(fx[j], pa_xxy[j], pa_xxyz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xxyz_xxyy[j] = ovlvecfunc::fvec_xxyz_xxyy_s_0(fx[j], pa_xxyz[j], pa_xxz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_xxyz_xxyz[j] = ovlvecfunc::fvec_xxyz_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (13) = (65,70)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_xy, pa_xyz, pa_xz, pa_y, pa_yz, pa_z, \
                                     pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, \
                                     pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     s_0_0, t_xxyz_xxzz, t_xxyz_xyyy, t_xxyz_xyyz, t_xxyz_xyzz, t_xxyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_xxzz[j] = ovlvecfunc::fvec_xxyz_xxzz_s_0(fx[j], pa_xxy[j], pa_xxyz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xxyz_xyyy[j] = ovlvecfunc::fvec_xxyz_xyyy_s_0(fx[j], pa_xxyz[j], pa_xxz[j], pa_xyz[j], pa_xz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]);

                t_xxyz_xyyz[j] = ovlvecfunc::fvec_xxyz_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xxyz_xyzz[j] = ovlvecfunc::fvec_xxyz_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_xy[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xxyz_xzzz[j] = ovlvecfunc::fvec_xxyz_xzzz_s_0(fx[j], pa_xxy[j], pa_xxyz[j], pa_xy[j], pa_xyz[j], pa_y[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (14) = (70,75)

            #pragma omp simd aligned(fx, pa_xx, pa_xxy, pa_xxyz, pa_xxz, pa_y, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     s_0_0, t_xxyz_yyyy, t_xxyz_yyyz, t_xxyz_yyzz, t_xxyz_yzzz, t_xxyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxyz_yyyy[j] = ovlvecfunc::fvec_xxyz_yyyy_s_0(fx[j], pa_xxyz[j], pa_xxz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_xxyz_yyyz[j] = ovlvecfunc::fvec_xxyz_yyyz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xxyz_yyzz[j] = ovlvecfunc::fvec_xxyz_yyzz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xxyz_yzzz[j] = ovlvecfunc::fvec_xxyz_yzzz_s_0(fx[j], pa_xx[j], pa_xxy[j], pa_xxyz[j], pa_xxz[j], pa_y[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);

                t_xxyz_zzzz[j] = ovlvecfunc::fvec_xxyz_zzzz_s_0(fx[j], pa_xxy[j], pa_xxyz[j], pa_y[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (15) = (75,80)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xxzz_xxxx, t_xxzz_xxxy, t_xxzz_xxxz, \
                                     t_xxzz_xxyy, t_xxzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xxxx[j] = ovlvecfunc::fvec_xxzz_xxxx_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xxzz_xxxy[j] = ovlvecfunc::fvec_xxzz_xxxy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xxzz_xxxz[j] = ovlvecfunc::fvec_xxzz_xxxz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xxzz_xxyy[j] = ovlvecfunc::fvec_xxzz_xxyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], s_0_0[j]);

                t_xxzz_xxyz[j] = ovlvecfunc::fvec_xxzz_xxyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xxy[j], pb_xxyz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]);
            }

            // Batch of Integrals (16) = (80,85)

            #pragma omp simd aligned(fx, pa_x, pa_xx, pa_xxz, pa_xxzz, pa_xz, pa_xzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xxzz_xxzz, \
                                     t_xxzz_xyyy, t_xxzz_xyyz, t_xxzz_xyzz, t_xxzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_xxzz[j] = ovlvecfunc::fvec_xxzz_xxzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xxzz_xyyy[j] = ovlvecfunc::fvec_xxzz_xyyy_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxzz[j], pa_xzz[j], pa_zz[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], s_0_0[j]);

                t_xxzz_xyyz[j] = ovlvecfunc::fvec_xxzz_xyyz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]);

                t_xxzz_xyzz[j] = ovlvecfunc::fvec_xxzz_xyzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]);

                t_xxzz_xzzz[j] = ovlvecfunc::fvec_xxzz_xzzz_s_0(fx[j], pa_x[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_xz[j], pa_xzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (17) = (85,90)

            #pragma omp simd aligned(fx, pa_xx, pa_xxz, pa_xxzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, \
                                     pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, \
                                     t_xxzz_yyyy, t_xxzz_yyyz, t_xxzz_yyzz, t_xxzz_yzzz, t_xxzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xxzz_yyyy[j] = ovlvecfunc::fvec_xxzz_yyyy_s_0(fx[j], pa_xx[j], pa_xxzz[j], pa_zz[j], pb_yy[j], pb_yyyy[j], s_0_0[j]);

                t_xxzz_yyyz[j] = ovlvecfunc::fvec_xxzz_yyyz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], s_0_0[j]);

                t_xxzz_yyzz[j] = ovlvecfunc::fvec_xxzz_yyzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xxzz_yzzz[j] = ovlvecfunc::fvec_xxzz_yzzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], s_0_0[j]);

                t_xxzz_zzzz[j] = ovlvecfunc::fvec_xxzz_zzzz_s_0(fx[j], pa_xx[j], pa_xxz[j], pa_xxzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (18) = (90,95)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xyyy_xxxx, t_xyyy_xxxy, t_xyyy_xxxz, \
                                     t_xyyy_xxyy, t_xyyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xxxx[j] = ovlvecfunc::fvec_xyyy_xxxx_s_0(fx[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xyyy_xxxy[j] = ovlvecfunc::fvec_xyyy_xxxy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xyyy_xxxz[j] = ovlvecfunc::fvec_xyyy_xxxz_s_0(fx[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xyyy_xxyy[j] = ovlvecfunc::fvec_xyyy_xxyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_xyyy_xxyz[j] = ovlvecfunc::fvec_xyyy_xxyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xxyz[j], pb_xxz[j], pb_xyz[j], pb_xz[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (19) = (95,100)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pa_y, pa_yy, pa_yyy, pb_x, pb_xx, pb_xxzz, pb_xy, \
                                     pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xyyy_xxzz, t_xyyy_xyyy, \
                                     t_xyyy_xyyz, t_xyyy_xyzz, t_xyyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_xxzz[j] = ovlvecfunc::fvec_xyyy_xxzz_s_0(fx[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_x[j], pb_xx[j], pb_xxzz[j], pb_xzz[j], pb_zz[j], s_0_0[j]);

                t_xyyy_xyyy[j] = ovlvecfunc::fvec_xyyy_xyyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]);

                t_xyyy_xyyz[j] = ovlvecfunc::fvec_xyyy_xyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xyyy_xyzz[j] = ovlvecfunc::fvec_xyyy_xyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pa_y[j], pa_yy[j], pa_yyy[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], pb_y[j], pb_yzz[j], pb_zz[j], s_0_0[j]);

                t_xyyy_xzzz[j] = ovlvecfunc::fvec_xyyy_xzzz_s_0(fx[j], pa_xy[j], pa_xyyy[j], pa_y[j], pa_yyy[j], pb_xz[j], pb_xzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (20) = (100,105)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xyyy_yyyy, \
                                     t_xyyy_yyyz, t_xyyy_yyzz, t_xyyy_yzzz, t_xyyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyy_yyyy[j] = ovlvecfunc::fvec_xyyy_yyyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_xyyy_yyyz[j] = ovlvecfunc::fvec_xyyy_yyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xyyy_yyzz[j] = ovlvecfunc::fvec_xyyy_yyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], s_0_0[j]);

                t_xyyy_yzzz[j] = ovlvecfunc::fvec_xyyy_yzzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]);

                t_xyyy_zzzz[j] = ovlvecfunc::fvec_xyyy_zzzz_s_0(fx[j], pa_xy[j], pa_xyyy[j], pb_zz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (21) = (105,110)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, \
                                     pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, \
                                     pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xyyz_xxxx, \
                                     t_xyyz_xxxy, t_xyyz_xxxz, t_xyyz_xxyy, t_xyyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xxxx[j] = ovlvecfunc::fvec_xyyz_xxxx_s_0(fx[j], pa_xyyz[j], pa_xz[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xyyz_xxxy[j] = ovlvecfunc::fvec_xyyz_xxxy_s_0(fx[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xyyz_xxxz[j] = ovlvecfunc::fvec_xyyz_xxxz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xyyz_xxyy[j] = ovlvecfunc::fvec_xyyz_xxyy_s_0(fx[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_xyyz_xxyz[j] = ovlvecfunc::fvec_xyyz_xxyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (22) = (110,115)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pa_y, pa_yy, pa_yyz, pa_yz, pa_z, \
                                     pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, \
                                     pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     s_0_0, t_xyyz_xxzz, t_xyyz_xyyy, t_xyyz_xyyz, t_xyyz_xyzz, t_xyyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_xxzz[j] = ovlvecfunc::fvec_xyyz_xxzz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xyyz_xyyy[j] = ovlvecfunc::fvec_xyyz_xyyy_s_0(fx[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]);

                t_xyyz_xyyz[j] = ovlvecfunc::fvec_xyyz_xyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xyyz_xyzz[j] = ovlvecfunc::fvec_xyyz_xyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xyyz_xzzz[j] = ovlvecfunc::fvec_xyyz_xzzz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pa_yy[j], pa_yyz[j], pa_z[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (23) = (115,120)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyy, pa_xyyz, pa_xyz, pa_xz, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, \
                                     t_xyyz_yyyy, t_xyyz_yyyz, t_xyyz_yyzz, t_xyyz_yzzz, t_xyyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyyz_yyyy[j] = ovlvecfunc::fvec_xyyz_yyyy_s_0(fx[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_xyyz_yyyz[j] = ovlvecfunc::fvec_xyyz_yyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xyyz_yyzz[j] = ovlvecfunc::fvec_xyyz_yyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xyyz_yzzz[j] = ovlvecfunc::fvec_xyyz_yzzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyy[j], pa_xyyz[j], pa_xyz[j], pa_xz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);

                t_xyyz_zzzz[j] = ovlvecfunc::fvec_xyyz_zzzz_s_0(fx[j], pa_x[j], pa_xyy[j], pa_xyyz[j], pa_xz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (24) = (120,125)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, \
                                     pb_xy, pb_xyy, pb_xyz, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xyzz_xxxx, \
                                     t_xyzz_xxxy, t_xyzz_xxxz, t_xyzz_xxyy, t_xyzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xxxx[j] = ovlvecfunc::fvec_xyzz_xxxx_s_0(fx[j], pa_xy[j], pa_xyzz[j], pa_y[j], pa_yzz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xyzz_xxxy[j] = ovlvecfunc::fvec_xyzz_xxxy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xyzz_xxxz[j] = ovlvecfunc::fvec_xyzz_xxxz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xyzz_xxyy[j] = ovlvecfunc::fvec_xyzz_xxyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_xy[j], pb_xyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_xyzz_xxyz[j] = ovlvecfunc::fvec_xyzz_xxyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_xy[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (25) = (125,130)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pa_y, pa_yz, pa_yzz, pa_z, pa_zz, \
                                     pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, \
                                     pb_xz, pb_xzz, pb_xzzz, pb_y, pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, \
                                     s_0_0, t_xyzz_xxzz, t_xyzz_xyyy, t_xyzz_xyyz, t_xyzz_xyzz, t_xyzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_xxzz[j] = ovlvecfunc::fvec_xyzz_xxzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xyzz_xyyy[j] = ovlvecfunc::fvec_xyzz_xyyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pa_y[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], s_0_0[j]);

                t_xyzz_xyyz[j] = ovlvecfunc::fvec_xyzz_xyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xyzz_xyzz[j] = ovlvecfunc::fvec_xyzz_xyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xyzz_xzzz[j] = ovlvecfunc::fvec_xyzz_xzzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_y[j], pa_yz[j], pa_yzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (26) = (130,135)

            #pragma omp simd aligned(fx, pa_x, pa_xy, pa_xyz, pa_xyzz, pa_xz, pa_xzz, pb_y, pb_yy, pb_yyy, pb_yyyy, \
                                     pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, \
                                     t_xyzz_yyyy, t_xyzz_yyyz, t_xyzz_yyzz, t_xyzz_yzzz, t_xyzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xyzz_yyyy[j] = ovlvecfunc::fvec_xyzz_yyyy_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyzz[j], pa_xzz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_xyzz_yyyz[j] = ovlvecfunc::fvec_xyzz_yyyz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_xyzz_yyzz[j] = ovlvecfunc::fvec_xyzz_yyzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xyzz_yzzz[j] = ovlvecfunc::fvec_xyzz_yzzz_s_0(fx[j], pa_x[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pa_xz[j], pa_xzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);

                t_xyzz_zzzz[j] = ovlvecfunc::fvec_xyzz_zzzz_s_0(fx[j], pa_xy[j], pa_xyz[j], pa_xyzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (27) = (135,140)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xyy, pb_xyz, \
                                     pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, t_xzzz_xxxx, t_xzzz_xxxy, t_xzzz_xxxz, \
                                     t_xzzz_xxyy, t_xzzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xxxx[j] = ovlvecfunc::fvec_xzzz_xxxx_s_0(fx[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxx[j], s_0_0[j]);

                t_xzzz_xxxy[j] = ovlvecfunc::fvec_xzzz_xxxy_s_0(fx[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_xxxy[j], pb_xxy[j], pb_xy[j], pb_y[j], s_0_0[j]);

                t_xzzz_xxxz[j] = ovlvecfunc::fvec_xzzz_xxxz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xx[j], pb_xxx[j], pb_xxxz[j], pb_xxz[j], pb_xz[j], pb_z[j], s_0_0[j]);

                t_xzzz_xxyy[j] = ovlvecfunc::fvec_xzzz_xxyy_s_0(fx[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xx[j], pb_xxyy[j], pb_xyy[j], pb_yy[j], s_0_0[j]);

                t_xzzz_xxyz[j] = ovlvecfunc::fvec_xzzz_xxyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xxy[j], pb_xxyz[j], pb_xy[j], pb_xyz[j], pb_y[j], pb_yz[j], s_0_0[j]);
            }

            // Batch of Integrals (28) = (140,145)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_y, \
                                     pb_yy, pb_yyy, pb_yyz, pb_yz, pb_yzz, pb_z, pb_zz, pb_zzz, s_0_0, t_xzzz_xxzz, \
                                     t_xzzz_xyyy, t_xzzz_xyyz, t_xzzz_xyzz, t_xzzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_xxzz[j] = ovlvecfunc::fvec_xzzz_xxzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_xz[j], pb_xzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xzzz_xyyy[j] = ovlvecfunc::fvec_xzzz_xyyy_s_0(fx[j], pa_xz[j], pa_xzzz[j], pa_z[j], pa_zzz[j], pb_xy[j], pb_xyyy[j], pb_y[j], pb_yyy[j], s_0_0[j]);

                t_xzzz_xyyz[j] = ovlvecfunc::fvec_xzzz_xyyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], pb_yy[j], pb_yyz[j], pb_z[j], s_0_0[j]);

                t_xzzz_xyzz[j] = ovlvecfunc::fvec_xzzz_xyzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_y[j], pb_yz[j], pb_yzz[j], s_0_0[j]);

                t_xzzz_xzzz[j] = ovlvecfunc::fvec_xzzz_xzzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);
            }

            // Batch of Integrals (29) = (145,150)

            #pragma omp simd aligned(fx, pa_x, pa_xz, pa_xzz, pa_xzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_xzzz_yyyy, \
                                     t_xzzz_yyyz, t_xzzz_yyzz, t_xzzz_yzzz, t_xzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_xzzz_yyyy[j] = ovlvecfunc::fvec_xzzz_yyyy_s_0(fx[j], pa_xz[j], pa_xzzz[j], pb_yy[j], pb_yyyy[j], s_0_0[j]);

                t_xzzz_yyyz[j] = ovlvecfunc::fvec_xzzz_yyyz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], s_0_0[j]);

                t_xzzz_yyzz[j] = ovlvecfunc::fvec_xzzz_yyzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_xzzz_yzzz[j] = ovlvecfunc::fvec_xzzz_yzzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], s_0_0[j]);

                t_xzzz_zzzz[j] = ovlvecfunc::fvec_xzzz_zzzz_s_0(fx[j], pa_x[j], pa_xz[j], pa_xzz[j], pa_xzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (30) = (150,155)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, pb_z, s_0_0, \
                                     t_yyyy_xxxx, t_yyyy_xxxy, t_yyyy_xxxz, t_yyyy_xxyy, t_yyyy_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xxxx[j] = ovlvecfunc::fvec_yyyy_xxxx_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_xx[j], pb_xxxx[j], s_0_0[j]);

                t_yyyy_xxxy[j] = ovlvecfunc::fvec_yyyy_xxxy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], s_0_0[j]);

                t_yyyy_xxxz[j] = ovlvecfunc::fvec_yyyy_xxxz_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_xxxz[j], pb_xz[j], s_0_0[j]);

                t_yyyy_xxyy[j] = ovlvecfunc::fvec_yyyy_xxyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_yyyy_xxyz[j] = ovlvecfunc::fvec_yyyy_xxyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_xxyz[j], pb_xxz[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (31) = (155,160)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_x, pb_xx, pb_xxzz, pb_xy, pb_xyy, pb_xyyy, \
                                     pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_zz, s_0_0, t_yyyy_xxzz, \
                                     t_yyyy_xyyy, t_yyyy_xyyz, t_yyyy_xyzz, t_yyyy_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_xxzz[j] = ovlvecfunc::fvec_yyyy_xxzz_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_xx[j], pb_xxzz[j], pb_zz[j], s_0_0[j]);

                t_yyyy_xyyy[j] = ovlvecfunc::fvec_yyyy_xyyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], s_0_0[j]);

                t_yyyy_xyyz[j] = ovlvecfunc::fvec_yyyy_xyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], s_0_0[j]);

                t_yyyy_xyzz[j] = ovlvecfunc::fvec_yyyy_xyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_x[j], pb_xy[j], pb_xyzz[j], pb_xzz[j], s_0_0[j]);

                t_yyyy_xzzz[j] = ovlvecfunc::fvec_yyyy_xzzz_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_xz[j], pb_xzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (32) = (160,165)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyy, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_yyyy_yyyy, \
                                     t_yyyy_yyyz, t_yyyy_yyzz, t_yyyy_yzzz, t_yyyy_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyy_yyyy[j] = ovlvecfunc::fvec_yyyy_yyyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_yyyy_yyyz[j] = ovlvecfunc::fvec_yyyy_yyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_yyyy_yyzz[j] = ovlvecfunc::fvec_yyyy_yyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_y[j], pb_yy[j], pb_yyzz[j], pb_yzz[j], pb_zz[j], s_0_0[j]);

                t_yyyy_yzzz[j] = ovlvecfunc::fvec_yyyy_yzzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyy[j], pb_yz[j], pb_yzzz[j], pb_z[j], pb_zzz[j], s_0_0[j]);

                t_yyyy_zzzz[j] = ovlvecfunc::fvec_yyyy_zzzz_s_0(fx[j], pa_yy[j], pa_yyyy[j], pb_zz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (33) = (165,170)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, s_0_0, t_yyyz_xxxx, t_yyyz_xxxy, t_yyyz_xxxz, t_yyyz_xxyy, \
                                     t_yyyz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xxxx[j] = ovlvecfunc::fvec_yyyz_xxxx_s_0(fx[j], pa_yyyz[j], pa_yz[j], pb_xx[j], pb_xxxx[j], s_0_0[j]);

                t_yyyz_xxxy[j] = ovlvecfunc::fvec_yyyz_xxxy_s_0(fx[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], s_0_0[j]);

                t_yyyz_xxxz[j] = ovlvecfunc::fvec_yyyz_xxxz_s_0(fx[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], s_0_0[j]);

                t_yyyz_xxyy[j] = ovlvecfunc::fvec_yyyz_xxyy_s_0(fx[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_yyyz_xxyz[j] = ovlvecfunc::fvec_yyyz_xxyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (34) = (170,175)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, \
                                     pb_zz, s_0_0, t_yyyz_xxzz, t_yyyz_xyyy, t_yyyz_xyyz, t_yyyz_xyzz, t_yyyz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_xxzz[j] = ovlvecfunc::fvec_yyyz_xxzz_s_0(fx[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_yyyz_xyyy[j] = ovlvecfunc::fvec_yyyz_xyyy_s_0(fx[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], s_0_0[j]);

                t_yyyz_xyyz[j] = ovlvecfunc::fvec_yyyz_xyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], s_0_0[j]);

                t_yyyz_xyzz[j] = ovlvecfunc::fvec_yyyz_xyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], s_0_0[j]);

                t_yyyz_xzzz[j] = ovlvecfunc::fvec_yyyz_xzzz_s_0(fx[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (35) = (175,180)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyy, pa_yyyz, pa_yyz, pa_yz, pa_z, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     s_0_0, t_yyyz_yyyy, t_yyyz_yyyz, t_yyyz_yyzz, t_yyyz_yzzz, t_yyyz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyyz_yyyy[j] = ovlvecfunc::fvec_yyyz_yyyy_s_0(fx[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_yyyz_yyyz[j] = ovlvecfunc::fvec_yyyz_yyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_yyyz_yyzz[j] = ovlvecfunc::fvec_yyyz_yyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_yyyz_yzzz[j] = ovlvecfunc::fvec_yyyz_yzzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyy[j], pa_yyyz[j], pa_yyz[j], pa_yz[j], pa_z[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);

                t_yyyz_zzzz[j] = ovlvecfunc::fvec_yyyz_zzzz_s_0(fx[j], pa_y[j], pa_yyy[j], pa_yyyz[j], pa_yz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (36) = (180,185)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, s_0_0, t_yyzz_xxxx, t_yyzz_xxxy, t_yyzz_xxxz, t_yyzz_xxyy, \
                                     t_yyzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xxxx[j] = ovlvecfunc::fvec_yyzz_xxxx_s_0(fx[j], pa_yy[j], pa_yyzz[j], pa_zz[j], pb_xx[j], pb_xxxx[j], s_0_0[j]);

                t_yyzz_xxxy[j] = ovlvecfunc::fvec_yyzz_xxxy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], s_0_0[j]);

                t_yyzz_xxxz[j] = ovlvecfunc::fvec_yyzz_xxxz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], s_0_0[j]);

                t_yyzz_xxyy[j] = ovlvecfunc::fvec_yyzz_xxyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_yyzz_xxyz[j] = ovlvecfunc::fvec_yyzz_xxyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (37) = (185,190)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, \
                                     pb_zz, s_0_0, t_yyzz_xxzz, t_yyzz_xyyy, t_yyzz_xyyz, t_yyzz_xyzz, t_yyzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_xxzz[j] = ovlvecfunc::fvec_yyzz_xxzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_yyzz_xyyy[j] = ovlvecfunc::fvec_yyzz_xyyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], s_0_0[j]);

                t_yyzz_xyyz[j] = ovlvecfunc::fvec_yyzz_xyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], s_0_0[j]);

                t_yyzz_xyzz[j] = ovlvecfunc::fvec_yyzz_xyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], s_0_0[j]);

                t_yyzz_xzzz[j] = ovlvecfunc::fvec_yyzz_xzzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (38) = (190,195)

            #pragma omp simd aligned(fx, pa_y, pa_yy, pa_yyz, pa_yyzz, pa_yz, pa_yzz, pa_z, pa_zz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     s_0_0, t_yyzz_yyyy, t_yyzz_yyyz, t_yyzz_yyzz, t_yyzz_yzzz, t_yyzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yyzz_yyyy[j] = ovlvecfunc::fvec_yyzz_yyyy_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyzz[j], pa_yzz[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_yyzz_yyyz[j] = ovlvecfunc::fvec_yyzz_yyyz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_yyzz_yyzz[j] = ovlvecfunc::fvec_yyzz_yyzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_yyzz_yzzz[j] = ovlvecfunc::fvec_yyzz_yzzz_s_0(fx[j], pa_y[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_yz[j], pa_yzz[j], pa_z[j], pa_zz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);

                t_yyzz_zzzz[j] = ovlvecfunc::fvec_yyzz_zzzz_s_0(fx[j], pa_yy[j], pa_yyz[j], pa_yyzz[j], pa_z[j], pa_zz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (39) = (195,200)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxx, \
                                     pb_xxxx, pb_xxxy, pb_xxxz, pb_xxy, pb_xxyy, pb_xxyz, pb_xxz, pb_xy, pb_xz, pb_y, pb_yy, \
                                     pb_yz, pb_z, s_0_0, t_yzzz_xxxx, t_yzzz_xxxy, t_yzzz_xxxz, t_yzzz_xxyy, \
                                     t_yzzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xxxx[j] = ovlvecfunc::fvec_yzzz_xxxx_s_0(fx[j], pa_yz[j], pa_yzzz[j], pb_xx[j], pb_xxxx[j], s_0_0[j]);

                t_yzzz_xxxy[j] = ovlvecfunc::fvec_yzzz_xxxy_s_0(fx[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xxx[j], pb_xxxy[j], pb_xy[j], s_0_0[j]);

                t_yzzz_xxxz[j] = ovlvecfunc::fvec_yzzz_xxxz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], s_0_0[j]);

                t_yzzz_xxyy[j] = ovlvecfunc::fvec_yzzz_xxyy_s_0(fx[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_xx[j], pb_xxy[j], pb_xxyy[j], pb_y[j], pb_yy[j], s_0_0[j]);

                t_yzzz_xxyz[j] = ovlvecfunc::fvec_yzzz_xxyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_xx[j], pb_xxy[j], pb_xxyz[j], pb_xxz[j], pb_y[j], pb_yz[j], pb_z[j], s_0_0[j]);
            }

            // Batch of Integrals (40) = (200,205)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_x, pb_xx, pb_xxz, \
                                     pb_xxzz, pb_xy, pb_xyy, pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, \
                                     pb_zz, s_0_0, t_yzzz_xxzz, t_yzzz_xyyy, t_yzzz_xyyz, t_yzzz_xyzz, t_yzzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_xxzz[j] = ovlvecfunc::fvec_yzzz_xxzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_yzzz_xyyy[j] = ovlvecfunc::fvec_yzzz_xyyy_s_0(fx[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyy[j], s_0_0[j]);

                t_yzzz_xyyz[j] = ovlvecfunc::fvec_yzzz_xyyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xy[j], pb_xyy[j], pb_xyyz[j], pb_xyz[j], pb_xz[j], s_0_0[j]);

                t_yzzz_xyzz[j] = ovlvecfunc::fvec_yzzz_xyzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_x[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], pb_xz[j], pb_xzz[j], s_0_0[j]);

                t_yzzz_xzzz[j] = ovlvecfunc::fvec_yzzz_xzzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (41) = (205,210)

            #pragma omp simd aligned(fx, pa_y, pa_yz, pa_yzz, pa_yzzz, pa_z, pa_zz, pa_zzz, pb_y, pb_yy, pb_yyy, \
                                     pb_yyyy, pb_yyyz, pb_yyz, pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, \
                                     s_0_0, t_yzzz_yyyy, t_yzzz_yyyz, t_yzzz_yyzz, t_yzzz_yzzz, t_yzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_yzzz_yyyy[j] = ovlvecfunc::fvec_yzzz_yyyy_s_0(fx[j], pa_yz[j], pa_yzzz[j], pa_z[j], pa_zzz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyy[j], s_0_0[j]);

                t_yzzz_yyyz[j] = ovlvecfunc::fvec_yzzz_yyyz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yy[j], pb_yyy[j], pb_yyyz[j], pb_yyz[j], pb_yz[j], pb_z[j], s_0_0[j]);

                t_yzzz_yyzz[j] = ovlvecfunc::fvec_yzzz_yyzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_yz[j], pb_yzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_yzzz_yzzz[j] = ovlvecfunc::fvec_yzzz_yzzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pa_z[j], pa_zz[j], pa_zzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], s_0_0[j]);

                t_yzzz_zzzz[j] = ovlvecfunc::fvec_yzzz_zzzz_s_0(fx[j], pa_y[j], pa_yz[j], pa_yzz[j], pa_yzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (42) = (210,215)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xx, pb_xxx, pb_xxxx, pb_xxxy, pb_xxxz, \
                                     pb_xxy, pb_xxyy, pb_xxyz, pb_xy, pb_xz, pb_y, pb_yy, pb_yz, s_0_0, t_zzzz_xxxx, \
                                     t_zzzz_xxxy, t_zzzz_xxxz, t_zzzz_xxyy, t_zzzz_xxyz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xxxx[j] = ovlvecfunc::fvec_zzzz_xxxx_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_xx[j], pb_xxxx[j], s_0_0[j]);

                t_zzzz_xxxy[j] = ovlvecfunc::fvec_zzzz_xxxy_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_xxxy[j], pb_xy[j], s_0_0[j]);

                t_zzzz_xxxz[j] = ovlvecfunc::fvec_zzzz_xxxz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_x[j], pb_xxx[j], pb_xxxz[j], pb_xz[j], s_0_0[j]);

                t_zzzz_xxyy[j] = ovlvecfunc::fvec_zzzz_xxyy_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_xx[j], pb_xxyy[j], pb_yy[j], s_0_0[j]);

                t_zzzz_xxyz[j] = ovlvecfunc::fvec_zzzz_xxyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_xxy[j], pb_xxyz[j], pb_y[j], pb_yz[j], s_0_0[j]);
            }

            // Batch of Integrals (43) = (215,220)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_x, pb_xx, pb_xxz, pb_xxzz, pb_xy, pb_xyy, \
                                     pb_xyyy, pb_xyyz, pb_xyz, pb_xyzz, pb_xz, pb_xzz, pb_xzzz, pb_z, pb_zz, s_0_0, \
                                     t_zzzz_xxzz, t_zzzz_xyyy, t_zzzz_xyyz, t_zzzz_xyzz, t_zzzz_xzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_xxzz[j] = ovlvecfunc::fvec_zzzz_xxzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_xx[j], pb_xxz[j], pb_xxzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_zzzz_xyyy[j] = ovlvecfunc::fvec_zzzz_xyyy_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_xy[j], pb_xyyy[j], s_0_0[j]);

                t_zzzz_xyyz[j] = ovlvecfunc::fvec_zzzz_xyyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_x[j], pb_xyy[j], pb_xyyz[j], pb_xz[j], s_0_0[j]);

                t_zzzz_xyzz[j] = ovlvecfunc::fvec_zzzz_xyzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_xy[j], pb_xyz[j], pb_xyzz[j], s_0_0[j]);

                t_zzzz_xzzz[j] = ovlvecfunc::fvec_zzzz_xzzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_x[j], pb_xz[j], pb_xzz[j], pb_xzzz[j], s_0_0[j]);
            }

            // Batch of Integrals (44) = (220,225)

            #pragma omp simd aligned(fx, pa_z, pa_zz, pa_zzz, pa_zzzz, pb_y, pb_yy, pb_yyy, pb_yyyy, pb_yyyz, pb_yyz, \
                                     pb_yyzz, pb_yz, pb_yzz, pb_yzzz, pb_z, pb_zz, pb_zzz, pb_zzzz, s_0_0, t_zzzz_yyyy, \
                                     t_zzzz_yyyz, t_zzzz_yyzz, t_zzzz_yzzz, t_zzzz_zzzz: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                t_zzzz_yyyy[j] = ovlvecfunc::fvec_zzzz_yyyy_s_0(fx[j], pa_zz[j], pa_zzzz[j], pb_yy[j], pb_yyyy[j], s_0_0[j]);

                t_zzzz_yyyz[j] = ovlvecfunc::fvec_zzzz_yyyz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_y[j], pb_yyy[j], pb_yyyz[j], pb_yz[j], s_0_0[j]);

                t_zzzz_yyzz[j] = ovlvecfunc::fvec_zzzz_yyzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_yy[j], pb_yyz[j], pb_yyzz[j], pb_z[j], pb_zz[j], s_0_0[j]);

                t_zzzz_yzzz[j] = ovlvecfunc::fvec_zzzz_yzzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_y[j], pb_yz[j], pb_yzz[j], pb_yzzz[j], s_0_0[j]);

                t_zzzz_zzzz[j] = ovlvecfunc::fvec_zzzz_zzzz_s_0(fx[j], pa_z[j], pa_zz[j], pa_zzz[j], pa_zzzz[j], pb_z[j], pb_zz[j], pb_zzz[j], pb_zzzz[j], s_0_0[j]);
            }

            idx++;
        }
    }


} // ovlrecfunc namespace

