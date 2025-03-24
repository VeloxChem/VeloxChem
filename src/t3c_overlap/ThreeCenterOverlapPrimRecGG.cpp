#include "ThreeCenterOverlapPrimRecGG.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_gg(CSimdArray<double>& pbuffer, 
                     const size_t idx_gg,
                     const size_t idx_dg,
                     const size_t idx_ff,
                     const size_t idx_fg,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

    // Set up components of auxiliary buffer : DG

    auto ts_xx_xxxx = pbuffer.data(idx_dg);

    auto ts_xx_xxxy = pbuffer.data(idx_dg + 1);

    auto ts_xx_xxxz = pbuffer.data(idx_dg + 2);

    auto ts_xx_xxyy = pbuffer.data(idx_dg + 3);

    auto ts_xx_xxyz = pbuffer.data(idx_dg + 4);

    auto ts_xx_xxzz = pbuffer.data(idx_dg + 5);

    auto ts_xx_xyyy = pbuffer.data(idx_dg + 6);

    auto ts_xx_xyyz = pbuffer.data(idx_dg + 7);

    auto ts_xx_xyzz = pbuffer.data(idx_dg + 8);

    auto ts_xx_xzzz = pbuffer.data(idx_dg + 9);

    auto ts_xx_yyyy = pbuffer.data(idx_dg + 10);

    auto ts_xx_yyyz = pbuffer.data(idx_dg + 11);

    auto ts_xx_yyzz = pbuffer.data(idx_dg + 12);

    auto ts_xx_yzzz = pbuffer.data(idx_dg + 13);

    auto ts_xx_zzzz = pbuffer.data(idx_dg + 14);

    auto ts_xy_yyyy = pbuffer.data(idx_dg + 25);

    auto ts_xy_yyyz = pbuffer.data(idx_dg + 26);

    auto ts_xy_yyzz = pbuffer.data(idx_dg + 27);

    auto ts_xy_yzzz = pbuffer.data(idx_dg + 28);

    auto ts_xz_yyyz = pbuffer.data(idx_dg + 41);

    auto ts_xz_yyzz = pbuffer.data(idx_dg + 42);

    auto ts_xz_yzzz = pbuffer.data(idx_dg + 43);

    auto ts_xz_zzzz = pbuffer.data(idx_dg + 44);

    auto ts_yy_xxxx = pbuffer.data(idx_dg + 45);

    auto ts_yy_xxxy = pbuffer.data(idx_dg + 46);

    auto ts_yy_xxxz = pbuffer.data(idx_dg + 47);

    auto ts_yy_xxyy = pbuffer.data(idx_dg + 48);

    auto ts_yy_xxyz = pbuffer.data(idx_dg + 49);

    auto ts_yy_xxzz = pbuffer.data(idx_dg + 50);

    auto ts_yy_xyyy = pbuffer.data(idx_dg + 51);

    auto ts_yy_xyyz = pbuffer.data(idx_dg + 52);

    auto ts_yy_xyzz = pbuffer.data(idx_dg + 53);

    auto ts_yy_xzzz = pbuffer.data(idx_dg + 54);

    auto ts_yy_yyyy = pbuffer.data(idx_dg + 55);

    auto ts_yy_yyyz = pbuffer.data(idx_dg + 56);

    auto ts_yy_yyzz = pbuffer.data(idx_dg + 57);

    auto ts_yy_yzzz = pbuffer.data(idx_dg + 58);

    auto ts_yy_zzzz = pbuffer.data(idx_dg + 59);

    auto ts_yz_xxxz = pbuffer.data(idx_dg + 62);

    auto ts_yz_xxzz = pbuffer.data(idx_dg + 65);

    auto ts_yz_xzzz = pbuffer.data(idx_dg + 69);

    auto ts_yz_yyyz = pbuffer.data(idx_dg + 71);

    auto ts_yz_yyzz = pbuffer.data(idx_dg + 72);

    auto ts_yz_yzzz = pbuffer.data(idx_dg + 73);

    auto ts_yz_zzzz = pbuffer.data(idx_dg + 74);

    auto ts_zz_xxxx = pbuffer.data(idx_dg + 75);

    auto ts_zz_xxxy = pbuffer.data(idx_dg + 76);

    auto ts_zz_xxxz = pbuffer.data(idx_dg + 77);

    auto ts_zz_xxyy = pbuffer.data(idx_dg + 78);

    auto ts_zz_xxyz = pbuffer.data(idx_dg + 79);

    auto ts_zz_xxzz = pbuffer.data(idx_dg + 80);

    auto ts_zz_xyyy = pbuffer.data(idx_dg + 81);

    auto ts_zz_xyyz = pbuffer.data(idx_dg + 82);

    auto ts_zz_xyzz = pbuffer.data(idx_dg + 83);

    auto ts_zz_xzzz = pbuffer.data(idx_dg + 84);

    auto ts_zz_yyyy = pbuffer.data(idx_dg + 85);

    auto ts_zz_yyyz = pbuffer.data(idx_dg + 86);

    auto ts_zz_yyzz = pbuffer.data(idx_dg + 87);

    auto ts_zz_yzzz = pbuffer.data(idx_dg + 88);

    auto ts_zz_zzzz = pbuffer.data(idx_dg + 89);

    // Set up components of auxiliary buffer : FF

    auto ts_xxx_xxx = pbuffer.data(idx_ff);

    auto ts_xxx_xxy = pbuffer.data(idx_ff + 1);

    auto ts_xxx_xxz = pbuffer.data(idx_ff + 2);

    auto ts_xxx_xyy = pbuffer.data(idx_ff + 3);

    auto ts_xxx_xyz = pbuffer.data(idx_ff + 4);

    auto ts_xxx_xzz = pbuffer.data(idx_ff + 5);

    auto ts_xxx_yyy = pbuffer.data(idx_ff + 6);

    auto ts_xxx_yyz = pbuffer.data(idx_ff + 7);

    auto ts_xxx_yzz = pbuffer.data(idx_ff + 8);

    auto ts_xxx_zzz = pbuffer.data(idx_ff + 9);

    auto ts_xxz_xxz = pbuffer.data(idx_ff + 22);

    auto ts_xxz_xyz = pbuffer.data(idx_ff + 24);

    auto ts_xxz_xzz = pbuffer.data(idx_ff + 25);

    auto ts_xyy_xxy = pbuffer.data(idx_ff + 31);

    auto ts_xyy_xyy = pbuffer.data(idx_ff + 33);

    auto ts_xyy_xyz = pbuffer.data(idx_ff + 34);

    auto ts_xyy_yyy = pbuffer.data(idx_ff + 36);

    auto ts_xyy_yyz = pbuffer.data(idx_ff + 37);

    auto ts_xyy_yzz = pbuffer.data(idx_ff + 38);

    auto ts_xzz_xxz = pbuffer.data(idx_ff + 52);

    auto ts_xzz_xyz = pbuffer.data(idx_ff + 54);

    auto ts_xzz_xzz = pbuffer.data(idx_ff + 55);

    auto ts_xzz_yyz = pbuffer.data(idx_ff + 57);

    auto ts_xzz_yzz = pbuffer.data(idx_ff + 58);

    auto ts_xzz_zzz = pbuffer.data(idx_ff + 59);

    auto ts_yyy_xxx = pbuffer.data(idx_ff + 60);

    auto ts_yyy_xxy = pbuffer.data(idx_ff + 61);

    auto ts_yyy_xxz = pbuffer.data(idx_ff + 62);

    auto ts_yyy_xyy = pbuffer.data(idx_ff + 63);

    auto ts_yyy_xyz = pbuffer.data(idx_ff + 64);

    auto ts_yyy_xzz = pbuffer.data(idx_ff + 65);

    auto ts_yyy_yyy = pbuffer.data(idx_ff + 66);

    auto ts_yyy_yyz = pbuffer.data(idx_ff + 67);

    auto ts_yyy_yzz = pbuffer.data(idx_ff + 68);

    auto ts_yyy_zzz = pbuffer.data(idx_ff + 69);

    auto ts_yyz_xxz = pbuffer.data(idx_ff + 72);

    auto ts_yyz_xyz = pbuffer.data(idx_ff + 74);

    auto ts_yyz_xzz = pbuffer.data(idx_ff + 75);

    auto ts_yyz_yyz = pbuffer.data(idx_ff + 77);

    auto ts_yyz_yzz = pbuffer.data(idx_ff + 78);

    auto ts_yyz_zzz = pbuffer.data(idx_ff + 79);

    auto ts_yzz_xxy = pbuffer.data(idx_ff + 81);

    auto ts_yzz_xxz = pbuffer.data(idx_ff + 82);

    auto ts_yzz_xyy = pbuffer.data(idx_ff + 83);

    auto ts_yzz_xyz = pbuffer.data(idx_ff + 84);

    auto ts_yzz_xzz = pbuffer.data(idx_ff + 85);

    auto ts_yzz_yyy = pbuffer.data(idx_ff + 86);

    auto ts_yzz_yyz = pbuffer.data(idx_ff + 87);

    auto ts_yzz_yzz = pbuffer.data(idx_ff + 88);

    auto ts_yzz_zzz = pbuffer.data(idx_ff + 89);

    auto ts_zzz_xxx = pbuffer.data(idx_ff + 90);

    auto ts_zzz_xxy = pbuffer.data(idx_ff + 91);

    auto ts_zzz_xxz = pbuffer.data(idx_ff + 92);

    auto ts_zzz_xyy = pbuffer.data(idx_ff + 93);

    auto ts_zzz_xyz = pbuffer.data(idx_ff + 94);

    auto ts_zzz_xzz = pbuffer.data(idx_ff + 95);

    auto ts_zzz_yyy = pbuffer.data(idx_ff + 96);

    auto ts_zzz_yyz = pbuffer.data(idx_ff + 97);

    auto ts_zzz_yzz = pbuffer.data(idx_ff + 98);

    auto ts_zzz_zzz = pbuffer.data(idx_ff + 99);

    // Set up components of auxiliary buffer : FG

    auto ts_xxx_xxxx = pbuffer.data(idx_fg);

    auto ts_xxx_xxxy = pbuffer.data(idx_fg + 1);

    auto ts_xxx_xxxz = pbuffer.data(idx_fg + 2);

    auto ts_xxx_xxyy = pbuffer.data(idx_fg + 3);

    auto ts_xxx_xxyz = pbuffer.data(idx_fg + 4);

    auto ts_xxx_xxzz = pbuffer.data(idx_fg + 5);

    auto ts_xxx_xyyy = pbuffer.data(idx_fg + 6);

    auto ts_xxx_xyyz = pbuffer.data(idx_fg + 7);

    auto ts_xxx_xyzz = pbuffer.data(idx_fg + 8);

    auto ts_xxx_xzzz = pbuffer.data(idx_fg + 9);

    auto ts_xxx_yyyy = pbuffer.data(idx_fg + 10);

    auto ts_xxx_yyyz = pbuffer.data(idx_fg + 11);

    auto ts_xxx_yyzz = pbuffer.data(idx_fg + 12);

    auto ts_xxx_yzzz = pbuffer.data(idx_fg + 13);

    auto ts_xxx_zzzz = pbuffer.data(idx_fg + 14);

    auto ts_xxy_xxxx = pbuffer.data(idx_fg + 15);

    auto ts_xxy_xxxy = pbuffer.data(idx_fg + 16);

    auto ts_xxy_xxxz = pbuffer.data(idx_fg + 17);

    auto ts_xxy_xxyy = pbuffer.data(idx_fg + 18);

    auto ts_xxy_xxzz = pbuffer.data(idx_fg + 20);

    auto ts_xxy_xyyy = pbuffer.data(idx_fg + 21);

    auto ts_xxy_xzzz = pbuffer.data(idx_fg + 24);

    auto ts_xxy_yyyy = pbuffer.data(idx_fg + 25);

    auto ts_xxy_yyyz = pbuffer.data(idx_fg + 26);

    auto ts_xxy_yyzz = pbuffer.data(idx_fg + 27);

    auto ts_xxy_yzzz = pbuffer.data(idx_fg + 28);

    auto ts_xxz_xxxx = pbuffer.data(idx_fg + 30);

    auto ts_xxz_xxxy = pbuffer.data(idx_fg + 31);

    auto ts_xxz_xxxz = pbuffer.data(idx_fg + 32);

    auto ts_xxz_xxyy = pbuffer.data(idx_fg + 33);

    auto ts_xxz_xxyz = pbuffer.data(idx_fg + 34);

    auto ts_xxz_xxzz = pbuffer.data(idx_fg + 35);

    auto ts_xxz_xyyy = pbuffer.data(idx_fg + 36);

    auto ts_xxz_xyyz = pbuffer.data(idx_fg + 37);

    auto ts_xxz_xyzz = pbuffer.data(idx_fg + 38);

    auto ts_xxz_xzzz = pbuffer.data(idx_fg + 39);

    auto ts_xxz_yyyz = pbuffer.data(idx_fg + 41);

    auto ts_xxz_yyzz = pbuffer.data(idx_fg + 42);

    auto ts_xxz_yzzz = pbuffer.data(idx_fg + 43);

    auto ts_xxz_zzzz = pbuffer.data(idx_fg + 44);

    auto ts_xyy_xxxx = pbuffer.data(idx_fg + 45);

    auto ts_xyy_xxxy = pbuffer.data(idx_fg + 46);

    auto ts_xyy_xxyy = pbuffer.data(idx_fg + 48);

    auto ts_xyy_xxyz = pbuffer.data(idx_fg + 49);

    auto ts_xyy_xyyy = pbuffer.data(idx_fg + 51);

    auto ts_xyy_xyyz = pbuffer.data(idx_fg + 52);

    auto ts_xyy_xyzz = pbuffer.data(idx_fg + 53);

    auto ts_xyy_yyyy = pbuffer.data(idx_fg + 55);

    auto ts_xyy_yyyz = pbuffer.data(idx_fg + 56);

    auto ts_xyy_yyzz = pbuffer.data(idx_fg + 57);

    auto ts_xyy_yzzz = pbuffer.data(idx_fg + 58);

    auto ts_xyy_zzzz = pbuffer.data(idx_fg + 59);

    auto ts_xyz_yyyz = pbuffer.data(idx_fg + 71);

    auto ts_xyz_yyzz = pbuffer.data(idx_fg + 72);

    auto ts_xyz_yzzz = pbuffer.data(idx_fg + 73);

    auto ts_xzz_xxxx = pbuffer.data(idx_fg + 75);

    auto ts_xzz_xxxz = pbuffer.data(idx_fg + 77);

    auto ts_xzz_xxyz = pbuffer.data(idx_fg + 79);

    auto ts_xzz_xxzz = pbuffer.data(idx_fg + 80);

    auto ts_xzz_xyyz = pbuffer.data(idx_fg + 82);

    auto ts_xzz_xyzz = pbuffer.data(idx_fg + 83);

    auto ts_xzz_xzzz = pbuffer.data(idx_fg + 84);

    auto ts_xzz_yyyy = pbuffer.data(idx_fg + 85);

    auto ts_xzz_yyyz = pbuffer.data(idx_fg + 86);

    auto ts_xzz_yyzz = pbuffer.data(idx_fg + 87);

    auto ts_xzz_yzzz = pbuffer.data(idx_fg + 88);

    auto ts_xzz_zzzz = pbuffer.data(idx_fg + 89);

    auto ts_yyy_xxxx = pbuffer.data(idx_fg + 90);

    auto ts_yyy_xxxy = pbuffer.data(idx_fg + 91);

    auto ts_yyy_xxxz = pbuffer.data(idx_fg + 92);

    auto ts_yyy_xxyy = pbuffer.data(idx_fg + 93);

    auto ts_yyy_xxyz = pbuffer.data(idx_fg + 94);

    auto ts_yyy_xxzz = pbuffer.data(idx_fg + 95);

    auto ts_yyy_xyyy = pbuffer.data(idx_fg + 96);

    auto ts_yyy_xyyz = pbuffer.data(idx_fg + 97);

    auto ts_yyy_xyzz = pbuffer.data(idx_fg + 98);

    auto ts_yyy_xzzz = pbuffer.data(idx_fg + 99);

    auto ts_yyy_yyyy = pbuffer.data(idx_fg + 100);

    auto ts_yyy_yyyz = pbuffer.data(idx_fg + 101);

    auto ts_yyy_yyzz = pbuffer.data(idx_fg + 102);

    auto ts_yyy_yzzz = pbuffer.data(idx_fg + 103);

    auto ts_yyy_zzzz = pbuffer.data(idx_fg + 104);

    auto ts_yyz_xxxy = pbuffer.data(idx_fg + 106);

    auto ts_yyz_xxxz = pbuffer.data(idx_fg + 107);

    auto ts_yyz_xxyy = pbuffer.data(idx_fg + 108);

    auto ts_yyz_xxyz = pbuffer.data(idx_fg + 109);

    auto ts_yyz_xxzz = pbuffer.data(idx_fg + 110);

    auto ts_yyz_xyyy = pbuffer.data(idx_fg + 111);

    auto ts_yyz_xyyz = pbuffer.data(idx_fg + 112);

    auto ts_yyz_xyzz = pbuffer.data(idx_fg + 113);

    auto ts_yyz_xzzz = pbuffer.data(idx_fg + 114);

    auto ts_yyz_yyyy = pbuffer.data(idx_fg + 115);

    auto ts_yyz_yyyz = pbuffer.data(idx_fg + 116);

    auto ts_yyz_yyzz = pbuffer.data(idx_fg + 117);

    auto ts_yyz_yzzz = pbuffer.data(idx_fg + 118);

    auto ts_yyz_zzzz = pbuffer.data(idx_fg + 119);

    auto ts_yzz_xxxx = pbuffer.data(idx_fg + 120);

    auto ts_yzz_xxxy = pbuffer.data(idx_fg + 121);

    auto ts_yzz_xxxz = pbuffer.data(idx_fg + 122);

    auto ts_yzz_xxyy = pbuffer.data(idx_fg + 123);

    auto ts_yzz_xxyz = pbuffer.data(idx_fg + 124);

    auto ts_yzz_xxzz = pbuffer.data(idx_fg + 125);

    auto ts_yzz_xyyy = pbuffer.data(idx_fg + 126);

    auto ts_yzz_xyyz = pbuffer.data(idx_fg + 127);

    auto ts_yzz_xyzz = pbuffer.data(idx_fg + 128);

    auto ts_yzz_xzzz = pbuffer.data(idx_fg + 129);

    auto ts_yzz_yyyy = pbuffer.data(idx_fg + 130);

    auto ts_yzz_yyyz = pbuffer.data(idx_fg + 131);

    auto ts_yzz_yyzz = pbuffer.data(idx_fg + 132);

    auto ts_yzz_yzzz = pbuffer.data(idx_fg + 133);

    auto ts_yzz_zzzz = pbuffer.data(idx_fg + 134);

    auto ts_zzz_xxxx = pbuffer.data(idx_fg + 135);

    auto ts_zzz_xxxy = pbuffer.data(idx_fg + 136);

    auto ts_zzz_xxxz = pbuffer.data(idx_fg + 137);

    auto ts_zzz_xxyy = pbuffer.data(idx_fg + 138);

    auto ts_zzz_xxyz = pbuffer.data(idx_fg + 139);

    auto ts_zzz_xxzz = pbuffer.data(idx_fg + 140);

    auto ts_zzz_xyyy = pbuffer.data(idx_fg + 141);

    auto ts_zzz_xyyz = pbuffer.data(idx_fg + 142);

    auto ts_zzz_xyzz = pbuffer.data(idx_fg + 143);

    auto ts_zzz_xzzz = pbuffer.data(idx_fg + 144);

    auto ts_zzz_yyyy = pbuffer.data(idx_fg + 145);

    auto ts_zzz_yyyz = pbuffer.data(idx_fg + 146);

    auto ts_zzz_yyzz = pbuffer.data(idx_fg + 147);

    auto ts_zzz_yzzz = pbuffer.data(idx_fg + 148);

    auto ts_zzz_zzzz = pbuffer.data(idx_fg + 149);

    // Set up 0-15 components of targeted buffer : GG

    auto ts_xxxx_xxxx = pbuffer.data(idx_gg);

    auto ts_xxxx_xxxy = pbuffer.data(idx_gg + 1);

    auto ts_xxxx_xxxz = pbuffer.data(idx_gg + 2);

    auto ts_xxxx_xxyy = pbuffer.data(idx_gg + 3);

    auto ts_xxxx_xxyz = pbuffer.data(idx_gg + 4);

    auto ts_xxxx_xxzz = pbuffer.data(idx_gg + 5);

    auto ts_xxxx_xyyy = pbuffer.data(idx_gg + 6);

    auto ts_xxxx_xyyz = pbuffer.data(idx_gg + 7);

    auto ts_xxxx_xyzz = pbuffer.data(idx_gg + 8);

    auto ts_xxxx_xzzz = pbuffer.data(idx_gg + 9);

    auto ts_xxxx_yyyy = pbuffer.data(idx_gg + 10);

    auto ts_xxxx_yyyz = pbuffer.data(idx_gg + 11);

    auto ts_xxxx_yyzz = pbuffer.data(idx_gg + 12);

    auto ts_xxxx_yzzz = pbuffer.data(idx_gg + 13);

    auto ts_xxxx_zzzz = pbuffer.data(idx_gg + 14);

    #pragma omp simd aligned(ga_x, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyy, ts_xxx_yyyy, ts_xxx_yyyz, ts_xxx_yyz, ts_xxx_yyzz, ts_xxx_yzz, ts_xxx_yzzz, ts_xxx_zzz, ts_xxx_zzzz, ts_xxxx_xxxx, ts_xxxx_xxxy, ts_xxxx_xxxz, ts_xxxx_xxyy, ts_xxxx_xxyz, ts_xxxx_xxzz, ts_xxxx_xyyy, ts_xxxx_xyyz, ts_xxxx_xyzz, ts_xxxx_xzzz, ts_xxxx_yyyy, ts_xxxx_yyyz, ts_xxxx_yyzz, ts_xxxx_yzzz, ts_xxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxx_xxxx[i] = 3.0 * ts_xx_xxxx[i] * gfe_0 + 4.0 * ts_xxx_xxx[i] * gfe_0 + ts_xxx_xxxx[i] * ga_x[i];

        ts_xxxx_xxxy[i] = 3.0 * ts_xx_xxxy[i] * gfe_0 + 3.0 * ts_xxx_xxy[i] * gfe_0 + ts_xxx_xxxy[i] * ga_x[i];

        ts_xxxx_xxxz[i] = 3.0 * ts_xx_xxxz[i] * gfe_0 + 3.0 * ts_xxx_xxz[i] * gfe_0 + ts_xxx_xxxz[i] * ga_x[i];

        ts_xxxx_xxyy[i] = 3.0 * ts_xx_xxyy[i] * gfe_0 + 2.0 * ts_xxx_xyy[i] * gfe_0 + ts_xxx_xxyy[i] * ga_x[i];

        ts_xxxx_xxyz[i] = 3.0 * ts_xx_xxyz[i] * gfe_0 + 2.0 * ts_xxx_xyz[i] * gfe_0 + ts_xxx_xxyz[i] * ga_x[i];

        ts_xxxx_xxzz[i] = 3.0 * ts_xx_xxzz[i] * gfe_0 + 2.0 * ts_xxx_xzz[i] * gfe_0 + ts_xxx_xxzz[i] * ga_x[i];

        ts_xxxx_xyyy[i] = 3.0 * ts_xx_xyyy[i] * gfe_0 + ts_xxx_yyy[i] * gfe_0 + ts_xxx_xyyy[i] * ga_x[i];

        ts_xxxx_xyyz[i] = 3.0 * ts_xx_xyyz[i] * gfe_0 + ts_xxx_yyz[i] * gfe_0 + ts_xxx_xyyz[i] * ga_x[i];

        ts_xxxx_xyzz[i] = 3.0 * ts_xx_xyzz[i] * gfe_0 + ts_xxx_yzz[i] * gfe_0 + ts_xxx_xyzz[i] * ga_x[i];

        ts_xxxx_xzzz[i] = 3.0 * ts_xx_xzzz[i] * gfe_0 + ts_xxx_zzz[i] * gfe_0 + ts_xxx_xzzz[i] * ga_x[i];

        ts_xxxx_yyyy[i] = 3.0 * ts_xx_yyyy[i] * gfe_0 + ts_xxx_yyyy[i] * ga_x[i];

        ts_xxxx_yyyz[i] = 3.0 * ts_xx_yyyz[i] * gfe_0 + ts_xxx_yyyz[i] * ga_x[i];

        ts_xxxx_yyzz[i] = 3.0 * ts_xx_yyzz[i] * gfe_0 + ts_xxx_yyzz[i] * ga_x[i];

        ts_xxxx_yzzz[i] = 3.0 * ts_xx_yzzz[i] * gfe_0 + ts_xxx_yzzz[i] * ga_x[i];

        ts_xxxx_zzzz[i] = 3.0 * ts_xx_zzzz[i] * gfe_0 + ts_xxx_zzzz[i] * ga_x[i];
    }

    // Set up 15-30 components of targeted buffer : GG

    auto ts_xxxy_xxxx = pbuffer.data(idx_gg + 15);

    auto ts_xxxy_xxxy = pbuffer.data(idx_gg + 16);

    auto ts_xxxy_xxxz = pbuffer.data(idx_gg + 17);

    auto ts_xxxy_xxyy = pbuffer.data(idx_gg + 18);

    auto ts_xxxy_xxyz = pbuffer.data(idx_gg + 19);

    auto ts_xxxy_xxzz = pbuffer.data(idx_gg + 20);

    auto ts_xxxy_xyyy = pbuffer.data(idx_gg + 21);

    auto ts_xxxy_xyyz = pbuffer.data(idx_gg + 22);

    auto ts_xxxy_xyzz = pbuffer.data(idx_gg + 23);

    auto ts_xxxy_xzzz = pbuffer.data(idx_gg + 24);

    auto ts_xxxy_yyyy = pbuffer.data(idx_gg + 25);

    auto ts_xxxy_yyyz = pbuffer.data(idx_gg + 26);

    auto ts_xxxy_yyzz = pbuffer.data(idx_gg + 27);

    auto ts_xxxy_yzzz = pbuffer.data(idx_gg + 28);

    auto ts_xxxy_zzzz = pbuffer.data(idx_gg + 29);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_zzzz, ts_xxxy_xxxx, ts_xxxy_xxxy, ts_xxxy_xxxz, ts_xxxy_xxyy, ts_xxxy_xxyz, ts_xxxy_xxzz, ts_xxxy_xyyy, ts_xxxy_xyyz, ts_xxxy_xyzz, ts_xxxy_xzzz, ts_xxxy_yyyy, ts_xxxy_yyyz, ts_xxxy_yyzz, ts_xxxy_yzzz, ts_xxxy_zzzz, ts_xxy_yyyy, ts_xxy_yyyz, ts_xxy_yyzz, ts_xxy_yzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxy_xxxx[i] = ts_xxx_xxxx[i] * ga_y[i];

        ts_xxxy_xxxy[i] = ts_xxx_xxx[i] * gfe_0 + ts_xxx_xxxy[i] * ga_y[i];

        ts_xxxy_xxxz[i] = ts_xxx_xxxz[i] * ga_y[i];

        ts_xxxy_xxyy[i] = 2.0 * ts_xxx_xxy[i] * gfe_0 + ts_xxx_xxyy[i] * ga_y[i];

        ts_xxxy_xxyz[i] = ts_xxx_xxz[i] * gfe_0 + ts_xxx_xxyz[i] * ga_y[i];

        ts_xxxy_xxzz[i] = ts_xxx_xxzz[i] * ga_y[i];

        ts_xxxy_xyyy[i] = 3.0 * ts_xxx_xyy[i] * gfe_0 + ts_xxx_xyyy[i] * ga_y[i];

        ts_xxxy_xyyz[i] = 2.0 * ts_xxx_xyz[i] * gfe_0 + ts_xxx_xyyz[i] * ga_y[i];

        ts_xxxy_xyzz[i] = ts_xxx_xzz[i] * gfe_0 + ts_xxx_xyzz[i] * ga_y[i];

        ts_xxxy_xzzz[i] = ts_xxx_xzzz[i] * ga_y[i];

        ts_xxxy_yyyy[i] = 2.0 * ts_xy_yyyy[i] * gfe_0 + ts_xxy_yyyy[i] * ga_x[i];

        ts_xxxy_yyyz[i] = 2.0 * ts_xy_yyyz[i] * gfe_0 + ts_xxy_yyyz[i] * ga_x[i];

        ts_xxxy_yyzz[i] = 2.0 * ts_xy_yyzz[i] * gfe_0 + ts_xxy_yyzz[i] * ga_x[i];

        ts_xxxy_yzzz[i] = 2.0 * ts_xy_yzzz[i] * gfe_0 + ts_xxy_yzzz[i] * ga_x[i];

        ts_xxxy_zzzz[i] = ts_xxx_zzzz[i] * ga_y[i];
    }

    // Set up 30-45 components of targeted buffer : GG

    auto ts_xxxz_xxxx = pbuffer.data(idx_gg + 30);

    auto ts_xxxz_xxxy = pbuffer.data(idx_gg + 31);

    auto ts_xxxz_xxxz = pbuffer.data(idx_gg + 32);

    auto ts_xxxz_xxyy = pbuffer.data(idx_gg + 33);

    auto ts_xxxz_xxyz = pbuffer.data(idx_gg + 34);

    auto ts_xxxz_xxzz = pbuffer.data(idx_gg + 35);

    auto ts_xxxz_xyyy = pbuffer.data(idx_gg + 36);

    auto ts_xxxz_xyyz = pbuffer.data(idx_gg + 37);

    auto ts_xxxz_xyzz = pbuffer.data(idx_gg + 38);

    auto ts_xxxz_xzzz = pbuffer.data(idx_gg + 39);

    auto ts_xxxz_yyyy = pbuffer.data(idx_gg + 40);

    auto ts_xxxz_yyyz = pbuffer.data(idx_gg + 41);

    auto ts_xxxz_yyzz = pbuffer.data(idx_gg + 42);

    auto ts_xxxz_yzzz = pbuffer.data(idx_gg + 43);

    auto ts_xxxz_zzzz = pbuffer.data(idx_gg + 44);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxx_xxx, ts_xxx_xxxx, ts_xxx_xxxy, ts_xxx_xxxz, ts_xxx_xxy, ts_xxx_xxyy, ts_xxx_xxyz, ts_xxx_xxz, ts_xxx_xxzz, ts_xxx_xyy, ts_xxx_xyyy, ts_xxx_xyyz, ts_xxx_xyz, ts_xxx_xyzz, ts_xxx_xzz, ts_xxx_xzzz, ts_xxx_yyyy, ts_xxxz_xxxx, ts_xxxz_xxxy, ts_xxxz_xxxz, ts_xxxz_xxyy, ts_xxxz_xxyz, ts_xxxz_xxzz, ts_xxxz_xyyy, ts_xxxz_xyyz, ts_xxxz_xyzz, ts_xxxz_xzzz, ts_xxxz_yyyy, ts_xxxz_yyyz, ts_xxxz_yyzz, ts_xxxz_yzzz, ts_xxxz_zzzz, ts_xxz_yyyz, ts_xxz_yyzz, ts_xxz_yzzz, ts_xxz_zzzz, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxz_xxxx[i] = ts_xxx_xxxx[i] * ga_z[i];

        ts_xxxz_xxxy[i] = ts_xxx_xxxy[i] * ga_z[i];

        ts_xxxz_xxxz[i] = ts_xxx_xxx[i] * gfe_0 + ts_xxx_xxxz[i] * ga_z[i];

        ts_xxxz_xxyy[i] = ts_xxx_xxyy[i] * ga_z[i];

        ts_xxxz_xxyz[i] = ts_xxx_xxy[i] * gfe_0 + ts_xxx_xxyz[i] * ga_z[i];

        ts_xxxz_xxzz[i] = 2.0 * ts_xxx_xxz[i] * gfe_0 + ts_xxx_xxzz[i] * ga_z[i];

        ts_xxxz_xyyy[i] = ts_xxx_xyyy[i] * ga_z[i];

        ts_xxxz_xyyz[i] = ts_xxx_xyy[i] * gfe_0 + ts_xxx_xyyz[i] * ga_z[i];

        ts_xxxz_xyzz[i] = 2.0 * ts_xxx_xyz[i] * gfe_0 + ts_xxx_xyzz[i] * ga_z[i];

        ts_xxxz_xzzz[i] = 3.0 * ts_xxx_xzz[i] * gfe_0 + ts_xxx_xzzz[i] * ga_z[i];

        ts_xxxz_yyyy[i] = ts_xxx_yyyy[i] * ga_z[i];

        ts_xxxz_yyyz[i] = 2.0 * ts_xz_yyyz[i] * gfe_0 + ts_xxz_yyyz[i] * ga_x[i];

        ts_xxxz_yyzz[i] = 2.0 * ts_xz_yyzz[i] * gfe_0 + ts_xxz_yyzz[i] * ga_x[i];

        ts_xxxz_yzzz[i] = 2.0 * ts_xz_yzzz[i] * gfe_0 + ts_xxz_yzzz[i] * ga_x[i];

        ts_xxxz_zzzz[i] = 2.0 * ts_xz_zzzz[i] * gfe_0 + ts_xxz_zzzz[i] * ga_x[i];
    }

    // Set up 45-60 components of targeted buffer : GG

    auto ts_xxyy_xxxx = pbuffer.data(idx_gg + 45);

    auto ts_xxyy_xxxy = pbuffer.data(idx_gg + 46);

    auto ts_xxyy_xxxz = pbuffer.data(idx_gg + 47);

    auto ts_xxyy_xxyy = pbuffer.data(idx_gg + 48);

    auto ts_xxyy_xxyz = pbuffer.data(idx_gg + 49);

    auto ts_xxyy_xxzz = pbuffer.data(idx_gg + 50);

    auto ts_xxyy_xyyy = pbuffer.data(idx_gg + 51);

    auto ts_xxyy_xyyz = pbuffer.data(idx_gg + 52);

    auto ts_xxyy_xyzz = pbuffer.data(idx_gg + 53);

    auto ts_xxyy_xzzz = pbuffer.data(idx_gg + 54);

    auto ts_xxyy_yyyy = pbuffer.data(idx_gg + 55);

    auto ts_xxyy_yyyz = pbuffer.data(idx_gg + 56);

    auto ts_xxyy_yyzz = pbuffer.data(idx_gg + 57);

    auto ts_xxyy_yzzz = pbuffer.data(idx_gg + 58);

    auto ts_xxyy_zzzz = pbuffer.data(idx_gg + 59);

    #pragma omp simd aligned(ga_x, ga_y, ts_xx_xxxx, ts_xx_xxxz, ts_xx_xxzz, ts_xx_xzzz, ts_xxy_xxxx, ts_xxy_xxxz, ts_xxy_xxzz, ts_xxy_xzzz, ts_xxyy_xxxx, ts_xxyy_xxxy, ts_xxyy_xxxz, ts_xxyy_xxyy, ts_xxyy_xxyz, ts_xxyy_xxzz, ts_xxyy_xyyy, ts_xxyy_xyyz, ts_xxyy_xyzz, ts_xxyy_xzzz, ts_xxyy_yyyy, ts_xxyy_yyyz, ts_xxyy_yyzz, ts_xxyy_yzzz, ts_xxyy_zzzz, ts_xyy_xxxy, ts_xyy_xxy, ts_xyy_xxyy, ts_xyy_xxyz, ts_xyy_xyy, ts_xyy_xyyy, ts_xyy_xyyz, ts_xyy_xyz, ts_xyy_xyzz, ts_xyy_yyy, ts_xyy_yyyy, ts_xyy_yyyz, ts_xyy_yyz, ts_xyy_yyzz, ts_xyy_yzz, ts_xyy_yzzz, ts_xyy_zzzz, ts_yy_xxxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyy_xxxx[i] = ts_xx_xxxx[i] * gfe_0 + ts_xxy_xxxx[i] * ga_y[i];

        ts_xxyy_xxxy[i] = ts_yy_xxxy[i] * gfe_0 + 3.0 * ts_xyy_xxy[i] * gfe_0 + ts_xyy_xxxy[i] * ga_x[i];

        ts_xxyy_xxxz[i] = ts_xx_xxxz[i] * gfe_0 + ts_xxy_xxxz[i] * ga_y[i];

        ts_xxyy_xxyy[i] = ts_yy_xxyy[i] * gfe_0 + 2.0 * ts_xyy_xyy[i] * gfe_0 + ts_xyy_xxyy[i] * ga_x[i];

        ts_xxyy_xxyz[i] = ts_yy_xxyz[i] * gfe_0 + 2.0 * ts_xyy_xyz[i] * gfe_0 + ts_xyy_xxyz[i] * ga_x[i];

        ts_xxyy_xxzz[i] = ts_xx_xxzz[i] * gfe_0 + ts_xxy_xxzz[i] * ga_y[i];

        ts_xxyy_xyyy[i] = ts_yy_xyyy[i] * gfe_0 + ts_xyy_yyy[i] * gfe_0 + ts_xyy_xyyy[i] * ga_x[i];

        ts_xxyy_xyyz[i] = ts_yy_xyyz[i] * gfe_0 + ts_xyy_yyz[i] * gfe_0 + ts_xyy_xyyz[i] * ga_x[i];

        ts_xxyy_xyzz[i] = ts_yy_xyzz[i] * gfe_0 + ts_xyy_yzz[i] * gfe_0 + ts_xyy_xyzz[i] * ga_x[i];

        ts_xxyy_xzzz[i] = ts_xx_xzzz[i] * gfe_0 + ts_xxy_xzzz[i] * ga_y[i];

        ts_xxyy_yyyy[i] = ts_yy_yyyy[i] * gfe_0 + ts_xyy_yyyy[i] * ga_x[i];

        ts_xxyy_yyyz[i] = ts_yy_yyyz[i] * gfe_0 + ts_xyy_yyyz[i] * ga_x[i];

        ts_xxyy_yyzz[i] = ts_yy_yyzz[i] * gfe_0 + ts_xyy_yyzz[i] * ga_x[i];

        ts_xxyy_yzzz[i] = ts_yy_yzzz[i] * gfe_0 + ts_xyy_yzzz[i] * ga_x[i];

        ts_xxyy_zzzz[i] = ts_yy_zzzz[i] * gfe_0 + ts_xyy_zzzz[i] * ga_x[i];
    }

    // Set up 60-75 components of targeted buffer : GG

    auto ts_xxyz_xxxx = pbuffer.data(idx_gg + 60);

    auto ts_xxyz_xxxy = pbuffer.data(idx_gg + 61);

    auto ts_xxyz_xxxz = pbuffer.data(idx_gg + 62);

    auto ts_xxyz_xxyy = pbuffer.data(idx_gg + 63);

    auto ts_xxyz_xxyz = pbuffer.data(idx_gg + 64);

    auto ts_xxyz_xxzz = pbuffer.data(idx_gg + 65);

    auto ts_xxyz_xyyy = pbuffer.data(idx_gg + 66);

    auto ts_xxyz_xyyz = pbuffer.data(idx_gg + 67);

    auto ts_xxyz_xyzz = pbuffer.data(idx_gg + 68);

    auto ts_xxyz_xzzz = pbuffer.data(idx_gg + 69);

    auto ts_xxyz_yyyy = pbuffer.data(idx_gg + 70);

    auto ts_xxyz_yyyz = pbuffer.data(idx_gg + 71);

    auto ts_xxyz_yyzz = pbuffer.data(idx_gg + 72);

    auto ts_xxyz_yzzz = pbuffer.data(idx_gg + 73);

    auto ts_xxyz_zzzz = pbuffer.data(idx_gg + 74);

    #pragma omp simd aligned(ga_x, ga_y, ga_z, ts_xxy_xxxy, ts_xxy_xxyy, ts_xxy_xyyy, ts_xxy_yyyy, ts_xxyz_xxxx, ts_xxyz_xxxy, ts_xxyz_xxxz, ts_xxyz_xxyy, ts_xxyz_xxyz, ts_xxyz_xxzz, ts_xxyz_xyyy, ts_xxyz_xyyz, ts_xxyz_xyzz, ts_xxyz_xzzz, ts_xxyz_yyyy, ts_xxyz_yyyz, ts_xxyz_yyzz, ts_xxyz_yzzz, ts_xxyz_zzzz, ts_xxz_xxxx, ts_xxz_xxxz, ts_xxz_xxyz, ts_xxz_xxz, ts_xxz_xxzz, ts_xxz_xyyz, ts_xxz_xyz, ts_xxz_xyzz, ts_xxz_xzz, ts_xxz_xzzz, ts_xxz_zzzz, ts_xyz_yyyz, ts_xyz_yyzz, ts_xyz_yzzz, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyz_xxxx[i] = ts_xxz_xxxx[i] * ga_y[i];

        ts_xxyz_xxxy[i] = ts_xxy_xxxy[i] * ga_z[i];

        ts_xxyz_xxxz[i] = ts_xxz_xxxz[i] * ga_y[i];

        ts_xxyz_xxyy[i] = ts_xxy_xxyy[i] * ga_z[i];

        ts_xxyz_xxyz[i] = ts_xxz_xxz[i] * gfe_0 + ts_xxz_xxyz[i] * ga_y[i];

        ts_xxyz_xxzz[i] = ts_xxz_xxzz[i] * ga_y[i];

        ts_xxyz_xyyy[i] = ts_xxy_xyyy[i] * ga_z[i];

        ts_xxyz_xyyz[i] = 2.0 * ts_xxz_xyz[i] * gfe_0 + ts_xxz_xyyz[i] * ga_y[i];

        ts_xxyz_xyzz[i] = ts_xxz_xzz[i] * gfe_0 + ts_xxz_xyzz[i] * ga_y[i];

        ts_xxyz_xzzz[i] = ts_xxz_xzzz[i] * ga_y[i];

        ts_xxyz_yyyy[i] = ts_xxy_yyyy[i] * ga_z[i];

        ts_xxyz_yyyz[i] = ts_yz_yyyz[i] * gfe_0 + ts_xyz_yyyz[i] * ga_x[i];

        ts_xxyz_yyzz[i] = ts_yz_yyzz[i] * gfe_0 + ts_xyz_yyzz[i] * ga_x[i];

        ts_xxyz_yzzz[i] = ts_yz_yzzz[i] * gfe_0 + ts_xyz_yzzz[i] * ga_x[i];

        ts_xxyz_zzzz[i] = ts_xxz_zzzz[i] * ga_y[i];
    }

    // Set up 75-90 components of targeted buffer : GG

    auto ts_xxzz_xxxx = pbuffer.data(idx_gg + 75);

    auto ts_xxzz_xxxy = pbuffer.data(idx_gg + 76);

    auto ts_xxzz_xxxz = pbuffer.data(idx_gg + 77);

    auto ts_xxzz_xxyy = pbuffer.data(idx_gg + 78);

    auto ts_xxzz_xxyz = pbuffer.data(idx_gg + 79);

    auto ts_xxzz_xxzz = pbuffer.data(idx_gg + 80);

    auto ts_xxzz_xyyy = pbuffer.data(idx_gg + 81);

    auto ts_xxzz_xyyz = pbuffer.data(idx_gg + 82);

    auto ts_xxzz_xyzz = pbuffer.data(idx_gg + 83);

    auto ts_xxzz_xzzz = pbuffer.data(idx_gg + 84);

    auto ts_xxzz_yyyy = pbuffer.data(idx_gg + 85);

    auto ts_xxzz_yyyz = pbuffer.data(idx_gg + 86);

    auto ts_xxzz_yyzz = pbuffer.data(idx_gg + 87);

    auto ts_xxzz_yzzz = pbuffer.data(idx_gg + 88);

    auto ts_xxzz_zzzz = pbuffer.data(idx_gg + 89);

    #pragma omp simd aligned(ga_x, ga_z, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxyy, ts_xx_xyyy, ts_xxz_xxxx, ts_xxz_xxxy, ts_xxz_xxyy, ts_xxz_xyyy, ts_xxzz_xxxx, ts_xxzz_xxxy, ts_xxzz_xxxz, ts_xxzz_xxyy, ts_xxzz_xxyz, ts_xxzz_xxzz, ts_xxzz_xyyy, ts_xxzz_xyyz, ts_xxzz_xyzz, ts_xxzz_xzzz, ts_xxzz_yyyy, ts_xxzz_yyyz, ts_xxzz_yyzz, ts_xxzz_yzzz, ts_xxzz_zzzz, ts_xzz_xxxz, ts_xzz_xxyz, ts_xzz_xxz, ts_xzz_xxzz, ts_xzz_xyyz, ts_xzz_xyz, ts_xzz_xyzz, ts_xzz_xzz, ts_xzz_xzzz, ts_xzz_yyyy, ts_xzz_yyyz, ts_xzz_yyz, ts_xzz_yyzz, ts_xzz_yzz, ts_xzz_yzzz, ts_xzz_zzz, ts_xzz_zzzz, ts_zz_xxxz, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxzz_xxxx[i] = ts_xx_xxxx[i] * gfe_0 + ts_xxz_xxxx[i] * ga_z[i];

        ts_xxzz_xxxy[i] = ts_xx_xxxy[i] * gfe_0 + ts_xxz_xxxy[i] * ga_z[i];

        ts_xxzz_xxxz[i] = ts_zz_xxxz[i] * gfe_0 + 3.0 * ts_xzz_xxz[i] * gfe_0 + ts_xzz_xxxz[i] * ga_x[i];

        ts_xxzz_xxyy[i] = ts_xx_xxyy[i] * gfe_0 + ts_xxz_xxyy[i] * ga_z[i];

        ts_xxzz_xxyz[i] = ts_zz_xxyz[i] * gfe_0 + 2.0 * ts_xzz_xyz[i] * gfe_0 + ts_xzz_xxyz[i] * ga_x[i];

        ts_xxzz_xxzz[i] = ts_zz_xxzz[i] * gfe_0 + 2.0 * ts_xzz_xzz[i] * gfe_0 + ts_xzz_xxzz[i] * ga_x[i];

        ts_xxzz_xyyy[i] = ts_xx_xyyy[i] * gfe_0 + ts_xxz_xyyy[i] * ga_z[i];

        ts_xxzz_xyyz[i] = ts_zz_xyyz[i] * gfe_0 + ts_xzz_yyz[i] * gfe_0 + ts_xzz_xyyz[i] * ga_x[i];

        ts_xxzz_xyzz[i] = ts_zz_xyzz[i] * gfe_0 + ts_xzz_yzz[i] * gfe_0 + ts_xzz_xyzz[i] * ga_x[i];

        ts_xxzz_xzzz[i] = ts_zz_xzzz[i] * gfe_0 + ts_xzz_zzz[i] * gfe_0 + ts_xzz_xzzz[i] * ga_x[i];

        ts_xxzz_yyyy[i] = ts_zz_yyyy[i] * gfe_0 + ts_xzz_yyyy[i] * ga_x[i];

        ts_xxzz_yyyz[i] = ts_zz_yyyz[i] * gfe_0 + ts_xzz_yyyz[i] * ga_x[i];

        ts_xxzz_yyzz[i] = ts_zz_yyzz[i] * gfe_0 + ts_xzz_yyzz[i] * ga_x[i];

        ts_xxzz_yzzz[i] = ts_zz_yzzz[i] * gfe_0 + ts_xzz_yzzz[i] * ga_x[i];

        ts_xxzz_zzzz[i] = ts_zz_zzzz[i] * gfe_0 + ts_xzz_zzzz[i] * ga_x[i];
    }

    // Set up 90-105 components of targeted buffer : GG

    auto ts_xyyy_xxxx = pbuffer.data(idx_gg + 90);

    auto ts_xyyy_xxxy = pbuffer.data(idx_gg + 91);

    auto ts_xyyy_xxxz = pbuffer.data(idx_gg + 92);

    auto ts_xyyy_xxyy = pbuffer.data(idx_gg + 93);

    auto ts_xyyy_xxyz = pbuffer.data(idx_gg + 94);

    auto ts_xyyy_xxzz = pbuffer.data(idx_gg + 95);

    auto ts_xyyy_xyyy = pbuffer.data(idx_gg + 96);

    auto ts_xyyy_xyyz = pbuffer.data(idx_gg + 97);

    auto ts_xyyy_xyzz = pbuffer.data(idx_gg + 98);

    auto ts_xyyy_xzzz = pbuffer.data(idx_gg + 99);

    auto ts_xyyy_yyyy = pbuffer.data(idx_gg + 100);

    auto ts_xyyy_yyyz = pbuffer.data(idx_gg + 101);

    auto ts_xyyy_yyzz = pbuffer.data(idx_gg + 102);

    auto ts_xyyy_yzzz = pbuffer.data(idx_gg + 103);

    auto ts_xyyy_zzzz = pbuffer.data(idx_gg + 104);

    #pragma omp simd aligned(ga_x, ts_xyyy_xxxx, ts_xyyy_xxxy, ts_xyyy_xxxz, ts_xyyy_xxyy, ts_xyyy_xxyz, ts_xyyy_xxzz, ts_xyyy_xyyy, ts_xyyy_xyyz, ts_xyyy_xyzz, ts_xyyy_xzzz, ts_xyyy_yyyy, ts_xyyy_yyyz, ts_xyyy_yyzz, ts_xyyy_yzzz, ts_xyyy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyy_xxxx[i] = 4.0 * ts_yyy_xxx[i] * gfe_0 + ts_yyy_xxxx[i] * ga_x[i];

        ts_xyyy_xxxy[i] = 3.0 * ts_yyy_xxy[i] * gfe_0 + ts_yyy_xxxy[i] * ga_x[i];

        ts_xyyy_xxxz[i] = 3.0 * ts_yyy_xxz[i] * gfe_0 + ts_yyy_xxxz[i] * ga_x[i];

        ts_xyyy_xxyy[i] = 2.0 * ts_yyy_xyy[i] * gfe_0 + ts_yyy_xxyy[i] * ga_x[i];

        ts_xyyy_xxyz[i] = 2.0 * ts_yyy_xyz[i] * gfe_0 + ts_yyy_xxyz[i] * ga_x[i];

        ts_xyyy_xxzz[i] = 2.0 * ts_yyy_xzz[i] * gfe_0 + ts_yyy_xxzz[i] * ga_x[i];

        ts_xyyy_xyyy[i] = ts_yyy_yyy[i] * gfe_0 + ts_yyy_xyyy[i] * ga_x[i];

        ts_xyyy_xyyz[i] = ts_yyy_yyz[i] * gfe_0 + ts_yyy_xyyz[i] * ga_x[i];

        ts_xyyy_xyzz[i] = ts_yyy_yzz[i] * gfe_0 + ts_yyy_xyzz[i] * ga_x[i];

        ts_xyyy_xzzz[i] = ts_yyy_zzz[i] * gfe_0 + ts_yyy_xzzz[i] * ga_x[i];

        ts_xyyy_yyyy[i] = ts_yyy_yyyy[i] * ga_x[i];

        ts_xyyy_yyyz[i] = ts_yyy_yyyz[i] * ga_x[i];

        ts_xyyy_yyzz[i] = ts_yyy_yyzz[i] * ga_x[i];

        ts_xyyy_yzzz[i] = ts_yyy_yzzz[i] * ga_x[i];

        ts_xyyy_zzzz[i] = ts_yyy_zzzz[i] * ga_x[i];
    }

    // Set up 105-120 components of targeted buffer : GG

    auto ts_xyyz_xxxx = pbuffer.data(idx_gg + 105);

    auto ts_xyyz_xxxy = pbuffer.data(idx_gg + 106);

    auto ts_xyyz_xxxz = pbuffer.data(idx_gg + 107);

    auto ts_xyyz_xxyy = pbuffer.data(idx_gg + 108);

    auto ts_xyyz_xxyz = pbuffer.data(idx_gg + 109);

    auto ts_xyyz_xxzz = pbuffer.data(idx_gg + 110);

    auto ts_xyyz_xyyy = pbuffer.data(idx_gg + 111);

    auto ts_xyyz_xyyz = pbuffer.data(idx_gg + 112);

    auto ts_xyyz_xyzz = pbuffer.data(idx_gg + 113);

    auto ts_xyyz_xzzz = pbuffer.data(idx_gg + 114);

    auto ts_xyyz_yyyy = pbuffer.data(idx_gg + 115);

    auto ts_xyyz_yyyz = pbuffer.data(idx_gg + 116);

    auto ts_xyyz_yyzz = pbuffer.data(idx_gg + 117);

    auto ts_xyyz_yzzz = pbuffer.data(idx_gg + 118);

    auto ts_xyyz_zzzz = pbuffer.data(idx_gg + 119);

    #pragma omp simd aligned(ga_x, ga_z, ts_xyy_xxxx, ts_xyy_xxxy, ts_xyy_xxyy, ts_xyy_xyyy, ts_xyyz_xxxx, ts_xyyz_xxxy, ts_xyyz_xxxz, ts_xyyz_xxyy, ts_xyyz_xxyz, ts_xyyz_xxzz, ts_xyyz_xyyy, ts_xyyz_xyyz, ts_xyyz_xyzz, ts_xyyz_xzzz, ts_xyyz_yyyy, ts_xyyz_yyyz, ts_xyyz_yyzz, ts_xyyz_yzzz, ts_xyyz_zzzz, ts_yyz_xxxz, ts_yyz_xxyz, ts_yyz_xxz, ts_yyz_xxzz, ts_yyz_xyyz, ts_yyz_xyz, ts_yyz_xyzz, ts_yyz_xzz, ts_yyz_xzzz, ts_yyz_yyyy, ts_yyz_yyyz, ts_yyz_yyz, ts_yyz_yyzz, ts_yyz_yzz, ts_yyz_yzzz, ts_yyz_zzz, ts_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyz_xxxx[i] = ts_xyy_xxxx[i] * ga_z[i];

        ts_xyyz_xxxy[i] = ts_xyy_xxxy[i] * ga_z[i];

        ts_xyyz_xxxz[i] = 3.0 * ts_yyz_xxz[i] * gfe_0 + ts_yyz_xxxz[i] * ga_x[i];

        ts_xyyz_xxyy[i] = ts_xyy_xxyy[i] * ga_z[i];

        ts_xyyz_xxyz[i] = 2.0 * ts_yyz_xyz[i] * gfe_0 + ts_yyz_xxyz[i] * ga_x[i];

        ts_xyyz_xxzz[i] = 2.0 * ts_yyz_xzz[i] * gfe_0 + ts_yyz_xxzz[i] * ga_x[i];

        ts_xyyz_xyyy[i] = ts_xyy_xyyy[i] * ga_z[i];

        ts_xyyz_xyyz[i] = ts_yyz_yyz[i] * gfe_0 + ts_yyz_xyyz[i] * ga_x[i];

        ts_xyyz_xyzz[i] = ts_yyz_yzz[i] * gfe_0 + ts_yyz_xyzz[i] * ga_x[i];

        ts_xyyz_xzzz[i] = ts_yyz_zzz[i] * gfe_0 + ts_yyz_xzzz[i] * ga_x[i];

        ts_xyyz_yyyy[i] = ts_yyz_yyyy[i] * ga_x[i];

        ts_xyyz_yyyz[i] = ts_yyz_yyyz[i] * ga_x[i];

        ts_xyyz_yyzz[i] = ts_yyz_yyzz[i] * ga_x[i];

        ts_xyyz_yzzz[i] = ts_yyz_yzzz[i] * ga_x[i];

        ts_xyyz_zzzz[i] = ts_yyz_zzzz[i] * ga_x[i];
    }

    // Set up 120-135 components of targeted buffer : GG

    auto ts_xyzz_xxxx = pbuffer.data(idx_gg + 120);

    auto ts_xyzz_xxxy = pbuffer.data(idx_gg + 121);

    auto ts_xyzz_xxxz = pbuffer.data(idx_gg + 122);

    auto ts_xyzz_xxyy = pbuffer.data(idx_gg + 123);

    auto ts_xyzz_xxyz = pbuffer.data(idx_gg + 124);

    auto ts_xyzz_xxzz = pbuffer.data(idx_gg + 125);

    auto ts_xyzz_xyyy = pbuffer.data(idx_gg + 126);

    auto ts_xyzz_xyyz = pbuffer.data(idx_gg + 127);

    auto ts_xyzz_xyzz = pbuffer.data(idx_gg + 128);

    auto ts_xyzz_xzzz = pbuffer.data(idx_gg + 129);

    auto ts_xyzz_yyyy = pbuffer.data(idx_gg + 130);

    auto ts_xyzz_yyyz = pbuffer.data(idx_gg + 131);

    auto ts_xyzz_yyzz = pbuffer.data(idx_gg + 132);

    auto ts_xyzz_yzzz = pbuffer.data(idx_gg + 133);

    auto ts_xyzz_zzzz = pbuffer.data(idx_gg + 134);

    #pragma omp simd aligned(ga_x, ga_y, ts_xyzz_xxxx, ts_xyzz_xxxy, ts_xyzz_xxxz, ts_xyzz_xxyy, ts_xyzz_xxyz, ts_xyzz_xxzz, ts_xyzz_xyyy, ts_xyzz_xyyz, ts_xyzz_xyzz, ts_xyzz_xzzz, ts_xyzz_yyyy, ts_xyzz_yyyz, ts_xyzz_yyzz, ts_xyzz_yzzz, ts_xyzz_zzzz, ts_xzz_xxxx, ts_xzz_xxxz, ts_xzz_xxzz, ts_xzz_xzzz, ts_yzz_xxxy, ts_yzz_xxy, ts_yzz_xxyy, ts_yzz_xxyz, ts_yzz_xyy, ts_yzz_xyyy, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_yyy, ts_yzz_yyyy, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyzz_xxxx[i] = ts_xzz_xxxx[i] * ga_y[i];

        ts_xyzz_xxxy[i] = 3.0 * ts_yzz_xxy[i] * gfe_0 + ts_yzz_xxxy[i] * ga_x[i];

        ts_xyzz_xxxz[i] = ts_xzz_xxxz[i] * ga_y[i];

        ts_xyzz_xxyy[i] = 2.0 * ts_yzz_xyy[i] * gfe_0 + ts_yzz_xxyy[i] * ga_x[i];

        ts_xyzz_xxyz[i] = 2.0 * ts_yzz_xyz[i] * gfe_0 + ts_yzz_xxyz[i] * ga_x[i];

        ts_xyzz_xxzz[i] = ts_xzz_xxzz[i] * ga_y[i];

        ts_xyzz_xyyy[i] = ts_yzz_yyy[i] * gfe_0 + ts_yzz_xyyy[i] * ga_x[i];

        ts_xyzz_xyyz[i] = ts_yzz_yyz[i] * gfe_0 + ts_yzz_xyyz[i] * ga_x[i];

        ts_xyzz_xyzz[i] = ts_yzz_yzz[i] * gfe_0 + ts_yzz_xyzz[i] * ga_x[i];

        ts_xyzz_xzzz[i] = ts_xzz_xzzz[i] * ga_y[i];

        ts_xyzz_yyyy[i] = ts_yzz_yyyy[i] * ga_x[i];

        ts_xyzz_yyyz[i] = ts_yzz_yyyz[i] * ga_x[i];

        ts_xyzz_yyzz[i] = ts_yzz_yyzz[i] * ga_x[i];

        ts_xyzz_yzzz[i] = ts_yzz_yzzz[i] * ga_x[i];

        ts_xyzz_zzzz[i] = ts_yzz_zzzz[i] * ga_x[i];
    }

    // Set up 135-150 components of targeted buffer : GG

    auto ts_xzzz_xxxx = pbuffer.data(idx_gg + 135);

    auto ts_xzzz_xxxy = pbuffer.data(idx_gg + 136);

    auto ts_xzzz_xxxz = pbuffer.data(idx_gg + 137);

    auto ts_xzzz_xxyy = pbuffer.data(idx_gg + 138);

    auto ts_xzzz_xxyz = pbuffer.data(idx_gg + 139);

    auto ts_xzzz_xxzz = pbuffer.data(idx_gg + 140);

    auto ts_xzzz_xyyy = pbuffer.data(idx_gg + 141);

    auto ts_xzzz_xyyz = pbuffer.data(idx_gg + 142);

    auto ts_xzzz_xyzz = pbuffer.data(idx_gg + 143);

    auto ts_xzzz_xzzz = pbuffer.data(idx_gg + 144);

    auto ts_xzzz_yyyy = pbuffer.data(idx_gg + 145);

    auto ts_xzzz_yyyz = pbuffer.data(idx_gg + 146);

    auto ts_xzzz_yyzz = pbuffer.data(idx_gg + 147);

    auto ts_xzzz_yzzz = pbuffer.data(idx_gg + 148);

    auto ts_xzzz_zzzz = pbuffer.data(idx_gg + 149);

    #pragma omp simd aligned(ga_x, ts_xzzz_xxxx, ts_xzzz_xxxy, ts_xzzz_xxxz, ts_xzzz_xxyy, ts_xzzz_xxyz, ts_xzzz_xxzz, ts_xzzz_xyyy, ts_xzzz_xyyz, ts_xzzz_xyzz, ts_xzzz_xzzz, ts_xzzz_yyyy, ts_xzzz_yyyz, ts_xzzz_yyzz, ts_xzzz_yzzz, ts_xzzz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xzzz_xxxx[i] = 4.0 * ts_zzz_xxx[i] * gfe_0 + ts_zzz_xxxx[i] * ga_x[i];

        ts_xzzz_xxxy[i] = 3.0 * ts_zzz_xxy[i] * gfe_0 + ts_zzz_xxxy[i] * ga_x[i];

        ts_xzzz_xxxz[i] = 3.0 * ts_zzz_xxz[i] * gfe_0 + ts_zzz_xxxz[i] * ga_x[i];

        ts_xzzz_xxyy[i] = 2.0 * ts_zzz_xyy[i] * gfe_0 + ts_zzz_xxyy[i] * ga_x[i];

        ts_xzzz_xxyz[i] = 2.0 * ts_zzz_xyz[i] * gfe_0 + ts_zzz_xxyz[i] * ga_x[i];

        ts_xzzz_xxzz[i] = 2.0 * ts_zzz_xzz[i] * gfe_0 + ts_zzz_xxzz[i] * ga_x[i];

        ts_xzzz_xyyy[i] = ts_zzz_yyy[i] * gfe_0 + ts_zzz_xyyy[i] * ga_x[i];

        ts_xzzz_xyyz[i] = ts_zzz_yyz[i] * gfe_0 + ts_zzz_xyyz[i] * ga_x[i];

        ts_xzzz_xyzz[i] = ts_zzz_yzz[i] * gfe_0 + ts_zzz_xyzz[i] * ga_x[i];

        ts_xzzz_xzzz[i] = ts_zzz_zzz[i] * gfe_0 + ts_zzz_xzzz[i] * ga_x[i];

        ts_xzzz_yyyy[i] = ts_zzz_yyyy[i] * ga_x[i];

        ts_xzzz_yyyz[i] = ts_zzz_yyyz[i] * ga_x[i];

        ts_xzzz_yyzz[i] = ts_zzz_yyzz[i] * ga_x[i];

        ts_xzzz_yzzz[i] = ts_zzz_yzzz[i] * ga_x[i];

        ts_xzzz_zzzz[i] = ts_zzz_zzzz[i] * ga_x[i];
    }

    // Set up 150-165 components of targeted buffer : GG

    auto ts_yyyy_xxxx = pbuffer.data(idx_gg + 150);

    auto ts_yyyy_xxxy = pbuffer.data(idx_gg + 151);

    auto ts_yyyy_xxxz = pbuffer.data(idx_gg + 152);

    auto ts_yyyy_xxyy = pbuffer.data(idx_gg + 153);

    auto ts_yyyy_xxyz = pbuffer.data(idx_gg + 154);

    auto ts_yyyy_xxzz = pbuffer.data(idx_gg + 155);

    auto ts_yyyy_xyyy = pbuffer.data(idx_gg + 156);

    auto ts_yyyy_xyyz = pbuffer.data(idx_gg + 157);

    auto ts_yyyy_xyzz = pbuffer.data(idx_gg + 158);

    auto ts_yyyy_xzzz = pbuffer.data(idx_gg + 159);

    auto ts_yyyy_yyyy = pbuffer.data(idx_gg + 160);

    auto ts_yyyy_yyyz = pbuffer.data(idx_gg + 161);

    auto ts_yyyy_yyzz = pbuffer.data(idx_gg + 162);

    auto ts_yyyy_yzzz = pbuffer.data(idx_gg + 163);

    auto ts_yyyy_zzzz = pbuffer.data(idx_gg + 164);

    #pragma omp simd aligned(ga_y, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, ts_yyy_xxx, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxxz, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xxz, ts_yyy_xxzz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_xzz, ts_yyy_xzzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyy_zzz, ts_yyy_zzzz, ts_yyyy_xxxx, ts_yyyy_xxxy, ts_yyyy_xxxz, ts_yyyy_xxyy, ts_yyyy_xxyz, ts_yyyy_xxzz, ts_yyyy_xyyy, ts_yyyy_xyyz, ts_yyyy_xyzz, ts_yyyy_xzzz, ts_yyyy_yyyy, ts_yyyy_yyyz, ts_yyyy_yyzz, ts_yyyy_yzzz, ts_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyy_xxxx[i] = 3.0 * ts_yy_xxxx[i] * gfe_0 + ts_yyy_xxxx[i] * ga_y[i];

        ts_yyyy_xxxy[i] = 3.0 * ts_yy_xxxy[i] * gfe_0 + ts_yyy_xxx[i] * gfe_0 + ts_yyy_xxxy[i] * ga_y[i];

        ts_yyyy_xxxz[i] = 3.0 * ts_yy_xxxz[i] * gfe_0 + ts_yyy_xxxz[i] * ga_y[i];

        ts_yyyy_xxyy[i] = 3.0 * ts_yy_xxyy[i] * gfe_0 + 2.0 * ts_yyy_xxy[i] * gfe_0 + ts_yyy_xxyy[i] * ga_y[i];

        ts_yyyy_xxyz[i] = 3.0 * ts_yy_xxyz[i] * gfe_0 + ts_yyy_xxz[i] * gfe_0 + ts_yyy_xxyz[i] * ga_y[i];

        ts_yyyy_xxzz[i] = 3.0 * ts_yy_xxzz[i] * gfe_0 + ts_yyy_xxzz[i] * ga_y[i];

        ts_yyyy_xyyy[i] = 3.0 * ts_yy_xyyy[i] * gfe_0 + 3.0 * ts_yyy_xyy[i] * gfe_0 + ts_yyy_xyyy[i] * ga_y[i];

        ts_yyyy_xyyz[i] = 3.0 * ts_yy_xyyz[i] * gfe_0 + 2.0 * ts_yyy_xyz[i] * gfe_0 + ts_yyy_xyyz[i] * ga_y[i];

        ts_yyyy_xyzz[i] = 3.0 * ts_yy_xyzz[i] * gfe_0 + ts_yyy_xzz[i] * gfe_0 + ts_yyy_xyzz[i] * ga_y[i];

        ts_yyyy_xzzz[i] = 3.0 * ts_yy_xzzz[i] * gfe_0 + ts_yyy_xzzz[i] * ga_y[i];

        ts_yyyy_yyyy[i] = 3.0 * ts_yy_yyyy[i] * gfe_0 + 4.0 * ts_yyy_yyy[i] * gfe_0 + ts_yyy_yyyy[i] * ga_y[i];

        ts_yyyy_yyyz[i] = 3.0 * ts_yy_yyyz[i] * gfe_0 + 3.0 * ts_yyy_yyz[i] * gfe_0 + ts_yyy_yyyz[i] * ga_y[i];

        ts_yyyy_yyzz[i] = 3.0 * ts_yy_yyzz[i] * gfe_0 + 2.0 * ts_yyy_yzz[i] * gfe_0 + ts_yyy_yyzz[i] * ga_y[i];

        ts_yyyy_yzzz[i] = 3.0 * ts_yy_yzzz[i] * gfe_0 + ts_yyy_zzz[i] * gfe_0 + ts_yyy_yzzz[i] * ga_y[i];

        ts_yyyy_zzzz[i] = 3.0 * ts_yy_zzzz[i] * gfe_0 + ts_yyy_zzzz[i] * ga_y[i];
    }

    // Set up 165-180 components of targeted buffer : GG

    auto ts_yyyz_xxxx = pbuffer.data(idx_gg + 165);

    auto ts_yyyz_xxxy = pbuffer.data(idx_gg + 166);

    auto ts_yyyz_xxxz = pbuffer.data(idx_gg + 167);

    auto ts_yyyz_xxyy = pbuffer.data(idx_gg + 168);

    auto ts_yyyz_xxyz = pbuffer.data(idx_gg + 169);

    auto ts_yyyz_xxzz = pbuffer.data(idx_gg + 170);

    auto ts_yyyz_xyyy = pbuffer.data(idx_gg + 171);

    auto ts_yyyz_xyyz = pbuffer.data(idx_gg + 172);

    auto ts_yyyz_xyzz = pbuffer.data(idx_gg + 173);

    auto ts_yyyz_xzzz = pbuffer.data(idx_gg + 174);

    auto ts_yyyz_yyyy = pbuffer.data(idx_gg + 175);

    auto ts_yyyz_yyyz = pbuffer.data(idx_gg + 176);

    auto ts_yyyz_yyzz = pbuffer.data(idx_gg + 177);

    auto ts_yyyz_yzzz = pbuffer.data(idx_gg + 178);

    auto ts_yyyz_zzzz = pbuffer.data(idx_gg + 179);

    #pragma omp simd aligned(ga_y, ga_z, ts_yyy_xxxx, ts_yyy_xxxy, ts_yyy_xxy, ts_yyy_xxyy, ts_yyy_xxyz, ts_yyy_xyy, ts_yyy_xyyy, ts_yyy_xyyz, ts_yyy_xyz, ts_yyy_xyzz, ts_yyy_yyy, ts_yyy_yyyy, ts_yyy_yyyz, ts_yyy_yyz, ts_yyy_yyzz, ts_yyy_yzz, ts_yyy_yzzz, ts_yyyz_xxxx, ts_yyyz_xxxy, ts_yyyz_xxxz, ts_yyyz_xxyy, ts_yyyz_xxyz, ts_yyyz_xxzz, ts_yyyz_xyyy, ts_yyyz_xyyz, ts_yyyz_xyzz, ts_yyyz_xzzz, ts_yyyz_yyyy, ts_yyyz_yyyz, ts_yyyz_yyzz, ts_yyyz_yzzz, ts_yyyz_zzzz, ts_yyz_xxxz, ts_yyz_xxzz, ts_yyz_xzzz, ts_yyz_zzzz, ts_yz_xxxz, ts_yz_xxzz, ts_yz_xzzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyz_xxxx[i] = ts_yyy_xxxx[i] * ga_z[i];

        ts_yyyz_xxxy[i] = ts_yyy_xxxy[i] * ga_z[i];

        ts_yyyz_xxxz[i] = 2.0 * ts_yz_xxxz[i] * gfe_0 + ts_yyz_xxxz[i] * ga_y[i];

        ts_yyyz_xxyy[i] = ts_yyy_xxyy[i] * ga_z[i];

        ts_yyyz_xxyz[i] = ts_yyy_xxy[i] * gfe_0 + ts_yyy_xxyz[i] * ga_z[i];

        ts_yyyz_xxzz[i] = 2.0 * ts_yz_xxzz[i] * gfe_0 + ts_yyz_xxzz[i] * ga_y[i];

        ts_yyyz_xyyy[i] = ts_yyy_xyyy[i] * ga_z[i];

        ts_yyyz_xyyz[i] = ts_yyy_xyy[i] * gfe_0 + ts_yyy_xyyz[i] * ga_z[i];

        ts_yyyz_xyzz[i] = 2.0 * ts_yyy_xyz[i] * gfe_0 + ts_yyy_xyzz[i] * ga_z[i];

        ts_yyyz_xzzz[i] = 2.0 * ts_yz_xzzz[i] * gfe_0 + ts_yyz_xzzz[i] * ga_y[i];

        ts_yyyz_yyyy[i] = ts_yyy_yyyy[i] * ga_z[i];

        ts_yyyz_yyyz[i] = ts_yyy_yyy[i] * gfe_0 + ts_yyy_yyyz[i] * ga_z[i];

        ts_yyyz_yyzz[i] = 2.0 * ts_yyy_yyz[i] * gfe_0 + ts_yyy_yyzz[i] * ga_z[i];

        ts_yyyz_yzzz[i] = 3.0 * ts_yyy_yzz[i] * gfe_0 + ts_yyy_yzzz[i] * ga_z[i];

        ts_yyyz_zzzz[i] = 2.0 * ts_yz_zzzz[i] * gfe_0 + ts_yyz_zzzz[i] * ga_y[i];
    }

    // Set up 180-195 components of targeted buffer : GG

    auto ts_yyzz_xxxx = pbuffer.data(idx_gg + 180);

    auto ts_yyzz_xxxy = pbuffer.data(idx_gg + 181);

    auto ts_yyzz_xxxz = pbuffer.data(idx_gg + 182);

    auto ts_yyzz_xxyy = pbuffer.data(idx_gg + 183);

    auto ts_yyzz_xxyz = pbuffer.data(idx_gg + 184);

    auto ts_yyzz_xxzz = pbuffer.data(idx_gg + 185);

    auto ts_yyzz_xyyy = pbuffer.data(idx_gg + 186);

    auto ts_yyzz_xyyz = pbuffer.data(idx_gg + 187);

    auto ts_yyzz_xyzz = pbuffer.data(idx_gg + 188);

    auto ts_yyzz_xzzz = pbuffer.data(idx_gg + 189);

    auto ts_yyzz_yyyy = pbuffer.data(idx_gg + 190);

    auto ts_yyzz_yyyz = pbuffer.data(idx_gg + 191);

    auto ts_yyzz_yyzz = pbuffer.data(idx_gg + 192);

    auto ts_yyzz_yzzz = pbuffer.data(idx_gg + 193);

    auto ts_yyzz_zzzz = pbuffer.data(idx_gg + 194);

    #pragma omp simd aligned(ga_y, ga_z, ts_yy_xxxy, ts_yy_xxyy, ts_yy_xyyy, ts_yy_yyyy, ts_yyz_xxxy, ts_yyz_xxyy, ts_yyz_xyyy, ts_yyz_yyyy, ts_yyzz_xxxx, ts_yyzz_xxxy, ts_yyzz_xxxz, ts_yyzz_xxyy, ts_yyzz_xxyz, ts_yyzz_xxzz, ts_yyzz_xyyy, ts_yyzz_xyyz, ts_yyzz_xyzz, ts_yyzz_xzzz, ts_yyzz_yyyy, ts_yyzz_yyyz, ts_yyzz_yyzz, ts_yyzz_yzzz, ts_yyzz_zzzz, ts_yzz_xxxx, ts_yzz_xxxz, ts_yzz_xxyz, ts_yzz_xxz, ts_yzz_xxzz, ts_yzz_xyyz, ts_yzz_xyz, ts_yzz_xyzz, ts_yzz_xzz, ts_yzz_xzzz, ts_yzz_yyyz, ts_yzz_yyz, ts_yzz_yyzz, ts_yzz_yzz, ts_yzz_yzzz, ts_yzz_zzz, ts_yzz_zzzz, ts_zz_xxxx, ts_zz_xxxz, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyzz_xxxx[i] = ts_zz_xxxx[i] * gfe_0 + ts_yzz_xxxx[i] * ga_y[i];

        ts_yyzz_xxxy[i] = ts_yy_xxxy[i] * gfe_0 + ts_yyz_xxxy[i] * ga_z[i];

        ts_yyzz_xxxz[i] = ts_zz_xxxz[i] * gfe_0 + ts_yzz_xxxz[i] * ga_y[i];

        ts_yyzz_xxyy[i] = ts_yy_xxyy[i] * gfe_0 + ts_yyz_xxyy[i] * ga_z[i];

        ts_yyzz_xxyz[i] = ts_zz_xxyz[i] * gfe_0 + ts_yzz_xxz[i] * gfe_0 + ts_yzz_xxyz[i] * ga_y[i];

        ts_yyzz_xxzz[i] = ts_zz_xxzz[i] * gfe_0 + ts_yzz_xxzz[i] * ga_y[i];

        ts_yyzz_xyyy[i] = ts_yy_xyyy[i] * gfe_0 + ts_yyz_xyyy[i] * ga_z[i];

        ts_yyzz_xyyz[i] = ts_zz_xyyz[i] * gfe_0 + 2.0 * ts_yzz_xyz[i] * gfe_0 + ts_yzz_xyyz[i] * ga_y[i];

        ts_yyzz_xyzz[i] = ts_zz_xyzz[i] * gfe_0 + ts_yzz_xzz[i] * gfe_0 + ts_yzz_xyzz[i] * ga_y[i];

        ts_yyzz_xzzz[i] = ts_zz_xzzz[i] * gfe_0 + ts_yzz_xzzz[i] * ga_y[i];

        ts_yyzz_yyyy[i] = ts_yy_yyyy[i] * gfe_0 + ts_yyz_yyyy[i] * ga_z[i];

        ts_yyzz_yyyz[i] = ts_zz_yyyz[i] * gfe_0 + 3.0 * ts_yzz_yyz[i] * gfe_0 + ts_yzz_yyyz[i] * ga_y[i];

        ts_yyzz_yyzz[i] = ts_zz_yyzz[i] * gfe_0 + 2.0 * ts_yzz_yzz[i] * gfe_0 + ts_yzz_yyzz[i] * ga_y[i];

        ts_yyzz_yzzz[i] = ts_zz_yzzz[i] * gfe_0 + ts_yzz_zzz[i] * gfe_0 + ts_yzz_yzzz[i] * ga_y[i];

        ts_yyzz_zzzz[i] = ts_zz_zzzz[i] * gfe_0 + ts_yzz_zzzz[i] * ga_y[i];
    }

    // Set up 195-210 components of targeted buffer : GG

    auto ts_yzzz_xxxx = pbuffer.data(idx_gg + 195);

    auto ts_yzzz_xxxy = pbuffer.data(idx_gg + 196);

    auto ts_yzzz_xxxz = pbuffer.data(idx_gg + 197);

    auto ts_yzzz_xxyy = pbuffer.data(idx_gg + 198);

    auto ts_yzzz_xxyz = pbuffer.data(idx_gg + 199);

    auto ts_yzzz_xxzz = pbuffer.data(idx_gg + 200);

    auto ts_yzzz_xyyy = pbuffer.data(idx_gg + 201);

    auto ts_yzzz_xyyz = pbuffer.data(idx_gg + 202);

    auto ts_yzzz_xyzz = pbuffer.data(idx_gg + 203);

    auto ts_yzzz_xzzz = pbuffer.data(idx_gg + 204);

    auto ts_yzzz_yyyy = pbuffer.data(idx_gg + 205);

    auto ts_yzzz_yyyz = pbuffer.data(idx_gg + 206);

    auto ts_yzzz_yyzz = pbuffer.data(idx_gg + 207);

    auto ts_yzzz_yzzz = pbuffer.data(idx_gg + 208);

    auto ts_yzzz_zzzz = pbuffer.data(idx_gg + 209);

    #pragma omp simd aligned(ga_y, ts_yzzz_xxxx, ts_yzzz_xxxy, ts_yzzz_xxxz, ts_yzzz_xxyy, ts_yzzz_xxyz, ts_yzzz_xxzz, ts_yzzz_xyyy, ts_yzzz_xyyz, ts_yzzz_xyzz, ts_yzzz_xzzz, ts_yzzz_yyyy, ts_yzzz_yyyz, ts_yzzz_yyzz, ts_yzzz_yzzz, ts_yzzz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yzzz_xxxx[i] = ts_zzz_xxxx[i] * ga_y[i];

        ts_yzzz_xxxy[i] = ts_zzz_xxx[i] * gfe_0 + ts_zzz_xxxy[i] * ga_y[i];

        ts_yzzz_xxxz[i] = ts_zzz_xxxz[i] * ga_y[i];

        ts_yzzz_xxyy[i] = 2.0 * ts_zzz_xxy[i] * gfe_0 + ts_zzz_xxyy[i] * ga_y[i];

        ts_yzzz_xxyz[i] = ts_zzz_xxz[i] * gfe_0 + ts_zzz_xxyz[i] * ga_y[i];

        ts_yzzz_xxzz[i] = ts_zzz_xxzz[i] * ga_y[i];

        ts_yzzz_xyyy[i] = 3.0 * ts_zzz_xyy[i] * gfe_0 + ts_zzz_xyyy[i] * ga_y[i];

        ts_yzzz_xyyz[i] = 2.0 * ts_zzz_xyz[i] * gfe_0 + ts_zzz_xyyz[i] * ga_y[i];

        ts_yzzz_xyzz[i] = ts_zzz_xzz[i] * gfe_0 + ts_zzz_xyzz[i] * ga_y[i];

        ts_yzzz_xzzz[i] = ts_zzz_xzzz[i] * ga_y[i];

        ts_yzzz_yyyy[i] = 4.0 * ts_zzz_yyy[i] * gfe_0 + ts_zzz_yyyy[i] * ga_y[i];

        ts_yzzz_yyyz[i] = 3.0 * ts_zzz_yyz[i] * gfe_0 + ts_zzz_yyyz[i] * ga_y[i];

        ts_yzzz_yyzz[i] = 2.0 * ts_zzz_yzz[i] * gfe_0 + ts_zzz_yyzz[i] * ga_y[i];

        ts_yzzz_yzzz[i] = ts_zzz_zzz[i] * gfe_0 + ts_zzz_yzzz[i] * ga_y[i];

        ts_yzzz_zzzz[i] = ts_zzz_zzzz[i] * ga_y[i];
    }

    // Set up 210-225 components of targeted buffer : GG

    auto ts_zzzz_xxxx = pbuffer.data(idx_gg + 210);

    auto ts_zzzz_xxxy = pbuffer.data(idx_gg + 211);

    auto ts_zzzz_xxxz = pbuffer.data(idx_gg + 212);

    auto ts_zzzz_xxyy = pbuffer.data(idx_gg + 213);

    auto ts_zzzz_xxyz = pbuffer.data(idx_gg + 214);

    auto ts_zzzz_xxzz = pbuffer.data(idx_gg + 215);

    auto ts_zzzz_xyyy = pbuffer.data(idx_gg + 216);

    auto ts_zzzz_xyyz = pbuffer.data(idx_gg + 217);

    auto ts_zzzz_xyzz = pbuffer.data(idx_gg + 218);

    auto ts_zzzz_xzzz = pbuffer.data(idx_gg + 219);

    auto ts_zzzz_yyyy = pbuffer.data(idx_gg + 220);

    auto ts_zzzz_yyyz = pbuffer.data(idx_gg + 221);

    auto ts_zzzz_yyzz = pbuffer.data(idx_gg + 222);

    auto ts_zzzz_yzzz = pbuffer.data(idx_gg + 223);

    auto ts_zzzz_zzzz = pbuffer.data(idx_gg + 224);

    #pragma omp simd aligned(ga_z, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, ts_zzz_xxx, ts_zzz_xxxx, ts_zzz_xxxy, ts_zzz_xxxz, ts_zzz_xxy, ts_zzz_xxyy, ts_zzz_xxyz, ts_zzz_xxz, ts_zzz_xxzz, ts_zzz_xyy, ts_zzz_xyyy, ts_zzz_xyyz, ts_zzz_xyz, ts_zzz_xyzz, ts_zzz_xzz, ts_zzz_xzzz, ts_zzz_yyy, ts_zzz_yyyy, ts_zzz_yyyz, ts_zzz_yyz, ts_zzz_yyzz, ts_zzz_yzz, ts_zzz_yzzz, ts_zzz_zzz, ts_zzz_zzzz, ts_zzzz_xxxx, ts_zzzz_xxxy, ts_zzzz_xxxz, ts_zzzz_xxyy, ts_zzzz_xxyz, ts_zzzz_xxzz, ts_zzzz_xyyy, ts_zzzz_xyyz, ts_zzzz_xyzz, ts_zzzz_xzzz, ts_zzzz_yyyy, ts_zzzz_yyyz, ts_zzzz_yyzz, ts_zzzz_yzzz, ts_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_zzzz_xxxx[i] = 3.0 * ts_zz_xxxx[i] * gfe_0 + ts_zzz_xxxx[i] * ga_z[i];

        ts_zzzz_xxxy[i] = 3.0 * ts_zz_xxxy[i] * gfe_0 + ts_zzz_xxxy[i] * ga_z[i];

        ts_zzzz_xxxz[i] = 3.0 * ts_zz_xxxz[i] * gfe_0 + ts_zzz_xxx[i] * gfe_0 + ts_zzz_xxxz[i] * ga_z[i];

        ts_zzzz_xxyy[i] = 3.0 * ts_zz_xxyy[i] * gfe_0 + ts_zzz_xxyy[i] * ga_z[i];

        ts_zzzz_xxyz[i] = 3.0 * ts_zz_xxyz[i] * gfe_0 + ts_zzz_xxy[i] * gfe_0 + ts_zzz_xxyz[i] * ga_z[i];

        ts_zzzz_xxzz[i] = 3.0 * ts_zz_xxzz[i] * gfe_0 + 2.0 * ts_zzz_xxz[i] * gfe_0 + ts_zzz_xxzz[i] * ga_z[i];

        ts_zzzz_xyyy[i] = 3.0 * ts_zz_xyyy[i] * gfe_0 + ts_zzz_xyyy[i] * ga_z[i];

        ts_zzzz_xyyz[i] = 3.0 * ts_zz_xyyz[i] * gfe_0 + ts_zzz_xyy[i] * gfe_0 + ts_zzz_xyyz[i] * ga_z[i];

        ts_zzzz_xyzz[i] = 3.0 * ts_zz_xyzz[i] * gfe_0 + 2.0 * ts_zzz_xyz[i] * gfe_0 + ts_zzz_xyzz[i] * ga_z[i];

        ts_zzzz_xzzz[i] = 3.0 * ts_zz_xzzz[i] * gfe_0 + 3.0 * ts_zzz_xzz[i] * gfe_0 + ts_zzz_xzzz[i] * ga_z[i];

        ts_zzzz_yyyy[i] = 3.0 * ts_zz_yyyy[i] * gfe_0 + ts_zzz_yyyy[i] * ga_z[i];

        ts_zzzz_yyyz[i] = 3.0 * ts_zz_yyyz[i] * gfe_0 + ts_zzz_yyy[i] * gfe_0 + ts_zzz_yyyz[i] * ga_z[i];

        ts_zzzz_yyzz[i] = 3.0 * ts_zz_yyzz[i] * gfe_0 + 2.0 * ts_zzz_yyz[i] * gfe_0 + ts_zzz_yyzz[i] * ga_z[i];

        ts_zzzz_yzzz[i] = 3.0 * ts_zz_yzzz[i] * gfe_0 + 3.0 * ts_zzz_yzz[i] * gfe_0 + ts_zzz_yzzz[i] * ga_z[i];

        ts_zzzz_zzzz[i] = 3.0 * ts_zz_zzzz[i] * gfe_0 + 4.0 * ts_zzz_zzz[i] * gfe_0 + ts_zzz_zzzz[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

