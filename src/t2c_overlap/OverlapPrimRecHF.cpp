#include "OverlapPrimRecHF.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_hf(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_hf,
                     const size_t              idx_ovl_ff,
                     const size_t              idx_ovl_gd,
                     const size_t              idx_ovl_gf,
                     const CSimdArray<double>& factors,
                     const size_t              idx_rpa,
                     const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : FF

    auto ts_xxx_xxx = pbuffer.data(idx_ovl_ff);

    auto ts_xxx_xxy = pbuffer.data(idx_ovl_ff + 1);

    auto ts_xxx_xxz = pbuffer.data(idx_ovl_ff + 2);

    auto ts_xxx_xyy = pbuffer.data(idx_ovl_ff + 3);

    auto ts_xxx_xyz = pbuffer.data(idx_ovl_ff + 4);

    auto ts_xxx_xzz = pbuffer.data(idx_ovl_ff + 5);

    auto ts_xxx_yyy = pbuffer.data(idx_ovl_ff + 6);

    auto ts_xxx_yyz = pbuffer.data(idx_ovl_ff + 7);

    auto ts_xxx_yzz = pbuffer.data(idx_ovl_ff + 8);

    auto ts_xxx_zzz = pbuffer.data(idx_ovl_ff + 9);

    auto ts_xxy_xxx = pbuffer.data(idx_ovl_ff + 10);

    auto ts_xxy_xxz = pbuffer.data(idx_ovl_ff + 12);

    auto ts_xxy_xzz = pbuffer.data(idx_ovl_ff + 15);

    auto ts_xxy_yyy = pbuffer.data(idx_ovl_ff + 16);

    auto ts_xxy_yyz = pbuffer.data(idx_ovl_ff + 17);

    auto ts_xxy_yzz = pbuffer.data(idx_ovl_ff + 18);

    auto ts_xxz_xxx = pbuffer.data(idx_ovl_ff + 20);

    auto ts_xxz_xxy = pbuffer.data(idx_ovl_ff + 21);

    auto ts_xxz_xxz = pbuffer.data(idx_ovl_ff + 22);

    auto ts_xxz_xyy = pbuffer.data(idx_ovl_ff + 23);

    auto ts_xxz_xzz = pbuffer.data(idx_ovl_ff + 25);

    auto ts_xxz_yyz = pbuffer.data(idx_ovl_ff + 27);

    auto ts_xxz_yzz = pbuffer.data(idx_ovl_ff + 28);

    auto ts_xxz_zzz = pbuffer.data(idx_ovl_ff + 29);

    auto ts_xyy_xxy = pbuffer.data(idx_ovl_ff + 31);

    auto ts_xyy_xyy = pbuffer.data(idx_ovl_ff + 33);

    auto ts_xyy_xyz = pbuffer.data(idx_ovl_ff + 34);

    auto ts_xyy_yyy = pbuffer.data(idx_ovl_ff + 36);

    auto ts_xyy_yyz = pbuffer.data(idx_ovl_ff + 37);

    auto ts_xyy_yzz = pbuffer.data(idx_ovl_ff + 38);

    auto ts_xyy_zzz = pbuffer.data(idx_ovl_ff + 39);

    auto ts_xyz_yyz = pbuffer.data(idx_ovl_ff + 47);

    auto ts_xyz_yzz = pbuffer.data(idx_ovl_ff + 48);

    auto ts_xzz_xxz = pbuffer.data(idx_ovl_ff + 52);

    auto ts_xzz_xyz = pbuffer.data(idx_ovl_ff + 54);

    auto ts_xzz_xzz = pbuffer.data(idx_ovl_ff + 55);

    auto ts_xzz_yyy = pbuffer.data(idx_ovl_ff + 56);

    auto ts_xzz_yyz = pbuffer.data(idx_ovl_ff + 57);

    auto ts_xzz_yzz = pbuffer.data(idx_ovl_ff + 58);

    auto ts_xzz_zzz = pbuffer.data(idx_ovl_ff + 59);

    auto ts_yyy_xxx = pbuffer.data(idx_ovl_ff + 60);

    auto ts_yyy_xxy = pbuffer.data(idx_ovl_ff + 61);

    auto ts_yyy_xxz = pbuffer.data(idx_ovl_ff + 62);

    auto ts_yyy_xyy = pbuffer.data(idx_ovl_ff + 63);

    auto ts_yyy_xyz = pbuffer.data(idx_ovl_ff + 64);

    auto ts_yyy_xzz = pbuffer.data(idx_ovl_ff + 65);

    auto ts_yyy_yyy = pbuffer.data(idx_ovl_ff + 66);

    auto ts_yyy_yyz = pbuffer.data(idx_ovl_ff + 67);

    auto ts_yyy_yzz = pbuffer.data(idx_ovl_ff + 68);

    auto ts_yyy_zzz = pbuffer.data(idx_ovl_ff + 69);

    auto ts_yyz_xxy = pbuffer.data(idx_ovl_ff + 71);

    auto ts_yyz_xxz = pbuffer.data(idx_ovl_ff + 72);

    auto ts_yyz_xyy = pbuffer.data(idx_ovl_ff + 73);

    auto ts_yyz_xzz = pbuffer.data(idx_ovl_ff + 75);

    auto ts_yyz_yyy = pbuffer.data(idx_ovl_ff + 76);

    auto ts_yyz_yyz = pbuffer.data(idx_ovl_ff + 77);

    auto ts_yyz_yzz = pbuffer.data(idx_ovl_ff + 78);

    auto ts_yyz_zzz = pbuffer.data(idx_ovl_ff + 79);

    auto ts_yzz_xxx = pbuffer.data(idx_ovl_ff + 80);

    auto ts_yzz_xxz = pbuffer.data(idx_ovl_ff + 82);

    auto ts_yzz_xyz = pbuffer.data(idx_ovl_ff + 84);

    auto ts_yzz_xzz = pbuffer.data(idx_ovl_ff + 85);

    auto ts_yzz_yyy = pbuffer.data(idx_ovl_ff + 86);

    auto ts_yzz_yyz = pbuffer.data(idx_ovl_ff + 87);

    auto ts_yzz_yzz = pbuffer.data(idx_ovl_ff + 88);

    auto ts_yzz_zzz = pbuffer.data(idx_ovl_ff + 89);

    auto ts_zzz_xxx = pbuffer.data(idx_ovl_ff + 90);

    auto ts_zzz_xxy = pbuffer.data(idx_ovl_ff + 91);

    auto ts_zzz_xxz = pbuffer.data(idx_ovl_ff + 92);

    auto ts_zzz_xyy = pbuffer.data(idx_ovl_ff + 93);

    auto ts_zzz_xyz = pbuffer.data(idx_ovl_ff + 94);

    auto ts_zzz_xzz = pbuffer.data(idx_ovl_ff + 95);

    auto ts_zzz_yyy = pbuffer.data(idx_ovl_ff + 96);

    auto ts_zzz_yyz = pbuffer.data(idx_ovl_ff + 97);

    auto ts_zzz_yzz = pbuffer.data(idx_ovl_ff + 98);

    auto ts_zzz_zzz = pbuffer.data(idx_ovl_ff + 99);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_ovl_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_ovl_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_ovl_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_ovl_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_ovl_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_ovl_gd + 5);

    auto ts_xxxz_xz = pbuffer.data(idx_ovl_gd + 14);

    auto ts_xxyy_xy = pbuffer.data(idx_ovl_gd + 19);

    auto ts_xxyy_yy = pbuffer.data(idx_ovl_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_ovl_gd + 22);

    auto ts_xxzz_xx = pbuffer.data(idx_ovl_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_ovl_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_ovl_gd + 32);

    auto ts_xxzz_yz = pbuffer.data(idx_ovl_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_ovl_gd + 35);

    auto ts_xyyy_xy = pbuffer.data(idx_ovl_gd + 37);

    auto ts_xyyy_yy = pbuffer.data(idx_ovl_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_ovl_gd + 40);

    auto ts_xzzz_xz = pbuffer.data(idx_ovl_gd + 56);

    auto ts_xzzz_yz = pbuffer.data(idx_ovl_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_ovl_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_ovl_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_ovl_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_ovl_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_ovl_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_ovl_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_ovl_gd + 65);

    auto ts_yyyz_xz = pbuffer.data(idx_ovl_gd + 68);

    auto ts_yyyz_yz = pbuffer.data(idx_ovl_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_ovl_gd + 71);

    auto ts_yyzz_xx = pbuffer.data(idx_ovl_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_ovl_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_ovl_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_ovl_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_ovl_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_ovl_gd + 77);

    auto ts_yzzz_xy = pbuffer.data(idx_ovl_gd + 79);

    auto ts_yzzz_xz = pbuffer.data(idx_ovl_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_ovl_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_ovl_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_ovl_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_ovl_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_ovl_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_ovl_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_ovl_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_ovl_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_ovl_gd + 89);

    // Set up components of auxiliary buffer : GF

    auto ts_xxxx_xxx = pbuffer.data(idx_ovl_gf);

    auto ts_xxxx_xxy = pbuffer.data(idx_ovl_gf + 1);

    auto ts_xxxx_xxz = pbuffer.data(idx_ovl_gf + 2);

    auto ts_xxxx_xyy = pbuffer.data(idx_ovl_gf + 3);

    auto ts_xxxx_xyz = pbuffer.data(idx_ovl_gf + 4);

    auto ts_xxxx_xzz = pbuffer.data(idx_ovl_gf + 5);

    auto ts_xxxx_yyy = pbuffer.data(idx_ovl_gf + 6);

    auto ts_xxxx_yyz = pbuffer.data(idx_ovl_gf + 7);

    auto ts_xxxx_yzz = pbuffer.data(idx_ovl_gf + 8);

    auto ts_xxxx_zzz = pbuffer.data(idx_ovl_gf + 9);

    auto ts_xxxy_xxx = pbuffer.data(idx_ovl_gf + 10);

    auto ts_xxxy_xxy = pbuffer.data(idx_ovl_gf + 11);

    auto ts_xxxy_xxz = pbuffer.data(idx_ovl_gf + 12);

    auto ts_xxxy_xyy = pbuffer.data(idx_ovl_gf + 13);

    auto ts_xxxy_xzz = pbuffer.data(idx_ovl_gf + 15);

    auto ts_xxxy_yyy = pbuffer.data(idx_ovl_gf + 16);

    auto ts_xxxy_yyz = pbuffer.data(idx_ovl_gf + 17);

    auto ts_xxxy_yzz = pbuffer.data(idx_ovl_gf + 18);

    auto ts_xxxz_xxx = pbuffer.data(idx_ovl_gf + 20);

    auto ts_xxxz_xxy = pbuffer.data(idx_ovl_gf + 21);

    auto ts_xxxz_xxz = pbuffer.data(idx_ovl_gf + 22);

    auto ts_xxxz_xyy = pbuffer.data(idx_ovl_gf + 23);

    auto ts_xxxz_xyz = pbuffer.data(idx_ovl_gf + 24);

    auto ts_xxxz_xzz = pbuffer.data(idx_ovl_gf + 25);

    auto ts_xxxz_yyz = pbuffer.data(idx_ovl_gf + 27);

    auto ts_xxxz_yzz = pbuffer.data(idx_ovl_gf + 28);

    auto ts_xxxz_zzz = pbuffer.data(idx_ovl_gf + 29);

    auto ts_xxyy_xxx = pbuffer.data(idx_ovl_gf + 30);

    auto ts_xxyy_xxy = pbuffer.data(idx_ovl_gf + 31);

    auto ts_xxyy_xxz = pbuffer.data(idx_ovl_gf + 32);

    auto ts_xxyy_xyy = pbuffer.data(idx_ovl_gf + 33);

    auto ts_xxyy_xyz = pbuffer.data(idx_ovl_gf + 34);

    auto ts_xxyy_xzz = pbuffer.data(idx_ovl_gf + 35);

    auto ts_xxyy_yyy = pbuffer.data(idx_ovl_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_ovl_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_ovl_gf + 38);

    auto ts_xxyy_zzz = pbuffer.data(idx_ovl_gf + 39);

    auto ts_xxyz_xxz = pbuffer.data(idx_ovl_gf + 42);

    auto ts_xxyz_xzz = pbuffer.data(idx_ovl_gf + 45);

    auto ts_xxyz_yyz = pbuffer.data(idx_ovl_gf + 47);

    auto ts_xxyz_yzz = pbuffer.data(idx_ovl_gf + 48);

    auto ts_xxzz_xxx = pbuffer.data(idx_ovl_gf + 50);

    auto ts_xxzz_xxy = pbuffer.data(idx_ovl_gf + 51);

    auto ts_xxzz_xxz = pbuffer.data(idx_ovl_gf + 52);

    auto ts_xxzz_xyy = pbuffer.data(idx_ovl_gf + 53);

    auto ts_xxzz_xyz = pbuffer.data(idx_ovl_gf + 54);

    auto ts_xxzz_xzz = pbuffer.data(idx_ovl_gf + 55);

    auto ts_xxzz_yyy = pbuffer.data(idx_ovl_gf + 56);

    auto ts_xxzz_yyz = pbuffer.data(idx_ovl_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_ovl_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_ovl_gf + 59);

    auto ts_xyyy_xxx = pbuffer.data(idx_ovl_gf + 60);

    auto ts_xyyy_xxy = pbuffer.data(idx_ovl_gf + 61);

    auto ts_xyyy_xyy = pbuffer.data(idx_ovl_gf + 63);

    auto ts_xyyy_xyz = pbuffer.data(idx_ovl_gf + 64);

    auto ts_xyyy_yyy = pbuffer.data(idx_ovl_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_ovl_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_ovl_gf + 68);

    auto ts_xyyy_zzz = pbuffer.data(idx_ovl_gf + 69);

    auto ts_xyyz_yyz = pbuffer.data(idx_ovl_gf + 77);

    auto ts_xyyz_yzz = pbuffer.data(idx_ovl_gf + 78);

    auto ts_xyyz_zzz = pbuffer.data(idx_ovl_gf + 79);

    auto ts_xyzz_yyy = pbuffer.data(idx_ovl_gf + 86);

    auto ts_xyzz_yyz = pbuffer.data(idx_ovl_gf + 87);

    auto ts_xyzz_yzz = pbuffer.data(idx_ovl_gf + 88);

    auto ts_xzzz_xxx = pbuffer.data(idx_ovl_gf + 90);

    auto ts_xzzz_xxz = pbuffer.data(idx_ovl_gf + 92);

    auto ts_xzzz_xyz = pbuffer.data(idx_ovl_gf + 94);

    auto ts_xzzz_xzz = pbuffer.data(idx_ovl_gf + 95);

    auto ts_xzzz_yyy = pbuffer.data(idx_ovl_gf + 96);

    auto ts_xzzz_yyz = pbuffer.data(idx_ovl_gf + 97);

    auto ts_xzzz_yzz = pbuffer.data(idx_ovl_gf + 98);

    auto ts_xzzz_zzz = pbuffer.data(idx_ovl_gf + 99);

    auto ts_yyyy_xxx = pbuffer.data(idx_ovl_gf + 100);

    auto ts_yyyy_xxy = pbuffer.data(idx_ovl_gf + 101);

    auto ts_yyyy_xxz = pbuffer.data(idx_ovl_gf + 102);

    auto ts_yyyy_xyy = pbuffer.data(idx_ovl_gf + 103);

    auto ts_yyyy_xyz = pbuffer.data(idx_ovl_gf + 104);

    auto ts_yyyy_xzz = pbuffer.data(idx_ovl_gf + 105);

    auto ts_yyyy_yyy = pbuffer.data(idx_ovl_gf + 106);

    auto ts_yyyy_yyz = pbuffer.data(idx_ovl_gf + 107);

    auto ts_yyyy_yzz = pbuffer.data(idx_ovl_gf + 108);

    auto ts_yyyy_zzz = pbuffer.data(idx_ovl_gf + 109);

    auto ts_yyyz_xxy = pbuffer.data(idx_ovl_gf + 111);

    auto ts_yyyz_xxz = pbuffer.data(idx_ovl_gf + 112);

    auto ts_yyyz_xyy = pbuffer.data(idx_ovl_gf + 113);

    auto ts_yyyz_xyz = pbuffer.data(idx_ovl_gf + 114);

    auto ts_yyyz_xzz = pbuffer.data(idx_ovl_gf + 115);

    auto ts_yyyz_yyy = pbuffer.data(idx_ovl_gf + 116);

    auto ts_yyyz_yyz = pbuffer.data(idx_ovl_gf + 117);

    auto ts_yyyz_yzz = pbuffer.data(idx_ovl_gf + 118);

    auto ts_yyyz_zzz = pbuffer.data(idx_ovl_gf + 119);

    auto ts_yyzz_xxx = pbuffer.data(idx_ovl_gf + 120);

    auto ts_yyzz_xxy = pbuffer.data(idx_ovl_gf + 121);

    auto ts_yyzz_xxz = pbuffer.data(idx_ovl_gf + 122);

    auto ts_yyzz_xyy = pbuffer.data(idx_ovl_gf + 123);

    auto ts_yyzz_xyz = pbuffer.data(idx_ovl_gf + 124);

    auto ts_yyzz_xzz = pbuffer.data(idx_ovl_gf + 125);

    auto ts_yyzz_yyy = pbuffer.data(idx_ovl_gf + 126);

    auto ts_yyzz_yyz = pbuffer.data(idx_ovl_gf + 127);

    auto ts_yyzz_yzz = pbuffer.data(idx_ovl_gf + 128);

    auto ts_yyzz_zzz = pbuffer.data(idx_ovl_gf + 129);

    auto ts_yzzz_xxx = pbuffer.data(idx_ovl_gf + 130);

    auto ts_yzzz_xxy = pbuffer.data(idx_ovl_gf + 131);

    auto ts_yzzz_xxz = pbuffer.data(idx_ovl_gf + 132);

    auto ts_yzzz_xyy = pbuffer.data(idx_ovl_gf + 133);

    auto ts_yzzz_xyz = pbuffer.data(idx_ovl_gf + 134);

    auto ts_yzzz_xzz = pbuffer.data(idx_ovl_gf + 135);

    auto ts_yzzz_yyy = pbuffer.data(idx_ovl_gf + 136);

    auto ts_yzzz_yyz = pbuffer.data(idx_ovl_gf + 137);

    auto ts_yzzz_yzz = pbuffer.data(idx_ovl_gf + 138);

    auto ts_yzzz_zzz = pbuffer.data(idx_ovl_gf + 139);

    auto ts_zzzz_xxx = pbuffer.data(idx_ovl_gf + 140);

    auto ts_zzzz_xxy = pbuffer.data(idx_ovl_gf + 141);

    auto ts_zzzz_xxz = pbuffer.data(idx_ovl_gf + 142);

    auto ts_zzzz_xyy = pbuffer.data(idx_ovl_gf + 143);

    auto ts_zzzz_xyz = pbuffer.data(idx_ovl_gf + 144);

    auto ts_zzzz_xzz = pbuffer.data(idx_ovl_gf + 145);

    auto ts_zzzz_yyy = pbuffer.data(idx_ovl_gf + 146);

    auto ts_zzzz_yyz = pbuffer.data(idx_ovl_gf + 147);

    auto ts_zzzz_yzz = pbuffer.data(idx_ovl_gf + 148);

    auto ts_zzzz_zzz = pbuffer.data(idx_ovl_gf + 149);

    // Set up 0-10 components of targeted buffer : HF

    auto ts_xxxxx_xxx = pbuffer.data(idx_ovl_hf);

    auto ts_xxxxx_xxy = pbuffer.data(idx_ovl_hf + 1);

    auto ts_xxxxx_xxz = pbuffer.data(idx_ovl_hf + 2);

    auto ts_xxxxx_xyy = pbuffer.data(idx_ovl_hf + 3);

    auto ts_xxxxx_xyz = pbuffer.data(idx_ovl_hf + 4);

    auto ts_xxxxx_xzz = pbuffer.data(idx_ovl_hf + 5);

    auto ts_xxxxx_yyy = pbuffer.data(idx_ovl_hf + 6);

    auto ts_xxxxx_yyz = pbuffer.data(idx_ovl_hf + 7);

    auto ts_xxxxx_yzz = pbuffer.data(idx_ovl_hf + 8);

    auto ts_xxxxx_zzz = pbuffer.data(idx_ovl_hf + 9);

#pragma omp simd aligned(pa_x,             \
                             ts_xxx_xxx,   \
                             ts_xxx_xxy,   \
                             ts_xxx_xxz,   \
                             ts_xxx_xyy,   \
                             ts_xxx_xyz,   \
                             ts_xxx_xzz,   \
                             ts_xxx_yyy,   \
                             ts_xxx_yyz,   \
                             ts_xxx_yzz,   \
                             ts_xxx_zzz,   \
                             ts_xxxx_xx,   \
                             ts_xxxx_xxx,  \
                             ts_xxxx_xxy,  \
                             ts_xxxx_xxz,  \
                             ts_xxxx_xy,   \
                             ts_xxxx_xyy,  \
                             ts_xxxx_xyz,  \
                             ts_xxxx_xz,   \
                             ts_xxxx_xzz,  \
                             ts_xxxx_yy,   \
                             ts_xxxx_yyy,  \
                             ts_xxxx_yyz,  \
                             ts_xxxx_yz,   \
                             ts_xxxx_yzz,  \
                             ts_xxxx_zz,   \
                             ts_xxxx_zzz,  \
                             ts_xxxxx_xxx, \
                             ts_xxxxx_xxy, \
                             ts_xxxxx_xxz, \
                             ts_xxxxx_xyy, \
                             ts_xxxxx_xyz, \
                             ts_xxxxx_xzz, \
                             ts_xxxxx_yyy, \
                             ts_xxxxx_yyz, \
                             ts_xxxxx_yzz, \
                             ts_xxxxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxx_xxx[i] = 4.0 * ts_xxx_xxx[i] * fe_0 + 3.0 * ts_xxxx_xx[i] * fe_0 + ts_xxxx_xxx[i] * pa_x[i];

        ts_xxxxx_xxy[i] = 4.0 * ts_xxx_xxy[i] * fe_0 + 2.0 * ts_xxxx_xy[i] * fe_0 + ts_xxxx_xxy[i] * pa_x[i];

        ts_xxxxx_xxz[i] = 4.0 * ts_xxx_xxz[i] * fe_0 + 2.0 * ts_xxxx_xz[i] * fe_0 + ts_xxxx_xxz[i] * pa_x[i];

        ts_xxxxx_xyy[i] = 4.0 * ts_xxx_xyy[i] * fe_0 + ts_xxxx_yy[i] * fe_0 + ts_xxxx_xyy[i] * pa_x[i];

        ts_xxxxx_xyz[i] = 4.0 * ts_xxx_xyz[i] * fe_0 + ts_xxxx_yz[i] * fe_0 + ts_xxxx_xyz[i] * pa_x[i];

        ts_xxxxx_xzz[i] = 4.0 * ts_xxx_xzz[i] * fe_0 + ts_xxxx_zz[i] * fe_0 + ts_xxxx_xzz[i] * pa_x[i];

        ts_xxxxx_yyy[i] = 4.0 * ts_xxx_yyy[i] * fe_0 + ts_xxxx_yyy[i] * pa_x[i];

        ts_xxxxx_yyz[i] = 4.0 * ts_xxx_yyz[i] * fe_0 + ts_xxxx_yyz[i] * pa_x[i];

        ts_xxxxx_yzz[i] = 4.0 * ts_xxx_yzz[i] * fe_0 + ts_xxxx_yzz[i] * pa_x[i];

        ts_xxxxx_zzz[i] = 4.0 * ts_xxx_zzz[i] * fe_0 + ts_xxxx_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : HF

    auto ts_xxxxy_xxx = pbuffer.data(idx_ovl_hf + 10);

    auto ts_xxxxy_xxy = pbuffer.data(idx_ovl_hf + 11);

    auto ts_xxxxy_xxz = pbuffer.data(idx_ovl_hf + 12);

    auto ts_xxxxy_xyy = pbuffer.data(idx_ovl_hf + 13);

    auto ts_xxxxy_xyz = pbuffer.data(idx_ovl_hf + 14);

    auto ts_xxxxy_xzz = pbuffer.data(idx_ovl_hf + 15);

    auto ts_xxxxy_yyy = pbuffer.data(idx_ovl_hf + 16);

    auto ts_xxxxy_yyz = pbuffer.data(idx_ovl_hf + 17);

    auto ts_xxxxy_yzz = pbuffer.data(idx_ovl_hf + 18);

    auto ts_xxxxy_zzz = pbuffer.data(idx_ovl_hf + 19);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             ts_xxxx_xx,   \
                             ts_xxxx_xxx,  \
                             ts_xxxx_xxy,  \
                             ts_xxxx_xxz,  \
                             ts_xxxx_xy,   \
                             ts_xxxx_xyy,  \
                             ts_xxxx_xyz,  \
                             ts_xxxx_xz,   \
                             ts_xxxx_xzz,  \
                             ts_xxxx_zzz,  \
                             ts_xxxxy_xxx, \
                             ts_xxxxy_xxy, \
                             ts_xxxxy_xxz, \
                             ts_xxxxy_xyy, \
                             ts_xxxxy_xyz, \
                             ts_xxxxy_xzz, \
                             ts_xxxxy_yyy, \
                             ts_xxxxy_yyz, \
                             ts_xxxxy_yzz, \
                             ts_xxxxy_zzz, \
                             ts_xxxy_yyy,  \
                             ts_xxxy_yyz,  \
                             ts_xxxy_yzz,  \
                             ts_xxy_yyy,   \
                             ts_xxy_yyz,   \
                             ts_xxy_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxy_xxx[i] = ts_xxxx_xxx[i] * pa_y[i];

        ts_xxxxy_xxy[i] = ts_xxxx_xx[i] * fe_0 + ts_xxxx_xxy[i] * pa_y[i];

        ts_xxxxy_xxz[i] = ts_xxxx_xxz[i] * pa_y[i];

        ts_xxxxy_xyy[i] = 2.0 * ts_xxxx_xy[i] * fe_0 + ts_xxxx_xyy[i] * pa_y[i];

        ts_xxxxy_xyz[i] = ts_xxxx_xz[i] * fe_0 + ts_xxxx_xyz[i] * pa_y[i];

        ts_xxxxy_xzz[i] = ts_xxxx_xzz[i] * pa_y[i];

        ts_xxxxy_yyy[i] = 3.0 * ts_xxy_yyy[i] * fe_0 + ts_xxxy_yyy[i] * pa_x[i];

        ts_xxxxy_yyz[i] = 3.0 * ts_xxy_yyz[i] * fe_0 + ts_xxxy_yyz[i] * pa_x[i];

        ts_xxxxy_yzz[i] = 3.0 * ts_xxy_yzz[i] * fe_0 + ts_xxxy_yzz[i] * pa_x[i];

        ts_xxxxy_zzz[i] = ts_xxxx_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : HF

    auto ts_xxxxz_xxx = pbuffer.data(idx_ovl_hf + 20);

    auto ts_xxxxz_xxy = pbuffer.data(idx_ovl_hf + 21);

    auto ts_xxxxz_xxz = pbuffer.data(idx_ovl_hf + 22);

    auto ts_xxxxz_xyy = pbuffer.data(idx_ovl_hf + 23);

    auto ts_xxxxz_xyz = pbuffer.data(idx_ovl_hf + 24);

    auto ts_xxxxz_xzz = pbuffer.data(idx_ovl_hf + 25);

    auto ts_xxxxz_yyy = pbuffer.data(idx_ovl_hf + 26);

    auto ts_xxxxz_yyz = pbuffer.data(idx_ovl_hf + 27);

    auto ts_xxxxz_yzz = pbuffer.data(idx_ovl_hf + 28);

    auto ts_xxxxz_zzz = pbuffer.data(idx_ovl_hf + 29);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             ts_xxxx_xx,   \
                             ts_xxxx_xxx,  \
                             ts_xxxx_xxy,  \
                             ts_xxxx_xxz,  \
                             ts_xxxx_xy,   \
                             ts_xxxx_xyy,  \
                             ts_xxxx_xyz,  \
                             ts_xxxx_xz,   \
                             ts_xxxx_xzz,  \
                             ts_xxxx_yyy,  \
                             ts_xxxxz_xxx, \
                             ts_xxxxz_xxy, \
                             ts_xxxxz_xxz, \
                             ts_xxxxz_xyy, \
                             ts_xxxxz_xyz, \
                             ts_xxxxz_xzz, \
                             ts_xxxxz_yyy, \
                             ts_xxxxz_yyz, \
                             ts_xxxxz_yzz, \
                             ts_xxxxz_zzz, \
                             ts_xxxz_yyz,  \
                             ts_xxxz_yzz,  \
                             ts_xxxz_zzz,  \
                             ts_xxz_yyz,   \
                             ts_xxz_yzz,   \
                             ts_xxz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxz_xxx[i] = ts_xxxx_xxx[i] * pa_z[i];

        ts_xxxxz_xxy[i] = ts_xxxx_xxy[i] * pa_z[i];

        ts_xxxxz_xxz[i] = ts_xxxx_xx[i] * fe_0 + ts_xxxx_xxz[i] * pa_z[i];

        ts_xxxxz_xyy[i] = ts_xxxx_xyy[i] * pa_z[i];

        ts_xxxxz_xyz[i] = ts_xxxx_xy[i] * fe_0 + ts_xxxx_xyz[i] * pa_z[i];

        ts_xxxxz_xzz[i] = 2.0 * ts_xxxx_xz[i] * fe_0 + ts_xxxx_xzz[i] * pa_z[i];

        ts_xxxxz_yyy[i] = ts_xxxx_yyy[i] * pa_z[i];

        ts_xxxxz_yyz[i] = 3.0 * ts_xxz_yyz[i] * fe_0 + ts_xxxz_yyz[i] * pa_x[i];

        ts_xxxxz_yzz[i] = 3.0 * ts_xxz_yzz[i] * fe_0 + ts_xxxz_yzz[i] * pa_x[i];

        ts_xxxxz_zzz[i] = 3.0 * ts_xxz_zzz[i] * fe_0 + ts_xxxz_zzz[i] * pa_x[i];
    }

    // Set up 30-40 components of targeted buffer : HF

    auto ts_xxxyy_xxx = pbuffer.data(idx_ovl_hf + 30);

    auto ts_xxxyy_xxy = pbuffer.data(idx_ovl_hf + 31);

    auto ts_xxxyy_xxz = pbuffer.data(idx_ovl_hf + 32);

    auto ts_xxxyy_xyy = pbuffer.data(idx_ovl_hf + 33);

    auto ts_xxxyy_xyz = pbuffer.data(idx_ovl_hf + 34);

    auto ts_xxxyy_xzz = pbuffer.data(idx_ovl_hf + 35);

    auto ts_xxxyy_yyy = pbuffer.data(idx_ovl_hf + 36);

    auto ts_xxxyy_yyz = pbuffer.data(idx_ovl_hf + 37);

    auto ts_xxxyy_yzz = pbuffer.data(idx_ovl_hf + 38);

    auto ts_xxxyy_zzz = pbuffer.data(idx_ovl_hf + 39);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             ts_xxx_xxx,   \
                             ts_xxx_xxz,   \
                             ts_xxx_xzz,   \
                             ts_xxxy_xxx,  \
                             ts_xxxy_xxz,  \
                             ts_xxxy_xzz,  \
                             ts_xxxyy_xxx, \
                             ts_xxxyy_xxy, \
                             ts_xxxyy_xxz, \
                             ts_xxxyy_xyy, \
                             ts_xxxyy_xyz, \
                             ts_xxxyy_xzz, \
                             ts_xxxyy_yyy, \
                             ts_xxxyy_yyz, \
                             ts_xxxyy_yzz, \
                             ts_xxxyy_zzz, \
                             ts_xxyy_xxy,  \
                             ts_xxyy_xy,   \
                             ts_xxyy_xyy,  \
                             ts_xxyy_xyz,  \
                             ts_xxyy_yy,   \
                             ts_xxyy_yyy,  \
                             ts_xxyy_yyz,  \
                             ts_xxyy_yz,   \
                             ts_xxyy_yzz,  \
                             ts_xxyy_zzz,  \
                             ts_xyy_xxy,   \
                             ts_xyy_xyy,   \
                             ts_xyy_xyz,   \
                             ts_xyy_yyy,   \
                             ts_xyy_yyz,   \
                             ts_xyy_yzz,   \
                             ts_xyy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyy_xxx[i] = ts_xxx_xxx[i] * fe_0 + ts_xxxy_xxx[i] * pa_y[i];

        ts_xxxyy_xxy[i] = 2.0 * ts_xyy_xxy[i] * fe_0 + 2.0 * ts_xxyy_xy[i] * fe_0 + ts_xxyy_xxy[i] * pa_x[i];

        ts_xxxyy_xxz[i] = ts_xxx_xxz[i] * fe_0 + ts_xxxy_xxz[i] * pa_y[i];

        ts_xxxyy_xyy[i] = 2.0 * ts_xyy_xyy[i] * fe_0 + ts_xxyy_yy[i] * fe_0 + ts_xxyy_xyy[i] * pa_x[i];

        ts_xxxyy_xyz[i] = 2.0 * ts_xyy_xyz[i] * fe_0 + ts_xxyy_yz[i] * fe_0 + ts_xxyy_xyz[i] * pa_x[i];

        ts_xxxyy_xzz[i] = ts_xxx_xzz[i] * fe_0 + ts_xxxy_xzz[i] * pa_y[i];

        ts_xxxyy_yyy[i] = 2.0 * ts_xyy_yyy[i] * fe_0 + ts_xxyy_yyy[i] * pa_x[i];

        ts_xxxyy_yyz[i] = 2.0 * ts_xyy_yyz[i] * fe_0 + ts_xxyy_yyz[i] * pa_x[i];

        ts_xxxyy_yzz[i] = 2.0 * ts_xyy_yzz[i] * fe_0 + ts_xxyy_yzz[i] * pa_x[i];

        ts_xxxyy_zzz[i] = 2.0 * ts_xyy_zzz[i] * fe_0 + ts_xxyy_zzz[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : HF

    auto ts_xxxyz_xxx = pbuffer.data(idx_ovl_hf + 40);

    auto ts_xxxyz_xxy = pbuffer.data(idx_ovl_hf + 41);

    auto ts_xxxyz_xxz = pbuffer.data(idx_ovl_hf + 42);

    auto ts_xxxyz_xyy = pbuffer.data(idx_ovl_hf + 43);

    auto ts_xxxyz_xyz = pbuffer.data(idx_ovl_hf + 44);

    auto ts_xxxyz_xzz = pbuffer.data(idx_ovl_hf + 45);

    auto ts_xxxyz_yyy = pbuffer.data(idx_ovl_hf + 46);

    auto ts_xxxyz_yyz = pbuffer.data(idx_ovl_hf + 47);

    auto ts_xxxyz_yzz = pbuffer.data(idx_ovl_hf + 48);

    auto ts_xxxyz_zzz = pbuffer.data(idx_ovl_hf + 49);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             ts_xxxy_xxy,  \
                             ts_xxxy_xyy,  \
                             ts_xxxy_yyy,  \
                             ts_xxxyz_xxx, \
                             ts_xxxyz_xxy, \
                             ts_xxxyz_xxz, \
                             ts_xxxyz_xyy, \
                             ts_xxxyz_xyz, \
                             ts_xxxyz_xzz, \
                             ts_xxxyz_yyy, \
                             ts_xxxyz_yyz, \
                             ts_xxxyz_yzz, \
                             ts_xxxyz_zzz, \
                             ts_xxxz_xxx,  \
                             ts_xxxz_xxz,  \
                             ts_xxxz_xyz,  \
                             ts_xxxz_xz,   \
                             ts_xxxz_xzz,  \
                             ts_xxxz_zzz,  \
                             ts_xxyz_yyz,  \
                             ts_xxyz_yzz,  \
                             ts_xyz_yyz,   \
                             ts_xyz_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyz_xxx[i] = ts_xxxz_xxx[i] * pa_y[i];

        ts_xxxyz_xxy[i] = ts_xxxy_xxy[i] * pa_z[i];

        ts_xxxyz_xxz[i] = ts_xxxz_xxz[i] * pa_y[i];

        ts_xxxyz_xyy[i] = ts_xxxy_xyy[i] * pa_z[i];

        ts_xxxyz_xyz[i] = ts_xxxz_xz[i] * fe_0 + ts_xxxz_xyz[i] * pa_y[i];

        ts_xxxyz_xzz[i] = ts_xxxz_xzz[i] * pa_y[i];

        ts_xxxyz_yyy[i] = ts_xxxy_yyy[i] * pa_z[i];

        ts_xxxyz_yyz[i] = 2.0 * ts_xyz_yyz[i] * fe_0 + ts_xxyz_yyz[i] * pa_x[i];

        ts_xxxyz_yzz[i] = 2.0 * ts_xyz_yzz[i] * fe_0 + ts_xxyz_yzz[i] * pa_x[i];

        ts_xxxyz_zzz[i] = ts_xxxz_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : HF

    auto ts_xxxzz_xxx = pbuffer.data(idx_ovl_hf + 50);

    auto ts_xxxzz_xxy = pbuffer.data(idx_ovl_hf + 51);

    auto ts_xxxzz_xxz = pbuffer.data(idx_ovl_hf + 52);

    auto ts_xxxzz_xyy = pbuffer.data(idx_ovl_hf + 53);

    auto ts_xxxzz_xyz = pbuffer.data(idx_ovl_hf + 54);

    auto ts_xxxzz_xzz = pbuffer.data(idx_ovl_hf + 55);

    auto ts_xxxzz_yyy = pbuffer.data(idx_ovl_hf + 56);

    auto ts_xxxzz_yyz = pbuffer.data(idx_ovl_hf + 57);

    auto ts_xxxzz_yzz = pbuffer.data(idx_ovl_hf + 58);

    auto ts_xxxzz_zzz = pbuffer.data(idx_ovl_hf + 59);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             ts_xxx_xxx,   \
                             ts_xxx_xxy,   \
                             ts_xxx_xyy,   \
                             ts_xxxz_xxx,  \
                             ts_xxxz_xxy,  \
                             ts_xxxz_xyy,  \
                             ts_xxxzz_xxx, \
                             ts_xxxzz_xxy, \
                             ts_xxxzz_xxz, \
                             ts_xxxzz_xyy, \
                             ts_xxxzz_xyz, \
                             ts_xxxzz_xzz, \
                             ts_xxxzz_yyy, \
                             ts_xxxzz_yyz, \
                             ts_xxxzz_yzz, \
                             ts_xxxzz_zzz, \
                             ts_xxzz_xxz,  \
                             ts_xxzz_xyz,  \
                             ts_xxzz_xz,   \
                             ts_xxzz_xzz,  \
                             ts_xxzz_yyy,  \
                             ts_xxzz_yyz,  \
                             ts_xxzz_yz,   \
                             ts_xxzz_yzz,  \
                             ts_xxzz_zz,   \
                             ts_xxzz_zzz,  \
                             ts_xzz_xxz,   \
                             ts_xzz_xyz,   \
                             ts_xzz_xzz,   \
                             ts_xzz_yyy,   \
                             ts_xzz_yyz,   \
                             ts_xzz_yzz,   \
                             ts_xzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzz_xxx[i] = ts_xxx_xxx[i] * fe_0 + ts_xxxz_xxx[i] * pa_z[i];

        ts_xxxzz_xxy[i] = ts_xxx_xxy[i] * fe_0 + ts_xxxz_xxy[i] * pa_z[i];

        ts_xxxzz_xxz[i] = 2.0 * ts_xzz_xxz[i] * fe_0 + 2.0 * ts_xxzz_xz[i] * fe_0 + ts_xxzz_xxz[i] * pa_x[i];

        ts_xxxzz_xyy[i] = ts_xxx_xyy[i] * fe_0 + ts_xxxz_xyy[i] * pa_z[i];

        ts_xxxzz_xyz[i] = 2.0 * ts_xzz_xyz[i] * fe_0 + ts_xxzz_yz[i] * fe_0 + ts_xxzz_xyz[i] * pa_x[i];

        ts_xxxzz_xzz[i] = 2.0 * ts_xzz_xzz[i] * fe_0 + ts_xxzz_zz[i] * fe_0 + ts_xxzz_xzz[i] * pa_x[i];

        ts_xxxzz_yyy[i] = 2.0 * ts_xzz_yyy[i] * fe_0 + ts_xxzz_yyy[i] * pa_x[i];

        ts_xxxzz_yyz[i] = 2.0 * ts_xzz_yyz[i] * fe_0 + ts_xxzz_yyz[i] * pa_x[i];

        ts_xxxzz_yzz[i] = 2.0 * ts_xzz_yzz[i] * fe_0 + ts_xxzz_yzz[i] * pa_x[i];

        ts_xxxzz_zzz[i] = 2.0 * ts_xzz_zzz[i] * fe_0 + ts_xxzz_zzz[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : HF

    auto ts_xxyyy_xxx = pbuffer.data(idx_ovl_hf + 60);

    auto ts_xxyyy_xxy = pbuffer.data(idx_ovl_hf + 61);

    auto ts_xxyyy_xxz = pbuffer.data(idx_ovl_hf + 62);

    auto ts_xxyyy_xyy = pbuffer.data(idx_ovl_hf + 63);

    auto ts_xxyyy_xyz = pbuffer.data(idx_ovl_hf + 64);

    auto ts_xxyyy_xzz = pbuffer.data(idx_ovl_hf + 65);

    auto ts_xxyyy_yyy = pbuffer.data(idx_ovl_hf + 66);

    auto ts_xxyyy_yyz = pbuffer.data(idx_ovl_hf + 67);

    auto ts_xxyyy_yzz = pbuffer.data(idx_ovl_hf + 68);

    auto ts_xxyyy_zzz = pbuffer.data(idx_ovl_hf + 69);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             ts_xxy_xxx,   \
                             ts_xxy_xxz,   \
                             ts_xxy_xzz,   \
                             ts_xxyy_xxx,  \
                             ts_xxyy_xxz,  \
                             ts_xxyy_xzz,  \
                             ts_xxyyy_xxx, \
                             ts_xxyyy_xxy, \
                             ts_xxyyy_xxz, \
                             ts_xxyyy_xyy, \
                             ts_xxyyy_xyz, \
                             ts_xxyyy_xzz, \
                             ts_xxyyy_yyy, \
                             ts_xxyyy_yyz, \
                             ts_xxyyy_yzz, \
                             ts_xxyyy_zzz, \
                             ts_xyyy_xxy,  \
                             ts_xyyy_xy,   \
                             ts_xyyy_xyy,  \
                             ts_xyyy_xyz,  \
                             ts_xyyy_yy,   \
                             ts_xyyy_yyy,  \
                             ts_xyyy_yyz,  \
                             ts_xyyy_yz,   \
                             ts_xyyy_yzz,  \
                             ts_xyyy_zzz,  \
                             ts_yyy_xxy,   \
                             ts_yyy_xyy,   \
                             ts_yyy_xyz,   \
                             ts_yyy_yyy,   \
                             ts_yyy_yyz,   \
                             ts_yyy_yzz,   \
                             ts_yyy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyy_xxx[i] = 2.0 * ts_xxy_xxx[i] * fe_0 + ts_xxyy_xxx[i] * pa_y[i];

        ts_xxyyy_xxy[i] = ts_yyy_xxy[i] * fe_0 + 2.0 * ts_xyyy_xy[i] * fe_0 + ts_xyyy_xxy[i] * pa_x[i];

        ts_xxyyy_xxz[i] = 2.0 * ts_xxy_xxz[i] * fe_0 + ts_xxyy_xxz[i] * pa_y[i];

        ts_xxyyy_xyy[i] = ts_yyy_xyy[i] * fe_0 + ts_xyyy_yy[i] * fe_0 + ts_xyyy_xyy[i] * pa_x[i];

        ts_xxyyy_xyz[i] = ts_yyy_xyz[i] * fe_0 + ts_xyyy_yz[i] * fe_0 + ts_xyyy_xyz[i] * pa_x[i];

        ts_xxyyy_xzz[i] = 2.0 * ts_xxy_xzz[i] * fe_0 + ts_xxyy_xzz[i] * pa_y[i];

        ts_xxyyy_yyy[i] = ts_yyy_yyy[i] * fe_0 + ts_xyyy_yyy[i] * pa_x[i];

        ts_xxyyy_yyz[i] = ts_yyy_yyz[i] * fe_0 + ts_xyyy_yyz[i] * pa_x[i];

        ts_xxyyy_yzz[i] = ts_yyy_yzz[i] * fe_0 + ts_xyyy_yzz[i] * pa_x[i];

        ts_xxyyy_zzz[i] = ts_yyy_zzz[i] * fe_0 + ts_xyyy_zzz[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : HF

    auto ts_xxyyz_xxx = pbuffer.data(idx_ovl_hf + 70);

    auto ts_xxyyz_xxy = pbuffer.data(idx_ovl_hf + 71);

    auto ts_xxyyz_xxz = pbuffer.data(idx_ovl_hf + 72);

    auto ts_xxyyz_xyy = pbuffer.data(idx_ovl_hf + 73);

    auto ts_xxyyz_xyz = pbuffer.data(idx_ovl_hf + 74);

    auto ts_xxyyz_xzz = pbuffer.data(idx_ovl_hf + 75);

    auto ts_xxyyz_yyy = pbuffer.data(idx_ovl_hf + 76);

    auto ts_xxyyz_yyz = pbuffer.data(idx_ovl_hf + 77);

    auto ts_xxyyz_yzz = pbuffer.data(idx_ovl_hf + 78);

    auto ts_xxyyz_zzz = pbuffer.data(idx_ovl_hf + 79);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             ts_xxyy_xxx,  \
                             ts_xxyy_xxy,  \
                             ts_xxyy_xy,   \
                             ts_xxyy_xyy,  \
                             ts_xxyy_xyz,  \
                             ts_xxyy_yyy,  \
                             ts_xxyyz_xxx, \
                             ts_xxyyz_xxy, \
                             ts_xxyyz_xxz, \
                             ts_xxyyz_xyy, \
                             ts_xxyyz_xyz, \
                             ts_xxyyz_xzz, \
                             ts_xxyyz_yyy, \
                             ts_xxyyz_yyz, \
                             ts_xxyyz_yzz, \
                             ts_xxyyz_zzz, \
                             ts_xxyz_xxz,  \
                             ts_xxyz_xzz,  \
                             ts_xxz_xxz,   \
                             ts_xxz_xzz,   \
                             ts_xyyz_yyz,  \
                             ts_xyyz_yzz,  \
                             ts_xyyz_zzz,  \
                             ts_yyz_yyz,   \
                             ts_yyz_yzz,   \
                             ts_yyz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyz_xxx[i] = ts_xxyy_xxx[i] * pa_z[i];

        ts_xxyyz_xxy[i] = ts_xxyy_xxy[i] * pa_z[i];

        ts_xxyyz_xxz[i] = ts_xxz_xxz[i] * fe_0 + ts_xxyz_xxz[i] * pa_y[i];

        ts_xxyyz_xyy[i] = ts_xxyy_xyy[i] * pa_z[i];

        ts_xxyyz_xyz[i] = ts_xxyy_xy[i] * fe_0 + ts_xxyy_xyz[i] * pa_z[i];

        ts_xxyyz_xzz[i] = ts_xxz_xzz[i] * fe_0 + ts_xxyz_xzz[i] * pa_y[i];

        ts_xxyyz_yyy[i] = ts_xxyy_yyy[i] * pa_z[i];

        ts_xxyyz_yyz[i] = ts_yyz_yyz[i] * fe_0 + ts_xyyz_yyz[i] * pa_x[i];

        ts_xxyyz_yzz[i] = ts_yyz_yzz[i] * fe_0 + ts_xyyz_yzz[i] * pa_x[i];

        ts_xxyyz_zzz[i] = ts_yyz_zzz[i] * fe_0 + ts_xyyz_zzz[i] * pa_x[i];
    }

    // Set up 80-90 components of targeted buffer : HF

    auto ts_xxyzz_xxx = pbuffer.data(idx_ovl_hf + 80);

    auto ts_xxyzz_xxy = pbuffer.data(idx_ovl_hf + 81);

    auto ts_xxyzz_xxz = pbuffer.data(idx_ovl_hf + 82);

    auto ts_xxyzz_xyy = pbuffer.data(idx_ovl_hf + 83);

    auto ts_xxyzz_xyz = pbuffer.data(idx_ovl_hf + 84);

    auto ts_xxyzz_xzz = pbuffer.data(idx_ovl_hf + 85);

    auto ts_xxyzz_yyy = pbuffer.data(idx_ovl_hf + 86);

    auto ts_xxyzz_yyz = pbuffer.data(idx_ovl_hf + 87);

    auto ts_xxyzz_yzz = pbuffer.data(idx_ovl_hf + 88);

    auto ts_xxyzz_zzz = pbuffer.data(idx_ovl_hf + 89);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             ts_xxyzz_xxx, \
                             ts_xxyzz_xxy, \
                             ts_xxyzz_xxz, \
                             ts_xxyzz_xyy, \
                             ts_xxyzz_xyz, \
                             ts_xxyzz_xzz, \
                             ts_xxyzz_yyy, \
                             ts_xxyzz_yyz, \
                             ts_xxyzz_yzz, \
                             ts_xxyzz_zzz, \
                             ts_xxzz_xx,   \
                             ts_xxzz_xxx,  \
                             ts_xxzz_xxy,  \
                             ts_xxzz_xxz,  \
                             ts_xxzz_xy,   \
                             ts_xxzz_xyy,  \
                             ts_xxzz_xyz,  \
                             ts_xxzz_xz,   \
                             ts_xxzz_xzz,  \
                             ts_xxzz_zzz,  \
                             ts_xyzz_yyy,  \
                             ts_xyzz_yyz,  \
                             ts_xyzz_yzz,  \
                             ts_yzz_yyy,   \
                             ts_yzz_yyz,   \
                             ts_yzz_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzz_xxx[i] = ts_xxzz_xxx[i] * pa_y[i];

        ts_xxyzz_xxy[i] = ts_xxzz_xx[i] * fe_0 + ts_xxzz_xxy[i] * pa_y[i];

        ts_xxyzz_xxz[i] = ts_xxzz_xxz[i] * pa_y[i];

        ts_xxyzz_xyy[i] = 2.0 * ts_xxzz_xy[i] * fe_0 + ts_xxzz_xyy[i] * pa_y[i];

        ts_xxyzz_xyz[i] = ts_xxzz_xz[i] * fe_0 + ts_xxzz_xyz[i] * pa_y[i];

        ts_xxyzz_xzz[i] = ts_xxzz_xzz[i] * pa_y[i];

        ts_xxyzz_yyy[i] = ts_yzz_yyy[i] * fe_0 + ts_xyzz_yyy[i] * pa_x[i];

        ts_xxyzz_yyz[i] = ts_yzz_yyz[i] * fe_0 + ts_xyzz_yyz[i] * pa_x[i];

        ts_xxyzz_yzz[i] = ts_yzz_yzz[i] * fe_0 + ts_xyzz_yzz[i] * pa_x[i];

        ts_xxyzz_zzz[i] = ts_xxzz_zzz[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : HF

    auto ts_xxzzz_xxx = pbuffer.data(idx_ovl_hf + 90);

    auto ts_xxzzz_xxy = pbuffer.data(idx_ovl_hf + 91);

    auto ts_xxzzz_xxz = pbuffer.data(idx_ovl_hf + 92);

    auto ts_xxzzz_xyy = pbuffer.data(idx_ovl_hf + 93);

    auto ts_xxzzz_xyz = pbuffer.data(idx_ovl_hf + 94);

    auto ts_xxzzz_xzz = pbuffer.data(idx_ovl_hf + 95);

    auto ts_xxzzz_yyy = pbuffer.data(idx_ovl_hf + 96);

    auto ts_xxzzz_yyz = pbuffer.data(idx_ovl_hf + 97);

    auto ts_xxzzz_yzz = pbuffer.data(idx_ovl_hf + 98);

    auto ts_xxzzz_zzz = pbuffer.data(idx_ovl_hf + 99);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             ts_xxz_xxx,   \
                             ts_xxz_xxy,   \
                             ts_xxz_xyy,   \
                             ts_xxzz_xxx,  \
                             ts_xxzz_xxy,  \
                             ts_xxzz_xyy,  \
                             ts_xxzzz_xxx, \
                             ts_xxzzz_xxy, \
                             ts_xxzzz_xxz, \
                             ts_xxzzz_xyy, \
                             ts_xxzzz_xyz, \
                             ts_xxzzz_xzz, \
                             ts_xxzzz_yyy, \
                             ts_xxzzz_yyz, \
                             ts_xxzzz_yzz, \
                             ts_xxzzz_zzz, \
                             ts_xzzz_xxz,  \
                             ts_xzzz_xyz,  \
                             ts_xzzz_xz,   \
                             ts_xzzz_xzz,  \
                             ts_xzzz_yyy,  \
                             ts_xzzz_yyz,  \
                             ts_xzzz_yz,   \
                             ts_xzzz_yzz,  \
                             ts_xzzz_zz,   \
                             ts_xzzz_zzz,  \
                             ts_zzz_xxz,   \
                             ts_zzz_xyz,   \
                             ts_zzz_xzz,   \
                             ts_zzz_yyy,   \
                             ts_zzz_yyz,   \
                             ts_zzz_yzz,   \
                             ts_zzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzz_xxx[i] = 2.0 * ts_xxz_xxx[i] * fe_0 + ts_xxzz_xxx[i] * pa_z[i];

        ts_xxzzz_xxy[i] = 2.0 * ts_xxz_xxy[i] * fe_0 + ts_xxzz_xxy[i] * pa_z[i];

        ts_xxzzz_xxz[i] = ts_zzz_xxz[i] * fe_0 + 2.0 * ts_xzzz_xz[i] * fe_0 + ts_xzzz_xxz[i] * pa_x[i];

        ts_xxzzz_xyy[i] = 2.0 * ts_xxz_xyy[i] * fe_0 + ts_xxzz_xyy[i] * pa_z[i];

        ts_xxzzz_xyz[i] = ts_zzz_xyz[i] * fe_0 + ts_xzzz_yz[i] * fe_0 + ts_xzzz_xyz[i] * pa_x[i];

        ts_xxzzz_xzz[i] = ts_zzz_xzz[i] * fe_0 + ts_xzzz_zz[i] * fe_0 + ts_xzzz_xzz[i] * pa_x[i];

        ts_xxzzz_yyy[i] = ts_zzz_yyy[i] * fe_0 + ts_xzzz_yyy[i] * pa_x[i];

        ts_xxzzz_yyz[i] = ts_zzz_yyz[i] * fe_0 + ts_xzzz_yyz[i] * pa_x[i];

        ts_xxzzz_yzz[i] = ts_zzz_yzz[i] * fe_0 + ts_xzzz_yzz[i] * pa_x[i];

        ts_xxzzz_zzz[i] = ts_zzz_zzz[i] * fe_0 + ts_xzzz_zzz[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : HF

    auto ts_xyyyy_xxx = pbuffer.data(idx_ovl_hf + 100);

    auto ts_xyyyy_xxy = pbuffer.data(idx_ovl_hf + 101);

    auto ts_xyyyy_xxz = pbuffer.data(idx_ovl_hf + 102);

    auto ts_xyyyy_xyy = pbuffer.data(idx_ovl_hf + 103);

    auto ts_xyyyy_xyz = pbuffer.data(idx_ovl_hf + 104);

    auto ts_xyyyy_xzz = pbuffer.data(idx_ovl_hf + 105);

    auto ts_xyyyy_yyy = pbuffer.data(idx_ovl_hf + 106);

    auto ts_xyyyy_yyz = pbuffer.data(idx_ovl_hf + 107);

    auto ts_xyyyy_yzz = pbuffer.data(idx_ovl_hf + 108);

    auto ts_xyyyy_zzz = pbuffer.data(idx_ovl_hf + 109);

#pragma omp simd aligned(pa_x,             \
                             ts_xyyyy_xxx, \
                             ts_xyyyy_xxy, \
                             ts_xyyyy_xxz, \
                             ts_xyyyy_xyy, \
                             ts_xyyyy_xyz, \
                             ts_xyyyy_xzz, \
                             ts_xyyyy_yyy, \
                             ts_xyyyy_yyz, \
                             ts_xyyyy_yzz, \
                             ts_xyyyy_zzz, \
                             ts_yyyy_xx,   \
                             ts_yyyy_xxx,  \
                             ts_yyyy_xxy,  \
                             ts_yyyy_xxz,  \
                             ts_yyyy_xy,   \
                             ts_yyyy_xyy,  \
                             ts_yyyy_xyz,  \
                             ts_yyyy_xz,   \
                             ts_yyyy_xzz,  \
                             ts_yyyy_yy,   \
                             ts_yyyy_yyy,  \
                             ts_yyyy_yyz,  \
                             ts_yyyy_yz,   \
                             ts_yyyy_yzz,  \
                             ts_yyyy_zz,   \
                             ts_yyyy_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyy_xxx[i] = 3.0 * ts_yyyy_xx[i] * fe_0 + ts_yyyy_xxx[i] * pa_x[i];

        ts_xyyyy_xxy[i] = 2.0 * ts_yyyy_xy[i] * fe_0 + ts_yyyy_xxy[i] * pa_x[i];

        ts_xyyyy_xxz[i] = 2.0 * ts_yyyy_xz[i] * fe_0 + ts_yyyy_xxz[i] * pa_x[i];

        ts_xyyyy_xyy[i] = ts_yyyy_yy[i] * fe_0 + ts_yyyy_xyy[i] * pa_x[i];

        ts_xyyyy_xyz[i] = ts_yyyy_yz[i] * fe_0 + ts_yyyy_xyz[i] * pa_x[i];

        ts_xyyyy_xzz[i] = ts_yyyy_zz[i] * fe_0 + ts_yyyy_xzz[i] * pa_x[i];

        ts_xyyyy_yyy[i] = ts_yyyy_yyy[i] * pa_x[i];

        ts_xyyyy_yyz[i] = ts_yyyy_yyz[i] * pa_x[i];

        ts_xyyyy_yzz[i] = ts_yyyy_yzz[i] * pa_x[i];

        ts_xyyyy_zzz[i] = ts_yyyy_zzz[i] * pa_x[i];
    }

    // Set up 110-120 components of targeted buffer : HF

    auto ts_xyyyz_xxx = pbuffer.data(idx_ovl_hf + 110);

    auto ts_xyyyz_xxy = pbuffer.data(idx_ovl_hf + 111);

    auto ts_xyyyz_xxz = pbuffer.data(idx_ovl_hf + 112);

    auto ts_xyyyz_xyy = pbuffer.data(idx_ovl_hf + 113);

    auto ts_xyyyz_xyz = pbuffer.data(idx_ovl_hf + 114);

    auto ts_xyyyz_xzz = pbuffer.data(idx_ovl_hf + 115);

    auto ts_xyyyz_yyy = pbuffer.data(idx_ovl_hf + 116);

    auto ts_xyyyz_yyz = pbuffer.data(idx_ovl_hf + 117);

    auto ts_xyyyz_yzz = pbuffer.data(idx_ovl_hf + 118);

    auto ts_xyyyz_zzz = pbuffer.data(idx_ovl_hf + 119);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             ts_xyyy_xxx,  \
                             ts_xyyy_xxy,  \
                             ts_xyyy_xyy,  \
                             ts_xyyyz_xxx, \
                             ts_xyyyz_xxy, \
                             ts_xyyyz_xxz, \
                             ts_xyyyz_xyy, \
                             ts_xyyyz_xyz, \
                             ts_xyyyz_xzz, \
                             ts_xyyyz_yyy, \
                             ts_xyyyz_yyz, \
                             ts_xyyyz_yzz, \
                             ts_xyyyz_zzz, \
                             ts_yyyz_xxz,  \
                             ts_yyyz_xyz,  \
                             ts_yyyz_xz,   \
                             ts_yyyz_xzz,  \
                             ts_yyyz_yyy,  \
                             ts_yyyz_yyz,  \
                             ts_yyyz_yz,   \
                             ts_yyyz_yzz,  \
                             ts_yyyz_zz,   \
                             ts_yyyz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyz_xxx[i] = ts_xyyy_xxx[i] * pa_z[i];

        ts_xyyyz_xxy[i] = ts_xyyy_xxy[i] * pa_z[i];

        ts_xyyyz_xxz[i] = 2.0 * ts_yyyz_xz[i] * fe_0 + ts_yyyz_xxz[i] * pa_x[i];

        ts_xyyyz_xyy[i] = ts_xyyy_xyy[i] * pa_z[i];

        ts_xyyyz_xyz[i] = ts_yyyz_yz[i] * fe_0 + ts_yyyz_xyz[i] * pa_x[i];

        ts_xyyyz_xzz[i] = ts_yyyz_zz[i] * fe_0 + ts_yyyz_xzz[i] * pa_x[i];

        ts_xyyyz_yyy[i] = ts_yyyz_yyy[i] * pa_x[i];

        ts_xyyyz_yyz[i] = ts_yyyz_yyz[i] * pa_x[i];

        ts_xyyyz_yzz[i] = ts_yyyz_yzz[i] * pa_x[i];

        ts_xyyyz_zzz[i] = ts_yyyz_zzz[i] * pa_x[i];
    }

    // Set up 120-130 components of targeted buffer : HF

    auto ts_xyyzz_xxx = pbuffer.data(idx_ovl_hf + 120);

    auto ts_xyyzz_xxy = pbuffer.data(idx_ovl_hf + 121);

    auto ts_xyyzz_xxz = pbuffer.data(idx_ovl_hf + 122);

    auto ts_xyyzz_xyy = pbuffer.data(idx_ovl_hf + 123);

    auto ts_xyyzz_xyz = pbuffer.data(idx_ovl_hf + 124);

    auto ts_xyyzz_xzz = pbuffer.data(idx_ovl_hf + 125);

    auto ts_xyyzz_yyy = pbuffer.data(idx_ovl_hf + 126);

    auto ts_xyyzz_yyz = pbuffer.data(idx_ovl_hf + 127);

    auto ts_xyyzz_yzz = pbuffer.data(idx_ovl_hf + 128);

    auto ts_xyyzz_zzz = pbuffer.data(idx_ovl_hf + 129);

#pragma omp simd aligned(pa_x,             \
                             ts_xyyzz_xxx, \
                             ts_xyyzz_xxy, \
                             ts_xyyzz_xxz, \
                             ts_xyyzz_xyy, \
                             ts_xyyzz_xyz, \
                             ts_xyyzz_xzz, \
                             ts_xyyzz_yyy, \
                             ts_xyyzz_yyz, \
                             ts_xyyzz_yzz, \
                             ts_xyyzz_zzz, \
                             ts_yyzz_xx,   \
                             ts_yyzz_xxx,  \
                             ts_yyzz_xxy,  \
                             ts_yyzz_xxz,  \
                             ts_yyzz_xy,   \
                             ts_yyzz_xyy,  \
                             ts_yyzz_xyz,  \
                             ts_yyzz_xz,   \
                             ts_yyzz_xzz,  \
                             ts_yyzz_yy,   \
                             ts_yyzz_yyy,  \
                             ts_yyzz_yyz,  \
                             ts_yyzz_yz,   \
                             ts_yyzz_yzz,  \
                             ts_yyzz_zz,   \
                             ts_yyzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzz_xxx[i] = 3.0 * ts_yyzz_xx[i] * fe_0 + ts_yyzz_xxx[i] * pa_x[i];

        ts_xyyzz_xxy[i] = 2.0 * ts_yyzz_xy[i] * fe_0 + ts_yyzz_xxy[i] * pa_x[i];

        ts_xyyzz_xxz[i] = 2.0 * ts_yyzz_xz[i] * fe_0 + ts_yyzz_xxz[i] * pa_x[i];

        ts_xyyzz_xyy[i] = ts_yyzz_yy[i] * fe_0 + ts_yyzz_xyy[i] * pa_x[i];

        ts_xyyzz_xyz[i] = ts_yyzz_yz[i] * fe_0 + ts_yyzz_xyz[i] * pa_x[i];

        ts_xyyzz_xzz[i] = ts_yyzz_zz[i] * fe_0 + ts_yyzz_xzz[i] * pa_x[i];

        ts_xyyzz_yyy[i] = ts_yyzz_yyy[i] * pa_x[i];

        ts_xyyzz_yyz[i] = ts_yyzz_yyz[i] * pa_x[i];

        ts_xyyzz_yzz[i] = ts_yyzz_yzz[i] * pa_x[i];

        ts_xyyzz_zzz[i] = ts_yyzz_zzz[i] * pa_x[i];
    }

    // Set up 130-140 components of targeted buffer : HF

    auto ts_xyzzz_xxx = pbuffer.data(idx_ovl_hf + 130);

    auto ts_xyzzz_xxy = pbuffer.data(idx_ovl_hf + 131);

    auto ts_xyzzz_xxz = pbuffer.data(idx_ovl_hf + 132);

    auto ts_xyzzz_xyy = pbuffer.data(idx_ovl_hf + 133);

    auto ts_xyzzz_xyz = pbuffer.data(idx_ovl_hf + 134);

    auto ts_xyzzz_xzz = pbuffer.data(idx_ovl_hf + 135);

    auto ts_xyzzz_yyy = pbuffer.data(idx_ovl_hf + 136);

    auto ts_xyzzz_yyz = pbuffer.data(idx_ovl_hf + 137);

    auto ts_xyzzz_yzz = pbuffer.data(idx_ovl_hf + 138);

    auto ts_xyzzz_zzz = pbuffer.data(idx_ovl_hf + 139);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             ts_xyzzz_xxx, \
                             ts_xyzzz_xxy, \
                             ts_xyzzz_xxz, \
                             ts_xyzzz_xyy, \
                             ts_xyzzz_xyz, \
                             ts_xyzzz_xzz, \
                             ts_xyzzz_yyy, \
                             ts_xyzzz_yyz, \
                             ts_xyzzz_yzz, \
                             ts_xyzzz_zzz, \
                             ts_xzzz_xxx,  \
                             ts_xzzz_xxz,  \
                             ts_xzzz_xzz,  \
                             ts_yzzz_xxy,  \
                             ts_yzzz_xy,   \
                             ts_yzzz_xyy,  \
                             ts_yzzz_xyz,  \
                             ts_yzzz_yy,   \
                             ts_yzzz_yyy,  \
                             ts_yzzz_yyz,  \
                             ts_yzzz_yz,   \
                             ts_yzzz_yzz,  \
                             ts_yzzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzzz_xxx[i] = ts_xzzz_xxx[i] * pa_y[i];

        ts_xyzzz_xxy[i] = 2.0 * ts_yzzz_xy[i] * fe_0 + ts_yzzz_xxy[i] * pa_x[i];

        ts_xyzzz_xxz[i] = ts_xzzz_xxz[i] * pa_y[i];

        ts_xyzzz_xyy[i] = ts_yzzz_yy[i] * fe_0 + ts_yzzz_xyy[i] * pa_x[i];

        ts_xyzzz_xyz[i] = ts_yzzz_yz[i] * fe_0 + ts_yzzz_xyz[i] * pa_x[i];

        ts_xyzzz_xzz[i] = ts_xzzz_xzz[i] * pa_y[i];

        ts_xyzzz_yyy[i] = ts_yzzz_yyy[i] * pa_x[i];

        ts_xyzzz_yyz[i] = ts_yzzz_yyz[i] * pa_x[i];

        ts_xyzzz_yzz[i] = ts_yzzz_yzz[i] * pa_x[i];

        ts_xyzzz_zzz[i] = ts_yzzz_zzz[i] * pa_x[i];
    }

    // Set up 140-150 components of targeted buffer : HF

    auto ts_xzzzz_xxx = pbuffer.data(idx_ovl_hf + 140);

    auto ts_xzzzz_xxy = pbuffer.data(idx_ovl_hf + 141);

    auto ts_xzzzz_xxz = pbuffer.data(idx_ovl_hf + 142);

    auto ts_xzzzz_xyy = pbuffer.data(idx_ovl_hf + 143);

    auto ts_xzzzz_xyz = pbuffer.data(idx_ovl_hf + 144);

    auto ts_xzzzz_xzz = pbuffer.data(idx_ovl_hf + 145);

    auto ts_xzzzz_yyy = pbuffer.data(idx_ovl_hf + 146);

    auto ts_xzzzz_yyz = pbuffer.data(idx_ovl_hf + 147);

    auto ts_xzzzz_yzz = pbuffer.data(idx_ovl_hf + 148);

    auto ts_xzzzz_zzz = pbuffer.data(idx_ovl_hf + 149);

#pragma omp simd aligned(pa_x,             \
                             ts_xzzzz_xxx, \
                             ts_xzzzz_xxy, \
                             ts_xzzzz_xxz, \
                             ts_xzzzz_xyy, \
                             ts_xzzzz_xyz, \
                             ts_xzzzz_xzz, \
                             ts_xzzzz_yyy, \
                             ts_xzzzz_yyz, \
                             ts_xzzzz_yzz, \
                             ts_xzzzz_zzz, \
                             ts_zzzz_xx,   \
                             ts_zzzz_xxx,  \
                             ts_zzzz_xxy,  \
                             ts_zzzz_xxz,  \
                             ts_zzzz_xy,   \
                             ts_zzzz_xyy,  \
                             ts_zzzz_xyz,  \
                             ts_zzzz_xz,   \
                             ts_zzzz_xzz,  \
                             ts_zzzz_yy,   \
                             ts_zzzz_yyy,  \
                             ts_zzzz_yyz,  \
                             ts_zzzz_yz,   \
                             ts_zzzz_yzz,  \
                             ts_zzzz_zz,   \
                             ts_zzzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzz_xxx[i] = 3.0 * ts_zzzz_xx[i] * fe_0 + ts_zzzz_xxx[i] * pa_x[i];

        ts_xzzzz_xxy[i] = 2.0 * ts_zzzz_xy[i] * fe_0 + ts_zzzz_xxy[i] * pa_x[i];

        ts_xzzzz_xxz[i] = 2.0 * ts_zzzz_xz[i] * fe_0 + ts_zzzz_xxz[i] * pa_x[i];

        ts_xzzzz_xyy[i] = ts_zzzz_yy[i] * fe_0 + ts_zzzz_xyy[i] * pa_x[i];

        ts_xzzzz_xyz[i] = ts_zzzz_yz[i] * fe_0 + ts_zzzz_xyz[i] * pa_x[i];

        ts_xzzzz_xzz[i] = ts_zzzz_zz[i] * fe_0 + ts_zzzz_xzz[i] * pa_x[i];

        ts_xzzzz_yyy[i] = ts_zzzz_yyy[i] * pa_x[i];

        ts_xzzzz_yyz[i] = ts_zzzz_yyz[i] * pa_x[i];

        ts_xzzzz_yzz[i] = ts_zzzz_yzz[i] * pa_x[i];

        ts_xzzzz_zzz[i] = ts_zzzz_zzz[i] * pa_x[i];
    }

    // Set up 150-160 components of targeted buffer : HF

    auto ts_yyyyy_xxx = pbuffer.data(idx_ovl_hf + 150);

    auto ts_yyyyy_xxy = pbuffer.data(idx_ovl_hf + 151);

    auto ts_yyyyy_xxz = pbuffer.data(idx_ovl_hf + 152);

    auto ts_yyyyy_xyy = pbuffer.data(idx_ovl_hf + 153);

    auto ts_yyyyy_xyz = pbuffer.data(idx_ovl_hf + 154);

    auto ts_yyyyy_xzz = pbuffer.data(idx_ovl_hf + 155);

    auto ts_yyyyy_yyy = pbuffer.data(idx_ovl_hf + 156);

    auto ts_yyyyy_yyz = pbuffer.data(idx_ovl_hf + 157);

    auto ts_yyyyy_yzz = pbuffer.data(idx_ovl_hf + 158);

    auto ts_yyyyy_zzz = pbuffer.data(idx_ovl_hf + 159);

#pragma omp simd aligned(pa_y,             \
                             ts_yyy_xxx,   \
                             ts_yyy_xxy,   \
                             ts_yyy_xxz,   \
                             ts_yyy_xyy,   \
                             ts_yyy_xyz,   \
                             ts_yyy_xzz,   \
                             ts_yyy_yyy,   \
                             ts_yyy_yyz,   \
                             ts_yyy_yzz,   \
                             ts_yyy_zzz,   \
                             ts_yyyy_xx,   \
                             ts_yyyy_xxx,  \
                             ts_yyyy_xxy,  \
                             ts_yyyy_xxz,  \
                             ts_yyyy_xy,   \
                             ts_yyyy_xyy,  \
                             ts_yyyy_xyz,  \
                             ts_yyyy_xz,   \
                             ts_yyyy_xzz,  \
                             ts_yyyy_yy,   \
                             ts_yyyy_yyy,  \
                             ts_yyyy_yyz,  \
                             ts_yyyy_yz,   \
                             ts_yyyy_yzz,  \
                             ts_yyyy_zz,   \
                             ts_yyyy_zzz,  \
                             ts_yyyyy_xxx, \
                             ts_yyyyy_xxy, \
                             ts_yyyyy_xxz, \
                             ts_yyyyy_xyy, \
                             ts_yyyyy_xyz, \
                             ts_yyyyy_xzz, \
                             ts_yyyyy_yyy, \
                             ts_yyyyy_yyz, \
                             ts_yyyyy_yzz, \
                             ts_yyyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyy_xxx[i] = 4.0 * ts_yyy_xxx[i] * fe_0 + ts_yyyy_xxx[i] * pa_y[i];

        ts_yyyyy_xxy[i] = 4.0 * ts_yyy_xxy[i] * fe_0 + ts_yyyy_xx[i] * fe_0 + ts_yyyy_xxy[i] * pa_y[i];

        ts_yyyyy_xxz[i] = 4.0 * ts_yyy_xxz[i] * fe_0 + ts_yyyy_xxz[i] * pa_y[i];

        ts_yyyyy_xyy[i] = 4.0 * ts_yyy_xyy[i] * fe_0 + 2.0 * ts_yyyy_xy[i] * fe_0 + ts_yyyy_xyy[i] * pa_y[i];

        ts_yyyyy_xyz[i] = 4.0 * ts_yyy_xyz[i] * fe_0 + ts_yyyy_xz[i] * fe_0 + ts_yyyy_xyz[i] * pa_y[i];

        ts_yyyyy_xzz[i] = 4.0 * ts_yyy_xzz[i] * fe_0 + ts_yyyy_xzz[i] * pa_y[i];

        ts_yyyyy_yyy[i] = 4.0 * ts_yyy_yyy[i] * fe_0 + 3.0 * ts_yyyy_yy[i] * fe_0 + ts_yyyy_yyy[i] * pa_y[i];

        ts_yyyyy_yyz[i] = 4.0 * ts_yyy_yyz[i] * fe_0 + 2.0 * ts_yyyy_yz[i] * fe_0 + ts_yyyy_yyz[i] * pa_y[i];

        ts_yyyyy_yzz[i] = 4.0 * ts_yyy_yzz[i] * fe_0 + ts_yyyy_zz[i] * fe_0 + ts_yyyy_yzz[i] * pa_y[i];

        ts_yyyyy_zzz[i] = 4.0 * ts_yyy_zzz[i] * fe_0 + ts_yyyy_zzz[i] * pa_y[i];
    }

    // Set up 160-170 components of targeted buffer : HF

    auto ts_yyyyz_xxx = pbuffer.data(idx_ovl_hf + 160);

    auto ts_yyyyz_xxy = pbuffer.data(idx_ovl_hf + 161);

    auto ts_yyyyz_xxz = pbuffer.data(idx_ovl_hf + 162);

    auto ts_yyyyz_xyy = pbuffer.data(idx_ovl_hf + 163);

    auto ts_yyyyz_xyz = pbuffer.data(idx_ovl_hf + 164);

    auto ts_yyyyz_xzz = pbuffer.data(idx_ovl_hf + 165);

    auto ts_yyyyz_yyy = pbuffer.data(idx_ovl_hf + 166);

    auto ts_yyyyz_yyz = pbuffer.data(idx_ovl_hf + 167);

    auto ts_yyyyz_yzz = pbuffer.data(idx_ovl_hf + 168);

    auto ts_yyyyz_zzz = pbuffer.data(idx_ovl_hf + 169);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             ts_yyyy_xxx,  \
                             ts_yyyy_xxy,  \
                             ts_yyyy_xy,   \
                             ts_yyyy_xyy,  \
                             ts_yyyy_xyz,  \
                             ts_yyyy_yy,   \
                             ts_yyyy_yyy,  \
                             ts_yyyy_yyz,  \
                             ts_yyyy_yz,   \
                             ts_yyyy_yzz,  \
                             ts_yyyyz_xxx, \
                             ts_yyyyz_xxy, \
                             ts_yyyyz_xxz, \
                             ts_yyyyz_xyy, \
                             ts_yyyyz_xyz, \
                             ts_yyyyz_xzz, \
                             ts_yyyyz_yyy, \
                             ts_yyyyz_yyz, \
                             ts_yyyyz_yzz, \
                             ts_yyyyz_zzz, \
                             ts_yyyz_xxz,  \
                             ts_yyyz_xzz,  \
                             ts_yyyz_zzz,  \
                             ts_yyz_xxz,   \
                             ts_yyz_xzz,   \
                             ts_yyz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyz_xxx[i] = ts_yyyy_xxx[i] * pa_z[i];

        ts_yyyyz_xxy[i] = ts_yyyy_xxy[i] * pa_z[i];

        ts_yyyyz_xxz[i] = 3.0 * ts_yyz_xxz[i] * fe_0 + ts_yyyz_xxz[i] * pa_y[i];

        ts_yyyyz_xyy[i] = ts_yyyy_xyy[i] * pa_z[i];

        ts_yyyyz_xyz[i] = ts_yyyy_xy[i] * fe_0 + ts_yyyy_xyz[i] * pa_z[i];

        ts_yyyyz_xzz[i] = 3.0 * ts_yyz_xzz[i] * fe_0 + ts_yyyz_xzz[i] * pa_y[i];

        ts_yyyyz_yyy[i] = ts_yyyy_yyy[i] * pa_z[i];

        ts_yyyyz_yyz[i] = ts_yyyy_yy[i] * fe_0 + ts_yyyy_yyz[i] * pa_z[i];

        ts_yyyyz_yzz[i] = 2.0 * ts_yyyy_yz[i] * fe_0 + ts_yyyy_yzz[i] * pa_z[i];

        ts_yyyyz_zzz[i] = 3.0 * ts_yyz_zzz[i] * fe_0 + ts_yyyz_zzz[i] * pa_y[i];
    }

    // Set up 170-180 components of targeted buffer : HF

    auto ts_yyyzz_xxx = pbuffer.data(idx_ovl_hf + 170);

    auto ts_yyyzz_xxy = pbuffer.data(idx_ovl_hf + 171);

    auto ts_yyyzz_xxz = pbuffer.data(idx_ovl_hf + 172);

    auto ts_yyyzz_xyy = pbuffer.data(idx_ovl_hf + 173);

    auto ts_yyyzz_xyz = pbuffer.data(idx_ovl_hf + 174);

    auto ts_yyyzz_xzz = pbuffer.data(idx_ovl_hf + 175);

    auto ts_yyyzz_yyy = pbuffer.data(idx_ovl_hf + 176);

    auto ts_yyyzz_yyz = pbuffer.data(idx_ovl_hf + 177);

    auto ts_yyyzz_yzz = pbuffer.data(idx_ovl_hf + 178);

    auto ts_yyyzz_zzz = pbuffer.data(idx_ovl_hf + 179);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             ts_yyy_xxy,   \
                             ts_yyy_xyy,   \
                             ts_yyy_yyy,   \
                             ts_yyyz_xxy,  \
                             ts_yyyz_xyy,  \
                             ts_yyyz_yyy,  \
                             ts_yyyzz_xxx, \
                             ts_yyyzz_xxy, \
                             ts_yyyzz_xxz, \
                             ts_yyyzz_xyy, \
                             ts_yyyzz_xyz, \
                             ts_yyyzz_xzz, \
                             ts_yyyzz_yyy, \
                             ts_yyyzz_yyz, \
                             ts_yyyzz_yzz, \
                             ts_yyyzz_zzz, \
                             ts_yyzz_xxx,  \
                             ts_yyzz_xxz,  \
                             ts_yyzz_xyz,  \
                             ts_yyzz_xz,   \
                             ts_yyzz_xzz,  \
                             ts_yyzz_yyz,  \
                             ts_yyzz_yz,   \
                             ts_yyzz_yzz,  \
                             ts_yyzz_zz,   \
                             ts_yyzz_zzz,  \
                             ts_yzz_xxx,   \
                             ts_yzz_xxz,   \
                             ts_yzz_xyz,   \
                             ts_yzz_xzz,   \
                             ts_yzz_yyz,   \
                             ts_yzz_yzz,   \
                             ts_yzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzz_xxx[i] = 2.0 * ts_yzz_xxx[i] * fe_0 + ts_yyzz_xxx[i] * pa_y[i];

        ts_yyyzz_xxy[i] = ts_yyy_xxy[i] * fe_0 + ts_yyyz_xxy[i] * pa_z[i];

        ts_yyyzz_xxz[i] = 2.0 * ts_yzz_xxz[i] * fe_0 + ts_yyzz_xxz[i] * pa_y[i];

        ts_yyyzz_xyy[i] = ts_yyy_xyy[i] * fe_0 + ts_yyyz_xyy[i] * pa_z[i];

        ts_yyyzz_xyz[i] = 2.0 * ts_yzz_xyz[i] * fe_0 + ts_yyzz_xz[i] * fe_0 + ts_yyzz_xyz[i] * pa_y[i];

        ts_yyyzz_xzz[i] = 2.0 * ts_yzz_xzz[i] * fe_0 + ts_yyzz_xzz[i] * pa_y[i];

        ts_yyyzz_yyy[i] = ts_yyy_yyy[i] * fe_0 + ts_yyyz_yyy[i] * pa_z[i];

        ts_yyyzz_yyz[i] = 2.0 * ts_yzz_yyz[i] * fe_0 + 2.0 * ts_yyzz_yz[i] * fe_0 + ts_yyzz_yyz[i] * pa_y[i];

        ts_yyyzz_yzz[i] = 2.0 * ts_yzz_yzz[i] * fe_0 + ts_yyzz_zz[i] * fe_0 + ts_yyzz_yzz[i] * pa_y[i];

        ts_yyyzz_zzz[i] = 2.0 * ts_yzz_zzz[i] * fe_0 + ts_yyzz_zzz[i] * pa_y[i];
    }

    // Set up 180-190 components of targeted buffer : HF

    auto ts_yyzzz_xxx = pbuffer.data(idx_ovl_hf + 180);

    auto ts_yyzzz_xxy = pbuffer.data(idx_ovl_hf + 181);

    auto ts_yyzzz_xxz = pbuffer.data(idx_ovl_hf + 182);

    auto ts_yyzzz_xyy = pbuffer.data(idx_ovl_hf + 183);

    auto ts_yyzzz_xyz = pbuffer.data(idx_ovl_hf + 184);

    auto ts_yyzzz_xzz = pbuffer.data(idx_ovl_hf + 185);

    auto ts_yyzzz_yyy = pbuffer.data(idx_ovl_hf + 186);

    auto ts_yyzzz_yyz = pbuffer.data(idx_ovl_hf + 187);

    auto ts_yyzzz_yzz = pbuffer.data(idx_ovl_hf + 188);

    auto ts_yyzzz_zzz = pbuffer.data(idx_ovl_hf + 189);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             ts_yyz_xxy,   \
                             ts_yyz_xyy,   \
                             ts_yyz_yyy,   \
                             ts_yyzz_xxy,  \
                             ts_yyzz_xyy,  \
                             ts_yyzz_yyy,  \
                             ts_yyzzz_xxx, \
                             ts_yyzzz_xxy, \
                             ts_yyzzz_xxz, \
                             ts_yyzzz_xyy, \
                             ts_yyzzz_xyz, \
                             ts_yyzzz_xzz, \
                             ts_yyzzz_yyy, \
                             ts_yyzzz_yyz, \
                             ts_yyzzz_yzz, \
                             ts_yyzzz_zzz, \
                             ts_yzzz_xxx,  \
                             ts_yzzz_xxz,  \
                             ts_yzzz_xyz,  \
                             ts_yzzz_xz,   \
                             ts_yzzz_xzz,  \
                             ts_yzzz_yyz,  \
                             ts_yzzz_yz,   \
                             ts_yzzz_yzz,  \
                             ts_yzzz_zz,   \
                             ts_yzzz_zzz,  \
                             ts_zzz_xxx,   \
                             ts_zzz_xxz,   \
                             ts_zzz_xyz,   \
                             ts_zzz_xzz,   \
                             ts_zzz_yyz,   \
                             ts_zzz_yzz,   \
                             ts_zzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzz_xxx[i] = ts_zzz_xxx[i] * fe_0 + ts_yzzz_xxx[i] * pa_y[i];

        ts_yyzzz_xxy[i] = 2.0 * ts_yyz_xxy[i] * fe_0 + ts_yyzz_xxy[i] * pa_z[i];

        ts_yyzzz_xxz[i] = ts_zzz_xxz[i] * fe_0 + ts_yzzz_xxz[i] * pa_y[i];

        ts_yyzzz_xyy[i] = 2.0 * ts_yyz_xyy[i] * fe_0 + ts_yyzz_xyy[i] * pa_z[i];

        ts_yyzzz_xyz[i] = ts_zzz_xyz[i] * fe_0 + ts_yzzz_xz[i] * fe_0 + ts_yzzz_xyz[i] * pa_y[i];

        ts_yyzzz_xzz[i] = ts_zzz_xzz[i] * fe_0 + ts_yzzz_xzz[i] * pa_y[i];

        ts_yyzzz_yyy[i] = 2.0 * ts_yyz_yyy[i] * fe_0 + ts_yyzz_yyy[i] * pa_z[i];

        ts_yyzzz_yyz[i] = ts_zzz_yyz[i] * fe_0 + 2.0 * ts_yzzz_yz[i] * fe_0 + ts_yzzz_yyz[i] * pa_y[i];

        ts_yyzzz_yzz[i] = ts_zzz_yzz[i] * fe_0 + ts_yzzz_zz[i] * fe_0 + ts_yzzz_yzz[i] * pa_y[i];

        ts_yyzzz_zzz[i] = ts_zzz_zzz[i] * fe_0 + ts_yzzz_zzz[i] * pa_y[i];
    }

    // Set up 190-200 components of targeted buffer : HF

    auto ts_yzzzz_xxx = pbuffer.data(idx_ovl_hf + 190);

    auto ts_yzzzz_xxy = pbuffer.data(idx_ovl_hf + 191);

    auto ts_yzzzz_xxz = pbuffer.data(idx_ovl_hf + 192);

    auto ts_yzzzz_xyy = pbuffer.data(idx_ovl_hf + 193);

    auto ts_yzzzz_xyz = pbuffer.data(idx_ovl_hf + 194);

    auto ts_yzzzz_xzz = pbuffer.data(idx_ovl_hf + 195);

    auto ts_yzzzz_yyy = pbuffer.data(idx_ovl_hf + 196);

    auto ts_yzzzz_yyz = pbuffer.data(idx_ovl_hf + 197);

    auto ts_yzzzz_yzz = pbuffer.data(idx_ovl_hf + 198);

    auto ts_yzzzz_zzz = pbuffer.data(idx_ovl_hf + 199);

#pragma omp simd aligned(pa_y,             \
                             ts_yzzzz_xxx, \
                             ts_yzzzz_xxy, \
                             ts_yzzzz_xxz, \
                             ts_yzzzz_xyy, \
                             ts_yzzzz_xyz, \
                             ts_yzzzz_xzz, \
                             ts_yzzzz_yyy, \
                             ts_yzzzz_yyz, \
                             ts_yzzzz_yzz, \
                             ts_yzzzz_zzz, \
                             ts_zzzz_xx,   \
                             ts_zzzz_xxx,  \
                             ts_zzzz_xxy,  \
                             ts_zzzz_xxz,  \
                             ts_zzzz_xy,   \
                             ts_zzzz_xyy,  \
                             ts_zzzz_xyz,  \
                             ts_zzzz_xz,   \
                             ts_zzzz_xzz,  \
                             ts_zzzz_yy,   \
                             ts_zzzz_yyy,  \
                             ts_zzzz_yyz,  \
                             ts_zzzz_yz,   \
                             ts_zzzz_yzz,  \
                             ts_zzzz_zz,   \
                             ts_zzzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzz_xxx[i] = ts_zzzz_xxx[i] * pa_y[i];

        ts_yzzzz_xxy[i] = ts_zzzz_xx[i] * fe_0 + ts_zzzz_xxy[i] * pa_y[i];

        ts_yzzzz_xxz[i] = ts_zzzz_xxz[i] * pa_y[i];

        ts_yzzzz_xyy[i] = 2.0 * ts_zzzz_xy[i] * fe_0 + ts_zzzz_xyy[i] * pa_y[i];

        ts_yzzzz_xyz[i] = ts_zzzz_xz[i] * fe_0 + ts_zzzz_xyz[i] * pa_y[i];

        ts_yzzzz_xzz[i] = ts_zzzz_xzz[i] * pa_y[i];

        ts_yzzzz_yyy[i] = 3.0 * ts_zzzz_yy[i] * fe_0 + ts_zzzz_yyy[i] * pa_y[i];

        ts_yzzzz_yyz[i] = 2.0 * ts_zzzz_yz[i] * fe_0 + ts_zzzz_yyz[i] * pa_y[i];

        ts_yzzzz_yzz[i] = ts_zzzz_zz[i] * fe_0 + ts_zzzz_yzz[i] * pa_y[i];

        ts_yzzzz_zzz[i] = ts_zzzz_zzz[i] * pa_y[i];
    }

    // Set up 200-210 components of targeted buffer : HF

    auto ts_zzzzz_xxx = pbuffer.data(idx_ovl_hf + 200);

    auto ts_zzzzz_xxy = pbuffer.data(idx_ovl_hf + 201);

    auto ts_zzzzz_xxz = pbuffer.data(idx_ovl_hf + 202);

    auto ts_zzzzz_xyy = pbuffer.data(idx_ovl_hf + 203);

    auto ts_zzzzz_xyz = pbuffer.data(idx_ovl_hf + 204);

    auto ts_zzzzz_xzz = pbuffer.data(idx_ovl_hf + 205);

    auto ts_zzzzz_yyy = pbuffer.data(idx_ovl_hf + 206);

    auto ts_zzzzz_yyz = pbuffer.data(idx_ovl_hf + 207);

    auto ts_zzzzz_yzz = pbuffer.data(idx_ovl_hf + 208);

    auto ts_zzzzz_zzz = pbuffer.data(idx_ovl_hf + 209);

#pragma omp simd aligned(pa_z,             \
                             ts_zzz_xxx,   \
                             ts_zzz_xxy,   \
                             ts_zzz_xxz,   \
                             ts_zzz_xyy,   \
                             ts_zzz_xyz,   \
                             ts_zzz_xzz,   \
                             ts_zzz_yyy,   \
                             ts_zzz_yyz,   \
                             ts_zzz_yzz,   \
                             ts_zzz_zzz,   \
                             ts_zzzz_xx,   \
                             ts_zzzz_xxx,  \
                             ts_zzzz_xxy,  \
                             ts_zzzz_xxz,  \
                             ts_zzzz_xy,   \
                             ts_zzzz_xyy,  \
                             ts_zzzz_xyz,  \
                             ts_zzzz_xz,   \
                             ts_zzzz_xzz,  \
                             ts_zzzz_yy,   \
                             ts_zzzz_yyy,  \
                             ts_zzzz_yyz,  \
                             ts_zzzz_yz,   \
                             ts_zzzz_yzz,  \
                             ts_zzzz_zz,   \
                             ts_zzzz_zzz,  \
                             ts_zzzzz_xxx, \
                             ts_zzzzz_xxy, \
                             ts_zzzzz_xxz, \
                             ts_zzzzz_xyy, \
                             ts_zzzzz_xyz, \
                             ts_zzzzz_xzz, \
                             ts_zzzzz_yyy, \
                             ts_zzzzz_yyz, \
                             ts_zzzzz_yzz, \
                             ts_zzzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzz_xxx[i] = 4.0 * ts_zzz_xxx[i] * fe_0 + ts_zzzz_xxx[i] * pa_z[i];

        ts_zzzzz_xxy[i] = 4.0 * ts_zzz_xxy[i] * fe_0 + ts_zzzz_xxy[i] * pa_z[i];

        ts_zzzzz_xxz[i] = 4.0 * ts_zzz_xxz[i] * fe_0 + ts_zzzz_xx[i] * fe_0 + ts_zzzz_xxz[i] * pa_z[i];

        ts_zzzzz_xyy[i] = 4.0 * ts_zzz_xyy[i] * fe_0 + ts_zzzz_xyy[i] * pa_z[i];

        ts_zzzzz_xyz[i] = 4.0 * ts_zzz_xyz[i] * fe_0 + ts_zzzz_xy[i] * fe_0 + ts_zzzz_xyz[i] * pa_z[i];

        ts_zzzzz_xzz[i] = 4.0 * ts_zzz_xzz[i] * fe_0 + 2.0 * ts_zzzz_xz[i] * fe_0 + ts_zzzz_xzz[i] * pa_z[i];

        ts_zzzzz_yyy[i] = 4.0 * ts_zzz_yyy[i] * fe_0 + ts_zzzz_yyy[i] * pa_z[i];

        ts_zzzzz_yyz[i] = 4.0 * ts_zzz_yyz[i] * fe_0 + ts_zzzz_yy[i] * fe_0 + ts_zzzz_yyz[i] * pa_z[i];

        ts_zzzzz_yzz[i] = 4.0 * ts_zzz_yzz[i] * fe_0 + 2.0 * ts_zzzz_yz[i] * fe_0 + ts_zzzz_yzz[i] * pa_z[i];

        ts_zzzzz_zzz[i] = 4.0 * ts_zzz_zzz[i] * fe_0 + 3.0 * ts_zzzz_zz[i] * fe_0 + ts_zzzz_zzz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
