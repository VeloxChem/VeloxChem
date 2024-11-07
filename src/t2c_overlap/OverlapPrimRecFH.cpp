#include "OverlapPrimRecFH.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_fh(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_fh,
                     const size_t              idx_ovl_ph,
                     const size_t              idx_ovl_dg,
                     const size_t              idx_ovl_dh,
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

    // Set up components of auxiliary buffer : PH

    auto ts_x_xxxxx = pbuffer.data(idx_ovl_ph);

    auto ts_x_xxxxy = pbuffer.data(idx_ovl_ph + 1);

    auto ts_x_xxxxz = pbuffer.data(idx_ovl_ph + 2);

    auto ts_x_xxxyy = pbuffer.data(idx_ovl_ph + 3);

    auto ts_x_xxxyz = pbuffer.data(idx_ovl_ph + 4);

    auto ts_x_xxxzz = pbuffer.data(idx_ovl_ph + 5);

    auto ts_x_xxyyy = pbuffer.data(idx_ovl_ph + 6);

    auto ts_x_xxyyz = pbuffer.data(idx_ovl_ph + 7);

    auto ts_x_xxyzz = pbuffer.data(idx_ovl_ph + 8);

    auto ts_x_xxzzz = pbuffer.data(idx_ovl_ph + 9);

    auto ts_x_xyyyy = pbuffer.data(idx_ovl_ph + 10);

    auto ts_x_xyyyz = pbuffer.data(idx_ovl_ph + 11);

    auto ts_x_xyyzz = pbuffer.data(idx_ovl_ph + 12);

    auto ts_x_xyzzz = pbuffer.data(idx_ovl_ph + 13);

    auto ts_x_xzzzz = pbuffer.data(idx_ovl_ph + 14);

    auto ts_x_yyyyy = pbuffer.data(idx_ovl_ph + 15);

    auto ts_x_yyyyz = pbuffer.data(idx_ovl_ph + 16);

    auto ts_x_yyyzz = pbuffer.data(idx_ovl_ph + 17);

    auto ts_x_yyzzz = pbuffer.data(idx_ovl_ph + 18);

    auto ts_x_yzzzz = pbuffer.data(idx_ovl_ph + 19);

    auto ts_x_zzzzz = pbuffer.data(idx_ovl_ph + 20);

    auto ts_y_xxxxx = pbuffer.data(idx_ovl_ph + 21);

    auto ts_y_xxxxy = pbuffer.data(idx_ovl_ph + 22);

    auto ts_y_xxxxz = pbuffer.data(idx_ovl_ph + 23);

    auto ts_y_xxxyy = pbuffer.data(idx_ovl_ph + 24);

    auto ts_y_xxxyz = pbuffer.data(idx_ovl_ph + 25);

    auto ts_y_xxxzz = pbuffer.data(idx_ovl_ph + 26);

    auto ts_y_xxyyy = pbuffer.data(idx_ovl_ph + 27);

    auto ts_y_xxyyz = pbuffer.data(idx_ovl_ph + 28);

    auto ts_y_xxyzz = pbuffer.data(idx_ovl_ph + 29);

    auto ts_y_xxzzz = pbuffer.data(idx_ovl_ph + 30);

    auto ts_y_xyyyy = pbuffer.data(idx_ovl_ph + 31);

    auto ts_y_xyyyz = pbuffer.data(idx_ovl_ph + 32);

    auto ts_y_xyyzz = pbuffer.data(idx_ovl_ph + 33);

    auto ts_y_xyzzz = pbuffer.data(idx_ovl_ph + 34);

    auto ts_y_xzzzz = pbuffer.data(idx_ovl_ph + 35);

    auto ts_y_yyyyy = pbuffer.data(idx_ovl_ph + 36);

    auto ts_y_yyyyz = pbuffer.data(idx_ovl_ph + 37);

    auto ts_y_yyyzz = pbuffer.data(idx_ovl_ph + 38);

    auto ts_y_yyzzz = pbuffer.data(idx_ovl_ph + 39);

    auto ts_y_yzzzz = pbuffer.data(idx_ovl_ph + 40);

    auto ts_y_zzzzz = pbuffer.data(idx_ovl_ph + 41);

    auto ts_z_xxxxx = pbuffer.data(idx_ovl_ph + 42);

    auto ts_z_xxxxy = pbuffer.data(idx_ovl_ph + 43);

    auto ts_z_xxxxz = pbuffer.data(idx_ovl_ph + 44);

    auto ts_z_xxxyy = pbuffer.data(idx_ovl_ph + 45);

    auto ts_z_xxxyz = pbuffer.data(idx_ovl_ph + 46);

    auto ts_z_xxxzz = pbuffer.data(idx_ovl_ph + 47);

    auto ts_z_xxyyy = pbuffer.data(idx_ovl_ph + 48);

    auto ts_z_xxyyz = pbuffer.data(idx_ovl_ph + 49);

    auto ts_z_xxyzz = pbuffer.data(idx_ovl_ph + 50);

    auto ts_z_xxzzz = pbuffer.data(idx_ovl_ph + 51);

    auto ts_z_xyyyy = pbuffer.data(idx_ovl_ph + 52);

    auto ts_z_xyyyz = pbuffer.data(idx_ovl_ph + 53);

    auto ts_z_xyyzz = pbuffer.data(idx_ovl_ph + 54);

    auto ts_z_xyzzz = pbuffer.data(idx_ovl_ph + 55);

    auto ts_z_xzzzz = pbuffer.data(idx_ovl_ph + 56);

    auto ts_z_yyyyy = pbuffer.data(idx_ovl_ph + 57);

    auto ts_z_yyyyz = pbuffer.data(idx_ovl_ph + 58);

    auto ts_z_yyyzz = pbuffer.data(idx_ovl_ph + 59);

    auto ts_z_yyzzz = pbuffer.data(idx_ovl_ph + 60);

    auto ts_z_yzzzz = pbuffer.data(idx_ovl_ph + 61);

    auto ts_z_zzzzz = pbuffer.data(idx_ovl_ph + 62);

    // Set up components of auxiliary buffer : DG

    auto ts_xx_xxxx = pbuffer.data(idx_ovl_dg);

    auto ts_xx_xxxy = pbuffer.data(idx_ovl_dg + 1);

    auto ts_xx_xxxz = pbuffer.data(idx_ovl_dg + 2);

    auto ts_xx_xxyy = pbuffer.data(idx_ovl_dg + 3);

    auto ts_xx_xxyz = pbuffer.data(idx_ovl_dg + 4);

    auto ts_xx_xxzz = pbuffer.data(idx_ovl_dg + 5);

    auto ts_xx_xyyy = pbuffer.data(idx_ovl_dg + 6);

    auto ts_xx_xyyz = pbuffer.data(idx_ovl_dg + 7);

    auto ts_xx_xyzz = pbuffer.data(idx_ovl_dg + 8);

    auto ts_xx_xzzz = pbuffer.data(idx_ovl_dg + 9);

    auto ts_xx_yyyy = pbuffer.data(idx_ovl_dg + 10);

    auto ts_xx_yyyz = pbuffer.data(idx_ovl_dg + 11);

    auto ts_xx_yyzz = pbuffer.data(idx_ovl_dg + 12);

    auto ts_xx_yzzz = pbuffer.data(idx_ovl_dg + 13);

    auto ts_xx_zzzz = pbuffer.data(idx_ovl_dg + 14);

    auto ts_yy_xxxx = pbuffer.data(idx_ovl_dg + 45);

    auto ts_yy_xxxy = pbuffer.data(idx_ovl_dg + 46);

    auto ts_yy_xxxz = pbuffer.data(idx_ovl_dg + 47);

    auto ts_yy_xxyy = pbuffer.data(idx_ovl_dg + 48);

    auto ts_yy_xxyz = pbuffer.data(idx_ovl_dg + 49);

    auto ts_yy_xxzz = pbuffer.data(idx_ovl_dg + 50);

    auto ts_yy_xyyy = pbuffer.data(idx_ovl_dg + 51);

    auto ts_yy_xyyz = pbuffer.data(idx_ovl_dg + 52);

    auto ts_yy_xyzz = pbuffer.data(idx_ovl_dg + 53);

    auto ts_yy_xzzz = pbuffer.data(idx_ovl_dg + 54);

    auto ts_yy_yyyy = pbuffer.data(idx_ovl_dg + 55);

    auto ts_yy_yyyz = pbuffer.data(idx_ovl_dg + 56);

    auto ts_yy_yyzz = pbuffer.data(idx_ovl_dg + 57);

    auto ts_yy_yzzz = pbuffer.data(idx_ovl_dg + 58);

    auto ts_yy_zzzz = pbuffer.data(idx_ovl_dg + 59);

    auto ts_yz_xxyz = pbuffer.data(idx_ovl_dg + 64);

    auto ts_yz_xyyz = pbuffer.data(idx_ovl_dg + 67);

    auto ts_yz_xyzz = pbuffer.data(idx_ovl_dg + 68);

    auto ts_yz_yyyz = pbuffer.data(idx_ovl_dg + 71);

    auto ts_yz_yyzz = pbuffer.data(idx_ovl_dg + 72);

    auto ts_yz_yzzz = pbuffer.data(idx_ovl_dg + 73);

    auto ts_zz_xxxx = pbuffer.data(idx_ovl_dg + 75);

    auto ts_zz_xxxy = pbuffer.data(idx_ovl_dg + 76);

    auto ts_zz_xxxz = pbuffer.data(idx_ovl_dg + 77);

    auto ts_zz_xxyy = pbuffer.data(idx_ovl_dg + 78);

    auto ts_zz_xxyz = pbuffer.data(idx_ovl_dg + 79);

    auto ts_zz_xxzz = pbuffer.data(idx_ovl_dg + 80);

    auto ts_zz_xyyy = pbuffer.data(idx_ovl_dg + 81);

    auto ts_zz_xyyz = pbuffer.data(idx_ovl_dg + 82);

    auto ts_zz_xyzz = pbuffer.data(idx_ovl_dg + 83);

    auto ts_zz_xzzz = pbuffer.data(idx_ovl_dg + 84);

    auto ts_zz_yyyy = pbuffer.data(idx_ovl_dg + 85);

    auto ts_zz_yyyz = pbuffer.data(idx_ovl_dg + 86);

    auto ts_zz_yyzz = pbuffer.data(idx_ovl_dg + 87);

    auto ts_zz_yzzz = pbuffer.data(idx_ovl_dg + 88);

    auto ts_zz_zzzz = pbuffer.data(idx_ovl_dg + 89);

    // Set up components of auxiliary buffer : DH

    auto ts_xx_xxxxx = pbuffer.data(idx_ovl_dh);

    auto ts_xx_xxxxy = pbuffer.data(idx_ovl_dh + 1);

    auto ts_xx_xxxxz = pbuffer.data(idx_ovl_dh + 2);

    auto ts_xx_xxxyy = pbuffer.data(idx_ovl_dh + 3);

    auto ts_xx_xxxyz = pbuffer.data(idx_ovl_dh + 4);

    auto ts_xx_xxxzz = pbuffer.data(idx_ovl_dh + 5);

    auto ts_xx_xxyyy = pbuffer.data(idx_ovl_dh + 6);

    auto ts_xx_xxyyz = pbuffer.data(idx_ovl_dh + 7);

    auto ts_xx_xxyzz = pbuffer.data(idx_ovl_dh + 8);

    auto ts_xx_xxzzz = pbuffer.data(idx_ovl_dh + 9);

    auto ts_xx_xyyyy = pbuffer.data(idx_ovl_dh + 10);

    auto ts_xx_xyyyz = pbuffer.data(idx_ovl_dh + 11);

    auto ts_xx_xyyzz = pbuffer.data(idx_ovl_dh + 12);

    auto ts_xx_xyzzz = pbuffer.data(idx_ovl_dh + 13);

    auto ts_xx_xzzzz = pbuffer.data(idx_ovl_dh + 14);

    auto ts_xx_yyyyy = pbuffer.data(idx_ovl_dh + 15);

    auto ts_xx_yyyyz = pbuffer.data(idx_ovl_dh + 16);

    auto ts_xx_yyyzz = pbuffer.data(idx_ovl_dh + 17);

    auto ts_xx_yyzzz = pbuffer.data(idx_ovl_dh + 18);

    auto ts_xx_yzzzz = pbuffer.data(idx_ovl_dh + 19);

    auto ts_xx_zzzzz = pbuffer.data(idx_ovl_dh + 20);

    auto ts_xy_xxxxy = pbuffer.data(idx_ovl_dh + 22);

    auto ts_xy_xxxyy = pbuffer.data(idx_ovl_dh + 24);

    auto ts_xy_xxyyy = pbuffer.data(idx_ovl_dh + 27);

    auto ts_xy_xyyyy = pbuffer.data(idx_ovl_dh + 31);

    auto ts_xy_yyyyy = pbuffer.data(idx_ovl_dh + 36);

    auto ts_xy_yyyyz = pbuffer.data(idx_ovl_dh + 37);

    auto ts_xy_yyyzz = pbuffer.data(idx_ovl_dh + 38);

    auto ts_xy_yyzzz = pbuffer.data(idx_ovl_dh + 39);

    auto ts_xy_yzzzz = pbuffer.data(idx_ovl_dh + 40);

    auto ts_xz_xxxxx = pbuffer.data(idx_ovl_dh + 42);

    auto ts_xz_xxxxz = pbuffer.data(idx_ovl_dh + 44);

    auto ts_xz_xxxzz = pbuffer.data(idx_ovl_dh + 47);

    auto ts_xz_xxzzz = pbuffer.data(idx_ovl_dh + 51);

    auto ts_xz_xzzzz = pbuffer.data(idx_ovl_dh + 56);

    auto ts_xz_yyyyz = pbuffer.data(idx_ovl_dh + 58);

    auto ts_xz_yyyzz = pbuffer.data(idx_ovl_dh + 59);

    auto ts_xz_yyzzz = pbuffer.data(idx_ovl_dh + 60);

    auto ts_xz_yzzzz = pbuffer.data(idx_ovl_dh + 61);

    auto ts_xz_zzzzz = pbuffer.data(idx_ovl_dh + 62);

    auto ts_yy_xxxxx = pbuffer.data(idx_ovl_dh + 63);

    auto ts_yy_xxxxy = pbuffer.data(idx_ovl_dh + 64);

    auto ts_yy_xxxxz = pbuffer.data(idx_ovl_dh + 65);

    auto ts_yy_xxxyy = pbuffer.data(idx_ovl_dh + 66);

    auto ts_yy_xxxyz = pbuffer.data(idx_ovl_dh + 67);

    auto ts_yy_xxxzz = pbuffer.data(idx_ovl_dh + 68);

    auto ts_yy_xxyyy = pbuffer.data(idx_ovl_dh + 69);

    auto ts_yy_xxyyz = pbuffer.data(idx_ovl_dh + 70);

    auto ts_yy_xxyzz = pbuffer.data(idx_ovl_dh + 71);

    auto ts_yy_xxzzz = pbuffer.data(idx_ovl_dh + 72);

    auto ts_yy_xyyyy = pbuffer.data(idx_ovl_dh + 73);

    auto ts_yy_xyyyz = pbuffer.data(idx_ovl_dh + 74);

    auto ts_yy_xyyzz = pbuffer.data(idx_ovl_dh + 75);

    auto ts_yy_xyzzz = pbuffer.data(idx_ovl_dh + 76);

    auto ts_yy_xzzzz = pbuffer.data(idx_ovl_dh + 77);

    auto ts_yy_yyyyy = pbuffer.data(idx_ovl_dh + 78);

    auto ts_yy_yyyyz = pbuffer.data(idx_ovl_dh + 79);

    auto ts_yy_yyyzz = pbuffer.data(idx_ovl_dh + 80);

    auto ts_yy_yyzzz = pbuffer.data(idx_ovl_dh + 81);

    auto ts_yy_yzzzz = pbuffer.data(idx_ovl_dh + 82);

    auto ts_yy_zzzzz = pbuffer.data(idx_ovl_dh + 83);

    auto ts_yz_xxxxz = pbuffer.data(idx_ovl_dh + 86);

    auto ts_yz_xxxyz = pbuffer.data(idx_ovl_dh + 88);

    auto ts_yz_xxxzz = pbuffer.data(idx_ovl_dh + 89);

    auto ts_yz_xxyyz = pbuffer.data(idx_ovl_dh + 91);

    auto ts_yz_xxyzz = pbuffer.data(idx_ovl_dh + 92);

    auto ts_yz_xxzzz = pbuffer.data(idx_ovl_dh + 93);

    auto ts_yz_xyyyz = pbuffer.data(idx_ovl_dh + 95);

    auto ts_yz_xyyzz = pbuffer.data(idx_ovl_dh + 96);

    auto ts_yz_xyzzz = pbuffer.data(idx_ovl_dh + 97);

    auto ts_yz_xzzzz = pbuffer.data(idx_ovl_dh + 98);

    auto ts_yz_yyyyy = pbuffer.data(idx_ovl_dh + 99);

    auto ts_yz_yyyyz = pbuffer.data(idx_ovl_dh + 100);

    auto ts_yz_yyyzz = pbuffer.data(idx_ovl_dh + 101);

    auto ts_yz_yyzzz = pbuffer.data(idx_ovl_dh + 102);

    auto ts_yz_yzzzz = pbuffer.data(idx_ovl_dh + 103);

    auto ts_yz_zzzzz = pbuffer.data(idx_ovl_dh + 104);

    auto ts_zz_xxxxx = pbuffer.data(idx_ovl_dh + 105);

    auto ts_zz_xxxxy = pbuffer.data(idx_ovl_dh + 106);

    auto ts_zz_xxxxz = pbuffer.data(idx_ovl_dh + 107);

    auto ts_zz_xxxyy = pbuffer.data(idx_ovl_dh + 108);

    auto ts_zz_xxxyz = pbuffer.data(idx_ovl_dh + 109);

    auto ts_zz_xxxzz = pbuffer.data(idx_ovl_dh + 110);

    auto ts_zz_xxyyy = pbuffer.data(idx_ovl_dh + 111);

    auto ts_zz_xxyyz = pbuffer.data(idx_ovl_dh + 112);

    auto ts_zz_xxyzz = pbuffer.data(idx_ovl_dh + 113);

    auto ts_zz_xxzzz = pbuffer.data(idx_ovl_dh + 114);

    auto ts_zz_xyyyy = pbuffer.data(idx_ovl_dh + 115);

    auto ts_zz_xyyyz = pbuffer.data(idx_ovl_dh + 116);

    auto ts_zz_xyyzz = pbuffer.data(idx_ovl_dh + 117);

    auto ts_zz_xyzzz = pbuffer.data(idx_ovl_dh + 118);

    auto ts_zz_xzzzz = pbuffer.data(idx_ovl_dh + 119);

    auto ts_zz_yyyyy = pbuffer.data(idx_ovl_dh + 120);

    auto ts_zz_yyyyz = pbuffer.data(idx_ovl_dh + 121);

    auto ts_zz_yyyzz = pbuffer.data(idx_ovl_dh + 122);

    auto ts_zz_yyzzz = pbuffer.data(idx_ovl_dh + 123);

    auto ts_zz_yzzzz = pbuffer.data(idx_ovl_dh + 124);

    auto ts_zz_zzzzz = pbuffer.data(idx_ovl_dh + 125);

    // Set up 0-21 components of targeted buffer : FH

    auto ts_xxx_xxxxx = pbuffer.data(idx_ovl_fh);

    auto ts_xxx_xxxxy = pbuffer.data(idx_ovl_fh + 1);

    auto ts_xxx_xxxxz = pbuffer.data(idx_ovl_fh + 2);

    auto ts_xxx_xxxyy = pbuffer.data(idx_ovl_fh + 3);

    auto ts_xxx_xxxyz = pbuffer.data(idx_ovl_fh + 4);

    auto ts_xxx_xxxzz = pbuffer.data(idx_ovl_fh + 5);

    auto ts_xxx_xxyyy = pbuffer.data(idx_ovl_fh + 6);

    auto ts_xxx_xxyyz = pbuffer.data(idx_ovl_fh + 7);

    auto ts_xxx_xxyzz = pbuffer.data(idx_ovl_fh + 8);

    auto ts_xxx_xxzzz = pbuffer.data(idx_ovl_fh + 9);

    auto ts_xxx_xyyyy = pbuffer.data(idx_ovl_fh + 10);

    auto ts_xxx_xyyyz = pbuffer.data(idx_ovl_fh + 11);

    auto ts_xxx_xyyzz = pbuffer.data(idx_ovl_fh + 12);

    auto ts_xxx_xyzzz = pbuffer.data(idx_ovl_fh + 13);

    auto ts_xxx_xzzzz = pbuffer.data(idx_ovl_fh + 14);

    auto ts_xxx_yyyyy = pbuffer.data(idx_ovl_fh + 15);

    auto ts_xxx_yyyyz = pbuffer.data(idx_ovl_fh + 16);

    auto ts_xxx_yyyzz = pbuffer.data(idx_ovl_fh + 17);

    auto ts_xxx_yyzzz = pbuffer.data(idx_ovl_fh + 18);

    auto ts_xxx_yzzzz = pbuffer.data(idx_ovl_fh + 19);

    auto ts_xxx_zzzzz = pbuffer.data(idx_ovl_fh + 20);

#pragma omp simd aligned(pa_x,             \
                             ts_x_xxxxx,   \
                             ts_x_xxxxy,   \
                             ts_x_xxxxz,   \
                             ts_x_xxxyy,   \
                             ts_x_xxxyz,   \
                             ts_x_xxxzz,   \
                             ts_x_xxyyy,   \
                             ts_x_xxyyz,   \
                             ts_x_xxyzz,   \
                             ts_x_xxzzz,   \
                             ts_x_xyyyy,   \
                             ts_x_xyyyz,   \
                             ts_x_xyyzz,   \
                             ts_x_xyzzz,   \
                             ts_x_xzzzz,   \
                             ts_x_yyyyy,   \
                             ts_x_yyyyz,   \
                             ts_x_yyyzz,   \
                             ts_x_yyzzz,   \
                             ts_x_yzzzz,   \
                             ts_x_zzzzz,   \
                             ts_xx_xxxx,   \
                             ts_xx_xxxxx,  \
                             ts_xx_xxxxy,  \
                             ts_xx_xxxxz,  \
                             ts_xx_xxxy,   \
                             ts_xx_xxxyy,  \
                             ts_xx_xxxyz,  \
                             ts_xx_xxxz,   \
                             ts_xx_xxxzz,  \
                             ts_xx_xxyy,   \
                             ts_xx_xxyyy,  \
                             ts_xx_xxyyz,  \
                             ts_xx_xxyz,   \
                             ts_xx_xxyzz,  \
                             ts_xx_xxzz,   \
                             ts_xx_xxzzz,  \
                             ts_xx_xyyy,   \
                             ts_xx_xyyyy,  \
                             ts_xx_xyyyz,  \
                             ts_xx_xyyz,   \
                             ts_xx_xyyzz,  \
                             ts_xx_xyzz,   \
                             ts_xx_xyzzz,  \
                             ts_xx_xzzz,   \
                             ts_xx_xzzzz,  \
                             ts_xx_yyyy,   \
                             ts_xx_yyyyy,  \
                             ts_xx_yyyyz,  \
                             ts_xx_yyyz,   \
                             ts_xx_yyyzz,  \
                             ts_xx_yyzz,   \
                             ts_xx_yyzzz,  \
                             ts_xx_yzzz,   \
                             ts_xx_yzzzz,  \
                             ts_xx_zzzz,   \
                             ts_xx_zzzzz,  \
                             ts_xxx_xxxxx, \
                             ts_xxx_xxxxy, \
                             ts_xxx_xxxxz, \
                             ts_xxx_xxxyy, \
                             ts_xxx_xxxyz, \
                             ts_xxx_xxxzz, \
                             ts_xxx_xxyyy, \
                             ts_xxx_xxyyz, \
                             ts_xxx_xxyzz, \
                             ts_xxx_xxzzz, \
                             ts_xxx_xyyyy, \
                             ts_xxx_xyyyz, \
                             ts_xxx_xyyzz, \
                             ts_xxx_xyzzz, \
                             ts_xxx_xzzzz, \
                             ts_xxx_yyyyy, \
                             ts_xxx_yyyyz, \
                             ts_xxx_yyyzz, \
                             ts_xxx_yyzzz, \
                             ts_xxx_yzzzz, \
                             ts_xxx_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxx_xxxxx[i] = 2.0 * ts_x_xxxxx[i] * fe_0 + 5.0 * ts_xx_xxxx[i] * fe_0 + ts_xx_xxxxx[i] * pa_x[i];

        ts_xxx_xxxxy[i] = 2.0 * ts_x_xxxxy[i] * fe_0 + 4.0 * ts_xx_xxxy[i] * fe_0 + ts_xx_xxxxy[i] * pa_x[i];

        ts_xxx_xxxxz[i] = 2.0 * ts_x_xxxxz[i] * fe_0 + 4.0 * ts_xx_xxxz[i] * fe_0 + ts_xx_xxxxz[i] * pa_x[i];

        ts_xxx_xxxyy[i] = 2.0 * ts_x_xxxyy[i] * fe_0 + 3.0 * ts_xx_xxyy[i] * fe_0 + ts_xx_xxxyy[i] * pa_x[i];

        ts_xxx_xxxyz[i] = 2.0 * ts_x_xxxyz[i] * fe_0 + 3.0 * ts_xx_xxyz[i] * fe_0 + ts_xx_xxxyz[i] * pa_x[i];

        ts_xxx_xxxzz[i] = 2.0 * ts_x_xxxzz[i] * fe_0 + 3.0 * ts_xx_xxzz[i] * fe_0 + ts_xx_xxxzz[i] * pa_x[i];

        ts_xxx_xxyyy[i] = 2.0 * ts_x_xxyyy[i] * fe_0 + 2.0 * ts_xx_xyyy[i] * fe_0 + ts_xx_xxyyy[i] * pa_x[i];

        ts_xxx_xxyyz[i] = 2.0 * ts_x_xxyyz[i] * fe_0 + 2.0 * ts_xx_xyyz[i] * fe_0 + ts_xx_xxyyz[i] * pa_x[i];

        ts_xxx_xxyzz[i] = 2.0 * ts_x_xxyzz[i] * fe_0 + 2.0 * ts_xx_xyzz[i] * fe_0 + ts_xx_xxyzz[i] * pa_x[i];

        ts_xxx_xxzzz[i] = 2.0 * ts_x_xxzzz[i] * fe_0 + 2.0 * ts_xx_xzzz[i] * fe_0 + ts_xx_xxzzz[i] * pa_x[i];

        ts_xxx_xyyyy[i] = 2.0 * ts_x_xyyyy[i] * fe_0 + ts_xx_yyyy[i] * fe_0 + ts_xx_xyyyy[i] * pa_x[i];

        ts_xxx_xyyyz[i] = 2.0 * ts_x_xyyyz[i] * fe_0 + ts_xx_yyyz[i] * fe_0 + ts_xx_xyyyz[i] * pa_x[i];

        ts_xxx_xyyzz[i] = 2.0 * ts_x_xyyzz[i] * fe_0 + ts_xx_yyzz[i] * fe_0 + ts_xx_xyyzz[i] * pa_x[i];

        ts_xxx_xyzzz[i] = 2.0 * ts_x_xyzzz[i] * fe_0 + ts_xx_yzzz[i] * fe_0 + ts_xx_xyzzz[i] * pa_x[i];

        ts_xxx_xzzzz[i] = 2.0 * ts_x_xzzzz[i] * fe_0 + ts_xx_zzzz[i] * fe_0 + ts_xx_xzzzz[i] * pa_x[i];

        ts_xxx_yyyyy[i] = 2.0 * ts_x_yyyyy[i] * fe_0 + ts_xx_yyyyy[i] * pa_x[i];

        ts_xxx_yyyyz[i] = 2.0 * ts_x_yyyyz[i] * fe_0 + ts_xx_yyyyz[i] * pa_x[i];

        ts_xxx_yyyzz[i] = 2.0 * ts_x_yyyzz[i] * fe_0 + ts_xx_yyyzz[i] * pa_x[i];

        ts_xxx_yyzzz[i] = 2.0 * ts_x_yyzzz[i] * fe_0 + ts_xx_yyzzz[i] * pa_x[i];

        ts_xxx_yzzzz[i] = 2.0 * ts_x_yzzzz[i] * fe_0 + ts_xx_yzzzz[i] * pa_x[i];

        ts_xxx_zzzzz[i] = 2.0 * ts_x_zzzzz[i] * fe_0 + ts_xx_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : FH

    auto ts_xxy_xxxxx = pbuffer.data(idx_ovl_fh + 21);

    auto ts_xxy_xxxxy = pbuffer.data(idx_ovl_fh + 22);

    auto ts_xxy_xxxxz = pbuffer.data(idx_ovl_fh + 23);

    auto ts_xxy_xxxyy = pbuffer.data(idx_ovl_fh + 24);

    auto ts_xxy_xxxyz = pbuffer.data(idx_ovl_fh + 25);

    auto ts_xxy_xxxzz = pbuffer.data(idx_ovl_fh + 26);

    auto ts_xxy_xxyyy = pbuffer.data(idx_ovl_fh + 27);

    auto ts_xxy_xxyyz = pbuffer.data(idx_ovl_fh + 28);

    auto ts_xxy_xxyzz = pbuffer.data(idx_ovl_fh + 29);

    auto ts_xxy_xxzzz = pbuffer.data(idx_ovl_fh + 30);

    auto ts_xxy_xyyyy = pbuffer.data(idx_ovl_fh + 31);

    auto ts_xxy_xyyyz = pbuffer.data(idx_ovl_fh + 32);

    auto ts_xxy_xyyzz = pbuffer.data(idx_ovl_fh + 33);

    auto ts_xxy_xyzzz = pbuffer.data(idx_ovl_fh + 34);

    auto ts_xxy_xzzzz = pbuffer.data(idx_ovl_fh + 35);

    auto ts_xxy_yyyyy = pbuffer.data(idx_ovl_fh + 36);

    auto ts_xxy_yyyyz = pbuffer.data(idx_ovl_fh + 37);

    auto ts_xxy_yyyzz = pbuffer.data(idx_ovl_fh + 38);

    auto ts_xxy_yyzzz = pbuffer.data(idx_ovl_fh + 39);

    auto ts_xxy_yzzzz = pbuffer.data(idx_ovl_fh + 40);

    auto ts_xxy_zzzzz = pbuffer.data(idx_ovl_fh + 41);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             ts_xx_xxxx,   \
                             ts_xx_xxxxx,  \
                             ts_xx_xxxxy,  \
                             ts_xx_xxxxz,  \
                             ts_xx_xxxy,   \
                             ts_xx_xxxyy,  \
                             ts_xx_xxxyz,  \
                             ts_xx_xxxz,   \
                             ts_xx_xxxzz,  \
                             ts_xx_xxyy,   \
                             ts_xx_xxyyy,  \
                             ts_xx_xxyyz,  \
                             ts_xx_xxyz,   \
                             ts_xx_xxyzz,  \
                             ts_xx_xxzz,   \
                             ts_xx_xxzzz,  \
                             ts_xx_xyyy,   \
                             ts_xx_xyyyy,  \
                             ts_xx_xyyyz,  \
                             ts_xx_xyyz,   \
                             ts_xx_xyyzz,  \
                             ts_xx_xyzz,   \
                             ts_xx_xyzzz,  \
                             ts_xx_xzzz,   \
                             ts_xx_xzzzz,  \
                             ts_xx_zzzzz,  \
                             ts_xxy_xxxxx, \
                             ts_xxy_xxxxy, \
                             ts_xxy_xxxxz, \
                             ts_xxy_xxxyy, \
                             ts_xxy_xxxyz, \
                             ts_xxy_xxxzz, \
                             ts_xxy_xxyyy, \
                             ts_xxy_xxyyz, \
                             ts_xxy_xxyzz, \
                             ts_xxy_xxzzz, \
                             ts_xxy_xyyyy, \
                             ts_xxy_xyyyz, \
                             ts_xxy_xyyzz, \
                             ts_xxy_xyzzz, \
                             ts_xxy_xzzzz, \
                             ts_xxy_yyyyy, \
                             ts_xxy_yyyyz, \
                             ts_xxy_yyyzz, \
                             ts_xxy_yyzzz, \
                             ts_xxy_yzzzz, \
                             ts_xxy_zzzzz, \
                             ts_xy_yyyyy,  \
                             ts_xy_yyyyz,  \
                             ts_xy_yyyzz,  \
                             ts_xy_yyzzz,  \
                             ts_xy_yzzzz,  \
                             ts_y_yyyyy,   \
                             ts_y_yyyyz,   \
                             ts_y_yyyzz,   \
                             ts_y_yyzzz,   \
                             ts_y_yzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxy_xxxxx[i] = ts_xx_xxxxx[i] * pa_y[i];

        ts_xxy_xxxxy[i] = ts_xx_xxxx[i] * fe_0 + ts_xx_xxxxy[i] * pa_y[i];

        ts_xxy_xxxxz[i] = ts_xx_xxxxz[i] * pa_y[i];

        ts_xxy_xxxyy[i] = 2.0 * ts_xx_xxxy[i] * fe_0 + ts_xx_xxxyy[i] * pa_y[i];

        ts_xxy_xxxyz[i] = ts_xx_xxxz[i] * fe_0 + ts_xx_xxxyz[i] * pa_y[i];

        ts_xxy_xxxzz[i] = ts_xx_xxxzz[i] * pa_y[i];

        ts_xxy_xxyyy[i] = 3.0 * ts_xx_xxyy[i] * fe_0 + ts_xx_xxyyy[i] * pa_y[i];

        ts_xxy_xxyyz[i] = 2.0 * ts_xx_xxyz[i] * fe_0 + ts_xx_xxyyz[i] * pa_y[i];

        ts_xxy_xxyzz[i] = ts_xx_xxzz[i] * fe_0 + ts_xx_xxyzz[i] * pa_y[i];

        ts_xxy_xxzzz[i] = ts_xx_xxzzz[i] * pa_y[i];

        ts_xxy_xyyyy[i] = 4.0 * ts_xx_xyyy[i] * fe_0 + ts_xx_xyyyy[i] * pa_y[i];

        ts_xxy_xyyyz[i] = 3.0 * ts_xx_xyyz[i] * fe_0 + ts_xx_xyyyz[i] * pa_y[i];

        ts_xxy_xyyzz[i] = 2.0 * ts_xx_xyzz[i] * fe_0 + ts_xx_xyyzz[i] * pa_y[i];

        ts_xxy_xyzzz[i] = ts_xx_xzzz[i] * fe_0 + ts_xx_xyzzz[i] * pa_y[i];

        ts_xxy_xzzzz[i] = ts_xx_xzzzz[i] * pa_y[i];

        ts_xxy_yyyyy[i] = ts_y_yyyyy[i] * fe_0 + ts_xy_yyyyy[i] * pa_x[i];

        ts_xxy_yyyyz[i] = ts_y_yyyyz[i] * fe_0 + ts_xy_yyyyz[i] * pa_x[i];

        ts_xxy_yyyzz[i] = ts_y_yyyzz[i] * fe_0 + ts_xy_yyyzz[i] * pa_x[i];

        ts_xxy_yyzzz[i] = ts_y_yyzzz[i] * fe_0 + ts_xy_yyzzz[i] * pa_x[i];

        ts_xxy_yzzzz[i] = ts_y_yzzzz[i] * fe_0 + ts_xy_yzzzz[i] * pa_x[i];

        ts_xxy_zzzzz[i] = ts_xx_zzzzz[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : FH

    auto ts_xxz_xxxxx = pbuffer.data(idx_ovl_fh + 42);

    auto ts_xxz_xxxxy = pbuffer.data(idx_ovl_fh + 43);

    auto ts_xxz_xxxxz = pbuffer.data(idx_ovl_fh + 44);

    auto ts_xxz_xxxyy = pbuffer.data(idx_ovl_fh + 45);

    auto ts_xxz_xxxyz = pbuffer.data(idx_ovl_fh + 46);

    auto ts_xxz_xxxzz = pbuffer.data(idx_ovl_fh + 47);

    auto ts_xxz_xxyyy = pbuffer.data(idx_ovl_fh + 48);

    auto ts_xxz_xxyyz = pbuffer.data(idx_ovl_fh + 49);

    auto ts_xxz_xxyzz = pbuffer.data(idx_ovl_fh + 50);

    auto ts_xxz_xxzzz = pbuffer.data(idx_ovl_fh + 51);

    auto ts_xxz_xyyyy = pbuffer.data(idx_ovl_fh + 52);

    auto ts_xxz_xyyyz = pbuffer.data(idx_ovl_fh + 53);

    auto ts_xxz_xyyzz = pbuffer.data(idx_ovl_fh + 54);

    auto ts_xxz_xyzzz = pbuffer.data(idx_ovl_fh + 55);

    auto ts_xxz_xzzzz = pbuffer.data(idx_ovl_fh + 56);

    auto ts_xxz_yyyyy = pbuffer.data(idx_ovl_fh + 57);

    auto ts_xxz_yyyyz = pbuffer.data(idx_ovl_fh + 58);

    auto ts_xxz_yyyzz = pbuffer.data(idx_ovl_fh + 59);

    auto ts_xxz_yyzzz = pbuffer.data(idx_ovl_fh + 60);

    auto ts_xxz_yzzzz = pbuffer.data(idx_ovl_fh + 61);

    auto ts_xxz_zzzzz = pbuffer.data(idx_ovl_fh + 62);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             ts_xx_xxxx,   \
                             ts_xx_xxxxx,  \
                             ts_xx_xxxxy,  \
                             ts_xx_xxxxz,  \
                             ts_xx_xxxy,   \
                             ts_xx_xxxyy,  \
                             ts_xx_xxxyz,  \
                             ts_xx_xxxz,   \
                             ts_xx_xxxzz,  \
                             ts_xx_xxyy,   \
                             ts_xx_xxyyy,  \
                             ts_xx_xxyyz,  \
                             ts_xx_xxyz,   \
                             ts_xx_xxyzz,  \
                             ts_xx_xxzz,   \
                             ts_xx_xxzzz,  \
                             ts_xx_xyyy,   \
                             ts_xx_xyyyy,  \
                             ts_xx_xyyyz,  \
                             ts_xx_xyyz,   \
                             ts_xx_xyyzz,  \
                             ts_xx_xyzz,   \
                             ts_xx_xyzzz,  \
                             ts_xx_xzzz,   \
                             ts_xx_xzzzz,  \
                             ts_xx_yyyyy,  \
                             ts_xxz_xxxxx, \
                             ts_xxz_xxxxy, \
                             ts_xxz_xxxxz, \
                             ts_xxz_xxxyy, \
                             ts_xxz_xxxyz, \
                             ts_xxz_xxxzz, \
                             ts_xxz_xxyyy, \
                             ts_xxz_xxyyz, \
                             ts_xxz_xxyzz, \
                             ts_xxz_xxzzz, \
                             ts_xxz_xyyyy, \
                             ts_xxz_xyyyz, \
                             ts_xxz_xyyzz, \
                             ts_xxz_xyzzz, \
                             ts_xxz_xzzzz, \
                             ts_xxz_yyyyy, \
                             ts_xxz_yyyyz, \
                             ts_xxz_yyyzz, \
                             ts_xxz_yyzzz, \
                             ts_xxz_yzzzz, \
                             ts_xxz_zzzzz, \
                             ts_xz_yyyyz,  \
                             ts_xz_yyyzz,  \
                             ts_xz_yyzzz,  \
                             ts_xz_yzzzz,  \
                             ts_xz_zzzzz,  \
                             ts_z_yyyyz,   \
                             ts_z_yyyzz,   \
                             ts_z_yyzzz,   \
                             ts_z_yzzzz,   \
                             ts_z_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxz_xxxxx[i] = ts_xx_xxxxx[i] * pa_z[i];

        ts_xxz_xxxxy[i] = ts_xx_xxxxy[i] * pa_z[i];

        ts_xxz_xxxxz[i] = ts_xx_xxxx[i] * fe_0 + ts_xx_xxxxz[i] * pa_z[i];

        ts_xxz_xxxyy[i] = ts_xx_xxxyy[i] * pa_z[i];

        ts_xxz_xxxyz[i] = ts_xx_xxxy[i] * fe_0 + ts_xx_xxxyz[i] * pa_z[i];

        ts_xxz_xxxzz[i] = 2.0 * ts_xx_xxxz[i] * fe_0 + ts_xx_xxxzz[i] * pa_z[i];

        ts_xxz_xxyyy[i] = ts_xx_xxyyy[i] * pa_z[i];

        ts_xxz_xxyyz[i] = ts_xx_xxyy[i] * fe_0 + ts_xx_xxyyz[i] * pa_z[i];

        ts_xxz_xxyzz[i] = 2.0 * ts_xx_xxyz[i] * fe_0 + ts_xx_xxyzz[i] * pa_z[i];

        ts_xxz_xxzzz[i] = 3.0 * ts_xx_xxzz[i] * fe_0 + ts_xx_xxzzz[i] * pa_z[i];

        ts_xxz_xyyyy[i] = ts_xx_xyyyy[i] * pa_z[i];

        ts_xxz_xyyyz[i] = ts_xx_xyyy[i] * fe_0 + ts_xx_xyyyz[i] * pa_z[i];

        ts_xxz_xyyzz[i] = 2.0 * ts_xx_xyyz[i] * fe_0 + ts_xx_xyyzz[i] * pa_z[i];

        ts_xxz_xyzzz[i] = 3.0 * ts_xx_xyzz[i] * fe_0 + ts_xx_xyzzz[i] * pa_z[i];

        ts_xxz_xzzzz[i] = 4.0 * ts_xx_xzzz[i] * fe_0 + ts_xx_xzzzz[i] * pa_z[i];

        ts_xxz_yyyyy[i] = ts_xx_yyyyy[i] * pa_z[i];

        ts_xxz_yyyyz[i] = ts_z_yyyyz[i] * fe_0 + ts_xz_yyyyz[i] * pa_x[i];

        ts_xxz_yyyzz[i] = ts_z_yyyzz[i] * fe_0 + ts_xz_yyyzz[i] * pa_x[i];

        ts_xxz_yyzzz[i] = ts_z_yyzzz[i] * fe_0 + ts_xz_yyzzz[i] * pa_x[i];

        ts_xxz_yzzzz[i] = ts_z_yzzzz[i] * fe_0 + ts_xz_yzzzz[i] * pa_x[i];

        ts_xxz_zzzzz[i] = ts_z_zzzzz[i] * fe_0 + ts_xz_zzzzz[i] * pa_x[i];
    }

    // Set up 63-84 components of targeted buffer : FH

    auto ts_xyy_xxxxx = pbuffer.data(idx_ovl_fh + 63);

    auto ts_xyy_xxxxy = pbuffer.data(idx_ovl_fh + 64);

    auto ts_xyy_xxxxz = pbuffer.data(idx_ovl_fh + 65);

    auto ts_xyy_xxxyy = pbuffer.data(idx_ovl_fh + 66);

    auto ts_xyy_xxxyz = pbuffer.data(idx_ovl_fh + 67);

    auto ts_xyy_xxxzz = pbuffer.data(idx_ovl_fh + 68);

    auto ts_xyy_xxyyy = pbuffer.data(idx_ovl_fh + 69);

    auto ts_xyy_xxyyz = pbuffer.data(idx_ovl_fh + 70);

    auto ts_xyy_xxyzz = pbuffer.data(idx_ovl_fh + 71);

    auto ts_xyy_xxzzz = pbuffer.data(idx_ovl_fh + 72);

    auto ts_xyy_xyyyy = pbuffer.data(idx_ovl_fh + 73);

    auto ts_xyy_xyyyz = pbuffer.data(idx_ovl_fh + 74);

    auto ts_xyy_xyyzz = pbuffer.data(idx_ovl_fh + 75);

    auto ts_xyy_xyzzz = pbuffer.data(idx_ovl_fh + 76);

    auto ts_xyy_xzzzz = pbuffer.data(idx_ovl_fh + 77);

    auto ts_xyy_yyyyy = pbuffer.data(idx_ovl_fh + 78);

    auto ts_xyy_yyyyz = pbuffer.data(idx_ovl_fh + 79);

    auto ts_xyy_yyyzz = pbuffer.data(idx_ovl_fh + 80);

    auto ts_xyy_yyzzz = pbuffer.data(idx_ovl_fh + 81);

    auto ts_xyy_yzzzz = pbuffer.data(idx_ovl_fh + 82);

    auto ts_xyy_zzzzz = pbuffer.data(idx_ovl_fh + 83);

#pragma omp simd aligned(pa_x,             \
                             ts_xyy_xxxxx, \
                             ts_xyy_xxxxy, \
                             ts_xyy_xxxxz, \
                             ts_xyy_xxxyy, \
                             ts_xyy_xxxyz, \
                             ts_xyy_xxxzz, \
                             ts_xyy_xxyyy, \
                             ts_xyy_xxyyz, \
                             ts_xyy_xxyzz, \
                             ts_xyy_xxzzz, \
                             ts_xyy_xyyyy, \
                             ts_xyy_xyyyz, \
                             ts_xyy_xyyzz, \
                             ts_xyy_xyzzz, \
                             ts_xyy_xzzzz, \
                             ts_xyy_yyyyy, \
                             ts_xyy_yyyyz, \
                             ts_xyy_yyyzz, \
                             ts_xyy_yyzzz, \
                             ts_xyy_yzzzz, \
                             ts_xyy_zzzzz, \
                             ts_yy_xxxx,   \
                             ts_yy_xxxxx,  \
                             ts_yy_xxxxy,  \
                             ts_yy_xxxxz,  \
                             ts_yy_xxxy,   \
                             ts_yy_xxxyy,  \
                             ts_yy_xxxyz,  \
                             ts_yy_xxxz,   \
                             ts_yy_xxxzz,  \
                             ts_yy_xxyy,   \
                             ts_yy_xxyyy,  \
                             ts_yy_xxyyz,  \
                             ts_yy_xxyz,   \
                             ts_yy_xxyzz,  \
                             ts_yy_xxzz,   \
                             ts_yy_xxzzz,  \
                             ts_yy_xyyy,   \
                             ts_yy_xyyyy,  \
                             ts_yy_xyyyz,  \
                             ts_yy_xyyz,   \
                             ts_yy_xyyzz,  \
                             ts_yy_xyzz,   \
                             ts_yy_xyzzz,  \
                             ts_yy_xzzz,   \
                             ts_yy_xzzzz,  \
                             ts_yy_yyyy,   \
                             ts_yy_yyyyy,  \
                             ts_yy_yyyyz,  \
                             ts_yy_yyyz,   \
                             ts_yy_yyyzz,  \
                             ts_yy_yyzz,   \
                             ts_yy_yyzzz,  \
                             ts_yy_yzzz,   \
                             ts_yy_yzzzz,  \
                             ts_yy_zzzz,   \
                             ts_yy_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyy_xxxxx[i] = 5.0 * ts_yy_xxxx[i] * fe_0 + ts_yy_xxxxx[i] * pa_x[i];

        ts_xyy_xxxxy[i] = 4.0 * ts_yy_xxxy[i] * fe_0 + ts_yy_xxxxy[i] * pa_x[i];

        ts_xyy_xxxxz[i] = 4.0 * ts_yy_xxxz[i] * fe_0 + ts_yy_xxxxz[i] * pa_x[i];

        ts_xyy_xxxyy[i] = 3.0 * ts_yy_xxyy[i] * fe_0 + ts_yy_xxxyy[i] * pa_x[i];

        ts_xyy_xxxyz[i] = 3.0 * ts_yy_xxyz[i] * fe_0 + ts_yy_xxxyz[i] * pa_x[i];

        ts_xyy_xxxzz[i] = 3.0 * ts_yy_xxzz[i] * fe_0 + ts_yy_xxxzz[i] * pa_x[i];

        ts_xyy_xxyyy[i] = 2.0 * ts_yy_xyyy[i] * fe_0 + ts_yy_xxyyy[i] * pa_x[i];

        ts_xyy_xxyyz[i] = 2.0 * ts_yy_xyyz[i] * fe_0 + ts_yy_xxyyz[i] * pa_x[i];

        ts_xyy_xxyzz[i] = 2.0 * ts_yy_xyzz[i] * fe_0 + ts_yy_xxyzz[i] * pa_x[i];

        ts_xyy_xxzzz[i] = 2.0 * ts_yy_xzzz[i] * fe_0 + ts_yy_xxzzz[i] * pa_x[i];

        ts_xyy_xyyyy[i] = ts_yy_yyyy[i] * fe_0 + ts_yy_xyyyy[i] * pa_x[i];

        ts_xyy_xyyyz[i] = ts_yy_yyyz[i] * fe_0 + ts_yy_xyyyz[i] * pa_x[i];

        ts_xyy_xyyzz[i] = ts_yy_yyzz[i] * fe_0 + ts_yy_xyyzz[i] * pa_x[i];

        ts_xyy_xyzzz[i] = ts_yy_yzzz[i] * fe_0 + ts_yy_xyzzz[i] * pa_x[i];

        ts_xyy_xzzzz[i] = ts_yy_zzzz[i] * fe_0 + ts_yy_xzzzz[i] * pa_x[i];

        ts_xyy_yyyyy[i] = ts_yy_yyyyy[i] * pa_x[i];

        ts_xyy_yyyyz[i] = ts_yy_yyyyz[i] * pa_x[i];

        ts_xyy_yyyzz[i] = ts_yy_yyyzz[i] * pa_x[i];

        ts_xyy_yyzzz[i] = ts_yy_yyzzz[i] * pa_x[i];

        ts_xyy_yzzzz[i] = ts_yy_yzzzz[i] * pa_x[i];

        ts_xyy_zzzzz[i] = ts_yy_zzzzz[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : FH

    auto ts_xyz_xxxxx = pbuffer.data(idx_ovl_fh + 84);

    auto ts_xyz_xxxxy = pbuffer.data(idx_ovl_fh + 85);

    auto ts_xyz_xxxxz = pbuffer.data(idx_ovl_fh + 86);

    auto ts_xyz_xxxyy = pbuffer.data(idx_ovl_fh + 87);

    auto ts_xyz_xxxyz = pbuffer.data(idx_ovl_fh + 88);

    auto ts_xyz_xxxzz = pbuffer.data(idx_ovl_fh + 89);

    auto ts_xyz_xxyyy = pbuffer.data(idx_ovl_fh + 90);

    auto ts_xyz_xxyyz = pbuffer.data(idx_ovl_fh + 91);

    auto ts_xyz_xxyzz = pbuffer.data(idx_ovl_fh + 92);

    auto ts_xyz_xxzzz = pbuffer.data(idx_ovl_fh + 93);

    auto ts_xyz_xyyyy = pbuffer.data(idx_ovl_fh + 94);

    auto ts_xyz_xyyyz = pbuffer.data(idx_ovl_fh + 95);

    auto ts_xyz_xyyzz = pbuffer.data(idx_ovl_fh + 96);

    auto ts_xyz_xyzzz = pbuffer.data(idx_ovl_fh + 97);

    auto ts_xyz_xzzzz = pbuffer.data(idx_ovl_fh + 98);

    auto ts_xyz_yyyyy = pbuffer.data(idx_ovl_fh + 99);

    auto ts_xyz_yyyyz = pbuffer.data(idx_ovl_fh + 100);

    auto ts_xyz_yyyzz = pbuffer.data(idx_ovl_fh + 101);

    auto ts_xyz_yyzzz = pbuffer.data(idx_ovl_fh + 102);

    auto ts_xyz_yzzzz = pbuffer.data(idx_ovl_fh + 103);

    auto ts_xyz_zzzzz = pbuffer.data(idx_ovl_fh + 104);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             ts_xy_xxxxy,  \
                             ts_xy_xxxyy,  \
                             ts_xy_xxyyy,  \
                             ts_xy_xyyyy,  \
                             ts_xyz_xxxxx, \
                             ts_xyz_xxxxy, \
                             ts_xyz_xxxxz, \
                             ts_xyz_xxxyy, \
                             ts_xyz_xxxyz, \
                             ts_xyz_xxxzz, \
                             ts_xyz_xxyyy, \
                             ts_xyz_xxyyz, \
                             ts_xyz_xxyzz, \
                             ts_xyz_xxzzz, \
                             ts_xyz_xyyyy, \
                             ts_xyz_xyyyz, \
                             ts_xyz_xyyzz, \
                             ts_xyz_xyzzz, \
                             ts_xyz_xzzzz, \
                             ts_xyz_yyyyy, \
                             ts_xyz_yyyyz, \
                             ts_xyz_yyyzz, \
                             ts_xyz_yyzzz, \
                             ts_xyz_yzzzz, \
                             ts_xyz_zzzzz, \
                             ts_xz_xxxxx,  \
                             ts_xz_xxxxz,  \
                             ts_xz_xxxzz,  \
                             ts_xz_xxzzz,  \
                             ts_xz_xzzzz,  \
                             ts_yz_xxxyz,  \
                             ts_yz_xxyyz,  \
                             ts_yz_xxyz,   \
                             ts_yz_xxyzz,  \
                             ts_yz_xyyyz,  \
                             ts_yz_xyyz,   \
                             ts_yz_xyyzz,  \
                             ts_yz_xyzz,   \
                             ts_yz_xyzzz,  \
                             ts_yz_yyyyy,  \
                             ts_yz_yyyyz,  \
                             ts_yz_yyyz,   \
                             ts_yz_yyyzz,  \
                             ts_yz_yyzz,   \
                             ts_yz_yyzzz,  \
                             ts_yz_yzzz,   \
                             ts_yz_yzzzz,  \
                             ts_yz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyz_xxxxx[i] = ts_xz_xxxxx[i] * pa_y[i];

        ts_xyz_xxxxy[i] = ts_xy_xxxxy[i] * pa_z[i];

        ts_xyz_xxxxz[i] = ts_xz_xxxxz[i] * pa_y[i];

        ts_xyz_xxxyy[i] = ts_xy_xxxyy[i] * pa_z[i];

        ts_xyz_xxxyz[i] = 3.0 * ts_yz_xxyz[i] * fe_0 + ts_yz_xxxyz[i] * pa_x[i];

        ts_xyz_xxxzz[i] = ts_xz_xxxzz[i] * pa_y[i];

        ts_xyz_xxyyy[i] = ts_xy_xxyyy[i] * pa_z[i];

        ts_xyz_xxyyz[i] = 2.0 * ts_yz_xyyz[i] * fe_0 + ts_yz_xxyyz[i] * pa_x[i];

        ts_xyz_xxyzz[i] = 2.0 * ts_yz_xyzz[i] * fe_0 + ts_yz_xxyzz[i] * pa_x[i];

        ts_xyz_xxzzz[i] = ts_xz_xxzzz[i] * pa_y[i];

        ts_xyz_xyyyy[i] = ts_xy_xyyyy[i] * pa_z[i];

        ts_xyz_xyyyz[i] = ts_yz_yyyz[i] * fe_0 + ts_yz_xyyyz[i] * pa_x[i];

        ts_xyz_xyyzz[i] = ts_yz_yyzz[i] * fe_0 + ts_yz_xyyzz[i] * pa_x[i];

        ts_xyz_xyzzz[i] = ts_yz_yzzz[i] * fe_0 + ts_yz_xyzzz[i] * pa_x[i];

        ts_xyz_xzzzz[i] = ts_xz_xzzzz[i] * pa_y[i];

        ts_xyz_yyyyy[i] = ts_yz_yyyyy[i] * pa_x[i];

        ts_xyz_yyyyz[i] = ts_yz_yyyyz[i] * pa_x[i];

        ts_xyz_yyyzz[i] = ts_yz_yyyzz[i] * pa_x[i];

        ts_xyz_yyzzz[i] = ts_yz_yyzzz[i] * pa_x[i];

        ts_xyz_yzzzz[i] = ts_yz_yzzzz[i] * pa_x[i];

        ts_xyz_zzzzz[i] = ts_yz_zzzzz[i] * pa_x[i];
    }

    // Set up 105-126 components of targeted buffer : FH

    auto ts_xzz_xxxxx = pbuffer.data(idx_ovl_fh + 105);

    auto ts_xzz_xxxxy = pbuffer.data(idx_ovl_fh + 106);

    auto ts_xzz_xxxxz = pbuffer.data(idx_ovl_fh + 107);

    auto ts_xzz_xxxyy = pbuffer.data(idx_ovl_fh + 108);

    auto ts_xzz_xxxyz = pbuffer.data(idx_ovl_fh + 109);

    auto ts_xzz_xxxzz = pbuffer.data(idx_ovl_fh + 110);

    auto ts_xzz_xxyyy = pbuffer.data(idx_ovl_fh + 111);

    auto ts_xzz_xxyyz = pbuffer.data(idx_ovl_fh + 112);

    auto ts_xzz_xxyzz = pbuffer.data(idx_ovl_fh + 113);

    auto ts_xzz_xxzzz = pbuffer.data(idx_ovl_fh + 114);

    auto ts_xzz_xyyyy = pbuffer.data(idx_ovl_fh + 115);

    auto ts_xzz_xyyyz = pbuffer.data(idx_ovl_fh + 116);

    auto ts_xzz_xyyzz = pbuffer.data(idx_ovl_fh + 117);

    auto ts_xzz_xyzzz = pbuffer.data(idx_ovl_fh + 118);

    auto ts_xzz_xzzzz = pbuffer.data(idx_ovl_fh + 119);

    auto ts_xzz_yyyyy = pbuffer.data(idx_ovl_fh + 120);

    auto ts_xzz_yyyyz = pbuffer.data(idx_ovl_fh + 121);

    auto ts_xzz_yyyzz = pbuffer.data(idx_ovl_fh + 122);

    auto ts_xzz_yyzzz = pbuffer.data(idx_ovl_fh + 123);

    auto ts_xzz_yzzzz = pbuffer.data(idx_ovl_fh + 124);

    auto ts_xzz_zzzzz = pbuffer.data(idx_ovl_fh + 125);

#pragma omp simd aligned(pa_x,             \
                             ts_xzz_xxxxx, \
                             ts_xzz_xxxxy, \
                             ts_xzz_xxxxz, \
                             ts_xzz_xxxyy, \
                             ts_xzz_xxxyz, \
                             ts_xzz_xxxzz, \
                             ts_xzz_xxyyy, \
                             ts_xzz_xxyyz, \
                             ts_xzz_xxyzz, \
                             ts_xzz_xxzzz, \
                             ts_xzz_xyyyy, \
                             ts_xzz_xyyyz, \
                             ts_xzz_xyyzz, \
                             ts_xzz_xyzzz, \
                             ts_xzz_xzzzz, \
                             ts_xzz_yyyyy, \
                             ts_xzz_yyyyz, \
                             ts_xzz_yyyzz, \
                             ts_xzz_yyzzz, \
                             ts_xzz_yzzzz, \
                             ts_xzz_zzzzz, \
                             ts_zz_xxxx,   \
                             ts_zz_xxxxx,  \
                             ts_zz_xxxxy,  \
                             ts_zz_xxxxz,  \
                             ts_zz_xxxy,   \
                             ts_zz_xxxyy,  \
                             ts_zz_xxxyz,  \
                             ts_zz_xxxz,   \
                             ts_zz_xxxzz,  \
                             ts_zz_xxyy,   \
                             ts_zz_xxyyy,  \
                             ts_zz_xxyyz,  \
                             ts_zz_xxyz,   \
                             ts_zz_xxyzz,  \
                             ts_zz_xxzz,   \
                             ts_zz_xxzzz,  \
                             ts_zz_xyyy,   \
                             ts_zz_xyyyy,  \
                             ts_zz_xyyyz,  \
                             ts_zz_xyyz,   \
                             ts_zz_xyyzz,  \
                             ts_zz_xyzz,   \
                             ts_zz_xyzzz,  \
                             ts_zz_xzzz,   \
                             ts_zz_xzzzz,  \
                             ts_zz_yyyy,   \
                             ts_zz_yyyyy,  \
                             ts_zz_yyyyz,  \
                             ts_zz_yyyz,   \
                             ts_zz_yyyzz,  \
                             ts_zz_yyzz,   \
                             ts_zz_yyzzz,  \
                             ts_zz_yzzz,   \
                             ts_zz_yzzzz,  \
                             ts_zz_zzzz,   \
                             ts_zz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzz_xxxxx[i] = 5.0 * ts_zz_xxxx[i] * fe_0 + ts_zz_xxxxx[i] * pa_x[i];

        ts_xzz_xxxxy[i] = 4.0 * ts_zz_xxxy[i] * fe_0 + ts_zz_xxxxy[i] * pa_x[i];

        ts_xzz_xxxxz[i] = 4.0 * ts_zz_xxxz[i] * fe_0 + ts_zz_xxxxz[i] * pa_x[i];

        ts_xzz_xxxyy[i] = 3.0 * ts_zz_xxyy[i] * fe_0 + ts_zz_xxxyy[i] * pa_x[i];

        ts_xzz_xxxyz[i] = 3.0 * ts_zz_xxyz[i] * fe_0 + ts_zz_xxxyz[i] * pa_x[i];

        ts_xzz_xxxzz[i] = 3.0 * ts_zz_xxzz[i] * fe_0 + ts_zz_xxxzz[i] * pa_x[i];

        ts_xzz_xxyyy[i] = 2.0 * ts_zz_xyyy[i] * fe_0 + ts_zz_xxyyy[i] * pa_x[i];

        ts_xzz_xxyyz[i] = 2.0 * ts_zz_xyyz[i] * fe_0 + ts_zz_xxyyz[i] * pa_x[i];

        ts_xzz_xxyzz[i] = 2.0 * ts_zz_xyzz[i] * fe_0 + ts_zz_xxyzz[i] * pa_x[i];

        ts_xzz_xxzzz[i] = 2.0 * ts_zz_xzzz[i] * fe_0 + ts_zz_xxzzz[i] * pa_x[i];

        ts_xzz_xyyyy[i] = ts_zz_yyyy[i] * fe_0 + ts_zz_xyyyy[i] * pa_x[i];

        ts_xzz_xyyyz[i] = ts_zz_yyyz[i] * fe_0 + ts_zz_xyyyz[i] * pa_x[i];

        ts_xzz_xyyzz[i] = ts_zz_yyzz[i] * fe_0 + ts_zz_xyyzz[i] * pa_x[i];

        ts_xzz_xyzzz[i] = ts_zz_yzzz[i] * fe_0 + ts_zz_xyzzz[i] * pa_x[i];

        ts_xzz_xzzzz[i] = ts_zz_zzzz[i] * fe_0 + ts_zz_xzzzz[i] * pa_x[i];

        ts_xzz_yyyyy[i] = ts_zz_yyyyy[i] * pa_x[i];

        ts_xzz_yyyyz[i] = ts_zz_yyyyz[i] * pa_x[i];

        ts_xzz_yyyzz[i] = ts_zz_yyyzz[i] * pa_x[i];

        ts_xzz_yyzzz[i] = ts_zz_yyzzz[i] * pa_x[i];

        ts_xzz_yzzzz[i] = ts_zz_yzzzz[i] * pa_x[i];

        ts_xzz_zzzzz[i] = ts_zz_zzzzz[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : FH

    auto ts_yyy_xxxxx = pbuffer.data(idx_ovl_fh + 126);

    auto ts_yyy_xxxxy = pbuffer.data(idx_ovl_fh + 127);

    auto ts_yyy_xxxxz = pbuffer.data(idx_ovl_fh + 128);

    auto ts_yyy_xxxyy = pbuffer.data(idx_ovl_fh + 129);

    auto ts_yyy_xxxyz = pbuffer.data(idx_ovl_fh + 130);

    auto ts_yyy_xxxzz = pbuffer.data(idx_ovl_fh + 131);

    auto ts_yyy_xxyyy = pbuffer.data(idx_ovl_fh + 132);

    auto ts_yyy_xxyyz = pbuffer.data(idx_ovl_fh + 133);

    auto ts_yyy_xxyzz = pbuffer.data(idx_ovl_fh + 134);

    auto ts_yyy_xxzzz = pbuffer.data(idx_ovl_fh + 135);

    auto ts_yyy_xyyyy = pbuffer.data(idx_ovl_fh + 136);

    auto ts_yyy_xyyyz = pbuffer.data(idx_ovl_fh + 137);

    auto ts_yyy_xyyzz = pbuffer.data(idx_ovl_fh + 138);

    auto ts_yyy_xyzzz = pbuffer.data(idx_ovl_fh + 139);

    auto ts_yyy_xzzzz = pbuffer.data(idx_ovl_fh + 140);

    auto ts_yyy_yyyyy = pbuffer.data(idx_ovl_fh + 141);

    auto ts_yyy_yyyyz = pbuffer.data(idx_ovl_fh + 142);

    auto ts_yyy_yyyzz = pbuffer.data(idx_ovl_fh + 143);

    auto ts_yyy_yyzzz = pbuffer.data(idx_ovl_fh + 144);

    auto ts_yyy_yzzzz = pbuffer.data(idx_ovl_fh + 145);

    auto ts_yyy_zzzzz = pbuffer.data(idx_ovl_fh + 146);

#pragma omp simd aligned(pa_y,             \
                             ts_y_xxxxx,   \
                             ts_y_xxxxy,   \
                             ts_y_xxxxz,   \
                             ts_y_xxxyy,   \
                             ts_y_xxxyz,   \
                             ts_y_xxxzz,   \
                             ts_y_xxyyy,   \
                             ts_y_xxyyz,   \
                             ts_y_xxyzz,   \
                             ts_y_xxzzz,   \
                             ts_y_xyyyy,   \
                             ts_y_xyyyz,   \
                             ts_y_xyyzz,   \
                             ts_y_xyzzz,   \
                             ts_y_xzzzz,   \
                             ts_y_yyyyy,   \
                             ts_y_yyyyz,   \
                             ts_y_yyyzz,   \
                             ts_y_yyzzz,   \
                             ts_y_yzzzz,   \
                             ts_y_zzzzz,   \
                             ts_yy_xxxx,   \
                             ts_yy_xxxxx,  \
                             ts_yy_xxxxy,  \
                             ts_yy_xxxxz,  \
                             ts_yy_xxxy,   \
                             ts_yy_xxxyy,  \
                             ts_yy_xxxyz,  \
                             ts_yy_xxxz,   \
                             ts_yy_xxxzz,  \
                             ts_yy_xxyy,   \
                             ts_yy_xxyyy,  \
                             ts_yy_xxyyz,  \
                             ts_yy_xxyz,   \
                             ts_yy_xxyzz,  \
                             ts_yy_xxzz,   \
                             ts_yy_xxzzz,  \
                             ts_yy_xyyy,   \
                             ts_yy_xyyyy,  \
                             ts_yy_xyyyz,  \
                             ts_yy_xyyz,   \
                             ts_yy_xyyzz,  \
                             ts_yy_xyzz,   \
                             ts_yy_xyzzz,  \
                             ts_yy_xzzz,   \
                             ts_yy_xzzzz,  \
                             ts_yy_yyyy,   \
                             ts_yy_yyyyy,  \
                             ts_yy_yyyyz,  \
                             ts_yy_yyyz,   \
                             ts_yy_yyyzz,  \
                             ts_yy_yyzz,   \
                             ts_yy_yyzzz,  \
                             ts_yy_yzzz,   \
                             ts_yy_yzzzz,  \
                             ts_yy_zzzz,   \
                             ts_yy_zzzzz,  \
                             ts_yyy_xxxxx, \
                             ts_yyy_xxxxy, \
                             ts_yyy_xxxxz, \
                             ts_yyy_xxxyy, \
                             ts_yyy_xxxyz, \
                             ts_yyy_xxxzz, \
                             ts_yyy_xxyyy, \
                             ts_yyy_xxyyz, \
                             ts_yyy_xxyzz, \
                             ts_yyy_xxzzz, \
                             ts_yyy_xyyyy, \
                             ts_yyy_xyyyz, \
                             ts_yyy_xyyzz, \
                             ts_yyy_xyzzz, \
                             ts_yyy_xzzzz, \
                             ts_yyy_yyyyy, \
                             ts_yyy_yyyyz, \
                             ts_yyy_yyyzz, \
                             ts_yyy_yyzzz, \
                             ts_yyy_yzzzz, \
                             ts_yyy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyy_xxxxx[i] = 2.0 * ts_y_xxxxx[i] * fe_0 + ts_yy_xxxxx[i] * pa_y[i];

        ts_yyy_xxxxy[i] = 2.0 * ts_y_xxxxy[i] * fe_0 + ts_yy_xxxx[i] * fe_0 + ts_yy_xxxxy[i] * pa_y[i];

        ts_yyy_xxxxz[i] = 2.0 * ts_y_xxxxz[i] * fe_0 + ts_yy_xxxxz[i] * pa_y[i];

        ts_yyy_xxxyy[i] = 2.0 * ts_y_xxxyy[i] * fe_0 + 2.0 * ts_yy_xxxy[i] * fe_0 + ts_yy_xxxyy[i] * pa_y[i];

        ts_yyy_xxxyz[i] = 2.0 * ts_y_xxxyz[i] * fe_0 + ts_yy_xxxz[i] * fe_0 + ts_yy_xxxyz[i] * pa_y[i];

        ts_yyy_xxxzz[i] = 2.0 * ts_y_xxxzz[i] * fe_0 + ts_yy_xxxzz[i] * pa_y[i];

        ts_yyy_xxyyy[i] = 2.0 * ts_y_xxyyy[i] * fe_0 + 3.0 * ts_yy_xxyy[i] * fe_0 + ts_yy_xxyyy[i] * pa_y[i];

        ts_yyy_xxyyz[i] = 2.0 * ts_y_xxyyz[i] * fe_0 + 2.0 * ts_yy_xxyz[i] * fe_0 + ts_yy_xxyyz[i] * pa_y[i];

        ts_yyy_xxyzz[i] = 2.0 * ts_y_xxyzz[i] * fe_0 + ts_yy_xxzz[i] * fe_0 + ts_yy_xxyzz[i] * pa_y[i];

        ts_yyy_xxzzz[i] = 2.0 * ts_y_xxzzz[i] * fe_0 + ts_yy_xxzzz[i] * pa_y[i];

        ts_yyy_xyyyy[i] = 2.0 * ts_y_xyyyy[i] * fe_0 + 4.0 * ts_yy_xyyy[i] * fe_0 + ts_yy_xyyyy[i] * pa_y[i];

        ts_yyy_xyyyz[i] = 2.0 * ts_y_xyyyz[i] * fe_0 + 3.0 * ts_yy_xyyz[i] * fe_0 + ts_yy_xyyyz[i] * pa_y[i];

        ts_yyy_xyyzz[i] = 2.0 * ts_y_xyyzz[i] * fe_0 + 2.0 * ts_yy_xyzz[i] * fe_0 + ts_yy_xyyzz[i] * pa_y[i];

        ts_yyy_xyzzz[i] = 2.0 * ts_y_xyzzz[i] * fe_0 + ts_yy_xzzz[i] * fe_0 + ts_yy_xyzzz[i] * pa_y[i];

        ts_yyy_xzzzz[i] = 2.0 * ts_y_xzzzz[i] * fe_0 + ts_yy_xzzzz[i] * pa_y[i];

        ts_yyy_yyyyy[i] = 2.0 * ts_y_yyyyy[i] * fe_0 + 5.0 * ts_yy_yyyy[i] * fe_0 + ts_yy_yyyyy[i] * pa_y[i];

        ts_yyy_yyyyz[i] = 2.0 * ts_y_yyyyz[i] * fe_0 + 4.0 * ts_yy_yyyz[i] * fe_0 + ts_yy_yyyyz[i] * pa_y[i];

        ts_yyy_yyyzz[i] = 2.0 * ts_y_yyyzz[i] * fe_0 + 3.0 * ts_yy_yyzz[i] * fe_0 + ts_yy_yyyzz[i] * pa_y[i];

        ts_yyy_yyzzz[i] = 2.0 * ts_y_yyzzz[i] * fe_0 + 2.0 * ts_yy_yzzz[i] * fe_0 + ts_yy_yyzzz[i] * pa_y[i];

        ts_yyy_yzzzz[i] = 2.0 * ts_y_yzzzz[i] * fe_0 + ts_yy_zzzz[i] * fe_0 + ts_yy_yzzzz[i] * pa_y[i];

        ts_yyy_zzzzz[i] = 2.0 * ts_y_zzzzz[i] * fe_0 + ts_yy_zzzzz[i] * pa_y[i];
    }

    // Set up 147-168 components of targeted buffer : FH

    auto ts_yyz_xxxxx = pbuffer.data(idx_ovl_fh + 147);

    auto ts_yyz_xxxxy = pbuffer.data(idx_ovl_fh + 148);

    auto ts_yyz_xxxxz = pbuffer.data(idx_ovl_fh + 149);

    auto ts_yyz_xxxyy = pbuffer.data(idx_ovl_fh + 150);

    auto ts_yyz_xxxyz = pbuffer.data(idx_ovl_fh + 151);

    auto ts_yyz_xxxzz = pbuffer.data(idx_ovl_fh + 152);

    auto ts_yyz_xxyyy = pbuffer.data(idx_ovl_fh + 153);

    auto ts_yyz_xxyyz = pbuffer.data(idx_ovl_fh + 154);

    auto ts_yyz_xxyzz = pbuffer.data(idx_ovl_fh + 155);

    auto ts_yyz_xxzzz = pbuffer.data(idx_ovl_fh + 156);

    auto ts_yyz_xyyyy = pbuffer.data(idx_ovl_fh + 157);

    auto ts_yyz_xyyyz = pbuffer.data(idx_ovl_fh + 158);

    auto ts_yyz_xyyzz = pbuffer.data(idx_ovl_fh + 159);

    auto ts_yyz_xyzzz = pbuffer.data(idx_ovl_fh + 160);

    auto ts_yyz_xzzzz = pbuffer.data(idx_ovl_fh + 161);

    auto ts_yyz_yyyyy = pbuffer.data(idx_ovl_fh + 162);

    auto ts_yyz_yyyyz = pbuffer.data(idx_ovl_fh + 163);

    auto ts_yyz_yyyzz = pbuffer.data(idx_ovl_fh + 164);

    auto ts_yyz_yyzzz = pbuffer.data(idx_ovl_fh + 165);

    auto ts_yyz_yzzzz = pbuffer.data(idx_ovl_fh + 166);

    auto ts_yyz_zzzzz = pbuffer.data(idx_ovl_fh + 167);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             ts_yy_xxxxx,  \
                             ts_yy_xxxxy,  \
                             ts_yy_xxxy,   \
                             ts_yy_xxxyy,  \
                             ts_yy_xxxyz,  \
                             ts_yy_xxyy,   \
                             ts_yy_xxyyy,  \
                             ts_yy_xxyyz,  \
                             ts_yy_xxyz,   \
                             ts_yy_xxyzz,  \
                             ts_yy_xyyy,   \
                             ts_yy_xyyyy,  \
                             ts_yy_xyyyz,  \
                             ts_yy_xyyz,   \
                             ts_yy_xyyzz,  \
                             ts_yy_xyzz,   \
                             ts_yy_xyzzz,  \
                             ts_yy_yyyy,   \
                             ts_yy_yyyyy,  \
                             ts_yy_yyyyz,  \
                             ts_yy_yyyz,   \
                             ts_yy_yyyzz,  \
                             ts_yy_yyzz,   \
                             ts_yy_yyzzz,  \
                             ts_yy_yzzz,   \
                             ts_yy_yzzzz,  \
                             ts_yyz_xxxxx, \
                             ts_yyz_xxxxy, \
                             ts_yyz_xxxxz, \
                             ts_yyz_xxxyy, \
                             ts_yyz_xxxyz, \
                             ts_yyz_xxxzz, \
                             ts_yyz_xxyyy, \
                             ts_yyz_xxyyz, \
                             ts_yyz_xxyzz, \
                             ts_yyz_xxzzz, \
                             ts_yyz_xyyyy, \
                             ts_yyz_xyyyz, \
                             ts_yyz_xyyzz, \
                             ts_yyz_xyzzz, \
                             ts_yyz_xzzzz, \
                             ts_yyz_yyyyy, \
                             ts_yyz_yyyyz, \
                             ts_yyz_yyyzz, \
                             ts_yyz_yyzzz, \
                             ts_yyz_yzzzz, \
                             ts_yyz_zzzzz, \
                             ts_yz_xxxxz,  \
                             ts_yz_xxxzz,  \
                             ts_yz_xxzzz,  \
                             ts_yz_xzzzz,  \
                             ts_yz_zzzzz,  \
                             ts_z_xxxxz,   \
                             ts_z_xxxzz,   \
                             ts_z_xxzzz,   \
                             ts_z_xzzzz,   \
                             ts_z_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyz_xxxxx[i] = ts_yy_xxxxx[i] * pa_z[i];

        ts_yyz_xxxxy[i] = ts_yy_xxxxy[i] * pa_z[i];

        ts_yyz_xxxxz[i] = ts_z_xxxxz[i] * fe_0 + ts_yz_xxxxz[i] * pa_y[i];

        ts_yyz_xxxyy[i] = ts_yy_xxxyy[i] * pa_z[i];

        ts_yyz_xxxyz[i] = ts_yy_xxxy[i] * fe_0 + ts_yy_xxxyz[i] * pa_z[i];

        ts_yyz_xxxzz[i] = ts_z_xxxzz[i] * fe_0 + ts_yz_xxxzz[i] * pa_y[i];

        ts_yyz_xxyyy[i] = ts_yy_xxyyy[i] * pa_z[i];

        ts_yyz_xxyyz[i] = ts_yy_xxyy[i] * fe_0 + ts_yy_xxyyz[i] * pa_z[i];

        ts_yyz_xxyzz[i] = 2.0 * ts_yy_xxyz[i] * fe_0 + ts_yy_xxyzz[i] * pa_z[i];

        ts_yyz_xxzzz[i] = ts_z_xxzzz[i] * fe_0 + ts_yz_xxzzz[i] * pa_y[i];

        ts_yyz_xyyyy[i] = ts_yy_xyyyy[i] * pa_z[i];

        ts_yyz_xyyyz[i] = ts_yy_xyyy[i] * fe_0 + ts_yy_xyyyz[i] * pa_z[i];

        ts_yyz_xyyzz[i] = 2.0 * ts_yy_xyyz[i] * fe_0 + ts_yy_xyyzz[i] * pa_z[i];

        ts_yyz_xyzzz[i] = 3.0 * ts_yy_xyzz[i] * fe_0 + ts_yy_xyzzz[i] * pa_z[i];

        ts_yyz_xzzzz[i] = ts_z_xzzzz[i] * fe_0 + ts_yz_xzzzz[i] * pa_y[i];

        ts_yyz_yyyyy[i] = ts_yy_yyyyy[i] * pa_z[i];

        ts_yyz_yyyyz[i] = ts_yy_yyyy[i] * fe_0 + ts_yy_yyyyz[i] * pa_z[i];

        ts_yyz_yyyzz[i] = 2.0 * ts_yy_yyyz[i] * fe_0 + ts_yy_yyyzz[i] * pa_z[i];

        ts_yyz_yyzzz[i] = 3.0 * ts_yy_yyzz[i] * fe_0 + ts_yy_yyzzz[i] * pa_z[i];

        ts_yyz_yzzzz[i] = 4.0 * ts_yy_yzzz[i] * fe_0 + ts_yy_yzzzz[i] * pa_z[i];

        ts_yyz_zzzzz[i] = ts_z_zzzzz[i] * fe_0 + ts_yz_zzzzz[i] * pa_y[i];
    }

    // Set up 168-189 components of targeted buffer : FH

    auto ts_yzz_xxxxx = pbuffer.data(idx_ovl_fh + 168);

    auto ts_yzz_xxxxy = pbuffer.data(idx_ovl_fh + 169);

    auto ts_yzz_xxxxz = pbuffer.data(idx_ovl_fh + 170);

    auto ts_yzz_xxxyy = pbuffer.data(idx_ovl_fh + 171);

    auto ts_yzz_xxxyz = pbuffer.data(idx_ovl_fh + 172);

    auto ts_yzz_xxxzz = pbuffer.data(idx_ovl_fh + 173);

    auto ts_yzz_xxyyy = pbuffer.data(idx_ovl_fh + 174);

    auto ts_yzz_xxyyz = pbuffer.data(idx_ovl_fh + 175);

    auto ts_yzz_xxyzz = pbuffer.data(idx_ovl_fh + 176);

    auto ts_yzz_xxzzz = pbuffer.data(idx_ovl_fh + 177);

    auto ts_yzz_xyyyy = pbuffer.data(idx_ovl_fh + 178);

    auto ts_yzz_xyyyz = pbuffer.data(idx_ovl_fh + 179);

    auto ts_yzz_xyyzz = pbuffer.data(idx_ovl_fh + 180);

    auto ts_yzz_xyzzz = pbuffer.data(idx_ovl_fh + 181);

    auto ts_yzz_xzzzz = pbuffer.data(idx_ovl_fh + 182);

    auto ts_yzz_yyyyy = pbuffer.data(idx_ovl_fh + 183);

    auto ts_yzz_yyyyz = pbuffer.data(idx_ovl_fh + 184);

    auto ts_yzz_yyyzz = pbuffer.data(idx_ovl_fh + 185);

    auto ts_yzz_yyzzz = pbuffer.data(idx_ovl_fh + 186);

    auto ts_yzz_yzzzz = pbuffer.data(idx_ovl_fh + 187);

    auto ts_yzz_zzzzz = pbuffer.data(idx_ovl_fh + 188);

#pragma omp simd aligned(pa_y,             \
                             ts_yzz_xxxxx, \
                             ts_yzz_xxxxy, \
                             ts_yzz_xxxxz, \
                             ts_yzz_xxxyy, \
                             ts_yzz_xxxyz, \
                             ts_yzz_xxxzz, \
                             ts_yzz_xxyyy, \
                             ts_yzz_xxyyz, \
                             ts_yzz_xxyzz, \
                             ts_yzz_xxzzz, \
                             ts_yzz_xyyyy, \
                             ts_yzz_xyyyz, \
                             ts_yzz_xyyzz, \
                             ts_yzz_xyzzz, \
                             ts_yzz_xzzzz, \
                             ts_yzz_yyyyy, \
                             ts_yzz_yyyyz, \
                             ts_yzz_yyyzz, \
                             ts_yzz_yyzzz, \
                             ts_yzz_yzzzz, \
                             ts_yzz_zzzzz, \
                             ts_zz_xxxx,   \
                             ts_zz_xxxxx,  \
                             ts_zz_xxxxy,  \
                             ts_zz_xxxxz,  \
                             ts_zz_xxxy,   \
                             ts_zz_xxxyy,  \
                             ts_zz_xxxyz,  \
                             ts_zz_xxxz,   \
                             ts_zz_xxxzz,  \
                             ts_zz_xxyy,   \
                             ts_zz_xxyyy,  \
                             ts_zz_xxyyz,  \
                             ts_zz_xxyz,   \
                             ts_zz_xxyzz,  \
                             ts_zz_xxzz,   \
                             ts_zz_xxzzz,  \
                             ts_zz_xyyy,   \
                             ts_zz_xyyyy,  \
                             ts_zz_xyyyz,  \
                             ts_zz_xyyz,   \
                             ts_zz_xyyzz,  \
                             ts_zz_xyzz,   \
                             ts_zz_xyzzz,  \
                             ts_zz_xzzz,   \
                             ts_zz_xzzzz,  \
                             ts_zz_yyyy,   \
                             ts_zz_yyyyy,  \
                             ts_zz_yyyyz,  \
                             ts_zz_yyyz,   \
                             ts_zz_yyyzz,  \
                             ts_zz_yyzz,   \
                             ts_zz_yyzzz,  \
                             ts_zz_yzzz,   \
                             ts_zz_yzzzz,  \
                             ts_zz_zzzz,   \
                             ts_zz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzz_xxxxx[i] = ts_zz_xxxxx[i] * pa_y[i];

        ts_yzz_xxxxy[i] = ts_zz_xxxx[i] * fe_0 + ts_zz_xxxxy[i] * pa_y[i];

        ts_yzz_xxxxz[i] = ts_zz_xxxxz[i] * pa_y[i];

        ts_yzz_xxxyy[i] = 2.0 * ts_zz_xxxy[i] * fe_0 + ts_zz_xxxyy[i] * pa_y[i];

        ts_yzz_xxxyz[i] = ts_zz_xxxz[i] * fe_0 + ts_zz_xxxyz[i] * pa_y[i];

        ts_yzz_xxxzz[i] = ts_zz_xxxzz[i] * pa_y[i];

        ts_yzz_xxyyy[i] = 3.0 * ts_zz_xxyy[i] * fe_0 + ts_zz_xxyyy[i] * pa_y[i];

        ts_yzz_xxyyz[i] = 2.0 * ts_zz_xxyz[i] * fe_0 + ts_zz_xxyyz[i] * pa_y[i];

        ts_yzz_xxyzz[i] = ts_zz_xxzz[i] * fe_0 + ts_zz_xxyzz[i] * pa_y[i];

        ts_yzz_xxzzz[i] = ts_zz_xxzzz[i] * pa_y[i];

        ts_yzz_xyyyy[i] = 4.0 * ts_zz_xyyy[i] * fe_0 + ts_zz_xyyyy[i] * pa_y[i];

        ts_yzz_xyyyz[i] = 3.0 * ts_zz_xyyz[i] * fe_0 + ts_zz_xyyyz[i] * pa_y[i];

        ts_yzz_xyyzz[i] = 2.0 * ts_zz_xyzz[i] * fe_0 + ts_zz_xyyzz[i] * pa_y[i];

        ts_yzz_xyzzz[i] = ts_zz_xzzz[i] * fe_0 + ts_zz_xyzzz[i] * pa_y[i];

        ts_yzz_xzzzz[i] = ts_zz_xzzzz[i] * pa_y[i];

        ts_yzz_yyyyy[i] = 5.0 * ts_zz_yyyy[i] * fe_0 + ts_zz_yyyyy[i] * pa_y[i];

        ts_yzz_yyyyz[i] = 4.0 * ts_zz_yyyz[i] * fe_0 + ts_zz_yyyyz[i] * pa_y[i];

        ts_yzz_yyyzz[i] = 3.0 * ts_zz_yyzz[i] * fe_0 + ts_zz_yyyzz[i] * pa_y[i];

        ts_yzz_yyzzz[i] = 2.0 * ts_zz_yzzz[i] * fe_0 + ts_zz_yyzzz[i] * pa_y[i];

        ts_yzz_yzzzz[i] = ts_zz_zzzz[i] * fe_0 + ts_zz_yzzzz[i] * pa_y[i];

        ts_yzz_zzzzz[i] = ts_zz_zzzzz[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : FH

    auto ts_zzz_xxxxx = pbuffer.data(idx_ovl_fh + 189);

    auto ts_zzz_xxxxy = pbuffer.data(idx_ovl_fh + 190);

    auto ts_zzz_xxxxz = pbuffer.data(idx_ovl_fh + 191);

    auto ts_zzz_xxxyy = pbuffer.data(idx_ovl_fh + 192);

    auto ts_zzz_xxxyz = pbuffer.data(idx_ovl_fh + 193);

    auto ts_zzz_xxxzz = pbuffer.data(idx_ovl_fh + 194);

    auto ts_zzz_xxyyy = pbuffer.data(idx_ovl_fh + 195);

    auto ts_zzz_xxyyz = pbuffer.data(idx_ovl_fh + 196);

    auto ts_zzz_xxyzz = pbuffer.data(idx_ovl_fh + 197);

    auto ts_zzz_xxzzz = pbuffer.data(idx_ovl_fh + 198);

    auto ts_zzz_xyyyy = pbuffer.data(idx_ovl_fh + 199);

    auto ts_zzz_xyyyz = pbuffer.data(idx_ovl_fh + 200);

    auto ts_zzz_xyyzz = pbuffer.data(idx_ovl_fh + 201);

    auto ts_zzz_xyzzz = pbuffer.data(idx_ovl_fh + 202);

    auto ts_zzz_xzzzz = pbuffer.data(idx_ovl_fh + 203);

    auto ts_zzz_yyyyy = pbuffer.data(idx_ovl_fh + 204);

    auto ts_zzz_yyyyz = pbuffer.data(idx_ovl_fh + 205);

    auto ts_zzz_yyyzz = pbuffer.data(idx_ovl_fh + 206);

    auto ts_zzz_yyzzz = pbuffer.data(idx_ovl_fh + 207);

    auto ts_zzz_yzzzz = pbuffer.data(idx_ovl_fh + 208);

    auto ts_zzz_zzzzz = pbuffer.data(idx_ovl_fh + 209);

#pragma omp simd aligned(pa_z,             \
                             ts_z_xxxxx,   \
                             ts_z_xxxxy,   \
                             ts_z_xxxxz,   \
                             ts_z_xxxyy,   \
                             ts_z_xxxyz,   \
                             ts_z_xxxzz,   \
                             ts_z_xxyyy,   \
                             ts_z_xxyyz,   \
                             ts_z_xxyzz,   \
                             ts_z_xxzzz,   \
                             ts_z_xyyyy,   \
                             ts_z_xyyyz,   \
                             ts_z_xyyzz,   \
                             ts_z_xyzzz,   \
                             ts_z_xzzzz,   \
                             ts_z_yyyyy,   \
                             ts_z_yyyyz,   \
                             ts_z_yyyzz,   \
                             ts_z_yyzzz,   \
                             ts_z_yzzzz,   \
                             ts_z_zzzzz,   \
                             ts_zz_xxxx,   \
                             ts_zz_xxxxx,  \
                             ts_zz_xxxxy,  \
                             ts_zz_xxxxz,  \
                             ts_zz_xxxy,   \
                             ts_zz_xxxyy,  \
                             ts_zz_xxxyz,  \
                             ts_zz_xxxz,   \
                             ts_zz_xxxzz,  \
                             ts_zz_xxyy,   \
                             ts_zz_xxyyy,  \
                             ts_zz_xxyyz,  \
                             ts_zz_xxyz,   \
                             ts_zz_xxyzz,  \
                             ts_zz_xxzz,   \
                             ts_zz_xxzzz,  \
                             ts_zz_xyyy,   \
                             ts_zz_xyyyy,  \
                             ts_zz_xyyyz,  \
                             ts_zz_xyyz,   \
                             ts_zz_xyyzz,  \
                             ts_zz_xyzz,   \
                             ts_zz_xyzzz,  \
                             ts_zz_xzzz,   \
                             ts_zz_xzzzz,  \
                             ts_zz_yyyy,   \
                             ts_zz_yyyyy,  \
                             ts_zz_yyyyz,  \
                             ts_zz_yyyz,   \
                             ts_zz_yyyzz,  \
                             ts_zz_yyzz,   \
                             ts_zz_yyzzz,  \
                             ts_zz_yzzz,   \
                             ts_zz_yzzzz,  \
                             ts_zz_zzzz,   \
                             ts_zz_zzzzz,  \
                             ts_zzz_xxxxx, \
                             ts_zzz_xxxxy, \
                             ts_zzz_xxxxz, \
                             ts_zzz_xxxyy, \
                             ts_zzz_xxxyz, \
                             ts_zzz_xxxzz, \
                             ts_zzz_xxyyy, \
                             ts_zzz_xxyyz, \
                             ts_zzz_xxyzz, \
                             ts_zzz_xxzzz, \
                             ts_zzz_xyyyy, \
                             ts_zzz_xyyyz, \
                             ts_zzz_xyyzz, \
                             ts_zzz_xyzzz, \
                             ts_zzz_xzzzz, \
                             ts_zzz_yyyyy, \
                             ts_zzz_yyyyz, \
                             ts_zzz_yyyzz, \
                             ts_zzz_yyzzz, \
                             ts_zzz_yzzzz, \
                             ts_zzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzz_xxxxx[i] = 2.0 * ts_z_xxxxx[i] * fe_0 + ts_zz_xxxxx[i] * pa_z[i];

        ts_zzz_xxxxy[i] = 2.0 * ts_z_xxxxy[i] * fe_0 + ts_zz_xxxxy[i] * pa_z[i];

        ts_zzz_xxxxz[i] = 2.0 * ts_z_xxxxz[i] * fe_0 + ts_zz_xxxx[i] * fe_0 + ts_zz_xxxxz[i] * pa_z[i];

        ts_zzz_xxxyy[i] = 2.0 * ts_z_xxxyy[i] * fe_0 + ts_zz_xxxyy[i] * pa_z[i];

        ts_zzz_xxxyz[i] = 2.0 * ts_z_xxxyz[i] * fe_0 + ts_zz_xxxy[i] * fe_0 + ts_zz_xxxyz[i] * pa_z[i];

        ts_zzz_xxxzz[i] = 2.0 * ts_z_xxxzz[i] * fe_0 + 2.0 * ts_zz_xxxz[i] * fe_0 + ts_zz_xxxzz[i] * pa_z[i];

        ts_zzz_xxyyy[i] = 2.0 * ts_z_xxyyy[i] * fe_0 + ts_zz_xxyyy[i] * pa_z[i];

        ts_zzz_xxyyz[i] = 2.0 * ts_z_xxyyz[i] * fe_0 + ts_zz_xxyy[i] * fe_0 + ts_zz_xxyyz[i] * pa_z[i];

        ts_zzz_xxyzz[i] = 2.0 * ts_z_xxyzz[i] * fe_0 + 2.0 * ts_zz_xxyz[i] * fe_0 + ts_zz_xxyzz[i] * pa_z[i];

        ts_zzz_xxzzz[i] = 2.0 * ts_z_xxzzz[i] * fe_0 + 3.0 * ts_zz_xxzz[i] * fe_0 + ts_zz_xxzzz[i] * pa_z[i];

        ts_zzz_xyyyy[i] = 2.0 * ts_z_xyyyy[i] * fe_0 + ts_zz_xyyyy[i] * pa_z[i];

        ts_zzz_xyyyz[i] = 2.0 * ts_z_xyyyz[i] * fe_0 + ts_zz_xyyy[i] * fe_0 + ts_zz_xyyyz[i] * pa_z[i];

        ts_zzz_xyyzz[i] = 2.0 * ts_z_xyyzz[i] * fe_0 + 2.0 * ts_zz_xyyz[i] * fe_0 + ts_zz_xyyzz[i] * pa_z[i];

        ts_zzz_xyzzz[i] = 2.0 * ts_z_xyzzz[i] * fe_0 + 3.0 * ts_zz_xyzz[i] * fe_0 + ts_zz_xyzzz[i] * pa_z[i];

        ts_zzz_xzzzz[i] = 2.0 * ts_z_xzzzz[i] * fe_0 + 4.0 * ts_zz_xzzz[i] * fe_0 + ts_zz_xzzzz[i] * pa_z[i];

        ts_zzz_yyyyy[i] = 2.0 * ts_z_yyyyy[i] * fe_0 + ts_zz_yyyyy[i] * pa_z[i];

        ts_zzz_yyyyz[i] = 2.0 * ts_z_yyyyz[i] * fe_0 + ts_zz_yyyy[i] * fe_0 + ts_zz_yyyyz[i] * pa_z[i];

        ts_zzz_yyyzz[i] = 2.0 * ts_z_yyyzz[i] * fe_0 + 2.0 * ts_zz_yyyz[i] * fe_0 + ts_zz_yyyzz[i] * pa_z[i];

        ts_zzz_yyzzz[i] = 2.0 * ts_z_yyzzz[i] * fe_0 + 3.0 * ts_zz_yyzz[i] * fe_0 + ts_zz_yyzzz[i] * pa_z[i];

        ts_zzz_yzzzz[i] = 2.0 * ts_z_yzzzz[i] * fe_0 + 4.0 * ts_zz_yzzz[i] * fe_0 + ts_zz_yzzzz[i] * pa_z[i];

        ts_zzz_zzzzz[i] = 2.0 * ts_z_zzzzz[i] * fe_0 + 5.0 * ts_zz_zzzz[i] * fe_0 + ts_zz_zzzzz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
