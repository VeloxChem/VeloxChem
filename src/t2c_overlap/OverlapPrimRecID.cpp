#include "OverlapPrimRecID.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_id(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_id,
                     const size_t idx_ovl_gd,
                     const size_t idx_ovl_hp,
                     const size_t idx_ovl_hd,
                     const CSimdArray<double>& factors,
                     const size_t idx_rpa,
                     const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_ovl_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_ovl_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_ovl_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_ovl_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_ovl_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_ovl_gd + 5);

    auto ts_xxxy_xx = pbuffer.data(idx_ovl_gd + 6);

    auto ts_xxxy_xz = pbuffer.data(idx_ovl_gd + 8);

    auto ts_xxxy_yy = pbuffer.data(idx_ovl_gd + 9);

    auto ts_xxxy_yz = pbuffer.data(idx_ovl_gd + 10);

    auto ts_xxxz_xx = pbuffer.data(idx_ovl_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_ovl_gd + 13);

    auto ts_xxxz_xz = pbuffer.data(idx_ovl_gd + 14);

    auto ts_xxxz_yz = pbuffer.data(idx_ovl_gd + 16);

    auto ts_xxxz_zz = pbuffer.data(idx_ovl_gd + 17);

    auto ts_xxyy_xx = pbuffer.data(idx_ovl_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_ovl_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_ovl_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_ovl_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_ovl_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_ovl_gd + 23);

    auto ts_xxyz_xz = pbuffer.data(idx_ovl_gd + 26);

    auto ts_xxyz_yz = pbuffer.data(idx_ovl_gd + 28);

    auto ts_xxzz_xx = pbuffer.data(idx_ovl_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_ovl_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_ovl_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_ovl_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_ovl_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_ovl_gd + 35);

    auto ts_xyyy_xy = pbuffer.data(idx_ovl_gd + 37);

    auto ts_xyyy_yy = pbuffer.data(idx_ovl_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_ovl_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_ovl_gd + 41);

    auto ts_xyyz_yz = pbuffer.data(idx_ovl_gd + 46);

    auto ts_xyyz_zz = pbuffer.data(idx_ovl_gd + 47);

    auto ts_xyzz_yy = pbuffer.data(idx_ovl_gd + 51);

    auto ts_xyzz_yz = pbuffer.data(idx_ovl_gd + 52);

    auto ts_xzzz_xz = pbuffer.data(idx_ovl_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_ovl_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_ovl_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_ovl_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_ovl_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_ovl_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_ovl_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_ovl_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_ovl_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_ovl_gd + 65);

    auto ts_yyyz_xy = pbuffer.data(idx_ovl_gd + 67);

    auto ts_yyyz_xz = pbuffer.data(idx_ovl_gd + 68);

    auto ts_yyyz_yy = pbuffer.data(idx_ovl_gd + 69);

    auto ts_yyyz_yz = pbuffer.data(idx_ovl_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_ovl_gd + 71);

    auto ts_yyzz_xx = pbuffer.data(idx_ovl_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_ovl_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_ovl_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_ovl_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_ovl_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_ovl_gd + 77);

    auto ts_yzzz_xx = pbuffer.data(idx_ovl_gd + 78);

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

    // Set up components of auxiliary buffer : HP

    auto ts_xxxxx_x = pbuffer.data(idx_ovl_hp);

    auto ts_xxxxx_y = pbuffer.data(idx_ovl_hp + 1);

    auto ts_xxxxx_z = pbuffer.data(idx_ovl_hp + 2);

    auto ts_xxxyy_y = pbuffer.data(idx_ovl_hp + 10);

    auto ts_xxxzz_x = pbuffer.data(idx_ovl_hp + 15);

    auto ts_xxxzz_z = pbuffer.data(idx_ovl_hp + 17);

    auto ts_xxyyy_y = pbuffer.data(idx_ovl_hp + 19);

    auto ts_xxzzz_x = pbuffer.data(idx_ovl_hp + 27);

    auto ts_xxzzz_z = pbuffer.data(idx_ovl_hp + 29);

    auto ts_xyyyy_y = pbuffer.data(idx_ovl_hp + 31);

    auto ts_xzzzz_z = pbuffer.data(idx_ovl_hp + 44);

    auto ts_yyyyy_x = pbuffer.data(idx_ovl_hp + 45);

    auto ts_yyyyy_y = pbuffer.data(idx_ovl_hp + 46);

    auto ts_yyyyy_z = pbuffer.data(idx_ovl_hp + 47);

    auto ts_yyyyz_z = pbuffer.data(idx_ovl_hp + 50);

    auto ts_yyyzz_x = pbuffer.data(idx_ovl_hp + 51);

    auto ts_yyyzz_y = pbuffer.data(idx_ovl_hp + 52);

    auto ts_yyyzz_z = pbuffer.data(idx_ovl_hp + 53);

    auto ts_yyzzz_x = pbuffer.data(idx_ovl_hp + 54);

    auto ts_yyzzz_y = pbuffer.data(idx_ovl_hp + 55);

    auto ts_yyzzz_z = pbuffer.data(idx_ovl_hp + 56);

    auto ts_yzzzz_y = pbuffer.data(idx_ovl_hp + 58);

    auto ts_yzzzz_z = pbuffer.data(idx_ovl_hp + 59);

    auto ts_zzzzz_x = pbuffer.data(idx_ovl_hp + 60);

    auto ts_zzzzz_y = pbuffer.data(idx_ovl_hp + 61);

    auto ts_zzzzz_z = pbuffer.data(idx_ovl_hp + 62);

    // Set up components of auxiliary buffer : HD

    auto ts_xxxxx_xx = pbuffer.data(idx_ovl_hd);

    auto ts_xxxxx_xy = pbuffer.data(idx_ovl_hd + 1);

    auto ts_xxxxx_xz = pbuffer.data(idx_ovl_hd + 2);

    auto ts_xxxxx_yy = pbuffer.data(idx_ovl_hd + 3);

    auto ts_xxxxx_yz = pbuffer.data(idx_ovl_hd + 4);

    auto ts_xxxxx_zz = pbuffer.data(idx_ovl_hd + 5);

    auto ts_xxxxy_xx = pbuffer.data(idx_ovl_hd + 6);

    auto ts_xxxxy_xy = pbuffer.data(idx_ovl_hd + 7);

    auto ts_xxxxy_xz = pbuffer.data(idx_ovl_hd + 8);

    auto ts_xxxxy_yy = pbuffer.data(idx_ovl_hd + 9);

    auto ts_xxxxy_yz = pbuffer.data(idx_ovl_hd + 10);

    auto ts_xxxxz_xx = pbuffer.data(idx_ovl_hd + 12);

    auto ts_xxxxz_xy = pbuffer.data(idx_ovl_hd + 13);

    auto ts_xxxxz_xz = pbuffer.data(idx_ovl_hd + 14);

    auto ts_xxxxz_yz = pbuffer.data(idx_ovl_hd + 16);

    auto ts_xxxxz_zz = pbuffer.data(idx_ovl_hd + 17);

    auto ts_xxxyy_xx = pbuffer.data(idx_ovl_hd + 18);

    auto ts_xxxyy_xy = pbuffer.data(idx_ovl_hd + 19);

    auto ts_xxxyy_xz = pbuffer.data(idx_ovl_hd + 20);

    auto ts_xxxyy_yy = pbuffer.data(idx_ovl_hd + 21);

    auto ts_xxxyy_yz = pbuffer.data(idx_ovl_hd + 22);

    auto ts_xxxyy_zz = pbuffer.data(idx_ovl_hd + 23);

    auto ts_xxxyz_xz = pbuffer.data(idx_ovl_hd + 26);

    auto ts_xxxyz_yz = pbuffer.data(idx_ovl_hd + 28);

    auto ts_xxxzz_xx = pbuffer.data(idx_ovl_hd + 30);

    auto ts_xxxzz_xy = pbuffer.data(idx_ovl_hd + 31);

    auto ts_xxxzz_xz = pbuffer.data(idx_ovl_hd + 32);

    auto ts_xxxzz_yy = pbuffer.data(idx_ovl_hd + 33);

    auto ts_xxxzz_yz = pbuffer.data(idx_ovl_hd + 34);

    auto ts_xxxzz_zz = pbuffer.data(idx_ovl_hd + 35);

    auto ts_xxyyy_xx = pbuffer.data(idx_ovl_hd + 36);

    auto ts_xxyyy_xy = pbuffer.data(idx_ovl_hd + 37);

    auto ts_xxyyy_xz = pbuffer.data(idx_ovl_hd + 38);

    auto ts_xxyyy_yy = pbuffer.data(idx_ovl_hd + 39);

    auto ts_xxyyy_yz = pbuffer.data(idx_ovl_hd + 40);

    auto ts_xxyyy_zz = pbuffer.data(idx_ovl_hd + 41);

    auto ts_xxyyz_xy = pbuffer.data(idx_ovl_hd + 43);

    auto ts_xxyyz_xz = pbuffer.data(idx_ovl_hd + 44);

    auto ts_xxyyz_yz = pbuffer.data(idx_ovl_hd + 46);

    auto ts_xxyyz_zz = pbuffer.data(idx_ovl_hd + 47);

    auto ts_xxyzz_xx = pbuffer.data(idx_ovl_hd + 48);

    auto ts_xxyzz_xz = pbuffer.data(idx_ovl_hd + 50);

    auto ts_xxyzz_yy = pbuffer.data(idx_ovl_hd + 51);

    auto ts_xxyzz_yz = pbuffer.data(idx_ovl_hd + 52);

    auto ts_xxzzz_xx = pbuffer.data(idx_ovl_hd + 54);

    auto ts_xxzzz_xy = pbuffer.data(idx_ovl_hd + 55);

    auto ts_xxzzz_xz = pbuffer.data(idx_ovl_hd + 56);

    auto ts_xxzzz_yy = pbuffer.data(idx_ovl_hd + 57);

    auto ts_xxzzz_yz = pbuffer.data(idx_ovl_hd + 58);

    auto ts_xxzzz_zz = pbuffer.data(idx_ovl_hd + 59);

    auto ts_xyyyy_xx = pbuffer.data(idx_ovl_hd + 60);

    auto ts_xyyyy_xy = pbuffer.data(idx_ovl_hd + 61);

    auto ts_xyyyy_yy = pbuffer.data(idx_ovl_hd + 63);

    auto ts_xyyyy_yz = pbuffer.data(idx_ovl_hd + 64);

    auto ts_xyyyy_zz = pbuffer.data(idx_ovl_hd + 65);

    auto ts_xyyyz_yz = pbuffer.data(idx_ovl_hd + 70);

    auto ts_xyyyz_zz = pbuffer.data(idx_ovl_hd + 71);

    auto ts_xyyzz_yy = pbuffer.data(idx_ovl_hd + 75);

    auto ts_xyyzz_yz = pbuffer.data(idx_ovl_hd + 76);

    auto ts_xyyzz_zz = pbuffer.data(idx_ovl_hd + 77);

    auto ts_xyzzz_yy = pbuffer.data(idx_ovl_hd + 81);

    auto ts_xyzzz_yz = pbuffer.data(idx_ovl_hd + 82);

    auto ts_xzzzz_xx = pbuffer.data(idx_ovl_hd + 84);

    auto ts_xzzzz_xz = pbuffer.data(idx_ovl_hd + 86);

    auto ts_xzzzz_yy = pbuffer.data(idx_ovl_hd + 87);

    auto ts_xzzzz_yz = pbuffer.data(idx_ovl_hd + 88);

    auto ts_xzzzz_zz = pbuffer.data(idx_ovl_hd + 89);

    auto ts_yyyyy_xx = pbuffer.data(idx_ovl_hd + 90);

    auto ts_yyyyy_xy = pbuffer.data(idx_ovl_hd + 91);

    auto ts_yyyyy_xz = pbuffer.data(idx_ovl_hd + 92);

    auto ts_yyyyy_yy = pbuffer.data(idx_ovl_hd + 93);

    auto ts_yyyyy_yz = pbuffer.data(idx_ovl_hd + 94);

    auto ts_yyyyy_zz = pbuffer.data(idx_ovl_hd + 95);

    auto ts_yyyyz_xy = pbuffer.data(idx_ovl_hd + 97);

    auto ts_yyyyz_xz = pbuffer.data(idx_ovl_hd + 98);

    auto ts_yyyyz_yy = pbuffer.data(idx_ovl_hd + 99);

    auto ts_yyyyz_yz = pbuffer.data(idx_ovl_hd + 100);

    auto ts_yyyyz_zz = pbuffer.data(idx_ovl_hd + 101);

    auto ts_yyyzz_xx = pbuffer.data(idx_ovl_hd + 102);

    auto ts_yyyzz_xy = pbuffer.data(idx_ovl_hd + 103);

    auto ts_yyyzz_xz = pbuffer.data(idx_ovl_hd + 104);

    auto ts_yyyzz_yy = pbuffer.data(idx_ovl_hd + 105);

    auto ts_yyyzz_yz = pbuffer.data(idx_ovl_hd + 106);

    auto ts_yyyzz_zz = pbuffer.data(idx_ovl_hd + 107);

    auto ts_yyzzz_xx = pbuffer.data(idx_ovl_hd + 108);

    auto ts_yyzzz_xy = pbuffer.data(idx_ovl_hd + 109);

    auto ts_yyzzz_xz = pbuffer.data(idx_ovl_hd + 110);

    auto ts_yyzzz_yy = pbuffer.data(idx_ovl_hd + 111);

    auto ts_yyzzz_yz = pbuffer.data(idx_ovl_hd + 112);

    auto ts_yyzzz_zz = pbuffer.data(idx_ovl_hd + 113);

    auto ts_yzzzz_xx = pbuffer.data(idx_ovl_hd + 114);

    auto ts_yzzzz_xy = pbuffer.data(idx_ovl_hd + 115);

    auto ts_yzzzz_xz = pbuffer.data(idx_ovl_hd + 116);

    auto ts_yzzzz_yy = pbuffer.data(idx_ovl_hd + 117);

    auto ts_yzzzz_yz = pbuffer.data(idx_ovl_hd + 118);

    auto ts_yzzzz_zz = pbuffer.data(idx_ovl_hd + 119);

    auto ts_zzzzz_xx = pbuffer.data(idx_ovl_hd + 120);

    auto ts_zzzzz_xy = pbuffer.data(idx_ovl_hd + 121);

    auto ts_zzzzz_xz = pbuffer.data(idx_ovl_hd + 122);

    auto ts_zzzzz_yy = pbuffer.data(idx_ovl_hd + 123);

    auto ts_zzzzz_yz = pbuffer.data(idx_ovl_hd + 124);

    auto ts_zzzzz_zz = pbuffer.data(idx_ovl_hd + 125);

    // Set up 0-6 components of targeted buffer : ID

    auto ts_xxxxxx_xx = pbuffer.data(idx_ovl_id);

    auto ts_xxxxxx_xy = pbuffer.data(idx_ovl_id + 1);

    auto ts_xxxxxx_xz = pbuffer.data(idx_ovl_id + 2);

    auto ts_xxxxxx_yy = pbuffer.data(idx_ovl_id + 3);

    auto ts_xxxxxx_yz = pbuffer.data(idx_ovl_id + 4);

    auto ts_xxxxxx_zz = pbuffer.data(idx_ovl_id + 5);

    #pragma omp simd aligned(pa_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_zz, ts_xxxxx_x, ts_xxxxx_xx, ts_xxxxx_xy, ts_xxxxx_xz, ts_xxxxx_y, ts_xxxxx_yy, ts_xxxxx_yz, ts_xxxxx_z, ts_xxxxx_zz, ts_xxxxxx_xx, ts_xxxxxx_xy, ts_xxxxxx_xz, ts_xxxxxx_yy, ts_xxxxxx_yz, ts_xxxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxx_xx[i] = 5.0 * ts_xxxx_xx[i] * fe_0 + 2.0 * ts_xxxxx_x[i] * fe_0 + ts_xxxxx_xx[i] * pa_x[i];

        ts_xxxxxx_xy[i] = 5.0 * ts_xxxx_xy[i] * fe_0 + ts_xxxxx_y[i] * fe_0 + ts_xxxxx_xy[i] * pa_x[i];

        ts_xxxxxx_xz[i] = 5.0 * ts_xxxx_xz[i] * fe_0 + ts_xxxxx_z[i] * fe_0 + ts_xxxxx_xz[i] * pa_x[i];

        ts_xxxxxx_yy[i] = 5.0 * ts_xxxx_yy[i] * fe_0 + ts_xxxxx_yy[i] * pa_x[i];

        ts_xxxxxx_yz[i] = 5.0 * ts_xxxx_yz[i] * fe_0 + ts_xxxxx_yz[i] * pa_x[i];

        ts_xxxxxx_zz[i] = 5.0 * ts_xxxx_zz[i] * fe_0 + ts_xxxxx_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : ID

    auto ts_xxxxxy_xx = pbuffer.data(idx_ovl_id + 6);

    auto ts_xxxxxy_xy = pbuffer.data(idx_ovl_id + 7);

    auto ts_xxxxxy_xz = pbuffer.data(idx_ovl_id + 8);

    auto ts_xxxxxy_yy = pbuffer.data(idx_ovl_id + 9);

    auto ts_xxxxxy_yz = pbuffer.data(idx_ovl_id + 10);

    auto ts_xxxxxy_zz = pbuffer.data(idx_ovl_id + 11);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxxx_x, ts_xxxxx_xx, ts_xxxxx_xy, ts_xxxxx_xz, ts_xxxxx_zz, ts_xxxxxy_xx, ts_xxxxxy_xy, ts_xxxxxy_xz, ts_xxxxxy_yy, ts_xxxxxy_yz, ts_xxxxxy_zz, ts_xxxxy_yy, ts_xxxxy_yz, ts_xxxy_yy, ts_xxxy_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxy_xx[i] = ts_xxxxx_xx[i] * pa_y[i];

        ts_xxxxxy_xy[i] = ts_xxxxx_x[i] * fe_0 + ts_xxxxx_xy[i] * pa_y[i];

        ts_xxxxxy_xz[i] = ts_xxxxx_xz[i] * pa_y[i];

        ts_xxxxxy_yy[i] = 4.0 * ts_xxxy_yy[i] * fe_0 + ts_xxxxy_yy[i] * pa_x[i];

        ts_xxxxxy_yz[i] = 4.0 * ts_xxxy_yz[i] * fe_0 + ts_xxxxy_yz[i] * pa_x[i];

        ts_xxxxxy_zz[i] = ts_xxxxx_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : ID

    auto ts_xxxxxz_xx = pbuffer.data(idx_ovl_id + 12);

    auto ts_xxxxxz_xy = pbuffer.data(idx_ovl_id + 13);

    auto ts_xxxxxz_xz = pbuffer.data(idx_ovl_id + 14);

    auto ts_xxxxxz_yy = pbuffer.data(idx_ovl_id + 15);

    auto ts_xxxxxz_yz = pbuffer.data(idx_ovl_id + 16);

    auto ts_xxxxxz_zz = pbuffer.data(idx_ovl_id + 17);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxxx_x, ts_xxxxx_xx, ts_xxxxx_xy, ts_xxxxx_xz, ts_xxxxx_yy, ts_xxxxxz_xx, ts_xxxxxz_xy, ts_xxxxxz_xz, ts_xxxxxz_yy, ts_xxxxxz_yz, ts_xxxxxz_zz, ts_xxxxz_yz, ts_xxxxz_zz, ts_xxxz_yz, ts_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxz_xx[i] = ts_xxxxx_xx[i] * pa_z[i];

        ts_xxxxxz_xy[i] = ts_xxxxx_xy[i] * pa_z[i];

        ts_xxxxxz_xz[i] = ts_xxxxx_x[i] * fe_0 + ts_xxxxx_xz[i] * pa_z[i];

        ts_xxxxxz_yy[i] = ts_xxxxx_yy[i] * pa_z[i];

        ts_xxxxxz_yz[i] = 4.0 * ts_xxxz_yz[i] * fe_0 + ts_xxxxz_yz[i] * pa_x[i];

        ts_xxxxxz_zz[i] = 4.0 * ts_xxxz_zz[i] * fe_0 + ts_xxxxz_zz[i] * pa_x[i];
    }

    // Set up 18-24 components of targeted buffer : ID

    auto ts_xxxxyy_xx = pbuffer.data(idx_ovl_id + 18);

    auto ts_xxxxyy_xy = pbuffer.data(idx_ovl_id + 19);

    auto ts_xxxxyy_xz = pbuffer.data(idx_ovl_id + 20);

    auto ts_xxxxyy_yy = pbuffer.data(idx_ovl_id + 21);

    auto ts_xxxxyy_yz = pbuffer.data(idx_ovl_id + 22);

    auto ts_xxxxyy_zz = pbuffer.data(idx_ovl_id + 23);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxx_xx, ts_xxxx_xz, ts_xxxxy_xx, ts_xxxxy_xz, ts_xxxxyy_xx, ts_xxxxyy_xy, ts_xxxxyy_xz, ts_xxxxyy_yy, ts_xxxxyy_yz, ts_xxxxyy_zz, ts_xxxyy_xy, ts_xxxyy_y, ts_xxxyy_yy, ts_xxxyy_yz, ts_xxxyy_zz, ts_xxyy_xy, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyy_xx[i] = ts_xxxx_xx[i] * fe_0 + ts_xxxxy_xx[i] * pa_y[i];

        ts_xxxxyy_xy[i] = 3.0 * ts_xxyy_xy[i] * fe_0 + ts_xxxyy_y[i] * fe_0 + ts_xxxyy_xy[i] * pa_x[i];

        ts_xxxxyy_xz[i] = ts_xxxx_xz[i] * fe_0 + ts_xxxxy_xz[i] * pa_y[i];

        ts_xxxxyy_yy[i] = 3.0 * ts_xxyy_yy[i] * fe_0 + ts_xxxyy_yy[i] * pa_x[i];

        ts_xxxxyy_yz[i] = 3.0 * ts_xxyy_yz[i] * fe_0 + ts_xxxyy_yz[i] * pa_x[i];

        ts_xxxxyy_zz[i] = 3.0 * ts_xxyy_zz[i] * fe_0 + ts_xxxyy_zz[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : ID

    auto ts_xxxxyz_xx = pbuffer.data(idx_ovl_id + 24);

    auto ts_xxxxyz_xy = pbuffer.data(idx_ovl_id + 25);

    auto ts_xxxxyz_xz = pbuffer.data(idx_ovl_id + 26);

    auto ts_xxxxyz_yy = pbuffer.data(idx_ovl_id + 27);

    auto ts_xxxxyz_yz = pbuffer.data(idx_ovl_id + 28);

    auto ts_xxxxyz_zz = pbuffer.data(idx_ovl_id + 29);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxxy_xy, ts_xxxxy_yy, ts_xxxxyz_xx, ts_xxxxyz_xy, ts_xxxxyz_xz, ts_xxxxyz_yy, ts_xxxxyz_yz, ts_xxxxyz_zz, ts_xxxxz_xx, ts_xxxxz_xz, ts_xxxxz_zz, ts_xxxyz_yz, ts_xxyz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyz_xx[i] = ts_xxxxz_xx[i] * pa_y[i];

        ts_xxxxyz_xy[i] = ts_xxxxy_xy[i] * pa_z[i];

        ts_xxxxyz_xz[i] = ts_xxxxz_xz[i] * pa_y[i];

        ts_xxxxyz_yy[i] = ts_xxxxy_yy[i] * pa_z[i];

        ts_xxxxyz_yz[i] = 3.0 * ts_xxyz_yz[i] * fe_0 + ts_xxxyz_yz[i] * pa_x[i];

        ts_xxxxyz_zz[i] = ts_xxxxz_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : ID

    auto ts_xxxxzz_xx = pbuffer.data(idx_ovl_id + 30);

    auto ts_xxxxzz_xy = pbuffer.data(idx_ovl_id + 31);

    auto ts_xxxxzz_xz = pbuffer.data(idx_ovl_id + 32);

    auto ts_xxxxzz_yy = pbuffer.data(idx_ovl_id + 33);

    auto ts_xxxxzz_yz = pbuffer.data(idx_ovl_id + 34);

    auto ts_xxxxzz_zz = pbuffer.data(idx_ovl_id + 35);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxx_xx, ts_xxxx_xy, ts_xxxxz_xx, ts_xxxxz_xy, ts_xxxxzz_xx, ts_xxxxzz_xy, ts_xxxxzz_xz, ts_xxxxzz_yy, ts_xxxxzz_yz, ts_xxxxzz_zz, ts_xxxzz_xz, ts_xxxzz_yy, ts_xxxzz_yz, ts_xxxzz_z, ts_xxxzz_zz, ts_xxzz_xz, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxzz_xx[i] = ts_xxxx_xx[i] * fe_0 + ts_xxxxz_xx[i] * pa_z[i];

        ts_xxxxzz_xy[i] = ts_xxxx_xy[i] * fe_0 + ts_xxxxz_xy[i] * pa_z[i];

        ts_xxxxzz_xz[i] = 3.0 * ts_xxzz_xz[i] * fe_0 + ts_xxxzz_z[i] * fe_0 + ts_xxxzz_xz[i] * pa_x[i];

        ts_xxxxzz_yy[i] = 3.0 * ts_xxzz_yy[i] * fe_0 + ts_xxxzz_yy[i] * pa_x[i];

        ts_xxxxzz_yz[i] = 3.0 * ts_xxzz_yz[i] * fe_0 + ts_xxxzz_yz[i] * pa_x[i];

        ts_xxxxzz_zz[i] = 3.0 * ts_xxzz_zz[i] * fe_0 + ts_xxxzz_zz[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : ID

    auto ts_xxxyyy_xx = pbuffer.data(idx_ovl_id + 36);

    auto ts_xxxyyy_xy = pbuffer.data(idx_ovl_id + 37);

    auto ts_xxxyyy_xz = pbuffer.data(idx_ovl_id + 38);

    auto ts_xxxyyy_yy = pbuffer.data(idx_ovl_id + 39);

    auto ts_xxxyyy_yz = pbuffer.data(idx_ovl_id + 40);

    auto ts_xxxyyy_zz = pbuffer.data(idx_ovl_id + 41);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxy_xx, ts_xxxy_xz, ts_xxxyy_xx, ts_xxxyy_xz, ts_xxxyyy_xx, ts_xxxyyy_xy, ts_xxxyyy_xz, ts_xxxyyy_yy, ts_xxxyyy_yz, ts_xxxyyy_zz, ts_xxyyy_xy, ts_xxyyy_y, ts_xxyyy_yy, ts_xxyyy_yz, ts_xxyyy_zz, ts_xyyy_xy, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyy_xx[i] = 2.0 * ts_xxxy_xx[i] * fe_0 + ts_xxxyy_xx[i] * pa_y[i];

        ts_xxxyyy_xy[i] = 2.0 * ts_xyyy_xy[i] * fe_0 + ts_xxyyy_y[i] * fe_0 + ts_xxyyy_xy[i] * pa_x[i];

        ts_xxxyyy_xz[i] = 2.0 * ts_xxxy_xz[i] * fe_0 + ts_xxxyy_xz[i] * pa_y[i];

        ts_xxxyyy_yy[i] = 2.0 * ts_xyyy_yy[i] * fe_0 + ts_xxyyy_yy[i] * pa_x[i];

        ts_xxxyyy_yz[i] = 2.0 * ts_xyyy_yz[i] * fe_0 + ts_xxyyy_yz[i] * pa_x[i];

        ts_xxxyyy_zz[i] = 2.0 * ts_xyyy_zz[i] * fe_0 + ts_xxyyy_zz[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : ID

    auto ts_xxxyyz_xx = pbuffer.data(idx_ovl_id + 42);

    auto ts_xxxyyz_xy = pbuffer.data(idx_ovl_id + 43);

    auto ts_xxxyyz_xz = pbuffer.data(idx_ovl_id + 44);

    auto ts_xxxyyz_yy = pbuffer.data(idx_ovl_id + 45);

    auto ts_xxxyyz_yz = pbuffer.data(idx_ovl_id + 46);

    auto ts_xxxyyz_zz = pbuffer.data(idx_ovl_id + 47);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxyy_xx, ts_xxxyy_xy, ts_xxxyy_yy, ts_xxxyyz_xx, ts_xxxyyz_xy, ts_xxxyyz_xz, ts_xxxyyz_yy, ts_xxxyyz_yz, ts_xxxyyz_zz, ts_xxxyz_xz, ts_xxxz_xz, ts_xxyyz_yz, ts_xxyyz_zz, ts_xyyz_yz, ts_xyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyz_xx[i] = ts_xxxyy_xx[i] * pa_z[i];

        ts_xxxyyz_xy[i] = ts_xxxyy_xy[i] * pa_z[i];

        ts_xxxyyz_xz[i] = ts_xxxz_xz[i] * fe_0 + ts_xxxyz_xz[i] * pa_y[i];

        ts_xxxyyz_yy[i] = ts_xxxyy_yy[i] * pa_z[i];

        ts_xxxyyz_yz[i] = 2.0 * ts_xyyz_yz[i] * fe_0 + ts_xxyyz_yz[i] * pa_x[i];

        ts_xxxyyz_zz[i] = 2.0 * ts_xyyz_zz[i] * fe_0 + ts_xxyyz_zz[i] * pa_x[i];
    }

    // Set up 48-54 components of targeted buffer : ID

    auto ts_xxxyzz_xx = pbuffer.data(idx_ovl_id + 48);

    auto ts_xxxyzz_xy = pbuffer.data(idx_ovl_id + 49);

    auto ts_xxxyzz_xz = pbuffer.data(idx_ovl_id + 50);

    auto ts_xxxyzz_yy = pbuffer.data(idx_ovl_id + 51);

    auto ts_xxxyzz_yz = pbuffer.data(idx_ovl_id + 52);

    auto ts_xxxyzz_zz = pbuffer.data(idx_ovl_id + 53);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxyzz_xx, ts_xxxyzz_xy, ts_xxxyzz_xz, ts_xxxyzz_yy, ts_xxxyzz_yz, ts_xxxyzz_zz, ts_xxxzz_x, ts_xxxzz_xx, ts_xxxzz_xy, ts_xxxzz_xz, ts_xxxzz_zz, ts_xxyzz_yy, ts_xxyzz_yz, ts_xyzz_yy, ts_xyzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyzz_xx[i] = ts_xxxzz_xx[i] * pa_y[i];

        ts_xxxyzz_xy[i] = ts_xxxzz_x[i] * fe_0 + ts_xxxzz_xy[i] * pa_y[i];

        ts_xxxyzz_xz[i] = ts_xxxzz_xz[i] * pa_y[i];

        ts_xxxyzz_yy[i] = 2.0 * ts_xyzz_yy[i] * fe_0 + ts_xxyzz_yy[i] * pa_x[i];

        ts_xxxyzz_yz[i] = 2.0 * ts_xyzz_yz[i] * fe_0 + ts_xxyzz_yz[i] * pa_x[i];

        ts_xxxyzz_zz[i] = ts_xxxzz_zz[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : ID

    auto ts_xxxzzz_xx = pbuffer.data(idx_ovl_id + 54);

    auto ts_xxxzzz_xy = pbuffer.data(idx_ovl_id + 55);

    auto ts_xxxzzz_xz = pbuffer.data(idx_ovl_id + 56);

    auto ts_xxxzzz_yy = pbuffer.data(idx_ovl_id + 57);

    auto ts_xxxzzz_yz = pbuffer.data(idx_ovl_id + 58);

    auto ts_xxxzzz_zz = pbuffer.data(idx_ovl_id + 59);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxz_xx, ts_xxxz_xy, ts_xxxzz_xx, ts_xxxzz_xy, ts_xxxzzz_xx, ts_xxxzzz_xy, ts_xxxzzz_xz, ts_xxxzzz_yy, ts_xxxzzz_yz, ts_xxxzzz_zz, ts_xxzzz_xz, ts_xxzzz_yy, ts_xxzzz_yz, ts_xxzzz_z, ts_xxzzz_zz, ts_xzzz_xz, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzzz_xx[i] = 2.0 * ts_xxxz_xx[i] * fe_0 + ts_xxxzz_xx[i] * pa_z[i];

        ts_xxxzzz_xy[i] = 2.0 * ts_xxxz_xy[i] * fe_0 + ts_xxxzz_xy[i] * pa_z[i];

        ts_xxxzzz_xz[i] = 2.0 * ts_xzzz_xz[i] * fe_0 + ts_xxzzz_z[i] * fe_0 + ts_xxzzz_xz[i] * pa_x[i];

        ts_xxxzzz_yy[i] = 2.0 * ts_xzzz_yy[i] * fe_0 + ts_xxzzz_yy[i] * pa_x[i];

        ts_xxxzzz_yz[i] = 2.0 * ts_xzzz_yz[i] * fe_0 + ts_xxzzz_yz[i] * pa_x[i];

        ts_xxxzzz_zz[i] = 2.0 * ts_xzzz_zz[i] * fe_0 + ts_xxzzz_zz[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : ID

    auto ts_xxyyyy_xx = pbuffer.data(idx_ovl_id + 60);

    auto ts_xxyyyy_xy = pbuffer.data(idx_ovl_id + 61);

    auto ts_xxyyyy_xz = pbuffer.data(idx_ovl_id + 62);

    auto ts_xxyyyy_yy = pbuffer.data(idx_ovl_id + 63);

    auto ts_xxyyyy_yz = pbuffer.data(idx_ovl_id + 64);

    auto ts_xxyyyy_zz = pbuffer.data(idx_ovl_id + 65);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyy_xx, ts_xxyy_xz, ts_xxyyy_xx, ts_xxyyy_xz, ts_xxyyyy_xx, ts_xxyyyy_xy, ts_xxyyyy_xz, ts_xxyyyy_yy, ts_xxyyyy_yz, ts_xxyyyy_zz, ts_xyyyy_xy, ts_xyyyy_y, ts_xyyyy_yy, ts_xyyyy_yz, ts_xyyyy_zz, ts_yyyy_xy, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyy_xx[i] = 3.0 * ts_xxyy_xx[i] * fe_0 + ts_xxyyy_xx[i] * pa_y[i];

        ts_xxyyyy_xy[i] = ts_yyyy_xy[i] * fe_0 + ts_xyyyy_y[i] * fe_0 + ts_xyyyy_xy[i] * pa_x[i];

        ts_xxyyyy_xz[i] = 3.0 * ts_xxyy_xz[i] * fe_0 + ts_xxyyy_xz[i] * pa_y[i];

        ts_xxyyyy_yy[i] = ts_yyyy_yy[i] * fe_0 + ts_xyyyy_yy[i] * pa_x[i];

        ts_xxyyyy_yz[i] = ts_yyyy_yz[i] * fe_0 + ts_xyyyy_yz[i] * pa_x[i];

        ts_xxyyyy_zz[i] = ts_yyyy_zz[i] * fe_0 + ts_xyyyy_zz[i] * pa_x[i];
    }

    // Set up 66-72 components of targeted buffer : ID

    auto ts_xxyyyz_xx = pbuffer.data(idx_ovl_id + 66);

    auto ts_xxyyyz_xy = pbuffer.data(idx_ovl_id + 67);

    auto ts_xxyyyz_xz = pbuffer.data(idx_ovl_id + 68);

    auto ts_xxyyyz_yy = pbuffer.data(idx_ovl_id + 69);

    auto ts_xxyyyz_yz = pbuffer.data(idx_ovl_id + 70);

    auto ts_xxyyyz_zz = pbuffer.data(idx_ovl_id + 71);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyyy_xx, ts_xxyyy_xy, ts_xxyyy_yy, ts_xxyyyz_xx, ts_xxyyyz_xy, ts_xxyyyz_xz, ts_xxyyyz_yy, ts_xxyyyz_yz, ts_xxyyyz_zz, ts_xxyyz_xz, ts_xxyz_xz, ts_xyyyz_yz, ts_xyyyz_zz, ts_yyyz_yz, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyz_xx[i] = ts_xxyyy_xx[i] * pa_z[i];

        ts_xxyyyz_xy[i] = ts_xxyyy_xy[i] * pa_z[i];

        ts_xxyyyz_xz[i] = 2.0 * ts_xxyz_xz[i] * fe_0 + ts_xxyyz_xz[i] * pa_y[i];

        ts_xxyyyz_yy[i] = ts_xxyyy_yy[i] * pa_z[i];

        ts_xxyyyz_yz[i] = ts_yyyz_yz[i] * fe_0 + ts_xyyyz_yz[i] * pa_x[i];

        ts_xxyyyz_zz[i] = ts_yyyz_zz[i] * fe_0 + ts_xyyyz_zz[i] * pa_x[i];
    }

    // Set up 72-78 components of targeted buffer : ID

    auto ts_xxyyzz_xx = pbuffer.data(idx_ovl_id + 72);

    auto ts_xxyyzz_xy = pbuffer.data(idx_ovl_id + 73);

    auto ts_xxyyzz_xz = pbuffer.data(idx_ovl_id + 74);

    auto ts_xxyyzz_yy = pbuffer.data(idx_ovl_id + 75);

    auto ts_xxyyzz_yz = pbuffer.data(idx_ovl_id + 76);

    auto ts_xxyyzz_zz = pbuffer.data(idx_ovl_id + 77);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyy_xy, ts_xxyyz_xy, ts_xxyyzz_xx, ts_xxyyzz_xy, ts_xxyyzz_xz, ts_xxyyzz_yy, ts_xxyyzz_yz, ts_xxyyzz_zz, ts_xxyzz_xx, ts_xxyzz_xz, ts_xxzz_xx, ts_xxzz_xz, ts_xyyzz_yy, ts_xyyzz_yz, ts_xyyzz_zz, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyzz_xx[i] = ts_xxzz_xx[i] * fe_0 + ts_xxyzz_xx[i] * pa_y[i];

        ts_xxyyzz_xy[i] = ts_xxyy_xy[i] * fe_0 + ts_xxyyz_xy[i] * pa_z[i];

        ts_xxyyzz_xz[i] = ts_xxzz_xz[i] * fe_0 + ts_xxyzz_xz[i] * pa_y[i];

        ts_xxyyzz_yy[i] = ts_yyzz_yy[i] * fe_0 + ts_xyyzz_yy[i] * pa_x[i];

        ts_xxyyzz_yz[i] = ts_yyzz_yz[i] * fe_0 + ts_xyyzz_yz[i] * pa_x[i];

        ts_xxyyzz_zz[i] = ts_yyzz_zz[i] * fe_0 + ts_xyyzz_zz[i] * pa_x[i];
    }

    // Set up 78-84 components of targeted buffer : ID

    auto ts_xxyzzz_xx = pbuffer.data(idx_ovl_id + 78);

    auto ts_xxyzzz_xy = pbuffer.data(idx_ovl_id + 79);

    auto ts_xxyzzz_xz = pbuffer.data(idx_ovl_id + 80);

    auto ts_xxyzzz_yy = pbuffer.data(idx_ovl_id + 81);

    auto ts_xxyzzz_yz = pbuffer.data(idx_ovl_id + 82);

    auto ts_xxyzzz_zz = pbuffer.data(idx_ovl_id + 83);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyzzz_xx, ts_xxyzzz_xy, ts_xxyzzz_xz, ts_xxyzzz_yy, ts_xxyzzz_yz, ts_xxyzzz_zz, ts_xxzzz_x, ts_xxzzz_xx, ts_xxzzz_xy, ts_xxzzz_xz, ts_xxzzz_zz, ts_xyzzz_yy, ts_xyzzz_yz, ts_yzzz_yy, ts_yzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzzz_xx[i] = ts_xxzzz_xx[i] * pa_y[i];

        ts_xxyzzz_xy[i] = ts_xxzzz_x[i] * fe_0 + ts_xxzzz_xy[i] * pa_y[i];

        ts_xxyzzz_xz[i] = ts_xxzzz_xz[i] * pa_y[i];

        ts_xxyzzz_yy[i] = ts_yzzz_yy[i] * fe_0 + ts_xyzzz_yy[i] * pa_x[i];

        ts_xxyzzz_yz[i] = ts_yzzz_yz[i] * fe_0 + ts_xyzzz_yz[i] * pa_x[i];

        ts_xxyzzz_zz[i] = ts_xxzzz_zz[i] * pa_y[i];
    }

    // Set up 84-90 components of targeted buffer : ID

    auto ts_xxzzzz_xx = pbuffer.data(idx_ovl_id + 84);

    auto ts_xxzzzz_xy = pbuffer.data(idx_ovl_id + 85);

    auto ts_xxzzzz_xz = pbuffer.data(idx_ovl_id + 86);

    auto ts_xxzzzz_yy = pbuffer.data(idx_ovl_id + 87);

    auto ts_xxzzzz_yz = pbuffer.data(idx_ovl_id + 88);

    auto ts_xxzzzz_zz = pbuffer.data(idx_ovl_id + 89);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxzz_xx, ts_xxzz_xy, ts_xxzzz_xx, ts_xxzzz_xy, ts_xxzzzz_xx, ts_xxzzzz_xy, ts_xxzzzz_xz, ts_xxzzzz_yy, ts_xxzzzz_yz, ts_xxzzzz_zz, ts_xzzzz_xz, ts_xzzzz_yy, ts_xzzzz_yz, ts_xzzzz_z, ts_xzzzz_zz, ts_zzzz_xz, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzzz_xx[i] = 3.0 * ts_xxzz_xx[i] * fe_0 + ts_xxzzz_xx[i] * pa_z[i];

        ts_xxzzzz_xy[i] = 3.0 * ts_xxzz_xy[i] * fe_0 + ts_xxzzz_xy[i] * pa_z[i];

        ts_xxzzzz_xz[i] = ts_zzzz_xz[i] * fe_0 + ts_xzzzz_z[i] * fe_0 + ts_xzzzz_xz[i] * pa_x[i];

        ts_xxzzzz_yy[i] = ts_zzzz_yy[i] * fe_0 + ts_xzzzz_yy[i] * pa_x[i];

        ts_xxzzzz_yz[i] = ts_zzzz_yz[i] * fe_0 + ts_xzzzz_yz[i] * pa_x[i];

        ts_xxzzzz_zz[i] = ts_zzzz_zz[i] * fe_0 + ts_xzzzz_zz[i] * pa_x[i];
    }

    // Set up 90-96 components of targeted buffer : ID

    auto ts_xyyyyy_xx = pbuffer.data(idx_ovl_id + 90);

    auto ts_xyyyyy_xy = pbuffer.data(idx_ovl_id + 91);

    auto ts_xyyyyy_xz = pbuffer.data(idx_ovl_id + 92);

    auto ts_xyyyyy_yy = pbuffer.data(idx_ovl_id + 93);

    auto ts_xyyyyy_yz = pbuffer.data(idx_ovl_id + 94);

    auto ts_xyyyyy_zz = pbuffer.data(idx_ovl_id + 95);

    #pragma omp simd aligned(pa_x, ts_xyyyyy_xx, ts_xyyyyy_xy, ts_xyyyyy_xz, ts_xyyyyy_yy, ts_xyyyyy_yz, ts_xyyyyy_zz, ts_yyyyy_x, ts_yyyyy_xx, ts_yyyyy_xy, ts_yyyyy_xz, ts_yyyyy_y, ts_yyyyy_yy, ts_yyyyy_yz, ts_yyyyy_z, ts_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyy_xx[i] = 2.0 * ts_yyyyy_x[i] * fe_0 + ts_yyyyy_xx[i] * pa_x[i];

        ts_xyyyyy_xy[i] = ts_yyyyy_y[i] * fe_0 + ts_yyyyy_xy[i] * pa_x[i];

        ts_xyyyyy_xz[i] = ts_yyyyy_z[i] * fe_0 + ts_yyyyy_xz[i] * pa_x[i];

        ts_xyyyyy_yy[i] = ts_yyyyy_yy[i] * pa_x[i];

        ts_xyyyyy_yz[i] = ts_yyyyy_yz[i] * pa_x[i];

        ts_xyyyyy_zz[i] = ts_yyyyy_zz[i] * pa_x[i];
    }

    // Set up 96-102 components of targeted buffer : ID

    auto ts_xyyyyz_xx = pbuffer.data(idx_ovl_id + 96);

    auto ts_xyyyyz_xy = pbuffer.data(idx_ovl_id + 97);

    auto ts_xyyyyz_xz = pbuffer.data(idx_ovl_id + 98);

    auto ts_xyyyyz_yy = pbuffer.data(idx_ovl_id + 99);

    auto ts_xyyyyz_yz = pbuffer.data(idx_ovl_id + 100);

    auto ts_xyyyyz_zz = pbuffer.data(idx_ovl_id + 101);

    #pragma omp simd aligned(pa_x, pa_z, ts_xyyyy_xx, ts_xyyyy_xy, ts_xyyyyz_xx, ts_xyyyyz_xy, ts_xyyyyz_xz, ts_xyyyyz_yy, ts_xyyyyz_yz, ts_xyyyyz_zz, ts_yyyyz_xz, ts_yyyyz_yy, ts_yyyyz_yz, ts_yyyyz_z, ts_yyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyz_xx[i] = ts_xyyyy_xx[i] * pa_z[i];

        ts_xyyyyz_xy[i] = ts_xyyyy_xy[i] * pa_z[i];

        ts_xyyyyz_xz[i] = ts_yyyyz_z[i] * fe_0 + ts_yyyyz_xz[i] * pa_x[i];

        ts_xyyyyz_yy[i] = ts_yyyyz_yy[i] * pa_x[i];

        ts_xyyyyz_yz[i] = ts_yyyyz_yz[i] * pa_x[i];

        ts_xyyyyz_zz[i] = ts_yyyyz_zz[i] * pa_x[i];
    }

    // Set up 102-108 components of targeted buffer : ID

    auto ts_xyyyzz_xx = pbuffer.data(idx_ovl_id + 102);

    auto ts_xyyyzz_xy = pbuffer.data(idx_ovl_id + 103);

    auto ts_xyyyzz_xz = pbuffer.data(idx_ovl_id + 104);

    auto ts_xyyyzz_yy = pbuffer.data(idx_ovl_id + 105);

    auto ts_xyyyzz_yz = pbuffer.data(idx_ovl_id + 106);

    auto ts_xyyyzz_zz = pbuffer.data(idx_ovl_id + 107);

    #pragma omp simd aligned(pa_x, ts_xyyyzz_xx, ts_xyyyzz_xy, ts_xyyyzz_xz, ts_xyyyzz_yy, ts_xyyyzz_yz, ts_xyyyzz_zz, ts_yyyzz_x, ts_yyyzz_xx, ts_yyyzz_xy, ts_yyyzz_xz, ts_yyyzz_y, ts_yyyzz_yy, ts_yyyzz_yz, ts_yyyzz_z, ts_yyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyzz_xx[i] = 2.0 * ts_yyyzz_x[i] * fe_0 + ts_yyyzz_xx[i] * pa_x[i];

        ts_xyyyzz_xy[i] = ts_yyyzz_y[i] * fe_0 + ts_yyyzz_xy[i] * pa_x[i];

        ts_xyyyzz_xz[i] = ts_yyyzz_z[i] * fe_0 + ts_yyyzz_xz[i] * pa_x[i];

        ts_xyyyzz_yy[i] = ts_yyyzz_yy[i] * pa_x[i];

        ts_xyyyzz_yz[i] = ts_yyyzz_yz[i] * pa_x[i];

        ts_xyyyzz_zz[i] = ts_yyyzz_zz[i] * pa_x[i];
    }

    // Set up 108-114 components of targeted buffer : ID

    auto ts_xyyzzz_xx = pbuffer.data(idx_ovl_id + 108);

    auto ts_xyyzzz_xy = pbuffer.data(idx_ovl_id + 109);

    auto ts_xyyzzz_xz = pbuffer.data(idx_ovl_id + 110);

    auto ts_xyyzzz_yy = pbuffer.data(idx_ovl_id + 111);

    auto ts_xyyzzz_yz = pbuffer.data(idx_ovl_id + 112);

    auto ts_xyyzzz_zz = pbuffer.data(idx_ovl_id + 113);

    #pragma omp simd aligned(pa_x, ts_xyyzzz_xx, ts_xyyzzz_xy, ts_xyyzzz_xz, ts_xyyzzz_yy, ts_xyyzzz_yz, ts_xyyzzz_zz, ts_yyzzz_x, ts_yyzzz_xx, ts_yyzzz_xy, ts_yyzzz_xz, ts_yyzzz_y, ts_yyzzz_yy, ts_yyzzz_yz, ts_yyzzz_z, ts_yyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzzz_xx[i] = 2.0 * ts_yyzzz_x[i] * fe_0 + ts_yyzzz_xx[i] * pa_x[i];

        ts_xyyzzz_xy[i] = ts_yyzzz_y[i] * fe_0 + ts_yyzzz_xy[i] * pa_x[i];

        ts_xyyzzz_xz[i] = ts_yyzzz_z[i] * fe_0 + ts_yyzzz_xz[i] * pa_x[i];

        ts_xyyzzz_yy[i] = ts_yyzzz_yy[i] * pa_x[i];

        ts_xyyzzz_yz[i] = ts_yyzzz_yz[i] * pa_x[i];

        ts_xyyzzz_zz[i] = ts_yyzzz_zz[i] * pa_x[i];
    }

    // Set up 114-120 components of targeted buffer : ID

    auto ts_xyzzzz_xx = pbuffer.data(idx_ovl_id + 114);

    auto ts_xyzzzz_xy = pbuffer.data(idx_ovl_id + 115);

    auto ts_xyzzzz_xz = pbuffer.data(idx_ovl_id + 116);

    auto ts_xyzzzz_yy = pbuffer.data(idx_ovl_id + 117);

    auto ts_xyzzzz_yz = pbuffer.data(idx_ovl_id + 118);

    auto ts_xyzzzz_zz = pbuffer.data(idx_ovl_id + 119);

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzzzz_xx, ts_xyzzzz_xy, ts_xyzzzz_xz, ts_xyzzzz_yy, ts_xyzzzz_yz, ts_xyzzzz_zz, ts_xzzzz_xx, ts_xzzzz_xz, ts_yzzzz_xy, ts_yzzzz_y, ts_yzzzz_yy, ts_yzzzz_yz, ts_yzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzzzz_xx[i] = ts_xzzzz_xx[i] * pa_y[i];

        ts_xyzzzz_xy[i] = ts_yzzzz_y[i] * fe_0 + ts_yzzzz_xy[i] * pa_x[i];

        ts_xyzzzz_xz[i] = ts_xzzzz_xz[i] * pa_y[i];

        ts_xyzzzz_yy[i] = ts_yzzzz_yy[i] * pa_x[i];

        ts_xyzzzz_yz[i] = ts_yzzzz_yz[i] * pa_x[i];

        ts_xyzzzz_zz[i] = ts_yzzzz_zz[i] * pa_x[i];
    }

    // Set up 120-126 components of targeted buffer : ID

    auto ts_xzzzzz_xx = pbuffer.data(idx_ovl_id + 120);

    auto ts_xzzzzz_xy = pbuffer.data(idx_ovl_id + 121);

    auto ts_xzzzzz_xz = pbuffer.data(idx_ovl_id + 122);

    auto ts_xzzzzz_yy = pbuffer.data(idx_ovl_id + 123);

    auto ts_xzzzzz_yz = pbuffer.data(idx_ovl_id + 124);

    auto ts_xzzzzz_zz = pbuffer.data(idx_ovl_id + 125);

    #pragma omp simd aligned(pa_x, ts_xzzzzz_xx, ts_xzzzzz_xy, ts_xzzzzz_xz, ts_xzzzzz_yy, ts_xzzzzz_yz, ts_xzzzzz_zz, ts_zzzzz_x, ts_zzzzz_xx, ts_zzzzz_xy, ts_zzzzz_xz, ts_zzzzz_y, ts_zzzzz_yy, ts_zzzzz_yz, ts_zzzzz_z, ts_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzzz_xx[i] = 2.0 * ts_zzzzz_x[i] * fe_0 + ts_zzzzz_xx[i] * pa_x[i];

        ts_xzzzzz_xy[i] = ts_zzzzz_y[i] * fe_0 + ts_zzzzz_xy[i] * pa_x[i];

        ts_xzzzzz_xz[i] = ts_zzzzz_z[i] * fe_0 + ts_zzzzz_xz[i] * pa_x[i];

        ts_xzzzzz_yy[i] = ts_zzzzz_yy[i] * pa_x[i];

        ts_xzzzzz_yz[i] = ts_zzzzz_yz[i] * pa_x[i];

        ts_xzzzzz_zz[i] = ts_zzzzz_zz[i] * pa_x[i];
    }

    // Set up 126-132 components of targeted buffer : ID

    auto ts_yyyyyy_xx = pbuffer.data(idx_ovl_id + 126);

    auto ts_yyyyyy_xy = pbuffer.data(idx_ovl_id + 127);

    auto ts_yyyyyy_xz = pbuffer.data(idx_ovl_id + 128);

    auto ts_yyyyyy_yy = pbuffer.data(idx_ovl_id + 129);

    auto ts_yyyyyy_yz = pbuffer.data(idx_ovl_id + 130);

    auto ts_yyyyyy_zz = pbuffer.data(idx_ovl_id + 131);

    #pragma omp simd aligned(pa_y, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_zz, ts_yyyyy_x, ts_yyyyy_xx, ts_yyyyy_xy, ts_yyyyy_xz, ts_yyyyy_y, ts_yyyyy_yy, ts_yyyyy_yz, ts_yyyyy_z, ts_yyyyy_zz, ts_yyyyyy_xx, ts_yyyyyy_xy, ts_yyyyyy_xz, ts_yyyyyy_yy, ts_yyyyyy_yz, ts_yyyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyy_xx[i] = 5.0 * ts_yyyy_xx[i] * fe_0 + ts_yyyyy_xx[i] * pa_y[i];

        ts_yyyyyy_xy[i] = 5.0 * ts_yyyy_xy[i] * fe_0 + ts_yyyyy_x[i] * fe_0 + ts_yyyyy_xy[i] * pa_y[i];

        ts_yyyyyy_xz[i] = 5.0 * ts_yyyy_xz[i] * fe_0 + ts_yyyyy_xz[i] * pa_y[i];

        ts_yyyyyy_yy[i] = 5.0 * ts_yyyy_yy[i] * fe_0 + 2.0 * ts_yyyyy_y[i] * fe_0 + ts_yyyyy_yy[i] * pa_y[i];

        ts_yyyyyy_yz[i] = 5.0 * ts_yyyy_yz[i] * fe_0 + ts_yyyyy_z[i] * fe_0 + ts_yyyyy_yz[i] * pa_y[i];

        ts_yyyyyy_zz[i] = 5.0 * ts_yyyy_zz[i] * fe_0 + ts_yyyyy_zz[i] * pa_y[i];
    }

    // Set up 132-138 components of targeted buffer : ID

    auto ts_yyyyyz_xx = pbuffer.data(idx_ovl_id + 132);

    auto ts_yyyyyz_xy = pbuffer.data(idx_ovl_id + 133);

    auto ts_yyyyyz_xz = pbuffer.data(idx_ovl_id + 134);

    auto ts_yyyyyz_yy = pbuffer.data(idx_ovl_id + 135);

    auto ts_yyyyyz_yz = pbuffer.data(idx_ovl_id + 136);

    auto ts_yyyyyz_zz = pbuffer.data(idx_ovl_id + 137);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyyy_xx, ts_yyyyy_xy, ts_yyyyy_y, ts_yyyyy_yy, ts_yyyyy_yz, ts_yyyyyz_xx, ts_yyyyyz_xy, ts_yyyyyz_xz, ts_yyyyyz_yy, ts_yyyyyz_yz, ts_yyyyyz_zz, ts_yyyyz_xz, ts_yyyyz_zz, ts_yyyz_xz, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyz_xx[i] = ts_yyyyy_xx[i] * pa_z[i];

        ts_yyyyyz_xy[i] = ts_yyyyy_xy[i] * pa_z[i];

        ts_yyyyyz_xz[i] = 4.0 * ts_yyyz_xz[i] * fe_0 + ts_yyyyz_xz[i] * pa_y[i];

        ts_yyyyyz_yy[i] = ts_yyyyy_yy[i] * pa_z[i];

        ts_yyyyyz_yz[i] = ts_yyyyy_y[i] * fe_0 + ts_yyyyy_yz[i] * pa_z[i];

        ts_yyyyyz_zz[i] = 4.0 * ts_yyyz_zz[i] * fe_0 + ts_yyyyz_zz[i] * pa_y[i];
    }

    // Set up 138-144 components of targeted buffer : ID

    auto ts_yyyyzz_xx = pbuffer.data(idx_ovl_id + 138);

    auto ts_yyyyzz_xy = pbuffer.data(idx_ovl_id + 139);

    auto ts_yyyyzz_xz = pbuffer.data(idx_ovl_id + 140);

    auto ts_yyyyzz_yy = pbuffer.data(idx_ovl_id + 141);

    auto ts_yyyyzz_yz = pbuffer.data(idx_ovl_id + 142);

    auto ts_yyyyzz_zz = pbuffer.data(idx_ovl_id + 143);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyy_xy, ts_yyyy_yy, ts_yyyyz_xy, ts_yyyyz_yy, ts_yyyyzz_xx, ts_yyyyzz_xy, ts_yyyyzz_xz, ts_yyyyzz_yy, ts_yyyyzz_yz, ts_yyyyzz_zz, ts_yyyzz_xx, ts_yyyzz_xz, ts_yyyzz_yz, ts_yyyzz_z, ts_yyyzz_zz, ts_yyzz_xx, ts_yyzz_xz, ts_yyzz_yz, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyzz_xx[i] = 3.0 * ts_yyzz_xx[i] * fe_0 + ts_yyyzz_xx[i] * pa_y[i];

        ts_yyyyzz_xy[i] = ts_yyyy_xy[i] * fe_0 + ts_yyyyz_xy[i] * pa_z[i];

        ts_yyyyzz_xz[i] = 3.0 * ts_yyzz_xz[i] * fe_0 + ts_yyyzz_xz[i] * pa_y[i];

        ts_yyyyzz_yy[i] = ts_yyyy_yy[i] * fe_0 + ts_yyyyz_yy[i] * pa_z[i];

        ts_yyyyzz_yz[i] = 3.0 * ts_yyzz_yz[i] * fe_0 + ts_yyyzz_z[i] * fe_0 + ts_yyyzz_yz[i] * pa_y[i];

        ts_yyyyzz_zz[i] = 3.0 * ts_yyzz_zz[i] * fe_0 + ts_yyyzz_zz[i] * pa_y[i];
    }

    // Set up 144-150 components of targeted buffer : ID

    auto ts_yyyzzz_xx = pbuffer.data(idx_ovl_id + 144);

    auto ts_yyyzzz_xy = pbuffer.data(idx_ovl_id + 145);

    auto ts_yyyzzz_xz = pbuffer.data(idx_ovl_id + 146);

    auto ts_yyyzzz_yy = pbuffer.data(idx_ovl_id + 147);

    auto ts_yyyzzz_yz = pbuffer.data(idx_ovl_id + 148);

    auto ts_yyyzzz_zz = pbuffer.data(idx_ovl_id + 149);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyz_xy, ts_yyyz_yy, ts_yyyzz_xy, ts_yyyzz_yy, ts_yyyzzz_xx, ts_yyyzzz_xy, ts_yyyzzz_xz, ts_yyyzzz_yy, ts_yyyzzz_yz, ts_yyyzzz_zz, ts_yyzzz_xx, ts_yyzzz_xz, ts_yyzzz_yz, ts_yyzzz_z, ts_yyzzz_zz, ts_yzzz_xx, ts_yzzz_xz, ts_yzzz_yz, ts_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzzz_xx[i] = 2.0 * ts_yzzz_xx[i] * fe_0 + ts_yyzzz_xx[i] * pa_y[i];

        ts_yyyzzz_xy[i] = 2.0 * ts_yyyz_xy[i] * fe_0 + ts_yyyzz_xy[i] * pa_z[i];

        ts_yyyzzz_xz[i] = 2.0 * ts_yzzz_xz[i] * fe_0 + ts_yyzzz_xz[i] * pa_y[i];

        ts_yyyzzz_yy[i] = 2.0 * ts_yyyz_yy[i] * fe_0 + ts_yyyzz_yy[i] * pa_z[i];

        ts_yyyzzz_yz[i] = 2.0 * ts_yzzz_yz[i] * fe_0 + ts_yyzzz_z[i] * fe_0 + ts_yyzzz_yz[i] * pa_y[i];

        ts_yyyzzz_zz[i] = 2.0 * ts_yzzz_zz[i] * fe_0 + ts_yyzzz_zz[i] * pa_y[i];
    }

    // Set up 150-156 components of targeted buffer : ID

    auto ts_yyzzzz_xx = pbuffer.data(idx_ovl_id + 150);

    auto ts_yyzzzz_xy = pbuffer.data(idx_ovl_id + 151);

    auto ts_yyzzzz_xz = pbuffer.data(idx_ovl_id + 152);

    auto ts_yyzzzz_yy = pbuffer.data(idx_ovl_id + 153);

    auto ts_yyzzzz_yz = pbuffer.data(idx_ovl_id + 154);

    auto ts_yyzzzz_zz = pbuffer.data(idx_ovl_id + 155);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyzz_xy, ts_yyzz_yy, ts_yyzzz_xy, ts_yyzzz_yy, ts_yyzzzz_xx, ts_yyzzzz_xy, ts_yyzzzz_xz, ts_yyzzzz_yy, ts_yyzzzz_yz, ts_yyzzzz_zz, ts_yzzzz_xx, ts_yzzzz_xz, ts_yzzzz_yz, ts_yzzzz_z, ts_yzzzz_zz, ts_zzzz_xx, ts_zzzz_xz, ts_zzzz_yz, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzzz_xx[i] = ts_zzzz_xx[i] * fe_0 + ts_yzzzz_xx[i] * pa_y[i];

        ts_yyzzzz_xy[i] = 3.0 * ts_yyzz_xy[i] * fe_0 + ts_yyzzz_xy[i] * pa_z[i];

        ts_yyzzzz_xz[i] = ts_zzzz_xz[i] * fe_0 + ts_yzzzz_xz[i] * pa_y[i];

        ts_yyzzzz_yy[i] = 3.0 * ts_yyzz_yy[i] * fe_0 + ts_yyzzz_yy[i] * pa_z[i];

        ts_yyzzzz_yz[i] = ts_zzzz_yz[i] * fe_0 + ts_yzzzz_z[i] * fe_0 + ts_yzzzz_yz[i] * pa_y[i];

        ts_yyzzzz_zz[i] = ts_zzzz_zz[i] * fe_0 + ts_yzzzz_zz[i] * pa_y[i];
    }

    // Set up 156-162 components of targeted buffer : ID

    auto ts_yzzzzz_xx = pbuffer.data(idx_ovl_id + 156);

    auto ts_yzzzzz_xy = pbuffer.data(idx_ovl_id + 157);

    auto ts_yzzzzz_xz = pbuffer.data(idx_ovl_id + 158);

    auto ts_yzzzzz_yy = pbuffer.data(idx_ovl_id + 159);

    auto ts_yzzzzz_yz = pbuffer.data(idx_ovl_id + 160);

    auto ts_yzzzzz_zz = pbuffer.data(idx_ovl_id + 161);

    #pragma omp simd aligned(pa_y, ts_yzzzzz_xx, ts_yzzzzz_xy, ts_yzzzzz_xz, ts_yzzzzz_yy, ts_yzzzzz_yz, ts_yzzzzz_zz, ts_zzzzz_x, ts_zzzzz_xx, ts_zzzzz_xy, ts_zzzzz_xz, ts_zzzzz_y, ts_zzzzz_yy, ts_zzzzz_yz, ts_zzzzz_z, ts_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzzz_xx[i] = ts_zzzzz_xx[i] * pa_y[i];

        ts_yzzzzz_xy[i] = ts_zzzzz_x[i] * fe_0 + ts_zzzzz_xy[i] * pa_y[i];

        ts_yzzzzz_xz[i] = ts_zzzzz_xz[i] * pa_y[i];

        ts_yzzzzz_yy[i] = 2.0 * ts_zzzzz_y[i] * fe_0 + ts_zzzzz_yy[i] * pa_y[i];

        ts_yzzzzz_yz[i] = ts_zzzzz_z[i] * fe_0 + ts_zzzzz_yz[i] * pa_y[i];

        ts_yzzzzz_zz[i] = ts_zzzzz_zz[i] * pa_y[i];
    }

    // Set up 162-168 components of targeted buffer : ID

    auto ts_zzzzzz_xx = pbuffer.data(idx_ovl_id + 162);

    auto ts_zzzzzz_xy = pbuffer.data(idx_ovl_id + 163);

    auto ts_zzzzzz_xz = pbuffer.data(idx_ovl_id + 164);

    auto ts_zzzzzz_yy = pbuffer.data(idx_ovl_id + 165);

    auto ts_zzzzzz_yz = pbuffer.data(idx_ovl_id + 166);

    auto ts_zzzzzz_zz = pbuffer.data(idx_ovl_id + 167);

    #pragma omp simd aligned(pa_z, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_zz, ts_zzzzz_x, ts_zzzzz_xx, ts_zzzzz_xy, ts_zzzzz_xz, ts_zzzzz_y, ts_zzzzz_yy, ts_zzzzz_yz, ts_zzzzz_z, ts_zzzzz_zz, ts_zzzzzz_xx, ts_zzzzzz_xy, ts_zzzzzz_xz, ts_zzzzzz_yy, ts_zzzzzz_yz, ts_zzzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzzz_xx[i] = 5.0 * ts_zzzz_xx[i] * fe_0 + ts_zzzzz_xx[i] * pa_z[i];

        ts_zzzzzz_xy[i] = 5.0 * ts_zzzz_xy[i] * fe_0 + ts_zzzzz_xy[i] * pa_z[i];

        ts_zzzzzz_xz[i] = 5.0 * ts_zzzz_xz[i] * fe_0 + ts_zzzzz_x[i] * fe_0 + ts_zzzzz_xz[i] * pa_z[i];

        ts_zzzzzz_yy[i] = 5.0 * ts_zzzz_yy[i] * fe_0 + ts_zzzzz_yy[i] * pa_z[i];

        ts_zzzzzz_yz[i] = 5.0 * ts_zzzz_yz[i] * fe_0 + ts_zzzzz_y[i] * fe_0 + ts_zzzzz_yz[i] * pa_z[i];

        ts_zzzzzz_zz[i] = 5.0 * ts_zzzz_zz[i] * fe_0 + 2.0 * ts_zzzzz_z[i] * fe_0 + ts_zzzzz_zz[i] * pa_z[i];
    }

}

} // ovlrec namespace

