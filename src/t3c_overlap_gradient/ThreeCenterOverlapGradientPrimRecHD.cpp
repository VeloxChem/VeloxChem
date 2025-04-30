#include "ThreeCenterOverlapGradientPrimRecHD.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_hd(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_hd,
                              const size_t idx_gd,
                              const size_t idx_hp,
                              const size_t idx_hd,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_gd + 5);

    auto ts_xxxy_xx = pbuffer.data(idx_gd + 6);

    auto ts_xxxy_xy = pbuffer.data(idx_gd + 7);

    auto ts_xxxy_xz = pbuffer.data(idx_gd + 8);

    auto ts_xxxy_yy = pbuffer.data(idx_gd + 9);

    auto ts_xxxy_yz = pbuffer.data(idx_gd + 10);

    auto ts_xxxy_zz = pbuffer.data(idx_gd + 11);

    auto ts_xxxz_xx = pbuffer.data(idx_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_gd + 13);

    auto ts_xxxz_xz = pbuffer.data(idx_gd + 14);

    auto ts_xxxz_yy = pbuffer.data(idx_gd + 15);

    auto ts_xxxz_yz = pbuffer.data(idx_gd + 16);

    auto ts_xxxz_zz = pbuffer.data(idx_gd + 17);

    auto ts_xxyy_xx = pbuffer.data(idx_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_gd + 23);

    auto ts_xxyz_xx = pbuffer.data(idx_gd + 24);

    auto ts_xxyz_xy = pbuffer.data(idx_gd + 25);

    auto ts_xxyz_xz = pbuffer.data(idx_gd + 26);

    auto ts_xxyz_yy = pbuffer.data(idx_gd + 27);

    auto ts_xxyz_yz = pbuffer.data(idx_gd + 28);

    auto ts_xxyz_zz = pbuffer.data(idx_gd + 29);

    auto ts_xxzz_xx = pbuffer.data(idx_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_gd + 35);

    auto ts_xyyy_xx = pbuffer.data(idx_gd + 36);

    auto ts_xyyy_xy = pbuffer.data(idx_gd + 37);

    auto ts_xyyy_xz = pbuffer.data(idx_gd + 38);

    auto ts_xyyy_yy = pbuffer.data(idx_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_gd + 41);

    auto ts_xyyz_xx = pbuffer.data(idx_gd + 42);

    auto ts_xyyz_xy = pbuffer.data(idx_gd + 43);

    auto ts_xyyz_xz = pbuffer.data(idx_gd + 44);

    auto ts_xyyz_yy = pbuffer.data(idx_gd + 45);

    auto ts_xyyz_yz = pbuffer.data(idx_gd + 46);

    auto ts_xyyz_zz = pbuffer.data(idx_gd + 47);

    auto ts_xyzz_xx = pbuffer.data(idx_gd + 48);

    auto ts_xyzz_xy = pbuffer.data(idx_gd + 49);

    auto ts_xyzz_xz = pbuffer.data(idx_gd + 50);

    auto ts_xyzz_yy = pbuffer.data(idx_gd + 51);

    auto ts_xyzz_yz = pbuffer.data(idx_gd + 52);

    auto ts_xyzz_zz = pbuffer.data(idx_gd + 53);

    auto ts_xzzz_xx = pbuffer.data(idx_gd + 54);

    auto ts_xzzz_xy = pbuffer.data(idx_gd + 55);

    auto ts_xzzz_xz = pbuffer.data(idx_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_gd + 65);

    auto ts_yyyz_xx = pbuffer.data(idx_gd + 66);

    auto ts_yyyz_xy = pbuffer.data(idx_gd + 67);

    auto ts_yyyz_xz = pbuffer.data(idx_gd + 68);

    auto ts_yyyz_yy = pbuffer.data(idx_gd + 69);

    auto ts_yyyz_yz = pbuffer.data(idx_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_gd + 71);

    auto ts_yyzz_xx = pbuffer.data(idx_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_gd + 77);

    auto ts_yzzz_xx = pbuffer.data(idx_gd + 78);

    auto ts_yzzz_xy = pbuffer.data(idx_gd + 79);

    auto ts_yzzz_xz = pbuffer.data(idx_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_gd + 89);

    // Set up components of auxiliary buffer : HP

    auto ts_xxxxx_x = pbuffer.data(idx_hp);

    auto ts_xxxxx_y = pbuffer.data(idx_hp + 1);

    auto ts_xxxxx_z = pbuffer.data(idx_hp + 2);

    auto ts_xxxxy_x = pbuffer.data(idx_hp + 3);

    auto ts_xxxxy_y = pbuffer.data(idx_hp + 4);

    auto ts_xxxxy_z = pbuffer.data(idx_hp + 5);

    auto ts_xxxxz_x = pbuffer.data(idx_hp + 6);

    auto ts_xxxxz_y = pbuffer.data(idx_hp + 7);

    auto ts_xxxxz_z = pbuffer.data(idx_hp + 8);

    auto ts_xxxyy_x = pbuffer.data(idx_hp + 9);

    auto ts_xxxyy_y = pbuffer.data(idx_hp + 10);

    auto ts_xxxyy_z = pbuffer.data(idx_hp + 11);

    auto ts_xxxyz_x = pbuffer.data(idx_hp + 12);

    auto ts_xxxyz_y = pbuffer.data(idx_hp + 13);

    auto ts_xxxyz_z = pbuffer.data(idx_hp + 14);

    auto ts_xxxzz_x = pbuffer.data(idx_hp + 15);

    auto ts_xxxzz_y = pbuffer.data(idx_hp + 16);

    auto ts_xxxzz_z = pbuffer.data(idx_hp + 17);

    auto ts_xxyyy_x = pbuffer.data(idx_hp + 18);

    auto ts_xxyyy_y = pbuffer.data(idx_hp + 19);

    auto ts_xxyyy_z = pbuffer.data(idx_hp + 20);

    auto ts_xxyyz_x = pbuffer.data(idx_hp + 21);

    auto ts_xxyyz_y = pbuffer.data(idx_hp + 22);

    auto ts_xxyyz_z = pbuffer.data(idx_hp + 23);

    auto ts_xxyzz_x = pbuffer.data(idx_hp + 24);

    auto ts_xxyzz_y = pbuffer.data(idx_hp + 25);

    auto ts_xxyzz_z = pbuffer.data(idx_hp + 26);

    auto ts_xxzzz_x = pbuffer.data(idx_hp + 27);

    auto ts_xxzzz_y = pbuffer.data(idx_hp + 28);

    auto ts_xxzzz_z = pbuffer.data(idx_hp + 29);

    auto ts_xyyyy_x = pbuffer.data(idx_hp + 30);

    auto ts_xyyyy_y = pbuffer.data(idx_hp + 31);

    auto ts_xyyyy_z = pbuffer.data(idx_hp + 32);

    auto ts_xyyyz_x = pbuffer.data(idx_hp + 33);

    auto ts_xyyyz_y = pbuffer.data(idx_hp + 34);

    auto ts_xyyyz_z = pbuffer.data(idx_hp + 35);

    auto ts_xyyzz_x = pbuffer.data(idx_hp + 36);

    auto ts_xyyzz_y = pbuffer.data(idx_hp + 37);

    auto ts_xyyzz_z = pbuffer.data(idx_hp + 38);

    auto ts_xyzzz_x = pbuffer.data(idx_hp + 39);

    auto ts_xyzzz_y = pbuffer.data(idx_hp + 40);

    auto ts_xyzzz_z = pbuffer.data(idx_hp + 41);

    auto ts_xzzzz_x = pbuffer.data(idx_hp + 42);

    auto ts_xzzzz_y = pbuffer.data(idx_hp + 43);

    auto ts_xzzzz_z = pbuffer.data(idx_hp + 44);

    auto ts_yyyyy_x = pbuffer.data(idx_hp + 45);

    auto ts_yyyyy_y = pbuffer.data(idx_hp + 46);

    auto ts_yyyyy_z = pbuffer.data(idx_hp + 47);

    auto ts_yyyyz_x = pbuffer.data(idx_hp + 48);

    auto ts_yyyyz_y = pbuffer.data(idx_hp + 49);

    auto ts_yyyyz_z = pbuffer.data(idx_hp + 50);

    auto ts_yyyzz_x = pbuffer.data(idx_hp + 51);

    auto ts_yyyzz_y = pbuffer.data(idx_hp + 52);

    auto ts_yyyzz_z = pbuffer.data(idx_hp + 53);

    auto ts_yyzzz_x = pbuffer.data(idx_hp + 54);

    auto ts_yyzzz_y = pbuffer.data(idx_hp + 55);

    auto ts_yyzzz_z = pbuffer.data(idx_hp + 56);

    auto ts_yzzzz_x = pbuffer.data(idx_hp + 57);

    auto ts_yzzzz_y = pbuffer.data(idx_hp + 58);

    auto ts_yzzzz_z = pbuffer.data(idx_hp + 59);

    auto ts_zzzzz_x = pbuffer.data(idx_hp + 60);

    auto ts_zzzzz_y = pbuffer.data(idx_hp + 61);

    auto ts_zzzzz_z = pbuffer.data(idx_hp + 62);

    // Set up components of auxiliary buffer : HD

    auto ts_xxxxx_xx = pbuffer.data(idx_hd);

    auto ts_xxxxx_xy = pbuffer.data(idx_hd + 1);

    auto ts_xxxxx_xz = pbuffer.data(idx_hd + 2);

    auto ts_xxxxx_yy = pbuffer.data(idx_hd + 3);

    auto ts_xxxxx_yz = pbuffer.data(idx_hd + 4);

    auto ts_xxxxx_zz = pbuffer.data(idx_hd + 5);

    auto ts_xxxxy_xx = pbuffer.data(idx_hd + 6);

    auto ts_xxxxy_xy = pbuffer.data(idx_hd + 7);

    auto ts_xxxxy_xz = pbuffer.data(idx_hd + 8);

    auto ts_xxxxy_yy = pbuffer.data(idx_hd + 9);

    auto ts_xxxxy_yz = pbuffer.data(idx_hd + 10);

    auto ts_xxxxy_zz = pbuffer.data(idx_hd + 11);

    auto ts_xxxxz_xx = pbuffer.data(idx_hd + 12);

    auto ts_xxxxz_xy = pbuffer.data(idx_hd + 13);

    auto ts_xxxxz_xz = pbuffer.data(idx_hd + 14);

    auto ts_xxxxz_yy = pbuffer.data(idx_hd + 15);

    auto ts_xxxxz_yz = pbuffer.data(idx_hd + 16);

    auto ts_xxxxz_zz = pbuffer.data(idx_hd + 17);

    auto ts_xxxyy_xx = pbuffer.data(idx_hd + 18);

    auto ts_xxxyy_xy = pbuffer.data(idx_hd + 19);

    auto ts_xxxyy_xz = pbuffer.data(idx_hd + 20);

    auto ts_xxxyy_yy = pbuffer.data(idx_hd + 21);

    auto ts_xxxyy_yz = pbuffer.data(idx_hd + 22);

    auto ts_xxxyy_zz = pbuffer.data(idx_hd + 23);

    auto ts_xxxyz_xx = pbuffer.data(idx_hd + 24);

    auto ts_xxxyz_xy = pbuffer.data(idx_hd + 25);

    auto ts_xxxyz_xz = pbuffer.data(idx_hd + 26);

    auto ts_xxxyz_yy = pbuffer.data(idx_hd + 27);

    auto ts_xxxyz_yz = pbuffer.data(idx_hd + 28);

    auto ts_xxxyz_zz = pbuffer.data(idx_hd + 29);

    auto ts_xxxzz_xx = pbuffer.data(idx_hd + 30);

    auto ts_xxxzz_xy = pbuffer.data(idx_hd + 31);

    auto ts_xxxzz_xz = pbuffer.data(idx_hd + 32);

    auto ts_xxxzz_yy = pbuffer.data(idx_hd + 33);

    auto ts_xxxzz_yz = pbuffer.data(idx_hd + 34);

    auto ts_xxxzz_zz = pbuffer.data(idx_hd + 35);

    auto ts_xxyyy_xx = pbuffer.data(idx_hd + 36);

    auto ts_xxyyy_xy = pbuffer.data(idx_hd + 37);

    auto ts_xxyyy_xz = pbuffer.data(idx_hd + 38);

    auto ts_xxyyy_yy = pbuffer.data(idx_hd + 39);

    auto ts_xxyyy_yz = pbuffer.data(idx_hd + 40);

    auto ts_xxyyy_zz = pbuffer.data(idx_hd + 41);

    auto ts_xxyyz_xx = pbuffer.data(idx_hd + 42);

    auto ts_xxyyz_xy = pbuffer.data(idx_hd + 43);

    auto ts_xxyyz_xz = pbuffer.data(idx_hd + 44);

    auto ts_xxyyz_yy = pbuffer.data(idx_hd + 45);

    auto ts_xxyyz_yz = pbuffer.data(idx_hd + 46);

    auto ts_xxyyz_zz = pbuffer.data(idx_hd + 47);

    auto ts_xxyzz_xx = pbuffer.data(idx_hd + 48);

    auto ts_xxyzz_xy = pbuffer.data(idx_hd + 49);

    auto ts_xxyzz_xz = pbuffer.data(idx_hd + 50);

    auto ts_xxyzz_yy = pbuffer.data(idx_hd + 51);

    auto ts_xxyzz_yz = pbuffer.data(idx_hd + 52);

    auto ts_xxyzz_zz = pbuffer.data(idx_hd + 53);

    auto ts_xxzzz_xx = pbuffer.data(idx_hd + 54);

    auto ts_xxzzz_xy = pbuffer.data(idx_hd + 55);

    auto ts_xxzzz_xz = pbuffer.data(idx_hd + 56);

    auto ts_xxzzz_yy = pbuffer.data(idx_hd + 57);

    auto ts_xxzzz_yz = pbuffer.data(idx_hd + 58);

    auto ts_xxzzz_zz = pbuffer.data(idx_hd + 59);

    auto ts_xyyyy_xx = pbuffer.data(idx_hd + 60);

    auto ts_xyyyy_xy = pbuffer.data(idx_hd + 61);

    auto ts_xyyyy_xz = pbuffer.data(idx_hd + 62);

    auto ts_xyyyy_yy = pbuffer.data(idx_hd + 63);

    auto ts_xyyyy_yz = pbuffer.data(idx_hd + 64);

    auto ts_xyyyy_zz = pbuffer.data(idx_hd + 65);

    auto ts_xyyyz_xx = pbuffer.data(idx_hd + 66);

    auto ts_xyyyz_xy = pbuffer.data(idx_hd + 67);

    auto ts_xyyyz_xz = pbuffer.data(idx_hd + 68);

    auto ts_xyyyz_yy = pbuffer.data(idx_hd + 69);

    auto ts_xyyyz_yz = pbuffer.data(idx_hd + 70);

    auto ts_xyyyz_zz = pbuffer.data(idx_hd + 71);

    auto ts_xyyzz_xx = pbuffer.data(idx_hd + 72);

    auto ts_xyyzz_xy = pbuffer.data(idx_hd + 73);

    auto ts_xyyzz_xz = pbuffer.data(idx_hd + 74);

    auto ts_xyyzz_yy = pbuffer.data(idx_hd + 75);

    auto ts_xyyzz_yz = pbuffer.data(idx_hd + 76);

    auto ts_xyyzz_zz = pbuffer.data(idx_hd + 77);

    auto ts_xyzzz_xx = pbuffer.data(idx_hd + 78);

    auto ts_xyzzz_xy = pbuffer.data(idx_hd + 79);

    auto ts_xyzzz_xz = pbuffer.data(idx_hd + 80);

    auto ts_xyzzz_yy = pbuffer.data(idx_hd + 81);

    auto ts_xyzzz_yz = pbuffer.data(idx_hd + 82);

    auto ts_xyzzz_zz = pbuffer.data(idx_hd + 83);

    auto ts_xzzzz_xx = pbuffer.data(idx_hd + 84);

    auto ts_xzzzz_xy = pbuffer.data(idx_hd + 85);

    auto ts_xzzzz_xz = pbuffer.data(idx_hd + 86);

    auto ts_xzzzz_yy = pbuffer.data(idx_hd + 87);

    auto ts_xzzzz_yz = pbuffer.data(idx_hd + 88);

    auto ts_xzzzz_zz = pbuffer.data(idx_hd + 89);

    auto ts_yyyyy_xx = pbuffer.data(idx_hd + 90);

    auto ts_yyyyy_xy = pbuffer.data(idx_hd + 91);

    auto ts_yyyyy_xz = pbuffer.data(idx_hd + 92);

    auto ts_yyyyy_yy = pbuffer.data(idx_hd + 93);

    auto ts_yyyyy_yz = pbuffer.data(idx_hd + 94);

    auto ts_yyyyy_zz = pbuffer.data(idx_hd + 95);

    auto ts_yyyyz_xx = pbuffer.data(idx_hd + 96);

    auto ts_yyyyz_xy = pbuffer.data(idx_hd + 97);

    auto ts_yyyyz_xz = pbuffer.data(idx_hd + 98);

    auto ts_yyyyz_yy = pbuffer.data(idx_hd + 99);

    auto ts_yyyyz_yz = pbuffer.data(idx_hd + 100);

    auto ts_yyyyz_zz = pbuffer.data(idx_hd + 101);

    auto ts_yyyzz_xx = pbuffer.data(idx_hd + 102);

    auto ts_yyyzz_xy = pbuffer.data(idx_hd + 103);

    auto ts_yyyzz_xz = pbuffer.data(idx_hd + 104);

    auto ts_yyyzz_yy = pbuffer.data(idx_hd + 105);

    auto ts_yyyzz_yz = pbuffer.data(idx_hd + 106);

    auto ts_yyyzz_zz = pbuffer.data(idx_hd + 107);

    auto ts_yyzzz_xx = pbuffer.data(idx_hd + 108);

    auto ts_yyzzz_xy = pbuffer.data(idx_hd + 109);

    auto ts_yyzzz_xz = pbuffer.data(idx_hd + 110);

    auto ts_yyzzz_yy = pbuffer.data(idx_hd + 111);

    auto ts_yyzzz_yz = pbuffer.data(idx_hd + 112);

    auto ts_yyzzz_zz = pbuffer.data(idx_hd + 113);

    auto ts_yzzzz_xx = pbuffer.data(idx_hd + 114);

    auto ts_yzzzz_xy = pbuffer.data(idx_hd + 115);

    auto ts_yzzzz_xz = pbuffer.data(idx_hd + 116);

    auto ts_yzzzz_yy = pbuffer.data(idx_hd + 117);

    auto ts_yzzzz_yz = pbuffer.data(idx_hd + 118);

    auto ts_yzzzz_zz = pbuffer.data(idx_hd + 119);

    auto ts_zzzzz_xx = pbuffer.data(idx_hd + 120);

    auto ts_zzzzz_xy = pbuffer.data(idx_hd + 121);

    auto ts_zzzzz_xz = pbuffer.data(idx_hd + 122);

    auto ts_zzzzz_yy = pbuffer.data(idx_hd + 123);

    auto ts_zzzzz_yz = pbuffer.data(idx_hd + 124);

    auto ts_zzzzz_zz = pbuffer.data(idx_hd + 125);

    // Set up 0-6 components of targeted buffer : HD

    auto gs_x_xxxxx_xx = pbuffer.data(idx_g_hd);

    auto gs_x_xxxxx_xy = pbuffer.data(idx_g_hd + 1);

    auto gs_x_xxxxx_xz = pbuffer.data(idx_g_hd + 2);

    auto gs_x_xxxxx_yy = pbuffer.data(idx_g_hd + 3);

    auto gs_x_xxxxx_yz = pbuffer.data(idx_g_hd + 4);

    auto gs_x_xxxxx_zz = pbuffer.data(idx_g_hd + 5);

    #pragma omp simd aligned(gc_x, gs_x_xxxxx_xx, gs_x_xxxxx_xy, gs_x_xxxxx_xz, gs_x_xxxxx_yy, gs_x_xxxxx_yz, gs_x_xxxxx_zz, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_zz, ts_xxxxx_x, ts_xxxxx_xx, ts_xxxxx_xy, ts_xxxxx_xz, ts_xxxxx_y, ts_xxxxx_yy, ts_xxxxx_yz, ts_xxxxx_z, ts_xxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxx_xx[i] = 10.0 * ts_xxxx_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xx[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xy[i] = 10.0 * ts_xxxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_xz[i] = 10.0 * ts_xxxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yy[i] = 10.0 * ts_xxxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yy[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_yz[i] = 10.0 * ts_xxxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yz[i] * gc_x[i] * tce_0;

        gs_x_xxxxx_zz[i] = 10.0 * ts_xxxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 6-12 components of targeted buffer : HD

    auto gs_x_xxxxy_xx = pbuffer.data(idx_g_hd + 6);

    auto gs_x_xxxxy_xy = pbuffer.data(idx_g_hd + 7);

    auto gs_x_xxxxy_xz = pbuffer.data(idx_g_hd + 8);

    auto gs_x_xxxxy_yy = pbuffer.data(idx_g_hd + 9);

    auto gs_x_xxxxy_yz = pbuffer.data(idx_g_hd + 10);

    auto gs_x_xxxxy_zz = pbuffer.data(idx_g_hd + 11);

    #pragma omp simd aligned(gc_x, gs_x_xxxxy_xx, gs_x_xxxxy_xy, gs_x_xxxxy_xz, gs_x_xxxxy_yy, gs_x_xxxxy_yz, gs_x_xxxxy_zz, ts_xxxxy_x, ts_xxxxy_xx, ts_xxxxy_xy, ts_xxxxy_xz, ts_xxxxy_y, ts_xxxxy_yy, ts_xxxxy_yz, ts_xxxxy_z, ts_xxxxy_zz, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxy_xx[i] = 8.0 * ts_xxxy_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xx[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xy[i] = 8.0 * ts_xxxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_xz[i] = 8.0 * ts_xxxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yy[i] = 8.0 * ts_xxxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yy[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_yz[i] = 8.0 * ts_xxxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yz[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_zz[i] = 8.0 * ts_xxxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 12-18 components of targeted buffer : HD

    auto gs_x_xxxxz_xx = pbuffer.data(idx_g_hd + 12);

    auto gs_x_xxxxz_xy = pbuffer.data(idx_g_hd + 13);

    auto gs_x_xxxxz_xz = pbuffer.data(idx_g_hd + 14);

    auto gs_x_xxxxz_yy = pbuffer.data(idx_g_hd + 15);

    auto gs_x_xxxxz_yz = pbuffer.data(idx_g_hd + 16);

    auto gs_x_xxxxz_zz = pbuffer.data(idx_g_hd + 17);

    #pragma omp simd aligned(gc_x, gs_x_xxxxz_xx, gs_x_xxxxz_xy, gs_x_xxxxz_xz, gs_x_xxxxz_yy, gs_x_xxxxz_yz, gs_x_xxxxz_zz, ts_xxxxz_x, ts_xxxxz_xx, ts_xxxxz_xy, ts_xxxxz_xz, ts_xxxxz_y, ts_xxxxz_yy, ts_xxxxz_yz, ts_xxxxz_z, ts_xxxxz_zz, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxz_xx[i] = 8.0 * ts_xxxz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xy[i] = 8.0 * ts_xxxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_xz[i] = 8.0 * ts_xxxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yy[i] = 8.0 * ts_xxxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_yz[i] = 8.0 * ts_xxxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_zz[i] = 8.0 * ts_xxxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 18-24 components of targeted buffer : HD

    auto gs_x_xxxyy_xx = pbuffer.data(idx_g_hd + 18);

    auto gs_x_xxxyy_xy = pbuffer.data(idx_g_hd + 19);

    auto gs_x_xxxyy_xz = pbuffer.data(idx_g_hd + 20);

    auto gs_x_xxxyy_yy = pbuffer.data(idx_g_hd + 21);

    auto gs_x_xxxyy_yz = pbuffer.data(idx_g_hd + 22);

    auto gs_x_xxxyy_zz = pbuffer.data(idx_g_hd + 23);

    #pragma omp simd aligned(gc_x, gs_x_xxxyy_xx, gs_x_xxxyy_xy, gs_x_xxxyy_xz, gs_x_xxxyy_yy, gs_x_xxxyy_yz, gs_x_xxxyy_zz, ts_xxxyy_x, ts_xxxyy_xx, ts_xxxyy_xy, ts_xxxyy_xz, ts_xxxyy_y, ts_xxxyy_yy, ts_xxxyy_yz, ts_xxxyy_z, ts_xxxyy_zz, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyy_xx[i] = 6.0 * ts_xxyy_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xx[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xy[i] = 6.0 * ts_xxyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_xz[i] = 6.0 * ts_xxyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yy[i] = 6.0 * ts_xxyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yy[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_yz[i] = 6.0 * ts_xxyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yz[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_zz[i] = 6.0 * ts_xxyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 24-30 components of targeted buffer : HD

    auto gs_x_xxxyz_xx = pbuffer.data(idx_g_hd + 24);

    auto gs_x_xxxyz_xy = pbuffer.data(idx_g_hd + 25);

    auto gs_x_xxxyz_xz = pbuffer.data(idx_g_hd + 26);

    auto gs_x_xxxyz_yy = pbuffer.data(idx_g_hd + 27);

    auto gs_x_xxxyz_yz = pbuffer.data(idx_g_hd + 28);

    auto gs_x_xxxyz_zz = pbuffer.data(idx_g_hd + 29);

    #pragma omp simd aligned(gc_x, gs_x_xxxyz_xx, gs_x_xxxyz_xy, gs_x_xxxyz_xz, gs_x_xxxyz_yy, gs_x_xxxyz_yz, gs_x_xxxyz_zz, ts_xxxyz_x, ts_xxxyz_xx, ts_xxxyz_xy, ts_xxxyz_xz, ts_xxxyz_y, ts_xxxyz_yy, ts_xxxyz_yz, ts_xxxyz_z, ts_xxxyz_zz, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxyz_xx[i] = 6.0 * ts_xxyz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xy[i] = 6.0 * ts_xxyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_xz[i] = 6.0 * ts_xxyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yy[i] = 6.0 * ts_xxyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_yz[i] = 6.0 * ts_xxyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_zz[i] = 6.0 * ts_xxyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-36 components of targeted buffer : HD

    auto gs_x_xxxzz_xx = pbuffer.data(idx_g_hd + 30);

    auto gs_x_xxxzz_xy = pbuffer.data(idx_g_hd + 31);

    auto gs_x_xxxzz_xz = pbuffer.data(idx_g_hd + 32);

    auto gs_x_xxxzz_yy = pbuffer.data(idx_g_hd + 33);

    auto gs_x_xxxzz_yz = pbuffer.data(idx_g_hd + 34);

    auto gs_x_xxxzz_zz = pbuffer.data(idx_g_hd + 35);

    #pragma omp simd aligned(gc_x, gs_x_xxxzz_xx, gs_x_xxxzz_xy, gs_x_xxxzz_xz, gs_x_xxxzz_yy, gs_x_xxxzz_yz, gs_x_xxxzz_zz, ts_xxxzz_x, ts_xxxzz_xx, ts_xxxzz_xy, ts_xxxzz_xz, ts_xxxzz_y, ts_xxxzz_yy, ts_xxxzz_yz, ts_xxxzz_z, ts_xxxzz_zz, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxzz_xx[i] = 6.0 * ts_xxzz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xy[i] = 6.0 * ts_xxzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_xz[i] = 6.0 * ts_xxzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yy[i] = 6.0 * ts_xxzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_yz[i] = 6.0 * ts_xxzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_zz[i] = 6.0 * ts_xxzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 36-42 components of targeted buffer : HD

    auto gs_x_xxyyy_xx = pbuffer.data(idx_g_hd + 36);

    auto gs_x_xxyyy_xy = pbuffer.data(idx_g_hd + 37);

    auto gs_x_xxyyy_xz = pbuffer.data(idx_g_hd + 38);

    auto gs_x_xxyyy_yy = pbuffer.data(idx_g_hd + 39);

    auto gs_x_xxyyy_yz = pbuffer.data(idx_g_hd + 40);

    auto gs_x_xxyyy_zz = pbuffer.data(idx_g_hd + 41);

    #pragma omp simd aligned(gc_x, gs_x_xxyyy_xx, gs_x_xxyyy_xy, gs_x_xxyyy_xz, gs_x_xxyyy_yy, gs_x_xxyyy_yz, gs_x_xxyyy_zz, ts_xxyyy_x, ts_xxyyy_xx, ts_xxyyy_xy, ts_xxyyy_xz, ts_xxyyy_y, ts_xxyyy_yy, ts_xxyyy_yz, ts_xxyyy_z, ts_xxyyy_zz, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyy_xx[i] = 4.0 * ts_xyyy_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xx[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xy[i] = 4.0 * ts_xyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_xz[i] = 4.0 * ts_xyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yy[i] = 4.0 * ts_xyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yy[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_yz[i] = 4.0 * ts_xyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yz[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_zz[i] = 4.0 * ts_xyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 42-48 components of targeted buffer : HD

    auto gs_x_xxyyz_xx = pbuffer.data(idx_g_hd + 42);

    auto gs_x_xxyyz_xy = pbuffer.data(idx_g_hd + 43);

    auto gs_x_xxyyz_xz = pbuffer.data(idx_g_hd + 44);

    auto gs_x_xxyyz_yy = pbuffer.data(idx_g_hd + 45);

    auto gs_x_xxyyz_yz = pbuffer.data(idx_g_hd + 46);

    auto gs_x_xxyyz_zz = pbuffer.data(idx_g_hd + 47);

    #pragma omp simd aligned(gc_x, gs_x_xxyyz_xx, gs_x_xxyyz_xy, gs_x_xxyyz_xz, gs_x_xxyyz_yy, gs_x_xxyyz_yz, gs_x_xxyyz_zz, ts_xxyyz_x, ts_xxyyz_xx, ts_xxyyz_xy, ts_xxyyz_xz, ts_xxyyz_y, ts_xxyyz_yy, ts_xxyyz_yz, ts_xxyyz_z, ts_xxyyz_zz, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyyz_xx[i] = 4.0 * ts_xyyz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xy[i] = 4.0 * ts_xyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_xz[i] = 4.0 * ts_xyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yy[i] = 4.0 * ts_xyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_yz[i] = 4.0 * ts_xyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_zz[i] = 4.0 * ts_xyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 48-54 components of targeted buffer : HD

    auto gs_x_xxyzz_xx = pbuffer.data(idx_g_hd + 48);

    auto gs_x_xxyzz_xy = pbuffer.data(idx_g_hd + 49);

    auto gs_x_xxyzz_xz = pbuffer.data(idx_g_hd + 50);

    auto gs_x_xxyzz_yy = pbuffer.data(idx_g_hd + 51);

    auto gs_x_xxyzz_yz = pbuffer.data(idx_g_hd + 52);

    auto gs_x_xxyzz_zz = pbuffer.data(idx_g_hd + 53);

    #pragma omp simd aligned(gc_x, gs_x_xxyzz_xx, gs_x_xxyzz_xy, gs_x_xxyzz_xz, gs_x_xxyzz_yy, gs_x_xxyzz_yz, gs_x_xxyzz_zz, ts_xxyzz_x, ts_xxyzz_xx, ts_xxyzz_xy, ts_xxyzz_xz, ts_xxyzz_y, ts_xxyzz_yy, ts_xxyzz_yz, ts_xxyzz_z, ts_xxyzz_zz, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyzz_xx[i] = 4.0 * ts_xyzz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xy[i] = 4.0 * ts_xyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_xz[i] = 4.0 * ts_xyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yy[i] = 4.0 * ts_xyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_yz[i] = 4.0 * ts_xyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_zz[i] = 4.0 * ts_xyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 54-60 components of targeted buffer : HD

    auto gs_x_xxzzz_xx = pbuffer.data(idx_g_hd + 54);

    auto gs_x_xxzzz_xy = pbuffer.data(idx_g_hd + 55);

    auto gs_x_xxzzz_xz = pbuffer.data(idx_g_hd + 56);

    auto gs_x_xxzzz_yy = pbuffer.data(idx_g_hd + 57);

    auto gs_x_xxzzz_yz = pbuffer.data(idx_g_hd + 58);

    auto gs_x_xxzzz_zz = pbuffer.data(idx_g_hd + 59);

    #pragma omp simd aligned(gc_x, gs_x_xxzzz_xx, gs_x_xxzzz_xy, gs_x_xxzzz_xz, gs_x_xxzzz_yy, gs_x_xxzzz_yz, gs_x_xxzzz_zz, ts_xxzzz_x, ts_xxzzz_xx, ts_xxzzz_xy, ts_xxzzz_xz, ts_xxzzz_y, ts_xxzzz_yy, ts_xxzzz_yz, ts_xxzzz_z, ts_xxzzz_zz, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzzz_xx[i] = 4.0 * ts_xzzz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xy[i] = 4.0 * ts_xzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_xz[i] = 4.0 * ts_xzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yy[i] = 4.0 * ts_xzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_yz[i] = 4.0 * ts_xzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_zz[i] = 4.0 * ts_xzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-66 components of targeted buffer : HD

    auto gs_x_xyyyy_xx = pbuffer.data(idx_g_hd + 60);

    auto gs_x_xyyyy_xy = pbuffer.data(idx_g_hd + 61);

    auto gs_x_xyyyy_xz = pbuffer.data(idx_g_hd + 62);

    auto gs_x_xyyyy_yy = pbuffer.data(idx_g_hd + 63);

    auto gs_x_xyyyy_yz = pbuffer.data(idx_g_hd + 64);

    auto gs_x_xyyyy_zz = pbuffer.data(idx_g_hd + 65);

    #pragma omp simd aligned(gc_x, gs_x_xyyyy_xx, gs_x_xyyyy_xy, gs_x_xyyyy_xz, gs_x_xyyyy_yy, gs_x_xyyyy_yz, gs_x_xyyyy_zz, ts_xyyyy_x, ts_xyyyy_xx, ts_xyyyy_xy, ts_xyyyy_xz, ts_xyyyy_y, ts_xyyyy_yy, ts_xyyyy_yz, ts_xyyyy_z, ts_xyyyy_zz, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyy_xx[i] = 2.0 * ts_yyyy_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xx[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xy[i] = 2.0 * ts_yyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_xz[i] = 2.0 * ts_yyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yy[i] = 2.0 * ts_yyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yy[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_yz[i] = 2.0 * ts_yyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yz[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_zz[i] = 2.0 * ts_yyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 66-72 components of targeted buffer : HD

    auto gs_x_xyyyz_xx = pbuffer.data(idx_g_hd + 66);

    auto gs_x_xyyyz_xy = pbuffer.data(idx_g_hd + 67);

    auto gs_x_xyyyz_xz = pbuffer.data(idx_g_hd + 68);

    auto gs_x_xyyyz_yy = pbuffer.data(idx_g_hd + 69);

    auto gs_x_xyyyz_yz = pbuffer.data(idx_g_hd + 70);

    auto gs_x_xyyyz_zz = pbuffer.data(idx_g_hd + 71);

    #pragma omp simd aligned(gc_x, gs_x_xyyyz_xx, gs_x_xyyyz_xy, gs_x_xyyyz_xz, gs_x_xyyyz_yy, gs_x_xyyyz_yz, gs_x_xyyyz_zz, ts_xyyyz_x, ts_xyyyz_xx, ts_xyyyz_xy, ts_xyyyz_xz, ts_xyyyz_y, ts_xyyyz_yy, ts_xyyyz_yz, ts_xyyyz_z, ts_xyyyz_zz, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyyz_xx[i] = 2.0 * ts_yyyz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xx[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xy[i] = 2.0 * ts_yyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_xz[i] = 2.0 * ts_yyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yy[i] = 2.0 * ts_yyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yy[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_yz[i] = 2.0 * ts_yyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yz[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_zz[i] = 2.0 * ts_yyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 72-78 components of targeted buffer : HD

    auto gs_x_xyyzz_xx = pbuffer.data(idx_g_hd + 72);

    auto gs_x_xyyzz_xy = pbuffer.data(idx_g_hd + 73);

    auto gs_x_xyyzz_xz = pbuffer.data(idx_g_hd + 74);

    auto gs_x_xyyzz_yy = pbuffer.data(idx_g_hd + 75);

    auto gs_x_xyyzz_yz = pbuffer.data(idx_g_hd + 76);

    auto gs_x_xyyzz_zz = pbuffer.data(idx_g_hd + 77);

    #pragma omp simd aligned(gc_x, gs_x_xyyzz_xx, gs_x_xyyzz_xy, gs_x_xyyzz_xz, gs_x_xyyzz_yy, gs_x_xyyzz_yz, gs_x_xyyzz_zz, ts_xyyzz_x, ts_xyyzz_xx, ts_xyyzz_xy, ts_xyyzz_xz, ts_xyyzz_y, ts_xyyzz_yy, ts_xyyzz_yz, ts_xyyzz_z, ts_xyyzz_zz, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyzz_xx[i] = 2.0 * ts_yyzz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xy[i] = 2.0 * ts_yyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_xz[i] = 2.0 * ts_yyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yy[i] = 2.0 * ts_yyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_yz[i] = 2.0 * ts_yyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_zz[i] = 2.0 * ts_yyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 78-84 components of targeted buffer : HD

    auto gs_x_xyzzz_xx = pbuffer.data(idx_g_hd + 78);

    auto gs_x_xyzzz_xy = pbuffer.data(idx_g_hd + 79);

    auto gs_x_xyzzz_xz = pbuffer.data(idx_g_hd + 80);

    auto gs_x_xyzzz_yy = pbuffer.data(idx_g_hd + 81);

    auto gs_x_xyzzz_yz = pbuffer.data(idx_g_hd + 82);

    auto gs_x_xyzzz_zz = pbuffer.data(idx_g_hd + 83);

    #pragma omp simd aligned(gc_x, gs_x_xyzzz_xx, gs_x_xyzzz_xy, gs_x_xyzzz_xz, gs_x_xyzzz_yy, gs_x_xyzzz_yz, gs_x_xyzzz_zz, ts_xyzzz_x, ts_xyzzz_xx, ts_xyzzz_xy, ts_xyzzz_xz, ts_xyzzz_y, ts_xyzzz_yy, ts_xyzzz_yz, ts_xyzzz_z, ts_xyzzz_zz, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzzz_xx[i] = 2.0 * ts_yzzz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xy[i] = 2.0 * ts_yzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_xz[i] = 2.0 * ts_yzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yy[i] = 2.0 * ts_yzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_yz[i] = 2.0 * ts_yzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_zz[i] = 2.0 * ts_yzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 84-90 components of targeted buffer : HD

    auto gs_x_xzzzz_xx = pbuffer.data(idx_g_hd + 84);

    auto gs_x_xzzzz_xy = pbuffer.data(idx_g_hd + 85);

    auto gs_x_xzzzz_xz = pbuffer.data(idx_g_hd + 86);

    auto gs_x_xzzzz_yy = pbuffer.data(idx_g_hd + 87);

    auto gs_x_xzzzz_yz = pbuffer.data(idx_g_hd + 88);

    auto gs_x_xzzzz_zz = pbuffer.data(idx_g_hd + 89);

    #pragma omp simd aligned(gc_x, gs_x_xzzzz_xx, gs_x_xzzzz_xy, gs_x_xzzzz_xz, gs_x_xzzzz_yy, gs_x_xzzzz_yz, gs_x_xzzzz_zz, ts_xzzzz_x, ts_xzzzz_xx, ts_xzzzz_xy, ts_xzzzz_xz, ts_xzzzz_y, ts_xzzzz_yy, ts_xzzzz_yz, ts_xzzzz_z, ts_xzzzz_zz, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzzz_xx[i] = 2.0 * ts_zzzz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xy[i] = 2.0 * ts_zzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_xz[i] = 2.0 * ts_zzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yy[i] = 2.0 * ts_zzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_yz[i] = 2.0 * ts_zzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_zz[i] = 2.0 * ts_zzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 90-96 components of targeted buffer : HD

    auto gs_x_yyyyy_xx = pbuffer.data(idx_g_hd + 90);

    auto gs_x_yyyyy_xy = pbuffer.data(idx_g_hd + 91);

    auto gs_x_yyyyy_xz = pbuffer.data(idx_g_hd + 92);

    auto gs_x_yyyyy_yy = pbuffer.data(idx_g_hd + 93);

    auto gs_x_yyyyy_yz = pbuffer.data(idx_g_hd + 94);

    auto gs_x_yyyyy_zz = pbuffer.data(idx_g_hd + 95);

    #pragma omp simd aligned(gc_x, gs_x_yyyyy_xx, gs_x_yyyyy_xy, gs_x_yyyyy_xz, gs_x_yyyyy_yy, gs_x_yyyyy_yz, gs_x_yyyyy_zz, ts_yyyyy_x, ts_yyyyy_xx, ts_yyyyy_xy, ts_yyyyy_xz, ts_yyyyy_y, ts_yyyyy_yy, ts_yyyyy_yz, ts_yyyyy_z, ts_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyy_xx[i] = 4.0 * ts_yyyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xx[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xy[i] = 2.0 * ts_yyyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_xz[i] = 2.0 * ts_yyyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yy[i] = 2.0 * ts_yyyyy_yy[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_yz[i] = 2.0 * ts_yyyyy_yz[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_zz[i] = 2.0 * ts_yyyyy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 96-102 components of targeted buffer : HD

    auto gs_x_yyyyz_xx = pbuffer.data(idx_g_hd + 96);

    auto gs_x_yyyyz_xy = pbuffer.data(idx_g_hd + 97);

    auto gs_x_yyyyz_xz = pbuffer.data(idx_g_hd + 98);

    auto gs_x_yyyyz_yy = pbuffer.data(idx_g_hd + 99);

    auto gs_x_yyyyz_yz = pbuffer.data(idx_g_hd + 100);

    auto gs_x_yyyyz_zz = pbuffer.data(idx_g_hd + 101);

    #pragma omp simd aligned(gc_x, gs_x_yyyyz_xx, gs_x_yyyyz_xy, gs_x_yyyyz_xz, gs_x_yyyyz_yy, gs_x_yyyyz_yz, gs_x_yyyyz_zz, ts_yyyyz_x, ts_yyyyz_xx, ts_yyyyz_xy, ts_yyyyz_xz, ts_yyyyz_y, ts_yyyyz_yy, ts_yyyyz_yz, ts_yyyyz_z, ts_yyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyyz_xx[i] = 4.0 * ts_yyyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xx[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xy[i] = 2.0 * ts_yyyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_xz[i] = 2.0 * ts_yyyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yy[i] = 2.0 * ts_yyyyz_yy[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_yz[i] = 2.0 * ts_yyyyz_yz[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_zz[i] = 2.0 * ts_yyyyz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 102-108 components of targeted buffer : HD

    auto gs_x_yyyzz_xx = pbuffer.data(idx_g_hd + 102);

    auto gs_x_yyyzz_xy = pbuffer.data(idx_g_hd + 103);

    auto gs_x_yyyzz_xz = pbuffer.data(idx_g_hd + 104);

    auto gs_x_yyyzz_yy = pbuffer.data(idx_g_hd + 105);

    auto gs_x_yyyzz_yz = pbuffer.data(idx_g_hd + 106);

    auto gs_x_yyyzz_zz = pbuffer.data(idx_g_hd + 107);

    #pragma omp simd aligned(gc_x, gs_x_yyyzz_xx, gs_x_yyyzz_xy, gs_x_yyyzz_xz, gs_x_yyyzz_yy, gs_x_yyyzz_yz, gs_x_yyyzz_zz, ts_yyyzz_x, ts_yyyzz_xx, ts_yyyzz_xy, ts_yyyzz_xz, ts_yyyzz_y, ts_yyyzz_yy, ts_yyyzz_yz, ts_yyyzz_z, ts_yyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyzz_xx[i] = 4.0 * ts_yyyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xx[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xy[i] = 2.0 * ts_yyyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_xz[i] = 2.0 * ts_yyyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yy[i] = 2.0 * ts_yyyzz_yy[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_yz[i] = 2.0 * ts_yyyzz_yz[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_zz[i] = 2.0 * ts_yyyzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 108-114 components of targeted buffer : HD

    auto gs_x_yyzzz_xx = pbuffer.data(idx_g_hd + 108);

    auto gs_x_yyzzz_xy = pbuffer.data(idx_g_hd + 109);

    auto gs_x_yyzzz_xz = pbuffer.data(idx_g_hd + 110);

    auto gs_x_yyzzz_yy = pbuffer.data(idx_g_hd + 111);

    auto gs_x_yyzzz_yz = pbuffer.data(idx_g_hd + 112);

    auto gs_x_yyzzz_zz = pbuffer.data(idx_g_hd + 113);

    #pragma omp simd aligned(gc_x, gs_x_yyzzz_xx, gs_x_yyzzz_xy, gs_x_yyzzz_xz, gs_x_yyzzz_yy, gs_x_yyzzz_yz, gs_x_yyzzz_zz, ts_yyzzz_x, ts_yyzzz_xx, ts_yyzzz_xy, ts_yyzzz_xz, ts_yyzzz_y, ts_yyzzz_yy, ts_yyzzz_yz, ts_yyzzz_z, ts_yyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzzz_xx[i] = 4.0 * ts_yyzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xx[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xy[i] = 2.0 * ts_yyzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_xz[i] = 2.0 * ts_yyzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yy[i] = 2.0 * ts_yyzzz_yy[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_yz[i] = 2.0 * ts_yyzzz_yz[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_zz[i] = 2.0 * ts_yyzzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 114-120 components of targeted buffer : HD

    auto gs_x_yzzzz_xx = pbuffer.data(idx_g_hd + 114);

    auto gs_x_yzzzz_xy = pbuffer.data(idx_g_hd + 115);

    auto gs_x_yzzzz_xz = pbuffer.data(idx_g_hd + 116);

    auto gs_x_yzzzz_yy = pbuffer.data(idx_g_hd + 117);

    auto gs_x_yzzzz_yz = pbuffer.data(idx_g_hd + 118);

    auto gs_x_yzzzz_zz = pbuffer.data(idx_g_hd + 119);

    #pragma omp simd aligned(gc_x, gs_x_yzzzz_xx, gs_x_yzzzz_xy, gs_x_yzzzz_xz, gs_x_yzzzz_yy, gs_x_yzzzz_yz, gs_x_yzzzz_zz, ts_yzzzz_x, ts_yzzzz_xx, ts_yzzzz_xy, ts_yzzzz_xz, ts_yzzzz_y, ts_yzzzz_yy, ts_yzzzz_yz, ts_yzzzz_z, ts_yzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzzz_xx[i] = 4.0 * ts_yzzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xx[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xy[i] = 2.0 * ts_yzzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_xz[i] = 2.0 * ts_yzzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yy[i] = 2.0 * ts_yzzzz_yy[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_yz[i] = 2.0 * ts_yzzzz_yz[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_zz[i] = 2.0 * ts_yzzzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 120-126 components of targeted buffer : HD

    auto gs_x_zzzzz_xx = pbuffer.data(idx_g_hd + 120);

    auto gs_x_zzzzz_xy = pbuffer.data(idx_g_hd + 121);

    auto gs_x_zzzzz_xz = pbuffer.data(idx_g_hd + 122);

    auto gs_x_zzzzz_yy = pbuffer.data(idx_g_hd + 123);

    auto gs_x_zzzzz_yz = pbuffer.data(idx_g_hd + 124);

    auto gs_x_zzzzz_zz = pbuffer.data(idx_g_hd + 125);

    #pragma omp simd aligned(gc_x, gs_x_zzzzz_xx, gs_x_zzzzz_xy, gs_x_zzzzz_xz, gs_x_zzzzz_yy, gs_x_zzzzz_yz, gs_x_zzzzz_zz, ts_zzzzz_x, ts_zzzzz_xx, ts_zzzzz_xy, ts_zzzzz_xz, ts_zzzzz_y, ts_zzzzz_yy, ts_zzzzz_yz, ts_zzzzz_z, ts_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzzz_xx[i] = 4.0 * ts_zzzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xx[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xy[i] = 2.0 * ts_zzzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_xz[i] = 2.0 * ts_zzzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yy[i] = 2.0 * ts_zzzzz_yy[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_yz[i] = 2.0 * ts_zzzzz_yz[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_zz[i] = 2.0 * ts_zzzzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 126-132 components of targeted buffer : HD

    auto gs_y_xxxxx_xx = pbuffer.data(idx_g_hd + 126);

    auto gs_y_xxxxx_xy = pbuffer.data(idx_g_hd + 127);

    auto gs_y_xxxxx_xz = pbuffer.data(idx_g_hd + 128);

    auto gs_y_xxxxx_yy = pbuffer.data(idx_g_hd + 129);

    auto gs_y_xxxxx_yz = pbuffer.data(idx_g_hd + 130);

    auto gs_y_xxxxx_zz = pbuffer.data(idx_g_hd + 131);

    #pragma omp simd aligned(gc_y, gs_y_xxxxx_xx, gs_y_xxxxx_xy, gs_y_xxxxx_xz, gs_y_xxxxx_yy, gs_y_xxxxx_yz, gs_y_xxxxx_zz, ts_xxxxx_x, ts_xxxxx_xx, ts_xxxxx_xy, ts_xxxxx_xz, ts_xxxxx_y, ts_xxxxx_yy, ts_xxxxx_yz, ts_xxxxx_z, ts_xxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxx_xx[i] = 2.0 * ts_xxxxx_xx[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xy[i] = 2.0 * ts_xxxxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_xz[i] = 2.0 * ts_xxxxx_xz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yy[i] = 4.0 * ts_xxxxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yy[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_yz[i] = 2.0 * ts_xxxxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yz[i] * gc_y[i] * tce_0;

        gs_y_xxxxx_zz[i] = 2.0 * ts_xxxxx_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 132-138 components of targeted buffer : HD

    auto gs_y_xxxxy_xx = pbuffer.data(idx_g_hd + 132);

    auto gs_y_xxxxy_xy = pbuffer.data(idx_g_hd + 133);

    auto gs_y_xxxxy_xz = pbuffer.data(idx_g_hd + 134);

    auto gs_y_xxxxy_yy = pbuffer.data(idx_g_hd + 135);

    auto gs_y_xxxxy_yz = pbuffer.data(idx_g_hd + 136);

    auto gs_y_xxxxy_zz = pbuffer.data(idx_g_hd + 137);

    #pragma omp simd aligned(gc_y, gs_y_xxxxy_xx, gs_y_xxxxy_xy, gs_y_xxxxy_xz, gs_y_xxxxy_yy, gs_y_xxxxy_yz, gs_y_xxxxy_zz, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_zz, ts_xxxxy_x, ts_xxxxy_xx, ts_xxxxy_xy, ts_xxxxy_xz, ts_xxxxy_y, ts_xxxxy_yy, ts_xxxxy_yz, ts_xxxxy_z, ts_xxxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxy_xx[i] = 2.0 * ts_xxxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xx[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xy[i] = 2.0 * ts_xxxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_xz[i] = 2.0 * ts_xxxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yy[i] = 2.0 * ts_xxxx_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yy[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_yz[i] = 2.0 * ts_xxxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yz[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_zz[i] = 2.0 * ts_xxxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 138-144 components of targeted buffer : HD

    auto gs_y_xxxxz_xx = pbuffer.data(idx_g_hd + 138);

    auto gs_y_xxxxz_xy = pbuffer.data(idx_g_hd + 139);

    auto gs_y_xxxxz_xz = pbuffer.data(idx_g_hd + 140);

    auto gs_y_xxxxz_yy = pbuffer.data(idx_g_hd + 141);

    auto gs_y_xxxxz_yz = pbuffer.data(idx_g_hd + 142);

    auto gs_y_xxxxz_zz = pbuffer.data(idx_g_hd + 143);

    #pragma omp simd aligned(gc_y, gs_y_xxxxz_xx, gs_y_xxxxz_xy, gs_y_xxxxz_xz, gs_y_xxxxz_yy, gs_y_xxxxz_yz, gs_y_xxxxz_zz, ts_xxxxz_x, ts_xxxxz_xx, ts_xxxxz_xy, ts_xxxxz_xz, ts_xxxxz_y, ts_xxxxz_yy, ts_xxxxz_yz, ts_xxxxz_z, ts_xxxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxxz_xx[i] = 2.0 * ts_xxxxz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xy[i] = 2.0 * ts_xxxxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_xz[i] = 2.0 * ts_xxxxz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yy[i] = 4.0 * ts_xxxxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_yz[i] = 2.0 * ts_xxxxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_zz[i] = 2.0 * ts_xxxxz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 144-150 components of targeted buffer : HD

    auto gs_y_xxxyy_xx = pbuffer.data(idx_g_hd + 144);

    auto gs_y_xxxyy_xy = pbuffer.data(idx_g_hd + 145);

    auto gs_y_xxxyy_xz = pbuffer.data(idx_g_hd + 146);

    auto gs_y_xxxyy_yy = pbuffer.data(idx_g_hd + 147);

    auto gs_y_xxxyy_yz = pbuffer.data(idx_g_hd + 148);

    auto gs_y_xxxyy_zz = pbuffer.data(idx_g_hd + 149);

    #pragma omp simd aligned(gc_y, gs_y_xxxyy_xx, gs_y_xxxyy_xy, gs_y_xxxyy_xz, gs_y_xxxyy_yy, gs_y_xxxyy_yz, gs_y_xxxyy_zz, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_zz, ts_xxxyy_x, ts_xxxyy_xx, ts_xxxyy_xy, ts_xxxyy_xz, ts_xxxyy_y, ts_xxxyy_yy, ts_xxxyy_yz, ts_xxxyy_z, ts_xxxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyy_xx[i] = 4.0 * ts_xxxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xx[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xy[i] = 4.0 * ts_xxxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_xz[i] = 4.0 * ts_xxxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yy[i] = 4.0 * ts_xxxy_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yy[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_yz[i] = 4.0 * ts_xxxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yz[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_zz[i] = 4.0 * ts_xxxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 150-156 components of targeted buffer : HD

    auto gs_y_xxxyz_xx = pbuffer.data(idx_g_hd + 150);

    auto gs_y_xxxyz_xy = pbuffer.data(idx_g_hd + 151);

    auto gs_y_xxxyz_xz = pbuffer.data(idx_g_hd + 152);

    auto gs_y_xxxyz_yy = pbuffer.data(idx_g_hd + 153);

    auto gs_y_xxxyz_yz = pbuffer.data(idx_g_hd + 154);

    auto gs_y_xxxyz_zz = pbuffer.data(idx_g_hd + 155);

    #pragma omp simd aligned(gc_y, gs_y_xxxyz_xx, gs_y_xxxyz_xy, gs_y_xxxyz_xz, gs_y_xxxyz_yy, gs_y_xxxyz_yz, gs_y_xxxyz_zz, ts_xxxyz_x, ts_xxxyz_xx, ts_xxxyz_xy, ts_xxxyz_xz, ts_xxxyz_y, ts_xxxyz_yy, ts_xxxyz_yz, ts_xxxyz_z, ts_xxxyz_zz, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxyz_xx[i] = 2.0 * ts_xxxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xy[i] = 2.0 * ts_xxxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_xz[i] = 2.0 * ts_xxxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yy[i] = 2.0 * ts_xxxz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_yz[i] = 2.0 * ts_xxxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_zz[i] = 2.0 * ts_xxxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 156-162 components of targeted buffer : HD

    auto gs_y_xxxzz_xx = pbuffer.data(idx_g_hd + 156);

    auto gs_y_xxxzz_xy = pbuffer.data(idx_g_hd + 157);

    auto gs_y_xxxzz_xz = pbuffer.data(idx_g_hd + 158);

    auto gs_y_xxxzz_yy = pbuffer.data(idx_g_hd + 159);

    auto gs_y_xxxzz_yz = pbuffer.data(idx_g_hd + 160);

    auto gs_y_xxxzz_zz = pbuffer.data(idx_g_hd + 161);

    #pragma omp simd aligned(gc_y, gs_y_xxxzz_xx, gs_y_xxxzz_xy, gs_y_xxxzz_xz, gs_y_xxxzz_yy, gs_y_xxxzz_yz, gs_y_xxxzz_zz, ts_xxxzz_x, ts_xxxzz_xx, ts_xxxzz_xy, ts_xxxzz_xz, ts_xxxzz_y, ts_xxxzz_yy, ts_xxxzz_yz, ts_xxxzz_z, ts_xxxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxzz_xx[i] = 2.0 * ts_xxxzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xy[i] = 2.0 * ts_xxxzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_xz[i] = 2.0 * ts_xxxzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yy[i] = 4.0 * ts_xxxzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_yz[i] = 2.0 * ts_xxxzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_zz[i] = 2.0 * ts_xxxzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 162-168 components of targeted buffer : HD

    auto gs_y_xxyyy_xx = pbuffer.data(idx_g_hd + 162);

    auto gs_y_xxyyy_xy = pbuffer.data(idx_g_hd + 163);

    auto gs_y_xxyyy_xz = pbuffer.data(idx_g_hd + 164);

    auto gs_y_xxyyy_yy = pbuffer.data(idx_g_hd + 165);

    auto gs_y_xxyyy_yz = pbuffer.data(idx_g_hd + 166);

    auto gs_y_xxyyy_zz = pbuffer.data(idx_g_hd + 167);

    #pragma omp simd aligned(gc_y, gs_y_xxyyy_xx, gs_y_xxyyy_xy, gs_y_xxyyy_xz, gs_y_xxyyy_yy, gs_y_xxyyy_yz, gs_y_xxyyy_zz, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_zz, ts_xxyyy_x, ts_xxyyy_xx, ts_xxyyy_xy, ts_xxyyy_xz, ts_xxyyy_y, ts_xxyyy_yy, ts_xxyyy_yz, ts_xxyyy_z, ts_xxyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyy_xx[i] = 6.0 * ts_xxyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xx[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xy[i] = 6.0 * ts_xxyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_xz[i] = 6.0 * ts_xxyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yy[i] = 6.0 * ts_xxyy_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yy[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_yz[i] = 6.0 * ts_xxyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yz[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_zz[i] = 6.0 * ts_xxyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 168-174 components of targeted buffer : HD

    auto gs_y_xxyyz_xx = pbuffer.data(idx_g_hd + 168);

    auto gs_y_xxyyz_xy = pbuffer.data(idx_g_hd + 169);

    auto gs_y_xxyyz_xz = pbuffer.data(idx_g_hd + 170);

    auto gs_y_xxyyz_yy = pbuffer.data(idx_g_hd + 171);

    auto gs_y_xxyyz_yz = pbuffer.data(idx_g_hd + 172);

    auto gs_y_xxyyz_zz = pbuffer.data(idx_g_hd + 173);

    #pragma omp simd aligned(gc_y, gs_y_xxyyz_xx, gs_y_xxyyz_xy, gs_y_xxyyz_xz, gs_y_xxyyz_yy, gs_y_xxyyz_yz, gs_y_xxyyz_zz, ts_xxyyz_x, ts_xxyyz_xx, ts_xxyyz_xy, ts_xxyyz_xz, ts_xxyyz_y, ts_xxyyz_yy, ts_xxyyz_yz, ts_xxyyz_z, ts_xxyyz_zz, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyyz_xx[i] = 4.0 * ts_xxyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xy[i] = 4.0 * ts_xxyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_xz[i] = 4.0 * ts_xxyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yy[i] = 4.0 * ts_xxyz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_yz[i] = 4.0 * ts_xxyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_zz[i] = 4.0 * ts_xxyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 174-180 components of targeted buffer : HD

    auto gs_y_xxyzz_xx = pbuffer.data(idx_g_hd + 174);

    auto gs_y_xxyzz_xy = pbuffer.data(idx_g_hd + 175);

    auto gs_y_xxyzz_xz = pbuffer.data(idx_g_hd + 176);

    auto gs_y_xxyzz_yy = pbuffer.data(idx_g_hd + 177);

    auto gs_y_xxyzz_yz = pbuffer.data(idx_g_hd + 178);

    auto gs_y_xxyzz_zz = pbuffer.data(idx_g_hd + 179);

    #pragma omp simd aligned(gc_y, gs_y_xxyzz_xx, gs_y_xxyzz_xy, gs_y_xxyzz_xz, gs_y_xxyzz_yy, gs_y_xxyzz_yz, gs_y_xxyzz_zz, ts_xxyzz_x, ts_xxyzz_xx, ts_xxyzz_xy, ts_xxyzz_xz, ts_xxyzz_y, ts_xxyzz_yy, ts_xxyzz_yz, ts_xxyzz_z, ts_xxyzz_zz, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyzz_xx[i] = 2.0 * ts_xxzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xy[i] = 2.0 * ts_xxzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_xz[i] = 2.0 * ts_xxzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yy[i] = 2.0 * ts_xxzz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_yz[i] = 2.0 * ts_xxzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_zz[i] = 2.0 * ts_xxzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 180-186 components of targeted buffer : HD

    auto gs_y_xxzzz_xx = pbuffer.data(idx_g_hd + 180);

    auto gs_y_xxzzz_xy = pbuffer.data(idx_g_hd + 181);

    auto gs_y_xxzzz_xz = pbuffer.data(idx_g_hd + 182);

    auto gs_y_xxzzz_yy = pbuffer.data(idx_g_hd + 183);

    auto gs_y_xxzzz_yz = pbuffer.data(idx_g_hd + 184);

    auto gs_y_xxzzz_zz = pbuffer.data(idx_g_hd + 185);

    #pragma omp simd aligned(gc_y, gs_y_xxzzz_xx, gs_y_xxzzz_xy, gs_y_xxzzz_xz, gs_y_xxzzz_yy, gs_y_xxzzz_yz, gs_y_xxzzz_zz, ts_xxzzz_x, ts_xxzzz_xx, ts_xxzzz_xy, ts_xxzzz_xz, ts_xxzzz_y, ts_xxzzz_yy, ts_xxzzz_yz, ts_xxzzz_z, ts_xxzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzzz_xx[i] = 2.0 * ts_xxzzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xy[i] = 2.0 * ts_xxzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_xz[i] = 2.0 * ts_xxzzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yy[i] = 4.0 * ts_xxzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_yz[i] = 2.0 * ts_xxzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_zz[i] = 2.0 * ts_xxzzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 186-192 components of targeted buffer : HD

    auto gs_y_xyyyy_xx = pbuffer.data(idx_g_hd + 186);

    auto gs_y_xyyyy_xy = pbuffer.data(idx_g_hd + 187);

    auto gs_y_xyyyy_xz = pbuffer.data(idx_g_hd + 188);

    auto gs_y_xyyyy_yy = pbuffer.data(idx_g_hd + 189);

    auto gs_y_xyyyy_yz = pbuffer.data(idx_g_hd + 190);

    auto gs_y_xyyyy_zz = pbuffer.data(idx_g_hd + 191);

    #pragma omp simd aligned(gc_y, gs_y_xyyyy_xx, gs_y_xyyyy_xy, gs_y_xyyyy_xz, gs_y_xyyyy_yy, gs_y_xyyyy_yz, gs_y_xyyyy_zz, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_zz, ts_xyyyy_x, ts_xyyyy_xx, ts_xyyyy_xy, ts_xyyyy_xz, ts_xyyyy_y, ts_xyyyy_yy, ts_xyyyy_yz, ts_xyyyy_z, ts_xyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyy_xx[i] = 8.0 * ts_xyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xx[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xy[i] = 8.0 * ts_xyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_xz[i] = 8.0 * ts_xyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yy[i] = 8.0 * ts_xyyy_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yy[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_yz[i] = 8.0 * ts_xyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yz[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_zz[i] = 8.0 * ts_xyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 192-198 components of targeted buffer : HD

    auto gs_y_xyyyz_xx = pbuffer.data(idx_g_hd + 192);

    auto gs_y_xyyyz_xy = pbuffer.data(idx_g_hd + 193);

    auto gs_y_xyyyz_xz = pbuffer.data(idx_g_hd + 194);

    auto gs_y_xyyyz_yy = pbuffer.data(idx_g_hd + 195);

    auto gs_y_xyyyz_yz = pbuffer.data(idx_g_hd + 196);

    auto gs_y_xyyyz_zz = pbuffer.data(idx_g_hd + 197);

    #pragma omp simd aligned(gc_y, gs_y_xyyyz_xx, gs_y_xyyyz_xy, gs_y_xyyyz_xz, gs_y_xyyyz_yy, gs_y_xyyyz_yz, gs_y_xyyyz_zz, ts_xyyyz_x, ts_xyyyz_xx, ts_xyyyz_xy, ts_xyyyz_xz, ts_xyyyz_y, ts_xyyyz_yy, ts_xyyyz_yz, ts_xyyyz_z, ts_xyyyz_zz, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyyz_xx[i] = 6.0 * ts_xyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xx[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xy[i] = 6.0 * ts_xyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_xz[i] = 6.0 * ts_xyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yy[i] = 6.0 * ts_xyyz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yy[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_yz[i] = 6.0 * ts_xyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yz[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_zz[i] = 6.0 * ts_xyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 198-204 components of targeted buffer : HD

    auto gs_y_xyyzz_xx = pbuffer.data(idx_g_hd + 198);

    auto gs_y_xyyzz_xy = pbuffer.data(idx_g_hd + 199);

    auto gs_y_xyyzz_xz = pbuffer.data(idx_g_hd + 200);

    auto gs_y_xyyzz_yy = pbuffer.data(idx_g_hd + 201);

    auto gs_y_xyyzz_yz = pbuffer.data(idx_g_hd + 202);

    auto gs_y_xyyzz_zz = pbuffer.data(idx_g_hd + 203);

    #pragma omp simd aligned(gc_y, gs_y_xyyzz_xx, gs_y_xyyzz_xy, gs_y_xyyzz_xz, gs_y_xyyzz_yy, gs_y_xyyzz_yz, gs_y_xyyzz_zz, ts_xyyzz_x, ts_xyyzz_xx, ts_xyyzz_xy, ts_xyyzz_xz, ts_xyyzz_y, ts_xyyzz_yy, ts_xyyzz_yz, ts_xyyzz_z, ts_xyyzz_zz, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyzz_xx[i] = 4.0 * ts_xyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xy[i] = 4.0 * ts_xyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_xz[i] = 4.0 * ts_xyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yy[i] = 4.0 * ts_xyzz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_yz[i] = 4.0 * ts_xyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_zz[i] = 4.0 * ts_xyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 204-210 components of targeted buffer : HD

    auto gs_y_xyzzz_xx = pbuffer.data(idx_g_hd + 204);

    auto gs_y_xyzzz_xy = pbuffer.data(idx_g_hd + 205);

    auto gs_y_xyzzz_xz = pbuffer.data(idx_g_hd + 206);

    auto gs_y_xyzzz_yy = pbuffer.data(idx_g_hd + 207);

    auto gs_y_xyzzz_yz = pbuffer.data(idx_g_hd + 208);

    auto gs_y_xyzzz_zz = pbuffer.data(idx_g_hd + 209);

    #pragma omp simd aligned(gc_y, gs_y_xyzzz_xx, gs_y_xyzzz_xy, gs_y_xyzzz_xz, gs_y_xyzzz_yy, gs_y_xyzzz_yz, gs_y_xyzzz_zz, ts_xyzzz_x, ts_xyzzz_xx, ts_xyzzz_xy, ts_xyzzz_xz, ts_xyzzz_y, ts_xyzzz_yy, ts_xyzzz_yz, ts_xyzzz_z, ts_xyzzz_zz, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzzz_xx[i] = 2.0 * ts_xzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xy[i] = 2.0 * ts_xzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_xz[i] = 2.0 * ts_xzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yy[i] = 2.0 * ts_xzzz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_yz[i] = 2.0 * ts_xzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_zz[i] = 2.0 * ts_xzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 210-216 components of targeted buffer : HD

    auto gs_y_xzzzz_xx = pbuffer.data(idx_g_hd + 210);

    auto gs_y_xzzzz_xy = pbuffer.data(idx_g_hd + 211);

    auto gs_y_xzzzz_xz = pbuffer.data(idx_g_hd + 212);

    auto gs_y_xzzzz_yy = pbuffer.data(idx_g_hd + 213);

    auto gs_y_xzzzz_yz = pbuffer.data(idx_g_hd + 214);

    auto gs_y_xzzzz_zz = pbuffer.data(idx_g_hd + 215);

    #pragma omp simd aligned(gc_y, gs_y_xzzzz_xx, gs_y_xzzzz_xy, gs_y_xzzzz_xz, gs_y_xzzzz_yy, gs_y_xzzzz_yz, gs_y_xzzzz_zz, ts_xzzzz_x, ts_xzzzz_xx, ts_xzzzz_xy, ts_xzzzz_xz, ts_xzzzz_y, ts_xzzzz_yy, ts_xzzzz_yz, ts_xzzzz_z, ts_xzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzzz_xx[i] = 2.0 * ts_xzzzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xy[i] = 2.0 * ts_xzzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_xz[i] = 2.0 * ts_xzzzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yy[i] = 4.0 * ts_xzzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_yz[i] = 2.0 * ts_xzzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_zz[i] = 2.0 * ts_xzzzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 216-222 components of targeted buffer : HD

    auto gs_y_yyyyy_xx = pbuffer.data(idx_g_hd + 216);

    auto gs_y_yyyyy_xy = pbuffer.data(idx_g_hd + 217);

    auto gs_y_yyyyy_xz = pbuffer.data(idx_g_hd + 218);

    auto gs_y_yyyyy_yy = pbuffer.data(idx_g_hd + 219);

    auto gs_y_yyyyy_yz = pbuffer.data(idx_g_hd + 220);

    auto gs_y_yyyyy_zz = pbuffer.data(idx_g_hd + 221);

    #pragma omp simd aligned(gc_y, gs_y_yyyyy_xx, gs_y_yyyyy_xy, gs_y_yyyyy_xz, gs_y_yyyyy_yy, gs_y_yyyyy_yz, gs_y_yyyyy_zz, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_zz, ts_yyyyy_x, ts_yyyyy_xx, ts_yyyyy_xy, ts_yyyyy_xz, ts_yyyyy_y, ts_yyyyy_yy, ts_yyyyy_yz, ts_yyyyy_z, ts_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyy_xx[i] = 10.0 * ts_yyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xx[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xy[i] = 10.0 * ts_yyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_xz[i] = 10.0 * ts_yyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yy[i] = 10.0 * ts_yyyy_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yy[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_yz[i] = 10.0 * ts_yyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yz[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_zz[i] = 10.0 * ts_yyyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 222-228 components of targeted buffer : HD

    auto gs_y_yyyyz_xx = pbuffer.data(idx_g_hd + 222);

    auto gs_y_yyyyz_xy = pbuffer.data(idx_g_hd + 223);

    auto gs_y_yyyyz_xz = pbuffer.data(idx_g_hd + 224);

    auto gs_y_yyyyz_yy = pbuffer.data(idx_g_hd + 225);

    auto gs_y_yyyyz_yz = pbuffer.data(idx_g_hd + 226);

    auto gs_y_yyyyz_zz = pbuffer.data(idx_g_hd + 227);

    #pragma omp simd aligned(gc_y, gs_y_yyyyz_xx, gs_y_yyyyz_xy, gs_y_yyyyz_xz, gs_y_yyyyz_yy, gs_y_yyyyz_yz, gs_y_yyyyz_zz, ts_yyyyz_x, ts_yyyyz_xx, ts_yyyyz_xy, ts_yyyyz_xz, ts_yyyyz_y, ts_yyyyz_yy, ts_yyyyz_yz, ts_yyyyz_z, ts_yyyyz_zz, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyyz_xx[i] = 8.0 * ts_yyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xx[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xy[i] = 8.0 * ts_yyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_xz[i] = 8.0 * ts_yyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yy[i] = 8.0 * ts_yyyz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yy[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_yz[i] = 8.0 * ts_yyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yz[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_zz[i] = 8.0 * ts_yyyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 228-234 components of targeted buffer : HD

    auto gs_y_yyyzz_xx = pbuffer.data(idx_g_hd + 228);

    auto gs_y_yyyzz_xy = pbuffer.data(idx_g_hd + 229);

    auto gs_y_yyyzz_xz = pbuffer.data(idx_g_hd + 230);

    auto gs_y_yyyzz_yy = pbuffer.data(idx_g_hd + 231);

    auto gs_y_yyyzz_yz = pbuffer.data(idx_g_hd + 232);

    auto gs_y_yyyzz_zz = pbuffer.data(idx_g_hd + 233);

    #pragma omp simd aligned(gc_y, gs_y_yyyzz_xx, gs_y_yyyzz_xy, gs_y_yyyzz_xz, gs_y_yyyzz_yy, gs_y_yyyzz_yz, gs_y_yyyzz_zz, ts_yyyzz_x, ts_yyyzz_xx, ts_yyyzz_xy, ts_yyyzz_xz, ts_yyyzz_y, ts_yyyzz_yy, ts_yyyzz_yz, ts_yyyzz_z, ts_yyyzz_zz, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyzz_xx[i] = 6.0 * ts_yyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xx[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xy[i] = 6.0 * ts_yyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_xz[i] = 6.0 * ts_yyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yy[i] = 6.0 * ts_yyzz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yy[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_yz[i] = 6.0 * ts_yyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yz[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_zz[i] = 6.0 * ts_yyzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 234-240 components of targeted buffer : HD

    auto gs_y_yyzzz_xx = pbuffer.data(idx_g_hd + 234);

    auto gs_y_yyzzz_xy = pbuffer.data(idx_g_hd + 235);

    auto gs_y_yyzzz_xz = pbuffer.data(idx_g_hd + 236);

    auto gs_y_yyzzz_yy = pbuffer.data(idx_g_hd + 237);

    auto gs_y_yyzzz_yz = pbuffer.data(idx_g_hd + 238);

    auto gs_y_yyzzz_zz = pbuffer.data(idx_g_hd + 239);

    #pragma omp simd aligned(gc_y, gs_y_yyzzz_xx, gs_y_yyzzz_xy, gs_y_yyzzz_xz, gs_y_yyzzz_yy, gs_y_yyzzz_yz, gs_y_yyzzz_zz, ts_yyzzz_x, ts_yyzzz_xx, ts_yyzzz_xy, ts_yyzzz_xz, ts_yyzzz_y, ts_yyzzz_yy, ts_yyzzz_yz, ts_yyzzz_z, ts_yyzzz_zz, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzzz_xx[i] = 4.0 * ts_yzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xx[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xy[i] = 4.0 * ts_yzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_xz[i] = 4.0 * ts_yzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yy[i] = 4.0 * ts_yzzz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yy[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_yz[i] = 4.0 * ts_yzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yz[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_zz[i] = 4.0 * ts_yzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 240-246 components of targeted buffer : HD

    auto gs_y_yzzzz_xx = pbuffer.data(idx_g_hd + 240);

    auto gs_y_yzzzz_xy = pbuffer.data(idx_g_hd + 241);

    auto gs_y_yzzzz_xz = pbuffer.data(idx_g_hd + 242);

    auto gs_y_yzzzz_yy = pbuffer.data(idx_g_hd + 243);

    auto gs_y_yzzzz_yz = pbuffer.data(idx_g_hd + 244);

    auto gs_y_yzzzz_zz = pbuffer.data(idx_g_hd + 245);

    #pragma omp simd aligned(gc_y, gs_y_yzzzz_xx, gs_y_yzzzz_xy, gs_y_yzzzz_xz, gs_y_yzzzz_yy, gs_y_yzzzz_yz, gs_y_yzzzz_zz, ts_yzzzz_x, ts_yzzzz_xx, ts_yzzzz_xy, ts_yzzzz_xz, ts_yzzzz_y, ts_yzzzz_yy, ts_yzzzz_yz, ts_yzzzz_z, ts_yzzzz_zz, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzzz_xx[i] = 2.0 * ts_zzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xx[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xy[i] = 2.0 * ts_zzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_xz[i] = 2.0 * ts_zzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yy[i] = 2.0 * ts_zzzz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yy[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_yz[i] = 2.0 * ts_zzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yz[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_zz[i] = 2.0 * ts_zzzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 246-252 components of targeted buffer : HD

    auto gs_y_zzzzz_xx = pbuffer.data(idx_g_hd + 246);

    auto gs_y_zzzzz_xy = pbuffer.data(idx_g_hd + 247);

    auto gs_y_zzzzz_xz = pbuffer.data(idx_g_hd + 248);

    auto gs_y_zzzzz_yy = pbuffer.data(idx_g_hd + 249);

    auto gs_y_zzzzz_yz = pbuffer.data(idx_g_hd + 250);

    auto gs_y_zzzzz_zz = pbuffer.data(idx_g_hd + 251);

    #pragma omp simd aligned(gc_y, gs_y_zzzzz_xx, gs_y_zzzzz_xy, gs_y_zzzzz_xz, gs_y_zzzzz_yy, gs_y_zzzzz_yz, gs_y_zzzzz_zz, ts_zzzzz_x, ts_zzzzz_xx, ts_zzzzz_xy, ts_zzzzz_xz, ts_zzzzz_y, ts_zzzzz_yy, ts_zzzzz_yz, ts_zzzzz_z, ts_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzzz_xx[i] = 2.0 * ts_zzzzz_xx[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xy[i] = 2.0 * ts_zzzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_xz[i] = 2.0 * ts_zzzzz_xz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yy[i] = 4.0 * ts_zzzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yy[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_yz[i] = 2.0 * ts_zzzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yz[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_zz[i] = 2.0 * ts_zzzzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 252-258 components of targeted buffer : HD

    auto gs_z_xxxxx_xx = pbuffer.data(idx_g_hd + 252);

    auto gs_z_xxxxx_xy = pbuffer.data(idx_g_hd + 253);

    auto gs_z_xxxxx_xz = pbuffer.data(idx_g_hd + 254);

    auto gs_z_xxxxx_yy = pbuffer.data(idx_g_hd + 255);

    auto gs_z_xxxxx_yz = pbuffer.data(idx_g_hd + 256);

    auto gs_z_xxxxx_zz = pbuffer.data(idx_g_hd + 257);

    #pragma omp simd aligned(gc_z, gs_z_xxxxx_xx, gs_z_xxxxx_xy, gs_z_xxxxx_xz, gs_z_xxxxx_yy, gs_z_xxxxx_yz, gs_z_xxxxx_zz, ts_xxxxx_x, ts_xxxxx_xx, ts_xxxxx_xy, ts_xxxxx_xz, ts_xxxxx_y, ts_xxxxx_yy, ts_xxxxx_yz, ts_xxxxx_z, ts_xxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxx_xx[i] = 2.0 * ts_xxxxx_xx[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xy[i] = 2.0 * ts_xxxxx_xy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_xz[i] = 2.0 * ts_xxxxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_xz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yy[i] = 2.0 * ts_xxxxx_yy[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_yz[i] = 2.0 * ts_xxxxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_yz[i] * gc_z[i] * tce_0;

        gs_z_xxxxx_zz[i] = 4.0 * ts_xxxxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 258-264 components of targeted buffer : HD

    auto gs_z_xxxxy_xx = pbuffer.data(idx_g_hd + 258);

    auto gs_z_xxxxy_xy = pbuffer.data(idx_g_hd + 259);

    auto gs_z_xxxxy_xz = pbuffer.data(idx_g_hd + 260);

    auto gs_z_xxxxy_yy = pbuffer.data(idx_g_hd + 261);

    auto gs_z_xxxxy_yz = pbuffer.data(idx_g_hd + 262);

    auto gs_z_xxxxy_zz = pbuffer.data(idx_g_hd + 263);

    #pragma omp simd aligned(gc_z, gs_z_xxxxy_xx, gs_z_xxxxy_xy, gs_z_xxxxy_xz, gs_z_xxxxy_yy, gs_z_xxxxy_yz, gs_z_xxxxy_zz, ts_xxxxy_x, ts_xxxxy_xx, ts_xxxxy_xy, ts_xxxxy_xz, ts_xxxxy_y, ts_xxxxy_yy, ts_xxxxy_yz, ts_xxxxy_z, ts_xxxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxy_xx[i] = 2.0 * ts_xxxxy_xx[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xy[i] = 2.0 * ts_xxxxy_xy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_xz[i] = 2.0 * ts_xxxxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_xz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yy[i] = 2.0 * ts_xxxxy_yy[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_yz[i] = 2.0 * ts_xxxxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_yz[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_zz[i] = 4.0 * ts_xxxxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 264-270 components of targeted buffer : HD

    auto gs_z_xxxxz_xx = pbuffer.data(idx_g_hd + 264);

    auto gs_z_xxxxz_xy = pbuffer.data(idx_g_hd + 265);

    auto gs_z_xxxxz_xz = pbuffer.data(idx_g_hd + 266);

    auto gs_z_xxxxz_yy = pbuffer.data(idx_g_hd + 267);

    auto gs_z_xxxxz_yz = pbuffer.data(idx_g_hd + 268);

    auto gs_z_xxxxz_zz = pbuffer.data(idx_g_hd + 269);

    #pragma omp simd aligned(gc_z, gs_z_xxxxz_xx, gs_z_xxxxz_xy, gs_z_xxxxz_xz, gs_z_xxxxz_yy, gs_z_xxxxz_yz, gs_z_xxxxz_zz, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_zz, ts_xxxxz_x, ts_xxxxz_xx, ts_xxxxz_xy, ts_xxxxz_xz, ts_xxxxz_y, ts_xxxxz_yy, ts_xxxxz_yz, ts_xxxxz_z, ts_xxxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxxz_xx[i] = 2.0 * ts_xxxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xy[i] = 2.0 * ts_xxxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_xz[i] = 2.0 * ts_xxxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yy[i] = 2.0 * ts_xxxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_yz[i] = 2.0 * ts_xxxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_zz[i] = 2.0 * ts_xxxx_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 270-276 components of targeted buffer : HD

    auto gs_z_xxxyy_xx = pbuffer.data(idx_g_hd + 270);

    auto gs_z_xxxyy_xy = pbuffer.data(idx_g_hd + 271);

    auto gs_z_xxxyy_xz = pbuffer.data(idx_g_hd + 272);

    auto gs_z_xxxyy_yy = pbuffer.data(idx_g_hd + 273);

    auto gs_z_xxxyy_yz = pbuffer.data(idx_g_hd + 274);

    auto gs_z_xxxyy_zz = pbuffer.data(idx_g_hd + 275);

    #pragma omp simd aligned(gc_z, gs_z_xxxyy_xx, gs_z_xxxyy_xy, gs_z_xxxyy_xz, gs_z_xxxyy_yy, gs_z_xxxyy_yz, gs_z_xxxyy_zz, ts_xxxyy_x, ts_xxxyy_xx, ts_xxxyy_xy, ts_xxxyy_xz, ts_xxxyy_y, ts_xxxyy_yy, ts_xxxyy_yz, ts_xxxyy_z, ts_xxxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyy_xx[i] = 2.0 * ts_xxxyy_xx[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xy[i] = 2.0 * ts_xxxyy_xy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_xz[i] = 2.0 * ts_xxxyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_xz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yy[i] = 2.0 * ts_xxxyy_yy[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_yz[i] = 2.0 * ts_xxxyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_yz[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_zz[i] = 4.0 * ts_xxxyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 276-282 components of targeted buffer : HD

    auto gs_z_xxxyz_xx = pbuffer.data(idx_g_hd + 276);

    auto gs_z_xxxyz_xy = pbuffer.data(idx_g_hd + 277);

    auto gs_z_xxxyz_xz = pbuffer.data(idx_g_hd + 278);

    auto gs_z_xxxyz_yy = pbuffer.data(idx_g_hd + 279);

    auto gs_z_xxxyz_yz = pbuffer.data(idx_g_hd + 280);

    auto gs_z_xxxyz_zz = pbuffer.data(idx_g_hd + 281);

    #pragma omp simd aligned(gc_z, gs_z_xxxyz_xx, gs_z_xxxyz_xy, gs_z_xxxyz_xz, gs_z_xxxyz_yy, gs_z_xxxyz_yz, gs_z_xxxyz_zz, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_zz, ts_xxxyz_x, ts_xxxyz_xx, ts_xxxyz_xy, ts_xxxyz_xz, ts_xxxyz_y, ts_xxxyz_yy, ts_xxxyz_yz, ts_xxxyz_z, ts_xxxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxyz_xx[i] = 2.0 * ts_xxxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xy[i] = 2.0 * ts_xxxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_xz[i] = 2.0 * ts_xxxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yy[i] = 2.0 * ts_xxxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_yz[i] = 2.0 * ts_xxxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_zz[i] = 2.0 * ts_xxxy_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 282-288 components of targeted buffer : HD

    auto gs_z_xxxzz_xx = pbuffer.data(idx_g_hd + 282);

    auto gs_z_xxxzz_xy = pbuffer.data(idx_g_hd + 283);

    auto gs_z_xxxzz_xz = pbuffer.data(idx_g_hd + 284);

    auto gs_z_xxxzz_yy = pbuffer.data(idx_g_hd + 285);

    auto gs_z_xxxzz_yz = pbuffer.data(idx_g_hd + 286);

    auto gs_z_xxxzz_zz = pbuffer.data(idx_g_hd + 287);

    #pragma omp simd aligned(gc_z, gs_z_xxxzz_xx, gs_z_xxxzz_xy, gs_z_xxxzz_xz, gs_z_xxxzz_yy, gs_z_xxxzz_yz, gs_z_xxxzz_zz, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_zz, ts_xxxzz_x, ts_xxxzz_xx, ts_xxxzz_xy, ts_xxxzz_xz, ts_xxxzz_y, ts_xxxzz_yy, ts_xxxzz_yz, ts_xxxzz_z, ts_xxxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxzz_xx[i] = 4.0 * ts_xxxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xy[i] = 4.0 * ts_xxxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_xz[i] = 4.0 * ts_xxxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yy[i] = 4.0 * ts_xxxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_yz[i] = 4.0 * ts_xxxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_zz[i] = 4.0 * ts_xxxz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 288-294 components of targeted buffer : HD

    auto gs_z_xxyyy_xx = pbuffer.data(idx_g_hd + 288);

    auto gs_z_xxyyy_xy = pbuffer.data(idx_g_hd + 289);

    auto gs_z_xxyyy_xz = pbuffer.data(idx_g_hd + 290);

    auto gs_z_xxyyy_yy = pbuffer.data(idx_g_hd + 291);

    auto gs_z_xxyyy_yz = pbuffer.data(idx_g_hd + 292);

    auto gs_z_xxyyy_zz = pbuffer.data(idx_g_hd + 293);

    #pragma omp simd aligned(gc_z, gs_z_xxyyy_xx, gs_z_xxyyy_xy, gs_z_xxyyy_xz, gs_z_xxyyy_yy, gs_z_xxyyy_yz, gs_z_xxyyy_zz, ts_xxyyy_x, ts_xxyyy_xx, ts_xxyyy_xy, ts_xxyyy_xz, ts_xxyyy_y, ts_xxyyy_yy, ts_xxyyy_yz, ts_xxyyy_z, ts_xxyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyy_xx[i] = 2.0 * ts_xxyyy_xx[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xy[i] = 2.0 * ts_xxyyy_xy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_xz[i] = 2.0 * ts_xxyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_xz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yy[i] = 2.0 * ts_xxyyy_yy[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_yz[i] = 2.0 * ts_xxyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_yz[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_zz[i] = 4.0 * ts_xxyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 294-300 components of targeted buffer : HD

    auto gs_z_xxyyz_xx = pbuffer.data(idx_g_hd + 294);

    auto gs_z_xxyyz_xy = pbuffer.data(idx_g_hd + 295);

    auto gs_z_xxyyz_xz = pbuffer.data(idx_g_hd + 296);

    auto gs_z_xxyyz_yy = pbuffer.data(idx_g_hd + 297);

    auto gs_z_xxyyz_yz = pbuffer.data(idx_g_hd + 298);

    auto gs_z_xxyyz_zz = pbuffer.data(idx_g_hd + 299);

    #pragma omp simd aligned(gc_z, gs_z_xxyyz_xx, gs_z_xxyyz_xy, gs_z_xxyyz_xz, gs_z_xxyyz_yy, gs_z_xxyyz_yz, gs_z_xxyyz_zz, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_zz, ts_xxyyz_x, ts_xxyyz_xx, ts_xxyyz_xy, ts_xxyyz_xz, ts_xxyyz_y, ts_xxyyz_yy, ts_xxyyz_yz, ts_xxyyz_z, ts_xxyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyyz_xx[i] = 2.0 * ts_xxyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xy[i] = 2.0 * ts_xxyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_xz[i] = 2.0 * ts_xxyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yy[i] = 2.0 * ts_xxyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_yz[i] = 2.0 * ts_xxyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_zz[i] = 2.0 * ts_xxyy_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 300-306 components of targeted buffer : HD

    auto gs_z_xxyzz_xx = pbuffer.data(idx_g_hd + 300);

    auto gs_z_xxyzz_xy = pbuffer.data(idx_g_hd + 301);

    auto gs_z_xxyzz_xz = pbuffer.data(idx_g_hd + 302);

    auto gs_z_xxyzz_yy = pbuffer.data(idx_g_hd + 303);

    auto gs_z_xxyzz_yz = pbuffer.data(idx_g_hd + 304);

    auto gs_z_xxyzz_zz = pbuffer.data(idx_g_hd + 305);

    #pragma omp simd aligned(gc_z, gs_z_xxyzz_xx, gs_z_xxyzz_xy, gs_z_xxyzz_xz, gs_z_xxyzz_yy, gs_z_xxyzz_yz, gs_z_xxyzz_zz, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_zz, ts_xxyzz_x, ts_xxyzz_xx, ts_xxyzz_xy, ts_xxyzz_xz, ts_xxyzz_y, ts_xxyzz_yy, ts_xxyzz_yz, ts_xxyzz_z, ts_xxyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyzz_xx[i] = 4.0 * ts_xxyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xy[i] = 4.0 * ts_xxyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_xz[i] = 4.0 * ts_xxyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yy[i] = 4.0 * ts_xxyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_yz[i] = 4.0 * ts_xxyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_zz[i] = 4.0 * ts_xxyz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 306-312 components of targeted buffer : HD

    auto gs_z_xxzzz_xx = pbuffer.data(idx_g_hd + 306);

    auto gs_z_xxzzz_xy = pbuffer.data(idx_g_hd + 307);

    auto gs_z_xxzzz_xz = pbuffer.data(idx_g_hd + 308);

    auto gs_z_xxzzz_yy = pbuffer.data(idx_g_hd + 309);

    auto gs_z_xxzzz_yz = pbuffer.data(idx_g_hd + 310);

    auto gs_z_xxzzz_zz = pbuffer.data(idx_g_hd + 311);

    #pragma omp simd aligned(gc_z, gs_z_xxzzz_xx, gs_z_xxzzz_xy, gs_z_xxzzz_xz, gs_z_xxzzz_yy, gs_z_xxzzz_yz, gs_z_xxzzz_zz, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_zz, ts_xxzzz_x, ts_xxzzz_xx, ts_xxzzz_xy, ts_xxzzz_xz, ts_xxzzz_y, ts_xxzzz_yy, ts_xxzzz_yz, ts_xxzzz_z, ts_xxzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzzz_xx[i] = 6.0 * ts_xxzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xy[i] = 6.0 * ts_xxzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_xz[i] = 6.0 * ts_xxzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yy[i] = 6.0 * ts_xxzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_yz[i] = 6.0 * ts_xxzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_zz[i] = 6.0 * ts_xxzz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 312-318 components of targeted buffer : HD

    auto gs_z_xyyyy_xx = pbuffer.data(idx_g_hd + 312);

    auto gs_z_xyyyy_xy = pbuffer.data(idx_g_hd + 313);

    auto gs_z_xyyyy_xz = pbuffer.data(idx_g_hd + 314);

    auto gs_z_xyyyy_yy = pbuffer.data(idx_g_hd + 315);

    auto gs_z_xyyyy_yz = pbuffer.data(idx_g_hd + 316);

    auto gs_z_xyyyy_zz = pbuffer.data(idx_g_hd + 317);

    #pragma omp simd aligned(gc_z, gs_z_xyyyy_xx, gs_z_xyyyy_xy, gs_z_xyyyy_xz, gs_z_xyyyy_yy, gs_z_xyyyy_yz, gs_z_xyyyy_zz, ts_xyyyy_x, ts_xyyyy_xx, ts_xyyyy_xy, ts_xyyyy_xz, ts_xyyyy_y, ts_xyyyy_yy, ts_xyyyy_yz, ts_xyyyy_z, ts_xyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyy_xx[i] = 2.0 * ts_xyyyy_xx[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xy[i] = 2.0 * ts_xyyyy_xy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_xz[i] = 2.0 * ts_xyyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_xz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yy[i] = 2.0 * ts_xyyyy_yy[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_yz[i] = 2.0 * ts_xyyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_yz[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_zz[i] = 4.0 * ts_xyyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 318-324 components of targeted buffer : HD

    auto gs_z_xyyyz_xx = pbuffer.data(idx_g_hd + 318);

    auto gs_z_xyyyz_xy = pbuffer.data(idx_g_hd + 319);

    auto gs_z_xyyyz_xz = pbuffer.data(idx_g_hd + 320);

    auto gs_z_xyyyz_yy = pbuffer.data(idx_g_hd + 321);

    auto gs_z_xyyyz_yz = pbuffer.data(idx_g_hd + 322);

    auto gs_z_xyyyz_zz = pbuffer.data(idx_g_hd + 323);

    #pragma omp simd aligned(gc_z, gs_z_xyyyz_xx, gs_z_xyyyz_xy, gs_z_xyyyz_xz, gs_z_xyyyz_yy, gs_z_xyyyz_yz, gs_z_xyyyz_zz, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_zz, ts_xyyyz_x, ts_xyyyz_xx, ts_xyyyz_xy, ts_xyyyz_xz, ts_xyyyz_y, ts_xyyyz_yy, ts_xyyyz_yz, ts_xyyyz_z, ts_xyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyyz_xx[i] = 2.0 * ts_xyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xx[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xy[i] = 2.0 * ts_xyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_xz[i] = 2.0 * ts_xyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_xz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yy[i] = 2.0 * ts_xyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yy[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_yz[i] = 2.0 * ts_xyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_yz[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_zz[i] = 2.0 * ts_xyyy_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 324-330 components of targeted buffer : HD

    auto gs_z_xyyzz_xx = pbuffer.data(idx_g_hd + 324);

    auto gs_z_xyyzz_xy = pbuffer.data(idx_g_hd + 325);

    auto gs_z_xyyzz_xz = pbuffer.data(idx_g_hd + 326);

    auto gs_z_xyyzz_yy = pbuffer.data(idx_g_hd + 327);

    auto gs_z_xyyzz_yz = pbuffer.data(idx_g_hd + 328);

    auto gs_z_xyyzz_zz = pbuffer.data(idx_g_hd + 329);

    #pragma omp simd aligned(gc_z, gs_z_xyyzz_xx, gs_z_xyyzz_xy, gs_z_xyyzz_xz, gs_z_xyyzz_yy, gs_z_xyyzz_yz, gs_z_xyyzz_zz, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_zz, ts_xyyzz_x, ts_xyyzz_xx, ts_xyyzz_xy, ts_xyyzz_xz, ts_xyyzz_y, ts_xyyzz_yy, ts_xyyzz_yz, ts_xyyzz_z, ts_xyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyzz_xx[i] = 4.0 * ts_xyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xy[i] = 4.0 * ts_xyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_xz[i] = 4.0 * ts_xyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yy[i] = 4.0 * ts_xyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_yz[i] = 4.0 * ts_xyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_zz[i] = 4.0 * ts_xyyz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 330-336 components of targeted buffer : HD

    auto gs_z_xyzzz_xx = pbuffer.data(idx_g_hd + 330);

    auto gs_z_xyzzz_xy = pbuffer.data(idx_g_hd + 331);

    auto gs_z_xyzzz_xz = pbuffer.data(idx_g_hd + 332);

    auto gs_z_xyzzz_yy = pbuffer.data(idx_g_hd + 333);

    auto gs_z_xyzzz_yz = pbuffer.data(idx_g_hd + 334);

    auto gs_z_xyzzz_zz = pbuffer.data(idx_g_hd + 335);

    #pragma omp simd aligned(gc_z, gs_z_xyzzz_xx, gs_z_xyzzz_xy, gs_z_xyzzz_xz, gs_z_xyzzz_yy, gs_z_xyzzz_yz, gs_z_xyzzz_zz, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_zz, ts_xyzzz_x, ts_xyzzz_xx, ts_xyzzz_xy, ts_xyzzz_xz, ts_xyzzz_y, ts_xyzzz_yy, ts_xyzzz_yz, ts_xyzzz_z, ts_xyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzzz_xx[i] = 6.0 * ts_xyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xy[i] = 6.0 * ts_xyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_xz[i] = 6.0 * ts_xyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yy[i] = 6.0 * ts_xyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_yz[i] = 6.0 * ts_xyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_zz[i] = 6.0 * ts_xyzz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 336-342 components of targeted buffer : HD

    auto gs_z_xzzzz_xx = pbuffer.data(idx_g_hd + 336);

    auto gs_z_xzzzz_xy = pbuffer.data(idx_g_hd + 337);

    auto gs_z_xzzzz_xz = pbuffer.data(idx_g_hd + 338);

    auto gs_z_xzzzz_yy = pbuffer.data(idx_g_hd + 339);

    auto gs_z_xzzzz_yz = pbuffer.data(idx_g_hd + 340);

    auto gs_z_xzzzz_zz = pbuffer.data(idx_g_hd + 341);

    #pragma omp simd aligned(gc_z, gs_z_xzzzz_xx, gs_z_xzzzz_xy, gs_z_xzzzz_xz, gs_z_xzzzz_yy, gs_z_xzzzz_yz, gs_z_xzzzz_zz, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_zz, ts_xzzzz_x, ts_xzzzz_xx, ts_xzzzz_xy, ts_xzzzz_xz, ts_xzzzz_y, ts_xzzzz_yy, ts_xzzzz_yz, ts_xzzzz_z, ts_xzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzzz_xx[i] = 8.0 * ts_xzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xy[i] = 8.0 * ts_xzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_xz[i] = 8.0 * ts_xzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yy[i] = 8.0 * ts_xzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_yz[i] = 8.0 * ts_xzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_zz[i] = 8.0 * ts_xzzz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 342-348 components of targeted buffer : HD

    auto gs_z_yyyyy_xx = pbuffer.data(idx_g_hd + 342);

    auto gs_z_yyyyy_xy = pbuffer.data(idx_g_hd + 343);

    auto gs_z_yyyyy_xz = pbuffer.data(idx_g_hd + 344);

    auto gs_z_yyyyy_yy = pbuffer.data(idx_g_hd + 345);

    auto gs_z_yyyyy_yz = pbuffer.data(idx_g_hd + 346);

    auto gs_z_yyyyy_zz = pbuffer.data(idx_g_hd + 347);

    #pragma omp simd aligned(gc_z, gs_z_yyyyy_xx, gs_z_yyyyy_xy, gs_z_yyyyy_xz, gs_z_yyyyy_yy, gs_z_yyyyy_yz, gs_z_yyyyy_zz, ts_yyyyy_x, ts_yyyyy_xx, ts_yyyyy_xy, ts_yyyyy_xz, ts_yyyyy_y, ts_yyyyy_yy, ts_yyyyy_yz, ts_yyyyy_z, ts_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyy_xx[i] = 2.0 * ts_yyyyy_xx[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xy[i] = 2.0 * ts_yyyyy_xy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_xz[i] = 2.0 * ts_yyyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_xz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yy[i] = 2.0 * ts_yyyyy_yy[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_yz[i] = 2.0 * ts_yyyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_yz[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_zz[i] = 4.0 * ts_yyyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 348-354 components of targeted buffer : HD

    auto gs_z_yyyyz_xx = pbuffer.data(idx_g_hd + 348);

    auto gs_z_yyyyz_xy = pbuffer.data(idx_g_hd + 349);

    auto gs_z_yyyyz_xz = pbuffer.data(idx_g_hd + 350);

    auto gs_z_yyyyz_yy = pbuffer.data(idx_g_hd + 351);

    auto gs_z_yyyyz_yz = pbuffer.data(idx_g_hd + 352);

    auto gs_z_yyyyz_zz = pbuffer.data(idx_g_hd + 353);

    #pragma omp simd aligned(gc_z, gs_z_yyyyz_xx, gs_z_yyyyz_xy, gs_z_yyyyz_xz, gs_z_yyyyz_yy, gs_z_yyyyz_yz, gs_z_yyyyz_zz, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_zz, ts_yyyyz_x, ts_yyyyz_xx, ts_yyyyz_xy, ts_yyyyz_xz, ts_yyyyz_y, ts_yyyyz_yy, ts_yyyyz_yz, ts_yyyyz_z, ts_yyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyyz_xx[i] = 2.0 * ts_yyyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xx[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xy[i] = 2.0 * ts_yyyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_xz[i] = 2.0 * ts_yyyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_xz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yy[i] = 2.0 * ts_yyyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yy[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_yz[i] = 2.0 * ts_yyyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_yz[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_zz[i] = 2.0 * ts_yyyy_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 354-360 components of targeted buffer : HD

    auto gs_z_yyyzz_xx = pbuffer.data(idx_g_hd + 354);

    auto gs_z_yyyzz_xy = pbuffer.data(idx_g_hd + 355);

    auto gs_z_yyyzz_xz = pbuffer.data(idx_g_hd + 356);

    auto gs_z_yyyzz_yy = pbuffer.data(idx_g_hd + 357);

    auto gs_z_yyyzz_yz = pbuffer.data(idx_g_hd + 358);

    auto gs_z_yyyzz_zz = pbuffer.data(idx_g_hd + 359);

    #pragma omp simd aligned(gc_z, gs_z_yyyzz_xx, gs_z_yyyzz_xy, gs_z_yyyzz_xz, gs_z_yyyzz_yy, gs_z_yyyzz_yz, gs_z_yyyzz_zz, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_zz, ts_yyyzz_x, ts_yyyzz_xx, ts_yyyzz_xy, ts_yyyzz_xz, ts_yyyzz_y, ts_yyyzz_yy, ts_yyyzz_yz, ts_yyyzz_z, ts_yyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyzz_xx[i] = 4.0 * ts_yyyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xx[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xy[i] = 4.0 * ts_yyyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_xz[i] = 4.0 * ts_yyyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_xz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yy[i] = 4.0 * ts_yyyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yy[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_yz[i] = 4.0 * ts_yyyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_yz[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_zz[i] = 4.0 * ts_yyyz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 360-366 components of targeted buffer : HD

    auto gs_z_yyzzz_xx = pbuffer.data(idx_g_hd + 360);

    auto gs_z_yyzzz_xy = pbuffer.data(idx_g_hd + 361);

    auto gs_z_yyzzz_xz = pbuffer.data(idx_g_hd + 362);

    auto gs_z_yyzzz_yy = pbuffer.data(idx_g_hd + 363);

    auto gs_z_yyzzz_yz = pbuffer.data(idx_g_hd + 364);

    auto gs_z_yyzzz_zz = pbuffer.data(idx_g_hd + 365);

    #pragma omp simd aligned(gc_z, gs_z_yyzzz_xx, gs_z_yyzzz_xy, gs_z_yyzzz_xz, gs_z_yyzzz_yy, gs_z_yyzzz_yz, gs_z_yyzzz_zz, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_zz, ts_yyzzz_x, ts_yyzzz_xx, ts_yyzzz_xy, ts_yyzzz_xz, ts_yyzzz_y, ts_yyzzz_yy, ts_yyzzz_yz, ts_yyzzz_z, ts_yyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzzz_xx[i] = 6.0 * ts_yyzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xx[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xy[i] = 6.0 * ts_yyzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_xz[i] = 6.0 * ts_yyzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_xz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yy[i] = 6.0 * ts_yyzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yy[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_yz[i] = 6.0 * ts_yyzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_yz[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_zz[i] = 6.0 * ts_yyzz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 366-372 components of targeted buffer : HD

    auto gs_z_yzzzz_xx = pbuffer.data(idx_g_hd + 366);

    auto gs_z_yzzzz_xy = pbuffer.data(idx_g_hd + 367);

    auto gs_z_yzzzz_xz = pbuffer.data(idx_g_hd + 368);

    auto gs_z_yzzzz_yy = pbuffer.data(idx_g_hd + 369);

    auto gs_z_yzzzz_yz = pbuffer.data(idx_g_hd + 370);

    auto gs_z_yzzzz_zz = pbuffer.data(idx_g_hd + 371);

    #pragma omp simd aligned(gc_z, gs_z_yzzzz_xx, gs_z_yzzzz_xy, gs_z_yzzzz_xz, gs_z_yzzzz_yy, gs_z_yzzzz_yz, gs_z_yzzzz_zz, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_zz, ts_yzzzz_x, ts_yzzzz_xx, ts_yzzzz_xy, ts_yzzzz_xz, ts_yzzzz_y, ts_yzzzz_yy, ts_yzzzz_yz, ts_yzzzz_z, ts_yzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzzz_xx[i] = 8.0 * ts_yzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xx[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xy[i] = 8.0 * ts_yzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_xz[i] = 8.0 * ts_yzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_xz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yy[i] = 8.0 * ts_yzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yy[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_yz[i] = 8.0 * ts_yzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_yz[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_zz[i] = 8.0 * ts_yzzz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 372-378 components of targeted buffer : HD

    auto gs_z_zzzzz_xx = pbuffer.data(idx_g_hd + 372);

    auto gs_z_zzzzz_xy = pbuffer.data(idx_g_hd + 373);

    auto gs_z_zzzzz_xz = pbuffer.data(idx_g_hd + 374);

    auto gs_z_zzzzz_yy = pbuffer.data(idx_g_hd + 375);

    auto gs_z_zzzzz_yz = pbuffer.data(idx_g_hd + 376);

    auto gs_z_zzzzz_zz = pbuffer.data(idx_g_hd + 377);

    #pragma omp simd aligned(gc_z, gs_z_zzzzz_xx, gs_z_zzzzz_xy, gs_z_zzzzz_xz, gs_z_zzzzz_yy, gs_z_zzzzz_yz, gs_z_zzzzz_zz, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_zz, ts_zzzzz_x, ts_zzzzz_xx, ts_zzzzz_xy, ts_zzzzz_xz, ts_zzzzz_y, ts_zzzzz_yy, ts_zzzzz_yz, ts_zzzzz_z, ts_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzzz_xx[i] = 10.0 * ts_zzzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xx[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xy[i] = 10.0 * ts_zzzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_xz[i] = 10.0 * ts_zzzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_xz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yy[i] = 10.0 * ts_zzzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yy[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_yz[i] = 10.0 * ts_zzzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_yz[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_zz[i] = 10.0 * ts_zzzz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_zz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

