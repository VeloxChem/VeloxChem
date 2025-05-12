#include "ThreeCenterOverlapGradientPrimRecDH.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_dh(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_dh,
                              const size_t idx_ph,
                              const size_t idx_dg,
                              const size_t idx_dh,
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

    // Set up components of auxiliary buffer : PH

    auto ts_x_xxxxx = pbuffer.data(idx_ph);

    auto ts_x_xxxxy = pbuffer.data(idx_ph + 1);

    auto ts_x_xxxxz = pbuffer.data(idx_ph + 2);

    auto ts_x_xxxyy = pbuffer.data(idx_ph + 3);

    auto ts_x_xxxyz = pbuffer.data(idx_ph + 4);

    auto ts_x_xxxzz = pbuffer.data(idx_ph + 5);

    auto ts_x_xxyyy = pbuffer.data(idx_ph + 6);

    auto ts_x_xxyyz = pbuffer.data(idx_ph + 7);

    auto ts_x_xxyzz = pbuffer.data(idx_ph + 8);

    auto ts_x_xxzzz = pbuffer.data(idx_ph + 9);

    auto ts_x_xyyyy = pbuffer.data(idx_ph + 10);

    auto ts_x_xyyyz = pbuffer.data(idx_ph + 11);

    auto ts_x_xyyzz = pbuffer.data(idx_ph + 12);

    auto ts_x_xyzzz = pbuffer.data(idx_ph + 13);

    auto ts_x_xzzzz = pbuffer.data(idx_ph + 14);

    auto ts_x_yyyyy = pbuffer.data(idx_ph + 15);

    auto ts_x_yyyyz = pbuffer.data(idx_ph + 16);

    auto ts_x_yyyzz = pbuffer.data(idx_ph + 17);

    auto ts_x_yyzzz = pbuffer.data(idx_ph + 18);

    auto ts_x_yzzzz = pbuffer.data(idx_ph + 19);

    auto ts_x_zzzzz = pbuffer.data(idx_ph + 20);

    auto ts_y_xxxxx = pbuffer.data(idx_ph + 21);

    auto ts_y_xxxxy = pbuffer.data(idx_ph + 22);

    auto ts_y_xxxxz = pbuffer.data(idx_ph + 23);

    auto ts_y_xxxyy = pbuffer.data(idx_ph + 24);

    auto ts_y_xxxyz = pbuffer.data(idx_ph + 25);

    auto ts_y_xxxzz = pbuffer.data(idx_ph + 26);

    auto ts_y_xxyyy = pbuffer.data(idx_ph + 27);

    auto ts_y_xxyyz = pbuffer.data(idx_ph + 28);

    auto ts_y_xxyzz = pbuffer.data(idx_ph + 29);

    auto ts_y_xxzzz = pbuffer.data(idx_ph + 30);

    auto ts_y_xyyyy = pbuffer.data(idx_ph + 31);

    auto ts_y_xyyyz = pbuffer.data(idx_ph + 32);

    auto ts_y_xyyzz = pbuffer.data(idx_ph + 33);

    auto ts_y_xyzzz = pbuffer.data(idx_ph + 34);

    auto ts_y_xzzzz = pbuffer.data(idx_ph + 35);

    auto ts_y_yyyyy = pbuffer.data(idx_ph + 36);

    auto ts_y_yyyyz = pbuffer.data(idx_ph + 37);

    auto ts_y_yyyzz = pbuffer.data(idx_ph + 38);

    auto ts_y_yyzzz = pbuffer.data(idx_ph + 39);

    auto ts_y_yzzzz = pbuffer.data(idx_ph + 40);

    auto ts_y_zzzzz = pbuffer.data(idx_ph + 41);

    auto ts_z_xxxxx = pbuffer.data(idx_ph + 42);

    auto ts_z_xxxxy = pbuffer.data(idx_ph + 43);

    auto ts_z_xxxxz = pbuffer.data(idx_ph + 44);

    auto ts_z_xxxyy = pbuffer.data(idx_ph + 45);

    auto ts_z_xxxyz = pbuffer.data(idx_ph + 46);

    auto ts_z_xxxzz = pbuffer.data(idx_ph + 47);

    auto ts_z_xxyyy = pbuffer.data(idx_ph + 48);

    auto ts_z_xxyyz = pbuffer.data(idx_ph + 49);

    auto ts_z_xxyzz = pbuffer.data(idx_ph + 50);

    auto ts_z_xxzzz = pbuffer.data(idx_ph + 51);

    auto ts_z_xyyyy = pbuffer.data(idx_ph + 52);

    auto ts_z_xyyyz = pbuffer.data(idx_ph + 53);

    auto ts_z_xyyzz = pbuffer.data(idx_ph + 54);

    auto ts_z_xyzzz = pbuffer.data(idx_ph + 55);

    auto ts_z_xzzzz = pbuffer.data(idx_ph + 56);

    auto ts_z_yyyyy = pbuffer.data(idx_ph + 57);

    auto ts_z_yyyyz = pbuffer.data(idx_ph + 58);

    auto ts_z_yyyzz = pbuffer.data(idx_ph + 59);

    auto ts_z_yyzzz = pbuffer.data(idx_ph + 60);

    auto ts_z_yzzzz = pbuffer.data(idx_ph + 61);

    auto ts_z_zzzzz = pbuffer.data(idx_ph + 62);

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

    auto ts_xy_xxxx = pbuffer.data(idx_dg + 15);

    auto ts_xy_xxxy = pbuffer.data(idx_dg + 16);

    auto ts_xy_xxxz = pbuffer.data(idx_dg + 17);

    auto ts_xy_xxyy = pbuffer.data(idx_dg + 18);

    auto ts_xy_xxyz = pbuffer.data(idx_dg + 19);

    auto ts_xy_xxzz = pbuffer.data(idx_dg + 20);

    auto ts_xy_xyyy = pbuffer.data(idx_dg + 21);

    auto ts_xy_xyyz = pbuffer.data(idx_dg + 22);

    auto ts_xy_xyzz = pbuffer.data(idx_dg + 23);

    auto ts_xy_xzzz = pbuffer.data(idx_dg + 24);

    auto ts_xy_yyyy = pbuffer.data(idx_dg + 25);

    auto ts_xy_yyyz = pbuffer.data(idx_dg + 26);

    auto ts_xy_yyzz = pbuffer.data(idx_dg + 27);

    auto ts_xy_yzzz = pbuffer.data(idx_dg + 28);

    auto ts_xy_zzzz = pbuffer.data(idx_dg + 29);

    auto ts_xz_xxxx = pbuffer.data(idx_dg + 30);

    auto ts_xz_xxxy = pbuffer.data(idx_dg + 31);

    auto ts_xz_xxxz = pbuffer.data(idx_dg + 32);

    auto ts_xz_xxyy = pbuffer.data(idx_dg + 33);

    auto ts_xz_xxyz = pbuffer.data(idx_dg + 34);

    auto ts_xz_xxzz = pbuffer.data(idx_dg + 35);

    auto ts_xz_xyyy = pbuffer.data(idx_dg + 36);

    auto ts_xz_xyyz = pbuffer.data(idx_dg + 37);

    auto ts_xz_xyzz = pbuffer.data(idx_dg + 38);

    auto ts_xz_xzzz = pbuffer.data(idx_dg + 39);

    auto ts_xz_yyyy = pbuffer.data(idx_dg + 40);

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

    auto ts_yz_xxxx = pbuffer.data(idx_dg + 60);

    auto ts_yz_xxxy = pbuffer.data(idx_dg + 61);

    auto ts_yz_xxxz = pbuffer.data(idx_dg + 62);

    auto ts_yz_xxyy = pbuffer.data(idx_dg + 63);

    auto ts_yz_xxyz = pbuffer.data(idx_dg + 64);

    auto ts_yz_xxzz = pbuffer.data(idx_dg + 65);

    auto ts_yz_xyyy = pbuffer.data(idx_dg + 66);

    auto ts_yz_xyyz = pbuffer.data(idx_dg + 67);

    auto ts_yz_xyzz = pbuffer.data(idx_dg + 68);

    auto ts_yz_xzzz = pbuffer.data(idx_dg + 69);

    auto ts_yz_yyyy = pbuffer.data(idx_dg + 70);

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

    // Set up components of auxiliary buffer : DH

    auto ts_xx_xxxxx = pbuffer.data(idx_dh);

    auto ts_xx_xxxxy = pbuffer.data(idx_dh + 1);

    auto ts_xx_xxxxz = pbuffer.data(idx_dh + 2);

    auto ts_xx_xxxyy = pbuffer.data(idx_dh + 3);

    auto ts_xx_xxxyz = pbuffer.data(idx_dh + 4);

    auto ts_xx_xxxzz = pbuffer.data(idx_dh + 5);

    auto ts_xx_xxyyy = pbuffer.data(idx_dh + 6);

    auto ts_xx_xxyyz = pbuffer.data(idx_dh + 7);

    auto ts_xx_xxyzz = pbuffer.data(idx_dh + 8);

    auto ts_xx_xxzzz = pbuffer.data(idx_dh + 9);

    auto ts_xx_xyyyy = pbuffer.data(idx_dh + 10);

    auto ts_xx_xyyyz = pbuffer.data(idx_dh + 11);

    auto ts_xx_xyyzz = pbuffer.data(idx_dh + 12);

    auto ts_xx_xyzzz = pbuffer.data(idx_dh + 13);

    auto ts_xx_xzzzz = pbuffer.data(idx_dh + 14);

    auto ts_xx_yyyyy = pbuffer.data(idx_dh + 15);

    auto ts_xx_yyyyz = pbuffer.data(idx_dh + 16);

    auto ts_xx_yyyzz = pbuffer.data(idx_dh + 17);

    auto ts_xx_yyzzz = pbuffer.data(idx_dh + 18);

    auto ts_xx_yzzzz = pbuffer.data(idx_dh + 19);

    auto ts_xx_zzzzz = pbuffer.data(idx_dh + 20);

    auto ts_xy_xxxxx = pbuffer.data(idx_dh + 21);

    auto ts_xy_xxxxy = pbuffer.data(idx_dh + 22);

    auto ts_xy_xxxxz = pbuffer.data(idx_dh + 23);

    auto ts_xy_xxxyy = pbuffer.data(idx_dh + 24);

    auto ts_xy_xxxyz = pbuffer.data(idx_dh + 25);

    auto ts_xy_xxxzz = pbuffer.data(idx_dh + 26);

    auto ts_xy_xxyyy = pbuffer.data(idx_dh + 27);

    auto ts_xy_xxyyz = pbuffer.data(idx_dh + 28);

    auto ts_xy_xxyzz = pbuffer.data(idx_dh + 29);

    auto ts_xy_xxzzz = pbuffer.data(idx_dh + 30);

    auto ts_xy_xyyyy = pbuffer.data(idx_dh + 31);

    auto ts_xy_xyyyz = pbuffer.data(idx_dh + 32);

    auto ts_xy_xyyzz = pbuffer.data(idx_dh + 33);

    auto ts_xy_xyzzz = pbuffer.data(idx_dh + 34);

    auto ts_xy_xzzzz = pbuffer.data(idx_dh + 35);

    auto ts_xy_yyyyy = pbuffer.data(idx_dh + 36);

    auto ts_xy_yyyyz = pbuffer.data(idx_dh + 37);

    auto ts_xy_yyyzz = pbuffer.data(idx_dh + 38);

    auto ts_xy_yyzzz = pbuffer.data(idx_dh + 39);

    auto ts_xy_yzzzz = pbuffer.data(idx_dh + 40);

    auto ts_xy_zzzzz = pbuffer.data(idx_dh + 41);

    auto ts_xz_xxxxx = pbuffer.data(idx_dh + 42);

    auto ts_xz_xxxxy = pbuffer.data(idx_dh + 43);

    auto ts_xz_xxxxz = pbuffer.data(idx_dh + 44);

    auto ts_xz_xxxyy = pbuffer.data(idx_dh + 45);

    auto ts_xz_xxxyz = pbuffer.data(idx_dh + 46);

    auto ts_xz_xxxzz = pbuffer.data(idx_dh + 47);

    auto ts_xz_xxyyy = pbuffer.data(idx_dh + 48);

    auto ts_xz_xxyyz = pbuffer.data(idx_dh + 49);

    auto ts_xz_xxyzz = pbuffer.data(idx_dh + 50);

    auto ts_xz_xxzzz = pbuffer.data(idx_dh + 51);

    auto ts_xz_xyyyy = pbuffer.data(idx_dh + 52);

    auto ts_xz_xyyyz = pbuffer.data(idx_dh + 53);

    auto ts_xz_xyyzz = pbuffer.data(idx_dh + 54);

    auto ts_xz_xyzzz = pbuffer.data(idx_dh + 55);

    auto ts_xz_xzzzz = pbuffer.data(idx_dh + 56);

    auto ts_xz_yyyyy = pbuffer.data(idx_dh + 57);

    auto ts_xz_yyyyz = pbuffer.data(idx_dh + 58);

    auto ts_xz_yyyzz = pbuffer.data(idx_dh + 59);

    auto ts_xz_yyzzz = pbuffer.data(idx_dh + 60);

    auto ts_xz_yzzzz = pbuffer.data(idx_dh + 61);

    auto ts_xz_zzzzz = pbuffer.data(idx_dh + 62);

    auto ts_yy_xxxxx = pbuffer.data(idx_dh + 63);

    auto ts_yy_xxxxy = pbuffer.data(idx_dh + 64);

    auto ts_yy_xxxxz = pbuffer.data(idx_dh + 65);

    auto ts_yy_xxxyy = pbuffer.data(idx_dh + 66);

    auto ts_yy_xxxyz = pbuffer.data(idx_dh + 67);

    auto ts_yy_xxxzz = pbuffer.data(idx_dh + 68);

    auto ts_yy_xxyyy = pbuffer.data(idx_dh + 69);

    auto ts_yy_xxyyz = pbuffer.data(idx_dh + 70);

    auto ts_yy_xxyzz = pbuffer.data(idx_dh + 71);

    auto ts_yy_xxzzz = pbuffer.data(idx_dh + 72);

    auto ts_yy_xyyyy = pbuffer.data(idx_dh + 73);

    auto ts_yy_xyyyz = pbuffer.data(idx_dh + 74);

    auto ts_yy_xyyzz = pbuffer.data(idx_dh + 75);

    auto ts_yy_xyzzz = pbuffer.data(idx_dh + 76);

    auto ts_yy_xzzzz = pbuffer.data(idx_dh + 77);

    auto ts_yy_yyyyy = pbuffer.data(idx_dh + 78);

    auto ts_yy_yyyyz = pbuffer.data(idx_dh + 79);

    auto ts_yy_yyyzz = pbuffer.data(idx_dh + 80);

    auto ts_yy_yyzzz = pbuffer.data(idx_dh + 81);

    auto ts_yy_yzzzz = pbuffer.data(idx_dh + 82);

    auto ts_yy_zzzzz = pbuffer.data(idx_dh + 83);

    auto ts_yz_xxxxx = pbuffer.data(idx_dh + 84);

    auto ts_yz_xxxxy = pbuffer.data(idx_dh + 85);

    auto ts_yz_xxxxz = pbuffer.data(idx_dh + 86);

    auto ts_yz_xxxyy = pbuffer.data(idx_dh + 87);

    auto ts_yz_xxxyz = pbuffer.data(idx_dh + 88);

    auto ts_yz_xxxzz = pbuffer.data(idx_dh + 89);

    auto ts_yz_xxyyy = pbuffer.data(idx_dh + 90);

    auto ts_yz_xxyyz = pbuffer.data(idx_dh + 91);

    auto ts_yz_xxyzz = pbuffer.data(idx_dh + 92);

    auto ts_yz_xxzzz = pbuffer.data(idx_dh + 93);

    auto ts_yz_xyyyy = pbuffer.data(idx_dh + 94);

    auto ts_yz_xyyyz = pbuffer.data(idx_dh + 95);

    auto ts_yz_xyyzz = pbuffer.data(idx_dh + 96);

    auto ts_yz_xyzzz = pbuffer.data(idx_dh + 97);

    auto ts_yz_xzzzz = pbuffer.data(idx_dh + 98);

    auto ts_yz_yyyyy = pbuffer.data(idx_dh + 99);

    auto ts_yz_yyyyz = pbuffer.data(idx_dh + 100);

    auto ts_yz_yyyzz = pbuffer.data(idx_dh + 101);

    auto ts_yz_yyzzz = pbuffer.data(idx_dh + 102);

    auto ts_yz_yzzzz = pbuffer.data(idx_dh + 103);

    auto ts_yz_zzzzz = pbuffer.data(idx_dh + 104);

    auto ts_zz_xxxxx = pbuffer.data(idx_dh + 105);

    auto ts_zz_xxxxy = pbuffer.data(idx_dh + 106);

    auto ts_zz_xxxxz = pbuffer.data(idx_dh + 107);

    auto ts_zz_xxxyy = pbuffer.data(idx_dh + 108);

    auto ts_zz_xxxyz = pbuffer.data(idx_dh + 109);

    auto ts_zz_xxxzz = pbuffer.data(idx_dh + 110);

    auto ts_zz_xxyyy = pbuffer.data(idx_dh + 111);

    auto ts_zz_xxyyz = pbuffer.data(idx_dh + 112);

    auto ts_zz_xxyzz = pbuffer.data(idx_dh + 113);

    auto ts_zz_xxzzz = pbuffer.data(idx_dh + 114);

    auto ts_zz_xyyyy = pbuffer.data(idx_dh + 115);

    auto ts_zz_xyyyz = pbuffer.data(idx_dh + 116);

    auto ts_zz_xyyzz = pbuffer.data(idx_dh + 117);

    auto ts_zz_xyzzz = pbuffer.data(idx_dh + 118);

    auto ts_zz_xzzzz = pbuffer.data(idx_dh + 119);

    auto ts_zz_yyyyy = pbuffer.data(idx_dh + 120);

    auto ts_zz_yyyyz = pbuffer.data(idx_dh + 121);

    auto ts_zz_yyyzz = pbuffer.data(idx_dh + 122);

    auto ts_zz_yyzzz = pbuffer.data(idx_dh + 123);

    auto ts_zz_yzzzz = pbuffer.data(idx_dh + 124);

    auto ts_zz_zzzzz = pbuffer.data(idx_dh + 125);

    // Set up 0-21 components of targeted buffer : DH

    auto gs_x_xx_xxxxx = pbuffer.data(idx_g_dh);

    auto gs_x_xx_xxxxy = pbuffer.data(idx_g_dh + 1);

    auto gs_x_xx_xxxxz = pbuffer.data(idx_g_dh + 2);

    auto gs_x_xx_xxxyy = pbuffer.data(idx_g_dh + 3);

    auto gs_x_xx_xxxyz = pbuffer.data(idx_g_dh + 4);

    auto gs_x_xx_xxxzz = pbuffer.data(idx_g_dh + 5);

    auto gs_x_xx_xxyyy = pbuffer.data(idx_g_dh + 6);

    auto gs_x_xx_xxyyz = pbuffer.data(idx_g_dh + 7);

    auto gs_x_xx_xxyzz = pbuffer.data(idx_g_dh + 8);

    auto gs_x_xx_xxzzz = pbuffer.data(idx_g_dh + 9);

    auto gs_x_xx_xyyyy = pbuffer.data(idx_g_dh + 10);

    auto gs_x_xx_xyyyz = pbuffer.data(idx_g_dh + 11);

    auto gs_x_xx_xyyzz = pbuffer.data(idx_g_dh + 12);

    auto gs_x_xx_xyzzz = pbuffer.data(idx_g_dh + 13);

    auto gs_x_xx_xzzzz = pbuffer.data(idx_g_dh + 14);

    auto gs_x_xx_yyyyy = pbuffer.data(idx_g_dh + 15);

    auto gs_x_xx_yyyyz = pbuffer.data(idx_g_dh + 16);

    auto gs_x_xx_yyyzz = pbuffer.data(idx_g_dh + 17);

    auto gs_x_xx_yyzzz = pbuffer.data(idx_g_dh + 18);

    auto gs_x_xx_yzzzz = pbuffer.data(idx_g_dh + 19);

    auto gs_x_xx_zzzzz = pbuffer.data(idx_g_dh + 20);

    #pragma omp simd aligned(gc_x, gs_x_xx_xxxxx, gs_x_xx_xxxxy, gs_x_xx_xxxxz, gs_x_xx_xxxyy, gs_x_xx_xxxyz, gs_x_xx_xxxzz, gs_x_xx_xxyyy, gs_x_xx_xxyyz, gs_x_xx_xxyzz, gs_x_xx_xxzzz, gs_x_xx_xyyyy, gs_x_xx_xyyyz, gs_x_xx_xyyzz, gs_x_xx_xyzzz, gs_x_xx_xzzzz, gs_x_xx_yyyyy, gs_x_xx_yyyyz, gs_x_xx_yyyzz, gs_x_xx_yyzzz, gs_x_xx_yzzzz, gs_x_xx_zzzzz, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxxz, ts_x_xxxyy, ts_x_xxxyz, ts_x_xxxzz, ts_x_xxyyy, ts_x_xxyyz, ts_x_xxyzz, ts_x_xxzzz, ts_x_xyyyy, ts_x_xyyyz, ts_x_xyyzz, ts_x_xyzzz, ts_x_xzzzz, ts_x_yyyyy, ts_x_yyyyz, ts_x_yyyzz, ts_x_yyzzz, ts_x_yzzzz, ts_x_zzzzz, ts_xx_xxxx, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxxz, ts_xx_xxxy, ts_xx_xxxyy, ts_xx_xxxyz, ts_xx_xxxz, ts_xx_xxxzz, ts_xx_xxyy, ts_xx_xxyyy, ts_xx_xxyyz, ts_xx_xxyz, ts_xx_xxyzz, ts_xx_xxzz, ts_xx_xxzzz, ts_xx_xyyy, ts_xx_xyyyy, ts_xx_xyyyz, ts_xx_xyyz, ts_xx_xyyzz, ts_xx_xyzz, ts_xx_xyzzz, ts_xx_xzzz, ts_xx_xzzzz, ts_xx_yyyy, ts_xx_yyyyy, ts_xx_yyyyz, ts_xx_yyyz, ts_xx_yyyzz, ts_xx_yyzz, ts_xx_yyzzz, ts_xx_yzzz, ts_xx_yzzzz, ts_xx_zzzz, ts_xx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xx_xxxxx[i] = 4.0 * ts_x_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xx_xxxxy[i] = 4.0 * ts_x_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xx_xxxxz[i] = 4.0 * ts_x_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xx_xxxyy[i] = 4.0 * ts_x_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xx_xxxyz[i] = 4.0 * ts_x_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xx_xxxzz[i] = 4.0 * ts_x_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xx_xxyyy[i] = 4.0 * ts_x_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xx_xxyyz[i] = 4.0 * ts_x_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xx_xxyzz[i] = 4.0 * ts_x_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xx_xxzzz[i] = 4.0 * ts_x_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xx_xyyyy[i] = 4.0 * ts_x_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xx_xyyyz[i] = 4.0 * ts_x_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xx_xyyzz[i] = 4.0 * ts_x_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xx_xyzzz[i] = 4.0 * ts_x_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xx_xzzzz[i] = 4.0 * ts_x_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xx_yyyyy[i] = 4.0 * ts_x_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xx_yyyyz[i] = 4.0 * ts_x_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xx_yyyzz[i] = 4.0 * ts_x_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xx_yyzzz[i] = 4.0 * ts_x_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xx_yzzzz[i] = 4.0 * ts_x_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xx_zzzzz[i] = 4.0 * ts_x_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 21-42 components of targeted buffer : DH

    auto gs_x_xy_xxxxx = pbuffer.data(idx_g_dh + 21);

    auto gs_x_xy_xxxxy = pbuffer.data(idx_g_dh + 22);

    auto gs_x_xy_xxxxz = pbuffer.data(idx_g_dh + 23);

    auto gs_x_xy_xxxyy = pbuffer.data(idx_g_dh + 24);

    auto gs_x_xy_xxxyz = pbuffer.data(idx_g_dh + 25);

    auto gs_x_xy_xxxzz = pbuffer.data(idx_g_dh + 26);

    auto gs_x_xy_xxyyy = pbuffer.data(idx_g_dh + 27);

    auto gs_x_xy_xxyyz = pbuffer.data(idx_g_dh + 28);

    auto gs_x_xy_xxyzz = pbuffer.data(idx_g_dh + 29);

    auto gs_x_xy_xxzzz = pbuffer.data(idx_g_dh + 30);

    auto gs_x_xy_xyyyy = pbuffer.data(idx_g_dh + 31);

    auto gs_x_xy_xyyyz = pbuffer.data(idx_g_dh + 32);

    auto gs_x_xy_xyyzz = pbuffer.data(idx_g_dh + 33);

    auto gs_x_xy_xyzzz = pbuffer.data(idx_g_dh + 34);

    auto gs_x_xy_xzzzz = pbuffer.data(idx_g_dh + 35);

    auto gs_x_xy_yyyyy = pbuffer.data(idx_g_dh + 36);

    auto gs_x_xy_yyyyz = pbuffer.data(idx_g_dh + 37);

    auto gs_x_xy_yyyzz = pbuffer.data(idx_g_dh + 38);

    auto gs_x_xy_yyzzz = pbuffer.data(idx_g_dh + 39);

    auto gs_x_xy_yzzzz = pbuffer.data(idx_g_dh + 40);

    auto gs_x_xy_zzzzz = pbuffer.data(idx_g_dh + 41);

    #pragma omp simd aligned(gc_x, gs_x_xy_xxxxx, gs_x_xy_xxxxy, gs_x_xy_xxxxz, gs_x_xy_xxxyy, gs_x_xy_xxxyz, gs_x_xy_xxxzz, gs_x_xy_xxyyy, gs_x_xy_xxyyz, gs_x_xy_xxyzz, gs_x_xy_xxzzz, gs_x_xy_xyyyy, gs_x_xy_xyyyz, gs_x_xy_xyyzz, gs_x_xy_xyzzz, gs_x_xy_xzzzz, gs_x_xy_yyyyy, gs_x_xy_yyyyz, gs_x_xy_yyyzz, gs_x_xy_yyzzz, gs_x_xy_yzzzz, gs_x_xy_zzzzz, ts_xy_xxxx, ts_xy_xxxxx, ts_xy_xxxxy, ts_xy_xxxxz, ts_xy_xxxy, ts_xy_xxxyy, ts_xy_xxxyz, ts_xy_xxxz, ts_xy_xxxzz, ts_xy_xxyy, ts_xy_xxyyy, ts_xy_xxyyz, ts_xy_xxyz, ts_xy_xxyzz, ts_xy_xxzz, ts_xy_xxzzz, ts_xy_xyyy, ts_xy_xyyyy, ts_xy_xyyyz, ts_xy_xyyz, ts_xy_xyyzz, ts_xy_xyzz, ts_xy_xyzzz, ts_xy_xzzz, ts_xy_xzzzz, ts_xy_yyyy, ts_xy_yyyyy, ts_xy_yyyyz, ts_xy_yyyz, ts_xy_yyyzz, ts_xy_yyzz, ts_xy_yyzzz, ts_xy_yzzz, ts_xy_yzzzz, ts_xy_zzzz, ts_xy_zzzzz, ts_y_xxxxx, ts_y_xxxxy, ts_y_xxxxz, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxxzz, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyzz, ts_y_xxzzz, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyzz, ts_y_xyzzz, ts_y_xzzzz, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyzz, ts_y_yyzzz, ts_y_yzzzz, ts_y_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xy_xxxxx[i] = 2.0 * ts_y_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xy_xxxxy[i] = 2.0 * ts_y_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xy_xxxxz[i] = 2.0 * ts_y_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xy_xxxyy[i] = 2.0 * ts_y_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xy_xxxyz[i] = 2.0 * ts_y_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xy_xxxzz[i] = 2.0 * ts_y_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xy_xxyyy[i] = 2.0 * ts_y_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xy_xxyyz[i] = 2.0 * ts_y_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xy_xxyzz[i] = 2.0 * ts_y_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xy_xxzzz[i] = 2.0 * ts_y_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xy_xyyyy[i] = 2.0 * ts_y_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xy_xyyyz[i] = 2.0 * ts_y_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xy_xyyzz[i] = 2.0 * ts_y_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xy_xyzzz[i] = 2.0 * ts_y_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xy_xzzzz[i] = 2.0 * ts_y_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xy_yyyyy[i] = 2.0 * ts_y_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xy_yyyyz[i] = 2.0 * ts_y_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xy_yyyzz[i] = 2.0 * ts_y_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xy_yyzzz[i] = 2.0 * ts_y_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xy_yzzzz[i] = 2.0 * ts_y_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xy_zzzzz[i] = 2.0 * ts_y_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 42-63 components of targeted buffer : DH

    auto gs_x_xz_xxxxx = pbuffer.data(idx_g_dh + 42);

    auto gs_x_xz_xxxxy = pbuffer.data(idx_g_dh + 43);

    auto gs_x_xz_xxxxz = pbuffer.data(idx_g_dh + 44);

    auto gs_x_xz_xxxyy = pbuffer.data(idx_g_dh + 45);

    auto gs_x_xz_xxxyz = pbuffer.data(idx_g_dh + 46);

    auto gs_x_xz_xxxzz = pbuffer.data(idx_g_dh + 47);

    auto gs_x_xz_xxyyy = pbuffer.data(idx_g_dh + 48);

    auto gs_x_xz_xxyyz = pbuffer.data(idx_g_dh + 49);

    auto gs_x_xz_xxyzz = pbuffer.data(idx_g_dh + 50);

    auto gs_x_xz_xxzzz = pbuffer.data(idx_g_dh + 51);

    auto gs_x_xz_xyyyy = pbuffer.data(idx_g_dh + 52);

    auto gs_x_xz_xyyyz = pbuffer.data(idx_g_dh + 53);

    auto gs_x_xz_xyyzz = pbuffer.data(idx_g_dh + 54);

    auto gs_x_xz_xyzzz = pbuffer.data(idx_g_dh + 55);

    auto gs_x_xz_xzzzz = pbuffer.data(idx_g_dh + 56);

    auto gs_x_xz_yyyyy = pbuffer.data(idx_g_dh + 57);

    auto gs_x_xz_yyyyz = pbuffer.data(idx_g_dh + 58);

    auto gs_x_xz_yyyzz = pbuffer.data(idx_g_dh + 59);

    auto gs_x_xz_yyzzz = pbuffer.data(idx_g_dh + 60);

    auto gs_x_xz_yzzzz = pbuffer.data(idx_g_dh + 61);

    auto gs_x_xz_zzzzz = pbuffer.data(idx_g_dh + 62);

    #pragma omp simd aligned(gc_x, gs_x_xz_xxxxx, gs_x_xz_xxxxy, gs_x_xz_xxxxz, gs_x_xz_xxxyy, gs_x_xz_xxxyz, gs_x_xz_xxxzz, gs_x_xz_xxyyy, gs_x_xz_xxyyz, gs_x_xz_xxyzz, gs_x_xz_xxzzz, gs_x_xz_xyyyy, gs_x_xz_xyyyz, gs_x_xz_xyyzz, gs_x_xz_xyzzz, gs_x_xz_xzzzz, gs_x_xz_yyyyy, gs_x_xz_yyyyz, gs_x_xz_yyyzz, gs_x_xz_yyzzz, gs_x_xz_yzzzz, gs_x_xz_zzzzz, ts_xz_xxxx, ts_xz_xxxxx, ts_xz_xxxxy, ts_xz_xxxxz, ts_xz_xxxy, ts_xz_xxxyy, ts_xz_xxxyz, ts_xz_xxxz, ts_xz_xxxzz, ts_xz_xxyy, ts_xz_xxyyy, ts_xz_xxyyz, ts_xz_xxyz, ts_xz_xxyzz, ts_xz_xxzz, ts_xz_xxzzz, ts_xz_xyyy, ts_xz_xyyyy, ts_xz_xyyyz, ts_xz_xyyz, ts_xz_xyyzz, ts_xz_xyzz, ts_xz_xyzzz, ts_xz_xzzz, ts_xz_xzzzz, ts_xz_yyyy, ts_xz_yyyyy, ts_xz_yyyyz, ts_xz_yyyz, ts_xz_yyyzz, ts_xz_yyzz, ts_xz_yyzzz, ts_xz_yzzz, ts_xz_yzzzz, ts_xz_zzzz, ts_xz_zzzzz, ts_z_xxxxx, ts_z_xxxxy, ts_z_xxxxz, ts_z_xxxyy, ts_z_xxxyz, ts_z_xxxzz, ts_z_xxyyy, ts_z_xxyyz, ts_z_xxyzz, ts_z_xxzzz, ts_z_xyyyy, ts_z_xyyyz, ts_z_xyyzz, ts_z_xyzzz, ts_z_xzzzz, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyzz, ts_z_yyzzz, ts_z_yzzzz, ts_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xz_xxxxx[i] = 2.0 * ts_z_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xz_xxxxy[i] = 2.0 * ts_z_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xz_xxxxz[i] = 2.0 * ts_z_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xz_xxxyy[i] = 2.0 * ts_z_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xz_xxxyz[i] = 2.0 * ts_z_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xz_xxxzz[i] = 2.0 * ts_z_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xz_xxyyy[i] = 2.0 * ts_z_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xz_xxyyz[i] = 2.0 * ts_z_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xz_xxyzz[i] = 2.0 * ts_z_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xz_xxzzz[i] = 2.0 * ts_z_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xz_xyyyy[i] = 2.0 * ts_z_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xz_xyyyz[i] = 2.0 * ts_z_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xz_xyyzz[i] = 2.0 * ts_z_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xz_xyzzz[i] = 2.0 * ts_z_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xz_xzzzz[i] = 2.0 * ts_z_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xz_yyyyy[i] = 2.0 * ts_z_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xz_yyyyz[i] = 2.0 * ts_z_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xz_yyyzz[i] = 2.0 * ts_z_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xz_yyzzz[i] = 2.0 * ts_z_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xz_yzzzz[i] = 2.0 * ts_z_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xz_zzzzz[i] = 2.0 * ts_z_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 63-84 components of targeted buffer : DH

    auto gs_x_yy_xxxxx = pbuffer.data(idx_g_dh + 63);

    auto gs_x_yy_xxxxy = pbuffer.data(idx_g_dh + 64);

    auto gs_x_yy_xxxxz = pbuffer.data(idx_g_dh + 65);

    auto gs_x_yy_xxxyy = pbuffer.data(idx_g_dh + 66);

    auto gs_x_yy_xxxyz = pbuffer.data(idx_g_dh + 67);

    auto gs_x_yy_xxxzz = pbuffer.data(idx_g_dh + 68);

    auto gs_x_yy_xxyyy = pbuffer.data(idx_g_dh + 69);

    auto gs_x_yy_xxyyz = pbuffer.data(idx_g_dh + 70);

    auto gs_x_yy_xxyzz = pbuffer.data(idx_g_dh + 71);

    auto gs_x_yy_xxzzz = pbuffer.data(idx_g_dh + 72);

    auto gs_x_yy_xyyyy = pbuffer.data(idx_g_dh + 73);

    auto gs_x_yy_xyyyz = pbuffer.data(idx_g_dh + 74);

    auto gs_x_yy_xyyzz = pbuffer.data(idx_g_dh + 75);

    auto gs_x_yy_xyzzz = pbuffer.data(idx_g_dh + 76);

    auto gs_x_yy_xzzzz = pbuffer.data(idx_g_dh + 77);

    auto gs_x_yy_yyyyy = pbuffer.data(idx_g_dh + 78);

    auto gs_x_yy_yyyyz = pbuffer.data(idx_g_dh + 79);

    auto gs_x_yy_yyyzz = pbuffer.data(idx_g_dh + 80);

    auto gs_x_yy_yyzzz = pbuffer.data(idx_g_dh + 81);

    auto gs_x_yy_yzzzz = pbuffer.data(idx_g_dh + 82);

    auto gs_x_yy_zzzzz = pbuffer.data(idx_g_dh + 83);

    #pragma omp simd aligned(gc_x, gs_x_yy_xxxxx, gs_x_yy_xxxxy, gs_x_yy_xxxxz, gs_x_yy_xxxyy, gs_x_yy_xxxyz, gs_x_yy_xxxzz, gs_x_yy_xxyyy, gs_x_yy_xxyyz, gs_x_yy_xxyzz, gs_x_yy_xxzzz, gs_x_yy_xyyyy, gs_x_yy_xyyyz, gs_x_yy_xyyzz, gs_x_yy_xyzzz, gs_x_yy_xzzzz, gs_x_yy_yyyyy, gs_x_yy_yyyyz, gs_x_yy_yyyzz, gs_x_yy_yyzzz, gs_x_yy_yzzzz, gs_x_yy_zzzzz, ts_yy_xxxx, ts_yy_xxxxx, ts_yy_xxxxy, ts_yy_xxxxz, ts_yy_xxxy, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxxz, ts_yy_xxxzz, ts_yy_xxyy, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyz, ts_yy_xxyzz, ts_yy_xxzz, ts_yy_xxzzz, ts_yy_xyyy, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyz, ts_yy_xyyzz, ts_yy_xyzz, ts_yy_xyzzz, ts_yy_xzzz, ts_yy_xzzzz, ts_yy_yyyy, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyz, ts_yy_yyyzz, ts_yy_yyzz, ts_yy_yyzzz, ts_yy_yzzz, ts_yy_yzzzz, ts_yy_zzzz, ts_yy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yy_xxxxx[i] = 10.0 * ts_yy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yy_xxxxy[i] = 8.0 * ts_yy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yy_xxxxz[i] = 8.0 * ts_yy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yy_xxxyy[i] = 6.0 * ts_yy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yy_xxxyz[i] = 6.0 * ts_yy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yy_xxxzz[i] = 6.0 * ts_yy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yy_xxyyy[i] = 4.0 * ts_yy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yy_xxyyz[i] = 4.0 * ts_yy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yy_xxyzz[i] = 4.0 * ts_yy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yy_xxzzz[i] = 4.0 * ts_yy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yy_xyyyy[i] = 2.0 * ts_yy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yy_xyyyz[i] = 2.0 * ts_yy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yy_xyyzz[i] = 2.0 * ts_yy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yy_xyzzz[i] = 2.0 * ts_yy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yy_xzzzz[i] = 2.0 * ts_yy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yy_yyyyy[i] = 2.0 * ts_yy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yy_yyyyz[i] = 2.0 * ts_yy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yy_yyyzz[i] = 2.0 * ts_yy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yy_yyzzz[i] = 2.0 * ts_yy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yy_yzzzz[i] = 2.0 * ts_yy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yy_zzzzz[i] = 2.0 * ts_yy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 84-105 components of targeted buffer : DH

    auto gs_x_yz_xxxxx = pbuffer.data(idx_g_dh + 84);

    auto gs_x_yz_xxxxy = pbuffer.data(idx_g_dh + 85);

    auto gs_x_yz_xxxxz = pbuffer.data(idx_g_dh + 86);

    auto gs_x_yz_xxxyy = pbuffer.data(idx_g_dh + 87);

    auto gs_x_yz_xxxyz = pbuffer.data(idx_g_dh + 88);

    auto gs_x_yz_xxxzz = pbuffer.data(idx_g_dh + 89);

    auto gs_x_yz_xxyyy = pbuffer.data(idx_g_dh + 90);

    auto gs_x_yz_xxyyz = pbuffer.data(idx_g_dh + 91);

    auto gs_x_yz_xxyzz = pbuffer.data(idx_g_dh + 92);

    auto gs_x_yz_xxzzz = pbuffer.data(idx_g_dh + 93);

    auto gs_x_yz_xyyyy = pbuffer.data(idx_g_dh + 94);

    auto gs_x_yz_xyyyz = pbuffer.data(idx_g_dh + 95);

    auto gs_x_yz_xyyzz = pbuffer.data(idx_g_dh + 96);

    auto gs_x_yz_xyzzz = pbuffer.data(idx_g_dh + 97);

    auto gs_x_yz_xzzzz = pbuffer.data(idx_g_dh + 98);

    auto gs_x_yz_yyyyy = pbuffer.data(idx_g_dh + 99);

    auto gs_x_yz_yyyyz = pbuffer.data(idx_g_dh + 100);

    auto gs_x_yz_yyyzz = pbuffer.data(idx_g_dh + 101);

    auto gs_x_yz_yyzzz = pbuffer.data(idx_g_dh + 102);

    auto gs_x_yz_yzzzz = pbuffer.data(idx_g_dh + 103);

    auto gs_x_yz_zzzzz = pbuffer.data(idx_g_dh + 104);

    #pragma omp simd aligned(gc_x, gs_x_yz_xxxxx, gs_x_yz_xxxxy, gs_x_yz_xxxxz, gs_x_yz_xxxyy, gs_x_yz_xxxyz, gs_x_yz_xxxzz, gs_x_yz_xxyyy, gs_x_yz_xxyyz, gs_x_yz_xxyzz, gs_x_yz_xxzzz, gs_x_yz_xyyyy, gs_x_yz_xyyyz, gs_x_yz_xyyzz, gs_x_yz_xyzzz, gs_x_yz_xzzzz, gs_x_yz_yyyyy, gs_x_yz_yyyyz, gs_x_yz_yyyzz, gs_x_yz_yyzzz, gs_x_yz_yzzzz, gs_x_yz_zzzzz, ts_yz_xxxx, ts_yz_xxxxx, ts_yz_xxxxy, ts_yz_xxxxz, ts_yz_xxxy, ts_yz_xxxyy, ts_yz_xxxyz, ts_yz_xxxz, ts_yz_xxxzz, ts_yz_xxyy, ts_yz_xxyyy, ts_yz_xxyyz, ts_yz_xxyz, ts_yz_xxyzz, ts_yz_xxzz, ts_yz_xxzzz, ts_yz_xyyy, ts_yz_xyyyy, ts_yz_xyyyz, ts_yz_xyyz, ts_yz_xyyzz, ts_yz_xyzz, ts_yz_xyzzz, ts_yz_xzzz, ts_yz_xzzzz, ts_yz_yyyy, ts_yz_yyyyy, ts_yz_yyyyz, ts_yz_yyyz, ts_yz_yyyzz, ts_yz_yyzz, ts_yz_yyzzz, ts_yz_yzzz, ts_yz_yzzzz, ts_yz_zzzz, ts_yz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yz_xxxxx[i] = 10.0 * ts_yz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yz_xxxxy[i] = 8.0 * ts_yz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yz_xxxxz[i] = 8.0 * ts_yz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yz_xxxyy[i] = 6.0 * ts_yz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yz_xxxyz[i] = 6.0 * ts_yz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yz_xxxzz[i] = 6.0 * ts_yz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yz_xxyyy[i] = 4.0 * ts_yz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yz_xxyyz[i] = 4.0 * ts_yz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yz_xxyzz[i] = 4.0 * ts_yz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yz_xxzzz[i] = 4.0 * ts_yz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yz_xyyyy[i] = 2.0 * ts_yz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yz_xyyyz[i] = 2.0 * ts_yz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yz_xyyzz[i] = 2.0 * ts_yz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yz_xyzzz[i] = 2.0 * ts_yz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yz_xzzzz[i] = 2.0 * ts_yz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yz_yyyyy[i] = 2.0 * ts_yz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yz_yyyyz[i] = 2.0 * ts_yz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yz_yyyzz[i] = 2.0 * ts_yz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yz_yyzzz[i] = 2.0 * ts_yz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yz_yzzzz[i] = 2.0 * ts_yz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yz_zzzzz[i] = 2.0 * ts_yz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 105-126 components of targeted buffer : DH

    auto gs_x_zz_xxxxx = pbuffer.data(idx_g_dh + 105);

    auto gs_x_zz_xxxxy = pbuffer.data(idx_g_dh + 106);

    auto gs_x_zz_xxxxz = pbuffer.data(idx_g_dh + 107);

    auto gs_x_zz_xxxyy = pbuffer.data(idx_g_dh + 108);

    auto gs_x_zz_xxxyz = pbuffer.data(idx_g_dh + 109);

    auto gs_x_zz_xxxzz = pbuffer.data(idx_g_dh + 110);

    auto gs_x_zz_xxyyy = pbuffer.data(idx_g_dh + 111);

    auto gs_x_zz_xxyyz = pbuffer.data(idx_g_dh + 112);

    auto gs_x_zz_xxyzz = pbuffer.data(idx_g_dh + 113);

    auto gs_x_zz_xxzzz = pbuffer.data(idx_g_dh + 114);

    auto gs_x_zz_xyyyy = pbuffer.data(idx_g_dh + 115);

    auto gs_x_zz_xyyyz = pbuffer.data(idx_g_dh + 116);

    auto gs_x_zz_xyyzz = pbuffer.data(idx_g_dh + 117);

    auto gs_x_zz_xyzzz = pbuffer.data(idx_g_dh + 118);

    auto gs_x_zz_xzzzz = pbuffer.data(idx_g_dh + 119);

    auto gs_x_zz_yyyyy = pbuffer.data(idx_g_dh + 120);

    auto gs_x_zz_yyyyz = pbuffer.data(idx_g_dh + 121);

    auto gs_x_zz_yyyzz = pbuffer.data(idx_g_dh + 122);

    auto gs_x_zz_yyzzz = pbuffer.data(idx_g_dh + 123);

    auto gs_x_zz_yzzzz = pbuffer.data(idx_g_dh + 124);

    auto gs_x_zz_zzzzz = pbuffer.data(idx_g_dh + 125);

    #pragma omp simd aligned(gc_x, gs_x_zz_xxxxx, gs_x_zz_xxxxy, gs_x_zz_xxxxz, gs_x_zz_xxxyy, gs_x_zz_xxxyz, gs_x_zz_xxxzz, gs_x_zz_xxyyy, gs_x_zz_xxyyz, gs_x_zz_xxyzz, gs_x_zz_xxzzz, gs_x_zz_xyyyy, gs_x_zz_xyyyz, gs_x_zz_xyyzz, gs_x_zz_xyzzz, gs_x_zz_xzzzz, gs_x_zz_yyyyy, gs_x_zz_yyyyz, gs_x_zz_yyyzz, gs_x_zz_yyzzz, gs_x_zz_yzzzz, gs_x_zz_zzzzz, ts_zz_xxxx, ts_zz_xxxxx, ts_zz_xxxxy, ts_zz_xxxxz, ts_zz_xxxy, ts_zz_xxxyy, ts_zz_xxxyz, ts_zz_xxxz, ts_zz_xxxzz, ts_zz_xxyy, ts_zz_xxyyy, ts_zz_xxyyz, ts_zz_xxyz, ts_zz_xxyzz, ts_zz_xxzz, ts_zz_xxzzz, ts_zz_xyyy, ts_zz_xyyyy, ts_zz_xyyyz, ts_zz_xyyz, ts_zz_xyyzz, ts_zz_xyzz, ts_zz_xyzzz, ts_zz_xzzz, ts_zz_xzzzz, ts_zz_yyyy, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyz, ts_zz_yyyzz, ts_zz_yyzz, ts_zz_yyzzz, ts_zz_yzzz, ts_zz_yzzzz, ts_zz_zzzz, ts_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zz_xxxxx[i] = 10.0 * ts_zz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_zz_xxxxy[i] = 8.0 * ts_zz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_zz_xxxxz[i] = 8.0 * ts_zz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_zz_xxxyy[i] = 6.0 * ts_zz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_zz_xxxyz[i] = 6.0 * ts_zz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_zz_xxxzz[i] = 6.0 * ts_zz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_zz_xxyyy[i] = 4.0 * ts_zz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_zz_xxyyz[i] = 4.0 * ts_zz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_zz_xxyzz[i] = 4.0 * ts_zz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_zz_xxzzz[i] = 4.0 * ts_zz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_zz_xyyyy[i] = 2.0 * ts_zz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_zz_xyyyz[i] = 2.0 * ts_zz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_zz_xyyzz[i] = 2.0 * ts_zz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_zz_xyzzz[i] = 2.0 * ts_zz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_zz_xzzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_zz_yyyyy[i] = 2.0 * ts_zz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_zz_yyyyz[i] = 2.0 * ts_zz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_zz_yyyzz[i] = 2.0 * ts_zz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_zz_yyzzz[i] = 2.0 * ts_zz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_zz_yzzzz[i] = 2.0 * ts_zz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_zz_zzzzz[i] = 2.0 * ts_zz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 126-147 components of targeted buffer : DH

    auto gs_y_xx_xxxxx = pbuffer.data(idx_g_dh + 126);

    auto gs_y_xx_xxxxy = pbuffer.data(idx_g_dh + 127);

    auto gs_y_xx_xxxxz = pbuffer.data(idx_g_dh + 128);

    auto gs_y_xx_xxxyy = pbuffer.data(idx_g_dh + 129);

    auto gs_y_xx_xxxyz = pbuffer.data(idx_g_dh + 130);

    auto gs_y_xx_xxxzz = pbuffer.data(idx_g_dh + 131);

    auto gs_y_xx_xxyyy = pbuffer.data(idx_g_dh + 132);

    auto gs_y_xx_xxyyz = pbuffer.data(idx_g_dh + 133);

    auto gs_y_xx_xxyzz = pbuffer.data(idx_g_dh + 134);

    auto gs_y_xx_xxzzz = pbuffer.data(idx_g_dh + 135);

    auto gs_y_xx_xyyyy = pbuffer.data(idx_g_dh + 136);

    auto gs_y_xx_xyyyz = pbuffer.data(idx_g_dh + 137);

    auto gs_y_xx_xyyzz = pbuffer.data(idx_g_dh + 138);

    auto gs_y_xx_xyzzz = pbuffer.data(idx_g_dh + 139);

    auto gs_y_xx_xzzzz = pbuffer.data(idx_g_dh + 140);

    auto gs_y_xx_yyyyy = pbuffer.data(idx_g_dh + 141);

    auto gs_y_xx_yyyyz = pbuffer.data(idx_g_dh + 142);

    auto gs_y_xx_yyyzz = pbuffer.data(idx_g_dh + 143);

    auto gs_y_xx_yyzzz = pbuffer.data(idx_g_dh + 144);

    auto gs_y_xx_yzzzz = pbuffer.data(idx_g_dh + 145);

    auto gs_y_xx_zzzzz = pbuffer.data(idx_g_dh + 146);

    #pragma omp simd aligned(gc_y, gs_y_xx_xxxxx, gs_y_xx_xxxxy, gs_y_xx_xxxxz, gs_y_xx_xxxyy, gs_y_xx_xxxyz, gs_y_xx_xxxzz, gs_y_xx_xxyyy, gs_y_xx_xxyyz, gs_y_xx_xxyzz, gs_y_xx_xxzzz, gs_y_xx_xyyyy, gs_y_xx_xyyyz, gs_y_xx_xyyzz, gs_y_xx_xyzzz, gs_y_xx_xzzzz, gs_y_xx_yyyyy, gs_y_xx_yyyyz, gs_y_xx_yyyzz, gs_y_xx_yyzzz, gs_y_xx_yzzzz, gs_y_xx_zzzzz, ts_xx_xxxx, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxxz, ts_xx_xxxy, ts_xx_xxxyy, ts_xx_xxxyz, ts_xx_xxxz, ts_xx_xxxzz, ts_xx_xxyy, ts_xx_xxyyy, ts_xx_xxyyz, ts_xx_xxyz, ts_xx_xxyzz, ts_xx_xxzz, ts_xx_xxzzz, ts_xx_xyyy, ts_xx_xyyyy, ts_xx_xyyyz, ts_xx_xyyz, ts_xx_xyyzz, ts_xx_xyzz, ts_xx_xyzzz, ts_xx_xzzz, ts_xx_xzzzz, ts_xx_yyyy, ts_xx_yyyyy, ts_xx_yyyyz, ts_xx_yyyz, ts_xx_yyyzz, ts_xx_yyzz, ts_xx_yyzzz, ts_xx_yzzz, ts_xx_yzzzz, ts_xx_zzzz, ts_xx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xx_xxxxx[i] = 2.0 * ts_xx_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xx_xxxxy[i] = 2.0 * ts_xx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xx_xxxxz[i] = 2.0 * ts_xx_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xx_xxxyy[i] = 4.0 * ts_xx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xx_xxxyz[i] = 2.0 * ts_xx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xx_xxxzz[i] = 2.0 * ts_xx_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xx_xxyyy[i] = 6.0 * ts_xx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xx_xxyyz[i] = 4.0 * ts_xx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xx_xxyzz[i] = 2.0 * ts_xx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xx_xxzzz[i] = 2.0 * ts_xx_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xx_xyyyy[i] = 8.0 * ts_xx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xx_xyyyz[i] = 6.0 * ts_xx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xx_xyyzz[i] = 4.0 * ts_xx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xx_xyzzz[i] = 2.0 * ts_xx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xx_xzzzz[i] = 2.0 * ts_xx_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xx_yyyyy[i] = 10.0 * ts_xx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xx_yyyyz[i] = 8.0 * ts_xx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xx_yyyzz[i] = 6.0 * ts_xx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xx_yyzzz[i] = 4.0 * ts_xx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xx_yzzzz[i] = 2.0 * ts_xx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xx_zzzzz[i] = 2.0 * ts_xx_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 147-168 components of targeted buffer : DH

    auto gs_y_xy_xxxxx = pbuffer.data(idx_g_dh + 147);

    auto gs_y_xy_xxxxy = pbuffer.data(idx_g_dh + 148);

    auto gs_y_xy_xxxxz = pbuffer.data(idx_g_dh + 149);

    auto gs_y_xy_xxxyy = pbuffer.data(idx_g_dh + 150);

    auto gs_y_xy_xxxyz = pbuffer.data(idx_g_dh + 151);

    auto gs_y_xy_xxxzz = pbuffer.data(idx_g_dh + 152);

    auto gs_y_xy_xxyyy = pbuffer.data(idx_g_dh + 153);

    auto gs_y_xy_xxyyz = pbuffer.data(idx_g_dh + 154);

    auto gs_y_xy_xxyzz = pbuffer.data(idx_g_dh + 155);

    auto gs_y_xy_xxzzz = pbuffer.data(idx_g_dh + 156);

    auto gs_y_xy_xyyyy = pbuffer.data(idx_g_dh + 157);

    auto gs_y_xy_xyyyz = pbuffer.data(idx_g_dh + 158);

    auto gs_y_xy_xyyzz = pbuffer.data(idx_g_dh + 159);

    auto gs_y_xy_xyzzz = pbuffer.data(idx_g_dh + 160);

    auto gs_y_xy_xzzzz = pbuffer.data(idx_g_dh + 161);

    auto gs_y_xy_yyyyy = pbuffer.data(idx_g_dh + 162);

    auto gs_y_xy_yyyyz = pbuffer.data(idx_g_dh + 163);

    auto gs_y_xy_yyyzz = pbuffer.data(idx_g_dh + 164);

    auto gs_y_xy_yyzzz = pbuffer.data(idx_g_dh + 165);

    auto gs_y_xy_yzzzz = pbuffer.data(idx_g_dh + 166);

    auto gs_y_xy_zzzzz = pbuffer.data(idx_g_dh + 167);

    #pragma omp simd aligned(gc_y, gs_y_xy_xxxxx, gs_y_xy_xxxxy, gs_y_xy_xxxxz, gs_y_xy_xxxyy, gs_y_xy_xxxyz, gs_y_xy_xxxzz, gs_y_xy_xxyyy, gs_y_xy_xxyyz, gs_y_xy_xxyzz, gs_y_xy_xxzzz, gs_y_xy_xyyyy, gs_y_xy_xyyyz, gs_y_xy_xyyzz, gs_y_xy_xyzzz, gs_y_xy_xzzzz, gs_y_xy_yyyyy, gs_y_xy_yyyyz, gs_y_xy_yyyzz, gs_y_xy_yyzzz, gs_y_xy_yzzzz, gs_y_xy_zzzzz, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxxz, ts_x_xxxyy, ts_x_xxxyz, ts_x_xxxzz, ts_x_xxyyy, ts_x_xxyyz, ts_x_xxyzz, ts_x_xxzzz, ts_x_xyyyy, ts_x_xyyyz, ts_x_xyyzz, ts_x_xyzzz, ts_x_xzzzz, ts_x_yyyyy, ts_x_yyyyz, ts_x_yyyzz, ts_x_yyzzz, ts_x_yzzzz, ts_x_zzzzz, ts_xy_xxxx, ts_xy_xxxxx, ts_xy_xxxxy, ts_xy_xxxxz, ts_xy_xxxy, ts_xy_xxxyy, ts_xy_xxxyz, ts_xy_xxxz, ts_xy_xxxzz, ts_xy_xxyy, ts_xy_xxyyy, ts_xy_xxyyz, ts_xy_xxyz, ts_xy_xxyzz, ts_xy_xxzz, ts_xy_xxzzz, ts_xy_xyyy, ts_xy_xyyyy, ts_xy_xyyyz, ts_xy_xyyz, ts_xy_xyyzz, ts_xy_xyzz, ts_xy_xyzzz, ts_xy_xzzz, ts_xy_xzzzz, ts_xy_yyyy, ts_xy_yyyyy, ts_xy_yyyyz, ts_xy_yyyz, ts_xy_yyyzz, ts_xy_yyzz, ts_xy_yyzzz, ts_xy_yzzz, ts_xy_yzzzz, ts_xy_zzzz, ts_xy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xy_xxxxx[i] = 2.0 * ts_x_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xy_xxxxy[i] = 2.0 * ts_x_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xy_xxxxz[i] = 2.0 * ts_x_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xy_xxxyy[i] = 2.0 * ts_x_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xy_xxxyz[i] = 2.0 * ts_x_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xy_xxxzz[i] = 2.0 * ts_x_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xy_xxyyy[i] = 2.0 * ts_x_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xy_xxyyz[i] = 2.0 * ts_x_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xy_xxyzz[i] = 2.0 * ts_x_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xy_xxzzz[i] = 2.0 * ts_x_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xy_xyyyy[i] = 2.0 * ts_x_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xy_xyyyz[i] = 2.0 * ts_x_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xy_xyyzz[i] = 2.0 * ts_x_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xy_xyzzz[i] = 2.0 * ts_x_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xy_xzzzz[i] = 2.0 * ts_x_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xy_yyyyy[i] = 2.0 * ts_x_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xy_yyyyz[i] = 2.0 * ts_x_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xy_yyyzz[i] = 2.0 * ts_x_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xy_yyzzz[i] = 2.0 * ts_x_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xy_yzzzz[i] = 2.0 * ts_x_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xy_zzzzz[i] = 2.0 * ts_x_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 168-189 components of targeted buffer : DH

    auto gs_y_xz_xxxxx = pbuffer.data(idx_g_dh + 168);

    auto gs_y_xz_xxxxy = pbuffer.data(idx_g_dh + 169);

    auto gs_y_xz_xxxxz = pbuffer.data(idx_g_dh + 170);

    auto gs_y_xz_xxxyy = pbuffer.data(idx_g_dh + 171);

    auto gs_y_xz_xxxyz = pbuffer.data(idx_g_dh + 172);

    auto gs_y_xz_xxxzz = pbuffer.data(idx_g_dh + 173);

    auto gs_y_xz_xxyyy = pbuffer.data(idx_g_dh + 174);

    auto gs_y_xz_xxyyz = pbuffer.data(idx_g_dh + 175);

    auto gs_y_xz_xxyzz = pbuffer.data(idx_g_dh + 176);

    auto gs_y_xz_xxzzz = pbuffer.data(idx_g_dh + 177);

    auto gs_y_xz_xyyyy = pbuffer.data(idx_g_dh + 178);

    auto gs_y_xz_xyyyz = pbuffer.data(idx_g_dh + 179);

    auto gs_y_xz_xyyzz = pbuffer.data(idx_g_dh + 180);

    auto gs_y_xz_xyzzz = pbuffer.data(idx_g_dh + 181);

    auto gs_y_xz_xzzzz = pbuffer.data(idx_g_dh + 182);

    auto gs_y_xz_yyyyy = pbuffer.data(idx_g_dh + 183);

    auto gs_y_xz_yyyyz = pbuffer.data(idx_g_dh + 184);

    auto gs_y_xz_yyyzz = pbuffer.data(idx_g_dh + 185);

    auto gs_y_xz_yyzzz = pbuffer.data(idx_g_dh + 186);

    auto gs_y_xz_yzzzz = pbuffer.data(idx_g_dh + 187);

    auto gs_y_xz_zzzzz = pbuffer.data(idx_g_dh + 188);

    #pragma omp simd aligned(gc_y, gs_y_xz_xxxxx, gs_y_xz_xxxxy, gs_y_xz_xxxxz, gs_y_xz_xxxyy, gs_y_xz_xxxyz, gs_y_xz_xxxzz, gs_y_xz_xxyyy, gs_y_xz_xxyyz, gs_y_xz_xxyzz, gs_y_xz_xxzzz, gs_y_xz_xyyyy, gs_y_xz_xyyyz, gs_y_xz_xyyzz, gs_y_xz_xyzzz, gs_y_xz_xzzzz, gs_y_xz_yyyyy, gs_y_xz_yyyyz, gs_y_xz_yyyzz, gs_y_xz_yyzzz, gs_y_xz_yzzzz, gs_y_xz_zzzzz, ts_xz_xxxx, ts_xz_xxxxx, ts_xz_xxxxy, ts_xz_xxxxz, ts_xz_xxxy, ts_xz_xxxyy, ts_xz_xxxyz, ts_xz_xxxz, ts_xz_xxxzz, ts_xz_xxyy, ts_xz_xxyyy, ts_xz_xxyyz, ts_xz_xxyz, ts_xz_xxyzz, ts_xz_xxzz, ts_xz_xxzzz, ts_xz_xyyy, ts_xz_xyyyy, ts_xz_xyyyz, ts_xz_xyyz, ts_xz_xyyzz, ts_xz_xyzz, ts_xz_xyzzz, ts_xz_xzzz, ts_xz_xzzzz, ts_xz_yyyy, ts_xz_yyyyy, ts_xz_yyyyz, ts_xz_yyyz, ts_xz_yyyzz, ts_xz_yyzz, ts_xz_yyzzz, ts_xz_yzzz, ts_xz_yzzzz, ts_xz_zzzz, ts_xz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xz_xxxxx[i] = 2.0 * ts_xz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xz_xxxxy[i] = 2.0 * ts_xz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xz_xxxxz[i] = 2.0 * ts_xz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xz_xxxyy[i] = 4.0 * ts_xz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xz_xxxyz[i] = 2.0 * ts_xz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xz_xxxzz[i] = 2.0 * ts_xz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xz_xxyyy[i] = 6.0 * ts_xz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xz_xxyyz[i] = 4.0 * ts_xz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xz_xxyzz[i] = 2.0 * ts_xz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xz_xxzzz[i] = 2.0 * ts_xz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xz_xyyyy[i] = 8.0 * ts_xz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xz_xyyyz[i] = 6.0 * ts_xz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xz_xyyzz[i] = 4.0 * ts_xz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xz_xyzzz[i] = 2.0 * ts_xz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xz_xzzzz[i] = 2.0 * ts_xz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xz_yyyyy[i] = 10.0 * ts_xz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xz_yyyyz[i] = 8.0 * ts_xz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xz_yyyzz[i] = 6.0 * ts_xz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xz_yyzzz[i] = 4.0 * ts_xz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xz_yzzzz[i] = 2.0 * ts_xz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xz_zzzzz[i] = 2.0 * ts_xz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 189-210 components of targeted buffer : DH

    auto gs_y_yy_xxxxx = pbuffer.data(idx_g_dh + 189);

    auto gs_y_yy_xxxxy = pbuffer.data(idx_g_dh + 190);

    auto gs_y_yy_xxxxz = pbuffer.data(idx_g_dh + 191);

    auto gs_y_yy_xxxyy = pbuffer.data(idx_g_dh + 192);

    auto gs_y_yy_xxxyz = pbuffer.data(idx_g_dh + 193);

    auto gs_y_yy_xxxzz = pbuffer.data(idx_g_dh + 194);

    auto gs_y_yy_xxyyy = pbuffer.data(idx_g_dh + 195);

    auto gs_y_yy_xxyyz = pbuffer.data(idx_g_dh + 196);

    auto gs_y_yy_xxyzz = pbuffer.data(idx_g_dh + 197);

    auto gs_y_yy_xxzzz = pbuffer.data(idx_g_dh + 198);

    auto gs_y_yy_xyyyy = pbuffer.data(idx_g_dh + 199);

    auto gs_y_yy_xyyyz = pbuffer.data(idx_g_dh + 200);

    auto gs_y_yy_xyyzz = pbuffer.data(idx_g_dh + 201);

    auto gs_y_yy_xyzzz = pbuffer.data(idx_g_dh + 202);

    auto gs_y_yy_xzzzz = pbuffer.data(idx_g_dh + 203);

    auto gs_y_yy_yyyyy = pbuffer.data(idx_g_dh + 204);

    auto gs_y_yy_yyyyz = pbuffer.data(idx_g_dh + 205);

    auto gs_y_yy_yyyzz = pbuffer.data(idx_g_dh + 206);

    auto gs_y_yy_yyzzz = pbuffer.data(idx_g_dh + 207);

    auto gs_y_yy_yzzzz = pbuffer.data(idx_g_dh + 208);

    auto gs_y_yy_zzzzz = pbuffer.data(idx_g_dh + 209);

    #pragma omp simd aligned(gc_y, gs_y_yy_xxxxx, gs_y_yy_xxxxy, gs_y_yy_xxxxz, gs_y_yy_xxxyy, gs_y_yy_xxxyz, gs_y_yy_xxxzz, gs_y_yy_xxyyy, gs_y_yy_xxyyz, gs_y_yy_xxyzz, gs_y_yy_xxzzz, gs_y_yy_xyyyy, gs_y_yy_xyyyz, gs_y_yy_xyyzz, gs_y_yy_xyzzz, gs_y_yy_xzzzz, gs_y_yy_yyyyy, gs_y_yy_yyyyz, gs_y_yy_yyyzz, gs_y_yy_yyzzz, gs_y_yy_yzzzz, gs_y_yy_zzzzz, ts_y_xxxxx, ts_y_xxxxy, ts_y_xxxxz, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxxzz, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyzz, ts_y_xxzzz, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyzz, ts_y_xyzzz, ts_y_xzzzz, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyzz, ts_y_yyzzz, ts_y_yzzzz, ts_y_zzzzz, ts_yy_xxxx, ts_yy_xxxxx, ts_yy_xxxxy, ts_yy_xxxxz, ts_yy_xxxy, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxxz, ts_yy_xxxzz, ts_yy_xxyy, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyz, ts_yy_xxyzz, ts_yy_xxzz, ts_yy_xxzzz, ts_yy_xyyy, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyz, ts_yy_xyyzz, ts_yy_xyzz, ts_yy_xyzzz, ts_yy_xzzz, ts_yy_xzzzz, ts_yy_yyyy, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyz, ts_yy_yyyzz, ts_yy_yyzz, ts_yy_yyzzz, ts_yy_yzzz, ts_yy_yzzzz, ts_yy_zzzz, ts_yy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yy_xxxxx[i] = 4.0 * ts_y_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yy_xxxxy[i] = 4.0 * ts_y_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yy_xxxxz[i] = 4.0 * ts_y_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yy_xxxyy[i] = 4.0 * ts_y_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yy_xxxyz[i] = 4.0 * ts_y_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yy_xxxzz[i] = 4.0 * ts_y_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yy_xxyyy[i] = 4.0 * ts_y_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yy_xxyyz[i] = 4.0 * ts_y_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yy_xxyzz[i] = 4.0 * ts_y_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yy_xxzzz[i] = 4.0 * ts_y_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yy_xyyyy[i] = 4.0 * ts_y_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yy_xyyyz[i] = 4.0 * ts_y_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yy_xyyzz[i] = 4.0 * ts_y_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yy_xyzzz[i] = 4.0 * ts_y_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yy_xzzzz[i] = 4.0 * ts_y_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yy_yyyyy[i] = 4.0 * ts_y_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yy_yyyyz[i] = 4.0 * ts_y_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yy_yyyzz[i] = 4.0 * ts_y_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yy_yyzzz[i] = 4.0 * ts_y_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yy_yzzzz[i] = 4.0 * ts_y_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yy_zzzzz[i] = 4.0 * ts_y_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 210-231 components of targeted buffer : DH

    auto gs_y_yz_xxxxx = pbuffer.data(idx_g_dh + 210);

    auto gs_y_yz_xxxxy = pbuffer.data(idx_g_dh + 211);

    auto gs_y_yz_xxxxz = pbuffer.data(idx_g_dh + 212);

    auto gs_y_yz_xxxyy = pbuffer.data(idx_g_dh + 213);

    auto gs_y_yz_xxxyz = pbuffer.data(idx_g_dh + 214);

    auto gs_y_yz_xxxzz = pbuffer.data(idx_g_dh + 215);

    auto gs_y_yz_xxyyy = pbuffer.data(idx_g_dh + 216);

    auto gs_y_yz_xxyyz = pbuffer.data(idx_g_dh + 217);

    auto gs_y_yz_xxyzz = pbuffer.data(idx_g_dh + 218);

    auto gs_y_yz_xxzzz = pbuffer.data(idx_g_dh + 219);

    auto gs_y_yz_xyyyy = pbuffer.data(idx_g_dh + 220);

    auto gs_y_yz_xyyyz = pbuffer.data(idx_g_dh + 221);

    auto gs_y_yz_xyyzz = pbuffer.data(idx_g_dh + 222);

    auto gs_y_yz_xyzzz = pbuffer.data(idx_g_dh + 223);

    auto gs_y_yz_xzzzz = pbuffer.data(idx_g_dh + 224);

    auto gs_y_yz_yyyyy = pbuffer.data(idx_g_dh + 225);

    auto gs_y_yz_yyyyz = pbuffer.data(idx_g_dh + 226);

    auto gs_y_yz_yyyzz = pbuffer.data(idx_g_dh + 227);

    auto gs_y_yz_yyzzz = pbuffer.data(idx_g_dh + 228);

    auto gs_y_yz_yzzzz = pbuffer.data(idx_g_dh + 229);

    auto gs_y_yz_zzzzz = pbuffer.data(idx_g_dh + 230);

    #pragma omp simd aligned(gc_y, gs_y_yz_xxxxx, gs_y_yz_xxxxy, gs_y_yz_xxxxz, gs_y_yz_xxxyy, gs_y_yz_xxxyz, gs_y_yz_xxxzz, gs_y_yz_xxyyy, gs_y_yz_xxyyz, gs_y_yz_xxyzz, gs_y_yz_xxzzz, gs_y_yz_xyyyy, gs_y_yz_xyyyz, gs_y_yz_xyyzz, gs_y_yz_xyzzz, gs_y_yz_xzzzz, gs_y_yz_yyyyy, gs_y_yz_yyyyz, gs_y_yz_yyyzz, gs_y_yz_yyzzz, gs_y_yz_yzzzz, gs_y_yz_zzzzz, ts_yz_xxxx, ts_yz_xxxxx, ts_yz_xxxxy, ts_yz_xxxxz, ts_yz_xxxy, ts_yz_xxxyy, ts_yz_xxxyz, ts_yz_xxxz, ts_yz_xxxzz, ts_yz_xxyy, ts_yz_xxyyy, ts_yz_xxyyz, ts_yz_xxyz, ts_yz_xxyzz, ts_yz_xxzz, ts_yz_xxzzz, ts_yz_xyyy, ts_yz_xyyyy, ts_yz_xyyyz, ts_yz_xyyz, ts_yz_xyyzz, ts_yz_xyzz, ts_yz_xyzzz, ts_yz_xzzz, ts_yz_xzzzz, ts_yz_yyyy, ts_yz_yyyyy, ts_yz_yyyyz, ts_yz_yyyz, ts_yz_yyyzz, ts_yz_yyzz, ts_yz_yyzzz, ts_yz_yzzz, ts_yz_yzzzz, ts_yz_zzzz, ts_yz_zzzzz, ts_z_xxxxx, ts_z_xxxxy, ts_z_xxxxz, ts_z_xxxyy, ts_z_xxxyz, ts_z_xxxzz, ts_z_xxyyy, ts_z_xxyyz, ts_z_xxyzz, ts_z_xxzzz, ts_z_xyyyy, ts_z_xyyyz, ts_z_xyyzz, ts_z_xyzzz, ts_z_xzzzz, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyzz, ts_z_yyzzz, ts_z_yzzzz, ts_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yz_xxxxx[i] = 2.0 * ts_z_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yz_xxxxy[i] = 2.0 * ts_z_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yz_xxxxz[i] = 2.0 * ts_z_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yz_xxxyy[i] = 2.0 * ts_z_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yz_xxxyz[i] = 2.0 * ts_z_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yz_xxxzz[i] = 2.0 * ts_z_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yz_xxyyy[i] = 2.0 * ts_z_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yz_xxyyz[i] = 2.0 * ts_z_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yz_xxyzz[i] = 2.0 * ts_z_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yz_xxzzz[i] = 2.0 * ts_z_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yz_xyyyy[i] = 2.0 * ts_z_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yz_xyyyz[i] = 2.0 * ts_z_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yz_xyyzz[i] = 2.0 * ts_z_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yz_xyzzz[i] = 2.0 * ts_z_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yz_xzzzz[i] = 2.0 * ts_z_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yz_yyyyy[i] = 2.0 * ts_z_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yz_yyyyz[i] = 2.0 * ts_z_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yz_yyyzz[i] = 2.0 * ts_z_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yz_yyzzz[i] = 2.0 * ts_z_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yz_yzzzz[i] = 2.0 * ts_z_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yz_zzzzz[i] = 2.0 * ts_z_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 231-252 components of targeted buffer : DH

    auto gs_y_zz_xxxxx = pbuffer.data(idx_g_dh + 231);

    auto gs_y_zz_xxxxy = pbuffer.data(idx_g_dh + 232);

    auto gs_y_zz_xxxxz = pbuffer.data(idx_g_dh + 233);

    auto gs_y_zz_xxxyy = pbuffer.data(idx_g_dh + 234);

    auto gs_y_zz_xxxyz = pbuffer.data(idx_g_dh + 235);

    auto gs_y_zz_xxxzz = pbuffer.data(idx_g_dh + 236);

    auto gs_y_zz_xxyyy = pbuffer.data(idx_g_dh + 237);

    auto gs_y_zz_xxyyz = pbuffer.data(idx_g_dh + 238);

    auto gs_y_zz_xxyzz = pbuffer.data(idx_g_dh + 239);

    auto gs_y_zz_xxzzz = pbuffer.data(idx_g_dh + 240);

    auto gs_y_zz_xyyyy = pbuffer.data(idx_g_dh + 241);

    auto gs_y_zz_xyyyz = pbuffer.data(idx_g_dh + 242);

    auto gs_y_zz_xyyzz = pbuffer.data(idx_g_dh + 243);

    auto gs_y_zz_xyzzz = pbuffer.data(idx_g_dh + 244);

    auto gs_y_zz_xzzzz = pbuffer.data(idx_g_dh + 245);

    auto gs_y_zz_yyyyy = pbuffer.data(idx_g_dh + 246);

    auto gs_y_zz_yyyyz = pbuffer.data(idx_g_dh + 247);

    auto gs_y_zz_yyyzz = pbuffer.data(idx_g_dh + 248);

    auto gs_y_zz_yyzzz = pbuffer.data(idx_g_dh + 249);

    auto gs_y_zz_yzzzz = pbuffer.data(idx_g_dh + 250);

    auto gs_y_zz_zzzzz = pbuffer.data(idx_g_dh + 251);

    #pragma omp simd aligned(gc_y, gs_y_zz_xxxxx, gs_y_zz_xxxxy, gs_y_zz_xxxxz, gs_y_zz_xxxyy, gs_y_zz_xxxyz, gs_y_zz_xxxzz, gs_y_zz_xxyyy, gs_y_zz_xxyyz, gs_y_zz_xxyzz, gs_y_zz_xxzzz, gs_y_zz_xyyyy, gs_y_zz_xyyyz, gs_y_zz_xyyzz, gs_y_zz_xyzzz, gs_y_zz_xzzzz, gs_y_zz_yyyyy, gs_y_zz_yyyyz, gs_y_zz_yyyzz, gs_y_zz_yyzzz, gs_y_zz_yzzzz, gs_y_zz_zzzzz, ts_zz_xxxx, ts_zz_xxxxx, ts_zz_xxxxy, ts_zz_xxxxz, ts_zz_xxxy, ts_zz_xxxyy, ts_zz_xxxyz, ts_zz_xxxz, ts_zz_xxxzz, ts_zz_xxyy, ts_zz_xxyyy, ts_zz_xxyyz, ts_zz_xxyz, ts_zz_xxyzz, ts_zz_xxzz, ts_zz_xxzzz, ts_zz_xyyy, ts_zz_xyyyy, ts_zz_xyyyz, ts_zz_xyyz, ts_zz_xyyzz, ts_zz_xyzz, ts_zz_xyzzz, ts_zz_xzzz, ts_zz_xzzzz, ts_zz_yyyy, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyz, ts_zz_yyyzz, ts_zz_yyzz, ts_zz_yyzzz, ts_zz_yzzz, ts_zz_yzzzz, ts_zz_zzzz, ts_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zz_xxxxx[i] = 2.0 * ts_zz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_zz_xxxxy[i] = 2.0 * ts_zz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_zz_xxxxz[i] = 2.0 * ts_zz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_zz_xxxyy[i] = 4.0 * ts_zz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_zz_xxxyz[i] = 2.0 * ts_zz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_zz_xxxzz[i] = 2.0 * ts_zz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_zz_xxyyy[i] = 6.0 * ts_zz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_zz_xxyyz[i] = 4.0 * ts_zz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_zz_xxyzz[i] = 2.0 * ts_zz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_zz_xxzzz[i] = 2.0 * ts_zz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_zz_xyyyy[i] = 8.0 * ts_zz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_zz_xyyyz[i] = 6.0 * ts_zz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_zz_xyyzz[i] = 4.0 * ts_zz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_zz_xyzzz[i] = 2.0 * ts_zz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_zz_xzzzz[i] = 2.0 * ts_zz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_zz_yyyyy[i] = 10.0 * ts_zz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_zz_yyyyz[i] = 8.0 * ts_zz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_zz_yyyzz[i] = 6.0 * ts_zz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_zz_yyzzz[i] = 4.0 * ts_zz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_zz_yzzzz[i] = 2.0 * ts_zz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_zz_zzzzz[i] = 2.0 * ts_zz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 252-273 components of targeted buffer : DH

    auto gs_z_xx_xxxxx = pbuffer.data(idx_g_dh + 252);

    auto gs_z_xx_xxxxy = pbuffer.data(idx_g_dh + 253);

    auto gs_z_xx_xxxxz = pbuffer.data(idx_g_dh + 254);

    auto gs_z_xx_xxxyy = pbuffer.data(idx_g_dh + 255);

    auto gs_z_xx_xxxyz = pbuffer.data(idx_g_dh + 256);

    auto gs_z_xx_xxxzz = pbuffer.data(idx_g_dh + 257);

    auto gs_z_xx_xxyyy = pbuffer.data(idx_g_dh + 258);

    auto gs_z_xx_xxyyz = pbuffer.data(idx_g_dh + 259);

    auto gs_z_xx_xxyzz = pbuffer.data(idx_g_dh + 260);

    auto gs_z_xx_xxzzz = pbuffer.data(idx_g_dh + 261);

    auto gs_z_xx_xyyyy = pbuffer.data(idx_g_dh + 262);

    auto gs_z_xx_xyyyz = pbuffer.data(idx_g_dh + 263);

    auto gs_z_xx_xyyzz = pbuffer.data(idx_g_dh + 264);

    auto gs_z_xx_xyzzz = pbuffer.data(idx_g_dh + 265);

    auto gs_z_xx_xzzzz = pbuffer.data(idx_g_dh + 266);

    auto gs_z_xx_yyyyy = pbuffer.data(idx_g_dh + 267);

    auto gs_z_xx_yyyyz = pbuffer.data(idx_g_dh + 268);

    auto gs_z_xx_yyyzz = pbuffer.data(idx_g_dh + 269);

    auto gs_z_xx_yyzzz = pbuffer.data(idx_g_dh + 270);

    auto gs_z_xx_yzzzz = pbuffer.data(idx_g_dh + 271);

    auto gs_z_xx_zzzzz = pbuffer.data(idx_g_dh + 272);

    #pragma omp simd aligned(gc_z, gs_z_xx_xxxxx, gs_z_xx_xxxxy, gs_z_xx_xxxxz, gs_z_xx_xxxyy, gs_z_xx_xxxyz, gs_z_xx_xxxzz, gs_z_xx_xxyyy, gs_z_xx_xxyyz, gs_z_xx_xxyzz, gs_z_xx_xxzzz, gs_z_xx_xyyyy, gs_z_xx_xyyyz, gs_z_xx_xyyzz, gs_z_xx_xyzzz, gs_z_xx_xzzzz, gs_z_xx_yyyyy, gs_z_xx_yyyyz, gs_z_xx_yyyzz, gs_z_xx_yyzzz, gs_z_xx_yzzzz, gs_z_xx_zzzzz, ts_xx_xxxx, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxxz, ts_xx_xxxy, ts_xx_xxxyy, ts_xx_xxxyz, ts_xx_xxxz, ts_xx_xxxzz, ts_xx_xxyy, ts_xx_xxyyy, ts_xx_xxyyz, ts_xx_xxyz, ts_xx_xxyzz, ts_xx_xxzz, ts_xx_xxzzz, ts_xx_xyyy, ts_xx_xyyyy, ts_xx_xyyyz, ts_xx_xyyz, ts_xx_xyyzz, ts_xx_xyzz, ts_xx_xyzzz, ts_xx_xzzz, ts_xx_xzzzz, ts_xx_yyyy, ts_xx_yyyyy, ts_xx_yyyyz, ts_xx_yyyz, ts_xx_yyyzz, ts_xx_yyzz, ts_xx_yyzzz, ts_xx_yzzz, ts_xx_yzzzz, ts_xx_zzzz, ts_xx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xx_xxxxx[i] = 2.0 * ts_xx_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xx_xxxxy[i] = 2.0 * ts_xx_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xx_xxxxz[i] = 2.0 * ts_xx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xx_xxxyy[i] = 2.0 * ts_xx_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xx_xxxyz[i] = 2.0 * ts_xx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xx_xxxzz[i] = 4.0 * ts_xx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xx_xxyyy[i] = 2.0 * ts_xx_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xx_xxyyz[i] = 2.0 * ts_xx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xx_xxyzz[i] = 4.0 * ts_xx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xx_xxzzz[i] = 6.0 * ts_xx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xx_xyyyy[i] = 2.0 * ts_xx_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xx_xyyyz[i] = 2.0 * ts_xx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xx_xyyzz[i] = 4.0 * ts_xx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xx_xyzzz[i] = 6.0 * ts_xx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xx_xzzzz[i] = 8.0 * ts_xx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xx_yyyyy[i] = 2.0 * ts_xx_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xx_yyyyz[i] = 2.0 * ts_xx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xx_yyyzz[i] = 4.0 * ts_xx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xx_yyzzz[i] = 6.0 * ts_xx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xx_yzzzz[i] = 8.0 * ts_xx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xx_zzzzz[i] = 10.0 * ts_xx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 273-294 components of targeted buffer : DH

    auto gs_z_xy_xxxxx = pbuffer.data(idx_g_dh + 273);

    auto gs_z_xy_xxxxy = pbuffer.data(idx_g_dh + 274);

    auto gs_z_xy_xxxxz = pbuffer.data(idx_g_dh + 275);

    auto gs_z_xy_xxxyy = pbuffer.data(idx_g_dh + 276);

    auto gs_z_xy_xxxyz = pbuffer.data(idx_g_dh + 277);

    auto gs_z_xy_xxxzz = pbuffer.data(idx_g_dh + 278);

    auto gs_z_xy_xxyyy = pbuffer.data(idx_g_dh + 279);

    auto gs_z_xy_xxyyz = pbuffer.data(idx_g_dh + 280);

    auto gs_z_xy_xxyzz = pbuffer.data(idx_g_dh + 281);

    auto gs_z_xy_xxzzz = pbuffer.data(idx_g_dh + 282);

    auto gs_z_xy_xyyyy = pbuffer.data(idx_g_dh + 283);

    auto gs_z_xy_xyyyz = pbuffer.data(idx_g_dh + 284);

    auto gs_z_xy_xyyzz = pbuffer.data(idx_g_dh + 285);

    auto gs_z_xy_xyzzz = pbuffer.data(idx_g_dh + 286);

    auto gs_z_xy_xzzzz = pbuffer.data(idx_g_dh + 287);

    auto gs_z_xy_yyyyy = pbuffer.data(idx_g_dh + 288);

    auto gs_z_xy_yyyyz = pbuffer.data(idx_g_dh + 289);

    auto gs_z_xy_yyyzz = pbuffer.data(idx_g_dh + 290);

    auto gs_z_xy_yyzzz = pbuffer.data(idx_g_dh + 291);

    auto gs_z_xy_yzzzz = pbuffer.data(idx_g_dh + 292);

    auto gs_z_xy_zzzzz = pbuffer.data(idx_g_dh + 293);

    #pragma omp simd aligned(gc_z, gs_z_xy_xxxxx, gs_z_xy_xxxxy, gs_z_xy_xxxxz, gs_z_xy_xxxyy, gs_z_xy_xxxyz, gs_z_xy_xxxzz, gs_z_xy_xxyyy, gs_z_xy_xxyyz, gs_z_xy_xxyzz, gs_z_xy_xxzzz, gs_z_xy_xyyyy, gs_z_xy_xyyyz, gs_z_xy_xyyzz, gs_z_xy_xyzzz, gs_z_xy_xzzzz, gs_z_xy_yyyyy, gs_z_xy_yyyyz, gs_z_xy_yyyzz, gs_z_xy_yyzzz, gs_z_xy_yzzzz, gs_z_xy_zzzzz, ts_xy_xxxx, ts_xy_xxxxx, ts_xy_xxxxy, ts_xy_xxxxz, ts_xy_xxxy, ts_xy_xxxyy, ts_xy_xxxyz, ts_xy_xxxz, ts_xy_xxxzz, ts_xy_xxyy, ts_xy_xxyyy, ts_xy_xxyyz, ts_xy_xxyz, ts_xy_xxyzz, ts_xy_xxzz, ts_xy_xxzzz, ts_xy_xyyy, ts_xy_xyyyy, ts_xy_xyyyz, ts_xy_xyyz, ts_xy_xyyzz, ts_xy_xyzz, ts_xy_xyzzz, ts_xy_xzzz, ts_xy_xzzzz, ts_xy_yyyy, ts_xy_yyyyy, ts_xy_yyyyz, ts_xy_yyyz, ts_xy_yyyzz, ts_xy_yyzz, ts_xy_yyzzz, ts_xy_yzzz, ts_xy_yzzzz, ts_xy_zzzz, ts_xy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xy_xxxxx[i] = 2.0 * ts_xy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xy_xxxxy[i] = 2.0 * ts_xy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xy_xxxxz[i] = 2.0 * ts_xy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xy_xxxyy[i] = 2.0 * ts_xy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xy_xxxyz[i] = 2.0 * ts_xy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xy_xxxzz[i] = 4.0 * ts_xy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xy_xxyyy[i] = 2.0 * ts_xy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xy_xxyyz[i] = 2.0 * ts_xy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xy_xxyzz[i] = 4.0 * ts_xy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xy_xxzzz[i] = 6.0 * ts_xy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xy_xyyyy[i] = 2.0 * ts_xy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xy_xyyyz[i] = 2.0 * ts_xy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xy_xyyzz[i] = 4.0 * ts_xy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xy_xyzzz[i] = 6.0 * ts_xy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xy_xzzzz[i] = 8.0 * ts_xy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xy_yyyyy[i] = 2.0 * ts_xy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xy_yyyyz[i] = 2.0 * ts_xy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xy_yyyzz[i] = 4.0 * ts_xy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xy_yyzzz[i] = 6.0 * ts_xy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xy_yzzzz[i] = 8.0 * ts_xy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xy_zzzzz[i] = 10.0 * ts_xy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 294-315 components of targeted buffer : DH

    auto gs_z_xz_xxxxx = pbuffer.data(idx_g_dh + 294);

    auto gs_z_xz_xxxxy = pbuffer.data(idx_g_dh + 295);

    auto gs_z_xz_xxxxz = pbuffer.data(idx_g_dh + 296);

    auto gs_z_xz_xxxyy = pbuffer.data(idx_g_dh + 297);

    auto gs_z_xz_xxxyz = pbuffer.data(idx_g_dh + 298);

    auto gs_z_xz_xxxzz = pbuffer.data(idx_g_dh + 299);

    auto gs_z_xz_xxyyy = pbuffer.data(idx_g_dh + 300);

    auto gs_z_xz_xxyyz = pbuffer.data(idx_g_dh + 301);

    auto gs_z_xz_xxyzz = pbuffer.data(idx_g_dh + 302);

    auto gs_z_xz_xxzzz = pbuffer.data(idx_g_dh + 303);

    auto gs_z_xz_xyyyy = pbuffer.data(idx_g_dh + 304);

    auto gs_z_xz_xyyyz = pbuffer.data(idx_g_dh + 305);

    auto gs_z_xz_xyyzz = pbuffer.data(idx_g_dh + 306);

    auto gs_z_xz_xyzzz = pbuffer.data(idx_g_dh + 307);

    auto gs_z_xz_xzzzz = pbuffer.data(idx_g_dh + 308);

    auto gs_z_xz_yyyyy = pbuffer.data(idx_g_dh + 309);

    auto gs_z_xz_yyyyz = pbuffer.data(idx_g_dh + 310);

    auto gs_z_xz_yyyzz = pbuffer.data(idx_g_dh + 311);

    auto gs_z_xz_yyzzz = pbuffer.data(idx_g_dh + 312);

    auto gs_z_xz_yzzzz = pbuffer.data(idx_g_dh + 313);

    auto gs_z_xz_zzzzz = pbuffer.data(idx_g_dh + 314);

    #pragma omp simd aligned(gc_z, gs_z_xz_xxxxx, gs_z_xz_xxxxy, gs_z_xz_xxxxz, gs_z_xz_xxxyy, gs_z_xz_xxxyz, gs_z_xz_xxxzz, gs_z_xz_xxyyy, gs_z_xz_xxyyz, gs_z_xz_xxyzz, gs_z_xz_xxzzz, gs_z_xz_xyyyy, gs_z_xz_xyyyz, gs_z_xz_xyyzz, gs_z_xz_xyzzz, gs_z_xz_xzzzz, gs_z_xz_yyyyy, gs_z_xz_yyyyz, gs_z_xz_yyyzz, gs_z_xz_yyzzz, gs_z_xz_yzzzz, gs_z_xz_zzzzz, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxxz, ts_x_xxxyy, ts_x_xxxyz, ts_x_xxxzz, ts_x_xxyyy, ts_x_xxyyz, ts_x_xxyzz, ts_x_xxzzz, ts_x_xyyyy, ts_x_xyyyz, ts_x_xyyzz, ts_x_xyzzz, ts_x_xzzzz, ts_x_yyyyy, ts_x_yyyyz, ts_x_yyyzz, ts_x_yyzzz, ts_x_yzzzz, ts_x_zzzzz, ts_xz_xxxx, ts_xz_xxxxx, ts_xz_xxxxy, ts_xz_xxxxz, ts_xz_xxxy, ts_xz_xxxyy, ts_xz_xxxyz, ts_xz_xxxz, ts_xz_xxxzz, ts_xz_xxyy, ts_xz_xxyyy, ts_xz_xxyyz, ts_xz_xxyz, ts_xz_xxyzz, ts_xz_xxzz, ts_xz_xxzzz, ts_xz_xyyy, ts_xz_xyyyy, ts_xz_xyyyz, ts_xz_xyyz, ts_xz_xyyzz, ts_xz_xyzz, ts_xz_xyzzz, ts_xz_xzzz, ts_xz_xzzzz, ts_xz_yyyy, ts_xz_yyyyy, ts_xz_yyyyz, ts_xz_yyyz, ts_xz_yyyzz, ts_xz_yyzz, ts_xz_yyzzz, ts_xz_yzzz, ts_xz_yzzzz, ts_xz_zzzz, ts_xz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xz_xxxxx[i] = 2.0 * ts_x_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xz_xxxxy[i] = 2.0 * ts_x_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xz_xxxxz[i] = 2.0 * ts_x_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xz_xxxyy[i] = 2.0 * ts_x_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xz_xxxyz[i] = 2.0 * ts_x_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xz_xxxzz[i] = 2.0 * ts_x_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xz_xxyyy[i] = 2.0 * ts_x_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xz_xxyyz[i] = 2.0 * ts_x_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xz_xxyzz[i] = 2.0 * ts_x_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xz_xxzzz[i] = 2.0 * ts_x_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xz_xyyyy[i] = 2.0 * ts_x_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xz_xyyyz[i] = 2.0 * ts_x_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xz_xyyzz[i] = 2.0 * ts_x_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xz_xyzzz[i] = 2.0 * ts_x_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xz_xzzzz[i] = 2.0 * ts_x_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xz_yyyyy[i] = 2.0 * ts_x_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xz_yyyyz[i] = 2.0 * ts_x_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xz_yyyzz[i] = 2.0 * ts_x_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xz_yyzzz[i] = 2.0 * ts_x_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xz_yzzzz[i] = 2.0 * ts_x_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xz_zzzzz[i] = 2.0 * ts_x_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 315-336 components of targeted buffer : DH

    auto gs_z_yy_xxxxx = pbuffer.data(idx_g_dh + 315);

    auto gs_z_yy_xxxxy = pbuffer.data(idx_g_dh + 316);

    auto gs_z_yy_xxxxz = pbuffer.data(idx_g_dh + 317);

    auto gs_z_yy_xxxyy = pbuffer.data(idx_g_dh + 318);

    auto gs_z_yy_xxxyz = pbuffer.data(idx_g_dh + 319);

    auto gs_z_yy_xxxzz = pbuffer.data(idx_g_dh + 320);

    auto gs_z_yy_xxyyy = pbuffer.data(idx_g_dh + 321);

    auto gs_z_yy_xxyyz = pbuffer.data(idx_g_dh + 322);

    auto gs_z_yy_xxyzz = pbuffer.data(idx_g_dh + 323);

    auto gs_z_yy_xxzzz = pbuffer.data(idx_g_dh + 324);

    auto gs_z_yy_xyyyy = pbuffer.data(idx_g_dh + 325);

    auto gs_z_yy_xyyyz = pbuffer.data(idx_g_dh + 326);

    auto gs_z_yy_xyyzz = pbuffer.data(idx_g_dh + 327);

    auto gs_z_yy_xyzzz = pbuffer.data(idx_g_dh + 328);

    auto gs_z_yy_xzzzz = pbuffer.data(idx_g_dh + 329);

    auto gs_z_yy_yyyyy = pbuffer.data(idx_g_dh + 330);

    auto gs_z_yy_yyyyz = pbuffer.data(idx_g_dh + 331);

    auto gs_z_yy_yyyzz = pbuffer.data(idx_g_dh + 332);

    auto gs_z_yy_yyzzz = pbuffer.data(idx_g_dh + 333);

    auto gs_z_yy_yzzzz = pbuffer.data(idx_g_dh + 334);

    auto gs_z_yy_zzzzz = pbuffer.data(idx_g_dh + 335);

    #pragma omp simd aligned(gc_z, gs_z_yy_xxxxx, gs_z_yy_xxxxy, gs_z_yy_xxxxz, gs_z_yy_xxxyy, gs_z_yy_xxxyz, gs_z_yy_xxxzz, gs_z_yy_xxyyy, gs_z_yy_xxyyz, gs_z_yy_xxyzz, gs_z_yy_xxzzz, gs_z_yy_xyyyy, gs_z_yy_xyyyz, gs_z_yy_xyyzz, gs_z_yy_xyzzz, gs_z_yy_xzzzz, gs_z_yy_yyyyy, gs_z_yy_yyyyz, gs_z_yy_yyyzz, gs_z_yy_yyzzz, gs_z_yy_yzzzz, gs_z_yy_zzzzz, ts_yy_xxxx, ts_yy_xxxxx, ts_yy_xxxxy, ts_yy_xxxxz, ts_yy_xxxy, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxxz, ts_yy_xxxzz, ts_yy_xxyy, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyz, ts_yy_xxyzz, ts_yy_xxzz, ts_yy_xxzzz, ts_yy_xyyy, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyz, ts_yy_xyyzz, ts_yy_xyzz, ts_yy_xyzzz, ts_yy_xzzz, ts_yy_xzzzz, ts_yy_yyyy, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyz, ts_yy_yyyzz, ts_yy_yyzz, ts_yy_yyzzz, ts_yy_yzzz, ts_yy_yzzzz, ts_yy_zzzz, ts_yy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yy_xxxxx[i] = 2.0 * ts_yy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yy_xxxxy[i] = 2.0 * ts_yy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yy_xxxxz[i] = 2.0 * ts_yy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yy_xxxyy[i] = 2.0 * ts_yy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yy_xxxyz[i] = 2.0 * ts_yy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yy_xxxzz[i] = 4.0 * ts_yy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yy_xxyyy[i] = 2.0 * ts_yy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yy_xxyyz[i] = 2.0 * ts_yy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yy_xxyzz[i] = 4.0 * ts_yy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yy_xxzzz[i] = 6.0 * ts_yy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yy_xyyyy[i] = 2.0 * ts_yy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yy_xyyyz[i] = 2.0 * ts_yy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yy_xyyzz[i] = 4.0 * ts_yy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yy_xyzzz[i] = 6.0 * ts_yy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yy_xzzzz[i] = 8.0 * ts_yy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yy_yyyyy[i] = 2.0 * ts_yy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yy_yyyyz[i] = 2.0 * ts_yy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yy_yyyzz[i] = 4.0 * ts_yy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yy_yyzzz[i] = 6.0 * ts_yy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yy_yzzzz[i] = 8.0 * ts_yy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yy_zzzzz[i] = 10.0 * ts_yy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 336-357 components of targeted buffer : DH

    auto gs_z_yz_xxxxx = pbuffer.data(idx_g_dh + 336);

    auto gs_z_yz_xxxxy = pbuffer.data(idx_g_dh + 337);

    auto gs_z_yz_xxxxz = pbuffer.data(idx_g_dh + 338);

    auto gs_z_yz_xxxyy = pbuffer.data(idx_g_dh + 339);

    auto gs_z_yz_xxxyz = pbuffer.data(idx_g_dh + 340);

    auto gs_z_yz_xxxzz = pbuffer.data(idx_g_dh + 341);

    auto gs_z_yz_xxyyy = pbuffer.data(idx_g_dh + 342);

    auto gs_z_yz_xxyyz = pbuffer.data(idx_g_dh + 343);

    auto gs_z_yz_xxyzz = pbuffer.data(idx_g_dh + 344);

    auto gs_z_yz_xxzzz = pbuffer.data(idx_g_dh + 345);

    auto gs_z_yz_xyyyy = pbuffer.data(idx_g_dh + 346);

    auto gs_z_yz_xyyyz = pbuffer.data(idx_g_dh + 347);

    auto gs_z_yz_xyyzz = pbuffer.data(idx_g_dh + 348);

    auto gs_z_yz_xyzzz = pbuffer.data(idx_g_dh + 349);

    auto gs_z_yz_xzzzz = pbuffer.data(idx_g_dh + 350);

    auto gs_z_yz_yyyyy = pbuffer.data(idx_g_dh + 351);

    auto gs_z_yz_yyyyz = pbuffer.data(idx_g_dh + 352);

    auto gs_z_yz_yyyzz = pbuffer.data(idx_g_dh + 353);

    auto gs_z_yz_yyzzz = pbuffer.data(idx_g_dh + 354);

    auto gs_z_yz_yzzzz = pbuffer.data(idx_g_dh + 355);

    auto gs_z_yz_zzzzz = pbuffer.data(idx_g_dh + 356);

    #pragma omp simd aligned(gc_z, gs_z_yz_xxxxx, gs_z_yz_xxxxy, gs_z_yz_xxxxz, gs_z_yz_xxxyy, gs_z_yz_xxxyz, gs_z_yz_xxxzz, gs_z_yz_xxyyy, gs_z_yz_xxyyz, gs_z_yz_xxyzz, gs_z_yz_xxzzz, gs_z_yz_xyyyy, gs_z_yz_xyyyz, gs_z_yz_xyyzz, gs_z_yz_xyzzz, gs_z_yz_xzzzz, gs_z_yz_yyyyy, gs_z_yz_yyyyz, gs_z_yz_yyyzz, gs_z_yz_yyzzz, gs_z_yz_yzzzz, gs_z_yz_zzzzz, ts_y_xxxxx, ts_y_xxxxy, ts_y_xxxxz, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxxzz, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyzz, ts_y_xxzzz, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyzz, ts_y_xyzzz, ts_y_xzzzz, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyzz, ts_y_yyzzz, ts_y_yzzzz, ts_y_zzzzz, ts_yz_xxxx, ts_yz_xxxxx, ts_yz_xxxxy, ts_yz_xxxxz, ts_yz_xxxy, ts_yz_xxxyy, ts_yz_xxxyz, ts_yz_xxxz, ts_yz_xxxzz, ts_yz_xxyy, ts_yz_xxyyy, ts_yz_xxyyz, ts_yz_xxyz, ts_yz_xxyzz, ts_yz_xxzz, ts_yz_xxzzz, ts_yz_xyyy, ts_yz_xyyyy, ts_yz_xyyyz, ts_yz_xyyz, ts_yz_xyyzz, ts_yz_xyzz, ts_yz_xyzzz, ts_yz_xzzz, ts_yz_xzzzz, ts_yz_yyyy, ts_yz_yyyyy, ts_yz_yyyyz, ts_yz_yyyz, ts_yz_yyyzz, ts_yz_yyzz, ts_yz_yyzzz, ts_yz_yzzz, ts_yz_yzzzz, ts_yz_zzzz, ts_yz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yz_xxxxx[i] = 2.0 * ts_y_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yz_xxxxy[i] = 2.0 * ts_y_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yz_xxxxz[i] = 2.0 * ts_y_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yz_xxxyy[i] = 2.0 * ts_y_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yz_xxxyz[i] = 2.0 * ts_y_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yz_xxxzz[i] = 2.0 * ts_y_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yz_xxyyy[i] = 2.0 * ts_y_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yz_xxyyz[i] = 2.0 * ts_y_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yz_xxyzz[i] = 2.0 * ts_y_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yz_xxzzz[i] = 2.0 * ts_y_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yz_xyyyy[i] = 2.0 * ts_y_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yz_xyyyz[i] = 2.0 * ts_y_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yz_xyyzz[i] = 2.0 * ts_y_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yz_xyzzz[i] = 2.0 * ts_y_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yz_xzzzz[i] = 2.0 * ts_y_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yz_yyyyy[i] = 2.0 * ts_y_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yz_yyyyz[i] = 2.0 * ts_y_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yz_yyyzz[i] = 2.0 * ts_y_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yz_yyzzz[i] = 2.0 * ts_y_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yz_yzzzz[i] = 2.0 * ts_y_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yz_zzzzz[i] = 2.0 * ts_y_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 357-378 components of targeted buffer : DH

    auto gs_z_zz_xxxxx = pbuffer.data(idx_g_dh + 357);

    auto gs_z_zz_xxxxy = pbuffer.data(idx_g_dh + 358);

    auto gs_z_zz_xxxxz = pbuffer.data(idx_g_dh + 359);

    auto gs_z_zz_xxxyy = pbuffer.data(idx_g_dh + 360);

    auto gs_z_zz_xxxyz = pbuffer.data(idx_g_dh + 361);

    auto gs_z_zz_xxxzz = pbuffer.data(idx_g_dh + 362);

    auto gs_z_zz_xxyyy = pbuffer.data(idx_g_dh + 363);

    auto gs_z_zz_xxyyz = pbuffer.data(idx_g_dh + 364);

    auto gs_z_zz_xxyzz = pbuffer.data(idx_g_dh + 365);

    auto gs_z_zz_xxzzz = pbuffer.data(idx_g_dh + 366);

    auto gs_z_zz_xyyyy = pbuffer.data(idx_g_dh + 367);

    auto gs_z_zz_xyyyz = pbuffer.data(idx_g_dh + 368);

    auto gs_z_zz_xyyzz = pbuffer.data(idx_g_dh + 369);

    auto gs_z_zz_xyzzz = pbuffer.data(idx_g_dh + 370);

    auto gs_z_zz_xzzzz = pbuffer.data(idx_g_dh + 371);

    auto gs_z_zz_yyyyy = pbuffer.data(idx_g_dh + 372);

    auto gs_z_zz_yyyyz = pbuffer.data(idx_g_dh + 373);

    auto gs_z_zz_yyyzz = pbuffer.data(idx_g_dh + 374);

    auto gs_z_zz_yyzzz = pbuffer.data(idx_g_dh + 375);

    auto gs_z_zz_yzzzz = pbuffer.data(idx_g_dh + 376);

    auto gs_z_zz_zzzzz = pbuffer.data(idx_g_dh + 377);

    #pragma omp simd aligned(gc_z, gs_z_zz_xxxxx, gs_z_zz_xxxxy, gs_z_zz_xxxxz, gs_z_zz_xxxyy, gs_z_zz_xxxyz, gs_z_zz_xxxzz, gs_z_zz_xxyyy, gs_z_zz_xxyyz, gs_z_zz_xxyzz, gs_z_zz_xxzzz, gs_z_zz_xyyyy, gs_z_zz_xyyyz, gs_z_zz_xyyzz, gs_z_zz_xyzzz, gs_z_zz_xzzzz, gs_z_zz_yyyyy, gs_z_zz_yyyyz, gs_z_zz_yyyzz, gs_z_zz_yyzzz, gs_z_zz_yzzzz, gs_z_zz_zzzzz, ts_z_xxxxx, ts_z_xxxxy, ts_z_xxxxz, ts_z_xxxyy, ts_z_xxxyz, ts_z_xxxzz, ts_z_xxyyy, ts_z_xxyyz, ts_z_xxyzz, ts_z_xxzzz, ts_z_xyyyy, ts_z_xyyyz, ts_z_xyyzz, ts_z_xyzzz, ts_z_xzzzz, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyzz, ts_z_yyzzz, ts_z_yzzzz, ts_z_zzzzz, ts_zz_xxxx, ts_zz_xxxxx, ts_zz_xxxxy, ts_zz_xxxxz, ts_zz_xxxy, ts_zz_xxxyy, ts_zz_xxxyz, ts_zz_xxxz, ts_zz_xxxzz, ts_zz_xxyy, ts_zz_xxyyy, ts_zz_xxyyz, ts_zz_xxyz, ts_zz_xxyzz, ts_zz_xxzz, ts_zz_xxzzz, ts_zz_xyyy, ts_zz_xyyyy, ts_zz_xyyyz, ts_zz_xyyz, ts_zz_xyyzz, ts_zz_xyzz, ts_zz_xyzzz, ts_zz_xzzz, ts_zz_xzzzz, ts_zz_yyyy, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyz, ts_zz_yyyzz, ts_zz_yyzz, ts_zz_yyzzz, ts_zz_yzzz, ts_zz_yzzzz, ts_zz_zzzz, ts_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zz_xxxxx[i] = 4.0 * ts_z_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_zz_xxxxy[i] = 4.0 * ts_z_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_zz_xxxxz[i] = 4.0 * ts_z_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_zz_xxxyy[i] = 4.0 * ts_z_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_zz_xxxyz[i] = 4.0 * ts_z_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_zz_xxxzz[i] = 4.0 * ts_z_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_zz_xxyyy[i] = 4.0 * ts_z_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_zz_xxyyz[i] = 4.0 * ts_z_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_zz_xxyzz[i] = 4.0 * ts_z_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_zz_xxzzz[i] = 4.0 * ts_z_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_zz_xyyyy[i] = 4.0 * ts_z_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_zz_xyyyz[i] = 4.0 * ts_z_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_zz_xyyzz[i] = 4.0 * ts_z_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_zz_xyzzz[i] = 4.0 * ts_z_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_zz_xzzzz[i] = 4.0 * ts_z_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_zz_yyyyy[i] = 4.0 * ts_z_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_zz_yyyyz[i] = 4.0 * ts_z_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_zz_yyyzz[i] = 4.0 * ts_z_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_zz_yyzzz[i] = 4.0 * ts_z_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_zz_yzzzz[i] = 4.0 * ts_z_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_zz_zzzzz[i] = 4.0 * ts_z_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_zz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_zzzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

