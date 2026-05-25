#include "LocalCorePotentialPrimRecDH.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_dh(CSimdArray<double>& pbuffer, 
                                  const size_t idx_dh,
                                  const size_t idx_sh,
                                  const size_t idx_pg,
                                  const size_t idx_ph,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx = pbuffer.data(idx_sh);

    auto tg_0_xxxxy = pbuffer.data(idx_sh + 1);

    auto tg_0_xxxxz = pbuffer.data(idx_sh + 2);

    auto tg_0_xxxyy = pbuffer.data(idx_sh + 3);

    auto tg_0_xxxyz = pbuffer.data(idx_sh + 4);

    auto tg_0_xxxzz = pbuffer.data(idx_sh + 5);

    auto tg_0_xxyyy = pbuffer.data(idx_sh + 6);

    auto tg_0_xxyyz = pbuffer.data(idx_sh + 7);

    auto tg_0_xxyzz = pbuffer.data(idx_sh + 8);

    auto tg_0_xxzzz = pbuffer.data(idx_sh + 9);

    auto tg_0_xyyyy = pbuffer.data(idx_sh + 10);

    auto tg_0_xyyyz = pbuffer.data(idx_sh + 11);

    auto tg_0_xyyzz = pbuffer.data(idx_sh + 12);

    auto tg_0_xyzzz = pbuffer.data(idx_sh + 13);

    auto tg_0_xzzzz = pbuffer.data(idx_sh + 14);

    auto tg_0_yyyyy = pbuffer.data(idx_sh + 15);

    auto tg_0_yyyyz = pbuffer.data(idx_sh + 16);

    auto tg_0_yyyzz = pbuffer.data(idx_sh + 17);

    auto tg_0_yyzzz = pbuffer.data(idx_sh + 18);

    auto tg_0_yzzzz = pbuffer.data(idx_sh + 19);

    auto tg_0_zzzzz = pbuffer.data(idx_sh + 20);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx = pbuffer.data(idx_pg);

    auto tg_x_xxxy = pbuffer.data(idx_pg + 1);

    auto tg_x_xxxz = pbuffer.data(idx_pg + 2);

    auto tg_x_xxyy = pbuffer.data(idx_pg + 3);

    auto tg_x_xxyz = pbuffer.data(idx_pg + 4);

    auto tg_x_xxzz = pbuffer.data(idx_pg + 5);

    auto tg_x_xyyy = pbuffer.data(idx_pg + 6);

    auto tg_x_xyyz = pbuffer.data(idx_pg + 7);

    auto tg_x_xyzz = pbuffer.data(idx_pg + 8);

    auto tg_x_xzzz = pbuffer.data(idx_pg + 9);

    auto tg_x_yyyy = pbuffer.data(idx_pg + 10);

    auto tg_x_yyyz = pbuffer.data(idx_pg + 11);

    auto tg_x_yyzz = pbuffer.data(idx_pg + 12);

    auto tg_x_yzzz = pbuffer.data(idx_pg + 13);

    auto tg_x_zzzz = pbuffer.data(idx_pg + 14);

    auto tg_y_xxxx = pbuffer.data(idx_pg + 15);

    auto tg_y_xxxy = pbuffer.data(idx_pg + 16);

    auto tg_y_xxxz = pbuffer.data(idx_pg + 17);

    auto tg_y_xxyy = pbuffer.data(idx_pg + 18);

    auto tg_y_xxyz = pbuffer.data(idx_pg + 19);

    auto tg_y_xxzz = pbuffer.data(idx_pg + 20);

    auto tg_y_xyyy = pbuffer.data(idx_pg + 21);

    auto tg_y_xyyz = pbuffer.data(idx_pg + 22);

    auto tg_y_xyzz = pbuffer.data(idx_pg + 23);

    auto tg_y_xzzz = pbuffer.data(idx_pg + 24);

    auto tg_y_yyyy = pbuffer.data(idx_pg + 25);

    auto tg_y_yyyz = pbuffer.data(idx_pg + 26);

    auto tg_y_yyzz = pbuffer.data(idx_pg + 27);

    auto tg_y_yzzz = pbuffer.data(idx_pg + 28);

    auto tg_y_zzzz = pbuffer.data(idx_pg + 29);

    auto tg_z_xxxx = pbuffer.data(idx_pg + 30);

    auto tg_z_xxxy = pbuffer.data(idx_pg + 31);

    auto tg_z_xxxz = pbuffer.data(idx_pg + 32);

    auto tg_z_xxyy = pbuffer.data(idx_pg + 33);

    auto tg_z_xxyz = pbuffer.data(idx_pg + 34);

    auto tg_z_xxzz = pbuffer.data(idx_pg + 35);

    auto tg_z_xyyy = pbuffer.data(idx_pg + 36);

    auto tg_z_xyyz = pbuffer.data(idx_pg + 37);

    auto tg_z_xyzz = pbuffer.data(idx_pg + 38);

    auto tg_z_xzzz = pbuffer.data(idx_pg + 39);

    auto tg_z_yyyy = pbuffer.data(idx_pg + 40);

    auto tg_z_yyyz = pbuffer.data(idx_pg + 41);

    auto tg_z_yyzz = pbuffer.data(idx_pg + 42);

    auto tg_z_yzzz = pbuffer.data(idx_pg + 43);

    auto tg_z_zzzz = pbuffer.data(idx_pg + 44);

    // Set up components of auxiliary buffer : PH

    auto tg_x_xxxxx = pbuffer.data(idx_ph);

    auto tg_x_xxxxy = pbuffer.data(idx_ph + 1);

    auto tg_x_xxxxz = pbuffer.data(idx_ph + 2);

    auto tg_x_xxxyy = pbuffer.data(idx_ph + 3);

    auto tg_x_xxxyz = pbuffer.data(idx_ph + 4);

    auto tg_x_xxxzz = pbuffer.data(idx_ph + 5);

    auto tg_x_xxyyy = pbuffer.data(idx_ph + 6);

    auto tg_x_xxyyz = pbuffer.data(idx_ph + 7);

    auto tg_x_xxyzz = pbuffer.data(idx_ph + 8);

    auto tg_x_xxzzz = pbuffer.data(idx_ph + 9);

    auto tg_x_xyyyy = pbuffer.data(idx_ph + 10);

    auto tg_x_xyyyz = pbuffer.data(idx_ph + 11);

    auto tg_x_xyyzz = pbuffer.data(idx_ph + 12);

    auto tg_x_xyzzz = pbuffer.data(idx_ph + 13);

    auto tg_x_xzzzz = pbuffer.data(idx_ph + 14);

    auto tg_x_yyyyy = pbuffer.data(idx_ph + 15);

    auto tg_x_yyyyz = pbuffer.data(idx_ph + 16);

    auto tg_x_yyyzz = pbuffer.data(idx_ph + 17);

    auto tg_x_yyzzz = pbuffer.data(idx_ph + 18);

    auto tg_x_yzzzz = pbuffer.data(idx_ph + 19);

    auto tg_x_zzzzz = pbuffer.data(idx_ph + 20);

    auto tg_y_xxxxx = pbuffer.data(idx_ph + 21);

    auto tg_y_xxxxy = pbuffer.data(idx_ph + 22);

    auto tg_y_xxxxz = pbuffer.data(idx_ph + 23);

    auto tg_y_xxxyy = pbuffer.data(idx_ph + 24);

    auto tg_y_xxxyz = pbuffer.data(idx_ph + 25);

    auto tg_y_xxxzz = pbuffer.data(idx_ph + 26);

    auto tg_y_xxyyy = pbuffer.data(idx_ph + 27);

    auto tg_y_xxyyz = pbuffer.data(idx_ph + 28);

    auto tg_y_xxyzz = pbuffer.data(idx_ph + 29);

    auto tg_y_xxzzz = pbuffer.data(idx_ph + 30);

    auto tg_y_xyyyy = pbuffer.data(idx_ph + 31);

    auto tg_y_xyyyz = pbuffer.data(idx_ph + 32);

    auto tg_y_xyyzz = pbuffer.data(idx_ph + 33);

    auto tg_y_xyzzz = pbuffer.data(idx_ph + 34);

    auto tg_y_xzzzz = pbuffer.data(idx_ph + 35);

    auto tg_y_yyyyy = pbuffer.data(idx_ph + 36);

    auto tg_y_yyyyz = pbuffer.data(idx_ph + 37);

    auto tg_y_yyyzz = pbuffer.data(idx_ph + 38);

    auto tg_y_yyzzz = pbuffer.data(idx_ph + 39);

    auto tg_y_yzzzz = pbuffer.data(idx_ph + 40);

    auto tg_y_zzzzz = pbuffer.data(idx_ph + 41);

    auto tg_z_xxxxx = pbuffer.data(idx_ph + 42);

    auto tg_z_xxxxy = pbuffer.data(idx_ph + 43);

    auto tg_z_xxxxz = pbuffer.data(idx_ph + 44);

    auto tg_z_xxxyy = pbuffer.data(idx_ph + 45);

    auto tg_z_xxxyz = pbuffer.data(idx_ph + 46);

    auto tg_z_xxxzz = pbuffer.data(idx_ph + 47);

    auto tg_z_xxyyy = pbuffer.data(idx_ph + 48);

    auto tg_z_xxyyz = pbuffer.data(idx_ph + 49);

    auto tg_z_xxyzz = pbuffer.data(idx_ph + 50);

    auto tg_z_xxzzz = pbuffer.data(idx_ph + 51);

    auto tg_z_xyyyy = pbuffer.data(idx_ph + 52);

    auto tg_z_xyyyz = pbuffer.data(idx_ph + 53);

    auto tg_z_xyyzz = pbuffer.data(idx_ph + 54);

    auto tg_z_xyzzz = pbuffer.data(idx_ph + 55);

    auto tg_z_xzzzz = pbuffer.data(idx_ph + 56);

    auto tg_z_yyyyy = pbuffer.data(idx_ph + 57);

    auto tg_z_yyyyz = pbuffer.data(idx_ph + 58);

    auto tg_z_yyyzz = pbuffer.data(idx_ph + 59);

    auto tg_z_yyzzz = pbuffer.data(idx_ph + 60);

    auto tg_z_yzzzz = pbuffer.data(idx_ph + 61);

    auto tg_z_zzzzz = pbuffer.data(idx_ph + 62);

    // Set up components of targeted buffer : DH

    auto tg_xx_xxxxx = pbuffer.data(idx_dh);

    auto tg_xx_xxxxy = pbuffer.data(idx_dh + 1);

    auto tg_xx_xxxxz = pbuffer.data(idx_dh + 2);

    auto tg_xx_xxxyy = pbuffer.data(idx_dh + 3);

    auto tg_xx_xxxyz = pbuffer.data(idx_dh + 4);

    auto tg_xx_xxxzz = pbuffer.data(idx_dh + 5);

    auto tg_xx_xxyyy = pbuffer.data(idx_dh + 6);

    auto tg_xx_xxyyz = pbuffer.data(idx_dh + 7);

    auto tg_xx_xxyzz = pbuffer.data(idx_dh + 8);

    auto tg_xx_xxzzz = pbuffer.data(idx_dh + 9);

    auto tg_xx_xyyyy = pbuffer.data(idx_dh + 10);

    auto tg_xx_xyyyz = pbuffer.data(idx_dh + 11);

    auto tg_xx_xyyzz = pbuffer.data(idx_dh + 12);

    auto tg_xx_xyzzz = pbuffer.data(idx_dh + 13);

    auto tg_xx_xzzzz = pbuffer.data(idx_dh + 14);

    auto tg_xx_yyyyy = pbuffer.data(idx_dh + 15);

    auto tg_xx_yyyyz = pbuffer.data(idx_dh + 16);

    auto tg_xx_yyyzz = pbuffer.data(idx_dh + 17);

    auto tg_xx_yyzzz = pbuffer.data(idx_dh + 18);

    auto tg_xx_yzzzz = pbuffer.data(idx_dh + 19);

    auto tg_xx_zzzzz = pbuffer.data(idx_dh + 20);

    auto tg_xy_xxxxx = pbuffer.data(idx_dh + 21);

    auto tg_xy_xxxxy = pbuffer.data(idx_dh + 22);

    auto tg_xy_xxxxz = pbuffer.data(idx_dh + 23);

    auto tg_xy_xxxyy = pbuffer.data(idx_dh + 24);

    auto tg_xy_xxxyz = pbuffer.data(idx_dh + 25);

    auto tg_xy_xxxzz = pbuffer.data(idx_dh + 26);

    auto tg_xy_xxyyy = pbuffer.data(idx_dh + 27);

    auto tg_xy_xxyyz = pbuffer.data(idx_dh + 28);

    auto tg_xy_xxyzz = pbuffer.data(idx_dh + 29);

    auto tg_xy_xxzzz = pbuffer.data(idx_dh + 30);

    auto tg_xy_xyyyy = pbuffer.data(idx_dh + 31);

    auto tg_xy_xyyyz = pbuffer.data(idx_dh + 32);

    auto tg_xy_xyyzz = pbuffer.data(idx_dh + 33);

    auto tg_xy_xyzzz = pbuffer.data(idx_dh + 34);

    auto tg_xy_xzzzz = pbuffer.data(idx_dh + 35);

    auto tg_xy_yyyyy = pbuffer.data(idx_dh + 36);

    auto tg_xy_yyyyz = pbuffer.data(idx_dh + 37);

    auto tg_xy_yyyzz = pbuffer.data(idx_dh + 38);

    auto tg_xy_yyzzz = pbuffer.data(idx_dh + 39);

    auto tg_xy_yzzzz = pbuffer.data(idx_dh + 40);

    auto tg_xy_zzzzz = pbuffer.data(idx_dh + 41);

    auto tg_xz_xxxxx = pbuffer.data(idx_dh + 42);

    auto tg_xz_xxxxy = pbuffer.data(idx_dh + 43);

    auto tg_xz_xxxxz = pbuffer.data(idx_dh + 44);

    auto tg_xz_xxxyy = pbuffer.data(idx_dh + 45);

    auto tg_xz_xxxyz = pbuffer.data(idx_dh + 46);

    auto tg_xz_xxxzz = pbuffer.data(idx_dh + 47);

    auto tg_xz_xxyyy = pbuffer.data(idx_dh + 48);

    auto tg_xz_xxyyz = pbuffer.data(idx_dh + 49);

    auto tg_xz_xxyzz = pbuffer.data(idx_dh + 50);

    auto tg_xz_xxzzz = pbuffer.data(idx_dh + 51);

    auto tg_xz_xyyyy = pbuffer.data(idx_dh + 52);

    auto tg_xz_xyyyz = pbuffer.data(idx_dh + 53);

    auto tg_xz_xyyzz = pbuffer.data(idx_dh + 54);

    auto tg_xz_xyzzz = pbuffer.data(idx_dh + 55);

    auto tg_xz_xzzzz = pbuffer.data(idx_dh + 56);

    auto tg_xz_yyyyy = pbuffer.data(idx_dh + 57);

    auto tg_xz_yyyyz = pbuffer.data(idx_dh + 58);

    auto tg_xz_yyyzz = pbuffer.data(idx_dh + 59);

    auto tg_xz_yyzzz = pbuffer.data(idx_dh + 60);

    auto tg_xz_yzzzz = pbuffer.data(idx_dh + 61);

    auto tg_xz_zzzzz = pbuffer.data(idx_dh + 62);

    auto tg_yy_xxxxx = pbuffer.data(idx_dh + 63);

    auto tg_yy_xxxxy = pbuffer.data(idx_dh + 64);

    auto tg_yy_xxxxz = pbuffer.data(idx_dh + 65);

    auto tg_yy_xxxyy = pbuffer.data(idx_dh + 66);

    auto tg_yy_xxxyz = pbuffer.data(idx_dh + 67);

    auto tg_yy_xxxzz = pbuffer.data(idx_dh + 68);

    auto tg_yy_xxyyy = pbuffer.data(idx_dh + 69);

    auto tg_yy_xxyyz = pbuffer.data(idx_dh + 70);

    auto tg_yy_xxyzz = pbuffer.data(idx_dh + 71);

    auto tg_yy_xxzzz = pbuffer.data(idx_dh + 72);

    auto tg_yy_xyyyy = pbuffer.data(idx_dh + 73);

    auto tg_yy_xyyyz = pbuffer.data(idx_dh + 74);

    auto tg_yy_xyyzz = pbuffer.data(idx_dh + 75);

    auto tg_yy_xyzzz = pbuffer.data(idx_dh + 76);

    auto tg_yy_xzzzz = pbuffer.data(idx_dh + 77);

    auto tg_yy_yyyyy = pbuffer.data(idx_dh + 78);

    auto tg_yy_yyyyz = pbuffer.data(idx_dh + 79);

    auto tg_yy_yyyzz = pbuffer.data(idx_dh + 80);

    auto tg_yy_yyzzz = pbuffer.data(idx_dh + 81);

    auto tg_yy_yzzzz = pbuffer.data(idx_dh + 82);

    auto tg_yy_zzzzz = pbuffer.data(idx_dh + 83);

    auto tg_yz_xxxxx = pbuffer.data(idx_dh + 84);

    auto tg_yz_xxxxy = pbuffer.data(idx_dh + 85);

    auto tg_yz_xxxxz = pbuffer.data(idx_dh + 86);

    auto tg_yz_xxxyy = pbuffer.data(idx_dh + 87);

    auto tg_yz_xxxyz = pbuffer.data(idx_dh + 88);

    auto tg_yz_xxxzz = pbuffer.data(idx_dh + 89);

    auto tg_yz_xxyyy = pbuffer.data(idx_dh + 90);

    auto tg_yz_xxyyz = pbuffer.data(idx_dh + 91);

    auto tg_yz_xxyzz = pbuffer.data(idx_dh + 92);

    auto tg_yz_xxzzz = pbuffer.data(idx_dh + 93);

    auto tg_yz_xyyyy = pbuffer.data(idx_dh + 94);

    auto tg_yz_xyyyz = pbuffer.data(idx_dh + 95);

    auto tg_yz_xyyzz = pbuffer.data(idx_dh + 96);

    auto tg_yz_xyzzz = pbuffer.data(idx_dh + 97);

    auto tg_yz_xzzzz = pbuffer.data(idx_dh + 98);

    auto tg_yz_yyyyy = pbuffer.data(idx_dh + 99);

    auto tg_yz_yyyyz = pbuffer.data(idx_dh + 100);

    auto tg_yz_yyyzz = pbuffer.data(idx_dh + 101);

    auto tg_yz_yyzzz = pbuffer.data(idx_dh + 102);

    auto tg_yz_yzzzz = pbuffer.data(idx_dh + 103);

    auto tg_yz_zzzzz = pbuffer.data(idx_dh + 104);

    auto tg_zz_xxxxx = pbuffer.data(idx_dh + 105);

    auto tg_zz_xxxxy = pbuffer.data(idx_dh + 106);

    auto tg_zz_xxxxz = pbuffer.data(idx_dh + 107);

    auto tg_zz_xxxyy = pbuffer.data(idx_dh + 108);

    auto tg_zz_xxxyz = pbuffer.data(idx_dh + 109);

    auto tg_zz_xxxzz = pbuffer.data(idx_dh + 110);

    auto tg_zz_xxyyy = pbuffer.data(idx_dh + 111);

    auto tg_zz_xxyyz = pbuffer.data(idx_dh + 112);

    auto tg_zz_xxyzz = pbuffer.data(idx_dh + 113);

    auto tg_zz_xxzzz = pbuffer.data(idx_dh + 114);

    auto tg_zz_xyyyy = pbuffer.data(idx_dh + 115);

    auto tg_zz_xyyyz = pbuffer.data(idx_dh + 116);

    auto tg_zz_xyyzz = pbuffer.data(idx_dh + 117);

    auto tg_zz_xyzzz = pbuffer.data(idx_dh + 118);

    auto tg_zz_xzzzz = pbuffer.data(idx_dh + 119);

    auto tg_zz_yyyyy = pbuffer.data(idx_dh + 120);

    auto tg_zz_yyyyz = pbuffer.data(idx_dh + 121);

    auto tg_zz_yyyzz = pbuffer.data(idx_dh + 122);

    auto tg_zz_yyzzz = pbuffer.data(idx_dh + 123);

    auto tg_zz_yzzzz = pbuffer.data(idx_dh + 124);

    auto tg_zz_zzzzz = pbuffer.data(idx_dh + 125);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_xxxxx, tg_0_xxxxy, tg_0_xxxxz, tg_0_xxxyy, tg_0_xxxyz, tg_0_xxxzz, tg_0_xxyyy, tg_0_xxyyz, tg_0_xxyzz, tg_0_xxzzz, tg_0_xyyyy, tg_0_xyyyz, tg_0_xyyzz, tg_0_xyzzz, tg_0_xzzzz, tg_0_yyyyy, tg_0_yyyyz, tg_0_yyyzz, tg_0_yyzzz, tg_0_yzzzz, tg_0_zzzzz, tg_x_xxxx, tg_x_xxxxx, tg_x_xxxxy, tg_x_xxxxz, tg_x_xxxy, tg_x_xxxyy, tg_x_xxxyz, tg_x_xxxz, tg_x_xxxzz, tg_x_xxyy, tg_x_xxyyy, tg_x_xxyyz, tg_x_xxyz, tg_x_xxyzz, tg_x_xxzz, tg_x_xxzzz, tg_x_xyyy, tg_x_xyyyy, tg_x_xyyyz, tg_x_xyyz, tg_x_xyyzz, tg_x_xyzz, tg_x_xyzzz, tg_x_xzzz, tg_x_xzzzz, tg_x_yyyy, tg_x_yyyyy, tg_x_yyyyz, tg_x_yyyz, tg_x_yyyzz, tg_x_yyzz, tg_x_yyzzz, tg_x_yzzz, tg_x_yzzzz, tg_x_zzzz, tg_x_zzzzz, tg_xx_xxxxx, tg_xx_xxxxy, tg_xx_xxxxz, tg_xx_xxxyy, tg_xx_xxxyz, tg_xx_xxxzz, tg_xx_xxyyy, tg_xx_xxyyz, tg_xx_xxyzz, tg_xx_xxzzz, tg_xx_xyyyy, tg_xx_xyyyz, tg_xx_xyyzz, tg_xx_xyzzz, tg_xx_xzzzz, tg_xx_yyyyy, tg_xx_yyyyz, tg_xx_yyyzz, tg_xx_yyzzz, tg_xx_yzzzz, tg_xx_zzzzz, tg_xy_xxxxx, tg_xy_xxxxy, tg_xy_xxxxz, tg_xy_xxxyy, tg_xy_xxxyz, tg_xy_xxxzz, tg_xy_xxyyy, tg_xy_xxyyz, tg_xy_xxyzz, tg_xy_xxzzz, tg_xy_xyyyy, tg_xy_xyyyz, tg_xy_xyyzz, tg_xy_xyzzz, tg_xy_xzzzz, tg_xy_yyyyy, tg_xy_yyyyz, tg_xy_yyyzz, tg_xy_yyzzz, tg_xy_yzzzz, tg_xy_zzzzz, tg_xz_xxxxx, tg_xz_xxxxy, tg_xz_xxxxz, tg_xz_xxxyy, tg_xz_xxxyz, tg_xz_xxxzz, tg_xz_xxyyy, tg_xz_xxyyz, tg_xz_xxyzz, tg_xz_xxzzz, tg_xz_xyyyy, tg_xz_xyyyz, tg_xz_xyyzz, tg_xz_xyzzz, tg_xz_xzzzz, tg_xz_yyyyy, tg_xz_yyyyz, tg_xz_yyyzz, tg_xz_yyzzz, tg_xz_yzzzz, tg_xz_zzzzz, tg_y_xxxx, tg_y_xxxxx, tg_y_xxxxy, tg_y_xxxxz, tg_y_xxxy, tg_y_xxxyy, tg_y_xxxyz, tg_y_xxxz, tg_y_xxxzz, tg_y_xxyy, tg_y_xxyyy, tg_y_xxyyz, tg_y_xxyz, tg_y_xxyzz, tg_y_xxzz, tg_y_xxzzz, tg_y_xyyy, tg_y_xyyyy, tg_y_xyyyz, tg_y_xyyz, tg_y_xyyzz, tg_y_xyzz, tg_y_xyzzz, tg_y_xzzz, tg_y_xzzzz, tg_y_yyyy, tg_y_yyyyy, tg_y_yyyyz, tg_y_yyyz, tg_y_yyyzz, tg_y_yyzz, tg_y_yyzzz, tg_y_yzzz, tg_y_yzzzz, tg_y_zzzz, tg_y_zzzzz, tg_yy_xxxxx, tg_yy_xxxxy, tg_yy_xxxxz, tg_yy_xxxyy, tg_yy_xxxyz, tg_yy_xxxzz, tg_yy_xxyyy, tg_yy_xxyyz, tg_yy_xxyzz, tg_yy_xxzzz, tg_yy_xyyyy, tg_yy_xyyyz, tg_yy_xyyzz, tg_yy_xyzzz, tg_yy_xzzzz, tg_yy_yyyyy, tg_yy_yyyyz, tg_yy_yyyzz, tg_yy_yyzzz, tg_yy_yzzzz, tg_yy_zzzzz, tg_yz_xxxxx, tg_yz_xxxxy, tg_yz_xxxxz, tg_yz_xxxyy, tg_yz_xxxyz, tg_yz_xxxzz, tg_yz_xxyyy, tg_yz_xxyyz, tg_yz_xxyzz, tg_yz_xxzzz, tg_yz_xyyyy, tg_yz_xyyyz, tg_yz_xyyzz, tg_yz_xyzzz, tg_yz_xzzzz, tg_yz_yyyyy, tg_yz_yyyyz, tg_yz_yyyzz, tg_yz_yyzzz, tg_yz_yzzzz, tg_yz_zzzzz, tg_z_xxxx, tg_z_xxxxx, tg_z_xxxxy, tg_z_xxxxz, tg_z_xxxy, tg_z_xxxyy, tg_z_xxxyz, tg_z_xxxz, tg_z_xxxzz, tg_z_xxyy, tg_z_xxyyy, tg_z_xxyyz, tg_z_xxyz, tg_z_xxyzz, tg_z_xxzz, tg_z_xxzzz, tg_z_xyyy, tg_z_xyyyy, tg_z_xyyyz, tg_z_xyyz, tg_z_xyyzz, tg_z_xyzz, tg_z_xyzzz, tg_z_xzzz, tg_z_xzzzz, tg_z_yyyy, tg_z_yyyyy, tg_z_yyyyz, tg_z_yyyz, tg_z_yyyzz, tg_z_yyzz, tg_z_yyzzz, tg_z_yzzz, tg_z_yzzzz, tg_z_zzzz, tg_z_zzzzz, tg_zz_xxxxx, tg_zz_xxxxy, tg_zz_xxxxz, tg_zz_xxxyy, tg_zz_xxxyz, tg_zz_xxxzz, tg_zz_xxyyy, tg_zz_xxyyz, tg_zz_xxyzz, tg_zz_xxzzz, tg_zz_xyyyy, tg_zz_xyyyz, tg_zz_xyyzz, tg_zz_xyzzz, tg_zz_xzzzz, tg_zz_yyyyy, tg_zz_yyyyz, tg_zz_yyyzz, tg_zz_yyzzz, tg_zz_yzzzz, tg_zz_zzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xx_xxxxx[i] = tg_0_xxxxx[i] * fxi[i] + 5.0 * tg_x_xxxx[i] * fxi[i] + tg_x_xxxxx[i] * ra_x[i];

        tg_xx_xxxxy[i] = tg_0_xxxxy[i] * fxi[i] + 4.0 * tg_x_xxxy[i] * fxi[i] + tg_x_xxxxy[i] * ra_x[i];

        tg_xx_xxxxz[i] = tg_0_xxxxz[i] * fxi[i] + 4.0 * tg_x_xxxz[i] * fxi[i] + tg_x_xxxxz[i] * ra_x[i];

        tg_xx_xxxyy[i] = tg_0_xxxyy[i] * fxi[i] + 3.0 * tg_x_xxyy[i] * fxi[i] + tg_x_xxxyy[i] * ra_x[i];

        tg_xx_xxxyz[i] = tg_0_xxxyz[i] * fxi[i] + 3.0 * tg_x_xxyz[i] * fxi[i] + tg_x_xxxyz[i] * ra_x[i];

        tg_xx_xxxzz[i] = tg_0_xxxzz[i] * fxi[i] + 3.0 * tg_x_xxzz[i] * fxi[i] + tg_x_xxxzz[i] * ra_x[i];

        tg_xx_xxyyy[i] = tg_0_xxyyy[i] * fxi[i] + 2.0 * tg_x_xyyy[i] * fxi[i] + tg_x_xxyyy[i] * ra_x[i];

        tg_xx_xxyyz[i] = tg_0_xxyyz[i] * fxi[i] + 2.0 * tg_x_xyyz[i] * fxi[i] + tg_x_xxyyz[i] * ra_x[i];

        tg_xx_xxyzz[i] = tg_0_xxyzz[i] * fxi[i] + 2.0 * tg_x_xyzz[i] * fxi[i] + tg_x_xxyzz[i] * ra_x[i];

        tg_xx_xxzzz[i] = tg_0_xxzzz[i] * fxi[i] + 2.0 * tg_x_xzzz[i] * fxi[i] + tg_x_xxzzz[i] * ra_x[i];

        tg_xx_xyyyy[i] = tg_0_xyyyy[i] * fxi[i] + tg_x_yyyy[i] * fxi[i] + tg_x_xyyyy[i] * ra_x[i];

        tg_xx_xyyyz[i] = tg_0_xyyyz[i] * fxi[i] + tg_x_yyyz[i] * fxi[i] + tg_x_xyyyz[i] * ra_x[i];

        tg_xx_xyyzz[i] = tg_0_xyyzz[i] * fxi[i] + tg_x_yyzz[i] * fxi[i] + tg_x_xyyzz[i] * ra_x[i];

        tg_xx_xyzzz[i] = tg_0_xyzzz[i] * fxi[i] + tg_x_yzzz[i] * fxi[i] + tg_x_xyzzz[i] * ra_x[i];

        tg_xx_xzzzz[i] = tg_0_xzzzz[i] * fxi[i] + tg_x_zzzz[i] * fxi[i] + tg_x_xzzzz[i] * ra_x[i];

        tg_xx_yyyyy[i] = tg_0_yyyyy[i] * fxi[i] + tg_x_yyyyy[i] * ra_x[i];

        tg_xx_yyyyz[i] = tg_0_yyyyz[i] * fxi[i] + tg_x_yyyyz[i] * ra_x[i];

        tg_xx_yyyzz[i] = tg_0_yyyzz[i] * fxi[i] + tg_x_yyyzz[i] * ra_x[i];

        tg_xx_yyzzz[i] = tg_0_yyzzz[i] * fxi[i] + tg_x_yyzzz[i] * ra_x[i];

        tg_xx_yzzzz[i] = tg_0_yzzzz[i] * fxi[i] + tg_x_yzzzz[i] * ra_x[i];

        tg_xx_zzzzz[i] = tg_0_zzzzz[i] * fxi[i] + tg_x_zzzzz[i] * ra_x[i];

        tg_xy_xxxxx[i] = tg_x_xxxxx[i] * ra_y[i];

        tg_xy_xxxxy[i] = 4.0 * tg_y_xxxy[i] * fxi[i] + tg_y_xxxxy[i] * ra_x[i];

        tg_xy_xxxxz[i] = tg_x_xxxxz[i] * ra_y[i];

        tg_xy_xxxyy[i] = 3.0 * tg_y_xxyy[i] * fxi[i] + tg_y_xxxyy[i] * ra_x[i];

        tg_xy_xxxyz[i] = 3.0 * tg_y_xxyz[i] * fxi[i] + tg_y_xxxyz[i] * ra_x[i];

        tg_xy_xxxzz[i] = tg_x_xxxzz[i] * ra_y[i];

        tg_xy_xxyyy[i] = 2.0 * tg_y_xyyy[i] * fxi[i] + tg_y_xxyyy[i] * ra_x[i];

        tg_xy_xxyyz[i] = 2.0 * tg_y_xyyz[i] * fxi[i] + tg_y_xxyyz[i] * ra_x[i];

        tg_xy_xxyzz[i] = 2.0 * tg_y_xyzz[i] * fxi[i] + tg_y_xxyzz[i] * ra_x[i];

        tg_xy_xxzzz[i] = tg_x_xxzzz[i] * ra_y[i];

        tg_xy_xyyyy[i] = tg_y_yyyy[i] * fxi[i] + tg_y_xyyyy[i] * ra_x[i];

        tg_xy_xyyyz[i] = tg_y_yyyz[i] * fxi[i] + tg_y_xyyyz[i] * ra_x[i];

        tg_xy_xyyzz[i] = tg_y_yyzz[i] * fxi[i] + tg_y_xyyzz[i] * ra_x[i];

        tg_xy_xyzzz[i] = tg_y_yzzz[i] * fxi[i] + tg_y_xyzzz[i] * ra_x[i];

        tg_xy_xzzzz[i] = tg_x_xzzzz[i] * ra_y[i];

        tg_xy_yyyyy[i] = tg_y_yyyyy[i] * ra_x[i];

        tg_xy_yyyyz[i] = tg_y_yyyyz[i] * ra_x[i];

        tg_xy_yyyzz[i] = tg_y_yyyzz[i] * ra_x[i];

        tg_xy_yyzzz[i] = tg_y_yyzzz[i] * ra_x[i];

        tg_xy_yzzzz[i] = tg_y_yzzzz[i] * ra_x[i];

        tg_xy_zzzzz[i] = tg_y_zzzzz[i] * ra_x[i];

        tg_xz_xxxxx[i] = tg_x_xxxxx[i] * ra_z[i];

        tg_xz_xxxxy[i] = tg_x_xxxxy[i] * ra_z[i];

        tg_xz_xxxxz[i] = 4.0 * tg_z_xxxz[i] * fxi[i] + tg_z_xxxxz[i] * ra_x[i];

        tg_xz_xxxyy[i] = tg_x_xxxyy[i] * ra_z[i];

        tg_xz_xxxyz[i] = 3.0 * tg_z_xxyz[i] * fxi[i] + tg_z_xxxyz[i] * ra_x[i];

        tg_xz_xxxzz[i] = 3.0 * tg_z_xxzz[i] * fxi[i] + tg_z_xxxzz[i] * ra_x[i];

        tg_xz_xxyyy[i] = tg_x_xxyyy[i] * ra_z[i];

        tg_xz_xxyyz[i] = 2.0 * tg_z_xyyz[i] * fxi[i] + tg_z_xxyyz[i] * ra_x[i];

        tg_xz_xxyzz[i] = 2.0 * tg_z_xyzz[i] * fxi[i] + tg_z_xxyzz[i] * ra_x[i];

        tg_xz_xxzzz[i] = 2.0 * tg_z_xzzz[i] * fxi[i] + tg_z_xxzzz[i] * ra_x[i];

        tg_xz_xyyyy[i] = tg_x_xyyyy[i] * ra_z[i];

        tg_xz_xyyyz[i] = tg_z_yyyz[i] * fxi[i] + tg_z_xyyyz[i] * ra_x[i];

        tg_xz_xyyzz[i] = tg_z_yyzz[i] * fxi[i] + tg_z_xyyzz[i] * ra_x[i];

        tg_xz_xyzzz[i] = tg_z_yzzz[i] * fxi[i] + tg_z_xyzzz[i] * ra_x[i];

        tg_xz_xzzzz[i] = tg_z_zzzz[i] * fxi[i] + tg_z_xzzzz[i] * ra_x[i];

        tg_xz_yyyyy[i] = tg_z_yyyyy[i] * ra_x[i];

        tg_xz_yyyyz[i] = tg_z_yyyyz[i] * ra_x[i];

        tg_xz_yyyzz[i] = tg_z_yyyzz[i] * ra_x[i];

        tg_xz_yyzzz[i] = tg_z_yyzzz[i] * ra_x[i];

        tg_xz_yzzzz[i] = tg_z_yzzzz[i] * ra_x[i];

        tg_xz_zzzzz[i] = tg_z_zzzzz[i] * ra_x[i];

        tg_yy_xxxxx[i] = tg_0_xxxxx[i] * fxi[i] + tg_y_xxxxx[i] * ra_y[i];

        tg_yy_xxxxy[i] = tg_0_xxxxy[i] * fxi[i] + tg_y_xxxx[i] * fxi[i] + tg_y_xxxxy[i] * ra_y[i];

        tg_yy_xxxxz[i] = tg_0_xxxxz[i] * fxi[i] + tg_y_xxxxz[i] * ra_y[i];

        tg_yy_xxxyy[i] = tg_0_xxxyy[i] * fxi[i] + 2.0 * tg_y_xxxy[i] * fxi[i] + tg_y_xxxyy[i] * ra_y[i];

        tg_yy_xxxyz[i] = tg_0_xxxyz[i] * fxi[i] + tg_y_xxxz[i] * fxi[i] + tg_y_xxxyz[i] * ra_y[i];

        tg_yy_xxxzz[i] = tg_0_xxxzz[i] * fxi[i] + tg_y_xxxzz[i] * ra_y[i];

        tg_yy_xxyyy[i] = tg_0_xxyyy[i] * fxi[i] + 3.0 * tg_y_xxyy[i] * fxi[i] + tg_y_xxyyy[i] * ra_y[i];

        tg_yy_xxyyz[i] = tg_0_xxyyz[i] * fxi[i] + 2.0 * tg_y_xxyz[i] * fxi[i] + tg_y_xxyyz[i] * ra_y[i];

        tg_yy_xxyzz[i] = tg_0_xxyzz[i] * fxi[i] + tg_y_xxzz[i] * fxi[i] + tg_y_xxyzz[i] * ra_y[i];

        tg_yy_xxzzz[i] = tg_0_xxzzz[i] * fxi[i] + tg_y_xxzzz[i] * ra_y[i];

        tg_yy_xyyyy[i] = tg_0_xyyyy[i] * fxi[i] + 4.0 * tg_y_xyyy[i] * fxi[i] + tg_y_xyyyy[i] * ra_y[i];

        tg_yy_xyyyz[i] = tg_0_xyyyz[i] * fxi[i] + 3.0 * tg_y_xyyz[i] * fxi[i] + tg_y_xyyyz[i] * ra_y[i];

        tg_yy_xyyzz[i] = tg_0_xyyzz[i] * fxi[i] + 2.0 * tg_y_xyzz[i] * fxi[i] + tg_y_xyyzz[i] * ra_y[i];

        tg_yy_xyzzz[i] = tg_0_xyzzz[i] * fxi[i] + tg_y_xzzz[i] * fxi[i] + tg_y_xyzzz[i] * ra_y[i];

        tg_yy_xzzzz[i] = tg_0_xzzzz[i] * fxi[i] + tg_y_xzzzz[i] * ra_y[i];

        tg_yy_yyyyy[i] = tg_0_yyyyy[i] * fxi[i] + 5.0 * tg_y_yyyy[i] * fxi[i] + tg_y_yyyyy[i] * ra_y[i];

        tg_yy_yyyyz[i] = tg_0_yyyyz[i] * fxi[i] + 4.0 * tg_y_yyyz[i] * fxi[i] + tg_y_yyyyz[i] * ra_y[i];

        tg_yy_yyyzz[i] = tg_0_yyyzz[i] * fxi[i] + 3.0 * tg_y_yyzz[i] * fxi[i] + tg_y_yyyzz[i] * ra_y[i];

        tg_yy_yyzzz[i] = tg_0_yyzzz[i] * fxi[i] + 2.0 * tg_y_yzzz[i] * fxi[i] + tg_y_yyzzz[i] * ra_y[i];

        tg_yy_yzzzz[i] = tg_0_yzzzz[i] * fxi[i] + tg_y_zzzz[i] * fxi[i] + tg_y_yzzzz[i] * ra_y[i];

        tg_yy_zzzzz[i] = tg_0_zzzzz[i] * fxi[i] + tg_y_zzzzz[i] * ra_y[i];

        tg_yz_xxxxx[i] = tg_z_xxxxx[i] * ra_y[i];

        tg_yz_xxxxy[i] = tg_y_xxxxy[i] * ra_z[i];

        tg_yz_xxxxz[i] = tg_z_xxxxz[i] * ra_y[i];

        tg_yz_xxxyy[i] = tg_y_xxxyy[i] * ra_z[i];

        tg_yz_xxxyz[i] = tg_z_xxxz[i] * fxi[i] + tg_z_xxxyz[i] * ra_y[i];

        tg_yz_xxxzz[i] = tg_z_xxxzz[i] * ra_y[i];

        tg_yz_xxyyy[i] = tg_y_xxyyy[i] * ra_z[i];

        tg_yz_xxyyz[i] = 2.0 * tg_z_xxyz[i] * fxi[i] + tg_z_xxyyz[i] * ra_y[i];

        tg_yz_xxyzz[i] = tg_z_xxzz[i] * fxi[i] + tg_z_xxyzz[i] * ra_y[i];

        tg_yz_xxzzz[i] = tg_z_xxzzz[i] * ra_y[i];

        tg_yz_xyyyy[i] = tg_y_xyyyy[i] * ra_z[i];

        tg_yz_xyyyz[i] = 3.0 * tg_z_xyyz[i] * fxi[i] + tg_z_xyyyz[i] * ra_y[i];

        tg_yz_xyyzz[i] = 2.0 * tg_z_xyzz[i] * fxi[i] + tg_z_xyyzz[i] * ra_y[i];

        tg_yz_xyzzz[i] = tg_z_xzzz[i] * fxi[i] + tg_z_xyzzz[i] * ra_y[i];

        tg_yz_xzzzz[i] = tg_z_xzzzz[i] * ra_y[i];

        tg_yz_yyyyy[i] = tg_y_yyyyy[i] * ra_z[i];

        tg_yz_yyyyz[i] = 4.0 * tg_z_yyyz[i] * fxi[i] + tg_z_yyyyz[i] * ra_y[i];

        tg_yz_yyyzz[i] = 3.0 * tg_z_yyzz[i] * fxi[i] + tg_z_yyyzz[i] * ra_y[i];

        tg_yz_yyzzz[i] = 2.0 * tg_z_yzzz[i] * fxi[i] + tg_z_yyzzz[i] * ra_y[i];

        tg_yz_yzzzz[i] = tg_z_zzzz[i] * fxi[i] + tg_z_yzzzz[i] * ra_y[i];

        tg_yz_zzzzz[i] = tg_z_zzzzz[i] * ra_y[i];

        tg_zz_xxxxx[i] = tg_0_xxxxx[i] * fxi[i] + tg_z_xxxxx[i] * ra_z[i];

        tg_zz_xxxxy[i] = tg_0_xxxxy[i] * fxi[i] + tg_z_xxxxy[i] * ra_z[i];

        tg_zz_xxxxz[i] = tg_0_xxxxz[i] * fxi[i] + tg_z_xxxx[i] * fxi[i] + tg_z_xxxxz[i] * ra_z[i];

        tg_zz_xxxyy[i] = tg_0_xxxyy[i] * fxi[i] + tg_z_xxxyy[i] * ra_z[i];

        tg_zz_xxxyz[i] = tg_0_xxxyz[i] * fxi[i] + tg_z_xxxy[i] * fxi[i] + tg_z_xxxyz[i] * ra_z[i];

        tg_zz_xxxzz[i] = tg_0_xxxzz[i] * fxi[i] + 2.0 * tg_z_xxxz[i] * fxi[i] + tg_z_xxxzz[i] * ra_z[i];

        tg_zz_xxyyy[i] = tg_0_xxyyy[i] * fxi[i] + tg_z_xxyyy[i] * ra_z[i];

        tg_zz_xxyyz[i] = tg_0_xxyyz[i] * fxi[i] + tg_z_xxyy[i] * fxi[i] + tg_z_xxyyz[i] * ra_z[i];

        tg_zz_xxyzz[i] = tg_0_xxyzz[i] * fxi[i] + 2.0 * tg_z_xxyz[i] * fxi[i] + tg_z_xxyzz[i] * ra_z[i];

        tg_zz_xxzzz[i] = tg_0_xxzzz[i] * fxi[i] + 3.0 * tg_z_xxzz[i] * fxi[i] + tg_z_xxzzz[i] * ra_z[i];

        tg_zz_xyyyy[i] = tg_0_xyyyy[i] * fxi[i] + tg_z_xyyyy[i] * ra_z[i];

        tg_zz_xyyyz[i] = tg_0_xyyyz[i] * fxi[i] + tg_z_xyyy[i] * fxi[i] + tg_z_xyyyz[i] * ra_z[i];

        tg_zz_xyyzz[i] = tg_0_xyyzz[i] * fxi[i] + 2.0 * tg_z_xyyz[i] * fxi[i] + tg_z_xyyzz[i] * ra_z[i];

        tg_zz_xyzzz[i] = tg_0_xyzzz[i] * fxi[i] + 3.0 * tg_z_xyzz[i] * fxi[i] + tg_z_xyzzz[i] * ra_z[i];

        tg_zz_xzzzz[i] = tg_0_xzzzz[i] * fxi[i] + 4.0 * tg_z_xzzz[i] * fxi[i] + tg_z_xzzzz[i] * ra_z[i];

        tg_zz_yyyyy[i] = tg_0_yyyyy[i] * fxi[i] + tg_z_yyyyy[i] * ra_z[i];

        tg_zz_yyyyz[i] = tg_0_yyyyz[i] * fxi[i] + tg_z_yyyy[i] * fxi[i] + tg_z_yyyyz[i] * ra_z[i];

        tg_zz_yyyzz[i] = tg_0_yyyzz[i] * fxi[i] + 2.0 * tg_z_yyyz[i] * fxi[i] + tg_z_yyyzz[i] * ra_z[i];

        tg_zz_yyzzz[i] = tg_0_yyzzz[i] * fxi[i] + 3.0 * tg_z_yyzz[i] * fxi[i] + tg_z_yyzzz[i] * ra_z[i];

        tg_zz_yzzzz[i] = tg_0_yzzzz[i] * fxi[i] + 4.0 * tg_z_yzzz[i] * fxi[i] + tg_z_yzzzz[i] * ra_z[i];

        tg_zz_zzzzz[i] = tg_0_zzzzz[i] * fxi[i] + 5.0 * tg_z_zzzz[i] * fxi[i] + tg_z_zzzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

