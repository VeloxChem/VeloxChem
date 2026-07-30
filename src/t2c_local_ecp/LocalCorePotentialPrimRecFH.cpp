#include "LocalCorePotentialPrimRecFH.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_fh(CSimdArray<double>& pbuffer, 
                                  const size_t idx_fh,
                                  const size_t idx_ph,
                                  const size_t idx_dg,
                                  const size_t idx_dh,
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

    // Set up components of auxiliary buffer : DG

    auto tg_xx_xxxx = pbuffer.data(idx_dg);

    auto tg_xx_xxxy = pbuffer.data(idx_dg + 1);

    auto tg_xx_xxxz = pbuffer.data(idx_dg + 2);

    auto tg_xx_xxyy = pbuffer.data(idx_dg + 3);

    auto tg_xx_xxyz = pbuffer.data(idx_dg + 4);

    auto tg_xx_xxzz = pbuffer.data(idx_dg + 5);

    auto tg_xx_xyyy = pbuffer.data(idx_dg + 6);

    auto tg_xx_xyyz = pbuffer.data(idx_dg + 7);

    auto tg_xx_xyzz = pbuffer.data(idx_dg + 8);

    auto tg_xx_xzzz = pbuffer.data(idx_dg + 9);

    auto tg_xx_yyyy = pbuffer.data(idx_dg + 10);

    auto tg_xx_yyyz = pbuffer.data(idx_dg + 11);

    auto tg_xx_yyzz = pbuffer.data(idx_dg + 12);

    auto tg_xx_yzzz = pbuffer.data(idx_dg + 13);

    auto tg_xx_zzzz = pbuffer.data(idx_dg + 14);

    auto tg_yy_xxxx = pbuffer.data(idx_dg + 45);

    auto tg_yy_xxxy = pbuffer.data(idx_dg + 46);

    auto tg_yy_xxxz = pbuffer.data(idx_dg + 47);

    auto tg_yy_xxyy = pbuffer.data(idx_dg + 48);

    auto tg_yy_xxyz = pbuffer.data(idx_dg + 49);

    auto tg_yy_xxzz = pbuffer.data(idx_dg + 50);

    auto tg_yy_xyyy = pbuffer.data(idx_dg + 51);

    auto tg_yy_xyyz = pbuffer.data(idx_dg + 52);

    auto tg_yy_xyzz = pbuffer.data(idx_dg + 53);

    auto tg_yy_xzzz = pbuffer.data(idx_dg + 54);

    auto tg_yy_yyyy = pbuffer.data(idx_dg + 55);

    auto tg_yy_yyyz = pbuffer.data(idx_dg + 56);

    auto tg_yy_yyzz = pbuffer.data(idx_dg + 57);

    auto tg_yy_yzzz = pbuffer.data(idx_dg + 58);

    auto tg_yy_zzzz = pbuffer.data(idx_dg + 59);

    auto tg_yz_xxyz = pbuffer.data(idx_dg + 64);

    auto tg_yz_xyyz = pbuffer.data(idx_dg + 67);

    auto tg_yz_xyzz = pbuffer.data(idx_dg + 68);

    auto tg_yz_yyyz = pbuffer.data(idx_dg + 71);

    auto tg_yz_yyzz = pbuffer.data(idx_dg + 72);

    auto tg_yz_yzzz = pbuffer.data(idx_dg + 73);

    auto tg_zz_xxxx = pbuffer.data(idx_dg + 75);

    auto tg_zz_xxxy = pbuffer.data(idx_dg + 76);

    auto tg_zz_xxxz = pbuffer.data(idx_dg + 77);

    auto tg_zz_xxyy = pbuffer.data(idx_dg + 78);

    auto tg_zz_xxyz = pbuffer.data(idx_dg + 79);

    auto tg_zz_xxzz = pbuffer.data(idx_dg + 80);

    auto tg_zz_xyyy = pbuffer.data(idx_dg + 81);

    auto tg_zz_xyyz = pbuffer.data(idx_dg + 82);

    auto tg_zz_xyzz = pbuffer.data(idx_dg + 83);

    auto tg_zz_xzzz = pbuffer.data(idx_dg + 84);

    auto tg_zz_yyyy = pbuffer.data(idx_dg + 85);

    auto tg_zz_yyyz = pbuffer.data(idx_dg + 86);

    auto tg_zz_yyzz = pbuffer.data(idx_dg + 87);

    auto tg_zz_yzzz = pbuffer.data(idx_dg + 88);

    auto tg_zz_zzzz = pbuffer.data(idx_dg + 89);

    // Set up components of auxiliary buffer : DH

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

    auto tg_xy_xxxxy = pbuffer.data(idx_dh + 22);

    auto tg_xy_xxxyy = pbuffer.data(idx_dh + 24);

    auto tg_xy_xxyyy = pbuffer.data(idx_dh + 27);

    auto tg_xy_xyyyy = pbuffer.data(idx_dh + 31);

    auto tg_xy_yyyyy = pbuffer.data(idx_dh + 36);

    auto tg_xy_yyyyz = pbuffer.data(idx_dh + 37);

    auto tg_xy_yyyzz = pbuffer.data(idx_dh + 38);

    auto tg_xy_yyzzz = pbuffer.data(idx_dh + 39);

    auto tg_xy_yzzzz = pbuffer.data(idx_dh + 40);

    auto tg_xz_xxxxx = pbuffer.data(idx_dh + 42);

    auto tg_xz_xxxxz = pbuffer.data(idx_dh + 44);

    auto tg_xz_xxxzz = pbuffer.data(idx_dh + 47);

    auto tg_xz_xxzzz = pbuffer.data(idx_dh + 51);

    auto tg_xz_xzzzz = pbuffer.data(idx_dh + 56);

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

    auto tg_yz_xxxxz = pbuffer.data(idx_dh + 86);

    auto tg_yz_xxxyz = pbuffer.data(idx_dh + 88);

    auto tg_yz_xxxzz = pbuffer.data(idx_dh + 89);

    auto tg_yz_xxyyz = pbuffer.data(idx_dh + 91);

    auto tg_yz_xxyzz = pbuffer.data(idx_dh + 92);

    auto tg_yz_xxzzz = pbuffer.data(idx_dh + 93);

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

    // Set up components of targeted buffer : FH

    auto tg_xxx_xxxxx = pbuffer.data(idx_fh);

    auto tg_xxx_xxxxy = pbuffer.data(idx_fh + 1);

    auto tg_xxx_xxxxz = pbuffer.data(idx_fh + 2);

    auto tg_xxx_xxxyy = pbuffer.data(idx_fh + 3);

    auto tg_xxx_xxxyz = pbuffer.data(idx_fh + 4);

    auto tg_xxx_xxxzz = pbuffer.data(idx_fh + 5);

    auto tg_xxx_xxyyy = pbuffer.data(idx_fh + 6);

    auto tg_xxx_xxyyz = pbuffer.data(idx_fh + 7);

    auto tg_xxx_xxyzz = pbuffer.data(idx_fh + 8);

    auto tg_xxx_xxzzz = pbuffer.data(idx_fh + 9);

    auto tg_xxx_xyyyy = pbuffer.data(idx_fh + 10);

    auto tg_xxx_xyyyz = pbuffer.data(idx_fh + 11);

    auto tg_xxx_xyyzz = pbuffer.data(idx_fh + 12);

    auto tg_xxx_xyzzz = pbuffer.data(idx_fh + 13);

    auto tg_xxx_xzzzz = pbuffer.data(idx_fh + 14);

    auto tg_xxx_yyyyy = pbuffer.data(idx_fh + 15);

    auto tg_xxx_yyyyz = pbuffer.data(idx_fh + 16);

    auto tg_xxx_yyyzz = pbuffer.data(idx_fh + 17);

    auto tg_xxx_yyzzz = pbuffer.data(idx_fh + 18);

    auto tg_xxx_yzzzz = pbuffer.data(idx_fh + 19);

    auto tg_xxx_zzzzz = pbuffer.data(idx_fh + 20);

    auto tg_xxy_xxxxx = pbuffer.data(idx_fh + 21);

    auto tg_xxy_xxxxy = pbuffer.data(idx_fh + 22);

    auto tg_xxy_xxxxz = pbuffer.data(idx_fh + 23);

    auto tg_xxy_xxxyy = pbuffer.data(idx_fh + 24);

    auto tg_xxy_xxxyz = pbuffer.data(idx_fh + 25);

    auto tg_xxy_xxxzz = pbuffer.data(idx_fh + 26);

    auto tg_xxy_xxyyy = pbuffer.data(idx_fh + 27);

    auto tg_xxy_xxyyz = pbuffer.data(idx_fh + 28);

    auto tg_xxy_xxyzz = pbuffer.data(idx_fh + 29);

    auto tg_xxy_xxzzz = pbuffer.data(idx_fh + 30);

    auto tg_xxy_xyyyy = pbuffer.data(idx_fh + 31);

    auto tg_xxy_xyyyz = pbuffer.data(idx_fh + 32);

    auto tg_xxy_xyyzz = pbuffer.data(idx_fh + 33);

    auto tg_xxy_xyzzz = pbuffer.data(idx_fh + 34);

    auto tg_xxy_xzzzz = pbuffer.data(idx_fh + 35);

    auto tg_xxy_yyyyy = pbuffer.data(idx_fh + 36);

    auto tg_xxy_yyyyz = pbuffer.data(idx_fh + 37);

    auto tg_xxy_yyyzz = pbuffer.data(idx_fh + 38);

    auto tg_xxy_yyzzz = pbuffer.data(idx_fh + 39);

    auto tg_xxy_yzzzz = pbuffer.data(idx_fh + 40);

    auto tg_xxy_zzzzz = pbuffer.data(idx_fh + 41);

    auto tg_xxz_xxxxx = pbuffer.data(idx_fh + 42);

    auto tg_xxz_xxxxy = pbuffer.data(idx_fh + 43);

    auto tg_xxz_xxxxz = pbuffer.data(idx_fh + 44);

    auto tg_xxz_xxxyy = pbuffer.data(idx_fh + 45);

    auto tg_xxz_xxxyz = pbuffer.data(idx_fh + 46);

    auto tg_xxz_xxxzz = pbuffer.data(idx_fh + 47);

    auto tg_xxz_xxyyy = pbuffer.data(idx_fh + 48);

    auto tg_xxz_xxyyz = pbuffer.data(idx_fh + 49);

    auto tg_xxz_xxyzz = pbuffer.data(idx_fh + 50);

    auto tg_xxz_xxzzz = pbuffer.data(idx_fh + 51);

    auto tg_xxz_xyyyy = pbuffer.data(idx_fh + 52);

    auto tg_xxz_xyyyz = pbuffer.data(idx_fh + 53);

    auto tg_xxz_xyyzz = pbuffer.data(idx_fh + 54);

    auto tg_xxz_xyzzz = pbuffer.data(idx_fh + 55);

    auto tg_xxz_xzzzz = pbuffer.data(idx_fh + 56);

    auto tg_xxz_yyyyy = pbuffer.data(idx_fh + 57);

    auto tg_xxz_yyyyz = pbuffer.data(idx_fh + 58);

    auto tg_xxz_yyyzz = pbuffer.data(idx_fh + 59);

    auto tg_xxz_yyzzz = pbuffer.data(idx_fh + 60);

    auto tg_xxz_yzzzz = pbuffer.data(idx_fh + 61);

    auto tg_xxz_zzzzz = pbuffer.data(idx_fh + 62);

    auto tg_xyy_xxxxx = pbuffer.data(idx_fh + 63);

    auto tg_xyy_xxxxy = pbuffer.data(idx_fh + 64);

    auto tg_xyy_xxxxz = pbuffer.data(idx_fh + 65);

    auto tg_xyy_xxxyy = pbuffer.data(idx_fh + 66);

    auto tg_xyy_xxxyz = pbuffer.data(idx_fh + 67);

    auto tg_xyy_xxxzz = pbuffer.data(idx_fh + 68);

    auto tg_xyy_xxyyy = pbuffer.data(idx_fh + 69);

    auto tg_xyy_xxyyz = pbuffer.data(idx_fh + 70);

    auto tg_xyy_xxyzz = pbuffer.data(idx_fh + 71);

    auto tg_xyy_xxzzz = pbuffer.data(idx_fh + 72);

    auto tg_xyy_xyyyy = pbuffer.data(idx_fh + 73);

    auto tg_xyy_xyyyz = pbuffer.data(idx_fh + 74);

    auto tg_xyy_xyyzz = pbuffer.data(idx_fh + 75);

    auto tg_xyy_xyzzz = pbuffer.data(idx_fh + 76);

    auto tg_xyy_xzzzz = pbuffer.data(idx_fh + 77);

    auto tg_xyy_yyyyy = pbuffer.data(idx_fh + 78);

    auto tg_xyy_yyyyz = pbuffer.data(idx_fh + 79);

    auto tg_xyy_yyyzz = pbuffer.data(idx_fh + 80);

    auto tg_xyy_yyzzz = pbuffer.data(idx_fh + 81);

    auto tg_xyy_yzzzz = pbuffer.data(idx_fh + 82);

    auto tg_xyy_zzzzz = pbuffer.data(idx_fh + 83);

    auto tg_xyz_xxxxx = pbuffer.data(idx_fh + 84);

    auto tg_xyz_xxxxy = pbuffer.data(idx_fh + 85);

    auto tg_xyz_xxxxz = pbuffer.data(idx_fh + 86);

    auto tg_xyz_xxxyy = pbuffer.data(idx_fh + 87);

    auto tg_xyz_xxxyz = pbuffer.data(idx_fh + 88);

    auto tg_xyz_xxxzz = pbuffer.data(idx_fh + 89);

    auto tg_xyz_xxyyy = pbuffer.data(idx_fh + 90);

    auto tg_xyz_xxyyz = pbuffer.data(idx_fh + 91);

    auto tg_xyz_xxyzz = pbuffer.data(idx_fh + 92);

    auto tg_xyz_xxzzz = pbuffer.data(idx_fh + 93);

    auto tg_xyz_xyyyy = pbuffer.data(idx_fh + 94);

    auto tg_xyz_xyyyz = pbuffer.data(idx_fh + 95);

    auto tg_xyz_xyyzz = pbuffer.data(idx_fh + 96);

    auto tg_xyz_xyzzz = pbuffer.data(idx_fh + 97);

    auto tg_xyz_xzzzz = pbuffer.data(idx_fh + 98);

    auto tg_xyz_yyyyy = pbuffer.data(idx_fh + 99);

    auto tg_xyz_yyyyz = pbuffer.data(idx_fh + 100);

    auto tg_xyz_yyyzz = pbuffer.data(idx_fh + 101);

    auto tg_xyz_yyzzz = pbuffer.data(idx_fh + 102);

    auto tg_xyz_yzzzz = pbuffer.data(idx_fh + 103);

    auto tg_xyz_zzzzz = pbuffer.data(idx_fh + 104);

    auto tg_xzz_xxxxx = pbuffer.data(idx_fh + 105);

    auto tg_xzz_xxxxy = pbuffer.data(idx_fh + 106);

    auto tg_xzz_xxxxz = pbuffer.data(idx_fh + 107);

    auto tg_xzz_xxxyy = pbuffer.data(idx_fh + 108);

    auto tg_xzz_xxxyz = pbuffer.data(idx_fh + 109);

    auto tg_xzz_xxxzz = pbuffer.data(idx_fh + 110);

    auto tg_xzz_xxyyy = pbuffer.data(idx_fh + 111);

    auto tg_xzz_xxyyz = pbuffer.data(idx_fh + 112);

    auto tg_xzz_xxyzz = pbuffer.data(idx_fh + 113);

    auto tg_xzz_xxzzz = pbuffer.data(idx_fh + 114);

    auto tg_xzz_xyyyy = pbuffer.data(idx_fh + 115);

    auto tg_xzz_xyyyz = pbuffer.data(idx_fh + 116);

    auto tg_xzz_xyyzz = pbuffer.data(idx_fh + 117);

    auto tg_xzz_xyzzz = pbuffer.data(idx_fh + 118);

    auto tg_xzz_xzzzz = pbuffer.data(idx_fh + 119);

    auto tg_xzz_yyyyy = pbuffer.data(idx_fh + 120);

    auto tg_xzz_yyyyz = pbuffer.data(idx_fh + 121);

    auto tg_xzz_yyyzz = pbuffer.data(idx_fh + 122);

    auto tg_xzz_yyzzz = pbuffer.data(idx_fh + 123);

    auto tg_xzz_yzzzz = pbuffer.data(idx_fh + 124);

    auto tg_xzz_zzzzz = pbuffer.data(idx_fh + 125);

    auto tg_yyy_xxxxx = pbuffer.data(idx_fh + 126);

    auto tg_yyy_xxxxy = pbuffer.data(idx_fh + 127);

    auto tg_yyy_xxxxz = pbuffer.data(idx_fh + 128);

    auto tg_yyy_xxxyy = pbuffer.data(idx_fh + 129);

    auto tg_yyy_xxxyz = pbuffer.data(idx_fh + 130);

    auto tg_yyy_xxxzz = pbuffer.data(idx_fh + 131);

    auto tg_yyy_xxyyy = pbuffer.data(idx_fh + 132);

    auto tg_yyy_xxyyz = pbuffer.data(idx_fh + 133);

    auto tg_yyy_xxyzz = pbuffer.data(idx_fh + 134);

    auto tg_yyy_xxzzz = pbuffer.data(idx_fh + 135);

    auto tg_yyy_xyyyy = pbuffer.data(idx_fh + 136);

    auto tg_yyy_xyyyz = pbuffer.data(idx_fh + 137);

    auto tg_yyy_xyyzz = pbuffer.data(idx_fh + 138);

    auto tg_yyy_xyzzz = pbuffer.data(idx_fh + 139);

    auto tg_yyy_xzzzz = pbuffer.data(idx_fh + 140);

    auto tg_yyy_yyyyy = pbuffer.data(idx_fh + 141);

    auto tg_yyy_yyyyz = pbuffer.data(idx_fh + 142);

    auto tg_yyy_yyyzz = pbuffer.data(idx_fh + 143);

    auto tg_yyy_yyzzz = pbuffer.data(idx_fh + 144);

    auto tg_yyy_yzzzz = pbuffer.data(idx_fh + 145);

    auto tg_yyy_zzzzz = pbuffer.data(idx_fh + 146);

    auto tg_yyz_xxxxx = pbuffer.data(idx_fh + 147);

    auto tg_yyz_xxxxy = pbuffer.data(idx_fh + 148);

    auto tg_yyz_xxxxz = pbuffer.data(idx_fh + 149);

    auto tg_yyz_xxxyy = pbuffer.data(idx_fh + 150);

    auto tg_yyz_xxxyz = pbuffer.data(idx_fh + 151);

    auto tg_yyz_xxxzz = pbuffer.data(idx_fh + 152);

    auto tg_yyz_xxyyy = pbuffer.data(idx_fh + 153);

    auto tg_yyz_xxyyz = pbuffer.data(idx_fh + 154);

    auto tg_yyz_xxyzz = pbuffer.data(idx_fh + 155);

    auto tg_yyz_xxzzz = pbuffer.data(idx_fh + 156);

    auto tg_yyz_xyyyy = pbuffer.data(idx_fh + 157);

    auto tg_yyz_xyyyz = pbuffer.data(idx_fh + 158);

    auto tg_yyz_xyyzz = pbuffer.data(idx_fh + 159);

    auto tg_yyz_xyzzz = pbuffer.data(idx_fh + 160);

    auto tg_yyz_xzzzz = pbuffer.data(idx_fh + 161);

    auto tg_yyz_yyyyy = pbuffer.data(idx_fh + 162);

    auto tg_yyz_yyyyz = pbuffer.data(idx_fh + 163);

    auto tg_yyz_yyyzz = pbuffer.data(idx_fh + 164);

    auto tg_yyz_yyzzz = pbuffer.data(idx_fh + 165);

    auto tg_yyz_yzzzz = pbuffer.data(idx_fh + 166);

    auto tg_yyz_zzzzz = pbuffer.data(idx_fh + 167);

    auto tg_yzz_xxxxx = pbuffer.data(idx_fh + 168);

    auto tg_yzz_xxxxy = pbuffer.data(idx_fh + 169);

    auto tg_yzz_xxxxz = pbuffer.data(idx_fh + 170);

    auto tg_yzz_xxxyy = pbuffer.data(idx_fh + 171);

    auto tg_yzz_xxxyz = pbuffer.data(idx_fh + 172);

    auto tg_yzz_xxxzz = pbuffer.data(idx_fh + 173);

    auto tg_yzz_xxyyy = pbuffer.data(idx_fh + 174);

    auto tg_yzz_xxyyz = pbuffer.data(idx_fh + 175);

    auto tg_yzz_xxyzz = pbuffer.data(idx_fh + 176);

    auto tg_yzz_xxzzz = pbuffer.data(idx_fh + 177);

    auto tg_yzz_xyyyy = pbuffer.data(idx_fh + 178);

    auto tg_yzz_xyyyz = pbuffer.data(idx_fh + 179);

    auto tg_yzz_xyyzz = pbuffer.data(idx_fh + 180);

    auto tg_yzz_xyzzz = pbuffer.data(idx_fh + 181);

    auto tg_yzz_xzzzz = pbuffer.data(idx_fh + 182);

    auto tg_yzz_yyyyy = pbuffer.data(idx_fh + 183);

    auto tg_yzz_yyyyz = pbuffer.data(idx_fh + 184);

    auto tg_yzz_yyyzz = pbuffer.data(idx_fh + 185);

    auto tg_yzz_yyzzz = pbuffer.data(idx_fh + 186);

    auto tg_yzz_yzzzz = pbuffer.data(idx_fh + 187);

    auto tg_yzz_zzzzz = pbuffer.data(idx_fh + 188);

    auto tg_zzz_xxxxx = pbuffer.data(idx_fh + 189);

    auto tg_zzz_xxxxy = pbuffer.data(idx_fh + 190);

    auto tg_zzz_xxxxz = pbuffer.data(idx_fh + 191);

    auto tg_zzz_xxxyy = pbuffer.data(idx_fh + 192);

    auto tg_zzz_xxxyz = pbuffer.data(idx_fh + 193);

    auto tg_zzz_xxxzz = pbuffer.data(idx_fh + 194);

    auto tg_zzz_xxyyy = pbuffer.data(idx_fh + 195);

    auto tg_zzz_xxyyz = pbuffer.data(idx_fh + 196);

    auto tg_zzz_xxyzz = pbuffer.data(idx_fh + 197);

    auto tg_zzz_xxzzz = pbuffer.data(idx_fh + 198);

    auto tg_zzz_xyyyy = pbuffer.data(idx_fh + 199);

    auto tg_zzz_xyyyz = pbuffer.data(idx_fh + 200);

    auto tg_zzz_xyyzz = pbuffer.data(idx_fh + 201);

    auto tg_zzz_xyzzz = pbuffer.data(idx_fh + 202);

    auto tg_zzz_xzzzz = pbuffer.data(idx_fh + 203);

    auto tg_zzz_yyyyy = pbuffer.data(idx_fh + 204);

    auto tg_zzz_yyyyz = pbuffer.data(idx_fh + 205);

    auto tg_zzz_yyyzz = pbuffer.data(idx_fh + 206);

    auto tg_zzz_yyzzz = pbuffer.data(idx_fh + 207);

    auto tg_zzz_yzzzz = pbuffer.data(idx_fh + 208);

    auto tg_zzz_zzzzz = pbuffer.data(idx_fh + 209);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_x_xxxxx, tg_x_xxxxy, tg_x_xxxxz, tg_x_xxxyy, tg_x_xxxyz, tg_x_xxxzz, tg_x_xxyyy, tg_x_xxyyz, tg_x_xxyzz, tg_x_xxzzz, tg_x_xyyyy, tg_x_xyyyz, tg_x_xyyzz, tg_x_xyzzz, tg_x_xzzzz, tg_x_yyyyy, tg_x_yyyyz, tg_x_yyyzz, tg_x_yyzzz, tg_x_yzzzz, tg_x_zzzzz, tg_xx_xxxx, tg_xx_xxxxx, tg_xx_xxxxy, tg_xx_xxxxz, tg_xx_xxxy, tg_xx_xxxyy, tg_xx_xxxyz, tg_xx_xxxz, tg_xx_xxxzz, tg_xx_xxyy, tg_xx_xxyyy, tg_xx_xxyyz, tg_xx_xxyz, tg_xx_xxyzz, tg_xx_xxzz, tg_xx_xxzzz, tg_xx_xyyy, tg_xx_xyyyy, tg_xx_xyyyz, tg_xx_xyyz, tg_xx_xyyzz, tg_xx_xyzz, tg_xx_xyzzz, tg_xx_xzzz, tg_xx_xzzzz, tg_xx_yyyy, tg_xx_yyyyy, tg_xx_yyyyz, tg_xx_yyyz, tg_xx_yyyzz, tg_xx_yyzz, tg_xx_yyzzz, tg_xx_yzzz, tg_xx_yzzzz, tg_xx_zzzz, tg_xx_zzzzz, tg_xxx_xxxxx, tg_xxx_xxxxy, tg_xxx_xxxxz, tg_xxx_xxxyy, tg_xxx_xxxyz, tg_xxx_xxxzz, tg_xxx_xxyyy, tg_xxx_xxyyz, tg_xxx_xxyzz, tg_xxx_xxzzz, tg_xxx_xyyyy, tg_xxx_xyyyz, tg_xxx_xyyzz, tg_xxx_xyzzz, tg_xxx_xzzzz, tg_xxx_yyyyy, tg_xxx_yyyyz, tg_xxx_yyyzz, tg_xxx_yyzzz, tg_xxx_yzzzz, tg_xxx_zzzzz, tg_xxy_xxxxx, tg_xxy_xxxxy, tg_xxy_xxxxz, tg_xxy_xxxyy, tg_xxy_xxxyz, tg_xxy_xxxzz, tg_xxy_xxyyy, tg_xxy_xxyyz, tg_xxy_xxyzz, tg_xxy_xxzzz, tg_xxy_xyyyy, tg_xxy_xyyyz, tg_xxy_xyyzz, tg_xxy_xyzzz, tg_xxy_xzzzz, tg_xxy_yyyyy, tg_xxy_yyyyz, tg_xxy_yyyzz, tg_xxy_yyzzz, tg_xxy_yzzzz, tg_xxy_zzzzz, tg_xxz_xxxxx, tg_xxz_xxxxy, tg_xxz_xxxxz, tg_xxz_xxxyy, tg_xxz_xxxyz, tg_xxz_xxxzz, tg_xxz_xxyyy, tg_xxz_xxyyz, tg_xxz_xxyzz, tg_xxz_xxzzz, tg_xxz_xyyyy, tg_xxz_xyyyz, tg_xxz_xyyzz, tg_xxz_xyzzz, tg_xxz_xzzzz, tg_xxz_yyyyy, tg_xxz_yyyyz, tg_xxz_yyyzz, tg_xxz_yyzzz, tg_xxz_yzzzz, tg_xxz_zzzzz, tg_xy_xxxxy, tg_xy_xxxyy, tg_xy_xxyyy, tg_xy_xyyyy, tg_xy_yyyyy, tg_xy_yyyyz, tg_xy_yyyzz, tg_xy_yyzzz, tg_xy_yzzzz, tg_xyy_xxxxx, tg_xyy_xxxxy, tg_xyy_xxxxz, tg_xyy_xxxyy, tg_xyy_xxxyz, tg_xyy_xxxzz, tg_xyy_xxyyy, tg_xyy_xxyyz, tg_xyy_xxyzz, tg_xyy_xxzzz, tg_xyy_xyyyy, tg_xyy_xyyyz, tg_xyy_xyyzz, tg_xyy_xyzzz, tg_xyy_xzzzz, tg_xyy_yyyyy, tg_xyy_yyyyz, tg_xyy_yyyzz, tg_xyy_yyzzz, tg_xyy_yzzzz, tg_xyy_zzzzz, tg_xyz_xxxxx, tg_xyz_xxxxy, tg_xyz_xxxxz, tg_xyz_xxxyy, tg_xyz_xxxyz, tg_xyz_xxxzz, tg_xyz_xxyyy, tg_xyz_xxyyz, tg_xyz_xxyzz, tg_xyz_xxzzz, tg_xyz_xyyyy, tg_xyz_xyyyz, tg_xyz_xyyzz, tg_xyz_xyzzz, tg_xyz_xzzzz, tg_xyz_yyyyy, tg_xyz_yyyyz, tg_xyz_yyyzz, tg_xyz_yyzzz, tg_xyz_yzzzz, tg_xyz_zzzzz, tg_xz_xxxxx, tg_xz_xxxxz, tg_xz_xxxzz, tg_xz_xxzzz, tg_xz_xzzzz, tg_xz_yyyyz, tg_xz_yyyzz, tg_xz_yyzzz, tg_xz_yzzzz, tg_xz_zzzzz, tg_xzz_xxxxx, tg_xzz_xxxxy, tg_xzz_xxxxz, tg_xzz_xxxyy, tg_xzz_xxxyz, tg_xzz_xxxzz, tg_xzz_xxyyy, tg_xzz_xxyyz, tg_xzz_xxyzz, tg_xzz_xxzzz, tg_xzz_xyyyy, tg_xzz_xyyyz, tg_xzz_xyyzz, tg_xzz_xyzzz, tg_xzz_xzzzz, tg_xzz_yyyyy, tg_xzz_yyyyz, tg_xzz_yyyzz, tg_xzz_yyzzz, tg_xzz_yzzzz, tg_xzz_zzzzz, tg_y_xxxxx, tg_y_xxxxy, tg_y_xxxxz, tg_y_xxxyy, tg_y_xxxyz, tg_y_xxxzz, tg_y_xxyyy, tg_y_xxyyz, tg_y_xxyzz, tg_y_xxzzz, tg_y_xyyyy, tg_y_xyyyz, tg_y_xyyzz, tg_y_xyzzz, tg_y_xzzzz, tg_y_yyyyy, tg_y_yyyyz, tg_y_yyyzz, tg_y_yyzzz, tg_y_yzzzz, tg_y_zzzzz, tg_yy_xxxx, tg_yy_xxxxx, tg_yy_xxxxy, tg_yy_xxxxz, tg_yy_xxxy, tg_yy_xxxyy, tg_yy_xxxyz, tg_yy_xxxz, tg_yy_xxxzz, tg_yy_xxyy, tg_yy_xxyyy, tg_yy_xxyyz, tg_yy_xxyz, tg_yy_xxyzz, tg_yy_xxzz, tg_yy_xxzzz, tg_yy_xyyy, tg_yy_xyyyy, tg_yy_xyyyz, tg_yy_xyyz, tg_yy_xyyzz, tg_yy_xyzz, tg_yy_xyzzz, tg_yy_xzzz, tg_yy_xzzzz, tg_yy_yyyy, tg_yy_yyyyy, tg_yy_yyyyz, tg_yy_yyyz, tg_yy_yyyzz, tg_yy_yyzz, tg_yy_yyzzz, tg_yy_yzzz, tg_yy_yzzzz, tg_yy_zzzz, tg_yy_zzzzz, tg_yyy_xxxxx, tg_yyy_xxxxy, tg_yyy_xxxxz, tg_yyy_xxxyy, tg_yyy_xxxyz, tg_yyy_xxxzz, tg_yyy_xxyyy, tg_yyy_xxyyz, tg_yyy_xxyzz, tg_yyy_xxzzz, tg_yyy_xyyyy, tg_yyy_xyyyz, tg_yyy_xyyzz, tg_yyy_xyzzz, tg_yyy_xzzzz, tg_yyy_yyyyy, tg_yyy_yyyyz, tg_yyy_yyyzz, tg_yyy_yyzzz, tg_yyy_yzzzz, tg_yyy_zzzzz, tg_yyz_xxxxx, tg_yyz_xxxxy, tg_yyz_xxxxz, tg_yyz_xxxyy, tg_yyz_xxxyz, tg_yyz_xxxzz, tg_yyz_xxyyy, tg_yyz_xxyyz, tg_yyz_xxyzz, tg_yyz_xxzzz, tg_yyz_xyyyy, tg_yyz_xyyyz, tg_yyz_xyyzz, tg_yyz_xyzzz, tg_yyz_xzzzz, tg_yyz_yyyyy, tg_yyz_yyyyz, tg_yyz_yyyzz, tg_yyz_yyzzz, tg_yyz_yzzzz, tg_yyz_zzzzz, tg_yz_xxxxz, tg_yz_xxxyz, tg_yz_xxxzz, tg_yz_xxyyz, tg_yz_xxyz, tg_yz_xxyzz, tg_yz_xxzzz, tg_yz_xyyyz, tg_yz_xyyz, tg_yz_xyyzz, tg_yz_xyzz, tg_yz_xyzzz, tg_yz_xzzzz, tg_yz_yyyyy, tg_yz_yyyyz, tg_yz_yyyz, tg_yz_yyyzz, tg_yz_yyzz, tg_yz_yyzzz, tg_yz_yzzz, tg_yz_yzzzz, tg_yz_zzzzz, tg_yzz_xxxxx, tg_yzz_xxxxy, tg_yzz_xxxxz, tg_yzz_xxxyy, tg_yzz_xxxyz, tg_yzz_xxxzz, tg_yzz_xxyyy, tg_yzz_xxyyz, tg_yzz_xxyzz, tg_yzz_xxzzz, tg_yzz_xyyyy, tg_yzz_xyyyz, tg_yzz_xyyzz, tg_yzz_xyzzz, tg_yzz_xzzzz, tg_yzz_yyyyy, tg_yzz_yyyyz, tg_yzz_yyyzz, tg_yzz_yyzzz, tg_yzz_yzzzz, tg_yzz_zzzzz, tg_z_xxxxx, tg_z_xxxxy, tg_z_xxxxz, tg_z_xxxyy, tg_z_xxxyz, tg_z_xxxzz, tg_z_xxyyy, tg_z_xxyyz, tg_z_xxyzz, tg_z_xxzzz, tg_z_xyyyy, tg_z_xyyyz, tg_z_xyyzz, tg_z_xyzzz, tg_z_xzzzz, tg_z_yyyyy, tg_z_yyyyz, tg_z_yyyzz, tg_z_yyzzz, tg_z_yzzzz, tg_z_zzzzz, tg_zz_xxxx, tg_zz_xxxxx, tg_zz_xxxxy, tg_zz_xxxxz, tg_zz_xxxy, tg_zz_xxxyy, tg_zz_xxxyz, tg_zz_xxxz, tg_zz_xxxzz, tg_zz_xxyy, tg_zz_xxyyy, tg_zz_xxyyz, tg_zz_xxyz, tg_zz_xxyzz, tg_zz_xxzz, tg_zz_xxzzz, tg_zz_xyyy, tg_zz_xyyyy, tg_zz_xyyyz, tg_zz_xyyz, tg_zz_xyyzz, tg_zz_xyzz, tg_zz_xyzzz, tg_zz_xzzz, tg_zz_xzzzz, tg_zz_yyyy, tg_zz_yyyyy, tg_zz_yyyyz, tg_zz_yyyz, tg_zz_yyyzz, tg_zz_yyzz, tg_zz_yyzzz, tg_zz_yzzz, tg_zz_yzzzz, tg_zz_zzzz, tg_zz_zzzzz, tg_zzz_xxxxx, tg_zzz_xxxxy, tg_zzz_xxxxz, tg_zzz_xxxyy, tg_zzz_xxxyz, tg_zzz_xxxzz, tg_zzz_xxyyy, tg_zzz_xxyyz, tg_zzz_xxyzz, tg_zzz_xxzzz, tg_zzz_xyyyy, tg_zzz_xyyyz, tg_zzz_xyyzz, tg_zzz_xyzzz, tg_zzz_xzzzz, tg_zzz_yyyyy, tg_zzz_yyyyz, tg_zzz_yyyzz, tg_zzz_yyzzz, tg_zzz_yzzzz, tg_zzz_zzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxx_xxxxx[i] = 2.0 * tg_x_xxxxx[i] * fxi[i] + 5.0 * tg_xx_xxxx[i] * fxi[i] + tg_xx_xxxxx[i] * ra_x[i];

        tg_xxx_xxxxy[i] = 2.0 * tg_x_xxxxy[i] * fxi[i] + 4.0 * tg_xx_xxxy[i] * fxi[i] + tg_xx_xxxxy[i] * ra_x[i];

        tg_xxx_xxxxz[i] = 2.0 * tg_x_xxxxz[i] * fxi[i] + 4.0 * tg_xx_xxxz[i] * fxi[i] + tg_xx_xxxxz[i] * ra_x[i];

        tg_xxx_xxxyy[i] = 2.0 * tg_x_xxxyy[i] * fxi[i] + 3.0 * tg_xx_xxyy[i] * fxi[i] + tg_xx_xxxyy[i] * ra_x[i];

        tg_xxx_xxxyz[i] = 2.0 * tg_x_xxxyz[i] * fxi[i] + 3.0 * tg_xx_xxyz[i] * fxi[i] + tg_xx_xxxyz[i] * ra_x[i];

        tg_xxx_xxxzz[i] = 2.0 * tg_x_xxxzz[i] * fxi[i] + 3.0 * tg_xx_xxzz[i] * fxi[i] + tg_xx_xxxzz[i] * ra_x[i];

        tg_xxx_xxyyy[i] = 2.0 * tg_x_xxyyy[i] * fxi[i] + 2.0 * tg_xx_xyyy[i] * fxi[i] + tg_xx_xxyyy[i] * ra_x[i];

        tg_xxx_xxyyz[i] = 2.0 * tg_x_xxyyz[i] * fxi[i] + 2.0 * tg_xx_xyyz[i] * fxi[i] + tg_xx_xxyyz[i] * ra_x[i];

        tg_xxx_xxyzz[i] = 2.0 * tg_x_xxyzz[i] * fxi[i] + 2.0 * tg_xx_xyzz[i] * fxi[i] + tg_xx_xxyzz[i] * ra_x[i];

        tg_xxx_xxzzz[i] = 2.0 * tg_x_xxzzz[i] * fxi[i] + 2.0 * tg_xx_xzzz[i] * fxi[i] + tg_xx_xxzzz[i] * ra_x[i];

        tg_xxx_xyyyy[i] = 2.0 * tg_x_xyyyy[i] * fxi[i] + tg_xx_yyyy[i] * fxi[i] + tg_xx_xyyyy[i] * ra_x[i];

        tg_xxx_xyyyz[i] = 2.0 * tg_x_xyyyz[i] * fxi[i] + tg_xx_yyyz[i] * fxi[i] + tg_xx_xyyyz[i] * ra_x[i];

        tg_xxx_xyyzz[i] = 2.0 * tg_x_xyyzz[i] * fxi[i] + tg_xx_yyzz[i] * fxi[i] + tg_xx_xyyzz[i] * ra_x[i];

        tg_xxx_xyzzz[i] = 2.0 * tg_x_xyzzz[i] * fxi[i] + tg_xx_yzzz[i] * fxi[i] + tg_xx_xyzzz[i] * ra_x[i];

        tg_xxx_xzzzz[i] = 2.0 * tg_x_xzzzz[i] * fxi[i] + tg_xx_zzzz[i] * fxi[i] + tg_xx_xzzzz[i] * ra_x[i];

        tg_xxx_yyyyy[i] = 2.0 * tg_x_yyyyy[i] * fxi[i] + tg_xx_yyyyy[i] * ra_x[i];

        tg_xxx_yyyyz[i] = 2.0 * tg_x_yyyyz[i] * fxi[i] + tg_xx_yyyyz[i] * ra_x[i];

        tg_xxx_yyyzz[i] = 2.0 * tg_x_yyyzz[i] * fxi[i] + tg_xx_yyyzz[i] * ra_x[i];

        tg_xxx_yyzzz[i] = 2.0 * tg_x_yyzzz[i] * fxi[i] + tg_xx_yyzzz[i] * ra_x[i];

        tg_xxx_yzzzz[i] = 2.0 * tg_x_yzzzz[i] * fxi[i] + tg_xx_yzzzz[i] * ra_x[i];

        tg_xxx_zzzzz[i] = 2.0 * tg_x_zzzzz[i] * fxi[i] + tg_xx_zzzzz[i] * ra_x[i];

        tg_xxy_xxxxx[i] = tg_xx_xxxxx[i] * ra_y[i];

        tg_xxy_xxxxy[i] = tg_xx_xxxx[i] * fxi[i] + tg_xx_xxxxy[i] * ra_y[i];

        tg_xxy_xxxxz[i] = tg_xx_xxxxz[i] * ra_y[i];

        tg_xxy_xxxyy[i] = 2.0 * tg_xx_xxxy[i] * fxi[i] + tg_xx_xxxyy[i] * ra_y[i];

        tg_xxy_xxxyz[i] = tg_xx_xxxz[i] * fxi[i] + tg_xx_xxxyz[i] * ra_y[i];

        tg_xxy_xxxzz[i] = tg_xx_xxxzz[i] * ra_y[i];

        tg_xxy_xxyyy[i] = 3.0 * tg_xx_xxyy[i] * fxi[i] + tg_xx_xxyyy[i] * ra_y[i];

        tg_xxy_xxyyz[i] = 2.0 * tg_xx_xxyz[i] * fxi[i] + tg_xx_xxyyz[i] * ra_y[i];

        tg_xxy_xxyzz[i] = tg_xx_xxzz[i] * fxi[i] + tg_xx_xxyzz[i] * ra_y[i];

        tg_xxy_xxzzz[i] = tg_xx_xxzzz[i] * ra_y[i];

        tg_xxy_xyyyy[i] = 4.0 * tg_xx_xyyy[i] * fxi[i] + tg_xx_xyyyy[i] * ra_y[i];

        tg_xxy_xyyyz[i] = 3.0 * tg_xx_xyyz[i] * fxi[i] + tg_xx_xyyyz[i] * ra_y[i];

        tg_xxy_xyyzz[i] = 2.0 * tg_xx_xyzz[i] * fxi[i] + tg_xx_xyyzz[i] * ra_y[i];

        tg_xxy_xyzzz[i] = tg_xx_xzzz[i] * fxi[i] + tg_xx_xyzzz[i] * ra_y[i];

        tg_xxy_xzzzz[i] = tg_xx_xzzzz[i] * ra_y[i];

        tg_xxy_yyyyy[i] = tg_y_yyyyy[i] * fxi[i] + tg_xy_yyyyy[i] * ra_x[i];

        tg_xxy_yyyyz[i] = tg_y_yyyyz[i] * fxi[i] + tg_xy_yyyyz[i] * ra_x[i];

        tg_xxy_yyyzz[i] = tg_y_yyyzz[i] * fxi[i] + tg_xy_yyyzz[i] * ra_x[i];

        tg_xxy_yyzzz[i] = tg_y_yyzzz[i] * fxi[i] + tg_xy_yyzzz[i] * ra_x[i];

        tg_xxy_yzzzz[i] = tg_y_yzzzz[i] * fxi[i] + tg_xy_yzzzz[i] * ra_x[i];

        tg_xxy_zzzzz[i] = tg_xx_zzzzz[i] * ra_y[i];

        tg_xxz_xxxxx[i] = tg_xx_xxxxx[i] * ra_z[i];

        tg_xxz_xxxxy[i] = tg_xx_xxxxy[i] * ra_z[i];

        tg_xxz_xxxxz[i] = tg_xx_xxxx[i] * fxi[i] + tg_xx_xxxxz[i] * ra_z[i];

        tg_xxz_xxxyy[i] = tg_xx_xxxyy[i] * ra_z[i];

        tg_xxz_xxxyz[i] = tg_xx_xxxy[i] * fxi[i] + tg_xx_xxxyz[i] * ra_z[i];

        tg_xxz_xxxzz[i] = 2.0 * tg_xx_xxxz[i] * fxi[i] + tg_xx_xxxzz[i] * ra_z[i];

        tg_xxz_xxyyy[i] = tg_xx_xxyyy[i] * ra_z[i];

        tg_xxz_xxyyz[i] = tg_xx_xxyy[i] * fxi[i] + tg_xx_xxyyz[i] * ra_z[i];

        tg_xxz_xxyzz[i] = 2.0 * tg_xx_xxyz[i] * fxi[i] + tg_xx_xxyzz[i] * ra_z[i];

        tg_xxz_xxzzz[i] = 3.0 * tg_xx_xxzz[i] * fxi[i] + tg_xx_xxzzz[i] * ra_z[i];

        tg_xxz_xyyyy[i] = tg_xx_xyyyy[i] * ra_z[i];

        tg_xxz_xyyyz[i] = tg_xx_xyyy[i] * fxi[i] + tg_xx_xyyyz[i] * ra_z[i];

        tg_xxz_xyyzz[i] = 2.0 * tg_xx_xyyz[i] * fxi[i] + tg_xx_xyyzz[i] * ra_z[i];

        tg_xxz_xyzzz[i] = 3.0 * tg_xx_xyzz[i] * fxi[i] + tg_xx_xyzzz[i] * ra_z[i];

        tg_xxz_xzzzz[i] = 4.0 * tg_xx_xzzz[i] * fxi[i] + tg_xx_xzzzz[i] * ra_z[i];

        tg_xxz_yyyyy[i] = tg_xx_yyyyy[i] * ra_z[i];

        tg_xxz_yyyyz[i] = tg_z_yyyyz[i] * fxi[i] + tg_xz_yyyyz[i] * ra_x[i];

        tg_xxz_yyyzz[i] = tg_z_yyyzz[i] * fxi[i] + tg_xz_yyyzz[i] * ra_x[i];

        tg_xxz_yyzzz[i] = tg_z_yyzzz[i] * fxi[i] + tg_xz_yyzzz[i] * ra_x[i];

        tg_xxz_yzzzz[i] = tg_z_yzzzz[i] * fxi[i] + tg_xz_yzzzz[i] * ra_x[i];

        tg_xxz_zzzzz[i] = tg_z_zzzzz[i] * fxi[i] + tg_xz_zzzzz[i] * ra_x[i];

        tg_xyy_xxxxx[i] = 5.0 * tg_yy_xxxx[i] * fxi[i] + tg_yy_xxxxx[i] * ra_x[i];

        tg_xyy_xxxxy[i] = 4.0 * tg_yy_xxxy[i] * fxi[i] + tg_yy_xxxxy[i] * ra_x[i];

        tg_xyy_xxxxz[i] = 4.0 * tg_yy_xxxz[i] * fxi[i] + tg_yy_xxxxz[i] * ra_x[i];

        tg_xyy_xxxyy[i] = 3.0 * tg_yy_xxyy[i] * fxi[i] + tg_yy_xxxyy[i] * ra_x[i];

        tg_xyy_xxxyz[i] = 3.0 * tg_yy_xxyz[i] * fxi[i] + tg_yy_xxxyz[i] * ra_x[i];

        tg_xyy_xxxzz[i] = 3.0 * tg_yy_xxzz[i] * fxi[i] + tg_yy_xxxzz[i] * ra_x[i];

        tg_xyy_xxyyy[i] = 2.0 * tg_yy_xyyy[i] * fxi[i] + tg_yy_xxyyy[i] * ra_x[i];

        tg_xyy_xxyyz[i] = 2.0 * tg_yy_xyyz[i] * fxi[i] + tg_yy_xxyyz[i] * ra_x[i];

        tg_xyy_xxyzz[i] = 2.0 * tg_yy_xyzz[i] * fxi[i] + tg_yy_xxyzz[i] * ra_x[i];

        tg_xyy_xxzzz[i] = 2.0 * tg_yy_xzzz[i] * fxi[i] + tg_yy_xxzzz[i] * ra_x[i];

        tg_xyy_xyyyy[i] = tg_yy_yyyy[i] * fxi[i] + tg_yy_xyyyy[i] * ra_x[i];

        tg_xyy_xyyyz[i] = tg_yy_yyyz[i] * fxi[i] + tg_yy_xyyyz[i] * ra_x[i];

        tg_xyy_xyyzz[i] = tg_yy_yyzz[i] * fxi[i] + tg_yy_xyyzz[i] * ra_x[i];

        tg_xyy_xyzzz[i] = tg_yy_yzzz[i] * fxi[i] + tg_yy_xyzzz[i] * ra_x[i];

        tg_xyy_xzzzz[i] = tg_yy_zzzz[i] * fxi[i] + tg_yy_xzzzz[i] * ra_x[i];

        tg_xyy_yyyyy[i] = tg_yy_yyyyy[i] * ra_x[i];

        tg_xyy_yyyyz[i] = tg_yy_yyyyz[i] * ra_x[i];

        tg_xyy_yyyzz[i] = tg_yy_yyyzz[i] * ra_x[i];

        tg_xyy_yyzzz[i] = tg_yy_yyzzz[i] * ra_x[i];

        tg_xyy_yzzzz[i] = tg_yy_yzzzz[i] * ra_x[i];

        tg_xyy_zzzzz[i] = tg_yy_zzzzz[i] * ra_x[i];

        tg_xyz_xxxxx[i] = tg_xz_xxxxx[i] * ra_y[i];

        tg_xyz_xxxxy[i] = tg_xy_xxxxy[i] * ra_z[i];

        tg_xyz_xxxxz[i] = tg_xz_xxxxz[i] * ra_y[i];

        tg_xyz_xxxyy[i] = tg_xy_xxxyy[i] * ra_z[i];

        tg_xyz_xxxyz[i] = 3.0 * tg_yz_xxyz[i] * fxi[i] + tg_yz_xxxyz[i] * ra_x[i];

        tg_xyz_xxxzz[i] = tg_xz_xxxzz[i] * ra_y[i];

        tg_xyz_xxyyy[i] = tg_xy_xxyyy[i] * ra_z[i];

        tg_xyz_xxyyz[i] = 2.0 * tg_yz_xyyz[i] * fxi[i] + tg_yz_xxyyz[i] * ra_x[i];

        tg_xyz_xxyzz[i] = 2.0 * tg_yz_xyzz[i] * fxi[i] + tg_yz_xxyzz[i] * ra_x[i];

        tg_xyz_xxzzz[i] = tg_xz_xxzzz[i] * ra_y[i];

        tg_xyz_xyyyy[i] = tg_xy_xyyyy[i] * ra_z[i];

        tg_xyz_xyyyz[i] = tg_yz_yyyz[i] * fxi[i] + tg_yz_xyyyz[i] * ra_x[i];

        tg_xyz_xyyzz[i] = tg_yz_yyzz[i] * fxi[i] + tg_yz_xyyzz[i] * ra_x[i];

        tg_xyz_xyzzz[i] = tg_yz_yzzz[i] * fxi[i] + tg_yz_xyzzz[i] * ra_x[i];

        tg_xyz_xzzzz[i] = tg_xz_xzzzz[i] * ra_y[i];

        tg_xyz_yyyyy[i] = tg_yz_yyyyy[i] * ra_x[i];

        tg_xyz_yyyyz[i] = tg_yz_yyyyz[i] * ra_x[i];

        tg_xyz_yyyzz[i] = tg_yz_yyyzz[i] * ra_x[i];

        tg_xyz_yyzzz[i] = tg_yz_yyzzz[i] * ra_x[i];

        tg_xyz_yzzzz[i] = tg_yz_yzzzz[i] * ra_x[i];

        tg_xyz_zzzzz[i] = tg_yz_zzzzz[i] * ra_x[i];

        tg_xzz_xxxxx[i] = 5.0 * tg_zz_xxxx[i] * fxi[i] + tg_zz_xxxxx[i] * ra_x[i];

        tg_xzz_xxxxy[i] = 4.0 * tg_zz_xxxy[i] * fxi[i] + tg_zz_xxxxy[i] * ra_x[i];

        tg_xzz_xxxxz[i] = 4.0 * tg_zz_xxxz[i] * fxi[i] + tg_zz_xxxxz[i] * ra_x[i];

        tg_xzz_xxxyy[i] = 3.0 * tg_zz_xxyy[i] * fxi[i] + tg_zz_xxxyy[i] * ra_x[i];

        tg_xzz_xxxyz[i] = 3.0 * tg_zz_xxyz[i] * fxi[i] + tg_zz_xxxyz[i] * ra_x[i];

        tg_xzz_xxxzz[i] = 3.0 * tg_zz_xxzz[i] * fxi[i] + tg_zz_xxxzz[i] * ra_x[i];

        tg_xzz_xxyyy[i] = 2.0 * tg_zz_xyyy[i] * fxi[i] + tg_zz_xxyyy[i] * ra_x[i];

        tg_xzz_xxyyz[i] = 2.0 * tg_zz_xyyz[i] * fxi[i] + tg_zz_xxyyz[i] * ra_x[i];

        tg_xzz_xxyzz[i] = 2.0 * tg_zz_xyzz[i] * fxi[i] + tg_zz_xxyzz[i] * ra_x[i];

        tg_xzz_xxzzz[i] = 2.0 * tg_zz_xzzz[i] * fxi[i] + tg_zz_xxzzz[i] * ra_x[i];

        tg_xzz_xyyyy[i] = tg_zz_yyyy[i] * fxi[i] + tg_zz_xyyyy[i] * ra_x[i];

        tg_xzz_xyyyz[i] = tg_zz_yyyz[i] * fxi[i] + tg_zz_xyyyz[i] * ra_x[i];

        tg_xzz_xyyzz[i] = tg_zz_yyzz[i] * fxi[i] + tg_zz_xyyzz[i] * ra_x[i];

        tg_xzz_xyzzz[i] = tg_zz_yzzz[i] * fxi[i] + tg_zz_xyzzz[i] * ra_x[i];

        tg_xzz_xzzzz[i] = tg_zz_zzzz[i] * fxi[i] + tg_zz_xzzzz[i] * ra_x[i];

        tg_xzz_yyyyy[i] = tg_zz_yyyyy[i] * ra_x[i];

        tg_xzz_yyyyz[i] = tg_zz_yyyyz[i] * ra_x[i];

        tg_xzz_yyyzz[i] = tg_zz_yyyzz[i] * ra_x[i];

        tg_xzz_yyzzz[i] = tg_zz_yyzzz[i] * ra_x[i];

        tg_xzz_yzzzz[i] = tg_zz_yzzzz[i] * ra_x[i];

        tg_xzz_zzzzz[i] = tg_zz_zzzzz[i] * ra_x[i];

        tg_yyy_xxxxx[i] = 2.0 * tg_y_xxxxx[i] * fxi[i] + tg_yy_xxxxx[i] * ra_y[i];

        tg_yyy_xxxxy[i] = 2.0 * tg_y_xxxxy[i] * fxi[i] + tg_yy_xxxx[i] * fxi[i] + tg_yy_xxxxy[i] * ra_y[i];

        tg_yyy_xxxxz[i] = 2.0 * tg_y_xxxxz[i] * fxi[i] + tg_yy_xxxxz[i] * ra_y[i];

        tg_yyy_xxxyy[i] = 2.0 * tg_y_xxxyy[i] * fxi[i] + 2.0 * tg_yy_xxxy[i] * fxi[i] + tg_yy_xxxyy[i] * ra_y[i];

        tg_yyy_xxxyz[i] = 2.0 * tg_y_xxxyz[i] * fxi[i] + tg_yy_xxxz[i] * fxi[i] + tg_yy_xxxyz[i] * ra_y[i];

        tg_yyy_xxxzz[i] = 2.0 * tg_y_xxxzz[i] * fxi[i] + tg_yy_xxxzz[i] * ra_y[i];

        tg_yyy_xxyyy[i] = 2.0 * tg_y_xxyyy[i] * fxi[i] + 3.0 * tg_yy_xxyy[i] * fxi[i] + tg_yy_xxyyy[i] * ra_y[i];

        tg_yyy_xxyyz[i] = 2.0 * tg_y_xxyyz[i] * fxi[i] + 2.0 * tg_yy_xxyz[i] * fxi[i] + tg_yy_xxyyz[i] * ra_y[i];

        tg_yyy_xxyzz[i] = 2.0 * tg_y_xxyzz[i] * fxi[i] + tg_yy_xxzz[i] * fxi[i] + tg_yy_xxyzz[i] * ra_y[i];

        tg_yyy_xxzzz[i] = 2.0 * tg_y_xxzzz[i] * fxi[i] + tg_yy_xxzzz[i] * ra_y[i];

        tg_yyy_xyyyy[i] = 2.0 * tg_y_xyyyy[i] * fxi[i] + 4.0 * tg_yy_xyyy[i] * fxi[i] + tg_yy_xyyyy[i] * ra_y[i];

        tg_yyy_xyyyz[i] = 2.0 * tg_y_xyyyz[i] * fxi[i] + 3.0 * tg_yy_xyyz[i] * fxi[i] + tg_yy_xyyyz[i] * ra_y[i];

        tg_yyy_xyyzz[i] = 2.0 * tg_y_xyyzz[i] * fxi[i] + 2.0 * tg_yy_xyzz[i] * fxi[i] + tg_yy_xyyzz[i] * ra_y[i];

        tg_yyy_xyzzz[i] = 2.0 * tg_y_xyzzz[i] * fxi[i] + tg_yy_xzzz[i] * fxi[i] + tg_yy_xyzzz[i] * ra_y[i];

        tg_yyy_xzzzz[i] = 2.0 * tg_y_xzzzz[i] * fxi[i] + tg_yy_xzzzz[i] * ra_y[i];

        tg_yyy_yyyyy[i] = 2.0 * tg_y_yyyyy[i] * fxi[i] + 5.0 * tg_yy_yyyy[i] * fxi[i] + tg_yy_yyyyy[i] * ra_y[i];

        tg_yyy_yyyyz[i] = 2.0 * tg_y_yyyyz[i] * fxi[i] + 4.0 * tg_yy_yyyz[i] * fxi[i] + tg_yy_yyyyz[i] * ra_y[i];

        tg_yyy_yyyzz[i] = 2.0 * tg_y_yyyzz[i] * fxi[i] + 3.0 * tg_yy_yyzz[i] * fxi[i] + tg_yy_yyyzz[i] * ra_y[i];

        tg_yyy_yyzzz[i] = 2.0 * tg_y_yyzzz[i] * fxi[i] + 2.0 * tg_yy_yzzz[i] * fxi[i] + tg_yy_yyzzz[i] * ra_y[i];

        tg_yyy_yzzzz[i] = 2.0 * tg_y_yzzzz[i] * fxi[i] + tg_yy_zzzz[i] * fxi[i] + tg_yy_yzzzz[i] * ra_y[i];

        tg_yyy_zzzzz[i] = 2.0 * tg_y_zzzzz[i] * fxi[i] + tg_yy_zzzzz[i] * ra_y[i];

        tg_yyz_xxxxx[i] = tg_yy_xxxxx[i] * ra_z[i];

        tg_yyz_xxxxy[i] = tg_yy_xxxxy[i] * ra_z[i];

        tg_yyz_xxxxz[i] = tg_z_xxxxz[i] * fxi[i] + tg_yz_xxxxz[i] * ra_y[i];

        tg_yyz_xxxyy[i] = tg_yy_xxxyy[i] * ra_z[i];

        tg_yyz_xxxyz[i] = tg_yy_xxxy[i] * fxi[i] + tg_yy_xxxyz[i] * ra_z[i];

        tg_yyz_xxxzz[i] = tg_z_xxxzz[i] * fxi[i] + tg_yz_xxxzz[i] * ra_y[i];

        tg_yyz_xxyyy[i] = tg_yy_xxyyy[i] * ra_z[i];

        tg_yyz_xxyyz[i] = tg_yy_xxyy[i] * fxi[i] + tg_yy_xxyyz[i] * ra_z[i];

        tg_yyz_xxyzz[i] = 2.0 * tg_yy_xxyz[i] * fxi[i] + tg_yy_xxyzz[i] * ra_z[i];

        tg_yyz_xxzzz[i] = tg_z_xxzzz[i] * fxi[i] + tg_yz_xxzzz[i] * ra_y[i];

        tg_yyz_xyyyy[i] = tg_yy_xyyyy[i] * ra_z[i];

        tg_yyz_xyyyz[i] = tg_yy_xyyy[i] * fxi[i] + tg_yy_xyyyz[i] * ra_z[i];

        tg_yyz_xyyzz[i] = 2.0 * tg_yy_xyyz[i] * fxi[i] + tg_yy_xyyzz[i] * ra_z[i];

        tg_yyz_xyzzz[i] = 3.0 * tg_yy_xyzz[i] * fxi[i] + tg_yy_xyzzz[i] * ra_z[i];

        tg_yyz_xzzzz[i] = tg_z_xzzzz[i] * fxi[i] + tg_yz_xzzzz[i] * ra_y[i];

        tg_yyz_yyyyy[i] = tg_yy_yyyyy[i] * ra_z[i];

        tg_yyz_yyyyz[i] = tg_yy_yyyy[i] * fxi[i] + tg_yy_yyyyz[i] * ra_z[i];

        tg_yyz_yyyzz[i] = 2.0 * tg_yy_yyyz[i] * fxi[i] + tg_yy_yyyzz[i] * ra_z[i];

        tg_yyz_yyzzz[i] = 3.0 * tg_yy_yyzz[i] * fxi[i] + tg_yy_yyzzz[i] * ra_z[i];

        tg_yyz_yzzzz[i] = 4.0 * tg_yy_yzzz[i] * fxi[i] + tg_yy_yzzzz[i] * ra_z[i];

        tg_yyz_zzzzz[i] = tg_z_zzzzz[i] * fxi[i] + tg_yz_zzzzz[i] * ra_y[i];

        tg_yzz_xxxxx[i] = tg_zz_xxxxx[i] * ra_y[i];

        tg_yzz_xxxxy[i] = tg_zz_xxxx[i] * fxi[i] + tg_zz_xxxxy[i] * ra_y[i];

        tg_yzz_xxxxz[i] = tg_zz_xxxxz[i] * ra_y[i];

        tg_yzz_xxxyy[i] = 2.0 * tg_zz_xxxy[i] * fxi[i] + tg_zz_xxxyy[i] * ra_y[i];

        tg_yzz_xxxyz[i] = tg_zz_xxxz[i] * fxi[i] + tg_zz_xxxyz[i] * ra_y[i];

        tg_yzz_xxxzz[i] = tg_zz_xxxzz[i] * ra_y[i];

        tg_yzz_xxyyy[i] = 3.0 * tg_zz_xxyy[i] * fxi[i] + tg_zz_xxyyy[i] * ra_y[i];

        tg_yzz_xxyyz[i] = 2.0 * tg_zz_xxyz[i] * fxi[i] + tg_zz_xxyyz[i] * ra_y[i];

        tg_yzz_xxyzz[i] = tg_zz_xxzz[i] * fxi[i] + tg_zz_xxyzz[i] * ra_y[i];

        tg_yzz_xxzzz[i] = tg_zz_xxzzz[i] * ra_y[i];

        tg_yzz_xyyyy[i] = 4.0 * tg_zz_xyyy[i] * fxi[i] + tg_zz_xyyyy[i] * ra_y[i];

        tg_yzz_xyyyz[i] = 3.0 * tg_zz_xyyz[i] * fxi[i] + tg_zz_xyyyz[i] * ra_y[i];

        tg_yzz_xyyzz[i] = 2.0 * tg_zz_xyzz[i] * fxi[i] + tg_zz_xyyzz[i] * ra_y[i];

        tg_yzz_xyzzz[i] = tg_zz_xzzz[i] * fxi[i] + tg_zz_xyzzz[i] * ra_y[i];

        tg_yzz_xzzzz[i] = tg_zz_xzzzz[i] * ra_y[i];

        tg_yzz_yyyyy[i] = 5.0 * tg_zz_yyyy[i] * fxi[i] + tg_zz_yyyyy[i] * ra_y[i];

        tg_yzz_yyyyz[i] = 4.0 * tg_zz_yyyz[i] * fxi[i] + tg_zz_yyyyz[i] * ra_y[i];

        tg_yzz_yyyzz[i] = 3.0 * tg_zz_yyzz[i] * fxi[i] + tg_zz_yyyzz[i] * ra_y[i];

        tg_yzz_yyzzz[i] = 2.0 * tg_zz_yzzz[i] * fxi[i] + tg_zz_yyzzz[i] * ra_y[i];

        tg_yzz_yzzzz[i] = tg_zz_zzzz[i] * fxi[i] + tg_zz_yzzzz[i] * ra_y[i];

        tg_yzz_zzzzz[i] = tg_zz_zzzzz[i] * ra_y[i];

        tg_zzz_xxxxx[i] = 2.0 * tg_z_xxxxx[i] * fxi[i] + tg_zz_xxxxx[i] * ra_z[i];

        tg_zzz_xxxxy[i] = 2.0 * tg_z_xxxxy[i] * fxi[i] + tg_zz_xxxxy[i] * ra_z[i];

        tg_zzz_xxxxz[i] = 2.0 * tg_z_xxxxz[i] * fxi[i] + tg_zz_xxxx[i] * fxi[i] + tg_zz_xxxxz[i] * ra_z[i];

        tg_zzz_xxxyy[i] = 2.0 * tg_z_xxxyy[i] * fxi[i] + tg_zz_xxxyy[i] * ra_z[i];

        tg_zzz_xxxyz[i] = 2.0 * tg_z_xxxyz[i] * fxi[i] + tg_zz_xxxy[i] * fxi[i] + tg_zz_xxxyz[i] * ra_z[i];

        tg_zzz_xxxzz[i] = 2.0 * tg_z_xxxzz[i] * fxi[i] + 2.0 * tg_zz_xxxz[i] * fxi[i] + tg_zz_xxxzz[i] * ra_z[i];

        tg_zzz_xxyyy[i] = 2.0 * tg_z_xxyyy[i] * fxi[i] + tg_zz_xxyyy[i] * ra_z[i];

        tg_zzz_xxyyz[i] = 2.0 * tg_z_xxyyz[i] * fxi[i] + tg_zz_xxyy[i] * fxi[i] + tg_zz_xxyyz[i] * ra_z[i];

        tg_zzz_xxyzz[i] = 2.0 * tg_z_xxyzz[i] * fxi[i] + 2.0 * tg_zz_xxyz[i] * fxi[i] + tg_zz_xxyzz[i] * ra_z[i];

        tg_zzz_xxzzz[i] = 2.0 * tg_z_xxzzz[i] * fxi[i] + 3.0 * tg_zz_xxzz[i] * fxi[i] + tg_zz_xxzzz[i] * ra_z[i];

        tg_zzz_xyyyy[i] = 2.0 * tg_z_xyyyy[i] * fxi[i] + tg_zz_xyyyy[i] * ra_z[i];

        tg_zzz_xyyyz[i] = 2.0 * tg_z_xyyyz[i] * fxi[i] + tg_zz_xyyy[i] * fxi[i] + tg_zz_xyyyz[i] * ra_z[i];

        tg_zzz_xyyzz[i] = 2.0 * tg_z_xyyzz[i] * fxi[i] + 2.0 * tg_zz_xyyz[i] * fxi[i] + tg_zz_xyyzz[i] * ra_z[i];

        tg_zzz_xyzzz[i] = 2.0 * tg_z_xyzzz[i] * fxi[i] + 3.0 * tg_zz_xyzz[i] * fxi[i] + tg_zz_xyzzz[i] * ra_z[i];

        tg_zzz_xzzzz[i] = 2.0 * tg_z_xzzzz[i] * fxi[i] + 4.0 * tg_zz_xzzz[i] * fxi[i] + tg_zz_xzzzz[i] * ra_z[i];

        tg_zzz_yyyyy[i] = 2.0 * tg_z_yyyyy[i] * fxi[i] + tg_zz_yyyyy[i] * ra_z[i];

        tg_zzz_yyyyz[i] = 2.0 * tg_z_yyyyz[i] * fxi[i] + tg_zz_yyyy[i] * fxi[i] + tg_zz_yyyyz[i] * ra_z[i];

        tg_zzz_yyyzz[i] = 2.0 * tg_z_yyyzz[i] * fxi[i] + 2.0 * tg_zz_yyyz[i] * fxi[i] + tg_zz_yyyzz[i] * ra_z[i];

        tg_zzz_yyzzz[i] = 2.0 * tg_z_yyzzz[i] * fxi[i] + 3.0 * tg_zz_yyzz[i] * fxi[i] + tg_zz_yyzzz[i] * ra_z[i];

        tg_zzz_yzzzz[i] = 2.0 * tg_z_yzzzz[i] * fxi[i] + 4.0 * tg_zz_yzzz[i] * fxi[i] + tg_zz_yzzzz[i] * ra_z[i];

        tg_zzz_zzzzz[i] = 2.0 * tg_z_zzzzz[i] * fxi[i] + 5.0 * tg_zz_zzzz[i] * fxi[i] + tg_zz_zzzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

