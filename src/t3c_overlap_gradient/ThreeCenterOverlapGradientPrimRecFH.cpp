#include "ThreeCenterOverlapGradientPrimRecFH.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_fh(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_fh,
                              const size_t idx_dh,
                              const size_t idx_fg,
                              const size_t idx_fh,
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

    auto ts_xxy_xxyz = pbuffer.data(idx_fg + 19);

    auto ts_xxy_xxzz = pbuffer.data(idx_fg + 20);

    auto ts_xxy_xyyy = pbuffer.data(idx_fg + 21);

    auto ts_xxy_xyyz = pbuffer.data(idx_fg + 22);

    auto ts_xxy_xyzz = pbuffer.data(idx_fg + 23);

    auto ts_xxy_xzzz = pbuffer.data(idx_fg + 24);

    auto ts_xxy_yyyy = pbuffer.data(idx_fg + 25);

    auto ts_xxy_yyyz = pbuffer.data(idx_fg + 26);

    auto ts_xxy_yyzz = pbuffer.data(idx_fg + 27);

    auto ts_xxy_yzzz = pbuffer.data(idx_fg + 28);

    auto ts_xxy_zzzz = pbuffer.data(idx_fg + 29);

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

    auto ts_xxz_yyyy = pbuffer.data(idx_fg + 40);

    auto ts_xxz_yyyz = pbuffer.data(idx_fg + 41);

    auto ts_xxz_yyzz = pbuffer.data(idx_fg + 42);

    auto ts_xxz_yzzz = pbuffer.data(idx_fg + 43);

    auto ts_xxz_zzzz = pbuffer.data(idx_fg + 44);

    auto ts_xyy_xxxx = pbuffer.data(idx_fg + 45);

    auto ts_xyy_xxxy = pbuffer.data(idx_fg + 46);

    auto ts_xyy_xxxz = pbuffer.data(idx_fg + 47);

    auto ts_xyy_xxyy = pbuffer.data(idx_fg + 48);

    auto ts_xyy_xxyz = pbuffer.data(idx_fg + 49);

    auto ts_xyy_xxzz = pbuffer.data(idx_fg + 50);

    auto ts_xyy_xyyy = pbuffer.data(idx_fg + 51);

    auto ts_xyy_xyyz = pbuffer.data(idx_fg + 52);

    auto ts_xyy_xyzz = pbuffer.data(idx_fg + 53);

    auto ts_xyy_xzzz = pbuffer.data(idx_fg + 54);

    auto ts_xyy_yyyy = pbuffer.data(idx_fg + 55);

    auto ts_xyy_yyyz = pbuffer.data(idx_fg + 56);

    auto ts_xyy_yyzz = pbuffer.data(idx_fg + 57);

    auto ts_xyy_yzzz = pbuffer.data(idx_fg + 58);

    auto ts_xyy_zzzz = pbuffer.data(idx_fg + 59);

    auto ts_xyz_xxxx = pbuffer.data(idx_fg + 60);

    auto ts_xyz_xxxy = pbuffer.data(idx_fg + 61);

    auto ts_xyz_xxxz = pbuffer.data(idx_fg + 62);

    auto ts_xyz_xxyy = pbuffer.data(idx_fg + 63);

    auto ts_xyz_xxyz = pbuffer.data(idx_fg + 64);

    auto ts_xyz_xxzz = pbuffer.data(idx_fg + 65);

    auto ts_xyz_xyyy = pbuffer.data(idx_fg + 66);

    auto ts_xyz_xyyz = pbuffer.data(idx_fg + 67);

    auto ts_xyz_xyzz = pbuffer.data(idx_fg + 68);

    auto ts_xyz_xzzz = pbuffer.data(idx_fg + 69);

    auto ts_xyz_yyyy = pbuffer.data(idx_fg + 70);

    auto ts_xyz_yyyz = pbuffer.data(idx_fg + 71);

    auto ts_xyz_yyzz = pbuffer.data(idx_fg + 72);

    auto ts_xyz_yzzz = pbuffer.data(idx_fg + 73);

    auto ts_xyz_zzzz = pbuffer.data(idx_fg + 74);

    auto ts_xzz_xxxx = pbuffer.data(idx_fg + 75);

    auto ts_xzz_xxxy = pbuffer.data(idx_fg + 76);

    auto ts_xzz_xxxz = pbuffer.data(idx_fg + 77);

    auto ts_xzz_xxyy = pbuffer.data(idx_fg + 78);

    auto ts_xzz_xxyz = pbuffer.data(idx_fg + 79);

    auto ts_xzz_xxzz = pbuffer.data(idx_fg + 80);

    auto ts_xzz_xyyy = pbuffer.data(idx_fg + 81);

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

    auto ts_yyz_xxxx = pbuffer.data(idx_fg + 105);

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

    // Set up components of auxiliary buffer : FH

    auto ts_xxx_xxxxx = pbuffer.data(idx_fh);

    auto ts_xxx_xxxxy = pbuffer.data(idx_fh + 1);

    auto ts_xxx_xxxxz = pbuffer.data(idx_fh + 2);

    auto ts_xxx_xxxyy = pbuffer.data(idx_fh + 3);

    auto ts_xxx_xxxyz = pbuffer.data(idx_fh + 4);

    auto ts_xxx_xxxzz = pbuffer.data(idx_fh + 5);

    auto ts_xxx_xxyyy = pbuffer.data(idx_fh + 6);

    auto ts_xxx_xxyyz = pbuffer.data(idx_fh + 7);

    auto ts_xxx_xxyzz = pbuffer.data(idx_fh + 8);

    auto ts_xxx_xxzzz = pbuffer.data(idx_fh + 9);

    auto ts_xxx_xyyyy = pbuffer.data(idx_fh + 10);

    auto ts_xxx_xyyyz = pbuffer.data(idx_fh + 11);

    auto ts_xxx_xyyzz = pbuffer.data(idx_fh + 12);

    auto ts_xxx_xyzzz = pbuffer.data(idx_fh + 13);

    auto ts_xxx_xzzzz = pbuffer.data(idx_fh + 14);

    auto ts_xxx_yyyyy = pbuffer.data(idx_fh + 15);

    auto ts_xxx_yyyyz = pbuffer.data(idx_fh + 16);

    auto ts_xxx_yyyzz = pbuffer.data(idx_fh + 17);

    auto ts_xxx_yyzzz = pbuffer.data(idx_fh + 18);

    auto ts_xxx_yzzzz = pbuffer.data(idx_fh + 19);

    auto ts_xxx_zzzzz = pbuffer.data(idx_fh + 20);

    auto ts_xxy_xxxxx = pbuffer.data(idx_fh + 21);

    auto ts_xxy_xxxxy = pbuffer.data(idx_fh + 22);

    auto ts_xxy_xxxxz = pbuffer.data(idx_fh + 23);

    auto ts_xxy_xxxyy = pbuffer.data(idx_fh + 24);

    auto ts_xxy_xxxyz = pbuffer.data(idx_fh + 25);

    auto ts_xxy_xxxzz = pbuffer.data(idx_fh + 26);

    auto ts_xxy_xxyyy = pbuffer.data(idx_fh + 27);

    auto ts_xxy_xxyyz = pbuffer.data(idx_fh + 28);

    auto ts_xxy_xxyzz = pbuffer.data(idx_fh + 29);

    auto ts_xxy_xxzzz = pbuffer.data(idx_fh + 30);

    auto ts_xxy_xyyyy = pbuffer.data(idx_fh + 31);

    auto ts_xxy_xyyyz = pbuffer.data(idx_fh + 32);

    auto ts_xxy_xyyzz = pbuffer.data(idx_fh + 33);

    auto ts_xxy_xyzzz = pbuffer.data(idx_fh + 34);

    auto ts_xxy_xzzzz = pbuffer.data(idx_fh + 35);

    auto ts_xxy_yyyyy = pbuffer.data(idx_fh + 36);

    auto ts_xxy_yyyyz = pbuffer.data(idx_fh + 37);

    auto ts_xxy_yyyzz = pbuffer.data(idx_fh + 38);

    auto ts_xxy_yyzzz = pbuffer.data(idx_fh + 39);

    auto ts_xxy_yzzzz = pbuffer.data(idx_fh + 40);

    auto ts_xxy_zzzzz = pbuffer.data(idx_fh + 41);

    auto ts_xxz_xxxxx = pbuffer.data(idx_fh + 42);

    auto ts_xxz_xxxxy = pbuffer.data(idx_fh + 43);

    auto ts_xxz_xxxxz = pbuffer.data(idx_fh + 44);

    auto ts_xxz_xxxyy = pbuffer.data(idx_fh + 45);

    auto ts_xxz_xxxyz = pbuffer.data(idx_fh + 46);

    auto ts_xxz_xxxzz = pbuffer.data(idx_fh + 47);

    auto ts_xxz_xxyyy = pbuffer.data(idx_fh + 48);

    auto ts_xxz_xxyyz = pbuffer.data(idx_fh + 49);

    auto ts_xxz_xxyzz = pbuffer.data(idx_fh + 50);

    auto ts_xxz_xxzzz = pbuffer.data(idx_fh + 51);

    auto ts_xxz_xyyyy = pbuffer.data(idx_fh + 52);

    auto ts_xxz_xyyyz = pbuffer.data(idx_fh + 53);

    auto ts_xxz_xyyzz = pbuffer.data(idx_fh + 54);

    auto ts_xxz_xyzzz = pbuffer.data(idx_fh + 55);

    auto ts_xxz_xzzzz = pbuffer.data(idx_fh + 56);

    auto ts_xxz_yyyyy = pbuffer.data(idx_fh + 57);

    auto ts_xxz_yyyyz = pbuffer.data(idx_fh + 58);

    auto ts_xxz_yyyzz = pbuffer.data(idx_fh + 59);

    auto ts_xxz_yyzzz = pbuffer.data(idx_fh + 60);

    auto ts_xxz_yzzzz = pbuffer.data(idx_fh + 61);

    auto ts_xxz_zzzzz = pbuffer.data(idx_fh + 62);

    auto ts_xyy_xxxxx = pbuffer.data(idx_fh + 63);

    auto ts_xyy_xxxxy = pbuffer.data(idx_fh + 64);

    auto ts_xyy_xxxxz = pbuffer.data(idx_fh + 65);

    auto ts_xyy_xxxyy = pbuffer.data(idx_fh + 66);

    auto ts_xyy_xxxyz = pbuffer.data(idx_fh + 67);

    auto ts_xyy_xxxzz = pbuffer.data(idx_fh + 68);

    auto ts_xyy_xxyyy = pbuffer.data(idx_fh + 69);

    auto ts_xyy_xxyyz = pbuffer.data(idx_fh + 70);

    auto ts_xyy_xxyzz = pbuffer.data(idx_fh + 71);

    auto ts_xyy_xxzzz = pbuffer.data(idx_fh + 72);

    auto ts_xyy_xyyyy = pbuffer.data(idx_fh + 73);

    auto ts_xyy_xyyyz = pbuffer.data(idx_fh + 74);

    auto ts_xyy_xyyzz = pbuffer.data(idx_fh + 75);

    auto ts_xyy_xyzzz = pbuffer.data(idx_fh + 76);

    auto ts_xyy_xzzzz = pbuffer.data(idx_fh + 77);

    auto ts_xyy_yyyyy = pbuffer.data(idx_fh + 78);

    auto ts_xyy_yyyyz = pbuffer.data(idx_fh + 79);

    auto ts_xyy_yyyzz = pbuffer.data(idx_fh + 80);

    auto ts_xyy_yyzzz = pbuffer.data(idx_fh + 81);

    auto ts_xyy_yzzzz = pbuffer.data(idx_fh + 82);

    auto ts_xyy_zzzzz = pbuffer.data(idx_fh + 83);

    auto ts_xyz_xxxxx = pbuffer.data(idx_fh + 84);

    auto ts_xyz_xxxxy = pbuffer.data(idx_fh + 85);

    auto ts_xyz_xxxxz = pbuffer.data(idx_fh + 86);

    auto ts_xyz_xxxyy = pbuffer.data(idx_fh + 87);

    auto ts_xyz_xxxyz = pbuffer.data(idx_fh + 88);

    auto ts_xyz_xxxzz = pbuffer.data(idx_fh + 89);

    auto ts_xyz_xxyyy = pbuffer.data(idx_fh + 90);

    auto ts_xyz_xxyyz = pbuffer.data(idx_fh + 91);

    auto ts_xyz_xxyzz = pbuffer.data(idx_fh + 92);

    auto ts_xyz_xxzzz = pbuffer.data(idx_fh + 93);

    auto ts_xyz_xyyyy = pbuffer.data(idx_fh + 94);

    auto ts_xyz_xyyyz = pbuffer.data(idx_fh + 95);

    auto ts_xyz_xyyzz = pbuffer.data(idx_fh + 96);

    auto ts_xyz_xyzzz = pbuffer.data(idx_fh + 97);

    auto ts_xyz_xzzzz = pbuffer.data(idx_fh + 98);

    auto ts_xyz_yyyyy = pbuffer.data(idx_fh + 99);

    auto ts_xyz_yyyyz = pbuffer.data(idx_fh + 100);

    auto ts_xyz_yyyzz = pbuffer.data(idx_fh + 101);

    auto ts_xyz_yyzzz = pbuffer.data(idx_fh + 102);

    auto ts_xyz_yzzzz = pbuffer.data(idx_fh + 103);

    auto ts_xyz_zzzzz = pbuffer.data(idx_fh + 104);

    auto ts_xzz_xxxxx = pbuffer.data(idx_fh + 105);

    auto ts_xzz_xxxxy = pbuffer.data(idx_fh + 106);

    auto ts_xzz_xxxxz = pbuffer.data(idx_fh + 107);

    auto ts_xzz_xxxyy = pbuffer.data(idx_fh + 108);

    auto ts_xzz_xxxyz = pbuffer.data(idx_fh + 109);

    auto ts_xzz_xxxzz = pbuffer.data(idx_fh + 110);

    auto ts_xzz_xxyyy = pbuffer.data(idx_fh + 111);

    auto ts_xzz_xxyyz = pbuffer.data(idx_fh + 112);

    auto ts_xzz_xxyzz = pbuffer.data(idx_fh + 113);

    auto ts_xzz_xxzzz = pbuffer.data(idx_fh + 114);

    auto ts_xzz_xyyyy = pbuffer.data(idx_fh + 115);

    auto ts_xzz_xyyyz = pbuffer.data(idx_fh + 116);

    auto ts_xzz_xyyzz = pbuffer.data(idx_fh + 117);

    auto ts_xzz_xyzzz = pbuffer.data(idx_fh + 118);

    auto ts_xzz_xzzzz = pbuffer.data(idx_fh + 119);

    auto ts_xzz_yyyyy = pbuffer.data(idx_fh + 120);

    auto ts_xzz_yyyyz = pbuffer.data(idx_fh + 121);

    auto ts_xzz_yyyzz = pbuffer.data(idx_fh + 122);

    auto ts_xzz_yyzzz = pbuffer.data(idx_fh + 123);

    auto ts_xzz_yzzzz = pbuffer.data(idx_fh + 124);

    auto ts_xzz_zzzzz = pbuffer.data(idx_fh + 125);

    auto ts_yyy_xxxxx = pbuffer.data(idx_fh + 126);

    auto ts_yyy_xxxxy = pbuffer.data(idx_fh + 127);

    auto ts_yyy_xxxxz = pbuffer.data(idx_fh + 128);

    auto ts_yyy_xxxyy = pbuffer.data(idx_fh + 129);

    auto ts_yyy_xxxyz = pbuffer.data(idx_fh + 130);

    auto ts_yyy_xxxzz = pbuffer.data(idx_fh + 131);

    auto ts_yyy_xxyyy = pbuffer.data(idx_fh + 132);

    auto ts_yyy_xxyyz = pbuffer.data(idx_fh + 133);

    auto ts_yyy_xxyzz = pbuffer.data(idx_fh + 134);

    auto ts_yyy_xxzzz = pbuffer.data(idx_fh + 135);

    auto ts_yyy_xyyyy = pbuffer.data(idx_fh + 136);

    auto ts_yyy_xyyyz = pbuffer.data(idx_fh + 137);

    auto ts_yyy_xyyzz = pbuffer.data(idx_fh + 138);

    auto ts_yyy_xyzzz = pbuffer.data(idx_fh + 139);

    auto ts_yyy_xzzzz = pbuffer.data(idx_fh + 140);

    auto ts_yyy_yyyyy = pbuffer.data(idx_fh + 141);

    auto ts_yyy_yyyyz = pbuffer.data(idx_fh + 142);

    auto ts_yyy_yyyzz = pbuffer.data(idx_fh + 143);

    auto ts_yyy_yyzzz = pbuffer.data(idx_fh + 144);

    auto ts_yyy_yzzzz = pbuffer.data(idx_fh + 145);

    auto ts_yyy_zzzzz = pbuffer.data(idx_fh + 146);

    auto ts_yyz_xxxxx = pbuffer.data(idx_fh + 147);

    auto ts_yyz_xxxxy = pbuffer.data(idx_fh + 148);

    auto ts_yyz_xxxxz = pbuffer.data(idx_fh + 149);

    auto ts_yyz_xxxyy = pbuffer.data(idx_fh + 150);

    auto ts_yyz_xxxyz = pbuffer.data(idx_fh + 151);

    auto ts_yyz_xxxzz = pbuffer.data(idx_fh + 152);

    auto ts_yyz_xxyyy = pbuffer.data(idx_fh + 153);

    auto ts_yyz_xxyyz = pbuffer.data(idx_fh + 154);

    auto ts_yyz_xxyzz = pbuffer.data(idx_fh + 155);

    auto ts_yyz_xxzzz = pbuffer.data(idx_fh + 156);

    auto ts_yyz_xyyyy = pbuffer.data(idx_fh + 157);

    auto ts_yyz_xyyyz = pbuffer.data(idx_fh + 158);

    auto ts_yyz_xyyzz = pbuffer.data(idx_fh + 159);

    auto ts_yyz_xyzzz = pbuffer.data(idx_fh + 160);

    auto ts_yyz_xzzzz = pbuffer.data(idx_fh + 161);

    auto ts_yyz_yyyyy = pbuffer.data(idx_fh + 162);

    auto ts_yyz_yyyyz = pbuffer.data(idx_fh + 163);

    auto ts_yyz_yyyzz = pbuffer.data(idx_fh + 164);

    auto ts_yyz_yyzzz = pbuffer.data(idx_fh + 165);

    auto ts_yyz_yzzzz = pbuffer.data(idx_fh + 166);

    auto ts_yyz_zzzzz = pbuffer.data(idx_fh + 167);

    auto ts_yzz_xxxxx = pbuffer.data(idx_fh + 168);

    auto ts_yzz_xxxxy = pbuffer.data(idx_fh + 169);

    auto ts_yzz_xxxxz = pbuffer.data(idx_fh + 170);

    auto ts_yzz_xxxyy = pbuffer.data(idx_fh + 171);

    auto ts_yzz_xxxyz = pbuffer.data(idx_fh + 172);

    auto ts_yzz_xxxzz = pbuffer.data(idx_fh + 173);

    auto ts_yzz_xxyyy = pbuffer.data(idx_fh + 174);

    auto ts_yzz_xxyyz = pbuffer.data(idx_fh + 175);

    auto ts_yzz_xxyzz = pbuffer.data(idx_fh + 176);

    auto ts_yzz_xxzzz = pbuffer.data(idx_fh + 177);

    auto ts_yzz_xyyyy = pbuffer.data(idx_fh + 178);

    auto ts_yzz_xyyyz = pbuffer.data(idx_fh + 179);

    auto ts_yzz_xyyzz = pbuffer.data(idx_fh + 180);

    auto ts_yzz_xyzzz = pbuffer.data(idx_fh + 181);

    auto ts_yzz_xzzzz = pbuffer.data(idx_fh + 182);

    auto ts_yzz_yyyyy = pbuffer.data(idx_fh + 183);

    auto ts_yzz_yyyyz = pbuffer.data(idx_fh + 184);

    auto ts_yzz_yyyzz = pbuffer.data(idx_fh + 185);

    auto ts_yzz_yyzzz = pbuffer.data(idx_fh + 186);

    auto ts_yzz_yzzzz = pbuffer.data(idx_fh + 187);

    auto ts_yzz_zzzzz = pbuffer.data(idx_fh + 188);

    auto ts_zzz_xxxxx = pbuffer.data(idx_fh + 189);

    auto ts_zzz_xxxxy = pbuffer.data(idx_fh + 190);

    auto ts_zzz_xxxxz = pbuffer.data(idx_fh + 191);

    auto ts_zzz_xxxyy = pbuffer.data(idx_fh + 192);

    auto ts_zzz_xxxyz = pbuffer.data(idx_fh + 193);

    auto ts_zzz_xxxzz = pbuffer.data(idx_fh + 194);

    auto ts_zzz_xxyyy = pbuffer.data(idx_fh + 195);

    auto ts_zzz_xxyyz = pbuffer.data(idx_fh + 196);

    auto ts_zzz_xxyzz = pbuffer.data(idx_fh + 197);

    auto ts_zzz_xxzzz = pbuffer.data(idx_fh + 198);

    auto ts_zzz_xyyyy = pbuffer.data(idx_fh + 199);

    auto ts_zzz_xyyyz = pbuffer.data(idx_fh + 200);

    auto ts_zzz_xyyzz = pbuffer.data(idx_fh + 201);

    auto ts_zzz_xyzzz = pbuffer.data(idx_fh + 202);

    auto ts_zzz_xzzzz = pbuffer.data(idx_fh + 203);

    auto ts_zzz_yyyyy = pbuffer.data(idx_fh + 204);

    auto ts_zzz_yyyyz = pbuffer.data(idx_fh + 205);

    auto ts_zzz_yyyzz = pbuffer.data(idx_fh + 206);

    auto ts_zzz_yyzzz = pbuffer.data(idx_fh + 207);

    auto ts_zzz_yzzzz = pbuffer.data(idx_fh + 208);

    auto ts_zzz_zzzzz = pbuffer.data(idx_fh + 209);

    // Set up 0-21 components of targeted buffer : FH

    auto gs_x_xxx_xxxxx = pbuffer.data(idx_g_fh);

    auto gs_x_xxx_xxxxy = pbuffer.data(idx_g_fh + 1);

    auto gs_x_xxx_xxxxz = pbuffer.data(idx_g_fh + 2);

    auto gs_x_xxx_xxxyy = pbuffer.data(idx_g_fh + 3);

    auto gs_x_xxx_xxxyz = pbuffer.data(idx_g_fh + 4);

    auto gs_x_xxx_xxxzz = pbuffer.data(idx_g_fh + 5);

    auto gs_x_xxx_xxyyy = pbuffer.data(idx_g_fh + 6);

    auto gs_x_xxx_xxyyz = pbuffer.data(idx_g_fh + 7);

    auto gs_x_xxx_xxyzz = pbuffer.data(idx_g_fh + 8);

    auto gs_x_xxx_xxzzz = pbuffer.data(idx_g_fh + 9);

    auto gs_x_xxx_xyyyy = pbuffer.data(idx_g_fh + 10);

    auto gs_x_xxx_xyyyz = pbuffer.data(idx_g_fh + 11);

    auto gs_x_xxx_xyyzz = pbuffer.data(idx_g_fh + 12);

    auto gs_x_xxx_xyzzz = pbuffer.data(idx_g_fh + 13);

    auto gs_x_xxx_xzzzz = pbuffer.data(idx_g_fh + 14);

    auto gs_x_xxx_yyyyy = pbuffer.data(idx_g_fh + 15);

    auto gs_x_xxx_yyyyz = pbuffer.data(idx_g_fh + 16);

    auto gs_x_xxx_yyyzz = pbuffer.data(idx_g_fh + 17);

    auto gs_x_xxx_yyzzz = pbuffer.data(idx_g_fh + 18);

    auto gs_x_xxx_yzzzz = pbuffer.data(idx_g_fh + 19);

    auto gs_x_xxx_zzzzz = pbuffer.data(idx_g_fh + 20);

    #pragma omp simd aligned(gc_x, gs_x_xxx_xxxxx, gs_x_xxx_xxxxy, gs_x_xxx_xxxxz, gs_x_xxx_xxxyy, gs_x_xxx_xxxyz, gs_x_xxx_xxxzz, gs_x_xxx_xxyyy, gs_x_xxx_xxyyz, gs_x_xxx_xxyzz, gs_x_xxx_xxzzz, gs_x_xxx_xyyyy, gs_x_xxx_xyyyz, gs_x_xxx_xyyzz, gs_x_xxx_xyzzz, gs_x_xxx_xzzzz, gs_x_xxx_yyyyy, gs_x_xxx_yyyyz, gs_x_xxx_yyyzz, gs_x_xxx_yyzzz, gs_x_xxx_yzzzz, gs_x_xxx_zzzzz, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxxz, ts_xx_xxxyy, ts_xx_xxxyz, ts_xx_xxxzz, ts_xx_xxyyy, ts_xx_xxyyz, ts_xx_xxyzz, ts_xx_xxzzz, ts_xx_xyyyy, ts_xx_xyyyz, ts_xx_xyyzz, ts_xx_xyzzz, ts_xx_xzzzz, ts_xx_yyyyy, ts_xx_yyyyz, ts_xx_yyyzz, ts_xx_yyzzz, ts_xx_yzzzz, ts_xx_zzzzz, ts_xxx_xxxx, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxy, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxz, ts_xxx_xxxzz, ts_xxx_xxyy, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyz, ts_xxx_xxyzz, ts_xxx_xxzz, ts_xxx_xxzzz, ts_xxx_xyyy, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyz, ts_xxx_xyyzz, ts_xxx_xyzz, ts_xxx_xyzzz, ts_xxx_xzzz, ts_xxx_xzzzz, ts_xxx_yyyy, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyz, ts_xxx_yyyzz, ts_xxx_yyzz, ts_xxx_yyzzz, ts_xxx_yzzz, ts_xxx_yzzzz, ts_xxx_zzzz, ts_xxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxx_xxxxx[i] = 6.0 * ts_xx_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxxxy[i] = 6.0 * ts_xx_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxxxz[i] = 6.0 * ts_xx_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxxyy[i] = 6.0 * ts_xx_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxxyz[i] = 6.0 * ts_xx_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxxzz[i] = 6.0 * ts_xx_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxyyy[i] = 6.0 * ts_xx_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxyyz[i] = 6.0 * ts_xx_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxyzz[i] = 6.0 * ts_xx_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xxzzz[i] = 6.0 * ts_xx_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xyyyy[i] = 6.0 * ts_xx_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xyyyz[i] = 6.0 * ts_xx_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xyyzz[i] = 6.0 * ts_xx_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xyzzz[i] = 6.0 * ts_xx_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_xzzzz[i] = 6.0 * ts_xx_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yyyyy[i] = 6.0 * ts_xx_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxx_yyyyz[i] = 6.0 * ts_xx_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yyyzz[i] = 6.0 * ts_xx_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yyzzz[i] = 6.0 * ts_xx_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yzzzz[i] = 6.0 * ts_xx_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxx_zzzzz[i] = 6.0 * ts_xx_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 21-42 components of targeted buffer : FH

    auto gs_x_xxy_xxxxx = pbuffer.data(idx_g_fh + 21);

    auto gs_x_xxy_xxxxy = pbuffer.data(idx_g_fh + 22);

    auto gs_x_xxy_xxxxz = pbuffer.data(idx_g_fh + 23);

    auto gs_x_xxy_xxxyy = pbuffer.data(idx_g_fh + 24);

    auto gs_x_xxy_xxxyz = pbuffer.data(idx_g_fh + 25);

    auto gs_x_xxy_xxxzz = pbuffer.data(idx_g_fh + 26);

    auto gs_x_xxy_xxyyy = pbuffer.data(idx_g_fh + 27);

    auto gs_x_xxy_xxyyz = pbuffer.data(idx_g_fh + 28);

    auto gs_x_xxy_xxyzz = pbuffer.data(idx_g_fh + 29);

    auto gs_x_xxy_xxzzz = pbuffer.data(idx_g_fh + 30);

    auto gs_x_xxy_xyyyy = pbuffer.data(idx_g_fh + 31);

    auto gs_x_xxy_xyyyz = pbuffer.data(idx_g_fh + 32);

    auto gs_x_xxy_xyyzz = pbuffer.data(idx_g_fh + 33);

    auto gs_x_xxy_xyzzz = pbuffer.data(idx_g_fh + 34);

    auto gs_x_xxy_xzzzz = pbuffer.data(idx_g_fh + 35);

    auto gs_x_xxy_yyyyy = pbuffer.data(idx_g_fh + 36);

    auto gs_x_xxy_yyyyz = pbuffer.data(idx_g_fh + 37);

    auto gs_x_xxy_yyyzz = pbuffer.data(idx_g_fh + 38);

    auto gs_x_xxy_yyzzz = pbuffer.data(idx_g_fh + 39);

    auto gs_x_xxy_yzzzz = pbuffer.data(idx_g_fh + 40);

    auto gs_x_xxy_zzzzz = pbuffer.data(idx_g_fh + 41);

    #pragma omp simd aligned(gc_x, gs_x_xxy_xxxxx, gs_x_xxy_xxxxy, gs_x_xxy_xxxxz, gs_x_xxy_xxxyy, gs_x_xxy_xxxyz, gs_x_xxy_xxxzz, gs_x_xxy_xxyyy, gs_x_xxy_xxyyz, gs_x_xxy_xxyzz, gs_x_xxy_xxzzz, gs_x_xxy_xyyyy, gs_x_xxy_xyyyz, gs_x_xxy_xyyzz, gs_x_xxy_xyzzz, gs_x_xxy_xzzzz, gs_x_xxy_yyyyy, gs_x_xxy_yyyyz, gs_x_xxy_yyyzz, gs_x_xxy_yyzzz, gs_x_xxy_yzzzz, gs_x_xxy_zzzzz, ts_xxy_xxxx, ts_xxy_xxxxx, ts_xxy_xxxxy, ts_xxy_xxxxz, ts_xxy_xxxy, ts_xxy_xxxyy, ts_xxy_xxxyz, ts_xxy_xxxz, ts_xxy_xxxzz, ts_xxy_xxyy, ts_xxy_xxyyy, ts_xxy_xxyyz, ts_xxy_xxyz, ts_xxy_xxyzz, ts_xxy_xxzz, ts_xxy_xxzzz, ts_xxy_xyyy, ts_xxy_xyyyy, ts_xxy_xyyyz, ts_xxy_xyyz, ts_xxy_xyyzz, ts_xxy_xyzz, ts_xxy_xyzzz, ts_xxy_xzzz, ts_xxy_xzzzz, ts_xxy_yyyy, ts_xxy_yyyyy, ts_xxy_yyyyz, ts_xxy_yyyz, ts_xxy_yyyzz, ts_xxy_yyzz, ts_xxy_yyzzz, ts_xxy_yzzz, ts_xxy_yzzzz, ts_xxy_zzzz, ts_xxy_zzzzz, ts_xy_xxxxx, ts_xy_xxxxy, ts_xy_xxxxz, ts_xy_xxxyy, ts_xy_xxxyz, ts_xy_xxxzz, ts_xy_xxyyy, ts_xy_xxyyz, ts_xy_xxyzz, ts_xy_xxzzz, ts_xy_xyyyy, ts_xy_xyyyz, ts_xy_xyyzz, ts_xy_xyzzz, ts_xy_xzzzz, ts_xy_yyyyy, ts_xy_yyyyz, ts_xy_yyyzz, ts_xy_yyzzz, ts_xy_yzzzz, ts_xy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxy_xxxxx[i] = 4.0 * ts_xy_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxxxy[i] = 4.0 * ts_xy_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxxxz[i] = 4.0 * ts_xy_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxxyy[i] = 4.0 * ts_xy_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxxyz[i] = 4.0 * ts_xy_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxxzz[i] = 4.0 * ts_xy_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxyyy[i] = 4.0 * ts_xy_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxyyz[i] = 4.0 * ts_xy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxyzz[i] = 4.0 * ts_xy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xxzzz[i] = 4.0 * ts_xy_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xyyyy[i] = 4.0 * ts_xy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xyyyz[i] = 4.0 * ts_xy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xyyzz[i] = 4.0 * ts_xy_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xyzzz[i] = 4.0 * ts_xy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_xzzzz[i] = 4.0 * ts_xy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yyyyy[i] = 4.0 * ts_xy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxy_yyyyz[i] = 4.0 * ts_xy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yyyzz[i] = 4.0 * ts_xy_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yyzzz[i] = 4.0 * ts_xy_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yzzzz[i] = 4.0 * ts_xy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxy_zzzzz[i] = 4.0 * ts_xy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 42-63 components of targeted buffer : FH

    auto gs_x_xxz_xxxxx = pbuffer.data(idx_g_fh + 42);

    auto gs_x_xxz_xxxxy = pbuffer.data(idx_g_fh + 43);

    auto gs_x_xxz_xxxxz = pbuffer.data(idx_g_fh + 44);

    auto gs_x_xxz_xxxyy = pbuffer.data(idx_g_fh + 45);

    auto gs_x_xxz_xxxyz = pbuffer.data(idx_g_fh + 46);

    auto gs_x_xxz_xxxzz = pbuffer.data(idx_g_fh + 47);

    auto gs_x_xxz_xxyyy = pbuffer.data(idx_g_fh + 48);

    auto gs_x_xxz_xxyyz = pbuffer.data(idx_g_fh + 49);

    auto gs_x_xxz_xxyzz = pbuffer.data(idx_g_fh + 50);

    auto gs_x_xxz_xxzzz = pbuffer.data(idx_g_fh + 51);

    auto gs_x_xxz_xyyyy = pbuffer.data(idx_g_fh + 52);

    auto gs_x_xxz_xyyyz = pbuffer.data(idx_g_fh + 53);

    auto gs_x_xxz_xyyzz = pbuffer.data(idx_g_fh + 54);

    auto gs_x_xxz_xyzzz = pbuffer.data(idx_g_fh + 55);

    auto gs_x_xxz_xzzzz = pbuffer.data(idx_g_fh + 56);

    auto gs_x_xxz_yyyyy = pbuffer.data(idx_g_fh + 57);

    auto gs_x_xxz_yyyyz = pbuffer.data(idx_g_fh + 58);

    auto gs_x_xxz_yyyzz = pbuffer.data(idx_g_fh + 59);

    auto gs_x_xxz_yyzzz = pbuffer.data(idx_g_fh + 60);

    auto gs_x_xxz_yzzzz = pbuffer.data(idx_g_fh + 61);

    auto gs_x_xxz_zzzzz = pbuffer.data(idx_g_fh + 62);

    #pragma omp simd aligned(gc_x, gs_x_xxz_xxxxx, gs_x_xxz_xxxxy, gs_x_xxz_xxxxz, gs_x_xxz_xxxyy, gs_x_xxz_xxxyz, gs_x_xxz_xxxzz, gs_x_xxz_xxyyy, gs_x_xxz_xxyyz, gs_x_xxz_xxyzz, gs_x_xxz_xxzzz, gs_x_xxz_xyyyy, gs_x_xxz_xyyyz, gs_x_xxz_xyyzz, gs_x_xxz_xyzzz, gs_x_xxz_xzzzz, gs_x_xxz_yyyyy, gs_x_xxz_yyyyz, gs_x_xxz_yyyzz, gs_x_xxz_yyzzz, gs_x_xxz_yzzzz, gs_x_xxz_zzzzz, ts_xxz_xxxx, ts_xxz_xxxxx, ts_xxz_xxxxy, ts_xxz_xxxxz, ts_xxz_xxxy, ts_xxz_xxxyy, ts_xxz_xxxyz, ts_xxz_xxxz, ts_xxz_xxxzz, ts_xxz_xxyy, ts_xxz_xxyyy, ts_xxz_xxyyz, ts_xxz_xxyz, ts_xxz_xxyzz, ts_xxz_xxzz, ts_xxz_xxzzz, ts_xxz_xyyy, ts_xxz_xyyyy, ts_xxz_xyyyz, ts_xxz_xyyz, ts_xxz_xyyzz, ts_xxz_xyzz, ts_xxz_xyzzz, ts_xxz_xzzz, ts_xxz_xzzzz, ts_xxz_yyyy, ts_xxz_yyyyy, ts_xxz_yyyyz, ts_xxz_yyyz, ts_xxz_yyyzz, ts_xxz_yyzz, ts_xxz_yyzzz, ts_xxz_yzzz, ts_xxz_yzzzz, ts_xxz_zzzz, ts_xxz_zzzzz, ts_xz_xxxxx, ts_xz_xxxxy, ts_xz_xxxxz, ts_xz_xxxyy, ts_xz_xxxyz, ts_xz_xxxzz, ts_xz_xxyyy, ts_xz_xxyyz, ts_xz_xxyzz, ts_xz_xxzzz, ts_xz_xyyyy, ts_xz_xyyyz, ts_xz_xyyzz, ts_xz_xyzzz, ts_xz_xzzzz, ts_xz_yyyyy, ts_xz_yyyyz, ts_xz_yyyzz, ts_xz_yyzzz, ts_xz_yzzzz, ts_xz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxz_xxxxx[i] = 4.0 * ts_xz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxxxy[i] = 4.0 * ts_xz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxxxz[i] = 4.0 * ts_xz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxxyy[i] = 4.0 * ts_xz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxxyz[i] = 4.0 * ts_xz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxxzz[i] = 4.0 * ts_xz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxyyy[i] = 4.0 * ts_xz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxyyz[i] = 4.0 * ts_xz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxyzz[i] = 4.0 * ts_xz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xxzzz[i] = 4.0 * ts_xz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xyyyy[i] = 4.0 * ts_xz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xyyyz[i] = 4.0 * ts_xz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xyyzz[i] = 4.0 * ts_xz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xyzzz[i] = 4.0 * ts_xz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_xzzzz[i] = 4.0 * ts_xz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yyyyy[i] = 4.0 * ts_xz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xxz_yyyyz[i] = 4.0 * ts_xz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yyyzz[i] = 4.0 * ts_xz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yyzzz[i] = 4.0 * ts_xz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yzzzz[i] = 4.0 * ts_xz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xxz_zzzzz[i] = 4.0 * ts_xz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 63-84 components of targeted buffer : FH

    auto gs_x_xyy_xxxxx = pbuffer.data(idx_g_fh + 63);

    auto gs_x_xyy_xxxxy = pbuffer.data(idx_g_fh + 64);

    auto gs_x_xyy_xxxxz = pbuffer.data(idx_g_fh + 65);

    auto gs_x_xyy_xxxyy = pbuffer.data(idx_g_fh + 66);

    auto gs_x_xyy_xxxyz = pbuffer.data(idx_g_fh + 67);

    auto gs_x_xyy_xxxzz = pbuffer.data(idx_g_fh + 68);

    auto gs_x_xyy_xxyyy = pbuffer.data(idx_g_fh + 69);

    auto gs_x_xyy_xxyyz = pbuffer.data(idx_g_fh + 70);

    auto gs_x_xyy_xxyzz = pbuffer.data(idx_g_fh + 71);

    auto gs_x_xyy_xxzzz = pbuffer.data(idx_g_fh + 72);

    auto gs_x_xyy_xyyyy = pbuffer.data(idx_g_fh + 73);

    auto gs_x_xyy_xyyyz = pbuffer.data(idx_g_fh + 74);

    auto gs_x_xyy_xyyzz = pbuffer.data(idx_g_fh + 75);

    auto gs_x_xyy_xyzzz = pbuffer.data(idx_g_fh + 76);

    auto gs_x_xyy_xzzzz = pbuffer.data(idx_g_fh + 77);

    auto gs_x_xyy_yyyyy = pbuffer.data(idx_g_fh + 78);

    auto gs_x_xyy_yyyyz = pbuffer.data(idx_g_fh + 79);

    auto gs_x_xyy_yyyzz = pbuffer.data(idx_g_fh + 80);

    auto gs_x_xyy_yyzzz = pbuffer.data(idx_g_fh + 81);

    auto gs_x_xyy_yzzzz = pbuffer.data(idx_g_fh + 82);

    auto gs_x_xyy_zzzzz = pbuffer.data(idx_g_fh + 83);

    #pragma omp simd aligned(gc_x, gs_x_xyy_xxxxx, gs_x_xyy_xxxxy, gs_x_xyy_xxxxz, gs_x_xyy_xxxyy, gs_x_xyy_xxxyz, gs_x_xyy_xxxzz, gs_x_xyy_xxyyy, gs_x_xyy_xxyyz, gs_x_xyy_xxyzz, gs_x_xyy_xxzzz, gs_x_xyy_xyyyy, gs_x_xyy_xyyyz, gs_x_xyy_xyyzz, gs_x_xyy_xyzzz, gs_x_xyy_xzzzz, gs_x_xyy_yyyyy, gs_x_xyy_yyyyz, gs_x_xyy_yyyzz, gs_x_xyy_yyzzz, gs_x_xyy_yzzzz, gs_x_xyy_zzzzz, ts_xyy_xxxx, ts_xyy_xxxxx, ts_xyy_xxxxy, ts_xyy_xxxxz, ts_xyy_xxxy, ts_xyy_xxxyy, ts_xyy_xxxyz, ts_xyy_xxxz, ts_xyy_xxxzz, ts_xyy_xxyy, ts_xyy_xxyyy, ts_xyy_xxyyz, ts_xyy_xxyz, ts_xyy_xxyzz, ts_xyy_xxzz, ts_xyy_xxzzz, ts_xyy_xyyy, ts_xyy_xyyyy, ts_xyy_xyyyz, ts_xyy_xyyz, ts_xyy_xyyzz, ts_xyy_xyzz, ts_xyy_xyzzz, ts_xyy_xzzz, ts_xyy_xzzzz, ts_xyy_yyyy, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyz, ts_xyy_yyyzz, ts_xyy_yyzz, ts_xyy_yyzzz, ts_xyy_yzzz, ts_xyy_yzzzz, ts_xyy_zzzz, ts_xyy_zzzzz, ts_yy_xxxxx, ts_yy_xxxxy, ts_yy_xxxxz, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxxzz, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyzz, ts_yy_xxzzz, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyzz, ts_yy_xyzzz, ts_yy_xzzzz, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyzz, ts_yy_yyzzz, ts_yy_yzzzz, ts_yy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyy_xxxxx[i] = 2.0 * ts_yy_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxxxy[i] = 2.0 * ts_yy_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxxxz[i] = 2.0 * ts_yy_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxxyy[i] = 2.0 * ts_yy_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxxyz[i] = 2.0 * ts_yy_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxxzz[i] = 2.0 * ts_yy_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxyyy[i] = 2.0 * ts_yy_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxyyz[i] = 2.0 * ts_yy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxyzz[i] = 2.0 * ts_yy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xxzzz[i] = 2.0 * ts_yy_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xyyyy[i] = 2.0 * ts_yy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xyyyz[i] = 2.0 * ts_yy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xyyzz[i] = 2.0 * ts_yy_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xyzzz[i] = 2.0 * ts_yy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_xzzzz[i] = 2.0 * ts_yy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yyyyy[i] = 2.0 * ts_yy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyy_yyyyz[i] = 2.0 * ts_yy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yyyzz[i] = 2.0 * ts_yy_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yyzzz[i] = 2.0 * ts_yy_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yzzzz[i] = 2.0 * ts_yy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyy_zzzzz[i] = 2.0 * ts_yy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 84-105 components of targeted buffer : FH

    auto gs_x_xyz_xxxxx = pbuffer.data(idx_g_fh + 84);

    auto gs_x_xyz_xxxxy = pbuffer.data(idx_g_fh + 85);

    auto gs_x_xyz_xxxxz = pbuffer.data(idx_g_fh + 86);

    auto gs_x_xyz_xxxyy = pbuffer.data(idx_g_fh + 87);

    auto gs_x_xyz_xxxyz = pbuffer.data(idx_g_fh + 88);

    auto gs_x_xyz_xxxzz = pbuffer.data(idx_g_fh + 89);

    auto gs_x_xyz_xxyyy = pbuffer.data(idx_g_fh + 90);

    auto gs_x_xyz_xxyyz = pbuffer.data(idx_g_fh + 91);

    auto gs_x_xyz_xxyzz = pbuffer.data(idx_g_fh + 92);

    auto gs_x_xyz_xxzzz = pbuffer.data(idx_g_fh + 93);

    auto gs_x_xyz_xyyyy = pbuffer.data(idx_g_fh + 94);

    auto gs_x_xyz_xyyyz = pbuffer.data(idx_g_fh + 95);

    auto gs_x_xyz_xyyzz = pbuffer.data(idx_g_fh + 96);

    auto gs_x_xyz_xyzzz = pbuffer.data(idx_g_fh + 97);

    auto gs_x_xyz_xzzzz = pbuffer.data(idx_g_fh + 98);

    auto gs_x_xyz_yyyyy = pbuffer.data(idx_g_fh + 99);

    auto gs_x_xyz_yyyyz = pbuffer.data(idx_g_fh + 100);

    auto gs_x_xyz_yyyzz = pbuffer.data(idx_g_fh + 101);

    auto gs_x_xyz_yyzzz = pbuffer.data(idx_g_fh + 102);

    auto gs_x_xyz_yzzzz = pbuffer.data(idx_g_fh + 103);

    auto gs_x_xyz_zzzzz = pbuffer.data(idx_g_fh + 104);

    #pragma omp simd aligned(gc_x, gs_x_xyz_xxxxx, gs_x_xyz_xxxxy, gs_x_xyz_xxxxz, gs_x_xyz_xxxyy, gs_x_xyz_xxxyz, gs_x_xyz_xxxzz, gs_x_xyz_xxyyy, gs_x_xyz_xxyyz, gs_x_xyz_xxyzz, gs_x_xyz_xxzzz, gs_x_xyz_xyyyy, gs_x_xyz_xyyyz, gs_x_xyz_xyyzz, gs_x_xyz_xyzzz, gs_x_xyz_xzzzz, gs_x_xyz_yyyyy, gs_x_xyz_yyyyz, gs_x_xyz_yyyzz, gs_x_xyz_yyzzz, gs_x_xyz_yzzzz, gs_x_xyz_zzzzz, ts_xyz_xxxx, ts_xyz_xxxxx, ts_xyz_xxxxy, ts_xyz_xxxxz, ts_xyz_xxxy, ts_xyz_xxxyy, ts_xyz_xxxyz, ts_xyz_xxxz, ts_xyz_xxxzz, ts_xyz_xxyy, ts_xyz_xxyyy, ts_xyz_xxyyz, ts_xyz_xxyz, ts_xyz_xxyzz, ts_xyz_xxzz, ts_xyz_xxzzz, ts_xyz_xyyy, ts_xyz_xyyyy, ts_xyz_xyyyz, ts_xyz_xyyz, ts_xyz_xyyzz, ts_xyz_xyzz, ts_xyz_xyzzz, ts_xyz_xzzz, ts_xyz_xzzzz, ts_xyz_yyyy, ts_xyz_yyyyy, ts_xyz_yyyyz, ts_xyz_yyyz, ts_xyz_yyyzz, ts_xyz_yyzz, ts_xyz_yyzzz, ts_xyz_yzzz, ts_xyz_yzzzz, ts_xyz_zzzz, ts_xyz_zzzzz, ts_yz_xxxxx, ts_yz_xxxxy, ts_yz_xxxxz, ts_yz_xxxyy, ts_yz_xxxyz, ts_yz_xxxzz, ts_yz_xxyyy, ts_yz_xxyyz, ts_yz_xxyzz, ts_yz_xxzzz, ts_yz_xyyyy, ts_yz_xyyyz, ts_yz_xyyzz, ts_yz_xyzzz, ts_yz_xzzzz, ts_yz_yyyyy, ts_yz_yyyyz, ts_yz_yyyzz, ts_yz_yyzzz, ts_yz_yzzzz, ts_yz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyz_xxxxx[i] = 2.0 * ts_yz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxxxy[i] = 2.0 * ts_yz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxxxz[i] = 2.0 * ts_yz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxxyy[i] = 2.0 * ts_yz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxxyz[i] = 2.0 * ts_yz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxxzz[i] = 2.0 * ts_yz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxyyy[i] = 2.0 * ts_yz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxyyz[i] = 2.0 * ts_yz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxyzz[i] = 2.0 * ts_yz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xxzzz[i] = 2.0 * ts_yz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xyyyy[i] = 2.0 * ts_yz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xyyyz[i] = 2.0 * ts_yz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xyyzz[i] = 2.0 * ts_yz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xyzzz[i] = 2.0 * ts_yz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_xzzzz[i] = 2.0 * ts_yz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yyyyy[i] = 2.0 * ts_yz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xyz_yyyyz[i] = 2.0 * ts_yz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yyyzz[i] = 2.0 * ts_yz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yyzzz[i] = 2.0 * ts_yz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yzzzz[i] = 2.0 * ts_yz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xyz_zzzzz[i] = 2.0 * ts_yz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 105-126 components of targeted buffer : FH

    auto gs_x_xzz_xxxxx = pbuffer.data(idx_g_fh + 105);

    auto gs_x_xzz_xxxxy = pbuffer.data(idx_g_fh + 106);

    auto gs_x_xzz_xxxxz = pbuffer.data(idx_g_fh + 107);

    auto gs_x_xzz_xxxyy = pbuffer.data(idx_g_fh + 108);

    auto gs_x_xzz_xxxyz = pbuffer.data(idx_g_fh + 109);

    auto gs_x_xzz_xxxzz = pbuffer.data(idx_g_fh + 110);

    auto gs_x_xzz_xxyyy = pbuffer.data(idx_g_fh + 111);

    auto gs_x_xzz_xxyyz = pbuffer.data(idx_g_fh + 112);

    auto gs_x_xzz_xxyzz = pbuffer.data(idx_g_fh + 113);

    auto gs_x_xzz_xxzzz = pbuffer.data(idx_g_fh + 114);

    auto gs_x_xzz_xyyyy = pbuffer.data(idx_g_fh + 115);

    auto gs_x_xzz_xyyyz = pbuffer.data(idx_g_fh + 116);

    auto gs_x_xzz_xyyzz = pbuffer.data(idx_g_fh + 117);

    auto gs_x_xzz_xyzzz = pbuffer.data(idx_g_fh + 118);

    auto gs_x_xzz_xzzzz = pbuffer.data(idx_g_fh + 119);

    auto gs_x_xzz_yyyyy = pbuffer.data(idx_g_fh + 120);

    auto gs_x_xzz_yyyyz = pbuffer.data(idx_g_fh + 121);

    auto gs_x_xzz_yyyzz = pbuffer.data(idx_g_fh + 122);

    auto gs_x_xzz_yyzzz = pbuffer.data(idx_g_fh + 123);

    auto gs_x_xzz_yzzzz = pbuffer.data(idx_g_fh + 124);

    auto gs_x_xzz_zzzzz = pbuffer.data(idx_g_fh + 125);

    #pragma omp simd aligned(gc_x, gs_x_xzz_xxxxx, gs_x_xzz_xxxxy, gs_x_xzz_xxxxz, gs_x_xzz_xxxyy, gs_x_xzz_xxxyz, gs_x_xzz_xxxzz, gs_x_xzz_xxyyy, gs_x_xzz_xxyyz, gs_x_xzz_xxyzz, gs_x_xzz_xxzzz, gs_x_xzz_xyyyy, gs_x_xzz_xyyyz, gs_x_xzz_xyyzz, gs_x_xzz_xyzzz, gs_x_xzz_xzzzz, gs_x_xzz_yyyyy, gs_x_xzz_yyyyz, gs_x_xzz_yyyzz, gs_x_xzz_yyzzz, gs_x_xzz_yzzzz, gs_x_xzz_zzzzz, ts_xzz_xxxx, ts_xzz_xxxxx, ts_xzz_xxxxy, ts_xzz_xxxxz, ts_xzz_xxxy, ts_xzz_xxxyy, ts_xzz_xxxyz, ts_xzz_xxxz, ts_xzz_xxxzz, ts_xzz_xxyy, ts_xzz_xxyyy, ts_xzz_xxyyz, ts_xzz_xxyz, ts_xzz_xxyzz, ts_xzz_xxzz, ts_xzz_xxzzz, ts_xzz_xyyy, ts_xzz_xyyyy, ts_xzz_xyyyz, ts_xzz_xyyz, ts_xzz_xyyzz, ts_xzz_xyzz, ts_xzz_xyzzz, ts_xzz_xzzz, ts_xzz_xzzzz, ts_xzz_yyyy, ts_xzz_yyyyy, ts_xzz_yyyyz, ts_xzz_yyyz, ts_xzz_yyyzz, ts_xzz_yyzz, ts_xzz_yyzzz, ts_xzz_yzzz, ts_xzz_yzzzz, ts_xzz_zzzz, ts_xzz_zzzzz, ts_zz_xxxxx, ts_zz_xxxxy, ts_zz_xxxxz, ts_zz_xxxyy, ts_zz_xxxyz, ts_zz_xxxzz, ts_zz_xxyyy, ts_zz_xxyyz, ts_zz_xxyzz, ts_zz_xxzzz, ts_zz_xyyyy, ts_zz_xyyyz, ts_zz_xyyzz, ts_zz_xyzzz, ts_zz_xzzzz, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyzz, ts_zz_yyzzz, ts_zz_yzzzz, ts_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzz_xxxxx[i] = 2.0 * ts_zz_xxxxx[i] * gfe_0 * tce_0 + 10.0 * ts_xzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxxxy[i] = 2.0 * ts_zz_xxxxy[i] * gfe_0 * tce_0 + 8.0 * ts_xzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxxxz[i] = 2.0 * ts_zz_xxxxz[i] * gfe_0 * tce_0 + 8.0 * ts_xzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxxyy[i] = 2.0 * ts_zz_xxxyy[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxxyz[i] = 2.0 * ts_zz_xxxyz[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxxzz[i] = 2.0 * ts_zz_xxxzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxyyy[i] = 2.0 * ts_zz_xxyyy[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxyyz[i] = 2.0 * ts_zz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxyzz[i] = 2.0 * ts_zz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xxzzz[i] = 2.0 * ts_zz_xxzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xyyyy[i] = 2.0 * ts_zz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xyyyz[i] = 2.0 * ts_zz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xyyzz[i] = 2.0 * ts_zz_xyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xyzzz[i] = 2.0 * ts_zz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_xzzzz[i] = 2.0 * ts_zz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yyyyy[i] = 2.0 * ts_zz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_xzz_yyyyz[i] = 2.0 * ts_zz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yyyzz[i] = 2.0 * ts_zz_yyyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yyzzz[i] = 2.0 * ts_zz_yyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yzzzz[i] = 2.0 * ts_zz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_xzz_zzzzz[i] = 2.0 * ts_zz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 126-147 components of targeted buffer : FH

    auto gs_x_yyy_xxxxx = pbuffer.data(idx_g_fh + 126);

    auto gs_x_yyy_xxxxy = pbuffer.data(idx_g_fh + 127);

    auto gs_x_yyy_xxxxz = pbuffer.data(idx_g_fh + 128);

    auto gs_x_yyy_xxxyy = pbuffer.data(idx_g_fh + 129);

    auto gs_x_yyy_xxxyz = pbuffer.data(idx_g_fh + 130);

    auto gs_x_yyy_xxxzz = pbuffer.data(idx_g_fh + 131);

    auto gs_x_yyy_xxyyy = pbuffer.data(idx_g_fh + 132);

    auto gs_x_yyy_xxyyz = pbuffer.data(idx_g_fh + 133);

    auto gs_x_yyy_xxyzz = pbuffer.data(idx_g_fh + 134);

    auto gs_x_yyy_xxzzz = pbuffer.data(idx_g_fh + 135);

    auto gs_x_yyy_xyyyy = pbuffer.data(idx_g_fh + 136);

    auto gs_x_yyy_xyyyz = pbuffer.data(idx_g_fh + 137);

    auto gs_x_yyy_xyyzz = pbuffer.data(idx_g_fh + 138);

    auto gs_x_yyy_xyzzz = pbuffer.data(idx_g_fh + 139);

    auto gs_x_yyy_xzzzz = pbuffer.data(idx_g_fh + 140);

    auto gs_x_yyy_yyyyy = pbuffer.data(idx_g_fh + 141);

    auto gs_x_yyy_yyyyz = pbuffer.data(idx_g_fh + 142);

    auto gs_x_yyy_yyyzz = pbuffer.data(idx_g_fh + 143);

    auto gs_x_yyy_yyzzz = pbuffer.data(idx_g_fh + 144);

    auto gs_x_yyy_yzzzz = pbuffer.data(idx_g_fh + 145);

    auto gs_x_yyy_zzzzz = pbuffer.data(idx_g_fh + 146);

    #pragma omp simd aligned(gc_x, gs_x_yyy_xxxxx, gs_x_yyy_xxxxy, gs_x_yyy_xxxxz, gs_x_yyy_xxxyy, gs_x_yyy_xxxyz, gs_x_yyy_xxxzz, gs_x_yyy_xxyyy, gs_x_yyy_xxyyz, gs_x_yyy_xxyzz, gs_x_yyy_xxzzz, gs_x_yyy_xyyyy, gs_x_yyy_xyyyz, gs_x_yyy_xyyzz, gs_x_yyy_xyzzz, gs_x_yyy_xzzzz, gs_x_yyy_yyyyy, gs_x_yyy_yyyyz, gs_x_yyy_yyyzz, gs_x_yyy_yyzzz, gs_x_yyy_yzzzz, gs_x_yyy_zzzzz, ts_yyy_xxxx, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxy, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxz, ts_yyy_xxxzz, ts_yyy_xxyy, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyz, ts_yyy_xxyzz, ts_yyy_xxzz, ts_yyy_xxzzz, ts_yyy_xyyy, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyz, ts_yyy_xyyzz, ts_yyy_xyzz, ts_yyy_xyzzz, ts_yyy_xzzz, ts_yyy_xzzzz, ts_yyy_yyyy, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyz, ts_yyy_yyyzz, ts_yyy_yyzz, ts_yyy_yyzzz, ts_yyy_yzzz, ts_yyy_yzzzz, ts_yyy_zzzz, ts_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyy_xxxxx[i] = 10.0 * ts_yyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxxxy[i] = 8.0 * ts_yyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxxxz[i] = 8.0 * ts_yyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxxyy[i] = 6.0 * ts_yyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxxyz[i] = 6.0 * ts_yyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxxzz[i] = 6.0 * ts_yyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxyyy[i] = 4.0 * ts_yyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxyyz[i] = 4.0 * ts_yyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxyzz[i] = 4.0 * ts_yyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xxzzz[i] = 4.0 * ts_yyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xyyyy[i] = 2.0 * ts_yyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xyyyz[i] = 2.0 * ts_yyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xyyzz[i] = 2.0 * ts_yyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xyzzz[i] = 2.0 * ts_yyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_xzzzz[i] = 2.0 * ts_yyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yyyyy[i] = 2.0 * ts_yyy_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyy_yyyyz[i] = 2.0 * ts_yyy_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yyyzz[i] = 2.0 * ts_yyy_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yyzzz[i] = 2.0 * ts_yyy_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yzzzz[i] = 2.0 * ts_yyy_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyy_zzzzz[i] = 2.0 * ts_yyy_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 147-168 components of targeted buffer : FH

    auto gs_x_yyz_xxxxx = pbuffer.data(idx_g_fh + 147);

    auto gs_x_yyz_xxxxy = pbuffer.data(idx_g_fh + 148);

    auto gs_x_yyz_xxxxz = pbuffer.data(idx_g_fh + 149);

    auto gs_x_yyz_xxxyy = pbuffer.data(idx_g_fh + 150);

    auto gs_x_yyz_xxxyz = pbuffer.data(idx_g_fh + 151);

    auto gs_x_yyz_xxxzz = pbuffer.data(idx_g_fh + 152);

    auto gs_x_yyz_xxyyy = pbuffer.data(idx_g_fh + 153);

    auto gs_x_yyz_xxyyz = pbuffer.data(idx_g_fh + 154);

    auto gs_x_yyz_xxyzz = pbuffer.data(idx_g_fh + 155);

    auto gs_x_yyz_xxzzz = pbuffer.data(idx_g_fh + 156);

    auto gs_x_yyz_xyyyy = pbuffer.data(idx_g_fh + 157);

    auto gs_x_yyz_xyyyz = pbuffer.data(idx_g_fh + 158);

    auto gs_x_yyz_xyyzz = pbuffer.data(idx_g_fh + 159);

    auto gs_x_yyz_xyzzz = pbuffer.data(idx_g_fh + 160);

    auto gs_x_yyz_xzzzz = pbuffer.data(idx_g_fh + 161);

    auto gs_x_yyz_yyyyy = pbuffer.data(idx_g_fh + 162);

    auto gs_x_yyz_yyyyz = pbuffer.data(idx_g_fh + 163);

    auto gs_x_yyz_yyyzz = pbuffer.data(idx_g_fh + 164);

    auto gs_x_yyz_yyzzz = pbuffer.data(idx_g_fh + 165);

    auto gs_x_yyz_yzzzz = pbuffer.data(idx_g_fh + 166);

    auto gs_x_yyz_zzzzz = pbuffer.data(idx_g_fh + 167);

    #pragma omp simd aligned(gc_x, gs_x_yyz_xxxxx, gs_x_yyz_xxxxy, gs_x_yyz_xxxxz, gs_x_yyz_xxxyy, gs_x_yyz_xxxyz, gs_x_yyz_xxxzz, gs_x_yyz_xxyyy, gs_x_yyz_xxyyz, gs_x_yyz_xxyzz, gs_x_yyz_xxzzz, gs_x_yyz_xyyyy, gs_x_yyz_xyyyz, gs_x_yyz_xyyzz, gs_x_yyz_xyzzz, gs_x_yyz_xzzzz, gs_x_yyz_yyyyy, gs_x_yyz_yyyyz, gs_x_yyz_yyyzz, gs_x_yyz_yyzzz, gs_x_yyz_yzzzz, gs_x_yyz_zzzzz, ts_yyz_xxxx, ts_yyz_xxxxx, ts_yyz_xxxxy, ts_yyz_xxxxz, ts_yyz_xxxy, ts_yyz_xxxyy, ts_yyz_xxxyz, ts_yyz_xxxz, ts_yyz_xxxzz, ts_yyz_xxyy, ts_yyz_xxyyy, ts_yyz_xxyyz, ts_yyz_xxyz, ts_yyz_xxyzz, ts_yyz_xxzz, ts_yyz_xxzzz, ts_yyz_xyyy, ts_yyz_xyyyy, ts_yyz_xyyyz, ts_yyz_xyyz, ts_yyz_xyyzz, ts_yyz_xyzz, ts_yyz_xyzzz, ts_yyz_xzzz, ts_yyz_xzzzz, ts_yyz_yyyy, ts_yyz_yyyyy, ts_yyz_yyyyz, ts_yyz_yyyz, ts_yyz_yyyzz, ts_yyz_yyzz, ts_yyz_yyzzz, ts_yyz_yzzz, ts_yyz_yzzzz, ts_yyz_zzzz, ts_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyz_xxxxx[i] = 10.0 * ts_yyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxxxy[i] = 8.0 * ts_yyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxxxz[i] = 8.0 * ts_yyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxxyy[i] = 6.0 * ts_yyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxxyz[i] = 6.0 * ts_yyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxxzz[i] = 6.0 * ts_yyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxyyy[i] = 4.0 * ts_yyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxyyz[i] = 4.0 * ts_yyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxyzz[i] = 4.0 * ts_yyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xxzzz[i] = 4.0 * ts_yyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xyyyy[i] = 2.0 * ts_yyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xyyyz[i] = 2.0 * ts_yyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xyyzz[i] = 2.0 * ts_yyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xyzzz[i] = 2.0 * ts_yyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_xzzzz[i] = 2.0 * ts_yyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yyyyy[i] = 2.0 * ts_yyz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yyz_yyyyz[i] = 2.0 * ts_yyz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yyyzz[i] = 2.0 * ts_yyz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yyzzz[i] = 2.0 * ts_yyz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yzzzz[i] = 2.0 * ts_yyz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yyz_zzzzz[i] = 2.0 * ts_yyz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 168-189 components of targeted buffer : FH

    auto gs_x_yzz_xxxxx = pbuffer.data(idx_g_fh + 168);

    auto gs_x_yzz_xxxxy = pbuffer.data(idx_g_fh + 169);

    auto gs_x_yzz_xxxxz = pbuffer.data(idx_g_fh + 170);

    auto gs_x_yzz_xxxyy = pbuffer.data(idx_g_fh + 171);

    auto gs_x_yzz_xxxyz = pbuffer.data(idx_g_fh + 172);

    auto gs_x_yzz_xxxzz = pbuffer.data(idx_g_fh + 173);

    auto gs_x_yzz_xxyyy = pbuffer.data(idx_g_fh + 174);

    auto gs_x_yzz_xxyyz = pbuffer.data(idx_g_fh + 175);

    auto gs_x_yzz_xxyzz = pbuffer.data(idx_g_fh + 176);

    auto gs_x_yzz_xxzzz = pbuffer.data(idx_g_fh + 177);

    auto gs_x_yzz_xyyyy = pbuffer.data(idx_g_fh + 178);

    auto gs_x_yzz_xyyyz = pbuffer.data(idx_g_fh + 179);

    auto gs_x_yzz_xyyzz = pbuffer.data(idx_g_fh + 180);

    auto gs_x_yzz_xyzzz = pbuffer.data(idx_g_fh + 181);

    auto gs_x_yzz_xzzzz = pbuffer.data(idx_g_fh + 182);

    auto gs_x_yzz_yyyyy = pbuffer.data(idx_g_fh + 183);

    auto gs_x_yzz_yyyyz = pbuffer.data(idx_g_fh + 184);

    auto gs_x_yzz_yyyzz = pbuffer.data(idx_g_fh + 185);

    auto gs_x_yzz_yyzzz = pbuffer.data(idx_g_fh + 186);

    auto gs_x_yzz_yzzzz = pbuffer.data(idx_g_fh + 187);

    auto gs_x_yzz_zzzzz = pbuffer.data(idx_g_fh + 188);

    #pragma omp simd aligned(gc_x, gs_x_yzz_xxxxx, gs_x_yzz_xxxxy, gs_x_yzz_xxxxz, gs_x_yzz_xxxyy, gs_x_yzz_xxxyz, gs_x_yzz_xxxzz, gs_x_yzz_xxyyy, gs_x_yzz_xxyyz, gs_x_yzz_xxyzz, gs_x_yzz_xxzzz, gs_x_yzz_xyyyy, gs_x_yzz_xyyyz, gs_x_yzz_xyyzz, gs_x_yzz_xyzzz, gs_x_yzz_xzzzz, gs_x_yzz_yyyyy, gs_x_yzz_yyyyz, gs_x_yzz_yyyzz, gs_x_yzz_yyzzz, gs_x_yzz_yzzzz, gs_x_yzz_zzzzz, ts_yzz_xxxx, ts_yzz_xxxxx, ts_yzz_xxxxy, ts_yzz_xxxxz, ts_yzz_xxxy, ts_yzz_xxxyy, ts_yzz_xxxyz, ts_yzz_xxxz, ts_yzz_xxxzz, ts_yzz_xxyy, ts_yzz_xxyyy, ts_yzz_xxyyz, ts_yzz_xxyz, ts_yzz_xxyzz, ts_yzz_xxzz, ts_yzz_xxzzz, ts_yzz_xyyy, ts_yzz_xyyyy, ts_yzz_xyyyz, ts_yzz_xyyz, ts_yzz_xyyzz, ts_yzz_xyzz, ts_yzz_xyzzz, ts_yzz_xzzz, ts_yzz_xzzzz, ts_yzz_yyyy, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyz, ts_yzz_yyyzz, ts_yzz_yyzz, ts_yzz_yyzzz, ts_yzz_yzzz, ts_yzz_yzzzz, ts_yzz_zzzz, ts_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzz_xxxxx[i] = 10.0 * ts_yzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxxxy[i] = 8.0 * ts_yzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxxxz[i] = 8.0 * ts_yzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxxyy[i] = 6.0 * ts_yzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxxyz[i] = 6.0 * ts_yzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxxzz[i] = 6.0 * ts_yzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxyyy[i] = 4.0 * ts_yzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxyyz[i] = 4.0 * ts_yzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxyzz[i] = 4.0 * ts_yzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xxzzz[i] = 4.0 * ts_yzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xyyyy[i] = 2.0 * ts_yzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xyyyz[i] = 2.0 * ts_yzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xyyzz[i] = 2.0 * ts_yzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xyzzz[i] = 2.0 * ts_yzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_xzzzz[i] = 2.0 * ts_yzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yyyyy[i] = 2.0 * ts_yzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_yzz_yyyyz[i] = 2.0 * ts_yzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yyyzz[i] = 2.0 * ts_yzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yyzzz[i] = 2.0 * ts_yzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yzzzz[i] = 2.0 * ts_yzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_yzz_zzzzz[i] = 2.0 * ts_yzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 189-210 components of targeted buffer : FH

    auto gs_x_zzz_xxxxx = pbuffer.data(idx_g_fh + 189);

    auto gs_x_zzz_xxxxy = pbuffer.data(idx_g_fh + 190);

    auto gs_x_zzz_xxxxz = pbuffer.data(idx_g_fh + 191);

    auto gs_x_zzz_xxxyy = pbuffer.data(idx_g_fh + 192);

    auto gs_x_zzz_xxxyz = pbuffer.data(idx_g_fh + 193);

    auto gs_x_zzz_xxxzz = pbuffer.data(idx_g_fh + 194);

    auto gs_x_zzz_xxyyy = pbuffer.data(idx_g_fh + 195);

    auto gs_x_zzz_xxyyz = pbuffer.data(idx_g_fh + 196);

    auto gs_x_zzz_xxyzz = pbuffer.data(idx_g_fh + 197);

    auto gs_x_zzz_xxzzz = pbuffer.data(idx_g_fh + 198);

    auto gs_x_zzz_xyyyy = pbuffer.data(idx_g_fh + 199);

    auto gs_x_zzz_xyyyz = pbuffer.data(idx_g_fh + 200);

    auto gs_x_zzz_xyyzz = pbuffer.data(idx_g_fh + 201);

    auto gs_x_zzz_xyzzz = pbuffer.data(idx_g_fh + 202);

    auto gs_x_zzz_xzzzz = pbuffer.data(idx_g_fh + 203);

    auto gs_x_zzz_yyyyy = pbuffer.data(idx_g_fh + 204);

    auto gs_x_zzz_yyyyz = pbuffer.data(idx_g_fh + 205);

    auto gs_x_zzz_yyyzz = pbuffer.data(idx_g_fh + 206);

    auto gs_x_zzz_yyzzz = pbuffer.data(idx_g_fh + 207);

    auto gs_x_zzz_yzzzz = pbuffer.data(idx_g_fh + 208);

    auto gs_x_zzz_zzzzz = pbuffer.data(idx_g_fh + 209);

    #pragma omp simd aligned(gc_x, gs_x_zzz_xxxxx, gs_x_zzz_xxxxy, gs_x_zzz_xxxxz, gs_x_zzz_xxxyy, gs_x_zzz_xxxyz, gs_x_zzz_xxxzz, gs_x_zzz_xxyyy, gs_x_zzz_xxyyz, gs_x_zzz_xxyzz, gs_x_zzz_xxzzz, gs_x_zzz_xyyyy, gs_x_zzz_xyyyz, gs_x_zzz_xyyzz, gs_x_zzz_xyzzz, gs_x_zzz_xzzzz, gs_x_zzz_yyyyy, gs_x_zzz_yyyyz, gs_x_zzz_yyyzz, gs_x_zzz_yyzzz, gs_x_zzz_yzzzz, gs_x_zzz_zzzzz, ts_zzz_xxxx, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxy, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxz, ts_zzz_xxxzz, ts_zzz_xxyy, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyz, ts_zzz_xxyzz, ts_zzz_xxzz, ts_zzz_xxzzz, ts_zzz_xyyy, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyz, ts_zzz_xyyzz, ts_zzz_xyzz, ts_zzz_xyzzz, ts_zzz_xzzz, ts_zzz_xzzzz, ts_zzz_yyyy, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyz, ts_zzz_yyyzz, ts_zzz_yyzz, ts_zzz_yyzzz, ts_zzz_yzzz, ts_zzz_yzzzz, ts_zzz_zzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzz_xxxxx[i] = 10.0 * ts_zzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxxxy[i] = 8.0 * ts_zzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxxxz[i] = 8.0 * ts_zzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxxyy[i] = 6.0 * ts_zzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxxyz[i] = 6.0 * ts_zzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxxzz[i] = 6.0 * ts_zzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxyyy[i] = 4.0 * ts_zzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxyyz[i] = 4.0 * ts_zzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxyzz[i] = 4.0 * ts_zzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xxzzz[i] = 4.0 * ts_zzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xyyyy[i] = 2.0 * ts_zzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xyyyz[i] = 2.0 * ts_zzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xyyzz[i] = 2.0 * ts_zzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xyzzz[i] = 2.0 * ts_zzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_xzzzz[i] = 2.0 * ts_zzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yyyyy[i] = 2.0 * ts_zzz_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_zzz_yyyyz[i] = 2.0 * ts_zzz_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yyyzz[i] = 2.0 * ts_zzz_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yyzzz[i] = 2.0 * ts_zzz_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yzzzz[i] = 2.0 * ts_zzz_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_zzz_zzzzz[i] = 2.0 * ts_zzz_zzzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 210-231 components of targeted buffer : FH

    auto gs_y_xxx_xxxxx = pbuffer.data(idx_g_fh + 210);

    auto gs_y_xxx_xxxxy = pbuffer.data(idx_g_fh + 211);

    auto gs_y_xxx_xxxxz = pbuffer.data(idx_g_fh + 212);

    auto gs_y_xxx_xxxyy = pbuffer.data(idx_g_fh + 213);

    auto gs_y_xxx_xxxyz = pbuffer.data(idx_g_fh + 214);

    auto gs_y_xxx_xxxzz = pbuffer.data(idx_g_fh + 215);

    auto gs_y_xxx_xxyyy = pbuffer.data(idx_g_fh + 216);

    auto gs_y_xxx_xxyyz = pbuffer.data(idx_g_fh + 217);

    auto gs_y_xxx_xxyzz = pbuffer.data(idx_g_fh + 218);

    auto gs_y_xxx_xxzzz = pbuffer.data(idx_g_fh + 219);

    auto gs_y_xxx_xyyyy = pbuffer.data(idx_g_fh + 220);

    auto gs_y_xxx_xyyyz = pbuffer.data(idx_g_fh + 221);

    auto gs_y_xxx_xyyzz = pbuffer.data(idx_g_fh + 222);

    auto gs_y_xxx_xyzzz = pbuffer.data(idx_g_fh + 223);

    auto gs_y_xxx_xzzzz = pbuffer.data(idx_g_fh + 224);

    auto gs_y_xxx_yyyyy = pbuffer.data(idx_g_fh + 225);

    auto gs_y_xxx_yyyyz = pbuffer.data(idx_g_fh + 226);

    auto gs_y_xxx_yyyzz = pbuffer.data(idx_g_fh + 227);

    auto gs_y_xxx_yyzzz = pbuffer.data(idx_g_fh + 228);

    auto gs_y_xxx_yzzzz = pbuffer.data(idx_g_fh + 229);

    auto gs_y_xxx_zzzzz = pbuffer.data(idx_g_fh + 230);

    #pragma omp simd aligned(gc_y, gs_y_xxx_xxxxx, gs_y_xxx_xxxxy, gs_y_xxx_xxxxz, gs_y_xxx_xxxyy, gs_y_xxx_xxxyz, gs_y_xxx_xxxzz, gs_y_xxx_xxyyy, gs_y_xxx_xxyyz, gs_y_xxx_xxyzz, gs_y_xxx_xxzzz, gs_y_xxx_xyyyy, gs_y_xxx_xyyyz, gs_y_xxx_xyyzz, gs_y_xxx_xyzzz, gs_y_xxx_xzzzz, gs_y_xxx_yyyyy, gs_y_xxx_yyyyz, gs_y_xxx_yyyzz, gs_y_xxx_yyzzz, gs_y_xxx_yzzzz, gs_y_xxx_zzzzz, ts_xxx_xxxx, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxy, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxz, ts_xxx_xxxzz, ts_xxx_xxyy, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyz, ts_xxx_xxyzz, ts_xxx_xxzz, ts_xxx_xxzzz, ts_xxx_xyyy, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyz, ts_xxx_xyyzz, ts_xxx_xyzz, ts_xxx_xyzzz, ts_xxx_xzzz, ts_xxx_xzzzz, ts_xxx_yyyy, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyz, ts_xxx_yyyzz, ts_xxx_yyzz, ts_xxx_yyzzz, ts_xxx_yzzz, ts_xxx_yzzzz, ts_xxx_zzzz, ts_xxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxx_xxxxx[i] = 2.0 * ts_xxx_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxxxy[i] = 2.0 * ts_xxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxxxz[i] = 2.0 * ts_xxx_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxxyy[i] = 4.0 * ts_xxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxxyz[i] = 2.0 * ts_xxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxxzz[i] = 2.0 * ts_xxx_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxyyy[i] = 6.0 * ts_xxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxyyz[i] = 4.0 * ts_xxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxyzz[i] = 2.0 * ts_xxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xxzzz[i] = 2.0 * ts_xxx_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xyyyy[i] = 8.0 * ts_xxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xyyyz[i] = 6.0 * ts_xxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xyyzz[i] = 4.0 * ts_xxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xyzzz[i] = 2.0 * ts_xxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_xzzzz[i] = 2.0 * ts_xxx_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yyyyy[i] = 10.0 * ts_xxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxx_yyyyz[i] = 8.0 * ts_xxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yyyzz[i] = 6.0 * ts_xxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yyzzz[i] = 4.0 * ts_xxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yzzzz[i] = 2.0 * ts_xxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxx_zzzzz[i] = 2.0 * ts_xxx_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 231-252 components of targeted buffer : FH

    auto gs_y_xxy_xxxxx = pbuffer.data(idx_g_fh + 231);

    auto gs_y_xxy_xxxxy = pbuffer.data(idx_g_fh + 232);

    auto gs_y_xxy_xxxxz = pbuffer.data(idx_g_fh + 233);

    auto gs_y_xxy_xxxyy = pbuffer.data(idx_g_fh + 234);

    auto gs_y_xxy_xxxyz = pbuffer.data(idx_g_fh + 235);

    auto gs_y_xxy_xxxzz = pbuffer.data(idx_g_fh + 236);

    auto gs_y_xxy_xxyyy = pbuffer.data(idx_g_fh + 237);

    auto gs_y_xxy_xxyyz = pbuffer.data(idx_g_fh + 238);

    auto gs_y_xxy_xxyzz = pbuffer.data(idx_g_fh + 239);

    auto gs_y_xxy_xxzzz = pbuffer.data(idx_g_fh + 240);

    auto gs_y_xxy_xyyyy = pbuffer.data(idx_g_fh + 241);

    auto gs_y_xxy_xyyyz = pbuffer.data(idx_g_fh + 242);

    auto gs_y_xxy_xyyzz = pbuffer.data(idx_g_fh + 243);

    auto gs_y_xxy_xyzzz = pbuffer.data(idx_g_fh + 244);

    auto gs_y_xxy_xzzzz = pbuffer.data(idx_g_fh + 245);

    auto gs_y_xxy_yyyyy = pbuffer.data(idx_g_fh + 246);

    auto gs_y_xxy_yyyyz = pbuffer.data(idx_g_fh + 247);

    auto gs_y_xxy_yyyzz = pbuffer.data(idx_g_fh + 248);

    auto gs_y_xxy_yyzzz = pbuffer.data(idx_g_fh + 249);

    auto gs_y_xxy_yzzzz = pbuffer.data(idx_g_fh + 250);

    auto gs_y_xxy_zzzzz = pbuffer.data(idx_g_fh + 251);

    #pragma omp simd aligned(gc_y, gs_y_xxy_xxxxx, gs_y_xxy_xxxxy, gs_y_xxy_xxxxz, gs_y_xxy_xxxyy, gs_y_xxy_xxxyz, gs_y_xxy_xxxzz, gs_y_xxy_xxyyy, gs_y_xxy_xxyyz, gs_y_xxy_xxyzz, gs_y_xxy_xxzzz, gs_y_xxy_xyyyy, gs_y_xxy_xyyyz, gs_y_xxy_xyyzz, gs_y_xxy_xyzzz, gs_y_xxy_xzzzz, gs_y_xxy_yyyyy, gs_y_xxy_yyyyz, gs_y_xxy_yyyzz, gs_y_xxy_yyzzz, gs_y_xxy_yzzzz, gs_y_xxy_zzzzz, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxxz, ts_xx_xxxyy, ts_xx_xxxyz, ts_xx_xxxzz, ts_xx_xxyyy, ts_xx_xxyyz, ts_xx_xxyzz, ts_xx_xxzzz, ts_xx_xyyyy, ts_xx_xyyyz, ts_xx_xyyzz, ts_xx_xyzzz, ts_xx_xzzzz, ts_xx_yyyyy, ts_xx_yyyyz, ts_xx_yyyzz, ts_xx_yyzzz, ts_xx_yzzzz, ts_xx_zzzzz, ts_xxy_xxxx, ts_xxy_xxxxx, ts_xxy_xxxxy, ts_xxy_xxxxz, ts_xxy_xxxy, ts_xxy_xxxyy, ts_xxy_xxxyz, ts_xxy_xxxz, ts_xxy_xxxzz, ts_xxy_xxyy, ts_xxy_xxyyy, ts_xxy_xxyyz, ts_xxy_xxyz, ts_xxy_xxyzz, ts_xxy_xxzz, ts_xxy_xxzzz, ts_xxy_xyyy, ts_xxy_xyyyy, ts_xxy_xyyyz, ts_xxy_xyyz, ts_xxy_xyyzz, ts_xxy_xyzz, ts_xxy_xyzzz, ts_xxy_xzzz, ts_xxy_xzzzz, ts_xxy_yyyy, ts_xxy_yyyyy, ts_xxy_yyyyz, ts_xxy_yyyz, ts_xxy_yyyzz, ts_xxy_yyzz, ts_xxy_yyzzz, ts_xxy_yzzz, ts_xxy_yzzzz, ts_xxy_zzzz, ts_xxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxy_xxxxx[i] = 2.0 * ts_xx_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxxxy[i] = 2.0 * ts_xx_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxxxz[i] = 2.0 * ts_xx_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxxyy[i] = 2.0 * ts_xx_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxxyz[i] = 2.0 * ts_xx_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxxzz[i] = 2.0 * ts_xx_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxyyy[i] = 2.0 * ts_xx_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxyyz[i] = 2.0 * ts_xx_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxyzz[i] = 2.0 * ts_xx_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xxzzz[i] = 2.0 * ts_xx_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xyyyy[i] = 2.0 * ts_xx_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xyyyz[i] = 2.0 * ts_xx_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xyyzz[i] = 2.0 * ts_xx_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xyzzz[i] = 2.0 * ts_xx_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_xzzzz[i] = 2.0 * ts_xx_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yyyyy[i] = 2.0 * ts_xx_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxy_yyyyz[i] = 2.0 * ts_xx_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yyyzz[i] = 2.0 * ts_xx_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yyzzz[i] = 2.0 * ts_xx_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yzzzz[i] = 2.0 * ts_xx_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxy_zzzzz[i] = 2.0 * ts_xx_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 252-273 components of targeted buffer : FH

    auto gs_y_xxz_xxxxx = pbuffer.data(idx_g_fh + 252);

    auto gs_y_xxz_xxxxy = pbuffer.data(idx_g_fh + 253);

    auto gs_y_xxz_xxxxz = pbuffer.data(idx_g_fh + 254);

    auto gs_y_xxz_xxxyy = pbuffer.data(idx_g_fh + 255);

    auto gs_y_xxz_xxxyz = pbuffer.data(idx_g_fh + 256);

    auto gs_y_xxz_xxxzz = pbuffer.data(idx_g_fh + 257);

    auto gs_y_xxz_xxyyy = pbuffer.data(idx_g_fh + 258);

    auto gs_y_xxz_xxyyz = pbuffer.data(idx_g_fh + 259);

    auto gs_y_xxz_xxyzz = pbuffer.data(idx_g_fh + 260);

    auto gs_y_xxz_xxzzz = pbuffer.data(idx_g_fh + 261);

    auto gs_y_xxz_xyyyy = pbuffer.data(idx_g_fh + 262);

    auto gs_y_xxz_xyyyz = pbuffer.data(idx_g_fh + 263);

    auto gs_y_xxz_xyyzz = pbuffer.data(idx_g_fh + 264);

    auto gs_y_xxz_xyzzz = pbuffer.data(idx_g_fh + 265);

    auto gs_y_xxz_xzzzz = pbuffer.data(idx_g_fh + 266);

    auto gs_y_xxz_yyyyy = pbuffer.data(idx_g_fh + 267);

    auto gs_y_xxz_yyyyz = pbuffer.data(idx_g_fh + 268);

    auto gs_y_xxz_yyyzz = pbuffer.data(idx_g_fh + 269);

    auto gs_y_xxz_yyzzz = pbuffer.data(idx_g_fh + 270);

    auto gs_y_xxz_yzzzz = pbuffer.data(idx_g_fh + 271);

    auto gs_y_xxz_zzzzz = pbuffer.data(idx_g_fh + 272);

    #pragma omp simd aligned(gc_y, gs_y_xxz_xxxxx, gs_y_xxz_xxxxy, gs_y_xxz_xxxxz, gs_y_xxz_xxxyy, gs_y_xxz_xxxyz, gs_y_xxz_xxxzz, gs_y_xxz_xxyyy, gs_y_xxz_xxyyz, gs_y_xxz_xxyzz, gs_y_xxz_xxzzz, gs_y_xxz_xyyyy, gs_y_xxz_xyyyz, gs_y_xxz_xyyzz, gs_y_xxz_xyzzz, gs_y_xxz_xzzzz, gs_y_xxz_yyyyy, gs_y_xxz_yyyyz, gs_y_xxz_yyyzz, gs_y_xxz_yyzzz, gs_y_xxz_yzzzz, gs_y_xxz_zzzzz, ts_xxz_xxxx, ts_xxz_xxxxx, ts_xxz_xxxxy, ts_xxz_xxxxz, ts_xxz_xxxy, ts_xxz_xxxyy, ts_xxz_xxxyz, ts_xxz_xxxz, ts_xxz_xxxzz, ts_xxz_xxyy, ts_xxz_xxyyy, ts_xxz_xxyyz, ts_xxz_xxyz, ts_xxz_xxyzz, ts_xxz_xxzz, ts_xxz_xxzzz, ts_xxz_xyyy, ts_xxz_xyyyy, ts_xxz_xyyyz, ts_xxz_xyyz, ts_xxz_xyyzz, ts_xxz_xyzz, ts_xxz_xyzzz, ts_xxz_xzzz, ts_xxz_xzzzz, ts_xxz_yyyy, ts_xxz_yyyyy, ts_xxz_yyyyz, ts_xxz_yyyz, ts_xxz_yyyzz, ts_xxz_yyzz, ts_xxz_yyzzz, ts_xxz_yzzz, ts_xxz_yzzzz, ts_xxz_zzzz, ts_xxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxz_xxxxx[i] = 2.0 * ts_xxz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxxxy[i] = 2.0 * ts_xxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxxxz[i] = 2.0 * ts_xxz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxxyy[i] = 4.0 * ts_xxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxxyz[i] = 2.0 * ts_xxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxxzz[i] = 2.0 * ts_xxz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxyyy[i] = 6.0 * ts_xxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxyyz[i] = 4.0 * ts_xxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxyzz[i] = 2.0 * ts_xxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xxzzz[i] = 2.0 * ts_xxz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xyyyy[i] = 8.0 * ts_xxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xyyyz[i] = 6.0 * ts_xxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xyyzz[i] = 4.0 * ts_xxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xyzzz[i] = 2.0 * ts_xxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_xzzzz[i] = 2.0 * ts_xxz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yyyyy[i] = 10.0 * ts_xxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xxz_yyyyz[i] = 8.0 * ts_xxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yyyzz[i] = 6.0 * ts_xxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yyzzz[i] = 4.0 * ts_xxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yzzzz[i] = 2.0 * ts_xxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xxz_zzzzz[i] = 2.0 * ts_xxz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 273-294 components of targeted buffer : FH

    auto gs_y_xyy_xxxxx = pbuffer.data(idx_g_fh + 273);

    auto gs_y_xyy_xxxxy = pbuffer.data(idx_g_fh + 274);

    auto gs_y_xyy_xxxxz = pbuffer.data(idx_g_fh + 275);

    auto gs_y_xyy_xxxyy = pbuffer.data(idx_g_fh + 276);

    auto gs_y_xyy_xxxyz = pbuffer.data(idx_g_fh + 277);

    auto gs_y_xyy_xxxzz = pbuffer.data(idx_g_fh + 278);

    auto gs_y_xyy_xxyyy = pbuffer.data(idx_g_fh + 279);

    auto gs_y_xyy_xxyyz = pbuffer.data(idx_g_fh + 280);

    auto gs_y_xyy_xxyzz = pbuffer.data(idx_g_fh + 281);

    auto gs_y_xyy_xxzzz = pbuffer.data(idx_g_fh + 282);

    auto gs_y_xyy_xyyyy = pbuffer.data(idx_g_fh + 283);

    auto gs_y_xyy_xyyyz = pbuffer.data(idx_g_fh + 284);

    auto gs_y_xyy_xyyzz = pbuffer.data(idx_g_fh + 285);

    auto gs_y_xyy_xyzzz = pbuffer.data(idx_g_fh + 286);

    auto gs_y_xyy_xzzzz = pbuffer.data(idx_g_fh + 287);

    auto gs_y_xyy_yyyyy = pbuffer.data(idx_g_fh + 288);

    auto gs_y_xyy_yyyyz = pbuffer.data(idx_g_fh + 289);

    auto gs_y_xyy_yyyzz = pbuffer.data(idx_g_fh + 290);

    auto gs_y_xyy_yyzzz = pbuffer.data(idx_g_fh + 291);

    auto gs_y_xyy_yzzzz = pbuffer.data(idx_g_fh + 292);

    auto gs_y_xyy_zzzzz = pbuffer.data(idx_g_fh + 293);

    #pragma omp simd aligned(gc_y, gs_y_xyy_xxxxx, gs_y_xyy_xxxxy, gs_y_xyy_xxxxz, gs_y_xyy_xxxyy, gs_y_xyy_xxxyz, gs_y_xyy_xxxzz, gs_y_xyy_xxyyy, gs_y_xyy_xxyyz, gs_y_xyy_xxyzz, gs_y_xyy_xxzzz, gs_y_xyy_xyyyy, gs_y_xyy_xyyyz, gs_y_xyy_xyyzz, gs_y_xyy_xyzzz, gs_y_xyy_xzzzz, gs_y_xyy_yyyyy, gs_y_xyy_yyyyz, gs_y_xyy_yyyzz, gs_y_xyy_yyzzz, gs_y_xyy_yzzzz, gs_y_xyy_zzzzz, ts_xy_xxxxx, ts_xy_xxxxy, ts_xy_xxxxz, ts_xy_xxxyy, ts_xy_xxxyz, ts_xy_xxxzz, ts_xy_xxyyy, ts_xy_xxyyz, ts_xy_xxyzz, ts_xy_xxzzz, ts_xy_xyyyy, ts_xy_xyyyz, ts_xy_xyyzz, ts_xy_xyzzz, ts_xy_xzzzz, ts_xy_yyyyy, ts_xy_yyyyz, ts_xy_yyyzz, ts_xy_yyzzz, ts_xy_yzzzz, ts_xy_zzzzz, ts_xyy_xxxx, ts_xyy_xxxxx, ts_xyy_xxxxy, ts_xyy_xxxxz, ts_xyy_xxxy, ts_xyy_xxxyy, ts_xyy_xxxyz, ts_xyy_xxxz, ts_xyy_xxxzz, ts_xyy_xxyy, ts_xyy_xxyyy, ts_xyy_xxyyz, ts_xyy_xxyz, ts_xyy_xxyzz, ts_xyy_xxzz, ts_xyy_xxzzz, ts_xyy_xyyy, ts_xyy_xyyyy, ts_xyy_xyyyz, ts_xyy_xyyz, ts_xyy_xyyzz, ts_xyy_xyzz, ts_xyy_xyzzz, ts_xyy_xzzz, ts_xyy_xzzzz, ts_xyy_yyyy, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyz, ts_xyy_yyyzz, ts_xyy_yyzz, ts_xyy_yyzzz, ts_xyy_yzzz, ts_xyy_yzzzz, ts_xyy_zzzz, ts_xyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyy_xxxxx[i] = 4.0 * ts_xy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxxxy[i] = 4.0 * ts_xy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxxxz[i] = 4.0 * ts_xy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxxyy[i] = 4.0 * ts_xy_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxxyz[i] = 4.0 * ts_xy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxxzz[i] = 4.0 * ts_xy_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxyyy[i] = 4.0 * ts_xy_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxyyz[i] = 4.0 * ts_xy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxyzz[i] = 4.0 * ts_xy_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xxzzz[i] = 4.0 * ts_xy_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xyyyy[i] = 4.0 * ts_xy_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xyyyz[i] = 4.0 * ts_xy_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xyyzz[i] = 4.0 * ts_xy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xyzzz[i] = 4.0 * ts_xy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_xzzzz[i] = 4.0 * ts_xy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yyyyy[i] = 4.0 * ts_xy_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyy_yyyyz[i] = 4.0 * ts_xy_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yyyzz[i] = 4.0 * ts_xy_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yyzzz[i] = 4.0 * ts_xy_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yzzzz[i] = 4.0 * ts_xy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyy_zzzzz[i] = 4.0 * ts_xy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 294-315 components of targeted buffer : FH

    auto gs_y_xyz_xxxxx = pbuffer.data(idx_g_fh + 294);

    auto gs_y_xyz_xxxxy = pbuffer.data(idx_g_fh + 295);

    auto gs_y_xyz_xxxxz = pbuffer.data(idx_g_fh + 296);

    auto gs_y_xyz_xxxyy = pbuffer.data(idx_g_fh + 297);

    auto gs_y_xyz_xxxyz = pbuffer.data(idx_g_fh + 298);

    auto gs_y_xyz_xxxzz = pbuffer.data(idx_g_fh + 299);

    auto gs_y_xyz_xxyyy = pbuffer.data(idx_g_fh + 300);

    auto gs_y_xyz_xxyyz = pbuffer.data(idx_g_fh + 301);

    auto gs_y_xyz_xxyzz = pbuffer.data(idx_g_fh + 302);

    auto gs_y_xyz_xxzzz = pbuffer.data(idx_g_fh + 303);

    auto gs_y_xyz_xyyyy = pbuffer.data(idx_g_fh + 304);

    auto gs_y_xyz_xyyyz = pbuffer.data(idx_g_fh + 305);

    auto gs_y_xyz_xyyzz = pbuffer.data(idx_g_fh + 306);

    auto gs_y_xyz_xyzzz = pbuffer.data(idx_g_fh + 307);

    auto gs_y_xyz_xzzzz = pbuffer.data(idx_g_fh + 308);

    auto gs_y_xyz_yyyyy = pbuffer.data(idx_g_fh + 309);

    auto gs_y_xyz_yyyyz = pbuffer.data(idx_g_fh + 310);

    auto gs_y_xyz_yyyzz = pbuffer.data(idx_g_fh + 311);

    auto gs_y_xyz_yyzzz = pbuffer.data(idx_g_fh + 312);

    auto gs_y_xyz_yzzzz = pbuffer.data(idx_g_fh + 313);

    auto gs_y_xyz_zzzzz = pbuffer.data(idx_g_fh + 314);

    #pragma omp simd aligned(gc_y, gs_y_xyz_xxxxx, gs_y_xyz_xxxxy, gs_y_xyz_xxxxz, gs_y_xyz_xxxyy, gs_y_xyz_xxxyz, gs_y_xyz_xxxzz, gs_y_xyz_xxyyy, gs_y_xyz_xxyyz, gs_y_xyz_xxyzz, gs_y_xyz_xxzzz, gs_y_xyz_xyyyy, gs_y_xyz_xyyyz, gs_y_xyz_xyyzz, gs_y_xyz_xyzzz, gs_y_xyz_xzzzz, gs_y_xyz_yyyyy, gs_y_xyz_yyyyz, gs_y_xyz_yyyzz, gs_y_xyz_yyzzz, gs_y_xyz_yzzzz, gs_y_xyz_zzzzz, ts_xyz_xxxx, ts_xyz_xxxxx, ts_xyz_xxxxy, ts_xyz_xxxxz, ts_xyz_xxxy, ts_xyz_xxxyy, ts_xyz_xxxyz, ts_xyz_xxxz, ts_xyz_xxxzz, ts_xyz_xxyy, ts_xyz_xxyyy, ts_xyz_xxyyz, ts_xyz_xxyz, ts_xyz_xxyzz, ts_xyz_xxzz, ts_xyz_xxzzz, ts_xyz_xyyy, ts_xyz_xyyyy, ts_xyz_xyyyz, ts_xyz_xyyz, ts_xyz_xyyzz, ts_xyz_xyzz, ts_xyz_xyzzz, ts_xyz_xzzz, ts_xyz_xzzzz, ts_xyz_yyyy, ts_xyz_yyyyy, ts_xyz_yyyyz, ts_xyz_yyyz, ts_xyz_yyyzz, ts_xyz_yyzz, ts_xyz_yyzzz, ts_xyz_yzzz, ts_xyz_yzzzz, ts_xyz_zzzz, ts_xyz_zzzzz, ts_xz_xxxxx, ts_xz_xxxxy, ts_xz_xxxxz, ts_xz_xxxyy, ts_xz_xxxyz, ts_xz_xxxzz, ts_xz_xxyyy, ts_xz_xxyyz, ts_xz_xxyzz, ts_xz_xxzzz, ts_xz_xyyyy, ts_xz_xyyyz, ts_xz_xyyzz, ts_xz_xyzzz, ts_xz_xzzzz, ts_xz_yyyyy, ts_xz_yyyyz, ts_xz_yyyzz, ts_xz_yyzzz, ts_xz_yzzzz, ts_xz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyz_xxxxx[i] = 2.0 * ts_xz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxxxy[i] = 2.0 * ts_xz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxxxz[i] = 2.0 * ts_xz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxxyy[i] = 2.0 * ts_xz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxxyz[i] = 2.0 * ts_xz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxxzz[i] = 2.0 * ts_xz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxyyy[i] = 2.0 * ts_xz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxyyz[i] = 2.0 * ts_xz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxyzz[i] = 2.0 * ts_xz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xxzzz[i] = 2.0 * ts_xz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xyyyy[i] = 2.0 * ts_xz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xyyyz[i] = 2.0 * ts_xz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xyyzz[i] = 2.0 * ts_xz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xyzzz[i] = 2.0 * ts_xz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_xzzzz[i] = 2.0 * ts_xz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yyyyy[i] = 2.0 * ts_xz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_xyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xyz_yyyyz[i] = 2.0 * ts_xz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_xyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yyyzz[i] = 2.0 * ts_xz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yyzzz[i] = 2.0 * ts_xz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yzzzz[i] = 2.0 * ts_xz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xyz_zzzzz[i] = 2.0 * ts_xz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 315-336 components of targeted buffer : FH

    auto gs_y_xzz_xxxxx = pbuffer.data(idx_g_fh + 315);

    auto gs_y_xzz_xxxxy = pbuffer.data(idx_g_fh + 316);

    auto gs_y_xzz_xxxxz = pbuffer.data(idx_g_fh + 317);

    auto gs_y_xzz_xxxyy = pbuffer.data(idx_g_fh + 318);

    auto gs_y_xzz_xxxyz = pbuffer.data(idx_g_fh + 319);

    auto gs_y_xzz_xxxzz = pbuffer.data(idx_g_fh + 320);

    auto gs_y_xzz_xxyyy = pbuffer.data(idx_g_fh + 321);

    auto gs_y_xzz_xxyyz = pbuffer.data(idx_g_fh + 322);

    auto gs_y_xzz_xxyzz = pbuffer.data(idx_g_fh + 323);

    auto gs_y_xzz_xxzzz = pbuffer.data(idx_g_fh + 324);

    auto gs_y_xzz_xyyyy = pbuffer.data(idx_g_fh + 325);

    auto gs_y_xzz_xyyyz = pbuffer.data(idx_g_fh + 326);

    auto gs_y_xzz_xyyzz = pbuffer.data(idx_g_fh + 327);

    auto gs_y_xzz_xyzzz = pbuffer.data(idx_g_fh + 328);

    auto gs_y_xzz_xzzzz = pbuffer.data(idx_g_fh + 329);

    auto gs_y_xzz_yyyyy = pbuffer.data(idx_g_fh + 330);

    auto gs_y_xzz_yyyyz = pbuffer.data(idx_g_fh + 331);

    auto gs_y_xzz_yyyzz = pbuffer.data(idx_g_fh + 332);

    auto gs_y_xzz_yyzzz = pbuffer.data(idx_g_fh + 333);

    auto gs_y_xzz_yzzzz = pbuffer.data(idx_g_fh + 334);

    auto gs_y_xzz_zzzzz = pbuffer.data(idx_g_fh + 335);

    #pragma omp simd aligned(gc_y, gs_y_xzz_xxxxx, gs_y_xzz_xxxxy, gs_y_xzz_xxxxz, gs_y_xzz_xxxyy, gs_y_xzz_xxxyz, gs_y_xzz_xxxzz, gs_y_xzz_xxyyy, gs_y_xzz_xxyyz, gs_y_xzz_xxyzz, gs_y_xzz_xxzzz, gs_y_xzz_xyyyy, gs_y_xzz_xyyyz, gs_y_xzz_xyyzz, gs_y_xzz_xyzzz, gs_y_xzz_xzzzz, gs_y_xzz_yyyyy, gs_y_xzz_yyyyz, gs_y_xzz_yyyzz, gs_y_xzz_yyzzz, gs_y_xzz_yzzzz, gs_y_xzz_zzzzz, ts_xzz_xxxx, ts_xzz_xxxxx, ts_xzz_xxxxy, ts_xzz_xxxxz, ts_xzz_xxxy, ts_xzz_xxxyy, ts_xzz_xxxyz, ts_xzz_xxxz, ts_xzz_xxxzz, ts_xzz_xxyy, ts_xzz_xxyyy, ts_xzz_xxyyz, ts_xzz_xxyz, ts_xzz_xxyzz, ts_xzz_xxzz, ts_xzz_xxzzz, ts_xzz_xyyy, ts_xzz_xyyyy, ts_xzz_xyyyz, ts_xzz_xyyz, ts_xzz_xyyzz, ts_xzz_xyzz, ts_xzz_xyzzz, ts_xzz_xzzz, ts_xzz_xzzzz, ts_xzz_yyyy, ts_xzz_yyyyy, ts_xzz_yyyyz, ts_xzz_yyyz, ts_xzz_yyyzz, ts_xzz_yyzz, ts_xzz_yyzzz, ts_xzz_yzzz, ts_xzz_yzzzz, ts_xzz_zzzz, ts_xzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzz_xxxxx[i] = 2.0 * ts_xzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxxxy[i] = 2.0 * ts_xzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxxxz[i] = 2.0 * ts_xzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxxyy[i] = 4.0 * ts_xzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxxyz[i] = 2.0 * ts_xzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxxzz[i] = 2.0 * ts_xzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxyyy[i] = 6.0 * ts_xzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxyyz[i] = 4.0 * ts_xzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxyzz[i] = 2.0 * ts_xzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xxzzz[i] = 2.0 * ts_xzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xyyyy[i] = 8.0 * ts_xzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xyyyz[i] = 6.0 * ts_xzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xyyzz[i] = 4.0 * ts_xzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xyzzz[i] = 2.0 * ts_xzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_xzzzz[i] = 2.0 * ts_xzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yyyyy[i] = 10.0 * ts_xzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_xzz_yyyyz[i] = 8.0 * ts_xzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yyyzz[i] = 6.0 * ts_xzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yyzzz[i] = 4.0 * ts_xzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yzzzz[i] = 2.0 * ts_xzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_xzz_zzzzz[i] = 2.0 * ts_xzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 336-357 components of targeted buffer : FH

    auto gs_y_yyy_xxxxx = pbuffer.data(idx_g_fh + 336);

    auto gs_y_yyy_xxxxy = pbuffer.data(idx_g_fh + 337);

    auto gs_y_yyy_xxxxz = pbuffer.data(idx_g_fh + 338);

    auto gs_y_yyy_xxxyy = pbuffer.data(idx_g_fh + 339);

    auto gs_y_yyy_xxxyz = pbuffer.data(idx_g_fh + 340);

    auto gs_y_yyy_xxxzz = pbuffer.data(idx_g_fh + 341);

    auto gs_y_yyy_xxyyy = pbuffer.data(idx_g_fh + 342);

    auto gs_y_yyy_xxyyz = pbuffer.data(idx_g_fh + 343);

    auto gs_y_yyy_xxyzz = pbuffer.data(idx_g_fh + 344);

    auto gs_y_yyy_xxzzz = pbuffer.data(idx_g_fh + 345);

    auto gs_y_yyy_xyyyy = pbuffer.data(idx_g_fh + 346);

    auto gs_y_yyy_xyyyz = pbuffer.data(idx_g_fh + 347);

    auto gs_y_yyy_xyyzz = pbuffer.data(idx_g_fh + 348);

    auto gs_y_yyy_xyzzz = pbuffer.data(idx_g_fh + 349);

    auto gs_y_yyy_xzzzz = pbuffer.data(idx_g_fh + 350);

    auto gs_y_yyy_yyyyy = pbuffer.data(idx_g_fh + 351);

    auto gs_y_yyy_yyyyz = pbuffer.data(idx_g_fh + 352);

    auto gs_y_yyy_yyyzz = pbuffer.data(idx_g_fh + 353);

    auto gs_y_yyy_yyzzz = pbuffer.data(idx_g_fh + 354);

    auto gs_y_yyy_yzzzz = pbuffer.data(idx_g_fh + 355);

    auto gs_y_yyy_zzzzz = pbuffer.data(idx_g_fh + 356);

    #pragma omp simd aligned(gc_y, gs_y_yyy_xxxxx, gs_y_yyy_xxxxy, gs_y_yyy_xxxxz, gs_y_yyy_xxxyy, gs_y_yyy_xxxyz, gs_y_yyy_xxxzz, gs_y_yyy_xxyyy, gs_y_yyy_xxyyz, gs_y_yyy_xxyzz, gs_y_yyy_xxzzz, gs_y_yyy_xyyyy, gs_y_yyy_xyyyz, gs_y_yyy_xyyzz, gs_y_yyy_xyzzz, gs_y_yyy_xzzzz, gs_y_yyy_yyyyy, gs_y_yyy_yyyyz, gs_y_yyy_yyyzz, gs_y_yyy_yyzzz, gs_y_yyy_yzzzz, gs_y_yyy_zzzzz, ts_yy_xxxxx, ts_yy_xxxxy, ts_yy_xxxxz, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxxzz, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyzz, ts_yy_xxzzz, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyzz, ts_yy_xyzzz, ts_yy_xzzzz, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyzz, ts_yy_yyzzz, ts_yy_yzzzz, ts_yy_zzzzz, ts_yyy_xxxx, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxy, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxz, ts_yyy_xxxzz, ts_yyy_xxyy, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyz, ts_yyy_xxyzz, ts_yyy_xxzz, ts_yyy_xxzzz, ts_yyy_xyyy, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyz, ts_yyy_xyyzz, ts_yyy_xyzz, ts_yyy_xyzzz, ts_yyy_xzzz, ts_yyy_xzzzz, ts_yyy_yyyy, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyz, ts_yyy_yyyzz, ts_yyy_yyzz, ts_yyy_yyzzz, ts_yyy_yzzz, ts_yyy_yzzzz, ts_yyy_zzzz, ts_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyy_xxxxx[i] = 6.0 * ts_yy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxxxy[i] = 6.0 * ts_yy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxxxz[i] = 6.0 * ts_yy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxxyy[i] = 6.0 * ts_yy_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxxyz[i] = 6.0 * ts_yy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxxzz[i] = 6.0 * ts_yy_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxyyy[i] = 6.0 * ts_yy_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxyyz[i] = 6.0 * ts_yy_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxyzz[i] = 6.0 * ts_yy_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xxzzz[i] = 6.0 * ts_yy_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xyyyy[i] = 6.0 * ts_yy_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xyyyz[i] = 6.0 * ts_yy_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xyyzz[i] = 6.0 * ts_yy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xyzzz[i] = 6.0 * ts_yy_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_xzzzz[i] = 6.0 * ts_yy_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yyyyy[i] = 6.0 * ts_yy_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyy_yyyyz[i] = 6.0 * ts_yy_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yyyzz[i] = 6.0 * ts_yy_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yyzzz[i] = 6.0 * ts_yy_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yzzzz[i] = 6.0 * ts_yy_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyy_zzzzz[i] = 6.0 * ts_yy_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 357-378 components of targeted buffer : FH

    auto gs_y_yyz_xxxxx = pbuffer.data(idx_g_fh + 357);

    auto gs_y_yyz_xxxxy = pbuffer.data(idx_g_fh + 358);

    auto gs_y_yyz_xxxxz = pbuffer.data(idx_g_fh + 359);

    auto gs_y_yyz_xxxyy = pbuffer.data(idx_g_fh + 360);

    auto gs_y_yyz_xxxyz = pbuffer.data(idx_g_fh + 361);

    auto gs_y_yyz_xxxzz = pbuffer.data(idx_g_fh + 362);

    auto gs_y_yyz_xxyyy = pbuffer.data(idx_g_fh + 363);

    auto gs_y_yyz_xxyyz = pbuffer.data(idx_g_fh + 364);

    auto gs_y_yyz_xxyzz = pbuffer.data(idx_g_fh + 365);

    auto gs_y_yyz_xxzzz = pbuffer.data(idx_g_fh + 366);

    auto gs_y_yyz_xyyyy = pbuffer.data(idx_g_fh + 367);

    auto gs_y_yyz_xyyyz = pbuffer.data(idx_g_fh + 368);

    auto gs_y_yyz_xyyzz = pbuffer.data(idx_g_fh + 369);

    auto gs_y_yyz_xyzzz = pbuffer.data(idx_g_fh + 370);

    auto gs_y_yyz_xzzzz = pbuffer.data(idx_g_fh + 371);

    auto gs_y_yyz_yyyyy = pbuffer.data(idx_g_fh + 372);

    auto gs_y_yyz_yyyyz = pbuffer.data(idx_g_fh + 373);

    auto gs_y_yyz_yyyzz = pbuffer.data(idx_g_fh + 374);

    auto gs_y_yyz_yyzzz = pbuffer.data(idx_g_fh + 375);

    auto gs_y_yyz_yzzzz = pbuffer.data(idx_g_fh + 376);

    auto gs_y_yyz_zzzzz = pbuffer.data(idx_g_fh + 377);

    #pragma omp simd aligned(gc_y, gs_y_yyz_xxxxx, gs_y_yyz_xxxxy, gs_y_yyz_xxxxz, gs_y_yyz_xxxyy, gs_y_yyz_xxxyz, gs_y_yyz_xxxzz, gs_y_yyz_xxyyy, gs_y_yyz_xxyyz, gs_y_yyz_xxyzz, gs_y_yyz_xxzzz, gs_y_yyz_xyyyy, gs_y_yyz_xyyyz, gs_y_yyz_xyyzz, gs_y_yyz_xyzzz, gs_y_yyz_xzzzz, gs_y_yyz_yyyyy, gs_y_yyz_yyyyz, gs_y_yyz_yyyzz, gs_y_yyz_yyzzz, gs_y_yyz_yzzzz, gs_y_yyz_zzzzz, ts_yyz_xxxx, ts_yyz_xxxxx, ts_yyz_xxxxy, ts_yyz_xxxxz, ts_yyz_xxxy, ts_yyz_xxxyy, ts_yyz_xxxyz, ts_yyz_xxxz, ts_yyz_xxxzz, ts_yyz_xxyy, ts_yyz_xxyyy, ts_yyz_xxyyz, ts_yyz_xxyz, ts_yyz_xxyzz, ts_yyz_xxzz, ts_yyz_xxzzz, ts_yyz_xyyy, ts_yyz_xyyyy, ts_yyz_xyyyz, ts_yyz_xyyz, ts_yyz_xyyzz, ts_yyz_xyzz, ts_yyz_xyzzz, ts_yyz_xzzz, ts_yyz_xzzzz, ts_yyz_yyyy, ts_yyz_yyyyy, ts_yyz_yyyyz, ts_yyz_yyyz, ts_yyz_yyyzz, ts_yyz_yyzz, ts_yyz_yyzzz, ts_yyz_yzzz, ts_yyz_yzzzz, ts_yyz_zzzz, ts_yyz_zzzzz, ts_yz_xxxxx, ts_yz_xxxxy, ts_yz_xxxxz, ts_yz_xxxyy, ts_yz_xxxyz, ts_yz_xxxzz, ts_yz_xxyyy, ts_yz_xxyyz, ts_yz_xxyzz, ts_yz_xxzzz, ts_yz_xyyyy, ts_yz_xyyyz, ts_yz_xyyzz, ts_yz_xyzzz, ts_yz_xzzzz, ts_yz_yyyyy, ts_yz_yyyyz, ts_yz_yyyzz, ts_yz_yyzzz, ts_yz_yzzzz, ts_yz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyz_xxxxx[i] = 4.0 * ts_yz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxxxy[i] = 4.0 * ts_yz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxxxz[i] = 4.0 * ts_yz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxxyy[i] = 4.0 * ts_yz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxxyz[i] = 4.0 * ts_yz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxxzz[i] = 4.0 * ts_yz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxyyy[i] = 4.0 * ts_yz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxyyz[i] = 4.0 * ts_yz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxyzz[i] = 4.0 * ts_yz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xxzzz[i] = 4.0 * ts_yz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xyyyy[i] = 4.0 * ts_yz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xyyyz[i] = 4.0 * ts_yz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xyyzz[i] = 4.0 * ts_yz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xyzzz[i] = 4.0 * ts_yz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_xzzzz[i] = 4.0 * ts_yz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yyyyy[i] = 4.0 * ts_yz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yyz_yyyyz[i] = 4.0 * ts_yz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yyyzz[i] = 4.0 * ts_yz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yyzzz[i] = 4.0 * ts_yz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yzzzz[i] = 4.0 * ts_yz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yyz_zzzzz[i] = 4.0 * ts_yz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 378-399 components of targeted buffer : FH

    auto gs_y_yzz_xxxxx = pbuffer.data(idx_g_fh + 378);

    auto gs_y_yzz_xxxxy = pbuffer.data(idx_g_fh + 379);

    auto gs_y_yzz_xxxxz = pbuffer.data(idx_g_fh + 380);

    auto gs_y_yzz_xxxyy = pbuffer.data(idx_g_fh + 381);

    auto gs_y_yzz_xxxyz = pbuffer.data(idx_g_fh + 382);

    auto gs_y_yzz_xxxzz = pbuffer.data(idx_g_fh + 383);

    auto gs_y_yzz_xxyyy = pbuffer.data(idx_g_fh + 384);

    auto gs_y_yzz_xxyyz = pbuffer.data(idx_g_fh + 385);

    auto gs_y_yzz_xxyzz = pbuffer.data(idx_g_fh + 386);

    auto gs_y_yzz_xxzzz = pbuffer.data(idx_g_fh + 387);

    auto gs_y_yzz_xyyyy = pbuffer.data(idx_g_fh + 388);

    auto gs_y_yzz_xyyyz = pbuffer.data(idx_g_fh + 389);

    auto gs_y_yzz_xyyzz = pbuffer.data(idx_g_fh + 390);

    auto gs_y_yzz_xyzzz = pbuffer.data(idx_g_fh + 391);

    auto gs_y_yzz_xzzzz = pbuffer.data(idx_g_fh + 392);

    auto gs_y_yzz_yyyyy = pbuffer.data(idx_g_fh + 393);

    auto gs_y_yzz_yyyyz = pbuffer.data(idx_g_fh + 394);

    auto gs_y_yzz_yyyzz = pbuffer.data(idx_g_fh + 395);

    auto gs_y_yzz_yyzzz = pbuffer.data(idx_g_fh + 396);

    auto gs_y_yzz_yzzzz = pbuffer.data(idx_g_fh + 397);

    auto gs_y_yzz_zzzzz = pbuffer.data(idx_g_fh + 398);

    #pragma omp simd aligned(gc_y, gs_y_yzz_xxxxx, gs_y_yzz_xxxxy, gs_y_yzz_xxxxz, gs_y_yzz_xxxyy, gs_y_yzz_xxxyz, gs_y_yzz_xxxzz, gs_y_yzz_xxyyy, gs_y_yzz_xxyyz, gs_y_yzz_xxyzz, gs_y_yzz_xxzzz, gs_y_yzz_xyyyy, gs_y_yzz_xyyyz, gs_y_yzz_xyyzz, gs_y_yzz_xyzzz, gs_y_yzz_xzzzz, gs_y_yzz_yyyyy, gs_y_yzz_yyyyz, gs_y_yzz_yyyzz, gs_y_yzz_yyzzz, gs_y_yzz_yzzzz, gs_y_yzz_zzzzz, ts_yzz_xxxx, ts_yzz_xxxxx, ts_yzz_xxxxy, ts_yzz_xxxxz, ts_yzz_xxxy, ts_yzz_xxxyy, ts_yzz_xxxyz, ts_yzz_xxxz, ts_yzz_xxxzz, ts_yzz_xxyy, ts_yzz_xxyyy, ts_yzz_xxyyz, ts_yzz_xxyz, ts_yzz_xxyzz, ts_yzz_xxzz, ts_yzz_xxzzz, ts_yzz_xyyy, ts_yzz_xyyyy, ts_yzz_xyyyz, ts_yzz_xyyz, ts_yzz_xyyzz, ts_yzz_xyzz, ts_yzz_xyzzz, ts_yzz_xzzz, ts_yzz_xzzzz, ts_yzz_yyyy, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyz, ts_yzz_yyyzz, ts_yzz_yyzz, ts_yzz_yyzzz, ts_yzz_yzzz, ts_yzz_yzzzz, ts_yzz_zzzz, ts_yzz_zzzzz, ts_zz_xxxxx, ts_zz_xxxxy, ts_zz_xxxxz, ts_zz_xxxyy, ts_zz_xxxyz, ts_zz_xxxzz, ts_zz_xxyyy, ts_zz_xxyyz, ts_zz_xxyzz, ts_zz_xxzzz, ts_zz_xyyyy, ts_zz_xyyyz, ts_zz_xyyzz, ts_zz_xyzzz, ts_zz_xzzzz, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyzz, ts_zz_yyzzz, ts_zz_yzzzz, ts_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzz_xxxxx[i] = 2.0 * ts_zz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxxxy[i] = 2.0 * ts_zz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxxxz[i] = 2.0 * ts_zz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxxyy[i] = 2.0 * ts_zz_xxxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxxyz[i] = 2.0 * ts_zz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxxzz[i] = 2.0 * ts_zz_xxxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxyyy[i] = 2.0 * ts_zz_xxyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxyyz[i] = 2.0 * ts_zz_xxyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxyzz[i] = 2.0 * ts_zz_xxyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xxzzz[i] = 2.0 * ts_zz_xxzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xyyyy[i] = 2.0 * ts_zz_xyyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xyyyz[i] = 2.0 * ts_zz_xyyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xyyzz[i] = 2.0 * ts_zz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xyzzz[i] = 2.0 * ts_zz_xyzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_xzzzz[i] = 2.0 * ts_zz_xzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yyyyy[i] = 2.0 * ts_zz_yyyyy[i] * gfe_0 * tce_0 + 10.0 * ts_yzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_yzz_yyyyz[i] = 2.0 * ts_zz_yyyyz[i] * gfe_0 * tce_0 + 8.0 * ts_yzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yyyzz[i] = 2.0 * ts_zz_yyyzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yyzzz[i] = 2.0 * ts_zz_yyzzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yzzzz[i] = 2.0 * ts_zz_yzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_yzz_zzzzz[i] = 2.0 * ts_zz_zzzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 399-420 components of targeted buffer : FH

    auto gs_y_zzz_xxxxx = pbuffer.data(idx_g_fh + 399);

    auto gs_y_zzz_xxxxy = pbuffer.data(idx_g_fh + 400);

    auto gs_y_zzz_xxxxz = pbuffer.data(idx_g_fh + 401);

    auto gs_y_zzz_xxxyy = pbuffer.data(idx_g_fh + 402);

    auto gs_y_zzz_xxxyz = pbuffer.data(idx_g_fh + 403);

    auto gs_y_zzz_xxxzz = pbuffer.data(idx_g_fh + 404);

    auto gs_y_zzz_xxyyy = pbuffer.data(idx_g_fh + 405);

    auto gs_y_zzz_xxyyz = pbuffer.data(idx_g_fh + 406);

    auto gs_y_zzz_xxyzz = pbuffer.data(idx_g_fh + 407);

    auto gs_y_zzz_xxzzz = pbuffer.data(idx_g_fh + 408);

    auto gs_y_zzz_xyyyy = pbuffer.data(idx_g_fh + 409);

    auto gs_y_zzz_xyyyz = pbuffer.data(idx_g_fh + 410);

    auto gs_y_zzz_xyyzz = pbuffer.data(idx_g_fh + 411);

    auto gs_y_zzz_xyzzz = pbuffer.data(idx_g_fh + 412);

    auto gs_y_zzz_xzzzz = pbuffer.data(idx_g_fh + 413);

    auto gs_y_zzz_yyyyy = pbuffer.data(idx_g_fh + 414);

    auto gs_y_zzz_yyyyz = pbuffer.data(idx_g_fh + 415);

    auto gs_y_zzz_yyyzz = pbuffer.data(idx_g_fh + 416);

    auto gs_y_zzz_yyzzz = pbuffer.data(idx_g_fh + 417);

    auto gs_y_zzz_yzzzz = pbuffer.data(idx_g_fh + 418);

    auto gs_y_zzz_zzzzz = pbuffer.data(idx_g_fh + 419);

    #pragma omp simd aligned(gc_y, gs_y_zzz_xxxxx, gs_y_zzz_xxxxy, gs_y_zzz_xxxxz, gs_y_zzz_xxxyy, gs_y_zzz_xxxyz, gs_y_zzz_xxxzz, gs_y_zzz_xxyyy, gs_y_zzz_xxyyz, gs_y_zzz_xxyzz, gs_y_zzz_xxzzz, gs_y_zzz_xyyyy, gs_y_zzz_xyyyz, gs_y_zzz_xyyzz, gs_y_zzz_xyzzz, gs_y_zzz_xzzzz, gs_y_zzz_yyyyy, gs_y_zzz_yyyyz, gs_y_zzz_yyyzz, gs_y_zzz_yyzzz, gs_y_zzz_yzzzz, gs_y_zzz_zzzzz, ts_zzz_xxxx, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxy, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxz, ts_zzz_xxxzz, ts_zzz_xxyy, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyz, ts_zzz_xxyzz, ts_zzz_xxzz, ts_zzz_xxzzz, ts_zzz_xyyy, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyz, ts_zzz_xyyzz, ts_zzz_xyzz, ts_zzz_xyzzz, ts_zzz_xzzz, ts_zzz_xzzzz, ts_zzz_yyyy, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyz, ts_zzz_yyyzz, ts_zzz_yyzz, ts_zzz_yyzzz, ts_zzz_yzzz, ts_zzz_yzzzz, ts_zzz_zzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzz_xxxxx[i] = 2.0 * ts_zzz_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxxxy[i] = 2.0 * ts_zzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxxxz[i] = 2.0 * ts_zzz_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxxyy[i] = 4.0 * ts_zzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxxyz[i] = 2.0 * ts_zzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxxzz[i] = 2.0 * ts_zzz_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxyyy[i] = 6.0 * ts_zzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxyyz[i] = 4.0 * ts_zzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxyzz[i] = 2.0 * ts_zzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xxzzz[i] = 2.0 * ts_zzz_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xyyyy[i] = 8.0 * ts_zzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xyyyz[i] = 6.0 * ts_zzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xyyzz[i] = 4.0 * ts_zzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xyzzz[i] = 2.0 * ts_zzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_xzzzz[i] = 2.0 * ts_zzz_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yyyyy[i] = 10.0 * ts_zzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_zzz_yyyyz[i] = 8.0 * ts_zzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yyyzz[i] = 6.0 * ts_zzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yyzzz[i] = 4.0 * ts_zzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yzzzz[i] = 2.0 * ts_zzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_zzz_zzzzz[i] = 2.0 * ts_zzz_zzzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 420-441 components of targeted buffer : FH

    auto gs_z_xxx_xxxxx = pbuffer.data(idx_g_fh + 420);

    auto gs_z_xxx_xxxxy = pbuffer.data(idx_g_fh + 421);

    auto gs_z_xxx_xxxxz = pbuffer.data(idx_g_fh + 422);

    auto gs_z_xxx_xxxyy = pbuffer.data(idx_g_fh + 423);

    auto gs_z_xxx_xxxyz = pbuffer.data(idx_g_fh + 424);

    auto gs_z_xxx_xxxzz = pbuffer.data(idx_g_fh + 425);

    auto gs_z_xxx_xxyyy = pbuffer.data(idx_g_fh + 426);

    auto gs_z_xxx_xxyyz = pbuffer.data(idx_g_fh + 427);

    auto gs_z_xxx_xxyzz = pbuffer.data(idx_g_fh + 428);

    auto gs_z_xxx_xxzzz = pbuffer.data(idx_g_fh + 429);

    auto gs_z_xxx_xyyyy = pbuffer.data(idx_g_fh + 430);

    auto gs_z_xxx_xyyyz = pbuffer.data(idx_g_fh + 431);

    auto gs_z_xxx_xyyzz = pbuffer.data(idx_g_fh + 432);

    auto gs_z_xxx_xyzzz = pbuffer.data(idx_g_fh + 433);

    auto gs_z_xxx_xzzzz = pbuffer.data(idx_g_fh + 434);

    auto gs_z_xxx_yyyyy = pbuffer.data(idx_g_fh + 435);

    auto gs_z_xxx_yyyyz = pbuffer.data(idx_g_fh + 436);

    auto gs_z_xxx_yyyzz = pbuffer.data(idx_g_fh + 437);

    auto gs_z_xxx_yyzzz = pbuffer.data(idx_g_fh + 438);

    auto gs_z_xxx_yzzzz = pbuffer.data(idx_g_fh + 439);

    auto gs_z_xxx_zzzzz = pbuffer.data(idx_g_fh + 440);

    #pragma omp simd aligned(gc_z, gs_z_xxx_xxxxx, gs_z_xxx_xxxxy, gs_z_xxx_xxxxz, gs_z_xxx_xxxyy, gs_z_xxx_xxxyz, gs_z_xxx_xxxzz, gs_z_xxx_xxyyy, gs_z_xxx_xxyyz, gs_z_xxx_xxyzz, gs_z_xxx_xxzzz, gs_z_xxx_xyyyy, gs_z_xxx_xyyyz, gs_z_xxx_xyyzz, gs_z_xxx_xyzzz, gs_z_xxx_xzzzz, gs_z_xxx_yyyyy, gs_z_xxx_yyyyz, gs_z_xxx_yyyzz, gs_z_xxx_yyzzz, gs_z_xxx_yzzzz, gs_z_xxx_zzzzz, ts_xxx_xxxx, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxy, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxz, ts_xxx_xxxzz, ts_xxx_xxyy, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyz, ts_xxx_xxyzz, ts_xxx_xxzz, ts_xxx_xxzzz, ts_xxx_xyyy, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyz, ts_xxx_xyyzz, ts_xxx_xyzz, ts_xxx_xyzzz, ts_xxx_xzzz, ts_xxx_xzzzz, ts_xxx_yyyy, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyz, ts_xxx_yyyzz, ts_xxx_yyzz, ts_xxx_yyzzz, ts_xxx_yzzz, ts_xxx_yzzzz, ts_xxx_zzzz, ts_xxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxx_xxxxx[i] = 2.0 * ts_xxx_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxxxy[i] = 2.0 * ts_xxx_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxxxz[i] = 2.0 * ts_xxx_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxxyy[i] = 2.0 * ts_xxx_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxxyz[i] = 2.0 * ts_xxx_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxxzz[i] = 4.0 * ts_xxx_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxyyy[i] = 2.0 * ts_xxx_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxyyz[i] = 2.0 * ts_xxx_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxyzz[i] = 4.0 * ts_xxx_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xxzzz[i] = 6.0 * ts_xxx_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xyyyy[i] = 2.0 * ts_xxx_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xyyyz[i] = 2.0 * ts_xxx_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xyyzz[i] = 4.0 * ts_xxx_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xyzzz[i] = 6.0 * ts_xxx_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_xzzzz[i] = 8.0 * ts_xxx_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yyyyy[i] = 2.0 * ts_xxx_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxx_yyyyz[i] = 2.0 * ts_xxx_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yyyzz[i] = 4.0 * ts_xxx_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yyzzz[i] = 6.0 * ts_xxx_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yzzzz[i] = 8.0 * ts_xxx_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxx_zzzzz[i] = 10.0 * ts_xxx_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 441-462 components of targeted buffer : FH

    auto gs_z_xxy_xxxxx = pbuffer.data(idx_g_fh + 441);

    auto gs_z_xxy_xxxxy = pbuffer.data(idx_g_fh + 442);

    auto gs_z_xxy_xxxxz = pbuffer.data(idx_g_fh + 443);

    auto gs_z_xxy_xxxyy = pbuffer.data(idx_g_fh + 444);

    auto gs_z_xxy_xxxyz = pbuffer.data(idx_g_fh + 445);

    auto gs_z_xxy_xxxzz = pbuffer.data(idx_g_fh + 446);

    auto gs_z_xxy_xxyyy = pbuffer.data(idx_g_fh + 447);

    auto gs_z_xxy_xxyyz = pbuffer.data(idx_g_fh + 448);

    auto gs_z_xxy_xxyzz = pbuffer.data(idx_g_fh + 449);

    auto gs_z_xxy_xxzzz = pbuffer.data(idx_g_fh + 450);

    auto gs_z_xxy_xyyyy = pbuffer.data(idx_g_fh + 451);

    auto gs_z_xxy_xyyyz = pbuffer.data(idx_g_fh + 452);

    auto gs_z_xxy_xyyzz = pbuffer.data(idx_g_fh + 453);

    auto gs_z_xxy_xyzzz = pbuffer.data(idx_g_fh + 454);

    auto gs_z_xxy_xzzzz = pbuffer.data(idx_g_fh + 455);

    auto gs_z_xxy_yyyyy = pbuffer.data(idx_g_fh + 456);

    auto gs_z_xxy_yyyyz = pbuffer.data(idx_g_fh + 457);

    auto gs_z_xxy_yyyzz = pbuffer.data(idx_g_fh + 458);

    auto gs_z_xxy_yyzzz = pbuffer.data(idx_g_fh + 459);

    auto gs_z_xxy_yzzzz = pbuffer.data(idx_g_fh + 460);

    auto gs_z_xxy_zzzzz = pbuffer.data(idx_g_fh + 461);

    #pragma omp simd aligned(gc_z, gs_z_xxy_xxxxx, gs_z_xxy_xxxxy, gs_z_xxy_xxxxz, gs_z_xxy_xxxyy, gs_z_xxy_xxxyz, gs_z_xxy_xxxzz, gs_z_xxy_xxyyy, gs_z_xxy_xxyyz, gs_z_xxy_xxyzz, gs_z_xxy_xxzzz, gs_z_xxy_xyyyy, gs_z_xxy_xyyyz, gs_z_xxy_xyyzz, gs_z_xxy_xyzzz, gs_z_xxy_xzzzz, gs_z_xxy_yyyyy, gs_z_xxy_yyyyz, gs_z_xxy_yyyzz, gs_z_xxy_yyzzz, gs_z_xxy_yzzzz, gs_z_xxy_zzzzz, ts_xxy_xxxx, ts_xxy_xxxxx, ts_xxy_xxxxy, ts_xxy_xxxxz, ts_xxy_xxxy, ts_xxy_xxxyy, ts_xxy_xxxyz, ts_xxy_xxxz, ts_xxy_xxxzz, ts_xxy_xxyy, ts_xxy_xxyyy, ts_xxy_xxyyz, ts_xxy_xxyz, ts_xxy_xxyzz, ts_xxy_xxzz, ts_xxy_xxzzz, ts_xxy_xyyy, ts_xxy_xyyyy, ts_xxy_xyyyz, ts_xxy_xyyz, ts_xxy_xyyzz, ts_xxy_xyzz, ts_xxy_xyzzz, ts_xxy_xzzz, ts_xxy_xzzzz, ts_xxy_yyyy, ts_xxy_yyyyy, ts_xxy_yyyyz, ts_xxy_yyyz, ts_xxy_yyyzz, ts_xxy_yyzz, ts_xxy_yyzzz, ts_xxy_yzzz, ts_xxy_yzzzz, ts_xxy_zzzz, ts_xxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxy_xxxxx[i] = 2.0 * ts_xxy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxxxy[i] = 2.0 * ts_xxy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxxxz[i] = 2.0 * ts_xxy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxxyy[i] = 2.0 * ts_xxy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxxyz[i] = 2.0 * ts_xxy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxxzz[i] = 4.0 * ts_xxy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxyyy[i] = 2.0 * ts_xxy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxyyz[i] = 2.0 * ts_xxy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxyzz[i] = 4.0 * ts_xxy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xxzzz[i] = 6.0 * ts_xxy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xyyyy[i] = 2.0 * ts_xxy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xyyyz[i] = 2.0 * ts_xxy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xyyzz[i] = 4.0 * ts_xxy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xyzzz[i] = 6.0 * ts_xxy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_xzzzz[i] = 8.0 * ts_xxy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yyyyy[i] = 2.0 * ts_xxy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxy_yyyyz[i] = 2.0 * ts_xxy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yyyzz[i] = 4.0 * ts_xxy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yyzzz[i] = 6.0 * ts_xxy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yzzzz[i] = 8.0 * ts_xxy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxy_zzzzz[i] = 10.0 * ts_xxy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 462-483 components of targeted buffer : FH

    auto gs_z_xxz_xxxxx = pbuffer.data(idx_g_fh + 462);

    auto gs_z_xxz_xxxxy = pbuffer.data(idx_g_fh + 463);

    auto gs_z_xxz_xxxxz = pbuffer.data(idx_g_fh + 464);

    auto gs_z_xxz_xxxyy = pbuffer.data(idx_g_fh + 465);

    auto gs_z_xxz_xxxyz = pbuffer.data(idx_g_fh + 466);

    auto gs_z_xxz_xxxzz = pbuffer.data(idx_g_fh + 467);

    auto gs_z_xxz_xxyyy = pbuffer.data(idx_g_fh + 468);

    auto gs_z_xxz_xxyyz = pbuffer.data(idx_g_fh + 469);

    auto gs_z_xxz_xxyzz = pbuffer.data(idx_g_fh + 470);

    auto gs_z_xxz_xxzzz = pbuffer.data(idx_g_fh + 471);

    auto gs_z_xxz_xyyyy = pbuffer.data(idx_g_fh + 472);

    auto gs_z_xxz_xyyyz = pbuffer.data(idx_g_fh + 473);

    auto gs_z_xxz_xyyzz = pbuffer.data(idx_g_fh + 474);

    auto gs_z_xxz_xyzzz = pbuffer.data(idx_g_fh + 475);

    auto gs_z_xxz_xzzzz = pbuffer.data(idx_g_fh + 476);

    auto gs_z_xxz_yyyyy = pbuffer.data(idx_g_fh + 477);

    auto gs_z_xxz_yyyyz = pbuffer.data(idx_g_fh + 478);

    auto gs_z_xxz_yyyzz = pbuffer.data(idx_g_fh + 479);

    auto gs_z_xxz_yyzzz = pbuffer.data(idx_g_fh + 480);

    auto gs_z_xxz_yzzzz = pbuffer.data(idx_g_fh + 481);

    auto gs_z_xxz_zzzzz = pbuffer.data(idx_g_fh + 482);

    #pragma omp simd aligned(gc_z, gs_z_xxz_xxxxx, gs_z_xxz_xxxxy, gs_z_xxz_xxxxz, gs_z_xxz_xxxyy, gs_z_xxz_xxxyz, gs_z_xxz_xxxzz, gs_z_xxz_xxyyy, gs_z_xxz_xxyyz, gs_z_xxz_xxyzz, gs_z_xxz_xxzzz, gs_z_xxz_xyyyy, gs_z_xxz_xyyyz, gs_z_xxz_xyyzz, gs_z_xxz_xyzzz, gs_z_xxz_xzzzz, gs_z_xxz_yyyyy, gs_z_xxz_yyyyz, gs_z_xxz_yyyzz, gs_z_xxz_yyzzz, gs_z_xxz_yzzzz, gs_z_xxz_zzzzz, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxxz, ts_xx_xxxyy, ts_xx_xxxyz, ts_xx_xxxzz, ts_xx_xxyyy, ts_xx_xxyyz, ts_xx_xxyzz, ts_xx_xxzzz, ts_xx_xyyyy, ts_xx_xyyyz, ts_xx_xyyzz, ts_xx_xyzzz, ts_xx_xzzzz, ts_xx_yyyyy, ts_xx_yyyyz, ts_xx_yyyzz, ts_xx_yyzzz, ts_xx_yzzzz, ts_xx_zzzzz, ts_xxz_xxxx, ts_xxz_xxxxx, ts_xxz_xxxxy, ts_xxz_xxxxz, ts_xxz_xxxy, ts_xxz_xxxyy, ts_xxz_xxxyz, ts_xxz_xxxz, ts_xxz_xxxzz, ts_xxz_xxyy, ts_xxz_xxyyy, ts_xxz_xxyyz, ts_xxz_xxyz, ts_xxz_xxyzz, ts_xxz_xxzz, ts_xxz_xxzzz, ts_xxz_xyyy, ts_xxz_xyyyy, ts_xxz_xyyyz, ts_xxz_xyyz, ts_xxz_xyyzz, ts_xxz_xyzz, ts_xxz_xyzzz, ts_xxz_xzzz, ts_xxz_xzzzz, ts_xxz_yyyy, ts_xxz_yyyyy, ts_xxz_yyyyz, ts_xxz_yyyz, ts_xxz_yyyzz, ts_xxz_yyzz, ts_xxz_yyzzz, ts_xxz_yzzz, ts_xxz_yzzzz, ts_xxz_zzzz, ts_xxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxz_xxxxx[i] = 2.0 * ts_xx_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxxxy[i] = 2.0 * ts_xx_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxxxz[i] = 2.0 * ts_xx_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxxyy[i] = 2.0 * ts_xx_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxxyz[i] = 2.0 * ts_xx_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxxzz[i] = 2.0 * ts_xx_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxyyy[i] = 2.0 * ts_xx_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxyyz[i] = 2.0 * ts_xx_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxyzz[i] = 2.0 * ts_xx_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xxzzz[i] = 2.0 * ts_xx_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xyyyy[i] = 2.0 * ts_xx_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xyyyz[i] = 2.0 * ts_xx_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xyyzz[i] = 2.0 * ts_xx_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xyzzz[i] = 2.0 * ts_xx_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_xzzzz[i] = 2.0 * ts_xx_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yyyyy[i] = 2.0 * ts_xx_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xxz_yyyyz[i] = 2.0 * ts_xx_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yyyzz[i] = 2.0 * ts_xx_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yyzzz[i] = 2.0 * ts_xx_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xxz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yzzzz[i] = 2.0 * ts_xx_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xxz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xxz_zzzzz[i] = 2.0 * ts_xx_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xxz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 483-504 components of targeted buffer : FH

    auto gs_z_xyy_xxxxx = pbuffer.data(idx_g_fh + 483);

    auto gs_z_xyy_xxxxy = pbuffer.data(idx_g_fh + 484);

    auto gs_z_xyy_xxxxz = pbuffer.data(idx_g_fh + 485);

    auto gs_z_xyy_xxxyy = pbuffer.data(idx_g_fh + 486);

    auto gs_z_xyy_xxxyz = pbuffer.data(idx_g_fh + 487);

    auto gs_z_xyy_xxxzz = pbuffer.data(idx_g_fh + 488);

    auto gs_z_xyy_xxyyy = pbuffer.data(idx_g_fh + 489);

    auto gs_z_xyy_xxyyz = pbuffer.data(idx_g_fh + 490);

    auto gs_z_xyy_xxyzz = pbuffer.data(idx_g_fh + 491);

    auto gs_z_xyy_xxzzz = pbuffer.data(idx_g_fh + 492);

    auto gs_z_xyy_xyyyy = pbuffer.data(idx_g_fh + 493);

    auto gs_z_xyy_xyyyz = pbuffer.data(idx_g_fh + 494);

    auto gs_z_xyy_xyyzz = pbuffer.data(idx_g_fh + 495);

    auto gs_z_xyy_xyzzz = pbuffer.data(idx_g_fh + 496);

    auto gs_z_xyy_xzzzz = pbuffer.data(idx_g_fh + 497);

    auto gs_z_xyy_yyyyy = pbuffer.data(idx_g_fh + 498);

    auto gs_z_xyy_yyyyz = pbuffer.data(idx_g_fh + 499);

    auto gs_z_xyy_yyyzz = pbuffer.data(idx_g_fh + 500);

    auto gs_z_xyy_yyzzz = pbuffer.data(idx_g_fh + 501);

    auto gs_z_xyy_yzzzz = pbuffer.data(idx_g_fh + 502);

    auto gs_z_xyy_zzzzz = pbuffer.data(idx_g_fh + 503);

    #pragma omp simd aligned(gc_z, gs_z_xyy_xxxxx, gs_z_xyy_xxxxy, gs_z_xyy_xxxxz, gs_z_xyy_xxxyy, gs_z_xyy_xxxyz, gs_z_xyy_xxxzz, gs_z_xyy_xxyyy, gs_z_xyy_xxyyz, gs_z_xyy_xxyzz, gs_z_xyy_xxzzz, gs_z_xyy_xyyyy, gs_z_xyy_xyyyz, gs_z_xyy_xyyzz, gs_z_xyy_xyzzz, gs_z_xyy_xzzzz, gs_z_xyy_yyyyy, gs_z_xyy_yyyyz, gs_z_xyy_yyyzz, gs_z_xyy_yyzzz, gs_z_xyy_yzzzz, gs_z_xyy_zzzzz, ts_xyy_xxxx, ts_xyy_xxxxx, ts_xyy_xxxxy, ts_xyy_xxxxz, ts_xyy_xxxy, ts_xyy_xxxyy, ts_xyy_xxxyz, ts_xyy_xxxz, ts_xyy_xxxzz, ts_xyy_xxyy, ts_xyy_xxyyy, ts_xyy_xxyyz, ts_xyy_xxyz, ts_xyy_xxyzz, ts_xyy_xxzz, ts_xyy_xxzzz, ts_xyy_xyyy, ts_xyy_xyyyy, ts_xyy_xyyyz, ts_xyy_xyyz, ts_xyy_xyyzz, ts_xyy_xyzz, ts_xyy_xyzzz, ts_xyy_xzzz, ts_xyy_xzzzz, ts_xyy_yyyy, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyz, ts_xyy_yyyzz, ts_xyy_yyzz, ts_xyy_yyzzz, ts_xyy_yzzz, ts_xyy_yzzzz, ts_xyy_zzzz, ts_xyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyy_xxxxx[i] = 2.0 * ts_xyy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxxxy[i] = 2.0 * ts_xyy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxxxz[i] = 2.0 * ts_xyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxxyy[i] = 2.0 * ts_xyy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxxyz[i] = 2.0 * ts_xyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxxzz[i] = 4.0 * ts_xyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxyyy[i] = 2.0 * ts_xyy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxyyz[i] = 2.0 * ts_xyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxyzz[i] = 4.0 * ts_xyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xxzzz[i] = 6.0 * ts_xyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xyyyy[i] = 2.0 * ts_xyy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xyyyz[i] = 2.0 * ts_xyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xyyzz[i] = 4.0 * ts_xyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xyzzz[i] = 6.0 * ts_xyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_xzzzz[i] = 8.0 * ts_xyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yyyyy[i] = 2.0 * ts_xyy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyy_yyyyz[i] = 2.0 * ts_xyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yyyzz[i] = 4.0 * ts_xyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yyzzz[i] = 6.0 * ts_xyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yzzzz[i] = 8.0 * ts_xyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyy_zzzzz[i] = 10.0 * ts_xyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 504-525 components of targeted buffer : FH

    auto gs_z_xyz_xxxxx = pbuffer.data(idx_g_fh + 504);

    auto gs_z_xyz_xxxxy = pbuffer.data(idx_g_fh + 505);

    auto gs_z_xyz_xxxxz = pbuffer.data(idx_g_fh + 506);

    auto gs_z_xyz_xxxyy = pbuffer.data(idx_g_fh + 507);

    auto gs_z_xyz_xxxyz = pbuffer.data(idx_g_fh + 508);

    auto gs_z_xyz_xxxzz = pbuffer.data(idx_g_fh + 509);

    auto gs_z_xyz_xxyyy = pbuffer.data(idx_g_fh + 510);

    auto gs_z_xyz_xxyyz = pbuffer.data(idx_g_fh + 511);

    auto gs_z_xyz_xxyzz = pbuffer.data(idx_g_fh + 512);

    auto gs_z_xyz_xxzzz = pbuffer.data(idx_g_fh + 513);

    auto gs_z_xyz_xyyyy = pbuffer.data(idx_g_fh + 514);

    auto gs_z_xyz_xyyyz = pbuffer.data(idx_g_fh + 515);

    auto gs_z_xyz_xyyzz = pbuffer.data(idx_g_fh + 516);

    auto gs_z_xyz_xyzzz = pbuffer.data(idx_g_fh + 517);

    auto gs_z_xyz_xzzzz = pbuffer.data(idx_g_fh + 518);

    auto gs_z_xyz_yyyyy = pbuffer.data(idx_g_fh + 519);

    auto gs_z_xyz_yyyyz = pbuffer.data(idx_g_fh + 520);

    auto gs_z_xyz_yyyzz = pbuffer.data(idx_g_fh + 521);

    auto gs_z_xyz_yyzzz = pbuffer.data(idx_g_fh + 522);

    auto gs_z_xyz_yzzzz = pbuffer.data(idx_g_fh + 523);

    auto gs_z_xyz_zzzzz = pbuffer.data(idx_g_fh + 524);

    #pragma omp simd aligned(gc_z, gs_z_xyz_xxxxx, gs_z_xyz_xxxxy, gs_z_xyz_xxxxz, gs_z_xyz_xxxyy, gs_z_xyz_xxxyz, gs_z_xyz_xxxzz, gs_z_xyz_xxyyy, gs_z_xyz_xxyyz, gs_z_xyz_xxyzz, gs_z_xyz_xxzzz, gs_z_xyz_xyyyy, gs_z_xyz_xyyyz, gs_z_xyz_xyyzz, gs_z_xyz_xyzzz, gs_z_xyz_xzzzz, gs_z_xyz_yyyyy, gs_z_xyz_yyyyz, gs_z_xyz_yyyzz, gs_z_xyz_yyzzz, gs_z_xyz_yzzzz, gs_z_xyz_zzzzz, ts_xy_xxxxx, ts_xy_xxxxy, ts_xy_xxxxz, ts_xy_xxxyy, ts_xy_xxxyz, ts_xy_xxxzz, ts_xy_xxyyy, ts_xy_xxyyz, ts_xy_xxyzz, ts_xy_xxzzz, ts_xy_xyyyy, ts_xy_xyyyz, ts_xy_xyyzz, ts_xy_xyzzz, ts_xy_xzzzz, ts_xy_yyyyy, ts_xy_yyyyz, ts_xy_yyyzz, ts_xy_yyzzz, ts_xy_yzzzz, ts_xy_zzzzz, ts_xyz_xxxx, ts_xyz_xxxxx, ts_xyz_xxxxy, ts_xyz_xxxxz, ts_xyz_xxxy, ts_xyz_xxxyy, ts_xyz_xxxyz, ts_xyz_xxxz, ts_xyz_xxxzz, ts_xyz_xxyy, ts_xyz_xxyyy, ts_xyz_xxyyz, ts_xyz_xxyz, ts_xyz_xxyzz, ts_xyz_xxzz, ts_xyz_xxzzz, ts_xyz_xyyy, ts_xyz_xyyyy, ts_xyz_xyyyz, ts_xyz_xyyz, ts_xyz_xyyzz, ts_xyz_xyzz, ts_xyz_xyzzz, ts_xyz_xzzz, ts_xyz_xzzzz, ts_xyz_yyyy, ts_xyz_yyyyy, ts_xyz_yyyyz, ts_xyz_yyyz, ts_xyz_yyyzz, ts_xyz_yyzz, ts_xyz_yyzzz, ts_xyz_yzzz, ts_xyz_yzzzz, ts_xyz_zzzz, ts_xyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyz_xxxxx[i] = 2.0 * ts_xy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxxxy[i] = 2.0 * ts_xy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxxxz[i] = 2.0 * ts_xy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxxyy[i] = 2.0 * ts_xy_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxxyz[i] = 2.0 * ts_xy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxxzz[i] = 2.0 * ts_xy_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxyyy[i] = 2.0 * ts_xy_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxyyz[i] = 2.0 * ts_xy_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxyzz[i] = 2.0 * ts_xy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xxzzz[i] = 2.0 * ts_xy_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xyyyy[i] = 2.0 * ts_xy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xyyyz[i] = 2.0 * ts_xy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xyyzz[i] = 2.0 * ts_xy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xyzzz[i] = 2.0 * ts_xy_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_xzzzz[i] = 2.0 * ts_xy_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yyyyy[i] = 2.0 * ts_xy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xyz_yyyyz[i] = 2.0 * ts_xy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yyyzz[i] = 2.0 * ts_xy_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yyzzz[i] = 2.0 * ts_xy_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yzzzz[i] = 2.0 * ts_xy_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xyz_zzzzz[i] = 2.0 * ts_xy_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 525-546 components of targeted buffer : FH

    auto gs_z_xzz_xxxxx = pbuffer.data(idx_g_fh + 525);

    auto gs_z_xzz_xxxxy = pbuffer.data(idx_g_fh + 526);

    auto gs_z_xzz_xxxxz = pbuffer.data(idx_g_fh + 527);

    auto gs_z_xzz_xxxyy = pbuffer.data(idx_g_fh + 528);

    auto gs_z_xzz_xxxyz = pbuffer.data(idx_g_fh + 529);

    auto gs_z_xzz_xxxzz = pbuffer.data(idx_g_fh + 530);

    auto gs_z_xzz_xxyyy = pbuffer.data(idx_g_fh + 531);

    auto gs_z_xzz_xxyyz = pbuffer.data(idx_g_fh + 532);

    auto gs_z_xzz_xxyzz = pbuffer.data(idx_g_fh + 533);

    auto gs_z_xzz_xxzzz = pbuffer.data(idx_g_fh + 534);

    auto gs_z_xzz_xyyyy = pbuffer.data(idx_g_fh + 535);

    auto gs_z_xzz_xyyyz = pbuffer.data(idx_g_fh + 536);

    auto gs_z_xzz_xyyzz = pbuffer.data(idx_g_fh + 537);

    auto gs_z_xzz_xyzzz = pbuffer.data(idx_g_fh + 538);

    auto gs_z_xzz_xzzzz = pbuffer.data(idx_g_fh + 539);

    auto gs_z_xzz_yyyyy = pbuffer.data(idx_g_fh + 540);

    auto gs_z_xzz_yyyyz = pbuffer.data(idx_g_fh + 541);

    auto gs_z_xzz_yyyzz = pbuffer.data(idx_g_fh + 542);

    auto gs_z_xzz_yyzzz = pbuffer.data(idx_g_fh + 543);

    auto gs_z_xzz_yzzzz = pbuffer.data(idx_g_fh + 544);

    auto gs_z_xzz_zzzzz = pbuffer.data(idx_g_fh + 545);

    #pragma omp simd aligned(gc_z, gs_z_xzz_xxxxx, gs_z_xzz_xxxxy, gs_z_xzz_xxxxz, gs_z_xzz_xxxyy, gs_z_xzz_xxxyz, gs_z_xzz_xxxzz, gs_z_xzz_xxyyy, gs_z_xzz_xxyyz, gs_z_xzz_xxyzz, gs_z_xzz_xxzzz, gs_z_xzz_xyyyy, gs_z_xzz_xyyyz, gs_z_xzz_xyyzz, gs_z_xzz_xyzzz, gs_z_xzz_xzzzz, gs_z_xzz_yyyyy, gs_z_xzz_yyyyz, gs_z_xzz_yyyzz, gs_z_xzz_yyzzz, gs_z_xzz_yzzzz, gs_z_xzz_zzzzz, ts_xz_xxxxx, ts_xz_xxxxy, ts_xz_xxxxz, ts_xz_xxxyy, ts_xz_xxxyz, ts_xz_xxxzz, ts_xz_xxyyy, ts_xz_xxyyz, ts_xz_xxyzz, ts_xz_xxzzz, ts_xz_xyyyy, ts_xz_xyyyz, ts_xz_xyyzz, ts_xz_xyzzz, ts_xz_xzzzz, ts_xz_yyyyy, ts_xz_yyyyz, ts_xz_yyyzz, ts_xz_yyzzz, ts_xz_yzzzz, ts_xz_zzzzz, ts_xzz_xxxx, ts_xzz_xxxxx, ts_xzz_xxxxy, ts_xzz_xxxxz, ts_xzz_xxxy, ts_xzz_xxxyy, ts_xzz_xxxyz, ts_xzz_xxxz, ts_xzz_xxxzz, ts_xzz_xxyy, ts_xzz_xxyyy, ts_xzz_xxyyz, ts_xzz_xxyz, ts_xzz_xxyzz, ts_xzz_xxzz, ts_xzz_xxzzz, ts_xzz_xyyy, ts_xzz_xyyyy, ts_xzz_xyyyz, ts_xzz_xyyz, ts_xzz_xyyzz, ts_xzz_xyzz, ts_xzz_xyzzz, ts_xzz_xzzz, ts_xzz_xzzzz, ts_xzz_yyyy, ts_xzz_yyyyy, ts_xzz_yyyyz, ts_xzz_yyyz, ts_xzz_yyyzz, ts_xzz_yyzz, ts_xzz_yyzzz, ts_xzz_yzzz, ts_xzz_yzzzz, ts_xzz_zzzz, ts_xzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzz_xxxxx[i] = 4.0 * ts_xz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxxxy[i] = 4.0 * ts_xz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxxxz[i] = 4.0 * ts_xz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxxyy[i] = 4.0 * ts_xz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxxyz[i] = 4.0 * ts_xz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxxzz[i] = 4.0 * ts_xz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxyyy[i] = 4.0 * ts_xz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxyyz[i] = 4.0 * ts_xz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxyzz[i] = 4.0 * ts_xz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xxzzz[i] = 4.0 * ts_xz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xyyyy[i] = 4.0 * ts_xz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xyyyz[i] = 4.0 * ts_xz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xyyzz[i] = 4.0 * ts_xz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xyzzz[i] = 4.0 * ts_xz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_xzzzz[i] = 4.0 * ts_xz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yyyyy[i] = 4.0 * ts_xz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_xzz_yyyyz[i] = 4.0 * ts_xz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yyyzz[i] = 4.0 * ts_xz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yyzzz[i] = 4.0 * ts_xz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yzzzz[i] = 4.0 * ts_xz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_xzz_zzzzz[i] = 4.0 * ts_xz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_xzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 546-567 components of targeted buffer : FH

    auto gs_z_yyy_xxxxx = pbuffer.data(idx_g_fh + 546);

    auto gs_z_yyy_xxxxy = pbuffer.data(idx_g_fh + 547);

    auto gs_z_yyy_xxxxz = pbuffer.data(idx_g_fh + 548);

    auto gs_z_yyy_xxxyy = pbuffer.data(idx_g_fh + 549);

    auto gs_z_yyy_xxxyz = pbuffer.data(idx_g_fh + 550);

    auto gs_z_yyy_xxxzz = pbuffer.data(idx_g_fh + 551);

    auto gs_z_yyy_xxyyy = pbuffer.data(idx_g_fh + 552);

    auto gs_z_yyy_xxyyz = pbuffer.data(idx_g_fh + 553);

    auto gs_z_yyy_xxyzz = pbuffer.data(idx_g_fh + 554);

    auto gs_z_yyy_xxzzz = pbuffer.data(idx_g_fh + 555);

    auto gs_z_yyy_xyyyy = pbuffer.data(idx_g_fh + 556);

    auto gs_z_yyy_xyyyz = pbuffer.data(idx_g_fh + 557);

    auto gs_z_yyy_xyyzz = pbuffer.data(idx_g_fh + 558);

    auto gs_z_yyy_xyzzz = pbuffer.data(idx_g_fh + 559);

    auto gs_z_yyy_xzzzz = pbuffer.data(idx_g_fh + 560);

    auto gs_z_yyy_yyyyy = pbuffer.data(idx_g_fh + 561);

    auto gs_z_yyy_yyyyz = pbuffer.data(idx_g_fh + 562);

    auto gs_z_yyy_yyyzz = pbuffer.data(idx_g_fh + 563);

    auto gs_z_yyy_yyzzz = pbuffer.data(idx_g_fh + 564);

    auto gs_z_yyy_yzzzz = pbuffer.data(idx_g_fh + 565);

    auto gs_z_yyy_zzzzz = pbuffer.data(idx_g_fh + 566);

    #pragma omp simd aligned(gc_z, gs_z_yyy_xxxxx, gs_z_yyy_xxxxy, gs_z_yyy_xxxxz, gs_z_yyy_xxxyy, gs_z_yyy_xxxyz, gs_z_yyy_xxxzz, gs_z_yyy_xxyyy, gs_z_yyy_xxyyz, gs_z_yyy_xxyzz, gs_z_yyy_xxzzz, gs_z_yyy_xyyyy, gs_z_yyy_xyyyz, gs_z_yyy_xyyzz, gs_z_yyy_xyzzz, gs_z_yyy_xzzzz, gs_z_yyy_yyyyy, gs_z_yyy_yyyyz, gs_z_yyy_yyyzz, gs_z_yyy_yyzzz, gs_z_yyy_yzzzz, gs_z_yyy_zzzzz, ts_yyy_xxxx, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxy, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxz, ts_yyy_xxxzz, ts_yyy_xxyy, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyz, ts_yyy_xxyzz, ts_yyy_xxzz, ts_yyy_xxzzz, ts_yyy_xyyy, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyz, ts_yyy_xyyzz, ts_yyy_xyzz, ts_yyy_xyzzz, ts_yyy_xzzz, ts_yyy_xzzzz, ts_yyy_yyyy, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyz, ts_yyy_yyyzz, ts_yyy_yyzz, ts_yyy_yyzzz, ts_yyy_yzzz, ts_yyy_yzzzz, ts_yyy_zzzz, ts_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyy_xxxxx[i] = 2.0 * ts_yyy_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxxxy[i] = 2.0 * ts_yyy_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxxxz[i] = 2.0 * ts_yyy_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxxyy[i] = 2.0 * ts_yyy_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxxyz[i] = 2.0 * ts_yyy_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxxzz[i] = 4.0 * ts_yyy_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxyyy[i] = 2.0 * ts_yyy_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxyyz[i] = 2.0 * ts_yyy_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxyzz[i] = 4.0 * ts_yyy_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xxzzz[i] = 6.0 * ts_yyy_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xyyyy[i] = 2.0 * ts_yyy_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xyyyz[i] = 2.0 * ts_yyy_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xyyzz[i] = 4.0 * ts_yyy_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xyzzz[i] = 6.0 * ts_yyy_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_xzzzz[i] = 8.0 * ts_yyy_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yyyyy[i] = 2.0 * ts_yyy_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyy_yyyyz[i] = 2.0 * ts_yyy_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yyyzz[i] = 4.0 * ts_yyy_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yyzzz[i] = 6.0 * ts_yyy_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yzzzz[i] = 8.0 * ts_yyy_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyy_zzzzz[i] = 10.0 * ts_yyy_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 567-588 components of targeted buffer : FH

    auto gs_z_yyz_xxxxx = pbuffer.data(idx_g_fh + 567);

    auto gs_z_yyz_xxxxy = pbuffer.data(idx_g_fh + 568);

    auto gs_z_yyz_xxxxz = pbuffer.data(idx_g_fh + 569);

    auto gs_z_yyz_xxxyy = pbuffer.data(idx_g_fh + 570);

    auto gs_z_yyz_xxxyz = pbuffer.data(idx_g_fh + 571);

    auto gs_z_yyz_xxxzz = pbuffer.data(idx_g_fh + 572);

    auto gs_z_yyz_xxyyy = pbuffer.data(idx_g_fh + 573);

    auto gs_z_yyz_xxyyz = pbuffer.data(idx_g_fh + 574);

    auto gs_z_yyz_xxyzz = pbuffer.data(idx_g_fh + 575);

    auto gs_z_yyz_xxzzz = pbuffer.data(idx_g_fh + 576);

    auto gs_z_yyz_xyyyy = pbuffer.data(idx_g_fh + 577);

    auto gs_z_yyz_xyyyz = pbuffer.data(idx_g_fh + 578);

    auto gs_z_yyz_xyyzz = pbuffer.data(idx_g_fh + 579);

    auto gs_z_yyz_xyzzz = pbuffer.data(idx_g_fh + 580);

    auto gs_z_yyz_xzzzz = pbuffer.data(idx_g_fh + 581);

    auto gs_z_yyz_yyyyy = pbuffer.data(idx_g_fh + 582);

    auto gs_z_yyz_yyyyz = pbuffer.data(idx_g_fh + 583);

    auto gs_z_yyz_yyyzz = pbuffer.data(idx_g_fh + 584);

    auto gs_z_yyz_yyzzz = pbuffer.data(idx_g_fh + 585);

    auto gs_z_yyz_yzzzz = pbuffer.data(idx_g_fh + 586);

    auto gs_z_yyz_zzzzz = pbuffer.data(idx_g_fh + 587);

    #pragma omp simd aligned(gc_z, gs_z_yyz_xxxxx, gs_z_yyz_xxxxy, gs_z_yyz_xxxxz, gs_z_yyz_xxxyy, gs_z_yyz_xxxyz, gs_z_yyz_xxxzz, gs_z_yyz_xxyyy, gs_z_yyz_xxyyz, gs_z_yyz_xxyzz, gs_z_yyz_xxzzz, gs_z_yyz_xyyyy, gs_z_yyz_xyyyz, gs_z_yyz_xyyzz, gs_z_yyz_xyzzz, gs_z_yyz_xzzzz, gs_z_yyz_yyyyy, gs_z_yyz_yyyyz, gs_z_yyz_yyyzz, gs_z_yyz_yyzzz, gs_z_yyz_yzzzz, gs_z_yyz_zzzzz, ts_yy_xxxxx, ts_yy_xxxxy, ts_yy_xxxxz, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxxzz, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyzz, ts_yy_xxzzz, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyzz, ts_yy_xyzzz, ts_yy_xzzzz, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyzz, ts_yy_yyzzz, ts_yy_yzzzz, ts_yy_zzzzz, ts_yyz_xxxx, ts_yyz_xxxxx, ts_yyz_xxxxy, ts_yyz_xxxxz, ts_yyz_xxxy, ts_yyz_xxxyy, ts_yyz_xxxyz, ts_yyz_xxxz, ts_yyz_xxxzz, ts_yyz_xxyy, ts_yyz_xxyyy, ts_yyz_xxyyz, ts_yyz_xxyz, ts_yyz_xxyzz, ts_yyz_xxzz, ts_yyz_xxzzz, ts_yyz_xyyy, ts_yyz_xyyyy, ts_yyz_xyyyz, ts_yyz_xyyz, ts_yyz_xyyzz, ts_yyz_xyzz, ts_yyz_xyzzz, ts_yyz_xzzz, ts_yyz_xzzzz, ts_yyz_yyyy, ts_yyz_yyyyy, ts_yyz_yyyyz, ts_yyz_yyyz, ts_yyz_yyyzz, ts_yyz_yyzz, ts_yyz_yyzzz, ts_yyz_yzzz, ts_yyz_yzzzz, ts_yyz_zzzz, ts_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyz_xxxxx[i] = 2.0 * ts_yy_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxxxy[i] = 2.0 * ts_yy_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxxxz[i] = 2.0 * ts_yy_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxxyy[i] = 2.0 * ts_yy_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxxyz[i] = 2.0 * ts_yy_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxxzz[i] = 2.0 * ts_yy_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxyyy[i] = 2.0 * ts_yy_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxyyz[i] = 2.0 * ts_yy_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxyzz[i] = 2.0 * ts_yy_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xxzzz[i] = 2.0 * ts_yy_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xyyyy[i] = 2.0 * ts_yy_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xyyyz[i] = 2.0 * ts_yy_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xyyzz[i] = 2.0 * ts_yy_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xyzzz[i] = 2.0 * ts_yy_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_xzzzz[i] = 2.0 * ts_yy_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yyyyy[i] = 2.0 * ts_yy_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yyz_yyyyz[i] = 2.0 * ts_yy_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yyyzz[i] = 2.0 * ts_yy_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yyzzz[i] = 2.0 * ts_yy_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yyz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yzzzz[i] = 2.0 * ts_yy_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yyz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yyz_zzzzz[i] = 2.0 * ts_yy_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yyz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 588-609 components of targeted buffer : FH

    auto gs_z_yzz_xxxxx = pbuffer.data(idx_g_fh + 588);

    auto gs_z_yzz_xxxxy = pbuffer.data(idx_g_fh + 589);

    auto gs_z_yzz_xxxxz = pbuffer.data(idx_g_fh + 590);

    auto gs_z_yzz_xxxyy = pbuffer.data(idx_g_fh + 591);

    auto gs_z_yzz_xxxyz = pbuffer.data(idx_g_fh + 592);

    auto gs_z_yzz_xxxzz = pbuffer.data(idx_g_fh + 593);

    auto gs_z_yzz_xxyyy = pbuffer.data(idx_g_fh + 594);

    auto gs_z_yzz_xxyyz = pbuffer.data(idx_g_fh + 595);

    auto gs_z_yzz_xxyzz = pbuffer.data(idx_g_fh + 596);

    auto gs_z_yzz_xxzzz = pbuffer.data(idx_g_fh + 597);

    auto gs_z_yzz_xyyyy = pbuffer.data(idx_g_fh + 598);

    auto gs_z_yzz_xyyyz = pbuffer.data(idx_g_fh + 599);

    auto gs_z_yzz_xyyzz = pbuffer.data(idx_g_fh + 600);

    auto gs_z_yzz_xyzzz = pbuffer.data(idx_g_fh + 601);

    auto gs_z_yzz_xzzzz = pbuffer.data(idx_g_fh + 602);

    auto gs_z_yzz_yyyyy = pbuffer.data(idx_g_fh + 603);

    auto gs_z_yzz_yyyyz = pbuffer.data(idx_g_fh + 604);

    auto gs_z_yzz_yyyzz = pbuffer.data(idx_g_fh + 605);

    auto gs_z_yzz_yyzzz = pbuffer.data(idx_g_fh + 606);

    auto gs_z_yzz_yzzzz = pbuffer.data(idx_g_fh + 607);

    auto gs_z_yzz_zzzzz = pbuffer.data(idx_g_fh + 608);

    #pragma omp simd aligned(gc_z, gs_z_yzz_xxxxx, gs_z_yzz_xxxxy, gs_z_yzz_xxxxz, gs_z_yzz_xxxyy, gs_z_yzz_xxxyz, gs_z_yzz_xxxzz, gs_z_yzz_xxyyy, gs_z_yzz_xxyyz, gs_z_yzz_xxyzz, gs_z_yzz_xxzzz, gs_z_yzz_xyyyy, gs_z_yzz_xyyyz, gs_z_yzz_xyyzz, gs_z_yzz_xyzzz, gs_z_yzz_xzzzz, gs_z_yzz_yyyyy, gs_z_yzz_yyyyz, gs_z_yzz_yyyzz, gs_z_yzz_yyzzz, gs_z_yzz_yzzzz, gs_z_yzz_zzzzz, ts_yz_xxxxx, ts_yz_xxxxy, ts_yz_xxxxz, ts_yz_xxxyy, ts_yz_xxxyz, ts_yz_xxxzz, ts_yz_xxyyy, ts_yz_xxyyz, ts_yz_xxyzz, ts_yz_xxzzz, ts_yz_xyyyy, ts_yz_xyyyz, ts_yz_xyyzz, ts_yz_xyzzz, ts_yz_xzzzz, ts_yz_yyyyy, ts_yz_yyyyz, ts_yz_yyyzz, ts_yz_yyzzz, ts_yz_yzzzz, ts_yz_zzzzz, ts_yzz_xxxx, ts_yzz_xxxxx, ts_yzz_xxxxy, ts_yzz_xxxxz, ts_yzz_xxxy, ts_yzz_xxxyy, ts_yzz_xxxyz, ts_yzz_xxxz, ts_yzz_xxxzz, ts_yzz_xxyy, ts_yzz_xxyyy, ts_yzz_xxyyz, ts_yzz_xxyz, ts_yzz_xxyzz, ts_yzz_xxzz, ts_yzz_xxzzz, ts_yzz_xyyy, ts_yzz_xyyyy, ts_yzz_xyyyz, ts_yzz_xyyz, ts_yzz_xyyzz, ts_yzz_xyzz, ts_yzz_xyzzz, ts_yzz_xzzz, ts_yzz_xzzzz, ts_yzz_yyyy, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyz, ts_yzz_yyyzz, ts_yzz_yyzz, ts_yzz_yyzzz, ts_yzz_yzzz, ts_yzz_yzzzz, ts_yzz_zzzz, ts_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzz_xxxxx[i] = 4.0 * ts_yz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxxxy[i] = 4.0 * ts_yz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxxxz[i] = 4.0 * ts_yz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxxyy[i] = 4.0 * ts_yz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxxyz[i] = 4.0 * ts_yz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxxzz[i] = 4.0 * ts_yz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxyyy[i] = 4.0 * ts_yz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxyyz[i] = 4.0 * ts_yz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxyzz[i] = 4.0 * ts_yz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xxzzz[i] = 4.0 * ts_yz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xyyyy[i] = 4.0 * ts_yz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xyyyz[i] = 4.0 * ts_yz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xyyzz[i] = 4.0 * ts_yz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xyzzz[i] = 4.0 * ts_yz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_xzzzz[i] = 4.0 * ts_yz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yyyyy[i] = 4.0 * ts_yz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_yzz_yyyyz[i] = 4.0 * ts_yz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yyyzz[i] = 4.0 * ts_yz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yyzzz[i] = 4.0 * ts_yz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yzzzz[i] = 4.0 * ts_yz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_yzz_zzzzz[i] = 4.0 * ts_yz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_yzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zzzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 609-630 components of targeted buffer : FH

    auto gs_z_zzz_xxxxx = pbuffer.data(idx_g_fh + 609);

    auto gs_z_zzz_xxxxy = pbuffer.data(idx_g_fh + 610);

    auto gs_z_zzz_xxxxz = pbuffer.data(idx_g_fh + 611);

    auto gs_z_zzz_xxxyy = pbuffer.data(idx_g_fh + 612);

    auto gs_z_zzz_xxxyz = pbuffer.data(idx_g_fh + 613);

    auto gs_z_zzz_xxxzz = pbuffer.data(idx_g_fh + 614);

    auto gs_z_zzz_xxyyy = pbuffer.data(idx_g_fh + 615);

    auto gs_z_zzz_xxyyz = pbuffer.data(idx_g_fh + 616);

    auto gs_z_zzz_xxyzz = pbuffer.data(idx_g_fh + 617);

    auto gs_z_zzz_xxzzz = pbuffer.data(idx_g_fh + 618);

    auto gs_z_zzz_xyyyy = pbuffer.data(idx_g_fh + 619);

    auto gs_z_zzz_xyyyz = pbuffer.data(idx_g_fh + 620);

    auto gs_z_zzz_xyyzz = pbuffer.data(idx_g_fh + 621);

    auto gs_z_zzz_xyzzz = pbuffer.data(idx_g_fh + 622);

    auto gs_z_zzz_xzzzz = pbuffer.data(idx_g_fh + 623);

    auto gs_z_zzz_yyyyy = pbuffer.data(idx_g_fh + 624);

    auto gs_z_zzz_yyyyz = pbuffer.data(idx_g_fh + 625);

    auto gs_z_zzz_yyyzz = pbuffer.data(idx_g_fh + 626);

    auto gs_z_zzz_yyzzz = pbuffer.data(idx_g_fh + 627);

    auto gs_z_zzz_yzzzz = pbuffer.data(idx_g_fh + 628);

    auto gs_z_zzz_zzzzz = pbuffer.data(idx_g_fh + 629);

    #pragma omp simd aligned(gc_z, gs_z_zzz_xxxxx, gs_z_zzz_xxxxy, gs_z_zzz_xxxxz, gs_z_zzz_xxxyy, gs_z_zzz_xxxyz, gs_z_zzz_xxxzz, gs_z_zzz_xxyyy, gs_z_zzz_xxyyz, gs_z_zzz_xxyzz, gs_z_zzz_xxzzz, gs_z_zzz_xyyyy, gs_z_zzz_xyyyz, gs_z_zzz_xyyzz, gs_z_zzz_xyzzz, gs_z_zzz_xzzzz, gs_z_zzz_yyyyy, gs_z_zzz_yyyyz, gs_z_zzz_yyyzz, gs_z_zzz_yyzzz, gs_z_zzz_yzzzz, gs_z_zzz_zzzzz, ts_zz_xxxxx, ts_zz_xxxxy, ts_zz_xxxxz, ts_zz_xxxyy, ts_zz_xxxyz, ts_zz_xxxzz, ts_zz_xxyyy, ts_zz_xxyyz, ts_zz_xxyzz, ts_zz_xxzzz, ts_zz_xyyyy, ts_zz_xyyyz, ts_zz_xyyzz, ts_zz_xyzzz, ts_zz_xzzzz, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyzz, ts_zz_yyzzz, ts_zz_yzzzz, ts_zz_zzzzz, ts_zzz_xxxx, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxy, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxz, ts_zzz_xxxzz, ts_zzz_xxyy, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyz, ts_zzz_xxyzz, ts_zzz_xxzz, ts_zzz_xxzzz, ts_zzz_xyyy, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyz, ts_zzz_xyyzz, ts_zzz_xyzz, ts_zzz_xyzzz, ts_zzz_xzzz, ts_zzz_xzzzz, ts_zzz_yyyy, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyz, ts_zzz_yyyzz, ts_zzz_yyzz, ts_zzz_yyzzz, ts_zzz_yzzz, ts_zzz_yzzzz, ts_zzz_zzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzz_xxxxx[i] = 6.0 * ts_zz_xxxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxxxy[i] = 6.0 * ts_zz_xxxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxxxz[i] = 6.0 * ts_zz_xxxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxxyy[i] = 6.0 * ts_zz_xxxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxxyz[i] = 6.0 * ts_zz_xxxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxxzz[i] = 6.0 * ts_zz_xxxzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxyyy[i] = 6.0 * ts_zz_xxyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxyyz[i] = 6.0 * ts_zz_xxyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxyzz[i] = 6.0 * ts_zz_xxyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xxzzz[i] = 6.0 * ts_zz_xxzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzz_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xyyyy[i] = 6.0 * ts_zz_xyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xyyyz[i] = 6.0 * ts_zz_xyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xyyzz[i] = 6.0 * ts_zz_xyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xyzzz[i] = 6.0 * ts_zz_xyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzz_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_xzzzz[i] = 6.0 * ts_zz_xzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zzz_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yyyyy[i] = 6.0 * ts_zz_yyyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_zzz_yyyyz[i] = 6.0 * ts_zz_yyyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yyyzz[i] = 6.0 * ts_zz_yyyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yyzzz[i] = 6.0 * ts_zz_yyzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zzz_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yzzzz[i] = 6.0 * ts_zz_yzzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zzz_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_zzz_zzzzz[i] = 6.0 * ts_zz_zzzzz[i] * gfe_0 * tce_0 + 10.0 * ts_zzz_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_zzzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

