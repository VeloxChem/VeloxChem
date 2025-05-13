#include "ThreeCenterOverlapPrimRecGH.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_gh(CSimdArray<double>& pbuffer, 
                     const size_t idx_gh,
                     const size_t idx_dh,
                     const size_t idx_fg,
                     const size_t idx_fh,
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

    auto ts_xy_yyyyy = pbuffer.data(idx_dh + 36);

    auto ts_xy_yyyyz = pbuffer.data(idx_dh + 37);

    auto ts_xy_yyyzz = pbuffer.data(idx_dh + 38);

    auto ts_xy_yyzzz = pbuffer.data(idx_dh + 39);

    auto ts_xy_yzzzz = pbuffer.data(idx_dh + 40);

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

    auto ts_yz_xxxxz = pbuffer.data(idx_dh + 86);

    auto ts_yz_xxxzz = pbuffer.data(idx_dh + 89);

    auto ts_yz_xxzzz = pbuffer.data(idx_dh + 93);

    auto ts_yz_xzzzz = pbuffer.data(idx_dh + 98);

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

    auto ts_xxz_xxxz = pbuffer.data(idx_fg + 32);

    auto ts_xxz_xxyz = pbuffer.data(idx_fg + 34);

    auto ts_xxz_xxzz = pbuffer.data(idx_fg + 35);

    auto ts_xxz_xyyz = pbuffer.data(idx_fg + 37);

    auto ts_xxz_xyzz = pbuffer.data(idx_fg + 38);

    auto ts_xxz_xzzz = pbuffer.data(idx_fg + 39);

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

    auto ts_xzz_xxxz = pbuffer.data(idx_fg + 77);

    auto ts_xzz_xxyz = pbuffer.data(idx_fg + 79);

    auto ts_xzz_xxzz = pbuffer.data(idx_fg + 80);

    auto ts_xzz_xyyz = pbuffer.data(idx_fg + 82);

    auto ts_xzz_xyzz = pbuffer.data(idx_fg + 83);

    auto ts_xzz_xzzz = pbuffer.data(idx_fg + 84);

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

    auto ts_yyz_xxxz = pbuffer.data(idx_fg + 107);

    auto ts_yyz_xxyz = pbuffer.data(idx_fg + 109);

    auto ts_yyz_xxzz = pbuffer.data(idx_fg + 110);

    auto ts_yyz_xyyz = pbuffer.data(idx_fg + 112);

    auto ts_yyz_xyzz = pbuffer.data(idx_fg + 113);

    auto ts_yyz_xzzz = pbuffer.data(idx_fg + 114);

    auto ts_yyz_yyyz = pbuffer.data(idx_fg + 116);

    auto ts_yyz_yyzz = pbuffer.data(idx_fg + 117);

    auto ts_yyz_yzzz = pbuffer.data(idx_fg + 118);

    auto ts_yyz_zzzz = pbuffer.data(idx_fg + 119);

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

    auto ts_xxy_xxxzz = pbuffer.data(idx_fh + 26);

    auto ts_xxy_xxyyy = pbuffer.data(idx_fh + 27);

    auto ts_xxy_xxzzz = pbuffer.data(idx_fh + 30);

    auto ts_xxy_xyyyy = pbuffer.data(idx_fh + 31);

    auto ts_xxy_xzzzz = pbuffer.data(idx_fh + 35);

    auto ts_xxy_yyyyy = pbuffer.data(idx_fh + 36);

    auto ts_xxy_yyyyz = pbuffer.data(idx_fh + 37);

    auto ts_xxy_yyyzz = pbuffer.data(idx_fh + 38);

    auto ts_xxy_yyzzz = pbuffer.data(idx_fh + 39);

    auto ts_xxy_yzzzz = pbuffer.data(idx_fh + 40);

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

    auto ts_xxz_yyyyz = pbuffer.data(idx_fh + 58);

    auto ts_xxz_yyyzz = pbuffer.data(idx_fh + 59);

    auto ts_xxz_yyzzz = pbuffer.data(idx_fh + 60);

    auto ts_xxz_yzzzz = pbuffer.data(idx_fh + 61);

    auto ts_xxz_zzzzz = pbuffer.data(idx_fh + 62);

    auto ts_xyy_xxxxx = pbuffer.data(idx_fh + 63);

    auto ts_xyy_xxxxy = pbuffer.data(idx_fh + 64);

    auto ts_xyy_xxxyy = pbuffer.data(idx_fh + 66);

    auto ts_xyy_xxxyz = pbuffer.data(idx_fh + 67);

    auto ts_xyy_xxyyy = pbuffer.data(idx_fh + 69);

    auto ts_xyy_xxyyz = pbuffer.data(idx_fh + 70);

    auto ts_xyy_xxyzz = pbuffer.data(idx_fh + 71);

    auto ts_xyy_xyyyy = pbuffer.data(idx_fh + 73);

    auto ts_xyy_xyyyz = pbuffer.data(idx_fh + 74);

    auto ts_xyy_xyyzz = pbuffer.data(idx_fh + 75);

    auto ts_xyy_xyzzz = pbuffer.data(idx_fh + 76);

    auto ts_xyy_yyyyy = pbuffer.data(idx_fh + 78);

    auto ts_xyy_yyyyz = pbuffer.data(idx_fh + 79);

    auto ts_xyy_yyyzz = pbuffer.data(idx_fh + 80);

    auto ts_xyy_yyzzz = pbuffer.data(idx_fh + 81);

    auto ts_xyy_yzzzz = pbuffer.data(idx_fh + 82);

    auto ts_xyy_zzzzz = pbuffer.data(idx_fh + 83);

    auto ts_xyz_yyyyz = pbuffer.data(idx_fh + 100);

    auto ts_xyz_yyyzz = pbuffer.data(idx_fh + 101);

    auto ts_xyz_yyzzz = pbuffer.data(idx_fh + 102);

    auto ts_xyz_yzzzz = pbuffer.data(idx_fh + 103);

    auto ts_xzz_xxxxx = pbuffer.data(idx_fh + 105);

    auto ts_xzz_xxxxz = pbuffer.data(idx_fh + 107);

    auto ts_xzz_xxxyz = pbuffer.data(idx_fh + 109);

    auto ts_xzz_xxxzz = pbuffer.data(idx_fh + 110);

    auto ts_xzz_xxyyz = pbuffer.data(idx_fh + 112);

    auto ts_xzz_xxyzz = pbuffer.data(idx_fh + 113);

    auto ts_xzz_xxzzz = pbuffer.data(idx_fh + 114);

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

    // Set up 0-21 components of targeted buffer : GH

    auto ts_xxxx_xxxxx = pbuffer.data(idx_gh);

    auto ts_xxxx_xxxxy = pbuffer.data(idx_gh + 1);

    auto ts_xxxx_xxxxz = pbuffer.data(idx_gh + 2);

    auto ts_xxxx_xxxyy = pbuffer.data(idx_gh + 3);

    auto ts_xxxx_xxxyz = pbuffer.data(idx_gh + 4);

    auto ts_xxxx_xxxzz = pbuffer.data(idx_gh + 5);

    auto ts_xxxx_xxyyy = pbuffer.data(idx_gh + 6);

    auto ts_xxxx_xxyyz = pbuffer.data(idx_gh + 7);

    auto ts_xxxx_xxyzz = pbuffer.data(idx_gh + 8);

    auto ts_xxxx_xxzzz = pbuffer.data(idx_gh + 9);

    auto ts_xxxx_xyyyy = pbuffer.data(idx_gh + 10);

    auto ts_xxxx_xyyyz = pbuffer.data(idx_gh + 11);

    auto ts_xxxx_xyyzz = pbuffer.data(idx_gh + 12);

    auto ts_xxxx_xyzzz = pbuffer.data(idx_gh + 13);

    auto ts_xxxx_xzzzz = pbuffer.data(idx_gh + 14);

    auto ts_xxxx_yyyyy = pbuffer.data(idx_gh + 15);

    auto ts_xxxx_yyyyz = pbuffer.data(idx_gh + 16);

    auto ts_xxxx_yyyzz = pbuffer.data(idx_gh + 17);

    auto ts_xxxx_yyzzz = pbuffer.data(idx_gh + 18);

    auto ts_xxxx_yzzzz = pbuffer.data(idx_gh + 19);

    auto ts_xxxx_zzzzz = pbuffer.data(idx_gh + 20);

    #pragma omp simd aligned(ga_x, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxxz, ts_xx_xxxyy, ts_xx_xxxyz, ts_xx_xxxzz, ts_xx_xxyyy, ts_xx_xxyyz, ts_xx_xxyzz, ts_xx_xxzzz, ts_xx_xyyyy, ts_xx_xyyyz, ts_xx_xyyzz, ts_xx_xyzzz, ts_xx_xzzzz, ts_xx_yyyyy, ts_xx_yyyyz, ts_xx_yyyzz, ts_xx_yyzzz, ts_xx_yzzzz, ts_xx_zzzzz, ts_xxx_xxxx, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxy, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxz, ts_xxx_xxxzz, ts_xxx_xxyy, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyz, ts_xxx_xxyzz, ts_xxx_xxzz, ts_xxx_xxzzz, ts_xxx_xyyy, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyz, ts_xxx_xyyzz, ts_xxx_xyzz, ts_xxx_xyzzz, ts_xxx_xzzz, ts_xxx_xzzzz, ts_xxx_yyyy, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyz, ts_xxx_yyyzz, ts_xxx_yyzz, ts_xxx_yyzzz, ts_xxx_yzzz, ts_xxx_yzzzz, ts_xxx_zzzz, ts_xxx_zzzzz, ts_xxxx_xxxxx, ts_xxxx_xxxxy, ts_xxxx_xxxxz, ts_xxxx_xxxyy, ts_xxxx_xxxyz, ts_xxxx_xxxzz, ts_xxxx_xxyyy, ts_xxxx_xxyyz, ts_xxxx_xxyzz, ts_xxxx_xxzzz, ts_xxxx_xyyyy, ts_xxxx_xyyyz, ts_xxxx_xyyzz, ts_xxxx_xyzzz, ts_xxxx_xzzzz, ts_xxxx_yyyyy, ts_xxxx_yyyyz, ts_xxxx_yyyzz, ts_xxxx_yyzzz, ts_xxxx_yzzzz, ts_xxxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxx_xxxxx[i] = 3.0 * ts_xx_xxxxx[i] * gfe_0 + 5.0 * ts_xxx_xxxx[i] * gfe_0 + ts_xxx_xxxxx[i] * ga_x[i];

        ts_xxxx_xxxxy[i] = 3.0 * ts_xx_xxxxy[i] * gfe_0 + 4.0 * ts_xxx_xxxy[i] * gfe_0 + ts_xxx_xxxxy[i] * ga_x[i];

        ts_xxxx_xxxxz[i] = 3.0 * ts_xx_xxxxz[i] * gfe_0 + 4.0 * ts_xxx_xxxz[i] * gfe_0 + ts_xxx_xxxxz[i] * ga_x[i];

        ts_xxxx_xxxyy[i] = 3.0 * ts_xx_xxxyy[i] * gfe_0 + 3.0 * ts_xxx_xxyy[i] * gfe_0 + ts_xxx_xxxyy[i] * ga_x[i];

        ts_xxxx_xxxyz[i] = 3.0 * ts_xx_xxxyz[i] * gfe_0 + 3.0 * ts_xxx_xxyz[i] * gfe_0 + ts_xxx_xxxyz[i] * ga_x[i];

        ts_xxxx_xxxzz[i] = 3.0 * ts_xx_xxxzz[i] * gfe_0 + 3.0 * ts_xxx_xxzz[i] * gfe_0 + ts_xxx_xxxzz[i] * ga_x[i];

        ts_xxxx_xxyyy[i] = 3.0 * ts_xx_xxyyy[i] * gfe_0 + 2.0 * ts_xxx_xyyy[i] * gfe_0 + ts_xxx_xxyyy[i] * ga_x[i];

        ts_xxxx_xxyyz[i] = 3.0 * ts_xx_xxyyz[i] * gfe_0 + 2.0 * ts_xxx_xyyz[i] * gfe_0 + ts_xxx_xxyyz[i] * ga_x[i];

        ts_xxxx_xxyzz[i] = 3.0 * ts_xx_xxyzz[i] * gfe_0 + 2.0 * ts_xxx_xyzz[i] * gfe_0 + ts_xxx_xxyzz[i] * ga_x[i];

        ts_xxxx_xxzzz[i] = 3.0 * ts_xx_xxzzz[i] * gfe_0 + 2.0 * ts_xxx_xzzz[i] * gfe_0 + ts_xxx_xxzzz[i] * ga_x[i];

        ts_xxxx_xyyyy[i] = 3.0 * ts_xx_xyyyy[i] * gfe_0 + ts_xxx_yyyy[i] * gfe_0 + ts_xxx_xyyyy[i] * ga_x[i];

        ts_xxxx_xyyyz[i] = 3.0 * ts_xx_xyyyz[i] * gfe_0 + ts_xxx_yyyz[i] * gfe_0 + ts_xxx_xyyyz[i] * ga_x[i];

        ts_xxxx_xyyzz[i] = 3.0 * ts_xx_xyyzz[i] * gfe_0 + ts_xxx_yyzz[i] * gfe_0 + ts_xxx_xyyzz[i] * ga_x[i];

        ts_xxxx_xyzzz[i] = 3.0 * ts_xx_xyzzz[i] * gfe_0 + ts_xxx_yzzz[i] * gfe_0 + ts_xxx_xyzzz[i] * ga_x[i];

        ts_xxxx_xzzzz[i] = 3.0 * ts_xx_xzzzz[i] * gfe_0 + ts_xxx_zzzz[i] * gfe_0 + ts_xxx_xzzzz[i] * ga_x[i];

        ts_xxxx_yyyyy[i] = 3.0 * ts_xx_yyyyy[i] * gfe_0 + ts_xxx_yyyyy[i] * ga_x[i];

        ts_xxxx_yyyyz[i] = 3.0 * ts_xx_yyyyz[i] * gfe_0 + ts_xxx_yyyyz[i] * ga_x[i];

        ts_xxxx_yyyzz[i] = 3.0 * ts_xx_yyyzz[i] * gfe_0 + ts_xxx_yyyzz[i] * ga_x[i];

        ts_xxxx_yyzzz[i] = 3.0 * ts_xx_yyzzz[i] * gfe_0 + ts_xxx_yyzzz[i] * ga_x[i];

        ts_xxxx_yzzzz[i] = 3.0 * ts_xx_yzzzz[i] * gfe_0 + ts_xxx_yzzzz[i] * ga_x[i];

        ts_xxxx_zzzzz[i] = 3.0 * ts_xx_zzzzz[i] * gfe_0 + ts_xxx_zzzzz[i] * ga_x[i];
    }

    // Set up 21-42 components of targeted buffer : GH

    auto ts_xxxy_xxxxx = pbuffer.data(idx_gh + 21);

    auto ts_xxxy_xxxxy = pbuffer.data(idx_gh + 22);

    auto ts_xxxy_xxxxz = pbuffer.data(idx_gh + 23);

    auto ts_xxxy_xxxyy = pbuffer.data(idx_gh + 24);

    auto ts_xxxy_xxxyz = pbuffer.data(idx_gh + 25);

    auto ts_xxxy_xxxzz = pbuffer.data(idx_gh + 26);

    auto ts_xxxy_xxyyy = pbuffer.data(idx_gh + 27);

    auto ts_xxxy_xxyyz = pbuffer.data(idx_gh + 28);

    auto ts_xxxy_xxyzz = pbuffer.data(idx_gh + 29);

    auto ts_xxxy_xxzzz = pbuffer.data(idx_gh + 30);

    auto ts_xxxy_xyyyy = pbuffer.data(idx_gh + 31);

    auto ts_xxxy_xyyyz = pbuffer.data(idx_gh + 32);

    auto ts_xxxy_xyyzz = pbuffer.data(idx_gh + 33);

    auto ts_xxxy_xyzzz = pbuffer.data(idx_gh + 34);

    auto ts_xxxy_xzzzz = pbuffer.data(idx_gh + 35);

    auto ts_xxxy_yyyyy = pbuffer.data(idx_gh + 36);

    auto ts_xxxy_yyyyz = pbuffer.data(idx_gh + 37);

    auto ts_xxxy_yyyzz = pbuffer.data(idx_gh + 38);

    auto ts_xxxy_yyzzz = pbuffer.data(idx_gh + 39);

    auto ts_xxxy_yzzzz = pbuffer.data(idx_gh + 40);

    auto ts_xxxy_zzzzz = pbuffer.data(idx_gh + 41);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxx_xxxx, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxy, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxz, ts_xxx_xxxzz, ts_xxx_xxyy, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyz, ts_xxx_xxyzz, ts_xxx_xxzz, ts_xxx_xxzzz, ts_xxx_xyyy, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyz, ts_xxx_xyyzz, ts_xxx_xyzz, ts_xxx_xyzzz, ts_xxx_xzzz, ts_xxx_xzzzz, ts_xxx_zzzzz, ts_xxxy_xxxxx, ts_xxxy_xxxxy, ts_xxxy_xxxxz, ts_xxxy_xxxyy, ts_xxxy_xxxyz, ts_xxxy_xxxzz, ts_xxxy_xxyyy, ts_xxxy_xxyyz, ts_xxxy_xxyzz, ts_xxxy_xxzzz, ts_xxxy_xyyyy, ts_xxxy_xyyyz, ts_xxxy_xyyzz, ts_xxxy_xyzzz, ts_xxxy_xzzzz, ts_xxxy_yyyyy, ts_xxxy_yyyyz, ts_xxxy_yyyzz, ts_xxxy_yyzzz, ts_xxxy_yzzzz, ts_xxxy_zzzzz, ts_xxy_yyyyy, ts_xxy_yyyyz, ts_xxy_yyyzz, ts_xxy_yyzzz, ts_xxy_yzzzz, ts_xy_yyyyy, ts_xy_yyyyz, ts_xy_yyyzz, ts_xy_yyzzz, ts_xy_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxy_xxxxx[i] = ts_xxx_xxxxx[i] * ga_y[i];

        ts_xxxy_xxxxy[i] = ts_xxx_xxxx[i] * gfe_0 + ts_xxx_xxxxy[i] * ga_y[i];

        ts_xxxy_xxxxz[i] = ts_xxx_xxxxz[i] * ga_y[i];

        ts_xxxy_xxxyy[i] = 2.0 * ts_xxx_xxxy[i] * gfe_0 + ts_xxx_xxxyy[i] * ga_y[i];

        ts_xxxy_xxxyz[i] = ts_xxx_xxxz[i] * gfe_0 + ts_xxx_xxxyz[i] * ga_y[i];

        ts_xxxy_xxxzz[i] = ts_xxx_xxxzz[i] * ga_y[i];

        ts_xxxy_xxyyy[i] = 3.0 * ts_xxx_xxyy[i] * gfe_0 + ts_xxx_xxyyy[i] * ga_y[i];

        ts_xxxy_xxyyz[i] = 2.0 * ts_xxx_xxyz[i] * gfe_0 + ts_xxx_xxyyz[i] * ga_y[i];

        ts_xxxy_xxyzz[i] = ts_xxx_xxzz[i] * gfe_0 + ts_xxx_xxyzz[i] * ga_y[i];

        ts_xxxy_xxzzz[i] = ts_xxx_xxzzz[i] * ga_y[i];

        ts_xxxy_xyyyy[i] = 4.0 * ts_xxx_xyyy[i] * gfe_0 + ts_xxx_xyyyy[i] * ga_y[i];

        ts_xxxy_xyyyz[i] = 3.0 * ts_xxx_xyyz[i] * gfe_0 + ts_xxx_xyyyz[i] * ga_y[i];

        ts_xxxy_xyyzz[i] = 2.0 * ts_xxx_xyzz[i] * gfe_0 + ts_xxx_xyyzz[i] * ga_y[i];

        ts_xxxy_xyzzz[i] = ts_xxx_xzzz[i] * gfe_0 + ts_xxx_xyzzz[i] * ga_y[i];

        ts_xxxy_xzzzz[i] = ts_xxx_xzzzz[i] * ga_y[i];

        ts_xxxy_yyyyy[i] = 2.0 * ts_xy_yyyyy[i] * gfe_0 + ts_xxy_yyyyy[i] * ga_x[i];

        ts_xxxy_yyyyz[i] = 2.0 * ts_xy_yyyyz[i] * gfe_0 + ts_xxy_yyyyz[i] * ga_x[i];

        ts_xxxy_yyyzz[i] = 2.0 * ts_xy_yyyzz[i] * gfe_0 + ts_xxy_yyyzz[i] * ga_x[i];

        ts_xxxy_yyzzz[i] = 2.0 * ts_xy_yyzzz[i] * gfe_0 + ts_xxy_yyzzz[i] * ga_x[i];

        ts_xxxy_yzzzz[i] = 2.0 * ts_xy_yzzzz[i] * gfe_0 + ts_xxy_yzzzz[i] * ga_x[i];

        ts_xxxy_zzzzz[i] = ts_xxx_zzzzz[i] * ga_y[i];
    }

    // Set up 42-63 components of targeted buffer : GH

    auto ts_xxxz_xxxxx = pbuffer.data(idx_gh + 42);

    auto ts_xxxz_xxxxy = pbuffer.data(idx_gh + 43);

    auto ts_xxxz_xxxxz = pbuffer.data(idx_gh + 44);

    auto ts_xxxz_xxxyy = pbuffer.data(idx_gh + 45);

    auto ts_xxxz_xxxyz = pbuffer.data(idx_gh + 46);

    auto ts_xxxz_xxxzz = pbuffer.data(idx_gh + 47);

    auto ts_xxxz_xxyyy = pbuffer.data(idx_gh + 48);

    auto ts_xxxz_xxyyz = pbuffer.data(idx_gh + 49);

    auto ts_xxxz_xxyzz = pbuffer.data(idx_gh + 50);

    auto ts_xxxz_xxzzz = pbuffer.data(idx_gh + 51);

    auto ts_xxxz_xyyyy = pbuffer.data(idx_gh + 52);

    auto ts_xxxz_xyyyz = pbuffer.data(idx_gh + 53);

    auto ts_xxxz_xyyzz = pbuffer.data(idx_gh + 54);

    auto ts_xxxz_xyzzz = pbuffer.data(idx_gh + 55);

    auto ts_xxxz_xzzzz = pbuffer.data(idx_gh + 56);

    auto ts_xxxz_yyyyy = pbuffer.data(idx_gh + 57);

    auto ts_xxxz_yyyyz = pbuffer.data(idx_gh + 58);

    auto ts_xxxz_yyyzz = pbuffer.data(idx_gh + 59);

    auto ts_xxxz_yyzzz = pbuffer.data(idx_gh + 60);

    auto ts_xxxz_yzzzz = pbuffer.data(idx_gh + 61);

    auto ts_xxxz_zzzzz = pbuffer.data(idx_gh + 62);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxx_xxxx, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxy, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxz, ts_xxx_xxxzz, ts_xxx_xxyy, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyz, ts_xxx_xxyzz, ts_xxx_xxzz, ts_xxx_xxzzz, ts_xxx_xyyy, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyz, ts_xxx_xyyzz, ts_xxx_xyzz, ts_xxx_xyzzz, ts_xxx_xzzz, ts_xxx_xzzzz, ts_xxx_yyyyy, ts_xxxz_xxxxx, ts_xxxz_xxxxy, ts_xxxz_xxxxz, ts_xxxz_xxxyy, ts_xxxz_xxxyz, ts_xxxz_xxxzz, ts_xxxz_xxyyy, ts_xxxz_xxyyz, ts_xxxz_xxyzz, ts_xxxz_xxzzz, ts_xxxz_xyyyy, ts_xxxz_xyyyz, ts_xxxz_xyyzz, ts_xxxz_xyzzz, ts_xxxz_xzzzz, ts_xxxz_yyyyy, ts_xxxz_yyyyz, ts_xxxz_yyyzz, ts_xxxz_yyzzz, ts_xxxz_yzzzz, ts_xxxz_zzzzz, ts_xxz_yyyyz, ts_xxz_yyyzz, ts_xxz_yyzzz, ts_xxz_yzzzz, ts_xxz_zzzzz, ts_xz_yyyyz, ts_xz_yyyzz, ts_xz_yyzzz, ts_xz_yzzzz, ts_xz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxz_xxxxx[i] = ts_xxx_xxxxx[i] * ga_z[i];

        ts_xxxz_xxxxy[i] = ts_xxx_xxxxy[i] * ga_z[i];

        ts_xxxz_xxxxz[i] = ts_xxx_xxxx[i] * gfe_0 + ts_xxx_xxxxz[i] * ga_z[i];

        ts_xxxz_xxxyy[i] = ts_xxx_xxxyy[i] * ga_z[i];

        ts_xxxz_xxxyz[i] = ts_xxx_xxxy[i] * gfe_0 + ts_xxx_xxxyz[i] * ga_z[i];

        ts_xxxz_xxxzz[i] = 2.0 * ts_xxx_xxxz[i] * gfe_0 + ts_xxx_xxxzz[i] * ga_z[i];

        ts_xxxz_xxyyy[i] = ts_xxx_xxyyy[i] * ga_z[i];

        ts_xxxz_xxyyz[i] = ts_xxx_xxyy[i] * gfe_0 + ts_xxx_xxyyz[i] * ga_z[i];

        ts_xxxz_xxyzz[i] = 2.0 * ts_xxx_xxyz[i] * gfe_0 + ts_xxx_xxyzz[i] * ga_z[i];

        ts_xxxz_xxzzz[i] = 3.0 * ts_xxx_xxzz[i] * gfe_0 + ts_xxx_xxzzz[i] * ga_z[i];

        ts_xxxz_xyyyy[i] = ts_xxx_xyyyy[i] * ga_z[i];

        ts_xxxz_xyyyz[i] = ts_xxx_xyyy[i] * gfe_0 + ts_xxx_xyyyz[i] * ga_z[i];

        ts_xxxz_xyyzz[i] = 2.0 * ts_xxx_xyyz[i] * gfe_0 + ts_xxx_xyyzz[i] * ga_z[i];

        ts_xxxz_xyzzz[i] = 3.0 * ts_xxx_xyzz[i] * gfe_0 + ts_xxx_xyzzz[i] * ga_z[i];

        ts_xxxz_xzzzz[i] = 4.0 * ts_xxx_xzzz[i] * gfe_0 + ts_xxx_xzzzz[i] * ga_z[i];

        ts_xxxz_yyyyy[i] = ts_xxx_yyyyy[i] * ga_z[i];

        ts_xxxz_yyyyz[i] = 2.0 * ts_xz_yyyyz[i] * gfe_0 + ts_xxz_yyyyz[i] * ga_x[i];

        ts_xxxz_yyyzz[i] = 2.0 * ts_xz_yyyzz[i] * gfe_0 + ts_xxz_yyyzz[i] * ga_x[i];

        ts_xxxz_yyzzz[i] = 2.0 * ts_xz_yyzzz[i] * gfe_0 + ts_xxz_yyzzz[i] * ga_x[i];

        ts_xxxz_yzzzz[i] = 2.0 * ts_xz_yzzzz[i] * gfe_0 + ts_xxz_yzzzz[i] * ga_x[i];

        ts_xxxz_zzzzz[i] = 2.0 * ts_xz_zzzzz[i] * gfe_0 + ts_xxz_zzzzz[i] * ga_x[i];
    }

    // Set up 63-84 components of targeted buffer : GH

    auto ts_xxyy_xxxxx = pbuffer.data(idx_gh + 63);

    auto ts_xxyy_xxxxy = pbuffer.data(idx_gh + 64);

    auto ts_xxyy_xxxxz = pbuffer.data(idx_gh + 65);

    auto ts_xxyy_xxxyy = pbuffer.data(idx_gh + 66);

    auto ts_xxyy_xxxyz = pbuffer.data(idx_gh + 67);

    auto ts_xxyy_xxxzz = pbuffer.data(idx_gh + 68);

    auto ts_xxyy_xxyyy = pbuffer.data(idx_gh + 69);

    auto ts_xxyy_xxyyz = pbuffer.data(idx_gh + 70);

    auto ts_xxyy_xxyzz = pbuffer.data(idx_gh + 71);

    auto ts_xxyy_xxzzz = pbuffer.data(idx_gh + 72);

    auto ts_xxyy_xyyyy = pbuffer.data(idx_gh + 73);

    auto ts_xxyy_xyyyz = pbuffer.data(idx_gh + 74);

    auto ts_xxyy_xyyzz = pbuffer.data(idx_gh + 75);

    auto ts_xxyy_xyzzz = pbuffer.data(idx_gh + 76);

    auto ts_xxyy_xzzzz = pbuffer.data(idx_gh + 77);

    auto ts_xxyy_yyyyy = pbuffer.data(idx_gh + 78);

    auto ts_xxyy_yyyyz = pbuffer.data(idx_gh + 79);

    auto ts_xxyy_yyyzz = pbuffer.data(idx_gh + 80);

    auto ts_xxyy_yyzzz = pbuffer.data(idx_gh + 81);

    auto ts_xxyy_yzzzz = pbuffer.data(idx_gh + 82);

    auto ts_xxyy_zzzzz = pbuffer.data(idx_gh + 83);

    #pragma omp simd aligned(ga_x, ga_y, ts_xx_xxxxx, ts_xx_xxxxz, ts_xx_xxxzz, ts_xx_xxzzz, ts_xx_xzzzz, ts_xxy_xxxxx, ts_xxy_xxxxz, ts_xxy_xxxzz, ts_xxy_xxzzz, ts_xxy_xzzzz, ts_xxyy_xxxxx, ts_xxyy_xxxxy, ts_xxyy_xxxxz, ts_xxyy_xxxyy, ts_xxyy_xxxyz, ts_xxyy_xxxzz, ts_xxyy_xxyyy, ts_xxyy_xxyyz, ts_xxyy_xxyzz, ts_xxyy_xxzzz, ts_xxyy_xyyyy, ts_xxyy_xyyyz, ts_xxyy_xyyzz, ts_xxyy_xyzzz, ts_xxyy_xzzzz, ts_xxyy_yyyyy, ts_xxyy_yyyyz, ts_xxyy_yyyzz, ts_xxyy_yyzzz, ts_xxyy_yzzzz, ts_xxyy_zzzzz, ts_xyy_xxxxy, ts_xyy_xxxy, ts_xyy_xxxyy, ts_xyy_xxxyz, ts_xyy_xxyy, ts_xyy_xxyyy, ts_xyy_xxyyz, ts_xyy_xxyz, ts_xyy_xxyzz, ts_xyy_xyyy, ts_xyy_xyyyy, ts_xyy_xyyyz, ts_xyy_xyyz, ts_xyy_xyyzz, ts_xyy_xyzz, ts_xyy_xyzzz, ts_xyy_yyyy, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyz, ts_xyy_yyyzz, ts_xyy_yyzz, ts_xyy_yyzzz, ts_xyy_yzzz, ts_xyy_yzzzz, ts_xyy_zzzzz, ts_yy_xxxxy, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyzz, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyzz, ts_yy_xyzzz, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyzz, ts_yy_yyzzz, ts_yy_yzzzz, ts_yy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyy_xxxxx[i] = ts_xx_xxxxx[i] * gfe_0 + ts_xxy_xxxxx[i] * ga_y[i];

        ts_xxyy_xxxxy[i] = ts_yy_xxxxy[i] * gfe_0 + 4.0 * ts_xyy_xxxy[i] * gfe_0 + ts_xyy_xxxxy[i] * ga_x[i];

        ts_xxyy_xxxxz[i] = ts_xx_xxxxz[i] * gfe_0 + ts_xxy_xxxxz[i] * ga_y[i];

        ts_xxyy_xxxyy[i] = ts_yy_xxxyy[i] * gfe_0 + 3.0 * ts_xyy_xxyy[i] * gfe_0 + ts_xyy_xxxyy[i] * ga_x[i];

        ts_xxyy_xxxyz[i] = ts_yy_xxxyz[i] * gfe_0 + 3.0 * ts_xyy_xxyz[i] * gfe_0 + ts_xyy_xxxyz[i] * ga_x[i];

        ts_xxyy_xxxzz[i] = ts_xx_xxxzz[i] * gfe_0 + ts_xxy_xxxzz[i] * ga_y[i];

        ts_xxyy_xxyyy[i] = ts_yy_xxyyy[i] * gfe_0 + 2.0 * ts_xyy_xyyy[i] * gfe_0 + ts_xyy_xxyyy[i] * ga_x[i];

        ts_xxyy_xxyyz[i] = ts_yy_xxyyz[i] * gfe_0 + 2.0 * ts_xyy_xyyz[i] * gfe_0 + ts_xyy_xxyyz[i] * ga_x[i];

        ts_xxyy_xxyzz[i] = ts_yy_xxyzz[i] * gfe_0 + 2.0 * ts_xyy_xyzz[i] * gfe_0 + ts_xyy_xxyzz[i] * ga_x[i];

        ts_xxyy_xxzzz[i] = ts_xx_xxzzz[i] * gfe_0 + ts_xxy_xxzzz[i] * ga_y[i];

        ts_xxyy_xyyyy[i] = ts_yy_xyyyy[i] * gfe_0 + ts_xyy_yyyy[i] * gfe_0 + ts_xyy_xyyyy[i] * ga_x[i];

        ts_xxyy_xyyyz[i] = ts_yy_xyyyz[i] * gfe_0 + ts_xyy_yyyz[i] * gfe_0 + ts_xyy_xyyyz[i] * ga_x[i];

        ts_xxyy_xyyzz[i] = ts_yy_xyyzz[i] * gfe_0 + ts_xyy_yyzz[i] * gfe_0 + ts_xyy_xyyzz[i] * ga_x[i];

        ts_xxyy_xyzzz[i] = ts_yy_xyzzz[i] * gfe_0 + ts_xyy_yzzz[i] * gfe_0 + ts_xyy_xyzzz[i] * ga_x[i];

        ts_xxyy_xzzzz[i] = ts_xx_xzzzz[i] * gfe_0 + ts_xxy_xzzzz[i] * ga_y[i];

        ts_xxyy_yyyyy[i] = ts_yy_yyyyy[i] * gfe_0 + ts_xyy_yyyyy[i] * ga_x[i];

        ts_xxyy_yyyyz[i] = ts_yy_yyyyz[i] * gfe_0 + ts_xyy_yyyyz[i] * ga_x[i];

        ts_xxyy_yyyzz[i] = ts_yy_yyyzz[i] * gfe_0 + ts_xyy_yyyzz[i] * ga_x[i];

        ts_xxyy_yyzzz[i] = ts_yy_yyzzz[i] * gfe_0 + ts_xyy_yyzzz[i] * ga_x[i];

        ts_xxyy_yzzzz[i] = ts_yy_yzzzz[i] * gfe_0 + ts_xyy_yzzzz[i] * ga_x[i];

        ts_xxyy_zzzzz[i] = ts_yy_zzzzz[i] * gfe_0 + ts_xyy_zzzzz[i] * ga_x[i];
    }

    // Set up 84-105 components of targeted buffer : GH

    auto ts_xxyz_xxxxx = pbuffer.data(idx_gh + 84);

    auto ts_xxyz_xxxxy = pbuffer.data(idx_gh + 85);

    auto ts_xxyz_xxxxz = pbuffer.data(idx_gh + 86);

    auto ts_xxyz_xxxyy = pbuffer.data(idx_gh + 87);

    auto ts_xxyz_xxxyz = pbuffer.data(idx_gh + 88);

    auto ts_xxyz_xxxzz = pbuffer.data(idx_gh + 89);

    auto ts_xxyz_xxyyy = pbuffer.data(idx_gh + 90);

    auto ts_xxyz_xxyyz = pbuffer.data(idx_gh + 91);

    auto ts_xxyz_xxyzz = pbuffer.data(idx_gh + 92);

    auto ts_xxyz_xxzzz = pbuffer.data(idx_gh + 93);

    auto ts_xxyz_xyyyy = pbuffer.data(idx_gh + 94);

    auto ts_xxyz_xyyyz = pbuffer.data(idx_gh + 95);

    auto ts_xxyz_xyyzz = pbuffer.data(idx_gh + 96);

    auto ts_xxyz_xyzzz = pbuffer.data(idx_gh + 97);

    auto ts_xxyz_xzzzz = pbuffer.data(idx_gh + 98);

    auto ts_xxyz_yyyyy = pbuffer.data(idx_gh + 99);

    auto ts_xxyz_yyyyz = pbuffer.data(idx_gh + 100);

    auto ts_xxyz_yyyzz = pbuffer.data(idx_gh + 101);

    auto ts_xxyz_yyzzz = pbuffer.data(idx_gh + 102);

    auto ts_xxyz_yzzzz = pbuffer.data(idx_gh + 103);

    auto ts_xxyz_zzzzz = pbuffer.data(idx_gh + 104);

    #pragma omp simd aligned(ga_x, ga_y, ga_z, ts_xxy_xxxxy, ts_xxy_xxxyy, ts_xxy_xxyyy, ts_xxy_xyyyy, ts_xxy_yyyyy, ts_xxyz_xxxxx, ts_xxyz_xxxxy, ts_xxyz_xxxxz, ts_xxyz_xxxyy, ts_xxyz_xxxyz, ts_xxyz_xxxzz, ts_xxyz_xxyyy, ts_xxyz_xxyyz, ts_xxyz_xxyzz, ts_xxyz_xxzzz, ts_xxyz_xyyyy, ts_xxyz_xyyyz, ts_xxyz_xyyzz, ts_xxyz_xyzzz, ts_xxyz_xzzzz, ts_xxyz_yyyyy, ts_xxyz_yyyyz, ts_xxyz_yyyzz, ts_xxyz_yyzzz, ts_xxyz_yzzzz, ts_xxyz_zzzzz, ts_xxz_xxxxx, ts_xxz_xxxxz, ts_xxz_xxxyz, ts_xxz_xxxz, ts_xxz_xxxzz, ts_xxz_xxyyz, ts_xxz_xxyz, ts_xxz_xxyzz, ts_xxz_xxzz, ts_xxz_xxzzz, ts_xxz_xyyyz, ts_xxz_xyyz, ts_xxz_xyyzz, ts_xxz_xyzz, ts_xxz_xyzzz, ts_xxz_xzzz, ts_xxz_xzzzz, ts_xxz_zzzzz, ts_xyz_yyyyz, ts_xyz_yyyzz, ts_xyz_yyzzz, ts_xyz_yzzzz, ts_yz_yyyyz, ts_yz_yyyzz, ts_yz_yyzzz, ts_yz_yzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyz_xxxxx[i] = ts_xxz_xxxxx[i] * ga_y[i];

        ts_xxyz_xxxxy[i] = ts_xxy_xxxxy[i] * ga_z[i];

        ts_xxyz_xxxxz[i] = ts_xxz_xxxxz[i] * ga_y[i];

        ts_xxyz_xxxyy[i] = ts_xxy_xxxyy[i] * ga_z[i];

        ts_xxyz_xxxyz[i] = ts_xxz_xxxz[i] * gfe_0 + ts_xxz_xxxyz[i] * ga_y[i];

        ts_xxyz_xxxzz[i] = ts_xxz_xxxzz[i] * ga_y[i];

        ts_xxyz_xxyyy[i] = ts_xxy_xxyyy[i] * ga_z[i];

        ts_xxyz_xxyyz[i] = 2.0 * ts_xxz_xxyz[i] * gfe_0 + ts_xxz_xxyyz[i] * ga_y[i];

        ts_xxyz_xxyzz[i] = ts_xxz_xxzz[i] * gfe_0 + ts_xxz_xxyzz[i] * ga_y[i];

        ts_xxyz_xxzzz[i] = ts_xxz_xxzzz[i] * ga_y[i];

        ts_xxyz_xyyyy[i] = ts_xxy_xyyyy[i] * ga_z[i];

        ts_xxyz_xyyyz[i] = 3.0 * ts_xxz_xyyz[i] * gfe_0 + ts_xxz_xyyyz[i] * ga_y[i];

        ts_xxyz_xyyzz[i] = 2.0 * ts_xxz_xyzz[i] * gfe_0 + ts_xxz_xyyzz[i] * ga_y[i];

        ts_xxyz_xyzzz[i] = ts_xxz_xzzz[i] * gfe_0 + ts_xxz_xyzzz[i] * ga_y[i];

        ts_xxyz_xzzzz[i] = ts_xxz_xzzzz[i] * ga_y[i];

        ts_xxyz_yyyyy[i] = ts_xxy_yyyyy[i] * ga_z[i];

        ts_xxyz_yyyyz[i] = ts_yz_yyyyz[i] * gfe_0 + ts_xyz_yyyyz[i] * ga_x[i];

        ts_xxyz_yyyzz[i] = ts_yz_yyyzz[i] * gfe_0 + ts_xyz_yyyzz[i] * ga_x[i];

        ts_xxyz_yyzzz[i] = ts_yz_yyzzz[i] * gfe_0 + ts_xyz_yyzzz[i] * ga_x[i];

        ts_xxyz_yzzzz[i] = ts_yz_yzzzz[i] * gfe_0 + ts_xyz_yzzzz[i] * ga_x[i];

        ts_xxyz_zzzzz[i] = ts_xxz_zzzzz[i] * ga_y[i];
    }

    // Set up 105-126 components of targeted buffer : GH

    auto ts_xxzz_xxxxx = pbuffer.data(idx_gh + 105);

    auto ts_xxzz_xxxxy = pbuffer.data(idx_gh + 106);

    auto ts_xxzz_xxxxz = pbuffer.data(idx_gh + 107);

    auto ts_xxzz_xxxyy = pbuffer.data(idx_gh + 108);

    auto ts_xxzz_xxxyz = pbuffer.data(idx_gh + 109);

    auto ts_xxzz_xxxzz = pbuffer.data(idx_gh + 110);

    auto ts_xxzz_xxyyy = pbuffer.data(idx_gh + 111);

    auto ts_xxzz_xxyyz = pbuffer.data(idx_gh + 112);

    auto ts_xxzz_xxyzz = pbuffer.data(idx_gh + 113);

    auto ts_xxzz_xxzzz = pbuffer.data(idx_gh + 114);

    auto ts_xxzz_xyyyy = pbuffer.data(idx_gh + 115);

    auto ts_xxzz_xyyyz = pbuffer.data(idx_gh + 116);

    auto ts_xxzz_xyyzz = pbuffer.data(idx_gh + 117);

    auto ts_xxzz_xyzzz = pbuffer.data(idx_gh + 118);

    auto ts_xxzz_xzzzz = pbuffer.data(idx_gh + 119);

    auto ts_xxzz_yyyyy = pbuffer.data(idx_gh + 120);

    auto ts_xxzz_yyyyz = pbuffer.data(idx_gh + 121);

    auto ts_xxzz_yyyzz = pbuffer.data(idx_gh + 122);

    auto ts_xxzz_yyzzz = pbuffer.data(idx_gh + 123);

    auto ts_xxzz_yzzzz = pbuffer.data(idx_gh + 124);

    auto ts_xxzz_zzzzz = pbuffer.data(idx_gh + 125);

    #pragma omp simd aligned(ga_x, ga_z, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxyy, ts_xx_xxyyy, ts_xx_xyyyy, ts_xxz_xxxxx, ts_xxz_xxxxy, ts_xxz_xxxyy, ts_xxz_xxyyy, ts_xxz_xyyyy, ts_xxzz_xxxxx, ts_xxzz_xxxxy, ts_xxzz_xxxxz, ts_xxzz_xxxyy, ts_xxzz_xxxyz, ts_xxzz_xxxzz, ts_xxzz_xxyyy, ts_xxzz_xxyyz, ts_xxzz_xxyzz, ts_xxzz_xxzzz, ts_xxzz_xyyyy, ts_xxzz_xyyyz, ts_xxzz_xyyzz, ts_xxzz_xyzzz, ts_xxzz_xzzzz, ts_xxzz_yyyyy, ts_xxzz_yyyyz, ts_xxzz_yyyzz, ts_xxzz_yyzzz, ts_xxzz_yzzzz, ts_xxzz_zzzzz, ts_xzz_xxxxz, ts_xzz_xxxyz, ts_xzz_xxxz, ts_xzz_xxxzz, ts_xzz_xxyyz, ts_xzz_xxyz, ts_xzz_xxyzz, ts_xzz_xxzz, ts_xzz_xxzzz, ts_xzz_xyyyz, ts_xzz_xyyz, ts_xzz_xyyzz, ts_xzz_xyzz, ts_xzz_xyzzz, ts_xzz_xzzz, ts_xzz_xzzzz, ts_xzz_yyyyy, ts_xzz_yyyyz, ts_xzz_yyyz, ts_xzz_yyyzz, ts_xzz_yyzz, ts_xzz_yyzzz, ts_xzz_yzzz, ts_xzz_yzzzz, ts_xzz_zzzz, ts_xzz_zzzzz, ts_zz_xxxxz, ts_zz_xxxyz, ts_zz_xxxzz, ts_zz_xxyyz, ts_zz_xxyzz, ts_zz_xxzzz, ts_zz_xyyyz, ts_zz_xyyzz, ts_zz_xyzzz, ts_zz_xzzzz, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyzz, ts_zz_yyzzz, ts_zz_yzzzz, ts_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxzz_xxxxx[i] = ts_xx_xxxxx[i] * gfe_0 + ts_xxz_xxxxx[i] * ga_z[i];

        ts_xxzz_xxxxy[i] = ts_xx_xxxxy[i] * gfe_0 + ts_xxz_xxxxy[i] * ga_z[i];

        ts_xxzz_xxxxz[i] = ts_zz_xxxxz[i] * gfe_0 + 4.0 * ts_xzz_xxxz[i] * gfe_0 + ts_xzz_xxxxz[i] * ga_x[i];

        ts_xxzz_xxxyy[i] = ts_xx_xxxyy[i] * gfe_0 + ts_xxz_xxxyy[i] * ga_z[i];

        ts_xxzz_xxxyz[i] = ts_zz_xxxyz[i] * gfe_0 + 3.0 * ts_xzz_xxyz[i] * gfe_0 + ts_xzz_xxxyz[i] * ga_x[i];

        ts_xxzz_xxxzz[i] = ts_zz_xxxzz[i] * gfe_0 + 3.0 * ts_xzz_xxzz[i] * gfe_0 + ts_xzz_xxxzz[i] * ga_x[i];

        ts_xxzz_xxyyy[i] = ts_xx_xxyyy[i] * gfe_0 + ts_xxz_xxyyy[i] * ga_z[i];

        ts_xxzz_xxyyz[i] = ts_zz_xxyyz[i] * gfe_0 + 2.0 * ts_xzz_xyyz[i] * gfe_0 + ts_xzz_xxyyz[i] * ga_x[i];

        ts_xxzz_xxyzz[i] = ts_zz_xxyzz[i] * gfe_0 + 2.0 * ts_xzz_xyzz[i] * gfe_0 + ts_xzz_xxyzz[i] * ga_x[i];

        ts_xxzz_xxzzz[i] = ts_zz_xxzzz[i] * gfe_0 + 2.0 * ts_xzz_xzzz[i] * gfe_0 + ts_xzz_xxzzz[i] * ga_x[i];

        ts_xxzz_xyyyy[i] = ts_xx_xyyyy[i] * gfe_0 + ts_xxz_xyyyy[i] * ga_z[i];

        ts_xxzz_xyyyz[i] = ts_zz_xyyyz[i] * gfe_0 + ts_xzz_yyyz[i] * gfe_0 + ts_xzz_xyyyz[i] * ga_x[i];

        ts_xxzz_xyyzz[i] = ts_zz_xyyzz[i] * gfe_0 + ts_xzz_yyzz[i] * gfe_0 + ts_xzz_xyyzz[i] * ga_x[i];

        ts_xxzz_xyzzz[i] = ts_zz_xyzzz[i] * gfe_0 + ts_xzz_yzzz[i] * gfe_0 + ts_xzz_xyzzz[i] * ga_x[i];

        ts_xxzz_xzzzz[i] = ts_zz_xzzzz[i] * gfe_0 + ts_xzz_zzzz[i] * gfe_0 + ts_xzz_xzzzz[i] * ga_x[i];

        ts_xxzz_yyyyy[i] = ts_zz_yyyyy[i] * gfe_0 + ts_xzz_yyyyy[i] * ga_x[i];

        ts_xxzz_yyyyz[i] = ts_zz_yyyyz[i] * gfe_0 + ts_xzz_yyyyz[i] * ga_x[i];

        ts_xxzz_yyyzz[i] = ts_zz_yyyzz[i] * gfe_0 + ts_xzz_yyyzz[i] * ga_x[i];

        ts_xxzz_yyzzz[i] = ts_zz_yyzzz[i] * gfe_0 + ts_xzz_yyzzz[i] * ga_x[i];

        ts_xxzz_yzzzz[i] = ts_zz_yzzzz[i] * gfe_0 + ts_xzz_yzzzz[i] * ga_x[i];

        ts_xxzz_zzzzz[i] = ts_zz_zzzzz[i] * gfe_0 + ts_xzz_zzzzz[i] * ga_x[i];
    }

    // Set up 126-147 components of targeted buffer : GH

    auto ts_xyyy_xxxxx = pbuffer.data(idx_gh + 126);

    auto ts_xyyy_xxxxy = pbuffer.data(idx_gh + 127);

    auto ts_xyyy_xxxxz = pbuffer.data(idx_gh + 128);

    auto ts_xyyy_xxxyy = pbuffer.data(idx_gh + 129);

    auto ts_xyyy_xxxyz = pbuffer.data(idx_gh + 130);

    auto ts_xyyy_xxxzz = pbuffer.data(idx_gh + 131);

    auto ts_xyyy_xxyyy = pbuffer.data(idx_gh + 132);

    auto ts_xyyy_xxyyz = pbuffer.data(idx_gh + 133);

    auto ts_xyyy_xxyzz = pbuffer.data(idx_gh + 134);

    auto ts_xyyy_xxzzz = pbuffer.data(idx_gh + 135);

    auto ts_xyyy_xyyyy = pbuffer.data(idx_gh + 136);

    auto ts_xyyy_xyyyz = pbuffer.data(idx_gh + 137);

    auto ts_xyyy_xyyzz = pbuffer.data(idx_gh + 138);

    auto ts_xyyy_xyzzz = pbuffer.data(idx_gh + 139);

    auto ts_xyyy_xzzzz = pbuffer.data(idx_gh + 140);

    auto ts_xyyy_yyyyy = pbuffer.data(idx_gh + 141);

    auto ts_xyyy_yyyyz = pbuffer.data(idx_gh + 142);

    auto ts_xyyy_yyyzz = pbuffer.data(idx_gh + 143);

    auto ts_xyyy_yyzzz = pbuffer.data(idx_gh + 144);

    auto ts_xyyy_yzzzz = pbuffer.data(idx_gh + 145);

    auto ts_xyyy_zzzzz = pbuffer.data(idx_gh + 146);

    #pragma omp simd aligned(ga_x, ts_xyyy_xxxxx, ts_xyyy_xxxxy, ts_xyyy_xxxxz, ts_xyyy_xxxyy, ts_xyyy_xxxyz, ts_xyyy_xxxzz, ts_xyyy_xxyyy, ts_xyyy_xxyyz, ts_xyyy_xxyzz, ts_xyyy_xxzzz, ts_xyyy_xyyyy, ts_xyyy_xyyyz, ts_xyyy_xyyzz, ts_xyyy_xyzzz, ts_xyyy_xzzzz, ts_xyyy_yyyyy, ts_xyyy_yyyyz, ts_xyyy_yyyzz, ts_xyyy_yyzzz, ts_xyyy_yzzzz, ts_xyyy_zzzzz, ts_yyy_xxxx, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxy, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxz, ts_yyy_xxxzz, ts_yyy_xxyy, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyz, ts_yyy_xxyzz, ts_yyy_xxzz, ts_yyy_xxzzz, ts_yyy_xyyy, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyz, ts_yyy_xyyzz, ts_yyy_xyzz, ts_yyy_xyzzz, ts_yyy_xzzz, ts_yyy_xzzzz, ts_yyy_yyyy, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyz, ts_yyy_yyyzz, ts_yyy_yyzz, ts_yyy_yyzzz, ts_yyy_yzzz, ts_yyy_yzzzz, ts_yyy_zzzz, ts_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyy_xxxxx[i] = 5.0 * ts_yyy_xxxx[i] * gfe_0 + ts_yyy_xxxxx[i] * ga_x[i];

        ts_xyyy_xxxxy[i] = 4.0 * ts_yyy_xxxy[i] * gfe_0 + ts_yyy_xxxxy[i] * ga_x[i];

        ts_xyyy_xxxxz[i] = 4.0 * ts_yyy_xxxz[i] * gfe_0 + ts_yyy_xxxxz[i] * ga_x[i];

        ts_xyyy_xxxyy[i] = 3.0 * ts_yyy_xxyy[i] * gfe_0 + ts_yyy_xxxyy[i] * ga_x[i];

        ts_xyyy_xxxyz[i] = 3.0 * ts_yyy_xxyz[i] * gfe_0 + ts_yyy_xxxyz[i] * ga_x[i];

        ts_xyyy_xxxzz[i] = 3.0 * ts_yyy_xxzz[i] * gfe_0 + ts_yyy_xxxzz[i] * ga_x[i];

        ts_xyyy_xxyyy[i] = 2.0 * ts_yyy_xyyy[i] * gfe_0 + ts_yyy_xxyyy[i] * ga_x[i];

        ts_xyyy_xxyyz[i] = 2.0 * ts_yyy_xyyz[i] * gfe_0 + ts_yyy_xxyyz[i] * ga_x[i];

        ts_xyyy_xxyzz[i] = 2.0 * ts_yyy_xyzz[i] * gfe_0 + ts_yyy_xxyzz[i] * ga_x[i];

        ts_xyyy_xxzzz[i] = 2.0 * ts_yyy_xzzz[i] * gfe_0 + ts_yyy_xxzzz[i] * ga_x[i];

        ts_xyyy_xyyyy[i] = ts_yyy_yyyy[i] * gfe_0 + ts_yyy_xyyyy[i] * ga_x[i];

        ts_xyyy_xyyyz[i] = ts_yyy_yyyz[i] * gfe_0 + ts_yyy_xyyyz[i] * ga_x[i];

        ts_xyyy_xyyzz[i] = ts_yyy_yyzz[i] * gfe_0 + ts_yyy_xyyzz[i] * ga_x[i];

        ts_xyyy_xyzzz[i] = ts_yyy_yzzz[i] * gfe_0 + ts_yyy_xyzzz[i] * ga_x[i];

        ts_xyyy_xzzzz[i] = ts_yyy_zzzz[i] * gfe_0 + ts_yyy_xzzzz[i] * ga_x[i];

        ts_xyyy_yyyyy[i] = ts_yyy_yyyyy[i] * ga_x[i];

        ts_xyyy_yyyyz[i] = ts_yyy_yyyyz[i] * ga_x[i];

        ts_xyyy_yyyzz[i] = ts_yyy_yyyzz[i] * ga_x[i];

        ts_xyyy_yyzzz[i] = ts_yyy_yyzzz[i] * ga_x[i];

        ts_xyyy_yzzzz[i] = ts_yyy_yzzzz[i] * ga_x[i];

        ts_xyyy_zzzzz[i] = ts_yyy_zzzzz[i] * ga_x[i];
    }

    // Set up 147-168 components of targeted buffer : GH

    auto ts_xyyz_xxxxx = pbuffer.data(idx_gh + 147);

    auto ts_xyyz_xxxxy = pbuffer.data(idx_gh + 148);

    auto ts_xyyz_xxxxz = pbuffer.data(idx_gh + 149);

    auto ts_xyyz_xxxyy = pbuffer.data(idx_gh + 150);

    auto ts_xyyz_xxxyz = pbuffer.data(idx_gh + 151);

    auto ts_xyyz_xxxzz = pbuffer.data(idx_gh + 152);

    auto ts_xyyz_xxyyy = pbuffer.data(idx_gh + 153);

    auto ts_xyyz_xxyyz = pbuffer.data(idx_gh + 154);

    auto ts_xyyz_xxyzz = pbuffer.data(idx_gh + 155);

    auto ts_xyyz_xxzzz = pbuffer.data(idx_gh + 156);

    auto ts_xyyz_xyyyy = pbuffer.data(idx_gh + 157);

    auto ts_xyyz_xyyyz = pbuffer.data(idx_gh + 158);

    auto ts_xyyz_xyyzz = pbuffer.data(idx_gh + 159);

    auto ts_xyyz_xyzzz = pbuffer.data(idx_gh + 160);

    auto ts_xyyz_xzzzz = pbuffer.data(idx_gh + 161);

    auto ts_xyyz_yyyyy = pbuffer.data(idx_gh + 162);

    auto ts_xyyz_yyyyz = pbuffer.data(idx_gh + 163);

    auto ts_xyyz_yyyzz = pbuffer.data(idx_gh + 164);

    auto ts_xyyz_yyzzz = pbuffer.data(idx_gh + 165);

    auto ts_xyyz_yzzzz = pbuffer.data(idx_gh + 166);

    auto ts_xyyz_zzzzz = pbuffer.data(idx_gh + 167);

    #pragma omp simd aligned(ga_x, ga_z, ts_xyy_xxxxx, ts_xyy_xxxxy, ts_xyy_xxxyy, ts_xyy_xxyyy, ts_xyy_xyyyy, ts_xyyz_xxxxx, ts_xyyz_xxxxy, ts_xyyz_xxxxz, ts_xyyz_xxxyy, ts_xyyz_xxxyz, ts_xyyz_xxxzz, ts_xyyz_xxyyy, ts_xyyz_xxyyz, ts_xyyz_xxyzz, ts_xyyz_xxzzz, ts_xyyz_xyyyy, ts_xyyz_xyyyz, ts_xyyz_xyyzz, ts_xyyz_xyzzz, ts_xyyz_xzzzz, ts_xyyz_yyyyy, ts_xyyz_yyyyz, ts_xyyz_yyyzz, ts_xyyz_yyzzz, ts_xyyz_yzzzz, ts_xyyz_zzzzz, ts_yyz_xxxxz, ts_yyz_xxxyz, ts_yyz_xxxz, ts_yyz_xxxzz, ts_yyz_xxyyz, ts_yyz_xxyz, ts_yyz_xxyzz, ts_yyz_xxzz, ts_yyz_xxzzz, ts_yyz_xyyyz, ts_yyz_xyyz, ts_yyz_xyyzz, ts_yyz_xyzz, ts_yyz_xyzzz, ts_yyz_xzzz, ts_yyz_xzzzz, ts_yyz_yyyyy, ts_yyz_yyyyz, ts_yyz_yyyz, ts_yyz_yyyzz, ts_yyz_yyzz, ts_yyz_yyzzz, ts_yyz_yzzz, ts_yyz_yzzzz, ts_yyz_zzzz, ts_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyz_xxxxx[i] = ts_xyy_xxxxx[i] * ga_z[i];

        ts_xyyz_xxxxy[i] = ts_xyy_xxxxy[i] * ga_z[i];

        ts_xyyz_xxxxz[i] = 4.0 * ts_yyz_xxxz[i] * gfe_0 + ts_yyz_xxxxz[i] * ga_x[i];

        ts_xyyz_xxxyy[i] = ts_xyy_xxxyy[i] * ga_z[i];

        ts_xyyz_xxxyz[i] = 3.0 * ts_yyz_xxyz[i] * gfe_0 + ts_yyz_xxxyz[i] * ga_x[i];

        ts_xyyz_xxxzz[i] = 3.0 * ts_yyz_xxzz[i] * gfe_0 + ts_yyz_xxxzz[i] * ga_x[i];

        ts_xyyz_xxyyy[i] = ts_xyy_xxyyy[i] * ga_z[i];

        ts_xyyz_xxyyz[i] = 2.0 * ts_yyz_xyyz[i] * gfe_0 + ts_yyz_xxyyz[i] * ga_x[i];

        ts_xyyz_xxyzz[i] = 2.0 * ts_yyz_xyzz[i] * gfe_0 + ts_yyz_xxyzz[i] * ga_x[i];

        ts_xyyz_xxzzz[i] = 2.0 * ts_yyz_xzzz[i] * gfe_0 + ts_yyz_xxzzz[i] * ga_x[i];

        ts_xyyz_xyyyy[i] = ts_xyy_xyyyy[i] * ga_z[i];

        ts_xyyz_xyyyz[i] = ts_yyz_yyyz[i] * gfe_0 + ts_yyz_xyyyz[i] * ga_x[i];

        ts_xyyz_xyyzz[i] = ts_yyz_yyzz[i] * gfe_0 + ts_yyz_xyyzz[i] * ga_x[i];

        ts_xyyz_xyzzz[i] = ts_yyz_yzzz[i] * gfe_0 + ts_yyz_xyzzz[i] * ga_x[i];

        ts_xyyz_xzzzz[i] = ts_yyz_zzzz[i] * gfe_0 + ts_yyz_xzzzz[i] * ga_x[i];

        ts_xyyz_yyyyy[i] = ts_yyz_yyyyy[i] * ga_x[i];

        ts_xyyz_yyyyz[i] = ts_yyz_yyyyz[i] * ga_x[i];

        ts_xyyz_yyyzz[i] = ts_yyz_yyyzz[i] * ga_x[i];

        ts_xyyz_yyzzz[i] = ts_yyz_yyzzz[i] * ga_x[i];

        ts_xyyz_yzzzz[i] = ts_yyz_yzzzz[i] * ga_x[i];

        ts_xyyz_zzzzz[i] = ts_yyz_zzzzz[i] * ga_x[i];
    }

    // Set up 168-189 components of targeted buffer : GH

    auto ts_xyzz_xxxxx = pbuffer.data(idx_gh + 168);

    auto ts_xyzz_xxxxy = pbuffer.data(idx_gh + 169);

    auto ts_xyzz_xxxxz = pbuffer.data(idx_gh + 170);

    auto ts_xyzz_xxxyy = pbuffer.data(idx_gh + 171);

    auto ts_xyzz_xxxyz = pbuffer.data(idx_gh + 172);

    auto ts_xyzz_xxxzz = pbuffer.data(idx_gh + 173);

    auto ts_xyzz_xxyyy = pbuffer.data(idx_gh + 174);

    auto ts_xyzz_xxyyz = pbuffer.data(idx_gh + 175);

    auto ts_xyzz_xxyzz = pbuffer.data(idx_gh + 176);

    auto ts_xyzz_xxzzz = pbuffer.data(idx_gh + 177);

    auto ts_xyzz_xyyyy = pbuffer.data(idx_gh + 178);

    auto ts_xyzz_xyyyz = pbuffer.data(idx_gh + 179);

    auto ts_xyzz_xyyzz = pbuffer.data(idx_gh + 180);

    auto ts_xyzz_xyzzz = pbuffer.data(idx_gh + 181);

    auto ts_xyzz_xzzzz = pbuffer.data(idx_gh + 182);

    auto ts_xyzz_yyyyy = pbuffer.data(idx_gh + 183);

    auto ts_xyzz_yyyyz = pbuffer.data(idx_gh + 184);

    auto ts_xyzz_yyyzz = pbuffer.data(idx_gh + 185);

    auto ts_xyzz_yyzzz = pbuffer.data(idx_gh + 186);

    auto ts_xyzz_yzzzz = pbuffer.data(idx_gh + 187);

    auto ts_xyzz_zzzzz = pbuffer.data(idx_gh + 188);

    #pragma omp simd aligned(ga_x, ga_y, ts_xyzz_xxxxx, ts_xyzz_xxxxy, ts_xyzz_xxxxz, ts_xyzz_xxxyy, ts_xyzz_xxxyz, ts_xyzz_xxxzz, ts_xyzz_xxyyy, ts_xyzz_xxyyz, ts_xyzz_xxyzz, ts_xyzz_xxzzz, ts_xyzz_xyyyy, ts_xyzz_xyyyz, ts_xyzz_xyyzz, ts_xyzz_xyzzz, ts_xyzz_xzzzz, ts_xyzz_yyyyy, ts_xyzz_yyyyz, ts_xyzz_yyyzz, ts_xyzz_yyzzz, ts_xyzz_yzzzz, ts_xyzz_zzzzz, ts_xzz_xxxxx, ts_xzz_xxxxz, ts_xzz_xxxzz, ts_xzz_xxzzz, ts_xzz_xzzzz, ts_yzz_xxxxy, ts_yzz_xxxy, ts_yzz_xxxyy, ts_yzz_xxxyz, ts_yzz_xxyy, ts_yzz_xxyyy, ts_yzz_xxyyz, ts_yzz_xxyz, ts_yzz_xxyzz, ts_yzz_xyyy, ts_yzz_xyyyy, ts_yzz_xyyyz, ts_yzz_xyyz, ts_yzz_xyyzz, ts_yzz_xyzz, ts_yzz_xyzzz, ts_yzz_yyyy, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyz, ts_yzz_yyyzz, ts_yzz_yyzz, ts_yzz_yyzzz, ts_yzz_yzzz, ts_yzz_yzzzz, ts_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyzz_xxxxx[i] = ts_xzz_xxxxx[i] * ga_y[i];

        ts_xyzz_xxxxy[i] = 4.0 * ts_yzz_xxxy[i] * gfe_0 + ts_yzz_xxxxy[i] * ga_x[i];

        ts_xyzz_xxxxz[i] = ts_xzz_xxxxz[i] * ga_y[i];

        ts_xyzz_xxxyy[i] = 3.0 * ts_yzz_xxyy[i] * gfe_0 + ts_yzz_xxxyy[i] * ga_x[i];

        ts_xyzz_xxxyz[i] = 3.0 * ts_yzz_xxyz[i] * gfe_0 + ts_yzz_xxxyz[i] * ga_x[i];

        ts_xyzz_xxxzz[i] = ts_xzz_xxxzz[i] * ga_y[i];

        ts_xyzz_xxyyy[i] = 2.0 * ts_yzz_xyyy[i] * gfe_0 + ts_yzz_xxyyy[i] * ga_x[i];

        ts_xyzz_xxyyz[i] = 2.0 * ts_yzz_xyyz[i] * gfe_0 + ts_yzz_xxyyz[i] * ga_x[i];

        ts_xyzz_xxyzz[i] = 2.0 * ts_yzz_xyzz[i] * gfe_0 + ts_yzz_xxyzz[i] * ga_x[i];

        ts_xyzz_xxzzz[i] = ts_xzz_xxzzz[i] * ga_y[i];

        ts_xyzz_xyyyy[i] = ts_yzz_yyyy[i] * gfe_0 + ts_yzz_xyyyy[i] * ga_x[i];

        ts_xyzz_xyyyz[i] = ts_yzz_yyyz[i] * gfe_0 + ts_yzz_xyyyz[i] * ga_x[i];

        ts_xyzz_xyyzz[i] = ts_yzz_yyzz[i] * gfe_0 + ts_yzz_xyyzz[i] * ga_x[i];

        ts_xyzz_xyzzz[i] = ts_yzz_yzzz[i] * gfe_0 + ts_yzz_xyzzz[i] * ga_x[i];

        ts_xyzz_xzzzz[i] = ts_xzz_xzzzz[i] * ga_y[i];

        ts_xyzz_yyyyy[i] = ts_yzz_yyyyy[i] * ga_x[i];

        ts_xyzz_yyyyz[i] = ts_yzz_yyyyz[i] * ga_x[i];

        ts_xyzz_yyyzz[i] = ts_yzz_yyyzz[i] * ga_x[i];

        ts_xyzz_yyzzz[i] = ts_yzz_yyzzz[i] * ga_x[i];

        ts_xyzz_yzzzz[i] = ts_yzz_yzzzz[i] * ga_x[i];

        ts_xyzz_zzzzz[i] = ts_yzz_zzzzz[i] * ga_x[i];
    }

    // Set up 189-210 components of targeted buffer : GH

    auto ts_xzzz_xxxxx = pbuffer.data(idx_gh + 189);

    auto ts_xzzz_xxxxy = pbuffer.data(idx_gh + 190);

    auto ts_xzzz_xxxxz = pbuffer.data(idx_gh + 191);

    auto ts_xzzz_xxxyy = pbuffer.data(idx_gh + 192);

    auto ts_xzzz_xxxyz = pbuffer.data(idx_gh + 193);

    auto ts_xzzz_xxxzz = pbuffer.data(idx_gh + 194);

    auto ts_xzzz_xxyyy = pbuffer.data(idx_gh + 195);

    auto ts_xzzz_xxyyz = pbuffer.data(idx_gh + 196);

    auto ts_xzzz_xxyzz = pbuffer.data(idx_gh + 197);

    auto ts_xzzz_xxzzz = pbuffer.data(idx_gh + 198);

    auto ts_xzzz_xyyyy = pbuffer.data(idx_gh + 199);

    auto ts_xzzz_xyyyz = pbuffer.data(idx_gh + 200);

    auto ts_xzzz_xyyzz = pbuffer.data(idx_gh + 201);

    auto ts_xzzz_xyzzz = pbuffer.data(idx_gh + 202);

    auto ts_xzzz_xzzzz = pbuffer.data(idx_gh + 203);

    auto ts_xzzz_yyyyy = pbuffer.data(idx_gh + 204);

    auto ts_xzzz_yyyyz = pbuffer.data(idx_gh + 205);

    auto ts_xzzz_yyyzz = pbuffer.data(idx_gh + 206);

    auto ts_xzzz_yyzzz = pbuffer.data(idx_gh + 207);

    auto ts_xzzz_yzzzz = pbuffer.data(idx_gh + 208);

    auto ts_xzzz_zzzzz = pbuffer.data(idx_gh + 209);

    #pragma omp simd aligned(ga_x, ts_xzzz_xxxxx, ts_xzzz_xxxxy, ts_xzzz_xxxxz, ts_xzzz_xxxyy, ts_xzzz_xxxyz, ts_xzzz_xxxzz, ts_xzzz_xxyyy, ts_xzzz_xxyyz, ts_xzzz_xxyzz, ts_xzzz_xxzzz, ts_xzzz_xyyyy, ts_xzzz_xyyyz, ts_xzzz_xyyzz, ts_xzzz_xyzzz, ts_xzzz_xzzzz, ts_xzzz_yyyyy, ts_xzzz_yyyyz, ts_xzzz_yyyzz, ts_xzzz_yyzzz, ts_xzzz_yzzzz, ts_xzzz_zzzzz, ts_zzz_xxxx, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxy, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxz, ts_zzz_xxxzz, ts_zzz_xxyy, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyz, ts_zzz_xxyzz, ts_zzz_xxzz, ts_zzz_xxzzz, ts_zzz_xyyy, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyz, ts_zzz_xyyzz, ts_zzz_xyzz, ts_zzz_xyzzz, ts_zzz_xzzz, ts_zzz_xzzzz, ts_zzz_yyyy, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyz, ts_zzz_yyyzz, ts_zzz_yyzz, ts_zzz_yyzzz, ts_zzz_yzzz, ts_zzz_yzzzz, ts_zzz_zzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xzzz_xxxxx[i] = 5.0 * ts_zzz_xxxx[i] * gfe_0 + ts_zzz_xxxxx[i] * ga_x[i];

        ts_xzzz_xxxxy[i] = 4.0 * ts_zzz_xxxy[i] * gfe_0 + ts_zzz_xxxxy[i] * ga_x[i];

        ts_xzzz_xxxxz[i] = 4.0 * ts_zzz_xxxz[i] * gfe_0 + ts_zzz_xxxxz[i] * ga_x[i];

        ts_xzzz_xxxyy[i] = 3.0 * ts_zzz_xxyy[i] * gfe_0 + ts_zzz_xxxyy[i] * ga_x[i];

        ts_xzzz_xxxyz[i] = 3.0 * ts_zzz_xxyz[i] * gfe_0 + ts_zzz_xxxyz[i] * ga_x[i];

        ts_xzzz_xxxzz[i] = 3.0 * ts_zzz_xxzz[i] * gfe_0 + ts_zzz_xxxzz[i] * ga_x[i];

        ts_xzzz_xxyyy[i] = 2.0 * ts_zzz_xyyy[i] * gfe_0 + ts_zzz_xxyyy[i] * ga_x[i];

        ts_xzzz_xxyyz[i] = 2.0 * ts_zzz_xyyz[i] * gfe_0 + ts_zzz_xxyyz[i] * ga_x[i];

        ts_xzzz_xxyzz[i] = 2.0 * ts_zzz_xyzz[i] * gfe_0 + ts_zzz_xxyzz[i] * ga_x[i];

        ts_xzzz_xxzzz[i] = 2.0 * ts_zzz_xzzz[i] * gfe_0 + ts_zzz_xxzzz[i] * ga_x[i];

        ts_xzzz_xyyyy[i] = ts_zzz_yyyy[i] * gfe_0 + ts_zzz_xyyyy[i] * ga_x[i];

        ts_xzzz_xyyyz[i] = ts_zzz_yyyz[i] * gfe_0 + ts_zzz_xyyyz[i] * ga_x[i];

        ts_xzzz_xyyzz[i] = ts_zzz_yyzz[i] * gfe_0 + ts_zzz_xyyzz[i] * ga_x[i];

        ts_xzzz_xyzzz[i] = ts_zzz_yzzz[i] * gfe_0 + ts_zzz_xyzzz[i] * ga_x[i];

        ts_xzzz_xzzzz[i] = ts_zzz_zzzz[i] * gfe_0 + ts_zzz_xzzzz[i] * ga_x[i];

        ts_xzzz_yyyyy[i] = ts_zzz_yyyyy[i] * ga_x[i];

        ts_xzzz_yyyyz[i] = ts_zzz_yyyyz[i] * ga_x[i];

        ts_xzzz_yyyzz[i] = ts_zzz_yyyzz[i] * ga_x[i];

        ts_xzzz_yyzzz[i] = ts_zzz_yyzzz[i] * ga_x[i];

        ts_xzzz_yzzzz[i] = ts_zzz_yzzzz[i] * ga_x[i];

        ts_xzzz_zzzzz[i] = ts_zzz_zzzzz[i] * ga_x[i];
    }

    // Set up 210-231 components of targeted buffer : GH

    auto ts_yyyy_xxxxx = pbuffer.data(idx_gh + 210);

    auto ts_yyyy_xxxxy = pbuffer.data(idx_gh + 211);

    auto ts_yyyy_xxxxz = pbuffer.data(idx_gh + 212);

    auto ts_yyyy_xxxyy = pbuffer.data(idx_gh + 213);

    auto ts_yyyy_xxxyz = pbuffer.data(idx_gh + 214);

    auto ts_yyyy_xxxzz = pbuffer.data(idx_gh + 215);

    auto ts_yyyy_xxyyy = pbuffer.data(idx_gh + 216);

    auto ts_yyyy_xxyyz = pbuffer.data(idx_gh + 217);

    auto ts_yyyy_xxyzz = pbuffer.data(idx_gh + 218);

    auto ts_yyyy_xxzzz = pbuffer.data(idx_gh + 219);

    auto ts_yyyy_xyyyy = pbuffer.data(idx_gh + 220);

    auto ts_yyyy_xyyyz = pbuffer.data(idx_gh + 221);

    auto ts_yyyy_xyyzz = pbuffer.data(idx_gh + 222);

    auto ts_yyyy_xyzzz = pbuffer.data(idx_gh + 223);

    auto ts_yyyy_xzzzz = pbuffer.data(idx_gh + 224);

    auto ts_yyyy_yyyyy = pbuffer.data(idx_gh + 225);

    auto ts_yyyy_yyyyz = pbuffer.data(idx_gh + 226);

    auto ts_yyyy_yyyzz = pbuffer.data(idx_gh + 227);

    auto ts_yyyy_yyzzz = pbuffer.data(idx_gh + 228);

    auto ts_yyyy_yzzzz = pbuffer.data(idx_gh + 229);

    auto ts_yyyy_zzzzz = pbuffer.data(idx_gh + 230);

    #pragma omp simd aligned(ga_y, ts_yy_xxxxx, ts_yy_xxxxy, ts_yy_xxxxz, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxxzz, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyzz, ts_yy_xxzzz, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyzz, ts_yy_xyzzz, ts_yy_xzzzz, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyzz, ts_yy_yyzzz, ts_yy_yzzzz, ts_yy_zzzzz, ts_yyy_xxxx, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxy, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxz, ts_yyy_xxxzz, ts_yyy_xxyy, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyz, ts_yyy_xxyzz, ts_yyy_xxzz, ts_yyy_xxzzz, ts_yyy_xyyy, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyz, ts_yyy_xyyzz, ts_yyy_xyzz, ts_yyy_xyzzz, ts_yyy_xzzz, ts_yyy_xzzzz, ts_yyy_yyyy, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyz, ts_yyy_yyyzz, ts_yyy_yyzz, ts_yyy_yyzzz, ts_yyy_yzzz, ts_yyy_yzzzz, ts_yyy_zzzz, ts_yyy_zzzzz, ts_yyyy_xxxxx, ts_yyyy_xxxxy, ts_yyyy_xxxxz, ts_yyyy_xxxyy, ts_yyyy_xxxyz, ts_yyyy_xxxzz, ts_yyyy_xxyyy, ts_yyyy_xxyyz, ts_yyyy_xxyzz, ts_yyyy_xxzzz, ts_yyyy_xyyyy, ts_yyyy_xyyyz, ts_yyyy_xyyzz, ts_yyyy_xyzzz, ts_yyyy_xzzzz, ts_yyyy_yyyyy, ts_yyyy_yyyyz, ts_yyyy_yyyzz, ts_yyyy_yyzzz, ts_yyyy_yzzzz, ts_yyyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyy_xxxxx[i] = 3.0 * ts_yy_xxxxx[i] * gfe_0 + ts_yyy_xxxxx[i] * ga_y[i];

        ts_yyyy_xxxxy[i] = 3.0 * ts_yy_xxxxy[i] * gfe_0 + ts_yyy_xxxx[i] * gfe_0 + ts_yyy_xxxxy[i] * ga_y[i];

        ts_yyyy_xxxxz[i] = 3.0 * ts_yy_xxxxz[i] * gfe_0 + ts_yyy_xxxxz[i] * ga_y[i];

        ts_yyyy_xxxyy[i] = 3.0 * ts_yy_xxxyy[i] * gfe_0 + 2.0 * ts_yyy_xxxy[i] * gfe_0 + ts_yyy_xxxyy[i] * ga_y[i];

        ts_yyyy_xxxyz[i] = 3.0 * ts_yy_xxxyz[i] * gfe_0 + ts_yyy_xxxz[i] * gfe_0 + ts_yyy_xxxyz[i] * ga_y[i];

        ts_yyyy_xxxzz[i] = 3.0 * ts_yy_xxxzz[i] * gfe_0 + ts_yyy_xxxzz[i] * ga_y[i];

        ts_yyyy_xxyyy[i] = 3.0 * ts_yy_xxyyy[i] * gfe_0 + 3.0 * ts_yyy_xxyy[i] * gfe_0 + ts_yyy_xxyyy[i] * ga_y[i];

        ts_yyyy_xxyyz[i] = 3.0 * ts_yy_xxyyz[i] * gfe_0 + 2.0 * ts_yyy_xxyz[i] * gfe_0 + ts_yyy_xxyyz[i] * ga_y[i];

        ts_yyyy_xxyzz[i] = 3.0 * ts_yy_xxyzz[i] * gfe_0 + ts_yyy_xxzz[i] * gfe_0 + ts_yyy_xxyzz[i] * ga_y[i];

        ts_yyyy_xxzzz[i] = 3.0 * ts_yy_xxzzz[i] * gfe_0 + ts_yyy_xxzzz[i] * ga_y[i];

        ts_yyyy_xyyyy[i] = 3.0 * ts_yy_xyyyy[i] * gfe_0 + 4.0 * ts_yyy_xyyy[i] * gfe_0 + ts_yyy_xyyyy[i] * ga_y[i];

        ts_yyyy_xyyyz[i] = 3.0 * ts_yy_xyyyz[i] * gfe_0 + 3.0 * ts_yyy_xyyz[i] * gfe_0 + ts_yyy_xyyyz[i] * ga_y[i];

        ts_yyyy_xyyzz[i] = 3.0 * ts_yy_xyyzz[i] * gfe_0 + 2.0 * ts_yyy_xyzz[i] * gfe_0 + ts_yyy_xyyzz[i] * ga_y[i];

        ts_yyyy_xyzzz[i] = 3.0 * ts_yy_xyzzz[i] * gfe_0 + ts_yyy_xzzz[i] * gfe_0 + ts_yyy_xyzzz[i] * ga_y[i];

        ts_yyyy_xzzzz[i] = 3.0 * ts_yy_xzzzz[i] * gfe_0 + ts_yyy_xzzzz[i] * ga_y[i];

        ts_yyyy_yyyyy[i] = 3.0 * ts_yy_yyyyy[i] * gfe_0 + 5.0 * ts_yyy_yyyy[i] * gfe_0 + ts_yyy_yyyyy[i] * ga_y[i];

        ts_yyyy_yyyyz[i] = 3.0 * ts_yy_yyyyz[i] * gfe_0 + 4.0 * ts_yyy_yyyz[i] * gfe_0 + ts_yyy_yyyyz[i] * ga_y[i];

        ts_yyyy_yyyzz[i] = 3.0 * ts_yy_yyyzz[i] * gfe_0 + 3.0 * ts_yyy_yyzz[i] * gfe_0 + ts_yyy_yyyzz[i] * ga_y[i];

        ts_yyyy_yyzzz[i] = 3.0 * ts_yy_yyzzz[i] * gfe_0 + 2.0 * ts_yyy_yzzz[i] * gfe_0 + ts_yyy_yyzzz[i] * ga_y[i];

        ts_yyyy_yzzzz[i] = 3.0 * ts_yy_yzzzz[i] * gfe_0 + ts_yyy_zzzz[i] * gfe_0 + ts_yyy_yzzzz[i] * ga_y[i];

        ts_yyyy_zzzzz[i] = 3.0 * ts_yy_zzzzz[i] * gfe_0 + ts_yyy_zzzzz[i] * ga_y[i];
    }

    // Set up 231-252 components of targeted buffer : GH

    auto ts_yyyz_xxxxx = pbuffer.data(idx_gh + 231);

    auto ts_yyyz_xxxxy = pbuffer.data(idx_gh + 232);

    auto ts_yyyz_xxxxz = pbuffer.data(idx_gh + 233);

    auto ts_yyyz_xxxyy = pbuffer.data(idx_gh + 234);

    auto ts_yyyz_xxxyz = pbuffer.data(idx_gh + 235);

    auto ts_yyyz_xxxzz = pbuffer.data(idx_gh + 236);

    auto ts_yyyz_xxyyy = pbuffer.data(idx_gh + 237);

    auto ts_yyyz_xxyyz = pbuffer.data(idx_gh + 238);

    auto ts_yyyz_xxyzz = pbuffer.data(idx_gh + 239);

    auto ts_yyyz_xxzzz = pbuffer.data(idx_gh + 240);

    auto ts_yyyz_xyyyy = pbuffer.data(idx_gh + 241);

    auto ts_yyyz_xyyyz = pbuffer.data(idx_gh + 242);

    auto ts_yyyz_xyyzz = pbuffer.data(idx_gh + 243);

    auto ts_yyyz_xyzzz = pbuffer.data(idx_gh + 244);

    auto ts_yyyz_xzzzz = pbuffer.data(idx_gh + 245);

    auto ts_yyyz_yyyyy = pbuffer.data(idx_gh + 246);

    auto ts_yyyz_yyyyz = pbuffer.data(idx_gh + 247);

    auto ts_yyyz_yyyzz = pbuffer.data(idx_gh + 248);

    auto ts_yyyz_yyzzz = pbuffer.data(idx_gh + 249);

    auto ts_yyyz_yzzzz = pbuffer.data(idx_gh + 250);

    auto ts_yyyz_zzzzz = pbuffer.data(idx_gh + 251);

    #pragma omp simd aligned(ga_y, ga_z, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxy, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxyy, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyz, ts_yyy_xxyzz, ts_yyy_xyyy, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyz, ts_yyy_xyyzz, ts_yyy_xyzz, ts_yyy_xyzzz, ts_yyy_yyyy, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyz, ts_yyy_yyyzz, ts_yyy_yyzz, ts_yyy_yyzzz, ts_yyy_yzzz, ts_yyy_yzzzz, ts_yyyz_xxxxx, ts_yyyz_xxxxy, ts_yyyz_xxxxz, ts_yyyz_xxxyy, ts_yyyz_xxxyz, ts_yyyz_xxxzz, ts_yyyz_xxyyy, ts_yyyz_xxyyz, ts_yyyz_xxyzz, ts_yyyz_xxzzz, ts_yyyz_xyyyy, ts_yyyz_xyyyz, ts_yyyz_xyyzz, ts_yyyz_xyzzz, ts_yyyz_xzzzz, ts_yyyz_yyyyy, ts_yyyz_yyyyz, ts_yyyz_yyyzz, ts_yyyz_yyzzz, ts_yyyz_yzzzz, ts_yyyz_zzzzz, ts_yyz_xxxxz, ts_yyz_xxxzz, ts_yyz_xxzzz, ts_yyz_xzzzz, ts_yyz_zzzzz, ts_yz_xxxxz, ts_yz_xxxzz, ts_yz_xxzzz, ts_yz_xzzzz, ts_yz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyz_xxxxx[i] = ts_yyy_xxxxx[i] * ga_z[i];

        ts_yyyz_xxxxy[i] = ts_yyy_xxxxy[i] * ga_z[i];

        ts_yyyz_xxxxz[i] = 2.0 * ts_yz_xxxxz[i] * gfe_0 + ts_yyz_xxxxz[i] * ga_y[i];

        ts_yyyz_xxxyy[i] = ts_yyy_xxxyy[i] * ga_z[i];

        ts_yyyz_xxxyz[i] = ts_yyy_xxxy[i] * gfe_0 + ts_yyy_xxxyz[i] * ga_z[i];

        ts_yyyz_xxxzz[i] = 2.0 * ts_yz_xxxzz[i] * gfe_0 + ts_yyz_xxxzz[i] * ga_y[i];

        ts_yyyz_xxyyy[i] = ts_yyy_xxyyy[i] * ga_z[i];

        ts_yyyz_xxyyz[i] = ts_yyy_xxyy[i] * gfe_0 + ts_yyy_xxyyz[i] * ga_z[i];

        ts_yyyz_xxyzz[i] = 2.0 * ts_yyy_xxyz[i] * gfe_0 + ts_yyy_xxyzz[i] * ga_z[i];

        ts_yyyz_xxzzz[i] = 2.0 * ts_yz_xxzzz[i] * gfe_0 + ts_yyz_xxzzz[i] * ga_y[i];

        ts_yyyz_xyyyy[i] = ts_yyy_xyyyy[i] * ga_z[i];

        ts_yyyz_xyyyz[i] = ts_yyy_xyyy[i] * gfe_0 + ts_yyy_xyyyz[i] * ga_z[i];

        ts_yyyz_xyyzz[i] = 2.0 * ts_yyy_xyyz[i] * gfe_0 + ts_yyy_xyyzz[i] * ga_z[i];

        ts_yyyz_xyzzz[i] = 3.0 * ts_yyy_xyzz[i] * gfe_0 + ts_yyy_xyzzz[i] * ga_z[i];

        ts_yyyz_xzzzz[i] = 2.0 * ts_yz_xzzzz[i] * gfe_0 + ts_yyz_xzzzz[i] * ga_y[i];

        ts_yyyz_yyyyy[i] = ts_yyy_yyyyy[i] * ga_z[i];

        ts_yyyz_yyyyz[i] = ts_yyy_yyyy[i] * gfe_0 + ts_yyy_yyyyz[i] * ga_z[i];

        ts_yyyz_yyyzz[i] = 2.0 * ts_yyy_yyyz[i] * gfe_0 + ts_yyy_yyyzz[i] * ga_z[i];

        ts_yyyz_yyzzz[i] = 3.0 * ts_yyy_yyzz[i] * gfe_0 + ts_yyy_yyzzz[i] * ga_z[i];

        ts_yyyz_yzzzz[i] = 4.0 * ts_yyy_yzzz[i] * gfe_0 + ts_yyy_yzzzz[i] * ga_z[i];

        ts_yyyz_zzzzz[i] = 2.0 * ts_yz_zzzzz[i] * gfe_0 + ts_yyz_zzzzz[i] * ga_y[i];
    }

    // Set up 252-273 components of targeted buffer : GH

    auto ts_yyzz_xxxxx = pbuffer.data(idx_gh + 252);

    auto ts_yyzz_xxxxy = pbuffer.data(idx_gh + 253);

    auto ts_yyzz_xxxxz = pbuffer.data(idx_gh + 254);

    auto ts_yyzz_xxxyy = pbuffer.data(idx_gh + 255);

    auto ts_yyzz_xxxyz = pbuffer.data(idx_gh + 256);

    auto ts_yyzz_xxxzz = pbuffer.data(idx_gh + 257);

    auto ts_yyzz_xxyyy = pbuffer.data(idx_gh + 258);

    auto ts_yyzz_xxyyz = pbuffer.data(idx_gh + 259);

    auto ts_yyzz_xxyzz = pbuffer.data(idx_gh + 260);

    auto ts_yyzz_xxzzz = pbuffer.data(idx_gh + 261);

    auto ts_yyzz_xyyyy = pbuffer.data(idx_gh + 262);

    auto ts_yyzz_xyyyz = pbuffer.data(idx_gh + 263);

    auto ts_yyzz_xyyzz = pbuffer.data(idx_gh + 264);

    auto ts_yyzz_xyzzz = pbuffer.data(idx_gh + 265);

    auto ts_yyzz_xzzzz = pbuffer.data(idx_gh + 266);

    auto ts_yyzz_yyyyy = pbuffer.data(idx_gh + 267);

    auto ts_yyzz_yyyyz = pbuffer.data(idx_gh + 268);

    auto ts_yyzz_yyyzz = pbuffer.data(idx_gh + 269);

    auto ts_yyzz_yyzzz = pbuffer.data(idx_gh + 270);

    auto ts_yyzz_yzzzz = pbuffer.data(idx_gh + 271);

    auto ts_yyzz_zzzzz = pbuffer.data(idx_gh + 272);

    #pragma omp simd aligned(ga_y, ga_z, ts_yy_xxxxy, ts_yy_xxxyy, ts_yy_xxyyy, ts_yy_xyyyy, ts_yy_yyyyy, ts_yyz_xxxxy, ts_yyz_xxxyy, ts_yyz_xxyyy, ts_yyz_xyyyy, ts_yyz_yyyyy, ts_yyzz_xxxxx, ts_yyzz_xxxxy, ts_yyzz_xxxxz, ts_yyzz_xxxyy, ts_yyzz_xxxyz, ts_yyzz_xxxzz, ts_yyzz_xxyyy, ts_yyzz_xxyyz, ts_yyzz_xxyzz, ts_yyzz_xxzzz, ts_yyzz_xyyyy, ts_yyzz_xyyyz, ts_yyzz_xyyzz, ts_yyzz_xyzzz, ts_yyzz_xzzzz, ts_yyzz_yyyyy, ts_yyzz_yyyyz, ts_yyzz_yyyzz, ts_yyzz_yyzzz, ts_yyzz_yzzzz, ts_yyzz_zzzzz, ts_yzz_xxxxx, ts_yzz_xxxxz, ts_yzz_xxxyz, ts_yzz_xxxz, ts_yzz_xxxzz, ts_yzz_xxyyz, ts_yzz_xxyz, ts_yzz_xxyzz, ts_yzz_xxzz, ts_yzz_xxzzz, ts_yzz_xyyyz, ts_yzz_xyyz, ts_yzz_xyyzz, ts_yzz_xyzz, ts_yzz_xyzzz, ts_yzz_xzzz, ts_yzz_xzzzz, ts_yzz_yyyyz, ts_yzz_yyyz, ts_yzz_yyyzz, ts_yzz_yyzz, ts_yzz_yyzzz, ts_yzz_yzzz, ts_yzz_yzzzz, ts_yzz_zzzz, ts_yzz_zzzzz, ts_zz_xxxxx, ts_zz_xxxxz, ts_zz_xxxyz, ts_zz_xxxzz, ts_zz_xxyyz, ts_zz_xxyzz, ts_zz_xxzzz, ts_zz_xyyyz, ts_zz_xyyzz, ts_zz_xyzzz, ts_zz_xzzzz, ts_zz_yyyyz, ts_zz_yyyzz, ts_zz_yyzzz, ts_zz_yzzzz, ts_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyzz_xxxxx[i] = ts_zz_xxxxx[i] * gfe_0 + ts_yzz_xxxxx[i] * ga_y[i];

        ts_yyzz_xxxxy[i] = ts_yy_xxxxy[i] * gfe_0 + ts_yyz_xxxxy[i] * ga_z[i];

        ts_yyzz_xxxxz[i] = ts_zz_xxxxz[i] * gfe_0 + ts_yzz_xxxxz[i] * ga_y[i];

        ts_yyzz_xxxyy[i] = ts_yy_xxxyy[i] * gfe_0 + ts_yyz_xxxyy[i] * ga_z[i];

        ts_yyzz_xxxyz[i] = ts_zz_xxxyz[i] * gfe_0 + ts_yzz_xxxz[i] * gfe_0 + ts_yzz_xxxyz[i] * ga_y[i];

        ts_yyzz_xxxzz[i] = ts_zz_xxxzz[i] * gfe_0 + ts_yzz_xxxzz[i] * ga_y[i];

        ts_yyzz_xxyyy[i] = ts_yy_xxyyy[i] * gfe_0 + ts_yyz_xxyyy[i] * ga_z[i];

        ts_yyzz_xxyyz[i] = ts_zz_xxyyz[i] * gfe_0 + 2.0 * ts_yzz_xxyz[i] * gfe_0 + ts_yzz_xxyyz[i] * ga_y[i];

        ts_yyzz_xxyzz[i] = ts_zz_xxyzz[i] * gfe_0 + ts_yzz_xxzz[i] * gfe_0 + ts_yzz_xxyzz[i] * ga_y[i];

        ts_yyzz_xxzzz[i] = ts_zz_xxzzz[i] * gfe_0 + ts_yzz_xxzzz[i] * ga_y[i];

        ts_yyzz_xyyyy[i] = ts_yy_xyyyy[i] * gfe_0 + ts_yyz_xyyyy[i] * ga_z[i];

        ts_yyzz_xyyyz[i] = ts_zz_xyyyz[i] * gfe_0 + 3.0 * ts_yzz_xyyz[i] * gfe_0 + ts_yzz_xyyyz[i] * ga_y[i];

        ts_yyzz_xyyzz[i] = ts_zz_xyyzz[i] * gfe_0 + 2.0 * ts_yzz_xyzz[i] * gfe_0 + ts_yzz_xyyzz[i] * ga_y[i];

        ts_yyzz_xyzzz[i] = ts_zz_xyzzz[i] * gfe_0 + ts_yzz_xzzz[i] * gfe_0 + ts_yzz_xyzzz[i] * ga_y[i];

        ts_yyzz_xzzzz[i] = ts_zz_xzzzz[i] * gfe_0 + ts_yzz_xzzzz[i] * ga_y[i];

        ts_yyzz_yyyyy[i] = ts_yy_yyyyy[i] * gfe_0 + ts_yyz_yyyyy[i] * ga_z[i];

        ts_yyzz_yyyyz[i] = ts_zz_yyyyz[i] * gfe_0 + 4.0 * ts_yzz_yyyz[i] * gfe_0 + ts_yzz_yyyyz[i] * ga_y[i];

        ts_yyzz_yyyzz[i] = ts_zz_yyyzz[i] * gfe_0 + 3.0 * ts_yzz_yyzz[i] * gfe_0 + ts_yzz_yyyzz[i] * ga_y[i];

        ts_yyzz_yyzzz[i] = ts_zz_yyzzz[i] * gfe_0 + 2.0 * ts_yzz_yzzz[i] * gfe_0 + ts_yzz_yyzzz[i] * ga_y[i];

        ts_yyzz_yzzzz[i] = ts_zz_yzzzz[i] * gfe_0 + ts_yzz_zzzz[i] * gfe_0 + ts_yzz_yzzzz[i] * ga_y[i];

        ts_yyzz_zzzzz[i] = ts_zz_zzzzz[i] * gfe_0 + ts_yzz_zzzzz[i] * ga_y[i];
    }

    // Set up 273-294 components of targeted buffer : GH

    auto ts_yzzz_xxxxx = pbuffer.data(idx_gh + 273);

    auto ts_yzzz_xxxxy = pbuffer.data(idx_gh + 274);

    auto ts_yzzz_xxxxz = pbuffer.data(idx_gh + 275);

    auto ts_yzzz_xxxyy = pbuffer.data(idx_gh + 276);

    auto ts_yzzz_xxxyz = pbuffer.data(idx_gh + 277);

    auto ts_yzzz_xxxzz = pbuffer.data(idx_gh + 278);

    auto ts_yzzz_xxyyy = pbuffer.data(idx_gh + 279);

    auto ts_yzzz_xxyyz = pbuffer.data(idx_gh + 280);

    auto ts_yzzz_xxyzz = pbuffer.data(idx_gh + 281);

    auto ts_yzzz_xxzzz = pbuffer.data(idx_gh + 282);

    auto ts_yzzz_xyyyy = pbuffer.data(idx_gh + 283);

    auto ts_yzzz_xyyyz = pbuffer.data(idx_gh + 284);

    auto ts_yzzz_xyyzz = pbuffer.data(idx_gh + 285);

    auto ts_yzzz_xyzzz = pbuffer.data(idx_gh + 286);

    auto ts_yzzz_xzzzz = pbuffer.data(idx_gh + 287);

    auto ts_yzzz_yyyyy = pbuffer.data(idx_gh + 288);

    auto ts_yzzz_yyyyz = pbuffer.data(idx_gh + 289);

    auto ts_yzzz_yyyzz = pbuffer.data(idx_gh + 290);

    auto ts_yzzz_yyzzz = pbuffer.data(idx_gh + 291);

    auto ts_yzzz_yzzzz = pbuffer.data(idx_gh + 292);

    auto ts_yzzz_zzzzz = pbuffer.data(idx_gh + 293);

    #pragma omp simd aligned(ga_y, ts_yzzz_xxxxx, ts_yzzz_xxxxy, ts_yzzz_xxxxz, ts_yzzz_xxxyy, ts_yzzz_xxxyz, ts_yzzz_xxxzz, ts_yzzz_xxyyy, ts_yzzz_xxyyz, ts_yzzz_xxyzz, ts_yzzz_xxzzz, ts_yzzz_xyyyy, ts_yzzz_xyyyz, ts_yzzz_xyyzz, ts_yzzz_xyzzz, ts_yzzz_xzzzz, ts_yzzz_yyyyy, ts_yzzz_yyyyz, ts_yzzz_yyyzz, ts_yzzz_yyzzz, ts_yzzz_yzzzz, ts_yzzz_zzzzz, ts_zzz_xxxx, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxy, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxz, ts_zzz_xxxzz, ts_zzz_xxyy, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyz, ts_zzz_xxyzz, ts_zzz_xxzz, ts_zzz_xxzzz, ts_zzz_xyyy, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyz, ts_zzz_xyyzz, ts_zzz_xyzz, ts_zzz_xyzzz, ts_zzz_xzzz, ts_zzz_xzzzz, ts_zzz_yyyy, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyz, ts_zzz_yyyzz, ts_zzz_yyzz, ts_zzz_yyzzz, ts_zzz_yzzz, ts_zzz_yzzzz, ts_zzz_zzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yzzz_xxxxx[i] = ts_zzz_xxxxx[i] * ga_y[i];

        ts_yzzz_xxxxy[i] = ts_zzz_xxxx[i] * gfe_0 + ts_zzz_xxxxy[i] * ga_y[i];

        ts_yzzz_xxxxz[i] = ts_zzz_xxxxz[i] * ga_y[i];

        ts_yzzz_xxxyy[i] = 2.0 * ts_zzz_xxxy[i] * gfe_0 + ts_zzz_xxxyy[i] * ga_y[i];

        ts_yzzz_xxxyz[i] = ts_zzz_xxxz[i] * gfe_0 + ts_zzz_xxxyz[i] * ga_y[i];

        ts_yzzz_xxxzz[i] = ts_zzz_xxxzz[i] * ga_y[i];

        ts_yzzz_xxyyy[i] = 3.0 * ts_zzz_xxyy[i] * gfe_0 + ts_zzz_xxyyy[i] * ga_y[i];

        ts_yzzz_xxyyz[i] = 2.0 * ts_zzz_xxyz[i] * gfe_0 + ts_zzz_xxyyz[i] * ga_y[i];

        ts_yzzz_xxyzz[i] = ts_zzz_xxzz[i] * gfe_0 + ts_zzz_xxyzz[i] * ga_y[i];

        ts_yzzz_xxzzz[i] = ts_zzz_xxzzz[i] * ga_y[i];

        ts_yzzz_xyyyy[i] = 4.0 * ts_zzz_xyyy[i] * gfe_0 + ts_zzz_xyyyy[i] * ga_y[i];

        ts_yzzz_xyyyz[i] = 3.0 * ts_zzz_xyyz[i] * gfe_0 + ts_zzz_xyyyz[i] * ga_y[i];

        ts_yzzz_xyyzz[i] = 2.0 * ts_zzz_xyzz[i] * gfe_0 + ts_zzz_xyyzz[i] * ga_y[i];

        ts_yzzz_xyzzz[i] = ts_zzz_xzzz[i] * gfe_0 + ts_zzz_xyzzz[i] * ga_y[i];

        ts_yzzz_xzzzz[i] = ts_zzz_xzzzz[i] * ga_y[i];

        ts_yzzz_yyyyy[i] = 5.0 * ts_zzz_yyyy[i] * gfe_0 + ts_zzz_yyyyy[i] * ga_y[i];

        ts_yzzz_yyyyz[i] = 4.0 * ts_zzz_yyyz[i] * gfe_0 + ts_zzz_yyyyz[i] * ga_y[i];

        ts_yzzz_yyyzz[i] = 3.0 * ts_zzz_yyzz[i] * gfe_0 + ts_zzz_yyyzz[i] * ga_y[i];

        ts_yzzz_yyzzz[i] = 2.0 * ts_zzz_yzzz[i] * gfe_0 + ts_zzz_yyzzz[i] * ga_y[i];

        ts_yzzz_yzzzz[i] = ts_zzz_zzzz[i] * gfe_0 + ts_zzz_yzzzz[i] * ga_y[i];

        ts_yzzz_zzzzz[i] = ts_zzz_zzzzz[i] * ga_y[i];
    }

    // Set up 294-315 components of targeted buffer : GH

    auto ts_zzzz_xxxxx = pbuffer.data(idx_gh + 294);

    auto ts_zzzz_xxxxy = pbuffer.data(idx_gh + 295);

    auto ts_zzzz_xxxxz = pbuffer.data(idx_gh + 296);

    auto ts_zzzz_xxxyy = pbuffer.data(idx_gh + 297);

    auto ts_zzzz_xxxyz = pbuffer.data(idx_gh + 298);

    auto ts_zzzz_xxxzz = pbuffer.data(idx_gh + 299);

    auto ts_zzzz_xxyyy = pbuffer.data(idx_gh + 300);

    auto ts_zzzz_xxyyz = pbuffer.data(idx_gh + 301);

    auto ts_zzzz_xxyzz = pbuffer.data(idx_gh + 302);

    auto ts_zzzz_xxzzz = pbuffer.data(idx_gh + 303);

    auto ts_zzzz_xyyyy = pbuffer.data(idx_gh + 304);

    auto ts_zzzz_xyyyz = pbuffer.data(idx_gh + 305);

    auto ts_zzzz_xyyzz = pbuffer.data(idx_gh + 306);

    auto ts_zzzz_xyzzz = pbuffer.data(idx_gh + 307);

    auto ts_zzzz_xzzzz = pbuffer.data(idx_gh + 308);

    auto ts_zzzz_yyyyy = pbuffer.data(idx_gh + 309);

    auto ts_zzzz_yyyyz = pbuffer.data(idx_gh + 310);

    auto ts_zzzz_yyyzz = pbuffer.data(idx_gh + 311);

    auto ts_zzzz_yyzzz = pbuffer.data(idx_gh + 312);

    auto ts_zzzz_yzzzz = pbuffer.data(idx_gh + 313);

    auto ts_zzzz_zzzzz = pbuffer.data(idx_gh + 314);

    #pragma omp simd aligned(ga_z, ts_zz_xxxxx, ts_zz_xxxxy, ts_zz_xxxxz, ts_zz_xxxyy, ts_zz_xxxyz, ts_zz_xxxzz, ts_zz_xxyyy, ts_zz_xxyyz, ts_zz_xxyzz, ts_zz_xxzzz, ts_zz_xyyyy, ts_zz_xyyyz, ts_zz_xyyzz, ts_zz_xyzzz, ts_zz_xzzzz, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyzz, ts_zz_yyzzz, ts_zz_yzzzz, ts_zz_zzzzz, ts_zzz_xxxx, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxy, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxz, ts_zzz_xxxzz, ts_zzz_xxyy, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyz, ts_zzz_xxyzz, ts_zzz_xxzz, ts_zzz_xxzzz, ts_zzz_xyyy, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyz, ts_zzz_xyyzz, ts_zzz_xyzz, ts_zzz_xyzzz, ts_zzz_xzzz, ts_zzz_xzzzz, ts_zzz_yyyy, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyz, ts_zzz_yyyzz, ts_zzz_yyzz, ts_zzz_yyzzz, ts_zzz_yzzz, ts_zzz_yzzzz, ts_zzz_zzzz, ts_zzz_zzzzz, ts_zzzz_xxxxx, ts_zzzz_xxxxy, ts_zzzz_xxxxz, ts_zzzz_xxxyy, ts_zzzz_xxxyz, ts_zzzz_xxxzz, ts_zzzz_xxyyy, ts_zzzz_xxyyz, ts_zzzz_xxyzz, ts_zzzz_xxzzz, ts_zzzz_xyyyy, ts_zzzz_xyyyz, ts_zzzz_xyyzz, ts_zzzz_xyzzz, ts_zzzz_xzzzz, ts_zzzz_yyyyy, ts_zzzz_yyyyz, ts_zzzz_yyyzz, ts_zzzz_yyzzz, ts_zzzz_yzzzz, ts_zzzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_zzzz_xxxxx[i] = 3.0 * ts_zz_xxxxx[i] * gfe_0 + ts_zzz_xxxxx[i] * ga_z[i];

        ts_zzzz_xxxxy[i] = 3.0 * ts_zz_xxxxy[i] * gfe_0 + ts_zzz_xxxxy[i] * ga_z[i];

        ts_zzzz_xxxxz[i] = 3.0 * ts_zz_xxxxz[i] * gfe_0 + ts_zzz_xxxx[i] * gfe_0 + ts_zzz_xxxxz[i] * ga_z[i];

        ts_zzzz_xxxyy[i] = 3.0 * ts_zz_xxxyy[i] * gfe_0 + ts_zzz_xxxyy[i] * ga_z[i];

        ts_zzzz_xxxyz[i] = 3.0 * ts_zz_xxxyz[i] * gfe_0 + ts_zzz_xxxy[i] * gfe_0 + ts_zzz_xxxyz[i] * ga_z[i];

        ts_zzzz_xxxzz[i] = 3.0 * ts_zz_xxxzz[i] * gfe_0 + 2.0 * ts_zzz_xxxz[i] * gfe_0 + ts_zzz_xxxzz[i] * ga_z[i];

        ts_zzzz_xxyyy[i] = 3.0 * ts_zz_xxyyy[i] * gfe_0 + ts_zzz_xxyyy[i] * ga_z[i];

        ts_zzzz_xxyyz[i] = 3.0 * ts_zz_xxyyz[i] * gfe_0 + ts_zzz_xxyy[i] * gfe_0 + ts_zzz_xxyyz[i] * ga_z[i];

        ts_zzzz_xxyzz[i] = 3.0 * ts_zz_xxyzz[i] * gfe_0 + 2.0 * ts_zzz_xxyz[i] * gfe_0 + ts_zzz_xxyzz[i] * ga_z[i];

        ts_zzzz_xxzzz[i] = 3.0 * ts_zz_xxzzz[i] * gfe_0 + 3.0 * ts_zzz_xxzz[i] * gfe_0 + ts_zzz_xxzzz[i] * ga_z[i];

        ts_zzzz_xyyyy[i] = 3.0 * ts_zz_xyyyy[i] * gfe_0 + ts_zzz_xyyyy[i] * ga_z[i];

        ts_zzzz_xyyyz[i] = 3.0 * ts_zz_xyyyz[i] * gfe_0 + ts_zzz_xyyy[i] * gfe_0 + ts_zzz_xyyyz[i] * ga_z[i];

        ts_zzzz_xyyzz[i] = 3.0 * ts_zz_xyyzz[i] * gfe_0 + 2.0 * ts_zzz_xyyz[i] * gfe_0 + ts_zzz_xyyzz[i] * ga_z[i];

        ts_zzzz_xyzzz[i] = 3.0 * ts_zz_xyzzz[i] * gfe_0 + 3.0 * ts_zzz_xyzz[i] * gfe_0 + ts_zzz_xyzzz[i] * ga_z[i];

        ts_zzzz_xzzzz[i] = 3.0 * ts_zz_xzzzz[i] * gfe_0 + 4.0 * ts_zzz_xzzz[i] * gfe_0 + ts_zzz_xzzzz[i] * ga_z[i];

        ts_zzzz_yyyyy[i] = 3.0 * ts_zz_yyyyy[i] * gfe_0 + ts_zzz_yyyyy[i] * ga_z[i];

        ts_zzzz_yyyyz[i] = 3.0 * ts_zz_yyyyz[i] * gfe_0 + ts_zzz_yyyy[i] * gfe_0 + ts_zzz_yyyyz[i] * ga_z[i];

        ts_zzzz_yyyzz[i] = 3.0 * ts_zz_yyyzz[i] * gfe_0 + 2.0 * ts_zzz_yyyz[i] * gfe_0 + ts_zzz_yyyzz[i] * ga_z[i];

        ts_zzzz_yyzzz[i] = 3.0 * ts_zz_yyzzz[i] * gfe_0 + 3.0 * ts_zzz_yyzz[i] * gfe_0 + ts_zzz_yyzzz[i] * ga_z[i];

        ts_zzzz_yzzzz[i] = 3.0 * ts_zz_yzzzz[i] * gfe_0 + 4.0 * ts_zzz_yzzz[i] * gfe_0 + ts_zzz_yzzzz[i] * ga_z[i];

        ts_zzzz_zzzzz[i] = 3.0 * ts_zz_zzzzz[i] * gfe_0 + 5.0 * ts_zzz_zzzz[i] * gfe_0 + ts_zzz_zzzzz[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

