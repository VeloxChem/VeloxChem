#include "T2CHrrABRecFG.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_fg(CSimdArray<double>& cbuffer, 
            const size_t idx_fg,
            const size_t idx_dg,
            const size_t idx_dh,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : DG

    auto t_xx_xxxx = cbuffer.data(idx_dg);

    auto t_xx_xxxy = cbuffer.data(idx_dg + 1);

    auto t_xx_xxxz = cbuffer.data(idx_dg + 2);

    auto t_xx_xxyy = cbuffer.data(idx_dg + 3);

    auto t_xx_xxyz = cbuffer.data(idx_dg + 4);

    auto t_xx_xxzz = cbuffer.data(idx_dg + 5);

    auto t_xx_xyyy = cbuffer.data(idx_dg + 6);

    auto t_xx_xyyz = cbuffer.data(idx_dg + 7);

    auto t_xx_xyzz = cbuffer.data(idx_dg + 8);

    auto t_xx_xzzz = cbuffer.data(idx_dg + 9);

    auto t_xx_yyyy = cbuffer.data(idx_dg + 10);

    auto t_xx_yyyz = cbuffer.data(idx_dg + 11);

    auto t_xx_yyzz = cbuffer.data(idx_dg + 12);

    auto t_xx_yzzz = cbuffer.data(idx_dg + 13);

    auto t_xx_zzzz = cbuffer.data(idx_dg + 14);

    auto t_xy_xxxx = cbuffer.data(idx_dg + 15);

    auto t_xy_xxxy = cbuffer.data(idx_dg + 16);

    auto t_xy_xxxz = cbuffer.data(idx_dg + 17);

    auto t_xy_xxyy = cbuffer.data(idx_dg + 18);

    auto t_xy_xxyz = cbuffer.data(idx_dg + 19);

    auto t_xy_xxzz = cbuffer.data(idx_dg + 20);

    auto t_xy_xyyy = cbuffer.data(idx_dg + 21);

    auto t_xy_xyyz = cbuffer.data(idx_dg + 22);

    auto t_xy_xyzz = cbuffer.data(idx_dg + 23);

    auto t_xy_xzzz = cbuffer.data(idx_dg + 24);

    auto t_xy_yyyy = cbuffer.data(idx_dg + 25);

    auto t_xy_yyyz = cbuffer.data(idx_dg + 26);

    auto t_xy_yyzz = cbuffer.data(idx_dg + 27);

    auto t_xy_yzzz = cbuffer.data(idx_dg + 28);

    auto t_xy_zzzz = cbuffer.data(idx_dg + 29);

    auto t_xz_xxxx = cbuffer.data(idx_dg + 30);

    auto t_xz_xxxy = cbuffer.data(idx_dg + 31);

    auto t_xz_xxxz = cbuffer.data(idx_dg + 32);

    auto t_xz_xxyy = cbuffer.data(idx_dg + 33);

    auto t_xz_xxyz = cbuffer.data(idx_dg + 34);

    auto t_xz_xxzz = cbuffer.data(idx_dg + 35);

    auto t_xz_xyyy = cbuffer.data(idx_dg + 36);

    auto t_xz_xyyz = cbuffer.data(idx_dg + 37);

    auto t_xz_xyzz = cbuffer.data(idx_dg + 38);

    auto t_xz_xzzz = cbuffer.data(idx_dg + 39);

    auto t_xz_yyyy = cbuffer.data(idx_dg + 40);

    auto t_xz_yyyz = cbuffer.data(idx_dg + 41);

    auto t_xz_yyzz = cbuffer.data(idx_dg + 42);

    auto t_xz_yzzz = cbuffer.data(idx_dg + 43);

    auto t_xz_zzzz = cbuffer.data(idx_dg + 44);

    auto t_yy_xxxx = cbuffer.data(idx_dg + 45);

    auto t_yy_xxxy = cbuffer.data(idx_dg + 46);

    auto t_yy_xxxz = cbuffer.data(idx_dg + 47);

    auto t_yy_xxyy = cbuffer.data(idx_dg + 48);

    auto t_yy_xxyz = cbuffer.data(idx_dg + 49);

    auto t_yy_xxzz = cbuffer.data(idx_dg + 50);

    auto t_yy_xyyy = cbuffer.data(idx_dg + 51);

    auto t_yy_xyyz = cbuffer.data(idx_dg + 52);

    auto t_yy_xyzz = cbuffer.data(idx_dg + 53);

    auto t_yy_xzzz = cbuffer.data(idx_dg + 54);

    auto t_yy_yyyy = cbuffer.data(idx_dg + 55);

    auto t_yy_yyyz = cbuffer.data(idx_dg + 56);

    auto t_yy_yyzz = cbuffer.data(idx_dg + 57);

    auto t_yy_yzzz = cbuffer.data(idx_dg + 58);

    auto t_yy_zzzz = cbuffer.data(idx_dg + 59);

    auto t_yz_xxxx = cbuffer.data(idx_dg + 60);

    auto t_yz_xxxy = cbuffer.data(idx_dg + 61);

    auto t_yz_xxxz = cbuffer.data(idx_dg + 62);

    auto t_yz_xxyy = cbuffer.data(idx_dg + 63);

    auto t_yz_xxyz = cbuffer.data(idx_dg + 64);

    auto t_yz_xxzz = cbuffer.data(idx_dg + 65);

    auto t_yz_xyyy = cbuffer.data(idx_dg + 66);

    auto t_yz_xyyz = cbuffer.data(idx_dg + 67);

    auto t_yz_xyzz = cbuffer.data(idx_dg + 68);

    auto t_yz_xzzz = cbuffer.data(idx_dg + 69);

    auto t_yz_yyyy = cbuffer.data(idx_dg + 70);

    auto t_yz_yyyz = cbuffer.data(idx_dg + 71);

    auto t_yz_yyzz = cbuffer.data(idx_dg + 72);

    auto t_yz_yzzz = cbuffer.data(idx_dg + 73);

    auto t_yz_zzzz = cbuffer.data(idx_dg + 74);

    auto t_zz_xxxx = cbuffer.data(idx_dg + 75);

    auto t_zz_xxxy = cbuffer.data(idx_dg + 76);

    auto t_zz_xxxz = cbuffer.data(idx_dg + 77);

    auto t_zz_xxyy = cbuffer.data(idx_dg + 78);

    auto t_zz_xxyz = cbuffer.data(idx_dg + 79);

    auto t_zz_xxzz = cbuffer.data(idx_dg + 80);

    auto t_zz_xyyy = cbuffer.data(idx_dg + 81);

    auto t_zz_xyyz = cbuffer.data(idx_dg + 82);

    auto t_zz_xyzz = cbuffer.data(idx_dg + 83);

    auto t_zz_xzzz = cbuffer.data(idx_dg + 84);

    auto t_zz_yyyy = cbuffer.data(idx_dg + 85);

    auto t_zz_yyyz = cbuffer.data(idx_dg + 86);

    auto t_zz_yyzz = cbuffer.data(idx_dg + 87);

    auto t_zz_yzzz = cbuffer.data(idx_dg + 88);

    auto t_zz_zzzz = cbuffer.data(idx_dg + 89);

    // Set up components of auxiliary buffer : DH

    auto t_xx_xxxxx = cbuffer.data(idx_dh);

    auto t_xx_xxxxy = cbuffer.data(idx_dh + 1);

    auto t_xx_xxxxz = cbuffer.data(idx_dh + 2);

    auto t_xx_xxxyy = cbuffer.data(idx_dh + 3);

    auto t_xx_xxxyz = cbuffer.data(idx_dh + 4);

    auto t_xx_xxxzz = cbuffer.data(idx_dh + 5);

    auto t_xx_xxyyy = cbuffer.data(idx_dh + 6);

    auto t_xx_xxyyz = cbuffer.data(idx_dh + 7);

    auto t_xx_xxyzz = cbuffer.data(idx_dh + 8);

    auto t_xx_xxzzz = cbuffer.data(idx_dh + 9);

    auto t_xx_xyyyy = cbuffer.data(idx_dh + 10);

    auto t_xx_xyyyz = cbuffer.data(idx_dh + 11);

    auto t_xx_xyyzz = cbuffer.data(idx_dh + 12);

    auto t_xx_xyzzz = cbuffer.data(idx_dh + 13);

    auto t_xx_xzzzz = cbuffer.data(idx_dh + 14);

    auto t_xy_xxxxx = cbuffer.data(idx_dh + 21);

    auto t_xy_xxxxy = cbuffer.data(idx_dh + 22);

    auto t_xy_xxxxz = cbuffer.data(idx_dh + 23);

    auto t_xy_xxxyy = cbuffer.data(idx_dh + 24);

    auto t_xy_xxxyz = cbuffer.data(idx_dh + 25);

    auto t_xy_xxxzz = cbuffer.data(idx_dh + 26);

    auto t_xy_xxyyy = cbuffer.data(idx_dh + 27);

    auto t_xy_xxyyz = cbuffer.data(idx_dh + 28);

    auto t_xy_xxyzz = cbuffer.data(idx_dh + 29);

    auto t_xy_xxzzz = cbuffer.data(idx_dh + 30);

    auto t_xy_xyyyy = cbuffer.data(idx_dh + 31);

    auto t_xy_xyyyz = cbuffer.data(idx_dh + 32);

    auto t_xy_xyyzz = cbuffer.data(idx_dh + 33);

    auto t_xy_xyzzz = cbuffer.data(idx_dh + 34);

    auto t_xy_xzzzz = cbuffer.data(idx_dh + 35);

    auto t_xz_xxxxx = cbuffer.data(idx_dh + 42);

    auto t_xz_xxxxy = cbuffer.data(idx_dh + 43);

    auto t_xz_xxxxz = cbuffer.data(idx_dh + 44);

    auto t_xz_xxxyy = cbuffer.data(idx_dh + 45);

    auto t_xz_xxxyz = cbuffer.data(idx_dh + 46);

    auto t_xz_xxxzz = cbuffer.data(idx_dh + 47);

    auto t_xz_xxyyy = cbuffer.data(idx_dh + 48);

    auto t_xz_xxyyz = cbuffer.data(idx_dh + 49);

    auto t_xz_xxyzz = cbuffer.data(idx_dh + 50);

    auto t_xz_xxzzz = cbuffer.data(idx_dh + 51);

    auto t_xz_xyyyy = cbuffer.data(idx_dh + 52);

    auto t_xz_xyyyz = cbuffer.data(idx_dh + 53);

    auto t_xz_xyyzz = cbuffer.data(idx_dh + 54);

    auto t_xz_xyzzz = cbuffer.data(idx_dh + 55);

    auto t_xz_xzzzz = cbuffer.data(idx_dh + 56);

    auto t_yy_xxxxx = cbuffer.data(idx_dh + 63);

    auto t_yy_xxxxy = cbuffer.data(idx_dh + 64);

    auto t_yy_xxxxz = cbuffer.data(idx_dh + 65);

    auto t_yy_xxxyy = cbuffer.data(idx_dh + 66);

    auto t_yy_xxxyz = cbuffer.data(idx_dh + 67);

    auto t_yy_xxxzz = cbuffer.data(idx_dh + 68);

    auto t_yy_xxyyy = cbuffer.data(idx_dh + 69);

    auto t_yy_xxyyz = cbuffer.data(idx_dh + 70);

    auto t_yy_xxyzz = cbuffer.data(idx_dh + 71);

    auto t_yy_xxzzz = cbuffer.data(idx_dh + 72);

    auto t_yy_xyyyy = cbuffer.data(idx_dh + 73);

    auto t_yy_xyyyz = cbuffer.data(idx_dh + 74);

    auto t_yy_xyyzz = cbuffer.data(idx_dh + 75);

    auto t_yy_xyzzz = cbuffer.data(idx_dh + 76);

    auto t_yy_xzzzz = cbuffer.data(idx_dh + 77);

    auto t_yy_yyyyy = cbuffer.data(idx_dh + 78);

    auto t_yy_yyyyz = cbuffer.data(idx_dh + 79);

    auto t_yy_yyyzz = cbuffer.data(idx_dh + 80);

    auto t_yy_yyzzz = cbuffer.data(idx_dh + 81);

    auto t_yy_yzzzz = cbuffer.data(idx_dh + 82);

    auto t_yz_xxxxx = cbuffer.data(idx_dh + 84);

    auto t_yz_xxxxy = cbuffer.data(idx_dh + 85);

    auto t_yz_xxxxz = cbuffer.data(idx_dh + 86);

    auto t_yz_xxxyy = cbuffer.data(idx_dh + 87);

    auto t_yz_xxxyz = cbuffer.data(idx_dh + 88);

    auto t_yz_xxxzz = cbuffer.data(idx_dh + 89);

    auto t_yz_xxyyy = cbuffer.data(idx_dh + 90);

    auto t_yz_xxyyz = cbuffer.data(idx_dh + 91);

    auto t_yz_xxyzz = cbuffer.data(idx_dh + 92);

    auto t_yz_xxzzz = cbuffer.data(idx_dh + 93);

    auto t_yz_xyyyy = cbuffer.data(idx_dh + 94);

    auto t_yz_xyyyz = cbuffer.data(idx_dh + 95);

    auto t_yz_xyyzz = cbuffer.data(idx_dh + 96);

    auto t_yz_xyzzz = cbuffer.data(idx_dh + 97);

    auto t_yz_xzzzz = cbuffer.data(idx_dh + 98);

    auto t_yz_yyyyy = cbuffer.data(idx_dh + 99);

    auto t_yz_yyyyz = cbuffer.data(idx_dh + 100);

    auto t_yz_yyyzz = cbuffer.data(idx_dh + 101);

    auto t_yz_yyzzz = cbuffer.data(idx_dh + 102);

    auto t_yz_yzzzz = cbuffer.data(idx_dh + 103);

    auto t_zz_xxxxx = cbuffer.data(idx_dh + 105);

    auto t_zz_xxxxy = cbuffer.data(idx_dh + 106);

    auto t_zz_xxxxz = cbuffer.data(idx_dh + 107);

    auto t_zz_xxxyy = cbuffer.data(idx_dh + 108);

    auto t_zz_xxxyz = cbuffer.data(idx_dh + 109);

    auto t_zz_xxxzz = cbuffer.data(idx_dh + 110);

    auto t_zz_xxyyy = cbuffer.data(idx_dh + 111);

    auto t_zz_xxyyz = cbuffer.data(idx_dh + 112);

    auto t_zz_xxyzz = cbuffer.data(idx_dh + 113);

    auto t_zz_xxzzz = cbuffer.data(idx_dh + 114);

    auto t_zz_xyyyy = cbuffer.data(idx_dh + 115);

    auto t_zz_xyyyz = cbuffer.data(idx_dh + 116);

    auto t_zz_xyyzz = cbuffer.data(idx_dh + 117);

    auto t_zz_xyzzz = cbuffer.data(idx_dh + 118);

    auto t_zz_xzzzz = cbuffer.data(idx_dh + 119);

    auto t_zz_yyyyy = cbuffer.data(idx_dh + 120);

    auto t_zz_yyyyz = cbuffer.data(idx_dh + 121);

    auto t_zz_yyyzz = cbuffer.data(idx_dh + 122);

    auto t_zz_yyzzz = cbuffer.data(idx_dh + 123);

    auto t_zz_yzzzz = cbuffer.data(idx_dh + 124);

    auto t_zz_zzzzz = cbuffer.data(idx_dh + 125);

    // Set up components of targeted buffer : FG

    auto t_xxx_xxxx = cbuffer.data(idx_fg);

    auto t_xxx_xxxy = cbuffer.data(idx_fg + 1);

    auto t_xxx_xxxz = cbuffer.data(idx_fg + 2);

    auto t_xxx_xxyy = cbuffer.data(idx_fg + 3);

    auto t_xxx_xxyz = cbuffer.data(idx_fg + 4);

    auto t_xxx_xxzz = cbuffer.data(idx_fg + 5);

    auto t_xxx_xyyy = cbuffer.data(idx_fg + 6);

    auto t_xxx_xyyz = cbuffer.data(idx_fg + 7);

    auto t_xxx_xyzz = cbuffer.data(idx_fg + 8);

    auto t_xxx_xzzz = cbuffer.data(idx_fg + 9);

    auto t_xxx_yyyy = cbuffer.data(idx_fg + 10);

    auto t_xxx_yyyz = cbuffer.data(idx_fg + 11);

    auto t_xxx_yyzz = cbuffer.data(idx_fg + 12);

    auto t_xxx_yzzz = cbuffer.data(idx_fg + 13);

    auto t_xxx_zzzz = cbuffer.data(idx_fg + 14);

    auto t_xxy_xxxx = cbuffer.data(idx_fg + 15);

    auto t_xxy_xxxy = cbuffer.data(idx_fg + 16);

    auto t_xxy_xxxz = cbuffer.data(idx_fg + 17);

    auto t_xxy_xxyy = cbuffer.data(idx_fg + 18);

    auto t_xxy_xxyz = cbuffer.data(idx_fg + 19);

    auto t_xxy_xxzz = cbuffer.data(idx_fg + 20);

    auto t_xxy_xyyy = cbuffer.data(idx_fg + 21);

    auto t_xxy_xyyz = cbuffer.data(idx_fg + 22);

    auto t_xxy_xyzz = cbuffer.data(idx_fg + 23);

    auto t_xxy_xzzz = cbuffer.data(idx_fg + 24);

    auto t_xxy_yyyy = cbuffer.data(idx_fg + 25);

    auto t_xxy_yyyz = cbuffer.data(idx_fg + 26);

    auto t_xxy_yyzz = cbuffer.data(idx_fg + 27);

    auto t_xxy_yzzz = cbuffer.data(idx_fg + 28);

    auto t_xxy_zzzz = cbuffer.data(idx_fg + 29);

    auto t_xxz_xxxx = cbuffer.data(idx_fg + 30);

    auto t_xxz_xxxy = cbuffer.data(idx_fg + 31);

    auto t_xxz_xxxz = cbuffer.data(idx_fg + 32);

    auto t_xxz_xxyy = cbuffer.data(idx_fg + 33);

    auto t_xxz_xxyz = cbuffer.data(idx_fg + 34);

    auto t_xxz_xxzz = cbuffer.data(idx_fg + 35);

    auto t_xxz_xyyy = cbuffer.data(idx_fg + 36);

    auto t_xxz_xyyz = cbuffer.data(idx_fg + 37);

    auto t_xxz_xyzz = cbuffer.data(idx_fg + 38);

    auto t_xxz_xzzz = cbuffer.data(idx_fg + 39);

    auto t_xxz_yyyy = cbuffer.data(idx_fg + 40);

    auto t_xxz_yyyz = cbuffer.data(idx_fg + 41);

    auto t_xxz_yyzz = cbuffer.data(idx_fg + 42);

    auto t_xxz_yzzz = cbuffer.data(idx_fg + 43);

    auto t_xxz_zzzz = cbuffer.data(idx_fg + 44);

    auto t_xyy_xxxx = cbuffer.data(idx_fg + 45);

    auto t_xyy_xxxy = cbuffer.data(idx_fg + 46);

    auto t_xyy_xxxz = cbuffer.data(idx_fg + 47);

    auto t_xyy_xxyy = cbuffer.data(idx_fg + 48);

    auto t_xyy_xxyz = cbuffer.data(idx_fg + 49);

    auto t_xyy_xxzz = cbuffer.data(idx_fg + 50);

    auto t_xyy_xyyy = cbuffer.data(idx_fg + 51);

    auto t_xyy_xyyz = cbuffer.data(idx_fg + 52);

    auto t_xyy_xyzz = cbuffer.data(idx_fg + 53);

    auto t_xyy_xzzz = cbuffer.data(idx_fg + 54);

    auto t_xyy_yyyy = cbuffer.data(idx_fg + 55);

    auto t_xyy_yyyz = cbuffer.data(idx_fg + 56);

    auto t_xyy_yyzz = cbuffer.data(idx_fg + 57);

    auto t_xyy_yzzz = cbuffer.data(idx_fg + 58);

    auto t_xyy_zzzz = cbuffer.data(idx_fg + 59);

    auto t_xyz_xxxx = cbuffer.data(idx_fg + 60);

    auto t_xyz_xxxy = cbuffer.data(idx_fg + 61);

    auto t_xyz_xxxz = cbuffer.data(idx_fg + 62);

    auto t_xyz_xxyy = cbuffer.data(idx_fg + 63);

    auto t_xyz_xxyz = cbuffer.data(idx_fg + 64);

    auto t_xyz_xxzz = cbuffer.data(idx_fg + 65);

    auto t_xyz_xyyy = cbuffer.data(idx_fg + 66);

    auto t_xyz_xyyz = cbuffer.data(idx_fg + 67);

    auto t_xyz_xyzz = cbuffer.data(idx_fg + 68);

    auto t_xyz_xzzz = cbuffer.data(idx_fg + 69);

    auto t_xyz_yyyy = cbuffer.data(idx_fg + 70);

    auto t_xyz_yyyz = cbuffer.data(idx_fg + 71);

    auto t_xyz_yyzz = cbuffer.data(idx_fg + 72);

    auto t_xyz_yzzz = cbuffer.data(idx_fg + 73);

    auto t_xyz_zzzz = cbuffer.data(idx_fg + 74);

    auto t_xzz_xxxx = cbuffer.data(idx_fg + 75);

    auto t_xzz_xxxy = cbuffer.data(idx_fg + 76);

    auto t_xzz_xxxz = cbuffer.data(idx_fg + 77);

    auto t_xzz_xxyy = cbuffer.data(idx_fg + 78);

    auto t_xzz_xxyz = cbuffer.data(idx_fg + 79);

    auto t_xzz_xxzz = cbuffer.data(idx_fg + 80);

    auto t_xzz_xyyy = cbuffer.data(idx_fg + 81);

    auto t_xzz_xyyz = cbuffer.data(idx_fg + 82);

    auto t_xzz_xyzz = cbuffer.data(idx_fg + 83);

    auto t_xzz_xzzz = cbuffer.data(idx_fg + 84);

    auto t_xzz_yyyy = cbuffer.data(idx_fg + 85);

    auto t_xzz_yyyz = cbuffer.data(idx_fg + 86);

    auto t_xzz_yyzz = cbuffer.data(idx_fg + 87);

    auto t_xzz_yzzz = cbuffer.data(idx_fg + 88);

    auto t_xzz_zzzz = cbuffer.data(idx_fg + 89);

    auto t_yyy_xxxx = cbuffer.data(idx_fg + 90);

    auto t_yyy_xxxy = cbuffer.data(idx_fg + 91);

    auto t_yyy_xxxz = cbuffer.data(idx_fg + 92);

    auto t_yyy_xxyy = cbuffer.data(idx_fg + 93);

    auto t_yyy_xxyz = cbuffer.data(idx_fg + 94);

    auto t_yyy_xxzz = cbuffer.data(idx_fg + 95);

    auto t_yyy_xyyy = cbuffer.data(idx_fg + 96);

    auto t_yyy_xyyz = cbuffer.data(idx_fg + 97);

    auto t_yyy_xyzz = cbuffer.data(idx_fg + 98);

    auto t_yyy_xzzz = cbuffer.data(idx_fg + 99);

    auto t_yyy_yyyy = cbuffer.data(idx_fg + 100);

    auto t_yyy_yyyz = cbuffer.data(idx_fg + 101);

    auto t_yyy_yyzz = cbuffer.data(idx_fg + 102);

    auto t_yyy_yzzz = cbuffer.data(idx_fg + 103);

    auto t_yyy_zzzz = cbuffer.data(idx_fg + 104);

    auto t_yyz_xxxx = cbuffer.data(idx_fg + 105);

    auto t_yyz_xxxy = cbuffer.data(idx_fg + 106);

    auto t_yyz_xxxz = cbuffer.data(idx_fg + 107);

    auto t_yyz_xxyy = cbuffer.data(idx_fg + 108);

    auto t_yyz_xxyz = cbuffer.data(idx_fg + 109);

    auto t_yyz_xxzz = cbuffer.data(idx_fg + 110);

    auto t_yyz_xyyy = cbuffer.data(idx_fg + 111);

    auto t_yyz_xyyz = cbuffer.data(idx_fg + 112);

    auto t_yyz_xyzz = cbuffer.data(idx_fg + 113);

    auto t_yyz_xzzz = cbuffer.data(idx_fg + 114);

    auto t_yyz_yyyy = cbuffer.data(idx_fg + 115);

    auto t_yyz_yyyz = cbuffer.data(idx_fg + 116);

    auto t_yyz_yyzz = cbuffer.data(idx_fg + 117);

    auto t_yyz_yzzz = cbuffer.data(idx_fg + 118);

    auto t_yyz_zzzz = cbuffer.data(idx_fg + 119);

    auto t_yzz_xxxx = cbuffer.data(idx_fg + 120);

    auto t_yzz_xxxy = cbuffer.data(idx_fg + 121);

    auto t_yzz_xxxz = cbuffer.data(idx_fg + 122);

    auto t_yzz_xxyy = cbuffer.data(idx_fg + 123);

    auto t_yzz_xxyz = cbuffer.data(idx_fg + 124);

    auto t_yzz_xxzz = cbuffer.data(idx_fg + 125);

    auto t_yzz_xyyy = cbuffer.data(idx_fg + 126);

    auto t_yzz_xyyz = cbuffer.data(idx_fg + 127);

    auto t_yzz_xyzz = cbuffer.data(idx_fg + 128);

    auto t_yzz_xzzz = cbuffer.data(idx_fg + 129);

    auto t_yzz_yyyy = cbuffer.data(idx_fg + 130);

    auto t_yzz_yyyz = cbuffer.data(idx_fg + 131);

    auto t_yzz_yyzz = cbuffer.data(idx_fg + 132);

    auto t_yzz_yzzz = cbuffer.data(idx_fg + 133);

    auto t_yzz_zzzz = cbuffer.data(idx_fg + 134);

    auto t_zzz_xxxx = cbuffer.data(idx_fg + 135);

    auto t_zzz_xxxy = cbuffer.data(idx_fg + 136);

    auto t_zzz_xxxz = cbuffer.data(idx_fg + 137);

    auto t_zzz_xxyy = cbuffer.data(idx_fg + 138);

    auto t_zzz_xxyz = cbuffer.data(idx_fg + 139);

    auto t_zzz_xxzz = cbuffer.data(idx_fg + 140);

    auto t_zzz_xyyy = cbuffer.data(idx_fg + 141);

    auto t_zzz_xyyz = cbuffer.data(idx_fg + 142);

    auto t_zzz_xyzz = cbuffer.data(idx_fg + 143);

    auto t_zzz_xzzz = cbuffer.data(idx_fg + 144);

    auto t_zzz_yyyy = cbuffer.data(idx_fg + 145);

    auto t_zzz_yyyz = cbuffer.data(idx_fg + 146);

    auto t_zzz_yyzz = cbuffer.data(idx_fg + 147);

    auto t_zzz_yzzz = cbuffer.data(idx_fg + 148);

    auto t_zzz_zzzz = cbuffer.data(idx_fg + 149);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xx_xxxx, t_xx_xxxxx, t_xx_xxxxy, t_xx_xxxxz, t_xx_xxxy, t_xx_xxxyy, t_xx_xxxyz, t_xx_xxxz, t_xx_xxxzz, t_xx_xxyy, t_xx_xxyyy, t_xx_xxyyz, t_xx_xxyz, t_xx_xxyzz, t_xx_xxzz, t_xx_xxzzz, t_xx_xyyy, t_xx_xyyyy, t_xx_xyyyz, t_xx_xyyz, t_xx_xyyzz, t_xx_xyzz, t_xx_xyzzz, t_xx_xzzz, t_xx_xzzzz, t_xx_yyyy, t_xx_yyyz, t_xx_yyzz, t_xx_yzzz, t_xx_zzzz, t_xxx_xxxx, t_xxx_xxxy, t_xxx_xxxz, t_xxx_xxyy, t_xxx_xxyz, t_xxx_xxzz, t_xxx_xyyy, t_xxx_xyyz, t_xxx_xyzz, t_xxx_xzzz, t_xxx_yyyy, t_xxx_yyyz, t_xxx_yyzz, t_xxx_yzzz, t_xxx_zzzz, t_xxy_xxxx, t_xxy_xxxy, t_xxy_xxxz, t_xxy_xxyy, t_xxy_xxyz, t_xxy_xxzz, t_xxy_xyyy, t_xxy_xyyz, t_xxy_xyzz, t_xxy_xzzz, t_xxy_yyyy, t_xxy_yyyz, t_xxy_yyzz, t_xxy_yzzz, t_xxy_zzzz, t_xxz_xxxx, t_xxz_xxxy, t_xxz_xxxz, t_xxz_xxyy, t_xxz_xxyz, t_xxz_xxzz, t_xxz_xyyy, t_xxz_xyyz, t_xxz_xyzz, t_xxz_xzzz, t_xxz_yyyy, t_xxz_yyyz, t_xxz_yyzz, t_xxz_yzzz, t_xxz_zzzz, t_xy_xxxx, t_xy_xxxxx, t_xy_xxxxy, t_xy_xxxxz, t_xy_xxxy, t_xy_xxxyy, t_xy_xxxyz, t_xy_xxxz, t_xy_xxxzz, t_xy_xxyy, t_xy_xxyyy, t_xy_xxyyz, t_xy_xxyz, t_xy_xxyzz, t_xy_xxzz, t_xy_xxzzz, t_xy_xyyy, t_xy_xyyyy, t_xy_xyyyz, t_xy_xyyz, t_xy_xyyzz, t_xy_xyzz, t_xy_xyzzz, t_xy_xzzz, t_xy_xzzzz, t_xy_yyyy, t_xy_yyyz, t_xy_yyzz, t_xy_yzzz, t_xy_zzzz, t_xyy_xxxx, t_xyy_xxxy, t_xyy_xxxz, t_xyy_xxyy, t_xyy_xxyz, t_xyy_xxzz, t_xyy_xyyy, t_xyy_xyyz, t_xyy_xyzz, t_xyy_xzzz, t_xyy_yyyy, t_xyy_yyyz, t_xyy_yyzz, t_xyy_yzzz, t_xyy_zzzz, t_xyz_xxxx, t_xyz_xxxy, t_xyz_xxxz, t_xyz_xxyy, t_xyz_xxyz, t_xyz_xxzz, t_xyz_xyyy, t_xyz_xyyz, t_xyz_xyzz, t_xyz_xzzz, t_xyz_yyyy, t_xyz_yyyz, t_xyz_yyzz, t_xyz_yzzz, t_xyz_zzzz, t_xz_xxxx, t_xz_xxxxx, t_xz_xxxxy, t_xz_xxxxz, t_xz_xxxy, t_xz_xxxyy, t_xz_xxxyz, t_xz_xxxz, t_xz_xxxzz, t_xz_xxyy, t_xz_xxyyy, t_xz_xxyyz, t_xz_xxyz, t_xz_xxyzz, t_xz_xxzz, t_xz_xxzzz, t_xz_xyyy, t_xz_xyyyy, t_xz_xyyyz, t_xz_xyyz, t_xz_xyyzz, t_xz_xyzz, t_xz_xyzzz, t_xz_xzzz, t_xz_xzzzz, t_xz_yyyy, t_xz_yyyz, t_xz_yyzz, t_xz_yzzz, t_xz_zzzz, t_xzz_xxxx, t_xzz_xxxy, t_xzz_xxxz, t_xzz_xxyy, t_xzz_xxyz, t_xzz_xxzz, t_xzz_xyyy, t_xzz_xyyz, t_xzz_xyzz, t_xzz_xzzz, t_xzz_yyyy, t_xzz_yyyz, t_xzz_yyzz, t_xzz_yzzz, t_xzz_zzzz, t_yy_xxxx, t_yy_xxxxx, t_yy_xxxxy, t_yy_xxxxz, t_yy_xxxy, t_yy_xxxyy, t_yy_xxxyz, t_yy_xxxz, t_yy_xxxzz, t_yy_xxyy, t_yy_xxyyy, t_yy_xxyyz, t_yy_xxyz, t_yy_xxyzz, t_yy_xxzz, t_yy_xxzzz, t_yy_xyyy, t_yy_xyyyy, t_yy_xyyyz, t_yy_xyyz, t_yy_xyyzz, t_yy_xyzz, t_yy_xyzzz, t_yy_xzzz, t_yy_xzzzz, t_yy_yyyy, t_yy_yyyyy, t_yy_yyyyz, t_yy_yyyz, t_yy_yyyzz, t_yy_yyzz, t_yy_yyzzz, t_yy_yzzz, t_yy_yzzzz, t_yy_zzzz, t_yyy_xxxx, t_yyy_xxxy, t_yyy_xxxz, t_yyy_xxyy, t_yyy_xxyz, t_yyy_xxzz, t_yyy_xyyy, t_yyy_xyyz, t_yyy_xyzz, t_yyy_xzzz, t_yyy_yyyy, t_yyy_yyyz, t_yyy_yyzz, t_yyy_yzzz, t_yyy_zzzz, t_yyz_xxxx, t_yyz_xxxy, t_yyz_xxxz, t_yyz_xxyy, t_yyz_xxyz, t_yyz_xxzz, t_yyz_xyyy, t_yyz_xyyz, t_yyz_xyzz, t_yyz_xzzz, t_yyz_yyyy, t_yyz_yyyz, t_yyz_yyzz, t_yyz_yzzz, t_yyz_zzzz, t_yz_xxxx, t_yz_xxxxx, t_yz_xxxxy, t_yz_xxxxz, t_yz_xxxy, t_yz_xxxyy, t_yz_xxxyz, t_yz_xxxz, t_yz_xxxzz, t_yz_xxyy, t_yz_xxyyy, t_yz_xxyyz, t_yz_xxyz, t_yz_xxyzz, t_yz_xxzz, t_yz_xxzzz, t_yz_xyyy, t_yz_xyyyy, t_yz_xyyyz, t_yz_xyyz, t_yz_xyyzz, t_yz_xyzz, t_yz_xyzzz, t_yz_xzzz, t_yz_xzzzz, t_yz_yyyy, t_yz_yyyyy, t_yz_yyyyz, t_yz_yyyz, t_yz_yyyzz, t_yz_yyzz, t_yz_yyzzz, t_yz_yzzz, t_yz_yzzzz, t_yz_zzzz, t_yzz_xxxx, t_yzz_xxxy, t_yzz_xxxz, t_yzz_xxyy, t_yzz_xxyz, t_yzz_xxzz, t_yzz_xyyy, t_yzz_xyyz, t_yzz_xyzz, t_yzz_xzzz, t_yzz_yyyy, t_yzz_yyyz, t_yzz_yyzz, t_yzz_yzzz, t_yzz_zzzz, t_zz_xxxx, t_zz_xxxxx, t_zz_xxxxy, t_zz_xxxxz, t_zz_xxxy, t_zz_xxxyy, t_zz_xxxyz, t_zz_xxxz, t_zz_xxxzz, t_zz_xxyy, t_zz_xxyyy, t_zz_xxyyz, t_zz_xxyz, t_zz_xxyzz, t_zz_xxzz, t_zz_xxzzz, t_zz_xyyy, t_zz_xyyyy, t_zz_xyyyz, t_zz_xyyz, t_zz_xyyzz, t_zz_xyzz, t_zz_xyzzz, t_zz_xzzz, t_zz_xzzzz, t_zz_yyyy, t_zz_yyyyy, t_zz_yyyyz, t_zz_yyyz, t_zz_yyyzz, t_zz_yyzz, t_zz_yyzzz, t_zz_yzzz, t_zz_yzzzz, t_zz_zzzz, t_zz_zzzzz, t_zzz_xxxx, t_zzz_xxxy, t_zzz_xxxz, t_zzz_xxyy, t_zzz_xxyz, t_zzz_xxzz, t_zzz_xyyy, t_zzz_xyyz, t_zzz_xyzz, t_zzz_xzzz, t_zzz_yyyy, t_zzz_yyyz, t_zzz_yyzz, t_zzz_yzzz, t_zzz_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxx_xxxx[i] = -t_xx_xxxx[i] * ab_x[i] + t_xx_xxxxx[i];

        t_xxx_xxxy[i] = -t_xx_xxxy[i] * ab_x[i] + t_xx_xxxxy[i];

        t_xxx_xxxz[i] = -t_xx_xxxz[i] * ab_x[i] + t_xx_xxxxz[i];

        t_xxx_xxyy[i] = -t_xx_xxyy[i] * ab_x[i] + t_xx_xxxyy[i];

        t_xxx_xxyz[i] = -t_xx_xxyz[i] * ab_x[i] + t_xx_xxxyz[i];

        t_xxx_xxzz[i] = -t_xx_xxzz[i] * ab_x[i] + t_xx_xxxzz[i];

        t_xxx_xyyy[i] = -t_xx_xyyy[i] * ab_x[i] + t_xx_xxyyy[i];

        t_xxx_xyyz[i] = -t_xx_xyyz[i] * ab_x[i] + t_xx_xxyyz[i];

        t_xxx_xyzz[i] = -t_xx_xyzz[i] * ab_x[i] + t_xx_xxyzz[i];

        t_xxx_xzzz[i] = -t_xx_xzzz[i] * ab_x[i] + t_xx_xxzzz[i];

        t_xxx_yyyy[i] = -t_xx_yyyy[i] * ab_x[i] + t_xx_xyyyy[i];

        t_xxx_yyyz[i] = -t_xx_yyyz[i] * ab_x[i] + t_xx_xyyyz[i];

        t_xxx_yyzz[i] = -t_xx_yyzz[i] * ab_x[i] + t_xx_xyyzz[i];

        t_xxx_yzzz[i] = -t_xx_yzzz[i] * ab_x[i] + t_xx_xyzzz[i];

        t_xxx_zzzz[i] = -t_xx_zzzz[i] * ab_x[i] + t_xx_xzzzz[i];

        t_xxy_xxxx[i] = -t_xy_xxxx[i] * ab_x[i] + t_xy_xxxxx[i];

        t_xxy_xxxy[i] = -t_xy_xxxy[i] * ab_x[i] + t_xy_xxxxy[i];

        t_xxy_xxxz[i] = -t_xy_xxxz[i] * ab_x[i] + t_xy_xxxxz[i];

        t_xxy_xxyy[i] = -t_xy_xxyy[i] * ab_x[i] + t_xy_xxxyy[i];

        t_xxy_xxyz[i] = -t_xy_xxyz[i] * ab_x[i] + t_xy_xxxyz[i];

        t_xxy_xxzz[i] = -t_xy_xxzz[i] * ab_x[i] + t_xy_xxxzz[i];

        t_xxy_xyyy[i] = -t_xy_xyyy[i] * ab_x[i] + t_xy_xxyyy[i];

        t_xxy_xyyz[i] = -t_xy_xyyz[i] * ab_x[i] + t_xy_xxyyz[i];

        t_xxy_xyzz[i] = -t_xy_xyzz[i] * ab_x[i] + t_xy_xxyzz[i];

        t_xxy_xzzz[i] = -t_xy_xzzz[i] * ab_x[i] + t_xy_xxzzz[i];

        t_xxy_yyyy[i] = -t_xy_yyyy[i] * ab_x[i] + t_xy_xyyyy[i];

        t_xxy_yyyz[i] = -t_xy_yyyz[i] * ab_x[i] + t_xy_xyyyz[i];

        t_xxy_yyzz[i] = -t_xy_yyzz[i] * ab_x[i] + t_xy_xyyzz[i];

        t_xxy_yzzz[i] = -t_xy_yzzz[i] * ab_x[i] + t_xy_xyzzz[i];

        t_xxy_zzzz[i] = -t_xy_zzzz[i] * ab_x[i] + t_xy_xzzzz[i];

        t_xxz_xxxx[i] = -t_xz_xxxx[i] * ab_x[i] + t_xz_xxxxx[i];

        t_xxz_xxxy[i] = -t_xz_xxxy[i] * ab_x[i] + t_xz_xxxxy[i];

        t_xxz_xxxz[i] = -t_xz_xxxz[i] * ab_x[i] + t_xz_xxxxz[i];

        t_xxz_xxyy[i] = -t_xz_xxyy[i] * ab_x[i] + t_xz_xxxyy[i];

        t_xxz_xxyz[i] = -t_xz_xxyz[i] * ab_x[i] + t_xz_xxxyz[i];

        t_xxz_xxzz[i] = -t_xz_xxzz[i] * ab_x[i] + t_xz_xxxzz[i];

        t_xxz_xyyy[i] = -t_xz_xyyy[i] * ab_x[i] + t_xz_xxyyy[i];

        t_xxz_xyyz[i] = -t_xz_xyyz[i] * ab_x[i] + t_xz_xxyyz[i];

        t_xxz_xyzz[i] = -t_xz_xyzz[i] * ab_x[i] + t_xz_xxyzz[i];

        t_xxz_xzzz[i] = -t_xz_xzzz[i] * ab_x[i] + t_xz_xxzzz[i];

        t_xxz_yyyy[i] = -t_xz_yyyy[i] * ab_x[i] + t_xz_xyyyy[i];

        t_xxz_yyyz[i] = -t_xz_yyyz[i] * ab_x[i] + t_xz_xyyyz[i];

        t_xxz_yyzz[i] = -t_xz_yyzz[i] * ab_x[i] + t_xz_xyyzz[i];

        t_xxz_yzzz[i] = -t_xz_yzzz[i] * ab_x[i] + t_xz_xyzzz[i];

        t_xxz_zzzz[i] = -t_xz_zzzz[i] * ab_x[i] + t_xz_xzzzz[i];

        t_xyy_xxxx[i] = -t_yy_xxxx[i] * ab_x[i] + t_yy_xxxxx[i];

        t_xyy_xxxy[i] = -t_yy_xxxy[i] * ab_x[i] + t_yy_xxxxy[i];

        t_xyy_xxxz[i] = -t_yy_xxxz[i] * ab_x[i] + t_yy_xxxxz[i];

        t_xyy_xxyy[i] = -t_yy_xxyy[i] * ab_x[i] + t_yy_xxxyy[i];

        t_xyy_xxyz[i] = -t_yy_xxyz[i] * ab_x[i] + t_yy_xxxyz[i];

        t_xyy_xxzz[i] = -t_yy_xxzz[i] * ab_x[i] + t_yy_xxxzz[i];

        t_xyy_xyyy[i] = -t_yy_xyyy[i] * ab_x[i] + t_yy_xxyyy[i];

        t_xyy_xyyz[i] = -t_yy_xyyz[i] * ab_x[i] + t_yy_xxyyz[i];

        t_xyy_xyzz[i] = -t_yy_xyzz[i] * ab_x[i] + t_yy_xxyzz[i];

        t_xyy_xzzz[i] = -t_yy_xzzz[i] * ab_x[i] + t_yy_xxzzz[i];

        t_xyy_yyyy[i] = -t_yy_yyyy[i] * ab_x[i] + t_yy_xyyyy[i];

        t_xyy_yyyz[i] = -t_yy_yyyz[i] * ab_x[i] + t_yy_xyyyz[i];

        t_xyy_yyzz[i] = -t_yy_yyzz[i] * ab_x[i] + t_yy_xyyzz[i];

        t_xyy_yzzz[i] = -t_yy_yzzz[i] * ab_x[i] + t_yy_xyzzz[i];

        t_xyy_zzzz[i] = -t_yy_zzzz[i] * ab_x[i] + t_yy_xzzzz[i];

        t_xyz_xxxx[i] = -t_yz_xxxx[i] * ab_x[i] + t_yz_xxxxx[i];

        t_xyz_xxxy[i] = -t_yz_xxxy[i] * ab_x[i] + t_yz_xxxxy[i];

        t_xyz_xxxz[i] = -t_yz_xxxz[i] * ab_x[i] + t_yz_xxxxz[i];

        t_xyz_xxyy[i] = -t_yz_xxyy[i] * ab_x[i] + t_yz_xxxyy[i];

        t_xyz_xxyz[i] = -t_yz_xxyz[i] * ab_x[i] + t_yz_xxxyz[i];

        t_xyz_xxzz[i] = -t_yz_xxzz[i] * ab_x[i] + t_yz_xxxzz[i];

        t_xyz_xyyy[i] = -t_yz_xyyy[i] * ab_x[i] + t_yz_xxyyy[i];

        t_xyz_xyyz[i] = -t_yz_xyyz[i] * ab_x[i] + t_yz_xxyyz[i];

        t_xyz_xyzz[i] = -t_yz_xyzz[i] * ab_x[i] + t_yz_xxyzz[i];

        t_xyz_xzzz[i] = -t_yz_xzzz[i] * ab_x[i] + t_yz_xxzzz[i];

        t_xyz_yyyy[i] = -t_yz_yyyy[i] * ab_x[i] + t_yz_xyyyy[i];

        t_xyz_yyyz[i] = -t_yz_yyyz[i] * ab_x[i] + t_yz_xyyyz[i];

        t_xyz_yyzz[i] = -t_yz_yyzz[i] * ab_x[i] + t_yz_xyyzz[i];

        t_xyz_yzzz[i] = -t_yz_yzzz[i] * ab_x[i] + t_yz_xyzzz[i];

        t_xyz_zzzz[i] = -t_yz_zzzz[i] * ab_x[i] + t_yz_xzzzz[i];

        t_xzz_xxxx[i] = -t_zz_xxxx[i] * ab_x[i] + t_zz_xxxxx[i];

        t_xzz_xxxy[i] = -t_zz_xxxy[i] * ab_x[i] + t_zz_xxxxy[i];

        t_xzz_xxxz[i] = -t_zz_xxxz[i] * ab_x[i] + t_zz_xxxxz[i];

        t_xzz_xxyy[i] = -t_zz_xxyy[i] * ab_x[i] + t_zz_xxxyy[i];

        t_xzz_xxyz[i] = -t_zz_xxyz[i] * ab_x[i] + t_zz_xxxyz[i];

        t_xzz_xxzz[i] = -t_zz_xxzz[i] * ab_x[i] + t_zz_xxxzz[i];

        t_xzz_xyyy[i] = -t_zz_xyyy[i] * ab_x[i] + t_zz_xxyyy[i];

        t_xzz_xyyz[i] = -t_zz_xyyz[i] * ab_x[i] + t_zz_xxyyz[i];

        t_xzz_xyzz[i] = -t_zz_xyzz[i] * ab_x[i] + t_zz_xxyzz[i];

        t_xzz_xzzz[i] = -t_zz_xzzz[i] * ab_x[i] + t_zz_xxzzz[i];

        t_xzz_yyyy[i] = -t_zz_yyyy[i] * ab_x[i] + t_zz_xyyyy[i];

        t_xzz_yyyz[i] = -t_zz_yyyz[i] * ab_x[i] + t_zz_xyyyz[i];

        t_xzz_yyzz[i] = -t_zz_yyzz[i] * ab_x[i] + t_zz_xyyzz[i];

        t_xzz_yzzz[i] = -t_zz_yzzz[i] * ab_x[i] + t_zz_xyzzz[i];

        t_xzz_zzzz[i] = -t_zz_zzzz[i] * ab_x[i] + t_zz_xzzzz[i];

        t_yyy_xxxx[i] = -t_yy_xxxx[i] * ab_y[i] + t_yy_xxxxy[i];

        t_yyy_xxxy[i] = -t_yy_xxxy[i] * ab_y[i] + t_yy_xxxyy[i];

        t_yyy_xxxz[i] = -t_yy_xxxz[i] * ab_y[i] + t_yy_xxxyz[i];

        t_yyy_xxyy[i] = -t_yy_xxyy[i] * ab_y[i] + t_yy_xxyyy[i];

        t_yyy_xxyz[i] = -t_yy_xxyz[i] * ab_y[i] + t_yy_xxyyz[i];

        t_yyy_xxzz[i] = -t_yy_xxzz[i] * ab_y[i] + t_yy_xxyzz[i];

        t_yyy_xyyy[i] = -t_yy_xyyy[i] * ab_y[i] + t_yy_xyyyy[i];

        t_yyy_xyyz[i] = -t_yy_xyyz[i] * ab_y[i] + t_yy_xyyyz[i];

        t_yyy_xyzz[i] = -t_yy_xyzz[i] * ab_y[i] + t_yy_xyyzz[i];

        t_yyy_xzzz[i] = -t_yy_xzzz[i] * ab_y[i] + t_yy_xyzzz[i];

        t_yyy_yyyy[i] = -t_yy_yyyy[i] * ab_y[i] + t_yy_yyyyy[i];

        t_yyy_yyyz[i] = -t_yy_yyyz[i] * ab_y[i] + t_yy_yyyyz[i];

        t_yyy_yyzz[i] = -t_yy_yyzz[i] * ab_y[i] + t_yy_yyyzz[i];

        t_yyy_yzzz[i] = -t_yy_yzzz[i] * ab_y[i] + t_yy_yyzzz[i];

        t_yyy_zzzz[i] = -t_yy_zzzz[i] * ab_y[i] + t_yy_yzzzz[i];

        t_yyz_xxxx[i] = -t_yz_xxxx[i] * ab_y[i] + t_yz_xxxxy[i];

        t_yyz_xxxy[i] = -t_yz_xxxy[i] * ab_y[i] + t_yz_xxxyy[i];

        t_yyz_xxxz[i] = -t_yz_xxxz[i] * ab_y[i] + t_yz_xxxyz[i];

        t_yyz_xxyy[i] = -t_yz_xxyy[i] * ab_y[i] + t_yz_xxyyy[i];

        t_yyz_xxyz[i] = -t_yz_xxyz[i] * ab_y[i] + t_yz_xxyyz[i];

        t_yyz_xxzz[i] = -t_yz_xxzz[i] * ab_y[i] + t_yz_xxyzz[i];

        t_yyz_xyyy[i] = -t_yz_xyyy[i] * ab_y[i] + t_yz_xyyyy[i];

        t_yyz_xyyz[i] = -t_yz_xyyz[i] * ab_y[i] + t_yz_xyyyz[i];

        t_yyz_xyzz[i] = -t_yz_xyzz[i] * ab_y[i] + t_yz_xyyzz[i];

        t_yyz_xzzz[i] = -t_yz_xzzz[i] * ab_y[i] + t_yz_xyzzz[i];

        t_yyz_yyyy[i] = -t_yz_yyyy[i] * ab_y[i] + t_yz_yyyyy[i];

        t_yyz_yyyz[i] = -t_yz_yyyz[i] * ab_y[i] + t_yz_yyyyz[i];

        t_yyz_yyzz[i] = -t_yz_yyzz[i] * ab_y[i] + t_yz_yyyzz[i];

        t_yyz_yzzz[i] = -t_yz_yzzz[i] * ab_y[i] + t_yz_yyzzz[i];

        t_yyz_zzzz[i] = -t_yz_zzzz[i] * ab_y[i] + t_yz_yzzzz[i];

        t_yzz_xxxx[i] = -t_zz_xxxx[i] * ab_y[i] + t_zz_xxxxy[i];

        t_yzz_xxxy[i] = -t_zz_xxxy[i] * ab_y[i] + t_zz_xxxyy[i];

        t_yzz_xxxz[i] = -t_zz_xxxz[i] * ab_y[i] + t_zz_xxxyz[i];

        t_yzz_xxyy[i] = -t_zz_xxyy[i] * ab_y[i] + t_zz_xxyyy[i];

        t_yzz_xxyz[i] = -t_zz_xxyz[i] * ab_y[i] + t_zz_xxyyz[i];

        t_yzz_xxzz[i] = -t_zz_xxzz[i] * ab_y[i] + t_zz_xxyzz[i];

        t_yzz_xyyy[i] = -t_zz_xyyy[i] * ab_y[i] + t_zz_xyyyy[i];

        t_yzz_xyyz[i] = -t_zz_xyyz[i] * ab_y[i] + t_zz_xyyyz[i];

        t_yzz_xyzz[i] = -t_zz_xyzz[i] * ab_y[i] + t_zz_xyyzz[i];

        t_yzz_xzzz[i] = -t_zz_xzzz[i] * ab_y[i] + t_zz_xyzzz[i];

        t_yzz_yyyy[i] = -t_zz_yyyy[i] * ab_y[i] + t_zz_yyyyy[i];

        t_yzz_yyyz[i] = -t_zz_yyyz[i] * ab_y[i] + t_zz_yyyyz[i];

        t_yzz_yyzz[i] = -t_zz_yyzz[i] * ab_y[i] + t_zz_yyyzz[i];

        t_yzz_yzzz[i] = -t_zz_yzzz[i] * ab_y[i] + t_zz_yyzzz[i];

        t_yzz_zzzz[i] = -t_zz_zzzz[i] * ab_y[i] + t_zz_yzzzz[i];

        t_zzz_xxxx[i] = -t_zz_xxxx[i] * ab_z[i] + t_zz_xxxxz[i];

        t_zzz_xxxy[i] = -t_zz_xxxy[i] * ab_z[i] + t_zz_xxxyz[i];

        t_zzz_xxxz[i] = -t_zz_xxxz[i] * ab_z[i] + t_zz_xxxzz[i];

        t_zzz_xxyy[i] = -t_zz_xxyy[i] * ab_z[i] + t_zz_xxyyz[i];

        t_zzz_xxyz[i] = -t_zz_xxyz[i] * ab_z[i] + t_zz_xxyzz[i];

        t_zzz_xxzz[i] = -t_zz_xxzz[i] * ab_z[i] + t_zz_xxzzz[i];

        t_zzz_xyyy[i] = -t_zz_xyyy[i] * ab_z[i] + t_zz_xyyyz[i];

        t_zzz_xyyz[i] = -t_zz_xyyz[i] * ab_z[i] + t_zz_xyyzz[i];

        t_zzz_xyzz[i] = -t_zz_xyzz[i] * ab_z[i] + t_zz_xyzzz[i];

        t_zzz_xzzz[i] = -t_zz_xzzz[i] * ab_z[i] + t_zz_xzzzz[i];

        t_zzz_yyyy[i] = -t_zz_yyyy[i] * ab_z[i] + t_zz_yyyyz[i];

        t_zzz_yyyz[i] = -t_zz_yyyz[i] * ab_z[i] + t_zz_yyyzz[i];

        t_zzz_yyzz[i] = -t_zz_yyzz[i] * ab_z[i] + t_zz_yyzzz[i];

        t_zzz_yzzz[i] = -t_zz_yzzz[i] * ab_z[i] + t_zz_yzzzz[i];

        t_zzz_zzzz[i] = -t_zz_zzzz[i] * ab_z[i] + t_zz_zzzzz[i];
    }
}

} // t2chrr namespace

