#include "T2CHrrABRecGF.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_gf(CSimdArray<double>& cbuffer, 
            const size_t idx_gf,
            const size_t idx_gd,
            const size_t idx_hd,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : GD

    auto t_xxxx_xx = cbuffer.data(idx_gd);

    auto t_xxxx_xy = cbuffer.data(idx_gd + 1);

    auto t_xxxx_xz = cbuffer.data(idx_gd + 2);

    auto t_xxxx_yy = cbuffer.data(idx_gd + 3);

    auto t_xxxx_yz = cbuffer.data(idx_gd + 4);

    auto t_xxxx_zz = cbuffer.data(idx_gd + 5);

    auto t_xxxy_xx = cbuffer.data(idx_gd + 6);

    auto t_xxxy_xy = cbuffer.data(idx_gd + 7);

    auto t_xxxy_xz = cbuffer.data(idx_gd + 8);

    auto t_xxxy_yy = cbuffer.data(idx_gd + 9);

    auto t_xxxy_yz = cbuffer.data(idx_gd + 10);

    auto t_xxxy_zz = cbuffer.data(idx_gd + 11);

    auto t_xxxz_xx = cbuffer.data(idx_gd + 12);

    auto t_xxxz_xy = cbuffer.data(idx_gd + 13);

    auto t_xxxz_xz = cbuffer.data(idx_gd + 14);

    auto t_xxxz_yy = cbuffer.data(idx_gd + 15);

    auto t_xxxz_yz = cbuffer.data(idx_gd + 16);

    auto t_xxxz_zz = cbuffer.data(idx_gd + 17);

    auto t_xxyy_xx = cbuffer.data(idx_gd + 18);

    auto t_xxyy_xy = cbuffer.data(idx_gd + 19);

    auto t_xxyy_xz = cbuffer.data(idx_gd + 20);

    auto t_xxyy_yy = cbuffer.data(idx_gd + 21);

    auto t_xxyy_yz = cbuffer.data(idx_gd + 22);

    auto t_xxyy_zz = cbuffer.data(idx_gd + 23);

    auto t_xxyz_xx = cbuffer.data(idx_gd + 24);

    auto t_xxyz_xy = cbuffer.data(idx_gd + 25);

    auto t_xxyz_xz = cbuffer.data(idx_gd + 26);

    auto t_xxyz_yy = cbuffer.data(idx_gd + 27);

    auto t_xxyz_yz = cbuffer.data(idx_gd + 28);

    auto t_xxyz_zz = cbuffer.data(idx_gd + 29);

    auto t_xxzz_xx = cbuffer.data(idx_gd + 30);

    auto t_xxzz_xy = cbuffer.data(idx_gd + 31);

    auto t_xxzz_xz = cbuffer.data(idx_gd + 32);

    auto t_xxzz_yy = cbuffer.data(idx_gd + 33);

    auto t_xxzz_yz = cbuffer.data(idx_gd + 34);

    auto t_xxzz_zz = cbuffer.data(idx_gd + 35);

    auto t_xyyy_xx = cbuffer.data(idx_gd + 36);

    auto t_xyyy_xy = cbuffer.data(idx_gd + 37);

    auto t_xyyy_xz = cbuffer.data(idx_gd + 38);

    auto t_xyyy_yy = cbuffer.data(idx_gd + 39);

    auto t_xyyy_yz = cbuffer.data(idx_gd + 40);

    auto t_xyyy_zz = cbuffer.data(idx_gd + 41);

    auto t_xyyz_xx = cbuffer.data(idx_gd + 42);

    auto t_xyyz_xy = cbuffer.data(idx_gd + 43);

    auto t_xyyz_xz = cbuffer.data(idx_gd + 44);

    auto t_xyyz_yy = cbuffer.data(idx_gd + 45);

    auto t_xyyz_yz = cbuffer.data(idx_gd + 46);

    auto t_xyyz_zz = cbuffer.data(idx_gd + 47);

    auto t_xyzz_xx = cbuffer.data(idx_gd + 48);

    auto t_xyzz_xy = cbuffer.data(idx_gd + 49);

    auto t_xyzz_xz = cbuffer.data(idx_gd + 50);

    auto t_xyzz_yy = cbuffer.data(idx_gd + 51);

    auto t_xyzz_yz = cbuffer.data(idx_gd + 52);

    auto t_xyzz_zz = cbuffer.data(idx_gd + 53);

    auto t_xzzz_xx = cbuffer.data(idx_gd + 54);

    auto t_xzzz_xy = cbuffer.data(idx_gd + 55);

    auto t_xzzz_xz = cbuffer.data(idx_gd + 56);

    auto t_xzzz_yy = cbuffer.data(idx_gd + 57);

    auto t_xzzz_yz = cbuffer.data(idx_gd + 58);

    auto t_xzzz_zz = cbuffer.data(idx_gd + 59);

    auto t_yyyy_xx = cbuffer.data(idx_gd + 60);

    auto t_yyyy_xy = cbuffer.data(idx_gd + 61);

    auto t_yyyy_xz = cbuffer.data(idx_gd + 62);

    auto t_yyyy_yy = cbuffer.data(idx_gd + 63);

    auto t_yyyy_yz = cbuffer.data(idx_gd + 64);

    auto t_yyyy_zz = cbuffer.data(idx_gd + 65);

    auto t_yyyz_xx = cbuffer.data(idx_gd + 66);

    auto t_yyyz_xy = cbuffer.data(idx_gd + 67);

    auto t_yyyz_xz = cbuffer.data(idx_gd + 68);

    auto t_yyyz_yy = cbuffer.data(idx_gd + 69);

    auto t_yyyz_yz = cbuffer.data(idx_gd + 70);

    auto t_yyyz_zz = cbuffer.data(idx_gd + 71);

    auto t_yyzz_xx = cbuffer.data(idx_gd + 72);

    auto t_yyzz_xy = cbuffer.data(idx_gd + 73);

    auto t_yyzz_xz = cbuffer.data(idx_gd + 74);

    auto t_yyzz_yy = cbuffer.data(idx_gd + 75);

    auto t_yyzz_yz = cbuffer.data(idx_gd + 76);

    auto t_yyzz_zz = cbuffer.data(idx_gd + 77);

    auto t_yzzz_xx = cbuffer.data(idx_gd + 78);

    auto t_yzzz_xy = cbuffer.data(idx_gd + 79);

    auto t_yzzz_xz = cbuffer.data(idx_gd + 80);

    auto t_yzzz_yy = cbuffer.data(idx_gd + 81);

    auto t_yzzz_yz = cbuffer.data(idx_gd + 82);

    auto t_yzzz_zz = cbuffer.data(idx_gd + 83);

    auto t_zzzz_xx = cbuffer.data(idx_gd + 84);

    auto t_zzzz_xy = cbuffer.data(idx_gd + 85);

    auto t_zzzz_xz = cbuffer.data(idx_gd + 86);

    auto t_zzzz_yy = cbuffer.data(idx_gd + 87);

    auto t_zzzz_yz = cbuffer.data(idx_gd + 88);

    auto t_zzzz_zz = cbuffer.data(idx_gd + 89);

    // Set up components of auxiliary buffer : HD

    auto t_xxxxx_xx = cbuffer.data(idx_hd);

    auto t_xxxxx_xy = cbuffer.data(idx_hd + 1);

    auto t_xxxxx_xz = cbuffer.data(idx_hd + 2);

    auto t_xxxxx_yy = cbuffer.data(idx_hd + 3);

    auto t_xxxxx_yz = cbuffer.data(idx_hd + 4);

    auto t_xxxxx_zz = cbuffer.data(idx_hd + 5);

    auto t_xxxxy_xx = cbuffer.data(idx_hd + 6);

    auto t_xxxxy_xy = cbuffer.data(idx_hd + 7);

    auto t_xxxxy_xz = cbuffer.data(idx_hd + 8);

    auto t_xxxxy_yy = cbuffer.data(idx_hd + 9);

    auto t_xxxxy_yz = cbuffer.data(idx_hd + 10);

    auto t_xxxxy_zz = cbuffer.data(idx_hd + 11);

    auto t_xxxxz_xx = cbuffer.data(idx_hd + 12);

    auto t_xxxxz_xy = cbuffer.data(idx_hd + 13);

    auto t_xxxxz_xz = cbuffer.data(idx_hd + 14);

    auto t_xxxxz_yy = cbuffer.data(idx_hd + 15);

    auto t_xxxxz_yz = cbuffer.data(idx_hd + 16);

    auto t_xxxxz_zz = cbuffer.data(idx_hd + 17);

    auto t_xxxyy_xx = cbuffer.data(idx_hd + 18);

    auto t_xxxyy_xy = cbuffer.data(idx_hd + 19);

    auto t_xxxyy_xz = cbuffer.data(idx_hd + 20);

    auto t_xxxyy_yy = cbuffer.data(idx_hd + 21);

    auto t_xxxyy_yz = cbuffer.data(idx_hd + 22);

    auto t_xxxyy_zz = cbuffer.data(idx_hd + 23);

    auto t_xxxyz_xx = cbuffer.data(idx_hd + 24);

    auto t_xxxyz_xy = cbuffer.data(idx_hd + 25);

    auto t_xxxyz_xz = cbuffer.data(idx_hd + 26);

    auto t_xxxyz_yy = cbuffer.data(idx_hd + 27);

    auto t_xxxyz_yz = cbuffer.data(idx_hd + 28);

    auto t_xxxyz_zz = cbuffer.data(idx_hd + 29);

    auto t_xxxzz_xx = cbuffer.data(idx_hd + 30);

    auto t_xxxzz_xy = cbuffer.data(idx_hd + 31);

    auto t_xxxzz_xz = cbuffer.data(idx_hd + 32);

    auto t_xxxzz_yy = cbuffer.data(idx_hd + 33);

    auto t_xxxzz_yz = cbuffer.data(idx_hd + 34);

    auto t_xxxzz_zz = cbuffer.data(idx_hd + 35);

    auto t_xxyyy_xx = cbuffer.data(idx_hd + 36);

    auto t_xxyyy_xy = cbuffer.data(idx_hd + 37);

    auto t_xxyyy_xz = cbuffer.data(idx_hd + 38);

    auto t_xxyyy_yy = cbuffer.data(idx_hd + 39);

    auto t_xxyyy_yz = cbuffer.data(idx_hd + 40);

    auto t_xxyyy_zz = cbuffer.data(idx_hd + 41);

    auto t_xxyyz_xx = cbuffer.data(idx_hd + 42);

    auto t_xxyyz_xy = cbuffer.data(idx_hd + 43);

    auto t_xxyyz_xz = cbuffer.data(idx_hd + 44);

    auto t_xxyyz_yy = cbuffer.data(idx_hd + 45);

    auto t_xxyyz_yz = cbuffer.data(idx_hd + 46);

    auto t_xxyyz_zz = cbuffer.data(idx_hd + 47);

    auto t_xxyzz_xx = cbuffer.data(idx_hd + 48);

    auto t_xxyzz_xy = cbuffer.data(idx_hd + 49);

    auto t_xxyzz_xz = cbuffer.data(idx_hd + 50);

    auto t_xxyzz_yy = cbuffer.data(idx_hd + 51);

    auto t_xxyzz_yz = cbuffer.data(idx_hd + 52);

    auto t_xxyzz_zz = cbuffer.data(idx_hd + 53);

    auto t_xxzzz_xx = cbuffer.data(idx_hd + 54);

    auto t_xxzzz_xy = cbuffer.data(idx_hd + 55);

    auto t_xxzzz_xz = cbuffer.data(idx_hd + 56);

    auto t_xxzzz_yy = cbuffer.data(idx_hd + 57);

    auto t_xxzzz_yz = cbuffer.data(idx_hd + 58);

    auto t_xxzzz_zz = cbuffer.data(idx_hd + 59);

    auto t_xyyyy_xx = cbuffer.data(idx_hd + 60);

    auto t_xyyyy_xy = cbuffer.data(idx_hd + 61);

    auto t_xyyyy_xz = cbuffer.data(idx_hd + 62);

    auto t_xyyyy_yy = cbuffer.data(idx_hd + 63);

    auto t_xyyyy_yz = cbuffer.data(idx_hd + 64);

    auto t_xyyyy_zz = cbuffer.data(idx_hd + 65);

    auto t_xyyyz_xx = cbuffer.data(idx_hd + 66);

    auto t_xyyyz_xy = cbuffer.data(idx_hd + 67);

    auto t_xyyyz_xz = cbuffer.data(idx_hd + 68);

    auto t_xyyyz_yy = cbuffer.data(idx_hd + 69);

    auto t_xyyyz_yz = cbuffer.data(idx_hd + 70);

    auto t_xyyyz_zz = cbuffer.data(idx_hd + 71);

    auto t_xyyzz_xx = cbuffer.data(idx_hd + 72);

    auto t_xyyzz_xy = cbuffer.data(idx_hd + 73);

    auto t_xyyzz_xz = cbuffer.data(idx_hd + 74);

    auto t_xyyzz_yy = cbuffer.data(idx_hd + 75);

    auto t_xyyzz_yz = cbuffer.data(idx_hd + 76);

    auto t_xyyzz_zz = cbuffer.data(idx_hd + 77);

    auto t_xyzzz_xx = cbuffer.data(idx_hd + 78);

    auto t_xyzzz_xy = cbuffer.data(idx_hd + 79);

    auto t_xyzzz_xz = cbuffer.data(idx_hd + 80);

    auto t_xyzzz_yy = cbuffer.data(idx_hd + 81);

    auto t_xyzzz_yz = cbuffer.data(idx_hd + 82);

    auto t_xyzzz_zz = cbuffer.data(idx_hd + 83);

    auto t_xzzzz_xx = cbuffer.data(idx_hd + 84);

    auto t_xzzzz_xy = cbuffer.data(idx_hd + 85);

    auto t_xzzzz_xz = cbuffer.data(idx_hd + 86);

    auto t_xzzzz_yy = cbuffer.data(idx_hd + 87);

    auto t_xzzzz_yz = cbuffer.data(idx_hd + 88);

    auto t_xzzzz_zz = cbuffer.data(idx_hd + 89);

    auto t_yyyyy_yy = cbuffer.data(idx_hd + 93);

    auto t_yyyyy_yz = cbuffer.data(idx_hd + 94);

    auto t_yyyyy_zz = cbuffer.data(idx_hd + 95);

    auto t_yyyyz_yy = cbuffer.data(idx_hd + 99);

    auto t_yyyyz_yz = cbuffer.data(idx_hd + 100);

    auto t_yyyyz_zz = cbuffer.data(idx_hd + 101);

    auto t_yyyzz_yy = cbuffer.data(idx_hd + 105);

    auto t_yyyzz_yz = cbuffer.data(idx_hd + 106);

    auto t_yyyzz_zz = cbuffer.data(idx_hd + 107);

    auto t_yyzzz_yy = cbuffer.data(idx_hd + 111);

    auto t_yyzzz_yz = cbuffer.data(idx_hd + 112);

    auto t_yyzzz_zz = cbuffer.data(idx_hd + 113);

    auto t_yzzzz_yy = cbuffer.data(idx_hd + 117);

    auto t_yzzzz_yz = cbuffer.data(idx_hd + 118);

    auto t_yzzzz_zz = cbuffer.data(idx_hd + 119);

    auto t_zzzzz_zz = cbuffer.data(idx_hd + 125);

    // Set up components of targeted buffer : GF

    auto t_xxxx_xxx = cbuffer.data(idx_gf);

    auto t_xxxx_xxy = cbuffer.data(idx_gf + 1);

    auto t_xxxx_xxz = cbuffer.data(idx_gf + 2);

    auto t_xxxx_xyy = cbuffer.data(idx_gf + 3);

    auto t_xxxx_xyz = cbuffer.data(idx_gf + 4);

    auto t_xxxx_xzz = cbuffer.data(idx_gf + 5);

    auto t_xxxx_yyy = cbuffer.data(idx_gf + 6);

    auto t_xxxx_yyz = cbuffer.data(idx_gf + 7);

    auto t_xxxx_yzz = cbuffer.data(idx_gf + 8);

    auto t_xxxx_zzz = cbuffer.data(idx_gf + 9);

    auto t_xxxy_xxx = cbuffer.data(idx_gf + 10);

    auto t_xxxy_xxy = cbuffer.data(idx_gf + 11);

    auto t_xxxy_xxz = cbuffer.data(idx_gf + 12);

    auto t_xxxy_xyy = cbuffer.data(idx_gf + 13);

    auto t_xxxy_xyz = cbuffer.data(idx_gf + 14);

    auto t_xxxy_xzz = cbuffer.data(idx_gf + 15);

    auto t_xxxy_yyy = cbuffer.data(idx_gf + 16);

    auto t_xxxy_yyz = cbuffer.data(idx_gf + 17);

    auto t_xxxy_yzz = cbuffer.data(idx_gf + 18);

    auto t_xxxy_zzz = cbuffer.data(idx_gf + 19);

    auto t_xxxz_xxx = cbuffer.data(idx_gf + 20);

    auto t_xxxz_xxy = cbuffer.data(idx_gf + 21);

    auto t_xxxz_xxz = cbuffer.data(idx_gf + 22);

    auto t_xxxz_xyy = cbuffer.data(idx_gf + 23);

    auto t_xxxz_xyz = cbuffer.data(idx_gf + 24);

    auto t_xxxz_xzz = cbuffer.data(idx_gf + 25);

    auto t_xxxz_yyy = cbuffer.data(idx_gf + 26);

    auto t_xxxz_yyz = cbuffer.data(idx_gf + 27);

    auto t_xxxz_yzz = cbuffer.data(idx_gf + 28);

    auto t_xxxz_zzz = cbuffer.data(idx_gf + 29);

    auto t_xxyy_xxx = cbuffer.data(idx_gf + 30);

    auto t_xxyy_xxy = cbuffer.data(idx_gf + 31);

    auto t_xxyy_xxz = cbuffer.data(idx_gf + 32);

    auto t_xxyy_xyy = cbuffer.data(idx_gf + 33);

    auto t_xxyy_xyz = cbuffer.data(idx_gf + 34);

    auto t_xxyy_xzz = cbuffer.data(idx_gf + 35);

    auto t_xxyy_yyy = cbuffer.data(idx_gf + 36);

    auto t_xxyy_yyz = cbuffer.data(idx_gf + 37);

    auto t_xxyy_yzz = cbuffer.data(idx_gf + 38);

    auto t_xxyy_zzz = cbuffer.data(idx_gf + 39);

    auto t_xxyz_xxx = cbuffer.data(idx_gf + 40);

    auto t_xxyz_xxy = cbuffer.data(idx_gf + 41);

    auto t_xxyz_xxz = cbuffer.data(idx_gf + 42);

    auto t_xxyz_xyy = cbuffer.data(idx_gf + 43);

    auto t_xxyz_xyz = cbuffer.data(idx_gf + 44);

    auto t_xxyz_xzz = cbuffer.data(idx_gf + 45);

    auto t_xxyz_yyy = cbuffer.data(idx_gf + 46);

    auto t_xxyz_yyz = cbuffer.data(idx_gf + 47);

    auto t_xxyz_yzz = cbuffer.data(idx_gf + 48);

    auto t_xxyz_zzz = cbuffer.data(idx_gf + 49);

    auto t_xxzz_xxx = cbuffer.data(idx_gf + 50);

    auto t_xxzz_xxy = cbuffer.data(idx_gf + 51);

    auto t_xxzz_xxz = cbuffer.data(idx_gf + 52);

    auto t_xxzz_xyy = cbuffer.data(idx_gf + 53);

    auto t_xxzz_xyz = cbuffer.data(idx_gf + 54);

    auto t_xxzz_xzz = cbuffer.data(idx_gf + 55);

    auto t_xxzz_yyy = cbuffer.data(idx_gf + 56);

    auto t_xxzz_yyz = cbuffer.data(idx_gf + 57);

    auto t_xxzz_yzz = cbuffer.data(idx_gf + 58);

    auto t_xxzz_zzz = cbuffer.data(idx_gf + 59);

    auto t_xyyy_xxx = cbuffer.data(idx_gf + 60);

    auto t_xyyy_xxy = cbuffer.data(idx_gf + 61);

    auto t_xyyy_xxz = cbuffer.data(idx_gf + 62);

    auto t_xyyy_xyy = cbuffer.data(idx_gf + 63);

    auto t_xyyy_xyz = cbuffer.data(idx_gf + 64);

    auto t_xyyy_xzz = cbuffer.data(idx_gf + 65);

    auto t_xyyy_yyy = cbuffer.data(idx_gf + 66);

    auto t_xyyy_yyz = cbuffer.data(idx_gf + 67);

    auto t_xyyy_yzz = cbuffer.data(idx_gf + 68);

    auto t_xyyy_zzz = cbuffer.data(idx_gf + 69);

    auto t_xyyz_xxx = cbuffer.data(idx_gf + 70);

    auto t_xyyz_xxy = cbuffer.data(idx_gf + 71);

    auto t_xyyz_xxz = cbuffer.data(idx_gf + 72);

    auto t_xyyz_xyy = cbuffer.data(idx_gf + 73);

    auto t_xyyz_xyz = cbuffer.data(idx_gf + 74);

    auto t_xyyz_xzz = cbuffer.data(idx_gf + 75);

    auto t_xyyz_yyy = cbuffer.data(idx_gf + 76);

    auto t_xyyz_yyz = cbuffer.data(idx_gf + 77);

    auto t_xyyz_yzz = cbuffer.data(idx_gf + 78);

    auto t_xyyz_zzz = cbuffer.data(idx_gf + 79);

    auto t_xyzz_xxx = cbuffer.data(idx_gf + 80);

    auto t_xyzz_xxy = cbuffer.data(idx_gf + 81);

    auto t_xyzz_xxz = cbuffer.data(idx_gf + 82);

    auto t_xyzz_xyy = cbuffer.data(idx_gf + 83);

    auto t_xyzz_xyz = cbuffer.data(idx_gf + 84);

    auto t_xyzz_xzz = cbuffer.data(idx_gf + 85);

    auto t_xyzz_yyy = cbuffer.data(idx_gf + 86);

    auto t_xyzz_yyz = cbuffer.data(idx_gf + 87);

    auto t_xyzz_yzz = cbuffer.data(idx_gf + 88);

    auto t_xyzz_zzz = cbuffer.data(idx_gf + 89);

    auto t_xzzz_xxx = cbuffer.data(idx_gf + 90);

    auto t_xzzz_xxy = cbuffer.data(idx_gf + 91);

    auto t_xzzz_xxz = cbuffer.data(idx_gf + 92);

    auto t_xzzz_xyy = cbuffer.data(idx_gf + 93);

    auto t_xzzz_xyz = cbuffer.data(idx_gf + 94);

    auto t_xzzz_xzz = cbuffer.data(idx_gf + 95);

    auto t_xzzz_yyy = cbuffer.data(idx_gf + 96);

    auto t_xzzz_yyz = cbuffer.data(idx_gf + 97);

    auto t_xzzz_yzz = cbuffer.data(idx_gf + 98);

    auto t_xzzz_zzz = cbuffer.data(idx_gf + 99);

    auto t_yyyy_xxx = cbuffer.data(idx_gf + 100);

    auto t_yyyy_xxy = cbuffer.data(idx_gf + 101);

    auto t_yyyy_xxz = cbuffer.data(idx_gf + 102);

    auto t_yyyy_xyy = cbuffer.data(idx_gf + 103);

    auto t_yyyy_xyz = cbuffer.data(idx_gf + 104);

    auto t_yyyy_xzz = cbuffer.data(idx_gf + 105);

    auto t_yyyy_yyy = cbuffer.data(idx_gf + 106);

    auto t_yyyy_yyz = cbuffer.data(idx_gf + 107);

    auto t_yyyy_yzz = cbuffer.data(idx_gf + 108);

    auto t_yyyy_zzz = cbuffer.data(idx_gf + 109);

    auto t_yyyz_xxx = cbuffer.data(idx_gf + 110);

    auto t_yyyz_xxy = cbuffer.data(idx_gf + 111);

    auto t_yyyz_xxz = cbuffer.data(idx_gf + 112);

    auto t_yyyz_xyy = cbuffer.data(idx_gf + 113);

    auto t_yyyz_xyz = cbuffer.data(idx_gf + 114);

    auto t_yyyz_xzz = cbuffer.data(idx_gf + 115);

    auto t_yyyz_yyy = cbuffer.data(idx_gf + 116);

    auto t_yyyz_yyz = cbuffer.data(idx_gf + 117);

    auto t_yyyz_yzz = cbuffer.data(idx_gf + 118);

    auto t_yyyz_zzz = cbuffer.data(idx_gf + 119);

    auto t_yyzz_xxx = cbuffer.data(idx_gf + 120);

    auto t_yyzz_xxy = cbuffer.data(idx_gf + 121);

    auto t_yyzz_xxz = cbuffer.data(idx_gf + 122);

    auto t_yyzz_xyy = cbuffer.data(idx_gf + 123);

    auto t_yyzz_xyz = cbuffer.data(idx_gf + 124);

    auto t_yyzz_xzz = cbuffer.data(idx_gf + 125);

    auto t_yyzz_yyy = cbuffer.data(idx_gf + 126);

    auto t_yyzz_yyz = cbuffer.data(idx_gf + 127);

    auto t_yyzz_yzz = cbuffer.data(idx_gf + 128);

    auto t_yyzz_zzz = cbuffer.data(idx_gf + 129);

    auto t_yzzz_xxx = cbuffer.data(idx_gf + 130);

    auto t_yzzz_xxy = cbuffer.data(idx_gf + 131);

    auto t_yzzz_xxz = cbuffer.data(idx_gf + 132);

    auto t_yzzz_xyy = cbuffer.data(idx_gf + 133);

    auto t_yzzz_xyz = cbuffer.data(idx_gf + 134);

    auto t_yzzz_xzz = cbuffer.data(idx_gf + 135);

    auto t_yzzz_yyy = cbuffer.data(idx_gf + 136);

    auto t_yzzz_yyz = cbuffer.data(idx_gf + 137);

    auto t_yzzz_yzz = cbuffer.data(idx_gf + 138);

    auto t_yzzz_zzz = cbuffer.data(idx_gf + 139);

    auto t_zzzz_xxx = cbuffer.data(idx_gf + 140);

    auto t_zzzz_xxy = cbuffer.data(idx_gf + 141);

    auto t_zzzz_xxz = cbuffer.data(idx_gf + 142);

    auto t_zzzz_xyy = cbuffer.data(idx_gf + 143);

    auto t_zzzz_xyz = cbuffer.data(idx_gf + 144);

    auto t_zzzz_xzz = cbuffer.data(idx_gf + 145);

    auto t_zzzz_yyy = cbuffer.data(idx_gf + 146);

    auto t_zzzz_yyz = cbuffer.data(idx_gf + 147);

    auto t_zzzz_yzz = cbuffer.data(idx_gf + 148);

    auto t_zzzz_zzz = cbuffer.data(idx_gf + 149);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxxx_xx, t_xxxx_xxx, t_xxxx_xxy, t_xxxx_xxz, t_xxxx_xy, t_xxxx_xyy, t_xxxx_xyz, t_xxxx_xz, t_xxxx_xzz, t_xxxx_yy, t_xxxx_yyy, t_xxxx_yyz, t_xxxx_yz, t_xxxx_yzz, t_xxxx_zz, t_xxxx_zzz, t_xxxxx_xx, t_xxxxx_xy, t_xxxxx_xz, t_xxxxx_yy, t_xxxxx_yz, t_xxxxx_zz, t_xxxxy_xx, t_xxxxy_xy, t_xxxxy_xz, t_xxxxy_yy, t_xxxxy_yz, t_xxxxy_zz, t_xxxxz_xx, t_xxxxz_xy, t_xxxxz_xz, t_xxxxz_yy, t_xxxxz_yz, t_xxxxz_zz, t_xxxy_xx, t_xxxy_xxx, t_xxxy_xxy, t_xxxy_xxz, t_xxxy_xy, t_xxxy_xyy, t_xxxy_xyz, t_xxxy_xz, t_xxxy_xzz, t_xxxy_yy, t_xxxy_yyy, t_xxxy_yyz, t_xxxy_yz, t_xxxy_yzz, t_xxxy_zz, t_xxxy_zzz, t_xxxyy_xx, t_xxxyy_xy, t_xxxyy_xz, t_xxxyy_yy, t_xxxyy_yz, t_xxxyy_zz, t_xxxyz_xx, t_xxxyz_xy, t_xxxyz_xz, t_xxxyz_yy, t_xxxyz_yz, t_xxxyz_zz, t_xxxz_xx, t_xxxz_xxx, t_xxxz_xxy, t_xxxz_xxz, t_xxxz_xy, t_xxxz_xyy, t_xxxz_xyz, t_xxxz_xz, t_xxxz_xzz, t_xxxz_yy, t_xxxz_yyy, t_xxxz_yyz, t_xxxz_yz, t_xxxz_yzz, t_xxxz_zz, t_xxxz_zzz, t_xxxzz_xx, t_xxxzz_xy, t_xxxzz_xz, t_xxxzz_yy, t_xxxzz_yz, t_xxxzz_zz, t_xxyy_xx, t_xxyy_xxx, t_xxyy_xxy, t_xxyy_xxz, t_xxyy_xy, t_xxyy_xyy, t_xxyy_xyz, t_xxyy_xz, t_xxyy_xzz, t_xxyy_yy, t_xxyy_yyy, t_xxyy_yyz, t_xxyy_yz, t_xxyy_yzz, t_xxyy_zz, t_xxyy_zzz, t_xxyyy_xx, t_xxyyy_xy, t_xxyyy_xz, t_xxyyy_yy, t_xxyyy_yz, t_xxyyy_zz, t_xxyyz_xx, t_xxyyz_xy, t_xxyyz_xz, t_xxyyz_yy, t_xxyyz_yz, t_xxyyz_zz, t_xxyz_xx, t_xxyz_xxx, t_xxyz_xxy, t_xxyz_xxz, t_xxyz_xy, t_xxyz_xyy, t_xxyz_xyz, t_xxyz_xz, t_xxyz_xzz, t_xxyz_yy, t_xxyz_yyy, t_xxyz_yyz, t_xxyz_yz, t_xxyz_yzz, t_xxyz_zz, t_xxyz_zzz, t_xxyzz_xx, t_xxyzz_xy, t_xxyzz_xz, t_xxyzz_yy, t_xxyzz_yz, t_xxyzz_zz, t_xxzz_xx, t_xxzz_xxx, t_xxzz_xxy, t_xxzz_xxz, t_xxzz_xy, t_xxzz_xyy, t_xxzz_xyz, t_xxzz_xz, t_xxzz_xzz, t_xxzz_yy, t_xxzz_yyy, t_xxzz_yyz, t_xxzz_yz, t_xxzz_yzz, t_xxzz_zz, t_xxzz_zzz, t_xxzzz_xx, t_xxzzz_xy, t_xxzzz_xz, t_xxzzz_yy, t_xxzzz_yz, t_xxzzz_zz, t_xyyy_xx, t_xyyy_xxx, t_xyyy_xxy, t_xyyy_xxz, t_xyyy_xy, t_xyyy_xyy, t_xyyy_xyz, t_xyyy_xz, t_xyyy_xzz, t_xyyy_yy, t_xyyy_yyy, t_xyyy_yyz, t_xyyy_yz, t_xyyy_yzz, t_xyyy_zz, t_xyyy_zzz, t_xyyyy_xx, t_xyyyy_xy, t_xyyyy_xz, t_xyyyy_yy, t_xyyyy_yz, t_xyyyy_zz, t_xyyyz_xx, t_xyyyz_xy, t_xyyyz_xz, t_xyyyz_yy, t_xyyyz_yz, t_xyyyz_zz, t_xyyz_xx, t_xyyz_xxx, t_xyyz_xxy, t_xyyz_xxz, t_xyyz_xy, t_xyyz_xyy, t_xyyz_xyz, t_xyyz_xz, t_xyyz_xzz, t_xyyz_yy, t_xyyz_yyy, t_xyyz_yyz, t_xyyz_yz, t_xyyz_yzz, t_xyyz_zz, t_xyyz_zzz, t_xyyzz_xx, t_xyyzz_xy, t_xyyzz_xz, t_xyyzz_yy, t_xyyzz_yz, t_xyyzz_zz, t_xyzz_xx, t_xyzz_xxx, t_xyzz_xxy, t_xyzz_xxz, t_xyzz_xy, t_xyzz_xyy, t_xyzz_xyz, t_xyzz_xz, t_xyzz_xzz, t_xyzz_yy, t_xyzz_yyy, t_xyzz_yyz, t_xyzz_yz, t_xyzz_yzz, t_xyzz_zz, t_xyzz_zzz, t_xyzzz_xx, t_xyzzz_xy, t_xyzzz_xz, t_xyzzz_yy, t_xyzzz_yz, t_xyzzz_zz, t_xzzz_xx, t_xzzz_xxx, t_xzzz_xxy, t_xzzz_xxz, t_xzzz_xy, t_xzzz_xyy, t_xzzz_xyz, t_xzzz_xz, t_xzzz_xzz, t_xzzz_yy, t_xzzz_yyy, t_xzzz_yyz, t_xzzz_yz, t_xzzz_yzz, t_xzzz_zz, t_xzzz_zzz, t_xzzzz_xx, t_xzzzz_xy, t_xzzzz_xz, t_xzzzz_yy, t_xzzzz_yz, t_xzzzz_zz, t_yyyy_xx, t_yyyy_xxx, t_yyyy_xxy, t_yyyy_xxz, t_yyyy_xy, t_yyyy_xyy, t_yyyy_xyz, t_yyyy_xz, t_yyyy_xzz, t_yyyy_yy, t_yyyy_yyy, t_yyyy_yyz, t_yyyy_yz, t_yyyy_yzz, t_yyyy_zz, t_yyyy_zzz, t_yyyyy_yy, t_yyyyy_yz, t_yyyyy_zz, t_yyyyz_yy, t_yyyyz_yz, t_yyyyz_zz, t_yyyz_xx, t_yyyz_xxx, t_yyyz_xxy, t_yyyz_xxz, t_yyyz_xy, t_yyyz_xyy, t_yyyz_xyz, t_yyyz_xz, t_yyyz_xzz, t_yyyz_yy, t_yyyz_yyy, t_yyyz_yyz, t_yyyz_yz, t_yyyz_yzz, t_yyyz_zz, t_yyyz_zzz, t_yyyzz_yy, t_yyyzz_yz, t_yyyzz_zz, t_yyzz_xx, t_yyzz_xxx, t_yyzz_xxy, t_yyzz_xxz, t_yyzz_xy, t_yyzz_xyy, t_yyzz_xyz, t_yyzz_xz, t_yyzz_xzz, t_yyzz_yy, t_yyzz_yyy, t_yyzz_yyz, t_yyzz_yz, t_yyzz_yzz, t_yyzz_zz, t_yyzz_zzz, t_yyzzz_yy, t_yyzzz_yz, t_yyzzz_zz, t_yzzz_xx, t_yzzz_xxx, t_yzzz_xxy, t_yzzz_xxz, t_yzzz_xy, t_yzzz_xyy, t_yzzz_xyz, t_yzzz_xz, t_yzzz_xzz, t_yzzz_yy, t_yzzz_yyy, t_yzzz_yyz, t_yzzz_yz, t_yzzz_yzz, t_yzzz_zz, t_yzzz_zzz, t_yzzzz_yy, t_yzzzz_yz, t_yzzzz_zz, t_zzzz_xx, t_zzzz_xxx, t_zzzz_xxy, t_zzzz_xxz, t_zzzz_xy, t_zzzz_xyy, t_zzzz_xyz, t_zzzz_xz, t_zzzz_xzz, t_zzzz_yy, t_zzzz_yyy, t_zzzz_yyz, t_zzzz_yz, t_zzzz_yzz, t_zzzz_zz, t_zzzz_zzz, t_zzzzz_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxx_xxx[i] = t_xxxx_xx[i] * ab_x[i] + t_xxxxx_xx[i];

        t_xxxx_xxy[i] = t_xxxx_xy[i] * ab_x[i] + t_xxxxx_xy[i];

        t_xxxx_xxz[i] = t_xxxx_xz[i] * ab_x[i] + t_xxxxx_xz[i];

        t_xxxx_xyy[i] = t_xxxx_yy[i] * ab_x[i] + t_xxxxx_yy[i];

        t_xxxx_xyz[i] = t_xxxx_yz[i] * ab_x[i] + t_xxxxx_yz[i];

        t_xxxx_xzz[i] = t_xxxx_zz[i] * ab_x[i] + t_xxxxx_zz[i];

        t_xxxx_yyy[i] = t_xxxx_yy[i] * ab_y[i] + t_xxxxy_yy[i];

        t_xxxx_yyz[i] = t_xxxx_yz[i] * ab_y[i] + t_xxxxy_yz[i];

        t_xxxx_yzz[i] = t_xxxx_zz[i] * ab_y[i] + t_xxxxy_zz[i];

        t_xxxx_zzz[i] = t_xxxx_zz[i] * ab_z[i] + t_xxxxz_zz[i];

        t_xxxy_xxx[i] = t_xxxy_xx[i] * ab_x[i] + t_xxxxy_xx[i];

        t_xxxy_xxy[i] = t_xxxy_xy[i] * ab_x[i] + t_xxxxy_xy[i];

        t_xxxy_xxz[i] = t_xxxy_xz[i] * ab_x[i] + t_xxxxy_xz[i];

        t_xxxy_xyy[i] = t_xxxy_yy[i] * ab_x[i] + t_xxxxy_yy[i];

        t_xxxy_xyz[i] = t_xxxy_yz[i] * ab_x[i] + t_xxxxy_yz[i];

        t_xxxy_xzz[i] = t_xxxy_zz[i] * ab_x[i] + t_xxxxy_zz[i];

        t_xxxy_yyy[i] = t_xxxy_yy[i] * ab_y[i] + t_xxxyy_yy[i];

        t_xxxy_yyz[i] = t_xxxy_yz[i] * ab_y[i] + t_xxxyy_yz[i];

        t_xxxy_yzz[i] = t_xxxy_zz[i] * ab_y[i] + t_xxxyy_zz[i];

        t_xxxy_zzz[i] = t_xxxy_zz[i] * ab_z[i] + t_xxxyz_zz[i];

        t_xxxz_xxx[i] = t_xxxz_xx[i] * ab_x[i] + t_xxxxz_xx[i];

        t_xxxz_xxy[i] = t_xxxz_xy[i] * ab_x[i] + t_xxxxz_xy[i];

        t_xxxz_xxz[i] = t_xxxz_xz[i] * ab_x[i] + t_xxxxz_xz[i];

        t_xxxz_xyy[i] = t_xxxz_yy[i] * ab_x[i] + t_xxxxz_yy[i];

        t_xxxz_xyz[i] = t_xxxz_yz[i] * ab_x[i] + t_xxxxz_yz[i];

        t_xxxz_xzz[i] = t_xxxz_zz[i] * ab_x[i] + t_xxxxz_zz[i];

        t_xxxz_yyy[i] = t_xxxz_yy[i] * ab_y[i] + t_xxxyz_yy[i];

        t_xxxz_yyz[i] = t_xxxz_yz[i] * ab_y[i] + t_xxxyz_yz[i];

        t_xxxz_yzz[i] = t_xxxz_zz[i] * ab_y[i] + t_xxxyz_zz[i];

        t_xxxz_zzz[i] = t_xxxz_zz[i] * ab_z[i] + t_xxxzz_zz[i];

        t_xxyy_xxx[i] = t_xxyy_xx[i] * ab_x[i] + t_xxxyy_xx[i];

        t_xxyy_xxy[i] = t_xxyy_xy[i] * ab_x[i] + t_xxxyy_xy[i];

        t_xxyy_xxz[i] = t_xxyy_xz[i] * ab_x[i] + t_xxxyy_xz[i];

        t_xxyy_xyy[i] = t_xxyy_yy[i] * ab_x[i] + t_xxxyy_yy[i];

        t_xxyy_xyz[i] = t_xxyy_yz[i] * ab_x[i] + t_xxxyy_yz[i];

        t_xxyy_xzz[i] = t_xxyy_zz[i] * ab_x[i] + t_xxxyy_zz[i];

        t_xxyy_yyy[i] = t_xxyy_yy[i] * ab_y[i] + t_xxyyy_yy[i];

        t_xxyy_yyz[i] = t_xxyy_yz[i] * ab_y[i] + t_xxyyy_yz[i];

        t_xxyy_yzz[i] = t_xxyy_zz[i] * ab_y[i] + t_xxyyy_zz[i];

        t_xxyy_zzz[i] = t_xxyy_zz[i] * ab_z[i] + t_xxyyz_zz[i];

        t_xxyz_xxx[i] = t_xxyz_xx[i] * ab_x[i] + t_xxxyz_xx[i];

        t_xxyz_xxy[i] = t_xxyz_xy[i] * ab_x[i] + t_xxxyz_xy[i];

        t_xxyz_xxz[i] = t_xxyz_xz[i] * ab_x[i] + t_xxxyz_xz[i];

        t_xxyz_xyy[i] = t_xxyz_yy[i] * ab_x[i] + t_xxxyz_yy[i];

        t_xxyz_xyz[i] = t_xxyz_yz[i] * ab_x[i] + t_xxxyz_yz[i];

        t_xxyz_xzz[i] = t_xxyz_zz[i] * ab_x[i] + t_xxxyz_zz[i];

        t_xxyz_yyy[i] = t_xxyz_yy[i] * ab_y[i] + t_xxyyz_yy[i];

        t_xxyz_yyz[i] = t_xxyz_yz[i] * ab_y[i] + t_xxyyz_yz[i];

        t_xxyz_yzz[i] = t_xxyz_zz[i] * ab_y[i] + t_xxyyz_zz[i];

        t_xxyz_zzz[i] = t_xxyz_zz[i] * ab_z[i] + t_xxyzz_zz[i];

        t_xxzz_xxx[i] = t_xxzz_xx[i] * ab_x[i] + t_xxxzz_xx[i];

        t_xxzz_xxy[i] = t_xxzz_xy[i] * ab_x[i] + t_xxxzz_xy[i];

        t_xxzz_xxz[i] = t_xxzz_xz[i] * ab_x[i] + t_xxxzz_xz[i];

        t_xxzz_xyy[i] = t_xxzz_yy[i] * ab_x[i] + t_xxxzz_yy[i];

        t_xxzz_xyz[i] = t_xxzz_yz[i] * ab_x[i] + t_xxxzz_yz[i];

        t_xxzz_xzz[i] = t_xxzz_zz[i] * ab_x[i] + t_xxxzz_zz[i];

        t_xxzz_yyy[i] = t_xxzz_yy[i] * ab_y[i] + t_xxyzz_yy[i];

        t_xxzz_yyz[i] = t_xxzz_yz[i] * ab_y[i] + t_xxyzz_yz[i];

        t_xxzz_yzz[i] = t_xxzz_zz[i] * ab_y[i] + t_xxyzz_zz[i];

        t_xxzz_zzz[i] = t_xxzz_zz[i] * ab_z[i] + t_xxzzz_zz[i];

        t_xyyy_xxx[i] = t_xyyy_xx[i] * ab_x[i] + t_xxyyy_xx[i];

        t_xyyy_xxy[i] = t_xyyy_xy[i] * ab_x[i] + t_xxyyy_xy[i];

        t_xyyy_xxz[i] = t_xyyy_xz[i] * ab_x[i] + t_xxyyy_xz[i];

        t_xyyy_xyy[i] = t_xyyy_yy[i] * ab_x[i] + t_xxyyy_yy[i];

        t_xyyy_xyz[i] = t_xyyy_yz[i] * ab_x[i] + t_xxyyy_yz[i];

        t_xyyy_xzz[i] = t_xyyy_zz[i] * ab_x[i] + t_xxyyy_zz[i];

        t_xyyy_yyy[i] = t_xyyy_yy[i] * ab_y[i] + t_xyyyy_yy[i];

        t_xyyy_yyz[i] = t_xyyy_yz[i] * ab_y[i] + t_xyyyy_yz[i];

        t_xyyy_yzz[i] = t_xyyy_zz[i] * ab_y[i] + t_xyyyy_zz[i];

        t_xyyy_zzz[i] = t_xyyy_zz[i] * ab_z[i] + t_xyyyz_zz[i];

        t_xyyz_xxx[i] = t_xyyz_xx[i] * ab_x[i] + t_xxyyz_xx[i];

        t_xyyz_xxy[i] = t_xyyz_xy[i] * ab_x[i] + t_xxyyz_xy[i];

        t_xyyz_xxz[i] = t_xyyz_xz[i] * ab_x[i] + t_xxyyz_xz[i];

        t_xyyz_xyy[i] = t_xyyz_yy[i] * ab_x[i] + t_xxyyz_yy[i];

        t_xyyz_xyz[i] = t_xyyz_yz[i] * ab_x[i] + t_xxyyz_yz[i];

        t_xyyz_xzz[i] = t_xyyz_zz[i] * ab_x[i] + t_xxyyz_zz[i];

        t_xyyz_yyy[i] = t_xyyz_yy[i] * ab_y[i] + t_xyyyz_yy[i];

        t_xyyz_yyz[i] = t_xyyz_yz[i] * ab_y[i] + t_xyyyz_yz[i];

        t_xyyz_yzz[i] = t_xyyz_zz[i] * ab_y[i] + t_xyyyz_zz[i];

        t_xyyz_zzz[i] = t_xyyz_zz[i] * ab_z[i] + t_xyyzz_zz[i];

        t_xyzz_xxx[i] = t_xyzz_xx[i] * ab_x[i] + t_xxyzz_xx[i];

        t_xyzz_xxy[i] = t_xyzz_xy[i] * ab_x[i] + t_xxyzz_xy[i];

        t_xyzz_xxz[i] = t_xyzz_xz[i] * ab_x[i] + t_xxyzz_xz[i];

        t_xyzz_xyy[i] = t_xyzz_yy[i] * ab_x[i] + t_xxyzz_yy[i];

        t_xyzz_xyz[i] = t_xyzz_yz[i] * ab_x[i] + t_xxyzz_yz[i];

        t_xyzz_xzz[i] = t_xyzz_zz[i] * ab_x[i] + t_xxyzz_zz[i];

        t_xyzz_yyy[i] = t_xyzz_yy[i] * ab_y[i] + t_xyyzz_yy[i];

        t_xyzz_yyz[i] = t_xyzz_yz[i] * ab_y[i] + t_xyyzz_yz[i];

        t_xyzz_yzz[i] = t_xyzz_zz[i] * ab_y[i] + t_xyyzz_zz[i];

        t_xyzz_zzz[i] = t_xyzz_zz[i] * ab_z[i] + t_xyzzz_zz[i];

        t_xzzz_xxx[i] = t_xzzz_xx[i] * ab_x[i] + t_xxzzz_xx[i];

        t_xzzz_xxy[i] = t_xzzz_xy[i] * ab_x[i] + t_xxzzz_xy[i];

        t_xzzz_xxz[i] = t_xzzz_xz[i] * ab_x[i] + t_xxzzz_xz[i];

        t_xzzz_xyy[i] = t_xzzz_yy[i] * ab_x[i] + t_xxzzz_yy[i];

        t_xzzz_xyz[i] = t_xzzz_yz[i] * ab_x[i] + t_xxzzz_yz[i];

        t_xzzz_xzz[i] = t_xzzz_zz[i] * ab_x[i] + t_xxzzz_zz[i];

        t_xzzz_yyy[i] = t_xzzz_yy[i] * ab_y[i] + t_xyzzz_yy[i];

        t_xzzz_yyz[i] = t_xzzz_yz[i] * ab_y[i] + t_xyzzz_yz[i];

        t_xzzz_yzz[i] = t_xzzz_zz[i] * ab_y[i] + t_xyzzz_zz[i];

        t_xzzz_zzz[i] = t_xzzz_zz[i] * ab_z[i] + t_xzzzz_zz[i];

        t_yyyy_xxx[i] = t_yyyy_xx[i] * ab_x[i] + t_xyyyy_xx[i];

        t_yyyy_xxy[i] = t_yyyy_xy[i] * ab_x[i] + t_xyyyy_xy[i];

        t_yyyy_xxz[i] = t_yyyy_xz[i] * ab_x[i] + t_xyyyy_xz[i];

        t_yyyy_xyy[i] = t_yyyy_yy[i] * ab_x[i] + t_xyyyy_yy[i];

        t_yyyy_xyz[i] = t_yyyy_yz[i] * ab_x[i] + t_xyyyy_yz[i];

        t_yyyy_xzz[i] = t_yyyy_zz[i] * ab_x[i] + t_xyyyy_zz[i];

        t_yyyy_yyy[i] = t_yyyy_yy[i] * ab_y[i] + t_yyyyy_yy[i];

        t_yyyy_yyz[i] = t_yyyy_yz[i] * ab_y[i] + t_yyyyy_yz[i];

        t_yyyy_yzz[i] = t_yyyy_zz[i] * ab_y[i] + t_yyyyy_zz[i];

        t_yyyy_zzz[i] = t_yyyy_zz[i] * ab_z[i] + t_yyyyz_zz[i];

        t_yyyz_xxx[i] = t_yyyz_xx[i] * ab_x[i] + t_xyyyz_xx[i];

        t_yyyz_xxy[i] = t_yyyz_xy[i] * ab_x[i] + t_xyyyz_xy[i];

        t_yyyz_xxz[i] = t_yyyz_xz[i] * ab_x[i] + t_xyyyz_xz[i];

        t_yyyz_xyy[i] = t_yyyz_yy[i] * ab_x[i] + t_xyyyz_yy[i];

        t_yyyz_xyz[i] = t_yyyz_yz[i] * ab_x[i] + t_xyyyz_yz[i];

        t_yyyz_xzz[i] = t_yyyz_zz[i] * ab_x[i] + t_xyyyz_zz[i];

        t_yyyz_yyy[i] = t_yyyz_yy[i] * ab_y[i] + t_yyyyz_yy[i];

        t_yyyz_yyz[i] = t_yyyz_yz[i] * ab_y[i] + t_yyyyz_yz[i];

        t_yyyz_yzz[i] = t_yyyz_zz[i] * ab_y[i] + t_yyyyz_zz[i];

        t_yyyz_zzz[i] = t_yyyz_zz[i] * ab_z[i] + t_yyyzz_zz[i];

        t_yyzz_xxx[i] = t_yyzz_xx[i] * ab_x[i] + t_xyyzz_xx[i];

        t_yyzz_xxy[i] = t_yyzz_xy[i] * ab_x[i] + t_xyyzz_xy[i];

        t_yyzz_xxz[i] = t_yyzz_xz[i] * ab_x[i] + t_xyyzz_xz[i];

        t_yyzz_xyy[i] = t_yyzz_yy[i] * ab_x[i] + t_xyyzz_yy[i];

        t_yyzz_xyz[i] = t_yyzz_yz[i] * ab_x[i] + t_xyyzz_yz[i];

        t_yyzz_xzz[i] = t_yyzz_zz[i] * ab_x[i] + t_xyyzz_zz[i];

        t_yyzz_yyy[i] = t_yyzz_yy[i] * ab_y[i] + t_yyyzz_yy[i];

        t_yyzz_yyz[i] = t_yyzz_yz[i] * ab_y[i] + t_yyyzz_yz[i];

        t_yyzz_yzz[i] = t_yyzz_zz[i] * ab_y[i] + t_yyyzz_zz[i];

        t_yyzz_zzz[i] = t_yyzz_zz[i] * ab_z[i] + t_yyzzz_zz[i];

        t_yzzz_xxx[i] = t_yzzz_xx[i] * ab_x[i] + t_xyzzz_xx[i];

        t_yzzz_xxy[i] = t_yzzz_xy[i] * ab_x[i] + t_xyzzz_xy[i];

        t_yzzz_xxz[i] = t_yzzz_xz[i] * ab_x[i] + t_xyzzz_xz[i];

        t_yzzz_xyy[i] = t_yzzz_yy[i] * ab_x[i] + t_xyzzz_yy[i];

        t_yzzz_xyz[i] = t_yzzz_yz[i] * ab_x[i] + t_xyzzz_yz[i];

        t_yzzz_xzz[i] = t_yzzz_zz[i] * ab_x[i] + t_xyzzz_zz[i];

        t_yzzz_yyy[i] = t_yzzz_yy[i] * ab_y[i] + t_yyzzz_yy[i];

        t_yzzz_yyz[i] = t_yzzz_yz[i] * ab_y[i] + t_yyzzz_yz[i];

        t_yzzz_yzz[i] = t_yzzz_zz[i] * ab_y[i] + t_yyzzz_zz[i];

        t_yzzz_zzz[i] = t_yzzz_zz[i] * ab_z[i] + t_yzzzz_zz[i];

        t_zzzz_xxx[i] = t_zzzz_xx[i] * ab_x[i] + t_xzzzz_xx[i];

        t_zzzz_xxy[i] = t_zzzz_xy[i] * ab_x[i] + t_xzzzz_xy[i];

        t_zzzz_xxz[i] = t_zzzz_xz[i] * ab_x[i] + t_xzzzz_xz[i];

        t_zzzz_xyy[i] = t_zzzz_yy[i] * ab_x[i] + t_xzzzz_yy[i];

        t_zzzz_xyz[i] = t_zzzz_yz[i] * ab_x[i] + t_xzzzz_yz[i];

        t_zzzz_xzz[i] = t_zzzz_zz[i] * ab_x[i] + t_xzzzz_zz[i];

        t_zzzz_yyy[i] = t_zzzz_yy[i] * ab_y[i] + t_yzzzz_yy[i];

        t_zzzz_yyz[i] = t_zzzz_yz[i] * ab_y[i] + t_yzzzz_yz[i];

        t_zzzz_yzz[i] = t_zzzz_zz[i] * ab_y[i] + t_yzzzz_zz[i];

        t_zzzz_zzz[i] = t_zzzz_zz[i] * ab_z[i] + t_zzzzz_zz[i];
    }
}

} // t2chrr namespace

