#include "GeometricalDerivatives0X1ForGG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_gg(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_gg,
                       const int idx_op_gf,
                       const int idx_op_gh,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : GF

        auto to_xxxx_xxx = pbuffer.data(idx_op_gf + i * 150 + 0);

        auto to_xxxx_xxy = pbuffer.data(idx_op_gf + i * 150 + 1);

        auto to_xxxx_xxz = pbuffer.data(idx_op_gf + i * 150 + 2);

        auto to_xxxx_xyy = pbuffer.data(idx_op_gf + i * 150 + 3);

        auto to_xxxx_xyz = pbuffer.data(idx_op_gf + i * 150 + 4);

        auto to_xxxx_xzz = pbuffer.data(idx_op_gf + i * 150 + 5);

        auto to_xxxx_yyy = pbuffer.data(idx_op_gf + i * 150 + 6);

        auto to_xxxx_yyz = pbuffer.data(idx_op_gf + i * 150 + 7);

        auto to_xxxx_yzz = pbuffer.data(idx_op_gf + i * 150 + 8);

        auto to_xxxx_zzz = pbuffer.data(idx_op_gf + i * 150 + 9);

        auto to_xxxy_xxx = pbuffer.data(idx_op_gf + i * 150 + 10);

        auto to_xxxy_xxy = pbuffer.data(idx_op_gf + i * 150 + 11);

        auto to_xxxy_xxz = pbuffer.data(idx_op_gf + i * 150 + 12);

        auto to_xxxy_xyy = pbuffer.data(idx_op_gf + i * 150 + 13);

        auto to_xxxy_xyz = pbuffer.data(idx_op_gf + i * 150 + 14);

        auto to_xxxy_xzz = pbuffer.data(idx_op_gf + i * 150 + 15);

        auto to_xxxy_yyy = pbuffer.data(idx_op_gf + i * 150 + 16);

        auto to_xxxy_yyz = pbuffer.data(idx_op_gf + i * 150 + 17);

        auto to_xxxy_yzz = pbuffer.data(idx_op_gf + i * 150 + 18);

        auto to_xxxy_zzz = pbuffer.data(idx_op_gf + i * 150 + 19);

        auto to_xxxz_xxx = pbuffer.data(idx_op_gf + i * 150 + 20);

        auto to_xxxz_xxy = pbuffer.data(idx_op_gf + i * 150 + 21);

        auto to_xxxz_xxz = pbuffer.data(idx_op_gf + i * 150 + 22);

        auto to_xxxz_xyy = pbuffer.data(idx_op_gf + i * 150 + 23);

        auto to_xxxz_xyz = pbuffer.data(idx_op_gf + i * 150 + 24);

        auto to_xxxz_xzz = pbuffer.data(idx_op_gf + i * 150 + 25);

        auto to_xxxz_yyy = pbuffer.data(idx_op_gf + i * 150 + 26);

        auto to_xxxz_yyz = pbuffer.data(idx_op_gf + i * 150 + 27);

        auto to_xxxz_yzz = pbuffer.data(idx_op_gf + i * 150 + 28);

        auto to_xxxz_zzz = pbuffer.data(idx_op_gf + i * 150 + 29);

        auto to_xxyy_xxx = pbuffer.data(idx_op_gf + i * 150 + 30);

        auto to_xxyy_xxy = pbuffer.data(idx_op_gf + i * 150 + 31);

        auto to_xxyy_xxz = pbuffer.data(idx_op_gf + i * 150 + 32);

        auto to_xxyy_xyy = pbuffer.data(idx_op_gf + i * 150 + 33);

        auto to_xxyy_xyz = pbuffer.data(idx_op_gf + i * 150 + 34);

        auto to_xxyy_xzz = pbuffer.data(idx_op_gf + i * 150 + 35);

        auto to_xxyy_yyy = pbuffer.data(idx_op_gf + i * 150 + 36);

        auto to_xxyy_yyz = pbuffer.data(idx_op_gf + i * 150 + 37);

        auto to_xxyy_yzz = pbuffer.data(idx_op_gf + i * 150 + 38);

        auto to_xxyy_zzz = pbuffer.data(idx_op_gf + i * 150 + 39);

        auto to_xxyz_xxx = pbuffer.data(idx_op_gf + i * 150 + 40);

        auto to_xxyz_xxy = pbuffer.data(idx_op_gf + i * 150 + 41);

        auto to_xxyz_xxz = pbuffer.data(idx_op_gf + i * 150 + 42);

        auto to_xxyz_xyy = pbuffer.data(idx_op_gf + i * 150 + 43);

        auto to_xxyz_xyz = pbuffer.data(idx_op_gf + i * 150 + 44);

        auto to_xxyz_xzz = pbuffer.data(idx_op_gf + i * 150 + 45);

        auto to_xxyz_yyy = pbuffer.data(idx_op_gf + i * 150 + 46);

        auto to_xxyz_yyz = pbuffer.data(idx_op_gf + i * 150 + 47);

        auto to_xxyz_yzz = pbuffer.data(idx_op_gf + i * 150 + 48);

        auto to_xxyz_zzz = pbuffer.data(idx_op_gf + i * 150 + 49);

        auto to_xxzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 50);

        auto to_xxzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 51);

        auto to_xxzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 52);

        auto to_xxzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 53);

        auto to_xxzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 54);

        auto to_xxzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 55);

        auto to_xxzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 56);

        auto to_xxzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 57);

        auto to_xxzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 58);

        auto to_xxzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 59);

        auto to_xyyy_xxx = pbuffer.data(idx_op_gf + i * 150 + 60);

        auto to_xyyy_xxy = pbuffer.data(idx_op_gf + i * 150 + 61);

        auto to_xyyy_xxz = pbuffer.data(idx_op_gf + i * 150 + 62);

        auto to_xyyy_xyy = pbuffer.data(idx_op_gf + i * 150 + 63);

        auto to_xyyy_xyz = pbuffer.data(idx_op_gf + i * 150 + 64);

        auto to_xyyy_xzz = pbuffer.data(idx_op_gf + i * 150 + 65);

        auto to_xyyy_yyy = pbuffer.data(idx_op_gf + i * 150 + 66);

        auto to_xyyy_yyz = pbuffer.data(idx_op_gf + i * 150 + 67);

        auto to_xyyy_yzz = pbuffer.data(idx_op_gf + i * 150 + 68);

        auto to_xyyy_zzz = pbuffer.data(idx_op_gf + i * 150 + 69);

        auto to_xyyz_xxx = pbuffer.data(idx_op_gf + i * 150 + 70);

        auto to_xyyz_xxy = pbuffer.data(idx_op_gf + i * 150 + 71);

        auto to_xyyz_xxz = pbuffer.data(idx_op_gf + i * 150 + 72);

        auto to_xyyz_xyy = pbuffer.data(idx_op_gf + i * 150 + 73);

        auto to_xyyz_xyz = pbuffer.data(idx_op_gf + i * 150 + 74);

        auto to_xyyz_xzz = pbuffer.data(idx_op_gf + i * 150 + 75);

        auto to_xyyz_yyy = pbuffer.data(idx_op_gf + i * 150 + 76);

        auto to_xyyz_yyz = pbuffer.data(idx_op_gf + i * 150 + 77);

        auto to_xyyz_yzz = pbuffer.data(idx_op_gf + i * 150 + 78);

        auto to_xyyz_zzz = pbuffer.data(idx_op_gf + i * 150 + 79);

        auto to_xyzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 80);

        auto to_xyzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 81);

        auto to_xyzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 82);

        auto to_xyzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 83);

        auto to_xyzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 84);

        auto to_xyzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 85);

        auto to_xyzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 86);

        auto to_xyzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 87);

        auto to_xyzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 88);

        auto to_xyzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 89);

        auto to_xzzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 90);

        auto to_xzzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 91);

        auto to_xzzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 92);

        auto to_xzzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 93);

        auto to_xzzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 94);

        auto to_xzzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 95);

        auto to_xzzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 96);

        auto to_xzzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 97);

        auto to_xzzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 98);

        auto to_xzzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 99);

        auto to_yyyy_xxx = pbuffer.data(idx_op_gf + i * 150 + 100);

        auto to_yyyy_xxy = pbuffer.data(idx_op_gf + i * 150 + 101);

        auto to_yyyy_xxz = pbuffer.data(idx_op_gf + i * 150 + 102);

        auto to_yyyy_xyy = pbuffer.data(idx_op_gf + i * 150 + 103);

        auto to_yyyy_xyz = pbuffer.data(idx_op_gf + i * 150 + 104);

        auto to_yyyy_xzz = pbuffer.data(idx_op_gf + i * 150 + 105);

        auto to_yyyy_yyy = pbuffer.data(idx_op_gf + i * 150 + 106);

        auto to_yyyy_yyz = pbuffer.data(idx_op_gf + i * 150 + 107);

        auto to_yyyy_yzz = pbuffer.data(idx_op_gf + i * 150 + 108);

        auto to_yyyy_zzz = pbuffer.data(idx_op_gf + i * 150 + 109);

        auto to_yyyz_xxx = pbuffer.data(idx_op_gf + i * 150 + 110);

        auto to_yyyz_xxy = pbuffer.data(idx_op_gf + i * 150 + 111);

        auto to_yyyz_xxz = pbuffer.data(idx_op_gf + i * 150 + 112);

        auto to_yyyz_xyy = pbuffer.data(idx_op_gf + i * 150 + 113);

        auto to_yyyz_xyz = pbuffer.data(idx_op_gf + i * 150 + 114);

        auto to_yyyz_xzz = pbuffer.data(idx_op_gf + i * 150 + 115);

        auto to_yyyz_yyy = pbuffer.data(idx_op_gf + i * 150 + 116);

        auto to_yyyz_yyz = pbuffer.data(idx_op_gf + i * 150 + 117);

        auto to_yyyz_yzz = pbuffer.data(idx_op_gf + i * 150 + 118);

        auto to_yyyz_zzz = pbuffer.data(idx_op_gf + i * 150 + 119);

        auto to_yyzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 120);

        auto to_yyzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 121);

        auto to_yyzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 122);

        auto to_yyzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 123);

        auto to_yyzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 124);

        auto to_yyzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 125);

        auto to_yyzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 126);

        auto to_yyzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 127);

        auto to_yyzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 128);

        auto to_yyzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 129);

        auto to_yzzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 130);

        auto to_yzzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 131);

        auto to_yzzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 132);

        auto to_yzzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 133);

        auto to_yzzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 134);

        auto to_yzzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 135);

        auto to_yzzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 136);

        auto to_yzzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 137);

        auto to_yzzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 138);

        auto to_yzzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 139);

        auto to_zzzz_xxx = pbuffer.data(idx_op_gf + i * 150 + 140);

        auto to_zzzz_xxy = pbuffer.data(idx_op_gf + i * 150 + 141);

        auto to_zzzz_xxz = pbuffer.data(idx_op_gf + i * 150 + 142);

        auto to_zzzz_xyy = pbuffer.data(idx_op_gf + i * 150 + 143);

        auto to_zzzz_xyz = pbuffer.data(idx_op_gf + i * 150 + 144);

        auto to_zzzz_xzz = pbuffer.data(idx_op_gf + i * 150 + 145);

        auto to_zzzz_yyy = pbuffer.data(idx_op_gf + i * 150 + 146);

        auto to_zzzz_yyz = pbuffer.data(idx_op_gf + i * 150 + 147);

        auto to_zzzz_yzz = pbuffer.data(idx_op_gf + i * 150 + 148);

        auto to_zzzz_zzz = pbuffer.data(idx_op_gf + i * 150 + 149);

        // Set up components of auxiliary buffer : GH

        auto to_xxxx_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 0);

        auto to_xxxx_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 1);

        auto to_xxxx_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 2);

        auto to_xxxx_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 3);

        auto to_xxxx_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 4);

        auto to_xxxx_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 5);

        auto to_xxxx_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 6);

        auto to_xxxx_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 7);

        auto to_xxxx_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 8);

        auto to_xxxx_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 9);

        auto to_xxxx_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 10);

        auto to_xxxx_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 11);

        auto to_xxxx_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 12);

        auto to_xxxx_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 13);

        auto to_xxxx_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 14);

        auto to_xxxx_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 15);

        auto to_xxxx_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 16);

        auto to_xxxx_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 17);

        auto to_xxxx_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 18);

        auto to_xxxx_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 19);

        auto to_xxxx_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 20);

        auto to_xxxy_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 21);

        auto to_xxxy_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 22);

        auto to_xxxy_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 23);

        auto to_xxxy_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 24);

        auto to_xxxy_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 25);

        auto to_xxxy_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 26);

        auto to_xxxy_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 27);

        auto to_xxxy_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 28);

        auto to_xxxy_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 29);

        auto to_xxxy_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 30);

        auto to_xxxy_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 31);

        auto to_xxxy_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 32);

        auto to_xxxy_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 33);

        auto to_xxxy_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 34);

        auto to_xxxy_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 35);

        auto to_xxxy_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 36);

        auto to_xxxy_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 37);

        auto to_xxxy_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 38);

        auto to_xxxy_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 39);

        auto to_xxxy_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 40);

        auto to_xxxy_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 41);

        auto to_xxxz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 42);

        auto to_xxxz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 43);

        auto to_xxxz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 44);

        auto to_xxxz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 45);

        auto to_xxxz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 46);

        auto to_xxxz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 47);

        auto to_xxxz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 48);

        auto to_xxxz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 49);

        auto to_xxxz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 50);

        auto to_xxxz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 51);

        auto to_xxxz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 52);

        auto to_xxxz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 53);

        auto to_xxxz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 54);

        auto to_xxxz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 55);

        auto to_xxxz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 56);

        auto to_xxxz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 57);

        auto to_xxxz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 58);

        auto to_xxxz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 59);

        auto to_xxxz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 60);

        auto to_xxxz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 61);

        auto to_xxxz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 62);

        auto to_xxyy_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 63);

        auto to_xxyy_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 64);

        auto to_xxyy_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 65);

        auto to_xxyy_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 66);

        auto to_xxyy_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 67);

        auto to_xxyy_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 68);

        auto to_xxyy_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 69);

        auto to_xxyy_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 70);

        auto to_xxyy_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 71);

        auto to_xxyy_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 72);

        auto to_xxyy_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 73);

        auto to_xxyy_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 74);

        auto to_xxyy_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 75);

        auto to_xxyy_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 76);

        auto to_xxyy_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 77);

        auto to_xxyy_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 78);

        auto to_xxyy_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 79);

        auto to_xxyy_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 80);

        auto to_xxyy_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 81);

        auto to_xxyy_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 82);

        auto to_xxyy_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 83);

        auto to_xxyz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 84);

        auto to_xxyz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 85);

        auto to_xxyz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 86);

        auto to_xxyz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 87);

        auto to_xxyz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 88);

        auto to_xxyz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 89);

        auto to_xxyz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 90);

        auto to_xxyz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 91);

        auto to_xxyz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 92);

        auto to_xxyz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 93);

        auto to_xxyz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 94);

        auto to_xxyz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 95);

        auto to_xxyz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 96);

        auto to_xxyz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 97);

        auto to_xxyz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 98);

        auto to_xxyz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 99);

        auto to_xxyz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 100);

        auto to_xxyz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 101);

        auto to_xxyz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 102);

        auto to_xxyz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 103);

        auto to_xxyz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 104);

        auto to_xxzz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 105);

        auto to_xxzz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 106);

        auto to_xxzz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 107);

        auto to_xxzz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 108);

        auto to_xxzz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 109);

        auto to_xxzz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 110);

        auto to_xxzz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 111);

        auto to_xxzz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 112);

        auto to_xxzz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 113);

        auto to_xxzz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 114);

        auto to_xxzz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 115);

        auto to_xxzz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 116);

        auto to_xxzz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 117);

        auto to_xxzz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 118);

        auto to_xxzz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 119);

        auto to_xxzz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 120);

        auto to_xxzz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 121);

        auto to_xxzz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 122);

        auto to_xxzz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 123);

        auto to_xxzz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 124);

        auto to_xxzz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 125);

        auto to_xyyy_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 126);

        auto to_xyyy_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 127);

        auto to_xyyy_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 128);

        auto to_xyyy_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 129);

        auto to_xyyy_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 130);

        auto to_xyyy_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 131);

        auto to_xyyy_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 132);

        auto to_xyyy_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 133);

        auto to_xyyy_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 134);

        auto to_xyyy_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 135);

        auto to_xyyy_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 136);

        auto to_xyyy_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 137);

        auto to_xyyy_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 138);

        auto to_xyyy_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 139);

        auto to_xyyy_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 140);

        auto to_xyyy_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 141);

        auto to_xyyy_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 142);

        auto to_xyyy_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 143);

        auto to_xyyy_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 144);

        auto to_xyyy_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 145);

        auto to_xyyy_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 146);

        auto to_xyyz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 147);

        auto to_xyyz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 148);

        auto to_xyyz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 149);

        auto to_xyyz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 150);

        auto to_xyyz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 151);

        auto to_xyyz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 152);

        auto to_xyyz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 153);

        auto to_xyyz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 154);

        auto to_xyyz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 155);

        auto to_xyyz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 156);

        auto to_xyyz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 157);

        auto to_xyyz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 158);

        auto to_xyyz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 159);

        auto to_xyyz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 160);

        auto to_xyyz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 161);

        auto to_xyyz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 162);

        auto to_xyyz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 163);

        auto to_xyyz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 164);

        auto to_xyyz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 165);

        auto to_xyyz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 166);

        auto to_xyyz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 167);

        auto to_xyzz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 168);

        auto to_xyzz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 169);

        auto to_xyzz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 170);

        auto to_xyzz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 171);

        auto to_xyzz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 172);

        auto to_xyzz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 173);

        auto to_xyzz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 174);

        auto to_xyzz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 175);

        auto to_xyzz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 176);

        auto to_xyzz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 177);

        auto to_xyzz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 178);

        auto to_xyzz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 179);

        auto to_xyzz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 180);

        auto to_xyzz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 181);

        auto to_xyzz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 182);

        auto to_xyzz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 183);

        auto to_xyzz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 184);

        auto to_xyzz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 185);

        auto to_xyzz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 186);

        auto to_xyzz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 187);

        auto to_xyzz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 188);

        auto to_xzzz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 189);

        auto to_xzzz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 190);

        auto to_xzzz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 191);

        auto to_xzzz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 192);

        auto to_xzzz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 193);

        auto to_xzzz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 194);

        auto to_xzzz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 195);

        auto to_xzzz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 196);

        auto to_xzzz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 197);

        auto to_xzzz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 198);

        auto to_xzzz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 199);

        auto to_xzzz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 200);

        auto to_xzzz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 201);

        auto to_xzzz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 202);

        auto to_xzzz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 203);

        auto to_xzzz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 204);

        auto to_xzzz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 205);

        auto to_xzzz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 206);

        auto to_xzzz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 207);

        auto to_xzzz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 208);

        auto to_xzzz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 209);

        auto to_yyyy_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 210);

        auto to_yyyy_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 211);

        auto to_yyyy_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 212);

        auto to_yyyy_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 213);

        auto to_yyyy_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 214);

        auto to_yyyy_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 215);

        auto to_yyyy_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 216);

        auto to_yyyy_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 217);

        auto to_yyyy_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 218);

        auto to_yyyy_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 219);

        auto to_yyyy_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 220);

        auto to_yyyy_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 221);

        auto to_yyyy_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 222);

        auto to_yyyy_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 223);

        auto to_yyyy_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 224);

        auto to_yyyy_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 225);

        auto to_yyyy_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 226);

        auto to_yyyy_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 227);

        auto to_yyyy_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 228);

        auto to_yyyy_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 229);

        auto to_yyyy_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 230);

        auto to_yyyz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 231);

        auto to_yyyz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 232);

        auto to_yyyz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 233);

        auto to_yyyz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 234);

        auto to_yyyz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 235);

        auto to_yyyz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 236);

        auto to_yyyz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 237);

        auto to_yyyz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 238);

        auto to_yyyz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 239);

        auto to_yyyz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 240);

        auto to_yyyz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 241);

        auto to_yyyz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 242);

        auto to_yyyz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 243);

        auto to_yyyz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 244);

        auto to_yyyz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 245);

        auto to_yyyz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 246);

        auto to_yyyz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 247);

        auto to_yyyz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 248);

        auto to_yyyz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 249);

        auto to_yyyz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 250);

        auto to_yyyz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 251);

        auto to_yyzz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 252);

        auto to_yyzz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 253);

        auto to_yyzz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 254);

        auto to_yyzz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 255);

        auto to_yyzz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 256);

        auto to_yyzz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 257);

        auto to_yyzz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 258);

        auto to_yyzz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 259);

        auto to_yyzz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 260);

        auto to_yyzz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 261);

        auto to_yyzz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 262);

        auto to_yyzz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 263);

        auto to_yyzz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 264);

        auto to_yyzz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 265);

        auto to_yyzz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 266);

        auto to_yyzz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 267);

        auto to_yyzz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 268);

        auto to_yyzz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 269);

        auto to_yyzz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 270);

        auto to_yyzz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 271);

        auto to_yyzz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 272);

        auto to_yzzz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 273);

        auto to_yzzz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 274);

        auto to_yzzz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 275);

        auto to_yzzz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 276);

        auto to_yzzz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 277);

        auto to_yzzz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 278);

        auto to_yzzz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 279);

        auto to_yzzz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 280);

        auto to_yzzz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 281);

        auto to_yzzz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 282);

        auto to_yzzz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 283);

        auto to_yzzz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 284);

        auto to_yzzz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 285);

        auto to_yzzz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 286);

        auto to_yzzz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 287);

        auto to_yzzz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 288);

        auto to_yzzz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 289);

        auto to_yzzz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 290);

        auto to_yzzz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 291);

        auto to_yzzz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 292);

        auto to_yzzz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 293);

        auto to_zzzz_xxxxx = pbuffer.data(idx_op_gh + i * 315 + 294);

        auto to_zzzz_xxxxy = pbuffer.data(idx_op_gh + i * 315 + 295);

        auto to_zzzz_xxxxz = pbuffer.data(idx_op_gh + i * 315 + 296);

        auto to_zzzz_xxxyy = pbuffer.data(idx_op_gh + i * 315 + 297);

        auto to_zzzz_xxxyz = pbuffer.data(idx_op_gh + i * 315 + 298);

        auto to_zzzz_xxxzz = pbuffer.data(idx_op_gh + i * 315 + 299);

        auto to_zzzz_xxyyy = pbuffer.data(idx_op_gh + i * 315 + 300);

        auto to_zzzz_xxyyz = pbuffer.data(idx_op_gh + i * 315 + 301);

        auto to_zzzz_xxyzz = pbuffer.data(idx_op_gh + i * 315 + 302);

        auto to_zzzz_xxzzz = pbuffer.data(idx_op_gh + i * 315 + 303);

        auto to_zzzz_xyyyy = pbuffer.data(idx_op_gh + i * 315 + 304);

        auto to_zzzz_xyyyz = pbuffer.data(idx_op_gh + i * 315 + 305);

        auto to_zzzz_xyyzz = pbuffer.data(idx_op_gh + i * 315 + 306);

        auto to_zzzz_xyzzz = pbuffer.data(idx_op_gh + i * 315 + 307);

        auto to_zzzz_xzzzz = pbuffer.data(idx_op_gh + i * 315 + 308);

        auto to_zzzz_yyyyy = pbuffer.data(idx_op_gh + i * 315 + 309);

        auto to_zzzz_yyyyz = pbuffer.data(idx_op_gh + i * 315 + 310);

        auto to_zzzz_yyyzz = pbuffer.data(idx_op_gh + i * 315 + 311);

        auto to_zzzz_yyzzz = pbuffer.data(idx_op_gh + i * 315 + 312);

        auto to_zzzz_yzzzz = pbuffer.data(idx_op_gh + i * 315 + 313);

        auto to_zzzz_zzzzz = pbuffer.data(idx_op_gh + i * 315 + 314);

        // Set up 0-15 components of targeted buffer : GG

        auto to_0_x_xxxx_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 0);

        auto to_0_x_xxxx_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 1);

        auto to_0_x_xxxx_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 2);

        auto to_0_x_xxxx_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 3);

        auto to_0_x_xxxx_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 4);

        auto to_0_x_xxxx_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 5);

        auto to_0_x_xxxx_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 6);

        auto to_0_x_xxxx_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 7);

        auto to_0_x_xxxx_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 8);

        auto to_0_x_xxxx_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 9);

        auto to_0_x_xxxx_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 10);

        auto to_0_x_xxxx_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 11);

        auto to_0_x_xxxx_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 12);

        auto to_0_x_xxxx_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 13);

        auto to_0_x_xxxx_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_0_x_xxxx_xxxx, to_0_x_xxxx_xxxy, to_0_x_xxxx_xxxz, to_0_x_xxxx_xxyy, to_0_x_xxxx_xxyz, to_0_x_xxxx_xxzz, to_0_x_xxxx_xyyy, to_0_x_xxxx_xyyz, to_0_x_xxxx_xyzz, to_0_x_xxxx_xzzz, to_0_x_xxxx_yyyy, to_0_x_xxxx_yyyz, to_0_x_xxxx_yyzz, to_0_x_xxxx_yzzz, to_0_x_xxxx_zzzz, to_xxxx_xxx, to_xxxx_xxxxx, to_xxxx_xxxxy, to_xxxx_xxxxz, to_xxxx_xxxyy, to_xxxx_xxxyz, to_xxxx_xxxzz, to_xxxx_xxy, to_xxxx_xxyyy, to_xxxx_xxyyz, to_xxxx_xxyzz, to_xxxx_xxz, to_xxxx_xxzzz, to_xxxx_xyy, to_xxxx_xyyyy, to_xxxx_xyyyz, to_xxxx_xyyzz, to_xxxx_xyz, to_xxxx_xyzzz, to_xxxx_xzz, to_xxxx_xzzzz, to_xxxx_yyy, to_xxxx_yyz, to_xxxx_yzz, to_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxx_xxxx[k] = -4.0 * to_xxxx_xxx[k] + 2.0 * to_xxxx_xxxxx[k] * tke_0;

            to_0_x_xxxx_xxxy[k] = -3.0 * to_xxxx_xxy[k] + 2.0 * to_xxxx_xxxxy[k] * tke_0;

            to_0_x_xxxx_xxxz[k] = -3.0 * to_xxxx_xxz[k] + 2.0 * to_xxxx_xxxxz[k] * tke_0;

            to_0_x_xxxx_xxyy[k] = -2.0 * to_xxxx_xyy[k] + 2.0 * to_xxxx_xxxyy[k] * tke_0;

            to_0_x_xxxx_xxyz[k] = -2.0 * to_xxxx_xyz[k] + 2.0 * to_xxxx_xxxyz[k] * tke_0;

            to_0_x_xxxx_xxzz[k] = -2.0 * to_xxxx_xzz[k] + 2.0 * to_xxxx_xxxzz[k] * tke_0;

            to_0_x_xxxx_xyyy[k] = -to_xxxx_yyy[k] + 2.0 * to_xxxx_xxyyy[k] * tke_0;

            to_0_x_xxxx_xyyz[k] = -to_xxxx_yyz[k] + 2.0 * to_xxxx_xxyyz[k] * tke_0;

            to_0_x_xxxx_xyzz[k] = -to_xxxx_yzz[k] + 2.0 * to_xxxx_xxyzz[k] * tke_0;

            to_0_x_xxxx_xzzz[k] = -to_xxxx_zzz[k] + 2.0 * to_xxxx_xxzzz[k] * tke_0;

            to_0_x_xxxx_yyyy[k] = 2.0 * to_xxxx_xyyyy[k] * tke_0;

            to_0_x_xxxx_yyyz[k] = 2.0 * to_xxxx_xyyyz[k] * tke_0;

            to_0_x_xxxx_yyzz[k] = 2.0 * to_xxxx_xyyzz[k] * tke_0;

            to_0_x_xxxx_yzzz[k] = 2.0 * to_xxxx_xyzzz[k] * tke_0;

            to_0_x_xxxx_zzzz[k] = 2.0 * to_xxxx_xzzzz[k] * tke_0;
        }

        // Set up 15-30 components of targeted buffer : GG

        auto to_0_x_xxxy_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 15);

        auto to_0_x_xxxy_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 16);

        auto to_0_x_xxxy_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 17);

        auto to_0_x_xxxy_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 18);

        auto to_0_x_xxxy_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 19);

        auto to_0_x_xxxy_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 20);

        auto to_0_x_xxxy_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 21);

        auto to_0_x_xxxy_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 22);

        auto to_0_x_xxxy_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 23);

        auto to_0_x_xxxy_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 24);

        auto to_0_x_xxxy_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 25);

        auto to_0_x_xxxy_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 26);

        auto to_0_x_xxxy_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 27);

        auto to_0_x_xxxy_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 28);

        auto to_0_x_xxxy_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_0_x_xxxy_xxxx, to_0_x_xxxy_xxxy, to_0_x_xxxy_xxxz, to_0_x_xxxy_xxyy, to_0_x_xxxy_xxyz, to_0_x_xxxy_xxzz, to_0_x_xxxy_xyyy, to_0_x_xxxy_xyyz, to_0_x_xxxy_xyzz, to_0_x_xxxy_xzzz, to_0_x_xxxy_yyyy, to_0_x_xxxy_yyyz, to_0_x_xxxy_yyzz, to_0_x_xxxy_yzzz, to_0_x_xxxy_zzzz, to_xxxy_xxx, to_xxxy_xxxxx, to_xxxy_xxxxy, to_xxxy_xxxxz, to_xxxy_xxxyy, to_xxxy_xxxyz, to_xxxy_xxxzz, to_xxxy_xxy, to_xxxy_xxyyy, to_xxxy_xxyyz, to_xxxy_xxyzz, to_xxxy_xxz, to_xxxy_xxzzz, to_xxxy_xyy, to_xxxy_xyyyy, to_xxxy_xyyyz, to_xxxy_xyyzz, to_xxxy_xyz, to_xxxy_xyzzz, to_xxxy_xzz, to_xxxy_xzzzz, to_xxxy_yyy, to_xxxy_yyz, to_xxxy_yzz, to_xxxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxy_xxxx[k] = -4.0 * to_xxxy_xxx[k] + 2.0 * to_xxxy_xxxxx[k] * tke_0;

            to_0_x_xxxy_xxxy[k] = -3.0 * to_xxxy_xxy[k] + 2.0 * to_xxxy_xxxxy[k] * tke_0;

            to_0_x_xxxy_xxxz[k] = -3.0 * to_xxxy_xxz[k] + 2.0 * to_xxxy_xxxxz[k] * tke_0;

            to_0_x_xxxy_xxyy[k] = -2.0 * to_xxxy_xyy[k] + 2.0 * to_xxxy_xxxyy[k] * tke_0;

            to_0_x_xxxy_xxyz[k] = -2.0 * to_xxxy_xyz[k] + 2.0 * to_xxxy_xxxyz[k] * tke_0;

            to_0_x_xxxy_xxzz[k] = -2.0 * to_xxxy_xzz[k] + 2.0 * to_xxxy_xxxzz[k] * tke_0;

            to_0_x_xxxy_xyyy[k] = -to_xxxy_yyy[k] + 2.0 * to_xxxy_xxyyy[k] * tke_0;

            to_0_x_xxxy_xyyz[k] = -to_xxxy_yyz[k] + 2.0 * to_xxxy_xxyyz[k] * tke_0;

            to_0_x_xxxy_xyzz[k] = -to_xxxy_yzz[k] + 2.0 * to_xxxy_xxyzz[k] * tke_0;

            to_0_x_xxxy_xzzz[k] = -to_xxxy_zzz[k] + 2.0 * to_xxxy_xxzzz[k] * tke_0;

            to_0_x_xxxy_yyyy[k] = 2.0 * to_xxxy_xyyyy[k] * tke_0;

            to_0_x_xxxy_yyyz[k] = 2.0 * to_xxxy_xyyyz[k] * tke_0;

            to_0_x_xxxy_yyzz[k] = 2.0 * to_xxxy_xyyzz[k] * tke_0;

            to_0_x_xxxy_yzzz[k] = 2.0 * to_xxxy_xyzzz[k] * tke_0;

            to_0_x_xxxy_zzzz[k] = 2.0 * to_xxxy_xzzzz[k] * tke_0;
        }

        // Set up 30-45 components of targeted buffer : GG

        auto to_0_x_xxxz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 30);

        auto to_0_x_xxxz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 31);

        auto to_0_x_xxxz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 32);

        auto to_0_x_xxxz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 33);

        auto to_0_x_xxxz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 34);

        auto to_0_x_xxxz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 35);

        auto to_0_x_xxxz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 36);

        auto to_0_x_xxxz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 37);

        auto to_0_x_xxxz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 38);

        auto to_0_x_xxxz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 39);

        auto to_0_x_xxxz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 40);

        auto to_0_x_xxxz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 41);

        auto to_0_x_xxxz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 42);

        auto to_0_x_xxxz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 43);

        auto to_0_x_xxxz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_0_x_xxxz_xxxx, to_0_x_xxxz_xxxy, to_0_x_xxxz_xxxz, to_0_x_xxxz_xxyy, to_0_x_xxxz_xxyz, to_0_x_xxxz_xxzz, to_0_x_xxxz_xyyy, to_0_x_xxxz_xyyz, to_0_x_xxxz_xyzz, to_0_x_xxxz_xzzz, to_0_x_xxxz_yyyy, to_0_x_xxxz_yyyz, to_0_x_xxxz_yyzz, to_0_x_xxxz_yzzz, to_0_x_xxxz_zzzz, to_xxxz_xxx, to_xxxz_xxxxx, to_xxxz_xxxxy, to_xxxz_xxxxz, to_xxxz_xxxyy, to_xxxz_xxxyz, to_xxxz_xxxzz, to_xxxz_xxy, to_xxxz_xxyyy, to_xxxz_xxyyz, to_xxxz_xxyzz, to_xxxz_xxz, to_xxxz_xxzzz, to_xxxz_xyy, to_xxxz_xyyyy, to_xxxz_xyyyz, to_xxxz_xyyzz, to_xxxz_xyz, to_xxxz_xyzzz, to_xxxz_xzz, to_xxxz_xzzzz, to_xxxz_yyy, to_xxxz_yyz, to_xxxz_yzz, to_xxxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxz_xxxx[k] = -4.0 * to_xxxz_xxx[k] + 2.0 * to_xxxz_xxxxx[k] * tke_0;

            to_0_x_xxxz_xxxy[k] = -3.0 * to_xxxz_xxy[k] + 2.0 * to_xxxz_xxxxy[k] * tke_0;

            to_0_x_xxxz_xxxz[k] = -3.0 * to_xxxz_xxz[k] + 2.0 * to_xxxz_xxxxz[k] * tke_0;

            to_0_x_xxxz_xxyy[k] = -2.0 * to_xxxz_xyy[k] + 2.0 * to_xxxz_xxxyy[k] * tke_0;

            to_0_x_xxxz_xxyz[k] = -2.0 * to_xxxz_xyz[k] + 2.0 * to_xxxz_xxxyz[k] * tke_0;

            to_0_x_xxxz_xxzz[k] = -2.0 * to_xxxz_xzz[k] + 2.0 * to_xxxz_xxxzz[k] * tke_0;

            to_0_x_xxxz_xyyy[k] = -to_xxxz_yyy[k] + 2.0 * to_xxxz_xxyyy[k] * tke_0;

            to_0_x_xxxz_xyyz[k] = -to_xxxz_yyz[k] + 2.0 * to_xxxz_xxyyz[k] * tke_0;

            to_0_x_xxxz_xyzz[k] = -to_xxxz_yzz[k] + 2.0 * to_xxxz_xxyzz[k] * tke_0;

            to_0_x_xxxz_xzzz[k] = -to_xxxz_zzz[k] + 2.0 * to_xxxz_xxzzz[k] * tke_0;

            to_0_x_xxxz_yyyy[k] = 2.0 * to_xxxz_xyyyy[k] * tke_0;

            to_0_x_xxxz_yyyz[k] = 2.0 * to_xxxz_xyyyz[k] * tke_0;

            to_0_x_xxxz_yyzz[k] = 2.0 * to_xxxz_xyyzz[k] * tke_0;

            to_0_x_xxxz_yzzz[k] = 2.0 * to_xxxz_xyzzz[k] * tke_0;

            to_0_x_xxxz_zzzz[k] = 2.0 * to_xxxz_xzzzz[k] * tke_0;
        }

        // Set up 45-60 components of targeted buffer : GG

        auto to_0_x_xxyy_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 45);

        auto to_0_x_xxyy_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 46);

        auto to_0_x_xxyy_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 47);

        auto to_0_x_xxyy_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 48);

        auto to_0_x_xxyy_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 49);

        auto to_0_x_xxyy_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 50);

        auto to_0_x_xxyy_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 51);

        auto to_0_x_xxyy_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 52);

        auto to_0_x_xxyy_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 53);

        auto to_0_x_xxyy_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 54);

        auto to_0_x_xxyy_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 55);

        auto to_0_x_xxyy_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 56);

        auto to_0_x_xxyy_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 57);

        auto to_0_x_xxyy_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 58);

        auto to_0_x_xxyy_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_0_x_xxyy_xxxx, to_0_x_xxyy_xxxy, to_0_x_xxyy_xxxz, to_0_x_xxyy_xxyy, to_0_x_xxyy_xxyz, to_0_x_xxyy_xxzz, to_0_x_xxyy_xyyy, to_0_x_xxyy_xyyz, to_0_x_xxyy_xyzz, to_0_x_xxyy_xzzz, to_0_x_xxyy_yyyy, to_0_x_xxyy_yyyz, to_0_x_xxyy_yyzz, to_0_x_xxyy_yzzz, to_0_x_xxyy_zzzz, to_xxyy_xxx, to_xxyy_xxxxx, to_xxyy_xxxxy, to_xxyy_xxxxz, to_xxyy_xxxyy, to_xxyy_xxxyz, to_xxyy_xxxzz, to_xxyy_xxy, to_xxyy_xxyyy, to_xxyy_xxyyz, to_xxyy_xxyzz, to_xxyy_xxz, to_xxyy_xxzzz, to_xxyy_xyy, to_xxyy_xyyyy, to_xxyy_xyyyz, to_xxyy_xyyzz, to_xxyy_xyz, to_xxyy_xyzzz, to_xxyy_xzz, to_xxyy_xzzzz, to_xxyy_yyy, to_xxyy_yyz, to_xxyy_yzz, to_xxyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxyy_xxxx[k] = -4.0 * to_xxyy_xxx[k] + 2.0 * to_xxyy_xxxxx[k] * tke_0;

            to_0_x_xxyy_xxxy[k] = -3.0 * to_xxyy_xxy[k] + 2.0 * to_xxyy_xxxxy[k] * tke_0;

            to_0_x_xxyy_xxxz[k] = -3.0 * to_xxyy_xxz[k] + 2.0 * to_xxyy_xxxxz[k] * tke_0;

            to_0_x_xxyy_xxyy[k] = -2.0 * to_xxyy_xyy[k] + 2.0 * to_xxyy_xxxyy[k] * tke_0;

            to_0_x_xxyy_xxyz[k] = -2.0 * to_xxyy_xyz[k] + 2.0 * to_xxyy_xxxyz[k] * tke_0;

            to_0_x_xxyy_xxzz[k] = -2.0 * to_xxyy_xzz[k] + 2.0 * to_xxyy_xxxzz[k] * tke_0;

            to_0_x_xxyy_xyyy[k] = -to_xxyy_yyy[k] + 2.0 * to_xxyy_xxyyy[k] * tke_0;

            to_0_x_xxyy_xyyz[k] = -to_xxyy_yyz[k] + 2.0 * to_xxyy_xxyyz[k] * tke_0;

            to_0_x_xxyy_xyzz[k] = -to_xxyy_yzz[k] + 2.0 * to_xxyy_xxyzz[k] * tke_0;

            to_0_x_xxyy_xzzz[k] = -to_xxyy_zzz[k] + 2.0 * to_xxyy_xxzzz[k] * tke_0;

            to_0_x_xxyy_yyyy[k] = 2.0 * to_xxyy_xyyyy[k] * tke_0;

            to_0_x_xxyy_yyyz[k] = 2.0 * to_xxyy_xyyyz[k] * tke_0;

            to_0_x_xxyy_yyzz[k] = 2.0 * to_xxyy_xyyzz[k] * tke_0;

            to_0_x_xxyy_yzzz[k] = 2.0 * to_xxyy_xyzzz[k] * tke_0;

            to_0_x_xxyy_zzzz[k] = 2.0 * to_xxyy_xzzzz[k] * tke_0;
        }

        // Set up 60-75 components of targeted buffer : GG

        auto to_0_x_xxyz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 60);

        auto to_0_x_xxyz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 61);

        auto to_0_x_xxyz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 62);

        auto to_0_x_xxyz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 63);

        auto to_0_x_xxyz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 64);

        auto to_0_x_xxyz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 65);

        auto to_0_x_xxyz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 66);

        auto to_0_x_xxyz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 67);

        auto to_0_x_xxyz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 68);

        auto to_0_x_xxyz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 69);

        auto to_0_x_xxyz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 70);

        auto to_0_x_xxyz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 71);

        auto to_0_x_xxyz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 72);

        auto to_0_x_xxyz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 73);

        auto to_0_x_xxyz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_0_x_xxyz_xxxx, to_0_x_xxyz_xxxy, to_0_x_xxyz_xxxz, to_0_x_xxyz_xxyy, to_0_x_xxyz_xxyz, to_0_x_xxyz_xxzz, to_0_x_xxyz_xyyy, to_0_x_xxyz_xyyz, to_0_x_xxyz_xyzz, to_0_x_xxyz_xzzz, to_0_x_xxyz_yyyy, to_0_x_xxyz_yyyz, to_0_x_xxyz_yyzz, to_0_x_xxyz_yzzz, to_0_x_xxyz_zzzz, to_xxyz_xxx, to_xxyz_xxxxx, to_xxyz_xxxxy, to_xxyz_xxxxz, to_xxyz_xxxyy, to_xxyz_xxxyz, to_xxyz_xxxzz, to_xxyz_xxy, to_xxyz_xxyyy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xxzzz, to_xxyz_xyy, to_xxyz_xyyyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_xzzzz, to_xxyz_yyy, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxyz_xxxx[k] = -4.0 * to_xxyz_xxx[k] + 2.0 * to_xxyz_xxxxx[k] * tke_0;

            to_0_x_xxyz_xxxy[k] = -3.0 * to_xxyz_xxy[k] + 2.0 * to_xxyz_xxxxy[k] * tke_0;

            to_0_x_xxyz_xxxz[k] = -3.0 * to_xxyz_xxz[k] + 2.0 * to_xxyz_xxxxz[k] * tke_0;

            to_0_x_xxyz_xxyy[k] = -2.0 * to_xxyz_xyy[k] + 2.0 * to_xxyz_xxxyy[k] * tke_0;

            to_0_x_xxyz_xxyz[k] = -2.0 * to_xxyz_xyz[k] + 2.0 * to_xxyz_xxxyz[k] * tke_0;

            to_0_x_xxyz_xxzz[k] = -2.0 * to_xxyz_xzz[k] + 2.0 * to_xxyz_xxxzz[k] * tke_0;

            to_0_x_xxyz_xyyy[k] = -to_xxyz_yyy[k] + 2.0 * to_xxyz_xxyyy[k] * tke_0;

            to_0_x_xxyz_xyyz[k] = -to_xxyz_yyz[k] + 2.0 * to_xxyz_xxyyz[k] * tke_0;

            to_0_x_xxyz_xyzz[k] = -to_xxyz_yzz[k] + 2.0 * to_xxyz_xxyzz[k] * tke_0;

            to_0_x_xxyz_xzzz[k] = -to_xxyz_zzz[k] + 2.0 * to_xxyz_xxzzz[k] * tke_0;

            to_0_x_xxyz_yyyy[k] = 2.0 * to_xxyz_xyyyy[k] * tke_0;

            to_0_x_xxyz_yyyz[k] = 2.0 * to_xxyz_xyyyz[k] * tke_0;

            to_0_x_xxyz_yyzz[k] = 2.0 * to_xxyz_xyyzz[k] * tke_0;

            to_0_x_xxyz_yzzz[k] = 2.0 * to_xxyz_xyzzz[k] * tke_0;

            to_0_x_xxyz_zzzz[k] = 2.0 * to_xxyz_xzzzz[k] * tke_0;
        }

        // Set up 75-90 components of targeted buffer : GG

        auto to_0_x_xxzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 75);

        auto to_0_x_xxzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 76);

        auto to_0_x_xxzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 77);

        auto to_0_x_xxzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 78);

        auto to_0_x_xxzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 79);

        auto to_0_x_xxzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 80);

        auto to_0_x_xxzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 81);

        auto to_0_x_xxzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 82);

        auto to_0_x_xxzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 83);

        auto to_0_x_xxzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 84);

        auto to_0_x_xxzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 85);

        auto to_0_x_xxzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 86);

        auto to_0_x_xxzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 87);

        auto to_0_x_xxzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 88);

        auto to_0_x_xxzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_0_x_xxzz_xxxx, to_0_x_xxzz_xxxy, to_0_x_xxzz_xxxz, to_0_x_xxzz_xxyy, to_0_x_xxzz_xxyz, to_0_x_xxzz_xxzz, to_0_x_xxzz_xyyy, to_0_x_xxzz_xyyz, to_0_x_xxzz_xyzz, to_0_x_xxzz_xzzz, to_0_x_xxzz_yyyy, to_0_x_xxzz_yyyz, to_0_x_xxzz_yyzz, to_0_x_xxzz_yzzz, to_0_x_xxzz_zzzz, to_xxzz_xxx, to_xxzz_xxxxx, to_xxzz_xxxxy, to_xxzz_xxxxz, to_xxzz_xxxyy, to_xxzz_xxxyz, to_xxzz_xxxzz, to_xxzz_xxy, to_xxzz_xxyyy, to_xxzz_xxyyz, to_xxzz_xxyzz, to_xxzz_xxz, to_xxzz_xxzzz, to_xxzz_xyy, to_xxzz_xyyyy, to_xxzz_xyyyz, to_xxzz_xyyzz, to_xxzz_xyz, to_xxzz_xyzzz, to_xxzz_xzz, to_xxzz_xzzzz, to_xxzz_yyy, to_xxzz_yyz, to_xxzz_yzz, to_xxzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxzz_xxxx[k] = -4.0 * to_xxzz_xxx[k] + 2.0 * to_xxzz_xxxxx[k] * tke_0;

            to_0_x_xxzz_xxxy[k] = -3.0 * to_xxzz_xxy[k] + 2.0 * to_xxzz_xxxxy[k] * tke_0;

            to_0_x_xxzz_xxxz[k] = -3.0 * to_xxzz_xxz[k] + 2.0 * to_xxzz_xxxxz[k] * tke_0;

            to_0_x_xxzz_xxyy[k] = -2.0 * to_xxzz_xyy[k] + 2.0 * to_xxzz_xxxyy[k] * tke_0;

            to_0_x_xxzz_xxyz[k] = -2.0 * to_xxzz_xyz[k] + 2.0 * to_xxzz_xxxyz[k] * tke_0;

            to_0_x_xxzz_xxzz[k] = -2.0 * to_xxzz_xzz[k] + 2.0 * to_xxzz_xxxzz[k] * tke_0;

            to_0_x_xxzz_xyyy[k] = -to_xxzz_yyy[k] + 2.0 * to_xxzz_xxyyy[k] * tke_0;

            to_0_x_xxzz_xyyz[k] = -to_xxzz_yyz[k] + 2.0 * to_xxzz_xxyyz[k] * tke_0;

            to_0_x_xxzz_xyzz[k] = -to_xxzz_yzz[k] + 2.0 * to_xxzz_xxyzz[k] * tke_0;

            to_0_x_xxzz_xzzz[k] = -to_xxzz_zzz[k] + 2.0 * to_xxzz_xxzzz[k] * tke_0;

            to_0_x_xxzz_yyyy[k] = 2.0 * to_xxzz_xyyyy[k] * tke_0;

            to_0_x_xxzz_yyyz[k] = 2.0 * to_xxzz_xyyyz[k] * tke_0;

            to_0_x_xxzz_yyzz[k] = 2.0 * to_xxzz_xyyzz[k] * tke_0;

            to_0_x_xxzz_yzzz[k] = 2.0 * to_xxzz_xyzzz[k] * tke_0;

            to_0_x_xxzz_zzzz[k] = 2.0 * to_xxzz_xzzzz[k] * tke_0;
        }

        // Set up 90-105 components of targeted buffer : GG

        auto to_0_x_xyyy_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 90);

        auto to_0_x_xyyy_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 91);

        auto to_0_x_xyyy_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 92);

        auto to_0_x_xyyy_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 93);

        auto to_0_x_xyyy_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 94);

        auto to_0_x_xyyy_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 95);

        auto to_0_x_xyyy_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 96);

        auto to_0_x_xyyy_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 97);

        auto to_0_x_xyyy_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 98);

        auto to_0_x_xyyy_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 99);

        auto to_0_x_xyyy_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 100);

        auto to_0_x_xyyy_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 101);

        auto to_0_x_xyyy_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 102);

        auto to_0_x_xyyy_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 103);

        auto to_0_x_xyyy_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_0_x_xyyy_xxxx, to_0_x_xyyy_xxxy, to_0_x_xyyy_xxxz, to_0_x_xyyy_xxyy, to_0_x_xyyy_xxyz, to_0_x_xyyy_xxzz, to_0_x_xyyy_xyyy, to_0_x_xyyy_xyyz, to_0_x_xyyy_xyzz, to_0_x_xyyy_xzzz, to_0_x_xyyy_yyyy, to_0_x_xyyy_yyyz, to_0_x_xyyy_yyzz, to_0_x_xyyy_yzzz, to_0_x_xyyy_zzzz, to_xyyy_xxx, to_xyyy_xxxxx, to_xyyy_xxxxy, to_xyyy_xxxxz, to_xyyy_xxxyy, to_xyyy_xxxyz, to_xyyy_xxxzz, to_xyyy_xxy, to_xyyy_xxyyy, to_xyyy_xxyyz, to_xyyy_xxyzz, to_xyyy_xxz, to_xyyy_xxzzz, to_xyyy_xyy, to_xyyy_xyyyy, to_xyyy_xyyyz, to_xyyy_xyyzz, to_xyyy_xyz, to_xyyy_xyzzz, to_xyyy_xzz, to_xyyy_xzzzz, to_xyyy_yyy, to_xyyy_yyz, to_xyyy_yzz, to_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyyy_xxxx[k] = -4.0 * to_xyyy_xxx[k] + 2.0 * to_xyyy_xxxxx[k] * tke_0;

            to_0_x_xyyy_xxxy[k] = -3.0 * to_xyyy_xxy[k] + 2.0 * to_xyyy_xxxxy[k] * tke_0;

            to_0_x_xyyy_xxxz[k] = -3.0 * to_xyyy_xxz[k] + 2.0 * to_xyyy_xxxxz[k] * tke_0;

            to_0_x_xyyy_xxyy[k] = -2.0 * to_xyyy_xyy[k] + 2.0 * to_xyyy_xxxyy[k] * tke_0;

            to_0_x_xyyy_xxyz[k] = -2.0 * to_xyyy_xyz[k] + 2.0 * to_xyyy_xxxyz[k] * tke_0;

            to_0_x_xyyy_xxzz[k] = -2.0 * to_xyyy_xzz[k] + 2.0 * to_xyyy_xxxzz[k] * tke_0;

            to_0_x_xyyy_xyyy[k] = -to_xyyy_yyy[k] + 2.0 * to_xyyy_xxyyy[k] * tke_0;

            to_0_x_xyyy_xyyz[k] = -to_xyyy_yyz[k] + 2.0 * to_xyyy_xxyyz[k] * tke_0;

            to_0_x_xyyy_xyzz[k] = -to_xyyy_yzz[k] + 2.0 * to_xyyy_xxyzz[k] * tke_0;

            to_0_x_xyyy_xzzz[k] = -to_xyyy_zzz[k] + 2.0 * to_xyyy_xxzzz[k] * tke_0;

            to_0_x_xyyy_yyyy[k] = 2.0 * to_xyyy_xyyyy[k] * tke_0;

            to_0_x_xyyy_yyyz[k] = 2.0 * to_xyyy_xyyyz[k] * tke_0;

            to_0_x_xyyy_yyzz[k] = 2.0 * to_xyyy_xyyzz[k] * tke_0;

            to_0_x_xyyy_yzzz[k] = 2.0 * to_xyyy_xyzzz[k] * tke_0;

            to_0_x_xyyy_zzzz[k] = 2.0 * to_xyyy_xzzzz[k] * tke_0;
        }

        // Set up 105-120 components of targeted buffer : GG

        auto to_0_x_xyyz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 105);

        auto to_0_x_xyyz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 106);

        auto to_0_x_xyyz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 107);

        auto to_0_x_xyyz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 108);

        auto to_0_x_xyyz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 109);

        auto to_0_x_xyyz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 110);

        auto to_0_x_xyyz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 111);

        auto to_0_x_xyyz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 112);

        auto to_0_x_xyyz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 113);

        auto to_0_x_xyyz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 114);

        auto to_0_x_xyyz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 115);

        auto to_0_x_xyyz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 116);

        auto to_0_x_xyyz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 117);

        auto to_0_x_xyyz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 118);

        auto to_0_x_xyyz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_0_x_xyyz_xxxx, to_0_x_xyyz_xxxy, to_0_x_xyyz_xxxz, to_0_x_xyyz_xxyy, to_0_x_xyyz_xxyz, to_0_x_xyyz_xxzz, to_0_x_xyyz_xyyy, to_0_x_xyyz_xyyz, to_0_x_xyyz_xyzz, to_0_x_xyyz_xzzz, to_0_x_xyyz_yyyy, to_0_x_xyyz_yyyz, to_0_x_xyyz_yyzz, to_0_x_xyyz_yzzz, to_0_x_xyyz_zzzz, to_xyyz_xxx, to_xyyz_xxxxx, to_xyyz_xxxxy, to_xyyz_xxxxz, to_xyyz_xxxyy, to_xyyz_xxxyz, to_xyyz_xxxzz, to_xyyz_xxy, to_xyyz_xxyyy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xxzzz, to_xyyz_xyy, to_xyyz_xyyyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_xzzzz, to_xyyz_yyy, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyyz_xxxx[k] = -4.0 * to_xyyz_xxx[k] + 2.0 * to_xyyz_xxxxx[k] * tke_0;

            to_0_x_xyyz_xxxy[k] = -3.0 * to_xyyz_xxy[k] + 2.0 * to_xyyz_xxxxy[k] * tke_0;

            to_0_x_xyyz_xxxz[k] = -3.0 * to_xyyz_xxz[k] + 2.0 * to_xyyz_xxxxz[k] * tke_0;

            to_0_x_xyyz_xxyy[k] = -2.0 * to_xyyz_xyy[k] + 2.0 * to_xyyz_xxxyy[k] * tke_0;

            to_0_x_xyyz_xxyz[k] = -2.0 * to_xyyz_xyz[k] + 2.0 * to_xyyz_xxxyz[k] * tke_0;

            to_0_x_xyyz_xxzz[k] = -2.0 * to_xyyz_xzz[k] + 2.0 * to_xyyz_xxxzz[k] * tke_0;

            to_0_x_xyyz_xyyy[k] = -to_xyyz_yyy[k] + 2.0 * to_xyyz_xxyyy[k] * tke_0;

            to_0_x_xyyz_xyyz[k] = -to_xyyz_yyz[k] + 2.0 * to_xyyz_xxyyz[k] * tke_0;

            to_0_x_xyyz_xyzz[k] = -to_xyyz_yzz[k] + 2.0 * to_xyyz_xxyzz[k] * tke_0;

            to_0_x_xyyz_xzzz[k] = -to_xyyz_zzz[k] + 2.0 * to_xyyz_xxzzz[k] * tke_0;

            to_0_x_xyyz_yyyy[k] = 2.0 * to_xyyz_xyyyy[k] * tke_0;

            to_0_x_xyyz_yyyz[k] = 2.0 * to_xyyz_xyyyz[k] * tke_0;

            to_0_x_xyyz_yyzz[k] = 2.0 * to_xyyz_xyyzz[k] * tke_0;

            to_0_x_xyyz_yzzz[k] = 2.0 * to_xyyz_xyzzz[k] * tke_0;

            to_0_x_xyyz_zzzz[k] = 2.0 * to_xyyz_xzzzz[k] * tke_0;
        }

        // Set up 120-135 components of targeted buffer : GG

        auto to_0_x_xyzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 120);

        auto to_0_x_xyzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 121);

        auto to_0_x_xyzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 122);

        auto to_0_x_xyzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 123);

        auto to_0_x_xyzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 124);

        auto to_0_x_xyzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 125);

        auto to_0_x_xyzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 126);

        auto to_0_x_xyzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 127);

        auto to_0_x_xyzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 128);

        auto to_0_x_xyzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 129);

        auto to_0_x_xyzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 130);

        auto to_0_x_xyzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 131);

        auto to_0_x_xyzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 132);

        auto to_0_x_xyzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 133);

        auto to_0_x_xyzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_0_x_xyzz_xxxx, to_0_x_xyzz_xxxy, to_0_x_xyzz_xxxz, to_0_x_xyzz_xxyy, to_0_x_xyzz_xxyz, to_0_x_xyzz_xxzz, to_0_x_xyzz_xyyy, to_0_x_xyzz_xyyz, to_0_x_xyzz_xyzz, to_0_x_xyzz_xzzz, to_0_x_xyzz_yyyy, to_0_x_xyzz_yyyz, to_0_x_xyzz_yyzz, to_0_x_xyzz_yzzz, to_0_x_xyzz_zzzz, to_xyzz_xxx, to_xyzz_xxxxx, to_xyzz_xxxxy, to_xyzz_xxxxz, to_xyzz_xxxyy, to_xyzz_xxxyz, to_xyzz_xxxzz, to_xyzz_xxy, to_xyzz_xxyyy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xxzzz, to_xyzz_xyy, to_xyzz_xyyyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_xzzzz, to_xyzz_yyy, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyzz_xxxx[k] = -4.0 * to_xyzz_xxx[k] + 2.0 * to_xyzz_xxxxx[k] * tke_0;

            to_0_x_xyzz_xxxy[k] = -3.0 * to_xyzz_xxy[k] + 2.0 * to_xyzz_xxxxy[k] * tke_0;

            to_0_x_xyzz_xxxz[k] = -3.0 * to_xyzz_xxz[k] + 2.0 * to_xyzz_xxxxz[k] * tke_0;

            to_0_x_xyzz_xxyy[k] = -2.0 * to_xyzz_xyy[k] + 2.0 * to_xyzz_xxxyy[k] * tke_0;

            to_0_x_xyzz_xxyz[k] = -2.0 * to_xyzz_xyz[k] + 2.0 * to_xyzz_xxxyz[k] * tke_0;

            to_0_x_xyzz_xxzz[k] = -2.0 * to_xyzz_xzz[k] + 2.0 * to_xyzz_xxxzz[k] * tke_0;

            to_0_x_xyzz_xyyy[k] = -to_xyzz_yyy[k] + 2.0 * to_xyzz_xxyyy[k] * tke_0;

            to_0_x_xyzz_xyyz[k] = -to_xyzz_yyz[k] + 2.0 * to_xyzz_xxyyz[k] * tke_0;

            to_0_x_xyzz_xyzz[k] = -to_xyzz_yzz[k] + 2.0 * to_xyzz_xxyzz[k] * tke_0;

            to_0_x_xyzz_xzzz[k] = -to_xyzz_zzz[k] + 2.0 * to_xyzz_xxzzz[k] * tke_0;

            to_0_x_xyzz_yyyy[k] = 2.0 * to_xyzz_xyyyy[k] * tke_0;

            to_0_x_xyzz_yyyz[k] = 2.0 * to_xyzz_xyyyz[k] * tke_0;

            to_0_x_xyzz_yyzz[k] = 2.0 * to_xyzz_xyyzz[k] * tke_0;

            to_0_x_xyzz_yzzz[k] = 2.0 * to_xyzz_xyzzz[k] * tke_0;

            to_0_x_xyzz_zzzz[k] = 2.0 * to_xyzz_xzzzz[k] * tke_0;
        }

        // Set up 135-150 components of targeted buffer : GG

        auto to_0_x_xzzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 135);

        auto to_0_x_xzzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 136);

        auto to_0_x_xzzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 137);

        auto to_0_x_xzzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 138);

        auto to_0_x_xzzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 139);

        auto to_0_x_xzzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 140);

        auto to_0_x_xzzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 141);

        auto to_0_x_xzzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 142);

        auto to_0_x_xzzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 143);

        auto to_0_x_xzzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 144);

        auto to_0_x_xzzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 145);

        auto to_0_x_xzzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 146);

        auto to_0_x_xzzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 147);

        auto to_0_x_xzzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 148);

        auto to_0_x_xzzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_0_x_xzzz_xxxx, to_0_x_xzzz_xxxy, to_0_x_xzzz_xxxz, to_0_x_xzzz_xxyy, to_0_x_xzzz_xxyz, to_0_x_xzzz_xxzz, to_0_x_xzzz_xyyy, to_0_x_xzzz_xyyz, to_0_x_xzzz_xyzz, to_0_x_xzzz_xzzz, to_0_x_xzzz_yyyy, to_0_x_xzzz_yyyz, to_0_x_xzzz_yyzz, to_0_x_xzzz_yzzz, to_0_x_xzzz_zzzz, to_xzzz_xxx, to_xzzz_xxxxx, to_xzzz_xxxxy, to_xzzz_xxxxz, to_xzzz_xxxyy, to_xzzz_xxxyz, to_xzzz_xxxzz, to_xzzz_xxy, to_xzzz_xxyyy, to_xzzz_xxyyz, to_xzzz_xxyzz, to_xzzz_xxz, to_xzzz_xxzzz, to_xzzz_xyy, to_xzzz_xyyyy, to_xzzz_xyyyz, to_xzzz_xyyzz, to_xzzz_xyz, to_xzzz_xyzzz, to_xzzz_xzz, to_xzzz_xzzzz, to_xzzz_yyy, to_xzzz_yyz, to_xzzz_yzz, to_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xzzz_xxxx[k] = -4.0 * to_xzzz_xxx[k] + 2.0 * to_xzzz_xxxxx[k] * tke_0;

            to_0_x_xzzz_xxxy[k] = -3.0 * to_xzzz_xxy[k] + 2.0 * to_xzzz_xxxxy[k] * tke_0;

            to_0_x_xzzz_xxxz[k] = -3.0 * to_xzzz_xxz[k] + 2.0 * to_xzzz_xxxxz[k] * tke_0;

            to_0_x_xzzz_xxyy[k] = -2.0 * to_xzzz_xyy[k] + 2.0 * to_xzzz_xxxyy[k] * tke_0;

            to_0_x_xzzz_xxyz[k] = -2.0 * to_xzzz_xyz[k] + 2.0 * to_xzzz_xxxyz[k] * tke_0;

            to_0_x_xzzz_xxzz[k] = -2.0 * to_xzzz_xzz[k] + 2.0 * to_xzzz_xxxzz[k] * tke_0;

            to_0_x_xzzz_xyyy[k] = -to_xzzz_yyy[k] + 2.0 * to_xzzz_xxyyy[k] * tke_0;

            to_0_x_xzzz_xyyz[k] = -to_xzzz_yyz[k] + 2.0 * to_xzzz_xxyyz[k] * tke_0;

            to_0_x_xzzz_xyzz[k] = -to_xzzz_yzz[k] + 2.0 * to_xzzz_xxyzz[k] * tke_0;

            to_0_x_xzzz_xzzz[k] = -to_xzzz_zzz[k] + 2.0 * to_xzzz_xxzzz[k] * tke_0;

            to_0_x_xzzz_yyyy[k] = 2.0 * to_xzzz_xyyyy[k] * tke_0;

            to_0_x_xzzz_yyyz[k] = 2.0 * to_xzzz_xyyyz[k] * tke_0;

            to_0_x_xzzz_yyzz[k] = 2.0 * to_xzzz_xyyzz[k] * tke_0;

            to_0_x_xzzz_yzzz[k] = 2.0 * to_xzzz_xyzzz[k] * tke_0;

            to_0_x_xzzz_zzzz[k] = 2.0 * to_xzzz_xzzzz[k] * tke_0;
        }

        // Set up 150-165 components of targeted buffer : GG

        auto to_0_x_yyyy_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 150);

        auto to_0_x_yyyy_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 151);

        auto to_0_x_yyyy_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 152);

        auto to_0_x_yyyy_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 153);

        auto to_0_x_yyyy_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 154);

        auto to_0_x_yyyy_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 155);

        auto to_0_x_yyyy_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 156);

        auto to_0_x_yyyy_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 157);

        auto to_0_x_yyyy_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 158);

        auto to_0_x_yyyy_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 159);

        auto to_0_x_yyyy_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 160);

        auto to_0_x_yyyy_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 161);

        auto to_0_x_yyyy_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 162);

        auto to_0_x_yyyy_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 163);

        auto to_0_x_yyyy_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_0_x_yyyy_xxxx, to_0_x_yyyy_xxxy, to_0_x_yyyy_xxxz, to_0_x_yyyy_xxyy, to_0_x_yyyy_xxyz, to_0_x_yyyy_xxzz, to_0_x_yyyy_xyyy, to_0_x_yyyy_xyyz, to_0_x_yyyy_xyzz, to_0_x_yyyy_xzzz, to_0_x_yyyy_yyyy, to_0_x_yyyy_yyyz, to_0_x_yyyy_yyzz, to_0_x_yyyy_yzzz, to_0_x_yyyy_zzzz, to_yyyy_xxx, to_yyyy_xxxxx, to_yyyy_xxxxy, to_yyyy_xxxxz, to_yyyy_xxxyy, to_yyyy_xxxyz, to_yyyy_xxxzz, to_yyyy_xxy, to_yyyy_xxyyy, to_yyyy_xxyyz, to_yyyy_xxyzz, to_yyyy_xxz, to_yyyy_xxzzz, to_yyyy_xyy, to_yyyy_xyyyy, to_yyyy_xyyyz, to_yyyy_xyyzz, to_yyyy_xyz, to_yyyy_xyzzz, to_yyyy_xzz, to_yyyy_xzzzz, to_yyyy_yyy, to_yyyy_yyz, to_yyyy_yzz, to_yyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyyy_xxxx[k] = -4.0 * to_yyyy_xxx[k] + 2.0 * to_yyyy_xxxxx[k] * tke_0;

            to_0_x_yyyy_xxxy[k] = -3.0 * to_yyyy_xxy[k] + 2.0 * to_yyyy_xxxxy[k] * tke_0;

            to_0_x_yyyy_xxxz[k] = -3.0 * to_yyyy_xxz[k] + 2.0 * to_yyyy_xxxxz[k] * tke_0;

            to_0_x_yyyy_xxyy[k] = -2.0 * to_yyyy_xyy[k] + 2.0 * to_yyyy_xxxyy[k] * tke_0;

            to_0_x_yyyy_xxyz[k] = -2.0 * to_yyyy_xyz[k] + 2.0 * to_yyyy_xxxyz[k] * tke_0;

            to_0_x_yyyy_xxzz[k] = -2.0 * to_yyyy_xzz[k] + 2.0 * to_yyyy_xxxzz[k] * tke_0;

            to_0_x_yyyy_xyyy[k] = -to_yyyy_yyy[k] + 2.0 * to_yyyy_xxyyy[k] * tke_0;

            to_0_x_yyyy_xyyz[k] = -to_yyyy_yyz[k] + 2.0 * to_yyyy_xxyyz[k] * tke_0;

            to_0_x_yyyy_xyzz[k] = -to_yyyy_yzz[k] + 2.0 * to_yyyy_xxyzz[k] * tke_0;

            to_0_x_yyyy_xzzz[k] = -to_yyyy_zzz[k] + 2.0 * to_yyyy_xxzzz[k] * tke_0;

            to_0_x_yyyy_yyyy[k] = 2.0 * to_yyyy_xyyyy[k] * tke_0;

            to_0_x_yyyy_yyyz[k] = 2.0 * to_yyyy_xyyyz[k] * tke_0;

            to_0_x_yyyy_yyzz[k] = 2.0 * to_yyyy_xyyzz[k] * tke_0;

            to_0_x_yyyy_yzzz[k] = 2.0 * to_yyyy_xyzzz[k] * tke_0;

            to_0_x_yyyy_zzzz[k] = 2.0 * to_yyyy_xzzzz[k] * tke_0;
        }

        // Set up 165-180 components of targeted buffer : GG

        auto to_0_x_yyyz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 165);

        auto to_0_x_yyyz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 166);

        auto to_0_x_yyyz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 167);

        auto to_0_x_yyyz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 168);

        auto to_0_x_yyyz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 169);

        auto to_0_x_yyyz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 170);

        auto to_0_x_yyyz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 171);

        auto to_0_x_yyyz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 172);

        auto to_0_x_yyyz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 173);

        auto to_0_x_yyyz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 174);

        auto to_0_x_yyyz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 175);

        auto to_0_x_yyyz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 176);

        auto to_0_x_yyyz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 177);

        auto to_0_x_yyyz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 178);

        auto to_0_x_yyyz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_0_x_yyyz_xxxx, to_0_x_yyyz_xxxy, to_0_x_yyyz_xxxz, to_0_x_yyyz_xxyy, to_0_x_yyyz_xxyz, to_0_x_yyyz_xxzz, to_0_x_yyyz_xyyy, to_0_x_yyyz_xyyz, to_0_x_yyyz_xyzz, to_0_x_yyyz_xzzz, to_0_x_yyyz_yyyy, to_0_x_yyyz_yyyz, to_0_x_yyyz_yyzz, to_0_x_yyyz_yzzz, to_0_x_yyyz_zzzz, to_yyyz_xxx, to_yyyz_xxxxx, to_yyyz_xxxxy, to_yyyz_xxxxz, to_yyyz_xxxyy, to_yyyz_xxxyz, to_yyyz_xxxzz, to_yyyz_xxy, to_yyyz_xxyyy, to_yyyz_xxyyz, to_yyyz_xxyzz, to_yyyz_xxz, to_yyyz_xxzzz, to_yyyz_xyy, to_yyyz_xyyyy, to_yyyz_xyyyz, to_yyyz_xyyzz, to_yyyz_xyz, to_yyyz_xyzzz, to_yyyz_xzz, to_yyyz_xzzzz, to_yyyz_yyy, to_yyyz_yyz, to_yyyz_yzz, to_yyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyyz_xxxx[k] = -4.0 * to_yyyz_xxx[k] + 2.0 * to_yyyz_xxxxx[k] * tke_0;

            to_0_x_yyyz_xxxy[k] = -3.0 * to_yyyz_xxy[k] + 2.0 * to_yyyz_xxxxy[k] * tke_0;

            to_0_x_yyyz_xxxz[k] = -3.0 * to_yyyz_xxz[k] + 2.0 * to_yyyz_xxxxz[k] * tke_0;

            to_0_x_yyyz_xxyy[k] = -2.0 * to_yyyz_xyy[k] + 2.0 * to_yyyz_xxxyy[k] * tke_0;

            to_0_x_yyyz_xxyz[k] = -2.0 * to_yyyz_xyz[k] + 2.0 * to_yyyz_xxxyz[k] * tke_0;

            to_0_x_yyyz_xxzz[k] = -2.0 * to_yyyz_xzz[k] + 2.0 * to_yyyz_xxxzz[k] * tke_0;

            to_0_x_yyyz_xyyy[k] = -to_yyyz_yyy[k] + 2.0 * to_yyyz_xxyyy[k] * tke_0;

            to_0_x_yyyz_xyyz[k] = -to_yyyz_yyz[k] + 2.0 * to_yyyz_xxyyz[k] * tke_0;

            to_0_x_yyyz_xyzz[k] = -to_yyyz_yzz[k] + 2.0 * to_yyyz_xxyzz[k] * tke_0;

            to_0_x_yyyz_xzzz[k] = -to_yyyz_zzz[k] + 2.0 * to_yyyz_xxzzz[k] * tke_0;

            to_0_x_yyyz_yyyy[k] = 2.0 * to_yyyz_xyyyy[k] * tke_0;

            to_0_x_yyyz_yyyz[k] = 2.0 * to_yyyz_xyyyz[k] * tke_0;

            to_0_x_yyyz_yyzz[k] = 2.0 * to_yyyz_xyyzz[k] * tke_0;

            to_0_x_yyyz_yzzz[k] = 2.0 * to_yyyz_xyzzz[k] * tke_0;

            to_0_x_yyyz_zzzz[k] = 2.0 * to_yyyz_xzzzz[k] * tke_0;
        }

        // Set up 180-195 components of targeted buffer : GG

        auto to_0_x_yyzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 180);

        auto to_0_x_yyzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 181);

        auto to_0_x_yyzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 182);

        auto to_0_x_yyzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 183);

        auto to_0_x_yyzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 184);

        auto to_0_x_yyzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 185);

        auto to_0_x_yyzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 186);

        auto to_0_x_yyzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 187);

        auto to_0_x_yyzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 188);

        auto to_0_x_yyzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 189);

        auto to_0_x_yyzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 190);

        auto to_0_x_yyzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 191);

        auto to_0_x_yyzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 192);

        auto to_0_x_yyzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 193);

        auto to_0_x_yyzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_0_x_yyzz_xxxx, to_0_x_yyzz_xxxy, to_0_x_yyzz_xxxz, to_0_x_yyzz_xxyy, to_0_x_yyzz_xxyz, to_0_x_yyzz_xxzz, to_0_x_yyzz_xyyy, to_0_x_yyzz_xyyz, to_0_x_yyzz_xyzz, to_0_x_yyzz_xzzz, to_0_x_yyzz_yyyy, to_0_x_yyzz_yyyz, to_0_x_yyzz_yyzz, to_0_x_yyzz_yzzz, to_0_x_yyzz_zzzz, to_yyzz_xxx, to_yyzz_xxxxx, to_yyzz_xxxxy, to_yyzz_xxxxz, to_yyzz_xxxyy, to_yyzz_xxxyz, to_yyzz_xxxzz, to_yyzz_xxy, to_yyzz_xxyyy, to_yyzz_xxyyz, to_yyzz_xxyzz, to_yyzz_xxz, to_yyzz_xxzzz, to_yyzz_xyy, to_yyzz_xyyyy, to_yyzz_xyyyz, to_yyzz_xyyzz, to_yyzz_xyz, to_yyzz_xyzzz, to_yyzz_xzz, to_yyzz_xzzzz, to_yyzz_yyy, to_yyzz_yyz, to_yyzz_yzz, to_yyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyzz_xxxx[k] = -4.0 * to_yyzz_xxx[k] + 2.0 * to_yyzz_xxxxx[k] * tke_0;

            to_0_x_yyzz_xxxy[k] = -3.0 * to_yyzz_xxy[k] + 2.0 * to_yyzz_xxxxy[k] * tke_0;

            to_0_x_yyzz_xxxz[k] = -3.0 * to_yyzz_xxz[k] + 2.0 * to_yyzz_xxxxz[k] * tke_0;

            to_0_x_yyzz_xxyy[k] = -2.0 * to_yyzz_xyy[k] + 2.0 * to_yyzz_xxxyy[k] * tke_0;

            to_0_x_yyzz_xxyz[k] = -2.0 * to_yyzz_xyz[k] + 2.0 * to_yyzz_xxxyz[k] * tke_0;

            to_0_x_yyzz_xxzz[k] = -2.0 * to_yyzz_xzz[k] + 2.0 * to_yyzz_xxxzz[k] * tke_0;

            to_0_x_yyzz_xyyy[k] = -to_yyzz_yyy[k] + 2.0 * to_yyzz_xxyyy[k] * tke_0;

            to_0_x_yyzz_xyyz[k] = -to_yyzz_yyz[k] + 2.0 * to_yyzz_xxyyz[k] * tke_0;

            to_0_x_yyzz_xyzz[k] = -to_yyzz_yzz[k] + 2.0 * to_yyzz_xxyzz[k] * tke_0;

            to_0_x_yyzz_xzzz[k] = -to_yyzz_zzz[k] + 2.0 * to_yyzz_xxzzz[k] * tke_0;

            to_0_x_yyzz_yyyy[k] = 2.0 * to_yyzz_xyyyy[k] * tke_0;

            to_0_x_yyzz_yyyz[k] = 2.0 * to_yyzz_xyyyz[k] * tke_0;

            to_0_x_yyzz_yyzz[k] = 2.0 * to_yyzz_xyyzz[k] * tke_0;

            to_0_x_yyzz_yzzz[k] = 2.0 * to_yyzz_xyzzz[k] * tke_0;

            to_0_x_yyzz_zzzz[k] = 2.0 * to_yyzz_xzzzz[k] * tke_0;
        }

        // Set up 195-210 components of targeted buffer : GG

        auto to_0_x_yzzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 195);

        auto to_0_x_yzzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 196);

        auto to_0_x_yzzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 197);

        auto to_0_x_yzzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 198);

        auto to_0_x_yzzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 199);

        auto to_0_x_yzzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 200);

        auto to_0_x_yzzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 201);

        auto to_0_x_yzzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 202);

        auto to_0_x_yzzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 203);

        auto to_0_x_yzzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 204);

        auto to_0_x_yzzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 205);

        auto to_0_x_yzzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 206);

        auto to_0_x_yzzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 207);

        auto to_0_x_yzzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 208);

        auto to_0_x_yzzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_0_x_yzzz_xxxx, to_0_x_yzzz_xxxy, to_0_x_yzzz_xxxz, to_0_x_yzzz_xxyy, to_0_x_yzzz_xxyz, to_0_x_yzzz_xxzz, to_0_x_yzzz_xyyy, to_0_x_yzzz_xyyz, to_0_x_yzzz_xyzz, to_0_x_yzzz_xzzz, to_0_x_yzzz_yyyy, to_0_x_yzzz_yyyz, to_0_x_yzzz_yyzz, to_0_x_yzzz_yzzz, to_0_x_yzzz_zzzz, to_yzzz_xxx, to_yzzz_xxxxx, to_yzzz_xxxxy, to_yzzz_xxxxz, to_yzzz_xxxyy, to_yzzz_xxxyz, to_yzzz_xxxzz, to_yzzz_xxy, to_yzzz_xxyyy, to_yzzz_xxyyz, to_yzzz_xxyzz, to_yzzz_xxz, to_yzzz_xxzzz, to_yzzz_xyy, to_yzzz_xyyyy, to_yzzz_xyyyz, to_yzzz_xyyzz, to_yzzz_xyz, to_yzzz_xyzzz, to_yzzz_xzz, to_yzzz_xzzzz, to_yzzz_yyy, to_yzzz_yyz, to_yzzz_yzz, to_yzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yzzz_xxxx[k] = -4.0 * to_yzzz_xxx[k] + 2.0 * to_yzzz_xxxxx[k] * tke_0;

            to_0_x_yzzz_xxxy[k] = -3.0 * to_yzzz_xxy[k] + 2.0 * to_yzzz_xxxxy[k] * tke_0;

            to_0_x_yzzz_xxxz[k] = -3.0 * to_yzzz_xxz[k] + 2.0 * to_yzzz_xxxxz[k] * tke_0;

            to_0_x_yzzz_xxyy[k] = -2.0 * to_yzzz_xyy[k] + 2.0 * to_yzzz_xxxyy[k] * tke_0;

            to_0_x_yzzz_xxyz[k] = -2.0 * to_yzzz_xyz[k] + 2.0 * to_yzzz_xxxyz[k] * tke_0;

            to_0_x_yzzz_xxzz[k] = -2.0 * to_yzzz_xzz[k] + 2.0 * to_yzzz_xxxzz[k] * tke_0;

            to_0_x_yzzz_xyyy[k] = -to_yzzz_yyy[k] + 2.0 * to_yzzz_xxyyy[k] * tke_0;

            to_0_x_yzzz_xyyz[k] = -to_yzzz_yyz[k] + 2.0 * to_yzzz_xxyyz[k] * tke_0;

            to_0_x_yzzz_xyzz[k] = -to_yzzz_yzz[k] + 2.0 * to_yzzz_xxyzz[k] * tke_0;

            to_0_x_yzzz_xzzz[k] = -to_yzzz_zzz[k] + 2.0 * to_yzzz_xxzzz[k] * tke_0;

            to_0_x_yzzz_yyyy[k] = 2.0 * to_yzzz_xyyyy[k] * tke_0;

            to_0_x_yzzz_yyyz[k] = 2.0 * to_yzzz_xyyyz[k] * tke_0;

            to_0_x_yzzz_yyzz[k] = 2.0 * to_yzzz_xyyzz[k] * tke_0;

            to_0_x_yzzz_yzzz[k] = 2.0 * to_yzzz_xyzzz[k] * tke_0;

            to_0_x_yzzz_zzzz[k] = 2.0 * to_yzzz_xzzzz[k] * tke_0;
        }

        // Set up 210-225 components of targeted buffer : GG

        auto to_0_x_zzzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 210);

        auto to_0_x_zzzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 211);

        auto to_0_x_zzzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 212);

        auto to_0_x_zzzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 213);

        auto to_0_x_zzzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 214);

        auto to_0_x_zzzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 215);

        auto to_0_x_zzzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 216);

        auto to_0_x_zzzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 217);

        auto to_0_x_zzzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 218);

        auto to_0_x_zzzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 219);

        auto to_0_x_zzzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 220);

        auto to_0_x_zzzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 221);

        auto to_0_x_zzzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 222);

        auto to_0_x_zzzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 223);

        auto to_0_x_zzzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 0 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_0_x_zzzz_xxxx, to_0_x_zzzz_xxxy, to_0_x_zzzz_xxxz, to_0_x_zzzz_xxyy, to_0_x_zzzz_xxyz, to_0_x_zzzz_xxzz, to_0_x_zzzz_xyyy, to_0_x_zzzz_xyyz, to_0_x_zzzz_xyzz, to_0_x_zzzz_xzzz, to_0_x_zzzz_yyyy, to_0_x_zzzz_yyyz, to_0_x_zzzz_yyzz, to_0_x_zzzz_yzzz, to_0_x_zzzz_zzzz, to_zzzz_xxx, to_zzzz_xxxxx, to_zzzz_xxxxy, to_zzzz_xxxxz, to_zzzz_xxxyy, to_zzzz_xxxyz, to_zzzz_xxxzz, to_zzzz_xxy, to_zzzz_xxyyy, to_zzzz_xxyyz, to_zzzz_xxyzz, to_zzzz_xxz, to_zzzz_xxzzz, to_zzzz_xyy, to_zzzz_xyyyy, to_zzzz_xyyyz, to_zzzz_xyyzz, to_zzzz_xyz, to_zzzz_xyzzz, to_zzzz_xzz, to_zzzz_xzzzz, to_zzzz_yyy, to_zzzz_yyz, to_zzzz_yzz, to_zzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zzzz_xxxx[k] = -4.0 * to_zzzz_xxx[k] + 2.0 * to_zzzz_xxxxx[k] * tke_0;

            to_0_x_zzzz_xxxy[k] = -3.0 * to_zzzz_xxy[k] + 2.0 * to_zzzz_xxxxy[k] * tke_0;

            to_0_x_zzzz_xxxz[k] = -3.0 * to_zzzz_xxz[k] + 2.0 * to_zzzz_xxxxz[k] * tke_0;

            to_0_x_zzzz_xxyy[k] = -2.0 * to_zzzz_xyy[k] + 2.0 * to_zzzz_xxxyy[k] * tke_0;

            to_0_x_zzzz_xxyz[k] = -2.0 * to_zzzz_xyz[k] + 2.0 * to_zzzz_xxxyz[k] * tke_0;

            to_0_x_zzzz_xxzz[k] = -2.0 * to_zzzz_xzz[k] + 2.0 * to_zzzz_xxxzz[k] * tke_0;

            to_0_x_zzzz_xyyy[k] = -to_zzzz_yyy[k] + 2.0 * to_zzzz_xxyyy[k] * tke_0;

            to_0_x_zzzz_xyyz[k] = -to_zzzz_yyz[k] + 2.0 * to_zzzz_xxyyz[k] * tke_0;

            to_0_x_zzzz_xyzz[k] = -to_zzzz_yzz[k] + 2.0 * to_zzzz_xxyzz[k] * tke_0;

            to_0_x_zzzz_xzzz[k] = -to_zzzz_zzz[k] + 2.0 * to_zzzz_xxzzz[k] * tke_0;

            to_0_x_zzzz_yyyy[k] = 2.0 * to_zzzz_xyyyy[k] * tke_0;

            to_0_x_zzzz_yyyz[k] = 2.0 * to_zzzz_xyyyz[k] * tke_0;

            to_0_x_zzzz_yyzz[k] = 2.0 * to_zzzz_xyyzz[k] * tke_0;

            to_0_x_zzzz_yzzz[k] = 2.0 * to_zzzz_xyzzz[k] * tke_0;

            to_0_x_zzzz_zzzz[k] = 2.0 * to_zzzz_xzzzz[k] * tke_0;
        }

        // Set up 225-240 components of targeted buffer : GG

        auto to_0_y_xxxx_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 0);

        auto to_0_y_xxxx_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 1);

        auto to_0_y_xxxx_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 2);

        auto to_0_y_xxxx_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 3);

        auto to_0_y_xxxx_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 4);

        auto to_0_y_xxxx_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 5);

        auto to_0_y_xxxx_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 6);

        auto to_0_y_xxxx_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 7);

        auto to_0_y_xxxx_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 8);

        auto to_0_y_xxxx_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 9);

        auto to_0_y_xxxx_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 10);

        auto to_0_y_xxxx_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 11);

        auto to_0_y_xxxx_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 12);

        auto to_0_y_xxxx_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 13);

        auto to_0_y_xxxx_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_0_y_xxxx_xxxx, to_0_y_xxxx_xxxy, to_0_y_xxxx_xxxz, to_0_y_xxxx_xxyy, to_0_y_xxxx_xxyz, to_0_y_xxxx_xxzz, to_0_y_xxxx_xyyy, to_0_y_xxxx_xyyz, to_0_y_xxxx_xyzz, to_0_y_xxxx_xzzz, to_0_y_xxxx_yyyy, to_0_y_xxxx_yyyz, to_0_y_xxxx_yyzz, to_0_y_xxxx_yzzz, to_0_y_xxxx_zzzz, to_xxxx_xxx, to_xxxx_xxxxy, to_xxxx_xxxyy, to_xxxx_xxxyz, to_xxxx_xxy, to_xxxx_xxyyy, to_xxxx_xxyyz, to_xxxx_xxyzz, to_xxxx_xxz, to_xxxx_xyy, to_xxxx_xyyyy, to_xxxx_xyyyz, to_xxxx_xyyzz, to_xxxx_xyz, to_xxxx_xyzzz, to_xxxx_xzz, to_xxxx_yyy, to_xxxx_yyyyy, to_xxxx_yyyyz, to_xxxx_yyyzz, to_xxxx_yyz, to_xxxx_yyzzz, to_xxxx_yzz, to_xxxx_yzzzz, to_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxx_xxxx[k] = 2.0 * to_xxxx_xxxxy[k] * tke_0;

            to_0_y_xxxx_xxxy[k] = -to_xxxx_xxx[k] + 2.0 * to_xxxx_xxxyy[k] * tke_0;

            to_0_y_xxxx_xxxz[k] = 2.0 * to_xxxx_xxxyz[k] * tke_0;

            to_0_y_xxxx_xxyy[k] = -2.0 * to_xxxx_xxy[k] + 2.0 * to_xxxx_xxyyy[k] * tke_0;

            to_0_y_xxxx_xxyz[k] = -to_xxxx_xxz[k] + 2.0 * to_xxxx_xxyyz[k] * tke_0;

            to_0_y_xxxx_xxzz[k] = 2.0 * to_xxxx_xxyzz[k] * tke_0;

            to_0_y_xxxx_xyyy[k] = -3.0 * to_xxxx_xyy[k] + 2.0 * to_xxxx_xyyyy[k] * tke_0;

            to_0_y_xxxx_xyyz[k] = -2.0 * to_xxxx_xyz[k] + 2.0 * to_xxxx_xyyyz[k] * tke_0;

            to_0_y_xxxx_xyzz[k] = -to_xxxx_xzz[k] + 2.0 * to_xxxx_xyyzz[k] * tke_0;

            to_0_y_xxxx_xzzz[k] = 2.0 * to_xxxx_xyzzz[k] * tke_0;

            to_0_y_xxxx_yyyy[k] = -4.0 * to_xxxx_yyy[k] + 2.0 * to_xxxx_yyyyy[k] * tke_0;

            to_0_y_xxxx_yyyz[k] = -3.0 * to_xxxx_yyz[k] + 2.0 * to_xxxx_yyyyz[k] * tke_0;

            to_0_y_xxxx_yyzz[k] = -2.0 * to_xxxx_yzz[k] + 2.0 * to_xxxx_yyyzz[k] * tke_0;

            to_0_y_xxxx_yzzz[k] = -to_xxxx_zzz[k] + 2.0 * to_xxxx_yyzzz[k] * tke_0;

            to_0_y_xxxx_zzzz[k] = 2.0 * to_xxxx_yzzzz[k] * tke_0;
        }

        // Set up 240-255 components of targeted buffer : GG

        auto to_0_y_xxxy_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 15);

        auto to_0_y_xxxy_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 16);

        auto to_0_y_xxxy_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 17);

        auto to_0_y_xxxy_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 18);

        auto to_0_y_xxxy_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 19);

        auto to_0_y_xxxy_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 20);

        auto to_0_y_xxxy_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 21);

        auto to_0_y_xxxy_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 22);

        auto to_0_y_xxxy_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 23);

        auto to_0_y_xxxy_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 24);

        auto to_0_y_xxxy_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 25);

        auto to_0_y_xxxy_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 26);

        auto to_0_y_xxxy_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 27);

        auto to_0_y_xxxy_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 28);

        auto to_0_y_xxxy_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_0_y_xxxy_xxxx, to_0_y_xxxy_xxxy, to_0_y_xxxy_xxxz, to_0_y_xxxy_xxyy, to_0_y_xxxy_xxyz, to_0_y_xxxy_xxzz, to_0_y_xxxy_xyyy, to_0_y_xxxy_xyyz, to_0_y_xxxy_xyzz, to_0_y_xxxy_xzzz, to_0_y_xxxy_yyyy, to_0_y_xxxy_yyyz, to_0_y_xxxy_yyzz, to_0_y_xxxy_yzzz, to_0_y_xxxy_zzzz, to_xxxy_xxx, to_xxxy_xxxxy, to_xxxy_xxxyy, to_xxxy_xxxyz, to_xxxy_xxy, to_xxxy_xxyyy, to_xxxy_xxyyz, to_xxxy_xxyzz, to_xxxy_xxz, to_xxxy_xyy, to_xxxy_xyyyy, to_xxxy_xyyyz, to_xxxy_xyyzz, to_xxxy_xyz, to_xxxy_xyzzz, to_xxxy_xzz, to_xxxy_yyy, to_xxxy_yyyyy, to_xxxy_yyyyz, to_xxxy_yyyzz, to_xxxy_yyz, to_xxxy_yyzzz, to_xxxy_yzz, to_xxxy_yzzzz, to_xxxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxy_xxxx[k] = 2.0 * to_xxxy_xxxxy[k] * tke_0;

            to_0_y_xxxy_xxxy[k] = -to_xxxy_xxx[k] + 2.0 * to_xxxy_xxxyy[k] * tke_0;

            to_0_y_xxxy_xxxz[k] = 2.0 * to_xxxy_xxxyz[k] * tke_0;

            to_0_y_xxxy_xxyy[k] = -2.0 * to_xxxy_xxy[k] + 2.0 * to_xxxy_xxyyy[k] * tke_0;

            to_0_y_xxxy_xxyz[k] = -to_xxxy_xxz[k] + 2.0 * to_xxxy_xxyyz[k] * tke_0;

            to_0_y_xxxy_xxzz[k] = 2.0 * to_xxxy_xxyzz[k] * tke_0;

            to_0_y_xxxy_xyyy[k] = -3.0 * to_xxxy_xyy[k] + 2.0 * to_xxxy_xyyyy[k] * tke_0;

            to_0_y_xxxy_xyyz[k] = -2.0 * to_xxxy_xyz[k] + 2.0 * to_xxxy_xyyyz[k] * tke_0;

            to_0_y_xxxy_xyzz[k] = -to_xxxy_xzz[k] + 2.0 * to_xxxy_xyyzz[k] * tke_0;

            to_0_y_xxxy_xzzz[k] = 2.0 * to_xxxy_xyzzz[k] * tke_0;

            to_0_y_xxxy_yyyy[k] = -4.0 * to_xxxy_yyy[k] + 2.0 * to_xxxy_yyyyy[k] * tke_0;

            to_0_y_xxxy_yyyz[k] = -3.0 * to_xxxy_yyz[k] + 2.0 * to_xxxy_yyyyz[k] * tke_0;

            to_0_y_xxxy_yyzz[k] = -2.0 * to_xxxy_yzz[k] + 2.0 * to_xxxy_yyyzz[k] * tke_0;

            to_0_y_xxxy_yzzz[k] = -to_xxxy_zzz[k] + 2.0 * to_xxxy_yyzzz[k] * tke_0;

            to_0_y_xxxy_zzzz[k] = 2.0 * to_xxxy_yzzzz[k] * tke_0;
        }

        // Set up 255-270 components of targeted buffer : GG

        auto to_0_y_xxxz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 30);

        auto to_0_y_xxxz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 31);

        auto to_0_y_xxxz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 32);

        auto to_0_y_xxxz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 33);

        auto to_0_y_xxxz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 34);

        auto to_0_y_xxxz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 35);

        auto to_0_y_xxxz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 36);

        auto to_0_y_xxxz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 37);

        auto to_0_y_xxxz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 38);

        auto to_0_y_xxxz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 39);

        auto to_0_y_xxxz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 40);

        auto to_0_y_xxxz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 41);

        auto to_0_y_xxxz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 42);

        auto to_0_y_xxxz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 43);

        auto to_0_y_xxxz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_0_y_xxxz_xxxx, to_0_y_xxxz_xxxy, to_0_y_xxxz_xxxz, to_0_y_xxxz_xxyy, to_0_y_xxxz_xxyz, to_0_y_xxxz_xxzz, to_0_y_xxxz_xyyy, to_0_y_xxxz_xyyz, to_0_y_xxxz_xyzz, to_0_y_xxxz_xzzz, to_0_y_xxxz_yyyy, to_0_y_xxxz_yyyz, to_0_y_xxxz_yyzz, to_0_y_xxxz_yzzz, to_0_y_xxxz_zzzz, to_xxxz_xxx, to_xxxz_xxxxy, to_xxxz_xxxyy, to_xxxz_xxxyz, to_xxxz_xxy, to_xxxz_xxyyy, to_xxxz_xxyyz, to_xxxz_xxyzz, to_xxxz_xxz, to_xxxz_xyy, to_xxxz_xyyyy, to_xxxz_xyyyz, to_xxxz_xyyzz, to_xxxz_xyz, to_xxxz_xyzzz, to_xxxz_xzz, to_xxxz_yyy, to_xxxz_yyyyy, to_xxxz_yyyyz, to_xxxz_yyyzz, to_xxxz_yyz, to_xxxz_yyzzz, to_xxxz_yzz, to_xxxz_yzzzz, to_xxxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxz_xxxx[k] = 2.0 * to_xxxz_xxxxy[k] * tke_0;

            to_0_y_xxxz_xxxy[k] = -to_xxxz_xxx[k] + 2.0 * to_xxxz_xxxyy[k] * tke_0;

            to_0_y_xxxz_xxxz[k] = 2.0 * to_xxxz_xxxyz[k] * tke_0;

            to_0_y_xxxz_xxyy[k] = -2.0 * to_xxxz_xxy[k] + 2.0 * to_xxxz_xxyyy[k] * tke_0;

            to_0_y_xxxz_xxyz[k] = -to_xxxz_xxz[k] + 2.0 * to_xxxz_xxyyz[k] * tke_0;

            to_0_y_xxxz_xxzz[k] = 2.0 * to_xxxz_xxyzz[k] * tke_0;

            to_0_y_xxxz_xyyy[k] = -3.0 * to_xxxz_xyy[k] + 2.0 * to_xxxz_xyyyy[k] * tke_0;

            to_0_y_xxxz_xyyz[k] = -2.0 * to_xxxz_xyz[k] + 2.0 * to_xxxz_xyyyz[k] * tke_0;

            to_0_y_xxxz_xyzz[k] = -to_xxxz_xzz[k] + 2.0 * to_xxxz_xyyzz[k] * tke_0;

            to_0_y_xxxz_xzzz[k] = 2.0 * to_xxxz_xyzzz[k] * tke_0;

            to_0_y_xxxz_yyyy[k] = -4.0 * to_xxxz_yyy[k] + 2.0 * to_xxxz_yyyyy[k] * tke_0;

            to_0_y_xxxz_yyyz[k] = -3.0 * to_xxxz_yyz[k] + 2.0 * to_xxxz_yyyyz[k] * tke_0;

            to_0_y_xxxz_yyzz[k] = -2.0 * to_xxxz_yzz[k] + 2.0 * to_xxxz_yyyzz[k] * tke_0;

            to_0_y_xxxz_yzzz[k] = -to_xxxz_zzz[k] + 2.0 * to_xxxz_yyzzz[k] * tke_0;

            to_0_y_xxxz_zzzz[k] = 2.0 * to_xxxz_yzzzz[k] * tke_0;
        }

        // Set up 270-285 components of targeted buffer : GG

        auto to_0_y_xxyy_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 45);

        auto to_0_y_xxyy_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 46);

        auto to_0_y_xxyy_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 47);

        auto to_0_y_xxyy_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 48);

        auto to_0_y_xxyy_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 49);

        auto to_0_y_xxyy_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 50);

        auto to_0_y_xxyy_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 51);

        auto to_0_y_xxyy_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 52);

        auto to_0_y_xxyy_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 53);

        auto to_0_y_xxyy_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 54);

        auto to_0_y_xxyy_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 55);

        auto to_0_y_xxyy_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 56);

        auto to_0_y_xxyy_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 57);

        auto to_0_y_xxyy_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 58);

        auto to_0_y_xxyy_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_0_y_xxyy_xxxx, to_0_y_xxyy_xxxy, to_0_y_xxyy_xxxz, to_0_y_xxyy_xxyy, to_0_y_xxyy_xxyz, to_0_y_xxyy_xxzz, to_0_y_xxyy_xyyy, to_0_y_xxyy_xyyz, to_0_y_xxyy_xyzz, to_0_y_xxyy_xzzz, to_0_y_xxyy_yyyy, to_0_y_xxyy_yyyz, to_0_y_xxyy_yyzz, to_0_y_xxyy_yzzz, to_0_y_xxyy_zzzz, to_xxyy_xxx, to_xxyy_xxxxy, to_xxyy_xxxyy, to_xxyy_xxxyz, to_xxyy_xxy, to_xxyy_xxyyy, to_xxyy_xxyyz, to_xxyy_xxyzz, to_xxyy_xxz, to_xxyy_xyy, to_xxyy_xyyyy, to_xxyy_xyyyz, to_xxyy_xyyzz, to_xxyy_xyz, to_xxyy_xyzzz, to_xxyy_xzz, to_xxyy_yyy, to_xxyy_yyyyy, to_xxyy_yyyyz, to_xxyy_yyyzz, to_xxyy_yyz, to_xxyy_yyzzz, to_xxyy_yzz, to_xxyy_yzzzz, to_xxyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxyy_xxxx[k] = 2.0 * to_xxyy_xxxxy[k] * tke_0;

            to_0_y_xxyy_xxxy[k] = -to_xxyy_xxx[k] + 2.0 * to_xxyy_xxxyy[k] * tke_0;

            to_0_y_xxyy_xxxz[k] = 2.0 * to_xxyy_xxxyz[k] * tke_0;

            to_0_y_xxyy_xxyy[k] = -2.0 * to_xxyy_xxy[k] + 2.0 * to_xxyy_xxyyy[k] * tke_0;

            to_0_y_xxyy_xxyz[k] = -to_xxyy_xxz[k] + 2.0 * to_xxyy_xxyyz[k] * tke_0;

            to_0_y_xxyy_xxzz[k] = 2.0 * to_xxyy_xxyzz[k] * tke_0;

            to_0_y_xxyy_xyyy[k] = -3.0 * to_xxyy_xyy[k] + 2.0 * to_xxyy_xyyyy[k] * tke_0;

            to_0_y_xxyy_xyyz[k] = -2.0 * to_xxyy_xyz[k] + 2.0 * to_xxyy_xyyyz[k] * tke_0;

            to_0_y_xxyy_xyzz[k] = -to_xxyy_xzz[k] + 2.0 * to_xxyy_xyyzz[k] * tke_0;

            to_0_y_xxyy_xzzz[k] = 2.0 * to_xxyy_xyzzz[k] * tke_0;

            to_0_y_xxyy_yyyy[k] = -4.0 * to_xxyy_yyy[k] + 2.0 * to_xxyy_yyyyy[k] * tke_0;

            to_0_y_xxyy_yyyz[k] = -3.0 * to_xxyy_yyz[k] + 2.0 * to_xxyy_yyyyz[k] * tke_0;

            to_0_y_xxyy_yyzz[k] = -2.0 * to_xxyy_yzz[k] + 2.0 * to_xxyy_yyyzz[k] * tke_0;

            to_0_y_xxyy_yzzz[k] = -to_xxyy_zzz[k] + 2.0 * to_xxyy_yyzzz[k] * tke_0;

            to_0_y_xxyy_zzzz[k] = 2.0 * to_xxyy_yzzzz[k] * tke_0;
        }

        // Set up 285-300 components of targeted buffer : GG

        auto to_0_y_xxyz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 60);

        auto to_0_y_xxyz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 61);

        auto to_0_y_xxyz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 62);

        auto to_0_y_xxyz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 63);

        auto to_0_y_xxyz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 64);

        auto to_0_y_xxyz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 65);

        auto to_0_y_xxyz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 66);

        auto to_0_y_xxyz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 67);

        auto to_0_y_xxyz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 68);

        auto to_0_y_xxyz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 69);

        auto to_0_y_xxyz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 70);

        auto to_0_y_xxyz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 71);

        auto to_0_y_xxyz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 72);

        auto to_0_y_xxyz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 73);

        auto to_0_y_xxyz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_0_y_xxyz_xxxx, to_0_y_xxyz_xxxy, to_0_y_xxyz_xxxz, to_0_y_xxyz_xxyy, to_0_y_xxyz_xxyz, to_0_y_xxyz_xxzz, to_0_y_xxyz_xyyy, to_0_y_xxyz_xyyz, to_0_y_xxyz_xyzz, to_0_y_xxyz_xzzz, to_0_y_xxyz_yyyy, to_0_y_xxyz_yyyz, to_0_y_xxyz_yyzz, to_0_y_xxyz_yzzz, to_0_y_xxyz_zzzz, to_xxyz_xxx, to_xxyz_xxxxy, to_xxyz_xxxyy, to_xxyz_xxxyz, to_xxyz_xxy, to_xxyz_xxyyy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xyy, to_xxyz_xyyyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_yyy, to_xxyz_yyyyy, to_xxyz_yyyyz, to_xxyz_yyyzz, to_xxyz_yyz, to_xxyz_yyzzz, to_xxyz_yzz, to_xxyz_yzzzz, to_xxyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxyz_xxxx[k] = 2.0 * to_xxyz_xxxxy[k] * tke_0;

            to_0_y_xxyz_xxxy[k] = -to_xxyz_xxx[k] + 2.0 * to_xxyz_xxxyy[k] * tke_0;

            to_0_y_xxyz_xxxz[k] = 2.0 * to_xxyz_xxxyz[k] * tke_0;

            to_0_y_xxyz_xxyy[k] = -2.0 * to_xxyz_xxy[k] + 2.0 * to_xxyz_xxyyy[k] * tke_0;

            to_0_y_xxyz_xxyz[k] = -to_xxyz_xxz[k] + 2.0 * to_xxyz_xxyyz[k] * tke_0;

            to_0_y_xxyz_xxzz[k] = 2.0 * to_xxyz_xxyzz[k] * tke_0;

            to_0_y_xxyz_xyyy[k] = -3.0 * to_xxyz_xyy[k] + 2.0 * to_xxyz_xyyyy[k] * tke_0;

            to_0_y_xxyz_xyyz[k] = -2.0 * to_xxyz_xyz[k] + 2.0 * to_xxyz_xyyyz[k] * tke_0;

            to_0_y_xxyz_xyzz[k] = -to_xxyz_xzz[k] + 2.0 * to_xxyz_xyyzz[k] * tke_0;

            to_0_y_xxyz_xzzz[k] = 2.0 * to_xxyz_xyzzz[k] * tke_0;

            to_0_y_xxyz_yyyy[k] = -4.0 * to_xxyz_yyy[k] + 2.0 * to_xxyz_yyyyy[k] * tke_0;

            to_0_y_xxyz_yyyz[k] = -3.0 * to_xxyz_yyz[k] + 2.0 * to_xxyz_yyyyz[k] * tke_0;

            to_0_y_xxyz_yyzz[k] = -2.0 * to_xxyz_yzz[k] + 2.0 * to_xxyz_yyyzz[k] * tke_0;

            to_0_y_xxyz_yzzz[k] = -to_xxyz_zzz[k] + 2.0 * to_xxyz_yyzzz[k] * tke_0;

            to_0_y_xxyz_zzzz[k] = 2.0 * to_xxyz_yzzzz[k] * tke_0;
        }

        // Set up 300-315 components of targeted buffer : GG

        auto to_0_y_xxzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 75);

        auto to_0_y_xxzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 76);

        auto to_0_y_xxzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 77);

        auto to_0_y_xxzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 78);

        auto to_0_y_xxzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 79);

        auto to_0_y_xxzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 80);

        auto to_0_y_xxzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 81);

        auto to_0_y_xxzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 82);

        auto to_0_y_xxzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 83);

        auto to_0_y_xxzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 84);

        auto to_0_y_xxzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 85);

        auto to_0_y_xxzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 86);

        auto to_0_y_xxzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 87);

        auto to_0_y_xxzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 88);

        auto to_0_y_xxzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_0_y_xxzz_xxxx, to_0_y_xxzz_xxxy, to_0_y_xxzz_xxxz, to_0_y_xxzz_xxyy, to_0_y_xxzz_xxyz, to_0_y_xxzz_xxzz, to_0_y_xxzz_xyyy, to_0_y_xxzz_xyyz, to_0_y_xxzz_xyzz, to_0_y_xxzz_xzzz, to_0_y_xxzz_yyyy, to_0_y_xxzz_yyyz, to_0_y_xxzz_yyzz, to_0_y_xxzz_yzzz, to_0_y_xxzz_zzzz, to_xxzz_xxx, to_xxzz_xxxxy, to_xxzz_xxxyy, to_xxzz_xxxyz, to_xxzz_xxy, to_xxzz_xxyyy, to_xxzz_xxyyz, to_xxzz_xxyzz, to_xxzz_xxz, to_xxzz_xyy, to_xxzz_xyyyy, to_xxzz_xyyyz, to_xxzz_xyyzz, to_xxzz_xyz, to_xxzz_xyzzz, to_xxzz_xzz, to_xxzz_yyy, to_xxzz_yyyyy, to_xxzz_yyyyz, to_xxzz_yyyzz, to_xxzz_yyz, to_xxzz_yyzzz, to_xxzz_yzz, to_xxzz_yzzzz, to_xxzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxzz_xxxx[k] = 2.0 * to_xxzz_xxxxy[k] * tke_0;

            to_0_y_xxzz_xxxy[k] = -to_xxzz_xxx[k] + 2.0 * to_xxzz_xxxyy[k] * tke_0;

            to_0_y_xxzz_xxxz[k] = 2.0 * to_xxzz_xxxyz[k] * tke_0;

            to_0_y_xxzz_xxyy[k] = -2.0 * to_xxzz_xxy[k] + 2.0 * to_xxzz_xxyyy[k] * tke_0;

            to_0_y_xxzz_xxyz[k] = -to_xxzz_xxz[k] + 2.0 * to_xxzz_xxyyz[k] * tke_0;

            to_0_y_xxzz_xxzz[k] = 2.0 * to_xxzz_xxyzz[k] * tke_0;

            to_0_y_xxzz_xyyy[k] = -3.0 * to_xxzz_xyy[k] + 2.0 * to_xxzz_xyyyy[k] * tke_0;

            to_0_y_xxzz_xyyz[k] = -2.0 * to_xxzz_xyz[k] + 2.0 * to_xxzz_xyyyz[k] * tke_0;

            to_0_y_xxzz_xyzz[k] = -to_xxzz_xzz[k] + 2.0 * to_xxzz_xyyzz[k] * tke_0;

            to_0_y_xxzz_xzzz[k] = 2.0 * to_xxzz_xyzzz[k] * tke_0;

            to_0_y_xxzz_yyyy[k] = -4.0 * to_xxzz_yyy[k] + 2.0 * to_xxzz_yyyyy[k] * tke_0;

            to_0_y_xxzz_yyyz[k] = -3.0 * to_xxzz_yyz[k] + 2.0 * to_xxzz_yyyyz[k] * tke_0;

            to_0_y_xxzz_yyzz[k] = -2.0 * to_xxzz_yzz[k] + 2.0 * to_xxzz_yyyzz[k] * tke_0;

            to_0_y_xxzz_yzzz[k] = -to_xxzz_zzz[k] + 2.0 * to_xxzz_yyzzz[k] * tke_0;

            to_0_y_xxzz_zzzz[k] = 2.0 * to_xxzz_yzzzz[k] * tke_0;
        }

        // Set up 315-330 components of targeted buffer : GG

        auto to_0_y_xyyy_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 90);

        auto to_0_y_xyyy_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 91);

        auto to_0_y_xyyy_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 92);

        auto to_0_y_xyyy_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 93);

        auto to_0_y_xyyy_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 94);

        auto to_0_y_xyyy_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 95);

        auto to_0_y_xyyy_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 96);

        auto to_0_y_xyyy_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 97);

        auto to_0_y_xyyy_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 98);

        auto to_0_y_xyyy_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 99);

        auto to_0_y_xyyy_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 100);

        auto to_0_y_xyyy_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 101);

        auto to_0_y_xyyy_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 102);

        auto to_0_y_xyyy_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 103);

        auto to_0_y_xyyy_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_0_y_xyyy_xxxx, to_0_y_xyyy_xxxy, to_0_y_xyyy_xxxz, to_0_y_xyyy_xxyy, to_0_y_xyyy_xxyz, to_0_y_xyyy_xxzz, to_0_y_xyyy_xyyy, to_0_y_xyyy_xyyz, to_0_y_xyyy_xyzz, to_0_y_xyyy_xzzz, to_0_y_xyyy_yyyy, to_0_y_xyyy_yyyz, to_0_y_xyyy_yyzz, to_0_y_xyyy_yzzz, to_0_y_xyyy_zzzz, to_xyyy_xxx, to_xyyy_xxxxy, to_xyyy_xxxyy, to_xyyy_xxxyz, to_xyyy_xxy, to_xyyy_xxyyy, to_xyyy_xxyyz, to_xyyy_xxyzz, to_xyyy_xxz, to_xyyy_xyy, to_xyyy_xyyyy, to_xyyy_xyyyz, to_xyyy_xyyzz, to_xyyy_xyz, to_xyyy_xyzzz, to_xyyy_xzz, to_xyyy_yyy, to_xyyy_yyyyy, to_xyyy_yyyyz, to_xyyy_yyyzz, to_xyyy_yyz, to_xyyy_yyzzz, to_xyyy_yzz, to_xyyy_yzzzz, to_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyyy_xxxx[k] = 2.0 * to_xyyy_xxxxy[k] * tke_0;

            to_0_y_xyyy_xxxy[k] = -to_xyyy_xxx[k] + 2.0 * to_xyyy_xxxyy[k] * tke_0;

            to_0_y_xyyy_xxxz[k] = 2.0 * to_xyyy_xxxyz[k] * tke_0;

            to_0_y_xyyy_xxyy[k] = -2.0 * to_xyyy_xxy[k] + 2.0 * to_xyyy_xxyyy[k] * tke_0;

            to_0_y_xyyy_xxyz[k] = -to_xyyy_xxz[k] + 2.0 * to_xyyy_xxyyz[k] * tke_0;

            to_0_y_xyyy_xxzz[k] = 2.0 * to_xyyy_xxyzz[k] * tke_0;

            to_0_y_xyyy_xyyy[k] = -3.0 * to_xyyy_xyy[k] + 2.0 * to_xyyy_xyyyy[k] * tke_0;

            to_0_y_xyyy_xyyz[k] = -2.0 * to_xyyy_xyz[k] + 2.0 * to_xyyy_xyyyz[k] * tke_0;

            to_0_y_xyyy_xyzz[k] = -to_xyyy_xzz[k] + 2.0 * to_xyyy_xyyzz[k] * tke_0;

            to_0_y_xyyy_xzzz[k] = 2.0 * to_xyyy_xyzzz[k] * tke_0;

            to_0_y_xyyy_yyyy[k] = -4.0 * to_xyyy_yyy[k] + 2.0 * to_xyyy_yyyyy[k] * tke_0;

            to_0_y_xyyy_yyyz[k] = -3.0 * to_xyyy_yyz[k] + 2.0 * to_xyyy_yyyyz[k] * tke_0;

            to_0_y_xyyy_yyzz[k] = -2.0 * to_xyyy_yzz[k] + 2.0 * to_xyyy_yyyzz[k] * tke_0;

            to_0_y_xyyy_yzzz[k] = -to_xyyy_zzz[k] + 2.0 * to_xyyy_yyzzz[k] * tke_0;

            to_0_y_xyyy_zzzz[k] = 2.0 * to_xyyy_yzzzz[k] * tke_0;
        }

        // Set up 330-345 components of targeted buffer : GG

        auto to_0_y_xyyz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 105);

        auto to_0_y_xyyz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 106);

        auto to_0_y_xyyz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 107);

        auto to_0_y_xyyz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 108);

        auto to_0_y_xyyz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 109);

        auto to_0_y_xyyz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 110);

        auto to_0_y_xyyz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 111);

        auto to_0_y_xyyz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 112);

        auto to_0_y_xyyz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 113);

        auto to_0_y_xyyz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 114);

        auto to_0_y_xyyz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 115);

        auto to_0_y_xyyz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 116);

        auto to_0_y_xyyz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 117);

        auto to_0_y_xyyz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 118);

        auto to_0_y_xyyz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_0_y_xyyz_xxxx, to_0_y_xyyz_xxxy, to_0_y_xyyz_xxxz, to_0_y_xyyz_xxyy, to_0_y_xyyz_xxyz, to_0_y_xyyz_xxzz, to_0_y_xyyz_xyyy, to_0_y_xyyz_xyyz, to_0_y_xyyz_xyzz, to_0_y_xyyz_xzzz, to_0_y_xyyz_yyyy, to_0_y_xyyz_yyyz, to_0_y_xyyz_yyzz, to_0_y_xyyz_yzzz, to_0_y_xyyz_zzzz, to_xyyz_xxx, to_xyyz_xxxxy, to_xyyz_xxxyy, to_xyyz_xxxyz, to_xyyz_xxy, to_xyyz_xxyyy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xyy, to_xyyz_xyyyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_yyy, to_xyyz_yyyyy, to_xyyz_yyyyz, to_xyyz_yyyzz, to_xyyz_yyz, to_xyyz_yyzzz, to_xyyz_yzz, to_xyyz_yzzzz, to_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyyz_xxxx[k] = 2.0 * to_xyyz_xxxxy[k] * tke_0;

            to_0_y_xyyz_xxxy[k] = -to_xyyz_xxx[k] + 2.0 * to_xyyz_xxxyy[k] * tke_0;

            to_0_y_xyyz_xxxz[k] = 2.0 * to_xyyz_xxxyz[k] * tke_0;

            to_0_y_xyyz_xxyy[k] = -2.0 * to_xyyz_xxy[k] + 2.0 * to_xyyz_xxyyy[k] * tke_0;

            to_0_y_xyyz_xxyz[k] = -to_xyyz_xxz[k] + 2.0 * to_xyyz_xxyyz[k] * tke_0;

            to_0_y_xyyz_xxzz[k] = 2.0 * to_xyyz_xxyzz[k] * tke_0;

            to_0_y_xyyz_xyyy[k] = -3.0 * to_xyyz_xyy[k] + 2.0 * to_xyyz_xyyyy[k] * tke_0;

            to_0_y_xyyz_xyyz[k] = -2.0 * to_xyyz_xyz[k] + 2.0 * to_xyyz_xyyyz[k] * tke_0;

            to_0_y_xyyz_xyzz[k] = -to_xyyz_xzz[k] + 2.0 * to_xyyz_xyyzz[k] * tke_0;

            to_0_y_xyyz_xzzz[k] = 2.0 * to_xyyz_xyzzz[k] * tke_0;

            to_0_y_xyyz_yyyy[k] = -4.0 * to_xyyz_yyy[k] + 2.0 * to_xyyz_yyyyy[k] * tke_0;

            to_0_y_xyyz_yyyz[k] = -3.0 * to_xyyz_yyz[k] + 2.0 * to_xyyz_yyyyz[k] * tke_0;

            to_0_y_xyyz_yyzz[k] = -2.0 * to_xyyz_yzz[k] + 2.0 * to_xyyz_yyyzz[k] * tke_0;

            to_0_y_xyyz_yzzz[k] = -to_xyyz_zzz[k] + 2.0 * to_xyyz_yyzzz[k] * tke_0;

            to_0_y_xyyz_zzzz[k] = 2.0 * to_xyyz_yzzzz[k] * tke_0;
        }

        // Set up 345-360 components of targeted buffer : GG

        auto to_0_y_xyzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 120);

        auto to_0_y_xyzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 121);

        auto to_0_y_xyzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 122);

        auto to_0_y_xyzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 123);

        auto to_0_y_xyzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 124);

        auto to_0_y_xyzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 125);

        auto to_0_y_xyzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 126);

        auto to_0_y_xyzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 127);

        auto to_0_y_xyzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 128);

        auto to_0_y_xyzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 129);

        auto to_0_y_xyzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 130);

        auto to_0_y_xyzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 131);

        auto to_0_y_xyzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 132);

        auto to_0_y_xyzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 133);

        auto to_0_y_xyzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_0_y_xyzz_xxxx, to_0_y_xyzz_xxxy, to_0_y_xyzz_xxxz, to_0_y_xyzz_xxyy, to_0_y_xyzz_xxyz, to_0_y_xyzz_xxzz, to_0_y_xyzz_xyyy, to_0_y_xyzz_xyyz, to_0_y_xyzz_xyzz, to_0_y_xyzz_xzzz, to_0_y_xyzz_yyyy, to_0_y_xyzz_yyyz, to_0_y_xyzz_yyzz, to_0_y_xyzz_yzzz, to_0_y_xyzz_zzzz, to_xyzz_xxx, to_xyzz_xxxxy, to_xyzz_xxxyy, to_xyzz_xxxyz, to_xyzz_xxy, to_xyzz_xxyyy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xyy, to_xyzz_xyyyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_yyy, to_xyzz_yyyyy, to_xyzz_yyyyz, to_xyzz_yyyzz, to_xyzz_yyz, to_xyzz_yyzzz, to_xyzz_yzz, to_xyzz_yzzzz, to_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyzz_xxxx[k] = 2.0 * to_xyzz_xxxxy[k] * tke_0;

            to_0_y_xyzz_xxxy[k] = -to_xyzz_xxx[k] + 2.0 * to_xyzz_xxxyy[k] * tke_0;

            to_0_y_xyzz_xxxz[k] = 2.0 * to_xyzz_xxxyz[k] * tke_0;

            to_0_y_xyzz_xxyy[k] = -2.0 * to_xyzz_xxy[k] + 2.0 * to_xyzz_xxyyy[k] * tke_0;

            to_0_y_xyzz_xxyz[k] = -to_xyzz_xxz[k] + 2.0 * to_xyzz_xxyyz[k] * tke_0;

            to_0_y_xyzz_xxzz[k] = 2.0 * to_xyzz_xxyzz[k] * tke_0;

            to_0_y_xyzz_xyyy[k] = -3.0 * to_xyzz_xyy[k] + 2.0 * to_xyzz_xyyyy[k] * tke_0;

            to_0_y_xyzz_xyyz[k] = -2.0 * to_xyzz_xyz[k] + 2.0 * to_xyzz_xyyyz[k] * tke_0;

            to_0_y_xyzz_xyzz[k] = -to_xyzz_xzz[k] + 2.0 * to_xyzz_xyyzz[k] * tke_0;

            to_0_y_xyzz_xzzz[k] = 2.0 * to_xyzz_xyzzz[k] * tke_0;

            to_0_y_xyzz_yyyy[k] = -4.0 * to_xyzz_yyy[k] + 2.0 * to_xyzz_yyyyy[k] * tke_0;

            to_0_y_xyzz_yyyz[k] = -3.0 * to_xyzz_yyz[k] + 2.0 * to_xyzz_yyyyz[k] * tke_0;

            to_0_y_xyzz_yyzz[k] = -2.0 * to_xyzz_yzz[k] + 2.0 * to_xyzz_yyyzz[k] * tke_0;

            to_0_y_xyzz_yzzz[k] = -to_xyzz_zzz[k] + 2.0 * to_xyzz_yyzzz[k] * tke_0;

            to_0_y_xyzz_zzzz[k] = 2.0 * to_xyzz_yzzzz[k] * tke_0;
        }

        // Set up 360-375 components of targeted buffer : GG

        auto to_0_y_xzzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 135);

        auto to_0_y_xzzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 136);

        auto to_0_y_xzzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 137);

        auto to_0_y_xzzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 138);

        auto to_0_y_xzzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 139);

        auto to_0_y_xzzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 140);

        auto to_0_y_xzzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 141);

        auto to_0_y_xzzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 142);

        auto to_0_y_xzzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 143);

        auto to_0_y_xzzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 144);

        auto to_0_y_xzzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 145);

        auto to_0_y_xzzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 146);

        auto to_0_y_xzzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 147);

        auto to_0_y_xzzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 148);

        auto to_0_y_xzzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_0_y_xzzz_xxxx, to_0_y_xzzz_xxxy, to_0_y_xzzz_xxxz, to_0_y_xzzz_xxyy, to_0_y_xzzz_xxyz, to_0_y_xzzz_xxzz, to_0_y_xzzz_xyyy, to_0_y_xzzz_xyyz, to_0_y_xzzz_xyzz, to_0_y_xzzz_xzzz, to_0_y_xzzz_yyyy, to_0_y_xzzz_yyyz, to_0_y_xzzz_yyzz, to_0_y_xzzz_yzzz, to_0_y_xzzz_zzzz, to_xzzz_xxx, to_xzzz_xxxxy, to_xzzz_xxxyy, to_xzzz_xxxyz, to_xzzz_xxy, to_xzzz_xxyyy, to_xzzz_xxyyz, to_xzzz_xxyzz, to_xzzz_xxz, to_xzzz_xyy, to_xzzz_xyyyy, to_xzzz_xyyyz, to_xzzz_xyyzz, to_xzzz_xyz, to_xzzz_xyzzz, to_xzzz_xzz, to_xzzz_yyy, to_xzzz_yyyyy, to_xzzz_yyyyz, to_xzzz_yyyzz, to_xzzz_yyz, to_xzzz_yyzzz, to_xzzz_yzz, to_xzzz_yzzzz, to_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xzzz_xxxx[k] = 2.0 * to_xzzz_xxxxy[k] * tke_0;

            to_0_y_xzzz_xxxy[k] = -to_xzzz_xxx[k] + 2.0 * to_xzzz_xxxyy[k] * tke_0;

            to_0_y_xzzz_xxxz[k] = 2.0 * to_xzzz_xxxyz[k] * tke_0;

            to_0_y_xzzz_xxyy[k] = -2.0 * to_xzzz_xxy[k] + 2.0 * to_xzzz_xxyyy[k] * tke_0;

            to_0_y_xzzz_xxyz[k] = -to_xzzz_xxz[k] + 2.0 * to_xzzz_xxyyz[k] * tke_0;

            to_0_y_xzzz_xxzz[k] = 2.0 * to_xzzz_xxyzz[k] * tke_0;

            to_0_y_xzzz_xyyy[k] = -3.0 * to_xzzz_xyy[k] + 2.0 * to_xzzz_xyyyy[k] * tke_0;

            to_0_y_xzzz_xyyz[k] = -2.0 * to_xzzz_xyz[k] + 2.0 * to_xzzz_xyyyz[k] * tke_0;

            to_0_y_xzzz_xyzz[k] = -to_xzzz_xzz[k] + 2.0 * to_xzzz_xyyzz[k] * tke_0;

            to_0_y_xzzz_xzzz[k] = 2.0 * to_xzzz_xyzzz[k] * tke_0;

            to_0_y_xzzz_yyyy[k] = -4.0 * to_xzzz_yyy[k] + 2.0 * to_xzzz_yyyyy[k] * tke_0;

            to_0_y_xzzz_yyyz[k] = -3.0 * to_xzzz_yyz[k] + 2.0 * to_xzzz_yyyyz[k] * tke_0;

            to_0_y_xzzz_yyzz[k] = -2.0 * to_xzzz_yzz[k] + 2.0 * to_xzzz_yyyzz[k] * tke_0;

            to_0_y_xzzz_yzzz[k] = -to_xzzz_zzz[k] + 2.0 * to_xzzz_yyzzz[k] * tke_0;

            to_0_y_xzzz_zzzz[k] = 2.0 * to_xzzz_yzzzz[k] * tke_0;
        }

        // Set up 375-390 components of targeted buffer : GG

        auto to_0_y_yyyy_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 150);

        auto to_0_y_yyyy_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 151);

        auto to_0_y_yyyy_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 152);

        auto to_0_y_yyyy_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 153);

        auto to_0_y_yyyy_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 154);

        auto to_0_y_yyyy_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 155);

        auto to_0_y_yyyy_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 156);

        auto to_0_y_yyyy_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 157);

        auto to_0_y_yyyy_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 158);

        auto to_0_y_yyyy_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 159);

        auto to_0_y_yyyy_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 160);

        auto to_0_y_yyyy_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 161);

        auto to_0_y_yyyy_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 162);

        auto to_0_y_yyyy_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 163);

        auto to_0_y_yyyy_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_0_y_yyyy_xxxx, to_0_y_yyyy_xxxy, to_0_y_yyyy_xxxz, to_0_y_yyyy_xxyy, to_0_y_yyyy_xxyz, to_0_y_yyyy_xxzz, to_0_y_yyyy_xyyy, to_0_y_yyyy_xyyz, to_0_y_yyyy_xyzz, to_0_y_yyyy_xzzz, to_0_y_yyyy_yyyy, to_0_y_yyyy_yyyz, to_0_y_yyyy_yyzz, to_0_y_yyyy_yzzz, to_0_y_yyyy_zzzz, to_yyyy_xxx, to_yyyy_xxxxy, to_yyyy_xxxyy, to_yyyy_xxxyz, to_yyyy_xxy, to_yyyy_xxyyy, to_yyyy_xxyyz, to_yyyy_xxyzz, to_yyyy_xxz, to_yyyy_xyy, to_yyyy_xyyyy, to_yyyy_xyyyz, to_yyyy_xyyzz, to_yyyy_xyz, to_yyyy_xyzzz, to_yyyy_xzz, to_yyyy_yyy, to_yyyy_yyyyy, to_yyyy_yyyyz, to_yyyy_yyyzz, to_yyyy_yyz, to_yyyy_yyzzz, to_yyyy_yzz, to_yyyy_yzzzz, to_yyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyyy_xxxx[k] = 2.0 * to_yyyy_xxxxy[k] * tke_0;

            to_0_y_yyyy_xxxy[k] = -to_yyyy_xxx[k] + 2.0 * to_yyyy_xxxyy[k] * tke_0;

            to_0_y_yyyy_xxxz[k] = 2.0 * to_yyyy_xxxyz[k] * tke_0;

            to_0_y_yyyy_xxyy[k] = -2.0 * to_yyyy_xxy[k] + 2.0 * to_yyyy_xxyyy[k] * tke_0;

            to_0_y_yyyy_xxyz[k] = -to_yyyy_xxz[k] + 2.0 * to_yyyy_xxyyz[k] * tke_0;

            to_0_y_yyyy_xxzz[k] = 2.0 * to_yyyy_xxyzz[k] * tke_0;

            to_0_y_yyyy_xyyy[k] = -3.0 * to_yyyy_xyy[k] + 2.0 * to_yyyy_xyyyy[k] * tke_0;

            to_0_y_yyyy_xyyz[k] = -2.0 * to_yyyy_xyz[k] + 2.0 * to_yyyy_xyyyz[k] * tke_0;

            to_0_y_yyyy_xyzz[k] = -to_yyyy_xzz[k] + 2.0 * to_yyyy_xyyzz[k] * tke_0;

            to_0_y_yyyy_xzzz[k] = 2.0 * to_yyyy_xyzzz[k] * tke_0;

            to_0_y_yyyy_yyyy[k] = -4.0 * to_yyyy_yyy[k] + 2.0 * to_yyyy_yyyyy[k] * tke_0;

            to_0_y_yyyy_yyyz[k] = -3.0 * to_yyyy_yyz[k] + 2.0 * to_yyyy_yyyyz[k] * tke_0;

            to_0_y_yyyy_yyzz[k] = -2.0 * to_yyyy_yzz[k] + 2.0 * to_yyyy_yyyzz[k] * tke_0;

            to_0_y_yyyy_yzzz[k] = -to_yyyy_zzz[k] + 2.0 * to_yyyy_yyzzz[k] * tke_0;

            to_0_y_yyyy_zzzz[k] = 2.0 * to_yyyy_yzzzz[k] * tke_0;
        }

        // Set up 390-405 components of targeted buffer : GG

        auto to_0_y_yyyz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 165);

        auto to_0_y_yyyz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 166);

        auto to_0_y_yyyz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 167);

        auto to_0_y_yyyz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 168);

        auto to_0_y_yyyz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 169);

        auto to_0_y_yyyz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 170);

        auto to_0_y_yyyz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 171);

        auto to_0_y_yyyz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 172);

        auto to_0_y_yyyz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 173);

        auto to_0_y_yyyz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 174);

        auto to_0_y_yyyz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 175);

        auto to_0_y_yyyz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 176);

        auto to_0_y_yyyz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 177);

        auto to_0_y_yyyz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 178);

        auto to_0_y_yyyz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_0_y_yyyz_xxxx, to_0_y_yyyz_xxxy, to_0_y_yyyz_xxxz, to_0_y_yyyz_xxyy, to_0_y_yyyz_xxyz, to_0_y_yyyz_xxzz, to_0_y_yyyz_xyyy, to_0_y_yyyz_xyyz, to_0_y_yyyz_xyzz, to_0_y_yyyz_xzzz, to_0_y_yyyz_yyyy, to_0_y_yyyz_yyyz, to_0_y_yyyz_yyzz, to_0_y_yyyz_yzzz, to_0_y_yyyz_zzzz, to_yyyz_xxx, to_yyyz_xxxxy, to_yyyz_xxxyy, to_yyyz_xxxyz, to_yyyz_xxy, to_yyyz_xxyyy, to_yyyz_xxyyz, to_yyyz_xxyzz, to_yyyz_xxz, to_yyyz_xyy, to_yyyz_xyyyy, to_yyyz_xyyyz, to_yyyz_xyyzz, to_yyyz_xyz, to_yyyz_xyzzz, to_yyyz_xzz, to_yyyz_yyy, to_yyyz_yyyyy, to_yyyz_yyyyz, to_yyyz_yyyzz, to_yyyz_yyz, to_yyyz_yyzzz, to_yyyz_yzz, to_yyyz_yzzzz, to_yyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyyz_xxxx[k] = 2.0 * to_yyyz_xxxxy[k] * tke_0;

            to_0_y_yyyz_xxxy[k] = -to_yyyz_xxx[k] + 2.0 * to_yyyz_xxxyy[k] * tke_0;

            to_0_y_yyyz_xxxz[k] = 2.0 * to_yyyz_xxxyz[k] * tke_0;

            to_0_y_yyyz_xxyy[k] = -2.0 * to_yyyz_xxy[k] + 2.0 * to_yyyz_xxyyy[k] * tke_0;

            to_0_y_yyyz_xxyz[k] = -to_yyyz_xxz[k] + 2.0 * to_yyyz_xxyyz[k] * tke_0;

            to_0_y_yyyz_xxzz[k] = 2.0 * to_yyyz_xxyzz[k] * tke_0;

            to_0_y_yyyz_xyyy[k] = -3.0 * to_yyyz_xyy[k] + 2.0 * to_yyyz_xyyyy[k] * tke_0;

            to_0_y_yyyz_xyyz[k] = -2.0 * to_yyyz_xyz[k] + 2.0 * to_yyyz_xyyyz[k] * tke_0;

            to_0_y_yyyz_xyzz[k] = -to_yyyz_xzz[k] + 2.0 * to_yyyz_xyyzz[k] * tke_0;

            to_0_y_yyyz_xzzz[k] = 2.0 * to_yyyz_xyzzz[k] * tke_0;

            to_0_y_yyyz_yyyy[k] = -4.0 * to_yyyz_yyy[k] + 2.0 * to_yyyz_yyyyy[k] * tke_0;

            to_0_y_yyyz_yyyz[k] = -3.0 * to_yyyz_yyz[k] + 2.0 * to_yyyz_yyyyz[k] * tke_0;

            to_0_y_yyyz_yyzz[k] = -2.0 * to_yyyz_yzz[k] + 2.0 * to_yyyz_yyyzz[k] * tke_0;

            to_0_y_yyyz_yzzz[k] = -to_yyyz_zzz[k] + 2.0 * to_yyyz_yyzzz[k] * tke_0;

            to_0_y_yyyz_zzzz[k] = 2.0 * to_yyyz_yzzzz[k] * tke_0;
        }

        // Set up 405-420 components of targeted buffer : GG

        auto to_0_y_yyzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 180);

        auto to_0_y_yyzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 181);

        auto to_0_y_yyzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 182);

        auto to_0_y_yyzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 183);

        auto to_0_y_yyzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 184);

        auto to_0_y_yyzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 185);

        auto to_0_y_yyzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 186);

        auto to_0_y_yyzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 187);

        auto to_0_y_yyzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 188);

        auto to_0_y_yyzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 189);

        auto to_0_y_yyzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 190);

        auto to_0_y_yyzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 191);

        auto to_0_y_yyzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 192);

        auto to_0_y_yyzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 193);

        auto to_0_y_yyzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_0_y_yyzz_xxxx, to_0_y_yyzz_xxxy, to_0_y_yyzz_xxxz, to_0_y_yyzz_xxyy, to_0_y_yyzz_xxyz, to_0_y_yyzz_xxzz, to_0_y_yyzz_xyyy, to_0_y_yyzz_xyyz, to_0_y_yyzz_xyzz, to_0_y_yyzz_xzzz, to_0_y_yyzz_yyyy, to_0_y_yyzz_yyyz, to_0_y_yyzz_yyzz, to_0_y_yyzz_yzzz, to_0_y_yyzz_zzzz, to_yyzz_xxx, to_yyzz_xxxxy, to_yyzz_xxxyy, to_yyzz_xxxyz, to_yyzz_xxy, to_yyzz_xxyyy, to_yyzz_xxyyz, to_yyzz_xxyzz, to_yyzz_xxz, to_yyzz_xyy, to_yyzz_xyyyy, to_yyzz_xyyyz, to_yyzz_xyyzz, to_yyzz_xyz, to_yyzz_xyzzz, to_yyzz_xzz, to_yyzz_yyy, to_yyzz_yyyyy, to_yyzz_yyyyz, to_yyzz_yyyzz, to_yyzz_yyz, to_yyzz_yyzzz, to_yyzz_yzz, to_yyzz_yzzzz, to_yyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyzz_xxxx[k] = 2.0 * to_yyzz_xxxxy[k] * tke_0;

            to_0_y_yyzz_xxxy[k] = -to_yyzz_xxx[k] + 2.0 * to_yyzz_xxxyy[k] * tke_0;

            to_0_y_yyzz_xxxz[k] = 2.0 * to_yyzz_xxxyz[k] * tke_0;

            to_0_y_yyzz_xxyy[k] = -2.0 * to_yyzz_xxy[k] + 2.0 * to_yyzz_xxyyy[k] * tke_0;

            to_0_y_yyzz_xxyz[k] = -to_yyzz_xxz[k] + 2.0 * to_yyzz_xxyyz[k] * tke_0;

            to_0_y_yyzz_xxzz[k] = 2.0 * to_yyzz_xxyzz[k] * tke_0;

            to_0_y_yyzz_xyyy[k] = -3.0 * to_yyzz_xyy[k] + 2.0 * to_yyzz_xyyyy[k] * tke_0;

            to_0_y_yyzz_xyyz[k] = -2.0 * to_yyzz_xyz[k] + 2.0 * to_yyzz_xyyyz[k] * tke_0;

            to_0_y_yyzz_xyzz[k] = -to_yyzz_xzz[k] + 2.0 * to_yyzz_xyyzz[k] * tke_0;

            to_0_y_yyzz_xzzz[k] = 2.0 * to_yyzz_xyzzz[k] * tke_0;

            to_0_y_yyzz_yyyy[k] = -4.0 * to_yyzz_yyy[k] + 2.0 * to_yyzz_yyyyy[k] * tke_0;

            to_0_y_yyzz_yyyz[k] = -3.0 * to_yyzz_yyz[k] + 2.0 * to_yyzz_yyyyz[k] * tke_0;

            to_0_y_yyzz_yyzz[k] = -2.0 * to_yyzz_yzz[k] + 2.0 * to_yyzz_yyyzz[k] * tke_0;

            to_0_y_yyzz_yzzz[k] = -to_yyzz_zzz[k] + 2.0 * to_yyzz_yyzzz[k] * tke_0;

            to_0_y_yyzz_zzzz[k] = 2.0 * to_yyzz_yzzzz[k] * tke_0;
        }

        // Set up 420-435 components of targeted buffer : GG

        auto to_0_y_yzzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 195);

        auto to_0_y_yzzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 196);

        auto to_0_y_yzzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 197);

        auto to_0_y_yzzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 198);

        auto to_0_y_yzzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 199);

        auto to_0_y_yzzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 200);

        auto to_0_y_yzzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 201);

        auto to_0_y_yzzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 202);

        auto to_0_y_yzzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 203);

        auto to_0_y_yzzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 204);

        auto to_0_y_yzzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 205);

        auto to_0_y_yzzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 206);

        auto to_0_y_yzzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 207);

        auto to_0_y_yzzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 208);

        auto to_0_y_yzzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_0_y_yzzz_xxxx, to_0_y_yzzz_xxxy, to_0_y_yzzz_xxxz, to_0_y_yzzz_xxyy, to_0_y_yzzz_xxyz, to_0_y_yzzz_xxzz, to_0_y_yzzz_xyyy, to_0_y_yzzz_xyyz, to_0_y_yzzz_xyzz, to_0_y_yzzz_xzzz, to_0_y_yzzz_yyyy, to_0_y_yzzz_yyyz, to_0_y_yzzz_yyzz, to_0_y_yzzz_yzzz, to_0_y_yzzz_zzzz, to_yzzz_xxx, to_yzzz_xxxxy, to_yzzz_xxxyy, to_yzzz_xxxyz, to_yzzz_xxy, to_yzzz_xxyyy, to_yzzz_xxyyz, to_yzzz_xxyzz, to_yzzz_xxz, to_yzzz_xyy, to_yzzz_xyyyy, to_yzzz_xyyyz, to_yzzz_xyyzz, to_yzzz_xyz, to_yzzz_xyzzz, to_yzzz_xzz, to_yzzz_yyy, to_yzzz_yyyyy, to_yzzz_yyyyz, to_yzzz_yyyzz, to_yzzz_yyz, to_yzzz_yyzzz, to_yzzz_yzz, to_yzzz_yzzzz, to_yzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yzzz_xxxx[k] = 2.0 * to_yzzz_xxxxy[k] * tke_0;

            to_0_y_yzzz_xxxy[k] = -to_yzzz_xxx[k] + 2.0 * to_yzzz_xxxyy[k] * tke_0;

            to_0_y_yzzz_xxxz[k] = 2.0 * to_yzzz_xxxyz[k] * tke_0;

            to_0_y_yzzz_xxyy[k] = -2.0 * to_yzzz_xxy[k] + 2.0 * to_yzzz_xxyyy[k] * tke_0;

            to_0_y_yzzz_xxyz[k] = -to_yzzz_xxz[k] + 2.0 * to_yzzz_xxyyz[k] * tke_0;

            to_0_y_yzzz_xxzz[k] = 2.0 * to_yzzz_xxyzz[k] * tke_0;

            to_0_y_yzzz_xyyy[k] = -3.0 * to_yzzz_xyy[k] + 2.0 * to_yzzz_xyyyy[k] * tke_0;

            to_0_y_yzzz_xyyz[k] = -2.0 * to_yzzz_xyz[k] + 2.0 * to_yzzz_xyyyz[k] * tke_0;

            to_0_y_yzzz_xyzz[k] = -to_yzzz_xzz[k] + 2.0 * to_yzzz_xyyzz[k] * tke_0;

            to_0_y_yzzz_xzzz[k] = 2.0 * to_yzzz_xyzzz[k] * tke_0;

            to_0_y_yzzz_yyyy[k] = -4.0 * to_yzzz_yyy[k] + 2.0 * to_yzzz_yyyyy[k] * tke_0;

            to_0_y_yzzz_yyyz[k] = -3.0 * to_yzzz_yyz[k] + 2.0 * to_yzzz_yyyyz[k] * tke_0;

            to_0_y_yzzz_yyzz[k] = -2.0 * to_yzzz_yzz[k] + 2.0 * to_yzzz_yyyzz[k] * tke_0;

            to_0_y_yzzz_yzzz[k] = -to_yzzz_zzz[k] + 2.0 * to_yzzz_yyzzz[k] * tke_0;

            to_0_y_yzzz_zzzz[k] = 2.0 * to_yzzz_yzzzz[k] * tke_0;
        }

        // Set up 435-450 components of targeted buffer : GG

        auto to_0_y_zzzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 210);

        auto to_0_y_zzzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 211);

        auto to_0_y_zzzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 212);

        auto to_0_y_zzzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 213);

        auto to_0_y_zzzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 214);

        auto to_0_y_zzzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 215);

        auto to_0_y_zzzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 216);

        auto to_0_y_zzzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 217);

        auto to_0_y_zzzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 218);

        auto to_0_y_zzzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 219);

        auto to_0_y_zzzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 220);

        auto to_0_y_zzzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 221);

        auto to_0_y_zzzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 222);

        auto to_0_y_zzzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 223);

        auto to_0_y_zzzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 1 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_0_y_zzzz_xxxx, to_0_y_zzzz_xxxy, to_0_y_zzzz_xxxz, to_0_y_zzzz_xxyy, to_0_y_zzzz_xxyz, to_0_y_zzzz_xxzz, to_0_y_zzzz_xyyy, to_0_y_zzzz_xyyz, to_0_y_zzzz_xyzz, to_0_y_zzzz_xzzz, to_0_y_zzzz_yyyy, to_0_y_zzzz_yyyz, to_0_y_zzzz_yyzz, to_0_y_zzzz_yzzz, to_0_y_zzzz_zzzz, to_zzzz_xxx, to_zzzz_xxxxy, to_zzzz_xxxyy, to_zzzz_xxxyz, to_zzzz_xxy, to_zzzz_xxyyy, to_zzzz_xxyyz, to_zzzz_xxyzz, to_zzzz_xxz, to_zzzz_xyy, to_zzzz_xyyyy, to_zzzz_xyyyz, to_zzzz_xyyzz, to_zzzz_xyz, to_zzzz_xyzzz, to_zzzz_xzz, to_zzzz_yyy, to_zzzz_yyyyy, to_zzzz_yyyyz, to_zzzz_yyyzz, to_zzzz_yyz, to_zzzz_yyzzz, to_zzzz_yzz, to_zzzz_yzzzz, to_zzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zzzz_xxxx[k] = 2.0 * to_zzzz_xxxxy[k] * tke_0;

            to_0_y_zzzz_xxxy[k] = -to_zzzz_xxx[k] + 2.0 * to_zzzz_xxxyy[k] * tke_0;

            to_0_y_zzzz_xxxz[k] = 2.0 * to_zzzz_xxxyz[k] * tke_0;

            to_0_y_zzzz_xxyy[k] = -2.0 * to_zzzz_xxy[k] + 2.0 * to_zzzz_xxyyy[k] * tke_0;

            to_0_y_zzzz_xxyz[k] = -to_zzzz_xxz[k] + 2.0 * to_zzzz_xxyyz[k] * tke_0;

            to_0_y_zzzz_xxzz[k] = 2.0 * to_zzzz_xxyzz[k] * tke_0;

            to_0_y_zzzz_xyyy[k] = -3.0 * to_zzzz_xyy[k] + 2.0 * to_zzzz_xyyyy[k] * tke_0;

            to_0_y_zzzz_xyyz[k] = -2.0 * to_zzzz_xyz[k] + 2.0 * to_zzzz_xyyyz[k] * tke_0;

            to_0_y_zzzz_xyzz[k] = -to_zzzz_xzz[k] + 2.0 * to_zzzz_xyyzz[k] * tke_0;

            to_0_y_zzzz_xzzz[k] = 2.0 * to_zzzz_xyzzz[k] * tke_0;

            to_0_y_zzzz_yyyy[k] = -4.0 * to_zzzz_yyy[k] + 2.0 * to_zzzz_yyyyy[k] * tke_0;

            to_0_y_zzzz_yyyz[k] = -3.0 * to_zzzz_yyz[k] + 2.0 * to_zzzz_yyyyz[k] * tke_0;

            to_0_y_zzzz_yyzz[k] = -2.0 * to_zzzz_yzz[k] + 2.0 * to_zzzz_yyyzz[k] * tke_0;

            to_0_y_zzzz_yzzz[k] = -to_zzzz_zzz[k] + 2.0 * to_zzzz_yyzzz[k] * tke_0;

            to_0_y_zzzz_zzzz[k] = 2.0 * to_zzzz_yzzzz[k] * tke_0;
        }

        // Set up 450-465 components of targeted buffer : GG

        auto to_0_z_xxxx_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 0);

        auto to_0_z_xxxx_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 1);

        auto to_0_z_xxxx_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 2);

        auto to_0_z_xxxx_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 3);

        auto to_0_z_xxxx_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 4);

        auto to_0_z_xxxx_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 5);

        auto to_0_z_xxxx_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 6);

        auto to_0_z_xxxx_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 7);

        auto to_0_z_xxxx_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 8);

        auto to_0_z_xxxx_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 9);

        auto to_0_z_xxxx_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 10);

        auto to_0_z_xxxx_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 11);

        auto to_0_z_xxxx_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 12);

        auto to_0_z_xxxx_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 13);

        auto to_0_z_xxxx_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 14);

        #pragma omp simd aligned(to_0_z_xxxx_xxxx, to_0_z_xxxx_xxxy, to_0_z_xxxx_xxxz, to_0_z_xxxx_xxyy, to_0_z_xxxx_xxyz, to_0_z_xxxx_xxzz, to_0_z_xxxx_xyyy, to_0_z_xxxx_xyyz, to_0_z_xxxx_xyzz, to_0_z_xxxx_xzzz, to_0_z_xxxx_yyyy, to_0_z_xxxx_yyyz, to_0_z_xxxx_yyzz, to_0_z_xxxx_yzzz, to_0_z_xxxx_zzzz, to_xxxx_xxx, to_xxxx_xxxxz, to_xxxx_xxxyz, to_xxxx_xxxzz, to_xxxx_xxy, to_xxxx_xxyyz, to_xxxx_xxyzz, to_xxxx_xxz, to_xxxx_xxzzz, to_xxxx_xyy, to_xxxx_xyyyz, to_xxxx_xyyzz, to_xxxx_xyz, to_xxxx_xyzzz, to_xxxx_xzz, to_xxxx_xzzzz, to_xxxx_yyy, to_xxxx_yyyyz, to_xxxx_yyyzz, to_xxxx_yyz, to_xxxx_yyzzz, to_xxxx_yzz, to_xxxx_yzzzz, to_xxxx_zzz, to_xxxx_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxx_xxxx[k] = 2.0 * to_xxxx_xxxxz[k] * tke_0;

            to_0_z_xxxx_xxxy[k] = 2.0 * to_xxxx_xxxyz[k] * tke_0;

            to_0_z_xxxx_xxxz[k] = -to_xxxx_xxx[k] + 2.0 * to_xxxx_xxxzz[k] * tke_0;

            to_0_z_xxxx_xxyy[k] = 2.0 * to_xxxx_xxyyz[k] * tke_0;

            to_0_z_xxxx_xxyz[k] = -to_xxxx_xxy[k] + 2.0 * to_xxxx_xxyzz[k] * tke_0;

            to_0_z_xxxx_xxzz[k] = -2.0 * to_xxxx_xxz[k] + 2.0 * to_xxxx_xxzzz[k] * tke_0;

            to_0_z_xxxx_xyyy[k] = 2.0 * to_xxxx_xyyyz[k] * tke_0;

            to_0_z_xxxx_xyyz[k] = -to_xxxx_xyy[k] + 2.0 * to_xxxx_xyyzz[k] * tke_0;

            to_0_z_xxxx_xyzz[k] = -2.0 * to_xxxx_xyz[k] + 2.0 * to_xxxx_xyzzz[k] * tke_0;

            to_0_z_xxxx_xzzz[k] = -3.0 * to_xxxx_xzz[k] + 2.0 * to_xxxx_xzzzz[k] * tke_0;

            to_0_z_xxxx_yyyy[k] = 2.0 * to_xxxx_yyyyz[k] * tke_0;

            to_0_z_xxxx_yyyz[k] = -to_xxxx_yyy[k] + 2.0 * to_xxxx_yyyzz[k] * tke_0;

            to_0_z_xxxx_yyzz[k] = -2.0 * to_xxxx_yyz[k] + 2.0 * to_xxxx_yyzzz[k] * tke_0;

            to_0_z_xxxx_yzzz[k] = -3.0 * to_xxxx_yzz[k] + 2.0 * to_xxxx_yzzzz[k] * tke_0;

            to_0_z_xxxx_zzzz[k] = -4.0 * to_xxxx_zzz[k] + 2.0 * to_xxxx_zzzzz[k] * tke_0;
        }

        // Set up 465-480 components of targeted buffer : GG

        auto to_0_z_xxxy_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 15);

        auto to_0_z_xxxy_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 16);

        auto to_0_z_xxxy_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 17);

        auto to_0_z_xxxy_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 18);

        auto to_0_z_xxxy_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 19);

        auto to_0_z_xxxy_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 20);

        auto to_0_z_xxxy_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 21);

        auto to_0_z_xxxy_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 22);

        auto to_0_z_xxxy_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 23);

        auto to_0_z_xxxy_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 24);

        auto to_0_z_xxxy_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 25);

        auto to_0_z_xxxy_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 26);

        auto to_0_z_xxxy_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 27);

        auto to_0_z_xxxy_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 28);

        auto to_0_z_xxxy_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 29);

        #pragma omp simd aligned(to_0_z_xxxy_xxxx, to_0_z_xxxy_xxxy, to_0_z_xxxy_xxxz, to_0_z_xxxy_xxyy, to_0_z_xxxy_xxyz, to_0_z_xxxy_xxzz, to_0_z_xxxy_xyyy, to_0_z_xxxy_xyyz, to_0_z_xxxy_xyzz, to_0_z_xxxy_xzzz, to_0_z_xxxy_yyyy, to_0_z_xxxy_yyyz, to_0_z_xxxy_yyzz, to_0_z_xxxy_yzzz, to_0_z_xxxy_zzzz, to_xxxy_xxx, to_xxxy_xxxxz, to_xxxy_xxxyz, to_xxxy_xxxzz, to_xxxy_xxy, to_xxxy_xxyyz, to_xxxy_xxyzz, to_xxxy_xxz, to_xxxy_xxzzz, to_xxxy_xyy, to_xxxy_xyyyz, to_xxxy_xyyzz, to_xxxy_xyz, to_xxxy_xyzzz, to_xxxy_xzz, to_xxxy_xzzzz, to_xxxy_yyy, to_xxxy_yyyyz, to_xxxy_yyyzz, to_xxxy_yyz, to_xxxy_yyzzz, to_xxxy_yzz, to_xxxy_yzzzz, to_xxxy_zzz, to_xxxy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxy_xxxx[k] = 2.0 * to_xxxy_xxxxz[k] * tke_0;

            to_0_z_xxxy_xxxy[k] = 2.0 * to_xxxy_xxxyz[k] * tke_0;

            to_0_z_xxxy_xxxz[k] = -to_xxxy_xxx[k] + 2.0 * to_xxxy_xxxzz[k] * tke_0;

            to_0_z_xxxy_xxyy[k] = 2.0 * to_xxxy_xxyyz[k] * tke_0;

            to_0_z_xxxy_xxyz[k] = -to_xxxy_xxy[k] + 2.0 * to_xxxy_xxyzz[k] * tke_0;

            to_0_z_xxxy_xxzz[k] = -2.0 * to_xxxy_xxz[k] + 2.0 * to_xxxy_xxzzz[k] * tke_0;

            to_0_z_xxxy_xyyy[k] = 2.0 * to_xxxy_xyyyz[k] * tke_0;

            to_0_z_xxxy_xyyz[k] = -to_xxxy_xyy[k] + 2.0 * to_xxxy_xyyzz[k] * tke_0;

            to_0_z_xxxy_xyzz[k] = -2.0 * to_xxxy_xyz[k] + 2.0 * to_xxxy_xyzzz[k] * tke_0;

            to_0_z_xxxy_xzzz[k] = -3.0 * to_xxxy_xzz[k] + 2.0 * to_xxxy_xzzzz[k] * tke_0;

            to_0_z_xxxy_yyyy[k] = 2.0 * to_xxxy_yyyyz[k] * tke_0;

            to_0_z_xxxy_yyyz[k] = -to_xxxy_yyy[k] + 2.0 * to_xxxy_yyyzz[k] * tke_0;

            to_0_z_xxxy_yyzz[k] = -2.0 * to_xxxy_yyz[k] + 2.0 * to_xxxy_yyzzz[k] * tke_0;

            to_0_z_xxxy_yzzz[k] = -3.0 * to_xxxy_yzz[k] + 2.0 * to_xxxy_yzzzz[k] * tke_0;

            to_0_z_xxxy_zzzz[k] = -4.0 * to_xxxy_zzz[k] + 2.0 * to_xxxy_zzzzz[k] * tke_0;
        }

        // Set up 480-495 components of targeted buffer : GG

        auto to_0_z_xxxz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 30);

        auto to_0_z_xxxz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 31);

        auto to_0_z_xxxz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 32);

        auto to_0_z_xxxz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 33);

        auto to_0_z_xxxz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 34);

        auto to_0_z_xxxz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 35);

        auto to_0_z_xxxz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 36);

        auto to_0_z_xxxz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 37);

        auto to_0_z_xxxz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 38);

        auto to_0_z_xxxz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 39);

        auto to_0_z_xxxz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 40);

        auto to_0_z_xxxz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 41);

        auto to_0_z_xxxz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 42);

        auto to_0_z_xxxz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 43);

        auto to_0_z_xxxz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 44);

        #pragma omp simd aligned(to_0_z_xxxz_xxxx, to_0_z_xxxz_xxxy, to_0_z_xxxz_xxxz, to_0_z_xxxz_xxyy, to_0_z_xxxz_xxyz, to_0_z_xxxz_xxzz, to_0_z_xxxz_xyyy, to_0_z_xxxz_xyyz, to_0_z_xxxz_xyzz, to_0_z_xxxz_xzzz, to_0_z_xxxz_yyyy, to_0_z_xxxz_yyyz, to_0_z_xxxz_yyzz, to_0_z_xxxz_yzzz, to_0_z_xxxz_zzzz, to_xxxz_xxx, to_xxxz_xxxxz, to_xxxz_xxxyz, to_xxxz_xxxzz, to_xxxz_xxy, to_xxxz_xxyyz, to_xxxz_xxyzz, to_xxxz_xxz, to_xxxz_xxzzz, to_xxxz_xyy, to_xxxz_xyyyz, to_xxxz_xyyzz, to_xxxz_xyz, to_xxxz_xyzzz, to_xxxz_xzz, to_xxxz_xzzzz, to_xxxz_yyy, to_xxxz_yyyyz, to_xxxz_yyyzz, to_xxxz_yyz, to_xxxz_yyzzz, to_xxxz_yzz, to_xxxz_yzzzz, to_xxxz_zzz, to_xxxz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxz_xxxx[k] = 2.0 * to_xxxz_xxxxz[k] * tke_0;

            to_0_z_xxxz_xxxy[k] = 2.0 * to_xxxz_xxxyz[k] * tke_0;

            to_0_z_xxxz_xxxz[k] = -to_xxxz_xxx[k] + 2.0 * to_xxxz_xxxzz[k] * tke_0;

            to_0_z_xxxz_xxyy[k] = 2.0 * to_xxxz_xxyyz[k] * tke_0;

            to_0_z_xxxz_xxyz[k] = -to_xxxz_xxy[k] + 2.0 * to_xxxz_xxyzz[k] * tke_0;

            to_0_z_xxxz_xxzz[k] = -2.0 * to_xxxz_xxz[k] + 2.0 * to_xxxz_xxzzz[k] * tke_0;

            to_0_z_xxxz_xyyy[k] = 2.0 * to_xxxz_xyyyz[k] * tke_0;

            to_0_z_xxxz_xyyz[k] = -to_xxxz_xyy[k] + 2.0 * to_xxxz_xyyzz[k] * tke_0;

            to_0_z_xxxz_xyzz[k] = -2.0 * to_xxxz_xyz[k] + 2.0 * to_xxxz_xyzzz[k] * tke_0;

            to_0_z_xxxz_xzzz[k] = -3.0 * to_xxxz_xzz[k] + 2.0 * to_xxxz_xzzzz[k] * tke_0;

            to_0_z_xxxz_yyyy[k] = 2.0 * to_xxxz_yyyyz[k] * tke_0;

            to_0_z_xxxz_yyyz[k] = -to_xxxz_yyy[k] + 2.0 * to_xxxz_yyyzz[k] * tke_0;

            to_0_z_xxxz_yyzz[k] = -2.0 * to_xxxz_yyz[k] + 2.0 * to_xxxz_yyzzz[k] * tke_0;

            to_0_z_xxxz_yzzz[k] = -3.0 * to_xxxz_yzz[k] + 2.0 * to_xxxz_yzzzz[k] * tke_0;

            to_0_z_xxxz_zzzz[k] = -4.0 * to_xxxz_zzz[k] + 2.0 * to_xxxz_zzzzz[k] * tke_0;
        }

        // Set up 495-510 components of targeted buffer : GG

        auto to_0_z_xxyy_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 45);

        auto to_0_z_xxyy_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 46);

        auto to_0_z_xxyy_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 47);

        auto to_0_z_xxyy_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 48);

        auto to_0_z_xxyy_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 49);

        auto to_0_z_xxyy_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 50);

        auto to_0_z_xxyy_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 51);

        auto to_0_z_xxyy_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 52);

        auto to_0_z_xxyy_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 53);

        auto to_0_z_xxyy_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 54);

        auto to_0_z_xxyy_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 55);

        auto to_0_z_xxyy_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 56);

        auto to_0_z_xxyy_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 57);

        auto to_0_z_xxyy_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 58);

        auto to_0_z_xxyy_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 59);

        #pragma omp simd aligned(to_0_z_xxyy_xxxx, to_0_z_xxyy_xxxy, to_0_z_xxyy_xxxz, to_0_z_xxyy_xxyy, to_0_z_xxyy_xxyz, to_0_z_xxyy_xxzz, to_0_z_xxyy_xyyy, to_0_z_xxyy_xyyz, to_0_z_xxyy_xyzz, to_0_z_xxyy_xzzz, to_0_z_xxyy_yyyy, to_0_z_xxyy_yyyz, to_0_z_xxyy_yyzz, to_0_z_xxyy_yzzz, to_0_z_xxyy_zzzz, to_xxyy_xxx, to_xxyy_xxxxz, to_xxyy_xxxyz, to_xxyy_xxxzz, to_xxyy_xxy, to_xxyy_xxyyz, to_xxyy_xxyzz, to_xxyy_xxz, to_xxyy_xxzzz, to_xxyy_xyy, to_xxyy_xyyyz, to_xxyy_xyyzz, to_xxyy_xyz, to_xxyy_xyzzz, to_xxyy_xzz, to_xxyy_xzzzz, to_xxyy_yyy, to_xxyy_yyyyz, to_xxyy_yyyzz, to_xxyy_yyz, to_xxyy_yyzzz, to_xxyy_yzz, to_xxyy_yzzzz, to_xxyy_zzz, to_xxyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxyy_xxxx[k] = 2.0 * to_xxyy_xxxxz[k] * tke_0;

            to_0_z_xxyy_xxxy[k] = 2.0 * to_xxyy_xxxyz[k] * tke_0;

            to_0_z_xxyy_xxxz[k] = -to_xxyy_xxx[k] + 2.0 * to_xxyy_xxxzz[k] * tke_0;

            to_0_z_xxyy_xxyy[k] = 2.0 * to_xxyy_xxyyz[k] * tke_0;

            to_0_z_xxyy_xxyz[k] = -to_xxyy_xxy[k] + 2.0 * to_xxyy_xxyzz[k] * tke_0;

            to_0_z_xxyy_xxzz[k] = -2.0 * to_xxyy_xxz[k] + 2.0 * to_xxyy_xxzzz[k] * tke_0;

            to_0_z_xxyy_xyyy[k] = 2.0 * to_xxyy_xyyyz[k] * tke_0;

            to_0_z_xxyy_xyyz[k] = -to_xxyy_xyy[k] + 2.0 * to_xxyy_xyyzz[k] * tke_0;

            to_0_z_xxyy_xyzz[k] = -2.0 * to_xxyy_xyz[k] + 2.0 * to_xxyy_xyzzz[k] * tke_0;

            to_0_z_xxyy_xzzz[k] = -3.0 * to_xxyy_xzz[k] + 2.0 * to_xxyy_xzzzz[k] * tke_0;

            to_0_z_xxyy_yyyy[k] = 2.0 * to_xxyy_yyyyz[k] * tke_0;

            to_0_z_xxyy_yyyz[k] = -to_xxyy_yyy[k] + 2.0 * to_xxyy_yyyzz[k] * tke_0;

            to_0_z_xxyy_yyzz[k] = -2.0 * to_xxyy_yyz[k] + 2.0 * to_xxyy_yyzzz[k] * tke_0;

            to_0_z_xxyy_yzzz[k] = -3.0 * to_xxyy_yzz[k] + 2.0 * to_xxyy_yzzzz[k] * tke_0;

            to_0_z_xxyy_zzzz[k] = -4.0 * to_xxyy_zzz[k] + 2.0 * to_xxyy_zzzzz[k] * tke_0;
        }

        // Set up 510-525 components of targeted buffer : GG

        auto to_0_z_xxyz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 60);

        auto to_0_z_xxyz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 61);

        auto to_0_z_xxyz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 62);

        auto to_0_z_xxyz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 63);

        auto to_0_z_xxyz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 64);

        auto to_0_z_xxyz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 65);

        auto to_0_z_xxyz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 66);

        auto to_0_z_xxyz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 67);

        auto to_0_z_xxyz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 68);

        auto to_0_z_xxyz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 69);

        auto to_0_z_xxyz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 70);

        auto to_0_z_xxyz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 71);

        auto to_0_z_xxyz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 72);

        auto to_0_z_xxyz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 73);

        auto to_0_z_xxyz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 74);

        #pragma omp simd aligned(to_0_z_xxyz_xxxx, to_0_z_xxyz_xxxy, to_0_z_xxyz_xxxz, to_0_z_xxyz_xxyy, to_0_z_xxyz_xxyz, to_0_z_xxyz_xxzz, to_0_z_xxyz_xyyy, to_0_z_xxyz_xyyz, to_0_z_xxyz_xyzz, to_0_z_xxyz_xzzz, to_0_z_xxyz_yyyy, to_0_z_xxyz_yyyz, to_0_z_xxyz_yyzz, to_0_z_xxyz_yzzz, to_0_z_xxyz_zzzz, to_xxyz_xxx, to_xxyz_xxxxz, to_xxyz_xxxyz, to_xxyz_xxxzz, to_xxyz_xxy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xxzzz, to_xxyz_xyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_xzzzz, to_xxyz_yyy, to_xxyz_yyyyz, to_xxyz_yyyzz, to_xxyz_yyz, to_xxyz_yyzzz, to_xxyz_yzz, to_xxyz_yzzzz, to_xxyz_zzz, to_xxyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxyz_xxxx[k] = 2.0 * to_xxyz_xxxxz[k] * tke_0;

            to_0_z_xxyz_xxxy[k] = 2.0 * to_xxyz_xxxyz[k] * tke_0;

            to_0_z_xxyz_xxxz[k] = -to_xxyz_xxx[k] + 2.0 * to_xxyz_xxxzz[k] * tke_0;

            to_0_z_xxyz_xxyy[k] = 2.0 * to_xxyz_xxyyz[k] * tke_0;

            to_0_z_xxyz_xxyz[k] = -to_xxyz_xxy[k] + 2.0 * to_xxyz_xxyzz[k] * tke_0;

            to_0_z_xxyz_xxzz[k] = -2.0 * to_xxyz_xxz[k] + 2.0 * to_xxyz_xxzzz[k] * tke_0;

            to_0_z_xxyz_xyyy[k] = 2.0 * to_xxyz_xyyyz[k] * tke_0;

            to_0_z_xxyz_xyyz[k] = -to_xxyz_xyy[k] + 2.0 * to_xxyz_xyyzz[k] * tke_0;

            to_0_z_xxyz_xyzz[k] = -2.0 * to_xxyz_xyz[k] + 2.0 * to_xxyz_xyzzz[k] * tke_0;

            to_0_z_xxyz_xzzz[k] = -3.0 * to_xxyz_xzz[k] + 2.0 * to_xxyz_xzzzz[k] * tke_0;

            to_0_z_xxyz_yyyy[k] = 2.0 * to_xxyz_yyyyz[k] * tke_0;

            to_0_z_xxyz_yyyz[k] = -to_xxyz_yyy[k] + 2.0 * to_xxyz_yyyzz[k] * tke_0;

            to_0_z_xxyz_yyzz[k] = -2.0 * to_xxyz_yyz[k] + 2.0 * to_xxyz_yyzzz[k] * tke_0;

            to_0_z_xxyz_yzzz[k] = -3.0 * to_xxyz_yzz[k] + 2.0 * to_xxyz_yzzzz[k] * tke_0;

            to_0_z_xxyz_zzzz[k] = -4.0 * to_xxyz_zzz[k] + 2.0 * to_xxyz_zzzzz[k] * tke_0;
        }

        // Set up 525-540 components of targeted buffer : GG

        auto to_0_z_xxzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 75);

        auto to_0_z_xxzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 76);

        auto to_0_z_xxzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 77);

        auto to_0_z_xxzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 78);

        auto to_0_z_xxzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 79);

        auto to_0_z_xxzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 80);

        auto to_0_z_xxzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 81);

        auto to_0_z_xxzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 82);

        auto to_0_z_xxzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 83);

        auto to_0_z_xxzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 84);

        auto to_0_z_xxzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 85);

        auto to_0_z_xxzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 86);

        auto to_0_z_xxzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 87);

        auto to_0_z_xxzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 88);

        auto to_0_z_xxzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 89);

        #pragma omp simd aligned(to_0_z_xxzz_xxxx, to_0_z_xxzz_xxxy, to_0_z_xxzz_xxxz, to_0_z_xxzz_xxyy, to_0_z_xxzz_xxyz, to_0_z_xxzz_xxzz, to_0_z_xxzz_xyyy, to_0_z_xxzz_xyyz, to_0_z_xxzz_xyzz, to_0_z_xxzz_xzzz, to_0_z_xxzz_yyyy, to_0_z_xxzz_yyyz, to_0_z_xxzz_yyzz, to_0_z_xxzz_yzzz, to_0_z_xxzz_zzzz, to_xxzz_xxx, to_xxzz_xxxxz, to_xxzz_xxxyz, to_xxzz_xxxzz, to_xxzz_xxy, to_xxzz_xxyyz, to_xxzz_xxyzz, to_xxzz_xxz, to_xxzz_xxzzz, to_xxzz_xyy, to_xxzz_xyyyz, to_xxzz_xyyzz, to_xxzz_xyz, to_xxzz_xyzzz, to_xxzz_xzz, to_xxzz_xzzzz, to_xxzz_yyy, to_xxzz_yyyyz, to_xxzz_yyyzz, to_xxzz_yyz, to_xxzz_yyzzz, to_xxzz_yzz, to_xxzz_yzzzz, to_xxzz_zzz, to_xxzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxzz_xxxx[k] = 2.0 * to_xxzz_xxxxz[k] * tke_0;

            to_0_z_xxzz_xxxy[k] = 2.0 * to_xxzz_xxxyz[k] * tke_0;

            to_0_z_xxzz_xxxz[k] = -to_xxzz_xxx[k] + 2.0 * to_xxzz_xxxzz[k] * tke_0;

            to_0_z_xxzz_xxyy[k] = 2.0 * to_xxzz_xxyyz[k] * tke_0;

            to_0_z_xxzz_xxyz[k] = -to_xxzz_xxy[k] + 2.0 * to_xxzz_xxyzz[k] * tke_0;

            to_0_z_xxzz_xxzz[k] = -2.0 * to_xxzz_xxz[k] + 2.0 * to_xxzz_xxzzz[k] * tke_0;

            to_0_z_xxzz_xyyy[k] = 2.0 * to_xxzz_xyyyz[k] * tke_0;

            to_0_z_xxzz_xyyz[k] = -to_xxzz_xyy[k] + 2.0 * to_xxzz_xyyzz[k] * tke_0;

            to_0_z_xxzz_xyzz[k] = -2.0 * to_xxzz_xyz[k] + 2.0 * to_xxzz_xyzzz[k] * tke_0;

            to_0_z_xxzz_xzzz[k] = -3.0 * to_xxzz_xzz[k] + 2.0 * to_xxzz_xzzzz[k] * tke_0;

            to_0_z_xxzz_yyyy[k] = 2.0 * to_xxzz_yyyyz[k] * tke_0;

            to_0_z_xxzz_yyyz[k] = -to_xxzz_yyy[k] + 2.0 * to_xxzz_yyyzz[k] * tke_0;

            to_0_z_xxzz_yyzz[k] = -2.0 * to_xxzz_yyz[k] + 2.0 * to_xxzz_yyzzz[k] * tke_0;

            to_0_z_xxzz_yzzz[k] = -3.0 * to_xxzz_yzz[k] + 2.0 * to_xxzz_yzzzz[k] * tke_0;

            to_0_z_xxzz_zzzz[k] = -4.0 * to_xxzz_zzz[k] + 2.0 * to_xxzz_zzzzz[k] * tke_0;
        }

        // Set up 540-555 components of targeted buffer : GG

        auto to_0_z_xyyy_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 90);

        auto to_0_z_xyyy_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 91);

        auto to_0_z_xyyy_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 92);

        auto to_0_z_xyyy_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 93);

        auto to_0_z_xyyy_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 94);

        auto to_0_z_xyyy_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 95);

        auto to_0_z_xyyy_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 96);

        auto to_0_z_xyyy_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 97);

        auto to_0_z_xyyy_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 98);

        auto to_0_z_xyyy_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 99);

        auto to_0_z_xyyy_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 100);

        auto to_0_z_xyyy_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 101);

        auto to_0_z_xyyy_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 102);

        auto to_0_z_xyyy_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 103);

        auto to_0_z_xyyy_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 104);

        #pragma omp simd aligned(to_0_z_xyyy_xxxx, to_0_z_xyyy_xxxy, to_0_z_xyyy_xxxz, to_0_z_xyyy_xxyy, to_0_z_xyyy_xxyz, to_0_z_xyyy_xxzz, to_0_z_xyyy_xyyy, to_0_z_xyyy_xyyz, to_0_z_xyyy_xyzz, to_0_z_xyyy_xzzz, to_0_z_xyyy_yyyy, to_0_z_xyyy_yyyz, to_0_z_xyyy_yyzz, to_0_z_xyyy_yzzz, to_0_z_xyyy_zzzz, to_xyyy_xxx, to_xyyy_xxxxz, to_xyyy_xxxyz, to_xyyy_xxxzz, to_xyyy_xxy, to_xyyy_xxyyz, to_xyyy_xxyzz, to_xyyy_xxz, to_xyyy_xxzzz, to_xyyy_xyy, to_xyyy_xyyyz, to_xyyy_xyyzz, to_xyyy_xyz, to_xyyy_xyzzz, to_xyyy_xzz, to_xyyy_xzzzz, to_xyyy_yyy, to_xyyy_yyyyz, to_xyyy_yyyzz, to_xyyy_yyz, to_xyyy_yyzzz, to_xyyy_yzz, to_xyyy_yzzzz, to_xyyy_zzz, to_xyyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyyy_xxxx[k] = 2.0 * to_xyyy_xxxxz[k] * tke_0;

            to_0_z_xyyy_xxxy[k] = 2.0 * to_xyyy_xxxyz[k] * tke_0;

            to_0_z_xyyy_xxxz[k] = -to_xyyy_xxx[k] + 2.0 * to_xyyy_xxxzz[k] * tke_0;

            to_0_z_xyyy_xxyy[k] = 2.0 * to_xyyy_xxyyz[k] * tke_0;

            to_0_z_xyyy_xxyz[k] = -to_xyyy_xxy[k] + 2.0 * to_xyyy_xxyzz[k] * tke_0;

            to_0_z_xyyy_xxzz[k] = -2.0 * to_xyyy_xxz[k] + 2.0 * to_xyyy_xxzzz[k] * tke_0;

            to_0_z_xyyy_xyyy[k] = 2.0 * to_xyyy_xyyyz[k] * tke_0;

            to_0_z_xyyy_xyyz[k] = -to_xyyy_xyy[k] + 2.0 * to_xyyy_xyyzz[k] * tke_0;

            to_0_z_xyyy_xyzz[k] = -2.0 * to_xyyy_xyz[k] + 2.0 * to_xyyy_xyzzz[k] * tke_0;

            to_0_z_xyyy_xzzz[k] = -3.0 * to_xyyy_xzz[k] + 2.0 * to_xyyy_xzzzz[k] * tke_0;

            to_0_z_xyyy_yyyy[k] = 2.0 * to_xyyy_yyyyz[k] * tke_0;

            to_0_z_xyyy_yyyz[k] = -to_xyyy_yyy[k] + 2.0 * to_xyyy_yyyzz[k] * tke_0;

            to_0_z_xyyy_yyzz[k] = -2.0 * to_xyyy_yyz[k] + 2.0 * to_xyyy_yyzzz[k] * tke_0;

            to_0_z_xyyy_yzzz[k] = -3.0 * to_xyyy_yzz[k] + 2.0 * to_xyyy_yzzzz[k] * tke_0;

            to_0_z_xyyy_zzzz[k] = -4.0 * to_xyyy_zzz[k] + 2.0 * to_xyyy_zzzzz[k] * tke_0;
        }

        // Set up 555-570 components of targeted buffer : GG

        auto to_0_z_xyyz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 105);

        auto to_0_z_xyyz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 106);

        auto to_0_z_xyyz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 107);

        auto to_0_z_xyyz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 108);

        auto to_0_z_xyyz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 109);

        auto to_0_z_xyyz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 110);

        auto to_0_z_xyyz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 111);

        auto to_0_z_xyyz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 112);

        auto to_0_z_xyyz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 113);

        auto to_0_z_xyyz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 114);

        auto to_0_z_xyyz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 115);

        auto to_0_z_xyyz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 116);

        auto to_0_z_xyyz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 117);

        auto to_0_z_xyyz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 118);

        auto to_0_z_xyyz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 119);

        #pragma omp simd aligned(to_0_z_xyyz_xxxx, to_0_z_xyyz_xxxy, to_0_z_xyyz_xxxz, to_0_z_xyyz_xxyy, to_0_z_xyyz_xxyz, to_0_z_xyyz_xxzz, to_0_z_xyyz_xyyy, to_0_z_xyyz_xyyz, to_0_z_xyyz_xyzz, to_0_z_xyyz_xzzz, to_0_z_xyyz_yyyy, to_0_z_xyyz_yyyz, to_0_z_xyyz_yyzz, to_0_z_xyyz_yzzz, to_0_z_xyyz_zzzz, to_xyyz_xxx, to_xyyz_xxxxz, to_xyyz_xxxyz, to_xyyz_xxxzz, to_xyyz_xxy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xxzzz, to_xyyz_xyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_xzzzz, to_xyyz_yyy, to_xyyz_yyyyz, to_xyyz_yyyzz, to_xyyz_yyz, to_xyyz_yyzzz, to_xyyz_yzz, to_xyyz_yzzzz, to_xyyz_zzz, to_xyyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyyz_xxxx[k] = 2.0 * to_xyyz_xxxxz[k] * tke_0;

            to_0_z_xyyz_xxxy[k] = 2.0 * to_xyyz_xxxyz[k] * tke_0;

            to_0_z_xyyz_xxxz[k] = -to_xyyz_xxx[k] + 2.0 * to_xyyz_xxxzz[k] * tke_0;

            to_0_z_xyyz_xxyy[k] = 2.0 * to_xyyz_xxyyz[k] * tke_0;

            to_0_z_xyyz_xxyz[k] = -to_xyyz_xxy[k] + 2.0 * to_xyyz_xxyzz[k] * tke_0;

            to_0_z_xyyz_xxzz[k] = -2.0 * to_xyyz_xxz[k] + 2.0 * to_xyyz_xxzzz[k] * tke_0;

            to_0_z_xyyz_xyyy[k] = 2.0 * to_xyyz_xyyyz[k] * tke_0;

            to_0_z_xyyz_xyyz[k] = -to_xyyz_xyy[k] + 2.0 * to_xyyz_xyyzz[k] * tke_0;

            to_0_z_xyyz_xyzz[k] = -2.0 * to_xyyz_xyz[k] + 2.0 * to_xyyz_xyzzz[k] * tke_0;

            to_0_z_xyyz_xzzz[k] = -3.0 * to_xyyz_xzz[k] + 2.0 * to_xyyz_xzzzz[k] * tke_0;

            to_0_z_xyyz_yyyy[k] = 2.0 * to_xyyz_yyyyz[k] * tke_0;

            to_0_z_xyyz_yyyz[k] = -to_xyyz_yyy[k] + 2.0 * to_xyyz_yyyzz[k] * tke_0;

            to_0_z_xyyz_yyzz[k] = -2.0 * to_xyyz_yyz[k] + 2.0 * to_xyyz_yyzzz[k] * tke_0;

            to_0_z_xyyz_yzzz[k] = -3.0 * to_xyyz_yzz[k] + 2.0 * to_xyyz_yzzzz[k] * tke_0;

            to_0_z_xyyz_zzzz[k] = -4.0 * to_xyyz_zzz[k] + 2.0 * to_xyyz_zzzzz[k] * tke_0;
        }

        // Set up 570-585 components of targeted buffer : GG

        auto to_0_z_xyzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 120);

        auto to_0_z_xyzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 121);

        auto to_0_z_xyzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 122);

        auto to_0_z_xyzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 123);

        auto to_0_z_xyzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 124);

        auto to_0_z_xyzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 125);

        auto to_0_z_xyzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 126);

        auto to_0_z_xyzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 127);

        auto to_0_z_xyzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 128);

        auto to_0_z_xyzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 129);

        auto to_0_z_xyzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 130);

        auto to_0_z_xyzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 131);

        auto to_0_z_xyzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 132);

        auto to_0_z_xyzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 133);

        auto to_0_z_xyzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 134);

        #pragma omp simd aligned(to_0_z_xyzz_xxxx, to_0_z_xyzz_xxxy, to_0_z_xyzz_xxxz, to_0_z_xyzz_xxyy, to_0_z_xyzz_xxyz, to_0_z_xyzz_xxzz, to_0_z_xyzz_xyyy, to_0_z_xyzz_xyyz, to_0_z_xyzz_xyzz, to_0_z_xyzz_xzzz, to_0_z_xyzz_yyyy, to_0_z_xyzz_yyyz, to_0_z_xyzz_yyzz, to_0_z_xyzz_yzzz, to_0_z_xyzz_zzzz, to_xyzz_xxx, to_xyzz_xxxxz, to_xyzz_xxxyz, to_xyzz_xxxzz, to_xyzz_xxy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xxzzz, to_xyzz_xyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_xzzzz, to_xyzz_yyy, to_xyzz_yyyyz, to_xyzz_yyyzz, to_xyzz_yyz, to_xyzz_yyzzz, to_xyzz_yzz, to_xyzz_yzzzz, to_xyzz_zzz, to_xyzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyzz_xxxx[k] = 2.0 * to_xyzz_xxxxz[k] * tke_0;

            to_0_z_xyzz_xxxy[k] = 2.0 * to_xyzz_xxxyz[k] * tke_0;

            to_0_z_xyzz_xxxz[k] = -to_xyzz_xxx[k] + 2.0 * to_xyzz_xxxzz[k] * tke_0;

            to_0_z_xyzz_xxyy[k] = 2.0 * to_xyzz_xxyyz[k] * tke_0;

            to_0_z_xyzz_xxyz[k] = -to_xyzz_xxy[k] + 2.0 * to_xyzz_xxyzz[k] * tke_0;

            to_0_z_xyzz_xxzz[k] = -2.0 * to_xyzz_xxz[k] + 2.0 * to_xyzz_xxzzz[k] * tke_0;

            to_0_z_xyzz_xyyy[k] = 2.0 * to_xyzz_xyyyz[k] * tke_0;

            to_0_z_xyzz_xyyz[k] = -to_xyzz_xyy[k] + 2.0 * to_xyzz_xyyzz[k] * tke_0;

            to_0_z_xyzz_xyzz[k] = -2.0 * to_xyzz_xyz[k] + 2.0 * to_xyzz_xyzzz[k] * tke_0;

            to_0_z_xyzz_xzzz[k] = -3.0 * to_xyzz_xzz[k] + 2.0 * to_xyzz_xzzzz[k] * tke_0;

            to_0_z_xyzz_yyyy[k] = 2.0 * to_xyzz_yyyyz[k] * tke_0;

            to_0_z_xyzz_yyyz[k] = -to_xyzz_yyy[k] + 2.0 * to_xyzz_yyyzz[k] * tke_0;

            to_0_z_xyzz_yyzz[k] = -2.0 * to_xyzz_yyz[k] + 2.0 * to_xyzz_yyzzz[k] * tke_0;

            to_0_z_xyzz_yzzz[k] = -3.0 * to_xyzz_yzz[k] + 2.0 * to_xyzz_yzzzz[k] * tke_0;

            to_0_z_xyzz_zzzz[k] = -4.0 * to_xyzz_zzz[k] + 2.0 * to_xyzz_zzzzz[k] * tke_0;
        }

        // Set up 585-600 components of targeted buffer : GG

        auto to_0_z_xzzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 135);

        auto to_0_z_xzzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 136);

        auto to_0_z_xzzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 137);

        auto to_0_z_xzzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 138);

        auto to_0_z_xzzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 139);

        auto to_0_z_xzzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 140);

        auto to_0_z_xzzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 141);

        auto to_0_z_xzzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 142);

        auto to_0_z_xzzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 143);

        auto to_0_z_xzzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 144);

        auto to_0_z_xzzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 145);

        auto to_0_z_xzzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 146);

        auto to_0_z_xzzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 147);

        auto to_0_z_xzzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 148);

        auto to_0_z_xzzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 149);

        #pragma omp simd aligned(to_0_z_xzzz_xxxx, to_0_z_xzzz_xxxy, to_0_z_xzzz_xxxz, to_0_z_xzzz_xxyy, to_0_z_xzzz_xxyz, to_0_z_xzzz_xxzz, to_0_z_xzzz_xyyy, to_0_z_xzzz_xyyz, to_0_z_xzzz_xyzz, to_0_z_xzzz_xzzz, to_0_z_xzzz_yyyy, to_0_z_xzzz_yyyz, to_0_z_xzzz_yyzz, to_0_z_xzzz_yzzz, to_0_z_xzzz_zzzz, to_xzzz_xxx, to_xzzz_xxxxz, to_xzzz_xxxyz, to_xzzz_xxxzz, to_xzzz_xxy, to_xzzz_xxyyz, to_xzzz_xxyzz, to_xzzz_xxz, to_xzzz_xxzzz, to_xzzz_xyy, to_xzzz_xyyyz, to_xzzz_xyyzz, to_xzzz_xyz, to_xzzz_xyzzz, to_xzzz_xzz, to_xzzz_xzzzz, to_xzzz_yyy, to_xzzz_yyyyz, to_xzzz_yyyzz, to_xzzz_yyz, to_xzzz_yyzzz, to_xzzz_yzz, to_xzzz_yzzzz, to_xzzz_zzz, to_xzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xzzz_xxxx[k] = 2.0 * to_xzzz_xxxxz[k] * tke_0;

            to_0_z_xzzz_xxxy[k] = 2.0 * to_xzzz_xxxyz[k] * tke_0;

            to_0_z_xzzz_xxxz[k] = -to_xzzz_xxx[k] + 2.0 * to_xzzz_xxxzz[k] * tke_0;

            to_0_z_xzzz_xxyy[k] = 2.0 * to_xzzz_xxyyz[k] * tke_0;

            to_0_z_xzzz_xxyz[k] = -to_xzzz_xxy[k] + 2.0 * to_xzzz_xxyzz[k] * tke_0;

            to_0_z_xzzz_xxzz[k] = -2.0 * to_xzzz_xxz[k] + 2.0 * to_xzzz_xxzzz[k] * tke_0;

            to_0_z_xzzz_xyyy[k] = 2.0 * to_xzzz_xyyyz[k] * tke_0;

            to_0_z_xzzz_xyyz[k] = -to_xzzz_xyy[k] + 2.0 * to_xzzz_xyyzz[k] * tke_0;

            to_0_z_xzzz_xyzz[k] = -2.0 * to_xzzz_xyz[k] + 2.0 * to_xzzz_xyzzz[k] * tke_0;

            to_0_z_xzzz_xzzz[k] = -3.0 * to_xzzz_xzz[k] + 2.0 * to_xzzz_xzzzz[k] * tke_0;

            to_0_z_xzzz_yyyy[k] = 2.0 * to_xzzz_yyyyz[k] * tke_0;

            to_0_z_xzzz_yyyz[k] = -to_xzzz_yyy[k] + 2.0 * to_xzzz_yyyzz[k] * tke_0;

            to_0_z_xzzz_yyzz[k] = -2.0 * to_xzzz_yyz[k] + 2.0 * to_xzzz_yyzzz[k] * tke_0;

            to_0_z_xzzz_yzzz[k] = -3.0 * to_xzzz_yzz[k] + 2.0 * to_xzzz_yzzzz[k] * tke_0;

            to_0_z_xzzz_zzzz[k] = -4.0 * to_xzzz_zzz[k] + 2.0 * to_xzzz_zzzzz[k] * tke_0;
        }

        // Set up 600-615 components of targeted buffer : GG

        auto to_0_z_yyyy_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 150);

        auto to_0_z_yyyy_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 151);

        auto to_0_z_yyyy_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 152);

        auto to_0_z_yyyy_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 153);

        auto to_0_z_yyyy_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 154);

        auto to_0_z_yyyy_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 155);

        auto to_0_z_yyyy_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 156);

        auto to_0_z_yyyy_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 157);

        auto to_0_z_yyyy_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 158);

        auto to_0_z_yyyy_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 159);

        auto to_0_z_yyyy_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 160);

        auto to_0_z_yyyy_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 161);

        auto to_0_z_yyyy_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 162);

        auto to_0_z_yyyy_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 163);

        auto to_0_z_yyyy_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 164);

        #pragma omp simd aligned(to_0_z_yyyy_xxxx, to_0_z_yyyy_xxxy, to_0_z_yyyy_xxxz, to_0_z_yyyy_xxyy, to_0_z_yyyy_xxyz, to_0_z_yyyy_xxzz, to_0_z_yyyy_xyyy, to_0_z_yyyy_xyyz, to_0_z_yyyy_xyzz, to_0_z_yyyy_xzzz, to_0_z_yyyy_yyyy, to_0_z_yyyy_yyyz, to_0_z_yyyy_yyzz, to_0_z_yyyy_yzzz, to_0_z_yyyy_zzzz, to_yyyy_xxx, to_yyyy_xxxxz, to_yyyy_xxxyz, to_yyyy_xxxzz, to_yyyy_xxy, to_yyyy_xxyyz, to_yyyy_xxyzz, to_yyyy_xxz, to_yyyy_xxzzz, to_yyyy_xyy, to_yyyy_xyyyz, to_yyyy_xyyzz, to_yyyy_xyz, to_yyyy_xyzzz, to_yyyy_xzz, to_yyyy_xzzzz, to_yyyy_yyy, to_yyyy_yyyyz, to_yyyy_yyyzz, to_yyyy_yyz, to_yyyy_yyzzz, to_yyyy_yzz, to_yyyy_yzzzz, to_yyyy_zzz, to_yyyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyyy_xxxx[k] = 2.0 * to_yyyy_xxxxz[k] * tke_0;

            to_0_z_yyyy_xxxy[k] = 2.0 * to_yyyy_xxxyz[k] * tke_0;

            to_0_z_yyyy_xxxz[k] = -to_yyyy_xxx[k] + 2.0 * to_yyyy_xxxzz[k] * tke_0;

            to_0_z_yyyy_xxyy[k] = 2.0 * to_yyyy_xxyyz[k] * tke_0;

            to_0_z_yyyy_xxyz[k] = -to_yyyy_xxy[k] + 2.0 * to_yyyy_xxyzz[k] * tke_0;

            to_0_z_yyyy_xxzz[k] = -2.0 * to_yyyy_xxz[k] + 2.0 * to_yyyy_xxzzz[k] * tke_0;

            to_0_z_yyyy_xyyy[k] = 2.0 * to_yyyy_xyyyz[k] * tke_0;

            to_0_z_yyyy_xyyz[k] = -to_yyyy_xyy[k] + 2.0 * to_yyyy_xyyzz[k] * tke_0;

            to_0_z_yyyy_xyzz[k] = -2.0 * to_yyyy_xyz[k] + 2.0 * to_yyyy_xyzzz[k] * tke_0;

            to_0_z_yyyy_xzzz[k] = -3.0 * to_yyyy_xzz[k] + 2.0 * to_yyyy_xzzzz[k] * tke_0;

            to_0_z_yyyy_yyyy[k] = 2.0 * to_yyyy_yyyyz[k] * tke_0;

            to_0_z_yyyy_yyyz[k] = -to_yyyy_yyy[k] + 2.0 * to_yyyy_yyyzz[k] * tke_0;

            to_0_z_yyyy_yyzz[k] = -2.0 * to_yyyy_yyz[k] + 2.0 * to_yyyy_yyzzz[k] * tke_0;

            to_0_z_yyyy_yzzz[k] = -3.0 * to_yyyy_yzz[k] + 2.0 * to_yyyy_yzzzz[k] * tke_0;

            to_0_z_yyyy_zzzz[k] = -4.0 * to_yyyy_zzz[k] + 2.0 * to_yyyy_zzzzz[k] * tke_0;
        }

        // Set up 615-630 components of targeted buffer : GG

        auto to_0_z_yyyz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 165);

        auto to_0_z_yyyz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 166);

        auto to_0_z_yyyz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 167);

        auto to_0_z_yyyz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 168);

        auto to_0_z_yyyz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 169);

        auto to_0_z_yyyz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 170);

        auto to_0_z_yyyz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 171);

        auto to_0_z_yyyz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 172);

        auto to_0_z_yyyz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 173);

        auto to_0_z_yyyz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 174);

        auto to_0_z_yyyz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 175);

        auto to_0_z_yyyz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 176);

        auto to_0_z_yyyz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 177);

        auto to_0_z_yyyz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 178);

        auto to_0_z_yyyz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 179);

        #pragma omp simd aligned(to_0_z_yyyz_xxxx, to_0_z_yyyz_xxxy, to_0_z_yyyz_xxxz, to_0_z_yyyz_xxyy, to_0_z_yyyz_xxyz, to_0_z_yyyz_xxzz, to_0_z_yyyz_xyyy, to_0_z_yyyz_xyyz, to_0_z_yyyz_xyzz, to_0_z_yyyz_xzzz, to_0_z_yyyz_yyyy, to_0_z_yyyz_yyyz, to_0_z_yyyz_yyzz, to_0_z_yyyz_yzzz, to_0_z_yyyz_zzzz, to_yyyz_xxx, to_yyyz_xxxxz, to_yyyz_xxxyz, to_yyyz_xxxzz, to_yyyz_xxy, to_yyyz_xxyyz, to_yyyz_xxyzz, to_yyyz_xxz, to_yyyz_xxzzz, to_yyyz_xyy, to_yyyz_xyyyz, to_yyyz_xyyzz, to_yyyz_xyz, to_yyyz_xyzzz, to_yyyz_xzz, to_yyyz_xzzzz, to_yyyz_yyy, to_yyyz_yyyyz, to_yyyz_yyyzz, to_yyyz_yyz, to_yyyz_yyzzz, to_yyyz_yzz, to_yyyz_yzzzz, to_yyyz_zzz, to_yyyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyyz_xxxx[k] = 2.0 * to_yyyz_xxxxz[k] * tke_0;

            to_0_z_yyyz_xxxy[k] = 2.0 * to_yyyz_xxxyz[k] * tke_0;

            to_0_z_yyyz_xxxz[k] = -to_yyyz_xxx[k] + 2.0 * to_yyyz_xxxzz[k] * tke_0;

            to_0_z_yyyz_xxyy[k] = 2.0 * to_yyyz_xxyyz[k] * tke_0;

            to_0_z_yyyz_xxyz[k] = -to_yyyz_xxy[k] + 2.0 * to_yyyz_xxyzz[k] * tke_0;

            to_0_z_yyyz_xxzz[k] = -2.0 * to_yyyz_xxz[k] + 2.0 * to_yyyz_xxzzz[k] * tke_0;

            to_0_z_yyyz_xyyy[k] = 2.0 * to_yyyz_xyyyz[k] * tke_0;

            to_0_z_yyyz_xyyz[k] = -to_yyyz_xyy[k] + 2.0 * to_yyyz_xyyzz[k] * tke_0;

            to_0_z_yyyz_xyzz[k] = -2.0 * to_yyyz_xyz[k] + 2.0 * to_yyyz_xyzzz[k] * tke_0;

            to_0_z_yyyz_xzzz[k] = -3.0 * to_yyyz_xzz[k] + 2.0 * to_yyyz_xzzzz[k] * tke_0;

            to_0_z_yyyz_yyyy[k] = 2.0 * to_yyyz_yyyyz[k] * tke_0;

            to_0_z_yyyz_yyyz[k] = -to_yyyz_yyy[k] + 2.0 * to_yyyz_yyyzz[k] * tke_0;

            to_0_z_yyyz_yyzz[k] = -2.0 * to_yyyz_yyz[k] + 2.0 * to_yyyz_yyzzz[k] * tke_0;

            to_0_z_yyyz_yzzz[k] = -3.0 * to_yyyz_yzz[k] + 2.0 * to_yyyz_yzzzz[k] * tke_0;

            to_0_z_yyyz_zzzz[k] = -4.0 * to_yyyz_zzz[k] + 2.0 * to_yyyz_zzzzz[k] * tke_0;
        }

        // Set up 630-645 components of targeted buffer : GG

        auto to_0_z_yyzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 180);

        auto to_0_z_yyzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 181);

        auto to_0_z_yyzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 182);

        auto to_0_z_yyzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 183);

        auto to_0_z_yyzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 184);

        auto to_0_z_yyzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 185);

        auto to_0_z_yyzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 186);

        auto to_0_z_yyzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 187);

        auto to_0_z_yyzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 188);

        auto to_0_z_yyzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 189);

        auto to_0_z_yyzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 190);

        auto to_0_z_yyzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 191);

        auto to_0_z_yyzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 192);

        auto to_0_z_yyzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 193);

        auto to_0_z_yyzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 194);

        #pragma omp simd aligned(to_0_z_yyzz_xxxx, to_0_z_yyzz_xxxy, to_0_z_yyzz_xxxz, to_0_z_yyzz_xxyy, to_0_z_yyzz_xxyz, to_0_z_yyzz_xxzz, to_0_z_yyzz_xyyy, to_0_z_yyzz_xyyz, to_0_z_yyzz_xyzz, to_0_z_yyzz_xzzz, to_0_z_yyzz_yyyy, to_0_z_yyzz_yyyz, to_0_z_yyzz_yyzz, to_0_z_yyzz_yzzz, to_0_z_yyzz_zzzz, to_yyzz_xxx, to_yyzz_xxxxz, to_yyzz_xxxyz, to_yyzz_xxxzz, to_yyzz_xxy, to_yyzz_xxyyz, to_yyzz_xxyzz, to_yyzz_xxz, to_yyzz_xxzzz, to_yyzz_xyy, to_yyzz_xyyyz, to_yyzz_xyyzz, to_yyzz_xyz, to_yyzz_xyzzz, to_yyzz_xzz, to_yyzz_xzzzz, to_yyzz_yyy, to_yyzz_yyyyz, to_yyzz_yyyzz, to_yyzz_yyz, to_yyzz_yyzzz, to_yyzz_yzz, to_yyzz_yzzzz, to_yyzz_zzz, to_yyzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyzz_xxxx[k] = 2.0 * to_yyzz_xxxxz[k] * tke_0;

            to_0_z_yyzz_xxxy[k] = 2.0 * to_yyzz_xxxyz[k] * tke_0;

            to_0_z_yyzz_xxxz[k] = -to_yyzz_xxx[k] + 2.0 * to_yyzz_xxxzz[k] * tke_0;

            to_0_z_yyzz_xxyy[k] = 2.0 * to_yyzz_xxyyz[k] * tke_0;

            to_0_z_yyzz_xxyz[k] = -to_yyzz_xxy[k] + 2.0 * to_yyzz_xxyzz[k] * tke_0;

            to_0_z_yyzz_xxzz[k] = -2.0 * to_yyzz_xxz[k] + 2.0 * to_yyzz_xxzzz[k] * tke_0;

            to_0_z_yyzz_xyyy[k] = 2.0 * to_yyzz_xyyyz[k] * tke_0;

            to_0_z_yyzz_xyyz[k] = -to_yyzz_xyy[k] + 2.0 * to_yyzz_xyyzz[k] * tke_0;

            to_0_z_yyzz_xyzz[k] = -2.0 * to_yyzz_xyz[k] + 2.0 * to_yyzz_xyzzz[k] * tke_0;

            to_0_z_yyzz_xzzz[k] = -3.0 * to_yyzz_xzz[k] + 2.0 * to_yyzz_xzzzz[k] * tke_0;

            to_0_z_yyzz_yyyy[k] = 2.0 * to_yyzz_yyyyz[k] * tke_0;

            to_0_z_yyzz_yyyz[k] = -to_yyzz_yyy[k] + 2.0 * to_yyzz_yyyzz[k] * tke_0;

            to_0_z_yyzz_yyzz[k] = -2.0 * to_yyzz_yyz[k] + 2.0 * to_yyzz_yyzzz[k] * tke_0;

            to_0_z_yyzz_yzzz[k] = -3.0 * to_yyzz_yzz[k] + 2.0 * to_yyzz_yzzzz[k] * tke_0;

            to_0_z_yyzz_zzzz[k] = -4.0 * to_yyzz_zzz[k] + 2.0 * to_yyzz_zzzzz[k] * tke_0;
        }

        // Set up 645-660 components of targeted buffer : GG

        auto to_0_z_yzzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 195);

        auto to_0_z_yzzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 196);

        auto to_0_z_yzzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 197);

        auto to_0_z_yzzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 198);

        auto to_0_z_yzzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 199);

        auto to_0_z_yzzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 200);

        auto to_0_z_yzzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 201);

        auto to_0_z_yzzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 202);

        auto to_0_z_yzzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 203);

        auto to_0_z_yzzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 204);

        auto to_0_z_yzzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 205);

        auto to_0_z_yzzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 206);

        auto to_0_z_yzzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 207);

        auto to_0_z_yzzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 208);

        auto to_0_z_yzzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 209);

        #pragma omp simd aligned(to_0_z_yzzz_xxxx, to_0_z_yzzz_xxxy, to_0_z_yzzz_xxxz, to_0_z_yzzz_xxyy, to_0_z_yzzz_xxyz, to_0_z_yzzz_xxzz, to_0_z_yzzz_xyyy, to_0_z_yzzz_xyyz, to_0_z_yzzz_xyzz, to_0_z_yzzz_xzzz, to_0_z_yzzz_yyyy, to_0_z_yzzz_yyyz, to_0_z_yzzz_yyzz, to_0_z_yzzz_yzzz, to_0_z_yzzz_zzzz, to_yzzz_xxx, to_yzzz_xxxxz, to_yzzz_xxxyz, to_yzzz_xxxzz, to_yzzz_xxy, to_yzzz_xxyyz, to_yzzz_xxyzz, to_yzzz_xxz, to_yzzz_xxzzz, to_yzzz_xyy, to_yzzz_xyyyz, to_yzzz_xyyzz, to_yzzz_xyz, to_yzzz_xyzzz, to_yzzz_xzz, to_yzzz_xzzzz, to_yzzz_yyy, to_yzzz_yyyyz, to_yzzz_yyyzz, to_yzzz_yyz, to_yzzz_yyzzz, to_yzzz_yzz, to_yzzz_yzzzz, to_yzzz_zzz, to_yzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yzzz_xxxx[k] = 2.0 * to_yzzz_xxxxz[k] * tke_0;

            to_0_z_yzzz_xxxy[k] = 2.0 * to_yzzz_xxxyz[k] * tke_0;

            to_0_z_yzzz_xxxz[k] = -to_yzzz_xxx[k] + 2.0 * to_yzzz_xxxzz[k] * tke_0;

            to_0_z_yzzz_xxyy[k] = 2.0 * to_yzzz_xxyyz[k] * tke_0;

            to_0_z_yzzz_xxyz[k] = -to_yzzz_xxy[k] + 2.0 * to_yzzz_xxyzz[k] * tke_0;

            to_0_z_yzzz_xxzz[k] = -2.0 * to_yzzz_xxz[k] + 2.0 * to_yzzz_xxzzz[k] * tke_0;

            to_0_z_yzzz_xyyy[k] = 2.0 * to_yzzz_xyyyz[k] * tke_0;

            to_0_z_yzzz_xyyz[k] = -to_yzzz_xyy[k] + 2.0 * to_yzzz_xyyzz[k] * tke_0;

            to_0_z_yzzz_xyzz[k] = -2.0 * to_yzzz_xyz[k] + 2.0 * to_yzzz_xyzzz[k] * tke_0;

            to_0_z_yzzz_xzzz[k] = -3.0 * to_yzzz_xzz[k] + 2.0 * to_yzzz_xzzzz[k] * tke_0;

            to_0_z_yzzz_yyyy[k] = 2.0 * to_yzzz_yyyyz[k] * tke_0;

            to_0_z_yzzz_yyyz[k] = -to_yzzz_yyy[k] + 2.0 * to_yzzz_yyyzz[k] * tke_0;

            to_0_z_yzzz_yyzz[k] = -2.0 * to_yzzz_yyz[k] + 2.0 * to_yzzz_yyzzz[k] * tke_0;

            to_0_z_yzzz_yzzz[k] = -3.0 * to_yzzz_yzz[k] + 2.0 * to_yzzz_yzzzz[k] * tke_0;

            to_0_z_yzzz_zzzz[k] = -4.0 * to_yzzz_zzz[k] + 2.0 * to_yzzz_zzzzz[k] * tke_0;
        }

        // Set up 660-675 components of targeted buffer : GG

        auto to_0_z_zzzz_xxxx = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 210);

        auto to_0_z_zzzz_xxxy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 211);

        auto to_0_z_zzzz_xxxz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 212);

        auto to_0_z_zzzz_xxyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 213);

        auto to_0_z_zzzz_xxyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 214);

        auto to_0_z_zzzz_xxzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 215);

        auto to_0_z_zzzz_xyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 216);

        auto to_0_z_zzzz_xyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 217);

        auto to_0_z_zzzz_xyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 218);

        auto to_0_z_zzzz_xzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 219);

        auto to_0_z_zzzz_yyyy = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 220);

        auto to_0_z_zzzz_yyyz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 221);

        auto to_0_z_zzzz_yyzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 222);

        auto to_0_z_zzzz_yzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 223);

        auto to_0_z_zzzz_zzzz = pbuffer.data(idx_op_geom_001_gg + 2 * op_comps * 225 + i * 225 + 224);

        #pragma omp simd aligned(to_0_z_zzzz_xxxx, to_0_z_zzzz_xxxy, to_0_z_zzzz_xxxz, to_0_z_zzzz_xxyy, to_0_z_zzzz_xxyz, to_0_z_zzzz_xxzz, to_0_z_zzzz_xyyy, to_0_z_zzzz_xyyz, to_0_z_zzzz_xyzz, to_0_z_zzzz_xzzz, to_0_z_zzzz_yyyy, to_0_z_zzzz_yyyz, to_0_z_zzzz_yyzz, to_0_z_zzzz_yzzz, to_0_z_zzzz_zzzz, to_zzzz_xxx, to_zzzz_xxxxz, to_zzzz_xxxyz, to_zzzz_xxxzz, to_zzzz_xxy, to_zzzz_xxyyz, to_zzzz_xxyzz, to_zzzz_xxz, to_zzzz_xxzzz, to_zzzz_xyy, to_zzzz_xyyyz, to_zzzz_xyyzz, to_zzzz_xyz, to_zzzz_xyzzz, to_zzzz_xzz, to_zzzz_xzzzz, to_zzzz_yyy, to_zzzz_yyyyz, to_zzzz_yyyzz, to_zzzz_yyz, to_zzzz_yyzzz, to_zzzz_yzz, to_zzzz_yzzzz, to_zzzz_zzz, to_zzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zzzz_xxxx[k] = 2.0 * to_zzzz_xxxxz[k] * tke_0;

            to_0_z_zzzz_xxxy[k] = 2.0 * to_zzzz_xxxyz[k] * tke_0;

            to_0_z_zzzz_xxxz[k] = -to_zzzz_xxx[k] + 2.0 * to_zzzz_xxxzz[k] * tke_0;

            to_0_z_zzzz_xxyy[k] = 2.0 * to_zzzz_xxyyz[k] * tke_0;

            to_0_z_zzzz_xxyz[k] = -to_zzzz_xxy[k] + 2.0 * to_zzzz_xxyzz[k] * tke_0;

            to_0_z_zzzz_xxzz[k] = -2.0 * to_zzzz_xxz[k] + 2.0 * to_zzzz_xxzzz[k] * tke_0;

            to_0_z_zzzz_xyyy[k] = 2.0 * to_zzzz_xyyyz[k] * tke_0;

            to_0_z_zzzz_xyyz[k] = -to_zzzz_xyy[k] + 2.0 * to_zzzz_xyyzz[k] * tke_0;

            to_0_z_zzzz_xyzz[k] = -2.0 * to_zzzz_xyz[k] + 2.0 * to_zzzz_xyzzz[k] * tke_0;

            to_0_z_zzzz_xzzz[k] = -3.0 * to_zzzz_xzz[k] + 2.0 * to_zzzz_xzzzz[k] * tke_0;

            to_0_z_zzzz_yyyy[k] = 2.0 * to_zzzz_yyyyz[k] * tke_0;

            to_0_z_zzzz_yyyz[k] = -to_zzzz_yyy[k] + 2.0 * to_zzzz_yyyzz[k] * tke_0;

            to_0_z_zzzz_yyzz[k] = -2.0 * to_zzzz_yyz[k] + 2.0 * to_zzzz_yyzzz[k] * tke_0;

            to_0_z_zzzz_yzzz[k] = -3.0 * to_zzzz_yzz[k] + 2.0 * to_zzzz_yzzzz[k] * tke_0;

            to_0_z_zzzz_zzzz[k] = -4.0 * to_zzzz_zzz[k] + 2.0 * to_zzzz_zzzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

