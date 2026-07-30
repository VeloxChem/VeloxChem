#include "GeometricalDerivatives010ForGG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_gg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_gg,
                         const int idx_op_fg,
                         const int idx_op_gf,
                         const int idx_op_gh,
                         const int idx_op_hg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : FG

    auto tr_xxx_xxxx = pbuffer.data(idx_op_fg);

    auto tr_xxx_xxxy = pbuffer.data(idx_op_fg + 1);

    auto tr_xxx_xxxz = pbuffer.data(idx_op_fg + 2);

    auto tr_xxx_xxyy = pbuffer.data(idx_op_fg + 3);

    auto tr_xxx_xxyz = pbuffer.data(idx_op_fg + 4);

    auto tr_xxx_xxzz = pbuffer.data(idx_op_fg + 5);

    auto tr_xxx_xyyy = pbuffer.data(idx_op_fg + 6);

    auto tr_xxx_xyyz = pbuffer.data(idx_op_fg + 7);

    auto tr_xxx_xyzz = pbuffer.data(idx_op_fg + 8);

    auto tr_xxx_xzzz = pbuffer.data(idx_op_fg + 9);

    auto tr_xxx_yyyy = pbuffer.data(idx_op_fg + 10);

    auto tr_xxx_yyyz = pbuffer.data(idx_op_fg + 11);

    auto tr_xxx_yyzz = pbuffer.data(idx_op_fg + 12);

    auto tr_xxx_yzzz = pbuffer.data(idx_op_fg + 13);

    auto tr_xxx_zzzz = pbuffer.data(idx_op_fg + 14);

    auto tr_xxy_xxxx = pbuffer.data(idx_op_fg + 15);

    auto tr_xxy_xxxy = pbuffer.data(idx_op_fg + 16);

    auto tr_xxy_xxxz = pbuffer.data(idx_op_fg + 17);

    auto tr_xxy_xxyy = pbuffer.data(idx_op_fg + 18);

    auto tr_xxy_xxyz = pbuffer.data(idx_op_fg + 19);

    auto tr_xxy_xxzz = pbuffer.data(idx_op_fg + 20);

    auto tr_xxy_xyyy = pbuffer.data(idx_op_fg + 21);

    auto tr_xxy_xyyz = pbuffer.data(idx_op_fg + 22);

    auto tr_xxy_xyzz = pbuffer.data(idx_op_fg + 23);

    auto tr_xxy_xzzz = pbuffer.data(idx_op_fg + 24);

    auto tr_xxy_yyyy = pbuffer.data(idx_op_fg + 25);

    auto tr_xxy_yyyz = pbuffer.data(idx_op_fg + 26);

    auto tr_xxy_yyzz = pbuffer.data(idx_op_fg + 27);

    auto tr_xxy_yzzz = pbuffer.data(idx_op_fg + 28);

    auto tr_xxy_zzzz = pbuffer.data(idx_op_fg + 29);

    auto tr_xxz_xxxx = pbuffer.data(idx_op_fg + 30);

    auto tr_xxz_xxxy = pbuffer.data(idx_op_fg + 31);

    auto tr_xxz_xxxz = pbuffer.data(idx_op_fg + 32);

    auto tr_xxz_xxyy = pbuffer.data(idx_op_fg + 33);

    auto tr_xxz_xxyz = pbuffer.data(idx_op_fg + 34);

    auto tr_xxz_xxzz = pbuffer.data(idx_op_fg + 35);

    auto tr_xxz_xyyy = pbuffer.data(idx_op_fg + 36);

    auto tr_xxz_xyyz = pbuffer.data(idx_op_fg + 37);

    auto tr_xxz_xyzz = pbuffer.data(idx_op_fg + 38);

    auto tr_xxz_xzzz = pbuffer.data(idx_op_fg + 39);

    auto tr_xxz_yyyy = pbuffer.data(idx_op_fg + 40);

    auto tr_xxz_yyyz = pbuffer.data(idx_op_fg + 41);

    auto tr_xxz_yyzz = pbuffer.data(idx_op_fg + 42);

    auto tr_xxz_yzzz = pbuffer.data(idx_op_fg + 43);

    auto tr_xxz_zzzz = pbuffer.data(idx_op_fg + 44);

    auto tr_xyy_xxxx = pbuffer.data(idx_op_fg + 45);

    auto tr_xyy_xxxy = pbuffer.data(idx_op_fg + 46);

    auto tr_xyy_xxxz = pbuffer.data(idx_op_fg + 47);

    auto tr_xyy_xxyy = pbuffer.data(idx_op_fg + 48);

    auto tr_xyy_xxyz = pbuffer.data(idx_op_fg + 49);

    auto tr_xyy_xxzz = pbuffer.data(idx_op_fg + 50);

    auto tr_xyy_xyyy = pbuffer.data(idx_op_fg + 51);

    auto tr_xyy_xyyz = pbuffer.data(idx_op_fg + 52);

    auto tr_xyy_xyzz = pbuffer.data(idx_op_fg + 53);

    auto tr_xyy_xzzz = pbuffer.data(idx_op_fg + 54);

    auto tr_xyy_yyyy = pbuffer.data(idx_op_fg + 55);

    auto tr_xyy_yyyz = pbuffer.data(idx_op_fg + 56);

    auto tr_xyy_yyzz = pbuffer.data(idx_op_fg + 57);

    auto tr_xyy_yzzz = pbuffer.data(idx_op_fg + 58);

    auto tr_xyy_zzzz = pbuffer.data(idx_op_fg + 59);

    auto tr_xyz_xxxx = pbuffer.data(idx_op_fg + 60);

    auto tr_xyz_xxxy = pbuffer.data(idx_op_fg + 61);

    auto tr_xyz_xxxz = pbuffer.data(idx_op_fg + 62);

    auto tr_xyz_xxyy = pbuffer.data(idx_op_fg + 63);

    auto tr_xyz_xxyz = pbuffer.data(idx_op_fg + 64);

    auto tr_xyz_xxzz = pbuffer.data(idx_op_fg + 65);

    auto tr_xyz_xyyy = pbuffer.data(idx_op_fg + 66);

    auto tr_xyz_xyyz = pbuffer.data(idx_op_fg + 67);

    auto tr_xyz_xyzz = pbuffer.data(idx_op_fg + 68);

    auto tr_xyz_xzzz = pbuffer.data(idx_op_fg + 69);

    auto tr_xyz_yyyy = pbuffer.data(idx_op_fg + 70);

    auto tr_xyz_yyyz = pbuffer.data(idx_op_fg + 71);

    auto tr_xyz_yyzz = pbuffer.data(idx_op_fg + 72);

    auto tr_xyz_yzzz = pbuffer.data(idx_op_fg + 73);

    auto tr_xyz_zzzz = pbuffer.data(idx_op_fg + 74);

    auto tr_xzz_xxxx = pbuffer.data(idx_op_fg + 75);

    auto tr_xzz_xxxy = pbuffer.data(idx_op_fg + 76);

    auto tr_xzz_xxxz = pbuffer.data(idx_op_fg + 77);

    auto tr_xzz_xxyy = pbuffer.data(idx_op_fg + 78);

    auto tr_xzz_xxyz = pbuffer.data(idx_op_fg + 79);

    auto tr_xzz_xxzz = pbuffer.data(idx_op_fg + 80);

    auto tr_xzz_xyyy = pbuffer.data(idx_op_fg + 81);

    auto tr_xzz_xyyz = pbuffer.data(idx_op_fg + 82);

    auto tr_xzz_xyzz = pbuffer.data(idx_op_fg + 83);

    auto tr_xzz_xzzz = pbuffer.data(idx_op_fg + 84);

    auto tr_xzz_yyyy = pbuffer.data(idx_op_fg + 85);

    auto tr_xzz_yyyz = pbuffer.data(idx_op_fg + 86);

    auto tr_xzz_yyzz = pbuffer.data(idx_op_fg + 87);

    auto tr_xzz_yzzz = pbuffer.data(idx_op_fg + 88);

    auto tr_xzz_zzzz = pbuffer.data(idx_op_fg + 89);

    auto tr_yyy_xxxx = pbuffer.data(idx_op_fg + 90);

    auto tr_yyy_xxxy = pbuffer.data(idx_op_fg + 91);

    auto tr_yyy_xxxz = pbuffer.data(idx_op_fg + 92);

    auto tr_yyy_xxyy = pbuffer.data(idx_op_fg + 93);

    auto tr_yyy_xxyz = pbuffer.data(idx_op_fg + 94);

    auto tr_yyy_xxzz = pbuffer.data(idx_op_fg + 95);

    auto tr_yyy_xyyy = pbuffer.data(idx_op_fg + 96);

    auto tr_yyy_xyyz = pbuffer.data(idx_op_fg + 97);

    auto tr_yyy_xyzz = pbuffer.data(idx_op_fg + 98);

    auto tr_yyy_xzzz = pbuffer.data(idx_op_fg + 99);

    auto tr_yyy_yyyy = pbuffer.data(idx_op_fg + 100);

    auto tr_yyy_yyyz = pbuffer.data(idx_op_fg + 101);

    auto tr_yyy_yyzz = pbuffer.data(idx_op_fg + 102);

    auto tr_yyy_yzzz = pbuffer.data(idx_op_fg + 103);

    auto tr_yyy_zzzz = pbuffer.data(idx_op_fg + 104);

    auto tr_yyz_xxxx = pbuffer.data(idx_op_fg + 105);

    auto tr_yyz_xxxy = pbuffer.data(idx_op_fg + 106);

    auto tr_yyz_xxxz = pbuffer.data(idx_op_fg + 107);

    auto tr_yyz_xxyy = pbuffer.data(idx_op_fg + 108);

    auto tr_yyz_xxyz = pbuffer.data(idx_op_fg + 109);

    auto tr_yyz_xxzz = pbuffer.data(idx_op_fg + 110);

    auto tr_yyz_xyyy = pbuffer.data(idx_op_fg + 111);

    auto tr_yyz_xyyz = pbuffer.data(idx_op_fg + 112);

    auto tr_yyz_xyzz = pbuffer.data(idx_op_fg + 113);

    auto tr_yyz_xzzz = pbuffer.data(idx_op_fg + 114);

    auto tr_yyz_yyyy = pbuffer.data(idx_op_fg + 115);

    auto tr_yyz_yyyz = pbuffer.data(idx_op_fg + 116);

    auto tr_yyz_yyzz = pbuffer.data(idx_op_fg + 117);

    auto tr_yyz_yzzz = pbuffer.data(idx_op_fg + 118);

    auto tr_yyz_zzzz = pbuffer.data(idx_op_fg + 119);

    auto tr_yzz_xxxx = pbuffer.data(idx_op_fg + 120);

    auto tr_yzz_xxxy = pbuffer.data(idx_op_fg + 121);

    auto tr_yzz_xxxz = pbuffer.data(idx_op_fg + 122);

    auto tr_yzz_xxyy = pbuffer.data(idx_op_fg + 123);

    auto tr_yzz_xxyz = pbuffer.data(idx_op_fg + 124);

    auto tr_yzz_xxzz = pbuffer.data(idx_op_fg + 125);

    auto tr_yzz_xyyy = pbuffer.data(idx_op_fg + 126);

    auto tr_yzz_xyyz = pbuffer.data(idx_op_fg + 127);

    auto tr_yzz_xyzz = pbuffer.data(idx_op_fg + 128);

    auto tr_yzz_xzzz = pbuffer.data(idx_op_fg + 129);

    auto tr_yzz_yyyy = pbuffer.data(idx_op_fg + 130);

    auto tr_yzz_yyyz = pbuffer.data(idx_op_fg + 131);

    auto tr_yzz_yyzz = pbuffer.data(idx_op_fg + 132);

    auto tr_yzz_yzzz = pbuffer.data(idx_op_fg + 133);

    auto tr_yzz_zzzz = pbuffer.data(idx_op_fg + 134);

    auto tr_zzz_xxxx = pbuffer.data(idx_op_fg + 135);

    auto tr_zzz_xxxy = pbuffer.data(idx_op_fg + 136);

    auto tr_zzz_xxxz = pbuffer.data(idx_op_fg + 137);

    auto tr_zzz_xxyy = pbuffer.data(idx_op_fg + 138);

    auto tr_zzz_xxyz = pbuffer.data(idx_op_fg + 139);

    auto tr_zzz_xxzz = pbuffer.data(idx_op_fg + 140);

    auto tr_zzz_xyyy = pbuffer.data(idx_op_fg + 141);

    auto tr_zzz_xyyz = pbuffer.data(idx_op_fg + 142);

    auto tr_zzz_xyzz = pbuffer.data(idx_op_fg + 143);

    auto tr_zzz_xzzz = pbuffer.data(idx_op_fg + 144);

    auto tr_zzz_yyyy = pbuffer.data(idx_op_fg + 145);

    auto tr_zzz_yyyz = pbuffer.data(idx_op_fg + 146);

    auto tr_zzz_yyzz = pbuffer.data(idx_op_fg + 147);

    auto tr_zzz_yzzz = pbuffer.data(idx_op_fg + 148);

    auto tr_zzz_zzzz = pbuffer.data(idx_op_fg + 149);

    // Set up components of auxiliary buffer : GF

    auto tr_xxxx_xxx = pbuffer.data(idx_op_gf);

    auto tr_xxxx_xxy = pbuffer.data(idx_op_gf + 1);

    auto tr_xxxx_xxz = pbuffer.data(idx_op_gf + 2);

    auto tr_xxxx_xyy = pbuffer.data(idx_op_gf + 3);

    auto tr_xxxx_xyz = pbuffer.data(idx_op_gf + 4);

    auto tr_xxxx_xzz = pbuffer.data(idx_op_gf + 5);

    auto tr_xxxx_yyy = pbuffer.data(idx_op_gf + 6);

    auto tr_xxxx_yyz = pbuffer.data(idx_op_gf + 7);

    auto tr_xxxx_yzz = pbuffer.data(idx_op_gf + 8);

    auto tr_xxxx_zzz = pbuffer.data(idx_op_gf + 9);

    auto tr_xxxy_xxx = pbuffer.data(idx_op_gf + 10);

    auto tr_xxxy_xxy = pbuffer.data(idx_op_gf + 11);

    auto tr_xxxy_xxz = pbuffer.data(idx_op_gf + 12);

    auto tr_xxxy_xyy = pbuffer.data(idx_op_gf + 13);

    auto tr_xxxy_xyz = pbuffer.data(idx_op_gf + 14);

    auto tr_xxxy_xzz = pbuffer.data(idx_op_gf + 15);

    auto tr_xxxy_yyy = pbuffer.data(idx_op_gf + 16);

    auto tr_xxxy_yyz = pbuffer.data(idx_op_gf + 17);

    auto tr_xxxy_yzz = pbuffer.data(idx_op_gf + 18);

    auto tr_xxxy_zzz = pbuffer.data(idx_op_gf + 19);

    auto tr_xxxz_xxx = pbuffer.data(idx_op_gf + 20);

    auto tr_xxxz_xxy = pbuffer.data(idx_op_gf + 21);

    auto tr_xxxz_xxz = pbuffer.data(idx_op_gf + 22);

    auto tr_xxxz_xyy = pbuffer.data(idx_op_gf + 23);

    auto tr_xxxz_xyz = pbuffer.data(idx_op_gf + 24);

    auto tr_xxxz_xzz = pbuffer.data(idx_op_gf + 25);

    auto tr_xxxz_yyy = pbuffer.data(idx_op_gf + 26);

    auto tr_xxxz_yyz = pbuffer.data(idx_op_gf + 27);

    auto tr_xxxz_yzz = pbuffer.data(idx_op_gf + 28);

    auto tr_xxxz_zzz = pbuffer.data(idx_op_gf + 29);

    auto tr_xxyy_xxx = pbuffer.data(idx_op_gf + 30);

    auto tr_xxyy_xxy = pbuffer.data(idx_op_gf + 31);

    auto tr_xxyy_xxz = pbuffer.data(idx_op_gf + 32);

    auto tr_xxyy_xyy = pbuffer.data(idx_op_gf + 33);

    auto tr_xxyy_xyz = pbuffer.data(idx_op_gf + 34);

    auto tr_xxyy_xzz = pbuffer.data(idx_op_gf + 35);

    auto tr_xxyy_yyy = pbuffer.data(idx_op_gf + 36);

    auto tr_xxyy_yyz = pbuffer.data(idx_op_gf + 37);

    auto tr_xxyy_yzz = pbuffer.data(idx_op_gf + 38);

    auto tr_xxyy_zzz = pbuffer.data(idx_op_gf + 39);

    auto tr_xxyz_xxx = pbuffer.data(idx_op_gf + 40);

    auto tr_xxyz_xxy = pbuffer.data(idx_op_gf + 41);

    auto tr_xxyz_xxz = pbuffer.data(idx_op_gf + 42);

    auto tr_xxyz_xyy = pbuffer.data(idx_op_gf + 43);

    auto tr_xxyz_xyz = pbuffer.data(idx_op_gf + 44);

    auto tr_xxyz_xzz = pbuffer.data(idx_op_gf + 45);

    auto tr_xxyz_yyy = pbuffer.data(idx_op_gf + 46);

    auto tr_xxyz_yyz = pbuffer.data(idx_op_gf + 47);

    auto tr_xxyz_yzz = pbuffer.data(idx_op_gf + 48);

    auto tr_xxyz_zzz = pbuffer.data(idx_op_gf + 49);

    auto tr_xxzz_xxx = pbuffer.data(idx_op_gf + 50);

    auto tr_xxzz_xxy = pbuffer.data(idx_op_gf + 51);

    auto tr_xxzz_xxz = pbuffer.data(idx_op_gf + 52);

    auto tr_xxzz_xyy = pbuffer.data(idx_op_gf + 53);

    auto tr_xxzz_xyz = pbuffer.data(idx_op_gf + 54);

    auto tr_xxzz_xzz = pbuffer.data(idx_op_gf + 55);

    auto tr_xxzz_yyy = pbuffer.data(idx_op_gf + 56);

    auto tr_xxzz_yyz = pbuffer.data(idx_op_gf + 57);

    auto tr_xxzz_yzz = pbuffer.data(idx_op_gf + 58);

    auto tr_xxzz_zzz = pbuffer.data(idx_op_gf + 59);

    auto tr_xyyy_xxx = pbuffer.data(idx_op_gf + 60);

    auto tr_xyyy_xxy = pbuffer.data(idx_op_gf + 61);

    auto tr_xyyy_xxz = pbuffer.data(idx_op_gf + 62);

    auto tr_xyyy_xyy = pbuffer.data(idx_op_gf + 63);

    auto tr_xyyy_xyz = pbuffer.data(idx_op_gf + 64);

    auto tr_xyyy_xzz = pbuffer.data(idx_op_gf + 65);

    auto tr_xyyy_yyy = pbuffer.data(idx_op_gf + 66);

    auto tr_xyyy_yyz = pbuffer.data(idx_op_gf + 67);

    auto tr_xyyy_yzz = pbuffer.data(idx_op_gf + 68);

    auto tr_xyyy_zzz = pbuffer.data(idx_op_gf + 69);

    auto tr_xyyz_xxx = pbuffer.data(idx_op_gf + 70);

    auto tr_xyyz_xxy = pbuffer.data(idx_op_gf + 71);

    auto tr_xyyz_xxz = pbuffer.data(idx_op_gf + 72);

    auto tr_xyyz_xyy = pbuffer.data(idx_op_gf + 73);

    auto tr_xyyz_xyz = pbuffer.data(idx_op_gf + 74);

    auto tr_xyyz_xzz = pbuffer.data(idx_op_gf + 75);

    auto tr_xyyz_yyy = pbuffer.data(idx_op_gf + 76);

    auto tr_xyyz_yyz = pbuffer.data(idx_op_gf + 77);

    auto tr_xyyz_yzz = pbuffer.data(idx_op_gf + 78);

    auto tr_xyyz_zzz = pbuffer.data(idx_op_gf + 79);

    auto tr_xyzz_xxx = pbuffer.data(idx_op_gf + 80);

    auto tr_xyzz_xxy = pbuffer.data(idx_op_gf + 81);

    auto tr_xyzz_xxz = pbuffer.data(idx_op_gf + 82);

    auto tr_xyzz_xyy = pbuffer.data(idx_op_gf + 83);

    auto tr_xyzz_xyz = pbuffer.data(idx_op_gf + 84);

    auto tr_xyzz_xzz = pbuffer.data(idx_op_gf + 85);

    auto tr_xyzz_yyy = pbuffer.data(idx_op_gf + 86);

    auto tr_xyzz_yyz = pbuffer.data(idx_op_gf + 87);

    auto tr_xyzz_yzz = pbuffer.data(idx_op_gf + 88);

    auto tr_xyzz_zzz = pbuffer.data(idx_op_gf + 89);

    auto tr_xzzz_xxx = pbuffer.data(idx_op_gf + 90);

    auto tr_xzzz_xxy = pbuffer.data(idx_op_gf + 91);

    auto tr_xzzz_xxz = pbuffer.data(idx_op_gf + 92);

    auto tr_xzzz_xyy = pbuffer.data(idx_op_gf + 93);

    auto tr_xzzz_xyz = pbuffer.data(idx_op_gf + 94);

    auto tr_xzzz_xzz = pbuffer.data(idx_op_gf + 95);

    auto tr_xzzz_yyy = pbuffer.data(idx_op_gf + 96);

    auto tr_xzzz_yyz = pbuffer.data(idx_op_gf + 97);

    auto tr_xzzz_yzz = pbuffer.data(idx_op_gf + 98);

    auto tr_xzzz_zzz = pbuffer.data(idx_op_gf + 99);

    auto tr_yyyy_xxx = pbuffer.data(idx_op_gf + 100);

    auto tr_yyyy_xxy = pbuffer.data(idx_op_gf + 101);

    auto tr_yyyy_xxz = pbuffer.data(idx_op_gf + 102);

    auto tr_yyyy_xyy = pbuffer.data(idx_op_gf + 103);

    auto tr_yyyy_xyz = pbuffer.data(idx_op_gf + 104);

    auto tr_yyyy_xzz = pbuffer.data(idx_op_gf + 105);

    auto tr_yyyy_yyy = pbuffer.data(idx_op_gf + 106);

    auto tr_yyyy_yyz = pbuffer.data(idx_op_gf + 107);

    auto tr_yyyy_yzz = pbuffer.data(idx_op_gf + 108);

    auto tr_yyyy_zzz = pbuffer.data(idx_op_gf + 109);

    auto tr_yyyz_xxx = pbuffer.data(idx_op_gf + 110);

    auto tr_yyyz_xxy = pbuffer.data(idx_op_gf + 111);

    auto tr_yyyz_xxz = pbuffer.data(idx_op_gf + 112);

    auto tr_yyyz_xyy = pbuffer.data(idx_op_gf + 113);

    auto tr_yyyz_xyz = pbuffer.data(idx_op_gf + 114);

    auto tr_yyyz_xzz = pbuffer.data(idx_op_gf + 115);

    auto tr_yyyz_yyy = pbuffer.data(idx_op_gf + 116);

    auto tr_yyyz_yyz = pbuffer.data(idx_op_gf + 117);

    auto tr_yyyz_yzz = pbuffer.data(idx_op_gf + 118);

    auto tr_yyyz_zzz = pbuffer.data(idx_op_gf + 119);

    auto tr_yyzz_xxx = pbuffer.data(idx_op_gf + 120);

    auto tr_yyzz_xxy = pbuffer.data(idx_op_gf + 121);

    auto tr_yyzz_xxz = pbuffer.data(idx_op_gf + 122);

    auto tr_yyzz_xyy = pbuffer.data(idx_op_gf + 123);

    auto tr_yyzz_xyz = pbuffer.data(idx_op_gf + 124);

    auto tr_yyzz_xzz = pbuffer.data(idx_op_gf + 125);

    auto tr_yyzz_yyy = pbuffer.data(idx_op_gf + 126);

    auto tr_yyzz_yyz = pbuffer.data(idx_op_gf + 127);

    auto tr_yyzz_yzz = pbuffer.data(idx_op_gf + 128);

    auto tr_yyzz_zzz = pbuffer.data(idx_op_gf + 129);

    auto tr_yzzz_xxx = pbuffer.data(idx_op_gf + 130);

    auto tr_yzzz_xxy = pbuffer.data(idx_op_gf + 131);

    auto tr_yzzz_xxz = pbuffer.data(idx_op_gf + 132);

    auto tr_yzzz_xyy = pbuffer.data(idx_op_gf + 133);

    auto tr_yzzz_xyz = pbuffer.data(idx_op_gf + 134);

    auto tr_yzzz_xzz = pbuffer.data(idx_op_gf + 135);

    auto tr_yzzz_yyy = pbuffer.data(idx_op_gf + 136);

    auto tr_yzzz_yyz = pbuffer.data(idx_op_gf + 137);

    auto tr_yzzz_yzz = pbuffer.data(idx_op_gf + 138);

    auto tr_yzzz_zzz = pbuffer.data(idx_op_gf + 139);

    auto tr_zzzz_xxx = pbuffer.data(idx_op_gf + 140);

    auto tr_zzzz_xxy = pbuffer.data(idx_op_gf + 141);

    auto tr_zzzz_xxz = pbuffer.data(idx_op_gf + 142);

    auto tr_zzzz_xyy = pbuffer.data(idx_op_gf + 143);

    auto tr_zzzz_xyz = pbuffer.data(idx_op_gf + 144);

    auto tr_zzzz_xzz = pbuffer.data(idx_op_gf + 145);

    auto tr_zzzz_yyy = pbuffer.data(idx_op_gf + 146);

    auto tr_zzzz_yyz = pbuffer.data(idx_op_gf + 147);

    auto tr_zzzz_yzz = pbuffer.data(idx_op_gf + 148);

    auto tr_zzzz_zzz = pbuffer.data(idx_op_gf + 149);

    // Set up components of auxiliary buffer : GH

    auto tr_xxxx_xxxxx = pbuffer.data(idx_op_gh);

    auto tr_xxxx_xxxxy = pbuffer.data(idx_op_gh + 1);

    auto tr_xxxx_xxxxz = pbuffer.data(idx_op_gh + 2);

    auto tr_xxxx_xxxyy = pbuffer.data(idx_op_gh + 3);

    auto tr_xxxx_xxxyz = pbuffer.data(idx_op_gh + 4);

    auto tr_xxxx_xxxzz = pbuffer.data(idx_op_gh + 5);

    auto tr_xxxx_xxyyy = pbuffer.data(idx_op_gh + 6);

    auto tr_xxxx_xxyyz = pbuffer.data(idx_op_gh + 7);

    auto tr_xxxx_xxyzz = pbuffer.data(idx_op_gh + 8);

    auto tr_xxxx_xxzzz = pbuffer.data(idx_op_gh + 9);

    auto tr_xxxx_xyyyy = pbuffer.data(idx_op_gh + 10);

    auto tr_xxxx_xyyyz = pbuffer.data(idx_op_gh + 11);

    auto tr_xxxx_xyyzz = pbuffer.data(idx_op_gh + 12);

    auto tr_xxxx_xyzzz = pbuffer.data(idx_op_gh + 13);

    auto tr_xxxx_xzzzz = pbuffer.data(idx_op_gh + 14);

    auto tr_xxxx_yyyyy = pbuffer.data(idx_op_gh + 15);

    auto tr_xxxx_yyyyz = pbuffer.data(idx_op_gh + 16);

    auto tr_xxxx_yyyzz = pbuffer.data(idx_op_gh + 17);

    auto tr_xxxx_yyzzz = pbuffer.data(idx_op_gh + 18);

    auto tr_xxxx_yzzzz = pbuffer.data(idx_op_gh + 19);

    auto tr_xxxx_zzzzz = pbuffer.data(idx_op_gh + 20);

    auto tr_xxxy_xxxxx = pbuffer.data(idx_op_gh + 21);

    auto tr_xxxy_xxxxy = pbuffer.data(idx_op_gh + 22);

    auto tr_xxxy_xxxxz = pbuffer.data(idx_op_gh + 23);

    auto tr_xxxy_xxxyy = pbuffer.data(idx_op_gh + 24);

    auto tr_xxxy_xxxyz = pbuffer.data(idx_op_gh + 25);

    auto tr_xxxy_xxxzz = pbuffer.data(idx_op_gh + 26);

    auto tr_xxxy_xxyyy = pbuffer.data(idx_op_gh + 27);

    auto tr_xxxy_xxyyz = pbuffer.data(idx_op_gh + 28);

    auto tr_xxxy_xxyzz = pbuffer.data(idx_op_gh + 29);

    auto tr_xxxy_xxzzz = pbuffer.data(idx_op_gh + 30);

    auto tr_xxxy_xyyyy = pbuffer.data(idx_op_gh + 31);

    auto tr_xxxy_xyyyz = pbuffer.data(idx_op_gh + 32);

    auto tr_xxxy_xyyzz = pbuffer.data(idx_op_gh + 33);

    auto tr_xxxy_xyzzz = pbuffer.data(idx_op_gh + 34);

    auto tr_xxxy_xzzzz = pbuffer.data(idx_op_gh + 35);

    auto tr_xxxy_yyyyy = pbuffer.data(idx_op_gh + 36);

    auto tr_xxxy_yyyyz = pbuffer.data(idx_op_gh + 37);

    auto tr_xxxy_yyyzz = pbuffer.data(idx_op_gh + 38);

    auto tr_xxxy_yyzzz = pbuffer.data(idx_op_gh + 39);

    auto tr_xxxy_yzzzz = pbuffer.data(idx_op_gh + 40);

    auto tr_xxxy_zzzzz = pbuffer.data(idx_op_gh + 41);

    auto tr_xxxz_xxxxx = pbuffer.data(idx_op_gh + 42);

    auto tr_xxxz_xxxxy = pbuffer.data(idx_op_gh + 43);

    auto tr_xxxz_xxxxz = pbuffer.data(idx_op_gh + 44);

    auto tr_xxxz_xxxyy = pbuffer.data(idx_op_gh + 45);

    auto tr_xxxz_xxxyz = pbuffer.data(idx_op_gh + 46);

    auto tr_xxxz_xxxzz = pbuffer.data(idx_op_gh + 47);

    auto tr_xxxz_xxyyy = pbuffer.data(idx_op_gh + 48);

    auto tr_xxxz_xxyyz = pbuffer.data(idx_op_gh + 49);

    auto tr_xxxz_xxyzz = pbuffer.data(idx_op_gh + 50);

    auto tr_xxxz_xxzzz = pbuffer.data(idx_op_gh + 51);

    auto tr_xxxz_xyyyy = pbuffer.data(idx_op_gh + 52);

    auto tr_xxxz_xyyyz = pbuffer.data(idx_op_gh + 53);

    auto tr_xxxz_xyyzz = pbuffer.data(idx_op_gh + 54);

    auto tr_xxxz_xyzzz = pbuffer.data(idx_op_gh + 55);

    auto tr_xxxz_xzzzz = pbuffer.data(idx_op_gh + 56);

    auto tr_xxxz_yyyyy = pbuffer.data(idx_op_gh + 57);

    auto tr_xxxz_yyyyz = pbuffer.data(idx_op_gh + 58);

    auto tr_xxxz_yyyzz = pbuffer.data(idx_op_gh + 59);

    auto tr_xxxz_yyzzz = pbuffer.data(idx_op_gh + 60);

    auto tr_xxxz_yzzzz = pbuffer.data(idx_op_gh + 61);

    auto tr_xxxz_zzzzz = pbuffer.data(idx_op_gh + 62);

    auto tr_xxyy_xxxxx = pbuffer.data(idx_op_gh + 63);

    auto tr_xxyy_xxxxy = pbuffer.data(idx_op_gh + 64);

    auto tr_xxyy_xxxxz = pbuffer.data(idx_op_gh + 65);

    auto tr_xxyy_xxxyy = pbuffer.data(idx_op_gh + 66);

    auto tr_xxyy_xxxyz = pbuffer.data(idx_op_gh + 67);

    auto tr_xxyy_xxxzz = pbuffer.data(idx_op_gh + 68);

    auto tr_xxyy_xxyyy = pbuffer.data(idx_op_gh + 69);

    auto tr_xxyy_xxyyz = pbuffer.data(idx_op_gh + 70);

    auto tr_xxyy_xxyzz = pbuffer.data(idx_op_gh + 71);

    auto tr_xxyy_xxzzz = pbuffer.data(idx_op_gh + 72);

    auto tr_xxyy_xyyyy = pbuffer.data(idx_op_gh + 73);

    auto tr_xxyy_xyyyz = pbuffer.data(idx_op_gh + 74);

    auto tr_xxyy_xyyzz = pbuffer.data(idx_op_gh + 75);

    auto tr_xxyy_xyzzz = pbuffer.data(idx_op_gh + 76);

    auto tr_xxyy_xzzzz = pbuffer.data(idx_op_gh + 77);

    auto tr_xxyy_yyyyy = pbuffer.data(idx_op_gh + 78);

    auto tr_xxyy_yyyyz = pbuffer.data(idx_op_gh + 79);

    auto tr_xxyy_yyyzz = pbuffer.data(idx_op_gh + 80);

    auto tr_xxyy_yyzzz = pbuffer.data(idx_op_gh + 81);

    auto tr_xxyy_yzzzz = pbuffer.data(idx_op_gh + 82);

    auto tr_xxyy_zzzzz = pbuffer.data(idx_op_gh + 83);

    auto tr_xxyz_xxxxx = pbuffer.data(idx_op_gh + 84);

    auto tr_xxyz_xxxxy = pbuffer.data(idx_op_gh + 85);

    auto tr_xxyz_xxxxz = pbuffer.data(idx_op_gh + 86);

    auto tr_xxyz_xxxyy = pbuffer.data(idx_op_gh + 87);

    auto tr_xxyz_xxxyz = pbuffer.data(idx_op_gh + 88);

    auto tr_xxyz_xxxzz = pbuffer.data(idx_op_gh + 89);

    auto tr_xxyz_xxyyy = pbuffer.data(idx_op_gh + 90);

    auto tr_xxyz_xxyyz = pbuffer.data(idx_op_gh + 91);

    auto tr_xxyz_xxyzz = pbuffer.data(idx_op_gh + 92);

    auto tr_xxyz_xxzzz = pbuffer.data(idx_op_gh + 93);

    auto tr_xxyz_xyyyy = pbuffer.data(idx_op_gh + 94);

    auto tr_xxyz_xyyyz = pbuffer.data(idx_op_gh + 95);

    auto tr_xxyz_xyyzz = pbuffer.data(idx_op_gh + 96);

    auto tr_xxyz_xyzzz = pbuffer.data(idx_op_gh + 97);

    auto tr_xxyz_xzzzz = pbuffer.data(idx_op_gh + 98);

    auto tr_xxyz_yyyyy = pbuffer.data(idx_op_gh + 99);

    auto tr_xxyz_yyyyz = pbuffer.data(idx_op_gh + 100);

    auto tr_xxyz_yyyzz = pbuffer.data(idx_op_gh + 101);

    auto tr_xxyz_yyzzz = pbuffer.data(idx_op_gh + 102);

    auto tr_xxyz_yzzzz = pbuffer.data(idx_op_gh + 103);

    auto tr_xxyz_zzzzz = pbuffer.data(idx_op_gh + 104);

    auto tr_xxzz_xxxxx = pbuffer.data(idx_op_gh + 105);

    auto tr_xxzz_xxxxy = pbuffer.data(idx_op_gh + 106);

    auto tr_xxzz_xxxxz = pbuffer.data(idx_op_gh + 107);

    auto tr_xxzz_xxxyy = pbuffer.data(idx_op_gh + 108);

    auto tr_xxzz_xxxyz = pbuffer.data(idx_op_gh + 109);

    auto tr_xxzz_xxxzz = pbuffer.data(idx_op_gh + 110);

    auto tr_xxzz_xxyyy = pbuffer.data(idx_op_gh + 111);

    auto tr_xxzz_xxyyz = pbuffer.data(idx_op_gh + 112);

    auto tr_xxzz_xxyzz = pbuffer.data(idx_op_gh + 113);

    auto tr_xxzz_xxzzz = pbuffer.data(idx_op_gh + 114);

    auto tr_xxzz_xyyyy = pbuffer.data(idx_op_gh + 115);

    auto tr_xxzz_xyyyz = pbuffer.data(idx_op_gh + 116);

    auto tr_xxzz_xyyzz = pbuffer.data(idx_op_gh + 117);

    auto tr_xxzz_xyzzz = pbuffer.data(idx_op_gh + 118);

    auto tr_xxzz_xzzzz = pbuffer.data(idx_op_gh + 119);

    auto tr_xxzz_yyyyy = pbuffer.data(idx_op_gh + 120);

    auto tr_xxzz_yyyyz = pbuffer.data(idx_op_gh + 121);

    auto tr_xxzz_yyyzz = pbuffer.data(idx_op_gh + 122);

    auto tr_xxzz_yyzzz = pbuffer.data(idx_op_gh + 123);

    auto tr_xxzz_yzzzz = pbuffer.data(idx_op_gh + 124);

    auto tr_xxzz_zzzzz = pbuffer.data(idx_op_gh + 125);

    auto tr_xyyy_xxxxx = pbuffer.data(idx_op_gh + 126);

    auto tr_xyyy_xxxxy = pbuffer.data(idx_op_gh + 127);

    auto tr_xyyy_xxxxz = pbuffer.data(idx_op_gh + 128);

    auto tr_xyyy_xxxyy = pbuffer.data(idx_op_gh + 129);

    auto tr_xyyy_xxxyz = pbuffer.data(idx_op_gh + 130);

    auto tr_xyyy_xxxzz = pbuffer.data(idx_op_gh + 131);

    auto tr_xyyy_xxyyy = pbuffer.data(idx_op_gh + 132);

    auto tr_xyyy_xxyyz = pbuffer.data(idx_op_gh + 133);

    auto tr_xyyy_xxyzz = pbuffer.data(idx_op_gh + 134);

    auto tr_xyyy_xxzzz = pbuffer.data(idx_op_gh + 135);

    auto tr_xyyy_xyyyy = pbuffer.data(idx_op_gh + 136);

    auto tr_xyyy_xyyyz = pbuffer.data(idx_op_gh + 137);

    auto tr_xyyy_xyyzz = pbuffer.data(idx_op_gh + 138);

    auto tr_xyyy_xyzzz = pbuffer.data(idx_op_gh + 139);

    auto tr_xyyy_xzzzz = pbuffer.data(idx_op_gh + 140);

    auto tr_xyyy_yyyyy = pbuffer.data(idx_op_gh + 141);

    auto tr_xyyy_yyyyz = pbuffer.data(idx_op_gh + 142);

    auto tr_xyyy_yyyzz = pbuffer.data(idx_op_gh + 143);

    auto tr_xyyy_yyzzz = pbuffer.data(idx_op_gh + 144);

    auto tr_xyyy_yzzzz = pbuffer.data(idx_op_gh + 145);

    auto tr_xyyy_zzzzz = pbuffer.data(idx_op_gh + 146);

    auto tr_xyyz_xxxxx = pbuffer.data(idx_op_gh + 147);

    auto tr_xyyz_xxxxy = pbuffer.data(idx_op_gh + 148);

    auto tr_xyyz_xxxxz = pbuffer.data(idx_op_gh + 149);

    auto tr_xyyz_xxxyy = pbuffer.data(idx_op_gh + 150);

    auto tr_xyyz_xxxyz = pbuffer.data(idx_op_gh + 151);

    auto tr_xyyz_xxxzz = pbuffer.data(idx_op_gh + 152);

    auto tr_xyyz_xxyyy = pbuffer.data(idx_op_gh + 153);

    auto tr_xyyz_xxyyz = pbuffer.data(idx_op_gh + 154);

    auto tr_xyyz_xxyzz = pbuffer.data(idx_op_gh + 155);

    auto tr_xyyz_xxzzz = pbuffer.data(idx_op_gh + 156);

    auto tr_xyyz_xyyyy = pbuffer.data(idx_op_gh + 157);

    auto tr_xyyz_xyyyz = pbuffer.data(idx_op_gh + 158);

    auto tr_xyyz_xyyzz = pbuffer.data(idx_op_gh + 159);

    auto tr_xyyz_xyzzz = pbuffer.data(idx_op_gh + 160);

    auto tr_xyyz_xzzzz = pbuffer.data(idx_op_gh + 161);

    auto tr_xyyz_yyyyy = pbuffer.data(idx_op_gh + 162);

    auto tr_xyyz_yyyyz = pbuffer.data(idx_op_gh + 163);

    auto tr_xyyz_yyyzz = pbuffer.data(idx_op_gh + 164);

    auto tr_xyyz_yyzzz = pbuffer.data(idx_op_gh + 165);

    auto tr_xyyz_yzzzz = pbuffer.data(idx_op_gh + 166);

    auto tr_xyyz_zzzzz = pbuffer.data(idx_op_gh + 167);

    auto tr_xyzz_xxxxx = pbuffer.data(idx_op_gh + 168);

    auto tr_xyzz_xxxxy = pbuffer.data(idx_op_gh + 169);

    auto tr_xyzz_xxxxz = pbuffer.data(idx_op_gh + 170);

    auto tr_xyzz_xxxyy = pbuffer.data(idx_op_gh + 171);

    auto tr_xyzz_xxxyz = pbuffer.data(idx_op_gh + 172);

    auto tr_xyzz_xxxzz = pbuffer.data(idx_op_gh + 173);

    auto tr_xyzz_xxyyy = pbuffer.data(idx_op_gh + 174);

    auto tr_xyzz_xxyyz = pbuffer.data(idx_op_gh + 175);

    auto tr_xyzz_xxyzz = pbuffer.data(idx_op_gh + 176);

    auto tr_xyzz_xxzzz = pbuffer.data(idx_op_gh + 177);

    auto tr_xyzz_xyyyy = pbuffer.data(idx_op_gh + 178);

    auto tr_xyzz_xyyyz = pbuffer.data(idx_op_gh + 179);

    auto tr_xyzz_xyyzz = pbuffer.data(idx_op_gh + 180);

    auto tr_xyzz_xyzzz = pbuffer.data(idx_op_gh + 181);

    auto tr_xyzz_xzzzz = pbuffer.data(idx_op_gh + 182);

    auto tr_xyzz_yyyyy = pbuffer.data(idx_op_gh + 183);

    auto tr_xyzz_yyyyz = pbuffer.data(idx_op_gh + 184);

    auto tr_xyzz_yyyzz = pbuffer.data(idx_op_gh + 185);

    auto tr_xyzz_yyzzz = pbuffer.data(idx_op_gh + 186);

    auto tr_xyzz_yzzzz = pbuffer.data(idx_op_gh + 187);

    auto tr_xyzz_zzzzz = pbuffer.data(idx_op_gh + 188);

    auto tr_xzzz_xxxxx = pbuffer.data(idx_op_gh + 189);

    auto tr_xzzz_xxxxy = pbuffer.data(idx_op_gh + 190);

    auto tr_xzzz_xxxxz = pbuffer.data(idx_op_gh + 191);

    auto tr_xzzz_xxxyy = pbuffer.data(idx_op_gh + 192);

    auto tr_xzzz_xxxyz = pbuffer.data(idx_op_gh + 193);

    auto tr_xzzz_xxxzz = pbuffer.data(idx_op_gh + 194);

    auto tr_xzzz_xxyyy = pbuffer.data(idx_op_gh + 195);

    auto tr_xzzz_xxyyz = pbuffer.data(idx_op_gh + 196);

    auto tr_xzzz_xxyzz = pbuffer.data(idx_op_gh + 197);

    auto tr_xzzz_xxzzz = pbuffer.data(idx_op_gh + 198);

    auto tr_xzzz_xyyyy = pbuffer.data(idx_op_gh + 199);

    auto tr_xzzz_xyyyz = pbuffer.data(idx_op_gh + 200);

    auto tr_xzzz_xyyzz = pbuffer.data(idx_op_gh + 201);

    auto tr_xzzz_xyzzz = pbuffer.data(idx_op_gh + 202);

    auto tr_xzzz_xzzzz = pbuffer.data(idx_op_gh + 203);

    auto tr_xzzz_yyyyy = pbuffer.data(idx_op_gh + 204);

    auto tr_xzzz_yyyyz = pbuffer.data(idx_op_gh + 205);

    auto tr_xzzz_yyyzz = pbuffer.data(idx_op_gh + 206);

    auto tr_xzzz_yyzzz = pbuffer.data(idx_op_gh + 207);

    auto tr_xzzz_yzzzz = pbuffer.data(idx_op_gh + 208);

    auto tr_xzzz_zzzzz = pbuffer.data(idx_op_gh + 209);

    auto tr_yyyy_xxxxx = pbuffer.data(idx_op_gh + 210);

    auto tr_yyyy_xxxxy = pbuffer.data(idx_op_gh + 211);

    auto tr_yyyy_xxxxz = pbuffer.data(idx_op_gh + 212);

    auto tr_yyyy_xxxyy = pbuffer.data(idx_op_gh + 213);

    auto tr_yyyy_xxxyz = pbuffer.data(idx_op_gh + 214);

    auto tr_yyyy_xxxzz = pbuffer.data(idx_op_gh + 215);

    auto tr_yyyy_xxyyy = pbuffer.data(idx_op_gh + 216);

    auto tr_yyyy_xxyyz = pbuffer.data(idx_op_gh + 217);

    auto tr_yyyy_xxyzz = pbuffer.data(idx_op_gh + 218);

    auto tr_yyyy_xxzzz = pbuffer.data(idx_op_gh + 219);

    auto tr_yyyy_xyyyy = pbuffer.data(idx_op_gh + 220);

    auto tr_yyyy_xyyyz = pbuffer.data(idx_op_gh + 221);

    auto tr_yyyy_xyyzz = pbuffer.data(idx_op_gh + 222);

    auto tr_yyyy_xyzzz = pbuffer.data(idx_op_gh + 223);

    auto tr_yyyy_xzzzz = pbuffer.data(idx_op_gh + 224);

    auto tr_yyyy_yyyyy = pbuffer.data(idx_op_gh + 225);

    auto tr_yyyy_yyyyz = pbuffer.data(idx_op_gh + 226);

    auto tr_yyyy_yyyzz = pbuffer.data(idx_op_gh + 227);

    auto tr_yyyy_yyzzz = pbuffer.data(idx_op_gh + 228);

    auto tr_yyyy_yzzzz = pbuffer.data(idx_op_gh + 229);

    auto tr_yyyy_zzzzz = pbuffer.data(idx_op_gh + 230);

    auto tr_yyyz_xxxxx = pbuffer.data(idx_op_gh + 231);

    auto tr_yyyz_xxxxy = pbuffer.data(idx_op_gh + 232);

    auto tr_yyyz_xxxxz = pbuffer.data(idx_op_gh + 233);

    auto tr_yyyz_xxxyy = pbuffer.data(idx_op_gh + 234);

    auto tr_yyyz_xxxyz = pbuffer.data(idx_op_gh + 235);

    auto tr_yyyz_xxxzz = pbuffer.data(idx_op_gh + 236);

    auto tr_yyyz_xxyyy = pbuffer.data(idx_op_gh + 237);

    auto tr_yyyz_xxyyz = pbuffer.data(idx_op_gh + 238);

    auto tr_yyyz_xxyzz = pbuffer.data(idx_op_gh + 239);

    auto tr_yyyz_xxzzz = pbuffer.data(idx_op_gh + 240);

    auto tr_yyyz_xyyyy = pbuffer.data(idx_op_gh + 241);

    auto tr_yyyz_xyyyz = pbuffer.data(idx_op_gh + 242);

    auto tr_yyyz_xyyzz = pbuffer.data(idx_op_gh + 243);

    auto tr_yyyz_xyzzz = pbuffer.data(idx_op_gh + 244);

    auto tr_yyyz_xzzzz = pbuffer.data(idx_op_gh + 245);

    auto tr_yyyz_yyyyy = pbuffer.data(idx_op_gh + 246);

    auto tr_yyyz_yyyyz = pbuffer.data(idx_op_gh + 247);

    auto tr_yyyz_yyyzz = pbuffer.data(idx_op_gh + 248);

    auto tr_yyyz_yyzzz = pbuffer.data(idx_op_gh + 249);

    auto tr_yyyz_yzzzz = pbuffer.data(idx_op_gh + 250);

    auto tr_yyyz_zzzzz = pbuffer.data(idx_op_gh + 251);

    auto tr_yyzz_xxxxx = pbuffer.data(idx_op_gh + 252);

    auto tr_yyzz_xxxxy = pbuffer.data(idx_op_gh + 253);

    auto tr_yyzz_xxxxz = pbuffer.data(idx_op_gh + 254);

    auto tr_yyzz_xxxyy = pbuffer.data(idx_op_gh + 255);

    auto tr_yyzz_xxxyz = pbuffer.data(idx_op_gh + 256);

    auto tr_yyzz_xxxzz = pbuffer.data(idx_op_gh + 257);

    auto tr_yyzz_xxyyy = pbuffer.data(idx_op_gh + 258);

    auto tr_yyzz_xxyyz = pbuffer.data(idx_op_gh + 259);

    auto tr_yyzz_xxyzz = pbuffer.data(idx_op_gh + 260);

    auto tr_yyzz_xxzzz = pbuffer.data(idx_op_gh + 261);

    auto tr_yyzz_xyyyy = pbuffer.data(idx_op_gh + 262);

    auto tr_yyzz_xyyyz = pbuffer.data(idx_op_gh + 263);

    auto tr_yyzz_xyyzz = pbuffer.data(idx_op_gh + 264);

    auto tr_yyzz_xyzzz = pbuffer.data(idx_op_gh + 265);

    auto tr_yyzz_xzzzz = pbuffer.data(idx_op_gh + 266);

    auto tr_yyzz_yyyyy = pbuffer.data(idx_op_gh + 267);

    auto tr_yyzz_yyyyz = pbuffer.data(idx_op_gh + 268);

    auto tr_yyzz_yyyzz = pbuffer.data(idx_op_gh + 269);

    auto tr_yyzz_yyzzz = pbuffer.data(idx_op_gh + 270);

    auto tr_yyzz_yzzzz = pbuffer.data(idx_op_gh + 271);

    auto tr_yyzz_zzzzz = pbuffer.data(idx_op_gh + 272);

    auto tr_yzzz_xxxxx = pbuffer.data(idx_op_gh + 273);

    auto tr_yzzz_xxxxy = pbuffer.data(idx_op_gh + 274);

    auto tr_yzzz_xxxxz = pbuffer.data(idx_op_gh + 275);

    auto tr_yzzz_xxxyy = pbuffer.data(idx_op_gh + 276);

    auto tr_yzzz_xxxyz = pbuffer.data(idx_op_gh + 277);

    auto tr_yzzz_xxxzz = pbuffer.data(idx_op_gh + 278);

    auto tr_yzzz_xxyyy = pbuffer.data(idx_op_gh + 279);

    auto tr_yzzz_xxyyz = pbuffer.data(idx_op_gh + 280);

    auto tr_yzzz_xxyzz = pbuffer.data(idx_op_gh + 281);

    auto tr_yzzz_xxzzz = pbuffer.data(idx_op_gh + 282);

    auto tr_yzzz_xyyyy = pbuffer.data(idx_op_gh + 283);

    auto tr_yzzz_xyyyz = pbuffer.data(idx_op_gh + 284);

    auto tr_yzzz_xyyzz = pbuffer.data(idx_op_gh + 285);

    auto tr_yzzz_xyzzz = pbuffer.data(idx_op_gh + 286);

    auto tr_yzzz_xzzzz = pbuffer.data(idx_op_gh + 287);

    auto tr_yzzz_yyyyy = pbuffer.data(idx_op_gh + 288);

    auto tr_yzzz_yyyyz = pbuffer.data(idx_op_gh + 289);

    auto tr_yzzz_yyyzz = pbuffer.data(idx_op_gh + 290);

    auto tr_yzzz_yyzzz = pbuffer.data(idx_op_gh + 291);

    auto tr_yzzz_yzzzz = pbuffer.data(idx_op_gh + 292);

    auto tr_yzzz_zzzzz = pbuffer.data(idx_op_gh + 293);

    auto tr_zzzz_xxxxx = pbuffer.data(idx_op_gh + 294);

    auto tr_zzzz_xxxxy = pbuffer.data(idx_op_gh + 295);

    auto tr_zzzz_xxxxz = pbuffer.data(idx_op_gh + 296);

    auto tr_zzzz_xxxyy = pbuffer.data(idx_op_gh + 297);

    auto tr_zzzz_xxxyz = pbuffer.data(idx_op_gh + 298);

    auto tr_zzzz_xxxzz = pbuffer.data(idx_op_gh + 299);

    auto tr_zzzz_xxyyy = pbuffer.data(idx_op_gh + 300);

    auto tr_zzzz_xxyyz = pbuffer.data(idx_op_gh + 301);

    auto tr_zzzz_xxyzz = pbuffer.data(idx_op_gh + 302);

    auto tr_zzzz_xxzzz = pbuffer.data(idx_op_gh + 303);

    auto tr_zzzz_xyyyy = pbuffer.data(idx_op_gh + 304);

    auto tr_zzzz_xyyyz = pbuffer.data(idx_op_gh + 305);

    auto tr_zzzz_xyyzz = pbuffer.data(idx_op_gh + 306);

    auto tr_zzzz_xyzzz = pbuffer.data(idx_op_gh + 307);

    auto tr_zzzz_xzzzz = pbuffer.data(idx_op_gh + 308);

    auto tr_zzzz_yyyyy = pbuffer.data(idx_op_gh + 309);

    auto tr_zzzz_yyyyz = pbuffer.data(idx_op_gh + 310);

    auto tr_zzzz_yyyzz = pbuffer.data(idx_op_gh + 311);

    auto tr_zzzz_yyzzz = pbuffer.data(idx_op_gh + 312);

    auto tr_zzzz_yzzzz = pbuffer.data(idx_op_gh + 313);

    auto tr_zzzz_zzzzz = pbuffer.data(idx_op_gh + 314);

    // Set up components of auxiliary buffer : HG

    auto tr_xxxxx_xxxx = pbuffer.data(idx_op_hg);

    auto tr_xxxxx_xxxy = pbuffer.data(idx_op_hg + 1);

    auto tr_xxxxx_xxxz = pbuffer.data(idx_op_hg + 2);

    auto tr_xxxxx_xxyy = pbuffer.data(idx_op_hg + 3);

    auto tr_xxxxx_xxyz = pbuffer.data(idx_op_hg + 4);

    auto tr_xxxxx_xxzz = pbuffer.data(idx_op_hg + 5);

    auto tr_xxxxx_xyyy = pbuffer.data(idx_op_hg + 6);

    auto tr_xxxxx_xyyz = pbuffer.data(idx_op_hg + 7);

    auto tr_xxxxx_xyzz = pbuffer.data(idx_op_hg + 8);

    auto tr_xxxxx_xzzz = pbuffer.data(idx_op_hg + 9);

    auto tr_xxxxx_yyyy = pbuffer.data(idx_op_hg + 10);

    auto tr_xxxxx_yyyz = pbuffer.data(idx_op_hg + 11);

    auto tr_xxxxx_yyzz = pbuffer.data(idx_op_hg + 12);

    auto tr_xxxxx_yzzz = pbuffer.data(idx_op_hg + 13);

    auto tr_xxxxx_zzzz = pbuffer.data(idx_op_hg + 14);

    auto tr_xxxxy_xxxx = pbuffer.data(idx_op_hg + 15);

    auto tr_xxxxy_xxxy = pbuffer.data(idx_op_hg + 16);

    auto tr_xxxxy_xxxz = pbuffer.data(idx_op_hg + 17);

    auto tr_xxxxy_xxyy = pbuffer.data(idx_op_hg + 18);

    auto tr_xxxxy_xxyz = pbuffer.data(idx_op_hg + 19);

    auto tr_xxxxy_xxzz = pbuffer.data(idx_op_hg + 20);

    auto tr_xxxxy_xyyy = pbuffer.data(idx_op_hg + 21);

    auto tr_xxxxy_xyyz = pbuffer.data(idx_op_hg + 22);

    auto tr_xxxxy_xyzz = pbuffer.data(idx_op_hg + 23);

    auto tr_xxxxy_xzzz = pbuffer.data(idx_op_hg + 24);

    auto tr_xxxxy_yyyy = pbuffer.data(idx_op_hg + 25);

    auto tr_xxxxy_yyyz = pbuffer.data(idx_op_hg + 26);

    auto tr_xxxxy_yyzz = pbuffer.data(idx_op_hg + 27);

    auto tr_xxxxy_yzzz = pbuffer.data(idx_op_hg + 28);

    auto tr_xxxxy_zzzz = pbuffer.data(idx_op_hg + 29);

    auto tr_xxxxz_xxxx = pbuffer.data(idx_op_hg + 30);

    auto tr_xxxxz_xxxy = pbuffer.data(idx_op_hg + 31);

    auto tr_xxxxz_xxxz = pbuffer.data(idx_op_hg + 32);

    auto tr_xxxxz_xxyy = pbuffer.data(idx_op_hg + 33);

    auto tr_xxxxz_xxyz = pbuffer.data(idx_op_hg + 34);

    auto tr_xxxxz_xxzz = pbuffer.data(idx_op_hg + 35);

    auto tr_xxxxz_xyyy = pbuffer.data(idx_op_hg + 36);

    auto tr_xxxxz_xyyz = pbuffer.data(idx_op_hg + 37);

    auto tr_xxxxz_xyzz = pbuffer.data(idx_op_hg + 38);

    auto tr_xxxxz_xzzz = pbuffer.data(idx_op_hg + 39);

    auto tr_xxxxz_yyyy = pbuffer.data(idx_op_hg + 40);

    auto tr_xxxxz_yyyz = pbuffer.data(idx_op_hg + 41);

    auto tr_xxxxz_yyzz = pbuffer.data(idx_op_hg + 42);

    auto tr_xxxxz_yzzz = pbuffer.data(idx_op_hg + 43);

    auto tr_xxxxz_zzzz = pbuffer.data(idx_op_hg + 44);

    auto tr_xxxyy_xxxx = pbuffer.data(idx_op_hg + 45);

    auto tr_xxxyy_xxxy = pbuffer.data(idx_op_hg + 46);

    auto tr_xxxyy_xxxz = pbuffer.data(idx_op_hg + 47);

    auto tr_xxxyy_xxyy = pbuffer.data(idx_op_hg + 48);

    auto tr_xxxyy_xxyz = pbuffer.data(idx_op_hg + 49);

    auto tr_xxxyy_xxzz = pbuffer.data(idx_op_hg + 50);

    auto tr_xxxyy_xyyy = pbuffer.data(idx_op_hg + 51);

    auto tr_xxxyy_xyyz = pbuffer.data(idx_op_hg + 52);

    auto tr_xxxyy_xyzz = pbuffer.data(idx_op_hg + 53);

    auto tr_xxxyy_xzzz = pbuffer.data(idx_op_hg + 54);

    auto tr_xxxyy_yyyy = pbuffer.data(idx_op_hg + 55);

    auto tr_xxxyy_yyyz = pbuffer.data(idx_op_hg + 56);

    auto tr_xxxyy_yyzz = pbuffer.data(idx_op_hg + 57);

    auto tr_xxxyy_yzzz = pbuffer.data(idx_op_hg + 58);

    auto tr_xxxyy_zzzz = pbuffer.data(idx_op_hg + 59);

    auto tr_xxxyz_xxxx = pbuffer.data(idx_op_hg + 60);

    auto tr_xxxyz_xxxy = pbuffer.data(idx_op_hg + 61);

    auto tr_xxxyz_xxxz = pbuffer.data(idx_op_hg + 62);

    auto tr_xxxyz_xxyy = pbuffer.data(idx_op_hg + 63);

    auto tr_xxxyz_xxyz = pbuffer.data(idx_op_hg + 64);

    auto tr_xxxyz_xxzz = pbuffer.data(idx_op_hg + 65);

    auto tr_xxxyz_xyyy = pbuffer.data(idx_op_hg + 66);

    auto tr_xxxyz_xyyz = pbuffer.data(idx_op_hg + 67);

    auto tr_xxxyz_xyzz = pbuffer.data(idx_op_hg + 68);

    auto tr_xxxyz_xzzz = pbuffer.data(idx_op_hg + 69);

    auto tr_xxxyz_yyyy = pbuffer.data(idx_op_hg + 70);

    auto tr_xxxyz_yyyz = pbuffer.data(idx_op_hg + 71);

    auto tr_xxxyz_yyzz = pbuffer.data(idx_op_hg + 72);

    auto tr_xxxyz_yzzz = pbuffer.data(idx_op_hg + 73);

    auto tr_xxxyz_zzzz = pbuffer.data(idx_op_hg + 74);

    auto tr_xxxzz_xxxx = pbuffer.data(idx_op_hg + 75);

    auto tr_xxxzz_xxxy = pbuffer.data(idx_op_hg + 76);

    auto tr_xxxzz_xxxz = pbuffer.data(idx_op_hg + 77);

    auto tr_xxxzz_xxyy = pbuffer.data(idx_op_hg + 78);

    auto tr_xxxzz_xxyz = pbuffer.data(idx_op_hg + 79);

    auto tr_xxxzz_xxzz = pbuffer.data(idx_op_hg + 80);

    auto tr_xxxzz_xyyy = pbuffer.data(idx_op_hg + 81);

    auto tr_xxxzz_xyyz = pbuffer.data(idx_op_hg + 82);

    auto tr_xxxzz_xyzz = pbuffer.data(idx_op_hg + 83);

    auto tr_xxxzz_xzzz = pbuffer.data(idx_op_hg + 84);

    auto tr_xxxzz_yyyy = pbuffer.data(idx_op_hg + 85);

    auto tr_xxxzz_yyyz = pbuffer.data(idx_op_hg + 86);

    auto tr_xxxzz_yyzz = pbuffer.data(idx_op_hg + 87);

    auto tr_xxxzz_yzzz = pbuffer.data(idx_op_hg + 88);

    auto tr_xxxzz_zzzz = pbuffer.data(idx_op_hg + 89);

    auto tr_xxyyy_xxxx = pbuffer.data(idx_op_hg + 90);

    auto tr_xxyyy_xxxy = pbuffer.data(idx_op_hg + 91);

    auto tr_xxyyy_xxxz = pbuffer.data(idx_op_hg + 92);

    auto tr_xxyyy_xxyy = pbuffer.data(idx_op_hg + 93);

    auto tr_xxyyy_xxyz = pbuffer.data(idx_op_hg + 94);

    auto tr_xxyyy_xxzz = pbuffer.data(idx_op_hg + 95);

    auto tr_xxyyy_xyyy = pbuffer.data(idx_op_hg + 96);

    auto tr_xxyyy_xyyz = pbuffer.data(idx_op_hg + 97);

    auto tr_xxyyy_xyzz = pbuffer.data(idx_op_hg + 98);

    auto tr_xxyyy_xzzz = pbuffer.data(idx_op_hg + 99);

    auto tr_xxyyy_yyyy = pbuffer.data(idx_op_hg + 100);

    auto tr_xxyyy_yyyz = pbuffer.data(idx_op_hg + 101);

    auto tr_xxyyy_yyzz = pbuffer.data(idx_op_hg + 102);

    auto tr_xxyyy_yzzz = pbuffer.data(idx_op_hg + 103);

    auto tr_xxyyy_zzzz = pbuffer.data(idx_op_hg + 104);

    auto tr_xxyyz_xxxx = pbuffer.data(idx_op_hg + 105);

    auto tr_xxyyz_xxxy = pbuffer.data(idx_op_hg + 106);

    auto tr_xxyyz_xxxz = pbuffer.data(idx_op_hg + 107);

    auto tr_xxyyz_xxyy = pbuffer.data(idx_op_hg + 108);

    auto tr_xxyyz_xxyz = pbuffer.data(idx_op_hg + 109);

    auto tr_xxyyz_xxzz = pbuffer.data(idx_op_hg + 110);

    auto tr_xxyyz_xyyy = pbuffer.data(idx_op_hg + 111);

    auto tr_xxyyz_xyyz = pbuffer.data(idx_op_hg + 112);

    auto tr_xxyyz_xyzz = pbuffer.data(idx_op_hg + 113);

    auto tr_xxyyz_xzzz = pbuffer.data(idx_op_hg + 114);

    auto tr_xxyyz_yyyy = pbuffer.data(idx_op_hg + 115);

    auto tr_xxyyz_yyyz = pbuffer.data(idx_op_hg + 116);

    auto tr_xxyyz_yyzz = pbuffer.data(idx_op_hg + 117);

    auto tr_xxyyz_yzzz = pbuffer.data(idx_op_hg + 118);

    auto tr_xxyyz_zzzz = pbuffer.data(idx_op_hg + 119);

    auto tr_xxyzz_xxxx = pbuffer.data(idx_op_hg + 120);

    auto tr_xxyzz_xxxy = pbuffer.data(idx_op_hg + 121);

    auto tr_xxyzz_xxxz = pbuffer.data(idx_op_hg + 122);

    auto tr_xxyzz_xxyy = pbuffer.data(idx_op_hg + 123);

    auto tr_xxyzz_xxyz = pbuffer.data(idx_op_hg + 124);

    auto tr_xxyzz_xxzz = pbuffer.data(idx_op_hg + 125);

    auto tr_xxyzz_xyyy = pbuffer.data(idx_op_hg + 126);

    auto tr_xxyzz_xyyz = pbuffer.data(idx_op_hg + 127);

    auto tr_xxyzz_xyzz = pbuffer.data(idx_op_hg + 128);

    auto tr_xxyzz_xzzz = pbuffer.data(idx_op_hg + 129);

    auto tr_xxyzz_yyyy = pbuffer.data(idx_op_hg + 130);

    auto tr_xxyzz_yyyz = pbuffer.data(idx_op_hg + 131);

    auto tr_xxyzz_yyzz = pbuffer.data(idx_op_hg + 132);

    auto tr_xxyzz_yzzz = pbuffer.data(idx_op_hg + 133);

    auto tr_xxyzz_zzzz = pbuffer.data(idx_op_hg + 134);

    auto tr_xxzzz_xxxx = pbuffer.data(idx_op_hg + 135);

    auto tr_xxzzz_xxxy = pbuffer.data(idx_op_hg + 136);

    auto tr_xxzzz_xxxz = pbuffer.data(idx_op_hg + 137);

    auto tr_xxzzz_xxyy = pbuffer.data(idx_op_hg + 138);

    auto tr_xxzzz_xxyz = pbuffer.data(idx_op_hg + 139);

    auto tr_xxzzz_xxzz = pbuffer.data(idx_op_hg + 140);

    auto tr_xxzzz_xyyy = pbuffer.data(idx_op_hg + 141);

    auto tr_xxzzz_xyyz = pbuffer.data(idx_op_hg + 142);

    auto tr_xxzzz_xyzz = pbuffer.data(idx_op_hg + 143);

    auto tr_xxzzz_xzzz = pbuffer.data(idx_op_hg + 144);

    auto tr_xxzzz_yyyy = pbuffer.data(idx_op_hg + 145);

    auto tr_xxzzz_yyyz = pbuffer.data(idx_op_hg + 146);

    auto tr_xxzzz_yyzz = pbuffer.data(idx_op_hg + 147);

    auto tr_xxzzz_yzzz = pbuffer.data(idx_op_hg + 148);

    auto tr_xxzzz_zzzz = pbuffer.data(idx_op_hg + 149);

    auto tr_xyyyy_xxxx = pbuffer.data(idx_op_hg + 150);

    auto tr_xyyyy_xxxy = pbuffer.data(idx_op_hg + 151);

    auto tr_xyyyy_xxxz = pbuffer.data(idx_op_hg + 152);

    auto tr_xyyyy_xxyy = pbuffer.data(idx_op_hg + 153);

    auto tr_xyyyy_xxyz = pbuffer.data(idx_op_hg + 154);

    auto tr_xyyyy_xxzz = pbuffer.data(idx_op_hg + 155);

    auto tr_xyyyy_xyyy = pbuffer.data(idx_op_hg + 156);

    auto tr_xyyyy_xyyz = pbuffer.data(idx_op_hg + 157);

    auto tr_xyyyy_xyzz = pbuffer.data(idx_op_hg + 158);

    auto tr_xyyyy_xzzz = pbuffer.data(idx_op_hg + 159);

    auto tr_xyyyy_yyyy = pbuffer.data(idx_op_hg + 160);

    auto tr_xyyyy_yyyz = pbuffer.data(idx_op_hg + 161);

    auto tr_xyyyy_yyzz = pbuffer.data(idx_op_hg + 162);

    auto tr_xyyyy_yzzz = pbuffer.data(idx_op_hg + 163);

    auto tr_xyyyy_zzzz = pbuffer.data(idx_op_hg + 164);

    auto tr_xyyyz_xxxx = pbuffer.data(idx_op_hg + 165);

    auto tr_xyyyz_xxxy = pbuffer.data(idx_op_hg + 166);

    auto tr_xyyyz_xxxz = pbuffer.data(idx_op_hg + 167);

    auto tr_xyyyz_xxyy = pbuffer.data(idx_op_hg + 168);

    auto tr_xyyyz_xxyz = pbuffer.data(idx_op_hg + 169);

    auto tr_xyyyz_xxzz = pbuffer.data(idx_op_hg + 170);

    auto tr_xyyyz_xyyy = pbuffer.data(idx_op_hg + 171);

    auto tr_xyyyz_xyyz = pbuffer.data(idx_op_hg + 172);

    auto tr_xyyyz_xyzz = pbuffer.data(idx_op_hg + 173);

    auto tr_xyyyz_xzzz = pbuffer.data(idx_op_hg + 174);

    auto tr_xyyyz_yyyy = pbuffer.data(idx_op_hg + 175);

    auto tr_xyyyz_yyyz = pbuffer.data(idx_op_hg + 176);

    auto tr_xyyyz_yyzz = pbuffer.data(idx_op_hg + 177);

    auto tr_xyyyz_yzzz = pbuffer.data(idx_op_hg + 178);

    auto tr_xyyyz_zzzz = pbuffer.data(idx_op_hg + 179);

    auto tr_xyyzz_xxxx = pbuffer.data(idx_op_hg + 180);

    auto tr_xyyzz_xxxy = pbuffer.data(idx_op_hg + 181);

    auto tr_xyyzz_xxxz = pbuffer.data(idx_op_hg + 182);

    auto tr_xyyzz_xxyy = pbuffer.data(idx_op_hg + 183);

    auto tr_xyyzz_xxyz = pbuffer.data(idx_op_hg + 184);

    auto tr_xyyzz_xxzz = pbuffer.data(idx_op_hg + 185);

    auto tr_xyyzz_xyyy = pbuffer.data(idx_op_hg + 186);

    auto tr_xyyzz_xyyz = pbuffer.data(idx_op_hg + 187);

    auto tr_xyyzz_xyzz = pbuffer.data(idx_op_hg + 188);

    auto tr_xyyzz_xzzz = pbuffer.data(idx_op_hg + 189);

    auto tr_xyyzz_yyyy = pbuffer.data(idx_op_hg + 190);

    auto tr_xyyzz_yyyz = pbuffer.data(idx_op_hg + 191);

    auto tr_xyyzz_yyzz = pbuffer.data(idx_op_hg + 192);

    auto tr_xyyzz_yzzz = pbuffer.data(idx_op_hg + 193);

    auto tr_xyyzz_zzzz = pbuffer.data(idx_op_hg + 194);

    auto tr_xyzzz_xxxx = pbuffer.data(idx_op_hg + 195);

    auto tr_xyzzz_xxxy = pbuffer.data(idx_op_hg + 196);

    auto tr_xyzzz_xxxz = pbuffer.data(idx_op_hg + 197);

    auto tr_xyzzz_xxyy = pbuffer.data(idx_op_hg + 198);

    auto tr_xyzzz_xxyz = pbuffer.data(idx_op_hg + 199);

    auto tr_xyzzz_xxzz = pbuffer.data(idx_op_hg + 200);

    auto tr_xyzzz_xyyy = pbuffer.data(idx_op_hg + 201);

    auto tr_xyzzz_xyyz = pbuffer.data(idx_op_hg + 202);

    auto tr_xyzzz_xyzz = pbuffer.data(idx_op_hg + 203);

    auto tr_xyzzz_xzzz = pbuffer.data(idx_op_hg + 204);

    auto tr_xyzzz_yyyy = pbuffer.data(idx_op_hg + 205);

    auto tr_xyzzz_yyyz = pbuffer.data(idx_op_hg + 206);

    auto tr_xyzzz_yyzz = pbuffer.data(idx_op_hg + 207);

    auto tr_xyzzz_yzzz = pbuffer.data(idx_op_hg + 208);

    auto tr_xyzzz_zzzz = pbuffer.data(idx_op_hg + 209);

    auto tr_xzzzz_xxxx = pbuffer.data(idx_op_hg + 210);

    auto tr_xzzzz_xxxy = pbuffer.data(idx_op_hg + 211);

    auto tr_xzzzz_xxxz = pbuffer.data(idx_op_hg + 212);

    auto tr_xzzzz_xxyy = pbuffer.data(idx_op_hg + 213);

    auto tr_xzzzz_xxyz = pbuffer.data(idx_op_hg + 214);

    auto tr_xzzzz_xxzz = pbuffer.data(idx_op_hg + 215);

    auto tr_xzzzz_xyyy = pbuffer.data(idx_op_hg + 216);

    auto tr_xzzzz_xyyz = pbuffer.data(idx_op_hg + 217);

    auto tr_xzzzz_xyzz = pbuffer.data(idx_op_hg + 218);

    auto tr_xzzzz_xzzz = pbuffer.data(idx_op_hg + 219);

    auto tr_xzzzz_yyyy = pbuffer.data(idx_op_hg + 220);

    auto tr_xzzzz_yyyz = pbuffer.data(idx_op_hg + 221);

    auto tr_xzzzz_yyzz = pbuffer.data(idx_op_hg + 222);

    auto tr_xzzzz_yzzz = pbuffer.data(idx_op_hg + 223);

    auto tr_xzzzz_zzzz = pbuffer.data(idx_op_hg + 224);

    auto tr_yyyyy_xxxx = pbuffer.data(idx_op_hg + 225);

    auto tr_yyyyy_xxxy = pbuffer.data(idx_op_hg + 226);

    auto tr_yyyyy_xxxz = pbuffer.data(idx_op_hg + 227);

    auto tr_yyyyy_xxyy = pbuffer.data(idx_op_hg + 228);

    auto tr_yyyyy_xxyz = pbuffer.data(idx_op_hg + 229);

    auto tr_yyyyy_xxzz = pbuffer.data(idx_op_hg + 230);

    auto tr_yyyyy_xyyy = pbuffer.data(idx_op_hg + 231);

    auto tr_yyyyy_xyyz = pbuffer.data(idx_op_hg + 232);

    auto tr_yyyyy_xyzz = pbuffer.data(idx_op_hg + 233);

    auto tr_yyyyy_xzzz = pbuffer.data(idx_op_hg + 234);

    auto tr_yyyyy_yyyy = pbuffer.data(idx_op_hg + 235);

    auto tr_yyyyy_yyyz = pbuffer.data(idx_op_hg + 236);

    auto tr_yyyyy_yyzz = pbuffer.data(idx_op_hg + 237);

    auto tr_yyyyy_yzzz = pbuffer.data(idx_op_hg + 238);

    auto tr_yyyyy_zzzz = pbuffer.data(idx_op_hg + 239);

    auto tr_yyyyz_xxxx = pbuffer.data(idx_op_hg + 240);

    auto tr_yyyyz_xxxy = pbuffer.data(idx_op_hg + 241);

    auto tr_yyyyz_xxxz = pbuffer.data(idx_op_hg + 242);

    auto tr_yyyyz_xxyy = pbuffer.data(idx_op_hg + 243);

    auto tr_yyyyz_xxyz = pbuffer.data(idx_op_hg + 244);

    auto tr_yyyyz_xxzz = pbuffer.data(idx_op_hg + 245);

    auto tr_yyyyz_xyyy = pbuffer.data(idx_op_hg + 246);

    auto tr_yyyyz_xyyz = pbuffer.data(idx_op_hg + 247);

    auto tr_yyyyz_xyzz = pbuffer.data(idx_op_hg + 248);

    auto tr_yyyyz_xzzz = pbuffer.data(idx_op_hg + 249);

    auto tr_yyyyz_yyyy = pbuffer.data(idx_op_hg + 250);

    auto tr_yyyyz_yyyz = pbuffer.data(idx_op_hg + 251);

    auto tr_yyyyz_yyzz = pbuffer.data(idx_op_hg + 252);

    auto tr_yyyyz_yzzz = pbuffer.data(idx_op_hg + 253);

    auto tr_yyyyz_zzzz = pbuffer.data(idx_op_hg + 254);

    auto tr_yyyzz_xxxx = pbuffer.data(idx_op_hg + 255);

    auto tr_yyyzz_xxxy = pbuffer.data(idx_op_hg + 256);

    auto tr_yyyzz_xxxz = pbuffer.data(idx_op_hg + 257);

    auto tr_yyyzz_xxyy = pbuffer.data(idx_op_hg + 258);

    auto tr_yyyzz_xxyz = pbuffer.data(idx_op_hg + 259);

    auto tr_yyyzz_xxzz = pbuffer.data(idx_op_hg + 260);

    auto tr_yyyzz_xyyy = pbuffer.data(idx_op_hg + 261);

    auto tr_yyyzz_xyyz = pbuffer.data(idx_op_hg + 262);

    auto tr_yyyzz_xyzz = pbuffer.data(idx_op_hg + 263);

    auto tr_yyyzz_xzzz = pbuffer.data(idx_op_hg + 264);

    auto tr_yyyzz_yyyy = pbuffer.data(idx_op_hg + 265);

    auto tr_yyyzz_yyyz = pbuffer.data(idx_op_hg + 266);

    auto tr_yyyzz_yyzz = pbuffer.data(idx_op_hg + 267);

    auto tr_yyyzz_yzzz = pbuffer.data(idx_op_hg + 268);

    auto tr_yyyzz_zzzz = pbuffer.data(idx_op_hg + 269);

    auto tr_yyzzz_xxxx = pbuffer.data(idx_op_hg + 270);

    auto tr_yyzzz_xxxy = pbuffer.data(idx_op_hg + 271);

    auto tr_yyzzz_xxxz = pbuffer.data(idx_op_hg + 272);

    auto tr_yyzzz_xxyy = pbuffer.data(idx_op_hg + 273);

    auto tr_yyzzz_xxyz = pbuffer.data(idx_op_hg + 274);

    auto tr_yyzzz_xxzz = pbuffer.data(idx_op_hg + 275);

    auto tr_yyzzz_xyyy = pbuffer.data(idx_op_hg + 276);

    auto tr_yyzzz_xyyz = pbuffer.data(idx_op_hg + 277);

    auto tr_yyzzz_xyzz = pbuffer.data(idx_op_hg + 278);

    auto tr_yyzzz_xzzz = pbuffer.data(idx_op_hg + 279);

    auto tr_yyzzz_yyyy = pbuffer.data(idx_op_hg + 280);

    auto tr_yyzzz_yyyz = pbuffer.data(idx_op_hg + 281);

    auto tr_yyzzz_yyzz = pbuffer.data(idx_op_hg + 282);

    auto tr_yyzzz_yzzz = pbuffer.data(idx_op_hg + 283);

    auto tr_yyzzz_zzzz = pbuffer.data(idx_op_hg + 284);

    auto tr_yzzzz_xxxx = pbuffer.data(idx_op_hg + 285);

    auto tr_yzzzz_xxxy = pbuffer.data(idx_op_hg + 286);

    auto tr_yzzzz_xxxz = pbuffer.data(idx_op_hg + 287);

    auto tr_yzzzz_xxyy = pbuffer.data(idx_op_hg + 288);

    auto tr_yzzzz_xxyz = pbuffer.data(idx_op_hg + 289);

    auto tr_yzzzz_xxzz = pbuffer.data(idx_op_hg + 290);

    auto tr_yzzzz_xyyy = pbuffer.data(idx_op_hg + 291);

    auto tr_yzzzz_xyyz = pbuffer.data(idx_op_hg + 292);

    auto tr_yzzzz_xyzz = pbuffer.data(idx_op_hg + 293);

    auto tr_yzzzz_xzzz = pbuffer.data(idx_op_hg + 294);

    auto tr_yzzzz_yyyy = pbuffer.data(idx_op_hg + 295);

    auto tr_yzzzz_yyyz = pbuffer.data(idx_op_hg + 296);

    auto tr_yzzzz_yyzz = pbuffer.data(idx_op_hg + 297);

    auto tr_yzzzz_yzzz = pbuffer.data(idx_op_hg + 298);

    auto tr_yzzzz_zzzz = pbuffer.data(idx_op_hg + 299);

    auto tr_zzzzz_xxxx = pbuffer.data(idx_op_hg + 300);

    auto tr_zzzzz_xxxy = pbuffer.data(idx_op_hg + 301);

    auto tr_zzzzz_xxxz = pbuffer.data(idx_op_hg + 302);

    auto tr_zzzzz_xxyy = pbuffer.data(idx_op_hg + 303);

    auto tr_zzzzz_xxyz = pbuffer.data(idx_op_hg + 304);

    auto tr_zzzzz_xxzz = pbuffer.data(idx_op_hg + 305);

    auto tr_zzzzz_xyyy = pbuffer.data(idx_op_hg + 306);

    auto tr_zzzzz_xyyz = pbuffer.data(idx_op_hg + 307);

    auto tr_zzzzz_xyzz = pbuffer.data(idx_op_hg + 308);

    auto tr_zzzzz_xzzz = pbuffer.data(idx_op_hg + 309);

    auto tr_zzzzz_yyyy = pbuffer.data(idx_op_hg + 310);

    auto tr_zzzzz_yyyz = pbuffer.data(idx_op_hg + 311);

    auto tr_zzzzz_yyzz = pbuffer.data(idx_op_hg + 312);

    auto tr_zzzzz_yzzz = pbuffer.data(idx_op_hg + 313);

    auto tr_zzzzz_zzzz = pbuffer.data(idx_op_hg + 314);

    // Set up 0-15 components of targeted buffer : GG

    auto tr_0_0_x_xxxx_xxxx = pbuffer.data(idx_op_geom_010_gg);

    auto tr_0_0_x_xxxx_xxxy = pbuffer.data(idx_op_geom_010_gg + 1);

    auto tr_0_0_x_xxxx_xxxz = pbuffer.data(idx_op_geom_010_gg + 2);

    auto tr_0_0_x_xxxx_xxyy = pbuffer.data(idx_op_geom_010_gg + 3);

    auto tr_0_0_x_xxxx_xxyz = pbuffer.data(idx_op_geom_010_gg + 4);

    auto tr_0_0_x_xxxx_xxzz = pbuffer.data(idx_op_geom_010_gg + 5);

    auto tr_0_0_x_xxxx_xyyy = pbuffer.data(idx_op_geom_010_gg + 6);

    auto tr_0_0_x_xxxx_xyyz = pbuffer.data(idx_op_geom_010_gg + 7);

    auto tr_0_0_x_xxxx_xyzz = pbuffer.data(idx_op_geom_010_gg + 8);

    auto tr_0_0_x_xxxx_xzzz = pbuffer.data(idx_op_geom_010_gg + 9);

    auto tr_0_0_x_xxxx_yyyy = pbuffer.data(idx_op_geom_010_gg + 10);

    auto tr_0_0_x_xxxx_yyyz = pbuffer.data(idx_op_geom_010_gg + 11);

    auto tr_0_0_x_xxxx_yyzz = pbuffer.data(idx_op_geom_010_gg + 12);

    auto tr_0_0_x_xxxx_yzzz = pbuffer.data(idx_op_geom_010_gg + 13);

    auto tr_0_0_x_xxxx_zzzz = pbuffer.data(idx_op_geom_010_gg + 14);

    #pragma omp simd aligned(tr_0_0_x_xxxx_xxxx, tr_0_0_x_xxxx_xxxy, tr_0_0_x_xxxx_xxxz, tr_0_0_x_xxxx_xxyy, tr_0_0_x_xxxx_xxyz, tr_0_0_x_xxxx_xxzz, tr_0_0_x_xxxx_xyyy, tr_0_0_x_xxxx_xyyz, tr_0_0_x_xxxx_xyzz, tr_0_0_x_xxxx_xzzz, tr_0_0_x_xxxx_yyyy, tr_0_0_x_xxxx_yyyz, tr_0_0_x_xxxx_yyzz, tr_0_0_x_xxxx_yzzz, tr_0_0_x_xxxx_zzzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxx_xxx, tr_xxxx_xxxxx, tr_xxxx_xxxxy, tr_xxxx_xxxxz, tr_xxxx_xxxyy, tr_xxxx_xxxyz, tr_xxxx_xxxzz, tr_xxxx_xxy, tr_xxxx_xxyyy, tr_xxxx_xxyyz, tr_xxxx_xxyzz, tr_xxxx_xxz, tr_xxxx_xxzzz, tr_xxxx_xyy, tr_xxxx_xyyyy, tr_xxxx_xyyyz, tr_xxxx_xyyzz, tr_xxxx_xyz, tr_xxxx_xyzzz, tr_xxxx_xzz, tr_xxxx_xzzzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, tr_xxxxx_xxxx, tr_xxxxx_xxxy, tr_xxxxx_xxxz, tr_xxxxx_xxyy, tr_xxxxx_xxyz, tr_xxxxx_xxzz, tr_xxxxx_xyyy, tr_xxxxx_xyyz, tr_xxxxx_xyzz, tr_xxxxx_xzzz, tr_xxxxx_yyyy, tr_xxxxx_yyyz, tr_xxxxx_yyzz, tr_xxxxx_yzzz, tr_xxxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxx_xxxx[i] = 2.0 * tr_xxxxx_xxxx[i] * tbe_0 + 2.0 * tr_xxxx_xxxxx[i] * tke_0 - 4.0 * tr_xxx_xxxx[i] - 4.0 * tr_xxxx_xxx[i];

        tr_0_0_x_xxxx_xxxy[i] = 2.0 * tr_xxxxx_xxxy[i] * tbe_0 + 2.0 * tr_xxxx_xxxxy[i] * tke_0 - 4.0 * tr_xxx_xxxy[i] - 3.0 * tr_xxxx_xxy[i];

        tr_0_0_x_xxxx_xxxz[i] = 2.0 * tr_xxxxx_xxxz[i] * tbe_0 + 2.0 * tr_xxxx_xxxxz[i] * tke_0 - 4.0 * tr_xxx_xxxz[i] - 3.0 * tr_xxxx_xxz[i];

        tr_0_0_x_xxxx_xxyy[i] = 2.0 * tr_xxxxx_xxyy[i] * tbe_0 + 2.0 * tr_xxxx_xxxyy[i] * tke_0 - 4.0 * tr_xxx_xxyy[i] - 2.0 * tr_xxxx_xyy[i];

        tr_0_0_x_xxxx_xxyz[i] = 2.0 * tr_xxxxx_xxyz[i] * tbe_0 + 2.0 * tr_xxxx_xxxyz[i] * tke_0 - 4.0 * tr_xxx_xxyz[i] - 2.0 * tr_xxxx_xyz[i];

        tr_0_0_x_xxxx_xxzz[i] = 2.0 * tr_xxxxx_xxzz[i] * tbe_0 + 2.0 * tr_xxxx_xxxzz[i] * tke_0 - 4.0 * tr_xxx_xxzz[i] - 2.0 * tr_xxxx_xzz[i];

        tr_0_0_x_xxxx_xyyy[i] = 2.0 * tr_xxxxx_xyyy[i] * tbe_0 + 2.0 * tr_xxxx_xxyyy[i] * tke_0 - 4.0 * tr_xxx_xyyy[i] - tr_xxxx_yyy[i];

        tr_0_0_x_xxxx_xyyz[i] = 2.0 * tr_xxxxx_xyyz[i] * tbe_0 + 2.0 * tr_xxxx_xxyyz[i] * tke_0 - 4.0 * tr_xxx_xyyz[i] - tr_xxxx_yyz[i];

        tr_0_0_x_xxxx_xyzz[i] = 2.0 * tr_xxxxx_xyzz[i] * tbe_0 + 2.0 * tr_xxxx_xxyzz[i] * tke_0 - 4.0 * tr_xxx_xyzz[i] - tr_xxxx_yzz[i];

        tr_0_0_x_xxxx_xzzz[i] = 2.0 * tr_xxxxx_xzzz[i] * tbe_0 + 2.0 * tr_xxxx_xxzzz[i] * tke_0 - 4.0 * tr_xxx_xzzz[i] - tr_xxxx_zzz[i];

        tr_0_0_x_xxxx_yyyy[i] = 2.0 * tr_xxxxx_yyyy[i] * tbe_0 + 2.0 * tr_xxxx_xyyyy[i] * tke_0 - 4.0 * tr_xxx_yyyy[i];

        tr_0_0_x_xxxx_yyyz[i] = 2.0 * tr_xxxxx_yyyz[i] * tbe_0 + 2.0 * tr_xxxx_xyyyz[i] * tke_0 - 4.0 * tr_xxx_yyyz[i];

        tr_0_0_x_xxxx_yyzz[i] = 2.0 * tr_xxxxx_yyzz[i] * tbe_0 + 2.0 * tr_xxxx_xyyzz[i] * tke_0 - 4.0 * tr_xxx_yyzz[i];

        tr_0_0_x_xxxx_yzzz[i] = 2.0 * tr_xxxxx_yzzz[i] * tbe_0 + 2.0 * tr_xxxx_xyzzz[i] * tke_0 - 4.0 * tr_xxx_yzzz[i];

        tr_0_0_x_xxxx_zzzz[i] = 2.0 * tr_xxxxx_zzzz[i] * tbe_0 + 2.0 * tr_xxxx_xzzzz[i] * tke_0 - 4.0 * tr_xxx_zzzz[i];
    }

    // Set up 15-30 components of targeted buffer : GG

    auto tr_0_0_x_xxxy_xxxx = pbuffer.data(idx_op_geom_010_gg + 15);

    auto tr_0_0_x_xxxy_xxxy = pbuffer.data(idx_op_geom_010_gg + 16);

    auto tr_0_0_x_xxxy_xxxz = pbuffer.data(idx_op_geom_010_gg + 17);

    auto tr_0_0_x_xxxy_xxyy = pbuffer.data(idx_op_geom_010_gg + 18);

    auto tr_0_0_x_xxxy_xxyz = pbuffer.data(idx_op_geom_010_gg + 19);

    auto tr_0_0_x_xxxy_xxzz = pbuffer.data(idx_op_geom_010_gg + 20);

    auto tr_0_0_x_xxxy_xyyy = pbuffer.data(idx_op_geom_010_gg + 21);

    auto tr_0_0_x_xxxy_xyyz = pbuffer.data(idx_op_geom_010_gg + 22);

    auto tr_0_0_x_xxxy_xyzz = pbuffer.data(idx_op_geom_010_gg + 23);

    auto tr_0_0_x_xxxy_xzzz = pbuffer.data(idx_op_geom_010_gg + 24);

    auto tr_0_0_x_xxxy_yyyy = pbuffer.data(idx_op_geom_010_gg + 25);

    auto tr_0_0_x_xxxy_yyyz = pbuffer.data(idx_op_geom_010_gg + 26);

    auto tr_0_0_x_xxxy_yyzz = pbuffer.data(idx_op_geom_010_gg + 27);

    auto tr_0_0_x_xxxy_yzzz = pbuffer.data(idx_op_geom_010_gg + 28);

    auto tr_0_0_x_xxxy_zzzz = pbuffer.data(idx_op_geom_010_gg + 29);

    #pragma omp simd aligned(tr_0_0_x_xxxy_xxxx, tr_0_0_x_xxxy_xxxy, tr_0_0_x_xxxy_xxxz, tr_0_0_x_xxxy_xxyy, tr_0_0_x_xxxy_xxyz, tr_0_0_x_xxxy_xxzz, tr_0_0_x_xxxy_xyyy, tr_0_0_x_xxxy_xyyz, tr_0_0_x_xxxy_xyzz, tr_0_0_x_xxxy_xzzz, tr_0_0_x_xxxy_yyyy, tr_0_0_x_xxxy_yyyz, tr_0_0_x_xxxy_yyzz, tr_0_0_x_xxxy_yzzz, tr_0_0_x_xxxy_zzzz, tr_xxxxy_xxxx, tr_xxxxy_xxxy, tr_xxxxy_xxxz, tr_xxxxy_xxyy, tr_xxxxy_xxyz, tr_xxxxy_xxzz, tr_xxxxy_xyyy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xzzz, tr_xxxxy_yyyy, tr_xxxxy_yyyz, tr_xxxxy_yyzz, tr_xxxxy_yzzz, tr_xxxxy_zzzz, tr_xxxy_xxx, tr_xxxy_xxxxx, tr_xxxy_xxxxy, tr_xxxy_xxxxz, tr_xxxy_xxxyy, tr_xxxy_xxxyz, tr_xxxy_xxxzz, tr_xxxy_xxy, tr_xxxy_xxyyy, tr_xxxy_xxyyz, tr_xxxy_xxyzz, tr_xxxy_xxz, tr_xxxy_xxzzz, tr_xxxy_xyy, tr_xxxy_xyyyy, tr_xxxy_xyyyz, tr_xxxy_xyyzz, tr_xxxy_xyz, tr_xxxy_xyzzz, tr_xxxy_xzz, tr_xxxy_xzzzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxy_xxxx[i] = 2.0 * tr_xxxxy_xxxx[i] * tbe_0 + 2.0 * tr_xxxy_xxxxx[i] * tke_0 - 3.0 * tr_xxy_xxxx[i] - 4.0 * tr_xxxy_xxx[i];

        tr_0_0_x_xxxy_xxxy[i] = 2.0 * tr_xxxxy_xxxy[i] * tbe_0 + 2.0 * tr_xxxy_xxxxy[i] * tke_0 - 3.0 * tr_xxy_xxxy[i] - 3.0 * tr_xxxy_xxy[i];

        tr_0_0_x_xxxy_xxxz[i] = 2.0 * tr_xxxxy_xxxz[i] * tbe_0 + 2.0 * tr_xxxy_xxxxz[i] * tke_0 - 3.0 * tr_xxy_xxxz[i] - 3.0 * tr_xxxy_xxz[i];

        tr_0_0_x_xxxy_xxyy[i] = 2.0 * tr_xxxxy_xxyy[i] * tbe_0 + 2.0 * tr_xxxy_xxxyy[i] * tke_0 - 3.0 * tr_xxy_xxyy[i] - 2.0 * tr_xxxy_xyy[i];

        tr_0_0_x_xxxy_xxyz[i] = 2.0 * tr_xxxxy_xxyz[i] * tbe_0 + 2.0 * tr_xxxy_xxxyz[i] * tke_0 - 3.0 * tr_xxy_xxyz[i] - 2.0 * tr_xxxy_xyz[i];

        tr_0_0_x_xxxy_xxzz[i] = 2.0 * tr_xxxxy_xxzz[i] * tbe_0 + 2.0 * tr_xxxy_xxxzz[i] * tke_0 - 3.0 * tr_xxy_xxzz[i] - 2.0 * tr_xxxy_xzz[i];

        tr_0_0_x_xxxy_xyyy[i] = 2.0 * tr_xxxxy_xyyy[i] * tbe_0 + 2.0 * tr_xxxy_xxyyy[i] * tke_0 - 3.0 * tr_xxy_xyyy[i] - tr_xxxy_yyy[i];

        tr_0_0_x_xxxy_xyyz[i] = 2.0 * tr_xxxxy_xyyz[i] * tbe_0 + 2.0 * tr_xxxy_xxyyz[i] * tke_0 - 3.0 * tr_xxy_xyyz[i] - tr_xxxy_yyz[i];

        tr_0_0_x_xxxy_xyzz[i] = 2.0 * tr_xxxxy_xyzz[i] * tbe_0 + 2.0 * tr_xxxy_xxyzz[i] * tke_0 - 3.0 * tr_xxy_xyzz[i] - tr_xxxy_yzz[i];

        tr_0_0_x_xxxy_xzzz[i] = 2.0 * tr_xxxxy_xzzz[i] * tbe_0 + 2.0 * tr_xxxy_xxzzz[i] * tke_0 - 3.0 * tr_xxy_xzzz[i] - tr_xxxy_zzz[i];

        tr_0_0_x_xxxy_yyyy[i] = 2.0 * tr_xxxxy_yyyy[i] * tbe_0 + 2.0 * tr_xxxy_xyyyy[i] * tke_0 - 3.0 * tr_xxy_yyyy[i];

        tr_0_0_x_xxxy_yyyz[i] = 2.0 * tr_xxxxy_yyyz[i] * tbe_0 + 2.0 * tr_xxxy_xyyyz[i] * tke_0 - 3.0 * tr_xxy_yyyz[i];

        tr_0_0_x_xxxy_yyzz[i] = 2.0 * tr_xxxxy_yyzz[i] * tbe_0 + 2.0 * tr_xxxy_xyyzz[i] * tke_0 - 3.0 * tr_xxy_yyzz[i];

        tr_0_0_x_xxxy_yzzz[i] = 2.0 * tr_xxxxy_yzzz[i] * tbe_0 + 2.0 * tr_xxxy_xyzzz[i] * tke_0 - 3.0 * tr_xxy_yzzz[i];

        tr_0_0_x_xxxy_zzzz[i] = 2.0 * tr_xxxxy_zzzz[i] * tbe_0 + 2.0 * tr_xxxy_xzzzz[i] * tke_0 - 3.0 * tr_xxy_zzzz[i];
    }

    // Set up 30-45 components of targeted buffer : GG

    auto tr_0_0_x_xxxz_xxxx = pbuffer.data(idx_op_geom_010_gg + 30);

    auto tr_0_0_x_xxxz_xxxy = pbuffer.data(idx_op_geom_010_gg + 31);

    auto tr_0_0_x_xxxz_xxxz = pbuffer.data(idx_op_geom_010_gg + 32);

    auto tr_0_0_x_xxxz_xxyy = pbuffer.data(idx_op_geom_010_gg + 33);

    auto tr_0_0_x_xxxz_xxyz = pbuffer.data(idx_op_geom_010_gg + 34);

    auto tr_0_0_x_xxxz_xxzz = pbuffer.data(idx_op_geom_010_gg + 35);

    auto tr_0_0_x_xxxz_xyyy = pbuffer.data(idx_op_geom_010_gg + 36);

    auto tr_0_0_x_xxxz_xyyz = pbuffer.data(idx_op_geom_010_gg + 37);

    auto tr_0_0_x_xxxz_xyzz = pbuffer.data(idx_op_geom_010_gg + 38);

    auto tr_0_0_x_xxxz_xzzz = pbuffer.data(idx_op_geom_010_gg + 39);

    auto tr_0_0_x_xxxz_yyyy = pbuffer.data(idx_op_geom_010_gg + 40);

    auto tr_0_0_x_xxxz_yyyz = pbuffer.data(idx_op_geom_010_gg + 41);

    auto tr_0_0_x_xxxz_yyzz = pbuffer.data(idx_op_geom_010_gg + 42);

    auto tr_0_0_x_xxxz_yzzz = pbuffer.data(idx_op_geom_010_gg + 43);

    auto tr_0_0_x_xxxz_zzzz = pbuffer.data(idx_op_geom_010_gg + 44);

    #pragma omp simd aligned(tr_0_0_x_xxxz_xxxx, tr_0_0_x_xxxz_xxxy, tr_0_0_x_xxxz_xxxz, tr_0_0_x_xxxz_xxyy, tr_0_0_x_xxxz_xxyz, tr_0_0_x_xxxz_xxzz, tr_0_0_x_xxxz_xyyy, tr_0_0_x_xxxz_xyyz, tr_0_0_x_xxxz_xyzz, tr_0_0_x_xxxz_xzzz, tr_0_0_x_xxxz_yyyy, tr_0_0_x_xxxz_yyyz, tr_0_0_x_xxxz_yyzz, tr_0_0_x_xxxz_yzzz, tr_0_0_x_xxxz_zzzz, tr_xxxxz_xxxx, tr_xxxxz_xxxy, tr_xxxxz_xxxz, tr_xxxxz_xxyy, tr_xxxxz_xxyz, tr_xxxxz_xxzz, tr_xxxxz_xyyy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xzzz, tr_xxxxz_yyyy, tr_xxxxz_yyyz, tr_xxxxz_yyzz, tr_xxxxz_yzzz, tr_xxxxz_zzzz, tr_xxxz_xxx, tr_xxxz_xxxxx, tr_xxxz_xxxxy, tr_xxxz_xxxxz, tr_xxxz_xxxyy, tr_xxxz_xxxyz, tr_xxxz_xxxzz, tr_xxxz_xxy, tr_xxxz_xxyyy, tr_xxxz_xxyyz, tr_xxxz_xxyzz, tr_xxxz_xxz, tr_xxxz_xxzzz, tr_xxxz_xyy, tr_xxxz_xyyyy, tr_xxxz_xyyyz, tr_xxxz_xyyzz, tr_xxxz_xyz, tr_xxxz_xyzzz, tr_xxxz_xzz, tr_xxxz_xzzzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxz_xxxx[i] = 2.0 * tr_xxxxz_xxxx[i] * tbe_0 + 2.0 * tr_xxxz_xxxxx[i] * tke_0 - 3.0 * tr_xxz_xxxx[i] - 4.0 * tr_xxxz_xxx[i];

        tr_0_0_x_xxxz_xxxy[i] = 2.0 * tr_xxxxz_xxxy[i] * tbe_0 + 2.0 * tr_xxxz_xxxxy[i] * tke_0 - 3.0 * tr_xxz_xxxy[i] - 3.0 * tr_xxxz_xxy[i];

        tr_0_0_x_xxxz_xxxz[i] = 2.0 * tr_xxxxz_xxxz[i] * tbe_0 + 2.0 * tr_xxxz_xxxxz[i] * tke_0 - 3.0 * tr_xxz_xxxz[i] - 3.0 * tr_xxxz_xxz[i];

        tr_0_0_x_xxxz_xxyy[i] = 2.0 * tr_xxxxz_xxyy[i] * tbe_0 + 2.0 * tr_xxxz_xxxyy[i] * tke_0 - 3.0 * tr_xxz_xxyy[i] - 2.0 * tr_xxxz_xyy[i];

        tr_0_0_x_xxxz_xxyz[i] = 2.0 * tr_xxxxz_xxyz[i] * tbe_0 + 2.0 * tr_xxxz_xxxyz[i] * tke_0 - 3.0 * tr_xxz_xxyz[i] - 2.0 * tr_xxxz_xyz[i];

        tr_0_0_x_xxxz_xxzz[i] = 2.0 * tr_xxxxz_xxzz[i] * tbe_0 + 2.0 * tr_xxxz_xxxzz[i] * tke_0 - 3.0 * tr_xxz_xxzz[i] - 2.0 * tr_xxxz_xzz[i];

        tr_0_0_x_xxxz_xyyy[i] = 2.0 * tr_xxxxz_xyyy[i] * tbe_0 + 2.0 * tr_xxxz_xxyyy[i] * tke_0 - 3.0 * tr_xxz_xyyy[i] - tr_xxxz_yyy[i];

        tr_0_0_x_xxxz_xyyz[i] = 2.0 * tr_xxxxz_xyyz[i] * tbe_0 + 2.0 * tr_xxxz_xxyyz[i] * tke_0 - 3.0 * tr_xxz_xyyz[i] - tr_xxxz_yyz[i];

        tr_0_0_x_xxxz_xyzz[i] = 2.0 * tr_xxxxz_xyzz[i] * tbe_0 + 2.0 * tr_xxxz_xxyzz[i] * tke_0 - 3.0 * tr_xxz_xyzz[i] - tr_xxxz_yzz[i];

        tr_0_0_x_xxxz_xzzz[i] = 2.0 * tr_xxxxz_xzzz[i] * tbe_0 + 2.0 * tr_xxxz_xxzzz[i] * tke_0 - 3.0 * tr_xxz_xzzz[i] - tr_xxxz_zzz[i];

        tr_0_0_x_xxxz_yyyy[i] = 2.0 * tr_xxxxz_yyyy[i] * tbe_0 + 2.0 * tr_xxxz_xyyyy[i] * tke_0 - 3.0 * tr_xxz_yyyy[i];

        tr_0_0_x_xxxz_yyyz[i] = 2.0 * tr_xxxxz_yyyz[i] * tbe_0 + 2.0 * tr_xxxz_xyyyz[i] * tke_0 - 3.0 * tr_xxz_yyyz[i];

        tr_0_0_x_xxxz_yyzz[i] = 2.0 * tr_xxxxz_yyzz[i] * tbe_0 + 2.0 * tr_xxxz_xyyzz[i] * tke_0 - 3.0 * tr_xxz_yyzz[i];

        tr_0_0_x_xxxz_yzzz[i] = 2.0 * tr_xxxxz_yzzz[i] * tbe_0 + 2.0 * tr_xxxz_xyzzz[i] * tke_0 - 3.0 * tr_xxz_yzzz[i];

        tr_0_0_x_xxxz_zzzz[i] = 2.0 * tr_xxxxz_zzzz[i] * tbe_0 + 2.0 * tr_xxxz_xzzzz[i] * tke_0 - 3.0 * tr_xxz_zzzz[i];
    }

    // Set up 45-60 components of targeted buffer : GG

    auto tr_0_0_x_xxyy_xxxx = pbuffer.data(idx_op_geom_010_gg + 45);

    auto tr_0_0_x_xxyy_xxxy = pbuffer.data(idx_op_geom_010_gg + 46);

    auto tr_0_0_x_xxyy_xxxz = pbuffer.data(idx_op_geom_010_gg + 47);

    auto tr_0_0_x_xxyy_xxyy = pbuffer.data(idx_op_geom_010_gg + 48);

    auto tr_0_0_x_xxyy_xxyz = pbuffer.data(idx_op_geom_010_gg + 49);

    auto tr_0_0_x_xxyy_xxzz = pbuffer.data(idx_op_geom_010_gg + 50);

    auto tr_0_0_x_xxyy_xyyy = pbuffer.data(idx_op_geom_010_gg + 51);

    auto tr_0_0_x_xxyy_xyyz = pbuffer.data(idx_op_geom_010_gg + 52);

    auto tr_0_0_x_xxyy_xyzz = pbuffer.data(idx_op_geom_010_gg + 53);

    auto tr_0_0_x_xxyy_xzzz = pbuffer.data(idx_op_geom_010_gg + 54);

    auto tr_0_0_x_xxyy_yyyy = pbuffer.data(idx_op_geom_010_gg + 55);

    auto tr_0_0_x_xxyy_yyyz = pbuffer.data(idx_op_geom_010_gg + 56);

    auto tr_0_0_x_xxyy_yyzz = pbuffer.data(idx_op_geom_010_gg + 57);

    auto tr_0_0_x_xxyy_yzzz = pbuffer.data(idx_op_geom_010_gg + 58);

    auto tr_0_0_x_xxyy_zzzz = pbuffer.data(idx_op_geom_010_gg + 59);

    #pragma omp simd aligned(tr_0_0_x_xxyy_xxxx, tr_0_0_x_xxyy_xxxy, tr_0_0_x_xxyy_xxxz, tr_0_0_x_xxyy_xxyy, tr_0_0_x_xxyy_xxyz, tr_0_0_x_xxyy_xxzz, tr_0_0_x_xxyy_xyyy, tr_0_0_x_xxyy_xyyz, tr_0_0_x_xxyy_xyzz, tr_0_0_x_xxyy_xzzz, tr_0_0_x_xxyy_yyyy, tr_0_0_x_xxyy_yyyz, tr_0_0_x_xxyy_yyzz, tr_0_0_x_xxyy_yzzz, tr_0_0_x_xxyy_zzzz, tr_xxxyy_xxxx, tr_xxxyy_xxxy, tr_xxxyy_xxxz, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xzzz, tr_xxxyy_yyyy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yzzz, tr_xxxyy_zzzz, tr_xxyy_xxx, tr_xxyy_xxxxx, tr_xxyy_xxxxy, tr_xxyy_xxxxz, tr_xxyy_xxxyy, tr_xxyy_xxxyz, tr_xxyy_xxxzz, tr_xxyy_xxy, tr_xxyy_xxyyy, tr_xxyy_xxyyz, tr_xxyy_xxyzz, tr_xxyy_xxz, tr_xxyy_xxzzz, tr_xxyy_xyy, tr_xxyy_xyyyy, tr_xxyy_xyyyz, tr_xxyy_xyyzz, tr_xxyy_xyz, tr_xxyy_xyzzz, tr_xxyy_xzz, tr_xxyy_xzzzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxyy_xxxx[i] = 2.0 * tr_xxxyy_xxxx[i] * tbe_0 + 2.0 * tr_xxyy_xxxxx[i] * tke_0 - 2.0 * tr_xyy_xxxx[i] - 4.0 * tr_xxyy_xxx[i];

        tr_0_0_x_xxyy_xxxy[i] = 2.0 * tr_xxxyy_xxxy[i] * tbe_0 + 2.0 * tr_xxyy_xxxxy[i] * tke_0 - 2.0 * tr_xyy_xxxy[i] - 3.0 * tr_xxyy_xxy[i];

        tr_0_0_x_xxyy_xxxz[i] = 2.0 * tr_xxxyy_xxxz[i] * tbe_0 + 2.0 * tr_xxyy_xxxxz[i] * tke_0 - 2.0 * tr_xyy_xxxz[i] - 3.0 * tr_xxyy_xxz[i];

        tr_0_0_x_xxyy_xxyy[i] = 2.0 * tr_xxxyy_xxyy[i] * tbe_0 + 2.0 * tr_xxyy_xxxyy[i] * tke_0 - 2.0 * tr_xyy_xxyy[i] - 2.0 * tr_xxyy_xyy[i];

        tr_0_0_x_xxyy_xxyz[i] = 2.0 * tr_xxxyy_xxyz[i] * tbe_0 + 2.0 * tr_xxyy_xxxyz[i] * tke_0 - 2.0 * tr_xyy_xxyz[i] - 2.0 * tr_xxyy_xyz[i];

        tr_0_0_x_xxyy_xxzz[i] = 2.0 * tr_xxxyy_xxzz[i] * tbe_0 + 2.0 * tr_xxyy_xxxzz[i] * tke_0 - 2.0 * tr_xyy_xxzz[i] - 2.0 * tr_xxyy_xzz[i];

        tr_0_0_x_xxyy_xyyy[i] = 2.0 * tr_xxxyy_xyyy[i] * tbe_0 + 2.0 * tr_xxyy_xxyyy[i] * tke_0 - 2.0 * tr_xyy_xyyy[i] - tr_xxyy_yyy[i];

        tr_0_0_x_xxyy_xyyz[i] = 2.0 * tr_xxxyy_xyyz[i] * tbe_0 + 2.0 * tr_xxyy_xxyyz[i] * tke_0 - 2.0 * tr_xyy_xyyz[i] - tr_xxyy_yyz[i];

        tr_0_0_x_xxyy_xyzz[i] = 2.0 * tr_xxxyy_xyzz[i] * tbe_0 + 2.0 * tr_xxyy_xxyzz[i] * tke_0 - 2.0 * tr_xyy_xyzz[i] - tr_xxyy_yzz[i];

        tr_0_0_x_xxyy_xzzz[i] = 2.0 * tr_xxxyy_xzzz[i] * tbe_0 + 2.0 * tr_xxyy_xxzzz[i] * tke_0 - 2.0 * tr_xyy_xzzz[i] - tr_xxyy_zzz[i];

        tr_0_0_x_xxyy_yyyy[i] = 2.0 * tr_xxxyy_yyyy[i] * tbe_0 + 2.0 * tr_xxyy_xyyyy[i] * tke_0 - 2.0 * tr_xyy_yyyy[i];

        tr_0_0_x_xxyy_yyyz[i] = 2.0 * tr_xxxyy_yyyz[i] * tbe_0 + 2.0 * tr_xxyy_xyyyz[i] * tke_0 - 2.0 * tr_xyy_yyyz[i];

        tr_0_0_x_xxyy_yyzz[i] = 2.0 * tr_xxxyy_yyzz[i] * tbe_0 + 2.0 * tr_xxyy_xyyzz[i] * tke_0 - 2.0 * tr_xyy_yyzz[i];

        tr_0_0_x_xxyy_yzzz[i] = 2.0 * tr_xxxyy_yzzz[i] * tbe_0 + 2.0 * tr_xxyy_xyzzz[i] * tke_0 - 2.0 * tr_xyy_yzzz[i];

        tr_0_0_x_xxyy_zzzz[i] = 2.0 * tr_xxxyy_zzzz[i] * tbe_0 + 2.0 * tr_xxyy_xzzzz[i] * tke_0 - 2.0 * tr_xyy_zzzz[i];
    }

    // Set up 60-75 components of targeted buffer : GG

    auto tr_0_0_x_xxyz_xxxx = pbuffer.data(idx_op_geom_010_gg + 60);

    auto tr_0_0_x_xxyz_xxxy = pbuffer.data(idx_op_geom_010_gg + 61);

    auto tr_0_0_x_xxyz_xxxz = pbuffer.data(idx_op_geom_010_gg + 62);

    auto tr_0_0_x_xxyz_xxyy = pbuffer.data(idx_op_geom_010_gg + 63);

    auto tr_0_0_x_xxyz_xxyz = pbuffer.data(idx_op_geom_010_gg + 64);

    auto tr_0_0_x_xxyz_xxzz = pbuffer.data(idx_op_geom_010_gg + 65);

    auto tr_0_0_x_xxyz_xyyy = pbuffer.data(idx_op_geom_010_gg + 66);

    auto tr_0_0_x_xxyz_xyyz = pbuffer.data(idx_op_geom_010_gg + 67);

    auto tr_0_0_x_xxyz_xyzz = pbuffer.data(idx_op_geom_010_gg + 68);

    auto tr_0_0_x_xxyz_xzzz = pbuffer.data(idx_op_geom_010_gg + 69);

    auto tr_0_0_x_xxyz_yyyy = pbuffer.data(idx_op_geom_010_gg + 70);

    auto tr_0_0_x_xxyz_yyyz = pbuffer.data(idx_op_geom_010_gg + 71);

    auto tr_0_0_x_xxyz_yyzz = pbuffer.data(idx_op_geom_010_gg + 72);

    auto tr_0_0_x_xxyz_yzzz = pbuffer.data(idx_op_geom_010_gg + 73);

    auto tr_0_0_x_xxyz_zzzz = pbuffer.data(idx_op_geom_010_gg + 74);

    #pragma omp simd aligned(tr_0_0_x_xxyz_xxxx, tr_0_0_x_xxyz_xxxy, tr_0_0_x_xxyz_xxxz, tr_0_0_x_xxyz_xxyy, tr_0_0_x_xxyz_xxyz, tr_0_0_x_xxyz_xxzz, tr_0_0_x_xxyz_xyyy, tr_0_0_x_xxyz_xyyz, tr_0_0_x_xxyz_xyzz, tr_0_0_x_xxyz_xzzz, tr_0_0_x_xxyz_yyyy, tr_0_0_x_xxyz_yyyz, tr_0_0_x_xxyz_yyzz, tr_0_0_x_xxyz_yzzz, tr_0_0_x_xxyz_zzzz, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxx, tr_xxyz_xxxxy, tr_xxyz_xxxxz, tr_xxyz_xxxyy, tr_xxyz_xxxyz, tr_xxyz_xxxzz, tr_xxyz_xxy, tr_xxyz_xxyyy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xxzzz, tr_xxyz_xyy, tr_xxyz_xyyyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_xzzzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxyz_xxxx[i] = 2.0 * tr_xxxyz_xxxx[i] * tbe_0 + 2.0 * tr_xxyz_xxxxx[i] * tke_0 - 2.0 * tr_xyz_xxxx[i] - 4.0 * tr_xxyz_xxx[i];

        tr_0_0_x_xxyz_xxxy[i] = 2.0 * tr_xxxyz_xxxy[i] * tbe_0 + 2.0 * tr_xxyz_xxxxy[i] * tke_0 - 2.0 * tr_xyz_xxxy[i] - 3.0 * tr_xxyz_xxy[i];

        tr_0_0_x_xxyz_xxxz[i] = 2.0 * tr_xxxyz_xxxz[i] * tbe_0 + 2.0 * tr_xxyz_xxxxz[i] * tke_0 - 2.0 * tr_xyz_xxxz[i] - 3.0 * tr_xxyz_xxz[i];

        tr_0_0_x_xxyz_xxyy[i] = 2.0 * tr_xxxyz_xxyy[i] * tbe_0 + 2.0 * tr_xxyz_xxxyy[i] * tke_0 - 2.0 * tr_xyz_xxyy[i] - 2.0 * tr_xxyz_xyy[i];

        tr_0_0_x_xxyz_xxyz[i] = 2.0 * tr_xxxyz_xxyz[i] * tbe_0 + 2.0 * tr_xxyz_xxxyz[i] * tke_0 - 2.0 * tr_xyz_xxyz[i] - 2.0 * tr_xxyz_xyz[i];

        tr_0_0_x_xxyz_xxzz[i] = 2.0 * tr_xxxyz_xxzz[i] * tbe_0 + 2.0 * tr_xxyz_xxxzz[i] * tke_0 - 2.0 * tr_xyz_xxzz[i] - 2.0 * tr_xxyz_xzz[i];

        tr_0_0_x_xxyz_xyyy[i] = 2.0 * tr_xxxyz_xyyy[i] * tbe_0 + 2.0 * tr_xxyz_xxyyy[i] * tke_0 - 2.0 * tr_xyz_xyyy[i] - tr_xxyz_yyy[i];

        tr_0_0_x_xxyz_xyyz[i] = 2.0 * tr_xxxyz_xyyz[i] * tbe_0 + 2.0 * tr_xxyz_xxyyz[i] * tke_0 - 2.0 * tr_xyz_xyyz[i] - tr_xxyz_yyz[i];

        tr_0_0_x_xxyz_xyzz[i] = 2.0 * tr_xxxyz_xyzz[i] * tbe_0 + 2.0 * tr_xxyz_xxyzz[i] * tke_0 - 2.0 * tr_xyz_xyzz[i] - tr_xxyz_yzz[i];

        tr_0_0_x_xxyz_xzzz[i] = 2.0 * tr_xxxyz_xzzz[i] * tbe_0 + 2.0 * tr_xxyz_xxzzz[i] * tke_0 - 2.0 * tr_xyz_xzzz[i] - tr_xxyz_zzz[i];

        tr_0_0_x_xxyz_yyyy[i] = 2.0 * tr_xxxyz_yyyy[i] * tbe_0 + 2.0 * tr_xxyz_xyyyy[i] * tke_0 - 2.0 * tr_xyz_yyyy[i];

        tr_0_0_x_xxyz_yyyz[i] = 2.0 * tr_xxxyz_yyyz[i] * tbe_0 + 2.0 * tr_xxyz_xyyyz[i] * tke_0 - 2.0 * tr_xyz_yyyz[i];

        tr_0_0_x_xxyz_yyzz[i] = 2.0 * tr_xxxyz_yyzz[i] * tbe_0 + 2.0 * tr_xxyz_xyyzz[i] * tke_0 - 2.0 * tr_xyz_yyzz[i];

        tr_0_0_x_xxyz_yzzz[i] = 2.0 * tr_xxxyz_yzzz[i] * tbe_0 + 2.0 * tr_xxyz_xyzzz[i] * tke_0 - 2.0 * tr_xyz_yzzz[i];

        tr_0_0_x_xxyz_zzzz[i] = 2.0 * tr_xxxyz_zzzz[i] * tbe_0 + 2.0 * tr_xxyz_xzzzz[i] * tke_0 - 2.0 * tr_xyz_zzzz[i];
    }

    // Set up 75-90 components of targeted buffer : GG

    auto tr_0_0_x_xxzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 75);

    auto tr_0_0_x_xxzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 76);

    auto tr_0_0_x_xxzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 77);

    auto tr_0_0_x_xxzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 78);

    auto tr_0_0_x_xxzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 79);

    auto tr_0_0_x_xxzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 80);

    auto tr_0_0_x_xxzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 81);

    auto tr_0_0_x_xxzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 82);

    auto tr_0_0_x_xxzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 83);

    auto tr_0_0_x_xxzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 84);

    auto tr_0_0_x_xxzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 85);

    auto tr_0_0_x_xxzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 86);

    auto tr_0_0_x_xxzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 87);

    auto tr_0_0_x_xxzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 88);

    auto tr_0_0_x_xxzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 89);

    #pragma omp simd aligned(tr_0_0_x_xxzz_xxxx, tr_0_0_x_xxzz_xxxy, tr_0_0_x_xxzz_xxxz, tr_0_0_x_xxzz_xxyy, tr_0_0_x_xxzz_xxyz, tr_0_0_x_xxzz_xxzz, tr_0_0_x_xxzz_xyyy, tr_0_0_x_xxzz_xyyz, tr_0_0_x_xxzz_xyzz, tr_0_0_x_xxzz_xzzz, tr_0_0_x_xxzz_yyyy, tr_0_0_x_xxzz_yyyz, tr_0_0_x_xxzz_yyzz, tr_0_0_x_xxzz_yzzz, tr_0_0_x_xxzz_zzzz, tr_xxxzz_xxxx, tr_xxxzz_xxxy, tr_xxxzz_xxxz, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xzzz, tr_xxxzz_yyyy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yzzz, tr_xxxzz_zzzz, tr_xxzz_xxx, tr_xxzz_xxxxx, tr_xxzz_xxxxy, tr_xxzz_xxxxz, tr_xxzz_xxxyy, tr_xxzz_xxxyz, tr_xxzz_xxxzz, tr_xxzz_xxy, tr_xxzz_xxyyy, tr_xxzz_xxyyz, tr_xxzz_xxyzz, tr_xxzz_xxz, tr_xxzz_xxzzz, tr_xxzz_xyy, tr_xxzz_xyyyy, tr_xxzz_xyyyz, tr_xxzz_xyyzz, tr_xxzz_xyz, tr_xxzz_xyzzz, tr_xxzz_xzz, tr_xxzz_xzzzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxzz_xxxx[i] = 2.0 * tr_xxxzz_xxxx[i] * tbe_0 + 2.0 * tr_xxzz_xxxxx[i] * tke_0 - 2.0 * tr_xzz_xxxx[i] - 4.0 * tr_xxzz_xxx[i];

        tr_0_0_x_xxzz_xxxy[i] = 2.0 * tr_xxxzz_xxxy[i] * tbe_0 + 2.0 * tr_xxzz_xxxxy[i] * tke_0 - 2.0 * tr_xzz_xxxy[i] - 3.0 * tr_xxzz_xxy[i];

        tr_0_0_x_xxzz_xxxz[i] = 2.0 * tr_xxxzz_xxxz[i] * tbe_0 + 2.0 * tr_xxzz_xxxxz[i] * tke_0 - 2.0 * tr_xzz_xxxz[i] - 3.0 * tr_xxzz_xxz[i];

        tr_0_0_x_xxzz_xxyy[i] = 2.0 * tr_xxxzz_xxyy[i] * tbe_0 + 2.0 * tr_xxzz_xxxyy[i] * tke_0 - 2.0 * tr_xzz_xxyy[i] - 2.0 * tr_xxzz_xyy[i];

        tr_0_0_x_xxzz_xxyz[i] = 2.0 * tr_xxxzz_xxyz[i] * tbe_0 + 2.0 * tr_xxzz_xxxyz[i] * tke_0 - 2.0 * tr_xzz_xxyz[i] - 2.0 * tr_xxzz_xyz[i];

        tr_0_0_x_xxzz_xxzz[i] = 2.0 * tr_xxxzz_xxzz[i] * tbe_0 + 2.0 * tr_xxzz_xxxzz[i] * tke_0 - 2.0 * tr_xzz_xxzz[i] - 2.0 * tr_xxzz_xzz[i];

        tr_0_0_x_xxzz_xyyy[i] = 2.0 * tr_xxxzz_xyyy[i] * tbe_0 + 2.0 * tr_xxzz_xxyyy[i] * tke_0 - 2.0 * tr_xzz_xyyy[i] - tr_xxzz_yyy[i];

        tr_0_0_x_xxzz_xyyz[i] = 2.0 * tr_xxxzz_xyyz[i] * tbe_0 + 2.0 * tr_xxzz_xxyyz[i] * tke_0 - 2.0 * tr_xzz_xyyz[i] - tr_xxzz_yyz[i];

        tr_0_0_x_xxzz_xyzz[i] = 2.0 * tr_xxxzz_xyzz[i] * tbe_0 + 2.0 * tr_xxzz_xxyzz[i] * tke_0 - 2.0 * tr_xzz_xyzz[i] - tr_xxzz_yzz[i];

        tr_0_0_x_xxzz_xzzz[i] = 2.0 * tr_xxxzz_xzzz[i] * tbe_0 + 2.0 * tr_xxzz_xxzzz[i] * tke_0 - 2.0 * tr_xzz_xzzz[i] - tr_xxzz_zzz[i];

        tr_0_0_x_xxzz_yyyy[i] = 2.0 * tr_xxxzz_yyyy[i] * tbe_0 + 2.0 * tr_xxzz_xyyyy[i] * tke_0 - 2.0 * tr_xzz_yyyy[i];

        tr_0_0_x_xxzz_yyyz[i] = 2.0 * tr_xxxzz_yyyz[i] * tbe_0 + 2.0 * tr_xxzz_xyyyz[i] * tke_0 - 2.0 * tr_xzz_yyyz[i];

        tr_0_0_x_xxzz_yyzz[i] = 2.0 * tr_xxxzz_yyzz[i] * tbe_0 + 2.0 * tr_xxzz_xyyzz[i] * tke_0 - 2.0 * tr_xzz_yyzz[i];

        tr_0_0_x_xxzz_yzzz[i] = 2.0 * tr_xxxzz_yzzz[i] * tbe_0 + 2.0 * tr_xxzz_xyzzz[i] * tke_0 - 2.0 * tr_xzz_yzzz[i];

        tr_0_0_x_xxzz_zzzz[i] = 2.0 * tr_xxxzz_zzzz[i] * tbe_0 + 2.0 * tr_xxzz_xzzzz[i] * tke_0 - 2.0 * tr_xzz_zzzz[i];
    }

    // Set up 90-105 components of targeted buffer : GG

    auto tr_0_0_x_xyyy_xxxx = pbuffer.data(idx_op_geom_010_gg + 90);

    auto tr_0_0_x_xyyy_xxxy = pbuffer.data(idx_op_geom_010_gg + 91);

    auto tr_0_0_x_xyyy_xxxz = pbuffer.data(idx_op_geom_010_gg + 92);

    auto tr_0_0_x_xyyy_xxyy = pbuffer.data(idx_op_geom_010_gg + 93);

    auto tr_0_0_x_xyyy_xxyz = pbuffer.data(idx_op_geom_010_gg + 94);

    auto tr_0_0_x_xyyy_xxzz = pbuffer.data(idx_op_geom_010_gg + 95);

    auto tr_0_0_x_xyyy_xyyy = pbuffer.data(idx_op_geom_010_gg + 96);

    auto tr_0_0_x_xyyy_xyyz = pbuffer.data(idx_op_geom_010_gg + 97);

    auto tr_0_0_x_xyyy_xyzz = pbuffer.data(idx_op_geom_010_gg + 98);

    auto tr_0_0_x_xyyy_xzzz = pbuffer.data(idx_op_geom_010_gg + 99);

    auto tr_0_0_x_xyyy_yyyy = pbuffer.data(idx_op_geom_010_gg + 100);

    auto tr_0_0_x_xyyy_yyyz = pbuffer.data(idx_op_geom_010_gg + 101);

    auto tr_0_0_x_xyyy_yyzz = pbuffer.data(idx_op_geom_010_gg + 102);

    auto tr_0_0_x_xyyy_yzzz = pbuffer.data(idx_op_geom_010_gg + 103);

    auto tr_0_0_x_xyyy_zzzz = pbuffer.data(idx_op_geom_010_gg + 104);

    #pragma omp simd aligned(tr_0_0_x_xyyy_xxxx, tr_0_0_x_xyyy_xxxy, tr_0_0_x_xyyy_xxxz, tr_0_0_x_xyyy_xxyy, tr_0_0_x_xyyy_xxyz, tr_0_0_x_xyyy_xxzz, tr_0_0_x_xyyy_xyyy, tr_0_0_x_xyyy_xyyz, tr_0_0_x_xyyy_xyzz, tr_0_0_x_xyyy_xzzz, tr_0_0_x_xyyy_yyyy, tr_0_0_x_xyyy_yyyz, tr_0_0_x_xyyy_yyzz, tr_0_0_x_xyyy_yzzz, tr_0_0_x_xyyy_zzzz, tr_xxyyy_xxxx, tr_xxyyy_xxxy, tr_xxyyy_xxxz, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xzzz, tr_xxyyy_yyyy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yzzz, tr_xxyyy_zzzz, tr_xyyy_xxx, tr_xyyy_xxxxx, tr_xyyy_xxxxy, tr_xyyy_xxxxz, tr_xyyy_xxxyy, tr_xyyy_xxxyz, tr_xyyy_xxxzz, tr_xyyy_xxy, tr_xyyy_xxyyy, tr_xyyy_xxyyz, tr_xyyy_xxyzz, tr_xyyy_xxz, tr_xyyy_xxzzz, tr_xyyy_xyy, tr_xyyy_xyyyy, tr_xyyy_xyyyz, tr_xyyy_xyyzz, tr_xyyy_xyz, tr_xyyy_xyzzz, tr_xyyy_xzz, tr_xyyy_xzzzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyyy_xxxx[i] = 2.0 * tr_xxyyy_xxxx[i] * tbe_0 + 2.0 * tr_xyyy_xxxxx[i] * tke_0 - tr_yyy_xxxx[i] - 4.0 * tr_xyyy_xxx[i];

        tr_0_0_x_xyyy_xxxy[i] = 2.0 * tr_xxyyy_xxxy[i] * tbe_0 + 2.0 * tr_xyyy_xxxxy[i] * tke_0 - tr_yyy_xxxy[i] - 3.0 * tr_xyyy_xxy[i];

        tr_0_0_x_xyyy_xxxz[i] = 2.0 * tr_xxyyy_xxxz[i] * tbe_0 + 2.0 * tr_xyyy_xxxxz[i] * tke_0 - tr_yyy_xxxz[i] - 3.0 * tr_xyyy_xxz[i];

        tr_0_0_x_xyyy_xxyy[i] = 2.0 * tr_xxyyy_xxyy[i] * tbe_0 + 2.0 * tr_xyyy_xxxyy[i] * tke_0 - tr_yyy_xxyy[i] - 2.0 * tr_xyyy_xyy[i];

        tr_0_0_x_xyyy_xxyz[i] = 2.0 * tr_xxyyy_xxyz[i] * tbe_0 + 2.0 * tr_xyyy_xxxyz[i] * tke_0 - tr_yyy_xxyz[i] - 2.0 * tr_xyyy_xyz[i];

        tr_0_0_x_xyyy_xxzz[i] = 2.0 * tr_xxyyy_xxzz[i] * tbe_0 + 2.0 * tr_xyyy_xxxzz[i] * tke_0 - tr_yyy_xxzz[i] - 2.0 * tr_xyyy_xzz[i];

        tr_0_0_x_xyyy_xyyy[i] = 2.0 * tr_xxyyy_xyyy[i] * tbe_0 + 2.0 * tr_xyyy_xxyyy[i] * tke_0 - tr_yyy_xyyy[i] - tr_xyyy_yyy[i];

        tr_0_0_x_xyyy_xyyz[i] = 2.0 * tr_xxyyy_xyyz[i] * tbe_0 + 2.0 * tr_xyyy_xxyyz[i] * tke_0 - tr_yyy_xyyz[i] - tr_xyyy_yyz[i];

        tr_0_0_x_xyyy_xyzz[i] = 2.0 * tr_xxyyy_xyzz[i] * tbe_0 + 2.0 * tr_xyyy_xxyzz[i] * tke_0 - tr_yyy_xyzz[i] - tr_xyyy_yzz[i];

        tr_0_0_x_xyyy_xzzz[i] = 2.0 * tr_xxyyy_xzzz[i] * tbe_0 + 2.0 * tr_xyyy_xxzzz[i] * tke_0 - tr_yyy_xzzz[i] - tr_xyyy_zzz[i];

        tr_0_0_x_xyyy_yyyy[i] = 2.0 * tr_xxyyy_yyyy[i] * tbe_0 + 2.0 * tr_xyyy_xyyyy[i] * tke_0 - tr_yyy_yyyy[i];

        tr_0_0_x_xyyy_yyyz[i] = 2.0 * tr_xxyyy_yyyz[i] * tbe_0 + 2.0 * tr_xyyy_xyyyz[i] * tke_0 - tr_yyy_yyyz[i];

        tr_0_0_x_xyyy_yyzz[i] = 2.0 * tr_xxyyy_yyzz[i] * tbe_0 + 2.0 * tr_xyyy_xyyzz[i] * tke_0 - tr_yyy_yyzz[i];

        tr_0_0_x_xyyy_yzzz[i] = 2.0 * tr_xxyyy_yzzz[i] * tbe_0 + 2.0 * tr_xyyy_xyzzz[i] * tke_0 - tr_yyy_yzzz[i];

        tr_0_0_x_xyyy_zzzz[i] = 2.0 * tr_xxyyy_zzzz[i] * tbe_0 + 2.0 * tr_xyyy_xzzzz[i] * tke_0 - tr_yyy_zzzz[i];
    }

    // Set up 105-120 components of targeted buffer : GG

    auto tr_0_0_x_xyyz_xxxx = pbuffer.data(idx_op_geom_010_gg + 105);

    auto tr_0_0_x_xyyz_xxxy = pbuffer.data(idx_op_geom_010_gg + 106);

    auto tr_0_0_x_xyyz_xxxz = pbuffer.data(idx_op_geom_010_gg + 107);

    auto tr_0_0_x_xyyz_xxyy = pbuffer.data(idx_op_geom_010_gg + 108);

    auto tr_0_0_x_xyyz_xxyz = pbuffer.data(idx_op_geom_010_gg + 109);

    auto tr_0_0_x_xyyz_xxzz = pbuffer.data(idx_op_geom_010_gg + 110);

    auto tr_0_0_x_xyyz_xyyy = pbuffer.data(idx_op_geom_010_gg + 111);

    auto tr_0_0_x_xyyz_xyyz = pbuffer.data(idx_op_geom_010_gg + 112);

    auto tr_0_0_x_xyyz_xyzz = pbuffer.data(idx_op_geom_010_gg + 113);

    auto tr_0_0_x_xyyz_xzzz = pbuffer.data(idx_op_geom_010_gg + 114);

    auto tr_0_0_x_xyyz_yyyy = pbuffer.data(idx_op_geom_010_gg + 115);

    auto tr_0_0_x_xyyz_yyyz = pbuffer.data(idx_op_geom_010_gg + 116);

    auto tr_0_0_x_xyyz_yyzz = pbuffer.data(idx_op_geom_010_gg + 117);

    auto tr_0_0_x_xyyz_yzzz = pbuffer.data(idx_op_geom_010_gg + 118);

    auto tr_0_0_x_xyyz_zzzz = pbuffer.data(idx_op_geom_010_gg + 119);

    #pragma omp simd aligned(tr_0_0_x_xyyz_xxxx, tr_0_0_x_xyyz_xxxy, tr_0_0_x_xyyz_xxxz, tr_0_0_x_xyyz_xxyy, tr_0_0_x_xyyz_xxyz, tr_0_0_x_xyyz_xxzz, tr_0_0_x_xyyz_xyyy, tr_0_0_x_xyyz_xyyz, tr_0_0_x_xyyz_xyzz, tr_0_0_x_xyyz_xzzz, tr_0_0_x_xyyz_yyyy, tr_0_0_x_xyyz_yyyz, tr_0_0_x_xyyz_yyzz, tr_0_0_x_xyyz_yzzz, tr_0_0_x_xyyz_zzzz, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxx, tr_xyyz_xxxxy, tr_xyyz_xxxxz, tr_xyyz_xxxyy, tr_xyyz_xxxyz, tr_xyyz_xxxzz, tr_xyyz_xxy, tr_xyyz_xxyyy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xxzzz, tr_xyyz_xyy, tr_xyyz_xyyyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_xzzzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyyz_xxxx[i] = 2.0 * tr_xxyyz_xxxx[i] * tbe_0 + 2.0 * tr_xyyz_xxxxx[i] * tke_0 - tr_yyz_xxxx[i] - 4.0 * tr_xyyz_xxx[i];

        tr_0_0_x_xyyz_xxxy[i] = 2.0 * tr_xxyyz_xxxy[i] * tbe_0 + 2.0 * tr_xyyz_xxxxy[i] * tke_0 - tr_yyz_xxxy[i] - 3.0 * tr_xyyz_xxy[i];

        tr_0_0_x_xyyz_xxxz[i] = 2.0 * tr_xxyyz_xxxz[i] * tbe_0 + 2.0 * tr_xyyz_xxxxz[i] * tke_0 - tr_yyz_xxxz[i] - 3.0 * tr_xyyz_xxz[i];

        tr_0_0_x_xyyz_xxyy[i] = 2.0 * tr_xxyyz_xxyy[i] * tbe_0 + 2.0 * tr_xyyz_xxxyy[i] * tke_0 - tr_yyz_xxyy[i] - 2.0 * tr_xyyz_xyy[i];

        tr_0_0_x_xyyz_xxyz[i] = 2.0 * tr_xxyyz_xxyz[i] * tbe_0 + 2.0 * tr_xyyz_xxxyz[i] * tke_0 - tr_yyz_xxyz[i] - 2.0 * tr_xyyz_xyz[i];

        tr_0_0_x_xyyz_xxzz[i] = 2.0 * tr_xxyyz_xxzz[i] * tbe_0 + 2.0 * tr_xyyz_xxxzz[i] * tke_0 - tr_yyz_xxzz[i] - 2.0 * tr_xyyz_xzz[i];

        tr_0_0_x_xyyz_xyyy[i] = 2.0 * tr_xxyyz_xyyy[i] * tbe_0 + 2.0 * tr_xyyz_xxyyy[i] * tke_0 - tr_yyz_xyyy[i] - tr_xyyz_yyy[i];

        tr_0_0_x_xyyz_xyyz[i] = 2.0 * tr_xxyyz_xyyz[i] * tbe_0 + 2.0 * tr_xyyz_xxyyz[i] * tke_0 - tr_yyz_xyyz[i] - tr_xyyz_yyz[i];

        tr_0_0_x_xyyz_xyzz[i] = 2.0 * tr_xxyyz_xyzz[i] * tbe_0 + 2.0 * tr_xyyz_xxyzz[i] * tke_0 - tr_yyz_xyzz[i] - tr_xyyz_yzz[i];

        tr_0_0_x_xyyz_xzzz[i] = 2.0 * tr_xxyyz_xzzz[i] * tbe_0 + 2.0 * tr_xyyz_xxzzz[i] * tke_0 - tr_yyz_xzzz[i] - tr_xyyz_zzz[i];

        tr_0_0_x_xyyz_yyyy[i] = 2.0 * tr_xxyyz_yyyy[i] * tbe_0 + 2.0 * tr_xyyz_xyyyy[i] * tke_0 - tr_yyz_yyyy[i];

        tr_0_0_x_xyyz_yyyz[i] = 2.0 * tr_xxyyz_yyyz[i] * tbe_0 + 2.0 * tr_xyyz_xyyyz[i] * tke_0 - tr_yyz_yyyz[i];

        tr_0_0_x_xyyz_yyzz[i] = 2.0 * tr_xxyyz_yyzz[i] * tbe_0 + 2.0 * tr_xyyz_xyyzz[i] * tke_0 - tr_yyz_yyzz[i];

        tr_0_0_x_xyyz_yzzz[i] = 2.0 * tr_xxyyz_yzzz[i] * tbe_0 + 2.0 * tr_xyyz_xyzzz[i] * tke_0 - tr_yyz_yzzz[i];

        tr_0_0_x_xyyz_zzzz[i] = 2.0 * tr_xxyyz_zzzz[i] * tbe_0 + 2.0 * tr_xyyz_xzzzz[i] * tke_0 - tr_yyz_zzzz[i];
    }

    // Set up 120-135 components of targeted buffer : GG

    auto tr_0_0_x_xyzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 120);

    auto tr_0_0_x_xyzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 121);

    auto tr_0_0_x_xyzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 122);

    auto tr_0_0_x_xyzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 123);

    auto tr_0_0_x_xyzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 124);

    auto tr_0_0_x_xyzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 125);

    auto tr_0_0_x_xyzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 126);

    auto tr_0_0_x_xyzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 127);

    auto tr_0_0_x_xyzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 128);

    auto tr_0_0_x_xyzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 129);

    auto tr_0_0_x_xyzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 130);

    auto tr_0_0_x_xyzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 131);

    auto tr_0_0_x_xyzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 132);

    auto tr_0_0_x_xyzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 133);

    auto tr_0_0_x_xyzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 134);

    #pragma omp simd aligned(tr_0_0_x_xyzz_xxxx, tr_0_0_x_xyzz_xxxy, tr_0_0_x_xyzz_xxxz, tr_0_0_x_xyzz_xxyy, tr_0_0_x_xyzz_xxyz, tr_0_0_x_xyzz_xxzz, tr_0_0_x_xyzz_xyyy, tr_0_0_x_xyzz_xyyz, tr_0_0_x_xyzz_xyzz, tr_0_0_x_xyzz_xzzz, tr_0_0_x_xyzz_yyyy, tr_0_0_x_xyzz_yyyz, tr_0_0_x_xyzz_yyzz, tr_0_0_x_xyzz_yzzz, tr_0_0_x_xyzz_zzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxx, tr_xyzz_xxxxy, tr_xyzz_xxxxz, tr_xyzz_xxxyy, tr_xyzz_xxxyz, tr_xyzz_xxxzz, tr_xyzz_xxy, tr_xyzz_xxyyy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xxzzz, tr_xyzz_xyy, tr_xyzz_xyyyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_xzzzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyzz_xxxx[i] = 2.0 * tr_xxyzz_xxxx[i] * tbe_0 + 2.0 * tr_xyzz_xxxxx[i] * tke_0 - tr_yzz_xxxx[i] - 4.0 * tr_xyzz_xxx[i];

        tr_0_0_x_xyzz_xxxy[i] = 2.0 * tr_xxyzz_xxxy[i] * tbe_0 + 2.0 * tr_xyzz_xxxxy[i] * tke_0 - tr_yzz_xxxy[i] - 3.0 * tr_xyzz_xxy[i];

        tr_0_0_x_xyzz_xxxz[i] = 2.0 * tr_xxyzz_xxxz[i] * tbe_0 + 2.0 * tr_xyzz_xxxxz[i] * tke_0 - tr_yzz_xxxz[i] - 3.0 * tr_xyzz_xxz[i];

        tr_0_0_x_xyzz_xxyy[i] = 2.0 * tr_xxyzz_xxyy[i] * tbe_0 + 2.0 * tr_xyzz_xxxyy[i] * tke_0 - tr_yzz_xxyy[i] - 2.0 * tr_xyzz_xyy[i];

        tr_0_0_x_xyzz_xxyz[i] = 2.0 * tr_xxyzz_xxyz[i] * tbe_0 + 2.0 * tr_xyzz_xxxyz[i] * tke_0 - tr_yzz_xxyz[i] - 2.0 * tr_xyzz_xyz[i];

        tr_0_0_x_xyzz_xxzz[i] = 2.0 * tr_xxyzz_xxzz[i] * tbe_0 + 2.0 * tr_xyzz_xxxzz[i] * tke_0 - tr_yzz_xxzz[i] - 2.0 * tr_xyzz_xzz[i];

        tr_0_0_x_xyzz_xyyy[i] = 2.0 * tr_xxyzz_xyyy[i] * tbe_0 + 2.0 * tr_xyzz_xxyyy[i] * tke_0 - tr_yzz_xyyy[i] - tr_xyzz_yyy[i];

        tr_0_0_x_xyzz_xyyz[i] = 2.0 * tr_xxyzz_xyyz[i] * tbe_0 + 2.0 * tr_xyzz_xxyyz[i] * tke_0 - tr_yzz_xyyz[i] - tr_xyzz_yyz[i];

        tr_0_0_x_xyzz_xyzz[i] = 2.0 * tr_xxyzz_xyzz[i] * tbe_0 + 2.0 * tr_xyzz_xxyzz[i] * tke_0 - tr_yzz_xyzz[i] - tr_xyzz_yzz[i];

        tr_0_0_x_xyzz_xzzz[i] = 2.0 * tr_xxyzz_xzzz[i] * tbe_0 + 2.0 * tr_xyzz_xxzzz[i] * tke_0 - tr_yzz_xzzz[i] - tr_xyzz_zzz[i];

        tr_0_0_x_xyzz_yyyy[i] = 2.0 * tr_xxyzz_yyyy[i] * tbe_0 + 2.0 * tr_xyzz_xyyyy[i] * tke_0 - tr_yzz_yyyy[i];

        tr_0_0_x_xyzz_yyyz[i] = 2.0 * tr_xxyzz_yyyz[i] * tbe_0 + 2.0 * tr_xyzz_xyyyz[i] * tke_0 - tr_yzz_yyyz[i];

        tr_0_0_x_xyzz_yyzz[i] = 2.0 * tr_xxyzz_yyzz[i] * tbe_0 + 2.0 * tr_xyzz_xyyzz[i] * tke_0 - tr_yzz_yyzz[i];

        tr_0_0_x_xyzz_yzzz[i] = 2.0 * tr_xxyzz_yzzz[i] * tbe_0 + 2.0 * tr_xyzz_xyzzz[i] * tke_0 - tr_yzz_yzzz[i];

        tr_0_0_x_xyzz_zzzz[i] = 2.0 * tr_xxyzz_zzzz[i] * tbe_0 + 2.0 * tr_xyzz_xzzzz[i] * tke_0 - tr_yzz_zzzz[i];
    }

    // Set up 135-150 components of targeted buffer : GG

    auto tr_0_0_x_xzzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 135);

    auto tr_0_0_x_xzzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 136);

    auto tr_0_0_x_xzzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 137);

    auto tr_0_0_x_xzzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 138);

    auto tr_0_0_x_xzzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 139);

    auto tr_0_0_x_xzzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 140);

    auto tr_0_0_x_xzzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 141);

    auto tr_0_0_x_xzzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 142);

    auto tr_0_0_x_xzzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 143);

    auto tr_0_0_x_xzzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 144);

    auto tr_0_0_x_xzzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 145);

    auto tr_0_0_x_xzzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 146);

    auto tr_0_0_x_xzzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 147);

    auto tr_0_0_x_xzzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 148);

    auto tr_0_0_x_xzzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 149);

    #pragma omp simd aligned(tr_0_0_x_xzzz_xxxx, tr_0_0_x_xzzz_xxxy, tr_0_0_x_xzzz_xxxz, tr_0_0_x_xzzz_xxyy, tr_0_0_x_xzzz_xxyz, tr_0_0_x_xzzz_xxzz, tr_0_0_x_xzzz_xyyy, tr_0_0_x_xzzz_xyyz, tr_0_0_x_xzzz_xyzz, tr_0_0_x_xzzz_xzzz, tr_0_0_x_xzzz_yyyy, tr_0_0_x_xzzz_yyyz, tr_0_0_x_xzzz_yyzz, tr_0_0_x_xzzz_yzzz, tr_0_0_x_xzzz_zzzz, tr_xxzzz_xxxx, tr_xxzzz_xxxy, tr_xxzzz_xxxz, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xzzz, tr_xxzzz_yyyy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yzzz, tr_xxzzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxxxx, tr_xzzz_xxxxy, tr_xzzz_xxxxz, tr_xzzz_xxxyy, tr_xzzz_xxxyz, tr_xzzz_xxxzz, tr_xzzz_xxy, tr_xzzz_xxyyy, tr_xzzz_xxyyz, tr_xzzz_xxyzz, tr_xzzz_xxz, tr_xzzz_xxzzz, tr_xzzz_xyy, tr_xzzz_xyyyy, tr_xzzz_xyyyz, tr_xzzz_xyyzz, tr_xzzz_xyz, tr_xzzz_xyzzz, tr_xzzz_xzz, tr_xzzz_xzzzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xzzz_xxxx[i] = 2.0 * tr_xxzzz_xxxx[i] * tbe_0 + 2.0 * tr_xzzz_xxxxx[i] * tke_0 - tr_zzz_xxxx[i] - 4.0 * tr_xzzz_xxx[i];

        tr_0_0_x_xzzz_xxxy[i] = 2.0 * tr_xxzzz_xxxy[i] * tbe_0 + 2.0 * tr_xzzz_xxxxy[i] * tke_0 - tr_zzz_xxxy[i] - 3.0 * tr_xzzz_xxy[i];

        tr_0_0_x_xzzz_xxxz[i] = 2.0 * tr_xxzzz_xxxz[i] * tbe_0 + 2.0 * tr_xzzz_xxxxz[i] * tke_0 - tr_zzz_xxxz[i] - 3.0 * tr_xzzz_xxz[i];

        tr_0_0_x_xzzz_xxyy[i] = 2.0 * tr_xxzzz_xxyy[i] * tbe_0 + 2.0 * tr_xzzz_xxxyy[i] * tke_0 - tr_zzz_xxyy[i] - 2.0 * tr_xzzz_xyy[i];

        tr_0_0_x_xzzz_xxyz[i] = 2.0 * tr_xxzzz_xxyz[i] * tbe_0 + 2.0 * tr_xzzz_xxxyz[i] * tke_0 - tr_zzz_xxyz[i] - 2.0 * tr_xzzz_xyz[i];

        tr_0_0_x_xzzz_xxzz[i] = 2.0 * tr_xxzzz_xxzz[i] * tbe_0 + 2.0 * tr_xzzz_xxxzz[i] * tke_0 - tr_zzz_xxzz[i] - 2.0 * tr_xzzz_xzz[i];

        tr_0_0_x_xzzz_xyyy[i] = 2.0 * tr_xxzzz_xyyy[i] * tbe_0 + 2.0 * tr_xzzz_xxyyy[i] * tke_0 - tr_zzz_xyyy[i] - tr_xzzz_yyy[i];

        tr_0_0_x_xzzz_xyyz[i] = 2.0 * tr_xxzzz_xyyz[i] * tbe_0 + 2.0 * tr_xzzz_xxyyz[i] * tke_0 - tr_zzz_xyyz[i] - tr_xzzz_yyz[i];

        tr_0_0_x_xzzz_xyzz[i] = 2.0 * tr_xxzzz_xyzz[i] * tbe_0 + 2.0 * tr_xzzz_xxyzz[i] * tke_0 - tr_zzz_xyzz[i] - tr_xzzz_yzz[i];

        tr_0_0_x_xzzz_xzzz[i] = 2.0 * tr_xxzzz_xzzz[i] * tbe_0 + 2.0 * tr_xzzz_xxzzz[i] * tke_0 - tr_zzz_xzzz[i] - tr_xzzz_zzz[i];

        tr_0_0_x_xzzz_yyyy[i] = 2.0 * tr_xxzzz_yyyy[i] * tbe_0 + 2.0 * tr_xzzz_xyyyy[i] * tke_0 - tr_zzz_yyyy[i];

        tr_0_0_x_xzzz_yyyz[i] = 2.0 * tr_xxzzz_yyyz[i] * tbe_0 + 2.0 * tr_xzzz_xyyyz[i] * tke_0 - tr_zzz_yyyz[i];

        tr_0_0_x_xzzz_yyzz[i] = 2.0 * tr_xxzzz_yyzz[i] * tbe_0 + 2.0 * tr_xzzz_xyyzz[i] * tke_0 - tr_zzz_yyzz[i];

        tr_0_0_x_xzzz_yzzz[i] = 2.0 * tr_xxzzz_yzzz[i] * tbe_0 + 2.0 * tr_xzzz_xyzzz[i] * tke_0 - tr_zzz_yzzz[i];

        tr_0_0_x_xzzz_zzzz[i] = 2.0 * tr_xxzzz_zzzz[i] * tbe_0 + 2.0 * tr_xzzz_xzzzz[i] * tke_0 - tr_zzz_zzzz[i];
    }

    // Set up 150-165 components of targeted buffer : GG

    auto tr_0_0_x_yyyy_xxxx = pbuffer.data(idx_op_geom_010_gg + 150);

    auto tr_0_0_x_yyyy_xxxy = pbuffer.data(idx_op_geom_010_gg + 151);

    auto tr_0_0_x_yyyy_xxxz = pbuffer.data(idx_op_geom_010_gg + 152);

    auto tr_0_0_x_yyyy_xxyy = pbuffer.data(idx_op_geom_010_gg + 153);

    auto tr_0_0_x_yyyy_xxyz = pbuffer.data(idx_op_geom_010_gg + 154);

    auto tr_0_0_x_yyyy_xxzz = pbuffer.data(idx_op_geom_010_gg + 155);

    auto tr_0_0_x_yyyy_xyyy = pbuffer.data(idx_op_geom_010_gg + 156);

    auto tr_0_0_x_yyyy_xyyz = pbuffer.data(idx_op_geom_010_gg + 157);

    auto tr_0_0_x_yyyy_xyzz = pbuffer.data(idx_op_geom_010_gg + 158);

    auto tr_0_0_x_yyyy_xzzz = pbuffer.data(idx_op_geom_010_gg + 159);

    auto tr_0_0_x_yyyy_yyyy = pbuffer.data(idx_op_geom_010_gg + 160);

    auto tr_0_0_x_yyyy_yyyz = pbuffer.data(idx_op_geom_010_gg + 161);

    auto tr_0_0_x_yyyy_yyzz = pbuffer.data(idx_op_geom_010_gg + 162);

    auto tr_0_0_x_yyyy_yzzz = pbuffer.data(idx_op_geom_010_gg + 163);

    auto tr_0_0_x_yyyy_zzzz = pbuffer.data(idx_op_geom_010_gg + 164);

    #pragma omp simd aligned(tr_0_0_x_yyyy_xxxx, tr_0_0_x_yyyy_xxxy, tr_0_0_x_yyyy_xxxz, tr_0_0_x_yyyy_xxyy, tr_0_0_x_yyyy_xxyz, tr_0_0_x_yyyy_xxzz, tr_0_0_x_yyyy_xyyy, tr_0_0_x_yyyy_xyyz, tr_0_0_x_yyyy_xyzz, tr_0_0_x_yyyy_xzzz, tr_0_0_x_yyyy_yyyy, tr_0_0_x_yyyy_yyyz, tr_0_0_x_yyyy_yyzz, tr_0_0_x_yyyy_yzzz, tr_0_0_x_yyyy_zzzz, tr_xyyyy_xxxx, tr_xyyyy_xxxy, tr_xyyyy_xxxz, tr_xyyyy_xxyy, tr_xyyyy_xxyz, tr_xyyyy_xxzz, tr_xyyyy_xyyy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xzzz, tr_xyyyy_yyyy, tr_xyyyy_yyyz, tr_xyyyy_yyzz, tr_xyyyy_yzzz, tr_xyyyy_zzzz, tr_yyyy_xxx, tr_yyyy_xxxxx, tr_yyyy_xxxxy, tr_yyyy_xxxxz, tr_yyyy_xxxyy, tr_yyyy_xxxyz, tr_yyyy_xxxzz, tr_yyyy_xxy, tr_yyyy_xxyyy, tr_yyyy_xxyyz, tr_yyyy_xxyzz, tr_yyyy_xxz, tr_yyyy_xxzzz, tr_yyyy_xyy, tr_yyyy_xyyyy, tr_yyyy_xyyyz, tr_yyyy_xyyzz, tr_yyyy_xyz, tr_yyyy_xyzzz, tr_yyyy_xzz, tr_yyyy_xzzzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyyy_xxxx[i] = 2.0 * tr_xyyyy_xxxx[i] * tbe_0 + 2.0 * tr_yyyy_xxxxx[i] * tke_0 - 4.0 * tr_yyyy_xxx[i];

        tr_0_0_x_yyyy_xxxy[i] = 2.0 * tr_xyyyy_xxxy[i] * tbe_0 + 2.0 * tr_yyyy_xxxxy[i] * tke_0 - 3.0 * tr_yyyy_xxy[i];

        tr_0_0_x_yyyy_xxxz[i] = 2.0 * tr_xyyyy_xxxz[i] * tbe_0 + 2.0 * tr_yyyy_xxxxz[i] * tke_0 - 3.0 * tr_yyyy_xxz[i];

        tr_0_0_x_yyyy_xxyy[i] = 2.0 * tr_xyyyy_xxyy[i] * tbe_0 + 2.0 * tr_yyyy_xxxyy[i] * tke_0 - 2.0 * tr_yyyy_xyy[i];

        tr_0_0_x_yyyy_xxyz[i] = 2.0 * tr_xyyyy_xxyz[i] * tbe_0 + 2.0 * tr_yyyy_xxxyz[i] * tke_0 - 2.0 * tr_yyyy_xyz[i];

        tr_0_0_x_yyyy_xxzz[i] = 2.0 * tr_xyyyy_xxzz[i] * tbe_0 + 2.0 * tr_yyyy_xxxzz[i] * tke_0 - 2.0 * tr_yyyy_xzz[i];

        tr_0_0_x_yyyy_xyyy[i] = 2.0 * tr_xyyyy_xyyy[i] * tbe_0 + 2.0 * tr_yyyy_xxyyy[i] * tke_0 - tr_yyyy_yyy[i];

        tr_0_0_x_yyyy_xyyz[i] = 2.0 * tr_xyyyy_xyyz[i] * tbe_0 + 2.0 * tr_yyyy_xxyyz[i] * tke_0 - tr_yyyy_yyz[i];

        tr_0_0_x_yyyy_xyzz[i] = 2.0 * tr_xyyyy_xyzz[i] * tbe_0 + 2.0 * tr_yyyy_xxyzz[i] * tke_0 - tr_yyyy_yzz[i];

        tr_0_0_x_yyyy_xzzz[i] = 2.0 * tr_xyyyy_xzzz[i] * tbe_0 + 2.0 * tr_yyyy_xxzzz[i] * tke_0 - tr_yyyy_zzz[i];

        tr_0_0_x_yyyy_yyyy[i] = 2.0 * tr_xyyyy_yyyy[i] * tbe_0 + 2.0 * tr_yyyy_xyyyy[i] * tke_0;

        tr_0_0_x_yyyy_yyyz[i] = 2.0 * tr_xyyyy_yyyz[i] * tbe_0 + 2.0 * tr_yyyy_xyyyz[i] * tke_0;

        tr_0_0_x_yyyy_yyzz[i] = 2.0 * tr_xyyyy_yyzz[i] * tbe_0 + 2.0 * tr_yyyy_xyyzz[i] * tke_0;

        tr_0_0_x_yyyy_yzzz[i] = 2.0 * tr_xyyyy_yzzz[i] * tbe_0 + 2.0 * tr_yyyy_xyzzz[i] * tke_0;

        tr_0_0_x_yyyy_zzzz[i] = 2.0 * tr_xyyyy_zzzz[i] * tbe_0 + 2.0 * tr_yyyy_xzzzz[i] * tke_0;
    }

    // Set up 165-180 components of targeted buffer : GG

    auto tr_0_0_x_yyyz_xxxx = pbuffer.data(idx_op_geom_010_gg + 165);

    auto tr_0_0_x_yyyz_xxxy = pbuffer.data(idx_op_geom_010_gg + 166);

    auto tr_0_0_x_yyyz_xxxz = pbuffer.data(idx_op_geom_010_gg + 167);

    auto tr_0_0_x_yyyz_xxyy = pbuffer.data(idx_op_geom_010_gg + 168);

    auto tr_0_0_x_yyyz_xxyz = pbuffer.data(idx_op_geom_010_gg + 169);

    auto tr_0_0_x_yyyz_xxzz = pbuffer.data(idx_op_geom_010_gg + 170);

    auto tr_0_0_x_yyyz_xyyy = pbuffer.data(idx_op_geom_010_gg + 171);

    auto tr_0_0_x_yyyz_xyyz = pbuffer.data(idx_op_geom_010_gg + 172);

    auto tr_0_0_x_yyyz_xyzz = pbuffer.data(idx_op_geom_010_gg + 173);

    auto tr_0_0_x_yyyz_xzzz = pbuffer.data(idx_op_geom_010_gg + 174);

    auto tr_0_0_x_yyyz_yyyy = pbuffer.data(idx_op_geom_010_gg + 175);

    auto tr_0_0_x_yyyz_yyyz = pbuffer.data(idx_op_geom_010_gg + 176);

    auto tr_0_0_x_yyyz_yyzz = pbuffer.data(idx_op_geom_010_gg + 177);

    auto tr_0_0_x_yyyz_yzzz = pbuffer.data(idx_op_geom_010_gg + 178);

    auto tr_0_0_x_yyyz_zzzz = pbuffer.data(idx_op_geom_010_gg + 179);

    #pragma omp simd aligned(tr_0_0_x_yyyz_xxxx, tr_0_0_x_yyyz_xxxy, tr_0_0_x_yyyz_xxxz, tr_0_0_x_yyyz_xxyy, tr_0_0_x_yyyz_xxyz, tr_0_0_x_yyyz_xxzz, tr_0_0_x_yyyz_xyyy, tr_0_0_x_yyyz_xyyz, tr_0_0_x_yyyz_xyzz, tr_0_0_x_yyyz_xzzz, tr_0_0_x_yyyz_yyyy, tr_0_0_x_yyyz_yyyz, tr_0_0_x_yyyz_yyzz, tr_0_0_x_yyyz_yzzz, tr_0_0_x_yyyz_zzzz, tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, tr_yyyz_xxx, tr_yyyz_xxxxx, tr_yyyz_xxxxy, tr_yyyz_xxxxz, tr_yyyz_xxxyy, tr_yyyz_xxxyz, tr_yyyz_xxxzz, tr_yyyz_xxy, tr_yyyz_xxyyy, tr_yyyz_xxyyz, tr_yyyz_xxyzz, tr_yyyz_xxz, tr_yyyz_xxzzz, tr_yyyz_xyy, tr_yyyz_xyyyy, tr_yyyz_xyyyz, tr_yyyz_xyyzz, tr_yyyz_xyz, tr_yyyz_xyzzz, tr_yyyz_xzz, tr_yyyz_xzzzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyyz_xxxx[i] = 2.0 * tr_xyyyz_xxxx[i] * tbe_0 + 2.0 * tr_yyyz_xxxxx[i] * tke_0 - 4.0 * tr_yyyz_xxx[i];

        tr_0_0_x_yyyz_xxxy[i] = 2.0 * tr_xyyyz_xxxy[i] * tbe_0 + 2.0 * tr_yyyz_xxxxy[i] * tke_0 - 3.0 * tr_yyyz_xxy[i];

        tr_0_0_x_yyyz_xxxz[i] = 2.0 * tr_xyyyz_xxxz[i] * tbe_0 + 2.0 * tr_yyyz_xxxxz[i] * tke_0 - 3.0 * tr_yyyz_xxz[i];

        tr_0_0_x_yyyz_xxyy[i] = 2.0 * tr_xyyyz_xxyy[i] * tbe_0 + 2.0 * tr_yyyz_xxxyy[i] * tke_0 - 2.0 * tr_yyyz_xyy[i];

        tr_0_0_x_yyyz_xxyz[i] = 2.0 * tr_xyyyz_xxyz[i] * tbe_0 + 2.0 * tr_yyyz_xxxyz[i] * tke_0 - 2.0 * tr_yyyz_xyz[i];

        tr_0_0_x_yyyz_xxzz[i] = 2.0 * tr_xyyyz_xxzz[i] * tbe_0 + 2.0 * tr_yyyz_xxxzz[i] * tke_0 - 2.0 * tr_yyyz_xzz[i];

        tr_0_0_x_yyyz_xyyy[i] = 2.0 * tr_xyyyz_xyyy[i] * tbe_0 + 2.0 * tr_yyyz_xxyyy[i] * tke_0 - tr_yyyz_yyy[i];

        tr_0_0_x_yyyz_xyyz[i] = 2.0 * tr_xyyyz_xyyz[i] * tbe_0 + 2.0 * tr_yyyz_xxyyz[i] * tke_0 - tr_yyyz_yyz[i];

        tr_0_0_x_yyyz_xyzz[i] = 2.0 * tr_xyyyz_xyzz[i] * tbe_0 + 2.0 * tr_yyyz_xxyzz[i] * tke_0 - tr_yyyz_yzz[i];

        tr_0_0_x_yyyz_xzzz[i] = 2.0 * tr_xyyyz_xzzz[i] * tbe_0 + 2.0 * tr_yyyz_xxzzz[i] * tke_0 - tr_yyyz_zzz[i];

        tr_0_0_x_yyyz_yyyy[i] = 2.0 * tr_xyyyz_yyyy[i] * tbe_0 + 2.0 * tr_yyyz_xyyyy[i] * tke_0;

        tr_0_0_x_yyyz_yyyz[i] = 2.0 * tr_xyyyz_yyyz[i] * tbe_0 + 2.0 * tr_yyyz_xyyyz[i] * tke_0;

        tr_0_0_x_yyyz_yyzz[i] = 2.0 * tr_xyyyz_yyzz[i] * tbe_0 + 2.0 * tr_yyyz_xyyzz[i] * tke_0;

        tr_0_0_x_yyyz_yzzz[i] = 2.0 * tr_xyyyz_yzzz[i] * tbe_0 + 2.0 * tr_yyyz_xyzzz[i] * tke_0;

        tr_0_0_x_yyyz_zzzz[i] = 2.0 * tr_xyyyz_zzzz[i] * tbe_0 + 2.0 * tr_yyyz_xzzzz[i] * tke_0;
    }

    // Set up 180-195 components of targeted buffer : GG

    auto tr_0_0_x_yyzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 180);

    auto tr_0_0_x_yyzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 181);

    auto tr_0_0_x_yyzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 182);

    auto tr_0_0_x_yyzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 183);

    auto tr_0_0_x_yyzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 184);

    auto tr_0_0_x_yyzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 185);

    auto tr_0_0_x_yyzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 186);

    auto tr_0_0_x_yyzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 187);

    auto tr_0_0_x_yyzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 188);

    auto tr_0_0_x_yyzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 189);

    auto tr_0_0_x_yyzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 190);

    auto tr_0_0_x_yyzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 191);

    auto tr_0_0_x_yyzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 192);

    auto tr_0_0_x_yyzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 193);

    auto tr_0_0_x_yyzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 194);

    #pragma omp simd aligned(tr_0_0_x_yyzz_xxxx, tr_0_0_x_yyzz_xxxy, tr_0_0_x_yyzz_xxxz, tr_0_0_x_yyzz_xxyy, tr_0_0_x_yyzz_xxyz, tr_0_0_x_yyzz_xxzz, tr_0_0_x_yyzz_xyyy, tr_0_0_x_yyzz_xyyz, tr_0_0_x_yyzz_xyzz, tr_0_0_x_yyzz_xzzz, tr_0_0_x_yyzz_yyyy, tr_0_0_x_yyzz_yyyz, tr_0_0_x_yyzz_yyzz, tr_0_0_x_yyzz_yzzz, tr_0_0_x_yyzz_zzzz, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, tr_yyzz_xxx, tr_yyzz_xxxxx, tr_yyzz_xxxxy, tr_yyzz_xxxxz, tr_yyzz_xxxyy, tr_yyzz_xxxyz, tr_yyzz_xxxzz, tr_yyzz_xxy, tr_yyzz_xxyyy, tr_yyzz_xxyyz, tr_yyzz_xxyzz, tr_yyzz_xxz, tr_yyzz_xxzzz, tr_yyzz_xyy, tr_yyzz_xyyyy, tr_yyzz_xyyyz, tr_yyzz_xyyzz, tr_yyzz_xyz, tr_yyzz_xyzzz, tr_yyzz_xzz, tr_yyzz_xzzzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyzz_xxxx[i] = 2.0 * tr_xyyzz_xxxx[i] * tbe_0 + 2.0 * tr_yyzz_xxxxx[i] * tke_0 - 4.0 * tr_yyzz_xxx[i];

        tr_0_0_x_yyzz_xxxy[i] = 2.0 * tr_xyyzz_xxxy[i] * tbe_0 + 2.0 * tr_yyzz_xxxxy[i] * tke_0 - 3.0 * tr_yyzz_xxy[i];

        tr_0_0_x_yyzz_xxxz[i] = 2.0 * tr_xyyzz_xxxz[i] * tbe_0 + 2.0 * tr_yyzz_xxxxz[i] * tke_0 - 3.0 * tr_yyzz_xxz[i];

        tr_0_0_x_yyzz_xxyy[i] = 2.0 * tr_xyyzz_xxyy[i] * tbe_0 + 2.0 * tr_yyzz_xxxyy[i] * tke_0 - 2.0 * tr_yyzz_xyy[i];

        tr_0_0_x_yyzz_xxyz[i] = 2.0 * tr_xyyzz_xxyz[i] * tbe_0 + 2.0 * tr_yyzz_xxxyz[i] * tke_0 - 2.0 * tr_yyzz_xyz[i];

        tr_0_0_x_yyzz_xxzz[i] = 2.0 * tr_xyyzz_xxzz[i] * tbe_0 + 2.0 * tr_yyzz_xxxzz[i] * tke_0 - 2.0 * tr_yyzz_xzz[i];

        tr_0_0_x_yyzz_xyyy[i] = 2.0 * tr_xyyzz_xyyy[i] * tbe_0 + 2.0 * tr_yyzz_xxyyy[i] * tke_0 - tr_yyzz_yyy[i];

        tr_0_0_x_yyzz_xyyz[i] = 2.0 * tr_xyyzz_xyyz[i] * tbe_0 + 2.0 * tr_yyzz_xxyyz[i] * tke_0 - tr_yyzz_yyz[i];

        tr_0_0_x_yyzz_xyzz[i] = 2.0 * tr_xyyzz_xyzz[i] * tbe_0 + 2.0 * tr_yyzz_xxyzz[i] * tke_0 - tr_yyzz_yzz[i];

        tr_0_0_x_yyzz_xzzz[i] = 2.0 * tr_xyyzz_xzzz[i] * tbe_0 + 2.0 * tr_yyzz_xxzzz[i] * tke_0 - tr_yyzz_zzz[i];

        tr_0_0_x_yyzz_yyyy[i] = 2.0 * tr_xyyzz_yyyy[i] * tbe_0 + 2.0 * tr_yyzz_xyyyy[i] * tke_0;

        tr_0_0_x_yyzz_yyyz[i] = 2.0 * tr_xyyzz_yyyz[i] * tbe_0 + 2.0 * tr_yyzz_xyyyz[i] * tke_0;

        tr_0_0_x_yyzz_yyzz[i] = 2.0 * tr_xyyzz_yyzz[i] * tbe_0 + 2.0 * tr_yyzz_xyyzz[i] * tke_0;

        tr_0_0_x_yyzz_yzzz[i] = 2.0 * tr_xyyzz_yzzz[i] * tbe_0 + 2.0 * tr_yyzz_xyzzz[i] * tke_0;

        tr_0_0_x_yyzz_zzzz[i] = 2.0 * tr_xyyzz_zzzz[i] * tbe_0 + 2.0 * tr_yyzz_xzzzz[i] * tke_0;
    }

    // Set up 195-210 components of targeted buffer : GG

    auto tr_0_0_x_yzzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 195);

    auto tr_0_0_x_yzzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 196);

    auto tr_0_0_x_yzzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 197);

    auto tr_0_0_x_yzzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 198);

    auto tr_0_0_x_yzzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 199);

    auto tr_0_0_x_yzzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 200);

    auto tr_0_0_x_yzzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 201);

    auto tr_0_0_x_yzzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 202);

    auto tr_0_0_x_yzzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 203);

    auto tr_0_0_x_yzzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 204);

    auto tr_0_0_x_yzzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 205);

    auto tr_0_0_x_yzzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 206);

    auto tr_0_0_x_yzzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 207);

    auto tr_0_0_x_yzzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 208);

    auto tr_0_0_x_yzzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 209);

    #pragma omp simd aligned(tr_0_0_x_yzzz_xxxx, tr_0_0_x_yzzz_xxxy, tr_0_0_x_yzzz_xxxz, tr_0_0_x_yzzz_xxyy, tr_0_0_x_yzzz_xxyz, tr_0_0_x_yzzz_xxzz, tr_0_0_x_yzzz_xyyy, tr_0_0_x_yzzz_xyyz, tr_0_0_x_yzzz_xyzz, tr_0_0_x_yzzz_xzzz, tr_0_0_x_yzzz_yyyy, tr_0_0_x_yzzz_yyyz, tr_0_0_x_yzzz_yyzz, tr_0_0_x_yzzz_yzzz, tr_0_0_x_yzzz_zzzz, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxxxx, tr_yzzz_xxxxy, tr_yzzz_xxxxz, tr_yzzz_xxxyy, tr_yzzz_xxxyz, tr_yzzz_xxxzz, tr_yzzz_xxy, tr_yzzz_xxyyy, tr_yzzz_xxyyz, tr_yzzz_xxyzz, tr_yzzz_xxz, tr_yzzz_xxzzz, tr_yzzz_xyy, tr_yzzz_xyyyy, tr_yzzz_xyyyz, tr_yzzz_xyyzz, tr_yzzz_xyz, tr_yzzz_xyzzz, tr_yzzz_xzz, tr_yzzz_xzzzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yzzz_xxxx[i] = 2.0 * tr_xyzzz_xxxx[i] * tbe_0 + 2.0 * tr_yzzz_xxxxx[i] * tke_0 - 4.0 * tr_yzzz_xxx[i];

        tr_0_0_x_yzzz_xxxy[i] = 2.0 * tr_xyzzz_xxxy[i] * tbe_0 + 2.0 * tr_yzzz_xxxxy[i] * tke_0 - 3.0 * tr_yzzz_xxy[i];

        tr_0_0_x_yzzz_xxxz[i] = 2.0 * tr_xyzzz_xxxz[i] * tbe_0 + 2.0 * tr_yzzz_xxxxz[i] * tke_0 - 3.0 * tr_yzzz_xxz[i];

        tr_0_0_x_yzzz_xxyy[i] = 2.0 * tr_xyzzz_xxyy[i] * tbe_0 + 2.0 * tr_yzzz_xxxyy[i] * tke_0 - 2.0 * tr_yzzz_xyy[i];

        tr_0_0_x_yzzz_xxyz[i] = 2.0 * tr_xyzzz_xxyz[i] * tbe_0 + 2.0 * tr_yzzz_xxxyz[i] * tke_0 - 2.0 * tr_yzzz_xyz[i];

        tr_0_0_x_yzzz_xxzz[i] = 2.0 * tr_xyzzz_xxzz[i] * tbe_0 + 2.0 * tr_yzzz_xxxzz[i] * tke_0 - 2.0 * tr_yzzz_xzz[i];

        tr_0_0_x_yzzz_xyyy[i] = 2.0 * tr_xyzzz_xyyy[i] * tbe_0 + 2.0 * tr_yzzz_xxyyy[i] * tke_0 - tr_yzzz_yyy[i];

        tr_0_0_x_yzzz_xyyz[i] = 2.0 * tr_xyzzz_xyyz[i] * tbe_0 + 2.0 * tr_yzzz_xxyyz[i] * tke_0 - tr_yzzz_yyz[i];

        tr_0_0_x_yzzz_xyzz[i] = 2.0 * tr_xyzzz_xyzz[i] * tbe_0 + 2.0 * tr_yzzz_xxyzz[i] * tke_0 - tr_yzzz_yzz[i];

        tr_0_0_x_yzzz_xzzz[i] = 2.0 * tr_xyzzz_xzzz[i] * tbe_0 + 2.0 * tr_yzzz_xxzzz[i] * tke_0 - tr_yzzz_zzz[i];

        tr_0_0_x_yzzz_yyyy[i] = 2.0 * tr_xyzzz_yyyy[i] * tbe_0 + 2.0 * tr_yzzz_xyyyy[i] * tke_0;

        tr_0_0_x_yzzz_yyyz[i] = 2.0 * tr_xyzzz_yyyz[i] * tbe_0 + 2.0 * tr_yzzz_xyyyz[i] * tke_0;

        tr_0_0_x_yzzz_yyzz[i] = 2.0 * tr_xyzzz_yyzz[i] * tbe_0 + 2.0 * tr_yzzz_xyyzz[i] * tke_0;

        tr_0_0_x_yzzz_yzzz[i] = 2.0 * tr_xyzzz_yzzz[i] * tbe_0 + 2.0 * tr_yzzz_xyzzz[i] * tke_0;

        tr_0_0_x_yzzz_zzzz[i] = 2.0 * tr_xyzzz_zzzz[i] * tbe_0 + 2.0 * tr_yzzz_xzzzz[i] * tke_0;
    }

    // Set up 210-225 components of targeted buffer : GG

    auto tr_0_0_x_zzzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 210);

    auto tr_0_0_x_zzzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 211);

    auto tr_0_0_x_zzzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 212);

    auto tr_0_0_x_zzzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 213);

    auto tr_0_0_x_zzzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 214);

    auto tr_0_0_x_zzzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 215);

    auto tr_0_0_x_zzzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 216);

    auto tr_0_0_x_zzzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 217);

    auto tr_0_0_x_zzzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 218);

    auto tr_0_0_x_zzzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 219);

    auto tr_0_0_x_zzzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 220);

    auto tr_0_0_x_zzzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 221);

    auto tr_0_0_x_zzzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 222);

    auto tr_0_0_x_zzzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 223);

    auto tr_0_0_x_zzzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 224);

    #pragma omp simd aligned(tr_0_0_x_zzzz_xxxx, tr_0_0_x_zzzz_xxxy, tr_0_0_x_zzzz_xxxz, tr_0_0_x_zzzz_xxyy, tr_0_0_x_zzzz_xxyz, tr_0_0_x_zzzz_xxzz, tr_0_0_x_zzzz_xyyy, tr_0_0_x_zzzz_xyyz, tr_0_0_x_zzzz_xyzz, tr_0_0_x_zzzz_xzzz, tr_0_0_x_zzzz_yyyy, tr_0_0_x_zzzz_yyyz, tr_0_0_x_zzzz_yyzz, tr_0_0_x_zzzz_yzzz, tr_0_0_x_zzzz_zzzz, tr_xzzzz_xxxx, tr_xzzzz_xxxy, tr_xzzzz_xxxz, tr_xzzzz_xxyy, tr_xzzzz_xxyz, tr_xzzzz_xxzz, tr_xzzzz_xyyy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xzzz, tr_xzzzz_yyyy, tr_xzzzz_yyyz, tr_xzzzz_yyzz, tr_xzzzz_yzzz, tr_xzzzz_zzzz, tr_zzzz_xxx, tr_zzzz_xxxxx, tr_zzzz_xxxxy, tr_zzzz_xxxxz, tr_zzzz_xxxyy, tr_zzzz_xxxyz, tr_zzzz_xxxzz, tr_zzzz_xxy, tr_zzzz_xxyyy, tr_zzzz_xxyyz, tr_zzzz_xxyzz, tr_zzzz_xxz, tr_zzzz_xxzzz, tr_zzzz_xyy, tr_zzzz_xyyyy, tr_zzzz_xyyyz, tr_zzzz_xyyzz, tr_zzzz_xyz, tr_zzzz_xyzzz, tr_zzzz_xzz, tr_zzzz_xzzzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zzzz_xxxx[i] = 2.0 * tr_xzzzz_xxxx[i] * tbe_0 + 2.0 * tr_zzzz_xxxxx[i] * tke_0 - 4.0 * tr_zzzz_xxx[i];

        tr_0_0_x_zzzz_xxxy[i] = 2.0 * tr_xzzzz_xxxy[i] * tbe_0 + 2.0 * tr_zzzz_xxxxy[i] * tke_0 - 3.0 * tr_zzzz_xxy[i];

        tr_0_0_x_zzzz_xxxz[i] = 2.0 * tr_xzzzz_xxxz[i] * tbe_0 + 2.0 * tr_zzzz_xxxxz[i] * tke_0 - 3.0 * tr_zzzz_xxz[i];

        tr_0_0_x_zzzz_xxyy[i] = 2.0 * tr_xzzzz_xxyy[i] * tbe_0 + 2.0 * tr_zzzz_xxxyy[i] * tke_0 - 2.0 * tr_zzzz_xyy[i];

        tr_0_0_x_zzzz_xxyz[i] = 2.0 * tr_xzzzz_xxyz[i] * tbe_0 + 2.0 * tr_zzzz_xxxyz[i] * tke_0 - 2.0 * tr_zzzz_xyz[i];

        tr_0_0_x_zzzz_xxzz[i] = 2.0 * tr_xzzzz_xxzz[i] * tbe_0 + 2.0 * tr_zzzz_xxxzz[i] * tke_0 - 2.0 * tr_zzzz_xzz[i];

        tr_0_0_x_zzzz_xyyy[i] = 2.0 * tr_xzzzz_xyyy[i] * tbe_0 + 2.0 * tr_zzzz_xxyyy[i] * tke_0 - tr_zzzz_yyy[i];

        tr_0_0_x_zzzz_xyyz[i] = 2.0 * tr_xzzzz_xyyz[i] * tbe_0 + 2.0 * tr_zzzz_xxyyz[i] * tke_0 - tr_zzzz_yyz[i];

        tr_0_0_x_zzzz_xyzz[i] = 2.0 * tr_xzzzz_xyzz[i] * tbe_0 + 2.0 * tr_zzzz_xxyzz[i] * tke_0 - tr_zzzz_yzz[i];

        tr_0_0_x_zzzz_xzzz[i] = 2.0 * tr_xzzzz_xzzz[i] * tbe_0 + 2.0 * tr_zzzz_xxzzz[i] * tke_0 - tr_zzzz_zzz[i];

        tr_0_0_x_zzzz_yyyy[i] = 2.0 * tr_xzzzz_yyyy[i] * tbe_0 + 2.0 * tr_zzzz_xyyyy[i] * tke_0;

        tr_0_0_x_zzzz_yyyz[i] = 2.0 * tr_xzzzz_yyyz[i] * tbe_0 + 2.0 * tr_zzzz_xyyyz[i] * tke_0;

        tr_0_0_x_zzzz_yyzz[i] = 2.0 * tr_xzzzz_yyzz[i] * tbe_0 + 2.0 * tr_zzzz_xyyzz[i] * tke_0;

        tr_0_0_x_zzzz_yzzz[i] = 2.0 * tr_xzzzz_yzzz[i] * tbe_0 + 2.0 * tr_zzzz_xyzzz[i] * tke_0;

        tr_0_0_x_zzzz_zzzz[i] = 2.0 * tr_xzzzz_zzzz[i] * tbe_0 + 2.0 * tr_zzzz_xzzzz[i] * tke_0;
    }

    // Set up 225-240 components of targeted buffer : GG

    auto tr_0_0_y_xxxx_xxxx = pbuffer.data(idx_op_geom_010_gg + 225);

    auto tr_0_0_y_xxxx_xxxy = pbuffer.data(idx_op_geom_010_gg + 226);

    auto tr_0_0_y_xxxx_xxxz = pbuffer.data(idx_op_geom_010_gg + 227);

    auto tr_0_0_y_xxxx_xxyy = pbuffer.data(idx_op_geom_010_gg + 228);

    auto tr_0_0_y_xxxx_xxyz = pbuffer.data(idx_op_geom_010_gg + 229);

    auto tr_0_0_y_xxxx_xxzz = pbuffer.data(idx_op_geom_010_gg + 230);

    auto tr_0_0_y_xxxx_xyyy = pbuffer.data(idx_op_geom_010_gg + 231);

    auto tr_0_0_y_xxxx_xyyz = pbuffer.data(idx_op_geom_010_gg + 232);

    auto tr_0_0_y_xxxx_xyzz = pbuffer.data(idx_op_geom_010_gg + 233);

    auto tr_0_0_y_xxxx_xzzz = pbuffer.data(idx_op_geom_010_gg + 234);

    auto tr_0_0_y_xxxx_yyyy = pbuffer.data(idx_op_geom_010_gg + 235);

    auto tr_0_0_y_xxxx_yyyz = pbuffer.data(idx_op_geom_010_gg + 236);

    auto tr_0_0_y_xxxx_yyzz = pbuffer.data(idx_op_geom_010_gg + 237);

    auto tr_0_0_y_xxxx_yzzz = pbuffer.data(idx_op_geom_010_gg + 238);

    auto tr_0_0_y_xxxx_zzzz = pbuffer.data(idx_op_geom_010_gg + 239);

    #pragma omp simd aligned(tr_0_0_y_xxxx_xxxx, tr_0_0_y_xxxx_xxxy, tr_0_0_y_xxxx_xxxz, tr_0_0_y_xxxx_xxyy, tr_0_0_y_xxxx_xxyz, tr_0_0_y_xxxx_xxzz, tr_0_0_y_xxxx_xyyy, tr_0_0_y_xxxx_xyyz, tr_0_0_y_xxxx_xyzz, tr_0_0_y_xxxx_xzzz, tr_0_0_y_xxxx_yyyy, tr_0_0_y_xxxx_yyyz, tr_0_0_y_xxxx_yyzz, tr_0_0_y_xxxx_yzzz, tr_0_0_y_xxxx_zzzz, tr_xxxx_xxx, tr_xxxx_xxxxy, tr_xxxx_xxxyy, tr_xxxx_xxxyz, tr_xxxx_xxy, tr_xxxx_xxyyy, tr_xxxx_xxyyz, tr_xxxx_xxyzz, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyyyy, tr_xxxx_xyyyz, tr_xxxx_xyyzz, tr_xxxx_xyz, tr_xxxx_xyzzz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyyyy, tr_xxxx_yyyyz, tr_xxxx_yyyzz, tr_xxxx_yyz, tr_xxxx_yyzzz, tr_xxxx_yzz, tr_xxxx_yzzzz, tr_xxxx_zzz, tr_xxxxy_xxxx, tr_xxxxy_xxxy, tr_xxxxy_xxxz, tr_xxxxy_xxyy, tr_xxxxy_xxyz, tr_xxxxy_xxzz, tr_xxxxy_xyyy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xzzz, tr_xxxxy_yyyy, tr_xxxxy_yyyz, tr_xxxxy_yyzz, tr_xxxxy_yzzz, tr_xxxxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxx_xxxx[i] = 2.0 * tr_xxxxy_xxxx[i] * tbe_0 + 2.0 * tr_xxxx_xxxxy[i] * tke_0;

        tr_0_0_y_xxxx_xxxy[i] = 2.0 * tr_xxxxy_xxxy[i] * tbe_0 + 2.0 * tr_xxxx_xxxyy[i] * tke_0 - tr_xxxx_xxx[i];

        tr_0_0_y_xxxx_xxxz[i] = 2.0 * tr_xxxxy_xxxz[i] * tbe_0 + 2.0 * tr_xxxx_xxxyz[i] * tke_0;

        tr_0_0_y_xxxx_xxyy[i] = 2.0 * tr_xxxxy_xxyy[i] * tbe_0 + 2.0 * tr_xxxx_xxyyy[i] * tke_0 - 2.0 * tr_xxxx_xxy[i];

        tr_0_0_y_xxxx_xxyz[i] = 2.0 * tr_xxxxy_xxyz[i] * tbe_0 + 2.0 * tr_xxxx_xxyyz[i] * tke_0 - tr_xxxx_xxz[i];

        tr_0_0_y_xxxx_xxzz[i] = 2.0 * tr_xxxxy_xxzz[i] * tbe_0 + 2.0 * tr_xxxx_xxyzz[i] * tke_0;

        tr_0_0_y_xxxx_xyyy[i] = 2.0 * tr_xxxxy_xyyy[i] * tbe_0 + 2.0 * tr_xxxx_xyyyy[i] * tke_0 - 3.0 * tr_xxxx_xyy[i];

        tr_0_0_y_xxxx_xyyz[i] = 2.0 * tr_xxxxy_xyyz[i] * tbe_0 + 2.0 * tr_xxxx_xyyyz[i] * tke_0 - 2.0 * tr_xxxx_xyz[i];

        tr_0_0_y_xxxx_xyzz[i] = 2.0 * tr_xxxxy_xyzz[i] * tbe_0 + 2.0 * tr_xxxx_xyyzz[i] * tke_0 - tr_xxxx_xzz[i];

        tr_0_0_y_xxxx_xzzz[i] = 2.0 * tr_xxxxy_xzzz[i] * tbe_0 + 2.0 * tr_xxxx_xyzzz[i] * tke_0;

        tr_0_0_y_xxxx_yyyy[i] = 2.0 * tr_xxxxy_yyyy[i] * tbe_0 + 2.0 * tr_xxxx_yyyyy[i] * tke_0 - 4.0 * tr_xxxx_yyy[i];

        tr_0_0_y_xxxx_yyyz[i] = 2.0 * tr_xxxxy_yyyz[i] * tbe_0 + 2.0 * tr_xxxx_yyyyz[i] * tke_0 - 3.0 * tr_xxxx_yyz[i];

        tr_0_0_y_xxxx_yyzz[i] = 2.0 * tr_xxxxy_yyzz[i] * tbe_0 + 2.0 * tr_xxxx_yyyzz[i] * tke_0 - 2.0 * tr_xxxx_yzz[i];

        tr_0_0_y_xxxx_yzzz[i] = 2.0 * tr_xxxxy_yzzz[i] * tbe_0 + 2.0 * tr_xxxx_yyzzz[i] * tke_0 - tr_xxxx_zzz[i];

        tr_0_0_y_xxxx_zzzz[i] = 2.0 * tr_xxxxy_zzzz[i] * tbe_0 + 2.0 * tr_xxxx_yzzzz[i] * tke_0;
    }

    // Set up 240-255 components of targeted buffer : GG

    auto tr_0_0_y_xxxy_xxxx = pbuffer.data(idx_op_geom_010_gg + 240);

    auto tr_0_0_y_xxxy_xxxy = pbuffer.data(idx_op_geom_010_gg + 241);

    auto tr_0_0_y_xxxy_xxxz = pbuffer.data(idx_op_geom_010_gg + 242);

    auto tr_0_0_y_xxxy_xxyy = pbuffer.data(idx_op_geom_010_gg + 243);

    auto tr_0_0_y_xxxy_xxyz = pbuffer.data(idx_op_geom_010_gg + 244);

    auto tr_0_0_y_xxxy_xxzz = pbuffer.data(idx_op_geom_010_gg + 245);

    auto tr_0_0_y_xxxy_xyyy = pbuffer.data(idx_op_geom_010_gg + 246);

    auto tr_0_0_y_xxxy_xyyz = pbuffer.data(idx_op_geom_010_gg + 247);

    auto tr_0_0_y_xxxy_xyzz = pbuffer.data(idx_op_geom_010_gg + 248);

    auto tr_0_0_y_xxxy_xzzz = pbuffer.data(idx_op_geom_010_gg + 249);

    auto tr_0_0_y_xxxy_yyyy = pbuffer.data(idx_op_geom_010_gg + 250);

    auto tr_0_0_y_xxxy_yyyz = pbuffer.data(idx_op_geom_010_gg + 251);

    auto tr_0_0_y_xxxy_yyzz = pbuffer.data(idx_op_geom_010_gg + 252);

    auto tr_0_0_y_xxxy_yzzz = pbuffer.data(idx_op_geom_010_gg + 253);

    auto tr_0_0_y_xxxy_zzzz = pbuffer.data(idx_op_geom_010_gg + 254);

    #pragma omp simd aligned(tr_0_0_y_xxxy_xxxx, tr_0_0_y_xxxy_xxxy, tr_0_0_y_xxxy_xxxz, tr_0_0_y_xxxy_xxyy, tr_0_0_y_xxxy_xxyz, tr_0_0_y_xxxy_xxzz, tr_0_0_y_xxxy_xyyy, tr_0_0_y_xxxy_xyyz, tr_0_0_y_xxxy_xyzz, tr_0_0_y_xxxy_xzzz, tr_0_0_y_xxxy_yyyy, tr_0_0_y_xxxy_yyyz, tr_0_0_y_xxxy_yyzz, tr_0_0_y_xxxy_yzzz, tr_0_0_y_xxxy_zzzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxy_xxx, tr_xxxy_xxxxy, tr_xxxy_xxxyy, tr_xxxy_xxxyz, tr_xxxy_xxy, tr_xxxy_xxyyy, tr_xxxy_xxyyz, tr_xxxy_xxyzz, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyyyy, tr_xxxy_xyyyz, tr_xxxy_xyyzz, tr_xxxy_xyz, tr_xxxy_xyzzz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyyyy, tr_xxxy_yyyyz, tr_xxxy_yyyzz, tr_xxxy_yyz, tr_xxxy_yyzzz, tr_xxxy_yzz, tr_xxxy_yzzzz, tr_xxxy_zzz, tr_xxxyy_xxxx, tr_xxxyy_xxxy, tr_xxxyy_xxxz, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xzzz, tr_xxxyy_yyyy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yzzz, tr_xxxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxy_xxxx[i] = 2.0 * tr_xxxyy_xxxx[i] * tbe_0 + 2.0 * tr_xxxy_xxxxy[i] * tke_0 - tr_xxx_xxxx[i];

        tr_0_0_y_xxxy_xxxy[i] = 2.0 * tr_xxxyy_xxxy[i] * tbe_0 + 2.0 * tr_xxxy_xxxyy[i] * tke_0 - tr_xxx_xxxy[i] - tr_xxxy_xxx[i];

        tr_0_0_y_xxxy_xxxz[i] = 2.0 * tr_xxxyy_xxxz[i] * tbe_0 + 2.0 * tr_xxxy_xxxyz[i] * tke_0 - tr_xxx_xxxz[i];

        tr_0_0_y_xxxy_xxyy[i] = 2.0 * tr_xxxyy_xxyy[i] * tbe_0 + 2.0 * tr_xxxy_xxyyy[i] * tke_0 - tr_xxx_xxyy[i] - 2.0 * tr_xxxy_xxy[i];

        tr_0_0_y_xxxy_xxyz[i] = 2.0 * tr_xxxyy_xxyz[i] * tbe_0 + 2.0 * tr_xxxy_xxyyz[i] * tke_0 - tr_xxx_xxyz[i] - tr_xxxy_xxz[i];

        tr_0_0_y_xxxy_xxzz[i] = 2.0 * tr_xxxyy_xxzz[i] * tbe_0 + 2.0 * tr_xxxy_xxyzz[i] * tke_0 - tr_xxx_xxzz[i];

        tr_0_0_y_xxxy_xyyy[i] = 2.0 * tr_xxxyy_xyyy[i] * tbe_0 + 2.0 * tr_xxxy_xyyyy[i] * tke_0 - tr_xxx_xyyy[i] - 3.0 * tr_xxxy_xyy[i];

        tr_0_0_y_xxxy_xyyz[i] = 2.0 * tr_xxxyy_xyyz[i] * tbe_0 + 2.0 * tr_xxxy_xyyyz[i] * tke_0 - tr_xxx_xyyz[i] - 2.0 * tr_xxxy_xyz[i];

        tr_0_0_y_xxxy_xyzz[i] = 2.0 * tr_xxxyy_xyzz[i] * tbe_0 + 2.0 * tr_xxxy_xyyzz[i] * tke_0 - tr_xxx_xyzz[i] - tr_xxxy_xzz[i];

        tr_0_0_y_xxxy_xzzz[i] = 2.0 * tr_xxxyy_xzzz[i] * tbe_0 + 2.0 * tr_xxxy_xyzzz[i] * tke_0 - tr_xxx_xzzz[i];

        tr_0_0_y_xxxy_yyyy[i] = 2.0 * tr_xxxyy_yyyy[i] * tbe_0 + 2.0 * tr_xxxy_yyyyy[i] * tke_0 - tr_xxx_yyyy[i] - 4.0 * tr_xxxy_yyy[i];

        tr_0_0_y_xxxy_yyyz[i] = 2.0 * tr_xxxyy_yyyz[i] * tbe_0 + 2.0 * tr_xxxy_yyyyz[i] * tke_0 - tr_xxx_yyyz[i] - 3.0 * tr_xxxy_yyz[i];

        tr_0_0_y_xxxy_yyzz[i] = 2.0 * tr_xxxyy_yyzz[i] * tbe_0 + 2.0 * tr_xxxy_yyyzz[i] * tke_0 - tr_xxx_yyzz[i] - 2.0 * tr_xxxy_yzz[i];

        tr_0_0_y_xxxy_yzzz[i] = 2.0 * tr_xxxyy_yzzz[i] * tbe_0 + 2.0 * tr_xxxy_yyzzz[i] * tke_0 - tr_xxx_yzzz[i] - tr_xxxy_zzz[i];

        tr_0_0_y_xxxy_zzzz[i] = 2.0 * tr_xxxyy_zzzz[i] * tbe_0 + 2.0 * tr_xxxy_yzzzz[i] * tke_0 - tr_xxx_zzzz[i];
    }

    // Set up 255-270 components of targeted buffer : GG

    auto tr_0_0_y_xxxz_xxxx = pbuffer.data(idx_op_geom_010_gg + 255);

    auto tr_0_0_y_xxxz_xxxy = pbuffer.data(idx_op_geom_010_gg + 256);

    auto tr_0_0_y_xxxz_xxxz = pbuffer.data(idx_op_geom_010_gg + 257);

    auto tr_0_0_y_xxxz_xxyy = pbuffer.data(idx_op_geom_010_gg + 258);

    auto tr_0_0_y_xxxz_xxyz = pbuffer.data(idx_op_geom_010_gg + 259);

    auto tr_0_0_y_xxxz_xxzz = pbuffer.data(idx_op_geom_010_gg + 260);

    auto tr_0_0_y_xxxz_xyyy = pbuffer.data(idx_op_geom_010_gg + 261);

    auto tr_0_0_y_xxxz_xyyz = pbuffer.data(idx_op_geom_010_gg + 262);

    auto tr_0_0_y_xxxz_xyzz = pbuffer.data(idx_op_geom_010_gg + 263);

    auto tr_0_0_y_xxxz_xzzz = pbuffer.data(idx_op_geom_010_gg + 264);

    auto tr_0_0_y_xxxz_yyyy = pbuffer.data(idx_op_geom_010_gg + 265);

    auto tr_0_0_y_xxxz_yyyz = pbuffer.data(idx_op_geom_010_gg + 266);

    auto tr_0_0_y_xxxz_yyzz = pbuffer.data(idx_op_geom_010_gg + 267);

    auto tr_0_0_y_xxxz_yzzz = pbuffer.data(idx_op_geom_010_gg + 268);

    auto tr_0_0_y_xxxz_zzzz = pbuffer.data(idx_op_geom_010_gg + 269);

    #pragma omp simd aligned(tr_0_0_y_xxxz_xxxx, tr_0_0_y_xxxz_xxxy, tr_0_0_y_xxxz_xxxz, tr_0_0_y_xxxz_xxyy, tr_0_0_y_xxxz_xxyz, tr_0_0_y_xxxz_xxzz, tr_0_0_y_xxxz_xyyy, tr_0_0_y_xxxz_xyyz, tr_0_0_y_xxxz_xyzz, tr_0_0_y_xxxz_xzzz, tr_0_0_y_xxxz_yyyy, tr_0_0_y_xxxz_yyyz, tr_0_0_y_xxxz_yyzz, tr_0_0_y_xxxz_yzzz, tr_0_0_y_xxxz_zzzz, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, tr_xxxz_xxx, tr_xxxz_xxxxy, tr_xxxz_xxxyy, tr_xxxz_xxxyz, tr_xxxz_xxy, tr_xxxz_xxyyy, tr_xxxz_xxyyz, tr_xxxz_xxyzz, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyyyy, tr_xxxz_xyyyz, tr_xxxz_xyyzz, tr_xxxz_xyz, tr_xxxz_xyzzz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyyyy, tr_xxxz_yyyyz, tr_xxxz_yyyzz, tr_xxxz_yyz, tr_xxxz_yyzzz, tr_xxxz_yzz, tr_xxxz_yzzzz, tr_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxz_xxxx[i] = 2.0 * tr_xxxyz_xxxx[i] * tbe_0 + 2.0 * tr_xxxz_xxxxy[i] * tke_0;

        tr_0_0_y_xxxz_xxxy[i] = 2.0 * tr_xxxyz_xxxy[i] * tbe_0 + 2.0 * tr_xxxz_xxxyy[i] * tke_0 - tr_xxxz_xxx[i];

        tr_0_0_y_xxxz_xxxz[i] = 2.0 * tr_xxxyz_xxxz[i] * tbe_0 + 2.0 * tr_xxxz_xxxyz[i] * tke_0;

        tr_0_0_y_xxxz_xxyy[i] = 2.0 * tr_xxxyz_xxyy[i] * tbe_0 + 2.0 * tr_xxxz_xxyyy[i] * tke_0 - 2.0 * tr_xxxz_xxy[i];

        tr_0_0_y_xxxz_xxyz[i] = 2.0 * tr_xxxyz_xxyz[i] * tbe_0 + 2.0 * tr_xxxz_xxyyz[i] * tke_0 - tr_xxxz_xxz[i];

        tr_0_0_y_xxxz_xxzz[i] = 2.0 * tr_xxxyz_xxzz[i] * tbe_0 + 2.0 * tr_xxxz_xxyzz[i] * tke_0;

        tr_0_0_y_xxxz_xyyy[i] = 2.0 * tr_xxxyz_xyyy[i] * tbe_0 + 2.0 * tr_xxxz_xyyyy[i] * tke_0 - 3.0 * tr_xxxz_xyy[i];

        tr_0_0_y_xxxz_xyyz[i] = 2.0 * tr_xxxyz_xyyz[i] * tbe_0 + 2.0 * tr_xxxz_xyyyz[i] * tke_0 - 2.0 * tr_xxxz_xyz[i];

        tr_0_0_y_xxxz_xyzz[i] = 2.0 * tr_xxxyz_xyzz[i] * tbe_0 + 2.0 * tr_xxxz_xyyzz[i] * tke_0 - tr_xxxz_xzz[i];

        tr_0_0_y_xxxz_xzzz[i] = 2.0 * tr_xxxyz_xzzz[i] * tbe_0 + 2.0 * tr_xxxz_xyzzz[i] * tke_0;

        tr_0_0_y_xxxz_yyyy[i] = 2.0 * tr_xxxyz_yyyy[i] * tbe_0 + 2.0 * tr_xxxz_yyyyy[i] * tke_0 - 4.0 * tr_xxxz_yyy[i];

        tr_0_0_y_xxxz_yyyz[i] = 2.0 * tr_xxxyz_yyyz[i] * tbe_0 + 2.0 * tr_xxxz_yyyyz[i] * tke_0 - 3.0 * tr_xxxz_yyz[i];

        tr_0_0_y_xxxz_yyzz[i] = 2.0 * tr_xxxyz_yyzz[i] * tbe_0 + 2.0 * tr_xxxz_yyyzz[i] * tke_0 - 2.0 * tr_xxxz_yzz[i];

        tr_0_0_y_xxxz_yzzz[i] = 2.0 * tr_xxxyz_yzzz[i] * tbe_0 + 2.0 * tr_xxxz_yyzzz[i] * tke_0 - tr_xxxz_zzz[i];

        tr_0_0_y_xxxz_zzzz[i] = 2.0 * tr_xxxyz_zzzz[i] * tbe_0 + 2.0 * tr_xxxz_yzzzz[i] * tke_0;
    }

    // Set up 270-285 components of targeted buffer : GG

    auto tr_0_0_y_xxyy_xxxx = pbuffer.data(idx_op_geom_010_gg + 270);

    auto tr_0_0_y_xxyy_xxxy = pbuffer.data(idx_op_geom_010_gg + 271);

    auto tr_0_0_y_xxyy_xxxz = pbuffer.data(idx_op_geom_010_gg + 272);

    auto tr_0_0_y_xxyy_xxyy = pbuffer.data(idx_op_geom_010_gg + 273);

    auto tr_0_0_y_xxyy_xxyz = pbuffer.data(idx_op_geom_010_gg + 274);

    auto tr_0_0_y_xxyy_xxzz = pbuffer.data(idx_op_geom_010_gg + 275);

    auto tr_0_0_y_xxyy_xyyy = pbuffer.data(idx_op_geom_010_gg + 276);

    auto tr_0_0_y_xxyy_xyyz = pbuffer.data(idx_op_geom_010_gg + 277);

    auto tr_0_0_y_xxyy_xyzz = pbuffer.data(idx_op_geom_010_gg + 278);

    auto tr_0_0_y_xxyy_xzzz = pbuffer.data(idx_op_geom_010_gg + 279);

    auto tr_0_0_y_xxyy_yyyy = pbuffer.data(idx_op_geom_010_gg + 280);

    auto tr_0_0_y_xxyy_yyyz = pbuffer.data(idx_op_geom_010_gg + 281);

    auto tr_0_0_y_xxyy_yyzz = pbuffer.data(idx_op_geom_010_gg + 282);

    auto tr_0_0_y_xxyy_yzzz = pbuffer.data(idx_op_geom_010_gg + 283);

    auto tr_0_0_y_xxyy_zzzz = pbuffer.data(idx_op_geom_010_gg + 284);

    #pragma omp simd aligned(tr_0_0_y_xxyy_xxxx, tr_0_0_y_xxyy_xxxy, tr_0_0_y_xxyy_xxxz, tr_0_0_y_xxyy_xxyy, tr_0_0_y_xxyy_xxyz, tr_0_0_y_xxyy_xxzz, tr_0_0_y_xxyy_xyyy, tr_0_0_y_xxyy_xyyz, tr_0_0_y_xxyy_xyzz, tr_0_0_y_xxyy_xzzz, tr_0_0_y_xxyy_yyyy, tr_0_0_y_xxyy_yyyz, tr_0_0_y_xxyy_yyzz, tr_0_0_y_xxyy_yzzz, tr_0_0_y_xxyy_zzzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyy_xxx, tr_xxyy_xxxxy, tr_xxyy_xxxyy, tr_xxyy_xxxyz, tr_xxyy_xxy, tr_xxyy_xxyyy, tr_xxyy_xxyyz, tr_xxyy_xxyzz, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyyyy, tr_xxyy_xyyyz, tr_xxyy_xyyzz, tr_xxyy_xyz, tr_xxyy_xyzzz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyyyy, tr_xxyy_yyyyz, tr_xxyy_yyyzz, tr_xxyy_yyz, tr_xxyy_yyzzz, tr_xxyy_yzz, tr_xxyy_yzzzz, tr_xxyy_zzz, tr_xxyyy_xxxx, tr_xxyyy_xxxy, tr_xxyyy_xxxz, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xzzz, tr_xxyyy_yyyy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yzzz, tr_xxyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxyy_xxxx[i] = 2.0 * tr_xxyyy_xxxx[i] * tbe_0 + 2.0 * tr_xxyy_xxxxy[i] * tke_0 - 2.0 * tr_xxy_xxxx[i];

        tr_0_0_y_xxyy_xxxy[i] = 2.0 * tr_xxyyy_xxxy[i] * tbe_0 + 2.0 * tr_xxyy_xxxyy[i] * tke_0 - 2.0 * tr_xxy_xxxy[i] - tr_xxyy_xxx[i];

        tr_0_0_y_xxyy_xxxz[i] = 2.0 * tr_xxyyy_xxxz[i] * tbe_0 + 2.0 * tr_xxyy_xxxyz[i] * tke_0 - 2.0 * tr_xxy_xxxz[i];

        tr_0_0_y_xxyy_xxyy[i] = 2.0 * tr_xxyyy_xxyy[i] * tbe_0 + 2.0 * tr_xxyy_xxyyy[i] * tke_0 - 2.0 * tr_xxy_xxyy[i] - 2.0 * tr_xxyy_xxy[i];

        tr_0_0_y_xxyy_xxyz[i] = 2.0 * tr_xxyyy_xxyz[i] * tbe_0 + 2.0 * tr_xxyy_xxyyz[i] * tke_0 - 2.0 * tr_xxy_xxyz[i] - tr_xxyy_xxz[i];

        tr_0_0_y_xxyy_xxzz[i] = 2.0 * tr_xxyyy_xxzz[i] * tbe_0 + 2.0 * tr_xxyy_xxyzz[i] * tke_0 - 2.0 * tr_xxy_xxzz[i];

        tr_0_0_y_xxyy_xyyy[i] = 2.0 * tr_xxyyy_xyyy[i] * tbe_0 + 2.0 * tr_xxyy_xyyyy[i] * tke_0 - 2.0 * tr_xxy_xyyy[i] - 3.0 * tr_xxyy_xyy[i];

        tr_0_0_y_xxyy_xyyz[i] = 2.0 * tr_xxyyy_xyyz[i] * tbe_0 + 2.0 * tr_xxyy_xyyyz[i] * tke_0 - 2.0 * tr_xxy_xyyz[i] - 2.0 * tr_xxyy_xyz[i];

        tr_0_0_y_xxyy_xyzz[i] = 2.0 * tr_xxyyy_xyzz[i] * tbe_0 + 2.0 * tr_xxyy_xyyzz[i] * tke_0 - 2.0 * tr_xxy_xyzz[i] - tr_xxyy_xzz[i];

        tr_0_0_y_xxyy_xzzz[i] = 2.0 * tr_xxyyy_xzzz[i] * tbe_0 + 2.0 * tr_xxyy_xyzzz[i] * tke_0 - 2.0 * tr_xxy_xzzz[i];

        tr_0_0_y_xxyy_yyyy[i] = 2.0 * tr_xxyyy_yyyy[i] * tbe_0 + 2.0 * tr_xxyy_yyyyy[i] * tke_0 - 2.0 * tr_xxy_yyyy[i] - 4.0 * tr_xxyy_yyy[i];

        tr_0_0_y_xxyy_yyyz[i] = 2.0 * tr_xxyyy_yyyz[i] * tbe_0 + 2.0 * tr_xxyy_yyyyz[i] * tke_0 - 2.0 * tr_xxy_yyyz[i] - 3.0 * tr_xxyy_yyz[i];

        tr_0_0_y_xxyy_yyzz[i] = 2.0 * tr_xxyyy_yyzz[i] * tbe_0 + 2.0 * tr_xxyy_yyyzz[i] * tke_0 - 2.0 * tr_xxy_yyzz[i] - 2.0 * tr_xxyy_yzz[i];

        tr_0_0_y_xxyy_yzzz[i] = 2.0 * tr_xxyyy_yzzz[i] * tbe_0 + 2.0 * tr_xxyy_yyzzz[i] * tke_0 - 2.0 * tr_xxy_yzzz[i] - tr_xxyy_zzz[i];

        tr_0_0_y_xxyy_zzzz[i] = 2.0 * tr_xxyyy_zzzz[i] * tbe_0 + 2.0 * tr_xxyy_yzzzz[i] * tke_0 - 2.0 * tr_xxy_zzzz[i];
    }

    // Set up 285-300 components of targeted buffer : GG

    auto tr_0_0_y_xxyz_xxxx = pbuffer.data(idx_op_geom_010_gg + 285);

    auto tr_0_0_y_xxyz_xxxy = pbuffer.data(idx_op_geom_010_gg + 286);

    auto tr_0_0_y_xxyz_xxxz = pbuffer.data(idx_op_geom_010_gg + 287);

    auto tr_0_0_y_xxyz_xxyy = pbuffer.data(idx_op_geom_010_gg + 288);

    auto tr_0_0_y_xxyz_xxyz = pbuffer.data(idx_op_geom_010_gg + 289);

    auto tr_0_0_y_xxyz_xxzz = pbuffer.data(idx_op_geom_010_gg + 290);

    auto tr_0_0_y_xxyz_xyyy = pbuffer.data(idx_op_geom_010_gg + 291);

    auto tr_0_0_y_xxyz_xyyz = pbuffer.data(idx_op_geom_010_gg + 292);

    auto tr_0_0_y_xxyz_xyzz = pbuffer.data(idx_op_geom_010_gg + 293);

    auto tr_0_0_y_xxyz_xzzz = pbuffer.data(idx_op_geom_010_gg + 294);

    auto tr_0_0_y_xxyz_yyyy = pbuffer.data(idx_op_geom_010_gg + 295);

    auto tr_0_0_y_xxyz_yyyz = pbuffer.data(idx_op_geom_010_gg + 296);

    auto tr_0_0_y_xxyz_yyzz = pbuffer.data(idx_op_geom_010_gg + 297);

    auto tr_0_0_y_xxyz_yzzz = pbuffer.data(idx_op_geom_010_gg + 298);

    auto tr_0_0_y_xxyz_zzzz = pbuffer.data(idx_op_geom_010_gg + 299);

    #pragma omp simd aligned(tr_0_0_y_xxyz_xxxx, tr_0_0_y_xxyz_xxxy, tr_0_0_y_xxyz_xxxz, tr_0_0_y_xxyz_xxyy, tr_0_0_y_xxyz_xxyz, tr_0_0_y_xxyz_xxzz, tr_0_0_y_xxyz_xyyy, tr_0_0_y_xxyz_xyyz, tr_0_0_y_xxyz_xyzz, tr_0_0_y_xxyz_xzzz, tr_0_0_y_xxyz_yyyy, tr_0_0_y_xxyz_yyyz, tr_0_0_y_xxyz_yyzz, tr_0_0_y_xxyz_yzzz, tr_0_0_y_xxyz_zzzz, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxy, tr_xxyz_xxxyy, tr_xxyz_xxxyz, tr_xxyz_xxy, tr_xxyz_xxyyy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyyyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyyyy, tr_xxyz_yyyyz, tr_xxyz_yyyzz, tr_xxyz_yyz, tr_xxyz_yyzzz, tr_xxyz_yzz, tr_xxyz_yzzzz, tr_xxyz_zzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxyz_xxxx[i] = 2.0 * tr_xxyyz_xxxx[i] * tbe_0 + 2.0 * tr_xxyz_xxxxy[i] * tke_0 - tr_xxz_xxxx[i];

        tr_0_0_y_xxyz_xxxy[i] = 2.0 * tr_xxyyz_xxxy[i] * tbe_0 + 2.0 * tr_xxyz_xxxyy[i] * tke_0 - tr_xxz_xxxy[i] - tr_xxyz_xxx[i];

        tr_0_0_y_xxyz_xxxz[i] = 2.0 * tr_xxyyz_xxxz[i] * tbe_0 + 2.0 * tr_xxyz_xxxyz[i] * tke_0 - tr_xxz_xxxz[i];

        tr_0_0_y_xxyz_xxyy[i] = 2.0 * tr_xxyyz_xxyy[i] * tbe_0 + 2.0 * tr_xxyz_xxyyy[i] * tke_0 - tr_xxz_xxyy[i] - 2.0 * tr_xxyz_xxy[i];

        tr_0_0_y_xxyz_xxyz[i] = 2.0 * tr_xxyyz_xxyz[i] * tbe_0 + 2.0 * tr_xxyz_xxyyz[i] * tke_0 - tr_xxz_xxyz[i] - tr_xxyz_xxz[i];

        tr_0_0_y_xxyz_xxzz[i] = 2.0 * tr_xxyyz_xxzz[i] * tbe_0 + 2.0 * tr_xxyz_xxyzz[i] * tke_0 - tr_xxz_xxzz[i];

        tr_0_0_y_xxyz_xyyy[i] = 2.0 * tr_xxyyz_xyyy[i] * tbe_0 + 2.0 * tr_xxyz_xyyyy[i] * tke_0 - tr_xxz_xyyy[i] - 3.0 * tr_xxyz_xyy[i];

        tr_0_0_y_xxyz_xyyz[i] = 2.0 * tr_xxyyz_xyyz[i] * tbe_0 + 2.0 * tr_xxyz_xyyyz[i] * tke_0 - tr_xxz_xyyz[i] - 2.0 * tr_xxyz_xyz[i];

        tr_0_0_y_xxyz_xyzz[i] = 2.0 * tr_xxyyz_xyzz[i] * tbe_0 + 2.0 * tr_xxyz_xyyzz[i] * tke_0 - tr_xxz_xyzz[i] - tr_xxyz_xzz[i];

        tr_0_0_y_xxyz_xzzz[i] = 2.0 * tr_xxyyz_xzzz[i] * tbe_0 + 2.0 * tr_xxyz_xyzzz[i] * tke_0 - tr_xxz_xzzz[i];

        tr_0_0_y_xxyz_yyyy[i] = 2.0 * tr_xxyyz_yyyy[i] * tbe_0 + 2.0 * tr_xxyz_yyyyy[i] * tke_0 - tr_xxz_yyyy[i] - 4.0 * tr_xxyz_yyy[i];

        tr_0_0_y_xxyz_yyyz[i] = 2.0 * tr_xxyyz_yyyz[i] * tbe_0 + 2.0 * tr_xxyz_yyyyz[i] * tke_0 - tr_xxz_yyyz[i] - 3.0 * tr_xxyz_yyz[i];

        tr_0_0_y_xxyz_yyzz[i] = 2.0 * tr_xxyyz_yyzz[i] * tbe_0 + 2.0 * tr_xxyz_yyyzz[i] * tke_0 - tr_xxz_yyzz[i] - 2.0 * tr_xxyz_yzz[i];

        tr_0_0_y_xxyz_yzzz[i] = 2.0 * tr_xxyyz_yzzz[i] * tbe_0 + 2.0 * tr_xxyz_yyzzz[i] * tke_0 - tr_xxz_yzzz[i] - tr_xxyz_zzz[i];

        tr_0_0_y_xxyz_zzzz[i] = 2.0 * tr_xxyyz_zzzz[i] * tbe_0 + 2.0 * tr_xxyz_yzzzz[i] * tke_0 - tr_xxz_zzzz[i];
    }

    // Set up 300-315 components of targeted buffer : GG

    auto tr_0_0_y_xxzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 300);

    auto tr_0_0_y_xxzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 301);

    auto tr_0_0_y_xxzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 302);

    auto tr_0_0_y_xxzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 303);

    auto tr_0_0_y_xxzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 304);

    auto tr_0_0_y_xxzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 305);

    auto tr_0_0_y_xxzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 306);

    auto tr_0_0_y_xxzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 307);

    auto tr_0_0_y_xxzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 308);

    auto tr_0_0_y_xxzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 309);

    auto tr_0_0_y_xxzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 310);

    auto tr_0_0_y_xxzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 311);

    auto tr_0_0_y_xxzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 312);

    auto tr_0_0_y_xxzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 313);

    auto tr_0_0_y_xxzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 314);

    #pragma omp simd aligned(tr_0_0_y_xxzz_xxxx, tr_0_0_y_xxzz_xxxy, tr_0_0_y_xxzz_xxxz, tr_0_0_y_xxzz_xxyy, tr_0_0_y_xxzz_xxyz, tr_0_0_y_xxzz_xxzz, tr_0_0_y_xxzz_xyyy, tr_0_0_y_xxzz_xyyz, tr_0_0_y_xxzz_xyzz, tr_0_0_y_xxzz_xzzz, tr_0_0_y_xxzz_yyyy, tr_0_0_y_xxzz_yyyz, tr_0_0_y_xxzz_yyzz, tr_0_0_y_xxzz_yzzz, tr_0_0_y_xxzz_zzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_xxzz_xxx, tr_xxzz_xxxxy, tr_xxzz_xxxyy, tr_xxzz_xxxyz, tr_xxzz_xxy, tr_xxzz_xxyyy, tr_xxzz_xxyyz, tr_xxzz_xxyzz, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyyyy, tr_xxzz_xyyyz, tr_xxzz_xyyzz, tr_xxzz_xyz, tr_xxzz_xyzzz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyyyy, tr_xxzz_yyyyz, tr_xxzz_yyyzz, tr_xxzz_yyz, tr_xxzz_yyzzz, tr_xxzz_yzz, tr_xxzz_yzzzz, tr_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxzz_xxxx[i] = 2.0 * tr_xxyzz_xxxx[i] * tbe_0 + 2.0 * tr_xxzz_xxxxy[i] * tke_0;

        tr_0_0_y_xxzz_xxxy[i] = 2.0 * tr_xxyzz_xxxy[i] * tbe_0 + 2.0 * tr_xxzz_xxxyy[i] * tke_0 - tr_xxzz_xxx[i];

        tr_0_0_y_xxzz_xxxz[i] = 2.0 * tr_xxyzz_xxxz[i] * tbe_0 + 2.0 * tr_xxzz_xxxyz[i] * tke_0;

        tr_0_0_y_xxzz_xxyy[i] = 2.0 * tr_xxyzz_xxyy[i] * tbe_0 + 2.0 * tr_xxzz_xxyyy[i] * tke_0 - 2.0 * tr_xxzz_xxy[i];

        tr_0_0_y_xxzz_xxyz[i] = 2.0 * tr_xxyzz_xxyz[i] * tbe_0 + 2.0 * tr_xxzz_xxyyz[i] * tke_0 - tr_xxzz_xxz[i];

        tr_0_0_y_xxzz_xxzz[i] = 2.0 * tr_xxyzz_xxzz[i] * tbe_0 + 2.0 * tr_xxzz_xxyzz[i] * tke_0;

        tr_0_0_y_xxzz_xyyy[i] = 2.0 * tr_xxyzz_xyyy[i] * tbe_0 + 2.0 * tr_xxzz_xyyyy[i] * tke_0 - 3.0 * tr_xxzz_xyy[i];

        tr_0_0_y_xxzz_xyyz[i] = 2.0 * tr_xxyzz_xyyz[i] * tbe_0 + 2.0 * tr_xxzz_xyyyz[i] * tke_0 - 2.0 * tr_xxzz_xyz[i];

        tr_0_0_y_xxzz_xyzz[i] = 2.0 * tr_xxyzz_xyzz[i] * tbe_0 + 2.0 * tr_xxzz_xyyzz[i] * tke_0 - tr_xxzz_xzz[i];

        tr_0_0_y_xxzz_xzzz[i] = 2.0 * tr_xxyzz_xzzz[i] * tbe_0 + 2.0 * tr_xxzz_xyzzz[i] * tke_0;

        tr_0_0_y_xxzz_yyyy[i] = 2.0 * tr_xxyzz_yyyy[i] * tbe_0 + 2.0 * tr_xxzz_yyyyy[i] * tke_0 - 4.0 * tr_xxzz_yyy[i];

        tr_0_0_y_xxzz_yyyz[i] = 2.0 * tr_xxyzz_yyyz[i] * tbe_0 + 2.0 * tr_xxzz_yyyyz[i] * tke_0 - 3.0 * tr_xxzz_yyz[i];

        tr_0_0_y_xxzz_yyzz[i] = 2.0 * tr_xxyzz_yyzz[i] * tbe_0 + 2.0 * tr_xxzz_yyyzz[i] * tke_0 - 2.0 * tr_xxzz_yzz[i];

        tr_0_0_y_xxzz_yzzz[i] = 2.0 * tr_xxyzz_yzzz[i] * tbe_0 + 2.0 * tr_xxzz_yyzzz[i] * tke_0 - tr_xxzz_zzz[i];

        tr_0_0_y_xxzz_zzzz[i] = 2.0 * tr_xxyzz_zzzz[i] * tbe_0 + 2.0 * tr_xxzz_yzzzz[i] * tke_0;
    }

    // Set up 315-330 components of targeted buffer : GG

    auto tr_0_0_y_xyyy_xxxx = pbuffer.data(idx_op_geom_010_gg + 315);

    auto tr_0_0_y_xyyy_xxxy = pbuffer.data(idx_op_geom_010_gg + 316);

    auto tr_0_0_y_xyyy_xxxz = pbuffer.data(idx_op_geom_010_gg + 317);

    auto tr_0_0_y_xyyy_xxyy = pbuffer.data(idx_op_geom_010_gg + 318);

    auto tr_0_0_y_xyyy_xxyz = pbuffer.data(idx_op_geom_010_gg + 319);

    auto tr_0_0_y_xyyy_xxzz = pbuffer.data(idx_op_geom_010_gg + 320);

    auto tr_0_0_y_xyyy_xyyy = pbuffer.data(idx_op_geom_010_gg + 321);

    auto tr_0_0_y_xyyy_xyyz = pbuffer.data(idx_op_geom_010_gg + 322);

    auto tr_0_0_y_xyyy_xyzz = pbuffer.data(idx_op_geom_010_gg + 323);

    auto tr_0_0_y_xyyy_xzzz = pbuffer.data(idx_op_geom_010_gg + 324);

    auto tr_0_0_y_xyyy_yyyy = pbuffer.data(idx_op_geom_010_gg + 325);

    auto tr_0_0_y_xyyy_yyyz = pbuffer.data(idx_op_geom_010_gg + 326);

    auto tr_0_0_y_xyyy_yyzz = pbuffer.data(idx_op_geom_010_gg + 327);

    auto tr_0_0_y_xyyy_yzzz = pbuffer.data(idx_op_geom_010_gg + 328);

    auto tr_0_0_y_xyyy_zzzz = pbuffer.data(idx_op_geom_010_gg + 329);

    #pragma omp simd aligned(tr_0_0_y_xyyy_xxxx, tr_0_0_y_xyyy_xxxy, tr_0_0_y_xyyy_xxxz, tr_0_0_y_xyyy_xxyy, tr_0_0_y_xyyy_xxyz, tr_0_0_y_xyyy_xxzz, tr_0_0_y_xyyy_xyyy, tr_0_0_y_xyyy_xyyz, tr_0_0_y_xyyy_xyzz, tr_0_0_y_xyyy_xzzz, tr_0_0_y_xyyy_yyyy, tr_0_0_y_xyyy_yyyz, tr_0_0_y_xyyy_yyzz, tr_0_0_y_xyyy_yzzz, tr_0_0_y_xyyy_zzzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyy_xxx, tr_xyyy_xxxxy, tr_xyyy_xxxyy, tr_xyyy_xxxyz, tr_xyyy_xxy, tr_xyyy_xxyyy, tr_xyyy_xxyyz, tr_xyyy_xxyzz, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyyyy, tr_xyyy_xyyyz, tr_xyyy_xyyzz, tr_xyyy_xyz, tr_xyyy_xyzzz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyyyy, tr_xyyy_yyyyz, tr_xyyy_yyyzz, tr_xyyy_yyz, tr_xyyy_yyzzz, tr_xyyy_yzz, tr_xyyy_yzzzz, tr_xyyy_zzz, tr_xyyyy_xxxx, tr_xyyyy_xxxy, tr_xyyyy_xxxz, tr_xyyyy_xxyy, tr_xyyyy_xxyz, tr_xyyyy_xxzz, tr_xyyyy_xyyy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xzzz, tr_xyyyy_yyyy, tr_xyyyy_yyyz, tr_xyyyy_yyzz, tr_xyyyy_yzzz, tr_xyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyyy_xxxx[i] = 2.0 * tr_xyyyy_xxxx[i] * tbe_0 + 2.0 * tr_xyyy_xxxxy[i] * tke_0 - 3.0 * tr_xyy_xxxx[i];

        tr_0_0_y_xyyy_xxxy[i] = 2.0 * tr_xyyyy_xxxy[i] * tbe_0 + 2.0 * tr_xyyy_xxxyy[i] * tke_0 - 3.0 * tr_xyy_xxxy[i] - tr_xyyy_xxx[i];

        tr_0_0_y_xyyy_xxxz[i] = 2.0 * tr_xyyyy_xxxz[i] * tbe_0 + 2.0 * tr_xyyy_xxxyz[i] * tke_0 - 3.0 * tr_xyy_xxxz[i];

        tr_0_0_y_xyyy_xxyy[i] = 2.0 * tr_xyyyy_xxyy[i] * tbe_0 + 2.0 * tr_xyyy_xxyyy[i] * tke_0 - 3.0 * tr_xyy_xxyy[i] - 2.0 * tr_xyyy_xxy[i];

        tr_0_0_y_xyyy_xxyz[i] = 2.0 * tr_xyyyy_xxyz[i] * tbe_0 + 2.0 * tr_xyyy_xxyyz[i] * tke_0 - 3.0 * tr_xyy_xxyz[i] - tr_xyyy_xxz[i];

        tr_0_0_y_xyyy_xxzz[i] = 2.0 * tr_xyyyy_xxzz[i] * tbe_0 + 2.0 * tr_xyyy_xxyzz[i] * tke_0 - 3.0 * tr_xyy_xxzz[i];

        tr_0_0_y_xyyy_xyyy[i] = 2.0 * tr_xyyyy_xyyy[i] * tbe_0 + 2.0 * tr_xyyy_xyyyy[i] * tke_0 - 3.0 * tr_xyy_xyyy[i] - 3.0 * tr_xyyy_xyy[i];

        tr_0_0_y_xyyy_xyyz[i] = 2.0 * tr_xyyyy_xyyz[i] * tbe_0 + 2.0 * tr_xyyy_xyyyz[i] * tke_0 - 3.0 * tr_xyy_xyyz[i] - 2.0 * tr_xyyy_xyz[i];

        tr_0_0_y_xyyy_xyzz[i] = 2.0 * tr_xyyyy_xyzz[i] * tbe_0 + 2.0 * tr_xyyy_xyyzz[i] * tke_0 - 3.0 * tr_xyy_xyzz[i] - tr_xyyy_xzz[i];

        tr_0_0_y_xyyy_xzzz[i] = 2.0 * tr_xyyyy_xzzz[i] * tbe_0 + 2.0 * tr_xyyy_xyzzz[i] * tke_0 - 3.0 * tr_xyy_xzzz[i];

        tr_0_0_y_xyyy_yyyy[i] = 2.0 * tr_xyyyy_yyyy[i] * tbe_0 + 2.0 * tr_xyyy_yyyyy[i] * tke_0 - 3.0 * tr_xyy_yyyy[i] - 4.0 * tr_xyyy_yyy[i];

        tr_0_0_y_xyyy_yyyz[i] = 2.0 * tr_xyyyy_yyyz[i] * tbe_0 + 2.0 * tr_xyyy_yyyyz[i] * tke_0 - 3.0 * tr_xyy_yyyz[i] - 3.0 * tr_xyyy_yyz[i];

        tr_0_0_y_xyyy_yyzz[i] = 2.0 * tr_xyyyy_yyzz[i] * tbe_0 + 2.0 * tr_xyyy_yyyzz[i] * tke_0 - 3.0 * tr_xyy_yyzz[i] - 2.0 * tr_xyyy_yzz[i];

        tr_0_0_y_xyyy_yzzz[i] = 2.0 * tr_xyyyy_yzzz[i] * tbe_0 + 2.0 * tr_xyyy_yyzzz[i] * tke_0 - 3.0 * tr_xyy_yzzz[i] - tr_xyyy_zzz[i];

        tr_0_0_y_xyyy_zzzz[i] = 2.0 * tr_xyyyy_zzzz[i] * tbe_0 + 2.0 * tr_xyyy_yzzzz[i] * tke_0 - 3.0 * tr_xyy_zzzz[i];
    }

    // Set up 330-345 components of targeted buffer : GG

    auto tr_0_0_y_xyyz_xxxx = pbuffer.data(idx_op_geom_010_gg + 330);

    auto tr_0_0_y_xyyz_xxxy = pbuffer.data(idx_op_geom_010_gg + 331);

    auto tr_0_0_y_xyyz_xxxz = pbuffer.data(idx_op_geom_010_gg + 332);

    auto tr_0_0_y_xyyz_xxyy = pbuffer.data(idx_op_geom_010_gg + 333);

    auto tr_0_0_y_xyyz_xxyz = pbuffer.data(idx_op_geom_010_gg + 334);

    auto tr_0_0_y_xyyz_xxzz = pbuffer.data(idx_op_geom_010_gg + 335);

    auto tr_0_0_y_xyyz_xyyy = pbuffer.data(idx_op_geom_010_gg + 336);

    auto tr_0_0_y_xyyz_xyyz = pbuffer.data(idx_op_geom_010_gg + 337);

    auto tr_0_0_y_xyyz_xyzz = pbuffer.data(idx_op_geom_010_gg + 338);

    auto tr_0_0_y_xyyz_xzzz = pbuffer.data(idx_op_geom_010_gg + 339);

    auto tr_0_0_y_xyyz_yyyy = pbuffer.data(idx_op_geom_010_gg + 340);

    auto tr_0_0_y_xyyz_yyyz = pbuffer.data(idx_op_geom_010_gg + 341);

    auto tr_0_0_y_xyyz_yyzz = pbuffer.data(idx_op_geom_010_gg + 342);

    auto tr_0_0_y_xyyz_yzzz = pbuffer.data(idx_op_geom_010_gg + 343);

    auto tr_0_0_y_xyyz_zzzz = pbuffer.data(idx_op_geom_010_gg + 344);

    #pragma omp simd aligned(tr_0_0_y_xyyz_xxxx, tr_0_0_y_xyyz_xxxy, tr_0_0_y_xyyz_xxxz, tr_0_0_y_xyyz_xxyy, tr_0_0_y_xyyz_xxyz, tr_0_0_y_xyyz_xxzz, tr_0_0_y_xyyz_xyyy, tr_0_0_y_xyyz_xyyz, tr_0_0_y_xyyz_xyzz, tr_0_0_y_xyyz_xzzz, tr_0_0_y_xyyz_yyyy, tr_0_0_y_xyyz_yyyz, tr_0_0_y_xyyz_yyzz, tr_0_0_y_xyyz_yzzz, tr_0_0_y_xyyz_zzzz, tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxy, tr_xyyz_xxxyy, tr_xyyz_xxxyz, tr_xyyz_xxy, tr_xyyz_xxyyy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyyyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyyyy, tr_xyyz_yyyyz, tr_xyyz_yyyzz, tr_xyyz_yyz, tr_xyyz_yyzzz, tr_xyyz_yzz, tr_xyyz_yzzzz, tr_xyyz_zzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyyz_xxxx[i] = 2.0 * tr_xyyyz_xxxx[i] * tbe_0 + 2.0 * tr_xyyz_xxxxy[i] * tke_0 - 2.0 * tr_xyz_xxxx[i];

        tr_0_0_y_xyyz_xxxy[i] = 2.0 * tr_xyyyz_xxxy[i] * tbe_0 + 2.0 * tr_xyyz_xxxyy[i] * tke_0 - 2.0 * tr_xyz_xxxy[i] - tr_xyyz_xxx[i];

        tr_0_0_y_xyyz_xxxz[i] = 2.0 * tr_xyyyz_xxxz[i] * tbe_0 + 2.0 * tr_xyyz_xxxyz[i] * tke_0 - 2.0 * tr_xyz_xxxz[i];

        tr_0_0_y_xyyz_xxyy[i] = 2.0 * tr_xyyyz_xxyy[i] * tbe_0 + 2.0 * tr_xyyz_xxyyy[i] * tke_0 - 2.0 * tr_xyz_xxyy[i] - 2.0 * tr_xyyz_xxy[i];

        tr_0_0_y_xyyz_xxyz[i] = 2.0 * tr_xyyyz_xxyz[i] * tbe_0 + 2.0 * tr_xyyz_xxyyz[i] * tke_0 - 2.0 * tr_xyz_xxyz[i] - tr_xyyz_xxz[i];

        tr_0_0_y_xyyz_xxzz[i] = 2.0 * tr_xyyyz_xxzz[i] * tbe_0 + 2.0 * tr_xyyz_xxyzz[i] * tke_0 - 2.0 * tr_xyz_xxzz[i];

        tr_0_0_y_xyyz_xyyy[i] = 2.0 * tr_xyyyz_xyyy[i] * tbe_0 + 2.0 * tr_xyyz_xyyyy[i] * tke_0 - 2.0 * tr_xyz_xyyy[i] - 3.0 * tr_xyyz_xyy[i];

        tr_0_0_y_xyyz_xyyz[i] = 2.0 * tr_xyyyz_xyyz[i] * tbe_0 + 2.0 * tr_xyyz_xyyyz[i] * tke_0 - 2.0 * tr_xyz_xyyz[i] - 2.0 * tr_xyyz_xyz[i];

        tr_0_0_y_xyyz_xyzz[i] = 2.0 * tr_xyyyz_xyzz[i] * tbe_0 + 2.0 * tr_xyyz_xyyzz[i] * tke_0 - 2.0 * tr_xyz_xyzz[i] - tr_xyyz_xzz[i];

        tr_0_0_y_xyyz_xzzz[i] = 2.0 * tr_xyyyz_xzzz[i] * tbe_0 + 2.0 * tr_xyyz_xyzzz[i] * tke_0 - 2.0 * tr_xyz_xzzz[i];

        tr_0_0_y_xyyz_yyyy[i] = 2.0 * tr_xyyyz_yyyy[i] * tbe_0 + 2.0 * tr_xyyz_yyyyy[i] * tke_0 - 2.0 * tr_xyz_yyyy[i] - 4.0 * tr_xyyz_yyy[i];

        tr_0_0_y_xyyz_yyyz[i] = 2.0 * tr_xyyyz_yyyz[i] * tbe_0 + 2.0 * tr_xyyz_yyyyz[i] * tke_0 - 2.0 * tr_xyz_yyyz[i] - 3.0 * tr_xyyz_yyz[i];

        tr_0_0_y_xyyz_yyzz[i] = 2.0 * tr_xyyyz_yyzz[i] * tbe_0 + 2.0 * tr_xyyz_yyyzz[i] * tke_0 - 2.0 * tr_xyz_yyzz[i] - 2.0 * tr_xyyz_yzz[i];

        tr_0_0_y_xyyz_yzzz[i] = 2.0 * tr_xyyyz_yzzz[i] * tbe_0 + 2.0 * tr_xyyz_yyzzz[i] * tke_0 - 2.0 * tr_xyz_yzzz[i] - tr_xyyz_zzz[i];

        tr_0_0_y_xyyz_zzzz[i] = 2.0 * tr_xyyyz_zzzz[i] * tbe_0 + 2.0 * tr_xyyz_yzzzz[i] * tke_0 - 2.0 * tr_xyz_zzzz[i];
    }

    // Set up 345-360 components of targeted buffer : GG

    auto tr_0_0_y_xyzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 345);

    auto tr_0_0_y_xyzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 346);

    auto tr_0_0_y_xyzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 347);

    auto tr_0_0_y_xyzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 348);

    auto tr_0_0_y_xyzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 349);

    auto tr_0_0_y_xyzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 350);

    auto tr_0_0_y_xyzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 351);

    auto tr_0_0_y_xyzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 352);

    auto tr_0_0_y_xyzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 353);

    auto tr_0_0_y_xyzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 354);

    auto tr_0_0_y_xyzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 355);

    auto tr_0_0_y_xyzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 356);

    auto tr_0_0_y_xyzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 357);

    auto tr_0_0_y_xyzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 358);

    auto tr_0_0_y_xyzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 359);

    #pragma omp simd aligned(tr_0_0_y_xyzz_xxxx, tr_0_0_y_xyzz_xxxy, tr_0_0_y_xyzz_xxxz, tr_0_0_y_xyzz_xxyy, tr_0_0_y_xyzz_xxyz, tr_0_0_y_xyzz_xxzz, tr_0_0_y_xyzz_xyyy, tr_0_0_y_xyzz_xyyz, tr_0_0_y_xyzz_xyzz, tr_0_0_y_xyzz_xzzz, tr_0_0_y_xyzz_yyyy, tr_0_0_y_xyzz_yyyz, tr_0_0_y_xyzz_yyzz, tr_0_0_y_xyzz_yzzz, tr_0_0_y_xyzz_zzzz, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxy, tr_xyzz_xxxyy, tr_xyzz_xxxyz, tr_xyzz_xxy, tr_xyzz_xxyyy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyyyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyyyy, tr_xyzz_yyyyz, tr_xyzz_yyyzz, tr_xyzz_yyz, tr_xyzz_yyzzz, tr_xyzz_yzz, tr_xyzz_yzzzz, tr_xyzz_zzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyzz_xxxx[i] = 2.0 * tr_xyyzz_xxxx[i] * tbe_0 + 2.0 * tr_xyzz_xxxxy[i] * tke_0 - tr_xzz_xxxx[i];

        tr_0_0_y_xyzz_xxxy[i] = 2.0 * tr_xyyzz_xxxy[i] * tbe_0 + 2.0 * tr_xyzz_xxxyy[i] * tke_0 - tr_xzz_xxxy[i] - tr_xyzz_xxx[i];

        tr_0_0_y_xyzz_xxxz[i] = 2.0 * tr_xyyzz_xxxz[i] * tbe_0 + 2.0 * tr_xyzz_xxxyz[i] * tke_0 - tr_xzz_xxxz[i];

        tr_0_0_y_xyzz_xxyy[i] = 2.0 * tr_xyyzz_xxyy[i] * tbe_0 + 2.0 * tr_xyzz_xxyyy[i] * tke_0 - tr_xzz_xxyy[i] - 2.0 * tr_xyzz_xxy[i];

        tr_0_0_y_xyzz_xxyz[i] = 2.0 * tr_xyyzz_xxyz[i] * tbe_0 + 2.0 * tr_xyzz_xxyyz[i] * tke_0 - tr_xzz_xxyz[i] - tr_xyzz_xxz[i];

        tr_0_0_y_xyzz_xxzz[i] = 2.0 * tr_xyyzz_xxzz[i] * tbe_0 + 2.0 * tr_xyzz_xxyzz[i] * tke_0 - tr_xzz_xxzz[i];

        tr_0_0_y_xyzz_xyyy[i] = 2.0 * tr_xyyzz_xyyy[i] * tbe_0 + 2.0 * tr_xyzz_xyyyy[i] * tke_0 - tr_xzz_xyyy[i] - 3.0 * tr_xyzz_xyy[i];

        tr_0_0_y_xyzz_xyyz[i] = 2.0 * tr_xyyzz_xyyz[i] * tbe_0 + 2.0 * tr_xyzz_xyyyz[i] * tke_0 - tr_xzz_xyyz[i] - 2.0 * tr_xyzz_xyz[i];

        tr_0_0_y_xyzz_xyzz[i] = 2.0 * tr_xyyzz_xyzz[i] * tbe_0 + 2.0 * tr_xyzz_xyyzz[i] * tke_0 - tr_xzz_xyzz[i] - tr_xyzz_xzz[i];

        tr_0_0_y_xyzz_xzzz[i] = 2.0 * tr_xyyzz_xzzz[i] * tbe_0 + 2.0 * tr_xyzz_xyzzz[i] * tke_0 - tr_xzz_xzzz[i];

        tr_0_0_y_xyzz_yyyy[i] = 2.0 * tr_xyyzz_yyyy[i] * tbe_0 + 2.0 * tr_xyzz_yyyyy[i] * tke_0 - tr_xzz_yyyy[i] - 4.0 * tr_xyzz_yyy[i];

        tr_0_0_y_xyzz_yyyz[i] = 2.0 * tr_xyyzz_yyyz[i] * tbe_0 + 2.0 * tr_xyzz_yyyyz[i] * tke_0 - tr_xzz_yyyz[i] - 3.0 * tr_xyzz_yyz[i];

        tr_0_0_y_xyzz_yyzz[i] = 2.0 * tr_xyyzz_yyzz[i] * tbe_0 + 2.0 * tr_xyzz_yyyzz[i] * tke_0 - tr_xzz_yyzz[i] - 2.0 * tr_xyzz_yzz[i];

        tr_0_0_y_xyzz_yzzz[i] = 2.0 * tr_xyyzz_yzzz[i] * tbe_0 + 2.0 * tr_xyzz_yyzzz[i] * tke_0 - tr_xzz_yzzz[i] - tr_xyzz_zzz[i];

        tr_0_0_y_xyzz_zzzz[i] = 2.0 * tr_xyyzz_zzzz[i] * tbe_0 + 2.0 * tr_xyzz_yzzzz[i] * tke_0 - tr_xzz_zzzz[i];
    }

    // Set up 360-375 components of targeted buffer : GG

    auto tr_0_0_y_xzzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 360);

    auto tr_0_0_y_xzzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 361);

    auto tr_0_0_y_xzzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 362);

    auto tr_0_0_y_xzzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 363);

    auto tr_0_0_y_xzzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 364);

    auto tr_0_0_y_xzzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 365);

    auto tr_0_0_y_xzzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 366);

    auto tr_0_0_y_xzzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 367);

    auto tr_0_0_y_xzzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 368);

    auto tr_0_0_y_xzzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 369);

    auto tr_0_0_y_xzzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 370);

    auto tr_0_0_y_xzzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 371);

    auto tr_0_0_y_xzzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 372);

    auto tr_0_0_y_xzzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 373);

    auto tr_0_0_y_xzzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 374);

    #pragma omp simd aligned(tr_0_0_y_xzzz_xxxx, tr_0_0_y_xzzz_xxxy, tr_0_0_y_xzzz_xxxz, tr_0_0_y_xzzz_xxyy, tr_0_0_y_xzzz_xxyz, tr_0_0_y_xzzz_xxzz, tr_0_0_y_xzzz_xyyy, tr_0_0_y_xzzz_xyyz, tr_0_0_y_xzzz_xyzz, tr_0_0_y_xzzz_xzzz, tr_0_0_y_xzzz_yyyy, tr_0_0_y_xzzz_yyyz, tr_0_0_y_xzzz_yyzz, tr_0_0_y_xzzz_yzzz, tr_0_0_y_xzzz_zzzz, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxxxy, tr_xzzz_xxxyy, tr_xzzz_xxxyz, tr_xzzz_xxy, tr_xzzz_xxyyy, tr_xzzz_xxyyz, tr_xzzz_xxyzz, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyyyy, tr_xzzz_xyyyz, tr_xzzz_xyyzz, tr_xzzz_xyz, tr_xzzz_xyzzz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyyyy, tr_xzzz_yyyyz, tr_xzzz_yyyzz, tr_xzzz_yyz, tr_xzzz_yyzzz, tr_xzzz_yzz, tr_xzzz_yzzzz, tr_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xzzz_xxxx[i] = 2.0 * tr_xyzzz_xxxx[i] * tbe_0 + 2.0 * tr_xzzz_xxxxy[i] * tke_0;

        tr_0_0_y_xzzz_xxxy[i] = 2.0 * tr_xyzzz_xxxy[i] * tbe_0 + 2.0 * tr_xzzz_xxxyy[i] * tke_0 - tr_xzzz_xxx[i];

        tr_0_0_y_xzzz_xxxz[i] = 2.0 * tr_xyzzz_xxxz[i] * tbe_0 + 2.0 * tr_xzzz_xxxyz[i] * tke_0;

        tr_0_0_y_xzzz_xxyy[i] = 2.0 * tr_xyzzz_xxyy[i] * tbe_0 + 2.0 * tr_xzzz_xxyyy[i] * tke_0 - 2.0 * tr_xzzz_xxy[i];

        tr_0_0_y_xzzz_xxyz[i] = 2.0 * tr_xyzzz_xxyz[i] * tbe_0 + 2.0 * tr_xzzz_xxyyz[i] * tke_0 - tr_xzzz_xxz[i];

        tr_0_0_y_xzzz_xxzz[i] = 2.0 * tr_xyzzz_xxzz[i] * tbe_0 + 2.0 * tr_xzzz_xxyzz[i] * tke_0;

        tr_0_0_y_xzzz_xyyy[i] = 2.0 * tr_xyzzz_xyyy[i] * tbe_0 + 2.0 * tr_xzzz_xyyyy[i] * tke_0 - 3.0 * tr_xzzz_xyy[i];

        tr_0_0_y_xzzz_xyyz[i] = 2.0 * tr_xyzzz_xyyz[i] * tbe_0 + 2.0 * tr_xzzz_xyyyz[i] * tke_0 - 2.0 * tr_xzzz_xyz[i];

        tr_0_0_y_xzzz_xyzz[i] = 2.0 * tr_xyzzz_xyzz[i] * tbe_0 + 2.0 * tr_xzzz_xyyzz[i] * tke_0 - tr_xzzz_xzz[i];

        tr_0_0_y_xzzz_xzzz[i] = 2.0 * tr_xyzzz_xzzz[i] * tbe_0 + 2.0 * tr_xzzz_xyzzz[i] * tke_0;

        tr_0_0_y_xzzz_yyyy[i] = 2.0 * tr_xyzzz_yyyy[i] * tbe_0 + 2.0 * tr_xzzz_yyyyy[i] * tke_0 - 4.0 * tr_xzzz_yyy[i];

        tr_0_0_y_xzzz_yyyz[i] = 2.0 * tr_xyzzz_yyyz[i] * tbe_0 + 2.0 * tr_xzzz_yyyyz[i] * tke_0 - 3.0 * tr_xzzz_yyz[i];

        tr_0_0_y_xzzz_yyzz[i] = 2.0 * tr_xyzzz_yyzz[i] * tbe_0 + 2.0 * tr_xzzz_yyyzz[i] * tke_0 - 2.0 * tr_xzzz_yzz[i];

        tr_0_0_y_xzzz_yzzz[i] = 2.0 * tr_xyzzz_yzzz[i] * tbe_0 + 2.0 * tr_xzzz_yyzzz[i] * tke_0 - tr_xzzz_zzz[i];

        tr_0_0_y_xzzz_zzzz[i] = 2.0 * tr_xyzzz_zzzz[i] * tbe_0 + 2.0 * tr_xzzz_yzzzz[i] * tke_0;
    }

    // Set up 375-390 components of targeted buffer : GG

    auto tr_0_0_y_yyyy_xxxx = pbuffer.data(idx_op_geom_010_gg + 375);

    auto tr_0_0_y_yyyy_xxxy = pbuffer.data(idx_op_geom_010_gg + 376);

    auto tr_0_0_y_yyyy_xxxz = pbuffer.data(idx_op_geom_010_gg + 377);

    auto tr_0_0_y_yyyy_xxyy = pbuffer.data(idx_op_geom_010_gg + 378);

    auto tr_0_0_y_yyyy_xxyz = pbuffer.data(idx_op_geom_010_gg + 379);

    auto tr_0_0_y_yyyy_xxzz = pbuffer.data(idx_op_geom_010_gg + 380);

    auto tr_0_0_y_yyyy_xyyy = pbuffer.data(idx_op_geom_010_gg + 381);

    auto tr_0_0_y_yyyy_xyyz = pbuffer.data(idx_op_geom_010_gg + 382);

    auto tr_0_0_y_yyyy_xyzz = pbuffer.data(idx_op_geom_010_gg + 383);

    auto tr_0_0_y_yyyy_xzzz = pbuffer.data(idx_op_geom_010_gg + 384);

    auto tr_0_0_y_yyyy_yyyy = pbuffer.data(idx_op_geom_010_gg + 385);

    auto tr_0_0_y_yyyy_yyyz = pbuffer.data(idx_op_geom_010_gg + 386);

    auto tr_0_0_y_yyyy_yyzz = pbuffer.data(idx_op_geom_010_gg + 387);

    auto tr_0_0_y_yyyy_yzzz = pbuffer.data(idx_op_geom_010_gg + 388);

    auto tr_0_0_y_yyyy_zzzz = pbuffer.data(idx_op_geom_010_gg + 389);

    #pragma omp simd aligned(tr_0_0_y_yyyy_xxxx, tr_0_0_y_yyyy_xxxy, tr_0_0_y_yyyy_xxxz, tr_0_0_y_yyyy_xxyy, tr_0_0_y_yyyy_xxyz, tr_0_0_y_yyyy_xxzz, tr_0_0_y_yyyy_xyyy, tr_0_0_y_yyyy_xyyz, tr_0_0_y_yyyy_xyzz, tr_0_0_y_yyyy_xzzz, tr_0_0_y_yyyy_yyyy, tr_0_0_y_yyyy_yyyz, tr_0_0_y_yyyy_yyzz, tr_0_0_y_yyyy_yzzz, tr_0_0_y_yyyy_zzzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, tr_yyyy_xxx, tr_yyyy_xxxxy, tr_yyyy_xxxyy, tr_yyyy_xxxyz, tr_yyyy_xxy, tr_yyyy_xxyyy, tr_yyyy_xxyyz, tr_yyyy_xxyzz, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyyyy, tr_yyyy_xyyyz, tr_yyyy_xyyzz, tr_yyyy_xyz, tr_yyyy_xyzzz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyyyy, tr_yyyy_yyyyz, tr_yyyy_yyyzz, tr_yyyy_yyz, tr_yyyy_yyzzz, tr_yyyy_yzz, tr_yyyy_yzzzz, tr_yyyy_zzz, tr_yyyyy_xxxx, tr_yyyyy_xxxy, tr_yyyyy_xxxz, tr_yyyyy_xxyy, tr_yyyyy_xxyz, tr_yyyyy_xxzz, tr_yyyyy_xyyy, tr_yyyyy_xyyz, tr_yyyyy_xyzz, tr_yyyyy_xzzz, tr_yyyyy_yyyy, tr_yyyyy_yyyz, tr_yyyyy_yyzz, tr_yyyyy_yzzz, tr_yyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyyy_xxxx[i] = 2.0 * tr_yyyyy_xxxx[i] * tbe_0 + 2.0 * tr_yyyy_xxxxy[i] * tke_0 - 4.0 * tr_yyy_xxxx[i];

        tr_0_0_y_yyyy_xxxy[i] = 2.0 * tr_yyyyy_xxxy[i] * tbe_0 + 2.0 * tr_yyyy_xxxyy[i] * tke_0 - 4.0 * tr_yyy_xxxy[i] - tr_yyyy_xxx[i];

        tr_0_0_y_yyyy_xxxz[i] = 2.0 * tr_yyyyy_xxxz[i] * tbe_0 + 2.0 * tr_yyyy_xxxyz[i] * tke_0 - 4.0 * tr_yyy_xxxz[i];

        tr_0_0_y_yyyy_xxyy[i] = 2.0 * tr_yyyyy_xxyy[i] * tbe_0 + 2.0 * tr_yyyy_xxyyy[i] * tke_0 - 4.0 * tr_yyy_xxyy[i] - 2.0 * tr_yyyy_xxy[i];

        tr_0_0_y_yyyy_xxyz[i] = 2.0 * tr_yyyyy_xxyz[i] * tbe_0 + 2.0 * tr_yyyy_xxyyz[i] * tke_0 - 4.0 * tr_yyy_xxyz[i] - tr_yyyy_xxz[i];

        tr_0_0_y_yyyy_xxzz[i] = 2.0 * tr_yyyyy_xxzz[i] * tbe_0 + 2.0 * tr_yyyy_xxyzz[i] * tke_0 - 4.0 * tr_yyy_xxzz[i];

        tr_0_0_y_yyyy_xyyy[i] = 2.0 * tr_yyyyy_xyyy[i] * tbe_0 + 2.0 * tr_yyyy_xyyyy[i] * tke_0 - 4.0 * tr_yyy_xyyy[i] - 3.0 * tr_yyyy_xyy[i];

        tr_0_0_y_yyyy_xyyz[i] = 2.0 * tr_yyyyy_xyyz[i] * tbe_0 + 2.0 * tr_yyyy_xyyyz[i] * tke_0 - 4.0 * tr_yyy_xyyz[i] - 2.0 * tr_yyyy_xyz[i];

        tr_0_0_y_yyyy_xyzz[i] = 2.0 * tr_yyyyy_xyzz[i] * tbe_0 + 2.0 * tr_yyyy_xyyzz[i] * tke_0 - 4.0 * tr_yyy_xyzz[i] - tr_yyyy_xzz[i];

        tr_0_0_y_yyyy_xzzz[i] = 2.0 * tr_yyyyy_xzzz[i] * tbe_0 + 2.0 * tr_yyyy_xyzzz[i] * tke_0 - 4.0 * tr_yyy_xzzz[i];

        tr_0_0_y_yyyy_yyyy[i] = 2.0 * tr_yyyyy_yyyy[i] * tbe_0 + 2.0 * tr_yyyy_yyyyy[i] * tke_0 - 4.0 * tr_yyy_yyyy[i] - 4.0 * tr_yyyy_yyy[i];

        tr_0_0_y_yyyy_yyyz[i] = 2.0 * tr_yyyyy_yyyz[i] * tbe_0 + 2.0 * tr_yyyy_yyyyz[i] * tke_0 - 4.0 * tr_yyy_yyyz[i] - 3.0 * tr_yyyy_yyz[i];

        tr_0_0_y_yyyy_yyzz[i] = 2.0 * tr_yyyyy_yyzz[i] * tbe_0 + 2.0 * tr_yyyy_yyyzz[i] * tke_0 - 4.0 * tr_yyy_yyzz[i] - 2.0 * tr_yyyy_yzz[i];

        tr_0_0_y_yyyy_yzzz[i] = 2.0 * tr_yyyyy_yzzz[i] * tbe_0 + 2.0 * tr_yyyy_yyzzz[i] * tke_0 - 4.0 * tr_yyy_yzzz[i] - tr_yyyy_zzz[i];

        tr_0_0_y_yyyy_zzzz[i] = 2.0 * tr_yyyyy_zzzz[i] * tbe_0 + 2.0 * tr_yyyy_yzzzz[i] * tke_0 - 4.0 * tr_yyy_zzzz[i];
    }

    // Set up 390-405 components of targeted buffer : GG

    auto tr_0_0_y_yyyz_xxxx = pbuffer.data(idx_op_geom_010_gg + 390);

    auto tr_0_0_y_yyyz_xxxy = pbuffer.data(idx_op_geom_010_gg + 391);

    auto tr_0_0_y_yyyz_xxxz = pbuffer.data(idx_op_geom_010_gg + 392);

    auto tr_0_0_y_yyyz_xxyy = pbuffer.data(idx_op_geom_010_gg + 393);

    auto tr_0_0_y_yyyz_xxyz = pbuffer.data(idx_op_geom_010_gg + 394);

    auto tr_0_0_y_yyyz_xxzz = pbuffer.data(idx_op_geom_010_gg + 395);

    auto tr_0_0_y_yyyz_xyyy = pbuffer.data(idx_op_geom_010_gg + 396);

    auto tr_0_0_y_yyyz_xyyz = pbuffer.data(idx_op_geom_010_gg + 397);

    auto tr_0_0_y_yyyz_xyzz = pbuffer.data(idx_op_geom_010_gg + 398);

    auto tr_0_0_y_yyyz_xzzz = pbuffer.data(idx_op_geom_010_gg + 399);

    auto tr_0_0_y_yyyz_yyyy = pbuffer.data(idx_op_geom_010_gg + 400);

    auto tr_0_0_y_yyyz_yyyz = pbuffer.data(idx_op_geom_010_gg + 401);

    auto tr_0_0_y_yyyz_yyzz = pbuffer.data(idx_op_geom_010_gg + 402);

    auto tr_0_0_y_yyyz_yzzz = pbuffer.data(idx_op_geom_010_gg + 403);

    auto tr_0_0_y_yyyz_zzzz = pbuffer.data(idx_op_geom_010_gg + 404);

    #pragma omp simd aligned(tr_0_0_y_yyyz_xxxx, tr_0_0_y_yyyz_xxxy, tr_0_0_y_yyyz_xxxz, tr_0_0_y_yyyz_xxyy, tr_0_0_y_yyyz_xxyz, tr_0_0_y_yyyz_xxzz, tr_0_0_y_yyyz_xyyy, tr_0_0_y_yyyz_xyyz, tr_0_0_y_yyyz_xyzz, tr_0_0_y_yyyz_xzzz, tr_0_0_y_yyyz_yyyy, tr_0_0_y_yyyz_yyyz, tr_0_0_y_yyyz_yyzz, tr_0_0_y_yyyz_yzzz, tr_0_0_y_yyyz_zzzz, tr_yyyyz_xxxx, tr_yyyyz_xxxy, tr_yyyyz_xxxz, tr_yyyyz_xxyy, tr_yyyyz_xxyz, tr_yyyyz_xxzz, tr_yyyyz_xyyy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xzzz, tr_yyyyz_yyyy, tr_yyyyz_yyyz, tr_yyyyz_yyzz, tr_yyyyz_yzzz, tr_yyyyz_zzzz, tr_yyyz_xxx, tr_yyyz_xxxxy, tr_yyyz_xxxyy, tr_yyyz_xxxyz, tr_yyyz_xxy, tr_yyyz_xxyyy, tr_yyyz_xxyyz, tr_yyyz_xxyzz, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyyyy, tr_yyyz_xyyyz, tr_yyyz_xyyzz, tr_yyyz_xyz, tr_yyyz_xyzzz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyyyy, tr_yyyz_yyyyz, tr_yyyz_yyyzz, tr_yyyz_yyz, tr_yyyz_yyzzz, tr_yyyz_yzz, tr_yyyz_yzzzz, tr_yyyz_zzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyyz_xxxx[i] = 2.0 * tr_yyyyz_xxxx[i] * tbe_0 + 2.0 * tr_yyyz_xxxxy[i] * tke_0 - 3.0 * tr_yyz_xxxx[i];

        tr_0_0_y_yyyz_xxxy[i] = 2.0 * tr_yyyyz_xxxy[i] * tbe_0 + 2.0 * tr_yyyz_xxxyy[i] * tke_0 - 3.0 * tr_yyz_xxxy[i] - tr_yyyz_xxx[i];

        tr_0_0_y_yyyz_xxxz[i] = 2.0 * tr_yyyyz_xxxz[i] * tbe_0 + 2.0 * tr_yyyz_xxxyz[i] * tke_0 - 3.0 * tr_yyz_xxxz[i];

        tr_0_0_y_yyyz_xxyy[i] = 2.0 * tr_yyyyz_xxyy[i] * tbe_0 + 2.0 * tr_yyyz_xxyyy[i] * tke_0 - 3.0 * tr_yyz_xxyy[i] - 2.0 * tr_yyyz_xxy[i];

        tr_0_0_y_yyyz_xxyz[i] = 2.0 * tr_yyyyz_xxyz[i] * tbe_0 + 2.0 * tr_yyyz_xxyyz[i] * tke_0 - 3.0 * tr_yyz_xxyz[i] - tr_yyyz_xxz[i];

        tr_0_0_y_yyyz_xxzz[i] = 2.0 * tr_yyyyz_xxzz[i] * tbe_0 + 2.0 * tr_yyyz_xxyzz[i] * tke_0 - 3.0 * tr_yyz_xxzz[i];

        tr_0_0_y_yyyz_xyyy[i] = 2.0 * tr_yyyyz_xyyy[i] * tbe_0 + 2.0 * tr_yyyz_xyyyy[i] * tke_0 - 3.0 * tr_yyz_xyyy[i] - 3.0 * tr_yyyz_xyy[i];

        tr_0_0_y_yyyz_xyyz[i] = 2.0 * tr_yyyyz_xyyz[i] * tbe_0 + 2.0 * tr_yyyz_xyyyz[i] * tke_0 - 3.0 * tr_yyz_xyyz[i] - 2.0 * tr_yyyz_xyz[i];

        tr_0_0_y_yyyz_xyzz[i] = 2.0 * tr_yyyyz_xyzz[i] * tbe_0 + 2.0 * tr_yyyz_xyyzz[i] * tke_0 - 3.0 * tr_yyz_xyzz[i] - tr_yyyz_xzz[i];

        tr_0_0_y_yyyz_xzzz[i] = 2.0 * tr_yyyyz_xzzz[i] * tbe_0 + 2.0 * tr_yyyz_xyzzz[i] * tke_0 - 3.0 * tr_yyz_xzzz[i];

        tr_0_0_y_yyyz_yyyy[i] = 2.0 * tr_yyyyz_yyyy[i] * tbe_0 + 2.0 * tr_yyyz_yyyyy[i] * tke_0 - 3.0 * tr_yyz_yyyy[i] - 4.0 * tr_yyyz_yyy[i];

        tr_0_0_y_yyyz_yyyz[i] = 2.0 * tr_yyyyz_yyyz[i] * tbe_0 + 2.0 * tr_yyyz_yyyyz[i] * tke_0 - 3.0 * tr_yyz_yyyz[i] - 3.0 * tr_yyyz_yyz[i];

        tr_0_0_y_yyyz_yyzz[i] = 2.0 * tr_yyyyz_yyzz[i] * tbe_0 + 2.0 * tr_yyyz_yyyzz[i] * tke_0 - 3.0 * tr_yyz_yyzz[i] - 2.0 * tr_yyyz_yzz[i];

        tr_0_0_y_yyyz_yzzz[i] = 2.0 * tr_yyyyz_yzzz[i] * tbe_0 + 2.0 * tr_yyyz_yyzzz[i] * tke_0 - 3.0 * tr_yyz_yzzz[i] - tr_yyyz_zzz[i];

        tr_0_0_y_yyyz_zzzz[i] = 2.0 * tr_yyyyz_zzzz[i] * tbe_0 + 2.0 * tr_yyyz_yzzzz[i] * tke_0 - 3.0 * tr_yyz_zzzz[i];
    }

    // Set up 405-420 components of targeted buffer : GG

    auto tr_0_0_y_yyzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 405);

    auto tr_0_0_y_yyzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 406);

    auto tr_0_0_y_yyzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 407);

    auto tr_0_0_y_yyzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 408);

    auto tr_0_0_y_yyzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 409);

    auto tr_0_0_y_yyzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 410);

    auto tr_0_0_y_yyzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 411);

    auto tr_0_0_y_yyzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 412);

    auto tr_0_0_y_yyzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 413);

    auto tr_0_0_y_yyzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 414);

    auto tr_0_0_y_yyzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 415);

    auto tr_0_0_y_yyzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 416);

    auto tr_0_0_y_yyzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 417);

    auto tr_0_0_y_yyzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 418);

    auto tr_0_0_y_yyzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 419);

    #pragma omp simd aligned(tr_0_0_y_yyzz_xxxx, tr_0_0_y_yyzz_xxxy, tr_0_0_y_yyzz_xxxz, tr_0_0_y_yyzz_xxyy, tr_0_0_y_yyzz_xxyz, tr_0_0_y_yyzz_xxzz, tr_0_0_y_yyzz_xyyy, tr_0_0_y_yyzz_xyyz, tr_0_0_y_yyzz_xyzz, tr_0_0_y_yyzz_xzzz, tr_0_0_y_yyzz_yyyy, tr_0_0_y_yyzz_yyyz, tr_0_0_y_yyzz_yyzz, tr_0_0_y_yyzz_yzzz, tr_0_0_y_yyzz_zzzz, tr_yyyzz_xxxx, tr_yyyzz_xxxy, tr_yyyzz_xxxz, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xzzz, tr_yyyzz_yyyy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yzzz, tr_yyyzz_zzzz, tr_yyzz_xxx, tr_yyzz_xxxxy, tr_yyzz_xxxyy, tr_yyzz_xxxyz, tr_yyzz_xxy, tr_yyzz_xxyyy, tr_yyzz_xxyyz, tr_yyzz_xxyzz, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyyyy, tr_yyzz_xyyyz, tr_yyzz_xyyzz, tr_yyzz_xyz, tr_yyzz_xyzzz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyyyy, tr_yyzz_yyyyz, tr_yyzz_yyyzz, tr_yyzz_yyz, tr_yyzz_yyzzz, tr_yyzz_yzz, tr_yyzz_yzzzz, tr_yyzz_zzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyzz_xxxx[i] = 2.0 * tr_yyyzz_xxxx[i] * tbe_0 + 2.0 * tr_yyzz_xxxxy[i] * tke_0 - 2.0 * tr_yzz_xxxx[i];

        tr_0_0_y_yyzz_xxxy[i] = 2.0 * tr_yyyzz_xxxy[i] * tbe_0 + 2.0 * tr_yyzz_xxxyy[i] * tke_0 - 2.0 * tr_yzz_xxxy[i] - tr_yyzz_xxx[i];

        tr_0_0_y_yyzz_xxxz[i] = 2.0 * tr_yyyzz_xxxz[i] * tbe_0 + 2.0 * tr_yyzz_xxxyz[i] * tke_0 - 2.0 * tr_yzz_xxxz[i];

        tr_0_0_y_yyzz_xxyy[i] = 2.0 * tr_yyyzz_xxyy[i] * tbe_0 + 2.0 * tr_yyzz_xxyyy[i] * tke_0 - 2.0 * tr_yzz_xxyy[i] - 2.0 * tr_yyzz_xxy[i];

        tr_0_0_y_yyzz_xxyz[i] = 2.0 * tr_yyyzz_xxyz[i] * tbe_0 + 2.0 * tr_yyzz_xxyyz[i] * tke_0 - 2.0 * tr_yzz_xxyz[i] - tr_yyzz_xxz[i];

        tr_0_0_y_yyzz_xxzz[i] = 2.0 * tr_yyyzz_xxzz[i] * tbe_0 + 2.0 * tr_yyzz_xxyzz[i] * tke_0 - 2.0 * tr_yzz_xxzz[i];

        tr_0_0_y_yyzz_xyyy[i] = 2.0 * tr_yyyzz_xyyy[i] * tbe_0 + 2.0 * tr_yyzz_xyyyy[i] * tke_0 - 2.0 * tr_yzz_xyyy[i] - 3.0 * tr_yyzz_xyy[i];

        tr_0_0_y_yyzz_xyyz[i] = 2.0 * tr_yyyzz_xyyz[i] * tbe_0 + 2.0 * tr_yyzz_xyyyz[i] * tke_0 - 2.0 * tr_yzz_xyyz[i] - 2.0 * tr_yyzz_xyz[i];

        tr_0_0_y_yyzz_xyzz[i] = 2.0 * tr_yyyzz_xyzz[i] * tbe_0 + 2.0 * tr_yyzz_xyyzz[i] * tke_0 - 2.0 * tr_yzz_xyzz[i] - tr_yyzz_xzz[i];

        tr_0_0_y_yyzz_xzzz[i] = 2.0 * tr_yyyzz_xzzz[i] * tbe_0 + 2.0 * tr_yyzz_xyzzz[i] * tke_0 - 2.0 * tr_yzz_xzzz[i];

        tr_0_0_y_yyzz_yyyy[i] = 2.0 * tr_yyyzz_yyyy[i] * tbe_0 + 2.0 * tr_yyzz_yyyyy[i] * tke_0 - 2.0 * tr_yzz_yyyy[i] - 4.0 * tr_yyzz_yyy[i];

        tr_0_0_y_yyzz_yyyz[i] = 2.0 * tr_yyyzz_yyyz[i] * tbe_0 + 2.0 * tr_yyzz_yyyyz[i] * tke_0 - 2.0 * tr_yzz_yyyz[i] - 3.0 * tr_yyzz_yyz[i];

        tr_0_0_y_yyzz_yyzz[i] = 2.0 * tr_yyyzz_yyzz[i] * tbe_0 + 2.0 * tr_yyzz_yyyzz[i] * tke_0 - 2.0 * tr_yzz_yyzz[i] - 2.0 * tr_yyzz_yzz[i];

        tr_0_0_y_yyzz_yzzz[i] = 2.0 * tr_yyyzz_yzzz[i] * tbe_0 + 2.0 * tr_yyzz_yyzzz[i] * tke_0 - 2.0 * tr_yzz_yzzz[i] - tr_yyzz_zzz[i];

        tr_0_0_y_yyzz_zzzz[i] = 2.0 * tr_yyyzz_zzzz[i] * tbe_0 + 2.0 * tr_yyzz_yzzzz[i] * tke_0 - 2.0 * tr_yzz_zzzz[i];
    }

    // Set up 420-435 components of targeted buffer : GG

    auto tr_0_0_y_yzzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 420);

    auto tr_0_0_y_yzzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 421);

    auto tr_0_0_y_yzzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 422);

    auto tr_0_0_y_yzzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 423);

    auto tr_0_0_y_yzzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 424);

    auto tr_0_0_y_yzzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 425);

    auto tr_0_0_y_yzzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 426);

    auto tr_0_0_y_yzzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 427);

    auto tr_0_0_y_yzzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 428);

    auto tr_0_0_y_yzzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 429);

    auto tr_0_0_y_yzzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 430);

    auto tr_0_0_y_yzzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 431);

    auto tr_0_0_y_yzzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 432);

    auto tr_0_0_y_yzzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 433);

    auto tr_0_0_y_yzzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 434);

    #pragma omp simd aligned(tr_0_0_y_yzzz_xxxx, tr_0_0_y_yzzz_xxxy, tr_0_0_y_yzzz_xxxz, tr_0_0_y_yzzz_xxyy, tr_0_0_y_yzzz_xxyz, tr_0_0_y_yzzz_xxzz, tr_0_0_y_yzzz_xyyy, tr_0_0_y_yzzz_xyyz, tr_0_0_y_yzzz_xyzz, tr_0_0_y_yzzz_xzzz, tr_0_0_y_yzzz_yyyy, tr_0_0_y_yzzz_yyyz, tr_0_0_y_yzzz_yyzz, tr_0_0_y_yzzz_yzzz, tr_0_0_y_yzzz_zzzz, tr_yyzzz_xxxx, tr_yyzzz_xxxy, tr_yyzzz_xxxz, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xzzz, tr_yyzzz_yyyy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yzzz, tr_yyzzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxxxy, tr_yzzz_xxxyy, tr_yzzz_xxxyz, tr_yzzz_xxy, tr_yzzz_xxyyy, tr_yzzz_xxyyz, tr_yzzz_xxyzz, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyyyy, tr_yzzz_xyyyz, tr_yzzz_xyyzz, tr_yzzz_xyz, tr_yzzz_xyzzz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyyyy, tr_yzzz_yyyyz, tr_yzzz_yyyzz, tr_yzzz_yyz, tr_yzzz_yyzzz, tr_yzzz_yzz, tr_yzzz_yzzzz, tr_yzzz_zzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yzzz_xxxx[i] = 2.0 * tr_yyzzz_xxxx[i] * tbe_0 + 2.0 * tr_yzzz_xxxxy[i] * tke_0 - tr_zzz_xxxx[i];

        tr_0_0_y_yzzz_xxxy[i] = 2.0 * tr_yyzzz_xxxy[i] * tbe_0 + 2.0 * tr_yzzz_xxxyy[i] * tke_0 - tr_zzz_xxxy[i] - tr_yzzz_xxx[i];

        tr_0_0_y_yzzz_xxxz[i] = 2.0 * tr_yyzzz_xxxz[i] * tbe_0 + 2.0 * tr_yzzz_xxxyz[i] * tke_0 - tr_zzz_xxxz[i];

        tr_0_0_y_yzzz_xxyy[i] = 2.0 * tr_yyzzz_xxyy[i] * tbe_0 + 2.0 * tr_yzzz_xxyyy[i] * tke_0 - tr_zzz_xxyy[i] - 2.0 * tr_yzzz_xxy[i];

        tr_0_0_y_yzzz_xxyz[i] = 2.0 * tr_yyzzz_xxyz[i] * tbe_0 + 2.0 * tr_yzzz_xxyyz[i] * tke_0 - tr_zzz_xxyz[i] - tr_yzzz_xxz[i];

        tr_0_0_y_yzzz_xxzz[i] = 2.0 * tr_yyzzz_xxzz[i] * tbe_0 + 2.0 * tr_yzzz_xxyzz[i] * tke_0 - tr_zzz_xxzz[i];

        tr_0_0_y_yzzz_xyyy[i] = 2.0 * tr_yyzzz_xyyy[i] * tbe_0 + 2.0 * tr_yzzz_xyyyy[i] * tke_0 - tr_zzz_xyyy[i] - 3.0 * tr_yzzz_xyy[i];

        tr_0_0_y_yzzz_xyyz[i] = 2.0 * tr_yyzzz_xyyz[i] * tbe_0 + 2.0 * tr_yzzz_xyyyz[i] * tke_0 - tr_zzz_xyyz[i] - 2.0 * tr_yzzz_xyz[i];

        tr_0_0_y_yzzz_xyzz[i] = 2.0 * tr_yyzzz_xyzz[i] * tbe_0 + 2.0 * tr_yzzz_xyyzz[i] * tke_0 - tr_zzz_xyzz[i] - tr_yzzz_xzz[i];

        tr_0_0_y_yzzz_xzzz[i] = 2.0 * tr_yyzzz_xzzz[i] * tbe_0 + 2.0 * tr_yzzz_xyzzz[i] * tke_0 - tr_zzz_xzzz[i];

        tr_0_0_y_yzzz_yyyy[i] = 2.0 * tr_yyzzz_yyyy[i] * tbe_0 + 2.0 * tr_yzzz_yyyyy[i] * tke_0 - tr_zzz_yyyy[i] - 4.0 * tr_yzzz_yyy[i];

        tr_0_0_y_yzzz_yyyz[i] = 2.0 * tr_yyzzz_yyyz[i] * tbe_0 + 2.0 * tr_yzzz_yyyyz[i] * tke_0 - tr_zzz_yyyz[i] - 3.0 * tr_yzzz_yyz[i];

        tr_0_0_y_yzzz_yyzz[i] = 2.0 * tr_yyzzz_yyzz[i] * tbe_0 + 2.0 * tr_yzzz_yyyzz[i] * tke_0 - tr_zzz_yyzz[i] - 2.0 * tr_yzzz_yzz[i];

        tr_0_0_y_yzzz_yzzz[i] = 2.0 * tr_yyzzz_yzzz[i] * tbe_0 + 2.0 * tr_yzzz_yyzzz[i] * tke_0 - tr_zzz_yzzz[i] - tr_yzzz_zzz[i];

        tr_0_0_y_yzzz_zzzz[i] = 2.0 * tr_yyzzz_zzzz[i] * tbe_0 + 2.0 * tr_yzzz_yzzzz[i] * tke_0 - tr_zzz_zzzz[i];
    }

    // Set up 435-450 components of targeted buffer : GG

    auto tr_0_0_y_zzzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 435);

    auto tr_0_0_y_zzzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 436);

    auto tr_0_0_y_zzzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 437);

    auto tr_0_0_y_zzzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 438);

    auto tr_0_0_y_zzzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 439);

    auto tr_0_0_y_zzzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 440);

    auto tr_0_0_y_zzzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 441);

    auto tr_0_0_y_zzzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 442);

    auto tr_0_0_y_zzzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 443);

    auto tr_0_0_y_zzzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 444);

    auto tr_0_0_y_zzzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 445);

    auto tr_0_0_y_zzzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 446);

    auto tr_0_0_y_zzzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 447);

    auto tr_0_0_y_zzzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 448);

    auto tr_0_0_y_zzzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 449);

    #pragma omp simd aligned(tr_0_0_y_zzzz_xxxx, tr_0_0_y_zzzz_xxxy, tr_0_0_y_zzzz_xxxz, tr_0_0_y_zzzz_xxyy, tr_0_0_y_zzzz_xxyz, tr_0_0_y_zzzz_xxzz, tr_0_0_y_zzzz_xyyy, tr_0_0_y_zzzz_xyyz, tr_0_0_y_zzzz_xyzz, tr_0_0_y_zzzz_xzzz, tr_0_0_y_zzzz_yyyy, tr_0_0_y_zzzz_yyyz, tr_0_0_y_zzzz_yyzz, tr_0_0_y_zzzz_yzzz, tr_0_0_y_zzzz_zzzz, tr_yzzzz_xxxx, tr_yzzzz_xxxy, tr_yzzzz_xxxz, tr_yzzzz_xxyy, tr_yzzzz_xxyz, tr_yzzzz_xxzz, tr_yzzzz_xyyy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xzzz, tr_yzzzz_yyyy, tr_yzzzz_yyyz, tr_yzzzz_yyzz, tr_yzzzz_yzzz, tr_yzzzz_zzzz, tr_zzzz_xxx, tr_zzzz_xxxxy, tr_zzzz_xxxyy, tr_zzzz_xxxyz, tr_zzzz_xxy, tr_zzzz_xxyyy, tr_zzzz_xxyyz, tr_zzzz_xxyzz, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyyyy, tr_zzzz_xyyyz, tr_zzzz_xyyzz, tr_zzzz_xyz, tr_zzzz_xyzzz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyyyy, tr_zzzz_yyyyz, tr_zzzz_yyyzz, tr_zzzz_yyz, tr_zzzz_yyzzz, tr_zzzz_yzz, tr_zzzz_yzzzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zzzz_xxxx[i] = 2.0 * tr_yzzzz_xxxx[i] * tbe_0 + 2.0 * tr_zzzz_xxxxy[i] * tke_0;

        tr_0_0_y_zzzz_xxxy[i] = 2.0 * tr_yzzzz_xxxy[i] * tbe_0 + 2.0 * tr_zzzz_xxxyy[i] * tke_0 - tr_zzzz_xxx[i];

        tr_0_0_y_zzzz_xxxz[i] = 2.0 * tr_yzzzz_xxxz[i] * tbe_0 + 2.0 * tr_zzzz_xxxyz[i] * tke_0;

        tr_0_0_y_zzzz_xxyy[i] = 2.0 * tr_yzzzz_xxyy[i] * tbe_0 + 2.0 * tr_zzzz_xxyyy[i] * tke_0 - 2.0 * tr_zzzz_xxy[i];

        tr_0_0_y_zzzz_xxyz[i] = 2.0 * tr_yzzzz_xxyz[i] * tbe_0 + 2.0 * tr_zzzz_xxyyz[i] * tke_0 - tr_zzzz_xxz[i];

        tr_0_0_y_zzzz_xxzz[i] = 2.0 * tr_yzzzz_xxzz[i] * tbe_0 + 2.0 * tr_zzzz_xxyzz[i] * tke_0;

        tr_0_0_y_zzzz_xyyy[i] = 2.0 * tr_yzzzz_xyyy[i] * tbe_0 + 2.0 * tr_zzzz_xyyyy[i] * tke_0 - 3.0 * tr_zzzz_xyy[i];

        tr_0_0_y_zzzz_xyyz[i] = 2.0 * tr_yzzzz_xyyz[i] * tbe_0 + 2.0 * tr_zzzz_xyyyz[i] * tke_0 - 2.0 * tr_zzzz_xyz[i];

        tr_0_0_y_zzzz_xyzz[i] = 2.0 * tr_yzzzz_xyzz[i] * tbe_0 + 2.0 * tr_zzzz_xyyzz[i] * tke_0 - tr_zzzz_xzz[i];

        tr_0_0_y_zzzz_xzzz[i] = 2.0 * tr_yzzzz_xzzz[i] * tbe_0 + 2.0 * tr_zzzz_xyzzz[i] * tke_0;

        tr_0_0_y_zzzz_yyyy[i] = 2.0 * tr_yzzzz_yyyy[i] * tbe_0 + 2.0 * tr_zzzz_yyyyy[i] * tke_0 - 4.0 * tr_zzzz_yyy[i];

        tr_0_0_y_zzzz_yyyz[i] = 2.0 * tr_yzzzz_yyyz[i] * tbe_0 + 2.0 * tr_zzzz_yyyyz[i] * tke_0 - 3.0 * tr_zzzz_yyz[i];

        tr_0_0_y_zzzz_yyzz[i] = 2.0 * tr_yzzzz_yyzz[i] * tbe_0 + 2.0 * tr_zzzz_yyyzz[i] * tke_0 - 2.0 * tr_zzzz_yzz[i];

        tr_0_0_y_zzzz_yzzz[i] = 2.0 * tr_yzzzz_yzzz[i] * tbe_0 + 2.0 * tr_zzzz_yyzzz[i] * tke_0 - tr_zzzz_zzz[i];

        tr_0_0_y_zzzz_zzzz[i] = 2.0 * tr_yzzzz_zzzz[i] * tbe_0 + 2.0 * tr_zzzz_yzzzz[i] * tke_0;
    }

    // Set up 450-465 components of targeted buffer : GG

    auto tr_0_0_z_xxxx_xxxx = pbuffer.data(idx_op_geom_010_gg + 450);

    auto tr_0_0_z_xxxx_xxxy = pbuffer.data(idx_op_geom_010_gg + 451);

    auto tr_0_0_z_xxxx_xxxz = pbuffer.data(idx_op_geom_010_gg + 452);

    auto tr_0_0_z_xxxx_xxyy = pbuffer.data(idx_op_geom_010_gg + 453);

    auto tr_0_0_z_xxxx_xxyz = pbuffer.data(idx_op_geom_010_gg + 454);

    auto tr_0_0_z_xxxx_xxzz = pbuffer.data(idx_op_geom_010_gg + 455);

    auto tr_0_0_z_xxxx_xyyy = pbuffer.data(idx_op_geom_010_gg + 456);

    auto tr_0_0_z_xxxx_xyyz = pbuffer.data(idx_op_geom_010_gg + 457);

    auto tr_0_0_z_xxxx_xyzz = pbuffer.data(idx_op_geom_010_gg + 458);

    auto tr_0_0_z_xxxx_xzzz = pbuffer.data(idx_op_geom_010_gg + 459);

    auto tr_0_0_z_xxxx_yyyy = pbuffer.data(idx_op_geom_010_gg + 460);

    auto tr_0_0_z_xxxx_yyyz = pbuffer.data(idx_op_geom_010_gg + 461);

    auto tr_0_0_z_xxxx_yyzz = pbuffer.data(idx_op_geom_010_gg + 462);

    auto tr_0_0_z_xxxx_yzzz = pbuffer.data(idx_op_geom_010_gg + 463);

    auto tr_0_0_z_xxxx_zzzz = pbuffer.data(idx_op_geom_010_gg + 464);

    #pragma omp simd aligned(tr_0_0_z_xxxx_xxxx, tr_0_0_z_xxxx_xxxy, tr_0_0_z_xxxx_xxxz, tr_0_0_z_xxxx_xxyy, tr_0_0_z_xxxx_xxyz, tr_0_0_z_xxxx_xxzz, tr_0_0_z_xxxx_xyyy, tr_0_0_z_xxxx_xyyz, tr_0_0_z_xxxx_xyzz, tr_0_0_z_xxxx_xzzz, tr_0_0_z_xxxx_yyyy, tr_0_0_z_xxxx_yyyz, tr_0_0_z_xxxx_yyzz, tr_0_0_z_xxxx_yzzz, tr_0_0_z_xxxx_zzzz, tr_xxxx_xxx, tr_xxxx_xxxxz, tr_xxxx_xxxyz, tr_xxxx_xxxzz, tr_xxxx_xxy, tr_xxxx_xxyyz, tr_xxxx_xxyzz, tr_xxxx_xxz, tr_xxxx_xxzzz, tr_xxxx_xyy, tr_xxxx_xyyyz, tr_xxxx_xyyzz, tr_xxxx_xyz, tr_xxxx_xyzzz, tr_xxxx_xzz, tr_xxxx_xzzzz, tr_xxxx_yyy, tr_xxxx_yyyyz, tr_xxxx_yyyzz, tr_xxxx_yyz, tr_xxxx_yyzzz, tr_xxxx_yzz, tr_xxxx_yzzzz, tr_xxxx_zzz, tr_xxxx_zzzzz, tr_xxxxz_xxxx, tr_xxxxz_xxxy, tr_xxxxz_xxxz, tr_xxxxz_xxyy, tr_xxxxz_xxyz, tr_xxxxz_xxzz, tr_xxxxz_xyyy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xzzz, tr_xxxxz_yyyy, tr_xxxxz_yyyz, tr_xxxxz_yyzz, tr_xxxxz_yzzz, tr_xxxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxx_xxxx[i] = 2.0 * tr_xxxxz_xxxx[i] * tbe_0 + 2.0 * tr_xxxx_xxxxz[i] * tke_0;

        tr_0_0_z_xxxx_xxxy[i] = 2.0 * tr_xxxxz_xxxy[i] * tbe_0 + 2.0 * tr_xxxx_xxxyz[i] * tke_0;

        tr_0_0_z_xxxx_xxxz[i] = 2.0 * tr_xxxxz_xxxz[i] * tbe_0 + 2.0 * tr_xxxx_xxxzz[i] * tke_0 - tr_xxxx_xxx[i];

        tr_0_0_z_xxxx_xxyy[i] = 2.0 * tr_xxxxz_xxyy[i] * tbe_0 + 2.0 * tr_xxxx_xxyyz[i] * tke_0;

        tr_0_0_z_xxxx_xxyz[i] = 2.0 * tr_xxxxz_xxyz[i] * tbe_0 + 2.0 * tr_xxxx_xxyzz[i] * tke_0 - tr_xxxx_xxy[i];

        tr_0_0_z_xxxx_xxzz[i] = 2.0 * tr_xxxxz_xxzz[i] * tbe_0 + 2.0 * tr_xxxx_xxzzz[i] * tke_0 - 2.0 * tr_xxxx_xxz[i];

        tr_0_0_z_xxxx_xyyy[i] = 2.0 * tr_xxxxz_xyyy[i] * tbe_0 + 2.0 * tr_xxxx_xyyyz[i] * tke_0;

        tr_0_0_z_xxxx_xyyz[i] = 2.0 * tr_xxxxz_xyyz[i] * tbe_0 + 2.0 * tr_xxxx_xyyzz[i] * tke_0 - tr_xxxx_xyy[i];

        tr_0_0_z_xxxx_xyzz[i] = 2.0 * tr_xxxxz_xyzz[i] * tbe_0 + 2.0 * tr_xxxx_xyzzz[i] * tke_0 - 2.0 * tr_xxxx_xyz[i];

        tr_0_0_z_xxxx_xzzz[i] = 2.0 * tr_xxxxz_xzzz[i] * tbe_0 + 2.0 * tr_xxxx_xzzzz[i] * tke_0 - 3.0 * tr_xxxx_xzz[i];

        tr_0_0_z_xxxx_yyyy[i] = 2.0 * tr_xxxxz_yyyy[i] * tbe_0 + 2.0 * tr_xxxx_yyyyz[i] * tke_0;

        tr_0_0_z_xxxx_yyyz[i] = 2.0 * tr_xxxxz_yyyz[i] * tbe_0 + 2.0 * tr_xxxx_yyyzz[i] * tke_0 - tr_xxxx_yyy[i];

        tr_0_0_z_xxxx_yyzz[i] = 2.0 * tr_xxxxz_yyzz[i] * tbe_0 + 2.0 * tr_xxxx_yyzzz[i] * tke_0 - 2.0 * tr_xxxx_yyz[i];

        tr_0_0_z_xxxx_yzzz[i] = 2.0 * tr_xxxxz_yzzz[i] * tbe_0 + 2.0 * tr_xxxx_yzzzz[i] * tke_0 - 3.0 * tr_xxxx_yzz[i];

        tr_0_0_z_xxxx_zzzz[i] = 2.0 * tr_xxxxz_zzzz[i] * tbe_0 + 2.0 * tr_xxxx_zzzzz[i] * tke_0 - 4.0 * tr_xxxx_zzz[i];
    }

    // Set up 465-480 components of targeted buffer : GG

    auto tr_0_0_z_xxxy_xxxx = pbuffer.data(idx_op_geom_010_gg + 465);

    auto tr_0_0_z_xxxy_xxxy = pbuffer.data(idx_op_geom_010_gg + 466);

    auto tr_0_0_z_xxxy_xxxz = pbuffer.data(idx_op_geom_010_gg + 467);

    auto tr_0_0_z_xxxy_xxyy = pbuffer.data(idx_op_geom_010_gg + 468);

    auto tr_0_0_z_xxxy_xxyz = pbuffer.data(idx_op_geom_010_gg + 469);

    auto tr_0_0_z_xxxy_xxzz = pbuffer.data(idx_op_geom_010_gg + 470);

    auto tr_0_0_z_xxxy_xyyy = pbuffer.data(idx_op_geom_010_gg + 471);

    auto tr_0_0_z_xxxy_xyyz = pbuffer.data(idx_op_geom_010_gg + 472);

    auto tr_0_0_z_xxxy_xyzz = pbuffer.data(idx_op_geom_010_gg + 473);

    auto tr_0_0_z_xxxy_xzzz = pbuffer.data(idx_op_geom_010_gg + 474);

    auto tr_0_0_z_xxxy_yyyy = pbuffer.data(idx_op_geom_010_gg + 475);

    auto tr_0_0_z_xxxy_yyyz = pbuffer.data(idx_op_geom_010_gg + 476);

    auto tr_0_0_z_xxxy_yyzz = pbuffer.data(idx_op_geom_010_gg + 477);

    auto tr_0_0_z_xxxy_yzzz = pbuffer.data(idx_op_geom_010_gg + 478);

    auto tr_0_0_z_xxxy_zzzz = pbuffer.data(idx_op_geom_010_gg + 479);

    #pragma omp simd aligned(tr_0_0_z_xxxy_xxxx, tr_0_0_z_xxxy_xxxy, tr_0_0_z_xxxy_xxxz, tr_0_0_z_xxxy_xxyy, tr_0_0_z_xxxy_xxyz, tr_0_0_z_xxxy_xxzz, tr_0_0_z_xxxy_xyyy, tr_0_0_z_xxxy_xyyz, tr_0_0_z_xxxy_xyzz, tr_0_0_z_xxxy_xzzz, tr_0_0_z_xxxy_yyyy, tr_0_0_z_xxxy_yyyz, tr_0_0_z_xxxy_yyzz, tr_0_0_z_xxxy_yzzz, tr_0_0_z_xxxy_zzzz, tr_xxxy_xxx, tr_xxxy_xxxxz, tr_xxxy_xxxyz, tr_xxxy_xxxzz, tr_xxxy_xxy, tr_xxxy_xxyyz, tr_xxxy_xxyzz, tr_xxxy_xxz, tr_xxxy_xxzzz, tr_xxxy_xyy, tr_xxxy_xyyyz, tr_xxxy_xyyzz, tr_xxxy_xyz, tr_xxxy_xyzzz, tr_xxxy_xzz, tr_xxxy_xzzzz, tr_xxxy_yyy, tr_xxxy_yyyyz, tr_xxxy_yyyzz, tr_xxxy_yyz, tr_xxxy_yyzzz, tr_xxxy_yzz, tr_xxxy_yzzzz, tr_xxxy_zzz, tr_xxxy_zzzzz, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxy_xxxx[i] = 2.0 * tr_xxxyz_xxxx[i] * tbe_0 + 2.0 * tr_xxxy_xxxxz[i] * tke_0;

        tr_0_0_z_xxxy_xxxy[i] = 2.0 * tr_xxxyz_xxxy[i] * tbe_0 + 2.0 * tr_xxxy_xxxyz[i] * tke_0;

        tr_0_0_z_xxxy_xxxz[i] = 2.0 * tr_xxxyz_xxxz[i] * tbe_0 + 2.0 * tr_xxxy_xxxzz[i] * tke_0 - tr_xxxy_xxx[i];

        tr_0_0_z_xxxy_xxyy[i] = 2.0 * tr_xxxyz_xxyy[i] * tbe_0 + 2.0 * tr_xxxy_xxyyz[i] * tke_0;

        tr_0_0_z_xxxy_xxyz[i] = 2.0 * tr_xxxyz_xxyz[i] * tbe_0 + 2.0 * tr_xxxy_xxyzz[i] * tke_0 - tr_xxxy_xxy[i];

        tr_0_0_z_xxxy_xxzz[i] = 2.0 * tr_xxxyz_xxzz[i] * tbe_0 + 2.0 * tr_xxxy_xxzzz[i] * tke_0 - 2.0 * tr_xxxy_xxz[i];

        tr_0_0_z_xxxy_xyyy[i] = 2.0 * tr_xxxyz_xyyy[i] * tbe_0 + 2.0 * tr_xxxy_xyyyz[i] * tke_0;

        tr_0_0_z_xxxy_xyyz[i] = 2.0 * tr_xxxyz_xyyz[i] * tbe_0 + 2.0 * tr_xxxy_xyyzz[i] * tke_0 - tr_xxxy_xyy[i];

        tr_0_0_z_xxxy_xyzz[i] = 2.0 * tr_xxxyz_xyzz[i] * tbe_0 + 2.0 * tr_xxxy_xyzzz[i] * tke_0 - 2.0 * tr_xxxy_xyz[i];

        tr_0_0_z_xxxy_xzzz[i] = 2.0 * tr_xxxyz_xzzz[i] * tbe_0 + 2.0 * tr_xxxy_xzzzz[i] * tke_0 - 3.0 * tr_xxxy_xzz[i];

        tr_0_0_z_xxxy_yyyy[i] = 2.0 * tr_xxxyz_yyyy[i] * tbe_0 + 2.0 * tr_xxxy_yyyyz[i] * tke_0;

        tr_0_0_z_xxxy_yyyz[i] = 2.0 * tr_xxxyz_yyyz[i] * tbe_0 + 2.0 * tr_xxxy_yyyzz[i] * tke_0 - tr_xxxy_yyy[i];

        tr_0_0_z_xxxy_yyzz[i] = 2.0 * tr_xxxyz_yyzz[i] * tbe_0 + 2.0 * tr_xxxy_yyzzz[i] * tke_0 - 2.0 * tr_xxxy_yyz[i];

        tr_0_0_z_xxxy_yzzz[i] = 2.0 * tr_xxxyz_yzzz[i] * tbe_0 + 2.0 * tr_xxxy_yzzzz[i] * tke_0 - 3.0 * tr_xxxy_yzz[i];

        tr_0_0_z_xxxy_zzzz[i] = 2.0 * tr_xxxyz_zzzz[i] * tbe_0 + 2.0 * tr_xxxy_zzzzz[i] * tke_0 - 4.0 * tr_xxxy_zzz[i];
    }

    // Set up 480-495 components of targeted buffer : GG

    auto tr_0_0_z_xxxz_xxxx = pbuffer.data(idx_op_geom_010_gg + 480);

    auto tr_0_0_z_xxxz_xxxy = pbuffer.data(idx_op_geom_010_gg + 481);

    auto tr_0_0_z_xxxz_xxxz = pbuffer.data(idx_op_geom_010_gg + 482);

    auto tr_0_0_z_xxxz_xxyy = pbuffer.data(idx_op_geom_010_gg + 483);

    auto tr_0_0_z_xxxz_xxyz = pbuffer.data(idx_op_geom_010_gg + 484);

    auto tr_0_0_z_xxxz_xxzz = pbuffer.data(idx_op_geom_010_gg + 485);

    auto tr_0_0_z_xxxz_xyyy = pbuffer.data(idx_op_geom_010_gg + 486);

    auto tr_0_0_z_xxxz_xyyz = pbuffer.data(idx_op_geom_010_gg + 487);

    auto tr_0_0_z_xxxz_xyzz = pbuffer.data(idx_op_geom_010_gg + 488);

    auto tr_0_0_z_xxxz_xzzz = pbuffer.data(idx_op_geom_010_gg + 489);

    auto tr_0_0_z_xxxz_yyyy = pbuffer.data(idx_op_geom_010_gg + 490);

    auto tr_0_0_z_xxxz_yyyz = pbuffer.data(idx_op_geom_010_gg + 491);

    auto tr_0_0_z_xxxz_yyzz = pbuffer.data(idx_op_geom_010_gg + 492);

    auto tr_0_0_z_xxxz_yzzz = pbuffer.data(idx_op_geom_010_gg + 493);

    auto tr_0_0_z_xxxz_zzzz = pbuffer.data(idx_op_geom_010_gg + 494);

    #pragma omp simd aligned(tr_0_0_z_xxxz_xxxx, tr_0_0_z_xxxz_xxxy, tr_0_0_z_xxxz_xxxz, tr_0_0_z_xxxz_xxyy, tr_0_0_z_xxxz_xxyz, tr_0_0_z_xxxz_xxzz, tr_0_0_z_xxxz_xyyy, tr_0_0_z_xxxz_xyyz, tr_0_0_z_xxxz_xyzz, tr_0_0_z_xxxz_xzzz, tr_0_0_z_xxxz_yyyy, tr_0_0_z_xxxz_yyyz, tr_0_0_z_xxxz_yyzz, tr_0_0_z_xxxz_yzzz, tr_0_0_z_xxxz_zzzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxz_xxx, tr_xxxz_xxxxz, tr_xxxz_xxxyz, tr_xxxz_xxxzz, tr_xxxz_xxy, tr_xxxz_xxyyz, tr_xxxz_xxyzz, tr_xxxz_xxz, tr_xxxz_xxzzz, tr_xxxz_xyy, tr_xxxz_xyyyz, tr_xxxz_xyyzz, tr_xxxz_xyz, tr_xxxz_xyzzz, tr_xxxz_xzz, tr_xxxz_xzzzz, tr_xxxz_yyy, tr_xxxz_yyyyz, tr_xxxz_yyyzz, tr_xxxz_yyz, tr_xxxz_yyzzz, tr_xxxz_yzz, tr_xxxz_yzzzz, tr_xxxz_zzz, tr_xxxz_zzzzz, tr_xxxzz_xxxx, tr_xxxzz_xxxy, tr_xxxzz_xxxz, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xzzz, tr_xxxzz_yyyy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yzzz, tr_xxxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxz_xxxx[i] = 2.0 * tr_xxxzz_xxxx[i] * tbe_0 + 2.0 * tr_xxxz_xxxxz[i] * tke_0 - tr_xxx_xxxx[i];

        tr_0_0_z_xxxz_xxxy[i] = 2.0 * tr_xxxzz_xxxy[i] * tbe_0 + 2.0 * tr_xxxz_xxxyz[i] * tke_0 - tr_xxx_xxxy[i];

        tr_0_0_z_xxxz_xxxz[i] = 2.0 * tr_xxxzz_xxxz[i] * tbe_0 + 2.0 * tr_xxxz_xxxzz[i] * tke_0 - tr_xxx_xxxz[i] - tr_xxxz_xxx[i];

        tr_0_0_z_xxxz_xxyy[i] = 2.0 * tr_xxxzz_xxyy[i] * tbe_0 + 2.0 * tr_xxxz_xxyyz[i] * tke_0 - tr_xxx_xxyy[i];

        tr_0_0_z_xxxz_xxyz[i] = 2.0 * tr_xxxzz_xxyz[i] * tbe_0 + 2.0 * tr_xxxz_xxyzz[i] * tke_0 - tr_xxx_xxyz[i] - tr_xxxz_xxy[i];

        tr_0_0_z_xxxz_xxzz[i] = 2.0 * tr_xxxzz_xxzz[i] * tbe_0 + 2.0 * tr_xxxz_xxzzz[i] * tke_0 - tr_xxx_xxzz[i] - 2.0 * tr_xxxz_xxz[i];

        tr_0_0_z_xxxz_xyyy[i] = 2.0 * tr_xxxzz_xyyy[i] * tbe_0 + 2.0 * tr_xxxz_xyyyz[i] * tke_0 - tr_xxx_xyyy[i];

        tr_0_0_z_xxxz_xyyz[i] = 2.0 * tr_xxxzz_xyyz[i] * tbe_0 + 2.0 * tr_xxxz_xyyzz[i] * tke_0 - tr_xxx_xyyz[i] - tr_xxxz_xyy[i];

        tr_0_0_z_xxxz_xyzz[i] = 2.0 * tr_xxxzz_xyzz[i] * tbe_0 + 2.0 * tr_xxxz_xyzzz[i] * tke_0 - tr_xxx_xyzz[i] - 2.0 * tr_xxxz_xyz[i];

        tr_0_0_z_xxxz_xzzz[i] = 2.0 * tr_xxxzz_xzzz[i] * tbe_0 + 2.0 * tr_xxxz_xzzzz[i] * tke_0 - tr_xxx_xzzz[i] - 3.0 * tr_xxxz_xzz[i];

        tr_0_0_z_xxxz_yyyy[i] = 2.0 * tr_xxxzz_yyyy[i] * tbe_0 + 2.0 * tr_xxxz_yyyyz[i] * tke_0 - tr_xxx_yyyy[i];

        tr_0_0_z_xxxz_yyyz[i] = 2.0 * tr_xxxzz_yyyz[i] * tbe_0 + 2.0 * tr_xxxz_yyyzz[i] * tke_0 - tr_xxx_yyyz[i] - tr_xxxz_yyy[i];

        tr_0_0_z_xxxz_yyzz[i] = 2.0 * tr_xxxzz_yyzz[i] * tbe_0 + 2.0 * tr_xxxz_yyzzz[i] * tke_0 - tr_xxx_yyzz[i] - 2.0 * tr_xxxz_yyz[i];

        tr_0_0_z_xxxz_yzzz[i] = 2.0 * tr_xxxzz_yzzz[i] * tbe_0 + 2.0 * tr_xxxz_yzzzz[i] * tke_0 - tr_xxx_yzzz[i] - 3.0 * tr_xxxz_yzz[i];

        tr_0_0_z_xxxz_zzzz[i] = 2.0 * tr_xxxzz_zzzz[i] * tbe_0 + 2.0 * tr_xxxz_zzzzz[i] * tke_0 - tr_xxx_zzzz[i] - 4.0 * tr_xxxz_zzz[i];
    }

    // Set up 495-510 components of targeted buffer : GG

    auto tr_0_0_z_xxyy_xxxx = pbuffer.data(idx_op_geom_010_gg + 495);

    auto tr_0_0_z_xxyy_xxxy = pbuffer.data(idx_op_geom_010_gg + 496);

    auto tr_0_0_z_xxyy_xxxz = pbuffer.data(idx_op_geom_010_gg + 497);

    auto tr_0_0_z_xxyy_xxyy = pbuffer.data(idx_op_geom_010_gg + 498);

    auto tr_0_0_z_xxyy_xxyz = pbuffer.data(idx_op_geom_010_gg + 499);

    auto tr_0_0_z_xxyy_xxzz = pbuffer.data(idx_op_geom_010_gg + 500);

    auto tr_0_0_z_xxyy_xyyy = pbuffer.data(idx_op_geom_010_gg + 501);

    auto tr_0_0_z_xxyy_xyyz = pbuffer.data(idx_op_geom_010_gg + 502);

    auto tr_0_0_z_xxyy_xyzz = pbuffer.data(idx_op_geom_010_gg + 503);

    auto tr_0_0_z_xxyy_xzzz = pbuffer.data(idx_op_geom_010_gg + 504);

    auto tr_0_0_z_xxyy_yyyy = pbuffer.data(idx_op_geom_010_gg + 505);

    auto tr_0_0_z_xxyy_yyyz = pbuffer.data(idx_op_geom_010_gg + 506);

    auto tr_0_0_z_xxyy_yyzz = pbuffer.data(idx_op_geom_010_gg + 507);

    auto tr_0_0_z_xxyy_yzzz = pbuffer.data(idx_op_geom_010_gg + 508);

    auto tr_0_0_z_xxyy_zzzz = pbuffer.data(idx_op_geom_010_gg + 509);

    #pragma omp simd aligned(tr_0_0_z_xxyy_xxxx, tr_0_0_z_xxyy_xxxy, tr_0_0_z_xxyy_xxxz, tr_0_0_z_xxyy_xxyy, tr_0_0_z_xxyy_xxyz, tr_0_0_z_xxyy_xxzz, tr_0_0_z_xxyy_xyyy, tr_0_0_z_xxyy_xyyz, tr_0_0_z_xxyy_xyzz, tr_0_0_z_xxyy_xzzz, tr_0_0_z_xxyy_yyyy, tr_0_0_z_xxyy_yyyz, tr_0_0_z_xxyy_yyzz, tr_0_0_z_xxyy_yzzz, tr_0_0_z_xxyy_zzzz, tr_xxyy_xxx, tr_xxyy_xxxxz, tr_xxyy_xxxyz, tr_xxyy_xxxzz, tr_xxyy_xxy, tr_xxyy_xxyyz, tr_xxyy_xxyzz, tr_xxyy_xxz, tr_xxyy_xxzzz, tr_xxyy_xyy, tr_xxyy_xyyyz, tr_xxyy_xyyzz, tr_xxyy_xyz, tr_xxyy_xyzzz, tr_xxyy_xzz, tr_xxyy_xzzzz, tr_xxyy_yyy, tr_xxyy_yyyyz, tr_xxyy_yyyzz, tr_xxyy_yyz, tr_xxyy_yyzzz, tr_xxyy_yzz, tr_xxyy_yzzzz, tr_xxyy_zzz, tr_xxyy_zzzzz, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxyy_xxxx[i] = 2.0 * tr_xxyyz_xxxx[i] * tbe_0 + 2.0 * tr_xxyy_xxxxz[i] * tke_0;

        tr_0_0_z_xxyy_xxxy[i] = 2.0 * tr_xxyyz_xxxy[i] * tbe_0 + 2.0 * tr_xxyy_xxxyz[i] * tke_0;

        tr_0_0_z_xxyy_xxxz[i] = 2.0 * tr_xxyyz_xxxz[i] * tbe_0 + 2.0 * tr_xxyy_xxxzz[i] * tke_0 - tr_xxyy_xxx[i];

        tr_0_0_z_xxyy_xxyy[i] = 2.0 * tr_xxyyz_xxyy[i] * tbe_0 + 2.0 * tr_xxyy_xxyyz[i] * tke_0;

        tr_0_0_z_xxyy_xxyz[i] = 2.0 * tr_xxyyz_xxyz[i] * tbe_0 + 2.0 * tr_xxyy_xxyzz[i] * tke_0 - tr_xxyy_xxy[i];

        tr_0_0_z_xxyy_xxzz[i] = 2.0 * tr_xxyyz_xxzz[i] * tbe_0 + 2.0 * tr_xxyy_xxzzz[i] * tke_0 - 2.0 * tr_xxyy_xxz[i];

        tr_0_0_z_xxyy_xyyy[i] = 2.0 * tr_xxyyz_xyyy[i] * tbe_0 + 2.0 * tr_xxyy_xyyyz[i] * tke_0;

        tr_0_0_z_xxyy_xyyz[i] = 2.0 * tr_xxyyz_xyyz[i] * tbe_0 + 2.0 * tr_xxyy_xyyzz[i] * tke_0 - tr_xxyy_xyy[i];

        tr_0_0_z_xxyy_xyzz[i] = 2.0 * tr_xxyyz_xyzz[i] * tbe_0 + 2.0 * tr_xxyy_xyzzz[i] * tke_0 - 2.0 * tr_xxyy_xyz[i];

        tr_0_0_z_xxyy_xzzz[i] = 2.0 * tr_xxyyz_xzzz[i] * tbe_0 + 2.0 * tr_xxyy_xzzzz[i] * tke_0 - 3.0 * tr_xxyy_xzz[i];

        tr_0_0_z_xxyy_yyyy[i] = 2.0 * tr_xxyyz_yyyy[i] * tbe_0 + 2.0 * tr_xxyy_yyyyz[i] * tke_0;

        tr_0_0_z_xxyy_yyyz[i] = 2.0 * tr_xxyyz_yyyz[i] * tbe_0 + 2.0 * tr_xxyy_yyyzz[i] * tke_0 - tr_xxyy_yyy[i];

        tr_0_0_z_xxyy_yyzz[i] = 2.0 * tr_xxyyz_yyzz[i] * tbe_0 + 2.0 * tr_xxyy_yyzzz[i] * tke_0 - 2.0 * tr_xxyy_yyz[i];

        tr_0_0_z_xxyy_yzzz[i] = 2.0 * tr_xxyyz_yzzz[i] * tbe_0 + 2.0 * tr_xxyy_yzzzz[i] * tke_0 - 3.0 * tr_xxyy_yzz[i];

        tr_0_0_z_xxyy_zzzz[i] = 2.0 * tr_xxyyz_zzzz[i] * tbe_0 + 2.0 * tr_xxyy_zzzzz[i] * tke_0 - 4.0 * tr_xxyy_zzz[i];
    }

    // Set up 510-525 components of targeted buffer : GG

    auto tr_0_0_z_xxyz_xxxx = pbuffer.data(idx_op_geom_010_gg + 510);

    auto tr_0_0_z_xxyz_xxxy = pbuffer.data(idx_op_geom_010_gg + 511);

    auto tr_0_0_z_xxyz_xxxz = pbuffer.data(idx_op_geom_010_gg + 512);

    auto tr_0_0_z_xxyz_xxyy = pbuffer.data(idx_op_geom_010_gg + 513);

    auto tr_0_0_z_xxyz_xxyz = pbuffer.data(idx_op_geom_010_gg + 514);

    auto tr_0_0_z_xxyz_xxzz = pbuffer.data(idx_op_geom_010_gg + 515);

    auto tr_0_0_z_xxyz_xyyy = pbuffer.data(idx_op_geom_010_gg + 516);

    auto tr_0_0_z_xxyz_xyyz = pbuffer.data(idx_op_geom_010_gg + 517);

    auto tr_0_0_z_xxyz_xyzz = pbuffer.data(idx_op_geom_010_gg + 518);

    auto tr_0_0_z_xxyz_xzzz = pbuffer.data(idx_op_geom_010_gg + 519);

    auto tr_0_0_z_xxyz_yyyy = pbuffer.data(idx_op_geom_010_gg + 520);

    auto tr_0_0_z_xxyz_yyyz = pbuffer.data(idx_op_geom_010_gg + 521);

    auto tr_0_0_z_xxyz_yyzz = pbuffer.data(idx_op_geom_010_gg + 522);

    auto tr_0_0_z_xxyz_yzzz = pbuffer.data(idx_op_geom_010_gg + 523);

    auto tr_0_0_z_xxyz_zzzz = pbuffer.data(idx_op_geom_010_gg + 524);

    #pragma omp simd aligned(tr_0_0_z_xxyz_xxxx, tr_0_0_z_xxyz_xxxy, tr_0_0_z_xxyz_xxxz, tr_0_0_z_xxyz_xxyy, tr_0_0_z_xxyz_xxyz, tr_0_0_z_xxyz_xxzz, tr_0_0_z_xxyz_xyyy, tr_0_0_z_xxyz_xyyz, tr_0_0_z_xxyz_xyzz, tr_0_0_z_xxyz_xzzz, tr_0_0_z_xxyz_yyyy, tr_0_0_z_xxyz_yyyz, tr_0_0_z_xxyz_yyzz, tr_0_0_z_xxyz_yzzz, tr_0_0_z_xxyz_zzzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxz, tr_xxyz_xxxyz, tr_xxyz_xxxzz, tr_xxyz_xxy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xxzzz, tr_xxyz_xyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_xzzzz, tr_xxyz_yyy, tr_xxyz_yyyyz, tr_xxyz_yyyzz, tr_xxyz_yyz, tr_xxyz_yyzzz, tr_xxyz_yzz, tr_xxyz_yzzzz, tr_xxyz_zzz, tr_xxyz_zzzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxyz_xxxx[i] = 2.0 * tr_xxyzz_xxxx[i] * tbe_0 + 2.0 * tr_xxyz_xxxxz[i] * tke_0 - tr_xxy_xxxx[i];

        tr_0_0_z_xxyz_xxxy[i] = 2.0 * tr_xxyzz_xxxy[i] * tbe_0 + 2.0 * tr_xxyz_xxxyz[i] * tke_0 - tr_xxy_xxxy[i];

        tr_0_0_z_xxyz_xxxz[i] = 2.0 * tr_xxyzz_xxxz[i] * tbe_0 + 2.0 * tr_xxyz_xxxzz[i] * tke_0 - tr_xxy_xxxz[i] - tr_xxyz_xxx[i];

        tr_0_0_z_xxyz_xxyy[i] = 2.0 * tr_xxyzz_xxyy[i] * tbe_0 + 2.0 * tr_xxyz_xxyyz[i] * tke_0 - tr_xxy_xxyy[i];

        tr_0_0_z_xxyz_xxyz[i] = 2.0 * tr_xxyzz_xxyz[i] * tbe_0 + 2.0 * tr_xxyz_xxyzz[i] * tke_0 - tr_xxy_xxyz[i] - tr_xxyz_xxy[i];

        tr_0_0_z_xxyz_xxzz[i] = 2.0 * tr_xxyzz_xxzz[i] * tbe_0 + 2.0 * tr_xxyz_xxzzz[i] * tke_0 - tr_xxy_xxzz[i] - 2.0 * tr_xxyz_xxz[i];

        tr_0_0_z_xxyz_xyyy[i] = 2.0 * tr_xxyzz_xyyy[i] * tbe_0 + 2.0 * tr_xxyz_xyyyz[i] * tke_0 - tr_xxy_xyyy[i];

        tr_0_0_z_xxyz_xyyz[i] = 2.0 * tr_xxyzz_xyyz[i] * tbe_0 + 2.0 * tr_xxyz_xyyzz[i] * tke_0 - tr_xxy_xyyz[i] - tr_xxyz_xyy[i];

        tr_0_0_z_xxyz_xyzz[i] = 2.0 * tr_xxyzz_xyzz[i] * tbe_0 + 2.0 * tr_xxyz_xyzzz[i] * tke_0 - tr_xxy_xyzz[i] - 2.0 * tr_xxyz_xyz[i];

        tr_0_0_z_xxyz_xzzz[i] = 2.0 * tr_xxyzz_xzzz[i] * tbe_0 + 2.0 * tr_xxyz_xzzzz[i] * tke_0 - tr_xxy_xzzz[i] - 3.0 * tr_xxyz_xzz[i];

        tr_0_0_z_xxyz_yyyy[i] = 2.0 * tr_xxyzz_yyyy[i] * tbe_0 + 2.0 * tr_xxyz_yyyyz[i] * tke_0 - tr_xxy_yyyy[i];

        tr_0_0_z_xxyz_yyyz[i] = 2.0 * tr_xxyzz_yyyz[i] * tbe_0 + 2.0 * tr_xxyz_yyyzz[i] * tke_0 - tr_xxy_yyyz[i] - tr_xxyz_yyy[i];

        tr_0_0_z_xxyz_yyzz[i] = 2.0 * tr_xxyzz_yyzz[i] * tbe_0 + 2.0 * tr_xxyz_yyzzz[i] * tke_0 - tr_xxy_yyzz[i] - 2.0 * tr_xxyz_yyz[i];

        tr_0_0_z_xxyz_yzzz[i] = 2.0 * tr_xxyzz_yzzz[i] * tbe_0 + 2.0 * tr_xxyz_yzzzz[i] * tke_0 - tr_xxy_yzzz[i] - 3.0 * tr_xxyz_yzz[i];

        tr_0_0_z_xxyz_zzzz[i] = 2.0 * tr_xxyzz_zzzz[i] * tbe_0 + 2.0 * tr_xxyz_zzzzz[i] * tke_0 - tr_xxy_zzzz[i] - 4.0 * tr_xxyz_zzz[i];
    }

    // Set up 525-540 components of targeted buffer : GG

    auto tr_0_0_z_xxzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 525);

    auto tr_0_0_z_xxzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 526);

    auto tr_0_0_z_xxzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 527);

    auto tr_0_0_z_xxzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 528);

    auto tr_0_0_z_xxzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 529);

    auto tr_0_0_z_xxzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 530);

    auto tr_0_0_z_xxzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 531);

    auto tr_0_0_z_xxzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 532);

    auto tr_0_0_z_xxzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 533);

    auto tr_0_0_z_xxzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 534);

    auto tr_0_0_z_xxzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 535);

    auto tr_0_0_z_xxzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 536);

    auto tr_0_0_z_xxzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 537);

    auto tr_0_0_z_xxzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 538);

    auto tr_0_0_z_xxzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 539);

    #pragma omp simd aligned(tr_0_0_z_xxzz_xxxx, tr_0_0_z_xxzz_xxxy, tr_0_0_z_xxzz_xxxz, tr_0_0_z_xxzz_xxyy, tr_0_0_z_xxzz_xxyz, tr_0_0_z_xxzz_xxzz, tr_0_0_z_xxzz_xyyy, tr_0_0_z_xxzz_xyyz, tr_0_0_z_xxzz_xyzz, tr_0_0_z_xxzz_xzzz, tr_0_0_z_xxzz_yyyy, tr_0_0_z_xxzz_yyyz, tr_0_0_z_xxzz_yyzz, tr_0_0_z_xxzz_yzzz, tr_0_0_z_xxzz_zzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xxzz_xxx, tr_xxzz_xxxxz, tr_xxzz_xxxyz, tr_xxzz_xxxzz, tr_xxzz_xxy, tr_xxzz_xxyyz, tr_xxzz_xxyzz, tr_xxzz_xxz, tr_xxzz_xxzzz, tr_xxzz_xyy, tr_xxzz_xyyyz, tr_xxzz_xyyzz, tr_xxzz_xyz, tr_xxzz_xyzzz, tr_xxzz_xzz, tr_xxzz_xzzzz, tr_xxzz_yyy, tr_xxzz_yyyyz, tr_xxzz_yyyzz, tr_xxzz_yyz, tr_xxzz_yyzzz, tr_xxzz_yzz, tr_xxzz_yzzzz, tr_xxzz_zzz, tr_xxzz_zzzzz, tr_xxzzz_xxxx, tr_xxzzz_xxxy, tr_xxzzz_xxxz, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xzzz, tr_xxzzz_yyyy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yzzz, tr_xxzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxzz_xxxx[i] = 2.0 * tr_xxzzz_xxxx[i] * tbe_0 + 2.0 * tr_xxzz_xxxxz[i] * tke_0 - 2.0 * tr_xxz_xxxx[i];

        tr_0_0_z_xxzz_xxxy[i] = 2.0 * tr_xxzzz_xxxy[i] * tbe_0 + 2.0 * tr_xxzz_xxxyz[i] * tke_0 - 2.0 * tr_xxz_xxxy[i];

        tr_0_0_z_xxzz_xxxz[i] = 2.0 * tr_xxzzz_xxxz[i] * tbe_0 + 2.0 * tr_xxzz_xxxzz[i] * tke_0 - 2.0 * tr_xxz_xxxz[i] - tr_xxzz_xxx[i];

        tr_0_0_z_xxzz_xxyy[i] = 2.0 * tr_xxzzz_xxyy[i] * tbe_0 + 2.0 * tr_xxzz_xxyyz[i] * tke_0 - 2.0 * tr_xxz_xxyy[i];

        tr_0_0_z_xxzz_xxyz[i] = 2.0 * tr_xxzzz_xxyz[i] * tbe_0 + 2.0 * tr_xxzz_xxyzz[i] * tke_0 - 2.0 * tr_xxz_xxyz[i] - tr_xxzz_xxy[i];

        tr_0_0_z_xxzz_xxzz[i] = 2.0 * tr_xxzzz_xxzz[i] * tbe_0 + 2.0 * tr_xxzz_xxzzz[i] * tke_0 - 2.0 * tr_xxz_xxzz[i] - 2.0 * tr_xxzz_xxz[i];

        tr_0_0_z_xxzz_xyyy[i] = 2.0 * tr_xxzzz_xyyy[i] * tbe_0 + 2.0 * tr_xxzz_xyyyz[i] * tke_0 - 2.0 * tr_xxz_xyyy[i];

        tr_0_0_z_xxzz_xyyz[i] = 2.0 * tr_xxzzz_xyyz[i] * tbe_0 + 2.0 * tr_xxzz_xyyzz[i] * tke_0 - 2.0 * tr_xxz_xyyz[i] - tr_xxzz_xyy[i];

        tr_0_0_z_xxzz_xyzz[i] = 2.0 * tr_xxzzz_xyzz[i] * tbe_0 + 2.0 * tr_xxzz_xyzzz[i] * tke_0 - 2.0 * tr_xxz_xyzz[i] - 2.0 * tr_xxzz_xyz[i];

        tr_0_0_z_xxzz_xzzz[i] = 2.0 * tr_xxzzz_xzzz[i] * tbe_0 + 2.0 * tr_xxzz_xzzzz[i] * tke_0 - 2.0 * tr_xxz_xzzz[i] - 3.0 * tr_xxzz_xzz[i];

        tr_0_0_z_xxzz_yyyy[i] = 2.0 * tr_xxzzz_yyyy[i] * tbe_0 + 2.0 * tr_xxzz_yyyyz[i] * tke_0 - 2.0 * tr_xxz_yyyy[i];

        tr_0_0_z_xxzz_yyyz[i] = 2.0 * tr_xxzzz_yyyz[i] * tbe_0 + 2.0 * tr_xxzz_yyyzz[i] * tke_0 - 2.0 * tr_xxz_yyyz[i] - tr_xxzz_yyy[i];

        tr_0_0_z_xxzz_yyzz[i] = 2.0 * tr_xxzzz_yyzz[i] * tbe_0 + 2.0 * tr_xxzz_yyzzz[i] * tke_0 - 2.0 * tr_xxz_yyzz[i] - 2.0 * tr_xxzz_yyz[i];

        tr_0_0_z_xxzz_yzzz[i] = 2.0 * tr_xxzzz_yzzz[i] * tbe_0 + 2.0 * tr_xxzz_yzzzz[i] * tke_0 - 2.0 * tr_xxz_yzzz[i] - 3.0 * tr_xxzz_yzz[i];

        tr_0_0_z_xxzz_zzzz[i] = 2.0 * tr_xxzzz_zzzz[i] * tbe_0 + 2.0 * tr_xxzz_zzzzz[i] * tke_0 - 2.0 * tr_xxz_zzzz[i] - 4.0 * tr_xxzz_zzz[i];
    }

    // Set up 540-555 components of targeted buffer : GG

    auto tr_0_0_z_xyyy_xxxx = pbuffer.data(idx_op_geom_010_gg + 540);

    auto tr_0_0_z_xyyy_xxxy = pbuffer.data(idx_op_geom_010_gg + 541);

    auto tr_0_0_z_xyyy_xxxz = pbuffer.data(idx_op_geom_010_gg + 542);

    auto tr_0_0_z_xyyy_xxyy = pbuffer.data(idx_op_geom_010_gg + 543);

    auto tr_0_0_z_xyyy_xxyz = pbuffer.data(idx_op_geom_010_gg + 544);

    auto tr_0_0_z_xyyy_xxzz = pbuffer.data(idx_op_geom_010_gg + 545);

    auto tr_0_0_z_xyyy_xyyy = pbuffer.data(idx_op_geom_010_gg + 546);

    auto tr_0_0_z_xyyy_xyyz = pbuffer.data(idx_op_geom_010_gg + 547);

    auto tr_0_0_z_xyyy_xyzz = pbuffer.data(idx_op_geom_010_gg + 548);

    auto tr_0_0_z_xyyy_xzzz = pbuffer.data(idx_op_geom_010_gg + 549);

    auto tr_0_0_z_xyyy_yyyy = pbuffer.data(idx_op_geom_010_gg + 550);

    auto tr_0_0_z_xyyy_yyyz = pbuffer.data(idx_op_geom_010_gg + 551);

    auto tr_0_0_z_xyyy_yyzz = pbuffer.data(idx_op_geom_010_gg + 552);

    auto tr_0_0_z_xyyy_yzzz = pbuffer.data(idx_op_geom_010_gg + 553);

    auto tr_0_0_z_xyyy_zzzz = pbuffer.data(idx_op_geom_010_gg + 554);

    #pragma omp simd aligned(tr_0_0_z_xyyy_xxxx, tr_0_0_z_xyyy_xxxy, tr_0_0_z_xyyy_xxxz, tr_0_0_z_xyyy_xxyy, tr_0_0_z_xyyy_xxyz, tr_0_0_z_xyyy_xxzz, tr_0_0_z_xyyy_xyyy, tr_0_0_z_xyyy_xyyz, tr_0_0_z_xyyy_xyzz, tr_0_0_z_xyyy_xzzz, tr_0_0_z_xyyy_yyyy, tr_0_0_z_xyyy_yyyz, tr_0_0_z_xyyy_yyzz, tr_0_0_z_xyyy_yzzz, tr_0_0_z_xyyy_zzzz, tr_xyyy_xxx, tr_xyyy_xxxxz, tr_xyyy_xxxyz, tr_xyyy_xxxzz, tr_xyyy_xxy, tr_xyyy_xxyyz, tr_xyyy_xxyzz, tr_xyyy_xxz, tr_xyyy_xxzzz, tr_xyyy_xyy, tr_xyyy_xyyyz, tr_xyyy_xyyzz, tr_xyyy_xyz, tr_xyyy_xyzzz, tr_xyyy_xzz, tr_xyyy_xzzzz, tr_xyyy_yyy, tr_xyyy_yyyyz, tr_xyyy_yyyzz, tr_xyyy_yyz, tr_xyyy_yyzzz, tr_xyyy_yzz, tr_xyyy_yzzzz, tr_xyyy_zzz, tr_xyyy_zzzzz, tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyyy_xxxx[i] = 2.0 * tr_xyyyz_xxxx[i] * tbe_0 + 2.0 * tr_xyyy_xxxxz[i] * tke_0;

        tr_0_0_z_xyyy_xxxy[i] = 2.0 * tr_xyyyz_xxxy[i] * tbe_0 + 2.0 * tr_xyyy_xxxyz[i] * tke_0;

        tr_0_0_z_xyyy_xxxz[i] = 2.0 * tr_xyyyz_xxxz[i] * tbe_0 + 2.0 * tr_xyyy_xxxzz[i] * tke_0 - tr_xyyy_xxx[i];

        tr_0_0_z_xyyy_xxyy[i] = 2.0 * tr_xyyyz_xxyy[i] * tbe_0 + 2.0 * tr_xyyy_xxyyz[i] * tke_0;

        tr_0_0_z_xyyy_xxyz[i] = 2.0 * tr_xyyyz_xxyz[i] * tbe_0 + 2.0 * tr_xyyy_xxyzz[i] * tke_0 - tr_xyyy_xxy[i];

        tr_0_0_z_xyyy_xxzz[i] = 2.0 * tr_xyyyz_xxzz[i] * tbe_0 + 2.0 * tr_xyyy_xxzzz[i] * tke_0 - 2.0 * tr_xyyy_xxz[i];

        tr_0_0_z_xyyy_xyyy[i] = 2.0 * tr_xyyyz_xyyy[i] * tbe_0 + 2.0 * tr_xyyy_xyyyz[i] * tke_0;

        tr_0_0_z_xyyy_xyyz[i] = 2.0 * tr_xyyyz_xyyz[i] * tbe_0 + 2.0 * tr_xyyy_xyyzz[i] * tke_0 - tr_xyyy_xyy[i];

        tr_0_0_z_xyyy_xyzz[i] = 2.0 * tr_xyyyz_xyzz[i] * tbe_0 + 2.0 * tr_xyyy_xyzzz[i] * tke_0 - 2.0 * tr_xyyy_xyz[i];

        tr_0_0_z_xyyy_xzzz[i] = 2.0 * tr_xyyyz_xzzz[i] * tbe_0 + 2.0 * tr_xyyy_xzzzz[i] * tke_0 - 3.0 * tr_xyyy_xzz[i];

        tr_0_0_z_xyyy_yyyy[i] = 2.0 * tr_xyyyz_yyyy[i] * tbe_0 + 2.0 * tr_xyyy_yyyyz[i] * tke_0;

        tr_0_0_z_xyyy_yyyz[i] = 2.0 * tr_xyyyz_yyyz[i] * tbe_0 + 2.0 * tr_xyyy_yyyzz[i] * tke_0 - tr_xyyy_yyy[i];

        tr_0_0_z_xyyy_yyzz[i] = 2.0 * tr_xyyyz_yyzz[i] * tbe_0 + 2.0 * tr_xyyy_yyzzz[i] * tke_0 - 2.0 * tr_xyyy_yyz[i];

        tr_0_0_z_xyyy_yzzz[i] = 2.0 * tr_xyyyz_yzzz[i] * tbe_0 + 2.0 * tr_xyyy_yzzzz[i] * tke_0 - 3.0 * tr_xyyy_yzz[i];

        tr_0_0_z_xyyy_zzzz[i] = 2.0 * tr_xyyyz_zzzz[i] * tbe_0 + 2.0 * tr_xyyy_zzzzz[i] * tke_0 - 4.0 * tr_xyyy_zzz[i];
    }

    // Set up 555-570 components of targeted buffer : GG

    auto tr_0_0_z_xyyz_xxxx = pbuffer.data(idx_op_geom_010_gg + 555);

    auto tr_0_0_z_xyyz_xxxy = pbuffer.data(idx_op_geom_010_gg + 556);

    auto tr_0_0_z_xyyz_xxxz = pbuffer.data(idx_op_geom_010_gg + 557);

    auto tr_0_0_z_xyyz_xxyy = pbuffer.data(idx_op_geom_010_gg + 558);

    auto tr_0_0_z_xyyz_xxyz = pbuffer.data(idx_op_geom_010_gg + 559);

    auto tr_0_0_z_xyyz_xxzz = pbuffer.data(idx_op_geom_010_gg + 560);

    auto tr_0_0_z_xyyz_xyyy = pbuffer.data(idx_op_geom_010_gg + 561);

    auto tr_0_0_z_xyyz_xyyz = pbuffer.data(idx_op_geom_010_gg + 562);

    auto tr_0_0_z_xyyz_xyzz = pbuffer.data(idx_op_geom_010_gg + 563);

    auto tr_0_0_z_xyyz_xzzz = pbuffer.data(idx_op_geom_010_gg + 564);

    auto tr_0_0_z_xyyz_yyyy = pbuffer.data(idx_op_geom_010_gg + 565);

    auto tr_0_0_z_xyyz_yyyz = pbuffer.data(idx_op_geom_010_gg + 566);

    auto tr_0_0_z_xyyz_yyzz = pbuffer.data(idx_op_geom_010_gg + 567);

    auto tr_0_0_z_xyyz_yzzz = pbuffer.data(idx_op_geom_010_gg + 568);

    auto tr_0_0_z_xyyz_zzzz = pbuffer.data(idx_op_geom_010_gg + 569);

    #pragma omp simd aligned(tr_0_0_z_xyyz_xxxx, tr_0_0_z_xyyz_xxxy, tr_0_0_z_xyyz_xxxz, tr_0_0_z_xyyz_xxyy, tr_0_0_z_xyyz_xxyz, tr_0_0_z_xyyz_xxzz, tr_0_0_z_xyyz_xyyy, tr_0_0_z_xyyz_xyyz, tr_0_0_z_xyyz_xyzz, tr_0_0_z_xyyz_xzzz, tr_0_0_z_xyyz_yyyy, tr_0_0_z_xyyz_yyyz, tr_0_0_z_xyyz_yyzz, tr_0_0_z_xyyz_yzzz, tr_0_0_z_xyyz_zzzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxz, tr_xyyz_xxxyz, tr_xyyz_xxxzz, tr_xyyz_xxy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xxzzz, tr_xyyz_xyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_xzzzz, tr_xyyz_yyy, tr_xyyz_yyyyz, tr_xyyz_yyyzz, tr_xyyz_yyz, tr_xyyz_yyzzz, tr_xyyz_yzz, tr_xyyz_yzzzz, tr_xyyz_zzz, tr_xyyz_zzzzz, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyyz_xxxx[i] = 2.0 * tr_xyyzz_xxxx[i] * tbe_0 + 2.0 * tr_xyyz_xxxxz[i] * tke_0 - tr_xyy_xxxx[i];

        tr_0_0_z_xyyz_xxxy[i] = 2.0 * tr_xyyzz_xxxy[i] * tbe_0 + 2.0 * tr_xyyz_xxxyz[i] * tke_0 - tr_xyy_xxxy[i];

        tr_0_0_z_xyyz_xxxz[i] = 2.0 * tr_xyyzz_xxxz[i] * tbe_0 + 2.0 * tr_xyyz_xxxzz[i] * tke_0 - tr_xyy_xxxz[i] - tr_xyyz_xxx[i];

        tr_0_0_z_xyyz_xxyy[i] = 2.0 * tr_xyyzz_xxyy[i] * tbe_0 + 2.0 * tr_xyyz_xxyyz[i] * tke_0 - tr_xyy_xxyy[i];

        tr_0_0_z_xyyz_xxyz[i] = 2.0 * tr_xyyzz_xxyz[i] * tbe_0 + 2.0 * tr_xyyz_xxyzz[i] * tke_0 - tr_xyy_xxyz[i] - tr_xyyz_xxy[i];

        tr_0_0_z_xyyz_xxzz[i] = 2.0 * tr_xyyzz_xxzz[i] * tbe_0 + 2.0 * tr_xyyz_xxzzz[i] * tke_0 - tr_xyy_xxzz[i] - 2.0 * tr_xyyz_xxz[i];

        tr_0_0_z_xyyz_xyyy[i] = 2.0 * tr_xyyzz_xyyy[i] * tbe_0 + 2.0 * tr_xyyz_xyyyz[i] * tke_0 - tr_xyy_xyyy[i];

        tr_0_0_z_xyyz_xyyz[i] = 2.0 * tr_xyyzz_xyyz[i] * tbe_0 + 2.0 * tr_xyyz_xyyzz[i] * tke_0 - tr_xyy_xyyz[i] - tr_xyyz_xyy[i];

        tr_0_0_z_xyyz_xyzz[i] = 2.0 * tr_xyyzz_xyzz[i] * tbe_0 + 2.0 * tr_xyyz_xyzzz[i] * tke_0 - tr_xyy_xyzz[i] - 2.0 * tr_xyyz_xyz[i];

        tr_0_0_z_xyyz_xzzz[i] = 2.0 * tr_xyyzz_xzzz[i] * tbe_0 + 2.0 * tr_xyyz_xzzzz[i] * tke_0 - tr_xyy_xzzz[i] - 3.0 * tr_xyyz_xzz[i];

        tr_0_0_z_xyyz_yyyy[i] = 2.0 * tr_xyyzz_yyyy[i] * tbe_0 + 2.0 * tr_xyyz_yyyyz[i] * tke_0 - tr_xyy_yyyy[i];

        tr_0_0_z_xyyz_yyyz[i] = 2.0 * tr_xyyzz_yyyz[i] * tbe_0 + 2.0 * tr_xyyz_yyyzz[i] * tke_0 - tr_xyy_yyyz[i] - tr_xyyz_yyy[i];

        tr_0_0_z_xyyz_yyzz[i] = 2.0 * tr_xyyzz_yyzz[i] * tbe_0 + 2.0 * tr_xyyz_yyzzz[i] * tke_0 - tr_xyy_yyzz[i] - 2.0 * tr_xyyz_yyz[i];

        tr_0_0_z_xyyz_yzzz[i] = 2.0 * tr_xyyzz_yzzz[i] * tbe_0 + 2.0 * tr_xyyz_yzzzz[i] * tke_0 - tr_xyy_yzzz[i] - 3.0 * tr_xyyz_yzz[i];

        tr_0_0_z_xyyz_zzzz[i] = 2.0 * tr_xyyzz_zzzz[i] * tbe_0 + 2.0 * tr_xyyz_zzzzz[i] * tke_0 - tr_xyy_zzzz[i] - 4.0 * tr_xyyz_zzz[i];
    }

    // Set up 570-585 components of targeted buffer : GG

    auto tr_0_0_z_xyzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 570);

    auto tr_0_0_z_xyzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 571);

    auto tr_0_0_z_xyzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 572);

    auto tr_0_0_z_xyzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 573);

    auto tr_0_0_z_xyzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 574);

    auto tr_0_0_z_xyzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 575);

    auto tr_0_0_z_xyzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 576);

    auto tr_0_0_z_xyzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 577);

    auto tr_0_0_z_xyzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 578);

    auto tr_0_0_z_xyzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 579);

    auto tr_0_0_z_xyzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 580);

    auto tr_0_0_z_xyzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 581);

    auto tr_0_0_z_xyzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 582);

    auto tr_0_0_z_xyzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 583);

    auto tr_0_0_z_xyzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 584);

    #pragma omp simd aligned(tr_0_0_z_xyzz_xxxx, tr_0_0_z_xyzz_xxxy, tr_0_0_z_xyzz_xxxz, tr_0_0_z_xyzz_xxyy, tr_0_0_z_xyzz_xxyz, tr_0_0_z_xyzz_xxzz, tr_0_0_z_xyzz_xyyy, tr_0_0_z_xyzz_xyyz, tr_0_0_z_xyzz_xyzz, tr_0_0_z_xyzz_xzzz, tr_0_0_z_xyzz_yyyy, tr_0_0_z_xyzz_yyyz, tr_0_0_z_xyzz_yyzz, tr_0_0_z_xyzz_yzzz, tr_0_0_z_xyzz_zzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxz, tr_xyzz_xxxyz, tr_xyzz_xxxzz, tr_xyzz_xxy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xxzzz, tr_xyzz_xyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_xzzzz, tr_xyzz_yyy, tr_xyzz_yyyyz, tr_xyzz_yyyzz, tr_xyzz_yyz, tr_xyzz_yyzzz, tr_xyzz_yzz, tr_xyzz_yzzzz, tr_xyzz_zzz, tr_xyzz_zzzzz, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyzz_xxxx[i] = 2.0 * tr_xyzzz_xxxx[i] * tbe_0 + 2.0 * tr_xyzz_xxxxz[i] * tke_0 - 2.0 * tr_xyz_xxxx[i];

        tr_0_0_z_xyzz_xxxy[i] = 2.0 * tr_xyzzz_xxxy[i] * tbe_0 + 2.0 * tr_xyzz_xxxyz[i] * tke_0 - 2.0 * tr_xyz_xxxy[i];

        tr_0_0_z_xyzz_xxxz[i] = 2.0 * tr_xyzzz_xxxz[i] * tbe_0 + 2.0 * tr_xyzz_xxxzz[i] * tke_0 - 2.0 * tr_xyz_xxxz[i] - tr_xyzz_xxx[i];

        tr_0_0_z_xyzz_xxyy[i] = 2.0 * tr_xyzzz_xxyy[i] * tbe_0 + 2.0 * tr_xyzz_xxyyz[i] * tke_0 - 2.0 * tr_xyz_xxyy[i];

        tr_0_0_z_xyzz_xxyz[i] = 2.0 * tr_xyzzz_xxyz[i] * tbe_0 + 2.0 * tr_xyzz_xxyzz[i] * tke_0 - 2.0 * tr_xyz_xxyz[i] - tr_xyzz_xxy[i];

        tr_0_0_z_xyzz_xxzz[i] = 2.0 * tr_xyzzz_xxzz[i] * tbe_0 + 2.0 * tr_xyzz_xxzzz[i] * tke_0 - 2.0 * tr_xyz_xxzz[i] - 2.0 * tr_xyzz_xxz[i];

        tr_0_0_z_xyzz_xyyy[i] = 2.0 * tr_xyzzz_xyyy[i] * tbe_0 + 2.0 * tr_xyzz_xyyyz[i] * tke_0 - 2.0 * tr_xyz_xyyy[i];

        tr_0_0_z_xyzz_xyyz[i] = 2.0 * tr_xyzzz_xyyz[i] * tbe_0 + 2.0 * tr_xyzz_xyyzz[i] * tke_0 - 2.0 * tr_xyz_xyyz[i] - tr_xyzz_xyy[i];

        tr_0_0_z_xyzz_xyzz[i] = 2.0 * tr_xyzzz_xyzz[i] * tbe_0 + 2.0 * tr_xyzz_xyzzz[i] * tke_0 - 2.0 * tr_xyz_xyzz[i] - 2.0 * tr_xyzz_xyz[i];

        tr_0_0_z_xyzz_xzzz[i] = 2.0 * tr_xyzzz_xzzz[i] * tbe_0 + 2.0 * tr_xyzz_xzzzz[i] * tke_0 - 2.0 * tr_xyz_xzzz[i] - 3.0 * tr_xyzz_xzz[i];

        tr_0_0_z_xyzz_yyyy[i] = 2.0 * tr_xyzzz_yyyy[i] * tbe_0 + 2.0 * tr_xyzz_yyyyz[i] * tke_0 - 2.0 * tr_xyz_yyyy[i];

        tr_0_0_z_xyzz_yyyz[i] = 2.0 * tr_xyzzz_yyyz[i] * tbe_0 + 2.0 * tr_xyzz_yyyzz[i] * tke_0 - 2.0 * tr_xyz_yyyz[i] - tr_xyzz_yyy[i];

        tr_0_0_z_xyzz_yyzz[i] = 2.0 * tr_xyzzz_yyzz[i] * tbe_0 + 2.0 * tr_xyzz_yyzzz[i] * tke_0 - 2.0 * tr_xyz_yyzz[i] - 2.0 * tr_xyzz_yyz[i];

        tr_0_0_z_xyzz_yzzz[i] = 2.0 * tr_xyzzz_yzzz[i] * tbe_0 + 2.0 * tr_xyzz_yzzzz[i] * tke_0 - 2.0 * tr_xyz_yzzz[i] - 3.0 * tr_xyzz_yzz[i];

        tr_0_0_z_xyzz_zzzz[i] = 2.0 * tr_xyzzz_zzzz[i] * tbe_0 + 2.0 * tr_xyzz_zzzzz[i] * tke_0 - 2.0 * tr_xyz_zzzz[i] - 4.0 * tr_xyzz_zzz[i];
    }

    // Set up 585-600 components of targeted buffer : GG

    auto tr_0_0_z_xzzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 585);

    auto tr_0_0_z_xzzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 586);

    auto tr_0_0_z_xzzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 587);

    auto tr_0_0_z_xzzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 588);

    auto tr_0_0_z_xzzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 589);

    auto tr_0_0_z_xzzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 590);

    auto tr_0_0_z_xzzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 591);

    auto tr_0_0_z_xzzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 592);

    auto tr_0_0_z_xzzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 593);

    auto tr_0_0_z_xzzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 594);

    auto tr_0_0_z_xzzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 595);

    auto tr_0_0_z_xzzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 596);

    auto tr_0_0_z_xzzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 597);

    auto tr_0_0_z_xzzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 598);

    auto tr_0_0_z_xzzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 599);

    #pragma omp simd aligned(tr_0_0_z_xzzz_xxxx, tr_0_0_z_xzzz_xxxy, tr_0_0_z_xzzz_xxxz, tr_0_0_z_xzzz_xxyy, tr_0_0_z_xzzz_xxyz, tr_0_0_z_xzzz_xxzz, tr_0_0_z_xzzz_xyyy, tr_0_0_z_xzzz_xyyz, tr_0_0_z_xzzz_xyzz, tr_0_0_z_xzzz_xzzz, tr_0_0_z_xzzz_yyyy, tr_0_0_z_xzzz_yyyz, tr_0_0_z_xzzz_yyzz, tr_0_0_z_xzzz_yzzz, tr_0_0_z_xzzz_zzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxxxz, tr_xzzz_xxxyz, tr_xzzz_xxxzz, tr_xzzz_xxy, tr_xzzz_xxyyz, tr_xzzz_xxyzz, tr_xzzz_xxz, tr_xzzz_xxzzz, tr_xzzz_xyy, tr_xzzz_xyyyz, tr_xzzz_xyyzz, tr_xzzz_xyz, tr_xzzz_xyzzz, tr_xzzz_xzz, tr_xzzz_xzzzz, tr_xzzz_yyy, tr_xzzz_yyyyz, tr_xzzz_yyyzz, tr_xzzz_yyz, tr_xzzz_yyzzz, tr_xzzz_yzz, tr_xzzz_yzzzz, tr_xzzz_zzz, tr_xzzz_zzzzz, tr_xzzzz_xxxx, tr_xzzzz_xxxy, tr_xzzzz_xxxz, tr_xzzzz_xxyy, tr_xzzzz_xxyz, tr_xzzzz_xxzz, tr_xzzzz_xyyy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xzzz, tr_xzzzz_yyyy, tr_xzzzz_yyyz, tr_xzzzz_yyzz, tr_xzzzz_yzzz, tr_xzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xzzz_xxxx[i] = 2.0 * tr_xzzzz_xxxx[i] * tbe_0 + 2.0 * tr_xzzz_xxxxz[i] * tke_0 - 3.0 * tr_xzz_xxxx[i];

        tr_0_0_z_xzzz_xxxy[i] = 2.0 * tr_xzzzz_xxxy[i] * tbe_0 + 2.0 * tr_xzzz_xxxyz[i] * tke_0 - 3.0 * tr_xzz_xxxy[i];

        tr_0_0_z_xzzz_xxxz[i] = 2.0 * tr_xzzzz_xxxz[i] * tbe_0 + 2.0 * tr_xzzz_xxxzz[i] * tke_0 - 3.0 * tr_xzz_xxxz[i] - tr_xzzz_xxx[i];

        tr_0_0_z_xzzz_xxyy[i] = 2.0 * tr_xzzzz_xxyy[i] * tbe_0 + 2.0 * tr_xzzz_xxyyz[i] * tke_0 - 3.0 * tr_xzz_xxyy[i];

        tr_0_0_z_xzzz_xxyz[i] = 2.0 * tr_xzzzz_xxyz[i] * tbe_0 + 2.0 * tr_xzzz_xxyzz[i] * tke_0 - 3.0 * tr_xzz_xxyz[i] - tr_xzzz_xxy[i];

        tr_0_0_z_xzzz_xxzz[i] = 2.0 * tr_xzzzz_xxzz[i] * tbe_0 + 2.0 * tr_xzzz_xxzzz[i] * tke_0 - 3.0 * tr_xzz_xxzz[i] - 2.0 * tr_xzzz_xxz[i];

        tr_0_0_z_xzzz_xyyy[i] = 2.0 * tr_xzzzz_xyyy[i] * tbe_0 + 2.0 * tr_xzzz_xyyyz[i] * tke_0 - 3.0 * tr_xzz_xyyy[i];

        tr_0_0_z_xzzz_xyyz[i] = 2.0 * tr_xzzzz_xyyz[i] * tbe_0 + 2.0 * tr_xzzz_xyyzz[i] * tke_0 - 3.0 * tr_xzz_xyyz[i] - tr_xzzz_xyy[i];

        tr_0_0_z_xzzz_xyzz[i] = 2.0 * tr_xzzzz_xyzz[i] * tbe_0 + 2.0 * tr_xzzz_xyzzz[i] * tke_0 - 3.0 * tr_xzz_xyzz[i] - 2.0 * tr_xzzz_xyz[i];

        tr_0_0_z_xzzz_xzzz[i] = 2.0 * tr_xzzzz_xzzz[i] * tbe_0 + 2.0 * tr_xzzz_xzzzz[i] * tke_0 - 3.0 * tr_xzz_xzzz[i] - 3.0 * tr_xzzz_xzz[i];

        tr_0_0_z_xzzz_yyyy[i] = 2.0 * tr_xzzzz_yyyy[i] * tbe_0 + 2.0 * tr_xzzz_yyyyz[i] * tke_0 - 3.0 * tr_xzz_yyyy[i];

        tr_0_0_z_xzzz_yyyz[i] = 2.0 * tr_xzzzz_yyyz[i] * tbe_0 + 2.0 * tr_xzzz_yyyzz[i] * tke_0 - 3.0 * tr_xzz_yyyz[i] - tr_xzzz_yyy[i];

        tr_0_0_z_xzzz_yyzz[i] = 2.0 * tr_xzzzz_yyzz[i] * tbe_0 + 2.0 * tr_xzzz_yyzzz[i] * tke_0 - 3.0 * tr_xzz_yyzz[i] - 2.0 * tr_xzzz_yyz[i];

        tr_0_0_z_xzzz_yzzz[i] = 2.0 * tr_xzzzz_yzzz[i] * tbe_0 + 2.0 * tr_xzzz_yzzzz[i] * tke_0 - 3.0 * tr_xzz_yzzz[i] - 3.0 * tr_xzzz_yzz[i];

        tr_0_0_z_xzzz_zzzz[i] = 2.0 * tr_xzzzz_zzzz[i] * tbe_0 + 2.0 * tr_xzzz_zzzzz[i] * tke_0 - 3.0 * tr_xzz_zzzz[i] - 4.0 * tr_xzzz_zzz[i];
    }

    // Set up 600-615 components of targeted buffer : GG

    auto tr_0_0_z_yyyy_xxxx = pbuffer.data(idx_op_geom_010_gg + 600);

    auto tr_0_0_z_yyyy_xxxy = pbuffer.data(idx_op_geom_010_gg + 601);

    auto tr_0_0_z_yyyy_xxxz = pbuffer.data(idx_op_geom_010_gg + 602);

    auto tr_0_0_z_yyyy_xxyy = pbuffer.data(idx_op_geom_010_gg + 603);

    auto tr_0_0_z_yyyy_xxyz = pbuffer.data(idx_op_geom_010_gg + 604);

    auto tr_0_0_z_yyyy_xxzz = pbuffer.data(idx_op_geom_010_gg + 605);

    auto tr_0_0_z_yyyy_xyyy = pbuffer.data(idx_op_geom_010_gg + 606);

    auto tr_0_0_z_yyyy_xyyz = pbuffer.data(idx_op_geom_010_gg + 607);

    auto tr_0_0_z_yyyy_xyzz = pbuffer.data(idx_op_geom_010_gg + 608);

    auto tr_0_0_z_yyyy_xzzz = pbuffer.data(idx_op_geom_010_gg + 609);

    auto tr_0_0_z_yyyy_yyyy = pbuffer.data(idx_op_geom_010_gg + 610);

    auto tr_0_0_z_yyyy_yyyz = pbuffer.data(idx_op_geom_010_gg + 611);

    auto tr_0_0_z_yyyy_yyzz = pbuffer.data(idx_op_geom_010_gg + 612);

    auto tr_0_0_z_yyyy_yzzz = pbuffer.data(idx_op_geom_010_gg + 613);

    auto tr_0_0_z_yyyy_zzzz = pbuffer.data(idx_op_geom_010_gg + 614);

    #pragma omp simd aligned(tr_0_0_z_yyyy_xxxx, tr_0_0_z_yyyy_xxxy, tr_0_0_z_yyyy_xxxz, tr_0_0_z_yyyy_xxyy, tr_0_0_z_yyyy_xxyz, tr_0_0_z_yyyy_xxzz, tr_0_0_z_yyyy_xyyy, tr_0_0_z_yyyy_xyyz, tr_0_0_z_yyyy_xyzz, tr_0_0_z_yyyy_xzzz, tr_0_0_z_yyyy_yyyy, tr_0_0_z_yyyy_yyyz, tr_0_0_z_yyyy_yyzz, tr_0_0_z_yyyy_yzzz, tr_0_0_z_yyyy_zzzz, tr_yyyy_xxx, tr_yyyy_xxxxz, tr_yyyy_xxxyz, tr_yyyy_xxxzz, tr_yyyy_xxy, tr_yyyy_xxyyz, tr_yyyy_xxyzz, tr_yyyy_xxz, tr_yyyy_xxzzz, tr_yyyy_xyy, tr_yyyy_xyyyz, tr_yyyy_xyyzz, tr_yyyy_xyz, tr_yyyy_xyzzz, tr_yyyy_xzz, tr_yyyy_xzzzz, tr_yyyy_yyy, tr_yyyy_yyyyz, tr_yyyy_yyyzz, tr_yyyy_yyz, tr_yyyy_yyzzz, tr_yyyy_yzz, tr_yyyy_yzzzz, tr_yyyy_zzz, tr_yyyy_zzzzz, tr_yyyyz_xxxx, tr_yyyyz_xxxy, tr_yyyyz_xxxz, tr_yyyyz_xxyy, tr_yyyyz_xxyz, tr_yyyyz_xxzz, tr_yyyyz_xyyy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xzzz, tr_yyyyz_yyyy, tr_yyyyz_yyyz, tr_yyyyz_yyzz, tr_yyyyz_yzzz, tr_yyyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyyy_xxxx[i] = 2.0 * tr_yyyyz_xxxx[i] * tbe_0 + 2.0 * tr_yyyy_xxxxz[i] * tke_0;

        tr_0_0_z_yyyy_xxxy[i] = 2.0 * tr_yyyyz_xxxy[i] * tbe_0 + 2.0 * tr_yyyy_xxxyz[i] * tke_0;

        tr_0_0_z_yyyy_xxxz[i] = 2.0 * tr_yyyyz_xxxz[i] * tbe_0 + 2.0 * tr_yyyy_xxxzz[i] * tke_0 - tr_yyyy_xxx[i];

        tr_0_0_z_yyyy_xxyy[i] = 2.0 * tr_yyyyz_xxyy[i] * tbe_0 + 2.0 * tr_yyyy_xxyyz[i] * tke_0;

        tr_0_0_z_yyyy_xxyz[i] = 2.0 * tr_yyyyz_xxyz[i] * tbe_0 + 2.0 * tr_yyyy_xxyzz[i] * tke_0 - tr_yyyy_xxy[i];

        tr_0_0_z_yyyy_xxzz[i] = 2.0 * tr_yyyyz_xxzz[i] * tbe_0 + 2.0 * tr_yyyy_xxzzz[i] * tke_0 - 2.0 * tr_yyyy_xxz[i];

        tr_0_0_z_yyyy_xyyy[i] = 2.0 * tr_yyyyz_xyyy[i] * tbe_0 + 2.0 * tr_yyyy_xyyyz[i] * tke_0;

        tr_0_0_z_yyyy_xyyz[i] = 2.0 * tr_yyyyz_xyyz[i] * tbe_0 + 2.0 * tr_yyyy_xyyzz[i] * tke_0 - tr_yyyy_xyy[i];

        tr_0_0_z_yyyy_xyzz[i] = 2.0 * tr_yyyyz_xyzz[i] * tbe_0 + 2.0 * tr_yyyy_xyzzz[i] * tke_0 - 2.0 * tr_yyyy_xyz[i];

        tr_0_0_z_yyyy_xzzz[i] = 2.0 * tr_yyyyz_xzzz[i] * tbe_0 + 2.0 * tr_yyyy_xzzzz[i] * tke_0 - 3.0 * tr_yyyy_xzz[i];

        tr_0_0_z_yyyy_yyyy[i] = 2.0 * tr_yyyyz_yyyy[i] * tbe_0 + 2.0 * tr_yyyy_yyyyz[i] * tke_0;

        tr_0_0_z_yyyy_yyyz[i] = 2.0 * tr_yyyyz_yyyz[i] * tbe_0 + 2.0 * tr_yyyy_yyyzz[i] * tke_0 - tr_yyyy_yyy[i];

        tr_0_0_z_yyyy_yyzz[i] = 2.0 * tr_yyyyz_yyzz[i] * tbe_0 + 2.0 * tr_yyyy_yyzzz[i] * tke_0 - 2.0 * tr_yyyy_yyz[i];

        tr_0_0_z_yyyy_yzzz[i] = 2.0 * tr_yyyyz_yzzz[i] * tbe_0 + 2.0 * tr_yyyy_yzzzz[i] * tke_0 - 3.0 * tr_yyyy_yzz[i];

        tr_0_0_z_yyyy_zzzz[i] = 2.0 * tr_yyyyz_zzzz[i] * tbe_0 + 2.0 * tr_yyyy_zzzzz[i] * tke_0 - 4.0 * tr_yyyy_zzz[i];
    }

    // Set up 615-630 components of targeted buffer : GG

    auto tr_0_0_z_yyyz_xxxx = pbuffer.data(idx_op_geom_010_gg + 615);

    auto tr_0_0_z_yyyz_xxxy = pbuffer.data(idx_op_geom_010_gg + 616);

    auto tr_0_0_z_yyyz_xxxz = pbuffer.data(idx_op_geom_010_gg + 617);

    auto tr_0_0_z_yyyz_xxyy = pbuffer.data(idx_op_geom_010_gg + 618);

    auto tr_0_0_z_yyyz_xxyz = pbuffer.data(idx_op_geom_010_gg + 619);

    auto tr_0_0_z_yyyz_xxzz = pbuffer.data(idx_op_geom_010_gg + 620);

    auto tr_0_0_z_yyyz_xyyy = pbuffer.data(idx_op_geom_010_gg + 621);

    auto tr_0_0_z_yyyz_xyyz = pbuffer.data(idx_op_geom_010_gg + 622);

    auto tr_0_0_z_yyyz_xyzz = pbuffer.data(idx_op_geom_010_gg + 623);

    auto tr_0_0_z_yyyz_xzzz = pbuffer.data(idx_op_geom_010_gg + 624);

    auto tr_0_0_z_yyyz_yyyy = pbuffer.data(idx_op_geom_010_gg + 625);

    auto tr_0_0_z_yyyz_yyyz = pbuffer.data(idx_op_geom_010_gg + 626);

    auto tr_0_0_z_yyyz_yyzz = pbuffer.data(idx_op_geom_010_gg + 627);

    auto tr_0_0_z_yyyz_yzzz = pbuffer.data(idx_op_geom_010_gg + 628);

    auto tr_0_0_z_yyyz_zzzz = pbuffer.data(idx_op_geom_010_gg + 629);

    #pragma omp simd aligned(tr_0_0_z_yyyz_xxxx, tr_0_0_z_yyyz_xxxy, tr_0_0_z_yyyz_xxxz, tr_0_0_z_yyyz_xxyy, tr_0_0_z_yyyz_xxyz, tr_0_0_z_yyyz_xxzz, tr_0_0_z_yyyz_xyyy, tr_0_0_z_yyyz_xyyz, tr_0_0_z_yyyz_xyzz, tr_0_0_z_yyyz_xzzz, tr_0_0_z_yyyz_yyyy, tr_0_0_z_yyyz_yyyz, tr_0_0_z_yyyz_yyzz, tr_0_0_z_yyyz_yzzz, tr_0_0_z_yyyz_zzzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, tr_yyyz_xxx, tr_yyyz_xxxxz, tr_yyyz_xxxyz, tr_yyyz_xxxzz, tr_yyyz_xxy, tr_yyyz_xxyyz, tr_yyyz_xxyzz, tr_yyyz_xxz, tr_yyyz_xxzzz, tr_yyyz_xyy, tr_yyyz_xyyyz, tr_yyyz_xyyzz, tr_yyyz_xyz, tr_yyyz_xyzzz, tr_yyyz_xzz, tr_yyyz_xzzzz, tr_yyyz_yyy, tr_yyyz_yyyyz, tr_yyyz_yyyzz, tr_yyyz_yyz, tr_yyyz_yyzzz, tr_yyyz_yzz, tr_yyyz_yzzzz, tr_yyyz_zzz, tr_yyyz_zzzzz, tr_yyyzz_xxxx, tr_yyyzz_xxxy, tr_yyyzz_xxxz, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xzzz, tr_yyyzz_yyyy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yzzz, tr_yyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyyz_xxxx[i] = 2.0 * tr_yyyzz_xxxx[i] * tbe_0 + 2.0 * tr_yyyz_xxxxz[i] * tke_0 - tr_yyy_xxxx[i];

        tr_0_0_z_yyyz_xxxy[i] = 2.0 * tr_yyyzz_xxxy[i] * tbe_0 + 2.0 * tr_yyyz_xxxyz[i] * tke_0 - tr_yyy_xxxy[i];

        tr_0_0_z_yyyz_xxxz[i] = 2.0 * tr_yyyzz_xxxz[i] * tbe_0 + 2.0 * tr_yyyz_xxxzz[i] * tke_0 - tr_yyy_xxxz[i] - tr_yyyz_xxx[i];

        tr_0_0_z_yyyz_xxyy[i] = 2.0 * tr_yyyzz_xxyy[i] * tbe_0 + 2.0 * tr_yyyz_xxyyz[i] * tke_0 - tr_yyy_xxyy[i];

        tr_0_0_z_yyyz_xxyz[i] = 2.0 * tr_yyyzz_xxyz[i] * tbe_0 + 2.0 * tr_yyyz_xxyzz[i] * tke_0 - tr_yyy_xxyz[i] - tr_yyyz_xxy[i];

        tr_0_0_z_yyyz_xxzz[i] = 2.0 * tr_yyyzz_xxzz[i] * tbe_0 + 2.0 * tr_yyyz_xxzzz[i] * tke_0 - tr_yyy_xxzz[i] - 2.0 * tr_yyyz_xxz[i];

        tr_0_0_z_yyyz_xyyy[i] = 2.0 * tr_yyyzz_xyyy[i] * tbe_0 + 2.0 * tr_yyyz_xyyyz[i] * tke_0 - tr_yyy_xyyy[i];

        tr_0_0_z_yyyz_xyyz[i] = 2.0 * tr_yyyzz_xyyz[i] * tbe_0 + 2.0 * tr_yyyz_xyyzz[i] * tke_0 - tr_yyy_xyyz[i] - tr_yyyz_xyy[i];

        tr_0_0_z_yyyz_xyzz[i] = 2.0 * tr_yyyzz_xyzz[i] * tbe_0 + 2.0 * tr_yyyz_xyzzz[i] * tke_0 - tr_yyy_xyzz[i] - 2.0 * tr_yyyz_xyz[i];

        tr_0_0_z_yyyz_xzzz[i] = 2.0 * tr_yyyzz_xzzz[i] * tbe_0 + 2.0 * tr_yyyz_xzzzz[i] * tke_0 - tr_yyy_xzzz[i] - 3.0 * tr_yyyz_xzz[i];

        tr_0_0_z_yyyz_yyyy[i] = 2.0 * tr_yyyzz_yyyy[i] * tbe_0 + 2.0 * tr_yyyz_yyyyz[i] * tke_0 - tr_yyy_yyyy[i];

        tr_0_0_z_yyyz_yyyz[i] = 2.0 * tr_yyyzz_yyyz[i] * tbe_0 + 2.0 * tr_yyyz_yyyzz[i] * tke_0 - tr_yyy_yyyz[i] - tr_yyyz_yyy[i];

        tr_0_0_z_yyyz_yyzz[i] = 2.0 * tr_yyyzz_yyzz[i] * tbe_0 + 2.0 * tr_yyyz_yyzzz[i] * tke_0 - tr_yyy_yyzz[i] - 2.0 * tr_yyyz_yyz[i];

        tr_0_0_z_yyyz_yzzz[i] = 2.0 * tr_yyyzz_yzzz[i] * tbe_0 + 2.0 * tr_yyyz_yzzzz[i] * tke_0 - tr_yyy_yzzz[i] - 3.0 * tr_yyyz_yzz[i];

        tr_0_0_z_yyyz_zzzz[i] = 2.0 * tr_yyyzz_zzzz[i] * tbe_0 + 2.0 * tr_yyyz_zzzzz[i] * tke_0 - tr_yyy_zzzz[i] - 4.0 * tr_yyyz_zzz[i];
    }

    // Set up 630-645 components of targeted buffer : GG

    auto tr_0_0_z_yyzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 630);

    auto tr_0_0_z_yyzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 631);

    auto tr_0_0_z_yyzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 632);

    auto tr_0_0_z_yyzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 633);

    auto tr_0_0_z_yyzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 634);

    auto tr_0_0_z_yyzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 635);

    auto tr_0_0_z_yyzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 636);

    auto tr_0_0_z_yyzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 637);

    auto tr_0_0_z_yyzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 638);

    auto tr_0_0_z_yyzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 639);

    auto tr_0_0_z_yyzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 640);

    auto tr_0_0_z_yyzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 641);

    auto tr_0_0_z_yyzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 642);

    auto tr_0_0_z_yyzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 643);

    auto tr_0_0_z_yyzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 644);

    #pragma omp simd aligned(tr_0_0_z_yyzz_xxxx, tr_0_0_z_yyzz_xxxy, tr_0_0_z_yyzz_xxxz, tr_0_0_z_yyzz_xxyy, tr_0_0_z_yyzz_xxyz, tr_0_0_z_yyzz_xxzz, tr_0_0_z_yyzz_xyyy, tr_0_0_z_yyzz_xyyz, tr_0_0_z_yyzz_xyzz, tr_0_0_z_yyzz_xzzz, tr_0_0_z_yyzz_yyyy, tr_0_0_z_yyzz_yyyz, tr_0_0_z_yyzz_yyzz, tr_0_0_z_yyzz_yzzz, tr_0_0_z_yyzz_zzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yyzz_xxx, tr_yyzz_xxxxz, tr_yyzz_xxxyz, tr_yyzz_xxxzz, tr_yyzz_xxy, tr_yyzz_xxyyz, tr_yyzz_xxyzz, tr_yyzz_xxz, tr_yyzz_xxzzz, tr_yyzz_xyy, tr_yyzz_xyyyz, tr_yyzz_xyyzz, tr_yyzz_xyz, tr_yyzz_xyzzz, tr_yyzz_xzz, tr_yyzz_xzzzz, tr_yyzz_yyy, tr_yyzz_yyyyz, tr_yyzz_yyyzz, tr_yyzz_yyz, tr_yyzz_yyzzz, tr_yyzz_yzz, tr_yyzz_yzzzz, tr_yyzz_zzz, tr_yyzz_zzzzz, tr_yyzzz_xxxx, tr_yyzzz_xxxy, tr_yyzzz_xxxz, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xzzz, tr_yyzzz_yyyy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yzzz, tr_yyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyzz_xxxx[i] = 2.0 * tr_yyzzz_xxxx[i] * tbe_0 + 2.0 * tr_yyzz_xxxxz[i] * tke_0 - 2.0 * tr_yyz_xxxx[i];

        tr_0_0_z_yyzz_xxxy[i] = 2.0 * tr_yyzzz_xxxy[i] * tbe_0 + 2.0 * tr_yyzz_xxxyz[i] * tke_0 - 2.0 * tr_yyz_xxxy[i];

        tr_0_0_z_yyzz_xxxz[i] = 2.0 * tr_yyzzz_xxxz[i] * tbe_0 + 2.0 * tr_yyzz_xxxzz[i] * tke_0 - 2.0 * tr_yyz_xxxz[i] - tr_yyzz_xxx[i];

        tr_0_0_z_yyzz_xxyy[i] = 2.0 * tr_yyzzz_xxyy[i] * tbe_0 + 2.0 * tr_yyzz_xxyyz[i] * tke_0 - 2.0 * tr_yyz_xxyy[i];

        tr_0_0_z_yyzz_xxyz[i] = 2.0 * tr_yyzzz_xxyz[i] * tbe_0 + 2.0 * tr_yyzz_xxyzz[i] * tke_0 - 2.0 * tr_yyz_xxyz[i] - tr_yyzz_xxy[i];

        tr_0_0_z_yyzz_xxzz[i] = 2.0 * tr_yyzzz_xxzz[i] * tbe_0 + 2.0 * tr_yyzz_xxzzz[i] * tke_0 - 2.0 * tr_yyz_xxzz[i] - 2.0 * tr_yyzz_xxz[i];

        tr_0_0_z_yyzz_xyyy[i] = 2.0 * tr_yyzzz_xyyy[i] * tbe_0 + 2.0 * tr_yyzz_xyyyz[i] * tke_0 - 2.0 * tr_yyz_xyyy[i];

        tr_0_0_z_yyzz_xyyz[i] = 2.0 * tr_yyzzz_xyyz[i] * tbe_0 + 2.0 * tr_yyzz_xyyzz[i] * tke_0 - 2.0 * tr_yyz_xyyz[i] - tr_yyzz_xyy[i];

        tr_0_0_z_yyzz_xyzz[i] = 2.0 * tr_yyzzz_xyzz[i] * tbe_0 + 2.0 * tr_yyzz_xyzzz[i] * tke_0 - 2.0 * tr_yyz_xyzz[i] - 2.0 * tr_yyzz_xyz[i];

        tr_0_0_z_yyzz_xzzz[i] = 2.0 * tr_yyzzz_xzzz[i] * tbe_0 + 2.0 * tr_yyzz_xzzzz[i] * tke_0 - 2.0 * tr_yyz_xzzz[i] - 3.0 * tr_yyzz_xzz[i];

        tr_0_0_z_yyzz_yyyy[i] = 2.0 * tr_yyzzz_yyyy[i] * tbe_0 + 2.0 * tr_yyzz_yyyyz[i] * tke_0 - 2.0 * tr_yyz_yyyy[i];

        tr_0_0_z_yyzz_yyyz[i] = 2.0 * tr_yyzzz_yyyz[i] * tbe_0 + 2.0 * tr_yyzz_yyyzz[i] * tke_0 - 2.0 * tr_yyz_yyyz[i] - tr_yyzz_yyy[i];

        tr_0_0_z_yyzz_yyzz[i] = 2.0 * tr_yyzzz_yyzz[i] * tbe_0 + 2.0 * tr_yyzz_yyzzz[i] * tke_0 - 2.0 * tr_yyz_yyzz[i] - 2.0 * tr_yyzz_yyz[i];

        tr_0_0_z_yyzz_yzzz[i] = 2.0 * tr_yyzzz_yzzz[i] * tbe_0 + 2.0 * tr_yyzz_yzzzz[i] * tke_0 - 2.0 * tr_yyz_yzzz[i] - 3.0 * tr_yyzz_yzz[i];

        tr_0_0_z_yyzz_zzzz[i] = 2.0 * tr_yyzzz_zzzz[i] * tbe_0 + 2.0 * tr_yyzz_zzzzz[i] * tke_0 - 2.0 * tr_yyz_zzzz[i] - 4.0 * tr_yyzz_zzz[i];
    }

    // Set up 645-660 components of targeted buffer : GG

    auto tr_0_0_z_yzzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 645);

    auto tr_0_0_z_yzzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 646);

    auto tr_0_0_z_yzzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 647);

    auto tr_0_0_z_yzzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 648);

    auto tr_0_0_z_yzzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 649);

    auto tr_0_0_z_yzzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 650);

    auto tr_0_0_z_yzzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 651);

    auto tr_0_0_z_yzzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 652);

    auto tr_0_0_z_yzzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 653);

    auto tr_0_0_z_yzzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 654);

    auto tr_0_0_z_yzzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 655);

    auto tr_0_0_z_yzzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 656);

    auto tr_0_0_z_yzzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 657);

    auto tr_0_0_z_yzzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 658);

    auto tr_0_0_z_yzzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 659);

    #pragma omp simd aligned(tr_0_0_z_yzzz_xxxx, tr_0_0_z_yzzz_xxxy, tr_0_0_z_yzzz_xxxz, tr_0_0_z_yzzz_xxyy, tr_0_0_z_yzzz_xxyz, tr_0_0_z_yzzz_xxzz, tr_0_0_z_yzzz_xyyy, tr_0_0_z_yzzz_xyyz, tr_0_0_z_yzzz_xyzz, tr_0_0_z_yzzz_xzzz, tr_0_0_z_yzzz_yyyy, tr_0_0_z_yzzz_yyyz, tr_0_0_z_yzzz_yyzz, tr_0_0_z_yzzz_yzzz, tr_0_0_z_yzzz_zzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxxxz, tr_yzzz_xxxyz, tr_yzzz_xxxzz, tr_yzzz_xxy, tr_yzzz_xxyyz, tr_yzzz_xxyzz, tr_yzzz_xxz, tr_yzzz_xxzzz, tr_yzzz_xyy, tr_yzzz_xyyyz, tr_yzzz_xyyzz, tr_yzzz_xyz, tr_yzzz_xyzzz, tr_yzzz_xzz, tr_yzzz_xzzzz, tr_yzzz_yyy, tr_yzzz_yyyyz, tr_yzzz_yyyzz, tr_yzzz_yyz, tr_yzzz_yyzzz, tr_yzzz_yzz, tr_yzzz_yzzzz, tr_yzzz_zzz, tr_yzzz_zzzzz, tr_yzzzz_xxxx, tr_yzzzz_xxxy, tr_yzzzz_xxxz, tr_yzzzz_xxyy, tr_yzzzz_xxyz, tr_yzzzz_xxzz, tr_yzzzz_xyyy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xzzz, tr_yzzzz_yyyy, tr_yzzzz_yyyz, tr_yzzzz_yyzz, tr_yzzzz_yzzz, tr_yzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yzzz_xxxx[i] = 2.0 * tr_yzzzz_xxxx[i] * tbe_0 + 2.0 * tr_yzzz_xxxxz[i] * tke_0 - 3.0 * tr_yzz_xxxx[i];

        tr_0_0_z_yzzz_xxxy[i] = 2.0 * tr_yzzzz_xxxy[i] * tbe_0 + 2.0 * tr_yzzz_xxxyz[i] * tke_0 - 3.0 * tr_yzz_xxxy[i];

        tr_0_0_z_yzzz_xxxz[i] = 2.0 * tr_yzzzz_xxxz[i] * tbe_0 + 2.0 * tr_yzzz_xxxzz[i] * tke_0 - 3.0 * tr_yzz_xxxz[i] - tr_yzzz_xxx[i];

        tr_0_0_z_yzzz_xxyy[i] = 2.0 * tr_yzzzz_xxyy[i] * tbe_0 + 2.0 * tr_yzzz_xxyyz[i] * tke_0 - 3.0 * tr_yzz_xxyy[i];

        tr_0_0_z_yzzz_xxyz[i] = 2.0 * tr_yzzzz_xxyz[i] * tbe_0 + 2.0 * tr_yzzz_xxyzz[i] * tke_0 - 3.0 * tr_yzz_xxyz[i] - tr_yzzz_xxy[i];

        tr_0_0_z_yzzz_xxzz[i] = 2.0 * tr_yzzzz_xxzz[i] * tbe_0 + 2.0 * tr_yzzz_xxzzz[i] * tke_0 - 3.0 * tr_yzz_xxzz[i] - 2.0 * tr_yzzz_xxz[i];

        tr_0_0_z_yzzz_xyyy[i] = 2.0 * tr_yzzzz_xyyy[i] * tbe_0 + 2.0 * tr_yzzz_xyyyz[i] * tke_0 - 3.0 * tr_yzz_xyyy[i];

        tr_0_0_z_yzzz_xyyz[i] = 2.0 * tr_yzzzz_xyyz[i] * tbe_0 + 2.0 * tr_yzzz_xyyzz[i] * tke_0 - 3.0 * tr_yzz_xyyz[i] - tr_yzzz_xyy[i];

        tr_0_0_z_yzzz_xyzz[i] = 2.0 * tr_yzzzz_xyzz[i] * tbe_0 + 2.0 * tr_yzzz_xyzzz[i] * tke_0 - 3.0 * tr_yzz_xyzz[i] - 2.0 * tr_yzzz_xyz[i];

        tr_0_0_z_yzzz_xzzz[i] = 2.0 * tr_yzzzz_xzzz[i] * tbe_0 + 2.0 * tr_yzzz_xzzzz[i] * tke_0 - 3.0 * tr_yzz_xzzz[i] - 3.0 * tr_yzzz_xzz[i];

        tr_0_0_z_yzzz_yyyy[i] = 2.0 * tr_yzzzz_yyyy[i] * tbe_0 + 2.0 * tr_yzzz_yyyyz[i] * tke_0 - 3.0 * tr_yzz_yyyy[i];

        tr_0_0_z_yzzz_yyyz[i] = 2.0 * tr_yzzzz_yyyz[i] * tbe_0 + 2.0 * tr_yzzz_yyyzz[i] * tke_0 - 3.0 * tr_yzz_yyyz[i] - tr_yzzz_yyy[i];

        tr_0_0_z_yzzz_yyzz[i] = 2.0 * tr_yzzzz_yyzz[i] * tbe_0 + 2.0 * tr_yzzz_yyzzz[i] * tke_0 - 3.0 * tr_yzz_yyzz[i] - 2.0 * tr_yzzz_yyz[i];

        tr_0_0_z_yzzz_yzzz[i] = 2.0 * tr_yzzzz_yzzz[i] * tbe_0 + 2.0 * tr_yzzz_yzzzz[i] * tke_0 - 3.0 * tr_yzz_yzzz[i] - 3.0 * tr_yzzz_yzz[i];

        tr_0_0_z_yzzz_zzzz[i] = 2.0 * tr_yzzzz_zzzz[i] * tbe_0 + 2.0 * tr_yzzz_zzzzz[i] * tke_0 - 3.0 * tr_yzz_zzzz[i] - 4.0 * tr_yzzz_zzz[i];
    }

    // Set up 660-675 components of targeted buffer : GG

    auto tr_0_0_z_zzzz_xxxx = pbuffer.data(idx_op_geom_010_gg + 660);

    auto tr_0_0_z_zzzz_xxxy = pbuffer.data(idx_op_geom_010_gg + 661);

    auto tr_0_0_z_zzzz_xxxz = pbuffer.data(idx_op_geom_010_gg + 662);

    auto tr_0_0_z_zzzz_xxyy = pbuffer.data(idx_op_geom_010_gg + 663);

    auto tr_0_0_z_zzzz_xxyz = pbuffer.data(idx_op_geom_010_gg + 664);

    auto tr_0_0_z_zzzz_xxzz = pbuffer.data(idx_op_geom_010_gg + 665);

    auto tr_0_0_z_zzzz_xyyy = pbuffer.data(idx_op_geom_010_gg + 666);

    auto tr_0_0_z_zzzz_xyyz = pbuffer.data(idx_op_geom_010_gg + 667);

    auto tr_0_0_z_zzzz_xyzz = pbuffer.data(idx_op_geom_010_gg + 668);

    auto tr_0_0_z_zzzz_xzzz = pbuffer.data(idx_op_geom_010_gg + 669);

    auto tr_0_0_z_zzzz_yyyy = pbuffer.data(idx_op_geom_010_gg + 670);

    auto tr_0_0_z_zzzz_yyyz = pbuffer.data(idx_op_geom_010_gg + 671);

    auto tr_0_0_z_zzzz_yyzz = pbuffer.data(idx_op_geom_010_gg + 672);

    auto tr_0_0_z_zzzz_yzzz = pbuffer.data(idx_op_geom_010_gg + 673);

    auto tr_0_0_z_zzzz_zzzz = pbuffer.data(idx_op_geom_010_gg + 674);

    #pragma omp simd aligned(tr_0_0_z_zzzz_xxxx, tr_0_0_z_zzzz_xxxy, tr_0_0_z_zzzz_xxxz, tr_0_0_z_zzzz_xxyy, tr_0_0_z_zzzz_xxyz, tr_0_0_z_zzzz_xxzz, tr_0_0_z_zzzz_xyyy, tr_0_0_z_zzzz_xyyz, tr_0_0_z_zzzz_xyzz, tr_0_0_z_zzzz_xzzz, tr_0_0_z_zzzz_yyyy, tr_0_0_z_zzzz_yyyz, tr_0_0_z_zzzz_yyzz, tr_0_0_z_zzzz_yzzz, tr_0_0_z_zzzz_zzzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, tr_zzzz_xxx, tr_zzzz_xxxxz, tr_zzzz_xxxyz, tr_zzzz_xxxzz, tr_zzzz_xxy, tr_zzzz_xxyyz, tr_zzzz_xxyzz, tr_zzzz_xxz, tr_zzzz_xxzzz, tr_zzzz_xyy, tr_zzzz_xyyyz, tr_zzzz_xyyzz, tr_zzzz_xyz, tr_zzzz_xyzzz, tr_zzzz_xzz, tr_zzzz_xzzzz, tr_zzzz_yyy, tr_zzzz_yyyyz, tr_zzzz_yyyzz, tr_zzzz_yyz, tr_zzzz_yyzzz, tr_zzzz_yzz, tr_zzzz_yzzzz, tr_zzzz_zzz, tr_zzzz_zzzzz, tr_zzzzz_xxxx, tr_zzzzz_xxxy, tr_zzzzz_xxxz, tr_zzzzz_xxyy, tr_zzzzz_xxyz, tr_zzzzz_xxzz, tr_zzzzz_xyyy, tr_zzzzz_xyyz, tr_zzzzz_xyzz, tr_zzzzz_xzzz, tr_zzzzz_yyyy, tr_zzzzz_yyyz, tr_zzzzz_yyzz, tr_zzzzz_yzzz, tr_zzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zzzz_xxxx[i] = 2.0 * tr_zzzzz_xxxx[i] * tbe_0 + 2.0 * tr_zzzz_xxxxz[i] * tke_0 - 4.0 * tr_zzz_xxxx[i];

        tr_0_0_z_zzzz_xxxy[i] = 2.0 * tr_zzzzz_xxxy[i] * tbe_0 + 2.0 * tr_zzzz_xxxyz[i] * tke_0 - 4.0 * tr_zzz_xxxy[i];

        tr_0_0_z_zzzz_xxxz[i] = 2.0 * tr_zzzzz_xxxz[i] * tbe_0 + 2.0 * tr_zzzz_xxxzz[i] * tke_0 - 4.0 * tr_zzz_xxxz[i] - tr_zzzz_xxx[i];

        tr_0_0_z_zzzz_xxyy[i] = 2.0 * tr_zzzzz_xxyy[i] * tbe_0 + 2.0 * tr_zzzz_xxyyz[i] * tke_0 - 4.0 * tr_zzz_xxyy[i];

        tr_0_0_z_zzzz_xxyz[i] = 2.0 * tr_zzzzz_xxyz[i] * tbe_0 + 2.0 * tr_zzzz_xxyzz[i] * tke_0 - 4.0 * tr_zzz_xxyz[i] - tr_zzzz_xxy[i];

        tr_0_0_z_zzzz_xxzz[i] = 2.0 * tr_zzzzz_xxzz[i] * tbe_0 + 2.0 * tr_zzzz_xxzzz[i] * tke_0 - 4.0 * tr_zzz_xxzz[i] - 2.0 * tr_zzzz_xxz[i];

        tr_0_0_z_zzzz_xyyy[i] = 2.0 * tr_zzzzz_xyyy[i] * tbe_0 + 2.0 * tr_zzzz_xyyyz[i] * tke_0 - 4.0 * tr_zzz_xyyy[i];

        tr_0_0_z_zzzz_xyyz[i] = 2.0 * tr_zzzzz_xyyz[i] * tbe_0 + 2.0 * tr_zzzz_xyyzz[i] * tke_0 - 4.0 * tr_zzz_xyyz[i] - tr_zzzz_xyy[i];

        tr_0_0_z_zzzz_xyzz[i] = 2.0 * tr_zzzzz_xyzz[i] * tbe_0 + 2.0 * tr_zzzz_xyzzz[i] * tke_0 - 4.0 * tr_zzz_xyzz[i] - 2.0 * tr_zzzz_xyz[i];

        tr_0_0_z_zzzz_xzzz[i] = 2.0 * tr_zzzzz_xzzz[i] * tbe_0 + 2.0 * tr_zzzz_xzzzz[i] * tke_0 - 4.0 * tr_zzz_xzzz[i] - 3.0 * tr_zzzz_xzz[i];

        tr_0_0_z_zzzz_yyyy[i] = 2.0 * tr_zzzzz_yyyy[i] * tbe_0 + 2.0 * tr_zzzz_yyyyz[i] * tke_0 - 4.0 * tr_zzz_yyyy[i];

        tr_0_0_z_zzzz_yyyz[i] = 2.0 * tr_zzzzz_yyyz[i] * tbe_0 + 2.0 * tr_zzzz_yyyzz[i] * tke_0 - 4.0 * tr_zzz_yyyz[i] - tr_zzzz_yyy[i];

        tr_0_0_z_zzzz_yyzz[i] = 2.0 * tr_zzzzz_yyzz[i] * tbe_0 + 2.0 * tr_zzzz_yyzzz[i] * tke_0 - 4.0 * tr_zzz_yyzz[i] - 2.0 * tr_zzzz_yyz[i];

        tr_0_0_z_zzzz_yzzz[i] = 2.0 * tr_zzzzz_yzzz[i] * tbe_0 + 2.0 * tr_zzzz_yzzzz[i] * tke_0 - 4.0 * tr_zzz_yzzz[i] - 3.0 * tr_zzzz_yzz[i];

        tr_0_0_z_zzzz_zzzz[i] = 2.0 * tr_zzzzz_zzzz[i] * tbe_0 + 2.0 * tr_zzzz_zzzzz[i] * tke_0 - 4.0 * tr_zzz_zzzz[i] - 4.0 * tr_zzzz_zzz[i];
    }

}

} // t2cgeom namespace

