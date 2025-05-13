#include "GeometricalDerivatives0X1ForGF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_gf(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_gf,
                       const int idx_op_gd,
                       const int idx_op_gg,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : GD

        auto to_xxxx_xx = pbuffer.data(idx_op_gd + i * 90 + 0);

        auto to_xxxx_xy = pbuffer.data(idx_op_gd + i * 90 + 1);

        auto to_xxxx_xz = pbuffer.data(idx_op_gd + i * 90 + 2);

        auto to_xxxx_yy = pbuffer.data(idx_op_gd + i * 90 + 3);

        auto to_xxxx_yz = pbuffer.data(idx_op_gd + i * 90 + 4);

        auto to_xxxx_zz = pbuffer.data(idx_op_gd + i * 90 + 5);

        auto to_xxxy_xx = pbuffer.data(idx_op_gd + i * 90 + 6);

        auto to_xxxy_xy = pbuffer.data(idx_op_gd + i * 90 + 7);

        auto to_xxxy_xz = pbuffer.data(idx_op_gd + i * 90 + 8);

        auto to_xxxy_yy = pbuffer.data(idx_op_gd + i * 90 + 9);

        auto to_xxxy_yz = pbuffer.data(idx_op_gd + i * 90 + 10);

        auto to_xxxy_zz = pbuffer.data(idx_op_gd + i * 90 + 11);

        auto to_xxxz_xx = pbuffer.data(idx_op_gd + i * 90 + 12);

        auto to_xxxz_xy = pbuffer.data(idx_op_gd + i * 90 + 13);

        auto to_xxxz_xz = pbuffer.data(idx_op_gd + i * 90 + 14);

        auto to_xxxz_yy = pbuffer.data(idx_op_gd + i * 90 + 15);

        auto to_xxxz_yz = pbuffer.data(idx_op_gd + i * 90 + 16);

        auto to_xxxz_zz = pbuffer.data(idx_op_gd + i * 90 + 17);

        auto to_xxyy_xx = pbuffer.data(idx_op_gd + i * 90 + 18);

        auto to_xxyy_xy = pbuffer.data(idx_op_gd + i * 90 + 19);

        auto to_xxyy_xz = pbuffer.data(idx_op_gd + i * 90 + 20);

        auto to_xxyy_yy = pbuffer.data(idx_op_gd + i * 90 + 21);

        auto to_xxyy_yz = pbuffer.data(idx_op_gd + i * 90 + 22);

        auto to_xxyy_zz = pbuffer.data(idx_op_gd + i * 90 + 23);

        auto to_xxyz_xx = pbuffer.data(idx_op_gd + i * 90 + 24);

        auto to_xxyz_xy = pbuffer.data(idx_op_gd + i * 90 + 25);

        auto to_xxyz_xz = pbuffer.data(idx_op_gd + i * 90 + 26);

        auto to_xxyz_yy = pbuffer.data(idx_op_gd + i * 90 + 27);

        auto to_xxyz_yz = pbuffer.data(idx_op_gd + i * 90 + 28);

        auto to_xxyz_zz = pbuffer.data(idx_op_gd + i * 90 + 29);

        auto to_xxzz_xx = pbuffer.data(idx_op_gd + i * 90 + 30);

        auto to_xxzz_xy = pbuffer.data(idx_op_gd + i * 90 + 31);

        auto to_xxzz_xz = pbuffer.data(idx_op_gd + i * 90 + 32);

        auto to_xxzz_yy = pbuffer.data(idx_op_gd + i * 90 + 33);

        auto to_xxzz_yz = pbuffer.data(idx_op_gd + i * 90 + 34);

        auto to_xxzz_zz = pbuffer.data(idx_op_gd + i * 90 + 35);

        auto to_xyyy_xx = pbuffer.data(idx_op_gd + i * 90 + 36);

        auto to_xyyy_xy = pbuffer.data(idx_op_gd + i * 90 + 37);

        auto to_xyyy_xz = pbuffer.data(idx_op_gd + i * 90 + 38);

        auto to_xyyy_yy = pbuffer.data(idx_op_gd + i * 90 + 39);

        auto to_xyyy_yz = pbuffer.data(idx_op_gd + i * 90 + 40);

        auto to_xyyy_zz = pbuffer.data(idx_op_gd + i * 90 + 41);

        auto to_xyyz_xx = pbuffer.data(idx_op_gd + i * 90 + 42);

        auto to_xyyz_xy = pbuffer.data(idx_op_gd + i * 90 + 43);

        auto to_xyyz_xz = pbuffer.data(idx_op_gd + i * 90 + 44);

        auto to_xyyz_yy = pbuffer.data(idx_op_gd + i * 90 + 45);

        auto to_xyyz_yz = pbuffer.data(idx_op_gd + i * 90 + 46);

        auto to_xyyz_zz = pbuffer.data(idx_op_gd + i * 90 + 47);

        auto to_xyzz_xx = pbuffer.data(idx_op_gd + i * 90 + 48);

        auto to_xyzz_xy = pbuffer.data(idx_op_gd + i * 90 + 49);

        auto to_xyzz_xz = pbuffer.data(idx_op_gd + i * 90 + 50);

        auto to_xyzz_yy = pbuffer.data(idx_op_gd + i * 90 + 51);

        auto to_xyzz_yz = pbuffer.data(idx_op_gd + i * 90 + 52);

        auto to_xyzz_zz = pbuffer.data(idx_op_gd + i * 90 + 53);

        auto to_xzzz_xx = pbuffer.data(idx_op_gd + i * 90 + 54);

        auto to_xzzz_xy = pbuffer.data(idx_op_gd + i * 90 + 55);

        auto to_xzzz_xz = pbuffer.data(idx_op_gd + i * 90 + 56);

        auto to_xzzz_yy = pbuffer.data(idx_op_gd + i * 90 + 57);

        auto to_xzzz_yz = pbuffer.data(idx_op_gd + i * 90 + 58);

        auto to_xzzz_zz = pbuffer.data(idx_op_gd + i * 90 + 59);

        auto to_yyyy_xx = pbuffer.data(idx_op_gd + i * 90 + 60);

        auto to_yyyy_xy = pbuffer.data(idx_op_gd + i * 90 + 61);

        auto to_yyyy_xz = pbuffer.data(idx_op_gd + i * 90 + 62);

        auto to_yyyy_yy = pbuffer.data(idx_op_gd + i * 90 + 63);

        auto to_yyyy_yz = pbuffer.data(idx_op_gd + i * 90 + 64);

        auto to_yyyy_zz = pbuffer.data(idx_op_gd + i * 90 + 65);

        auto to_yyyz_xx = pbuffer.data(idx_op_gd + i * 90 + 66);

        auto to_yyyz_xy = pbuffer.data(idx_op_gd + i * 90 + 67);

        auto to_yyyz_xz = pbuffer.data(idx_op_gd + i * 90 + 68);

        auto to_yyyz_yy = pbuffer.data(idx_op_gd + i * 90 + 69);

        auto to_yyyz_yz = pbuffer.data(idx_op_gd + i * 90 + 70);

        auto to_yyyz_zz = pbuffer.data(idx_op_gd + i * 90 + 71);

        auto to_yyzz_xx = pbuffer.data(idx_op_gd + i * 90 + 72);

        auto to_yyzz_xy = pbuffer.data(idx_op_gd + i * 90 + 73);

        auto to_yyzz_xz = pbuffer.data(idx_op_gd + i * 90 + 74);

        auto to_yyzz_yy = pbuffer.data(idx_op_gd + i * 90 + 75);

        auto to_yyzz_yz = pbuffer.data(idx_op_gd + i * 90 + 76);

        auto to_yyzz_zz = pbuffer.data(idx_op_gd + i * 90 + 77);

        auto to_yzzz_xx = pbuffer.data(idx_op_gd + i * 90 + 78);

        auto to_yzzz_xy = pbuffer.data(idx_op_gd + i * 90 + 79);

        auto to_yzzz_xz = pbuffer.data(idx_op_gd + i * 90 + 80);

        auto to_yzzz_yy = pbuffer.data(idx_op_gd + i * 90 + 81);

        auto to_yzzz_yz = pbuffer.data(idx_op_gd + i * 90 + 82);

        auto to_yzzz_zz = pbuffer.data(idx_op_gd + i * 90 + 83);

        auto to_zzzz_xx = pbuffer.data(idx_op_gd + i * 90 + 84);

        auto to_zzzz_xy = pbuffer.data(idx_op_gd + i * 90 + 85);

        auto to_zzzz_xz = pbuffer.data(idx_op_gd + i * 90 + 86);

        auto to_zzzz_yy = pbuffer.data(idx_op_gd + i * 90 + 87);

        auto to_zzzz_yz = pbuffer.data(idx_op_gd + i * 90 + 88);

        auto to_zzzz_zz = pbuffer.data(idx_op_gd + i * 90 + 89);

        // Set up components of auxiliary buffer : GG

        auto to_xxxx_xxxx = pbuffer.data(idx_op_gg + i * 225 + 0);

        auto to_xxxx_xxxy = pbuffer.data(idx_op_gg + i * 225 + 1);

        auto to_xxxx_xxxz = pbuffer.data(idx_op_gg + i * 225 + 2);

        auto to_xxxx_xxyy = pbuffer.data(idx_op_gg + i * 225 + 3);

        auto to_xxxx_xxyz = pbuffer.data(idx_op_gg + i * 225 + 4);

        auto to_xxxx_xxzz = pbuffer.data(idx_op_gg + i * 225 + 5);

        auto to_xxxx_xyyy = pbuffer.data(idx_op_gg + i * 225 + 6);

        auto to_xxxx_xyyz = pbuffer.data(idx_op_gg + i * 225 + 7);

        auto to_xxxx_xyzz = pbuffer.data(idx_op_gg + i * 225 + 8);

        auto to_xxxx_xzzz = pbuffer.data(idx_op_gg + i * 225 + 9);

        auto to_xxxx_yyyy = pbuffer.data(idx_op_gg + i * 225 + 10);

        auto to_xxxx_yyyz = pbuffer.data(idx_op_gg + i * 225 + 11);

        auto to_xxxx_yyzz = pbuffer.data(idx_op_gg + i * 225 + 12);

        auto to_xxxx_yzzz = pbuffer.data(idx_op_gg + i * 225 + 13);

        auto to_xxxx_zzzz = pbuffer.data(idx_op_gg + i * 225 + 14);

        auto to_xxxy_xxxx = pbuffer.data(idx_op_gg + i * 225 + 15);

        auto to_xxxy_xxxy = pbuffer.data(idx_op_gg + i * 225 + 16);

        auto to_xxxy_xxxz = pbuffer.data(idx_op_gg + i * 225 + 17);

        auto to_xxxy_xxyy = pbuffer.data(idx_op_gg + i * 225 + 18);

        auto to_xxxy_xxyz = pbuffer.data(idx_op_gg + i * 225 + 19);

        auto to_xxxy_xxzz = pbuffer.data(idx_op_gg + i * 225 + 20);

        auto to_xxxy_xyyy = pbuffer.data(idx_op_gg + i * 225 + 21);

        auto to_xxxy_xyyz = pbuffer.data(idx_op_gg + i * 225 + 22);

        auto to_xxxy_xyzz = pbuffer.data(idx_op_gg + i * 225 + 23);

        auto to_xxxy_xzzz = pbuffer.data(idx_op_gg + i * 225 + 24);

        auto to_xxxy_yyyy = pbuffer.data(idx_op_gg + i * 225 + 25);

        auto to_xxxy_yyyz = pbuffer.data(idx_op_gg + i * 225 + 26);

        auto to_xxxy_yyzz = pbuffer.data(idx_op_gg + i * 225 + 27);

        auto to_xxxy_yzzz = pbuffer.data(idx_op_gg + i * 225 + 28);

        auto to_xxxy_zzzz = pbuffer.data(idx_op_gg + i * 225 + 29);

        auto to_xxxz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 30);

        auto to_xxxz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 31);

        auto to_xxxz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 32);

        auto to_xxxz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 33);

        auto to_xxxz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 34);

        auto to_xxxz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 35);

        auto to_xxxz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 36);

        auto to_xxxz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 37);

        auto to_xxxz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 38);

        auto to_xxxz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 39);

        auto to_xxxz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 40);

        auto to_xxxz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 41);

        auto to_xxxz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 42);

        auto to_xxxz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 43);

        auto to_xxxz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 44);

        auto to_xxyy_xxxx = pbuffer.data(idx_op_gg + i * 225 + 45);

        auto to_xxyy_xxxy = pbuffer.data(idx_op_gg + i * 225 + 46);

        auto to_xxyy_xxxz = pbuffer.data(idx_op_gg + i * 225 + 47);

        auto to_xxyy_xxyy = pbuffer.data(idx_op_gg + i * 225 + 48);

        auto to_xxyy_xxyz = pbuffer.data(idx_op_gg + i * 225 + 49);

        auto to_xxyy_xxzz = pbuffer.data(idx_op_gg + i * 225 + 50);

        auto to_xxyy_xyyy = pbuffer.data(idx_op_gg + i * 225 + 51);

        auto to_xxyy_xyyz = pbuffer.data(idx_op_gg + i * 225 + 52);

        auto to_xxyy_xyzz = pbuffer.data(idx_op_gg + i * 225 + 53);

        auto to_xxyy_xzzz = pbuffer.data(idx_op_gg + i * 225 + 54);

        auto to_xxyy_yyyy = pbuffer.data(idx_op_gg + i * 225 + 55);

        auto to_xxyy_yyyz = pbuffer.data(idx_op_gg + i * 225 + 56);

        auto to_xxyy_yyzz = pbuffer.data(idx_op_gg + i * 225 + 57);

        auto to_xxyy_yzzz = pbuffer.data(idx_op_gg + i * 225 + 58);

        auto to_xxyy_zzzz = pbuffer.data(idx_op_gg + i * 225 + 59);

        auto to_xxyz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 60);

        auto to_xxyz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 61);

        auto to_xxyz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 62);

        auto to_xxyz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 63);

        auto to_xxyz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 64);

        auto to_xxyz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 65);

        auto to_xxyz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 66);

        auto to_xxyz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 67);

        auto to_xxyz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 68);

        auto to_xxyz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 69);

        auto to_xxyz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 70);

        auto to_xxyz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 71);

        auto to_xxyz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 72);

        auto to_xxyz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 73);

        auto to_xxyz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 74);

        auto to_xxzz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 75);

        auto to_xxzz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 76);

        auto to_xxzz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 77);

        auto to_xxzz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 78);

        auto to_xxzz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 79);

        auto to_xxzz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 80);

        auto to_xxzz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 81);

        auto to_xxzz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 82);

        auto to_xxzz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 83);

        auto to_xxzz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 84);

        auto to_xxzz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 85);

        auto to_xxzz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 86);

        auto to_xxzz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 87);

        auto to_xxzz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 88);

        auto to_xxzz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 89);

        auto to_xyyy_xxxx = pbuffer.data(idx_op_gg + i * 225 + 90);

        auto to_xyyy_xxxy = pbuffer.data(idx_op_gg + i * 225 + 91);

        auto to_xyyy_xxxz = pbuffer.data(idx_op_gg + i * 225 + 92);

        auto to_xyyy_xxyy = pbuffer.data(idx_op_gg + i * 225 + 93);

        auto to_xyyy_xxyz = pbuffer.data(idx_op_gg + i * 225 + 94);

        auto to_xyyy_xxzz = pbuffer.data(idx_op_gg + i * 225 + 95);

        auto to_xyyy_xyyy = pbuffer.data(idx_op_gg + i * 225 + 96);

        auto to_xyyy_xyyz = pbuffer.data(idx_op_gg + i * 225 + 97);

        auto to_xyyy_xyzz = pbuffer.data(idx_op_gg + i * 225 + 98);

        auto to_xyyy_xzzz = pbuffer.data(idx_op_gg + i * 225 + 99);

        auto to_xyyy_yyyy = pbuffer.data(idx_op_gg + i * 225 + 100);

        auto to_xyyy_yyyz = pbuffer.data(idx_op_gg + i * 225 + 101);

        auto to_xyyy_yyzz = pbuffer.data(idx_op_gg + i * 225 + 102);

        auto to_xyyy_yzzz = pbuffer.data(idx_op_gg + i * 225 + 103);

        auto to_xyyy_zzzz = pbuffer.data(idx_op_gg + i * 225 + 104);

        auto to_xyyz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 105);

        auto to_xyyz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 106);

        auto to_xyyz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 107);

        auto to_xyyz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 108);

        auto to_xyyz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 109);

        auto to_xyyz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 110);

        auto to_xyyz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 111);

        auto to_xyyz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 112);

        auto to_xyyz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 113);

        auto to_xyyz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 114);

        auto to_xyyz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 115);

        auto to_xyyz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 116);

        auto to_xyyz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 117);

        auto to_xyyz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 118);

        auto to_xyyz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 119);

        auto to_xyzz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 120);

        auto to_xyzz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 121);

        auto to_xyzz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 122);

        auto to_xyzz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 123);

        auto to_xyzz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 124);

        auto to_xyzz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 125);

        auto to_xyzz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 126);

        auto to_xyzz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 127);

        auto to_xyzz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 128);

        auto to_xyzz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 129);

        auto to_xyzz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 130);

        auto to_xyzz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 131);

        auto to_xyzz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 132);

        auto to_xyzz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 133);

        auto to_xyzz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 134);

        auto to_xzzz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 135);

        auto to_xzzz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 136);

        auto to_xzzz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 137);

        auto to_xzzz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 138);

        auto to_xzzz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 139);

        auto to_xzzz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 140);

        auto to_xzzz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 141);

        auto to_xzzz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 142);

        auto to_xzzz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 143);

        auto to_xzzz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 144);

        auto to_xzzz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 145);

        auto to_xzzz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 146);

        auto to_xzzz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 147);

        auto to_xzzz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 148);

        auto to_xzzz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 149);

        auto to_yyyy_xxxx = pbuffer.data(idx_op_gg + i * 225 + 150);

        auto to_yyyy_xxxy = pbuffer.data(idx_op_gg + i * 225 + 151);

        auto to_yyyy_xxxz = pbuffer.data(idx_op_gg + i * 225 + 152);

        auto to_yyyy_xxyy = pbuffer.data(idx_op_gg + i * 225 + 153);

        auto to_yyyy_xxyz = pbuffer.data(idx_op_gg + i * 225 + 154);

        auto to_yyyy_xxzz = pbuffer.data(idx_op_gg + i * 225 + 155);

        auto to_yyyy_xyyy = pbuffer.data(idx_op_gg + i * 225 + 156);

        auto to_yyyy_xyyz = pbuffer.data(idx_op_gg + i * 225 + 157);

        auto to_yyyy_xyzz = pbuffer.data(idx_op_gg + i * 225 + 158);

        auto to_yyyy_xzzz = pbuffer.data(idx_op_gg + i * 225 + 159);

        auto to_yyyy_yyyy = pbuffer.data(idx_op_gg + i * 225 + 160);

        auto to_yyyy_yyyz = pbuffer.data(idx_op_gg + i * 225 + 161);

        auto to_yyyy_yyzz = pbuffer.data(idx_op_gg + i * 225 + 162);

        auto to_yyyy_yzzz = pbuffer.data(idx_op_gg + i * 225 + 163);

        auto to_yyyy_zzzz = pbuffer.data(idx_op_gg + i * 225 + 164);

        auto to_yyyz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 165);

        auto to_yyyz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 166);

        auto to_yyyz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 167);

        auto to_yyyz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 168);

        auto to_yyyz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 169);

        auto to_yyyz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 170);

        auto to_yyyz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 171);

        auto to_yyyz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 172);

        auto to_yyyz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 173);

        auto to_yyyz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 174);

        auto to_yyyz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 175);

        auto to_yyyz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 176);

        auto to_yyyz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 177);

        auto to_yyyz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 178);

        auto to_yyyz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 179);

        auto to_yyzz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 180);

        auto to_yyzz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 181);

        auto to_yyzz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 182);

        auto to_yyzz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 183);

        auto to_yyzz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 184);

        auto to_yyzz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 185);

        auto to_yyzz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 186);

        auto to_yyzz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 187);

        auto to_yyzz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 188);

        auto to_yyzz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 189);

        auto to_yyzz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 190);

        auto to_yyzz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 191);

        auto to_yyzz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 192);

        auto to_yyzz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 193);

        auto to_yyzz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 194);

        auto to_yzzz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 195);

        auto to_yzzz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 196);

        auto to_yzzz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 197);

        auto to_yzzz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 198);

        auto to_yzzz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 199);

        auto to_yzzz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 200);

        auto to_yzzz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 201);

        auto to_yzzz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 202);

        auto to_yzzz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 203);

        auto to_yzzz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 204);

        auto to_yzzz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 205);

        auto to_yzzz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 206);

        auto to_yzzz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 207);

        auto to_yzzz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 208);

        auto to_yzzz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 209);

        auto to_zzzz_xxxx = pbuffer.data(idx_op_gg + i * 225 + 210);

        auto to_zzzz_xxxy = pbuffer.data(idx_op_gg + i * 225 + 211);

        auto to_zzzz_xxxz = pbuffer.data(idx_op_gg + i * 225 + 212);

        auto to_zzzz_xxyy = pbuffer.data(idx_op_gg + i * 225 + 213);

        auto to_zzzz_xxyz = pbuffer.data(idx_op_gg + i * 225 + 214);

        auto to_zzzz_xxzz = pbuffer.data(idx_op_gg + i * 225 + 215);

        auto to_zzzz_xyyy = pbuffer.data(idx_op_gg + i * 225 + 216);

        auto to_zzzz_xyyz = pbuffer.data(idx_op_gg + i * 225 + 217);

        auto to_zzzz_xyzz = pbuffer.data(idx_op_gg + i * 225 + 218);

        auto to_zzzz_xzzz = pbuffer.data(idx_op_gg + i * 225 + 219);

        auto to_zzzz_yyyy = pbuffer.data(idx_op_gg + i * 225 + 220);

        auto to_zzzz_yyyz = pbuffer.data(idx_op_gg + i * 225 + 221);

        auto to_zzzz_yyzz = pbuffer.data(idx_op_gg + i * 225 + 222);

        auto to_zzzz_yzzz = pbuffer.data(idx_op_gg + i * 225 + 223);

        auto to_zzzz_zzzz = pbuffer.data(idx_op_gg + i * 225 + 224);

        // Set up 0-10 components of targeted buffer : GF

        auto to_0_x_xxxx_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 0);

        auto to_0_x_xxxx_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 1);

        auto to_0_x_xxxx_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 2);

        auto to_0_x_xxxx_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 3);

        auto to_0_x_xxxx_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 4);

        auto to_0_x_xxxx_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 5);

        auto to_0_x_xxxx_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 6);

        auto to_0_x_xxxx_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 7);

        auto to_0_x_xxxx_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 8);

        auto to_0_x_xxxx_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_0_x_xxxx_xxx, to_0_x_xxxx_xxy, to_0_x_xxxx_xxz, to_0_x_xxxx_xyy, to_0_x_xxxx_xyz, to_0_x_xxxx_xzz, to_0_x_xxxx_yyy, to_0_x_xxxx_yyz, to_0_x_xxxx_yzz, to_0_x_xxxx_zzz, to_xxxx_xx, to_xxxx_xxxx, to_xxxx_xxxy, to_xxxx_xxxz, to_xxxx_xxyy, to_xxxx_xxyz, to_xxxx_xxzz, to_xxxx_xy, to_xxxx_xyyy, to_xxxx_xyyz, to_xxxx_xyzz, to_xxxx_xz, to_xxxx_xzzz, to_xxxx_yy, to_xxxx_yz, to_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxx_xxx[k] = -3.0 * to_xxxx_xx[k] + 2.0 * to_xxxx_xxxx[k] * tke_0;

            to_0_x_xxxx_xxy[k] = -2.0 * to_xxxx_xy[k] + 2.0 * to_xxxx_xxxy[k] * tke_0;

            to_0_x_xxxx_xxz[k] = -2.0 * to_xxxx_xz[k] + 2.0 * to_xxxx_xxxz[k] * tke_0;

            to_0_x_xxxx_xyy[k] = -to_xxxx_yy[k] + 2.0 * to_xxxx_xxyy[k] * tke_0;

            to_0_x_xxxx_xyz[k] = -to_xxxx_yz[k] + 2.0 * to_xxxx_xxyz[k] * tke_0;

            to_0_x_xxxx_xzz[k] = -to_xxxx_zz[k] + 2.0 * to_xxxx_xxzz[k] * tke_0;

            to_0_x_xxxx_yyy[k] = 2.0 * to_xxxx_xyyy[k] * tke_0;

            to_0_x_xxxx_yyz[k] = 2.0 * to_xxxx_xyyz[k] * tke_0;

            to_0_x_xxxx_yzz[k] = 2.0 * to_xxxx_xyzz[k] * tke_0;

            to_0_x_xxxx_zzz[k] = 2.0 * to_xxxx_xzzz[k] * tke_0;
        }

        // Set up 10-20 components of targeted buffer : GF

        auto to_0_x_xxxy_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 10);

        auto to_0_x_xxxy_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 11);

        auto to_0_x_xxxy_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 12);

        auto to_0_x_xxxy_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 13);

        auto to_0_x_xxxy_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 14);

        auto to_0_x_xxxy_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 15);

        auto to_0_x_xxxy_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 16);

        auto to_0_x_xxxy_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 17);

        auto to_0_x_xxxy_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 18);

        auto to_0_x_xxxy_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_0_x_xxxy_xxx, to_0_x_xxxy_xxy, to_0_x_xxxy_xxz, to_0_x_xxxy_xyy, to_0_x_xxxy_xyz, to_0_x_xxxy_xzz, to_0_x_xxxy_yyy, to_0_x_xxxy_yyz, to_0_x_xxxy_yzz, to_0_x_xxxy_zzz, to_xxxy_xx, to_xxxy_xxxx, to_xxxy_xxxy, to_xxxy_xxxz, to_xxxy_xxyy, to_xxxy_xxyz, to_xxxy_xxzz, to_xxxy_xy, to_xxxy_xyyy, to_xxxy_xyyz, to_xxxy_xyzz, to_xxxy_xz, to_xxxy_xzzz, to_xxxy_yy, to_xxxy_yz, to_xxxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxy_xxx[k] = -3.0 * to_xxxy_xx[k] + 2.0 * to_xxxy_xxxx[k] * tke_0;

            to_0_x_xxxy_xxy[k] = -2.0 * to_xxxy_xy[k] + 2.0 * to_xxxy_xxxy[k] * tke_0;

            to_0_x_xxxy_xxz[k] = -2.0 * to_xxxy_xz[k] + 2.0 * to_xxxy_xxxz[k] * tke_0;

            to_0_x_xxxy_xyy[k] = -to_xxxy_yy[k] + 2.0 * to_xxxy_xxyy[k] * tke_0;

            to_0_x_xxxy_xyz[k] = -to_xxxy_yz[k] + 2.0 * to_xxxy_xxyz[k] * tke_0;

            to_0_x_xxxy_xzz[k] = -to_xxxy_zz[k] + 2.0 * to_xxxy_xxzz[k] * tke_0;

            to_0_x_xxxy_yyy[k] = 2.0 * to_xxxy_xyyy[k] * tke_0;

            to_0_x_xxxy_yyz[k] = 2.0 * to_xxxy_xyyz[k] * tke_0;

            to_0_x_xxxy_yzz[k] = 2.0 * to_xxxy_xyzz[k] * tke_0;

            to_0_x_xxxy_zzz[k] = 2.0 * to_xxxy_xzzz[k] * tke_0;
        }

        // Set up 20-30 components of targeted buffer : GF

        auto to_0_x_xxxz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 20);

        auto to_0_x_xxxz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 21);

        auto to_0_x_xxxz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 22);

        auto to_0_x_xxxz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 23);

        auto to_0_x_xxxz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 24);

        auto to_0_x_xxxz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 25);

        auto to_0_x_xxxz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 26);

        auto to_0_x_xxxz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 27);

        auto to_0_x_xxxz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 28);

        auto to_0_x_xxxz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_0_x_xxxz_xxx, to_0_x_xxxz_xxy, to_0_x_xxxz_xxz, to_0_x_xxxz_xyy, to_0_x_xxxz_xyz, to_0_x_xxxz_xzz, to_0_x_xxxz_yyy, to_0_x_xxxz_yyz, to_0_x_xxxz_yzz, to_0_x_xxxz_zzz, to_xxxz_xx, to_xxxz_xxxx, to_xxxz_xxxy, to_xxxz_xxxz, to_xxxz_xxyy, to_xxxz_xxyz, to_xxxz_xxzz, to_xxxz_xy, to_xxxz_xyyy, to_xxxz_xyyz, to_xxxz_xyzz, to_xxxz_xz, to_xxxz_xzzz, to_xxxz_yy, to_xxxz_yz, to_xxxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxz_xxx[k] = -3.0 * to_xxxz_xx[k] + 2.0 * to_xxxz_xxxx[k] * tke_0;

            to_0_x_xxxz_xxy[k] = -2.0 * to_xxxz_xy[k] + 2.0 * to_xxxz_xxxy[k] * tke_0;

            to_0_x_xxxz_xxz[k] = -2.0 * to_xxxz_xz[k] + 2.0 * to_xxxz_xxxz[k] * tke_0;

            to_0_x_xxxz_xyy[k] = -to_xxxz_yy[k] + 2.0 * to_xxxz_xxyy[k] * tke_0;

            to_0_x_xxxz_xyz[k] = -to_xxxz_yz[k] + 2.0 * to_xxxz_xxyz[k] * tke_0;

            to_0_x_xxxz_xzz[k] = -to_xxxz_zz[k] + 2.0 * to_xxxz_xxzz[k] * tke_0;

            to_0_x_xxxz_yyy[k] = 2.0 * to_xxxz_xyyy[k] * tke_0;

            to_0_x_xxxz_yyz[k] = 2.0 * to_xxxz_xyyz[k] * tke_0;

            to_0_x_xxxz_yzz[k] = 2.0 * to_xxxz_xyzz[k] * tke_0;

            to_0_x_xxxz_zzz[k] = 2.0 * to_xxxz_xzzz[k] * tke_0;
        }

        // Set up 30-40 components of targeted buffer : GF

        auto to_0_x_xxyy_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 30);

        auto to_0_x_xxyy_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 31);

        auto to_0_x_xxyy_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 32);

        auto to_0_x_xxyy_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 33);

        auto to_0_x_xxyy_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 34);

        auto to_0_x_xxyy_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 35);

        auto to_0_x_xxyy_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 36);

        auto to_0_x_xxyy_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 37);

        auto to_0_x_xxyy_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 38);

        auto to_0_x_xxyy_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_0_x_xxyy_xxx, to_0_x_xxyy_xxy, to_0_x_xxyy_xxz, to_0_x_xxyy_xyy, to_0_x_xxyy_xyz, to_0_x_xxyy_xzz, to_0_x_xxyy_yyy, to_0_x_xxyy_yyz, to_0_x_xxyy_yzz, to_0_x_xxyy_zzz, to_xxyy_xx, to_xxyy_xxxx, to_xxyy_xxxy, to_xxyy_xxxz, to_xxyy_xxyy, to_xxyy_xxyz, to_xxyy_xxzz, to_xxyy_xy, to_xxyy_xyyy, to_xxyy_xyyz, to_xxyy_xyzz, to_xxyy_xz, to_xxyy_xzzz, to_xxyy_yy, to_xxyy_yz, to_xxyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxyy_xxx[k] = -3.0 * to_xxyy_xx[k] + 2.0 * to_xxyy_xxxx[k] * tke_0;

            to_0_x_xxyy_xxy[k] = -2.0 * to_xxyy_xy[k] + 2.0 * to_xxyy_xxxy[k] * tke_0;

            to_0_x_xxyy_xxz[k] = -2.0 * to_xxyy_xz[k] + 2.0 * to_xxyy_xxxz[k] * tke_0;

            to_0_x_xxyy_xyy[k] = -to_xxyy_yy[k] + 2.0 * to_xxyy_xxyy[k] * tke_0;

            to_0_x_xxyy_xyz[k] = -to_xxyy_yz[k] + 2.0 * to_xxyy_xxyz[k] * tke_0;

            to_0_x_xxyy_xzz[k] = -to_xxyy_zz[k] + 2.0 * to_xxyy_xxzz[k] * tke_0;

            to_0_x_xxyy_yyy[k] = 2.0 * to_xxyy_xyyy[k] * tke_0;

            to_0_x_xxyy_yyz[k] = 2.0 * to_xxyy_xyyz[k] * tke_0;

            to_0_x_xxyy_yzz[k] = 2.0 * to_xxyy_xyzz[k] * tke_0;

            to_0_x_xxyy_zzz[k] = 2.0 * to_xxyy_xzzz[k] * tke_0;
        }

        // Set up 40-50 components of targeted buffer : GF

        auto to_0_x_xxyz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 40);

        auto to_0_x_xxyz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 41);

        auto to_0_x_xxyz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 42);

        auto to_0_x_xxyz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 43);

        auto to_0_x_xxyz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 44);

        auto to_0_x_xxyz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 45);

        auto to_0_x_xxyz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 46);

        auto to_0_x_xxyz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 47);

        auto to_0_x_xxyz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 48);

        auto to_0_x_xxyz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_0_x_xxyz_xxx, to_0_x_xxyz_xxy, to_0_x_xxyz_xxz, to_0_x_xxyz_xyy, to_0_x_xxyz_xyz, to_0_x_xxyz_xzz, to_0_x_xxyz_yyy, to_0_x_xxyz_yyz, to_0_x_xxyz_yzz, to_0_x_xxyz_zzz, to_xxyz_xx, to_xxyz_xxxx, to_xxyz_xxxy, to_xxyz_xxxz, to_xxyz_xxyy, to_xxyz_xxyz, to_xxyz_xxzz, to_xxyz_xy, to_xxyz_xyyy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_xzzz, to_xxyz_yy, to_xxyz_yz, to_xxyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxyz_xxx[k] = -3.0 * to_xxyz_xx[k] + 2.0 * to_xxyz_xxxx[k] * tke_0;

            to_0_x_xxyz_xxy[k] = -2.0 * to_xxyz_xy[k] + 2.0 * to_xxyz_xxxy[k] * tke_0;

            to_0_x_xxyz_xxz[k] = -2.0 * to_xxyz_xz[k] + 2.0 * to_xxyz_xxxz[k] * tke_0;

            to_0_x_xxyz_xyy[k] = -to_xxyz_yy[k] + 2.0 * to_xxyz_xxyy[k] * tke_0;

            to_0_x_xxyz_xyz[k] = -to_xxyz_yz[k] + 2.0 * to_xxyz_xxyz[k] * tke_0;

            to_0_x_xxyz_xzz[k] = -to_xxyz_zz[k] + 2.0 * to_xxyz_xxzz[k] * tke_0;

            to_0_x_xxyz_yyy[k] = 2.0 * to_xxyz_xyyy[k] * tke_0;

            to_0_x_xxyz_yyz[k] = 2.0 * to_xxyz_xyyz[k] * tke_0;

            to_0_x_xxyz_yzz[k] = 2.0 * to_xxyz_xyzz[k] * tke_0;

            to_0_x_xxyz_zzz[k] = 2.0 * to_xxyz_xzzz[k] * tke_0;
        }

        // Set up 50-60 components of targeted buffer : GF

        auto to_0_x_xxzz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 50);

        auto to_0_x_xxzz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 51);

        auto to_0_x_xxzz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 52);

        auto to_0_x_xxzz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 53);

        auto to_0_x_xxzz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 54);

        auto to_0_x_xxzz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 55);

        auto to_0_x_xxzz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 56);

        auto to_0_x_xxzz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 57);

        auto to_0_x_xxzz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 58);

        auto to_0_x_xxzz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_0_x_xxzz_xxx, to_0_x_xxzz_xxy, to_0_x_xxzz_xxz, to_0_x_xxzz_xyy, to_0_x_xxzz_xyz, to_0_x_xxzz_xzz, to_0_x_xxzz_yyy, to_0_x_xxzz_yyz, to_0_x_xxzz_yzz, to_0_x_xxzz_zzz, to_xxzz_xx, to_xxzz_xxxx, to_xxzz_xxxy, to_xxzz_xxxz, to_xxzz_xxyy, to_xxzz_xxyz, to_xxzz_xxzz, to_xxzz_xy, to_xxzz_xyyy, to_xxzz_xyyz, to_xxzz_xyzz, to_xxzz_xz, to_xxzz_xzzz, to_xxzz_yy, to_xxzz_yz, to_xxzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxzz_xxx[k] = -3.0 * to_xxzz_xx[k] + 2.0 * to_xxzz_xxxx[k] * tke_0;

            to_0_x_xxzz_xxy[k] = -2.0 * to_xxzz_xy[k] + 2.0 * to_xxzz_xxxy[k] * tke_0;

            to_0_x_xxzz_xxz[k] = -2.0 * to_xxzz_xz[k] + 2.0 * to_xxzz_xxxz[k] * tke_0;

            to_0_x_xxzz_xyy[k] = -to_xxzz_yy[k] + 2.0 * to_xxzz_xxyy[k] * tke_0;

            to_0_x_xxzz_xyz[k] = -to_xxzz_yz[k] + 2.0 * to_xxzz_xxyz[k] * tke_0;

            to_0_x_xxzz_xzz[k] = -to_xxzz_zz[k] + 2.0 * to_xxzz_xxzz[k] * tke_0;

            to_0_x_xxzz_yyy[k] = 2.0 * to_xxzz_xyyy[k] * tke_0;

            to_0_x_xxzz_yyz[k] = 2.0 * to_xxzz_xyyz[k] * tke_0;

            to_0_x_xxzz_yzz[k] = 2.0 * to_xxzz_xyzz[k] * tke_0;

            to_0_x_xxzz_zzz[k] = 2.0 * to_xxzz_xzzz[k] * tke_0;
        }

        // Set up 60-70 components of targeted buffer : GF

        auto to_0_x_xyyy_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 60);

        auto to_0_x_xyyy_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 61);

        auto to_0_x_xyyy_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 62);

        auto to_0_x_xyyy_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 63);

        auto to_0_x_xyyy_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 64);

        auto to_0_x_xyyy_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 65);

        auto to_0_x_xyyy_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 66);

        auto to_0_x_xyyy_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 67);

        auto to_0_x_xyyy_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 68);

        auto to_0_x_xyyy_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_0_x_xyyy_xxx, to_0_x_xyyy_xxy, to_0_x_xyyy_xxz, to_0_x_xyyy_xyy, to_0_x_xyyy_xyz, to_0_x_xyyy_xzz, to_0_x_xyyy_yyy, to_0_x_xyyy_yyz, to_0_x_xyyy_yzz, to_0_x_xyyy_zzz, to_xyyy_xx, to_xyyy_xxxx, to_xyyy_xxxy, to_xyyy_xxxz, to_xyyy_xxyy, to_xyyy_xxyz, to_xyyy_xxzz, to_xyyy_xy, to_xyyy_xyyy, to_xyyy_xyyz, to_xyyy_xyzz, to_xyyy_xz, to_xyyy_xzzz, to_xyyy_yy, to_xyyy_yz, to_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyyy_xxx[k] = -3.0 * to_xyyy_xx[k] + 2.0 * to_xyyy_xxxx[k] * tke_0;

            to_0_x_xyyy_xxy[k] = -2.0 * to_xyyy_xy[k] + 2.0 * to_xyyy_xxxy[k] * tke_0;

            to_0_x_xyyy_xxz[k] = -2.0 * to_xyyy_xz[k] + 2.0 * to_xyyy_xxxz[k] * tke_0;

            to_0_x_xyyy_xyy[k] = -to_xyyy_yy[k] + 2.0 * to_xyyy_xxyy[k] * tke_0;

            to_0_x_xyyy_xyz[k] = -to_xyyy_yz[k] + 2.0 * to_xyyy_xxyz[k] * tke_0;

            to_0_x_xyyy_xzz[k] = -to_xyyy_zz[k] + 2.0 * to_xyyy_xxzz[k] * tke_0;

            to_0_x_xyyy_yyy[k] = 2.0 * to_xyyy_xyyy[k] * tke_0;

            to_0_x_xyyy_yyz[k] = 2.0 * to_xyyy_xyyz[k] * tke_0;

            to_0_x_xyyy_yzz[k] = 2.0 * to_xyyy_xyzz[k] * tke_0;

            to_0_x_xyyy_zzz[k] = 2.0 * to_xyyy_xzzz[k] * tke_0;
        }

        // Set up 70-80 components of targeted buffer : GF

        auto to_0_x_xyyz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 70);

        auto to_0_x_xyyz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 71);

        auto to_0_x_xyyz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 72);

        auto to_0_x_xyyz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 73);

        auto to_0_x_xyyz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 74);

        auto to_0_x_xyyz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 75);

        auto to_0_x_xyyz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 76);

        auto to_0_x_xyyz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 77);

        auto to_0_x_xyyz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 78);

        auto to_0_x_xyyz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_0_x_xyyz_xxx, to_0_x_xyyz_xxy, to_0_x_xyyz_xxz, to_0_x_xyyz_xyy, to_0_x_xyyz_xyz, to_0_x_xyyz_xzz, to_0_x_xyyz_yyy, to_0_x_xyyz_yyz, to_0_x_xyyz_yzz, to_0_x_xyyz_zzz, to_xyyz_xx, to_xyyz_xxxx, to_xyyz_xxxy, to_xyyz_xxxz, to_xyyz_xxyy, to_xyyz_xxyz, to_xyyz_xxzz, to_xyyz_xy, to_xyyz_xyyy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_xzzz, to_xyyz_yy, to_xyyz_yz, to_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyyz_xxx[k] = -3.0 * to_xyyz_xx[k] + 2.0 * to_xyyz_xxxx[k] * tke_0;

            to_0_x_xyyz_xxy[k] = -2.0 * to_xyyz_xy[k] + 2.0 * to_xyyz_xxxy[k] * tke_0;

            to_0_x_xyyz_xxz[k] = -2.0 * to_xyyz_xz[k] + 2.0 * to_xyyz_xxxz[k] * tke_0;

            to_0_x_xyyz_xyy[k] = -to_xyyz_yy[k] + 2.0 * to_xyyz_xxyy[k] * tke_0;

            to_0_x_xyyz_xyz[k] = -to_xyyz_yz[k] + 2.0 * to_xyyz_xxyz[k] * tke_0;

            to_0_x_xyyz_xzz[k] = -to_xyyz_zz[k] + 2.0 * to_xyyz_xxzz[k] * tke_0;

            to_0_x_xyyz_yyy[k] = 2.0 * to_xyyz_xyyy[k] * tke_0;

            to_0_x_xyyz_yyz[k] = 2.0 * to_xyyz_xyyz[k] * tke_0;

            to_0_x_xyyz_yzz[k] = 2.0 * to_xyyz_xyzz[k] * tke_0;

            to_0_x_xyyz_zzz[k] = 2.0 * to_xyyz_xzzz[k] * tke_0;
        }

        // Set up 80-90 components of targeted buffer : GF

        auto to_0_x_xyzz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 80);

        auto to_0_x_xyzz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 81);

        auto to_0_x_xyzz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 82);

        auto to_0_x_xyzz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 83);

        auto to_0_x_xyzz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 84);

        auto to_0_x_xyzz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 85);

        auto to_0_x_xyzz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 86);

        auto to_0_x_xyzz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 87);

        auto to_0_x_xyzz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 88);

        auto to_0_x_xyzz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_0_x_xyzz_xxx, to_0_x_xyzz_xxy, to_0_x_xyzz_xxz, to_0_x_xyzz_xyy, to_0_x_xyzz_xyz, to_0_x_xyzz_xzz, to_0_x_xyzz_yyy, to_0_x_xyzz_yyz, to_0_x_xyzz_yzz, to_0_x_xyzz_zzz, to_xyzz_xx, to_xyzz_xxxx, to_xyzz_xxxy, to_xyzz_xxxz, to_xyzz_xxyy, to_xyzz_xxyz, to_xyzz_xxzz, to_xyzz_xy, to_xyzz_xyyy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_xzzz, to_xyzz_yy, to_xyzz_yz, to_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyzz_xxx[k] = -3.0 * to_xyzz_xx[k] + 2.0 * to_xyzz_xxxx[k] * tke_0;

            to_0_x_xyzz_xxy[k] = -2.0 * to_xyzz_xy[k] + 2.0 * to_xyzz_xxxy[k] * tke_0;

            to_0_x_xyzz_xxz[k] = -2.0 * to_xyzz_xz[k] + 2.0 * to_xyzz_xxxz[k] * tke_0;

            to_0_x_xyzz_xyy[k] = -to_xyzz_yy[k] + 2.0 * to_xyzz_xxyy[k] * tke_0;

            to_0_x_xyzz_xyz[k] = -to_xyzz_yz[k] + 2.0 * to_xyzz_xxyz[k] * tke_0;

            to_0_x_xyzz_xzz[k] = -to_xyzz_zz[k] + 2.0 * to_xyzz_xxzz[k] * tke_0;

            to_0_x_xyzz_yyy[k] = 2.0 * to_xyzz_xyyy[k] * tke_0;

            to_0_x_xyzz_yyz[k] = 2.0 * to_xyzz_xyyz[k] * tke_0;

            to_0_x_xyzz_yzz[k] = 2.0 * to_xyzz_xyzz[k] * tke_0;

            to_0_x_xyzz_zzz[k] = 2.0 * to_xyzz_xzzz[k] * tke_0;
        }

        // Set up 90-100 components of targeted buffer : GF

        auto to_0_x_xzzz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 90);

        auto to_0_x_xzzz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 91);

        auto to_0_x_xzzz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 92);

        auto to_0_x_xzzz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 93);

        auto to_0_x_xzzz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 94);

        auto to_0_x_xzzz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 95);

        auto to_0_x_xzzz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 96);

        auto to_0_x_xzzz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 97);

        auto to_0_x_xzzz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 98);

        auto to_0_x_xzzz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_0_x_xzzz_xxx, to_0_x_xzzz_xxy, to_0_x_xzzz_xxz, to_0_x_xzzz_xyy, to_0_x_xzzz_xyz, to_0_x_xzzz_xzz, to_0_x_xzzz_yyy, to_0_x_xzzz_yyz, to_0_x_xzzz_yzz, to_0_x_xzzz_zzz, to_xzzz_xx, to_xzzz_xxxx, to_xzzz_xxxy, to_xzzz_xxxz, to_xzzz_xxyy, to_xzzz_xxyz, to_xzzz_xxzz, to_xzzz_xy, to_xzzz_xyyy, to_xzzz_xyyz, to_xzzz_xyzz, to_xzzz_xz, to_xzzz_xzzz, to_xzzz_yy, to_xzzz_yz, to_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xzzz_xxx[k] = -3.0 * to_xzzz_xx[k] + 2.0 * to_xzzz_xxxx[k] * tke_0;

            to_0_x_xzzz_xxy[k] = -2.0 * to_xzzz_xy[k] + 2.0 * to_xzzz_xxxy[k] * tke_0;

            to_0_x_xzzz_xxz[k] = -2.0 * to_xzzz_xz[k] + 2.0 * to_xzzz_xxxz[k] * tke_0;

            to_0_x_xzzz_xyy[k] = -to_xzzz_yy[k] + 2.0 * to_xzzz_xxyy[k] * tke_0;

            to_0_x_xzzz_xyz[k] = -to_xzzz_yz[k] + 2.0 * to_xzzz_xxyz[k] * tke_0;

            to_0_x_xzzz_xzz[k] = -to_xzzz_zz[k] + 2.0 * to_xzzz_xxzz[k] * tke_0;

            to_0_x_xzzz_yyy[k] = 2.0 * to_xzzz_xyyy[k] * tke_0;

            to_0_x_xzzz_yyz[k] = 2.0 * to_xzzz_xyyz[k] * tke_0;

            to_0_x_xzzz_yzz[k] = 2.0 * to_xzzz_xyzz[k] * tke_0;

            to_0_x_xzzz_zzz[k] = 2.0 * to_xzzz_xzzz[k] * tke_0;
        }

        // Set up 100-110 components of targeted buffer : GF

        auto to_0_x_yyyy_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 100);

        auto to_0_x_yyyy_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 101);

        auto to_0_x_yyyy_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 102);

        auto to_0_x_yyyy_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 103);

        auto to_0_x_yyyy_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 104);

        auto to_0_x_yyyy_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 105);

        auto to_0_x_yyyy_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 106);

        auto to_0_x_yyyy_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 107);

        auto to_0_x_yyyy_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 108);

        auto to_0_x_yyyy_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_0_x_yyyy_xxx, to_0_x_yyyy_xxy, to_0_x_yyyy_xxz, to_0_x_yyyy_xyy, to_0_x_yyyy_xyz, to_0_x_yyyy_xzz, to_0_x_yyyy_yyy, to_0_x_yyyy_yyz, to_0_x_yyyy_yzz, to_0_x_yyyy_zzz, to_yyyy_xx, to_yyyy_xxxx, to_yyyy_xxxy, to_yyyy_xxxz, to_yyyy_xxyy, to_yyyy_xxyz, to_yyyy_xxzz, to_yyyy_xy, to_yyyy_xyyy, to_yyyy_xyyz, to_yyyy_xyzz, to_yyyy_xz, to_yyyy_xzzz, to_yyyy_yy, to_yyyy_yz, to_yyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyyy_xxx[k] = -3.0 * to_yyyy_xx[k] + 2.0 * to_yyyy_xxxx[k] * tke_0;

            to_0_x_yyyy_xxy[k] = -2.0 * to_yyyy_xy[k] + 2.0 * to_yyyy_xxxy[k] * tke_0;

            to_0_x_yyyy_xxz[k] = -2.0 * to_yyyy_xz[k] + 2.0 * to_yyyy_xxxz[k] * tke_0;

            to_0_x_yyyy_xyy[k] = -to_yyyy_yy[k] + 2.0 * to_yyyy_xxyy[k] * tke_0;

            to_0_x_yyyy_xyz[k] = -to_yyyy_yz[k] + 2.0 * to_yyyy_xxyz[k] * tke_0;

            to_0_x_yyyy_xzz[k] = -to_yyyy_zz[k] + 2.0 * to_yyyy_xxzz[k] * tke_0;

            to_0_x_yyyy_yyy[k] = 2.0 * to_yyyy_xyyy[k] * tke_0;

            to_0_x_yyyy_yyz[k] = 2.0 * to_yyyy_xyyz[k] * tke_0;

            to_0_x_yyyy_yzz[k] = 2.0 * to_yyyy_xyzz[k] * tke_0;

            to_0_x_yyyy_zzz[k] = 2.0 * to_yyyy_xzzz[k] * tke_0;
        }

        // Set up 110-120 components of targeted buffer : GF

        auto to_0_x_yyyz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 110);

        auto to_0_x_yyyz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 111);

        auto to_0_x_yyyz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 112);

        auto to_0_x_yyyz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 113);

        auto to_0_x_yyyz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 114);

        auto to_0_x_yyyz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 115);

        auto to_0_x_yyyz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 116);

        auto to_0_x_yyyz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 117);

        auto to_0_x_yyyz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 118);

        auto to_0_x_yyyz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_0_x_yyyz_xxx, to_0_x_yyyz_xxy, to_0_x_yyyz_xxz, to_0_x_yyyz_xyy, to_0_x_yyyz_xyz, to_0_x_yyyz_xzz, to_0_x_yyyz_yyy, to_0_x_yyyz_yyz, to_0_x_yyyz_yzz, to_0_x_yyyz_zzz, to_yyyz_xx, to_yyyz_xxxx, to_yyyz_xxxy, to_yyyz_xxxz, to_yyyz_xxyy, to_yyyz_xxyz, to_yyyz_xxzz, to_yyyz_xy, to_yyyz_xyyy, to_yyyz_xyyz, to_yyyz_xyzz, to_yyyz_xz, to_yyyz_xzzz, to_yyyz_yy, to_yyyz_yz, to_yyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyyz_xxx[k] = -3.0 * to_yyyz_xx[k] + 2.0 * to_yyyz_xxxx[k] * tke_0;

            to_0_x_yyyz_xxy[k] = -2.0 * to_yyyz_xy[k] + 2.0 * to_yyyz_xxxy[k] * tke_0;

            to_0_x_yyyz_xxz[k] = -2.0 * to_yyyz_xz[k] + 2.0 * to_yyyz_xxxz[k] * tke_0;

            to_0_x_yyyz_xyy[k] = -to_yyyz_yy[k] + 2.0 * to_yyyz_xxyy[k] * tke_0;

            to_0_x_yyyz_xyz[k] = -to_yyyz_yz[k] + 2.0 * to_yyyz_xxyz[k] * tke_0;

            to_0_x_yyyz_xzz[k] = -to_yyyz_zz[k] + 2.0 * to_yyyz_xxzz[k] * tke_0;

            to_0_x_yyyz_yyy[k] = 2.0 * to_yyyz_xyyy[k] * tke_0;

            to_0_x_yyyz_yyz[k] = 2.0 * to_yyyz_xyyz[k] * tke_0;

            to_0_x_yyyz_yzz[k] = 2.0 * to_yyyz_xyzz[k] * tke_0;

            to_0_x_yyyz_zzz[k] = 2.0 * to_yyyz_xzzz[k] * tke_0;
        }

        // Set up 120-130 components of targeted buffer : GF

        auto to_0_x_yyzz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 120);

        auto to_0_x_yyzz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 121);

        auto to_0_x_yyzz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 122);

        auto to_0_x_yyzz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 123);

        auto to_0_x_yyzz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 124);

        auto to_0_x_yyzz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 125);

        auto to_0_x_yyzz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 126);

        auto to_0_x_yyzz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 127);

        auto to_0_x_yyzz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 128);

        auto to_0_x_yyzz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_0_x_yyzz_xxx, to_0_x_yyzz_xxy, to_0_x_yyzz_xxz, to_0_x_yyzz_xyy, to_0_x_yyzz_xyz, to_0_x_yyzz_xzz, to_0_x_yyzz_yyy, to_0_x_yyzz_yyz, to_0_x_yyzz_yzz, to_0_x_yyzz_zzz, to_yyzz_xx, to_yyzz_xxxx, to_yyzz_xxxy, to_yyzz_xxxz, to_yyzz_xxyy, to_yyzz_xxyz, to_yyzz_xxzz, to_yyzz_xy, to_yyzz_xyyy, to_yyzz_xyyz, to_yyzz_xyzz, to_yyzz_xz, to_yyzz_xzzz, to_yyzz_yy, to_yyzz_yz, to_yyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyzz_xxx[k] = -3.0 * to_yyzz_xx[k] + 2.0 * to_yyzz_xxxx[k] * tke_0;

            to_0_x_yyzz_xxy[k] = -2.0 * to_yyzz_xy[k] + 2.0 * to_yyzz_xxxy[k] * tke_0;

            to_0_x_yyzz_xxz[k] = -2.0 * to_yyzz_xz[k] + 2.0 * to_yyzz_xxxz[k] * tke_0;

            to_0_x_yyzz_xyy[k] = -to_yyzz_yy[k] + 2.0 * to_yyzz_xxyy[k] * tke_0;

            to_0_x_yyzz_xyz[k] = -to_yyzz_yz[k] + 2.0 * to_yyzz_xxyz[k] * tke_0;

            to_0_x_yyzz_xzz[k] = -to_yyzz_zz[k] + 2.0 * to_yyzz_xxzz[k] * tke_0;

            to_0_x_yyzz_yyy[k] = 2.0 * to_yyzz_xyyy[k] * tke_0;

            to_0_x_yyzz_yyz[k] = 2.0 * to_yyzz_xyyz[k] * tke_0;

            to_0_x_yyzz_yzz[k] = 2.0 * to_yyzz_xyzz[k] * tke_0;

            to_0_x_yyzz_zzz[k] = 2.0 * to_yyzz_xzzz[k] * tke_0;
        }

        // Set up 130-140 components of targeted buffer : GF

        auto to_0_x_yzzz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 130);

        auto to_0_x_yzzz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 131);

        auto to_0_x_yzzz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 132);

        auto to_0_x_yzzz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 133);

        auto to_0_x_yzzz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 134);

        auto to_0_x_yzzz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 135);

        auto to_0_x_yzzz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 136);

        auto to_0_x_yzzz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 137);

        auto to_0_x_yzzz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 138);

        auto to_0_x_yzzz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_0_x_yzzz_xxx, to_0_x_yzzz_xxy, to_0_x_yzzz_xxz, to_0_x_yzzz_xyy, to_0_x_yzzz_xyz, to_0_x_yzzz_xzz, to_0_x_yzzz_yyy, to_0_x_yzzz_yyz, to_0_x_yzzz_yzz, to_0_x_yzzz_zzz, to_yzzz_xx, to_yzzz_xxxx, to_yzzz_xxxy, to_yzzz_xxxz, to_yzzz_xxyy, to_yzzz_xxyz, to_yzzz_xxzz, to_yzzz_xy, to_yzzz_xyyy, to_yzzz_xyyz, to_yzzz_xyzz, to_yzzz_xz, to_yzzz_xzzz, to_yzzz_yy, to_yzzz_yz, to_yzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yzzz_xxx[k] = -3.0 * to_yzzz_xx[k] + 2.0 * to_yzzz_xxxx[k] * tke_0;

            to_0_x_yzzz_xxy[k] = -2.0 * to_yzzz_xy[k] + 2.0 * to_yzzz_xxxy[k] * tke_0;

            to_0_x_yzzz_xxz[k] = -2.0 * to_yzzz_xz[k] + 2.0 * to_yzzz_xxxz[k] * tke_0;

            to_0_x_yzzz_xyy[k] = -to_yzzz_yy[k] + 2.0 * to_yzzz_xxyy[k] * tke_0;

            to_0_x_yzzz_xyz[k] = -to_yzzz_yz[k] + 2.0 * to_yzzz_xxyz[k] * tke_0;

            to_0_x_yzzz_xzz[k] = -to_yzzz_zz[k] + 2.0 * to_yzzz_xxzz[k] * tke_0;

            to_0_x_yzzz_yyy[k] = 2.0 * to_yzzz_xyyy[k] * tke_0;

            to_0_x_yzzz_yyz[k] = 2.0 * to_yzzz_xyyz[k] * tke_0;

            to_0_x_yzzz_yzz[k] = 2.0 * to_yzzz_xyzz[k] * tke_0;

            to_0_x_yzzz_zzz[k] = 2.0 * to_yzzz_xzzz[k] * tke_0;
        }

        // Set up 140-150 components of targeted buffer : GF

        auto to_0_x_zzzz_xxx = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 140);

        auto to_0_x_zzzz_xxy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 141);

        auto to_0_x_zzzz_xxz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 142);

        auto to_0_x_zzzz_xyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 143);

        auto to_0_x_zzzz_xyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 144);

        auto to_0_x_zzzz_xzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 145);

        auto to_0_x_zzzz_yyy = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 146);

        auto to_0_x_zzzz_yyz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 147);

        auto to_0_x_zzzz_yzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 148);

        auto to_0_x_zzzz_zzz = pbuffer.data(idx_op_geom_001_gf + 0 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_0_x_zzzz_xxx, to_0_x_zzzz_xxy, to_0_x_zzzz_xxz, to_0_x_zzzz_xyy, to_0_x_zzzz_xyz, to_0_x_zzzz_xzz, to_0_x_zzzz_yyy, to_0_x_zzzz_yyz, to_0_x_zzzz_yzz, to_0_x_zzzz_zzz, to_zzzz_xx, to_zzzz_xxxx, to_zzzz_xxxy, to_zzzz_xxxz, to_zzzz_xxyy, to_zzzz_xxyz, to_zzzz_xxzz, to_zzzz_xy, to_zzzz_xyyy, to_zzzz_xyyz, to_zzzz_xyzz, to_zzzz_xz, to_zzzz_xzzz, to_zzzz_yy, to_zzzz_yz, to_zzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zzzz_xxx[k] = -3.0 * to_zzzz_xx[k] + 2.0 * to_zzzz_xxxx[k] * tke_0;

            to_0_x_zzzz_xxy[k] = -2.0 * to_zzzz_xy[k] + 2.0 * to_zzzz_xxxy[k] * tke_0;

            to_0_x_zzzz_xxz[k] = -2.0 * to_zzzz_xz[k] + 2.0 * to_zzzz_xxxz[k] * tke_0;

            to_0_x_zzzz_xyy[k] = -to_zzzz_yy[k] + 2.0 * to_zzzz_xxyy[k] * tke_0;

            to_0_x_zzzz_xyz[k] = -to_zzzz_yz[k] + 2.0 * to_zzzz_xxyz[k] * tke_0;

            to_0_x_zzzz_xzz[k] = -to_zzzz_zz[k] + 2.0 * to_zzzz_xxzz[k] * tke_0;

            to_0_x_zzzz_yyy[k] = 2.0 * to_zzzz_xyyy[k] * tke_0;

            to_0_x_zzzz_yyz[k] = 2.0 * to_zzzz_xyyz[k] * tke_0;

            to_0_x_zzzz_yzz[k] = 2.0 * to_zzzz_xyzz[k] * tke_0;

            to_0_x_zzzz_zzz[k] = 2.0 * to_zzzz_xzzz[k] * tke_0;
        }

        // Set up 150-160 components of targeted buffer : GF

        auto to_0_y_xxxx_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 0);

        auto to_0_y_xxxx_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 1);

        auto to_0_y_xxxx_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 2);

        auto to_0_y_xxxx_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 3);

        auto to_0_y_xxxx_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 4);

        auto to_0_y_xxxx_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 5);

        auto to_0_y_xxxx_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 6);

        auto to_0_y_xxxx_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 7);

        auto to_0_y_xxxx_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 8);

        auto to_0_y_xxxx_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_0_y_xxxx_xxx, to_0_y_xxxx_xxy, to_0_y_xxxx_xxz, to_0_y_xxxx_xyy, to_0_y_xxxx_xyz, to_0_y_xxxx_xzz, to_0_y_xxxx_yyy, to_0_y_xxxx_yyz, to_0_y_xxxx_yzz, to_0_y_xxxx_zzz, to_xxxx_xx, to_xxxx_xxxy, to_xxxx_xxyy, to_xxxx_xxyz, to_xxxx_xy, to_xxxx_xyyy, to_xxxx_xyyz, to_xxxx_xyzz, to_xxxx_xz, to_xxxx_yy, to_xxxx_yyyy, to_xxxx_yyyz, to_xxxx_yyzz, to_xxxx_yz, to_xxxx_yzzz, to_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxx_xxx[k] = 2.0 * to_xxxx_xxxy[k] * tke_0;

            to_0_y_xxxx_xxy[k] = -to_xxxx_xx[k] + 2.0 * to_xxxx_xxyy[k] * tke_0;

            to_0_y_xxxx_xxz[k] = 2.0 * to_xxxx_xxyz[k] * tke_0;

            to_0_y_xxxx_xyy[k] = -2.0 * to_xxxx_xy[k] + 2.0 * to_xxxx_xyyy[k] * tke_0;

            to_0_y_xxxx_xyz[k] = -to_xxxx_xz[k] + 2.0 * to_xxxx_xyyz[k] * tke_0;

            to_0_y_xxxx_xzz[k] = 2.0 * to_xxxx_xyzz[k] * tke_0;

            to_0_y_xxxx_yyy[k] = -3.0 * to_xxxx_yy[k] + 2.0 * to_xxxx_yyyy[k] * tke_0;

            to_0_y_xxxx_yyz[k] = -2.0 * to_xxxx_yz[k] + 2.0 * to_xxxx_yyyz[k] * tke_0;

            to_0_y_xxxx_yzz[k] = -to_xxxx_zz[k] + 2.0 * to_xxxx_yyzz[k] * tke_0;

            to_0_y_xxxx_zzz[k] = 2.0 * to_xxxx_yzzz[k] * tke_0;
        }

        // Set up 160-170 components of targeted buffer : GF

        auto to_0_y_xxxy_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 10);

        auto to_0_y_xxxy_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 11);

        auto to_0_y_xxxy_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 12);

        auto to_0_y_xxxy_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 13);

        auto to_0_y_xxxy_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 14);

        auto to_0_y_xxxy_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 15);

        auto to_0_y_xxxy_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 16);

        auto to_0_y_xxxy_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 17);

        auto to_0_y_xxxy_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 18);

        auto to_0_y_xxxy_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_0_y_xxxy_xxx, to_0_y_xxxy_xxy, to_0_y_xxxy_xxz, to_0_y_xxxy_xyy, to_0_y_xxxy_xyz, to_0_y_xxxy_xzz, to_0_y_xxxy_yyy, to_0_y_xxxy_yyz, to_0_y_xxxy_yzz, to_0_y_xxxy_zzz, to_xxxy_xx, to_xxxy_xxxy, to_xxxy_xxyy, to_xxxy_xxyz, to_xxxy_xy, to_xxxy_xyyy, to_xxxy_xyyz, to_xxxy_xyzz, to_xxxy_xz, to_xxxy_yy, to_xxxy_yyyy, to_xxxy_yyyz, to_xxxy_yyzz, to_xxxy_yz, to_xxxy_yzzz, to_xxxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxy_xxx[k] = 2.0 * to_xxxy_xxxy[k] * tke_0;

            to_0_y_xxxy_xxy[k] = -to_xxxy_xx[k] + 2.0 * to_xxxy_xxyy[k] * tke_0;

            to_0_y_xxxy_xxz[k] = 2.0 * to_xxxy_xxyz[k] * tke_0;

            to_0_y_xxxy_xyy[k] = -2.0 * to_xxxy_xy[k] + 2.0 * to_xxxy_xyyy[k] * tke_0;

            to_0_y_xxxy_xyz[k] = -to_xxxy_xz[k] + 2.0 * to_xxxy_xyyz[k] * tke_0;

            to_0_y_xxxy_xzz[k] = 2.0 * to_xxxy_xyzz[k] * tke_0;

            to_0_y_xxxy_yyy[k] = -3.0 * to_xxxy_yy[k] + 2.0 * to_xxxy_yyyy[k] * tke_0;

            to_0_y_xxxy_yyz[k] = -2.0 * to_xxxy_yz[k] + 2.0 * to_xxxy_yyyz[k] * tke_0;

            to_0_y_xxxy_yzz[k] = -to_xxxy_zz[k] + 2.0 * to_xxxy_yyzz[k] * tke_0;

            to_0_y_xxxy_zzz[k] = 2.0 * to_xxxy_yzzz[k] * tke_0;
        }

        // Set up 170-180 components of targeted buffer : GF

        auto to_0_y_xxxz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 20);

        auto to_0_y_xxxz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 21);

        auto to_0_y_xxxz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 22);

        auto to_0_y_xxxz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 23);

        auto to_0_y_xxxz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 24);

        auto to_0_y_xxxz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 25);

        auto to_0_y_xxxz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 26);

        auto to_0_y_xxxz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 27);

        auto to_0_y_xxxz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 28);

        auto to_0_y_xxxz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_0_y_xxxz_xxx, to_0_y_xxxz_xxy, to_0_y_xxxz_xxz, to_0_y_xxxz_xyy, to_0_y_xxxz_xyz, to_0_y_xxxz_xzz, to_0_y_xxxz_yyy, to_0_y_xxxz_yyz, to_0_y_xxxz_yzz, to_0_y_xxxz_zzz, to_xxxz_xx, to_xxxz_xxxy, to_xxxz_xxyy, to_xxxz_xxyz, to_xxxz_xy, to_xxxz_xyyy, to_xxxz_xyyz, to_xxxz_xyzz, to_xxxz_xz, to_xxxz_yy, to_xxxz_yyyy, to_xxxz_yyyz, to_xxxz_yyzz, to_xxxz_yz, to_xxxz_yzzz, to_xxxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxz_xxx[k] = 2.0 * to_xxxz_xxxy[k] * tke_0;

            to_0_y_xxxz_xxy[k] = -to_xxxz_xx[k] + 2.0 * to_xxxz_xxyy[k] * tke_0;

            to_0_y_xxxz_xxz[k] = 2.0 * to_xxxz_xxyz[k] * tke_0;

            to_0_y_xxxz_xyy[k] = -2.0 * to_xxxz_xy[k] + 2.0 * to_xxxz_xyyy[k] * tke_0;

            to_0_y_xxxz_xyz[k] = -to_xxxz_xz[k] + 2.0 * to_xxxz_xyyz[k] * tke_0;

            to_0_y_xxxz_xzz[k] = 2.0 * to_xxxz_xyzz[k] * tke_0;

            to_0_y_xxxz_yyy[k] = -3.0 * to_xxxz_yy[k] + 2.0 * to_xxxz_yyyy[k] * tke_0;

            to_0_y_xxxz_yyz[k] = -2.0 * to_xxxz_yz[k] + 2.0 * to_xxxz_yyyz[k] * tke_0;

            to_0_y_xxxz_yzz[k] = -to_xxxz_zz[k] + 2.0 * to_xxxz_yyzz[k] * tke_0;

            to_0_y_xxxz_zzz[k] = 2.0 * to_xxxz_yzzz[k] * tke_0;
        }

        // Set up 180-190 components of targeted buffer : GF

        auto to_0_y_xxyy_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 30);

        auto to_0_y_xxyy_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 31);

        auto to_0_y_xxyy_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 32);

        auto to_0_y_xxyy_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 33);

        auto to_0_y_xxyy_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 34);

        auto to_0_y_xxyy_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 35);

        auto to_0_y_xxyy_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 36);

        auto to_0_y_xxyy_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 37);

        auto to_0_y_xxyy_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 38);

        auto to_0_y_xxyy_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_0_y_xxyy_xxx, to_0_y_xxyy_xxy, to_0_y_xxyy_xxz, to_0_y_xxyy_xyy, to_0_y_xxyy_xyz, to_0_y_xxyy_xzz, to_0_y_xxyy_yyy, to_0_y_xxyy_yyz, to_0_y_xxyy_yzz, to_0_y_xxyy_zzz, to_xxyy_xx, to_xxyy_xxxy, to_xxyy_xxyy, to_xxyy_xxyz, to_xxyy_xy, to_xxyy_xyyy, to_xxyy_xyyz, to_xxyy_xyzz, to_xxyy_xz, to_xxyy_yy, to_xxyy_yyyy, to_xxyy_yyyz, to_xxyy_yyzz, to_xxyy_yz, to_xxyy_yzzz, to_xxyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxyy_xxx[k] = 2.0 * to_xxyy_xxxy[k] * tke_0;

            to_0_y_xxyy_xxy[k] = -to_xxyy_xx[k] + 2.0 * to_xxyy_xxyy[k] * tke_0;

            to_0_y_xxyy_xxz[k] = 2.0 * to_xxyy_xxyz[k] * tke_0;

            to_0_y_xxyy_xyy[k] = -2.0 * to_xxyy_xy[k] + 2.0 * to_xxyy_xyyy[k] * tke_0;

            to_0_y_xxyy_xyz[k] = -to_xxyy_xz[k] + 2.0 * to_xxyy_xyyz[k] * tke_0;

            to_0_y_xxyy_xzz[k] = 2.0 * to_xxyy_xyzz[k] * tke_0;

            to_0_y_xxyy_yyy[k] = -3.0 * to_xxyy_yy[k] + 2.0 * to_xxyy_yyyy[k] * tke_0;

            to_0_y_xxyy_yyz[k] = -2.0 * to_xxyy_yz[k] + 2.0 * to_xxyy_yyyz[k] * tke_0;

            to_0_y_xxyy_yzz[k] = -to_xxyy_zz[k] + 2.0 * to_xxyy_yyzz[k] * tke_0;

            to_0_y_xxyy_zzz[k] = 2.0 * to_xxyy_yzzz[k] * tke_0;
        }

        // Set up 190-200 components of targeted buffer : GF

        auto to_0_y_xxyz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 40);

        auto to_0_y_xxyz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 41);

        auto to_0_y_xxyz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 42);

        auto to_0_y_xxyz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 43);

        auto to_0_y_xxyz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 44);

        auto to_0_y_xxyz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 45);

        auto to_0_y_xxyz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 46);

        auto to_0_y_xxyz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 47);

        auto to_0_y_xxyz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 48);

        auto to_0_y_xxyz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_0_y_xxyz_xxx, to_0_y_xxyz_xxy, to_0_y_xxyz_xxz, to_0_y_xxyz_xyy, to_0_y_xxyz_xyz, to_0_y_xxyz_xzz, to_0_y_xxyz_yyy, to_0_y_xxyz_yyz, to_0_y_xxyz_yzz, to_0_y_xxyz_zzz, to_xxyz_xx, to_xxyz_xxxy, to_xxyz_xxyy, to_xxyz_xxyz, to_xxyz_xy, to_xxyz_xyyy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_yy, to_xxyz_yyyy, to_xxyz_yyyz, to_xxyz_yyzz, to_xxyz_yz, to_xxyz_yzzz, to_xxyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxyz_xxx[k] = 2.0 * to_xxyz_xxxy[k] * tke_0;

            to_0_y_xxyz_xxy[k] = -to_xxyz_xx[k] + 2.0 * to_xxyz_xxyy[k] * tke_0;

            to_0_y_xxyz_xxz[k] = 2.0 * to_xxyz_xxyz[k] * tke_0;

            to_0_y_xxyz_xyy[k] = -2.0 * to_xxyz_xy[k] + 2.0 * to_xxyz_xyyy[k] * tke_0;

            to_0_y_xxyz_xyz[k] = -to_xxyz_xz[k] + 2.0 * to_xxyz_xyyz[k] * tke_0;

            to_0_y_xxyz_xzz[k] = 2.0 * to_xxyz_xyzz[k] * tke_0;

            to_0_y_xxyz_yyy[k] = -3.0 * to_xxyz_yy[k] + 2.0 * to_xxyz_yyyy[k] * tke_0;

            to_0_y_xxyz_yyz[k] = -2.0 * to_xxyz_yz[k] + 2.0 * to_xxyz_yyyz[k] * tke_0;

            to_0_y_xxyz_yzz[k] = -to_xxyz_zz[k] + 2.0 * to_xxyz_yyzz[k] * tke_0;

            to_0_y_xxyz_zzz[k] = 2.0 * to_xxyz_yzzz[k] * tke_0;
        }

        // Set up 200-210 components of targeted buffer : GF

        auto to_0_y_xxzz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 50);

        auto to_0_y_xxzz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 51);

        auto to_0_y_xxzz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 52);

        auto to_0_y_xxzz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 53);

        auto to_0_y_xxzz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 54);

        auto to_0_y_xxzz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 55);

        auto to_0_y_xxzz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 56);

        auto to_0_y_xxzz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 57);

        auto to_0_y_xxzz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 58);

        auto to_0_y_xxzz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_0_y_xxzz_xxx, to_0_y_xxzz_xxy, to_0_y_xxzz_xxz, to_0_y_xxzz_xyy, to_0_y_xxzz_xyz, to_0_y_xxzz_xzz, to_0_y_xxzz_yyy, to_0_y_xxzz_yyz, to_0_y_xxzz_yzz, to_0_y_xxzz_zzz, to_xxzz_xx, to_xxzz_xxxy, to_xxzz_xxyy, to_xxzz_xxyz, to_xxzz_xy, to_xxzz_xyyy, to_xxzz_xyyz, to_xxzz_xyzz, to_xxzz_xz, to_xxzz_yy, to_xxzz_yyyy, to_xxzz_yyyz, to_xxzz_yyzz, to_xxzz_yz, to_xxzz_yzzz, to_xxzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxzz_xxx[k] = 2.0 * to_xxzz_xxxy[k] * tke_0;

            to_0_y_xxzz_xxy[k] = -to_xxzz_xx[k] + 2.0 * to_xxzz_xxyy[k] * tke_0;

            to_0_y_xxzz_xxz[k] = 2.0 * to_xxzz_xxyz[k] * tke_0;

            to_0_y_xxzz_xyy[k] = -2.0 * to_xxzz_xy[k] + 2.0 * to_xxzz_xyyy[k] * tke_0;

            to_0_y_xxzz_xyz[k] = -to_xxzz_xz[k] + 2.0 * to_xxzz_xyyz[k] * tke_0;

            to_0_y_xxzz_xzz[k] = 2.0 * to_xxzz_xyzz[k] * tke_0;

            to_0_y_xxzz_yyy[k] = -3.0 * to_xxzz_yy[k] + 2.0 * to_xxzz_yyyy[k] * tke_0;

            to_0_y_xxzz_yyz[k] = -2.0 * to_xxzz_yz[k] + 2.0 * to_xxzz_yyyz[k] * tke_0;

            to_0_y_xxzz_yzz[k] = -to_xxzz_zz[k] + 2.0 * to_xxzz_yyzz[k] * tke_0;

            to_0_y_xxzz_zzz[k] = 2.0 * to_xxzz_yzzz[k] * tke_0;
        }

        // Set up 210-220 components of targeted buffer : GF

        auto to_0_y_xyyy_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 60);

        auto to_0_y_xyyy_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 61);

        auto to_0_y_xyyy_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 62);

        auto to_0_y_xyyy_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 63);

        auto to_0_y_xyyy_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 64);

        auto to_0_y_xyyy_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 65);

        auto to_0_y_xyyy_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 66);

        auto to_0_y_xyyy_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 67);

        auto to_0_y_xyyy_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 68);

        auto to_0_y_xyyy_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_0_y_xyyy_xxx, to_0_y_xyyy_xxy, to_0_y_xyyy_xxz, to_0_y_xyyy_xyy, to_0_y_xyyy_xyz, to_0_y_xyyy_xzz, to_0_y_xyyy_yyy, to_0_y_xyyy_yyz, to_0_y_xyyy_yzz, to_0_y_xyyy_zzz, to_xyyy_xx, to_xyyy_xxxy, to_xyyy_xxyy, to_xyyy_xxyz, to_xyyy_xy, to_xyyy_xyyy, to_xyyy_xyyz, to_xyyy_xyzz, to_xyyy_xz, to_xyyy_yy, to_xyyy_yyyy, to_xyyy_yyyz, to_xyyy_yyzz, to_xyyy_yz, to_xyyy_yzzz, to_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyyy_xxx[k] = 2.0 * to_xyyy_xxxy[k] * tke_0;

            to_0_y_xyyy_xxy[k] = -to_xyyy_xx[k] + 2.0 * to_xyyy_xxyy[k] * tke_0;

            to_0_y_xyyy_xxz[k] = 2.0 * to_xyyy_xxyz[k] * tke_0;

            to_0_y_xyyy_xyy[k] = -2.0 * to_xyyy_xy[k] + 2.0 * to_xyyy_xyyy[k] * tke_0;

            to_0_y_xyyy_xyz[k] = -to_xyyy_xz[k] + 2.0 * to_xyyy_xyyz[k] * tke_0;

            to_0_y_xyyy_xzz[k] = 2.0 * to_xyyy_xyzz[k] * tke_0;

            to_0_y_xyyy_yyy[k] = -3.0 * to_xyyy_yy[k] + 2.0 * to_xyyy_yyyy[k] * tke_0;

            to_0_y_xyyy_yyz[k] = -2.0 * to_xyyy_yz[k] + 2.0 * to_xyyy_yyyz[k] * tke_0;

            to_0_y_xyyy_yzz[k] = -to_xyyy_zz[k] + 2.0 * to_xyyy_yyzz[k] * tke_0;

            to_0_y_xyyy_zzz[k] = 2.0 * to_xyyy_yzzz[k] * tke_0;
        }

        // Set up 220-230 components of targeted buffer : GF

        auto to_0_y_xyyz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 70);

        auto to_0_y_xyyz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 71);

        auto to_0_y_xyyz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 72);

        auto to_0_y_xyyz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 73);

        auto to_0_y_xyyz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 74);

        auto to_0_y_xyyz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 75);

        auto to_0_y_xyyz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 76);

        auto to_0_y_xyyz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 77);

        auto to_0_y_xyyz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 78);

        auto to_0_y_xyyz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_0_y_xyyz_xxx, to_0_y_xyyz_xxy, to_0_y_xyyz_xxz, to_0_y_xyyz_xyy, to_0_y_xyyz_xyz, to_0_y_xyyz_xzz, to_0_y_xyyz_yyy, to_0_y_xyyz_yyz, to_0_y_xyyz_yzz, to_0_y_xyyz_zzz, to_xyyz_xx, to_xyyz_xxxy, to_xyyz_xxyy, to_xyyz_xxyz, to_xyyz_xy, to_xyyz_xyyy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_yy, to_xyyz_yyyy, to_xyyz_yyyz, to_xyyz_yyzz, to_xyyz_yz, to_xyyz_yzzz, to_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyyz_xxx[k] = 2.0 * to_xyyz_xxxy[k] * tke_0;

            to_0_y_xyyz_xxy[k] = -to_xyyz_xx[k] + 2.0 * to_xyyz_xxyy[k] * tke_0;

            to_0_y_xyyz_xxz[k] = 2.0 * to_xyyz_xxyz[k] * tke_0;

            to_0_y_xyyz_xyy[k] = -2.0 * to_xyyz_xy[k] + 2.0 * to_xyyz_xyyy[k] * tke_0;

            to_0_y_xyyz_xyz[k] = -to_xyyz_xz[k] + 2.0 * to_xyyz_xyyz[k] * tke_0;

            to_0_y_xyyz_xzz[k] = 2.0 * to_xyyz_xyzz[k] * tke_0;

            to_0_y_xyyz_yyy[k] = -3.0 * to_xyyz_yy[k] + 2.0 * to_xyyz_yyyy[k] * tke_0;

            to_0_y_xyyz_yyz[k] = -2.0 * to_xyyz_yz[k] + 2.0 * to_xyyz_yyyz[k] * tke_0;

            to_0_y_xyyz_yzz[k] = -to_xyyz_zz[k] + 2.0 * to_xyyz_yyzz[k] * tke_0;

            to_0_y_xyyz_zzz[k] = 2.0 * to_xyyz_yzzz[k] * tke_0;
        }

        // Set up 230-240 components of targeted buffer : GF

        auto to_0_y_xyzz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 80);

        auto to_0_y_xyzz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 81);

        auto to_0_y_xyzz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 82);

        auto to_0_y_xyzz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 83);

        auto to_0_y_xyzz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 84);

        auto to_0_y_xyzz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 85);

        auto to_0_y_xyzz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 86);

        auto to_0_y_xyzz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 87);

        auto to_0_y_xyzz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 88);

        auto to_0_y_xyzz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_0_y_xyzz_xxx, to_0_y_xyzz_xxy, to_0_y_xyzz_xxz, to_0_y_xyzz_xyy, to_0_y_xyzz_xyz, to_0_y_xyzz_xzz, to_0_y_xyzz_yyy, to_0_y_xyzz_yyz, to_0_y_xyzz_yzz, to_0_y_xyzz_zzz, to_xyzz_xx, to_xyzz_xxxy, to_xyzz_xxyy, to_xyzz_xxyz, to_xyzz_xy, to_xyzz_xyyy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_yy, to_xyzz_yyyy, to_xyzz_yyyz, to_xyzz_yyzz, to_xyzz_yz, to_xyzz_yzzz, to_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyzz_xxx[k] = 2.0 * to_xyzz_xxxy[k] * tke_0;

            to_0_y_xyzz_xxy[k] = -to_xyzz_xx[k] + 2.0 * to_xyzz_xxyy[k] * tke_0;

            to_0_y_xyzz_xxz[k] = 2.0 * to_xyzz_xxyz[k] * tke_0;

            to_0_y_xyzz_xyy[k] = -2.0 * to_xyzz_xy[k] + 2.0 * to_xyzz_xyyy[k] * tke_0;

            to_0_y_xyzz_xyz[k] = -to_xyzz_xz[k] + 2.0 * to_xyzz_xyyz[k] * tke_0;

            to_0_y_xyzz_xzz[k] = 2.0 * to_xyzz_xyzz[k] * tke_0;

            to_0_y_xyzz_yyy[k] = -3.0 * to_xyzz_yy[k] + 2.0 * to_xyzz_yyyy[k] * tke_0;

            to_0_y_xyzz_yyz[k] = -2.0 * to_xyzz_yz[k] + 2.0 * to_xyzz_yyyz[k] * tke_0;

            to_0_y_xyzz_yzz[k] = -to_xyzz_zz[k] + 2.0 * to_xyzz_yyzz[k] * tke_0;

            to_0_y_xyzz_zzz[k] = 2.0 * to_xyzz_yzzz[k] * tke_0;
        }

        // Set up 240-250 components of targeted buffer : GF

        auto to_0_y_xzzz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 90);

        auto to_0_y_xzzz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 91);

        auto to_0_y_xzzz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 92);

        auto to_0_y_xzzz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 93);

        auto to_0_y_xzzz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 94);

        auto to_0_y_xzzz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 95);

        auto to_0_y_xzzz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 96);

        auto to_0_y_xzzz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 97);

        auto to_0_y_xzzz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 98);

        auto to_0_y_xzzz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_0_y_xzzz_xxx, to_0_y_xzzz_xxy, to_0_y_xzzz_xxz, to_0_y_xzzz_xyy, to_0_y_xzzz_xyz, to_0_y_xzzz_xzz, to_0_y_xzzz_yyy, to_0_y_xzzz_yyz, to_0_y_xzzz_yzz, to_0_y_xzzz_zzz, to_xzzz_xx, to_xzzz_xxxy, to_xzzz_xxyy, to_xzzz_xxyz, to_xzzz_xy, to_xzzz_xyyy, to_xzzz_xyyz, to_xzzz_xyzz, to_xzzz_xz, to_xzzz_yy, to_xzzz_yyyy, to_xzzz_yyyz, to_xzzz_yyzz, to_xzzz_yz, to_xzzz_yzzz, to_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xzzz_xxx[k] = 2.0 * to_xzzz_xxxy[k] * tke_0;

            to_0_y_xzzz_xxy[k] = -to_xzzz_xx[k] + 2.0 * to_xzzz_xxyy[k] * tke_0;

            to_0_y_xzzz_xxz[k] = 2.0 * to_xzzz_xxyz[k] * tke_0;

            to_0_y_xzzz_xyy[k] = -2.0 * to_xzzz_xy[k] + 2.0 * to_xzzz_xyyy[k] * tke_0;

            to_0_y_xzzz_xyz[k] = -to_xzzz_xz[k] + 2.0 * to_xzzz_xyyz[k] * tke_0;

            to_0_y_xzzz_xzz[k] = 2.0 * to_xzzz_xyzz[k] * tke_0;

            to_0_y_xzzz_yyy[k] = -3.0 * to_xzzz_yy[k] + 2.0 * to_xzzz_yyyy[k] * tke_0;

            to_0_y_xzzz_yyz[k] = -2.0 * to_xzzz_yz[k] + 2.0 * to_xzzz_yyyz[k] * tke_0;

            to_0_y_xzzz_yzz[k] = -to_xzzz_zz[k] + 2.0 * to_xzzz_yyzz[k] * tke_0;

            to_0_y_xzzz_zzz[k] = 2.0 * to_xzzz_yzzz[k] * tke_0;
        }

        // Set up 250-260 components of targeted buffer : GF

        auto to_0_y_yyyy_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 100);

        auto to_0_y_yyyy_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 101);

        auto to_0_y_yyyy_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 102);

        auto to_0_y_yyyy_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 103);

        auto to_0_y_yyyy_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 104);

        auto to_0_y_yyyy_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 105);

        auto to_0_y_yyyy_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 106);

        auto to_0_y_yyyy_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 107);

        auto to_0_y_yyyy_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 108);

        auto to_0_y_yyyy_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_0_y_yyyy_xxx, to_0_y_yyyy_xxy, to_0_y_yyyy_xxz, to_0_y_yyyy_xyy, to_0_y_yyyy_xyz, to_0_y_yyyy_xzz, to_0_y_yyyy_yyy, to_0_y_yyyy_yyz, to_0_y_yyyy_yzz, to_0_y_yyyy_zzz, to_yyyy_xx, to_yyyy_xxxy, to_yyyy_xxyy, to_yyyy_xxyz, to_yyyy_xy, to_yyyy_xyyy, to_yyyy_xyyz, to_yyyy_xyzz, to_yyyy_xz, to_yyyy_yy, to_yyyy_yyyy, to_yyyy_yyyz, to_yyyy_yyzz, to_yyyy_yz, to_yyyy_yzzz, to_yyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyyy_xxx[k] = 2.0 * to_yyyy_xxxy[k] * tke_0;

            to_0_y_yyyy_xxy[k] = -to_yyyy_xx[k] + 2.0 * to_yyyy_xxyy[k] * tke_0;

            to_0_y_yyyy_xxz[k] = 2.0 * to_yyyy_xxyz[k] * tke_0;

            to_0_y_yyyy_xyy[k] = -2.0 * to_yyyy_xy[k] + 2.0 * to_yyyy_xyyy[k] * tke_0;

            to_0_y_yyyy_xyz[k] = -to_yyyy_xz[k] + 2.0 * to_yyyy_xyyz[k] * tke_0;

            to_0_y_yyyy_xzz[k] = 2.0 * to_yyyy_xyzz[k] * tke_0;

            to_0_y_yyyy_yyy[k] = -3.0 * to_yyyy_yy[k] + 2.0 * to_yyyy_yyyy[k] * tke_0;

            to_0_y_yyyy_yyz[k] = -2.0 * to_yyyy_yz[k] + 2.0 * to_yyyy_yyyz[k] * tke_0;

            to_0_y_yyyy_yzz[k] = -to_yyyy_zz[k] + 2.0 * to_yyyy_yyzz[k] * tke_0;

            to_0_y_yyyy_zzz[k] = 2.0 * to_yyyy_yzzz[k] * tke_0;
        }

        // Set up 260-270 components of targeted buffer : GF

        auto to_0_y_yyyz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 110);

        auto to_0_y_yyyz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 111);

        auto to_0_y_yyyz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 112);

        auto to_0_y_yyyz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 113);

        auto to_0_y_yyyz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 114);

        auto to_0_y_yyyz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 115);

        auto to_0_y_yyyz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 116);

        auto to_0_y_yyyz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 117);

        auto to_0_y_yyyz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 118);

        auto to_0_y_yyyz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_0_y_yyyz_xxx, to_0_y_yyyz_xxy, to_0_y_yyyz_xxz, to_0_y_yyyz_xyy, to_0_y_yyyz_xyz, to_0_y_yyyz_xzz, to_0_y_yyyz_yyy, to_0_y_yyyz_yyz, to_0_y_yyyz_yzz, to_0_y_yyyz_zzz, to_yyyz_xx, to_yyyz_xxxy, to_yyyz_xxyy, to_yyyz_xxyz, to_yyyz_xy, to_yyyz_xyyy, to_yyyz_xyyz, to_yyyz_xyzz, to_yyyz_xz, to_yyyz_yy, to_yyyz_yyyy, to_yyyz_yyyz, to_yyyz_yyzz, to_yyyz_yz, to_yyyz_yzzz, to_yyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyyz_xxx[k] = 2.0 * to_yyyz_xxxy[k] * tke_0;

            to_0_y_yyyz_xxy[k] = -to_yyyz_xx[k] + 2.0 * to_yyyz_xxyy[k] * tke_0;

            to_0_y_yyyz_xxz[k] = 2.0 * to_yyyz_xxyz[k] * tke_0;

            to_0_y_yyyz_xyy[k] = -2.0 * to_yyyz_xy[k] + 2.0 * to_yyyz_xyyy[k] * tke_0;

            to_0_y_yyyz_xyz[k] = -to_yyyz_xz[k] + 2.0 * to_yyyz_xyyz[k] * tke_0;

            to_0_y_yyyz_xzz[k] = 2.0 * to_yyyz_xyzz[k] * tke_0;

            to_0_y_yyyz_yyy[k] = -3.0 * to_yyyz_yy[k] + 2.0 * to_yyyz_yyyy[k] * tke_0;

            to_0_y_yyyz_yyz[k] = -2.0 * to_yyyz_yz[k] + 2.0 * to_yyyz_yyyz[k] * tke_0;

            to_0_y_yyyz_yzz[k] = -to_yyyz_zz[k] + 2.0 * to_yyyz_yyzz[k] * tke_0;

            to_0_y_yyyz_zzz[k] = 2.0 * to_yyyz_yzzz[k] * tke_0;
        }

        // Set up 270-280 components of targeted buffer : GF

        auto to_0_y_yyzz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 120);

        auto to_0_y_yyzz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 121);

        auto to_0_y_yyzz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 122);

        auto to_0_y_yyzz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 123);

        auto to_0_y_yyzz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 124);

        auto to_0_y_yyzz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 125);

        auto to_0_y_yyzz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 126);

        auto to_0_y_yyzz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 127);

        auto to_0_y_yyzz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 128);

        auto to_0_y_yyzz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_0_y_yyzz_xxx, to_0_y_yyzz_xxy, to_0_y_yyzz_xxz, to_0_y_yyzz_xyy, to_0_y_yyzz_xyz, to_0_y_yyzz_xzz, to_0_y_yyzz_yyy, to_0_y_yyzz_yyz, to_0_y_yyzz_yzz, to_0_y_yyzz_zzz, to_yyzz_xx, to_yyzz_xxxy, to_yyzz_xxyy, to_yyzz_xxyz, to_yyzz_xy, to_yyzz_xyyy, to_yyzz_xyyz, to_yyzz_xyzz, to_yyzz_xz, to_yyzz_yy, to_yyzz_yyyy, to_yyzz_yyyz, to_yyzz_yyzz, to_yyzz_yz, to_yyzz_yzzz, to_yyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyzz_xxx[k] = 2.0 * to_yyzz_xxxy[k] * tke_0;

            to_0_y_yyzz_xxy[k] = -to_yyzz_xx[k] + 2.0 * to_yyzz_xxyy[k] * tke_0;

            to_0_y_yyzz_xxz[k] = 2.0 * to_yyzz_xxyz[k] * tke_0;

            to_0_y_yyzz_xyy[k] = -2.0 * to_yyzz_xy[k] + 2.0 * to_yyzz_xyyy[k] * tke_0;

            to_0_y_yyzz_xyz[k] = -to_yyzz_xz[k] + 2.0 * to_yyzz_xyyz[k] * tke_0;

            to_0_y_yyzz_xzz[k] = 2.0 * to_yyzz_xyzz[k] * tke_0;

            to_0_y_yyzz_yyy[k] = -3.0 * to_yyzz_yy[k] + 2.0 * to_yyzz_yyyy[k] * tke_0;

            to_0_y_yyzz_yyz[k] = -2.0 * to_yyzz_yz[k] + 2.0 * to_yyzz_yyyz[k] * tke_0;

            to_0_y_yyzz_yzz[k] = -to_yyzz_zz[k] + 2.0 * to_yyzz_yyzz[k] * tke_0;

            to_0_y_yyzz_zzz[k] = 2.0 * to_yyzz_yzzz[k] * tke_0;
        }

        // Set up 280-290 components of targeted buffer : GF

        auto to_0_y_yzzz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 130);

        auto to_0_y_yzzz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 131);

        auto to_0_y_yzzz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 132);

        auto to_0_y_yzzz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 133);

        auto to_0_y_yzzz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 134);

        auto to_0_y_yzzz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 135);

        auto to_0_y_yzzz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 136);

        auto to_0_y_yzzz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 137);

        auto to_0_y_yzzz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 138);

        auto to_0_y_yzzz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_0_y_yzzz_xxx, to_0_y_yzzz_xxy, to_0_y_yzzz_xxz, to_0_y_yzzz_xyy, to_0_y_yzzz_xyz, to_0_y_yzzz_xzz, to_0_y_yzzz_yyy, to_0_y_yzzz_yyz, to_0_y_yzzz_yzz, to_0_y_yzzz_zzz, to_yzzz_xx, to_yzzz_xxxy, to_yzzz_xxyy, to_yzzz_xxyz, to_yzzz_xy, to_yzzz_xyyy, to_yzzz_xyyz, to_yzzz_xyzz, to_yzzz_xz, to_yzzz_yy, to_yzzz_yyyy, to_yzzz_yyyz, to_yzzz_yyzz, to_yzzz_yz, to_yzzz_yzzz, to_yzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yzzz_xxx[k] = 2.0 * to_yzzz_xxxy[k] * tke_0;

            to_0_y_yzzz_xxy[k] = -to_yzzz_xx[k] + 2.0 * to_yzzz_xxyy[k] * tke_0;

            to_0_y_yzzz_xxz[k] = 2.0 * to_yzzz_xxyz[k] * tke_0;

            to_0_y_yzzz_xyy[k] = -2.0 * to_yzzz_xy[k] + 2.0 * to_yzzz_xyyy[k] * tke_0;

            to_0_y_yzzz_xyz[k] = -to_yzzz_xz[k] + 2.0 * to_yzzz_xyyz[k] * tke_0;

            to_0_y_yzzz_xzz[k] = 2.0 * to_yzzz_xyzz[k] * tke_0;

            to_0_y_yzzz_yyy[k] = -3.0 * to_yzzz_yy[k] + 2.0 * to_yzzz_yyyy[k] * tke_0;

            to_0_y_yzzz_yyz[k] = -2.0 * to_yzzz_yz[k] + 2.0 * to_yzzz_yyyz[k] * tke_0;

            to_0_y_yzzz_yzz[k] = -to_yzzz_zz[k] + 2.0 * to_yzzz_yyzz[k] * tke_0;

            to_0_y_yzzz_zzz[k] = 2.0 * to_yzzz_yzzz[k] * tke_0;
        }

        // Set up 290-300 components of targeted buffer : GF

        auto to_0_y_zzzz_xxx = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 140);

        auto to_0_y_zzzz_xxy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 141);

        auto to_0_y_zzzz_xxz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 142);

        auto to_0_y_zzzz_xyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 143);

        auto to_0_y_zzzz_xyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 144);

        auto to_0_y_zzzz_xzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 145);

        auto to_0_y_zzzz_yyy = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 146);

        auto to_0_y_zzzz_yyz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 147);

        auto to_0_y_zzzz_yzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 148);

        auto to_0_y_zzzz_zzz = pbuffer.data(idx_op_geom_001_gf + 1 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_0_y_zzzz_xxx, to_0_y_zzzz_xxy, to_0_y_zzzz_xxz, to_0_y_zzzz_xyy, to_0_y_zzzz_xyz, to_0_y_zzzz_xzz, to_0_y_zzzz_yyy, to_0_y_zzzz_yyz, to_0_y_zzzz_yzz, to_0_y_zzzz_zzz, to_zzzz_xx, to_zzzz_xxxy, to_zzzz_xxyy, to_zzzz_xxyz, to_zzzz_xy, to_zzzz_xyyy, to_zzzz_xyyz, to_zzzz_xyzz, to_zzzz_xz, to_zzzz_yy, to_zzzz_yyyy, to_zzzz_yyyz, to_zzzz_yyzz, to_zzzz_yz, to_zzzz_yzzz, to_zzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zzzz_xxx[k] = 2.0 * to_zzzz_xxxy[k] * tke_0;

            to_0_y_zzzz_xxy[k] = -to_zzzz_xx[k] + 2.0 * to_zzzz_xxyy[k] * tke_0;

            to_0_y_zzzz_xxz[k] = 2.0 * to_zzzz_xxyz[k] * tke_0;

            to_0_y_zzzz_xyy[k] = -2.0 * to_zzzz_xy[k] + 2.0 * to_zzzz_xyyy[k] * tke_0;

            to_0_y_zzzz_xyz[k] = -to_zzzz_xz[k] + 2.0 * to_zzzz_xyyz[k] * tke_0;

            to_0_y_zzzz_xzz[k] = 2.0 * to_zzzz_xyzz[k] * tke_0;

            to_0_y_zzzz_yyy[k] = -3.0 * to_zzzz_yy[k] + 2.0 * to_zzzz_yyyy[k] * tke_0;

            to_0_y_zzzz_yyz[k] = -2.0 * to_zzzz_yz[k] + 2.0 * to_zzzz_yyyz[k] * tke_0;

            to_0_y_zzzz_yzz[k] = -to_zzzz_zz[k] + 2.0 * to_zzzz_yyzz[k] * tke_0;

            to_0_y_zzzz_zzz[k] = 2.0 * to_zzzz_yzzz[k] * tke_0;
        }

        // Set up 300-310 components of targeted buffer : GF

        auto to_0_z_xxxx_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 0);

        auto to_0_z_xxxx_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 1);

        auto to_0_z_xxxx_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 2);

        auto to_0_z_xxxx_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 3);

        auto to_0_z_xxxx_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 4);

        auto to_0_z_xxxx_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 5);

        auto to_0_z_xxxx_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 6);

        auto to_0_z_xxxx_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 7);

        auto to_0_z_xxxx_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 8);

        auto to_0_z_xxxx_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 9);

        #pragma omp simd aligned(to_0_z_xxxx_xxx, to_0_z_xxxx_xxy, to_0_z_xxxx_xxz, to_0_z_xxxx_xyy, to_0_z_xxxx_xyz, to_0_z_xxxx_xzz, to_0_z_xxxx_yyy, to_0_z_xxxx_yyz, to_0_z_xxxx_yzz, to_0_z_xxxx_zzz, to_xxxx_xx, to_xxxx_xxxz, to_xxxx_xxyz, to_xxxx_xxzz, to_xxxx_xy, to_xxxx_xyyz, to_xxxx_xyzz, to_xxxx_xz, to_xxxx_xzzz, to_xxxx_yy, to_xxxx_yyyz, to_xxxx_yyzz, to_xxxx_yz, to_xxxx_yzzz, to_xxxx_zz, to_xxxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxx_xxx[k] = 2.0 * to_xxxx_xxxz[k] * tke_0;

            to_0_z_xxxx_xxy[k] = 2.0 * to_xxxx_xxyz[k] * tke_0;

            to_0_z_xxxx_xxz[k] = -to_xxxx_xx[k] + 2.0 * to_xxxx_xxzz[k] * tke_0;

            to_0_z_xxxx_xyy[k] = 2.0 * to_xxxx_xyyz[k] * tke_0;

            to_0_z_xxxx_xyz[k] = -to_xxxx_xy[k] + 2.0 * to_xxxx_xyzz[k] * tke_0;

            to_0_z_xxxx_xzz[k] = -2.0 * to_xxxx_xz[k] + 2.0 * to_xxxx_xzzz[k] * tke_0;

            to_0_z_xxxx_yyy[k] = 2.0 * to_xxxx_yyyz[k] * tke_0;

            to_0_z_xxxx_yyz[k] = -to_xxxx_yy[k] + 2.0 * to_xxxx_yyzz[k] * tke_0;

            to_0_z_xxxx_yzz[k] = -2.0 * to_xxxx_yz[k] + 2.0 * to_xxxx_yzzz[k] * tke_0;

            to_0_z_xxxx_zzz[k] = -3.0 * to_xxxx_zz[k] + 2.0 * to_xxxx_zzzz[k] * tke_0;
        }

        // Set up 310-320 components of targeted buffer : GF

        auto to_0_z_xxxy_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 10);

        auto to_0_z_xxxy_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 11);

        auto to_0_z_xxxy_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 12);

        auto to_0_z_xxxy_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 13);

        auto to_0_z_xxxy_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 14);

        auto to_0_z_xxxy_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 15);

        auto to_0_z_xxxy_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 16);

        auto to_0_z_xxxy_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 17);

        auto to_0_z_xxxy_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 18);

        auto to_0_z_xxxy_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 19);

        #pragma omp simd aligned(to_0_z_xxxy_xxx, to_0_z_xxxy_xxy, to_0_z_xxxy_xxz, to_0_z_xxxy_xyy, to_0_z_xxxy_xyz, to_0_z_xxxy_xzz, to_0_z_xxxy_yyy, to_0_z_xxxy_yyz, to_0_z_xxxy_yzz, to_0_z_xxxy_zzz, to_xxxy_xx, to_xxxy_xxxz, to_xxxy_xxyz, to_xxxy_xxzz, to_xxxy_xy, to_xxxy_xyyz, to_xxxy_xyzz, to_xxxy_xz, to_xxxy_xzzz, to_xxxy_yy, to_xxxy_yyyz, to_xxxy_yyzz, to_xxxy_yz, to_xxxy_yzzz, to_xxxy_zz, to_xxxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxy_xxx[k] = 2.0 * to_xxxy_xxxz[k] * tke_0;

            to_0_z_xxxy_xxy[k] = 2.0 * to_xxxy_xxyz[k] * tke_0;

            to_0_z_xxxy_xxz[k] = -to_xxxy_xx[k] + 2.0 * to_xxxy_xxzz[k] * tke_0;

            to_0_z_xxxy_xyy[k] = 2.0 * to_xxxy_xyyz[k] * tke_0;

            to_0_z_xxxy_xyz[k] = -to_xxxy_xy[k] + 2.0 * to_xxxy_xyzz[k] * tke_0;

            to_0_z_xxxy_xzz[k] = -2.0 * to_xxxy_xz[k] + 2.0 * to_xxxy_xzzz[k] * tke_0;

            to_0_z_xxxy_yyy[k] = 2.0 * to_xxxy_yyyz[k] * tke_0;

            to_0_z_xxxy_yyz[k] = -to_xxxy_yy[k] + 2.0 * to_xxxy_yyzz[k] * tke_0;

            to_0_z_xxxy_yzz[k] = -2.0 * to_xxxy_yz[k] + 2.0 * to_xxxy_yzzz[k] * tke_0;

            to_0_z_xxxy_zzz[k] = -3.0 * to_xxxy_zz[k] + 2.0 * to_xxxy_zzzz[k] * tke_0;
        }

        // Set up 320-330 components of targeted buffer : GF

        auto to_0_z_xxxz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 20);

        auto to_0_z_xxxz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 21);

        auto to_0_z_xxxz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 22);

        auto to_0_z_xxxz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 23);

        auto to_0_z_xxxz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 24);

        auto to_0_z_xxxz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 25);

        auto to_0_z_xxxz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 26);

        auto to_0_z_xxxz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 27);

        auto to_0_z_xxxz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 28);

        auto to_0_z_xxxz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_0_z_xxxz_xxx, to_0_z_xxxz_xxy, to_0_z_xxxz_xxz, to_0_z_xxxz_xyy, to_0_z_xxxz_xyz, to_0_z_xxxz_xzz, to_0_z_xxxz_yyy, to_0_z_xxxz_yyz, to_0_z_xxxz_yzz, to_0_z_xxxz_zzz, to_xxxz_xx, to_xxxz_xxxz, to_xxxz_xxyz, to_xxxz_xxzz, to_xxxz_xy, to_xxxz_xyyz, to_xxxz_xyzz, to_xxxz_xz, to_xxxz_xzzz, to_xxxz_yy, to_xxxz_yyyz, to_xxxz_yyzz, to_xxxz_yz, to_xxxz_yzzz, to_xxxz_zz, to_xxxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxz_xxx[k] = 2.0 * to_xxxz_xxxz[k] * tke_0;

            to_0_z_xxxz_xxy[k] = 2.0 * to_xxxz_xxyz[k] * tke_0;

            to_0_z_xxxz_xxz[k] = -to_xxxz_xx[k] + 2.0 * to_xxxz_xxzz[k] * tke_0;

            to_0_z_xxxz_xyy[k] = 2.0 * to_xxxz_xyyz[k] * tke_0;

            to_0_z_xxxz_xyz[k] = -to_xxxz_xy[k] + 2.0 * to_xxxz_xyzz[k] * tke_0;

            to_0_z_xxxz_xzz[k] = -2.0 * to_xxxz_xz[k] + 2.0 * to_xxxz_xzzz[k] * tke_0;

            to_0_z_xxxz_yyy[k] = 2.0 * to_xxxz_yyyz[k] * tke_0;

            to_0_z_xxxz_yyz[k] = -to_xxxz_yy[k] + 2.0 * to_xxxz_yyzz[k] * tke_0;

            to_0_z_xxxz_yzz[k] = -2.0 * to_xxxz_yz[k] + 2.0 * to_xxxz_yzzz[k] * tke_0;

            to_0_z_xxxz_zzz[k] = -3.0 * to_xxxz_zz[k] + 2.0 * to_xxxz_zzzz[k] * tke_0;
        }

        // Set up 330-340 components of targeted buffer : GF

        auto to_0_z_xxyy_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 30);

        auto to_0_z_xxyy_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 31);

        auto to_0_z_xxyy_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 32);

        auto to_0_z_xxyy_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 33);

        auto to_0_z_xxyy_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 34);

        auto to_0_z_xxyy_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 35);

        auto to_0_z_xxyy_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 36);

        auto to_0_z_xxyy_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 37);

        auto to_0_z_xxyy_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 38);

        auto to_0_z_xxyy_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 39);

        #pragma omp simd aligned(to_0_z_xxyy_xxx, to_0_z_xxyy_xxy, to_0_z_xxyy_xxz, to_0_z_xxyy_xyy, to_0_z_xxyy_xyz, to_0_z_xxyy_xzz, to_0_z_xxyy_yyy, to_0_z_xxyy_yyz, to_0_z_xxyy_yzz, to_0_z_xxyy_zzz, to_xxyy_xx, to_xxyy_xxxz, to_xxyy_xxyz, to_xxyy_xxzz, to_xxyy_xy, to_xxyy_xyyz, to_xxyy_xyzz, to_xxyy_xz, to_xxyy_xzzz, to_xxyy_yy, to_xxyy_yyyz, to_xxyy_yyzz, to_xxyy_yz, to_xxyy_yzzz, to_xxyy_zz, to_xxyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxyy_xxx[k] = 2.0 * to_xxyy_xxxz[k] * tke_0;

            to_0_z_xxyy_xxy[k] = 2.0 * to_xxyy_xxyz[k] * tke_0;

            to_0_z_xxyy_xxz[k] = -to_xxyy_xx[k] + 2.0 * to_xxyy_xxzz[k] * tke_0;

            to_0_z_xxyy_xyy[k] = 2.0 * to_xxyy_xyyz[k] * tke_0;

            to_0_z_xxyy_xyz[k] = -to_xxyy_xy[k] + 2.0 * to_xxyy_xyzz[k] * tke_0;

            to_0_z_xxyy_xzz[k] = -2.0 * to_xxyy_xz[k] + 2.0 * to_xxyy_xzzz[k] * tke_0;

            to_0_z_xxyy_yyy[k] = 2.0 * to_xxyy_yyyz[k] * tke_0;

            to_0_z_xxyy_yyz[k] = -to_xxyy_yy[k] + 2.0 * to_xxyy_yyzz[k] * tke_0;

            to_0_z_xxyy_yzz[k] = -2.0 * to_xxyy_yz[k] + 2.0 * to_xxyy_yzzz[k] * tke_0;

            to_0_z_xxyy_zzz[k] = -3.0 * to_xxyy_zz[k] + 2.0 * to_xxyy_zzzz[k] * tke_0;
        }

        // Set up 340-350 components of targeted buffer : GF

        auto to_0_z_xxyz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 40);

        auto to_0_z_xxyz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 41);

        auto to_0_z_xxyz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 42);

        auto to_0_z_xxyz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 43);

        auto to_0_z_xxyz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 44);

        auto to_0_z_xxyz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 45);

        auto to_0_z_xxyz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 46);

        auto to_0_z_xxyz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 47);

        auto to_0_z_xxyz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 48);

        auto to_0_z_xxyz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 49);

        #pragma omp simd aligned(to_0_z_xxyz_xxx, to_0_z_xxyz_xxy, to_0_z_xxyz_xxz, to_0_z_xxyz_xyy, to_0_z_xxyz_xyz, to_0_z_xxyz_xzz, to_0_z_xxyz_yyy, to_0_z_xxyz_yyz, to_0_z_xxyz_yzz, to_0_z_xxyz_zzz, to_xxyz_xx, to_xxyz_xxxz, to_xxyz_xxyz, to_xxyz_xxzz, to_xxyz_xy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_xzzz, to_xxyz_yy, to_xxyz_yyyz, to_xxyz_yyzz, to_xxyz_yz, to_xxyz_yzzz, to_xxyz_zz, to_xxyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxyz_xxx[k] = 2.0 * to_xxyz_xxxz[k] * tke_0;

            to_0_z_xxyz_xxy[k] = 2.0 * to_xxyz_xxyz[k] * tke_0;

            to_0_z_xxyz_xxz[k] = -to_xxyz_xx[k] + 2.0 * to_xxyz_xxzz[k] * tke_0;

            to_0_z_xxyz_xyy[k] = 2.0 * to_xxyz_xyyz[k] * tke_0;

            to_0_z_xxyz_xyz[k] = -to_xxyz_xy[k] + 2.0 * to_xxyz_xyzz[k] * tke_0;

            to_0_z_xxyz_xzz[k] = -2.0 * to_xxyz_xz[k] + 2.0 * to_xxyz_xzzz[k] * tke_0;

            to_0_z_xxyz_yyy[k] = 2.0 * to_xxyz_yyyz[k] * tke_0;

            to_0_z_xxyz_yyz[k] = -to_xxyz_yy[k] + 2.0 * to_xxyz_yyzz[k] * tke_0;

            to_0_z_xxyz_yzz[k] = -2.0 * to_xxyz_yz[k] + 2.0 * to_xxyz_yzzz[k] * tke_0;

            to_0_z_xxyz_zzz[k] = -3.0 * to_xxyz_zz[k] + 2.0 * to_xxyz_zzzz[k] * tke_0;
        }

        // Set up 350-360 components of targeted buffer : GF

        auto to_0_z_xxzz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 50);

        auto to_0_z_xxzz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 51);

        auto to_0_z_xxzz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 52);

        auto to_0_z_xxzz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 53);

        auto to_0_z_xxzz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 54);

        auto to_0_z_xxzz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 55);

        auto to_0_z_xxzz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 56);

        auto to_0_z_xxzz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 57);

        auto to_0_z_xxzz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 58);

        auto to_0_z_xxzz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_0_z_xxzz_xxx, to_0_z_xxzz_xxy, to_0_z_xxzz_xxz, to_0_z_xxzz_xyy, to_0_z_xxzz_xyz, to_0_z_xxzz_xzz, to_0_z_xxzz_yyy, to_0_z_xxzz_yyz, to_0_z_xxzz_yzz, to_0_z_xxzz_zzz, to_xxzz_xx, to_xxzz_xxxz, to_xxzz_xxyz, to_xxzz_xxzz, to_xxzz_xy, to_xxzz_xyyz, to_xxzz_xyzz, to_xxzz_xz, to_xxzz_xzzz, to_xxzz_yy, to_xxzz_yyyz, to_xxzz_yyzz, to_xxzz_yz, to_xxzz_yzzz, to_xxzz_zz, to_xxzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxzz_xxx[k] = 2.0 * to_xxzz_xxxz[k] * tke_0;

            to_0_z_xxzz_xxy[k] = 2.0 * to_xxzz_xxyz[k] * tke_0;

            to_0_z_xxzz_xxz[k] = -to_xxzz_xx[k] + 2.0 * to_xxzz_xxzz[k] * tke_0;

            to_0_z_xxzz_xyy[k] = 2.0 * to_xxzz_xyyz[k] * tke_0;

            to_0_z_xxzz_xyz[k] = -to_xxzz_xy[k] + 2.0 * to_xxzz_xyzz[k] * tke_0;

            to_0_z_xxzz_xzz[k] = -2.0 * to_xxzz_xz[k] + 2.0 * to_xxzz_xzzz[k] * tke_0;

            to_0_z_xxzz_yyy[k] = 2.0 * to_xxzz_yyyz[k] * tke_0;

            to_0_z_xxzz_yyz[k] = -to_xxzz_yy[k] + 2.0 * to_xxzz_yyzz[k] * tke_0;

            to_0_z_xxzz_yzz[k] = -2.0 * to_xxzz_yz[k] + 2.0 * to_xxzz_yzzz[k] * tke_0;

            to_0_z_xxzz_zzz[k] = -3.0 * to_xxzz_zz[k] + 2.0 * to_xxzz_zzzz[k] * tke_0;
        }

        // Set up 360-370 components of targeted buffer : GF

        auto to_0_z_xyyy_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 60);

        auto to_0_z_xyyy_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 61);

        auto to_0_z_xyyy_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 62);

        auto to_0_z_xyyy_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 63);

        auto to_0_z_xyyy_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 64);

        auto to_0_z_xyyy_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 65);

        auto to_0_z_xyyy_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 66);

        auto to_0_z_xyyy_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 67);

        auto to_0_z_xyyy_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 68);

        auto to_0_z_xyyy_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 69);

        #pragma omp simd aligned(to_0_z_xyyy_xxx, to_0_z_xyyy_xxy, to_0_z_xyyy_xxz, to_0_z_xyyy_xyy, to_0_z_xyyy_xyz, to_0_z_xyyy_xzz, to_0_z_xyyy_yyy, to_0_z_xyyy_yyz, to_0_z_xyyy_yzz, to_0_z_xyyy_zzz, to_xyyy_xx, to_xyyy_xxxz, to_xyyy_xxyz, to_xyyy_xxzz, to_xyyy_xy, to_xyyy_xyyz, to_xyyy_xyzz, to_xyyy_xz, to_xyyy_xzzz, to_xyyy_yy, to_xyyy_yyyz, to_xyyy_yyzz, to_xyyy_yz, to_xyyy_yzzz, to_xyyy_zz, to_xyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyyy_xxx[k] = 2.0 * to_xyyy_xxxz[k] * tke_0;

            to_0_z_xyyy_xxy[k] = 2.0 * to_xyyy_xxyz[k] * tke_0;

            to_0_z_xyyy_xxz[k] = -to_xyyy_xx[k] + 2.0 * to_xyyy_xxzz[k] * tke_0;

            to_0_z_xyyy_xyy[k] = 2.0 * to_xyyy_xyyz[k] * tke_0;

            to_0_z_xyyy_xyz[k] = -to_xyyy_xy[k] + 2.0 * to_xyyy_xyzz[k] * tke_0;

            to_0_z_xyyy_xzz[k] = -2.0 * to_xyyy_xz[k] + 2.0 * to_xyyy_xzzz[k] * tke_0;

            to_0_z_xyyy_yyy[k] = 2.0 * to_xyyy_yyyz[k] * tke_0;

            to_0_z_xyyy_yyz[k] = -to_xyyy_yy[k] + 2.0 * to_xyyy_yyzz[k] * tke_0;

            to_0_z_xyyy_yzz[k] = -2.0 * to_xyyy_yz[k] + 2.0 * to_xyyy_yzzz[k] * tke_0;

            to_0_z_xyyy_zzz[k] = -3.0 * to_xyyy_zz[k] + 2.0 * to_xyyy_zzzz[k] * tke_0;
        }

        // Set up 370-380 components of targeted buffer : GF

        auto to_0_z_xyyz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 70);

        auto to_0_z_xyyz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 71);

        auto to_0_z_xyyz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 72);

        auto to_0_z_xyyz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 73);

        auto to_0_z_xyyz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 74);

        auto to_0_z_xyyz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 75);

        auto to_0_z_xyyz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 76);

        auto to_0_z_xyyz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 77);

        auto to_0_z_xyyz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 78);

        auto to_0_z_xyyz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 79);

        #pragma omp simd aligned(to_0_z_xyyz_xxx, to_0_z_xyyz_xxy, to_0_z_xyyz_xxz, to_0_z_xyyz_xyy, to_0_z_xyyz_xyz, to_0_z_xyyz_xzz, to_0_z_xyyz_yyy, to_0_z_xyyz_yyz, to_0_z_xyyz_yzz, to_0_z_xyyz_zzz, to_xyyz_xx, to_xyyz_xxxz, to_xyyz_xxyz, to_xyyz_xxzz, to_xyyz_xy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_xzzz, to_xyyz_yy, to_xyyz_yyyz, to_xyyz_yyzz, to_xyyz_yz, to_xyyz_yzzz, to_xyyz_zz, to_xyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyyz_xxx[k] = 2.0 * to_xyyz_xxxz[k] * tke_0;

            to_0_z_xyyz_xxy[k] = 2.0 * to_xyyz_xxyz[k] * tke_0;

            to_0_z_xyyz_xxz[k] = -to_xyyz_xx[k] + 2.0 * to_xyyz_xxzz[k] * tke_0;

            to_0_z_xyyz_xyy[k] = 2.0 * to_xyyz_xyyz[k] * tke_0;

            to_0_z_xyyz_xyz[k] = -to_xyyz_xy[k] + 2.0 * to_xyyz_xyzz[k] * tke_0;

            to_0_z_xyyz_xzz[k] = -2.0 * to_xyyz_xz[k] + 2.0 * to_xyyz_xzzz[k] * tke_0;

            to_0_z_xyyz_yyy[k] = 2.0 * to_xyyz_yyyz[k] * tke_0;

            to_0_z_xyyz_yyz[k] = -to_xyyz_yy[k] + 2.0 * to_xyyz_yyzz[k] * tke_0;

            to_0_z_xyyz_yzz[k] = -2.0 * to_xyyz_yz[k] + 2.0 * to_xyyz_yzzz[k] * tke_0;

            to_0_z_xyyz_zzz[k] = -3.0 * to_xyyz_zz[k] + 2.0 * to_xyyz_zzzz[k] * tke_0;
        }

        // Set up 380-390 components of targeted buffer : GF

        auto to_0_z_xyzz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 80);

        auto to_0_z_xyzz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 81);

        auto to_0_z_xyzz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 82);

        auto to_0_z_xyzz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 83);

        auto to_0_z_xyzz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 84);

        auto to_0_z_xyzz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 85);

        auto to_0_z_xyzz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 86);

        auto to_0_z_xyzz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 87);

        auto to_0_z_xyzz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 88);

        auto to_0_z_xyzz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_0_z_xyzz_xxx, to_0_z_xyzz_xxy, to_0_z_xyzz_xxz, to_0_z_xyzz_xyy, to_0_z_xyzz_xyz, to_0_z_xyzz_xzz, to_0_z_xyzz_yyy, to_0_z_xyzz_yyz, to_0_z_xyzz_yzz, to_0_z_xyzz_zzz, to_xyzz_xx, to_xyzz_xxxz, to_xyzz_xxyz, to_xyzz_xxzz, to_xyzz_xy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_xzzz, to_xyzz_yy, to_xyzz_yyyz, to_xyzz_yyzz, to_xyzz_yz, to_xyzz_yzzz, to_xyzz_zz, to_xyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyzz_xxx[k] = 2.0 * to_xyzz_xxxz[k] * tke_0;

            to_0_z_xyzz_xxy[k] = 2.0 * to_xyzz_xxyz[k] * tke_0;

            to_0_z_xyzz_xxz[k] = -to_xyzz_xx[k] + 2.0 * to_xyzz_xxzz[k] * tke_0;

            to_0_z_xyzz_xyy[k] = 2.0 * to_xyzz_xyyz[k] * tke_0;

            to_0_z_xyzz_xyz[k] = -to_xyzz_xy[k] + 2.0 * to_xyzz_xyzz[k] * tke_0;

            to_0_z_xyzz_xzz[k] = -2.0 * to_xyzz_xz[k] + 2.0 * to_xyzz_xzzz[k] * tke_0;

            to_0_z_xyzz_yyy[k] = 2.0 * to_xyzz_yyyz[k] * tke_0;

            to_0_z_xyzz_yyz[k] = -to_xyzz_yy[k] + 2.0 * to_xyzz_yyzz[k] * tke_0;

            to_0_z_xyzz_yzz[k] = -2.0 * to_xyzz_yz[k] + 2.0 * to_xyzz_yzzz[k] * tke_0;

            to_0_z_xyzz_zzz[k] = -3.0 * to_xyzz_zz[k] + 2.0 * to_xyzz_zzzz[k] * tke_0;
        }

        // Set up 390-400 components of targeted buffer : GF

        auto to_0_z_xzzz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 90);

        auto to_0_z_xzzz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 91);

        auto to_0_z_xzzz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 92);

        auto to_0_z_xzzz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 93);

        auto to_0_z_xzzz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 94);

        auto to_0_z_xzzz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 95);

        auto to_0_z_xzzz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 96);

        auto to_0_z_xzzz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 97);

        auto to_0_z_xzzz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 98);

        auto to_0_z_xzzz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 99);

        #pragma omp simd aligned(to_0_z_xzzz_xxx, to_0_z_xzzz_xxy, to_0_z_xzzz_xxz, to_0_z_xzzz_xyy, to_0_z_xzzz_xyz, to_0_z_xzzz_xzz, to_0_z_xzzz_yyy, to_0_z_xzzz_yyz, to_0_z_xzzz_yzz, to_0_z_xzzz_zzz, to_xzzz_xx, to_xzzz_xxxz, to_xzzz_xxyz, to_xzzz_xxzz, to_xzzz_xy, to_xzzz_xyyz, to_xzzz_xyzz, to_xzzz_xz, to_xzzz_xzzz, to_xzzz_yy, to_xzzz_yyyz, to_xzzz_yyzz, to_xzzz_yz, to_xzzz_yzzz, to_xzzz_zz, to_xzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xzzz_xxx[k] = 2.0 * to_xzzz_xxxz[k] * tke_0;

            to_0_z_xzzz_xxy[k] = 2.0 * to_xzzz_xxyz[k] * tke_0;

            to_0_z_xzzz_xxz[k] = -to_xzzz_xx[k] + 2.0 * to_xzzz_xxzz[k] * tke_0;

            to_0_z_xzzz_xyy[k] = 2.0 * to_xzzz_xyyz[k] * tke_0;

            to_0_z_xzzz_xyz[k] = -to_xzzz_xy[k] + 2.0 * to_xzzz_xyzz[k] * tke_0;

            to_0_z_xzzz_xzz[k] = -2.0 * to_xzzz_xz[k] + 2.0 * to_xzzz_xzzz[k] * tke_0;

            to_0_z_xzzz_yyy[k] = 2.0 * to_xzzz_yyyz[k] * tke_0;

            to_0_z_xzzz_yyz[k] = -to_xzzz_yy[k] + 2.0 * to_xzzz_yyzz[k] * tke_0;

            to_0_z_xzzz_yzz[k] = -2.0 * to_xzzz_yz[k] + 2.0 * to_xzzz_yzzz[k] * tke_0;

            to_0_z_xzzz_zzz[k] = -3.0 * to_xzzz_zz[k] + 2.0 * to_xzzz_zzzz[k] * tke_0;
        }

        // Set up 400-410 components of targeted buffer : GF

        auto to_0_z_yyyy_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 100);

        auto to_0_z_yyyy_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 101);

        auto to_0_z_yyyy_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 102);

        auto to_0_z_yyyy_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 103);

        auto to_0_z_yyyy_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 104);

        auto to_0_z_yyyy_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 105);

        auto to_0_z_yyyy_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 106);

        auto to_0_z_yyyy_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 107);

        auto to_0_z_yyyy_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 108);

        auto to_0_z_yyyy_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 109);

        #pragma omp simd aligned(to_0_z_yyyy_xxx, to_0_z_yyyy_xxy, to_0_z_yyyy_xxz, to_0_z_yyyy_xyy, to_0_z_yyyy_xyz, to_0_z_yyyy_xzz, to_0_z_yyyy_yyy, to_0_z_yyyy_yyz, to_0_z_yyyy_yzz, to_0_z_yyyy_zzz, to_yyyy_xx, to_yyyy_xxxz, to_yyyy_xxyz, to_yyyy_xxzz, to_yyyy_xy, to_yyyy_xyyz, to_yyyy_xyzz, to_yyyy_xz, to_yyyy_xzzz, to_yyyy_yy, to_yyyy_yyyz, to_yyyy_yyzz, to_yyyy_yz, to_yyyy_yzzz, to_yyyy_zz, to_yyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyyy_xxx[k] = 2.0 * to_yyyy_xxxz[k] * tke_0;

            to_0_z_yyyy_xxy[k] = 2.0 * to_yyyy_xxyz[k] * tke_0;

            to_0_z_yyyy_xxz[k] = -to_yyyy_xx[k] + 2.0 * to_yyyy_xxzz[k] * tke_0;

            to_0_z_yyyy_xyy[k] = 2.0 * to_yyyy_xyyz[k] * tke_0;

            to_0_z_yyyy_xyz[k] = -to_yyyy_xy[k] + 2.0 * to_yyyy_xyzz[k] * tke_0;

            to_0_z_yyyy_xzz[k] = -2.0 * to_yyyy_xz[k] + 2.0 * to_yyyy_xzzz[k] * tke_0;

            to_0_z_yyyy_yyy[k] = 2.0 * to_yyyy_yyyz[k] * tke_0;

            to_0_z_yyyy_yyz[k] = -to_yyyy_yy[k] + 2.0 * to_yyyy_yyzz[k] * tke_0;

            to_0_z_yyyy_yzz[k] = -2.0 * to_yyyy_yz[k] + 2.0 * to_yyyy_yzzz[k] * tke_0;

            to_0_z_yyyy_zzz[k] = -3.0 * to_yyyy_zz[k] + 2.0 * to_yyyy_zzzz[k] * tke_0;
        }

        // Set up 410-420 components of targeted buffer : GF

        auto to_0_z_yyyz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 110);

        auto to_0_z_yyyz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 111);

        auto to_0_z_yyyz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 112);

        auto to_0_z_yyyz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 113);

        auto to_0_z_yyyz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 114);

        auto to_0_z_yyyz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 115);

        auto to_0_z_yyyz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 116);

        auto to_0_z_yyyz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 117);

        auto to_0_z_yyyz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 118);

        auto to_0_z_yyyz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_0_z_yyyz_xxx, to_0_z_yyyz_xxy, to_0_z_yyyz_xxz, to_0_z_yyyz_xyy, to_0_z_yyyz_xyz, to_0_z_yyyz_xzz, to_0_z_yyyz_yyy, to_0_z_yyyz_yyz, to_0_z_yyyz_yzz, to_0_z_yyyz_zzz, to_yyyz_xx, to_yyyz_xxxz, to_yyyz_xxyz, to_yyyz_xxzz, to_yyyz_xy, to_yyyz_xyyz, to_yyyz_xyzz, to_yyyz_xz, to_yyyz_xzzz, to_yyyz_yy, to_yyyz_yyyz, to_yyyz_yyzz, to_yyyz_yz, to_yyyz_yzzz, to_yyyz_zz, to_yyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyyz_xxx[k] = 2.0 * to_yyyz_xxxz[k] * tke_0;

            to_0_z_yyyz_xxy[k] = 2.0 * to_yyyz_xxyz[k] * tke_0;

            to_0_z_yyyz_xxz[k] = -to_yyyz_xx[k] + 2.0 * to_yyyz_xxzz[k] * tke_0;

            to_0_z_yyyz_xyy[k] = 2.0 * to_yyyz_xyyz[k] * tke_0;

            to_0_z_yyyz_xyz[k] = -to_yyyz_xy[k] + 2.0 * to_yyyz_xyzz[k] * tke_0;

            to_0_z_yyyz_xzz[k] = -2.0 * to_yyyz_xz[k] + 2.0 * to_yyyz_xzzz[k] * tke_0;

            to_0_z_yyyz_yyy[k] = 2.0 * to_yyyz_yyyz[k] * tke_0;

            to_0_z_yyyz_yyz[k] = -to_yyyz_yy[k] + 2.0 * to_yyyz_yyzz[k] * tke_0;

            to_0_z_yyyz_yzz[k] = -2.0 * to_yyyz_yz[k] + 2.0 * to_yyyz_yzzz[k] * tke_0;

            to_0_z_yyyz_zzz[k] = -3.0 * to_yyyz_zz[k] + 2.0 * to_yyyz_zzzz[k] * tke_0;
        }

        // Set up 420-430 components of targeted buffer : GF

        auto to_0_z_yyzz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 120);

        auto to_0_z_yyzz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 121);

        auto to_0_z_yyzz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 122);

        auto to_0_z_yyzz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 123);

        auto to_0_z_yyzz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 124);

        auto to_0_z_yyzz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 125);

        auto to_0_z_yyzz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 126);

        auto to_0_z_yyzz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 127);

        auto to_0_z_yyzz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 128);

        auto to_0_z_yyzz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 129);

        #pragma omp simd aligned(to_0_z_yyzz_xxx, to_0_z_yyzz_xxy, to_0_z_yyzz_xxz, to_0_z_yyzz_xyy, to_0_z_yyzz_xyz, to_0_z_yyzz_xzz, to_0_z_yyzz_yyy, to_0_z_yyzz_yyz, to_0_z_yyzz_yzz, to_0_z_yyzz_zzz, to_yyzz_xx, to_yyzz_xxxz, to_yyzz_xxyz, to_yyzz_xxzz, to_yyzz_xy, to_yyzz_xyyz, to_yyzz_xyzz, to_yyzz_xz, to_yyzz_xzzz, to_yyzz_yy, to_yyzz_yyyz, to_yyzz_yyzz, to_yyzz_yz, to_yyzz_yzzz, to_yyzz_zz, to_yyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyzz_xxx[k] = 2.0 * to_yyzz_xxxz[k] * tke_0;

            to_0_z_yyzz_xxy[k] = 2.0 * to_yyzz_xxyz[k] * tke_0;

            to_0_z_yyzz_xxz[k] = -to_yyzz_xx[k] + 2.0 * to_yyzz_xxzz[k] * tke_0;

            to_0_z_yyzz_xyy[k] = 2.0 * to_yyzz_xyyz[k] * tke_0;

            to_0_z_yyzz_xyz[k] = -to_yyzz_xy[k] + 2.0 * to_yyzz_xyzz[k] * tke_0;

            to_0_z_yyzz_xzz[k] = -2.0 * to_yyzz_xz[k] + 2.0 * to_yyzz_xzzz[k] * tke_0;

            to_0_z_yyzz_yyy[k] = 2.0 * to_yyzz_yyyz[k] * tke_0;

            to_0_z_yyzz_yyz[k] = -to_yyzz_yy[k] + 2.0 * to_yyzz_yyzz[k] * tke_0;

            to_0_z_yyzz_yzz[k] = -2.0 * to_yyzz_yz[k] + 2.0 * to_yyzz_yzzz[k] * tke_0;

            to_0_z_yyzz_zzz[k] = -3.0 * to_yyzz_zz[k] + 2.0 * to_yyzz_zzzz[k] * tke_0;
        }

        // Set up 430-440 components of targeted buffer : GF

        auto to_0_z_yzzz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 130);

        auto to_0_z_yzzz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 131);

        auto to_0_z_yzzz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 132);

        auto to_0_z_yzzz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 133);

        auto to_0_z_yzzz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 134);

        auto to_0_z_yzzz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 135);

        auto to_0_z_yzzz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 136);

        auto to_0_z_yzzz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 137);

        auto to_0_z_yzzz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 138);

        auto to_0_z_yzzz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 139);

        #pragma omp simd aligned(to_0_z_yzzz_xxx, to_0_z_yzzz_xxy, to_0_z_yzzz_xxz, to_0_z_yzzz_xyy, to_0_z_yzzz_xyz, to_0_z_yzzz_xzz, to_0_z_yzzz_yyy, to_0_z_yzzz_yyz, to_0_z_yzzz_yzz, to_0_z_yzzz_zzz, to_yzzz_xx, to_yzzz_xxxz, to_yzzz_xxyz, to_yzzz_xxzz, to_yzzz_xy, to_yzzz_xyyz, to_yzzz_xyzz, to_yzzz_xz, to_yzzz_xzzz, to_yzzz_yy, to_yzzz_yyyz, to_yzzz_yyzz, to_yzzz_yz, to_yzzz_yzzz, to_yzzz_zz, to_yzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yzzz_xxx[k] = 2.0 * to_yzzz_xxxz[k] * tke_0;

            to_0_z_yzzz_xxy[k] = 2.0 * to_yzzz_xxyz[k] * tke_0;

            to_0_z_yzzz_xxz[k] = -to_yzzz_xx[k] + 2.0 * to_yzzz_xxzz[k] * tke_0;

            to_0_z_yzzz_xyy[k] = 2.0 * to_yzzz_xyyz[k] * tke_0;

            to_0_z_yzzz_xyz[k] = -to_yzzz_xy[k] + 2.0 * to_yzzz_xyzz[k] * tke_0;

            to_0_z_yzzz_xzz[k] = -2.0 * to_yzzz_xz[k] + 2.0 * to_yzzz_xzzz[k] * tke_0;

            to_0_z_yzzz_yyy[k] = 2.0 * to_yzzz_yyyz[k] * tke_0;

            to_0_z_yzzz_yyz[k] = -to_yzzz_yy[k] + 2.0 * to_yzzz_yyzz[k] * tke_0;

            to_0_z_yzzz_yzz[k] = -2.0 * to_yzzz_yz[k] + 2.0 * to_yzzz_yzzz[k] * tke_0;

            to_0_z_yzzz_zzz[k] = -3.0 * to_yzzz_zz[k] + 2.0 * to_yzzz_zzzz[k] * tke_0;
        }

        // Set up 440-450 components of targeted buffer : GF

        auto to_0_z_zzzz_xxx = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 140);

        auto to_0_z_zzzz_xxy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 141);

        auto to_0_z_zzzz_xxz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 142);

        auto to_0_z_zzzz_xyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 143);

        auto to_0_z_zzzz_xyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 144);

        auto to_0_z_zzzz_xzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 145);

        auto to_0_z_zzzz_yyy = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 146);

        auto to_0_z_zzzz_yyz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 147);

        auto to_0_z_zzzz_yzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 148);

        auto to_0_z_zzzz_zzz = pbuffer.data(idx_op_geom_001_gf + 2 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_0_z_zzzz_xxx, to_0_z_zzzz_xxy, to_0_z_zzzz_xxz, to_0_z_zzzz_xyy, to_0_z_zzzz_xyz, to_0_z_zzzz_xzz, to_0_z_zzzz_yyy, to_0_z_zzzz_yyz, to_0_z_zzzz_yzz, to_0_z_zzzz_zzz, to_zzzz_xx, to_zzzz_xxxz, to_zzzz_xxyz, to_zzzz_xxzz, to_zzzz_xy, to_zzzz_xyyz, to_zzzz_xyzz, to_zzzz_xz, to_zzzz_xzzz, to_zzzz_yy, to_zzzz_yyyz, to_zzzz_yyzz, to_zzzz_yz, to_zzzz_yzzz, to_zzzz_zz, to_zzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zzzz_xxx[k] = 2.0 * to_zzzz_xxxz[k] * tke_0;

            to_0_z_zzzz_xxy[k] = 2.0 * to_zzzz_xxyz[k] * tke_0;

            to_0_z_zzzz_xxz[k] = -to_zzzz_xx[k] + 2.0 * to_zzzz_xxzz[k] * tke_0;

            to_0_z_zzzz_xyy[k] = 2.0 * to_zzzz_xyyz[k] * tke_0;

            to_0_z_zzzz_xyz[k] = -to_zzzz_xy[k] + 2.0 * to_zzzz_xyzz[k] * tke_0;

            to_0_z_zzzz_xzz[k] = -2.0 * to_zzzz_xz[k] + 2.0 * to_zzzz_xzzz[k] * tke_0;

            to_0_z_zzzz_yyy[k] = 2.0 * to_zzzz_yyyz[k] * tke_0;

            to_0_z_zzzz_yyz[k] = -to_zzzz_yy[k] + 2.0 * to_zzzz_yyzz[k] * tke_0;

            to_0_z_zzzz_yzz[k] = -2.0 * to_zzzz_yz[k] + 2.0 * to_zzzz_yzzz[k] * tke_0;

            to_0_z_zzzz_zzz[k] = -3.0 * to_zzzz_zz[k] + 2.0 * to_zzzz_zzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

