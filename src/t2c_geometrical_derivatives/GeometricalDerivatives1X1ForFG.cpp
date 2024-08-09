#include "GeometricalDerivatives1X1ForFG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_fg(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_fg,
                        const size_t idx_op_df,
                        const size_t idx_op_dh,
                        const size_t idx_op_gf,
                        const size_t idx_op_gh,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : DF

        auto to_xx_xxx = pbuffer.data(idx_op_df + i * 60 + 0);

        auto to_xx_xxy = pbuffer.data(idx_op_df + i * 60 + 1);

        auto to_xx_xxz = pbuffer.data(idx_op_df + i * 60 + 2);

        auto to_xx_xyy = pbuffer.data(idx_op_df + i * 60 + 3);

        auto to_xx_xyz = pbuffer.data(idx_op_df + i * 60 + 4);

        auto to_xx_xzz = pbuffer.data(idx_op_df + i * 60 + 5);

        auto to_xx_yyy = pbuffer.data(idx_op_df + i * 60 + 6);

        auto to_xx_yyz = pbuffer.data(idx_op_df + i * 60 + 7);

        auto to_xx_yzz = pbuffer.data(idx_op_df + i * 60 + 8);

        auto to_xx_zzz = pbuffer.data(idx_op_df + i * 60 + 9);

        auto to_xy_xxx = pbuffer.data(idx_op_df + i * 60 + 10);

        auto to_xy_xxy = pbuffer.data(idx_op_df + i * 60 + 11);

        auto to_xy_xxz = pbuffer.data(idx_op_df + i * 60 + 12);

        auto to_xy_xyy = pbuffer.data(idx_op_df + i * 60 + 13);

        auto to_xy_xyz = pbuffer.data(idx_op_df + i * 60 + 14);

        auto to_xy_xzz = pbuffer.data(idx_op_df + i * 60 + 15);

        auto to_xy_yyy = pbuffer.data(idx_op_df + i * 60 + 16);

        auto to_xy_yyz = pbuffer.data(idx_op_df + i * 60 + 17);

        auto to_xy_yzz = pbuffer.data(idx_op_df + i * 60 + 18);

        auto to_xy_zzz = pbuffer.data(idx_op_df + i * 60 + 19);

        auto to_xz_xxx = pbuffer.data(idx_op_df + i * 60 + 20);

        auto to_xz_xxy = pbuffer.data(idx_op_df + i * 60 + 21);

        auto to_xz_xxz = pbuffer.data(idx_op_df + i * 60 + 22);

        auto to_xz_xyy = pbuffer.data(idx_op_df + i * 60 + 23);

        auto to_xz_xyz = pbuffer.data(idx_op_df + i * 60 + 24);

        auto to_xz_xzz = pbuffer.data(idx_op_df + i * 60 + 25);

        auto to_xz_yyy = pbuffer.data(idx_op_df + i * 60 + 26);

        auto to_xz_yyz = pbuffer.data(idx_op_df + i * 60 + 27);

        auto to_xz_yzz = pbuffer.data(idx_op_df + i * 60 + 28);

        auto to_xz_zzz = pbuffer.data(idx_op_df + i * 60 + 29);

        auto to_yy_xxx = pbuffer.data(idx_op_df + i * 60 + 30);

        auto to_yy_xxy = pbuffer.data(idx_op_df + i * 60 + 31);

        auto to_yy_xxz = pbuffer.data(idx_op_df + i * 60 + 32);

        auto to_yy_xyy = pbuffer.data(idx_op_df + i * 60 + 33);

        auto to_yy_xyz = pbuffer.data(idx_op_df + i * 60 + 34);

        auto to_yy_xzz = pbuffer.data(idx_op_df + i * 60 + 35);

        auto to_yy_yyy = pbuffer.data(idx_op_df + i * 60 + 36);

        auto to_yy_yyz = pbuffer.data(idx_op_df + i * 60 + 37);

        auto to_yy_yzz = pbuffer.data(idx_op_df + i * 60 + 38);

        auto to_yy_zzz = pbuffer.data(idx_op_df + i * 60 + 39);

        auto to_yz_xxx = pbuffer.data(idx_op_df + i * 60 + 40);

        auto to_yz_xxy = pbuffer.data(idx_op_df + i * 60 + 41);

        auto to_yz_xxz = pbuffer.data(idx_op_df + i * 60 + 42);

        auto to_yz_xyy = pbuffer.data(idx_op_df + i * 60 + 43);

        auto to_yz_xyz = pbuffer.data(idx_op_df + i * 60 + 44);

        auto to_yz_xzz = pbuffer.data(idx_op_df + i * 60 + 45);

        auto to_yz_yyy = pbuffer.data(idx_op_df + i * 60 + 46);

        auto to_yz_yyz = pbuffer.data(idx_op_df + i * 60 + 47);

        auto to_yz_yzz = pbuffer.data(idx_op_df + i * 60 + 48);

        auto to_yz_zzz = pbuffer.data(idx_op_df + i * 60 + 49);

        auto to_zz_xxx = pbuffer.data(idx_op_df + i * 60 + 50);

        auto to_zz_xxy = pbuffer.data(idx_op_df + i * 60 + 51);

        auto to_zz_xxz = pbuffer.data(idx_op_df + i * 60 + 52);

        auto to_zz_xyy = pbuffer.data(idx_op_df + i * 60 + 53);

        auto to_zz_xyz = pbuffer.data(idx_op_df + i * 60 + 54);

        auto to_zz_xzz = pbuffer.data(idx_op_df + i * 60 + 55);

        auto to_zz_yyy = pbuffer.data(idx_op_df + i * 60 + 56);

        auto to_zz_yyz = pbuffer.data(idx_op_df + i * 60 + 57);

        auto to_zz_yzz = pbuffer.data(idx_op_df + i * 60 + 58);

        auto to_zz_zzz = pbuffer.data(idx_op_df + i * 60 + 59);

        // Set up components of auxiliary buffer : DH

        auto to_xx_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 0);

        auto to_xx_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 1);

        auto to_xx_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 2);

        auto to_xx_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 3);

        auto to_xx_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 4);

        auto to_xx_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 5);

        auto to_xx_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 6);

        auto to_xx_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 7);

        auto to_xx_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 8);

        auto to_xx_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 9);

        auto to_xx_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 10);

        auto to_xx_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 11);

        auto to_xx_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 12);

        auto to_xx_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 13);

        auto to_xx_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 14);

        auto to_xx_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 15);

        auto to_xx_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 16);

        auto to_xx_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 17);

        auto to_xx_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 18);

        auto to_xx_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 19);

        auto to_xx_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 20);

        auto to_xy_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 21);

        auto to_xy_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 22);

        auto to_xy_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 23);

        auto to_xy_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 24);

        auto to_xy_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 25);

        auto to_xy_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 26);

        auto to_xy_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 27);

        auto to_xy_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 28);

        auto to_xy_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 29);

        auto to_xy_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 30);

        auto to_xy_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 31);

        auto to_xy_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 32);

        auto to_xy_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 33);

        auto to_xy_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 34);

        auto to_xy_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 35);

        auto to_xy_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 36);

        auto to_xy_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 37);

        auto to_xy_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 38);

        auto to_xy_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 39);

        auto to_xy_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 40);

        auto to_xy_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 41);

        auto to_xz_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 42);

        auto to_xz_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 43);

        auto to_xz_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 44);

        auto to_xz_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 45);

        auto to_xz_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 46);

        auto to_xz_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 47);

        auto to_xz_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 48);

        auto to_xz_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 49);

        auto to_xz_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 50);

        auto to_xz_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 51);

        auto to_xz_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 52);

        auto to_xz_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 53);

        auto to_xz_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 54);

        auto to_xz_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 55);

        auto to_xz_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 56);

        auto to_xz_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 57);

        auto to_xz_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 58);

        auto to_xz_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 59);

        auto to_xz_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 60);

        auto to_xz_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 61);

        auto to_xz_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 62);

        auto to_yy_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 63);

        auto to_yy_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 64);

        auto to_yy_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 65);

        auto to_yy_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 66);

        auto to_yy_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 67);

        auto to_yy_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 68);

        auto to_yy_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 69);

        auto to_yy_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 70);

        auto to_yy_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 71);

        auto to_yy_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 72);

        auto to_yy_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 73);

        auto to_yy_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 74);

        auto to_yy_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 75);

        auto to_yy_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 76);

        auto to_yy_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 77);

        auto to_yy_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 78);

        auto to_yy_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 79);

        auto to_yy_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 80);

        auto to_yy_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 81);

        auto to_yy_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 82);

        auto to_yy_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 83);

        auto to_yz_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 84);

        auto to_yz_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 85);

        auto to_yz_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 86);

        auto to_yz_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 87);

        auto to_yz_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 88);

        auto to_yz_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 89);

        auto to_yz_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 90);

        auto to_yz_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 91);

        auto to_yz_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 92);

        auto to_yz_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 93);

        auto to_yz_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 94);

        auto to_yz_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 95);

        auto to_yz_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 96);

        auto to_yz_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 97);

        auto to_yz_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 98);

        auto to_yz_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 99);

        auto to_yz_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 100);

        auto to_yz_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 101);

        auto to_yz_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 102);

        auto to_yz_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 103);

        auto to_yz_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 104);

        auto to_zz_xxxxx = pbuffer.data(idx_op_dh + i * 126 + 105);

        auto to_zz_xxxxy = pbuffer.data(idx_op_dh + i * 126 + 106);

        auto to_zz_xxxxz = pbuffer.data(idx_op_dh + i * 126 + 107);

        auto to_zz_xxxyy = pbuffer.data(idx_op_dh + i * 126 + 108);

        auto to_zz_xxxyz = pbuffer.data(idx_op_dh + i * 126 + 109);

        auto to_zz_xxxzz = pbuffer.data(idx_op_dh + i * 126 + 110);

        auto to_zz_xxyyy = pbuffer.data(idx_op_dh + i * 126 + 111);

        auto to_zz_xxyyz = pbuffer.data(idx_op_dh + i * 126 + 112);

        auto to_zz_xxyzz = pbuffer.data(idx_op_dh + i * 126 + 113);

        auto to_zz_xxzzz = pbuffer.data(idx_op_dh + i * 126 + 114);

        auto to_zz_xyyyy = pbuffer.data(idx_op_dh + i * 126 + 115);

        auto to_zz_xyyyz = pbuffer.data(idx_op_dh + i * 126 + 116);

        auto to_zz_xyyzz = pbuffer.data(idx_op_dh + i * 126 + 117);

        auto to_zz_xyzzz = pbuffer.data(idx_op_dh + i * 126 + 118);

        auto to_zz_xzzzz = pbuffer.data(idx_op_dh + i * 126 + 119);

        auto to_zz_yyyyy = pbuffer.data(idx_op_dh + i * 126 + 120);

        auto to_zz_yyyyz = pbuffer.data(idx_op_dh + i * 126 + 121);

        auto to_zz_yyyzz = pbuffer.data(idx_op_dh + i * 126 + 122);

        auto to_zz_yyzzz = pbuffer.data(idx_op_dh + i * 126 + 123);

        auto to_zz_yzzzz = pbuffer.data(idx_op_dh + i * 126 + 124);

        auto to_zz_zzzzz = pbuffer.data(idx_op_dh + i * 126 + 125);

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

        // Set up 0-15 components of targeted buffer : FG

        auto to_x_x_xxx_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 0);

        auto to_x_x_xxx_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 1);

        auto to_x_x_xxx_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 2);

        auto to_x_x_xxx_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 3);

        auto to_x_x_xxx_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 4);

        auto to_x_x_xxx_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 5);

        auto to_x_x_xxx_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 6);

        auto to_x_x_xxx_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 7);

        auto to_x_x_xxx_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 8);

        auto to_x_x_xxx_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 9);

        auto to_x_x_xxx_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 10);

        auto to_x_x_xxx_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 11);

        auto to_x_x_xxx_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 12);

        auto to_x_x_xxx_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 13);

        auto to_x_x_xxx_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_x_x_xxx_xxxx, to_x_x_xxx_xxxy, to_x_x_xxx_xxxz, to_x_x_xxx_xxyy, to_x_x_xxx_xxyz, to_x_x_xxx_xxzz, to_x_x_xxx_xyyy, to_x_x_xxx_xyyz, to_x_x_xxx_xyzz, to_x_x_xxx_xzzz, to_x_x_xxx_yyyy, to_x_x_xxx_yyyz, to_x_x_xxx_yyzz, to_x_x_xxx_yzzz, to_x_x_xxx_zzzz, to_xx_xxx, to_xx_xxxxx, to_xx_xxxxy, to_xx_xxxxz, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyz, to_xx_yzz, to_xx_zzz, to_xxxx_xxx, to_xxxx_xxxxx, to_xxxx_xxxxy, to_xxxx_xxxxz, to_xxxx_xxxyy, to_xxxx_xxxyz, to_xxxx_xxxzz, to_xxxx_xxy, to_xxxx_xxyyy, to_xxxx_xxyyz, to_xxxx_xxyzz, to_xxxx_xxz, to_xxxx_xxzzz, to_xxxx_xyy, to_xxxx_xyyyy, to_xxxx_xyyyz, to_xxxx_xyyzz, to_xxxx_xyz, to_xxxx_xyzzz, to_xxxx_xzz, to_xxxx_xzzzz, to_xxxx_yyy, to_xxxx_yyz, to_xxxx_yzz, to_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxx_xxxx[k] = 12.0 * to_xx_xxx[k] - 6.0 * to_xx_xxxxx[k] * tke_0 - 8.0 * to_xxxx_xxx[k] * tbe_0 + 4.0 * to_xxxx_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xxx_xxxy[k] = 9.0 * to_xx_xxy[k] - 6.0 * to_xx_xxxxy[k] * tke_0 - 6.0 * to_xxxx_xxy[k] * tbe_0 + 4.0 * to_xxxx_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xxx_xxxz[k] = 9.0 * to_xx_xxz[k] - 6.0 * to_xx_xxxxz[k] * tke_0 - 6.0 * to_xxxx_xxz[k] * tbe_0 + 4.0 * to_xxxx_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xxx_xxyy[k] = 6.0 * to_xx_xyy[k] - 6.0 * to_xx_xxxyy[k] * tke_0 - 4.0 * to_xxxx_xyy[k] * tbe_0 + 4.0 * to_xxxx_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xxx_xxyz[k] = 6.0 * to_xx_xyz[k] - 6.0 * to_xx_xxxyz[k] * tke_0 - 4.0 * to_xxxx_xyz[k] * tbe_0 + 4.0 * to_xxxx_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xxx_xxzz[k] = 6.0 * to_xx_xzz[k] - 6.0 * to_xx_xxxzz[k] * tke_0 - 4.0 * to_xxxx_xzz[k] * tbe_0 + 4.0 * to_xxxx_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xxx_xyyy[k] = 3.0 * to_xx_yyy[k] - 6.0 * to_xx_xxyyy[k] * tke_0 - 2.0 * to_xxxx_yyy[k] * tbe_0 + 4.0 * to_xxxx_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xxx_xyyz[k] = 3.0 * to_xx_yyz[k] - 6.0 * to_xx_xxyyz[k] * tke_0 - 2.0 * to_xxxx_yyz[k] * tbe_0 + 4.0 * to_xxxx_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xxx_xyzz[k] = 3.0 * to_xx_yzz[k] - 6.0 * to_xx_xxyzz[k] * tke_0 - 2.0 * to_xxxx_yzz[k] * tbe_0 + 4.0 * to_xxxx_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xxx_xzzz[k] = 3.0 * to_xx_zzz[k] - 6.0 * to_xx_xxzzz[k] * tke_0 - 2.0 * to_xxxx_zzz[k] * tbe_0 + 4.0 * to_xxxx_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xxx_yyyy[k] = -6.0 * to_xx_xyyyy[k] * tke_0 + 4.0 * to_xxxx_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xxx_yyyz[k] = -6.0 * to_xx_xyyyz[k] * tke_0 + 4.0 * to_xxxx_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xxx_yyzz[k] = -6.0 * to_xx_xyyzz[k] * tke_0 + 4.0 * to_xxxx_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xxx_yzzz[k] = -6.0 * to_xx_xyzzz[k] * tke_0 + 4.0 * to_xxxx_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xxx_zzzz[k] = -6.0 * to_xx_xzzzz[k] * tke_0 + 4.0 * to_xxxx_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 15-30 components of targeted buffer : FG

        auto to_x_x_xxy_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 15);

        auto to_x_x_xxy_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 16);

        auto to_x_x_xxy_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 17);

        auto to_x_x_xxy_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 18);

        auto to_x_x_xxy_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 19);

        auto to_x_x_xxy_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 20);

        auto to_x_x_xxy_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 21);

        auto to_x_x_xxy_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 22);

        auto to_x_x_xxy_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 23);

        auto to_x_x_xxy_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 24);

        auto to_x_x_xxy_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 25);

        auto to_x_x_xxy_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 26);

        auto to_x_x_xxy_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 27);

        auto to_x_x_xxy_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 28);

        auto to_x_x_xxy_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_x_x_xxy_xxxx, to_x_x_xxy_xxxy, to_x_x_xxy_xxxz, to_x_x_xxy_xxyy, to_x_x_xxy_xxyz, to_x_x_xxy_xxzz, to_x_x_xxy_xyyy, to_x_x_xxy_xyyz, to_x_x_xxy_xyzz, to_x_x_xxy_xzzz, to_x_x_xxy_yyyy, to_x_x_xxy_yyyz, to_x_x_xxy_yyzz, to_x_x_xxy_yzzz, to_x_x_xxy_zzzz, to_xxxy_xxx, to_xxxy_xxxxx, to_xxxy_xxxxy, to_xxxy_xxxxz, to_xxxy_xxxyy, to_xxxy_xxxyz, to_xxxy_xxxzz, to_xxxy_xxy, to_xxxy_xxyyy, to_xxxy_xxyyz, to_xxxy_xxyzz, to_xxxy_xxz, to_xxxy_xxzzz, to_xxxy_xyy, to_xxxy_xyyyy, to_xxxy_xyyyz, to_xxxy_xyyzz, to_xxxy_xyz, to_xxxy_xyzzz, to_xxxy_xzz, to_xxxy_xzzzz, to_xxxy_yyy, to_xxxy_yyz, to_xxxy_yzz, to_xxxy_zzz, to_xy_xxx, to_xy_xxxxx, to_xy_xxxxy, to_xy_xxxxz, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxy_xxxx[k] = 8.0 * to_xy_xxx[k] - 4.0 * to_xy_xxxxx[k] * tke_0 - 8.0 * to_xxxy_xxx[k] * tbe_0 + 4.0 * to_xxxy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xxy_xxxy[k] = 6.0 * to_xy_xxy[k] - 4.0 * to_xy_xxxxy[k] * tke_0 - 6.0 * to_xxxy_xxy[k] * tbe_0 + 4.0 * to_xxxy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xxy_xxxz[k] = 6.0 * to_xy_xxz[k] - 4.0 * to_xy_xxxxz[k] * tke_0 - 6.0 * to_xxxy_xxz[k] * tbe_0 + 4.0 * to_xxxy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xxy_xxyy[k] = 4.0 * to_xy_xyy[k] - 4.0 * to_xy_xxxyy[k] * tke_0 - 4.0 * to_xxxy_xyy[k] * tbe_0 + 4.0 * to_xxxy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xxy_xxyz[k] = 4.0 * to_xy_xyz[k] - 4.0 * to_xy_xxxyz[k] * tke_0 - 4.0 * to_xxxy_xyz[k] * tbe_0 + 4.0 * to_xxxy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xxy_xxzz[k] = 4.0 * to_xy_xzz[k] - 4.0 * to_xy_xxxzz[k] * tke_0 - 4.0 * to_xxxy_xzz[k] * tbe_0 + 4.0 * to_xxxy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xxy_xyyy[k] = 2.0 * to_xy_yyy[k] - 4.0 * to_xy_xxyyy[k] * tke_0 - 2.0 * to_xxxy_yyy[k] * tbe_0 + 4.0 * to_xxxy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xxy_xyyz[k] = 2.0 * to_xy_yyz[k] - 4.0 * to_xy_xxyyz[k] * tke_0 - 2.0 * to_xxxy_yyz[k] * tbe_0 + 4.0 * to_xxxy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xxy_xyzz[k] = 2.0 * to_xy_yzz[k] - 4.0 * to_xy_xxyzz[k] * tke_0 - 2.0 * to_xxxy_yzz[k] * tbe_0 + 4.0 * to_xxxy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xxy_xzzz[k] = 2.0 * to_xy_zzz[k] - 4.0 * to_xy_xxzzz[k] * tke_0 - 2.0 * to_xxxy_zzz[k] * tbe_0 + 4.0 * to_xxxy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xxy_yyyy[k] = -4.0 * to_xy_xyyyy[k] * tke_0 + 4.0 * to_xxxy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xxy_yyyz[k] = -4.0 * to_xy_xyyyz[k] * tke_0 + 4.0 * to_xxxy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xxy_yyzz[k] = -4.0 * to_xy_xyyzz[k] * tke_0 + 4.0 * to_xxxy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xxy_yzzz[k] = -4.0 * to_xy_xyzzz[k] * tke_0 + 4.0 * to_xxxy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xxy_zzzz[k] = -4.0 * to_xy_xzzzz[k] * tke_0 + 4.0 * to_xxxy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-45 components of targeted buffer : FG

        auto to_x_x_xxz_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 30);

        auto to_x_x_xxz_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 31);

        auto to_x_x_xxz_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 32);

        auto to_x_x_xxz_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 33);

        auto to_x_x_xxz_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 34);

        auto to_x_x_xxz_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 35);

        auto to_x_x_xxz_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 36);

        auto to_x_x_xxz_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 37);

        auto to_x_x_xxz_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 38);

        auto to_x_x_xxz_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 39);

        auto to_x_x_xxz_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 40);

        auto to_x_x_xxz_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 41);

        auto to_x_x_xxz_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 42);

        auto to_x_x_xxz_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 43);

        auto to_x_x_xxz_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_x_x_xxz_xxxx, to_x_x_xxz_xxxy, to_x_x_xxz_xxxz, to_x_x_xxz_xxyy, to_x_x_xxz_xxyz, to_x_x_xxz_xxzz, to_x_x_xxz_xyyy, to_x_x_xxz_xyyz, to_x_x_xxz_xyzz, to_x_x_xxz_xzzz, to_x_x_xxz_yyyy, to_x_x_xxz_yyyz, to_x_x_xxz_yyzz, to_x_x_xxz_yzzz, to_x_x_xxz_zzzz, to_xxxz_xxx, to_xxxz_xxxxx, to_xxxz_xxxxy, to_xxxz_xxxxz, to_xxxz_xxxyy, to_xxxz_xxxyz, to_xxxz_xxxzz, to_xxxz_xxy, to_xxxz_xxyyy, to_xxxz_xxyyz, to_xxxz_xxyzz, to_xxxz_xxz, to_xxxz_xxzzz, to_xxxz_xyy, to_xxxz_xyyyy, to_xxxz_xyyyz, to_xxxz_xyyzz, to_xxxz_xyz, to_xxxz_xyzzz, to_xxxz_xzz, to_xxxz_xzzzz, to_xxxz_yyy, to_xxxz_yyz, to_xxxz_yzz, to_xxxz_zzz, to_xz_xxx, to_xz_xxxxx, to_xz_xxxxy, to_xz_xxxxz, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxz_xxxx[k] = 8.0 * to_xz_xxx[k] - 4.0 * to_xz_xxxxx[k] * tke_0 - 8.0 * to_xxxz_xxx[k] * tbe_0 + 4.0 * to_xxxz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xxz_xxxy[k] = 6.0 * to_xz_xxy[k] - 4.0 * to_xz_xxxxy[k] * tke_0 - 6.0 * to_xxxz_xxy[k] * tbe_0 + 4.0 * to_xxxz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xxz_xxxz[k] = 6.0 * to_xz_xxz[k] - 4.0 * to_xz_xxxxz[k] * tke_0 - 6.0 * to_xxxz_xxz[k] * tbe_0 + 4.0 * to_xxxz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xxz_xxyy[k] = 4.0 * to_xz_xyy[k] - 4.0 * to_xz_xxxyy[k] * tke_0 - 4.0 * to_xxxz_xyy[k] * tbe_0 + 4.0 * to_xxxz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xxz_xxyz[k] = 4.0 * to_xz_xyz[k] - 4.0 * to_xz_xxxyz[k] * tke_0 - 4.0 * to_xxxz_xyz[k] * tbe_0 + 4.0 * to_xxxz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xxz_xxzz[k] = 4.0 * to_xz_xzz[k] - 4.0 * to_xz_xxxzz[k] * tke_0 - 4.0 * to_xxxz_xzz[k] * tbe_0 + 4.0 * to_xxxz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xxz_xyyy[k] = 2.0 * to_xz_yyy[k] - 4.0 * to_xz_xxyyy[k] * tke_0 - 2.0 * to_xxxz_yyy[k] * tbe_0 + 4.0 * to_xxxz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xxz_xyyz[k] = 2.0 * to_xz_yyz[k] - 4.0 * to_xz_xxyyz[k] * tke_0 - 2.0 * to_xxxz_yyz[k] * tbe_0 + 4.0 * to_xxxz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xxz_xyzz[k] = 2.0 * to_xz_yzz[k] - 4.0 * to_xz_xxyzz[k] * tke_0 - 2.0 * to_xxxz_yzz[k] * tbe_0 + 4.0 * to_xxxz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xxz_xzzz[k] = 2.0 * to_xz_zzz[k] - 4.0 * to_xz_xxzzz[k] * tke_0 - 2.0 * to_xxxz_zzz[k] * tbe_0 + 4.0 * to_xxxz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xxz_yyyy[k] = -4.0 * to_xz_xyyyy[k] * tke_0 + 4.0 * to_xxxz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xxz_yyyz[k] = -4.0 * to_xz_xyyyz[k] * tke_0 + 4.0 * to_xxxz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xxz_yyzz[k] = -4.0 * to_xz_xyyzz[k] * tke_0 + 4.0 * to_xxxz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xxz_yzzz[k] = -4.0 * to_xz_xyzzz[k] * tke_0 + 4.0 * to_xxxz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xxz_zzzz[k] = -4.0 * to_xz_xzzzz[k] * tke_0 + 4.0 * to_xxxz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 45-60 components of targeted buffer : FG

        auto to_x_x_xyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 45);

        auto to_x_x_xyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 46);

        auto to_x_x_xyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 47);

        auto to_x_x_xyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 48);

        auto to_x_x_xyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 49);

        auto to_x_x_xyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 50);

        auto to_x_x_xyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 51);

        auto to_x_x_xyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 52);

        auto to_x_x_xyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 53);

        auto to_x_x_xyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 54);

        auto to_x_x_xyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 55);

        auto to_x_x_xyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 56);

        auto to_x_x_xyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 57);

        auto to_x_x_xyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 58);

        auto to_x_x_xyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_x_x_xyy_xxxx, to_x_x_xyy_xxxy, to_x_x_xyy_xxxz, to_x_x_xyy_xxyy, to_x_x_xyy_xxyz, to_x_x_xyy_xxzz, to_x_x_xyy_xyyy, to_x_x_xyy_xyyz, to_x_x_xyy_xyzz, to_x_x_xyy_xzzz, to_x_x_xyy_yyyy, to_x_x_xyy_yyyz, to_x_x_xyy_yyzz, to_x_x_xyy_yzzz, to_x_x_xyy_zzzz, to_xxyy_xxx, to_xxyy_xxxxx, to_xxyy_xxxxy, to_xxyy_xxxxz, to_xxyy_xxxyy, to_xxyy_xxxyz, to_xxyy_xxxzz, to_xxyy_xxy, to_xxyy_xxyyy, to_xxyy_xxyyz, to_xxyy_xxyzz, to_xxyy_xxz, to_xxyy_xxzzz, to_xxyy_xyy, to_xxyy_xyyyy, to_xxyy_xyyyz, to_xxyy_xyyzz, to_xxyy_xyz, to_xxyy_xyzzz, to_xxyy_xzz, to_xxyy_xzzzz, to_xxyy_yyy, to_xxyy_yyz, to_xxyy_yzz, to_xxyy_zzz, to_yy_xxx, to_yy_xxxxx, to_yy_xxxxy, to_yy_xxxxz, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyz, to_yy_yzz, to_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyy_xxxx[k] = 4.0 * to_yy_xxx[k] - 2.0 * to_yy_xxxxx[k] * tke_0 - 8.0 * to_xxyy_xxx[k] * tbe_0 + 4.0 * to_xxyy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xyy_xxxy[k] = 3.0 * to_yy_xxy[k] - 2.0 * to_yy_xxxxy[k] * tke_0 - 6.0 * to_xxyy_xxy[k] * tbe_0 + 4.0 * to_xxyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xyy_xxxz[k] = 3.0 * to_yy_xxz[k] - 2.0 * to_yy_xxxxz[k] * tke_0 - 6.0 * to_xxyy_xxz[k] * tbe_0 + 4.0 * to_xxyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xyy_xxyy[k] = 2.0 * to_yy_xyy[k] - 2.0 * to_yy_xxxyy[k] * tke_0 - 4.0 * to_xxyy_xyy[k] * tbe_0 + 4.0 * to_xxyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xyy_xxyz[k] = 2.0 * to_yy_xyz[k] - 2.0 * to_yy_xxxyz[k] * tke_0 - 4.0 * to_xxyy_xyz[k] * tbe_0 + 4.0 * to_xxyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xyy_xxzz[k] = 2.0 * to_yy_xzz[k] - 2.0 * to_yy_xxxzz[k] * tke_0 - 4.0 * to_xxyy_xzz[k] * tbe_0 + 4.0 * to_xxyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xyy_xyyy[k] = to_yy_yyy[k] - 2.0 * to_yy_xxyyy[k] * tke_0 - 2.0 * to_xxyy_yyy[k] * tbe_0 + 4.0 * to_xxyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xyy_xyyz[k] = to_yy_yyz[k] - 2.0 * to_yy_xxyyz[k] * tke_0 - 2.0 * to_xxyy_yyz[k] * tbe_0 + 4.0 * to_xxyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xyy_xyzz[k] = to_yy_yzz[k] - 2.0 * to_yy_xxyzz[k] * tke_0 - 2.0 * to_xxyy_yzz[k] * tbe_0 + 4.0 * to_xxyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xyy_xzzz[k] = to_yy_zzz[k] - 2.0 * to_yy_xxzzz[k] * tke_0 - 2.0 * to_xxyy_zzz[k] * tbe_0 + 4.0 * to_xxyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xyy_yyyy[k] = -2.0 * to_yy_xyyyy[k] * tke_0 + 4.0 * to_xxyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xyy_yyyz[k] = -2.0 * to_yy_xyyyz[k] * tke_0 + 4.0 * to_xxyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xyy_yyzz[k] = -2.0 * to_yy_xyyzz[k] * tke_0 + 4.0 * to_xxyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xyy_yzzz[k] = -2.0 * to_yy_xyzzz[k] * tke_0 + 4.0 * to_xxyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xyy_zzzz[k] = -2.0 * to_yy_xzzzz[k] * tke_0 + 4.0 * to_xxyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-75 components of targeted buffer : FG

        auto to_x_x_xyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 60);

        auto to_x_x_xyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 61);

        auto to_x_x_xyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 62);

        auto to_x_x_xyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 63);

        auto to_x_x_xyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 64);

        auto to_x_x_xyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 65);

        auto to_x_x_xyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 66);

        auto to_x_x_xyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 67);

        auto to_x_x_xyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 68);

        auto to_x_x_xyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 69);

        auto to_x_x_xyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 70);

        auto to_x_x_xyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 71);

        auto to_x_x_xyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 72);

        auto to_x_x_xyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 73);

        auto to_x_x_xyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_x_x_xyz_xxxx, to_x_x_xyz_xxxy, to_x_x_xyz_xxxz, to_x_x_xyz_xxyy, to_x_x_xyz_xxyz, to_x_x_xyz_xxzz, to_x_x_xyz_xyyy, to_x_x_xyz_xyyz, to_x_x_xyz_xyzz, to_x_x_xyz_xzzz, to_x_x_xyz_yyyy, to_x_x_xyz_yyyz, to_x_x_xyz_yyzz, to_x_x_xyz_yzzz, to_x_x_xyz_zzzz, to_xxyz_xxx, to_xxyz_xxxxx, to_xxyz_xxxxy, to_xxyz_xxxxz, to_xxyz_xxxyy, to_xxyz_xxxyz, to_xxyz_xxxzz, to_xxyz_xxy, to_xxyz_xxyyy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xxzzz, to_xxyz_xyy, to_xxyz_xyyyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_xzzzz, to_xxyz_yyy, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_zzz, to_yz_xxx, to_yz_xxxxx, to_yz_xxxxy, to_yz_xxxxz, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyz_xxxx[k] = 4.0 * to_yz_xxx[k] - 2.0 * to_yz_xxxxx[k] * tke_0 - 8.0 * to_xxyz_xxx[k] * tbe_0 + 4.0 * to_xxyz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xyz_xxxy[k] = 3.0 * to_yz_xxy[k] - 2.0 * to_yz_xxxxy[k] * tke_0 - 6.0 * to_xxyz_xxy[k] * tbe_0 + 4.0 * to_xxyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xyz_xxxz[k] = 3.0 * to_yz_xxz[k] - 2.0 * to_yz_xxxxz[k] * tke_0 - 6.0 * to_xxyz_xxz[k] * tbe_0 + 4.0 * to_xxyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xyz_xxyy[k] = 2.0 * to_yz_xyy[k] - 2.0 * to_yz_xxxyy[k] * tke_0 - 4.0 * to_xxyz_xyy[k] * tbe_0 + 4.0 * to_xxyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xyz_xxyz[k] = 2.0 * to_yz_xyz[k] - 2.0 * to_yz_xxxyz[k] * tke_0 - 4.0 * to_xxyz_xyz[k] * tbe_0 + 4.0 * to_xxyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xyz_xxzz[k] = 2.0 * to_yz_xzz[k] - 2.0 * to_yz_xxxzz[k] * tke_0 - 4.0 * to_xxyz_xzz[k] * tbe_0 + 4.0 * to_xxyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xyz_xyyy[k] = to_yz_yyy[k] - 2.0 * to_yz_xxyyy[k] * tke_0 - 2.0 * to_xxyz_yyy[k] * tbe_0 + 4.0 * to_xxyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xyz_xyyz[k] = to_yz_yyz[k] - 2.0 * to_yz_xxyyz[k] * tke_0 - 2.0 * to_xxyz_yyz[k] * tbe_0 + 4.0 * to_xxyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xyz_xyzz[k] = to_yz_yzz[k] - 2.0 * to_yz_xxyzz[k] * tke_0 - 2.0 * to_xxyz_yzz[k] * tbe_0 + 4.0 * to_xxyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xyz_xzzz[k] = to_yz_zzz[k] - 2.0 * to_yz_xxzzz[k] * tke_0 - 2.0 * to_xxyz_zzz[k] * tbe_0 + 4.0 * to_xxyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xyz_yyyy[k] = -2.0 * to_yz_xyyyy[k] * tke_0 + 4.0 * to_xxyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xyz_yyyz[k] = -2.0 * to_yz_xyyyz[k] * tke_0 + 4.0 * to_xxyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xyz_yyzz[k] = -2.0 * to_yz_xyyzz[k] * tke_0 + 4.0 * to_xxyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xyz_yzzz[k] = -2.0 * to_yz_xyzzz[k] * tke_0 + 4.0 * to_xxyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xyz_zzzz[k] = -2.0 * to_yz_xzzzz[k] * tke_0 + 4.0 * to_xxyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 75-90 components of targeted buffer : FG

        auto to_x_x_xzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 75);

        auto to_x_x_xzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 76);

        auto to_x_x_xzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 77);

        auto to_x_x_xzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 78);

        auto to_x_x_xzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 79);

        auto to_x_x_xzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 80);

        auto to_x_x_xzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 81);

        auto to_x_x_xzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 82);

        auto to_x_x_xzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 83);

        auto to_x_x_xzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 84);

        auto to_x_x_xzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 85);

        auto to_x_x_xzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 86);

        auto to_x_x_xzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 87);

        auto to_x_x_xzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 88);

        auto to_x_x_xzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_x_x_xzz_xxxx, to_x_x_xzz_xxxy, to_x_x_xzz_xxxz, to_x_x_xzz_xxyy, to_x_x_xzz_xxyz, to_x_x_xzz_xxzz, to_x_x_xzz_xyyy, to_x_x_xzz_xyyz, to_x_x_xzz_xyzz, to_x_x_xzz_xzzz, to_x_x_xzz_yyyy, to_x_x_xzz_yyyz, to_x_x_xzz_yyzz, to_x_x_xzz_yzzz, to_x_x_xzz_zzzz, to_xxzz_xxx, to_xxzz_xxxxx, to_xxzz_xxxxy, to_xxzz_xxxxz, to_xxzz_xxxyy, to_xxzz_xxxyz, to_xxzz_xxxzz, to_xxzz_xxy, to_xxzz_xxyyy, to_xxzz_xxyyz, to_xxzz_xxyzz, to_xxzz_xxz, to_xxzz_xxzzz, to_xxzz_xyy, to_xxzz_xyyyy, to_xxzz_xyyyz, to_xxzz_xyyzz, to_xxzz_xyz, to_xxzz_xyzzz, to_xxzz_xzz, to_xxzz_xzzzz, to_xxzz_yyy, to_xxzz_yyz, to_xxzz_yzz, to_xxzz_zzz, to_zz_xxx, to_zz_xxxxx, to_zz_xxxxy, to_zz_xxxxz, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyz, to_zz_yzz, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xzz_xxxx[k] = 4.0 * to_zz_xxx[k] - 2.0 * to_zz_xxxxx[k] * tke_0 - 8.0 * to_xxzz_xxx[k] * tbe_0 + 4.0 * to_xxzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_xzz_xxxy[k] = 3.0 * to_zz_xxy[k] - 2.0 * to_zz_xxxxy[k] * tke_0 - 6.0 * to_xxzz_xxy[k] * tbe_0 + 4.0 * to_xxzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_xzz_xxxz[k] = 3.0 * to_zz_xxz[k] - 2.0 * to_zz_xxxxz[k] * tke_0 - 6.0 * to_xxzz_xxz[k] * tbe_0 + 4.0 * to_xxzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_xzz_xxyy[k] = 2.0 * to_zz_xyy[k] - 2.0 * to_zz_xxxyy[k] * tke_0 - 4.0 * to_xxzz_xyy[k] * tbe_0 + 4.0 * to_xxzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_xzz_xxyz[k] = 2.0 * to_zz_xyz[k] - 2.0 * to_zz_xxxyz[k] * tke_0 - 4.0 * to_xxzz_xyz[k] * tbe_0 + 4.0 * to_xxzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_xzz_xxzz[k] = 2.0 * to_zz_xzz[k] - 2.0 * to_zz_xxxzz[k] * tke_0 - 4.0 * to_xxzz_xzz[k] * tbe_0 + 4.0 * to_xxzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_xzz_xyyy[k] = to_zz_yyy[k] - 2.0 * to_zz_xxyyy[k] * tke_0 - 2.0 * to_xxzz_yyy[k] * tbe_0 + 4.0 * to_xxzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_xzz_xyyz[k] = to_zz_yyz[k] - 2.0 * to_zz_xxyyz[k] * tke_0 - 2.0 * to_xxzz_yyz[k] * tbe_0 + 4.0 * to_xxzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_xzz_xyzz[k] = to_zz_yzz[k] - 2.0 * to_zz_xxyzz[k] * tke_0 - 2.0 * to_xxzz_yzz[k] * tbe_0 + 4.0 * to_xxzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_xzz_xzzz[k] = to_zz_zzz[k] - 2.0 * to_zz_xxzzz[k] * tke_0 - 2.0 * to_xxzz_zzz[k] * tbe_0 + 4.0 * to_xxzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_xzz_yyyy[k] = -2.0 * to_zz_xyyyy[k] * tke_0 + 4.0 * to_xxzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_xzz_yyyz[k] = -2.0 * to_zz_xyyyz[k] * tke_0 + 4.0 * to_xxzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_xzz_yyzz[k] = -2.0 * to_zz_xyyzz[k] * tke_0 + 4.0 * to_xxzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_xzz_yzzz[k] = -2.0 * to_zz_xyzzz[k] * tke_0 + 4.0 * to_xxzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_xzz_zzzz[k] = -2.0 * to_zz_xzzzz[k] * tke_0 + 4.0 * to_xxzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-105 components of targeted buffer : FG

        auto to_x_x_yyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 90);

        auto to_x_x_yyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 91);

        auto to_x_x_yyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 92);

        auto to_x_x_yyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 93);

        auto to_x_x_yyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 94);

        auto to_x_x_yyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 95);

        auto to_x_x_yyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 96);

        auto to_x_x_yyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 97);

        auto to_x_x_yyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 98);

        auto to_x_x_yyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 99);

        auto to_x_x_yyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 100);

        auto to_x_x_yyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 101);

        auto to_x_x_yyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 102);

        auto to_x_x_yyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 103);

        auto to_x_x_yyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_x_x_yyy_xxxx, to_x_x_yyy_xxxy, to_x_x_yyy_xxxz, to_x_x_yyy_xxyy, to_x_x_yyy_xxyz, to_x_x_yyy_xxzz, to_x_x_yyy_xyyy, to_x_x_yyy_xyyz, to_x_x_yyy_xyzz, to_x_x_yyy_xzzz, to_x_x_yyy_yyyy, to_x_x_yyy_yyyz, to_x_x_yyy_yyzz, to_x_x_yyy_yzzz, to_x_x_yyy_zzzz, to_xyyy_xxx, to_xyyy_xxxxx, to_xyyy_xxxxy, to_xyyy_xxxxz, to_xyyy_xxxyy, to_xyyy_xxxyz, to_xyyy_xxxzz, to_xyyy_xxy, to_xyyy_xxyyy, to_xyyy_xxyyz, to_xyyy_xxyzz, to_xyyy_xxz, to_xyyy_xxzzz, to_xyyy_xyy, to_xyyy_xyyyy, to_xyyy_xyyyz, to_xyyy_xyyzz, to_xyyy_xyz, to_xyyy_xyzzz, to_xyyy_xzz, to_xyyy_xzzzz, to_xyyy_yyy, to_xyyy_yyz, to_xyyy_yzz, to_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyy_xxxx[k] = -8.0 * to_xyyy_xxx[k] * tbe_0 + 4.0 * to_xyyy_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_yyy_xxxy[k] = -6.0 * to_xyyy_xxy[k] * tbe_0 + 4.0 * to_xyyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_yyy_xxxz[k] = -6.0 * to_xyyy_xxz[k] * tbe_0 + 4.0 * to_xyyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_yyy_xxyy[k] = -4.0 * to_xyyy_xyy[k] * tbe_0 + 4.0 * to_xyyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_yyy_xxyz[k] = -4.0 * to_xyyy_xyz[k] * tbe_0 + 4.0 * to_xyyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_yyy_xxzz[k] = -4.0 * to_xyyy_xzz[k] * tbe_0 + 4.0 * to_xyyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_yyy_xyyy[k] = -2.0 * to_xyyy_yyy[k] * tbe_0 + 4.0 * to_xyyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_yyy_xyyz[k] = -2.0 * to_xyyy_yyz[k] * tbe_0 + 4.0 * to_xyyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_yyy_xyzz[k] = -2.0 * to_xyyy_yzz[k] * tbe_0 + 4.0 * to_xyyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_yyy_xzzz[k] = -2.0 * to_xyyy_zzz[k] * tbe_0 + 4.0 * to_xyyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_yyy_yyyy[k] = 4.0 * to_xyyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_yyy_yyyz[k] = 4.0 * to_xyyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_yyy_yyzz[k] = 4.0 * to_xyyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_yyy_yzzz[k] = 4.0 * to_xyyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_yyy_zzzz[k] = 4.0 * to_xyyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 105-120 components of targeted buffer : FG

        auto to_x_x_yyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 105);

        auto to_x_x_yyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 106);

        auto to_x_x_yyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 107);

        auto to_x_x_yyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 108);

        auto to_x_x_yyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 109);

        auto to_x_x_yyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 110);

        auto to_x_x_yyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 111);

        auto to_x_x_yyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 112);

        auto to_x_x_yyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 113);

        auto to_x_x_yyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 114);

        auto to_x_x_yyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 115);

        auto to_x_x_yyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 116);

        auto to_x_x_yyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 117);

        auto to_x_x_yyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 118);

        auto to_x_x_yyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_x_x_yyz_xxxx, to_x_x_yyz_xxxy, to_x_x_yyz_xxxz, to_x_x_yyz_xxyy, to_x_x_yyz_xxyz, to_x_x_yyz_xxzz, to_x_x_yyz_xyyy, to_x_x_yyz_xyyz, to_x_x_yyz_xyzz, to_x_x_yyz_xzzz, to_x_x_yyz_yyyy, to_x_x_yyz_yyyz, to_x_x_yyz_yyzz, to_x_x_yyz_yzzz, to_x_x_yyz_zzzz, to_xyyz_xxx, to_xyyz_xxxxx, to_xyyz_xxxxy, to_xyyz_xxxxz, to_xyyz_xxxyy, to_xyyz_xxxyz, to_xyyz_xxxzz, to_xyyz_xxy, to_xyyz_xxyyy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xxzzz, to_xyyz_xyy, to_xyyz_xyyyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_xzzzz, to_xyyz_yyy, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyz_xxxx[k] = -8.0 * to_xyyz_xxx[k] * tbe_0 + 4.0 * to_xyyz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_yyz_xxxy[k] = -6.0 * to_xyyz_xxy[k] * tbe_0 + 4.0 * to_xyyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_yyz_xxxz[k] = -6.0 * to_xyyz_xxz[k] * tbe_0 + 4.0 * to_xyyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_yyz_xxyy[k] = -4.0 * to_xyyz_xyy[k] * tbe_0 + 4.0 * to_xyyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_yyz_xxyz[k] = -4.0 * to_xyyz_xyz[k] * tbe_0 + 4.0 * to_xyyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_yyz_xxzz[k] = -4.0 * to_xyyz_xzz[k] * tbe_0 + 4.0 * to_xyyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_yyz_xyyy[k] = -2.0 * to_xyyz_yyy[k] * tbe_0 + 4.0 * to_xyyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_yyz_xyyz[k] = -2.0 * to_xyyz_yyz[k] * tbe_0 + 4.0 * to_xyyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_yyz_xyzz[k] = -2.0 * to_xyyz_yzz[k] * tbe_0 + 4.0 * to_xyyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_yyz_xzzz[k] = -2.0 * to_xyyz_zzz[k] * tbe_0 + 4.0 * to_xyyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_yyz_yyyy[k] = 4.0 * to_xyyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_yyz_yyyz[k] = 4.0 * to_xyyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_yyz_yyzz[k] = 4.0 * to_xyyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_yyz_yzzz[k] = 4.0 * to_xyyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_yyz_zzzz[k] = 4.0 * to_xyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-135 components of targeted buffer : FG

        auto to_x_x_yzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 120);

        auto to_x_x_yzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 121);

        auto to_x_x_yzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 122);

        auto to_x_x_yzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 123);

        auto to_x_x_yzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 124);

        auto to_x_x_yzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 125);

        auto to_x_x_yzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 126);

        auto to_x_x_yzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 127);

        auto to_x_x_yzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 128);

        auto to_x_x_yzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 129);

        auto to_x_x_yzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 130);

        auto to_x_x_yzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 131);

        auto to_x_x_yzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 132);

        auto to_x_x_yzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 133);

        auto to_x_x_yzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_x_x_yzz_xxxx, to_x_x_yzz_xxxy, to_x_x_yzz_xxxz, to_x_x_yzz_xxyy, to_x_x_yzz_xxyz, to_x_x_yzz_xxzz, to_x_x_yzz_xyyy, to_x_x_yzz_xyyz, to_x_x_yzz_xyzz, to_x_x_yzz_xzzz, to_x_x_yzz_yyyy, to_x_x_yzz_yyyz, to_x_x_yzz_yyzz, to_x_x_yzz_yzzz, to_x_x_yzz_zzzz, to_xyzz_xxx, to_xyzz_xxxxx, to_xyzz_xxxxy, to_xyzz_xxxxz, to_xyzz_xxxyy, to_xyzz_xxxyz, to_xyzz_xxxzz, to_xyzz_xxy, to_xyzz_xxyyy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xxzzz, to_xyzz_xyy, to_xyzz_xyyyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_xzzzz, to_xyzz_yyy, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yzz_xxxx[k] = -8.0 * to_xyzz_xxx[k] * tbe_0 + 4.0 * to_xyzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_yzz_xxxy[k] = -6.0 * to_xyzz_xxy[k] * tbe_0 + 4.0 * to_xyzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_yzz_xxxz[k] = -6.0 * to_xyzz_xxz[k] * tbe_0 + 4.0 * to_xyzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_yzz_xxyy[k] = -4.0 * to_xyzz_xyy[k] * tbe_0 + 4.0 * to_xyzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_yzz_xxyz[k] = -4.0 * to_xyzz_xyz[k] * tbe_0 + 4.0 * to_xyzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_yzz_xxzz[k] = -4.0 * to_xyzz_xzz[k] * tbe_0 + 4.0 * to_xyzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_yzz_xyyy[k] = -2.0 * to_xyzz_yyy[k] * tbe_0 + 4.0 * to_xyzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_yzz_xyyz[k] = -2.0 * to_xyzz_yyz[k] * tbe_0 + 4.0 * to_xyzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_yzz_xyzz[k] = -2.0 * to_xyzz_yzz[k] * tbe_0 + 4.0 * to_xyzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_yzz_xzzz[k] = -2.0 * to_xyzz_zzz[k] * tbe_0 + 4.0 * to_xyzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_yzz_yyyy[k] = 4.0 * to_xyzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_yzz_yyyz[k] = 4.0 * to_xyzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_yzz_yyzz[k] = 4.0 * to_xyzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_yzz_yzzz[k] = 4.0 * to_xyzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_yzz_zzzz[k] = 4.0 * to_xyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 135-150 components of targeted buffer : FG

        auto to_x_x_zzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 135);

        auto to_x_x_zzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 136);

        auto to_x_x_zzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 137);

        auto to_x_x_zzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 138);

        auto to_x_x_zzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 139);

        auto to_x_x_zzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 140);

        auto to_x_x_zzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 141);

        auto to_x_x_zzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 142);

        auto to_x_x_zzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 143);

        auto to_x_x_zzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 144);

        auto to_x_x_zzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 145);

        auto to_x_x_zzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 146);

        auto to_x_x_zzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 147);

        auto to_x_x_zzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 148);

        auto to_x_x_zzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 0 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_x_x_zzz_xxxx, to_x_x_zzz_xxxy, to_x_x_zzz_xxxz, to_x_x_zzz_xxyy, to_x_x_zzz_xxyz, to_x_x_zzz_xxzz, to_x_x_zzz_xyyy, to_x_x_zzz_xyyz, to_x_x_zzz_xyzz, to_x_x_zzz_xzzz, to_x_x_zzz_yyyy, to_x_x_zzz_yyyz, to_x_x_zzz_yyzz, to_x_x_zzz_yzzz, to_x_x_zzz_zzzz, to_xzzz_xxx, to_xzzz_xxxxx, to_xzzz_xxxxy, to_xzzz_xxxxz, to_xzzz_xxxyy, to_xzzz_xxxyz, to_xzzz_xxxzz, to_xzzz_xxy, to_xzzz_xxyyy, to_xzzz_xxyyz, to_xzzz_xxyzz, to_xzzz_xxz, to_xzzz_xxzzz, to_xzzz_xyy, to_xzzz_xyyyy, to_xzzz_xyyyz, to_xzzz_xyyzz, to_xzzz_xyz, to_xzzz_xyzzz, to_xzzz_xzz, to_xzzz_xzzzz, to_xzzz_yyy, to_xzzz_yyz, to_xzzz_yzz, to_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zzz_xxxx[k] = -8.0 * to_xzzz_xxx[k] * tbe_0 + 4.0 * to_xzzz_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_zzz_xxxy[k] = -6.0 * to_xzzz_xxy[k] * tbe_0 + 4.0 * to_xzzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_zzz_xxxz[k] = -6.0 * to_xzzz_xxz[k] * tbe_0 + 4.0 * to_xzzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_zzz_xxyy[k] = -4.0 * to_xzzz_xyy[k] * tbe_0 + 4.0 * to_xzzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_zzz_xxyz[k] = -4.0 * to_xzzz_xyz[k] * tbe_0 + 4.0 * to_xzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_zzz_xxzz[k] = -4.0 * to_xzzz_xzz[k] * tbe_0 + 4.0 * to_xzzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_zzz_xyyy[k] = -2.0 * to_xzzz_yyy[k] * tbe_0 + 4.0 * to_xzzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_zzz_xyyz[k] = -2.0 * to_xzzz_yyz[k] * tbe_0 + 4.0 * to_xzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_zzz_xyzz[k] = -2.0 * to_xzzz_yzz[k] * tbe_0 + 4.0 * to_xzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_zzz_xzzz[k] = -2.0 * to_xzzz_zzz[k] * tbe_0 + 4.0 * to_xzzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_zzz_yyyy[k] = 4.0 * to_xzzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_zzz_yyyz[k] = 4.0 * to_xzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_zzz_yyzz[k] = 4.0 * to_xzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_zzz_yzzz[k] = 4.0 * to_xzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_zzz_zzzz[k] = 4.0 * to_xzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-165 components of targeted buffer : FG

        auto to_x_y_xxx_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 0);

        auto to_x_y_xxx_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 1);

        auto to_x_y_xxx_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 2);

        auto to_x_y_xxx_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 3);

        auto to_x_y_xxx_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 4);

        auto to_x_y_xxx_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 5);

        auto to_x_y_xxx_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 6);

        auto to_x_y_xxx_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 7);

        auto to_x_y_xxx_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 8);

        auto to_x_y_xxx_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 9);

        auto to_x_y_xxx_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 10);

        auto to_x_y_xxx_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 11);

        auto to_x_y_xxx_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 12);

        auto to_x_y_xxx_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 13);

        auto to_x_y_xxx_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_x_y_xxx_xxxx, to_x_y_xxx_xxxy, to_x_y_xxx_xxxz, to_x_y_xxx_xxyy, to_x_y_xxx_xxyz, to_x_y_xxx_xxzz, to_x_y_xxx_xyyy, to_x_y_xxx_xyyz, to_x_y_xxx_xyzz, to_x_y_xxx_xzzz, to_x_y_xxx_yyyy, to_x_y_xxx_yyyz, to_x_y_xxx_yyzz, to_x_y_xxx_yzzz, to_x_y_xxx_zzzz, to_xx_xxx, to_xx_xxxxy, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_yyy, to_xx_yyyyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, to_xxxx_xxx, to_xxxx_xxxxy, to_xxxx_xxxyy, to_xxxx_xxxyz, to_xxxx_xxy, to_xxxx_xxyyy, to_xxxx_xxyyz, to_xxxx_xxyzz, to_xxxx_xxz, to_xxxx_xyy, to_xxxx_xyyyy, to_xxxx_xyyyz, to_xxxx_xyyzz, to_xxxx_xyz, to_xxxx_xyzzz, to_xxxx_xzz, to_xxxx_yyy, to_xxxx_yyyyy, to_xxxx_yyyyz, to_xxxx_yyyzz, to_xxxx_yyz, to_xxxx_yyzzz, to_xxxx_yzz, to_xxxx_yzzzz, to_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxx_xxxx[k] = -6.0 * to_xx_xxxxy[k] * tke_0 + 4.0 * to_xxxx_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xxx_xxxy[k] = 3.0 * to_xx_xxx[k] - 6.0 * to_xx_xxxyy[k] * tke_0 - 2.0 * to_xxxx_xxx[k] * tbe_0 + 4.0 * to_xxxx_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xxx_xxxz[k] = -6.0 * to_xx_xxxyz[k] * tke_0 + 4.0 * to_xxxx_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xxx_xxyy[k] = 6.0 * to_xx_xxy[k] - 6.0 * to_xx_xxyyy[k] * tke_0 - 4.0 * to_xxxx_xxy[k] * tbe_0 + 4.0 * to_xxxx_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xxx_xxyz[k] = 3.0 * to_xx_xxz[k] - 6.0 * to_xx_xxyyz[k] * tke_0 - 2.0 * to_xxxx_xxz[k] * tbe_0 + 4.0 * to_xxxx_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xxx_xxzz[k] = -6.0 * to_xx_xxyzz[k] * tke_0 + 4.0 * to_xxxx_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xxx_xyyy[k] = 9.0 * to_xx_xyy[k] - 6.0 * to_xx_xyyyy[k] * tke_0 - 6.0 * to_xxxx_xyy[k] * tbe_0 + 4.0 * to_xxxx_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxx_xyyz[k] = 6.0 * to_xx_xyz[k] - 6.0 * to_xx_xyyyz[k] * tke_0 - 4.0 * to_xxxx_xyz[k] * tbe_0 + 4.0 * to_xxxx_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxx_xyzz[k] = 3.0 * to_xx_xzz[k] - 6.0 * to_xx_xyyzz[k] * tke_0 - 2.0 * to_xxxx_xzz[k] * tbe_0 + 4.0 * to_xxxx_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxx_xzzz[k] = -6.0 * to_xx_xyzzz[k] * tke_0 + 4.0 * to_xxxx_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxx_yyyy[k] = 12.0 * to_xx_yyy[k] - 6.0 * to_xx_yyyyy[k] * tke_0 - 8.0 * to_xxxx_yyy[k] * tbe_0 + 4.0 * to_xxxx_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxx_yyyz[k] = 9.0 * to_xx_yyz[k] - 6.0 * to_xx_yyyyz[k] * tke_0 - 6.0 * to_xxxx_yyz[k] * tbe_0 + 4.0 * to_xxxx_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxx_yyzz[k] = 6.0 * to_xx_yzz[k] - 6.0 * to_xx_yyyzz[k] * tke_0 - 4.0 * to_xxxx_yzz[k] * tbe_0 + 4.0 * to_xxxx_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxx_yzzz[k] = 3.0 * to_xx_zzz[k] - 6.0 * to_xx_yyzzz[k] * tke_0 - 2.0 * to_xxxx_zzz[k] * tbe_0 + 4.0 * to_xxxx_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxx_zzzz[k] = -6.0 * to_xx_yzzzz[k] * tke_0 + 4.0 * to_xxxx_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 165-180 components of targeted buffer : FG

        auto to_x_y_xxy_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 15);

        auto to_x_y_xxy_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 16);

        auto to_x_y_xxy_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 17);

        auto to_x_y_xxy_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 18);

        auto to_x_y_xxy_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 19);

        auto to_x_y_xxy_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 20);

        auto to_x_y_xxy_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 21);

        auto to_x_y_xxy_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 22);

        auto to_x_y_xxy_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 23);

        auto to_x_y_xxy_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 24);

        auto to_x_y_xxy_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 25);

        auto to_x_y_xxy_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 26);

        auto to_x_y_xxy_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 27);

        auto to_x_y_xxy_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 28);

        auto to_x_y_xxy_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_x_y_xxy_xxxx, to_x_y_xxy_xxxy, to_x_y_xxy_xxxz, to_x_y_xxy_xxyy, to_x_y_xxy_xxyz, to_x_y_xxy_xxzz, to_x_y_xxy_xyyy, to_x_y_xxy_xyyz, to_x_y_xxy_xyzz, to_x_y_xxy_xzzz, to_x_y_xxy_yyyy, to_x_y_xxy_yyyz, to_x_y_xxy_yyzz, to_x_y_xxy_yzzz, to_x_y_xxy_zzzz, to_xxxy_xxx, to_xxxy_xxxxy, to_xxxy_xxxyy, to_xxxy_xxxyz, to_xxxy_xxy, to_xxxy_xxyyy, to_xxxy_xxyyz, to_xxxy_xxyzz, to_xxxy_xxz, to_xxxy_xyy, to_xxxy_xyyyy, to_xxxy_xyyyz, to_xxxy_xyyzz, to_xxxy_xyz, to_xxxy_xyzzz, to_xxxy_xzz, to_xxxy_yyy, to_xxxy_yyyyy, to_xxxy_yyyyz, to_xxxy_yyyzz, to_xxxy_yyz, to_xxxy_yyzzz, to_xxxy_yzz, to_xxxy_yzzzz, to_xxxy_zzz, to_xy_xxx, to_xy_xxxxy, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_yyy, to_xy_yyyyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxy_xxxx[k] = -4.0 * to_xy_xxxxy[k] * tke_0 + 4.0 * to_xxxy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xxy_xxxy[k] = 2.0 * to_xy_xxx[k] - 4.0 * to_xy_xxxyy[k] * tke_0 - 2.0 * to_xxxy_xxx[k] * tbe_0 + 4.0 * to_xxxy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xxy_xxxz[k] = -4.0 * to_xy_xxxyz[k] * tke_0 + 4.0 * to_xxxy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xxy_xxyy[k] = 4.0 * to_xy_xxy[k] - 4.0 * to_xy_xxyyy[k] * tke_0 - 4.0 * to_xxxy_xxy[k] * tbe_0 + 4.0 * to_xxxy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xxy_xxyz[k] = 2.0 * to_xy_xxz[k] - 4.0 * to_xy_xxyyz[k] * tke_0 - 2.0 * to_xxxy_xxz[k] * tbe_0 + 4.0 * to_xxxy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xxy_xxzz[k] = -4.0 * to_xy_xxyzz[k] * tke_0 + 4.0 * to_xxxy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xxy_xyyy[k] = 6.0 * to_xy_xyy[k] - 4.0 * to_xy_xyyyy[k] * tke_0 - 6.0 * to_xxxy_xyy[k] * tbe_0 + 4.0 * to_xxxy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxy_xyyz[k] = 4.0 * to_xy_xyz[k] - 4.0 * to_xy_xyyyz[k] * tke_0 - 4.0 * to_xxxy_xyz[k] * tbe_0 + 4.0 * to_xxxy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxy_xyzz[k] = 2.0 * to_xy_xzz[k] - 4.0 * to_xy_xyyzz[k] * tke_0 - 2.0 * to_xxxy_xzz[k] * tbe_0 + 4.0 * to_xxxy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxy_xzzz[k] = -4.0 * to_xy_xyzzz[k] * tke_0 + 4.0 * to_xxxy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxy_yyyy[k] = 8.0 * to_xy_yyy[k] - 4.0 * to_xy_yyyyy[k] * tke_0 - 8.0 * to_xxxy_yyy[k] * tbe_0 + 4.0 * to_xxxy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxy_yyyz[k] = 6.0 * to_xy_yyz[k] - 4.0 * to_xy_yyyyz[k] * tke_0 - 6.0 * to_xxxy_yyz[k] * tbe_0 + 4.0 * to_xxxy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxy_yyzz[k] = 4.0 * to_xy_yzz[k] - 4.0 * to_xy_yyyzz[k] * tke_0 - 4.0 * to_xxxy_yzz[k] * tbe_0 + 4.0 * to_xxxy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxy_yzzz[k] = 2.0 * to_xy_zzz[k] - 4.0 * to_xy_yyzzz[k] * tke_0 - 2.0 * to_xxxy_zzz[k] * tbe_0 + 4.0 * to_xxxy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxy_zzzz[k] = -4.0 * to_xy_yzzzz[k] * tke_0 + 4.0 * to_xxxy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-195 components of targeted buffer : FG

        auto to_x_y_xxz_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 30);

        auto to_x_y_xxz_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 31);

        auto to_x_y_xxz_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 32);

        auto to_x_y_xxz_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 33);

        auto to_x_y_xxz_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 34);

        auto to_x_y_xxz_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 35);

        auto to_x_y_xxz_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 36);

        auto to_x_y_xxz_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 37);

        auto to_x_y_xxz_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 38);

        auto to_x_y_xxz_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 39);

        auto to_x_y_xxz_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 40);

        auto to_x_y_xxz_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 41);

        auto to_x_y_xxz_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 42);

        auto to_x_y_xxz_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 43);

        auto to_x_y_xxz_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_x_y_xxz_xxxx, to_x_y_xxz_xxxy, to_x_y_xxz_xxxz, to_x_y_xxz_xxyy, to_x_y_xxz_xxyz, to_x_y_xxz_xxzz, to_x_y_xxz_xyyy, to_x_y_xxz_xyyz, to_x_y_xxz_xyzz, to_x_y_xxz_xzzz, to_x_y_xxz_yyyy, to_x_y_xxz_yyyz, to_x_y_xxz_yyzz, to_x_y_xxz_yzzz, to_x_y_xxz_zzzz, to_xxxz_xxx, to_xxxz_xxxxy, to_xxxz_xxxyy, to_xxxz_xxxyz, to_xxxz_xxy, to_xxxz_xxyyy, to_xxxz_xxyyz, to_xxxz_xxyzz, to_xxxz_xxz, to_xxxz_xyy, to_xxxz_xyyyy, to_xxxz_xyyyz, to_xxxz_xyyzz, to_xxxz_xyz, to_xxxz_xyzzz, to_xxxz_xzz, to_xxxz_yyy, to_xxxz_yyyyy, to_xxxz_yyyyz, to_xxxz_yyyzz, to_xxxz_yyz, to_xxxz_yyzzz, to_xxxz_yzz, to_xxxz_yzzzz, to_xxxz_zzz, to_xz_xxx, to_xz_xxxxy, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_yyy, to_xz_yyyyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxz_xxxx[k] = -4.0 * to_xz_xxxxy[k] * tke_0 + 4.0 * to_xxxz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xxz_xxxy[k] = 2.0 * to_xz_xxx[k] - 4.0 * to_xz_xxxyy[k] * tke_0 - 2.0 * to_xxxz_xxx[k] * tbe_0 + 4.0 * to_xxxz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xxz_xxxz[k] = -4.0 * to_xz_xxxyz[k] * tke_0 + 4.0 * to_xxxz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xxz_xxyy[k] = 4.0 * to_xz_xxy[k] - 4.0 * to_xz_xxyyy[k] * tke_0 - 4.0 * to_xxxz_xxy[k] * tbe_0 + 4.0 * to_xxxz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xxz_xxyz[k] = 2.0 * to_xz_xxz[k] - 4.0 * to_xz_xxyyz[k] * tke_0 - 2.0 * to_xxxz_xxz[k] * tbe_0 + 4.0 * to_xxxz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xxz_xxzz[k] = -4.0 * to_xz_xxyzz[k] * tke_0 + 4.0 * to_xxxz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xxz_xyyy[k] = 6.0 * to_xz_xyy[k] - 4.0 * to_xz_xyyyy[k] * tke_0 - 6.0 * to_xxxz_xyy[k] * tbe_0 + 4.0 * to_xxxz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxz_xyyz[k] = 4.0 * to_xz_xyz[k] - 4.0 * to_xz_xyyyz[k] * tke_0 - 4.0 * to_xxxz_xyz[k] * tbe_0 + 4.0 * to_xxxz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxz_xyzz[k] = 2.0 * to_xz_xzz[k] - 4.0 * to_xz_xyyzz[k] * tke_0 - 2.0 * to_xxxz_xzz[k] * tbe_0 + 4.0 * to_xxxz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxz_xzzz[k] = -4.0 * to_xz_xyzzz[k] * tke_0 + 4.0 * to_xxxz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxz_yyyy[k] = 8.0 * to_xz_yyy[k] - 4.0 * to_xz_yyyyy[k] * tke_0 - 8.0 * to_xxxz_yyy[k] * tbe_0 + 4.0 * to_xxxz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xxz_yyyz[k] = 6.0 * to_xz_yyz[k] - 4.0 * to_xz_yyyyz[k] * tke_0 - 6.0 * to_xxxz_yyz[k] * tbe_0 + 4.0 * to_xxxz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xxz_yyzz[k] = 4.0 * to_xz_yzz[k] - 4.0 * to_xz_yyyzz[k] * tke_0 - 4.0 * to_xxxz_yzz[k] * tbe_0 + 4.0 * to_xxxz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xxz_yzzz[k] = 2.0 * to_xz_zzz[k] - 4.0 * to_xz_yyzzz[k] * tke_0 - 2.0 * to_xxxz_zzz[k] * tbe_0 + 4.0 * to_xxxz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xxz_zzzz[k] = -4.0 * to_xz_yzzzz[k] * tke_0 + 4.0 * to_xxxz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 195-210 components of targeted buffer : FG

        auto to_x_y_xyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 45);

        auto to_x_y_xyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 46);

        auto to_x_y_xyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 47);

        auto to_x_y_xyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 48);

        auto to_x_y_xyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 49);

        auto to_x_y_xyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 50);

        auto to_x_y_xyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 51);

        auto to_x_y_xyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 52);

        auto to_x_y_xyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 53);

        auto to_x_y_xyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 54);

        auto to_x_y_xyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 55);

        auto to_x_y_xyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 56);

        auto to_x_y_xyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 57);

        auto to_x_y_xyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 58);

        auto to_x_y_xyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_x_y_xyy_xxxx, to_x_y_xyy_xxxy, to_x_y_xyy_xxxz, to_x_y_xyy_xxyy, to_x_y_xyy_xxyz, to_x_y_xyy_xxzz, to_x_y_xyy_xyyy, to_x_y_xyy_xyyz, to_x_y_xyy_xyzz, to_x_y_xyy_xzzz, to_x_y_xyy_yyyy, to_x_y_xyy_yyyz, to_x_y_xyy_yyzz, to_x_y_xyy_yzzz, to_x_y_xyy_zzzz, to_xxyy_xxx, to_xxyy_xxxxy, to_xxyy_xxxyy, to_xxyy_xxxyz, to_xxyy_xxy, to_xxyy_xxyyy, to_xxyy_xxyyz, to_xxyy_xxyzz, to_xxyy_xxz, to_xxyy_xyy, to_xxyy_xyyyy, to_xxyy_xyyyz, to_xxyy_xyyzz, to_xxyy_xyz, to_xxyy_xyzzz, to_xxyy_xzz, to_xxyy_yyy, to_xxyy_yyyyy, to_xxyy_yyyyz, to_xxyy_yyyzz, to_xxyy_yyz, to_xxyy_yyzzz, to_xxyy_yzz, to_xxyy_yzzzz, to_xxyy_zzz, to_yy_xxx, to_yy_xxxxy, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_yyy, to_yy_yyyyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyy_xxxx[k] = -2.0 * to_yy_xxxxy[k] * tke_0 + 4.0 * to_xxyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xyy_xxxy[k] = to_yy_xxx[k] - 2.0 * to_yy_xxxyy[k] * tke_0 - 2.0 * to_xxyy_xxx[k] * tbe_0 + 4.0 * to_xxyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xyy_xxxz[k] = -2.0 * to_yy_xxxyz[k] * tke_0 + 4.0 * to_xxyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xyy_xxyy[k] = 2.0 * to_yy_xxy[k] - 2.0 * to_yy_xxyyy[k] * tke_0 - 4.0 * to_xxyy_xxy[k] * tbe_0 + 4.0 * to_xxyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xyy_xxyz[k] = to_yy_xxz[k] - 2.0 * to_yy_xxyyz[k] * tke_0 - 2.0 * to_xxyy_xxz[k] * tbe_0 + 4.0 * to_xxyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xyy_xxzz[k] = -2.0 * to_yy_xxyzz[k] * tke_0 + 4.0 * to_xxyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xyy_xyyy[k] = 3.0 * to_yy_xyy[k] - 2.0 * to_yy_xyyyy[k] * tke_0 - 6.0 * to_xxyy_xyy[k] * tbe_0 + 4.0 * to_xxyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyy_xyyz[k] = 2.0 * to_yy_xyz[k] - 2.0 * to_yy_xyyyz[k] * tke_0 - 4.0 * to_xxyy_xyz[k] * tbe_0 + 4.0 * to_xxyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyy_xyzz[k] = to_yy_xzz[k] - 2.0 * to_yy_xyyzz[k] * tke_0 - 2.0 * to_xxyy_xzz[k] * tbe_0 + 4.0 * to_xxyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyy_xzzz[k] = -2.0 * to_yy_xyzzz[k] * tke_0 + 4.0 * to_xxyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyy_yyyy[k] = 4.0 * to_yy_yyy[k] - 2.0 * to_yy_yyyyy[k] * tke_0 - 8.0 * to_xxyy_yyy[k] * tbe_0 + 4.0 * to_xxyy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyy_yyyz[k] = 3.0 * to_yy_yyz[k] - 2.0 * to_yy_yyyyz[k] * tke_0 - 6.0 * to_xxyy_yyz[k] * tbe_0 + 4.0 * to_xxyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyy_yyzz[k] = 2.0 * to_yy_yzz[k] - 2.0 * to_yy_yyyzz[k] * tke_0 - 4.0 * to_xxyy_yzz[k] * tbe_0 + 4.0 * to_xxyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyy_yzzz[k] = to_yy_zzz[k] - 2.0 * to_yy_yyzzz[k] * tke_0 - 2.0 * to_xxyy_zzz[k] * tbe_0 + 4.0 * to_xxyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyy_zzzz[k] = -2.0 * to_yy_yzzzz[k] * tke_0 + 4.0 * to_xxyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-225 components of targeted buffer : FG

        auto to_x_y_xyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 60);

        auto to_x_y_xyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 61);

        auto to_x_y_xyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 62);

        auto to_x_y_xyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 63);

        auto to_x_y_xyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 64);

        auto to_x_y_xyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 65);

        auto to_x_y_xyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 66);

        auto to_x_y_xyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 67);

        auto to_x_y_xyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 68);

        auto to_x_y_xyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 69);

        auto to_x_y_xyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 70);

        auto to_x_y_xyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 71);

        auto to_x_y_xyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 72);

        auto to_x_y_xyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 73);

        auto to_x_y_xyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_x_y_xyz_xxxx, to_x_y_xyz_xxxy, to_x_y_xyz_xxxz, to_x_y_xyz_xxyy, to_x_y_xyz_xxyz, to_x_y_xyz_xxzz, to_x_y_xyz_xyyy, to_x_y_xyz_xyyz, to_x_y_xyz_xyzz, to_x_y_xyz_xzzz, to_x_y_xyz_yyyy, to_x_y_xyz_yyyz, to_x_y_xyz_yyzz, to_x_y_xyz_yzzz, to_x_y_xyz_zzzz, to_xxyz_xxx, to_xxyz_xxxxy, to_xxyz_xxxyy, to_xxyz_xxxyz, to_xxyz_xxy, to_xxyz_xxyyy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xyy, to_xxyz_xyyyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_yyy, to_xxyz_yyyyy, to_xxyz_yyyyz, to_xxyz_yyyzz, to_xxyz_yyz, to_xxyz_yyzzz, to_xxyz_yzz, to_xxyz_yzzzz, to_xxyz_zzz, to_yz_xxx, to_yz_xxxxy, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_yyy, to_yz_yyyyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyz_xxxx[k] = -2.0 * to_yz_xxxxy[k] * tke_0 + 4.0 * to_xxyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xyz_xxxy[k] = to_yz_xxx[k] - 2.0 * to_yz_xxxyy[k] * tke_0 - 2.0 * to_xxyz_xxx[k] * tbe_0 + 4.0 * to_xxyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xyz_xxxz[k] = -2.0 * to_yz_xxxyz[k] * tke_0 + 4.0 * to_xxyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xyz_xxyy[k] = 2.0 * to_yz_xxy[k] - 2.0 * to_yz_xxyyy[k] * tke_0 - 4.0 * to_xxyz_xxy[k] * tbe_0 + 4.0 * to_xxyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xyz_xxyz[k] = to_yz_xxz[k] - 2.0 * to_yz_xxyyz[k] * tke_0 - 2.0 * to_xxyz_xxz[k] * tbe_0 + 4.0 * to_xxyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xyz_xxzz[k] = -2.0 * to_yz_xxyzz[k] * tke_0 + 4.0 * to_xxyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xyz_xyyy[k] = 3.0 * to_yz_xyy[k] - 2.0 * to_yz_xyyyy[k] * tke_0 - 6.0 * to_xxyz_xyy[k] * tbe_0 + 4.0 * to_xxyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyz_xyyz[k] = 2.0 * to_yz_xyz[k] - 2.0 * to_yz_xyyyz[k] * tke_0 - 4.0 * to_xxyz_xyz[k] * tbe_0 + 4.0 * to_xxyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyz_xyzz[k] = to_yz_xzz[k] - 2.0 * to_yz_xyyzz[k] * tke_0 - 2.0 * to_xxyz_xzz[k] * tbe_0 + 4.0 * to_xxyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyz_xzzz[k] = -2.0 * to_yz_xyzzz[k] * tke_0 + 4.0 * to_xxyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyz_yyyy[k] = 4.0 * to_yz_yyy[k] - 2.0 * to_yz_yyyyy[k] * tke_0 - 8.0 * to_xxyz_yyy[k] * tbe_0 + 4.0 * to_xxyz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xyz_yyyz[k] = 3.0 * to_yz_yyz[k] - 2.0 * to_yz_yyyyz[k] * tke_0 - 6.0 * to_xxyz_yyz[k] * tbe_0 + 4.0 * to_xxyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xyz_yyzz[k] = 2.0 * to_yz_yzz[k] - 2.0 * to_yz_yyyzz[k] * tke_0 - 4.0 * to_xxyz_yzz[k] * tbe_0 + 4.0 * to_xxyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xyz_yzzz[k] = to_yz_zzz[k] - 2.0 * to_yz_yyzzz[k] * tke_0 - 2.0 * to_xxyz_zzz[k] * tbe_0 + 4.0 * to_xxyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xyz_zzzz[k] = -2.0 * to_yz_yzzzz[k] * tke_0 + 4.0 * to_xxyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 225-240 components of targeted buffer : FG

        auto to_x_y_xzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 75);

        auto to_x_y_xzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 76);

        auto to_x_y_xzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 77);

        auto to_x_y_xzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 78);

        auto to_x_y_xzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 79);

        auto to_x_y_xzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 80);

        auto to_x_y_xzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 81);

        auto to_x_y_xzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 82);

        auto to_x_y_xzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 83);

        auto to_x_y_xzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 84);

        auto to_x_y_xzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 85);

        auto to_x_y_xzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 86);

        auto to_x_y_xzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 87);

        auto to_x_y_xzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 88);

        auto to_x_y_xzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_x_y_xzz_xxxx, to_x_y_xzz_xxxy, to_x_y_xzz_xxxz, to_x_y_xzz_xxyy, to_x_y_xzz_xxyz, to_x_y_xzz_xxzz, to_x_y_xzz_xyyy, to_x_y_xzz_xyyz, to_x_y_xzz_xyzz, to_x_y_xzz_xzzz, to_x_y_xzz_yyyy, to_x_y_xzz_yyyz, to_x_y_xzz_yyzz, to_x_y_xzz_yzzz, to_x_y_xzz_zzzz, to_xxzz_xxx, to_xxzz_xxxxy, to_xxzz_xxxyy, to_xxzz_xxxyz, to_xxzz_xxy, to_xxzz_xxyyy, to_xxzz_xxyyz, to_xxzz_xxyzz, to_xxzz_xxz, to_xxzz_xyy, to_xxzz_xyyyy, to_xxzz_xyyyz, to_xxzz_xyyzz, to_xxzz_xyz, to_xxzz_xyzzz, to_xxzz_xzz, to_xxzz_yyy, to_xxzz_yyyyy, to_xxzz_yyyyz, to_xxzz_yyyzz, to_xxzz_yyz, to_xxzz_yyzzz, to_xxzz_yzz, to_xxzz_yzzzz, to_xxzz_zzz, to_zz_xxx, to_zz_xxxxy, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_yyy, to_zz_yyyyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xzz_xxxx[k] = -2.0 * to_zz_xxxxy[k] * tke_0 + 4.0 * to_xxzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_xzz_xxxy[k] = to_zz_xxx[k] - 2.0 * to_zz_xxxyy[k] * tke_0 - 2.0 * to_xxzz_xxx[k] * tbe_0 + 4.0 * to_xxzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_xzz_xxxz[k] = -2.0 * to_zz_xxxyz[k] * tke_0 + 4.0 * to_xxzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_xzz_xxyy[k] = 2.0 * to_zz_xxy[k] - 2.0 * to_zz_xxyyy[k] * tke_0 - 4.0 * to_xxzz_xxy[k] * tbe_0 + 4.0 * to_xxzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_xzz_xxyz[k] = to_zz_xxz[k] - 2.0 * to_zz_xxyyz[k] * tke_0 - 2.0 * to_xxzz_xxz[k] * tbe_0 + 4.0 * to_xxzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_xzz_xxzz[k] = -2.0 * to_zz_xxyzz[k] * tke_0 + 4.0 * to_xxzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_xzz_xyyy[k] = 3.0 * to_zz_xyy[k] - 2.0 * to_zz_xyyyy[k] * tke_0 - 6.0 * to_xxzz_xyy[k] * tbe_0 + 4.0 * to_xxzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_xzz_xyyz[k] = 2.0 * to_zz_xyz[k] - 2.0 * to_zz_xyyyz[k] * tke_0 - 4.0 * to_xxzz_xyz[k] * tbe_0 + 4.0 * to_xxzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_xzz_xyzz[k] = to_zz_xzz[k] - 2.0 * to_zz_xyyzz[k] * tke_0 - 2.0 * to_xxzz_xzz[k] * tbe_0 + 4.0 * to_xxzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_xzz_xzzz[k] = -2.0 * to_zz_xyzzz[k] * tke_0 + 4.0 * to_xxzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_xzz_yyyy[k] = 4.0 * to_zz_yyy[k] - 2.0 * to_zz_yyyyy[k] * tke_0 - 8.0 * to_xxzz_yyy[k] * tbe_0 + 4.0 * to_xxzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_xzz_yyyz[k] = 3.0 * to_zz_yyz[k] - 2.0 * to_zz_yyyyz[k] * tke_0 - 6.0 * to_xxzz_yyz[k] * tbe_0 + 4.0 * to_xxzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_xzz_yyzz[k] = 2.0 * to_zz_yzz[k] - 2.0 * to_zz_yyyzz[k] * tke_0 - 4.0 * to_xxzz_yzz[k] * tbe_0 + 4.0 * to_xxzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_xzz_yzzz[k] = to_zz_zzz[k] - 2.0 * to_zz_yyzzz[k] * tke_0 - 2.0 * to_xxzz_zzz[k] * tbe_0 + 4.0 * to_xxzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_xzz_zzzz[k] = -2.0 * to_zz_yzzzz[k] * tke_0 + 4.0 * to_xxzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-255 components of targeted buffer : FG

        auto to_x_y_yyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 90);

        auto to_x_y_yyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 91);

        auto to_x_y_yyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 92);

        auto to_x_y_yyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 93);

        auto to_x_y_yyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 94);

        auto to_x_y_yyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 95);

        auto to_x_y_yyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 96);

        auto to_x_y_yyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 97);

        auto to_x_y_yyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 98);

        auto to_x_y_yyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 99);

        auto to_x_y_yyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 100);

        auto to_x_y_yyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 101);

        auto to_x_y_yyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 102);

        auto to_x_y_yyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 103);

        auto to_x_y_yyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_x_y_yyy_xxxx, to_x_y_yyy_xxxy, to_x_y_yyy_xxxz, to_x_y_yyy_xxyy, to_x_y_yyy_xxyz, to_x_y_yyy_xxzz, to_x_y_yyy_xyyy, to_x_y_yyy_xyyz, to_x_y_yyy_xyzz, to_x_y_yyy_xzzz, to_x_y_yyy_yyyy, to_x_y_yyy_yyyz, to_x_y_yyy_yyzz, to_x_y_yyy_yzzz, to_x_y_yyy_zzzz, to_xyyy_xxx, to_xyyy_xxxxy, to_xyyy_xxxyy, to_xyyy_xxxyz, to_xyyy_xxy, to_xyyy_xxyyy, to_xyyy_xxyyz, to_xyyy_xxyzz, to_xyyy_xxz, to_xyyy_xyy, to_xyyy_xyyyy, to_xyyy_xyyyz, to_xyyy_xyyzz, to_xyyy_xyz, to_xyyy_xyzzz, to_xyyy_xzz, to_xyyy_yyy, to_xyyy_yyyyy, to_xyyy_yyyyz, to_xyyy_yyyzz, to_xyyy_yyz, to_xyyy_yyzzz, to_xyyy_yzz, to_xyyy_yzzzz, to_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyy_xxxx[k] = 4.0 * to_xyyy_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_yyy_xxxy[k] = -2.0 * to_xyyy_xxx[k] * tbe_0 + 4.0 * to_xyyy_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_yyy_xxxz[k] = 4.0 * to_xyyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_yyy_xxyy[k] = -4.0 * to_xyyy_xxy[k] * tbe_0 + 4.0 * to_xyyy_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_yyy_xxyz[k] = -2.0 * to_xyyy_xxz[k] * tbe_0 + 4.0 * to_xyyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_yyy_xxzz[k] = 4.0 * to_xyyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_yyy_xyyy[k] = -6.0 * to_xyyy_xyy[k] * tbe_0 + 4.0 * to_xyyy_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyy_xyyz[k] = -4.0 * to_xyyy_xyz[k] * tbe_0 + 4.0 * to_xyyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyy_xyzz[k] = -2.0 * to_xyyy_xzz[k] * tbe_0 + 4.0 * to_xyyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyy_xzzz[k] = 4.0 * to_xyyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyy_yyyy[k] = -8.0 * to_xyyy_yyy[k] * tbe_0 + 4.0 * to_xyyy_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyy_yyyz[k] = -6.0 * to_xyyy_yyz[k] * tbe_0 + 4.0 * to_xyyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyy_yyzz[k] = -4.0 * to_xyyy_yzz[k] * tbe_0 + 4.0 * to_xyyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyy_yzzz[k] = -2.0 * to_xyyy_zzz[k] * tbe_0 + 4.0 * to_xyyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyy_zzzz[k] = 4.0 * to_xyyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 255-270 components of targeted buffer : FG

        auto to_x_y_yyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 105);

        auto to_x_y_yyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 106);

        auto to_x_y_yyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 107);

        auto to_x_y_yyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 108);

        auto to_x_y_yyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 109);

        auto to_x_y_yyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 110);

        auto to_x_y_yyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 111);

        auto to_x_y_yyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 112);

        auto to_x_y_yyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 113);

        auto to_x_y_yyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 114);

        auto to_x_y_yyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 115);

        auto to_x_y_yyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 116);

        auto to_x_y_yyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 117);

        auto to_x_y_yyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 118);

        auto to_x_y_yyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_x_y_yyz_xxxx, to_x_y_yyz_xxxy, to_x_y_yyz_xxxz, to_x_y_yyz_xxyy, to_x_y_yyz_xxyz, to_x_y_yyz_xxzz, to_x_y_yyz_xyyy, to_x_y_yyz_xyyz, to_x_y_yyz_xyzz, to_x_y_yyz_xzzz, to_x_y_yyz_yyyy, to_x_y_yyz_yyyz, to_x_y_yyz_yyzz, to_x_y_yyz_yzzz, to_x_y_yyz_zzzz, to_xyyz_xxx, to_xyyz_xxxxy, to_xyyz_xxxyy, to_xyyz_xxxyz, to_xyyz_xxy, to_xyyz_xxyyy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xyy, to_xyyz_xyyyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_yyy, to_xyyz_yyyyy, to_xyyz_yyyyz, to_xyyz_yyyzz, to_xyyz_yyz, to_xyyz_yyzzz, to_xyyz_yzz, to_xyyz_yzzzz, to_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyz_xxxx[k] = 4.0 * to_xyyz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_yyz_xxxy[k] = -2.0 * to_xyyz_xxx[k] * tbe_0 + 4.0 * to_xyyz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_yyz_xxxz[k] = 4.0 * to_xyyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_yyz_xxyy[k] = -4.0 * to_xyyz_xxy[k] * tbe_0 + 4.0 * to_xyyz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_yyz_xxyz[k] = -2.0 * to_xyyz_xxz[k] * tbe_0 + 4.0 * to_xyyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_yyz_xxzz[k] = 4.0 * to_xyyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_yyz_xyyy[k] = -6.0 * to_xyyz_xyy[k] * tbe_0 + 4.0 * to_xyyz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyz_xyyz[k] = -4.0 * to_xyyz_xyz[k] * tbe_0 + 4.0 * to_xyyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyz_xyzz[k] = -2.0 * to_xyyz_xzz[k] * tbe_0 + 4.0 * to_xyyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyz_xzzz[k] = 4.0 * to_xyyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyz_yyyy[k] = -8.0 * to_xyyz_yyy[k] * tbe_0 + 4.0 * to_xyyz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_yyz_yyyz[k] = -6.0 * to_xyyz_yyz[k] * tbe_0 + 4.0 * to_xyyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_yyz_yyzz[k] = -4.0 * to_xyyz_yzz[k] * tbe_0 + 4.0 * to_xyyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_yyz_yzzz[k] = -2.0 * to_xyyz_zzz[k] * tbe_0 + 4.0 * to_xyyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_yyz_zzzz[k] = 4.0 * to_xyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-285 components of targeted buffer : FG

        auto to_x_y_yzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 120);

        auto to_x_y_yzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 121);

        auto to_x_y_yzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 122);

        auto to_x_y_yzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 123);

        auto to_x_y_yzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 124);

        auto to_x_y_yzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 125);

        auto to_x_y_yzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 126);

        auto to_x_y_yzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 127);

        auto to_x_y_yzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 128);

        auto to_x_y_yzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 129);

        auto to_x_y_yzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 130);

        auto to_x_y_yzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 131);

        auto to_x_y_yzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 132);

        auto to_x_y_yzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 133);

        auto to_x_y_yzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_x_y_yzz_xxxx, to_x_y_yzz_xxxy, to_x_y_yzz_xxxz, to_x_y_yzz_xxyy, to_x_y_yzz_xxyz, to_x_y_yzz_xxzz, to_x_y_yzz_xyyy, to_x_y_yzz_xyyz, to_x_y_yzz_xyzz, to_x_y_yzz_xzzz, to_x_y_yzz_yyyy, to_x_y_yzz_yyyz, to_x_y_yzz_yyzz, to_x_y_yzz_yzzz, to_x_y_yzz_zzzz, to_xyzz_xxx, to_xyzz_xxxxy, to_xyzz_xxxyy, to_xyzz_xxxyz, to_xyzz_xxy, to_xyzz_xxyyy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xyy, to_xyzz_xyyyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_yyy, to_xyzz_yyyyy, to_xyzz_yyyyz, to_xyzz_yyyzz, to_xyzz_yyz, to_xyzz_yyzzz, to_xyzz_yzz, to_xyzz_yzzzz, to_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yzz_xxxx[k] = 4.0 * to_xyzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_yzz_xxxy[k] = -2.0 * to_xyzz_xxx[k] * tbe_0 + 4.0 * to_xyzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_yzz_xxxz[k] = 4.0 * to_xyzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_yzz_xxyy[k] = -4.0 * to_xyzz_xxy[k] * tbe_0 + 4.0 * to_xyzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_yzz_xxyz[k] = -2.0 * to_xyzz_xxz[k] * tbe_0 + 4.0 * to_xyzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_yzz_xxzz[k] = 4.0 * to_xyzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_yzz_xyyy[k] = -6.0 * to_xyzz_xyy[k] * tbe_0 + 4.0 * to_xyzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_yzz_xyyz[k] = -4.0 * to_xyzz_xyz[k] * tbe_0 + 4.0 * to_xyzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_yzz_xyzz[k] = -2.0 * to_xyzz_xzz[k] * tbe_0 + 4.0 * to_xyzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_yzz_xzzz[k] = 4.0 * to_xyzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_yzz_yyyy[k] = -8.0 * to_xyzz_yyy[k] * tbe_0 + 4.0 * to_xyzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_yzz_yyyz[k] = -6.0 * to_xyzz_yyz[k] * tbe_0 + 4.0 * to_xyzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_yzz_yyzz[k] = -4.0 * to_xyzz_yzz[k] * tbe_0 + 4.0 * to_xyzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_yzz_yzzz[k] = -2.0 * to_xyzz_zzz[k] * tbe_0 + 4.0 * to_xyzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_yzz_zzzz[k] = 4.0 * to_xyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 285-300 components of targeted buffer : FG

        auto to_x_y_zzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 135);

        auto to_x_y_zzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 136);

        auto to_x_y_zzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 137);

        auto to_x_y_zzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 138);

        auto to_x_y_zzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 139);

        auto to_x_y_zzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 140);

        auto to_x_y_zzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 141);

        auto to_x_y_zzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 142);

        auto to_x_y_zzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 143);

        auto to_x_y_zzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 144);

        auto to_x_y_zzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 145);

        auto to_x_y_zzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 146);

        auto to_x_y_zzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 147);

        auto to_x_y_zzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 148);

        auto to_x_y_zzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 1 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_x_y_zzz_xxxx, to_x_y_zzz_xxxy, to_x_y_zzz_xxxz, to_x_y_zzz_xxyy, to_x_y_zzz_xxyz, to_x_y_zzz_xxzz, to_x_y_zzz_xyyy, to_x_y_zzz_xyyz, to_x_y_zzz_xyzz, to_x_y_zzz_xzzz, to_x_y_zzz_yyyy, to_x_y_zzz_yyyz, to_x_y_zzz_yyzz, to_x_y_zzz_yzzz, to_x_y_zzz_zzzz, to_xzzz_xxx, to_xzzz_xxxxy, to_xzzz_xxxyy, to_xzzz_xxxyz, to_xzzz_xxy, to_xzzz_xxyyy, to_xzzz_xxyyz, to_xzzz_xxyzz, to_xzzz_xxz, to_xzzz_xyy, to_xzzz_xyyyy, to_xzzz_xyyyz, to_xzzz_xyyzz, to_xzzz_xyz, to_xzzz_xyzzz, to_xzzz_xzz, to_xzzz_yyy, to_xzzz_yyyyy, to_xzzz_yyyyz, to_xzzz_yyyzz, to_xzzz_yyz, to_xzzz_yyzzz, to_xzzz_yzz, to_xzzz_yzzzz, to_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zzz_xxxx[k] = 4.0 * to_xzzz_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_zzz_xxxy[k] = -2.0 * to_xzzz_xxx[k] * tbe_0 + 4.0 * to_xzzz_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_zzz_xxxz[k] = 4.0 * to_xzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_zzz_xxyy[k] = -4.0 * to_xzzz_xxy[k] * tbe_0 + 4.0 * to_xzzz_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_zzz_xxyz[k] = -2.0 * to_xzzz_xxz[k] * tbe_0 + 4.0 * to_xzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_zzz_xxzz[k] = 4.0 * to_xzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_zzz_xyyy[k] = -6.0 * to_xzzz_xyy[k] * tbe_0 + 4.0 * to_xzzz_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_zzz_xyyz[k] = -4.0 * to_xzzz_xyz[k] * tbe_0 + 4.0 * to_xzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_zzz_xyzz[k] = -2.0 * to_xzzz_xzz[k] * tbe_0 + 4.0 * to_xzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_zzz_xzzz[k] = 4.0 * to_xzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_zzz_yyyy[k] = -8.0 * to_xzzz_yyy[k] * tbe_0 + 4.0 * to_xzzz_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_zzz_yyyz[k] = -6.0 * to_xzzz_yyz[k] * tbe_0 + 4.0 * to_xzzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_zzz_yyzz[k] = -4.0 * to_xzzz_yzz[k] * tbe_0 + 4.0 * to_xzzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_zzz_yzzz[k] = -2.0 * to_xzzz_zzz[k] * tbe_0 + 4.0 * to_xzzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_zzz_zzzz[k] = 4.0 * to_xzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-315 components of targeted buffer : FG

        auto to_x_z_xxx_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 0);

        auto to_x_z_xxx_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 1);

        auto to_x_z_xxx_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 2);

        auto to_x_z_xxx_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 3);

        auto to_x_z_xxx_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 4);

        auto to_x_z_xxx_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 5);

        auto to_x_z_xxx_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 6);

        auto to_x_z_xxx_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 7);

        auto to_x_z_xxx_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 8);

        auto to_x_z_xxx_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 9);

        auto to_x_z_xxx_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 10);

        auto to_x_z_xxx_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 11);

        auto to_x_z_xxx_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 12);

        auto to_x_z_xxx_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 13);

        auto to_x_z_xxx_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_x_z_xxx_xxxx, to_x_z_xxx_xxxy, to_x_z_xxx_xxxz, to_x_z_xxx_xxyy, to_x_z_xxx_xxyz, to_x_z_xxx_xxzz, to_x_z_xxx_xyyy, to_x_z_xxx_xyyz, to_x_z_xxx_xyzz, to_x_z_xxx_xzzz, to_x_z_xxx_yyyy, to_x_z_xxx_yyyz, to_x_z_xxx_yyzz, to_x_z_xxx_yzzz, to_x_z_xxx_zzzz, to_xx_xxx, to_xx_xxxxz, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, to_xx_zzzzz, to_xxxx_xxx, to_xxxx_xxxxz, to_xxxx_xxxyz, to_xxxx_xxxzz, to_xxxx_xxy, to_xxxx_xxyyz, to_xxxx_xxyzz, to_xxxx_xxz, to_xxxx_xxzzz, to_xxxx_xyy, to_xxxx_xyyyz, to_xxxx_xyyzz, to_xxxx_xyz, to_xxxx_xyzzz, to_xxxx_xzz, to_xxxx_xzzzz, to_xxxx_yyy, to_xxxx_yyyyz, to_xxxx_yyyzz, to_xxxx_yyz, to_xxxx_yyzzz, to_xxxx_yzz, to_xxxx_yzzzz, to_xxxx_zzz, to_xxxx_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxx_xxxx[k] = -6.0 * to_xx_xxxxz[k] * tke_0 + 4.0 * to_xxxx_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xxxy[k] = -6.0 * to_xx_xxxyz[k] * tke_0 + 4.0 * to_xxxx_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xxxz[k] = 3.0 * to_xx_xxx[k] - 6.0 * to_xx_xxxzz[k] * tke_0 - 2.0 * to_xxxx_xxx[k] * tbe_0 + 4.0 * to_xxxx_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xxyy[k] = -6.0 * to_xx_xxyyz[k] * tke_0 + 4.0 * to_xxxx_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xxyz[k] = 3.0 * to_xx_xxy[k] - 6.0 * to_xx_xxyzz[k] * tke_0 - 2.0 * to_xxxx_xxy[k] * tbe_0 + 4.0 * to_xxxx_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xxzz[k] = 6.0 * to_xx_xxz[k] - 6.0 * to_xx_xxzzz[k] * tke_0 - 4.0 * to_xxxx_xxz[k] * tbe_0 + 4.0 * to_xxxx_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xyyy[k] = -6.0 * to_xx_xyyyz[k] * tke_0 + 4.0 * to_xxxx_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xyyz[k] = 3.0 * to_xx_xyy[k] - 6.0 * to_xx_xyyzz[k] * tke_0 - 2.0 * to_xxxx_xyy[k] * tbe_0 + 4.0 * to_xxxx_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xyzz[k] = 6.0 * to_xx_xyz[k] - 6.0 * to_xx_xyzzz[k] * tke_0 - 4.0 * to_xxxx_xyz[k] * tbe_0 + 4.0 * to_xxxx_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xzzz[k] = 9.0 * to_xx_xzz[k] - 6.0 * to_xx_xzzzz[k] * tke_0 - 6.0 * to_xxxx_xzz[k] * tbe_0 + 4.0 * to_xxxx_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_yyyy[k] = -6.0 * to_xx_yyyyz[k] * tke_0 + 4.0 * to_xxxx_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxx_yyyz[k] = 3.0 * to_xx_yyy[k] - 6.0 * to_xx_yyyzz[k] * tke_0 - 2.0 * to_xxxx_yyy[k] * tbe_0 + 4.0 * to_xxxx_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_yyzz[k] = 6.0 * to_xx_yyz[k] - 6.0 * to_xx_yyzzz[k] * tke_0 - 4.0 * to_xxxx_yyz[k] * tbe_0 + 4.0 * to_xxxx_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_yzzz[k] = 9.0 * to_xx_yzz[k] - 6.0 * to_xx_yzzzz[k] * tke_0 - 6.0 * to_xxxx_yzz[k] * tbe_0 + 4.0 * to_xxxx_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_zzzz[k] = 12.0 * to_xx_zzz[k] - 6.0 * to_xx_zzzzz[k] * tke_0 - 8.0 * to_xxxx_zzz[k] * tbe_0 + 4.0 * to_xxxx_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 315-330 components of targeted buffer : FG

        auto to_x_z_xxy_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 15);

        auto to_x_z_xxy_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 16);

        auto to_x_z_xxy_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 17);

        auto to_x_z_xxy_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 18);

        auto to_x_z_xxy_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 19);

        auto to_x_z_xxy_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 20);

        auto to_x_z_xxy_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 21);

        auto to_x_z_xxy_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 22);

        auto to_x_z_xxy_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 23);

        auto to_x_z_xxy_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 24);

        auto to_x_z_xxy_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 25);

        auto to_x_z_xxy_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 26);

        auto to_x_z_xxy_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 27);

        auto to_x_z_xxy_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 28);

        auto to_x_z_xxy_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_x_z_xxy_xxxx, to_x_z_xxy_xxxy, to_x_z_xxy_xxxz, to_x_z_xxy_xxyy, to_x_z_xxy_xxyz, to_x_z_xxy_xxzz, to_x_z_xxy_xyyy, to_x_z_xxy_xyyz, to_x_z_xxy_xyzz, to_x_z_xxy_xzzz, to_x_z_xxy_yyyy, to_x_z_xxy_yyyz, to_x_z_xxy_yyzz, to_x_z_xxy_yzzz, to_x_z_xxy_zzzz, to_xxxy_xxx, to_xxxy_xxxxz, to_xxxy_xxxyz, to_xxxy_xxxzz, to_xxxy_xxy, to_xxxy_xxyyz, to_xxxy_xxyzz, to_xxxy_xxz, to_xxxy_xxzzz, to_xxxy_xyy, to_xxxy_xyyyz, to_xxxy_xyyzz, to_xxxy_xyz, to_xxxy_xyzzz, to_xxxy_xzz, to_xxxy_xzzzz, to_xxxy_yyy, to_xxxy_yyyyz, to_xxxy_yyyzz, to_xxxy_yyz, to_xxxy_yyzzz, to_xxxy_yzz, to_xxxy_yzzzz, to_xxxy_zzz, to_xxxy_zzzzz, to_xy_xxx, to_xy_xxxxz, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, to_xy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxy_xxxx[k] = -4.0 * to_xy_xxxxz[k] * tke_0 + 4.0 * to_xxxy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xxxy[k] = -4.0 * to_xy_xxxyz[k] * tke_0 + 4.0 * to_xxxy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xxxz[k] = 2.0 * to_xy_xxx[k] - 4.0 * to_xy_xxxzz[k] * tke_0 - 2.0 * to_xxxy_xxx[k] * tbe_0 + 4.0 * to_xxxy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xxyy[k] = -4.0 * to_xy_xxyyz[k] * tke_0 + 4.0 * to_xxxy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xxyz[k] = 2.0 * to_xy_xxy[k] - 4.0 * to_xy_xxyzz[k] * tke_0 - 2.0 * to_xxxy_xxy[k] * tbe_0 + 4.0 * to_xxxy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xxzz[k] = 4.0 * to_xy_xxz[k] - 4.0 * to_xy_xxzzz[k] * tke_0 - 4.0 * to_xxxy_xxz[k] * tbe_0 + 4.0 * to_xxxy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xyyy[k] = -4.0 * to_xy_xyyyz[k] * tke_0 + 4.0 * to_xxxy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xyyz[k] = 2.0 * to_xy_xyy[k] - 4.0 * to_xy_xyyzz[k] * tke_0 - 2.0 * to_xxxy_xyy[k] * tbe_0 + 4.0 * to_xxxy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xyzz[k] = 4.0 * to_xy_xyz[k] - 4.0 * to_xy_xyzzz[k] * tke_0 - 4.0 * to_xxxy_xyz[k] * tbe_0 + 4.0 * to_xxxy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xzzz[k] = 6.0 * to_xy_xzz[k] - 4.0 * to_xy_xzzzz[k] * tke_0 - 6.0 * to_xxxy_xzz[k] * tbe_0 + 4.0 * to_xxxy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_yyyy[k] = -4.0 * to_xy_yyyyz[k] * tke_0 + 4.0 * to_xxxy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxy_yyyz[k] = 2.0 * to_xy_yyy[k] - 4.0 * to_xy_yyyzz[k] * tke_0 - 2.0 * to_xxxy_yyy[k] * tbe_0 + 4.0 * to_xxxy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_yyzz[k] = 4.0 * to_xy_yyz[k] - 4.0 * to_xy_yyzzz[k] * tke_0 - 4.0 * to_xxxy_yyz[k] * tbe_0 + 4.0 * to_xxxy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_yzzz[k] = 6.0 * to_xy_yzz[k] - 4.0 * to_xy_yzzzz[k] * tke_0 - 6.0 * to_xxxy_yzz[k] * tbe_0 + 4.0 * to_xxxy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_zzzz[k] = 8.0 * to_xy_zzz[k] - 4.0 * to_xy_zzzzz[k] * tke_0 - 8.0 * to_xxxy_zzz[k] * tbe_0 + 4.0 * to_xxxy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 330-345 components of targeted buffer : FG

        auto to_x_z_xxz_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 30);

        auto to_x_z_xxz_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 31);

        auto to_x_z_xxz_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 32);

        auto to_x_z_xxz_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 33);

        auto to_x_z_xxz_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 34);

        auto to_x_z_xxz_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 35);

        auto to_x_z_xxz_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 36);

        auto to_x_z_xxz_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 37);

        auto to_x_z_xxz_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 38);

        auto to_x_z_xxz_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 39);

        auto to_x_z_xxz_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 40);

        auto to_x_z_xxz_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 41);

        auto to_x_z_xxz_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 42);

        auto to_x_z_xxz_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 43);

        auto to_x_z_xxz_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_x_z_xxz_xxxx, to_x_z_xxz_xxxy, to_x_z_xxz_xxxz, to_x_z_xxz_xxyy, to_x_z_xxz_xxyz, to_x_z_xxz_xxzz, to_x_z_xxz_xyyy, to_x_z_xxz_xyyz, to_x_z_xxz_xyzz, to_x_z_xxz_xzzz, to_x_z_xxz_yyyy, to_x_z_xxz_yyyz, to_x_z_xxz_yyzz, to_x_z_xxz_yzzz, to_x_z_xxz_zzzz, to_xxxz_xxx, to_xxxz_xxxxz, to_xxxz_xxxyz, to_xxxz_xxxzz, to_xxxz_xxy, to_xxxz_xxyyz, to_xxxz_xxyzz, to_xxxz_xxz, to_xxxz_xxzzz, to_xxxz_xyy, to_xxxz_xyyyz, to_xxxz_xyyzz, to_xxxz_xyz, to_xxxz_xyzzz, to_xxxz_xzz, to_xxxz_xzzzz, to_xxxz_yyy, to_xxxz_yyyyz, to_xxxz_yyyzz, to_xxxz_yyz, to_xxxz_yyzzz, to_xxxz_yzz, to_xxxz_yzzzz, to_xxxz_zzz, to_xxxz_zzzzz, to_xz_xxx, to_xz_xxxxz, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, to_xz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxz_xxxx[k] = -4.0 * to_xz_xxxxz[k] * tke_0 + 4.0 * to_xxxz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xxxy[k] = -4.0 * to_xz_xxxyz[k] * tke_0 + 4.0 * to_xxxz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xxxz[k] = 2.0 * to_xz_xxx[k] - 4.0 * to_xz_xxxzz[k] * tke_0 - 2.0 * to_xxxz_xxx[k] * tbe_0 + 4.0 * to_xxxz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xxyy[k] = -4.0 * to_xz_xxyyz[k] * tke_0 + 4.0 * to_xxxz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xxyz[k] = 2.0 * to_xz_xxy[k] - 4.0 * to_xz_xxyzz[k] * tke_0 - 2.0 * to_xxxz_xxy[k] * tbe_0 + 4.0 * to_xxxz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xxzz[k] = 4.0 * to_xz_xxz[k] - 4.0 * to_xz_xxzzz[k] * tke_0 - 4.0 * to_xxxz_xxz[k] * tbe_0 + 4.0 * to_xxxz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xyyy[k] = -4.0 * to_xz_xyyyz[k] * tke_0 + 4.0 * to_xxxz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xyyz[k] = 2.0 * to_xz_xyy[k] - 4.0 * to_xz_xyyzz[k] * tke_0 - 2.0 * to_xxxz_xyy[k] * tbe_0 + 4.0 * to_xxxz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xyzz[k] = 4.0 * to_xz_xyz[k] - 4.0 * to_xz_xyzzz[k] * tke_0 - 4.0 * to_xxxz_xyz[k] * tbe_0 + 4.0 * to_xxxz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xzzz[k] = 6.0 * to_xz_xzz[k] - 4.0 * to_xz_xzzzz[k] * tke_0 - 6.0 * to_xxxz_xzz[k] * tbe_0 + 4.0 * to_xxxz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_yyyy[k] = -4.0 * to_xz_yyyyz[k] * tke_0 + 4.0 * to_xxxz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xxz_yyyz[k] = 2.0 * to_xz_yyy[k] - 4.0 * to_xz_yyyzz[k] * tke_0 - 2.0 * to_xxxz_yyy[k] * tbe_0 + 4.0 * to_xxxz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_yyzz[k] = 4.0 * to_xz_yyz[k] - 4.0 * to_xz_yyzzz[k] * tke_0 - 4.0 * to_xxxz_yyz[k] * tbe_0 + 4.0 * to_xxxz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_yzzz[k] = 6.0 * to_xz_yzz[k] - 4.0 * to_xz_yzzzz[k] * tke_0 - 6.0 * to_xxxz_yzz[k] * tbe_0 + 4.0 * to_xxxz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_zzzz[k] = 8.0 * to_xz_zzz[k] - 4.0 * to_xz_zzzzz[k] * tke_0 - 8.0 * to_xxxz_zzz[k] * tbe_0 + 4.0 * to_xxxz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 345-360 components of targeted buffer : FG

        auto to_x_z_xyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 45);

        auto to_x_z_xyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 46);

        auto to_x_z_xyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 47);

        auto to_x_z_xyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 48);

        auto to_x_z_xyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 49);

        auto to_x_z_xyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 50);

        auto to_x_z_xyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 51);

        auto to_x_z_xyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 52);

        auto to_x_z_xyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 53);

        auto to_x_z_xyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 54);

        auto to_x_z_xyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 55);

        auto to_x_z_xyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 56);

        auto to_x_z_xyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 57);

        auto to_x_z_xyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 58);

        auto to_x_z_xyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_x_z_xyy_xxxx, to_x_z_xyy_xxxy, to_x_z_xyy_xxxz, to_x_z_xyy_xxyy, to_x_z_xyy_xxyz, to_x_z_xyy_xxzz, to_x_z_xyy_xyyy, to_x_z_xyy_xyyz, to_x_z_xyy_xyzz, to_x_z_xyy_xzzz, to_x_z_xyy_yyyy, to_x_z_xyy_yyyz, to_x_z_xyy_yyzz, to_x_z_xyy_yzzz, to_x_z_xyy_zzzz, to_xxyy_xxx, to_xxyy_xxxxz, to_xxyy_xxxyz, to_xxyy_xxxzz, to_xxyy_xxy, to_xxyy_xxyyz, to_xxyy_xxyzz, to_xxyy_xxz, to_xxyy_xxzzz, to_xxyy_xyy, to_xxyy_xyyyz, to_xxyy_xyyzz, to_xxyy_xyz, to_xxyy_xyzzz, to_xxyy_xzz, to_xxyy_xzzzz, to_xxyy_yyy, to_xxyy_yyyyz, to_xxyy_yyyzz, to_xxyy_yyz, to_xxyy_yyzzz, to_xxyy_yzz, to_xxyy_yzzzz, to_xxyy_zzz, to_xxyy_zzzzz, to_yy_xxx, to_yy_xxxxz, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, to_yy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyy_xxxx[k] = -2.0 * to_yy_xxxxz[k] * tke_0 + 4.0 * to_xxyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xxxy[k] = -2.0 * to_yy_xxxyz[k] * tke_0 + 4.0 * to_xxyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xxxz[k] = to_yy_xxx[k] - 2.0 * to_yy_xxxzz[k] * tke_0 - 2.0 * to_xxyy_xxx[k] * tbe_0 + 4.0 * to_xxyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xxyy[k] = -2.0 * to_yy_xxyyz[k] * tke_0 + 4.0 * to_xxyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xxyz[k] = to_yy_xxy[k] - 2.0 * to_yy_xxyzz[k] * tke_0 - 2.0 * to_xxyy_xxy[k] * tbe_0 + 4.0 * to_xxyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xxzz[k] = 2.0 * to_yy_xxz[k] - 2.0 * to_yy_xxzzz[k] * tke_0 - 4.0 * to_xxyy_xxz[k] * tbe_0 + 4.0 * to_xxyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xyyy[k] = -2.0 * to_yy_xyyyz[k] * tke_0 + 4.0 * to_xxyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xyyz[k] = to_yy_xyy[k] - 2.0 * to_yy_xyyzz[k] * tke_0 - 2.0 * to_xxyy_xyy[k] * tbe_0 + 4.0 * to_xxyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xyzz[k] = 2.0 * to_yy_xyz[k] - 2.0 * to_yy_xyzzz[k] * tke_0 - 4.0 * to_xxyy_xyz[k] * tbe_0 + 4.0 * to_xxyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xzzz[k] = 3.0 * to_yy_xzz[k] - 2.0 * to_yy_xzzzz[k] * tke_0 - 6.0 * to_xxyy_xzz[k] * tbe_0 + 4.0 * to_xxyy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_yyyy[k] = -2.0 * to_yy_yyyyz[k] * tke_0 + 4.0 * to_xxyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyy_yyyz[k] = to_yy_yyy[k] - 2.0 * to_yy_yyyzz[k] * tke_0 - 2.0 * to_xxyy_yyy[k] * tbe_0 + 4.0 * to_xxyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_yyzz[k] = 2.0 * to_yy_yyz[k] - 2.0 * to_yy_yyzzz[k] * tke_0 - 4.0 * to_xxyy_yyz[k] * tbe_0 + 4.0 * to_xxyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_yzzz[k] = 3.0 * to_yy_yzz[k] - 2.0 * to_yy_yzzzz[k] * tke_0 - 6.0 * to_xxyy_yzz[k] * tbe_0 + 4.0 * to_xxyy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_zzzz[k] = 4.0 * to_yy_zzz[k] - 2.0 * to_yy_zzzzz[k] * tke_0 - 8.0 * to_xxyy_zzz[k] * tbe_0 + 4.0 * to_xxyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 360-375 components of targeted buffer : FG

        auto to_x_z_xyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 60);

        auto to_x_z_xyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 61);

        auto to_x_z_xyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 62);

        auto to_x_z_xyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 63);

        auto to_x_z_xyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 64);

        auto to_x_z_xyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 65);

        auto to_x_z_xyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 66);

        auto to_x_z_xyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 67);

        auto to_x_z_xyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 68);

        auto to_x_z_xyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 69);

        auto to_x_z_xyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 70);

        auto to_x_z_xyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 71);

        auto to_x_z_xyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 72);

        auto to_x_z_xyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 73);

        auto to_x_z_xyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_x_z_xyz_xxxx, to_x_z_xyz_xxxy, to_x_z_xyz_xxxz, to_x_z_xyz_xxyy, to_x_z_xyz_xxyz, to_x_z_xyz_xxzz, to_x_z_xyz_xyyy, to_x_z_xyz_xyyz, to_x_z_xyz_xyzz, to_x_z_xyz_xzzz, to_x_z_xyz_yyyy, to_x_z_xyz_yyyz, to_x_z_xyz_yyzz, to_x_z_xyz_yzzz, to_x_z_xyz_zzzz, to_xxyz_xxx, to_xxyz_xxxxz, to_xxyz_xxxyz, to_xxyz_xxxzz, to_xxyz_xxy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xxzzz, to_xxyz_xyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_xzzzz, to_xxyz_yyy, to_xxyz_yyyyz, to_xxyz_yyyzz, to_xxyz_yyz, to_xxyz_yyzzz, to_xxyz_yzz, to_xxyz_yzzzz, to_xxyz_zzz, to_xxyz_zzzzz, to_yz_xxx, to_yz_xxxxz, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, to_yz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyz_xxxx[k] = -2.0 * to_yz_xxxxz[k] * tke_0 + 4.0 * to_xxyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xxxy[k] = -2.0 * to_yz_xxxyz[k] * tke_0 + 4.0 * to_xxyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xxxz[k] = to_yz_xxx[k] - 2.0 * to_yz_xxxzz[k] * tke_0 - 2.0 * to_xxyz_xxx[k] * tbe_0 + 4.0 * to_xxyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xxyy[k] = -2.0 * to_yz_xxyyz[k] * tke_0 + 4.0 * to_xxyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xxyz[k] = to_yz_xxy[k] - 2.0 * to_yz_xxyzz[k] * tke_0 - 2.0 * to_xxyz_xxy[k] * tbe_0 + 4.0 * to_xxyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xxzz[k] = 2.0 * to_yz_xxz[k] - 2.0 * to_yz_xxzzz[k] * tke_0 - 4.0 * to_xxyz_xxz[k] * tbe_0 + 4.0 * to_xxyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xyyy[k] = -2.0 * to_yz_xyyyz[k] * tke_0 + 4.0 * to_xxyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xyyz[k] = to_yz_xyy[k] - 2.0 * to_yz_xyyzz[k] * tke_0 - 2.0 * to_xxyz_xyy[k] * tbe_0 + 4.0 * to_xxyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xyzz[k] = 2.0 * to_yz_xyz[k] - 2.0 * to_yz_xyzzz[k] * tke_0 - 4.0 * to_xxyz_xyz[k] * tbe_0 + 4.0 * to_xxyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xzzz[k] = 3.0 * to_yz_xzz[k] - 2.0 * to_yz_xzzzz[k] * tke_0 - 6.0 * to_xxyz_xzz[k] * tbe_0 + 4.0 * to_xxyz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_yyyy[k] = -2.0 * to_yz_yyyyz[k] * tke_0 + 4.0 * to_xxyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xyz_yyyz[k] = to_yz_yyy[k] - 2.0 * to_yz_yyyzz[k] * tke_0 - 2.0 * to_xxyz_yyy[k] * tbe_0 + 4.0 * to_xxyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_yyzz[k] = 2.0 * to_yz_yyz[k] - 2.0 * to_yz_yyzzz[k] * tke_0 - 4.0 * to_xxyz_yyz[k] * tbe_0 + 4.0 * to_xxyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_yzzz[k] = 3.0 * to_yz_yzz[k] - 2.0 * to_yz_yzzzz[k] * tke_0 - 6.0 * to_xxyz_yzz[k] * tbe_0 + 4.0 * to_xxyz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_zzzz[k] = 4.0 * to_yz_zzz[k] - 2.0 * to_yz_zzzzz[k] * tke_0 - 8.0 * to_xxyz_zzz[k] * tbe_0 + 4.0 * to_xxyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 375-390 components of targeted buffer : FG

        auto to_x_z_xzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 75);

        auto to_x_z_xzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 76);

        auto to_x_z_xzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 77);

        auto to_x_z_xzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 78);

        auto to_x_z_xzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 79);

        auto to_x_z_xzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 80);

        auto to_x_z_xzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 81);

        auto to_x_z_xzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 82);

        auto to_x_z_xzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 83);

        auto to_x_z_xzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 84);

        auto to_x_z_xzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 85);

        auto to_x_z_xzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 86);

        auto to_x_z_xzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 87);

        auto to_x_z_xzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 88);

        auto to_x_z_xzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_x_z_xzz_xxxx, to_x_z_xzz_xxxy, to_x_z_xzz_xxxz, to_x_z_xzz_xxyy, to_x_z_xzz_xxyz, to_x_z_xzz_xxzz, to_x_z_xzz_xyyy, to_x_z_xzz_xyyz, to_x_z_xzz_xyzz, to_x_z_xzz_xzzz, to_x_z_xzz_yyyy, to_x_z_xzz_yyyz, to_x_z_xzz_yyzz, to_x_z_xzz_yzzz, to_x_z_xzz_zzzz, to_xxzz_xxx, to_xxzz_xxxxz, to_xxzz_xxxyz, to_xxzz_xxxzz, to_xxzz_xxy, to_xxzz_xxyyz, to_xxzz_xxyzz, to_xxzz_xxz, to_xxzz_xxzzz, to_xxzz_xyy, to_xxzz_xyyyz, to_xxzz_xyyzz, to_xxzz_xyz, to_xxzz_xyzzz, to_xxzz_xzz, to_xxzz_xzzzz, to_xxzz_yyy, to_xxzz_yyyyz, to_xxzz_yyyzz, to_xxzz_yyz, to_xxzz_yyzzz, to_xxzz_yzz, to_xxzz_yzzzz, to_xxzz_zzz, to_xxzz_zzzzz, to_zz_xxx, to_zz_xxxxz, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, to_zz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xzz_xxxx[k] = -2.0 * to_zz_xxxxz[k] * tke_0 + 4.0 * to_xxzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xxxy[k] = -2.0 * to_zz_xxxyz[k] * tke_0 + 4.0 * to_xxzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xxxz[k] = to_zz_xxx[k] - 2.0 * to_zz_xxxzz[k] * tke_0 - 2.0 * to_xxzz_xxx[k] * tbe_0 + 4.0 * to_xxzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xxyy[k] = -2.0 * to_zz_xxyyz[k] * tke_0 + 4.0 * to_xxzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xxyz[k] = to_zz_xxy[k] - 2.0 * to_zz_xxyzz[k] * tke_0 - 2.0 * to_xxzz_xxy[k] * tbe_0 + 4.0 * to_xxzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xxzz[k] = 2.0 * to_zz_xxz[k] - 2.0 * to_zz_xxzzz[k] * tke_0 - 4.0 * to_xxzz_xxz[k] * tbe_0 + 4.0 * to_xxzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xyyy[k] = -2.0 * to_zz_xyyyz[k] * tke_0 + 4.0 * to_xxzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xyyz[k] = to_zz_xyy[k] - 2.0 * to_zz_xyyzz[k] * tke_0 - 2.0 * to_xxzz_xyy[k] * tbe_0 + 4.0 * to_xxzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xyzz[k] = 2.0 * to_zz_xyz[k] - 2.0 * to_zz_xyzzz[k] * tke_0 - 4.0 * to_xxzz_xyz[k] * tbe_0 + 4.0 * to_xxzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xzzz[k] = 3.0 * to_zz_xzz[k] - 2.0 * to_zz_xzzzz[k] * tke_0 - 6.0 * to_xxzz_xzz[k] * tbe_0 + 4.0 * to_xxzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_yyyy[k] = -2.0 * to_zz_yyyyz[k] * tke_0 + 4.0 * to_xxzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_xzz_yyyz[k] = to_zz_yyy[k] - 2.0 * to_zz_yyyzz[k] * tke_0 - 2.0 * to_xxzz_yyy[k] * tbe_0 + 4.0 * to_xxzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_yyzz[k] = 2.0 * to_zz_yyz[k] - 2.0 * to_zz_yyzzz[k] * tke_0 - 4.0 * to_xxzz_yyz[k] * tbe_0 + 4.0 * to_xxzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_yzzz[k] = 3.0 * to_zz_yzz[k] - 2.0 * to_zz_yzzzz[k] * tke_0 - 6.0 * to_xxzz_yzz[k] * tbe_0 + 4.0 * to_xxzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_zzzz[k] = 4.0 * to_zz_zzz[k] - 2.0 * to_zz_zzzzz[k] * tke_0 - 8.0 * to_xxzz_zzz[k] * tbe_0 + 4.0 * to_xxzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 390-405 components of targeted buffer : FG

        auto to_x_z_yyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 90);

        auto to_x_z_yyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 91);

        auto to_x_z_yyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 92);

        auto to_x_z_yyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 93);

        auto to_x_z_yyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 94);

        auto to_x_z_yyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 95);

        auto to_x_z_yyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 96);

        auto to_x_z_yyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 97);

        auto to_x_z_yyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 98);

        auto to_x_z_yyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 99);

        auto to_x_z_yyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 100);

        auto to_x_z_yyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 101);

        auto to_x_z_yyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 102);

        auto to_x_z_yyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 103);

        auto to_x_z_yyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_x_z_yyy_xxxx, to_x_z_yyy_xxxy, to_x_z_yyy_xxxz, to_x_z_yyy_xxyy, to_x_z_yyy_xxyz, to_x_z_yyy_xxzz, to_x_z_yyy_xyyy, to_x_z_yyy_xyyz, to_x_z_yyy_xyzz, to_x_z_yyy_xzzz, to_x_z_yyy_yyyy, to_x_z_yyy_yyyz, to_x_z_yyy_yyzz, to_x_z_yyy_yzzz, to_x_z_yyy_zzzz, to_xyyy_xxx, to_xyyy_xxxxz, to_xyyy_xxxyz, to_xyyy_xxxzz, to_xyyy_xxy, to_xyyy_xxyyz, to_xyyy_xxyzz, to_xyyy_xxz, to_xyyy_xxzzz, to_xyyy_xyy, to_xyyy_xyyyz, to_xyyy_xyyzz, to_xyyy_xyz, to_xyyy_xyzzz, to_xyyy_xzz, to_xyyy_xzzzz, to_xyyy_yyy, to_xyyy_yyyyz, to_xyyy_yyyzz, to_xyyy_yyz, to_xyyy_yyzzz, to_xyyy_yzz, to_xyyy_yzzzz, to_xyyy_zzz, to_xyyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyy_xxxx[k] = 4.0 * to_xyyy_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xxxy[k] = 4.0 * to_xyyy_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xxxz[k] = -2.0 * to_xyyy_xxx[k] * tbe_0 + 4.0 * to_xyyy_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xxyy[k] = 4.0 * to_xyyy_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xxyz[k] = -2.0 * to_xyyy_xxy[k] * tbe_0 + 4.0 * to_xyyy_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xxzz[k] = -4.0 * to_xyyy_xxz[k] * tbe_0 + 4.0 * to_xyyy_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xyyy[k] = 4.0 * to_xyyy_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xyyz[k] = -2.0 * to_xyyy_xyy[k] * tbe_0 + 4.0 * to_xyyy_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xyzz[k] = -4.0 * to_xyyy_xyz[k] * tbe_0 + 4.0 * to_xyyy_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xzzz[k] = -6.0 * to_xyyy_xzz[k] * tbe_0 + 4.0 * to_xyyy_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_yyyy[k] = 4.0 * to_xyyy_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyy_yyyz[k] = -2.0 * to_xyyy_yyy[k] * tbe_0 + 4.0 * to_xyyy_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_yyzz[k] = -4.0 * to_xyyy_yyz[k] * tbe_0 + 4.0 * to_xyyy_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_yzzz[k] = -6.0 * to_xyyy_yzz[k] * tbe_0 + 4.0 * to_xyyy_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_zzzz[k] = -8.0 * to_xyyy_zzz[k] * tbe_0 + 4.0 * to_xyyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 405-420 components of targeted buffer : FG

        auto to_x_z_yyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 105);

        auto to_x_z_yyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 106);

        auto to_x_z_yyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 107);

        auto to_x_z_yyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 108);

        auto to_x_z_yyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 109);

        auto to_x_z_yyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 110);

        auto to_x_z_yyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 111);

        auto to_x_z_yyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 112);

        auto to_x_z_yyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 113);

        auto to_x_z_yyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 114);

        auto to_x_z_yyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 115);

        auto to_x_z_yyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 116);

        auto to_x_z_yyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 117);

        auto to_x_z_yyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 118);

        auto to_x_z_yyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_x_z_yyz_xxxx, to_x_z_yyz_xxxy, to_x_z_yyz_xxxz, to_x_z_yyz_xxyy, to_x_z_yyz_xxyz, to_x_z_yyz_xxzz, to_x_z_yyz_xyyy, to_x_z_yyz_xyyz, to_x_z_yyz_xyzz, to_x_z_yyz_xzzz, to_x_z_yyz_yyyy, to_x_z_yyz_yyyz, to_x_z_yyz_yyzz, to_x_z_yyz_yzzz, to_x_z_yyz_zzzz, to_xyyz_xxx, to_xyyz_xxxxz, to_xyyz_xxxyz, to_xyyz_xxxzz, to_xyyz_xxy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xxzzz, to_xyyz_xyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_xzzzz, to_xyyz_yyy, to_xyyz_yyyyz, to_xyyz_yyyzz, to_xyyz_yyz, to_xyyz_yyzzz, to_xyyz_yzz, to_xyyz_yzzzz, to_xyyz_zzz, to_xyyz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyz_xxxx[k] = 4.0 * to_xyyz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xxxy[k] = 4.0 * to_xyyz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xxxz[k] = -2.0 * to_xyyz_xxx[k] * tbe_0 + 4.0 * to_xyyz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xxyy[k] = 4.0 * to_xyyz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xxyz[k] = -2.0 * to_xyyz_xxy[k] * tbe_0 + 4.0 * to_xyyz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xxzz[k] = -4.0 * to_xyyz_xxz[k] * tbe_0 + 4.0 * to_xyyz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xyyy[k] = 4.0 * to_xyyz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xyyz[k] = -2.0 * to_xyyz_xyy[k] * tbe_0 + 4.0 * to_xyyz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xyzz[k] = -4.0 * to_xyyz_xyz[k] * tbe_0 + 4.0 * to_xyyz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xzzz[k] = -6.0 * to_xyyz_xzz[k] * tbe_0 + 4.0 * to_xyyz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_yyyy[k] = 4.0 * to_xyyz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_yyz_yyyz[k] = -2.0 * to_xyyz_yyy[k] * tbe_0 + 4.0 * to_xyyz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_yyzz[k] = -4.0 * to_xyyz_yyz[k] * tbe_0 + 4.0 * to_xyyz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_yzzz[k] = -6.0 * to_xyyz_yzz[k] * tbe_0 + 4.0 * to_xyyz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_zzzz[k] = -8.0 * to_xyyz_zzz[k] * tbe_0 + 4.0 * to_xyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 420-435 components of targeted buffer : FG

        auto to_x_z_yzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 120);

        auto to_x_z_yzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 121);

        auto to_x_z_yzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 122);

        auto to_x_z_yzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 123);

        auto to_x_z_yzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 124);

        auto to_x_z_yzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 125);

        auto to_x_z_yzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 126);

        auto to_x_z_yzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 127);

        auto to_x_z_yzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 128);

        auto to_x_z_yzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 129);

        auto to_x_z_yzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 130);

        auto to_x_z_yzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 131);

        auto to_x_z_yzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 132);

        auto to_x_z_yzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 133);

        auto to_x_z_yzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_x_z_yzz_xxxx, to_x_z_yzz_xxxy, to_x_z_yzz_xxxz, to_x_z_yzz_xxyy, to_x_z_yzz_xxyz, to_x_z_yzz_xxzz, to_x_z_yzz_xyyy, to_x_z_yzz_xyyz, to_x_z_yzz_xyzz, to_x_z_yzz_xzzz, to_x_z_yzz_yyyy, to_x_z_yzz_yyyz, to_x_z_yzz_yyzz, to_x_z_yzz_yzzz, to_x_z_yzz_zzzz, to_xyzz_xxx, to_xyzz_xxxxz, to_xyzz_xxxyz, to_xyzz_xxxzz, to_xyzz_xxy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xxzzz, to_xyzz_xyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_xzzzz, to_xyzz_yyy, to_xyzz_yyyyz, to_xyzz_yyyzz, to_xyzz_yyz, to_xyzz_yyzzz, to_xyzz_yzz, to_xyzz_yzzzz, to_xyzz_zzz, to_xyzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yzz_xxxx[k] = 4.0 * to_xyzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xxxy[k] = 4.0 * to_xyzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xxxz[k] = -2.0 * to_xyzz_xxx[k] * tbe_0 + 4.0 * to_xyzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xxyy[k] = 4.0 * to_xyzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xxyz[k] = -2.0 * to_xyzz_xxy[k] * tbe_0 + 4.0 * to_xyzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xxzz[k] = -4.0 * to_xyzz_xxz[k] * tbe_0 + 4.0 * to_xyzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xyyy[k] = 4.0 * to_xyzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xyyz[k] = -2.0 * to_xyzz_xyy[k] * tbe_0 + 4.0 * to_xyzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xyzz[k] = -4.0 * to_xyzz_xyz[k] * tbe_0 + 4.0 * to_xyzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xzzz[k] = -6.0 * to_xyzz_xzz[k] * tbe_0 + 4.0 * to_xyzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_yyyy[k] = 4.0 * to_xyzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_yzz_yyyz[k] = -2.0 * to_xyzz_yyy[k] * tbe_0 + 4.0 * to_xyzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_yyzz[k] = -4.0 * to_xyzz_yyz[k] * tbe_0 + 4.0 * to_xyzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_yzzz[k] = -6.0 * to_xyzz_yzz[k] * tbe_0 + 4.0 * to_xyzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_zzzz[k] = -8.0 * to_xyzz_zzz[k] * tbe_0 + 4.0 * to_xyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 435-450 components of targeted buffer : FG

        auto to_x_z_zzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 135);

        auto to_x_z_zzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 136);

        auto to_x_z_zzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 137);

        auto to_x_z_zzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 138);

        auto to_x_z_zzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 139);

        auto to_x_z_zzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 140);

        auto to_x_z_zzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 141);

        auto to_x_z_zzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 142);

        auto to_x_z_zzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 143);

        auto to_x_z_zzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 144);

        auto to_x_z_zzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 145);

        auto to_x_z_zzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 146);

        auto to_x_z_zzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 147);

        auto to_x_z_zzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 148);

        auto to_x_z_zzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 2 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_x_z_zzz_xxxx, to_x_z_zzz_xxxy, to_x_z_zzz_xxxz, to_x_z_zzz_xxyy, to_x_z_zzz_xxyz, to_x_z_zzz_xxzz, to_x_z_zzz_xyyy, to_x_z_zzz_xyyz, to_x_z_zzz_xyzz, to_x_z_zzz_xzzz, to_x_z_zzz_yyyy, to_x_z_zzz_yyyz, to_x_z_zzz_yyzz, to_x_z_zzz_yzzz, to_x_z_zzz_zzzz, to_xzzz_xxx, to_xzzz_xxxxz, to_xzzz_xxxyz, to_xzzz_xxxzz, to_xzzz_xxy, to_xzzz_xxyyz, to_xzzz_xxyzz, to_xzzz_xxz, to_xzzz_xxzzz, to_xzzz_xyy, to_xzzz_xyyyz, to_xzzz_xyyzz, to_xzzz_xyz, to_xzzz_xyzzz, to_xzzz_xzz, to_xzzz_xzzzz, to_xzzz_yyy, to_xzzz_yyyyz, to_xzzz_yyyzz, to_xzzz_yyz, to_xzzz_yyzzz, to_xzzz_yzz, to_xzzz_yzzzz, to_xzzz_zzz, to_xzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zzz_xxxx[k] = 4.0 * to_xzzz_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xxxy[k] = 4.0 * to_xzzz_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xxxz[k] = -2.0 * to_xzzz_xxx[k] * tbe_0 + 4.0 * to_xzzz_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xxyy[k] = 4.0 * to_xzzz_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xxyz[k] = -2.0 * to_xzzz_xxy[k] * tbe_0 + 4.0 * to_xzzz_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xxzz[k] = -4.0 * to_xzzz_xxz[k] * tbe_0 + 4.0 * to_xzzz_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xyyy[k] = 4.0 * to_xzzz_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xyyz[k] = -2.0 * to_xzzz_xyy[k] * tbe_0 + 4.0 * to_xzzz_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xyzz[k] = -4.0 * to_xzzz_xyz[k] * tbe_0 + 4.0 * to_xzzz_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xzzz[k] = -6.0 * to_xzzz_xzz[k] * tbe_0 + 4.0 * to_xzzz_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_yyyy[k] = 4.0 * to_xzzz_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_zzz_yyyz[k] = -2.0 * to_xzzz_yyy[k] * tbe_0 + 4.0 * to_xzzz_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_yyzz[k] = -4.0 * to_xzzz_yyz[k] * tbe_0 + 4.0 * to_xzzz_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_yzzz[k] = -6.0 * to_xzzz_yzz[k] * tbe_0 + 4.0 * to_xzzz_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_zzzz[k] = -8.0 * to_xzzz_zzz[k] * tbe_0 + 4.0 * to_xzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 450-465 components of targeted buffer : FG

        auto to_y_x_xxx_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 0);

        auto to_y_x_xxx_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 1);

        auto to_y_x_xxx_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 2);

        auto to_y_x_xxx_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 3);

        auto to_y_x_xxx_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 4);

        auto to_y_x_xxx_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 5);

        auto to_y_x_xxx_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 6);

        auto to_y_x_xxx_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 7);

        auto to_y_x_xxx_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 8);

        auto to_y_x_xxx_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 9);

        auto to_y_x_xxx_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 10);

        auto to_y_x_xxx_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 11);

        auto to_y_x_xxx_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 12);

        auto to_y_x_xxx_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 13);

        auto to_y_x_xxx_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_xxxy_xxx, to_xxxy_xxxxx, to_xxxy_xxxxy, to_xxxy_xxxxz, to_xxxy_xxxyy, to_xxxy_xxxyz, to_xxxy_xxxzz, to_xxxy_xxy, to_xxxy_xxyyy, to_xxxy_xxyyz, to_xxxy_xxyzz, to_xxxy_xxz, to_xxxy_xxzzz, to_xxxy_xyy, to_xxxy_xyyyy, to_xxxy_xyyyz, to_xxxy_xyyzz, to_xxxy_xyz, to_xxxy_xyzzz, to_xxxy_xzz, to_xxxy_xzzzz, to_xxxy_yyy, to_xxxy_yyz, to_xxxy_yzz, to_xxxy_zzz, to_y_x_xxx_xxxx, to_y_x_xxx_xxxy, to_y_x_xxx_xxxz, to_y_x_xxx_xxyy, to_y_x_xxx_xxyz, to_y_x_xxx_xxzz, to_y_x_xxx_xyyy, to_y_x_xxx_xyyz, to_y_x_xxx_xyzz, to_y_x_xxx_xzzz, to_y_x_xxx_yyyy, to_y_x_xxx_yyyz, to_y_x_xxx_yyzz, to_y_x_xxx_yzzz, to_y_x_xxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxx_xxxx[k] = -8.0 * to_xxxy_xxx[k] * tbe_0 + 4.0 * to_xxxy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xxx_xxxy[k] = -6.0 * to_xxxy_xxy[k] * tbe_0 + 4.0 * to_xxxy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xxx_xxxz[k] = -6.0 * to_xxxy_xxz[k] * tbe_0 + 4.0 * to_xxxy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xxx_xxyy[k] = -4.0 * to_xxxy_xyy[k] * tbe_0 + 4.0 * to_xxxy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xxx_xxyz[k] = -4.0 * to_xxxy_xyz[k] * tbe_0 + 4.0 * to_xxxy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xxx_xxzz[k] = -4.0 * to_xxxy_xzz[k] * tbe_0 + 4.0 * to_xxxy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xxx_xyyy[k] = -2.0 * to_xxxy_yyy[k] * tbe_0 + 4.0 * to_xxxy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xxx_xyyz[k] = -2.0 * to_xxxy_yyz[k] * tbe_0 + 4.0 * to_xxxy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xxx_xyzz[k] = -2.0 * to_xxxy_yzz[k] * tbe_0 + 4.0 * to_xxxy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xxx_xzzz[k] = -2.0 * to_xxxy_zzz[k] * tbe_0 + 4.0 * to_xxxy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xxx_yyyy[k] = 4.0 * to_xxxy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xxx_yyyz[k] = 4.0 * to_xxxy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xxx_yyzz[k] = 4.0 * to_xxxy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xxx_yzzz[k] = 4.0 * to_xxxy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xxx_zzzz[k] = 4.0 * to_xxxy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 465-480 components of targeted buffer : FG

        auto to_y_x_xxy_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 15);

        auto to_y_x_xxy_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 16);

        auto to_y_x_xxy_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 17);

        auto to_y_x_xxy_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 18);

        auto to_y_x_xxy_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 19);

        auto to_y_x_xxy_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 20);

        auto to_y_x_xxy_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 21);

        auto to_y_x_xxy_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 22);

        auto to_y_x_xxy_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 23);

        auto to_y_x_xxy_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 24);

        auto to_y_x_xxy_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 25);

        auto to_y_x_xxy_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 26);

        auto to_y_x_xxy_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 27);

        auto to_y_x_xxy_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 28);

        auto to_y_x_xxy_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xx_xxx, to_xx_xxxxx, to_xx_xxxxy, to_xx_xxxxz, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyz, to_xx_yzz, to_xx_zzz, to_xxyy_xxx, to_xxyy_xxxxx, to_xxyy_xxxxy, to_xxyy_xxxxz, to_xxyy_xxxyy, to_xxyy_xxxyz, to_xxyy_xxxzz, to_xxyy_xxy, to_xxyy_xxyyy, to_xxyy_xxyyz, to_xxyy_xxyzz, to_xxyy_xxz, to_xxyy_xxzzz, to_xxyy_xyy, to_xxyy_xyyyy, to_xxyy_xyyyz, to_xxyy_xyyzz, to_xxyy_xyz, to_xxyy_xyzzz, to_xxyy_xzz, to_xxyy_xzzzz, to_xxyy_yyy, to_xxyy_yyz, to_xxyy_yzz, to_xxyy_zzz, to_y_x_xxy_xxxx, to_y_x_xxy_xxxy, to_y_x_xxy_xxxz, to_y_x_xxy_xxyy, to_y_x_xxy_xxyz, to_y_x_xxy_xxzz, to_y_x_xxy_xyyy, to_y_x_xxy_xyyz, to_y_x_xxy_xyzz, to_y_x_xxy_xzzz, to_y_x_xxy_yyyy, to_y_x_xxy_yyyz, to_y_x_xxy_yyzz, to_y_x_xxy_yzzz, to_y_x_xxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxy_xxxx[k] = 4.0 * to_xx_xxx[k] - 2.0 * to_xx_xxxxx[k] * tke_0 - 8.0 * to_xxyy_xxx[k] * tbe_0 + 4.0 * to_xxyy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xxy_xxxy[k] = 3.0 * to_xx_xxy[k] - 2.0 * to_xx_xxxxy[k] * tke_0 - 6.0 * to_xxyy_xxy[k] * tbe_0 + 4.0 * to_xxyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xxy_xxxz[k] = 3.0 * to_xx_xxz[k] - 2.0 * to_xx_xxxxz[k] * tke_0 - 6.0 * to_xxyy_xxz[k] * tbe_0 + 4.0 * to_xxyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xxy_xxyy[k] = 2.0 * to_xx_xyy[k] - 2.0 * to_xx_xxxyy[k] * tke_0 - 4.0 * to_xxyy_xyy[k] * tbe_0 + 4.0 * to_xxyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xxy_xxyz[k] = 2.0 * to_xx_xyz[k] - 2.0 * to_xx_xxxyz[k] * tke_0 - 4.0 * to_xxyy_xyz[k] * tbe_0 + 4.0 * to_xxyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xxy_xxzz[k] = 2.0 * to_xx_xzz[k] - 2.0 * to_xx_xxxzz[k] * tke_0 - 4.0 * to_xxyy_xzz[k] * tbe_0 + 4.0 * to_xxyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xxy_xyyy[k] = to_xx_yyy[k] - 2.0 * to_xx_xxyyy[k] * tke_0 - 2.0 * to_xxyy_yyy[k] * tbe_0 + 4.0 * to_xxyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xxy_xyyz[k] = to_xx_yyz[k] - 2.0 * to_xx_xxyyz[k] * tke_0 - 2.0 * to_xxyy_yyz[k] * tbe_0 + 4.0 * to_xxyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xxy_xyzz[k] = to_xx_yzz[k] - 2.0 * to_xx_xxyzz[k] * tke_0 - 2.0 * to_xxyy_yzz[k] * tbe_0 + 4.0 * to_xxyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xxy_xzzz[k] = to_xx_zzz[k] - 2.0 * to_xx_xxzzz[k] * tke_0 - 2.0 * to_xxyy_zzz[k] * tbe_0 + 4.0 * to_xxyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xxy_yyyy[k] = -2.0 * to_xx_xyyyy[k] * tke_0 + 4.0 * to_xxyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xxy_yyyz[k] = -2.0 * to_xx_xyyyz[k] * tke_0 + 4.0 * to_xxyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xxy_yyzz[k] = -2.0 * to_xx_xyyzz[k] * tke_0 + 4.0 * to_xxyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xxy_yzzz[k] = -2.0 * to_xx_xyzzz[k] * tke_0 + 4.0 * to_xxyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xxy_zzzz[k] = -2.0 * to_xx_xzzzz[k] * tke_0 + 4.0 * to_xxyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 480-495 components of targeted buffer : FG

        auto to_y_x_xxz_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 30);

        auto to_y_x_xxz_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 31);

        auto to_y_x_xxz_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 32);

        auto to_y_x_xxz_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 33);

        auto to_y_x_xxz_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 34);

        auto to_y_x_xxz_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 35);

        auto to_y_x_xxz_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 36);

        auto to_y_x_xxz_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 37);

        auto to_y_x_xxz_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 38);

        auto to_y_x_xxz_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 39);

        auto to_y_x_xxz_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 40);

        auto to_y_x_xxz_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 41);

        auto to_y_x_xxz_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 42);

        auto to_y_x_xxz_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 43);

        auto to_y_x_xxz_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_xxyz_xxx, to_xxyz_xxxxx, to_xxyz_xxxxy, to_xxyz_xxxxz, to_xxyz_xxxyy, to_xxyz_xxxyz, to_xxyz_xxxzz, to_xxyz_xxy, to_xxyz_xxyyy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xxzzz, to_xxyz_xyy, to_xxyz_xyyyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_xzzzz, to_xxyz_yyy, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_zzz, to_y_x_xxz_xxxx, to_y_x_xxz_xxxy, to_y_x_xxz_xxxz, to_y_x_xxz_xxyy, to_y_x_xxz_xxyz, to_y_x_xxz_xxzz, to_y_x_xxz_xyyy, to_y_x_xxz_xyyz, to_y_x_xxz_xyzz, to_y_x_xxz_xzzz, to_y_x_xxz_yyyy, to_y_x_xxz_yyyz, to_y_x_xxz_yyzz, to_y_x_xxz_yzzz, to_y_x_xxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxz_xxxx[k] = -8.0 * to_xxyz_xxx[k] * tbe_0 + 4.0 * to_xxyz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xxz_xxxy[k] = -6.0 * to_xxyz_xxy[k] * tbe_0 + 4.0 * to_xxyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xxz_xxxz[k] = -6.0 * to_xxyz_xxz[k] * tbe_0 + 4.0 * to_xxyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xxz_xxyy[k] = -4.0 * to_xxyz_xyy[k] * tbe_0 + 4.0 * to_xxyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xxz_xxyz[k] = -4.0 * to_xxyz_xyz[k] * tbe_0 + 4.0 * to_xxyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xxz_xxzz[k] = -4.0 * to_xxyz_xzz[k] * tbe_0 + 4.0 * to_xxyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xxz_xyyy[k] = -2.0 * to_xxyz_yyy[k] * tbe_0 + 4.0 * to_xxyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xxz_xyyz[k] = -2.0 * to_xxyz_yyz[k] * tbe_0 + 4.0 * to_xxyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xxz_xyzz[k] = -2.0 * to_xxyz_yzz[k] * tbe_0 + 4.0 * to_xxyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xxz_xzzz[k] = -2.0 * to_xxyz_zzz[k] * tbe_0 + 4.0 * to_xxyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xxz_yyyy[k] = 4.0 * to_xxyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xxz_yyyz[k] = 4.0 * to_xxyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xxz_yyzz[k] = 4.0 * to_xxyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xxz_yzzz[k] = 4.0 * to_xxyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xxz_zzzz[k] = 4.0 * to_xxyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 495-510 components of targeted buffer : FG

        auto to_y_x_xyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 45);

        auto to_y_x_xyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 46);

        auto to_y_x_xyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 47);

        auto to_y_x_xyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 48);

        auto to_y_x_xyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 49);

        auto to_y_x_xyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 50);

        auto to_y_x_xyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 51);

        auto to_y_x_xyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 52);

        auto to_y_x_xyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 53);

        auto to_y_x_xyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 54);

        auto to_y_x_xyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 55);

        auto to_y_x_xyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 56);

        auto to_y_x_xyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 57);

        auto to_y_x_xyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 58);

        auto to_y_x_xyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xy_xxx, to_xy_xxxxx, to_xy_xxxxy, to_xy_xxxxz, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_zzz, to_xyyy_xxx, to_xyyy_xxxxx, to_xyyy_xxxxy, to_xyyy_xxxxz, to_xyyy_xxxyy, to_xyyy_xxxyz, to_xyyy_xxxzz, to_xyyy_xxy, to_xyyy_xxyyy, to_xyyy_xxyyz, to_xyyy_xxyzz, to_xyyy_xxz, to_xyyy_xxzzz, to_xyyy_xyy, to_xyyy_xyyyy, to_xyyy_xyyyz, to_xyyy_xyyzz, to_xyyy_xyz, to_xyyy_xyzzz, to_xyyy_xzz, to_xyyy_xzzzz, to_xyyy_yyy, to_xyyy_yyz, to_xyyy_yzz, to_xyyy_zzz, to_y_x_xyy_xxxx, to_y_x_xyy_xxxy, to_y_x_xyy_xxxz, to_y_x_xyy_xxyy, to_y_x_xyy_xxyz, to_y_x_xyy_xxzz, to_y_x_xyy_xyyy, to_y_x_xyy_xyyz, to_y_x_xyy_xyzz, to_y_x_xyy_xzzz, to_y_x_xyy_yyyy, to_y_x_xyy_yyyz, to_y_x_xyy_yyzz, to_y_x_xyy_yzzz, to_y_x_xyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyy_xxxx[k] = 8.0 * to_xy_xxx[k] - 4.0 * to_xy_xxxxx[k] * tke_0 - 8.0 * to_xyyy_xxx[k] * tbe_0 + 4.0 * to_xyyy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xyy_xxxy[k] = 6.0 * to_xy_xxy[k] - 4.0 * to_xy_xxxxy[k] * tke_0 - 6.0 * to_xyyy_xxy[k] * tbe_0 + 4.0 * to_xyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xyy_xxxz[k] = 6.0 * to_xy_xxz[k] - 4.0 * to_xy_xxxxz[k] * tke_0 - 6.0 * to_xyyy_xxz[k] * tbe_0 + 4.0 * to_xyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xyy_xxyy[k] = 4.0 * to_xy_xyy[k] - 4.0 * to_xy_xxxyy[k] * tke_0 - 4.0 * to_xyyy_xyy[k] * tbe_0 + 4.0 * to_xyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xyy_xxyz[k] = 4.0 * to_xy_xyz[k] - 4.0 * to_xy_xxxyz[k] * tke_0 - 4.0 * to_xyyy_xyz[k] * tbe_0 + 4.0 * to_xyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xyy_xxzz[k] = 4.0 * to_xy_xzz[k] - 4.0 * to_xy_xxxzz[k] * tke_0 - 4.0 * to_xyyy_xzz[k] * tbe_0 + 4.0 * to_xyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xyy_xyyy[k] = 2.0 * to_xy_yyy[k] - 4.0 * to_xy_xxyyy[k] * tke_0 - 2.0 * to_xyyy_yyy[k] * tbe_0 + 4.0 * to_xyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xyy_xyyz[k] = 2.0 * to_xy_yyz[k] - 4.0 * to_xy_xxyyz[k] * tke_0 - 2.0 * to_xyyy_yyz[k] * tbe_0 + 4.0 * to_xyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xyy_xyzz[k] = 2.0 * to_xy_yzz[k] - 4.0 * to_xy_xxyzz[k] * tke_0 - 2.0 * to_xyyy_yzz[k] * tbe_0 + 4.0 * to_xyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xyy_xzzz[k] = 2.0 * to_xy_zzz[k] - 4.0 * to_xy_xxzzz[k] * tke_0 - 2.0 * to_xyyy_zzz[k] * tbe_0 + 4.0 * to_xyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xyy_yyyy[k] = -4.0 * to_xy_xyyyy[k] * tke_0 + 4.0 * to_xyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xyy_yyyz[k] = -4.0 * to_xy_xyyyz[k] * tke_0 + 4.0 * to_xyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xyy_yyzz[k] = -4.0 * to_xy_xyyzz[k] * tke_0 + 4.0 * to_xyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xyy_yzzz[k] = -4.0 * to_xy_xyzzz[k] * tke_0 + 4.0 * to_xyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xyy_zzzz[k] = -4.0 * to_xy_xzzzz[k] * tke_0 + 4.0 * to_xyyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 510-525 components of targeted buffer : FG

        auto to_y_x_xyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 60);

        auto to_y_x_xyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 61);

        auto to_y_x_xyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 62);

        auto to_y_x_xyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 63);

        auto to_y_x_xyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 64);

        auto to_y_x_xyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 65);

        auto to_y_x_xyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 66);

        auto to_y_x_xyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 67);

        auto to_y_x_xyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 68);

        auto to_y_x_xyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 69);

        auto to_y_x_xyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 70);

        auto to_y_x_xyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 71);

        auto to_y_x_xyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 72);

        auto to_y_x_xyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 73);

        auto to_y_x_xyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_xyyz_xxx, to_xyyz_xxxxx, to_xyyz_xxxxy, to_xyyz_xxxxz, to_xyyz_xxxyy, to_xyyz_xxxyz, to_xyyz_xxxzz, to_xyyz_xxy, to_xyyz_xxyyy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xxzzz, to_xyyz_xyy, to_xyyz_xyyyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_xzzzz, to_xyyz_yyy, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_zzz, to_xz_xxx, to_xz_xxxxx, to_xz_xxxxy, to_xz_xxxxz, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_zzz, to_y_x_xyz_xxxx, to_y_x_xyz_xxxy, to_y_x_xyz_xxxz, to_y_x_xyz_xxyy, to_y_x_xyz_xxyz, to_y_x_xyz_xxzz, to_y_x_xyz_xyyy, to_y_x_xyz_xyyz, to_y_x_xyz_xyzz, to_y_x_xyz_xzzz, to_y_x_xyz_yyyy, to_y_x_xyz_yyyz, to_y_x_xyz_yyzz, to_y_x_xyz_yzzz, to_y_x_xyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyz_xxxx[k] = 4.0 * to_xz_xxx[k] - 2.0 * to_xz_xxxxx[k] * tke_0 - 8.0 * to_xyyz_xxx[k] * tbe_0 + 4.0 * to_xyyz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xyz_xxxy[k] = 3.0 * to_xz_xxy[k] - 2.0 * to_xz_xxxxy[k] * tke_0 - 6.0 * to_xyyz_xxy[k] * tbe_0 + 4.0 * to_xyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xyz_xxxz[k] = 3.0 * to_xz_xxz[k] - 2.0 * to_xz_xxxxz[k] * tke_0 - 6.0 * to_xyyz_xxz[k] * tbe_0 + 4.0 * to_xyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xyz_xxyy[k] = 2.0 * to_xz_xyy[k] - 2.0 * to_xz_xxxyy[k] * tke_0 - 4.0 * to_xyyz_xyy[k] * tbe_0 + 4.0 * to_xyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xyz_xxyz[k] = 2.0 * to_xz_xyz[k] - 2.0 * to_xz_xxxyz[k] * tke_0 - 4.0 * to_xyyz_xyz[k] * tbe_0 + 4.0 * to_xyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xyz_xxzz[k] = 2.0 * to_xz_xzz[k] - 2.0 * to_xz_xxxzz[k] * tke_0 - 4.0 * to_xyyz_xzz[k] * tbe_0 + 4.0 * to_xyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xyz_xyyy[k] = to_xz_yyy[k] - 2.0 * to_xz_xxyyy[k] * tke_0 - 2.0 * to_xyyz_yyy[k] * tbe_0 + 4.0 * to_xyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xyz_xyyz[k] = to_xz_yyz[k] - 2.0 * to_xz_xxyyz[k] * tke_0 - 2.0 * to_xyyz_yyz[k] * tbe_0 + 4.0 * to_xyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xyz_xyzz[k] = to_xz_yzz[k] - 2.0 * to_xz_xxyzz[k] * tke_0 - 2.0 * to_xyyz_yzz[k] * tbe_0 + 4.0 * to_xyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xyz_xzzz[k] = to_xz_zzz[k] - 2.0 * to_xz_xxzzz[k] * tke_0 - 2.0 * to_xyyz_zzz[k] * tbe_0 + 4.0 * to_xyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xyz_yyyy[k] = -2.0 * to_xz_xyyyy[k] * tke_0 + 4.0 * to_xyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xyz_yyyz[k] = -2.0 * to_xz_xyyyz[k] * tke_0 + 4.0 * to_xyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xyz_yyzz[k] = -2.0 * to_xz_xyyzz[k] * tke_0 + 4.0 * to_xyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xyz_yzzz[k] = -2.0 * to_xz_xyzzz[k] * tke_0 + 4.0 * to_xyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xyz_zzzz[k] = -2.0 * to_xz_xzzzz[k] * tke_0 + 4.0 * to_xyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 525-540 components of targeted buffer : FG

        auto to_y_x_xzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 75);

        auto to_y_x_xzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 76);

        auto to_y_x_xzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 77);

        auto to_y_x_xzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 78);

        auto to_y_x_xzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 79);

        auto to_y_x_xzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 80);

        auto to_y_x_xzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 81);

        auto to_y_x_xzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 82);

        auto to_y_x_xzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 83);

        auto to_y_x_xzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 84);

        auto to_y_x_xzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 85);

        auto to_y_x_xzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 86);

        auto to_y_x_xzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 87);

        auto to_y_x_xzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 88);

        auto to_y_x_xzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xyzz_xxx, to_xyzz_xxxxx, to_xyzz_xxxxy, to_xyzz_xxxxz, to_xyzz_xxxyy, to_xyzz_xxxyz, to_xyzz_xxxzz, to_xyzz_xxy, to_xyzz_xxyyy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xxzzz, to_xyzz_xyy, to_xyzz_xyyyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_xzzzz, to_xyzz_yyy, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_zzz, to_y_x_xzz_xxxx, to_y_x_xzz_xxxy, to_y_x_xzz_xxxz, to_y_x_xzz_xxyy, to_y_x_xzz_xxyz, to_y_x_xzz_xxzz, to_y_x_xzz_xyyy, to_y_x_xzz_xyyz, to_y_x_xzz_xyzz, to_y_x_xzz_xzzz, to_y_x_xzz_yyyy, to_y_x_xzz_yyyz, to_y_x_xzz_yyzz, to_y_x_xzz_yzzz, to_y_x_xzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xzz_xxxx[k] = -8.0 * to_xyzz_xxx[k] * tbe_0 + 4.0 * to_xyzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_xzz_xxxy[k] = -6.0 * to_xyzz_xxy[k] * tbe_0 + 4.0 * to_xyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_xzz_xxxz[k] = -6.0 * to_xyzz_xxz[k] * tbe_0 + 4.0 * to_xyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_xzz_xxyy[k] = -4.0 * to_xyzz_xyy[k] * tbe_0 + 4.0 * to_xyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_xzz_xxyz[k] = -4.0 * to_xyzz_xyz[k] * tbe_0 + 4.0 * to_xyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_xzz_xxzz[k] = -4.0 * to_xyzz_xzz[k] * tbe_0 + 4.0 * to_xyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_xzz_xyyy[k] = -2.0 * to_xyzz_yyy[k] * tbe_0 + 4.0 * to_xyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_xzz_xyyz[k] = -2.0 * to_xyzz_yyz[k] * tbe_0 + 4.0 * to_xyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_xzz_xyzz[k] = -2.0 * to_xyzz_yzz[k] * tbe_0 + 4.0 * to_xyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_xzz_xzzz[k] = -2.0 * to_xyzz_zzz[k] * tbe_0 + 4.0 * to_xyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_xzz_yyyy[k] = 4.0 * to_xyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_xzz_yyyz[k] = 4.0 * to_xyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_xzz_yyzz[k] = 4.0 * to_xyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_xzz_yzzz[k] = 4.0 * to_xyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_xzz_zzzz[k] = 4.0 * to_xyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 540-555 components of targeted buffer : FG

        auto to_y_x_yyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 90);

        auto to_y_x_yyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 91);

        auto to_y_x_yyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 92);

        auto to_y_x_yyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 93);

        auto to_y_x_yyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 94);

        auto to_y_x_yyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 95);

        auto to_y_x_yyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 96);

        auto to_y_x_yyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 97);

        auto to_y_x_yyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 98);

        auto to_y_x_yyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 99);

        auto to_y_x_yyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 100);

        auto to_y_x_yyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 101);

        auto to_y_x_yyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 102);

        auto to_y_x_yyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 103);

        auto to_y_x_yyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_y_x_yyy_xxxx, to_y_x_yyy_xxxy, to_y_x_yyy_xxxz, to_y_x_yyy_xxyy, to_y_x_yyy_xxyz, to_y_x_yyy_xxzz, to_y_x_yyy_xyyy, to_y_x_yyy_xyyz, to_y_x_yyy_xyzz, to_y_x_yyy_xzzz, to_y_x_yyy_yyyy, to_y_x_yyy_yyyz, to_y_x_yyy_yyzz, to_y_x_yyy_yzzz, to_y_x_yyy_zzzz, to_yy_xxx, to_yy_xxxxx, to_yy_xxxxy, to_yy_xxxxz, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyz, to_yy_yzz, to_yy_zzz, to_yyyy_xxx, to_yyyy_xxxxx, to_yyyy_xxxxy, to_yyyy_xxxxz, to_yyyy_xxxyy, to_yyyy_xxxyz, to_yyyy_xxxzz, to_yyyy_xxy, to_yyyy_xxyyy, to_yyyy_xxyyz, to_yyyy_xxyzz, to_yyyy_xxz, to_yyyy_xxzzz, to_yyyy_xyy, to_yyyy_xyyyy, to_yyyy_xyyyz, to_yyyy_xyyzz, to_yyyy_xyz, to_yyyy_xyzzz, to_yyyy_xzz, to_yyyy_xzzzz, to_yyyy_yyy, to_yyyy_yyz, to_yyyy_yzz, to_yyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyy_xxxx[k] = 12.0 * to_yy_xxx[k] - 6.0 * to_yy_xxxxx[k] * tke_0 - 8.0 * to_yyyy_xxx[k] * tbe_0 + 4.0 * to_yyyy_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_yyy_xxxy[k] = 9.0 * to_yy_xxy[k] - 6.0 * to_yy_xxxxy[k] * tke_0 - 6.0 * to_yyyy_xxy[k] * tbe_0 + 4.0 * to_yyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_yyy_xxxz[k] = 9.0 * to_yy_xxz[k] - 6.0 * to_yy_xxxxz[k] * tke_0 - 6.0 * to_yyyy_xxz[k] * tbe_0 + 4.0 * to_yyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_yyy_xxyy[k] = 6.0 * to_yy_xyy[k] - 6.0 * to_yy_xxxyy[k] * tke_0 - 4.0 * to_yyyy_xyy[k] * tbe_0 + 4.0 * to_yyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_yyy_xxyz[k] = 6.0 * to_yy_xyz[k] - 6.0 * to_yy_xxxyz[k] * tke_0 - 4.0 * to_yyyy_xyz[k] * tbe_0 + 4.0 * to_yyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_yyy_xxzz[k] = 6.0 * to_yy_xzz[k] - 6.0 * to_yy_xxxzz[k] * tke_0 - 4.0 * to_yyyy_xzz[k] * tbe_0 + 4.0 * to_yyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_yyy_xyyy[k] = 3.0 * to_yy_yyy[k] - 6.0 * to_yy_xxyyy[k] * tke_0 - 2.0 * to_yyyy_yyy[k] * tbe_0 + 4.0 * to_yyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_yyy_xyyz[k] = 3.0 * to_yy_yyz[k] - 6.0 * to_yy_xxyyz[k] * tke_0 - 2.0 * to_yyyy_yyz[k] * tbe_0 + 4.0 * to_yyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_yyy_xyzz[k] = 3.0 * to_yy_yzz[k] - 6.0 * to_yy_xxyzz[k] * tke_0 - 2.0 * to_yyyy_yzz[k] * tbe_0 + 4.0 * to_yyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_yyy_xzzz[k] = 3.0 * to_yy_zzz[k] - 6.0 * to_yy_xxzzz[k] * tke_0 - 2.0 * to_yyyy_zzz[k] * tbe_0 + 4.0 * to_yyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_yyy_yyyy[k] = -6.0 * to_yy_xyyyy[k] * tke_0 + 4.0 * to_yyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_yyy_yyyz[k] = -6.0 * to_yy_xyyyz[k] * tke_0 + 4.0 * to_yyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_yyy_yyzz[k] = -6.0 * to_yy_xyyzz[k] * tke_0 + 4.0 * to_yyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_yyy_yzzz[k] = -6.0 * to_yy_xyzzz[k] * tke_0 + 4.0 * to_yyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_yyy_zzzz[k] = -6.0 * to_yy_xzzzz[k] * tke_0 + 4.0 * to_yyyy_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 555-570 components of targeted buffer : FG

        auto to_y_x_yyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 105);

        auto to_y_x_yyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 106);

        auto to_y_x_yyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 107);

        auto to_y_x_yyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 108);

        auto to_y_x_yyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 109);

        auto to_y_x_yyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 110);

        auto to_y_x_yyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 111);

        auto to_y_x_yyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 112);

        auto to_y_x_yyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 113);

        auto to_y_x_yyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 114);

        auto to_y_x_yyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 115);

        auto to_y_x_yyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 116);

        auto to_y_x_yyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 117);

        auto to_y_x_yyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 118);

        auto to_y_x_yyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_y_x_yyz_xxxx, to_y_x_yyz_xxxy, to_y_x_yyz_xxxz, to_y_x_yyz_xxyy, to_y_x_yyz_xxyz, to_y_x_yyz_xxzz, to_y_x_yyz_xyyy, to_y_x_yyz_xyyz, to_y_x_yyz_xyzz, to_y_x_yyz_xzzz, to_y_x_yyz_yyyy, to_y_x_yyz_yyyz, to_y_x_yyz_yyzz, to_y_x_yyz_yzzz, to_y_x_yyz_zzzz, to_yyyz_xxx, to_yyyz_xxxxx, to_yyyz_xxxxy, to_yyyz_xxxxz, to_yyyz_xxxyy, to_yyyz_xxxyz, to_yyyz_xxxzz, to_yyyz_xxy, to_yyyz_xxyyy, to_yyyz_xxyyz, to_yyyz_xxyzz, to_yyyz_xxz, to_yyyz_xxzzz, to_yyyz_xyy, to_yyyz_xyyyy, to_yyyz_xyyyz, to_yyyz_xyyzz, to_yyyz_xyz, to_yyyz_xyzzz, to_yyyz_xzz, to_yyyz_xzzzz, to_yyyz_yyy, to_yyyz_yyz, to_yyyz_yzz, to_yyyz_zzz, to_yz_xxx, to_yz_xxxxx, to_yz_xxxxy, to_yz_xxxxz, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyz_xxxx[k] = 8.0 * to_yz_xxx[k] - 4.0 * to_yz_xxxxx[k] * tke_0 - 8.0 * to_yyyz_xxx[k] * tbe_0 + 4.0 * to_yyyz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_yyz_xxxy[k] = 6.0 * to_yz_xxy[k] - 4.0 * to_yz_xxxxy[k] * tke_0 - 6.0 * to_yyyz_xxy[k] * tbe_0 + 4.0 * to_yyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_yyz_xxxz[k] = 6.0 * to_yz_xxz[k] - 4.0 * to_yz_xxxxz[k] * tke_0 - 6.0 * to_yyyz_xxz[k] * tbe_0 + 4.0 * to_yyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_yyz_xxyy[k] = 4.0 * to_yz_xyy[k] - 4.0 * to_yz_xxxyy[k] * tke_0 - 4.0 * to_yyyz_xyy[k] * tbe_0 + 4.0 * to_yyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_yyz_xxyz[k] = 4.0 * to_yz_xyz[k] - 4.0 * to_yz_xxxyz[k] * tke_0 - 4.0 * to_yyyz_xyz[k] * tbe_0 + 4.0 * to_yyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_yyz_xxzz[k] = 4.0 * to_yz_xzz[k] - 4.0 * to_yz_xxxzz[k] * tke_0 - 4.0 * to_yyyz_xzz[k] * tbe_0 + 4.0 * to_yyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_yyz_xyyy[k] = 2.0 * to_yz_yyy[k] - 4.0 * to_yz_xxyyy[k] * tke_0 - 2.0 * to_yyyz_yyy[k] * tbe_0 + 4.0 * to_yyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_yyz_xyyz[k] = 2.0 * to_yz_yyz[k] - 4.0 * to_yz_xxyyz[k] * tke_0 - 2.0 * to_yyyz_yyz[k] * tbe_0 + 4.0 * to_yyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_yyz_xyzz[k] = 2.0 * to_yz_yzz[k] - 4.0 * to_yz_xxyzz[k] * tke_0 - 2.0 * to_yyyz_yzz[k] * tbe_0 + 4.0 * to_yyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_yyz_xzzz[k] = 2.0 * to_yz_zzz[k] - 4.0 * to_yz_xxzzz[k] * tke_0 - 2.0 * to_yyyz_zzz[k] * tbe_0 + 4.0 * to_yyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_yyz_yyyy[k] = -4.0 * to_yz_xyyyy[k] * tke_0 + 4.0 * to_yyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_yyz_yyyz[k] = -4.0 * to_yz_xyyyz[k] * tke_0 + 4.0 * to_yyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_yyz_yyzz[k] = -4.0 * to_yz_xyyzz[k] * tke_0 + 4.0 * to_yyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_yyz_yzzz[k] = -4.0 * to_yz_xyzzz[k] * tke_0 + 4.0 * to_yyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_yyz_zzzz[k] = -4.0 * to_yz_xzzzz[k] * tke_0 + 4.0 * to_yyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 570-585 components of targeted buffer : FG

        auto to_y_x_yzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 120);

        auto to_y_x_yzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 121);

        auto to_y_x_yzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 122);

        auto to_y_x_yzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 123);

        auto to_y_x_yzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 124);

        auto to_y_x_yzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 125);

        auto to_y_x_yzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 126);

        auto to_y_x_yzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 127);

        auto to_y_x_yzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 128);

        auto to_y_x_yzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 129);

        auto to_y_x_yzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 130);

        auto to_y_x_yzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 131);

        auto to_y_x_yzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 132);

        auto to_y_x_yzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 133);

        auto to_y_x_yzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_y_x_yzz_xxxx, to_y_x_yzz_xxxy, to_y_x_yzz_xxxz, to_y_x_yzz_xxyy, to_y_x_yzz_xxyz, to_y_x_yzz_xxzz, to_y_x_yzz_xyyy, to_y_x_yzz_xyyz, to_y_x_yzz_xyzz, to_y_x_yzz_xzzz, to_y_x_yzz_yyyy, to_y_x_yzz_yyyz, to_y_x_yzz_yyzz, to_y_x_yzz_yzzz, to_y_x_yzz_zzzz, to_yyzz_xxx, to_yyzz_xxxxx, to_yyzz_xxxxy, to_yyzz_xxxxz, to_yyzz_xxxyy, to_yyzz_xxxyz, to_yyzz_xxxzz, to_yyzz_xxy, to_yyzz_xxyyy, to_yyzz_xxyyz, to_yyzz_xxyzz, to_yyzz_xxz, to_yyzz_xxzzz, to_yyzz_xyy, to_yyzz_xyyyy, to_yyzz_xyyyz, to_yyzz_xyyzz, to_yyzz_xyz, to_yyzz_xyzzz, to_yyzz_xzz, to_yyzz_xzzzz, to_yyzz_yyy, to_yyzz_yyz, to_yyzz_yzz, to_yyzz_zzz, to_zz_xxx, to_zz_xxxxx, to_zz_xxxxy, to_zz_xxxxz, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyz, to_zz_yzz, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yzz_xxxx[k] = 4.0 * to_zz_xxx[k] - 2.0 * to_zz_xxxxx[k] * tke_0 - 8.0 * to_yyzz_xxx[k] * tbe_0 + 4.0 * to_yyzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_yzz_xxxy[k] = 3.0 * to_zz_xxy[k] - 2.0 * to_zz_xxxxy[k] * tke_0 - 6.0 * to_yyzz_xxy[k] * tbe_0 + 4.0 * to_yyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_yzz_xxxz[k] = 3.0 * to_zz_xxz[k] - 2.0 * to_zz_xxxxz[k] * tke_0 - 6.0 * to_yyzz_xxz[k] * tbe_0 + 4.0 * to_yyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_yzz_xxyy[k] = 2.0 * to_zz_xyy[k] - 2.0 * to_zz_xxxyy[k] * tke_0 - 4.0 * to_yyzz_xyy[k] * tbe_0 + 4.0 * to_yyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_yzz_xxyz[k] = 2.0 * to_zz_xyz[k] - 2.0 * to_zz_xxxyz[k] * tke_0 - 4.0 * to_yyzz_xyz[k] * tbe_0 + 4.0 * to_yyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_yzz_xxzz[k] = 2.0 * to_zz_xzz[k] - 2.0 * to_zz_xxxzz[k] * tke_0 - 4.0 * to_yyzz_xzz[k] * tbe_0 + 4.0 * to_yyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_yzz_xyyy[k] = to_zz_yyy[k] - 2.0 * to_zz_xxyyy[k] * tke_0 - 2.0 * to_yyzz_yyy[k] * tbe_0 + 4.0 * to_yyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_yzz_xyyz[k] = to_zz_yyz[k] - 2.0 * to_zz_xxyyz[k] * tke_0 - 2.0 * to_yyzz_yyz[k] * tbe_0 + 4.0 * to_yyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_yzz_xyzz[k] = to_zz_yzz[k] - 2.0 * to_zz_xxyzz[k] * tke_0 - 2.0 * to_yyzz_yzz[k] * tbe_0 + 4.0 * to_yyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_yzz_xzzz[k] = to_zz_zzz[k] - 2.0 * to_zz_xxzzz[k] * tke_0 - 2.0 * to_yyzz_zzz[k] * tbe_0 + 4.0 * to_yyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_yzz_yyyy[k] = -2.0 * to_zz_xyyyy[k] * tke_0 + 4.0 * to_yyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_yzz_yyyz[k] = -2.0 * to_zz_xyyyz[k] * tke_0 + 4.0 * to_yyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_yzz_yyzz[k] = -2.0 * to_zz_xyyzz[k] * tke_0 + 4.0 * to_yyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_yzz_yzzz[k] = -2.0 * to_zz_xyzzz[k] * tke_0 + 4.0 * to_yyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_yzz_zzzz[k] = -2.0 * to_zz_xzzzz[k] * tke_0 + 4.0 * to_yyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 585-600 components of targeted buffer : FG

        auto to_y_x_zzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 135);

        auto to_y_x_zzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 136);

        auto to_y_x_zzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 137);

        auto to_y_x_zzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 138);

        auto to_y_x_zzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 139);

        auto to_y_x_zzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 140);

        auto to_y_x_zzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 141);

        auto to_y_x_zzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 142);

        auto to_y_x_zzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 143);

        auto to_y_x_zzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 144);

        auto to_y_x_zzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 145);

        auto to_y_x_zzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 146);

        auto to_y_x_zzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 147);

        auto to_y_x_zzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 148);

        auto to_y_x_zzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 3 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_y_x_zzz_xxxx, to_y_x_zzz_xxxy, to_y_x_zzz_xxxz, to_y_x_zzz_xxyy, to_y_x_zzz_xxyz, to_y_x_zzz_xxzz, to_y_x_zzz_xyyy, to_y_x_zzz_xyyz, to_y_x_zzz_xyzz, to_y_x_zzz_xzzz, to_y_x_zzz_yyyy, to_y_x_zzz_yyyz, to_y_x_zzz_yyzz, to_y_x_zzz_yzzz, to_y_x_zzz_zzzz, to_yzzz_xxx, to_yzzz_xxxxx, to_yzzz_xxxxy, to_yzzz_xxxxz, to_yzzz_xxxyy, to_yzzz_xxxyz, to_yzzz_xxxzz, to_yzzz_xxy, to_yzzz_xxyyy, to_yzzz_xxyyz, to_yzzz_xxyzz, to_yzzz_xxz, to_yzzz_xxzzz, to_yzzz_xyy, to_yzzz_xyyyy, to_yzzz_xyyyz, to_yzzz_xyyzz, to_yzzz_xyz, to_yzzz_xyzzz, to_yzzz_xzz, to_yzzz_xzzzz, to_yzzz_yyy, to_yzzz_yyz, to_yzzz_yzz, to_yzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zzz_xxxx[k] = -8.0 * to_yzzz_xxx[k] * tbe_0 + 4.0 * to_yzzz_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_zzz_xxxy[k] = -6.0 * to_yzzz_xxy[k] * tbe_0 + 4.0 * to_yzzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_zzz_xxxz[k] = -6.0 * to_yzzz_xxz[k] * tbe_0 + 4.0 * to_yzzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_zzz_xxyy[k] = -4.0 * to_yzzz_xyy[k] * tbe_0 + 4.0 * to_yzzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_zzz_xxyz[k] = -4.0 * to_yzzz_xyz[k] * tbe_0 + 4.0 * to_yzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_zzz_xxzz[k] = -4.0 * to_yzzz_xzz[k] * tbe_0 + 4.0 * to_yzzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_zzz_xyyy[k] = -2.0 * to_yzzz_yyy[k] * tbe_0 + 4.0 * to_yzzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_zzz_xyyz[k] = -2.0 * to_yzzz_yyz[k] * tbe_0 + 4.0 * to_yzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_zzz_xyzz[k] = -2.0 * to_yzzz_yzz[k] * tbe_0 + 4.0 * to_yzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_zzz_xzzz[k] = -2.0 * to_yzzz_zzz[k] * tbe_0 + 4.0 * to_yzzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_zzz_yyyy[k] = 4.0 * to_yzzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_zzz_yyyz[k] = 4.0 * to_yzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_zzz_yyzz[k] = 4.0 * to_yzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_zzz_yzzz[k] = 4.0 * to_yzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_zzz_zzzz[k] = 4.0 * to_yzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 600-615 components of targeted buffer : FG

        auto to_y_y_xxx_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 0);

        auto to_y_y_xxx_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 1);

        auto to_y_y_xxx_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 2);

        auto to_y_y_xxx_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 3);

        auto to_y_y_xxx_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 4);

        auto to_y_y_xxx_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 5);

        auto to_y_y_xxx_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 6);

        auto to_y_y_xxx_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 7);

        auto to_y_y_xxx_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 8);

        auto to_y_y_xxx_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 9);

        auto to_y_y_xxx_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 10);

        auto to_y_y_xxx_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 11);

        auto to_y_y_xxx_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 12);

        auto to_y_y_xxx_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 13);

        auto to_y_y_xxx_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_xxxy_xxx, to_xxxy_xxxxy, to_xxxy_xxxyy, to_xxxy_xxxyz, to_xxxy_xxy, to_xxxy_xxyyy, to_xxxy_xxyyz, to_xxxy_xxyzz, to_xxxy_xxz, to_xxxy_xyy, to_xxxy_xyyyy, to_xxxy_xyyyz, to_xxxy_xyyzz, to_xxxy_xyz, to_xxxy_xyzzz, to_xxxy_xzz, to_xxxy_yyy, to_xxxy_yyyyy, to_xxxy_yyyyz, to_xxxy_yyyzz, to_xxxy_yyz, to_xxxy_yyzzz, to_xxxy_yzz, to_xxxy_yzzzz, to_xxxy_zzz, to_y_y_xxx_xxxx, to_y_y_xxx_xxxy, to_y_y_xxx_xxxz, to_y_y_xxx_xxyy, to_y_y_xxx_xxyz, to_y_y_xxx_xxzz, to_y_y_xxx_xyyy, to_y_y_xxx_xyyz, to_y_y_xxx_xyzz, to_y_y_xxx_xzzz, to_y_y_xxx_yyyy, to_y_y_xxx_yyyz, to_y_y_xxx_yyzz, to_y_y_xxx_yzzz, to_y_y_xxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxx_xxxx[k] = 4.0 * to_xxxy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xxx_xxxy[k] = -2.0 * to_xxxy_xxx[k] * tbe_0 + 4.0 * to_xxxy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xxx_xxxz[k] = 4.0 * to_xxxy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xxx_xxyy[k] = -4.0 * to_xxxy_xxy[k] * tbe_0 + 4.0 * to_xxxy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xxx_xxyz[k] = -2.0 * to_xxxy_xxz[k] * tbe_0 + 4.0 * to_xxxy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xxx_xxzz[k] = 4.0 * to_xxxy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xxx_xyyy[k] = -6.0 * to_xxxy_xyy[k] * tbe_0 + 4.0 * to_xxxy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxx_xyyz[k] = -4.0 * to_xxxy_xyz[k] * tbe_0 + 4.0 * to_xxxy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxx_xyzz[k] = -2.0 * to_xxxy_xzz[k] * tbe_0 + 4.0 * to_xxxy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxx_xzzz[k] = 4.0 * to_xxxy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxx_yyyy[k] = -8.0 * to_xxxy_yyy[k] * tbe_0 + 4.0 * to_xxxy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxx_yyyz[k] = -6.0 * to_xxxy_yyz[k] * tbe_0 + 4.0 * to_xxxy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxx_yyzz[k] = -4.0 * to_xxxy_yzz[k] * tbe_0 + 4.0 * to_xxxy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxx_yzzz[k] = -2.0 * to_xxxy_zzz[k] * tbe_0 + 4.0 * to_xxxy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxx_zzzz[k] = 4.0 * to_xxxy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 615-630 components of targeted buffer : FG

        auto to_y_y_xxy_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 15);

        auto to_y_y_xxy_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 16);

        auto to_y_y_xxy_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 17);

        auto to_y_y_xxy_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 18);

        auto to_y_y_xxy_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 19);

        auto to_y_y_xxy_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 20);

        auto to_y_y_xxy_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 21);

        auto to_y_y_xxy_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 22);

        auto to_y_y_xxy_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 23);

        auto to_y_y_xxy_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 24);

        auto to_y_y_xxy_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 25);

        auto to_y_y_xxy_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 26);

        auto to_y_y_xxy_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 27);

        auto to_y_y_xxy_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 28);

        auto to_y_y_xxy_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xx_xxx, to_xx_xxxxy, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_yyy, to_xx_yyyyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, to_xxyy_xxx, to_xxyy_xxxxy, to_xxyy_xxxyy, to_xxyy_xxxyz, to_xxyy_xxy, to_xxyy_xxyyy, to_xxyy_xxyyz, to_xxyy_xxyzz, to_xxyy_xxz, to_xxyy_xyy, to_xxyy_xyyyy, to_xxyy_xyyyz, to_xxyy_xyyzz, to_xxyy_xyz, to_xxyy_xyzzz, to_xxyy_xzz, to_xxyy_yyy, to_xxyy_yyyyy, to_xxyy_yyyyz, to_xxyy_yyyzz, to_xxyy_yyz, to_xxyy_yyzzz, to_xxyy_yzz, to_xxyy_yzzzz, to_xxyy_zzz, to_y_y_xxy_xxxx, to_y_y_xxy_xxxy, to_y_y_xxy_xxxz, to_y_y_xxy_xxyy, to_y_y_xxy_xxyz, to_y_y_xxy_xxzz, to_y_y_xxy_xyyy, to_y_y_xxy_xyyz, to_y_y_xxy_xyzz, to_y_y_xxy_xzzz, to_y_y_xxy_yyyy, to_y_y_xxy_yyyz, to_y_y_xxy_yyzz, to_y_y_xxy_yzzz, to_y_y_xxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxy_xxxx[k] = -2.0 * to_xx_xxxxy[k] * tke_0 + 4.0 * to_xxyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xxy_xxxy[k] = to_xx_xxx[k] - 2.0 * to_xx_xxxyy[k] * tke_0 - 2.0 * to_xxyy_xxx[k] * tbe_0 + 4.0 * to_xxyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xxy_xxxz[k] = -2.0 * to_xx_xxxyz[k] * tke_0 + 4.0 * to_xxyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xxy_xxyy[k] = 2.0 * to_xx_xxy[k] - 2.0 * to_xx_xxyyy[k] * tke_0 - 4.0 * to_xxyy_xxy[k] * tbe_0 + 4.0 * to_xxyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xxy_xxyz[k] = to_xx_xxz[k] - 2.0 * to_xx_xxyyz[k] * tke_0 - 2.0 * to_xxyy_xxz[k] * tbe_0 + 4.0 * to_xxyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xxy_xxzz[k] = -2.0 * to_xx_xxyzz[k] * tke_0 + 4.0 * to_xxyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xxy_xyyy[k] = 3.0 * to_xx_xyy[k] - 2.0 * to_xx_xyyyy[k] * tke_0 - 6.0 * to_xxyy_xyy[k] * tbe_0 + 4.0 * to_xxyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxy_xyyz[k] = 2.0 * to_xx_xyz[k] - 2.0 * to_xx_xyyyz[k] * tke_0 - 4.0 * to_xxyy_xyz[k] * tbe_0 + 4.0 * to_xxyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxy_xyzz[k] = to_xx_xzz[k] - 2.0 * to_xx_xyyzz[k] * tke_0 - 2.0 * to_xxyy_xzz[k] * tbe_0 + 4.0 * to_xxyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxy_xzzz[k] = -2.0 * to_xx_xyzzz[k] * tke_0 + 4.0 * to_xxyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxy_yyyy[k] = 4.0 * to_xx_yyy[k] - 2.0 * to_xx_yyyyy[k] * tke_0 - 8.0 * to_xxyy_yyy[k] * tbe_0 + 4.0 * to_xxyy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxy_yyyz[k] = 3.0 * to_xx_yyz[k] - 2.0 * to_xx_yyyyz[k] * tke_0 - 6.0 * to_xxyy_yyz[k] * tbe_0 + 4.0 * to_xxyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxy_yyzz[k] = 2.0 * to_xx_yzz[k] - 2.0 * to_xx_yyyzz[k] * tke_0 - 4.0 * to_xxyy_yzz[k] * tbe_0 + 4.0 * to_xxyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxy_yzzz[k] = to_xx_zzz[k] - 2.0 * to_xx_yyzzz[k] * tke_0 - 2.0 * to_xxyy_zzz[k] * tbe_0 + 4.0 * to_xxyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxy_zzzz[k] = -2.0 * to_xx_yzzzz[k] * tke_0 + 4.0 * to_xxyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 630-645 components of targeted buffer : FG

        auto to_y_y_xxz_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 30);

        auto to_y_y_xxz_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 31);

        auto to_y_y_xxz_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 32);

        auto to_y_y_xxz_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 33);

        auto to_y_y_xxz_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 34);

        auto to_y_y_xxz_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 35);

        auto to_y_y_xxz_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 36);

        auto to_y_y_xxz_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 37);

        auto to_y_y_xxz_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 38);

        auto to_y_y_xxz_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 39);

        auto to_y_y_xxz_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 40);

        auto to_y_y_xxz_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 41);

        auto to_y_y_xxz_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 42);

        auto to_y_y_xxz_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 43);

        auto to_y_y_xxz_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_xxyz_xxx, to_xxyz_xxxxy, to_xxyz_xxxyy, to_xxyz_xxxyz, to_xxyz_xxy, to_xxyz_xxyyy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xyy, to_xxyz_xyyyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_yyy, to_xxyz_yyyyy, to_xxyz_yyyyz, to_xxyz_yyyzz, to_xxyz_yyz, to_xxyz_yyzzz, to_xxyz_yzz, to_xxyz_yzzzz, to_xxyz_zzz, to_y_y_xxz_xxxx, to_y_y_xxz_xxxy, to_y_y_xxz_xxxz, to_y_y_xxz_xxyy, to_y_y_xxz_xxyz, to_y_y_xxz_xxzz, to_y_y_xxz_xyyy, to_y_y_xxz_xyyz, to_y_y_xxz_xyzz, to_y_y_xxz_xzzz, to_y_y_xxz_yyyy, to_y_y_xxz_yyyz, to_y_y_xxz_yyzz, to_y_y_xxz_yzzz, to_y_y_xxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxz_xxxx[k] = 4.0 * to_xxyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xxz_xxxy[k] = -2.0 * to_xxyz_xxx[k] * tbe_0 + 4.0 * to_xxyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xxz_xxxz[k] = 4.0 * to_xxyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xxz_xxyy[k] = -4.0 * to_xxyz_xxy[k] * tbe_0 + 4.0 * to_xxyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xxz_xxyz[k] = -2.0 * to_xxyz_xxz[k] * tbe_0 + 4.0 * to_xxyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xxz_xxzz[k] = 4.0 * to_xxyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xxz_xyyy[k] = -6.0 * to_xxyz_xyy[k] * tbe_0 + 4.0 * to_xxyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxz_xyyz[k] = -4.0 * to_xxyz_xyz[k] * tbe_0 + 4.0 * to_xxyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxz_xyzz[k] = -2.0 * to_xxyz_xzz[k] * tbe_0 + 4.0 * to_xxyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxz_xzzz[k] = 4.0 * to_xxyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxz_yyyy[k] = -8.0 * to_xxyz_yyy[k] * tbe_0 + 4.0 * to_xxyz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xxz_yyyz[k] = -6.0 * to_xxyz_yyz[k] * tbe_0 + 4.0 * to_xxyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xxz_yyzz[k] = -4.0 * to_xxyz_yzz[k] * tbe_0 + 4.0 * to_xxyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xxz_yzzz[k] = -2.0 * to_xxyz_zzz[k] * tbe_0 + 4.0 * to_xxyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xxz_zzzz[k] = 4.0 * to_xxyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 645-660 components of targeted buffer : FG

        auto to_y_y_xyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 45);

        auto to_y_y_xyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 46);

        auto to_y_y_xyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 47);

        auto to_y_y_xyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 48);

        auto to_y_y_xyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 49);

        auto to_y_y_xyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 50);

        auto to_y_y_xyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 51);

        auto to_y_y_xyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 52);

        auto to_y_y_xyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 53);

        auto to_y_y_xyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 54);

        auto to_y_y_xyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 55);

        auto to_y_y_xyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 56);

        auto to_y_y_xyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 57);

        auto to_y_y_xyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 58);

        auto to_y_y_xyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xy_xxx, to_xy_xxxxy, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_yyy, to_xy_yyyyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, to_xyyy_xxx, to_xyyy_xxxxy, to_xyyy_xxxyy, to_xyyy_xxxyz, to_xyyy_xxy, to_xyyy_xxyyy, to_xyyy_xxyyz, to_xyyy_xxyzz, to_xyyy_xxz, to_xyyy_xyy, to_xyyy_xyyyy, to_xyyy_xyyyz, to_xyyy_xyyzz, to_xyyy_xyz, to_xyyy_xyzzz, to_xyyy_xzz, to_xyyy_yyy, to_xyyy_yyyyy, to_xyyy_yyyyz, to_xyyy_yyyzz, to_xyyy_yyz, to_xyyy_yyzzz, to_xyyy_yzz, to_xyyy_yzzzz, to_xyyy_zzz, to_y_y_xyy_xxxx, to_y_y_xyy_xxxy, to_y_y_xyy_xxxz, to_y_y_xyy_xxyy, to_y_y_xyy_xxyz, to_y_y_xyy_xxzz, to_y_y_xyy_xyyy, to_y_y_xyy_xyyz, to_y_y_xyy_xyzz, to_y_y_xyy_xzzz, to_y_y_xyy_yyyy, to_y_y_xyy_yyyz, to_y_y_xyy_yyzz, to_y_y_xyy_yzzz, to_y_y_xyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyy_xxxx[k] = -4.0 * to_xy_xxxxy[k] * tke_0 + 4.0 * to_xyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xyy_xxxy[k] = 2.0 * to_xy_xxx[k] - 4.0 * to_xy_xxxyy[k] * tke_0 - 2.0 * to_xyyy_xxx[k] * tbe_0 + 4.0 * to_xyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xyy_xxxz[k] = -4.0 * to_xy_xxxyz[k] * tke_0 + 4.0 * to_xyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xyy_xxyy[k] = 4.0 * to_xy_xxy[k] - 4.0 * to_xy_xxyyy[k] * tke_0 - 4.0 * to_xyyy_xxy[k] * tbe_0 + 4.0 * to_xyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xyy_xxyz[k] = 2.0 * to_xy_xxz[k] - 4.0 * to_xy_xxyyz[k] * tke_0 - 2.0 * to_xyyy_xxz[k] * tbe_0 + 4.0 * to_xyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xyy_xxzz[k] = -4.0 * to_xy_xxyzz[k] * tke_0 + 4.0 * to_xyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xyy_xyyy[k] = 6.0 * to_xy_xyy[k] - 4.0 * to_xy_xyyyy[k] * tke_0 - 6.0 * to_xyyy_xyy[k] * tbe_0 + 4.0 * to_xyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyy_xyyz[k] = 4.0 * to_xy_xyz[k] - 4.0 * to_xy_xyyyz[k] * tke_0 - 4.0 * to_xyyy_xyz[k] * tbe_0 + 4.0 * to_xyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyy_xyzz[k] = 2.0 * to_xy_xzz[k] - 4.0 * to_xy_xyyzz[k] * tke_0 - 2.0 * to_xyyy_xzz[k] * tbe_0 + 4.0 * to_xyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyy_xzzz[k] = -4.0 * to_xy_xyzzz[k] * tke_0 + 4.0 * to_xyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyy_yyyy[k] = 8.0 * to_xy_yyy[k] - 4.0 * to_xy_yyyyy[k] * tke_0 - 8.0 * to_xyyy_yyy[k] * tbe_0 + 4.0 * to_xyyy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyy_yyyz[k] = 6.0 * to_xy_yyz[k] - 4.0 * to_xy_yyyyz[k] * tke_0 - 6.0 * to_xyyy_yyz[k] * tbe_0 + 4.0 * to_xyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyy_yyzz[k] = 4.0 * to_xy_yzz[k] - 4.0 * to_xy_yyyzz[k] * tke_0 - 4.0 * to_xyyy_yzz[k] * tbe_0 + 4.0 * to_xyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyy_yzzz[k] = 2.0 * to_xy_zzz[k] - 4.0 * to_xy_yyzzz[k] * tke_0 - 2.0 * to_xyyy_zzz[k] * tbe_0 + 4.0 * to_xyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyy_zzzz[k] = -4.0 * to_xy_yzzzz[k] * tke_0 + 4.0 * to_xyyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 660-675 components of targeted buffer : FG

        auto to_y_y_xyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 60);

        auto to_y_y_xyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 61);

        auto to_y_y_xyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 62);

        auto to_y_y_xyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 63);

        auto to_y_y_xyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 64);

        auto to_y_y_xyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 65);

        auto to_y_y_xyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 66);

        auto to_y_y_xyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 67);

        auto to_y_y_xyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 68);

        auto to_y_y_xyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 69);

        auto to_y_y_xyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 70);

        auto to_y_y_xyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 71);

        auto to_y_y_xyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 72);

        auto to_y_y_xyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 73);

        auto to_y_y_xyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_xyyz_xxx, to_xyyz_xxxxy, to_xyyz_xxxyy, to_xyyz_xxxyz, to_xyyz_xxy, to_xyyz_xxyyy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xyy, to_xyyz_xyyyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_yyy, to_xyyz_yyyyy, to_xyyz_yyyyz, to_xyyz_yyyzz, to_xyyz_yyz, to_xyyz_yyzzz, to_xyyz_yzz, to_xyyz_yzzzz, to_xyyz_zzz, to_xz_xxx, to_xz_xxxxy, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_yyy, to_xz_yyyyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, to_y_y_xyz_xxxx, to_y_y_xyz_xxxy, to_y_y_xyz_xxxz, to_y_y_xyz_xxyy, to_y_y_xyz_xxyz, to_y_y_xyz_xxzz, to_y_y_xyz_xyyy, to_y_y_xyz_xyyz, to_y_y_xyz_xyzz, to_y_y_xyz_xzzz, to_y_y_xyz_yyyy, to_y_y_xyz_yyyz, to_y_y_xyz_yyzz, to_y_y_xyz_yzzz, to_y_y_xyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyz_xxxx[k] = -2.0 * to_xz_xxxxy[k] * tke_0 + 4.0 * to_xyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xyz_xxxy[k] = to_xz_xxx[k] - 2.0 * to_xz_xxxyy[k] * tke_0 - 2.0 * to_xyyz_xxx[k] * tbe_0 + 4.0 * to_xyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xyz_xxxz[k] = -2.0 * to_xz_xxxyz[k] * tke_0 + 4.0 * to_xyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xyz_xxyy[k] = 2.0 * to_xz_xxy[k] - 2.0 * to_xz_xxyyy[k] * tke_0 - 4.0 * to_xyyz_xxy[k] * tbe_0 + 4.0 * to_xyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xyz_xxyz[k] = to_xz_xxz[k] - 2.0 * to_xz_xxyyz[k] * tke_0 - 2.0 * to_xyyz_xxz[k] * tbe_0 + 4.0 * to_xyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xyz_xxzz[k] = -2.0 * to_xz_xxyzz[k] * tke_0 + 4.0 * to_xyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xyz_xyyy[k] = 3.0 * to_xz_xyy[k] - 2.0 * to_xz_xyyyy[k] * tke_0 - 6.0 * to_xyyz_xyy[k] * tbe_0 + 4.0 * to_xyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyz_xyyz[k] = 2.0 * to_xz_xyz[k] - 2.0 * to_xz_xyyyz[k] * tke_0 - 4.0 * to_xyyz_xyz[k] * tbe_0 + 4.0 * to_xyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyz_xyzz[k] = to_xz_xzz[k] - 2.0 * to_xz_xyyzz[k] * tke_0 - 2.0 * to_xyyz_xzz[k] * tbe_0 + 4.0 * to_xyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyz_xzzz[k] = -2.0 * to_xz_xyzzz[k] * tke_0 + 4.0 * to_xyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyz_yyyy[k] = 4.0 * to_xz_yyy[k] - 2.0 * to_xz_yyyyy[k] * tke_0 - 8.0 * to_xyyz_yyy[k] * tbe_0 + 4.0 * to_xyyz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xyz_yyyz[k] = 3.0 * to_xz_yyz[k] - 2.0 * to_xz_yyyyz[k] * tke_0 - 6.0 * to_xyyz_yyz[k] * tbe_0 + 4.0 * to_xyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xyz_yyzz[k] = 2.0 * to_xz_yzz[k] - 2.0 * to_xz_yyyzz[k] * tke_0 - 4.0 * to_xyyz_yzz[k] * tbe_0 + 4.0 * to_xyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xyz_yzzz[k] = to_xz_zzz[k] - 2.0 * to_xz_yyzzz[k] * tke_0 - 2.0 * to_xyyz_zzz[k] * tbe_0 + 4.0 * to_xyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xyz_zzzz[k] = -2.0 * to_xz_yzzzz[k] * tke_0 + 4.0 * to_xyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 675-690 components of targeted buffer : FG

        auto to_y_y_xzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 75);

        auto to_y_y_xzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 76);

        auto to_y_y_xzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 77);

        auto to_y_y_xzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 78);

        auto to_y_y_xzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 79);

        auto to_y_y_xzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 80);

        auto to_y_y_xzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 81);

        auto to_y_y_xzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 82);

        auto to_y_y_xzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 83);

        auto to_y_y_xzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 84);

        auto to_y_y_xzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 85);

        auto to_y_y_xzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 86);

        auto to_y_y_xzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 87);

        auto to_y_y_xzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 88);

        auto to_y_y_xzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xyzz_xxx, to_xyzz_xxxxy, to_xyzz_xxxyy, to_xyzz_xxxyz, to_xyzz_xxy, to_xyzz_xxyyy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xyy, to_xyzz_xyyyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_yyy, to_xyzz_yyyyy, to_xyzz_yyyyz, to_xyzz_yyyzz, to_xyzz_yyz, to_xyzz_yyzzz, to_xyzz_yzz, to_xyzz_yzzzz, to_xyzz_zzz, to_y_y_xzz_xxxx, to_y_y_xzz_xxxy, to_y_y_xzz_xxxz, to_y_y_xzz_xxyy, to_y_y_xzz_xxyz, to_y_y_xzz_xxzz, to_y_y_xzz_xyyy, to_y_y_xzz_xyyz, to_y_y_xzz_xyzz, to_y_y_xzz_xzzz, to_y_y_xzz_yyyy, to_y_y_xzz_yyyz, to_y_y_xzz_yyzz, to_y_y_xzz_yzzz, to_y_y_xzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xzz_xxxx[k] = 4.0 * to_xyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_xzz_xxxy[k] = -2.0 * to_xyzz_xxx[k] * tbe_0 + 4.0 * to_xyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_xzz_xxxz[k] = 4.0 * to_xyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_xzz_xxyy[k] = -4.0 * to_xyzz_xxy[k] * tbe_0 + 4.0 * to_xyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_xzz_xxyz[k] = -2.0 * to_xyzz_xxz[k] * tbe_0 + 4.0 * to_xyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_xzz_xxzz[k] = 4.0 * to_xyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_xzz_xyyy[k] = -6.0 * to_xyzz_xyy[k] * tbe_0 + 4.0 * to_xyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_xzz_xyyz[k] = -4.0 * to_xyzz_xyz[k] * tbe_0 + 4.0 * to_xyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_xzz_xyzz[k] = -2.0 * to_xyzz_xzz[k] * tbe_0 + 4.0 * to_xyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_xzz_xzzz[k] = 4.0 * to_xyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_xzz_yyyy[k] = -8.0 * to_xyzz_yyy[k] * tbe_0 + 4.0 * to_xyzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_xzz_yyyz[k] = -6.0 * to_xyzz_yyz[k] * tbe_0 + 4.0 * to_xyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_xzz_yyzz[k] = -4.0 * to_xyzz_yzz[k] * tbe_0 + 4.0 * to_xyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_xzz_yzzz[k] = -2.0 * to_xyzz_zzz[k] * tbe_0 + 4.0 * to_xyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_xzz_zzzz[k] = 4.0 * to_xyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 690-705 components of targeted buffer : FG

        auto to_y_y_yyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 90);

        auto to_y_y_yyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 91);

        auto to_y_y_yyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 92);

        auto to_y_y_yyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 93);

        auto to_y_y_yyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 94);

        auto to_y_y_yyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 95);

        auto to_y_y_yyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 96);

        auto to_y_y_yyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 97);

        auto to_y_y_yyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 98);

        auto to_y_y_yyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 99);

        auto to_y_y_yyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 100);

        auto to_y_y_yyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 101);

        auto to_y_y_yyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 102);

        auto to_y_y_yyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 103);

        auto to_y_y_yyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_y_y_yyy_xxxx, to_y_y_yyy_xxxy, to_y_y_yyy_xxxz, to_y_y_yyy_xxyy, to_y_y_yyy_xxyz, to_y_y_yyy_xxzz, to_y_y_yyy_xyyy, to_y_y_yyy_xyyz, to_y_y_yyy_xyzz, to_y_y_yyy_xzzz, to_y_y_yyy_yyyy, to_y_y_yyy_yyyz, to_y_y_yyy_yyzz, to_y_y_yyy_yzzz, to_y_y_yyy_zzzz, to_yy_xxx, to_yy_xxxxy, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_yyy, to_yy_yyyyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, to_yyyy_xxx, to_yyyy_xxxxy, to_yyyy_xxxyy, to_yyyy_xxxyz, to_yyyy_xxy, to_yyyy_xxyyy, to_yyyy_xxyyz, to_yyyy_xxyzz, to_yyyy_xxz, to_yyyy_xyy, to_yyyy_xyyyy, to_yyyy_xyyyz, to_yyyy_xyyzz, to_yyyy_xyz, to_yyyy_xyzzz, to_yyyy_xzz, to_yyyy_yyy, to_yyyy_yyyyy, to_yyyy_yyyyz, to_yyyy_yyyzz, to_yyyy_yyz, to_yyyy_yyzzz, to_yyyy_yzz, to_yyyy_yzzzz, to_yyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyy_xxxx[k] = -6.0 * to_yy_xxxxy[k] * tke_0 + 4.0 * to_yyyy_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_yyy_xxxy[k] = 3.0 * to_yy_xxx[k] - 6.0 * to_yy_xxxyy[k] * tke_0 - 2.0 * to_yyyy_xxx[k] * tbe_0 + 4.0 * to_yyyy_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_yyy_xxxz[k] = -6.0 * to_yy_xxxyz[k] * tke_0 + 4.0 * to_yyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_yyy_xxyy[k] = 6.0 * to_yy_xxy[k] - 6.0 * to_yy_xxyyy[k] * tke_0 - 4.0 * to_yyyy_xxy[k] * tbe_0 + 4.0 * to_yyyy_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_yyy_xxyz[k] = 3.0 * to_yy_xxz[k] - 6.0 * to_yy_xxyyz[k] * tke_0 - 2.0 * to_yyyy_xxz[k] * tbe_0 + 4.0 * to_yyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_yyy_xxzz[k] = -6.0 * to_yy_xxyzz[k] * tke_0 + 4.0 * to_yyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_yyy_xyyy[k] = 9.0 * to_yy_xyy[k] - 6.0 * to_yy_xyyyy[k] * tke_0 - 6.0 * to_yyyy_xyy[k] * tbe_0 + 4.0 * to_yyyy_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyy_xyyz[k] = 6.0 * to_yy_xyz[k] - 6.0 * to_yy_xyyyz[k] * tke_0 - 4.0 * to_yyyy_xyz[k] * tbe_0 + 4.0 * to_yyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyy_xyzz[k] = 3.0 * to_yy_xzz[k] - 6.0 * to_yy_xyyzz[k] * tke_0 - 2.0 * to_yyyy_xzz[k] * tbe_0 + 4.0 * to_yyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyy_xzzz[k] = -6.0 * to_yy_xyzzz[k] * tke_0 + 4.0 * to_yyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyy_yyyy[k] = 12.0 * to_yy_yyy[k] - 6.0 * to_yy_yyyyy[k] * tke_0 - 8.0 * to_yyyy_yyy[k] * tbe_0 + 4.0 * to_yyyy_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyy_yyyz[k] = 9.0 * to_yy_yyz[k] - 6.0 * to_yy_yyyyz[k] * tke_0 - 6.0 * to_yyyy_yyz[k] * tbe_0 + 4.0 * to_yyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyy_yyzz[k] = 6.0 * to_yy_yzz[k] - 6.0 * to_yy_yyyzz[k] * tke_0 - 4.0 * to_yyyy_yzz[k] * tbe_0 + 4.0 * to_yyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyy_yzzz[k] = 3.0 * to_yy_zzz[k] - 6.0 * to_yy_yyzzz[k] * tke_0 - 2.0 * to_yyyy_zzz[k] * tbe_0 + 4.0 * to_yyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyy_zzzz[k] = -6.0 * to_yy_yzzzz[k] * tke_0 + 4.0 * to_yyyy_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 705-720 components of targeted buffer : FG

        auto to_y_y_yyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 105);

        auto to_y_y_yyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 106);

        auto to_y_y_yyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 107);

        auto to_y_y_yyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 108);

        auto to_y_y_yyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 109);

        auto to_y_y_yyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 110);

        auto to_y_y_yyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 111);

        auto to_y_y_yyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 112);

        auto to_y_y_yyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 113);

        auto to_y_y_yyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 114);

        auto to_y_y_yyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 115);

        auto to_y_y_yyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 116);

        auto to_y_y_yyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 117);

        auto to_y_y_yyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 118);

        auto to_y_y_yyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_y_y_yyz_xxxx, to_y_y_yyz_xxxy, to_y_y_yyz_xxxz, to_y_y_yyz_xxyy, to_y_y_yyz_xxyz, to_y_y_yyz_xxzz, to_y_y_yyz_xyyy, to_y_y_yyz_xyyz, to_y_y_yyz_xyzz, to_y_y_yyz_xzzz, to_y_y_yyz_yyyy, to_y_y_yyz_yyyz, to_y_y_yyz_yyzz, to_y_y_yyz_yzzz, to_y_y_yyz_zzzz, to_yyyz_xxx, to_yyyz_xxxxy, to_yyyz_xxxyy, to_yyyz_xxxyz, to_yyyz_xxy, to_yyyz_xxyyy, to_yyyz_xxyyz, to_yyyz_xxyzz, to_yyyz_xxz, to_yyyz_xyy, to_yyyz_xyyyy, to_yyyz_xyyyz, to_yyyz_xyyzz, to_yyyz_xyz, to_yyyz_xyzzz, to_yyyz_xzz, to_yyyz_yyy, to_yyyz_yyyyy, to_yyyz_yyyyz, to_yyyz_yyyzz, to_yyyz_yyz, to_yyyz_yyzzz, to_yyyz_yzz, to_yyyz_yzzzz, to_yyyz_zzz, to_yz_xxx, to_yz_xxxxy, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_yyy, to_yz_yyyyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyz_xxxx[k] = -4.0 * to_yz_xxxxy[k] * tke_0 + 4.0 * to_yyyz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_yyz_xxxy[k] = 2.0 * to_yz_xxx[k] - 4.0 * to_yz_xxxyy[k] * tke_0 - 2.0 * to_yyyz_xxx[k] * tbe_0 + 4.0 * to_yyyz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_yyz_xxxz[k] = -4.0 * to_yz_xxxyz[k] * tke_0 + 4.0 * to_yyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_yyz_xxyy[k] = 4.0 * to_yz_xxy[k] - 4.0 * to_yz_xxyyy[k] * tke_0 - 4.0 * to_yyyz_xxy[k] * tbe_0 + 4.0 * to_yyyz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_yyz_xxyz[k] = 2.0 * to_yz_xxz[k] - 4.0 * to_yz_xxyyz[k] * tke_0 - 2.0 * to_yyyz_xxz[k] * tbe_0 + 4.0 * to_yyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_yyz_xxzz[k] = -4.0 * to_yz_xxyzz[k] * tke_0 + 4.0 * to_yyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_yyz_xyyy[k] = 6.0 * to_yz_xyy[k] - 4.0 * to_yz_xyyyy[k] * tke_0 - 6.0 * to_yyyz_xyy[k] * tbe_0 + 4.0 * to_yyyz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyz_xyyz[k] = 4.0 * to_yz_xyz[k] - 4.0 * to_yz_xyyyz[k] * tke_0 - 4.0 * to_yyyz_xyz[k] * tbe_0 + 4.0 * to_yyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyz_xyzz[k] = 2.0 * to_yz_xzz[k] - 4.0 * to_yz_xyyzz[k] * tke_0 - 2.0 * to_yyyz_xzz[k] * tbe_0 + 4.0 * to_yyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyz_xzzz[k] = -4.0 * to_yz_xyzzz[k] * tke_0 + 4.0 * to_yyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyz_yyyy[k] = 8.0 * to_yz_yyy[k] - 4.0 * to_yz_yyyyy[k] * tke_0 - 8.0 * to_yyyz_yyy[k] * tbe_0 + 4.0 * to_yyyz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_yyz_yyyz[k] = 6.0 * to_yz_yyz[k] - 4.0 * to_yz_yyyyz[k] * tke_0 - 6.0 * to_yyyz_yyz[k] * tbe_0 + 4.0 * to_yyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_yyz_yyzz[k] = 4.0 * to_yz_yzz[k] - 4.0 * to_yz_yyyzz[k] * tke_0 - 4.0 * to_yyyz_yzz[k] * tbe_0 + 4.0 * to_yyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_yyz_yzzz[k] = 2.0 * to_yz_zzz[k] - 4.0 * to_yz_yyzzz[k] * tke_0 - 2.0 * to_yyyz_zzz[k] * tbe_0 + 4.0 * to_yyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_yyz_zzzz[k] = -4.0 * to_yz_yzzzz[k] * tke_0 + 4.0 * to_yyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 720-735 components of targeted buffer : FG

        auto to_y_y_yzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 120);

        auto to_y_y_yzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 121);

        auto to_y_y_yzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 122);

        auto to_y_y_yzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 123);

        auto to_y_y_yzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 124);

        auto to_y_y_yzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 125);

        auto to_y_y_yzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 126);

        auto to_y_y_yzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 127);

        auto to_y_y_yzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 128);

        auto to_y_y_yzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 129);

        auto to_y_y_yzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 130);

        auto to_y_y_yzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 131);

        auto to_y_y_yzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 132);

        auto to_y_y_yzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 133);

        auto to_y_y_yzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_y_y_yzz_xxxx, to_y_y_yzz_xxxy, to_y_y_yzz_xxxz, to_y_y_yzz_xxyy, to_y_y_yzz_xxyz, to_y_y_yzz_xxzz, to_y_y_yzz_xyyy, to_y_y_yzz_xyyz, to_y_y_yzz_xyzz, to_y_y_yzz_xzzz, to_y_y_yzz_yyyy, to_y_y_yzz_yyyz, to_y_y_yzz_yyzz, to_y_y_yzz_yzzz, to_y_y_yzz_zzzz, to_yyzz_xxx, to_yyzz_xxxxy, to_yyzz_xxxyy, to_yyzz_xxxyz, to_yyzz_xxy, to_yyzz_xxyyy, to_yyzz_xxyyz, to_yyzz_xxyzz, to_yyzz_xxz, to_yyzz_xyy, to_yyzz_xyyyy, to_yyzz_xyyyz, to_yyzz_xyyzz, to_yyzz_xyz, to_yyzz_xyzzz, to_yyzz_xzz, to_yyzz_yyy, to_yyzz_yyyyy, to_yyzz_yyyyz, to_yyzz_yyyzz, to_yyzz_yyz, to_yyzz_yyzzz, to_yyzz_yzz, to_yyzz_yzzzz, to_yyzz_zzz, to_zz_xxx, to_zz_xxxxy, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_yyy, to_zz_yyyyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yzz_xxxx[k] = -2.0 * to_zz_xxxxy[k] * tke_0 + 4.0 * to_yyzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_yzz_xxxy[k] = to_zz_xxx[k] - 2.0 * to_zz_xxxyy[k] * tke_0 - 2.0 * to_yyzz_xxx[k] * tbe_0 + 4.0 * to_yyzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_yzz_xxxz[k] = -2.0 * to_zz_xxxyz[k] * tke_0 + 4.0 * to_yyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_yzz_xxyy[k] = 2.0 * to_zz_xxy[k] - 2.0 * to_zz_xxyyy[k] * tke_0 - 4.0 * to_yyzz_xxy[k] * tbe_0 + 4.0 * to_yyzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_yzz_xxyz[k] = to_zz_xxz[k] - 2.0 * to_zz_xxyyz[k] * tke_0 - 2.0 * to_yyzz_xxz[k] * tbe_0 + 4.0 * to_yyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_yzz_xxzz[k] = -2.0 * to_zz_xxyzz[k] * tke_0 + 4.0 * to_yyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_yzz_xyyy[k] = 3.0 * to_zz_xyy[k] - 2.0 * to_zz_xyyyy[k] * tke_0 - 6.0 * to_yyzz_xyy[k] * tbe_0 + 4.0 * to_yyzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_yzz_xyyz[k] = 2.0 * to_zz_xyz[k] - 2.0 * to_zz_xyyyz[k] * tke_0 - 4.0 * to_yyzz_xyz[k] * tbe_0 + 4.0 * to_yyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_yzz_xyzz[k] = to_zz_xzz[k] - 2.0 * to_zz_xyyzz[k] * tke_0 - 2.0 * to_yyzz_xzz[k] * tbe_0 + 4.0 * to_yyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_yzz_xzzz[k] = -2.0 * to_zz_xyzzz[k] * tke_0 + 4.0 * to_yyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_yzz_yyyy[k] = 4.0 * to_zz_yyy[k] - 2.0 * to_zz_yyyyy[k] * tke_0 - 8.0 * to_yyzz_yyy[k] * tbe_0 + 4.0 * to_yyzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_yzz_yyyz[k] = 3.0 * to_zz_yyz[k] - 2.0 * to_zz_yyyyz[k] * tke_0 - 6.0 * to_yyzz_yyz[k] * tbe_0 + 4.0 * to_yyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_yzz_yyzz[k] = 2.0 * to_zz_yzz[k] - 2.0 * to_zz_yyyzz[k] * tke_0 - 4.0 * to_yyzz_yzz[k] * tbe_0 + 4.0 * to_yyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_yzz_yzzz[k] = to_zz_zzz[k] - 2.0 * to_zz_yyzzz[k] * tke_0 - 2.0 * to_yyzz_zzz[k] * tbe_0 + 4.0 * to_yyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_yzz_zzzz[k] = -2.0 * to_zz_yzzzz[k] * tke_0 + 4.0 * to_yyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 735-750 components of targeted buffer : FG

        auto to_y_y_zzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 135);

        auto to_y_y_zzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 136);

        auto to_y_y_zzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 137);

        auto to_y_y_zzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 138);

        auto to_y_y_zzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 139);

        auto to_y_y_zzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 140);

        auto to_y_y_zzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 141);

        auto to_y_y_zzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 142);

        auto to_y_y_zzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 143);

        auto to_y_y_zzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 144);

        auto to_y_y_zzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 145);

        auto to_y_y_zzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 146);

        auto to_y_y_zzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 147);

        auto to_y_y_zzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 148);

        auto to_y_y_zzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 4 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_y_y_zzz_xxxx, to_y_y_zzz_xxxy, to_y_y_zzz_xxxz, to_y_y_zzz_xxyy, to_y_y_zzz_xxyz, to_y_y_zzz_xxzz, to_y_y_zzz_xyyy, to_y_y_zzz_xyyz, to_y_y_zzz_xyzz, to_y_y_zzz_xzzz, to_y_y_zzz_yyyy, to_y_y_zzz_yyyz, to_y_y_zzz_yyzz, to_y_y_zzz_yzzz, to_y_y_zzz_zzzz, to_yzzz_xxx, to_yzzz_xxxxy, to_yzzz_xxxyy, to_yzzz_xxxyz, to_yzzz_xxy, to_yzzz_xxyyy, to_yzzz_xxyyz, to_yzzz_xxyzz, to_yzzz_xxz, to_yzzz_xyy, to_yzzz_xyyyy, to_yzzz_xyyyz, to_yzzz_xyyzz, to_yzzz_xyz, to_yzzz_xyzzz, to_yzzz_xzz, to_yzzz_yyy, to_yzzz_yyyyy, to_yzzz_yyyyz, to_yzzz_yyyzz, to_yzzz_yyz, to_yzzz_yyzzz, to_yzzz_yzz, to_yzzz_yzzzz, to_yzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zzz_xxxx[k] = 4.0 * to_yzzz_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_zzz_xxxy[k] = -2.0 * to_yzzz_xxx[k] * tbe_0 + 4.0 * to_yzzz_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_zzz_xxxz[k] = 4.0 * to_yzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_zzz_xxyy[k] = -4.0 * to_yzzz_xxy[k] * tbe_0 + 4.0 * to_yzzz_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_zzz_xxyz[k] = -2.0 * to_yzzz_xxz[k] * tbe_0 + 4.0 * to_yzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_zzz_xxzz[k] = 4.0 * to_yzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_zzz_xyyy[k] = -6.0 * to_yzzz_xyy[k] * tbe_0 + 4.0 * to_yzzz_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_zzz_xyyz[k] = -4.0 * to_yzzz_xyz[k] * tbe_0 + 4.0 * to_yzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_zzz_xyzz[k] = -2.0 * to_yzzz_xzz[k] * tbe_0 + 4.0 * to_yzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_zzz_xzzz[k] = 4.0 * to_yzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_zzz_yyyy[k] = -8.0 * to_yzzz_yyy[k] * tbe_0 + 4.0 * to_yzzz_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_zzz_yyyz[k] = -6.0 * to_yzzz_yyz[k] * tbe_0 + 4.0 * to_yzzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_zzz_yyzz[k] = -4.0 * to_yzzz_yzz[k] * tbe_0 + 4.0 * to_yzzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_zzz_yzzz[k] = -2.0 * to_yzzz_zzz[k] * tbe_0 + 4.0 * to_yzzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_zzz_zzzz[k] = 4.0 * to_yzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 750-765 components of targeted buffer : FG

        auto to_y_z_xxx_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 0);

        auto to_y_z_xxx_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 1);

        auto to_y_z_xxx_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 2);

        auto to_y_z_xxx_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 3);

        auto to_y_z_xxx_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 4);

        auto to_y_z_xxx_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 5);

        auto to_y_z_xxx_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 6);

        auto to_y_z_xxx_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 7);

        auto to_y_z_xxx_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 8);

        auto to_y_z_xxx_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 9);

        auto to_y_z_xxx_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 10);

        auto to_y_z_xxx_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 11);

        auto to_y_z_xxx_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 12);

        auto to_y_z_xxx_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 13);

        auto to_y_z_xxx_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_xxxy_xxx, to_xxxy_xxxxz, to_xxxy_xxxyz, to_xxxy_xxxzz, to_xxxy_xxy, to_xxxy_xxyyz, to_xxxy_xxyzz, to_xxxy_xxz, to_xxxy_xxzzz, to_xxxy_xyy, to_xxxy_xyyyz, to_xxxy_xyyzz, to_xxxy_xyz, to_xxxy_xyzzz, to_xxxy_xzz, to_xxxy_xzzzz, to_xxxy_yyy, to_xxxy_yyyyz, to_xxxy_yyyzz, to_xxxy_yyz, to_xxxy_yyzzz, to_xxxy_yzz, to_xxxy_yzzzz, to_xxxy_zzz, to_xxxy_zzzzz, to_y_z_xxx_xxxx, to_y_z_xxx_xxxy, to_y_z_xxx_xxxz, to_y_z_xxx_xxyy, to_y_z_xxx_xxyz, to_y_z_xxx_xxzz, to_y_z_xxx_xyyy, to_y_z_xxx_xyyz, to_y_z_xxx_xyzz, to_y_z_xxx_xzzz, to_y_z_xxx_yyyy, to_y_z_xxx_yyyz, to_y_z_xxx_yyzz, to_y_z_xxx_yzzz, to_y_z_xxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxx_xxxx[k] = 4.0 * to_xxxy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xxxy[k] = 4.0 * to_xxxy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xxxz[k] = -2.0 * to_xxxy_xxx[k] * tbe_0 + 4.0 * to_xxxy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xxyy[k] = 4.0 * to_xxxy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xxyz[k] = -2.0 * to_xxxy_xxy[k] * tbe_0 + 4.0 * to_xxxy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xxzz[k] = -4.0 * to_xxxy_xxz[k] * tbe_0 + 4.0 * to_xxxy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xyyy[k] = 4.0 * to_xxxy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xyyz[k] = -2.0 * to_xxxy_xyy[k] * tbe_0 + 4.0 * to_xxxy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xyzz[k] = -4.0 * to_xxxy_xyz[k] * tbe_0 + 4.0 * to_xxxy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xzzz[k] = -6.0 * to_xxxy_xzz[k] * tbe_0 + 4.0 * to_xxxy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_yyyy[k] = 4.0 * to_xxxy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxx_yyyz[k] = -2.0 * to_xxxy_yyy[k] * tbe_0 + 4.0 * to_xxxy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_yyzz[k] = -4.0 * to_xxxy_yyz[k] * tbe_0 + 4.0 * to_xxxy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_yzzz[k] = -6.0 * to_xxxy_yzz[k] * tbe_0 + 4.0 * to_xxxy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_zzzz[k] = -8.0 * to_xxxy_zzz[k] * tbe_0 + 4.0 * to_xxxy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 765-780 components of targeted buffer : FG

        auto to_y_z_xxy_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 15);

        auto to_y_z_xxy_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 16);

        auto to_y_z_xxy_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 17);

        auto to_y_z_xxy_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 18);

        auto to_y_z_xxy_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 19);

        auto to_y_z_xxy_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 20);

        auto to_y_z_xxy_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 21);

        auto to_y_z_xxy_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 22);

        auto to_y_z_xxy_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 23);

        auto to_y_z_xxy_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 24);

        auto to_y_z_xxy_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 25);

        auto to_y_z_xxy_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 26);

        auto to_y_z_xxy_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 27);

        auto to_y_z_xxy_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 28);

        auto to_y_z_xxy_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xx_xxx, to_xx_xxxxz, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, to_xx_zzzzz, to_xxyy_xxx, to_xxyy_xxxxz, to_xxyy_xxxyz, to_xxyy_xxxzz, to_xxyy_xxy, to_xxyy_xxyyz, to_xxyy_xxyzz, to_xxyy_xxz, to_xxyy_xxzzz, to_xxyy_xyy, to_xxyy_xyyyz, to_xxyy_xyyzz, to_xxyy_xyz, to_xxyy_xyzzz, to_xxyy_xzz, to_xxyy_xzzzz, to_xxyy_yyy, to_xxyy_yyyyz, to_xxyy_yyyzz, to_xxyy_yyz, to_xxyy_yyzzz, to_xxyy_yzz, to_xxyy_yzzzz, to_xxyy_zzz, to_xxyy_zzzzz, to_y_z_xxy_xxxx, to_y_z_xxy_xxxy, to_y_z_xxy_xxxz, to_y_z_xxy_xxyy, to_y_z_xxy_xxyz, to_y_z_xxy_xxzz, to_y_z_xxy_xyyy, to_y_z_xxy_xyyz, to_y_z_xxy_xyzz, to_y_z_xxy_xzzz, to_y_z_xxy_yyyy, to_y_z_xxy_yyyz, to_y_z_xxy_yyzz, to_y_z_xxy_yzzz, to_y_z_xxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxy_xxxx[k] = -2.0 * to_xx_xxxxz[k] * tke_0 + 4.0 * to_xxyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xxxy[k] = -2.0 * to_xx_xxxyz[k] * tke_0 + 4.0 * to_xxyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xxxz[k] = to_xx_xxx[k] - 2.0 * to_xx_xxxzz[k] * tke_0 - 2.0 * to_xxyy_xxx[k] * tbe_0 + 4.0 * to_xxyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xxyy[k] = -2.0 * to_xx_xxyyz[k] * tke_0 + 4.0 * to_xxyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xxyz[k] = to_xx_xxy[k] - 2.0 * to_xx_xxyzz[k] * tke_0 - 2.0 * to_xxyy_xxy[k] * tbe_0 + 4.0 * to_xxyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xxzz[k] = 2.0 * to_xx_xxz[k] - 2.0 * to_xx_xxzzz[k] * tke_0 - 4.0 * to_xxyy_xxz[k] * tbe_0 + 4.0 * to_xxyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xyyy[k] = -2.0 * to_xx_xyyyz[k] * tke_0 + 4.0 * to_xxyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xyyz[k] = to_xx_xyy[k] - 2.0 * to_xx_xyyzz[k] * tke_0 - 2.0 * to_xxyy_xyy[k] * tbe_0 + 4.0 * to_xxyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xyzz[k] = 2.0 * to_xx_xyz[k] - 2.0 * to_xx_xyzzz[k] * tke_0 - 4.0 * to_xxyy_xyz[k] * tbe_0 + 4.0 * to_xxyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xzzz[k] = 3.0 * to_xx_xzz[k] - 2.0 * to_xx_xzzzz[k] * tke_0 - 6.0 * to_xxyy_xzz[k] * tbe_0 + 4.0 * to_xxyy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_yyyy[k] = -2.0 * to_xx_yyyyz[k] * tke_0 + 4.0 * to_xxyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxy_yyyz[k] = to_xx_yyy[k] - 2.0 * to_xx_yyyzz[k] * tke_0 - 2.0 * to_xxyy_yyy[k] * tbe_0 + 4.0 * to_xxyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_yyzz[k] = 2.0 * to_xx_yyz[k] - 2.0 * to_xx_yyzzz[k] * tke_0 - 4.0 * to_xxyy_yyz[k] * tbe_0 + 4.0 * to_xxyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_yzzz[k] = 3.0 * to_xx_yzz[k] - 2.0 * to_xx_yzzzz[k] * tke_0 - 6.0 * to_xxyy_yzz[k] * tbe_0 + 4.0 * to_xxyy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_zzzz[k] = 4.0 * to_xx_zzz[k] - 2.0 * to_xx_zzzzz[k] * tke_0 - 8.0 * to_xxyy_zzz[k] * tbe_0 + 4.0 * to_xxyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 780-795 components of targeted buffer : FG

        auto to_y_z_xxz_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 30);

        auto to_y_z_xxz_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 31);

        auto to_y_z_xxz_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 32);

        auto to_y_z_xxz_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 33);

        auto to_y_z_xxz_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 34);

        auto to_y_z_xxz_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 35);

        auto to_y_z_xxz_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 36);

        auto to_y_z_xxz_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 37);

        auto to_y_z_xxz_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 38);

        auto to_y_z_xxz_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 39);

        auto to_y_z_xxz_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 40);

        auto to_y_z_xxz_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 41);

        auto to_y_z_xxz_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 42);

        auto to_y_z_xxz_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 43);

        auto to_y_z_xxz_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_xxyz_xxx, to_xxyz_xxxxz, to_xxyz_xxxyz, to_xxyz_xxxzz, to_xxyz_xxy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xxzzz, to_xxyz_xyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_xzzzz, to_xxyz_yyy, to_xxyz_yyyyz, to_xxyz_yyyzz, to_xxyz_yyz, to_xxyz_yyzzz, to_xxyz_yzz, to_xxyz_yzzzz, to_xxyz_zzz, to_xxyz_zzzzz, to_y_z_xxz_xxxx, to_y_z_xxz_xxxy, to_y_z_xxz_xxxz, to_y_z_xxz_xxyy, to_y_z_xxz_xxyz, to_y_z_xxz_xxzz, to_y_z_xxz_xyyy, to_y_z_xxz_xyyz, to_y_z_xxz_xyzz, to_y_z_xxz_xzzz, to_y_z_xxz_yyyy, to_y_z_xxz_yyyz, to_y_z_xxz_yyzz, to_y_z_xxz_yzzz, to_y_z_xxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxz_xxxx[k] = 4.0 * to_xxyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xxxy[k] = 4.0 * to_xxyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xxxz[k] = -2.0 * to_xxyz_xxx[k] * tbe_0 + 4.0 * to_xxyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xxyy[k] = 4.0 * to_xxyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xxyz[k] = -2.0 * to_xxyz_xxy[k] * tbe_0 + 4.0 * to_xxyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xxzz[k] = -4.0 * to_xxyz_xxz[k] * tbe_0 + 4.0 * to_xxyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xyyy[k] = 4.0 * to_xxyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xyyz[k] = -2.0 * to_xxyz_xyy[k] * tbe_0 + 4.0 * to_xxyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xyzz[k] = -4.0 * to_xxyz_xyz[k] * tbe_0 + 4.0 * to_xxyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xzzz[k] = -6.0 * to_xxyz_xzz[k] * tbe_0 + 4.0 * to_xxyz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_yyyy[k] = 4.0 * to_xxyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xxz_yyyz[k] = -2.0 * to_xxyz_yyy[k] * tbe_0 + 4.0 * to_xxyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_yyzz[k] = -4.0 * to_xxyz_yyz[k] * tbe_0 + 4.0 * to_xxyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_yzzz[k] = -6.0 * to_xxyz_yzz[k] * tbe_0 + 4.0 * to_xxyz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_zzzz[k] = -8.0 * to_xxyz_zzz[k] * tbe_0 + 4.0 * to_xxyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 795-810 components of targeted buffer : FG

        auto to_y_z_xyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 45);

        auto to_y_z_xyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 46);

        auto to_y_z_xyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 47);

        auto to_y_z_xyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 48);

        auto to_y_z_xyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 49);

        auto to_y_z_xyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 50);

        auto to_y_z_xyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 51);

        auto to_y_z_xyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 52);

        auto to_y_z_xyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 53);

        auto to_y_z_xyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 54);

        auto to_y_z_xyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 55);

        auto to_y_z_xyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 56);

        auto to_y_z_xyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 57);

        auto to_y_z_xyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 58);

        auto to_y_z_xyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xy_xxx, to_xy_xxxxz, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, to_xy_zzzzz, to_xyyy_xxx, to_xyyy_xxxxz, to_xyyy_xxxyz, to_xyyy_xxxzz, to_xyyy_xxy, to_xyyy_xxyyz, to_xyyy_xxyzz, to_xyyy_xxz, to_xyyy_xxzzz, to_xyyy_xyy, to_xyyy_xyyyz, to_xyyy_xyyzz, to_xyyy_xyz, to_xyyy_xyzzz, to_xyyy_xzz, to_xyyy_xzzzz, to_xyyy_yyy, to_xyyy_yyyyz, to_xyyy_yyyzz, to_xyyy_yyz, to_xyyy_yyzzz, to_xyyy_yzz, to_xyyy_yzzzz, to_xyyy_zzz, to_xyyy_zzzzz, to_y_z_xyy_xxxx, to_y_z_xyy_xxxy, to_y_z_xyy_xxxz, to_y_z_xyy_xxyy, to_y_z_xyy_xxyz, to_y_z_xyy_xxzz, to_y_z_xyy_xyyy, to_y_z_xyy_xyyz, to_y_z_xyy_xyzz, to_y_z_xyy_xzzz, to_y_z_xyy_yyyy, to_y_z_xyy_yyyz, to_y_z_xyy_yyzz, to_y_z_xyy_yzzz, to_y_z_xyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyy_xxxx[k] = -4.0 * to_xy_xxxxz[k] * tke_0 + 4.0 * to_xyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xxxy[k] = -4.0 * to_xy_xxxyz[k] * tke_0 + 4.0 * to_xyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xxxz[k] = 2.0 * to_xy_xxx[k] - 4.0 * to_xy_xxxzz[k] * tke_0 - 2.0 * to_xyyy_xxx[k] * tbe_0 + 4.0 * to_xyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xxyy[k] = -4.0 * to_xy_xxyyz[k] * tke_0 + 4.0 * to_xyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xxyz[k] = 2.0 * to_xy_xxy[k] - 4.0 * to_xy_xxyzz[k] * tke_0 - 2.0 * to_xyyy_xxy[k] * tbe_0 + 4.0 * to_xyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xxzz[k] = 4.0 * to_xy_xxz[k] - 4.0 * to_xy_xxzzz[k] * tke_0 - 4.0 * to_xyyy_xxz[k] * tbe_0 + 4.0 * to_xyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xyyy[k] = -4.0 * to_xy_xyyyz[k] * tke_0 + 4.0 * to_xyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xyyz[k] = 2.0 * to_xy_xyy[k] - 4.0 * to_xy_xyyzz[k] * tke_0 - 2.0 * to_xyyy_xyy[k] * tbe_0 + 4.0 * to_xyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xyzz[k] = 4.0 * to_xy_xyz[k] - 4.0 * to_xy_xyzzz[k] * tke_0 - 4.0 * to_xyyy_xyz[k] * tbe_0 + 4.0 * to_xyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xzzz[k] = 6.0 * to_xy_xzz[k] - 4.0 * to_xy_xzzzz[k] * tke_0 - 6.0 * to_xyyy_xzz[k] * tbe_0 + 4.0 * to_xyyy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_yyyy[k] = -4.0 * to_xy_yyyyz[k] * tke_0 + 4.0 * to_xyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyy_yyyz[k] = 2.0 * to_xy_yyy[k] - 4.0 * to_xy_yyyzz[k] * tke_0 - 2.0 * to_xyyy_yyy[k] * tbe_0 + 4.0 * to_xyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_yyzz[k] = 4.0 * to_xy_yyz[k] - 4.0 * to_xy_yyzzz[k] * tke_0 - 4.0 * to_xyyy_yyz[k] * tbe_0 + 4.0 * to_xyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_yzzz[k] = 6.0 * to_xy_yzz[k] - 4.0 * to_xy_yzzzz[k] * tke_0 - 6.0 * to_xyyy_yzz[k] * tbe_0 + 4.0 * to_xyyy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_zzzz[k] = 8.0 * to_xy_zzz[k] - 4.0 * to_xy_zzzzz[k] * tke_0 - 8.0 * to_xyyy_zzz[k] * tbe_0 + 4.0 * to_xyyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 810-825 components of targeted buffer : FG

        auto to_y_z_xyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 60);

        auto to_y_z_xyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 61);

        auto to_y_z_xyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 62);

        auto to_y_z_xyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 63);

        auto to_y_z_xyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 64);

        auto to_y_z_xyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 65);

        auto to_y_z_xyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 66);

        auto to_y_z_xyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 67);

        auto to_y_z_xyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 68);

        auto to_y_z_xyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 69);

        auto to_y_z_xyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 70);

        auto to_y_z_xyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 71);

        auto to_y_z_xyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 72);

        auto to_y_z_xyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 73);

        auto to_y_z_xyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_xyyz_xxx, to_xyyz_xxxxz, to_xyyz_xxxyz, to_xyyz_xxxzz, to_xyyz_xxy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xxzzz, to_xyyz_xyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_xzzzz, to_xyyz_yyy, to_xyyz_yyyyz, to_xyyz_yyyzz, to_xyyz_yyz, to_xyyz_yyzzz, to_xyyz_yzz, to_xyyz_yzzzz, to_xyyz_zzz, to_xyyz_zzzzz, to_xz_xxx, to_xz_xxxxz, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, to_xz_zzzzz, to_y_z_xyz_xxxx, to_y_z_xyz_xxxy, to_y_z_xyz_xxxz, to_y_z_xyz_xxyy, to_y_z_xyz_xxyz, to_y_z_xyz_xxzz, to_y_z_xyz_xyyy, to_y_z_xyz_xyyz, to_y_z_xyz_xyzz, to_y_z_xyz_xzzz, to_y_z_xyz_yyyy, to_y_z_xyz_yyyz, to_y_z_xyz_yyzz, to_y_z_xyz_yzzz, to_y_z_xyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyz_xxxx[k] = -2.0 * to_xz_xxxxz[k] * tke_0 + 4.0 * to_xyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xxxy[k] = -2.0 * to_xz_xxxyz[k] * tke_0 + 4.0 * to_xyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xxxz[k] = to_xz_xxx[k] - 2.0 * to_xz_xxxzz[k] * tke_0 - 2.0 * to_xyyz_xxx[k] * tbe_0 + 4.0 * to_xyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xxyy[k] = -2.0 * to_xz_xxyyz[k] * tke_0 + 4.0 * to_xyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xxyz[k] = to_xz_xxy[k] - 2.0 * to_xz_xxyzz[k] * tke_0 - 2.0 * to_xyyz_xxy[k] * tbe_0 + 4.0 * to_xyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xxzz[k] = 2.0 * to_xz_xxz[k] - 2.0 * to_xz_xxzzz[k] * tke_0 - 4.0 * to_xyyz_xxz[k] * tbe_0 + 4.0 * to_xyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xyyy[k] = -2.0 * to_xz_xyyyz[k] * tke_0 + 4.0 * to_xyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xyyz[k] = to_xz_xyy[k] - 2.0 * to_xz_xyyzz[k] * tke_0 - 2.0 * to_xyyz_xyy[k] * tbe_0 + 4.0 * to_xyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xyzz[k] = 2.0 * to_xz_xyz[k] - 2.0 * to_xz_xyzzz[k] * tke_0 - 4.0 * to_xyyz_xyz[k] * tbe_0 + 4.0 * to_xyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xzzz[k] = 3.0 * to_xz_xzz[k] - 2.0 * to_xz_xzzzz[k] * tke_0 - 6.0 * to_xyyz_xzz[k] * tbe_0 + 4.0 * to_xyyz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_yyyy[k] = -2.0 * to_xz_yyyyz[k] * tke_0 + 4.0 * to_xyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xyz_yyyz[k] = to_xz_yyy[k] - 2.0 * to_xz_yyyzz[k] * tke_0 - 2.0 * to_xyyz_yyy[k] * tbe_0 + 4.0 * to_xyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_yyzz[k] = 2.0 * to_xz_yyz[k] - 2.0 * to_xz_yyzzz[k] * tke_0 - 4.0 * to_xyyz_yyz[k] * tbe_0 + 4.0 * to_xyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_yzzz[k] = 3.0 * to_xz_yzz[k] - 2.0 * to_xz_yzzzz[k] * tke_0 - 6.0 * to_xyyz_yzz[k] * tbe_0 + 4.0 * to_xyyz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_zzzz[k] = 4.0 * to_xz_zzz[k] - 2.0 * to_xz_zzzzz[k] * tke_0 - 8.0 * to_xyyz_zzz[k] * tbe_0 + 4.0 * to_xyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 825-840 components of targeted buffer : FG

        auto to_y_z_xzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 75);

        auto to_y_z_xzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 76);

        auto to_y_z_xzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 77);

        auto to_y_z_xzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 78);

        auto to_y_z_xzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 79);

        auto to_y_z_xzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 80);

        auto to_y_z_xzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 81);

        auto to_y_z_xzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 82);

        auto to_y_z_xzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 83);

        auto to_y_z_xzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 84);

        auto to_y_z_xzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 85);

        auto to_y_z_xzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 86);

        auto to_y_z_xzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 87);

        auto to_y_z_xzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 88);

        auto to_y_z_xzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xyzz_xxx, to_xyzz_xxxxz, to_xyzz_xxxyz, to_xyzz_xxxzz, to_xyzz_xxy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xxzzz, to_xyzz_xyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_xzzzz, to_xyzz_yyy, to_xyzz_yyyyz, to_xyzz_yyyzz, to_xyzz_yyz, to_xyzz_yyzzz, to_xyzz_yzz, to_xyzz_yzzzz, to_xyzz_zzz, to_xyzz_zzzzz, to_y_z_xzz_xxxx, to_y_z_xzz_xxxy, to_y_z_xzz_xxxz, to_y_z_xzz_xxyy, to_y_z_xzz_xxyz, to_y_z_xzz_xxzz, to_y_z_xzz_xyyy, to_y_z_xzz_xyyz, to_y_z_xzz_xyzz, to_y_z_xzz_xzzz, to_y_z_xzz_yyyy, to_y_z_xzz_yyyz, to_y_z_xzz_yyzz, to_y_z_xzz_yzzz, to_y_z_xzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xzz_xxxx[k] = 4.0 * to_xyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xxxy[k] = 4.0 * to_xyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xxxz[k] = -2.0 * to_xyzz_xxx[k] * tbe_0 + 4.0 * to_xyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xxyy[k] = 4.0 * to_xyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xxyz[k] = -2.0 * to_xyzz_xxy[k] * tbe_0 + 4.0 * to_xyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xxzz[k] = -4.0 * to_xyzz_xxz[k] * tbe_0 + 4.0 * to_xyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xyyy[k] = 4.0 * to_xyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xyyz[k] = -2.0 * to_xyzz_xyy[k] * tbe_0 + 4.0 * to_xyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xyzz[k] = -4.0 * to_xyzz_xyz[k] * tbe_0 + 4.0 * to_xyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xzzz[k] = -6.0 * to_xyzz_xzz[k] * tbe_0 + 4.0 * to_xyzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_yyyy[k] = 4.0 * to_xyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_xzz_yyyz[k] = -2.0 * to_xyzz_yyy[k] * tbe_0 + 4.0 * to_xyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_yyzz[k] = -4.0 * to_xyzz_yyz[k] * tbe_0 + 4.0 * to_xyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_yzzz[k] = -6.0 * to_xyzz_yzz[k] * tbe_0 + 4.0 * to_xyzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_zzzz[k] = -8.0 * to_xyzz_zzz[k] * tbe_0 + 4.0 * to_xyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 840-855 components of targeted buffer : FG

        auto to_y_z_yyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 90);

        auto to_y_z_yyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 91);

        auto to_y_z_yyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 92);

        auto to_y_z_yyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 93);

        auto to_y_z_yyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 94);

        auto to_y_z_yyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 95);

        auto to_y_z_yyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 96);

        auto to_y_z_yyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 97);

        auto to_y_z_yyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 98);

        auto to_y_z_yyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 99);

        auto to_y_z_yyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 100);

        auto to_y_z_yyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 101);

        auto to_y_z_yyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 102);

        auto to_y_z_yyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 103);

        auto to_y_z_yyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_y_z_yyy_xxxx, to_y_z_yyy_xxxy, to_y_z_yyy_xxxz, to_y_z_yyy_xxyy, to_y_z_yyy_xxyz, to_y_z_yyy_xxzz, to_y_z_yyy_xyyy, to_y_z_yyy_xyyz, to_y_z_yyy_xyzz, to_y_z_yyy_xzzz, to_y_z_yyy_yyyy, to_y_z_yyy_yyyz, to_y_z_yyy_yyzz, to_y_z_yyy_yzzz, to_y_z_yyy_zzzz, to_yy_xxx, to_yy_xxxxz, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, to_yy_zzzzz, to_yyyy_xxx, to_yyyy_xxxxz, to_yyyy_xxxyz, to_yyyy_xxxzz, to_yyyy_xxy, to_yyyy_xxyyz, to_yyyy_xxyzz, to_yyyy_xxz, to_yyyy_xxzzz, to_yyyy_xyy, to_yyyy_xyyyz, to_yyyy_xyyzz, to_yyyy_xyz, to_yyyy_xyzzz, to_yyyy_xzz, to_yyyy_xzzzz, to_yyyy_yyy, to_yyyy_yyyyz, to_yyyy_yyyzz, to_yyyy_yyz, to_yyyy_yyzzz, to_yyyy_yzz, to_yyyy_yzzzz, to_yyyy_zzz, to_yyyy_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyy_xxxx[k] = -6.0 * to_yy_xxxxz[k] * tke_0 + 4.0 * to_yyyy_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xxxy[k] = -6.0 * to_yy_xxxyz[k] * tke_0 + 4.0 * to_yyyy_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xxxz[k] = 3.0 * to_yy_xxx[k] - 6.0 * to_yy_xxxzz[k] * tke_0 - 2.0 * to_yyyy_xxx[k] * tbe_0 + 4.0 * to_yyyy_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xxyy[k] = -6.0 * to_yy_xxyyz[k] * tke_0 + 4.0 * to_yyyy_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xxyz[k] = 3.0 * to_yy_xxy[k] - 6.0 * to_yy_xxyzz[k] * tke_0 - 2.0 * to_yyyy_xxy[k] * tbe_0 + 4.0 * to_yyyy_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xxzz[k] = 6.0 * to_yy_xxz[k] - 6.0 * to_yy_xxzzz[k] * tke_0 - 4.0 * to_yyyy_xxz[k] * tbe_0 + 4.0 * to_yyyy_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xyyy[k] = -6.0 * to_yy_xyyyz[k] * tke_0 + 4.0 * to_yyyy_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xyyz[k] = 3.0 * to_yy_xyy[k] - 6.0 * to_yy_xyyzz[k] * tke_0 - 2.0 * to_yyyy_xyy[k] * tbe_0 + 4.0 * to_yyyy_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xyzz[k] = 6.0 * to_yy_xyz[k] - 6.0 * to_yy_xyzzz[k] * tke_0 - 4.0 * to_yyyy_xyz[k] * tbe_0 + 4.0 * to_yyyy_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xzzz[k] = 9.0 * to_yy_xzz[k] - 6.0 * to_yy_xzzzz[k] * tke_0 - 6.0 * to_yyyy_xzz[k] * tbe_0 + 4.0 * to_yyyy_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_yyyy[k] = -6.0 * to_yy_yyyyz[k] * tke_0 + 4.0 * to_yyyy_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyy_yyyz[k] = 3.0 * to_yy_yyy[k] - 6.0 * to_yy_yyyzz[k] * tke_0 - 2.0 * to_yyyy_yyy[k] * tbe_0 + 4.0 * to_yyyy_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_yyzz[k] = 6.0 * to_yy_yyz[k] - 6.0 * to_yy_yyzzz[k] * tke_0 - 4.0 * to_yyyy_yyz[k] * tbe_0 + 4.0 * to_yyyy_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_yzzz[k] = 9.0 * to_yy_yzz[k] - 6.0 * to_yy_yzzzz[k] * tke_0 - 6.0 * to_yyyy_yzz[k] * tbe_0 + 4.0 * to_yyyy_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_zzzz[k] = 12.0 * to_yy_zzz[k] - 6.0 * to_yy_zzzzz[k] * tke_0 - 8.0 * to_yyyy_zzz[k] * tbe_0 + 4.0 * to_yyyy_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 855-870 components of targeted buffer : FG

        auto to_y_z_yyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 105);

        auto to_y_z_yyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 106);

        auto to_y_z_yyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 107);

        auto to_y_z_yyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 108);

        auto to_y_z_yyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 109);

        auto to_y_z_yyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 110);

        auto to_y_z_yyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 111);

        auto to_y_z_yyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 112);

        auto to_y_z_yyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 113);

        auto to_y_z_yyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 114);

        auto to_y_z_yyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 115);

        auto to_y_z_yyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 116);

        auto to_y_z_yyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 117);

        auto to_y_z_yyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 118);

        auto to_y_z_yyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_y_z_yyz_xxxx, to_y_z_yyz_xxxy, to_y_z_yyz_xxxz, to_y_z_yyz_xxyy, to_y_z_yyz_xxyz, to_y_z_yyz_xxzz, to_y_z_yyz_xyyy, to_y_z_yyz_xyyz, to_y_z_yyz_xyzz, to_y_z_yyz_xzzz, to_y_z_yyz_yyyy, to_y_z_yyz_yyyz, to_y_z_yyz_yyzz, to_y_z_yyz_yzzz, to_y_z_yyz_zzzz, to_yyyz_xxx, to_yyyz_xxxxz, to_yyyz_xxxyz, to_yyyz_xxxzz, to_yyyz_xxy, to_yyyz_xxyyz, to_yyyz_xxyzz, to_yyyz_xxz, to_yyyz_xxzzz, to_yyyz_xyy, to_yyyz_xyyyz, to_yyyz_xyyzz, to_yyyz_xyz, to_yyyz_xyzzz, to_yyyz_xzz, to_yyyz_xzzzz, to_yyyz_yyy, to_yyyz_yyyyz, to_yyyz_yyyzz, to_yyyz_yyz, to_yyyz_yyzzz, to_yyyz_yzz, to_yyyz_yzzzz, to_yyyz_zzz, to_yyyz_zzzzz, to_yz_xxx, to_yz_xxxxz, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, to_yz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyz_xxxx[k] = -4.0 * to_yz_xxxxz[k] * tke_0 + 4.0 * to_yyyz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xxxy[k] = -4.0 * to_yz_xxxyz[k] * tke_0 + 4.0 * to_yyyz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xxxz[k] = 2.0 * to_yz_xxx[k] - 4.0 * to_yz_xxxzz[k] * tke_0 - 2.0 * to_yyyz_xxx[k] * tbe_0 + 4.0 * to_yyyz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xxyy[k] = -4.0 * to_yz_xxyyz[k] * tke_0 + 4.0 * to_yyyz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xxyz[k] = 2.0 * to_yz_xxy[k] - 4.0 * to_yz_xxyzz[k] * tke_0 - 2.0 * to_yyyz_xxy[k] * tbe_0 + 4.0 * to_yyyz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xxzz[k] = 4.0 * to_yz_xxz[k] - 4.0 * to_yz_xxzzz[k] * tke_0 - 4.0 * to_yyyz_xxz[k] * tbe_0 + 4.0 * to_yyyz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xyyy[k] = -4.0 * to_yz_xyyyz[k] * tke_0 + 4.0 * to_yyyz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xyyz[k] = 2.0 * to_yz_xyy[k] - 4.0 * to_yz_xyyzz[k] * tke_0 - 2.0 * to_yyyz_xyy[k] * tbe_0 + 4.0 * to_yyyz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xyzz[k] = 4.0 * to_yz_xyz[k] - 4.0 * to_yz_xyzzz[k] * tke_0 - 4.0 * to_yyyz_xyz[k] * tbe_0 + 4.0 * to_yyyz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xzzz[k] = 6.0 * to_yz_xzz[k] - 4.0 * to_yz_xzzzz[k] * tke_0 - 6.0 * to_yyyz_xzz[k] * tbe_0 + 4.0 * to_yyyz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_yyyy[k] = -4.0 * to_yz_yyyyz[k] * tke_0 + 4.0 * to_yyyz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_yyz_yyyz[k] = 2.0 * to_yz_yyy[k] - 4.0 * to_yz_yyyzz[k] * tke_0 - 2.0 * to_yyyz_yyy[k] * tbe_0 + 4.0 * to_yyyz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_yyzz[k] = 4.0 * to_yz_yyz[k] - 4.0 * to_yz_yyzzz[k] * tke_0 - 4.0 * to_yyyz_yyz[k] * tbe_0 + 4.0 * to_yyyz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_yzzz[k] = 6.0 * to_yz_yzz[k] - 4.0 * to_yz_yzzzz[k] * tke_0 - 6.0 * to_yyyz_yzz[k] * tbe_0 + 4.0 * to_yyyz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_zzzz[k] = 8.0 * to_yz_zzz[k] - 4.0 * to_yz_zzzzz[k] * tke_0 - 8.0 * to_yyyz_zzz[k] * tbe_0 + 4.0 * to_yyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 870-885 components of targeted buffer : FG

        auto to_y_z_yzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 120);

        auto to_y_z_yzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 121);

        auto to_y_z_yzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 122);

        auto to_y_z_yzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 123);

        auto to_y_z_yzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 124);

        auto to_y_z_yzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 125);

        auto to_y_z_yzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 126);

        auto to_y_z_yzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 127);

        auto to_y_z_yzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 128);

        auto to_y_z_yzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 129);

        auto to_y_z_yzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 130);

        auto to_y_z_yzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 131);

        auto to_y_z_yzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 132);

        auto to_y_z_yzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 133);

        auto to_y_z_yzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_y_z_yzz_xxxx, to_y_z_yzz_xxxy, to_y_z_yzz_xxxz, to_y_z_yzz_xxyy, to_y_z_yzz_xxyz, to_y_z_yzz_xxzz, to_y_z_yzz_xyyy, to_y_z_yzz_xyyz, to_y_z_yzz_xyzz, to_y_z_yzz_xzzz, to_y_z_yzz_yyyy, to_y_z_yzz_yyyz, to_y_z_yzz_yyzz, to_y_z_yzz_yzzz, to_y_z_yzz_zzzz, to_yyzz_xxx, to_yyzz_xxxxz, to_yyzz_xxxyz, to_yyzz_xxxzz, to_yyzz_xxy, to_yyzz_xxyyz, to_yyzz_xxyzz, to_yyzz_xxz, to_yyzz_xxzzz, to_yyzz_xyy, to_yyzz_xyyyz, to_yyzz_xyyzz, to_yyzz_xyz, to_yyzz_xyzzz, to_yyzz_xzz, to_yyzz_xzzzz, to_yyzz_yyy, to_yyzz_yyyyz, to_yyzz_yyyzz, to_yyzz_yyz, to_yyzz_yyzzz, to_yyzz_yzz, to_yyzz_yzzzz, to_yyzz_zzz, to_yyzz_zzzzz, to_zz_xxx, to_zz_xxxxz, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, to_zz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yzz_xxxx[k] = -2.0 * to_zz_xxxxz[k] * tke_0 + 4.0 * to_yyzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xxxy[k] = -2.0 * to_zz_xxxyz[k] * tke_0 + 4.0 * to_yyzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xxxz[k] = to_zz_xxx[k] - 2.0 * to_zz_xxxzz[k] * tke_0 - 2.0 * to_yyzz_xxx[k] * tbe_0 + 4.0 * to_yyzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xxyy[k] = -2.0 * to_zz_xxyyz[k] * tke_0 + 4.0 * to_yyzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xxyz[k] = to_zz_xxy[k] - 2.0 * to_zz_xxyzz[k] * tke_0 - 2.0 * to_yyzz_xxy[k] * tbe_0 + 4.0 * to_yyzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xxzz[k] = 2.0 * to_zz_xxz[k] - 2.0 * to_zz_xxzzz[k] * tke_0 - 4.0 * to_yyzz_xxz[k] * tbe_0 + 4.0 * to_yyzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xyyy[k] = -2.0 * to_zz_xyyyz[k] * tke_0 + 4.0 * to_yyzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xyyz[k] = to_zz_xyy[k] - 2.0 * to_zz_xyyzz[k] * tke_0 - 2.0 * to_yyzz_xyy[k] * tbe_0 + 4.0 * to_yyzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xyzz[k] = 2.0 * to_zz_xyz[k] - 2.0 * to_zz_xyzzz[k] * tke_0 - 4.0 * to_yyzz_xyz[k] * tbe_0 + 4.0 * to_yyzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xzzz[k] = 3.0 * to_zz_xzz[k] - 2.0 * to_zz_xzzzz[k] * tke_0 - 6.0 * to_yyzz_xzz[k] * tbe_0 + 4.0 * to_yyzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_yyyy[k] = -2.0 * to_zz_yyyyz[k] * tke_0 + 4.0 * to_yyzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_yzz_yyyz[k] = to_zz_yyy[k] - 2.0 * to_zz_yyyzz[k] * tke_0 - 2.0 * to_yyzz_yyy[k] * tbe_0 + 4.0 * to_yyzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_yyzz[k] = 2.0 * to_zz_yyz[k] - 2.0 * to_zz_yyzzz[k] * tke_0 - 4.0 * to_yyzz_yyz[k] * tbe_0 + 4.0 * to_yyzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_yzzz[k] = 3.0 * to_zz_yzz[k] - 2.0 * to_zz_yzzzz[k] * tke_0 - 6.0 * to_yyzz_yzz[k] * tbe_0 + 4.0 * to_yyzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_zzzz[k] = 4.0 * to_zz_zzz[k] - 2.0 * to_zz_zzzzz[k] * tke_0 - 8.0 * to_yyzz_zzz[k] * tbe_0 + 4.0 * to_yyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 885-900 components of targeted buffer : FG

        auto to_y_z_zzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 135);

        auto to_y_z_zzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 136);

        auto to_y_z_zzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 137);

        auto to_y_z_zzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 138);

        auto to_y_z_zzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 139);

        auto to_y_z_zzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 140);

        auto to_y_z_zzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 141);

        auto to_y_z_zzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 142);

        auto to_y_z_zzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 143);

        auto to_y_z_zzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 144);

        auto to_y_z_zzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 145);

        auto to_y_z_zzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 146);

        auto to_y_z_zzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 147);

        auto to_y_z_zzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 148);

        auto to_y_z_zzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 5 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_y_z_zzz_xxxx, to_y_z_zzz_xxxy, to_y_z_zzz_xxxz, to_y_z_zzz_xxyy, to_y_z_zzz_xxyz, to_y_z_zzz_xxzz, to_y_z_zzz_xyyy, to_y_z_zzz_xyyz, to_y_z_zzz_xyzz, to_y_z_zzz_xzzz, to_y_z_zzz_yyyy, to_y_z_zzz_yyyz, to_y_z_zzz_yyzz, to_y_z_zzz_yzzz, to_y_z_zzz_zzzz, to_yzzz_xxx, to_yzzz_xxxxz, to_yzzz_xxxyz, to_yzzz_xxxzz, to_yzzz_xxy, to_yzzz_xxyyz, to_yzzz_xxyzz, to_yzzz_xxz, to_yzzz_xxzzz, to_yzzz_xyy, to_yzzz_xyyyz, to_yzzz_xyyzz, to_yzzz_xyz, to_yzzz_xyzzz, to_yzzz_xzz, to_yzzz_xzzzz, to_yzzz_yyy, to_yzzz_yyyyz, to_yzzz_yyyzz, to_yzzz_yyz, to_yzzz_yyzzz, to_yzzz_yzz, to_yzzz_yzzzz, to_yzzz_zzz, to_yzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zzz_xxxx[k] = 4.0 * to_yzzz_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xxxy[k] = 4.0 * to_yzzz_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xxxz[k] = -2.0 * to_yzzz_xxx[k] * tbe_0 + 4.0 * to_yzzz_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xxyy[k] = 4.0 * to_yzzz_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xxyz[k] = -2.0 * to_yzzz_xxy[k] * tbe_0 + 4.0 * to_yzzz_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xxzz[k] = -4.0 * to_yzzz_xxz[k] * tbe_0 + 4.0 * to_yzzz_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xyyy[k] = 4.0 * to_yzzz_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xyyz[k] = -2.0 * to_yzzz_xyy[k] * tbe_0 + 4.0 * to_yzzz_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xyzz[k] = -4.0 * to_yzzz_xyz[k] * tbe_0 + 4.0 * to_yzzz_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xzzz[k] = -6.0 * to_yzzz_xzz[k] * tbe_0 + 4.0 * to_yzzz_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_yyyy[k] = 4.0 * to_yzzz_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_zzz_yyyz[k] = -2.0 * to_yzzz_yyy[k] * tbe_0 + 4.0 * to_yzzz_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_yyzz[k] = -4.0 * to_yzzz_yyz[k] * tbe_0 + 4.0 * to_yzzz_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_yzzz[k] = -6.0 * to_yzzz_yzz[k] * tbe_0 + 4.0 * to_yzzz_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_zzzz[k] = -8.0 * to_yzzz_zzz[k] * tbe_0 + 4.0 * to_yzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 900-915 components of targeted buffer : FG

        auto to_z_x_xxx_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 0);

        auto to_z_x_xxx_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 1);

        auto to_z_x_xxx_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 2);

        auto to_z_x_xxx_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 3);

        auto to_z_x_xxx_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 4);

        auto to_z_x_xxx_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 5);

        auto to_z_x_xxx_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 6);

        auto to_z_x_xxx_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 7);

        auto to_z_x_xxx_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 8);

        auto to_z_x_xxx_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 9);

        auto to_z_x_xxx_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 10);

        auto to_z_x_xxx_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 11);

        auto to_z_x_xxx_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 12);

        auto to_z_x_xxx_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 13);

        auto to_z_x_xxx_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_xxxz_xxx, to_xxxz_xxxxx, to_xxxz_xxxxy, to_xxxz_xxxxz, to_xxxz_xxxyy, to_xxxz_xxxyz, to_xxxz_xxxzz, to_xxxz_xxy, to_xxxz_xxyyy, to_xxxz_xxyyz, to_xxxz_xxyzz, to_xxxz_xxz, to_xxxz_xxzzz, to_xxxz_xyy, to_xxxz_xyyyy, to_xxxz_xyyyz, to_xxxz_xyyzz, to_xxxz_xyz, to_xxxz_xyzzz, to_xxxz_xzz, to_xxxz_xzzzz, to_xxxz_yyy, to_xxxz_yyz, to_xxxz_yzz, to_xxxz_zzz, to_z_x_xxx_xxxx, to_z_x_xxx_xxxy, to_z_x_xxx_xxxz, to_z_x_xxx_xxyy, to_z_x_xxx_xxyz, to_z_x_xxx_xxzz, to_z_x_xxx_xyyy, to_z_x_xxx_xyyz, to_z_x_xxx_xyzz, to_z_x_xxx_xzzz, to_z_x_xxx_yyyy, to_z_x_xxx_yyyz, to_z_x_xxx_yyzz, to_z_x_xxx_yzzz, to_z_x_xxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxx_xxxx[k] = -8.0 * to_xxxz_xxx[k] * tbe_0 + 4.0 * to_xxxz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xxx_xxxy[k] = -6.0 * to_xxxz_xxy[k] * tbe_0 + 4.0 * to_xxxz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xxx_xxxz[k] = -6.0 * to_xxxz_xxz[k] * tbe_0 + 4.0 * to_xxxz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xxx_xxyy[k] = -4.0 * to_xxxz_xyy[k] * tbe_0 + 4.0 * to_xxxz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xxx_xxyz[k] = -4.0 * to_xxxz_xyz[k] * tbe_0 + 4.0 * to_xxxz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xxx_xxzz[k] = -4.0 * to_xxxz_xzz[k] * tbe_0 + 4.0 * to_xxxz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xxx_xyyy[k] = -2.0 * to_xxxz_yyy[k] * tbe_0 + 4.0 * to_xxxz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xxx_xyyz[k] = -2.0 * to_xxxz_yyz[k] * tbe_0 + 4.0 * to_xxxz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xxx_xyzz[k] = -2.0 * to_xxxz_yzz[k] * tbe_0 + 4.0 * to_xxxz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xxx_xzzz[k] = -2.0 * to_xxxz_zzz[k] * tbe_0 + 4.0 * to_xxxz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xxx_yyyy[k] = 4.0 * to_xxxz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xxx_yyyz[k] = 4.0 * to_xxxz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xxx_yyzz[k] = 4.0 * to_xxxz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xxx_yzzz[k] = 4.0 * to_xxxz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xxx_zzzz[k] = 4.0 * to_xxxz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 915-930 components of targeted buffer : FG

        auto to_z_x_xxy_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 15);

        auto to_z_x_xxy_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 16);

        auto to_z_x_xxy_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 17);

        auto to_z_x_xxy_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 18);

        auto to_z_x_xxy_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 19);

        auto to_z_x_xxy_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 20);

        auto to_z_x_xxy_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 21);

        auto to_z_x_xxy_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 22);

        auto to_z_x_xxy_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 23);

        auto to_z_x_xxy_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 24);

        auto to_z_x_xxy_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 25);

        auto to_z_x_xxy_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 26);

        auto to_z_x_xxy_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 27);

        auto to_z_x_xxy_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 28);

        auto to_z_x_xxy_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xxyz_xxx, to_xxyz_xxxxx, to_xxyz_xxxxy, to_xxyz_xxxxz, to_xxyz_xxxyy, to_xxyz_xxxyz, to_xxyz_xxxzz, to_xxyz_xxy, to_xxyz_xxyyy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xxzzz, to_xxyz_xyy, to_xxyz_xyyyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_xzzzz, to_xxyz_yyy, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_zzz, to_z_x_xxy_xxxx, to_z_x_xxy_xxxy, to_z_x_xxy_xxxz, to_z_x_xxy_xxyy, to_z_x_xxy_xxyz, to_z_x_xxy_xxzz, to_z_x_xxy_xyyy, to_z_x_xxy_xyyz, to_z_x_xxy_xyzz, to_z_x_xxy_xzzz, to_z_x_xxy_yyyy, to_z_x_xxy_yyyz, to_z_x_xxy_yyzz, to_z_x_xxy_yzzz, to_z_x_xxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxy_xxxx[k] = -8.0 * to_xxyz_xxx[k] * tbe_0 + 4.0 * to_xxyz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xxy_xxxy[k] = -6.0 * to_xxyz_xxy[k] * tbe_0 + 4.0 * to_xxyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xxy_xxxz[k] = -6.0 * to_xxyz_xxz[k] * tbe_0 + 4.0 * to_xxyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xxy_xxyy[k] = -4.0 * to_xxyz_xyy[k] * tbe_0 + 4.0 * to_xxyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xxy_xxyz[k] = -4.0 * to_xxyz_xyz[k] * tbe_0 + 4.0 * to_xxyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xxy_xxzz[k] = -4.0 * to_xxyz_xzz[k] * tbe_0 + 4.0 * to_xxyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xxy_xyyy[k] = -2.0 * to_xxyz_yyy[k] * tbe_0 + 4.0 * to_xxyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xxy_xyyz[k] = -2.0 * to_xxyz_yyz[k] * tbe_0 + 4.0 * to_xxyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xxy_xyzz[k] = -2.0 * to_xxyz_yzz[k] * tbe_0 + 4.0 * to_xxyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xxy_xzzz[k] = -2.0 * to_xxyz_zzz[k] * tbe_0 + 4.0 * to_xxyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xxy_yyyy[k] = 4.0 * to_xxyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xxy_yyyz[k] = 4.0 * to_xxyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xxy_yyzz[k] = 4.0 * to_xxyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xxy_yzzz[k] = 4.0 * to_xxyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xxy_zzzz[k] = 4.0 * to_xxyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 930-945 components of targeted buffer : FG

        auto to_z_x_xxz_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 30);

        auto to_z_x_xxz_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 31);

        auto to_z_x_xxz_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 32);

        auto to_z_x_xxz_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 33);

        auto to_z_x_xxz_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 34);

        auto to_z_x_xxz_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 35);

        auto to_z_x_xxz_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 36);

        auto to_z_x_xxz_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 37);

        auto to_z_x_xxz_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 38);

        auto to_z_x_xxz_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 39);

        auto to_z_x_xxz_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 40);

        auto to_z_x_xxz_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 41);

        auto to_z_x_xxz_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 42);

        auto to_z_x_xxz_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 43);

        auto to_z_x_xxz_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_xx_xxx, to_xx_xxxxx, to_xx_xxxxy, to_xx_xxxxz, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyz, to_xx_yzz, to_xx_zzz, to_xxzz_xxx, to_xxzz_xxxxx, to_xxzz_xxxxy, to_xxzz_xxxxz, to_xxzz_xxxyy, to_xxzz_xxxyz, to_xxzz_xxxzz, to_xxzz_xxy, to_xxzz_xxyyy, to_xxzz_xxyyz, to_xxzz_xxyzz, to_xxzz_xxz, to_xxzz_xxzzz, to_xxzz_xyy, to_xxzz_xyyyy, to_xxzz_xyyyz, to_xxzz_xyyzz, to_xxzz_xyz, to_xxzz_xyzzz, to_xxzz_xzz, to_xxzz_xzzzz, to_xxzz_yyy, to_xxzz_yyz, to_xxzz_yzz, to_xxzz_zzz, to_z_x_xxz_xxxx, to_z_x_xxz_xxxy, to_z_x_xxz_xxxz, to_z_x_xxz_xxyy, to_z_x_xxz_xxyz, to_z_x_xxz_xxzz, to_z_x_xxz_xyyy, to_z_x_xxz_xyyz, to_z_x_xxz_xyzz, to_z_x_xxz_xzzz, to_z_x_xxz_yyyy, to_z_x_xxz_yyyz, to_z_x_xxz_yyzz, to_z_x_xxz_yzzz, to_z_x_xxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxz_xxxx[k] = 4.0 * to_xx_xxx[k] - 2.0 * to_xx_xxxxx[k] * tke_0 - 8.0 * to_xxzz_xxx[k] * tbe_0 + 4.0 * to_xxzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xxz_xxxy[k] = 3.0 * to_xx_xxy[k] - 2.0 * to_xx_xxxxy[k] * tke_0 - 6.0 * to_xxzz_xxy[k] * tbe_0 + 4.0 * to_xxzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xxz_xxxz[k] = 3.0 * to_xx_xxz[k] - 2.0 * to_xx_xxxxz[k] * tke_0 - 6.0 * to_xxzz_xxz[k] * tbe_0 + 4.0 * to_xxzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xxz_xxyy[k] = 2.0 * to_xx_xyy[k] - 2.0 * to_xx_xxxyy[k] * tke_0 - 4.0 * to_xxzz_xyy[k] * tbe_0 + 4.0 * to_xxzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xxz_xxyz[k] = 2.0 * to_xx_xyz[k] - 2.0 * to_xx_xxxyz[k] * tke_0 - 4.0 * to_xxzz_xyz[k] * tbe_0 + 4.0 * to_xxzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xxz_xxzz[k] = 2.0 * to_xx_xzz[k] - 2.0 * to_xx_xxxzz[k] * tke_0 - 4.0 * to_xxzz_xzz[k] * tbe_0 + 4.0 * to_xxzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xxz_xyyy[k] = to_xx_yyy[k] - 2.0 * to_xx_xxyyy[k] * tke_0 - 2.0 * to_xxzz_yyy[k] * tbe_0 + 4.0 * to_xxzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xxz_xyyz[k] = to_xx_yyz[k] - 2.0 * to_xx_xxyyz[k] * tke_0 - 2.0 * to_xxzz_yyz[k] * tbe_0 + 4.0 * to_xxzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xxz_xyzz[k] = to_xx_yzz[k] - 2.0 * to_xx_xxyzz[k] * tke_0 - 2.0 * to_xxzz_yzz[k] * tbe_0 + 4.0 * to_xxzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xxz_xzzz[k] = to_xx_zzz[k] - 2.0 * to_xx_xxzzz[k] * tke_0 - 2.0 * to_xxzz_zzz[k] * tbe_0 + 4.0 * to_xxzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xxz_yyyy[k] = -2.0 * to_xx_xyyyy[k] * tke_0 + 4.0 * to_xxzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xxz_yyyz[k] = -2.0 * to_xx_xyyyz[k] * tke_0 + 4.0 * to_xxzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xxz_yyzz[k] = -2.0 * to_xx_xyyzz[k] * tke_0 + 4.0 * to_xxzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xxz_yzzz[k] = -2.0 * to_xx_xyzzz[k] * tke_0 + 4.0 * to_xxzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xxz_zzzz[k] = -2.0 * to_xx_xzzzz[k] * tke_0 + 4.0 * to_xxzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 945-960 components of targeted buffer : FG

        auto to_z_x_xyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 45);

        auto to_z_x_xyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 46);

        auto to_z_x_xyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 47);

        auto to_z_x_xyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 48);

        auto to_z_x_xyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 49);

        auto to_z_x_xyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 50);

        auto to_z_x_xyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 51);

        auto to_z_x_xyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 52);

        auto to_z_x_xyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 53);

        auto to_z_x_xyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 54);

        auto to_z_x_xyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 55);

        auto to_z_x_xyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 56);

        auto to_z_x_xyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 57);

        auto to_z_x_xyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 58);

        auto to_z_x_xyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xyyz_xxx, to_xyyz_xxxxx, to_xyyz_xxxxy, to_xyyz_xxxxz, to_xyyz_xxxyy, to_xyyz_xxxyz, to_xyyz_xxxzz, to_xyyz_xxy, to_xyyz_xxyyy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xxzzz, to_xyyz_xyy, to_xyyz_xyyyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_xzzzz, to_xyyz_yyy, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_zzz, to_z_x_xyy_xxxx, to_z_x_xyy_xxxy, to_z_x_xyy_xxxz, to_z_x_xyy_xxyy, to_z_x_xyy_xxyz, to_z_x_xyy_xxzz, to_z_x_xyy_xyyy, to_z_x_xyy_xyyz, to_z_x_xyy_xyzz, to_z_x_xyy_xzzz, to_z_x_xyy_yyyy, to_z_x_xyy_yyyz, to_z_x_xyy_yyzz, to_z_x_xyy_yzzz, to_z_x_xyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyy_xxxx[k] = -8.0 * to_xyyz_xxx[k] * tbe_0 + 4.0 * to_xyyz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xyy_xxxy[k] = -6.0 * to_xyyz_xxy[k] * tbe_0 + 4.0 * to_xyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xyy_xxxz[k] = -6.0 * to_xyyz_xxz[k] * tbe_0 + 4.0 * to_xyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xyy_xxyy[k] = -4.0 * to_xyyz_xyy[k] * tbe_0 + 4.0 * to_xyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xyy_xxyz[k] = -4.0 * to_xyyz_xyz[k] * tbe_0 + 4.0 * to_xyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xyy_xxzz[k] = -4.0 * to_xyyz_xzz[k] * tbe_0 + 4.0 * to_xyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xyy_xyyy[k] = -2.0 * to_xyyz_yyy[k] * tbe_0 + 4.0 * to_xyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xyy_xyyz[k] = -2.0 * to_xyyz_yyz[k] * tbe_0 + 4.0 * to_xyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xyy_xyzz[k] = -2.0 * to_xyyz_yzz[k] * tbe_0 + 4.0 * to_xyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xyy_xzzz[k] = -2.0 * to_xyyz_zzz[k] * tbe_0 + 4.0 * to_xyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xyy_yyyy[k] = 4.0 * to_xyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xyy_yyyz[k] = 4.0 * to_xyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xyy_yyzz[k] = 4.0 * to_xyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xyy_yzzz[k] = 4.0 * to_xyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xyy_zzzz[k] = 4.0 * to_xyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 960-975 components of targeted buffer : FG

        auto to_z_x_xyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 60);

        auto to_z_x_xyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 61);

        auto to_z_x_xyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 62);

        auto to_z_x_xyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 63);

        auto to_z_x_xyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 64);

        auto to_z_x_xyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 65);

        auto to_z_x_xyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 66);

        auto to_z_x_xyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 67);

        auto to_z_x_xyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 68);

        auto to_z_x_xyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 69);

        auto to_z_x_xyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 70);

        auto to_z_x_xyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 71);

        auto to_z_x_xyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 72);

        auto to_z_x_xyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 73);

        auto to_z_x_xyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_xy_xxx, to_xy_xxxxx, to_xy_xxxxy, to_xy_xxxxz, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_zzz, to_xyzz_xxx, to_xyzz_xxxxx, to_xyzz_xxxxy, to_xyzz_xxxxz, to_xyzz_xxxyy, to_xyzz_xxxyz, to_xyzz_xxxzz, to_xyzz_xxy, to_xyzz_xxyyy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xxzzz, to_xyzz_xyy, to_xyzz_xyyyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_xzzzz, to_xyzz_yyy, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_zzz, to_z_x_xyz_xxxx, to_z_x_xyz_xxxy, to_z_x_xyz_xxxz, to_z_x_xyz_xxyy, to_z_x_xyz_xxyz, to_z_x_xyz_xxzz, to_z_x_xyz_xyyy, to_z_x_xyz_xyyz, to_z_x_xyz_xyzz, to_z_x_xyz_xzzz, to_z_x_xyz_yyyy, to_z_x_xyz_yyyz, to_z_x_xyz_yyzz, to_z_x_xyz_yzzz, to_z_x_xyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyz_xxxx[k] = 4.0 * to_xy_xxx[k] - 2.0 * to_xy_xxxxx[k] * tke_0 - 8.0 * to_xyzz_xxx[k] * tbe_0 + 4.0 * to_xyzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xyz_xxxy[k] = 3.0 * to_xy_xxy[k] - 2.0 * to_xy_xxxxy[k] * tke_0 - 6.0 * to_xyzz_xxy[k] * tbe_0 + 4.0 * to_xyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xyz_xxxz[k] = 3.0 * to_xy_xxz[k] - 2.0 * to_xy_xxxxz[k] * tke_0 - 6.0 * to_xyzz_xxz[k] * tbe_0 + 4.0 * to_xyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xyz_xxyy[k] = 2.0 * to_xy_xyy[k] - 2.0 * to_xy_xxxyy[k] * tke_0 - 4.0 * to_xyzz_xyy[k] * tbe_0 + 4.0 * to_xyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xyz_xxyz[k] = 2.0 * to_xy_xyz[k] - 2.0 * to_xy_xxxyz[k] * tke_0 - 4.0 * to_xyzz_xyz[k] * tbe_0 + 4.0 * to_xyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xyz_xxzz[k] = 2.0 * to_xy_xzz[k] - 2.0 * to_xy_xxxzz[k] * tke_0 - 4.0 * to_xyzz_xzz[k] * tbe_0 + 4.0 * to_xyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xyz_xyyy[k] = to_xy_yyy[k] - 2.0 * to_xy_xxyyy[k] * tke_0 - 2.0 * to_xyzz_yyy[k] * tbe_0 + 4.0 * to_xyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xyz_xyyz[k] = to_xy_yyz[k] - 2.0 * to_xy_xxyyz[k] * tke_0 - 2.0 * to_xyzz_yyz[k] * tbe_0 + 4.0 * to_xyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xyz_xyzz[k] = to_xy_yzz[k] - 2.0 * to_xy_xxyzz[k] * tke_0 - 2.0 * to_xyzz_yzz[k] * tbe_0 + 4.0 * to_xyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xyz_xzzz[k] = to_xy_zzz[k] - 2.0 * to_xy_xxzzz[k] * tke_0 - 2.0 * to_xyzz_zzz[k] * tbe_0 + 4.0 * to_xyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xyz_yyyy[k] = -2.0 * to_xy_xyyyy[k] * tke_0 + 4.0 * to_xyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xyz_yyyz[k] = -2.0 * to_xy_xyyyz[k] * tke_0 + 4.0 * to_xyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xyz_yyzz[k] = -2.0 * to_xy_xyyzz[k] * tke_0 + 4.0 * to_xyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xyz_yzzz[k] = -2.0 * to_xy_xyzzz[k] * tke_0 + 4.0 * to_xyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xyz_zzzz[k] = -2.0 * to_xy_xzzzz[k] * tke_0 + 4.0 * to_xyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 975-990 components of targeted buffer : FG

        auto to_z_x_xzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 75);

        auto to_z_x_xzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 76);

        auto to_z_x_xzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 77);

        auto to_z_x_xzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 78);

        auto to_z_x_xzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 79);

        auto to_z_x_xzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 80);

        auto to_z_x_xzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 81);

        auto to_z_x_xzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 82);

        auto to_z_x_xzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 83);

        auto to_z_x_xzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 84);

        auto to_z_x_xzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 85);

        auto to_z_x_xzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 86);

        auto to_z_x_xzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 87);

        auto to_z_x_xzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 88);

        auto to_z_x_xzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xz_xxx, to_xz_xxxxx, to_xz_xxxxy, to_xz_xxxxz, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_zzz, to_xzzz_xxx, to_xzzz_xxxxx, to_xzzz_xxxxy, to_xzzz_xxxxz, to_xzzz_xxxyy, to_xzzz_xxxyz, to_xzzz_xxxzz, to_xzzz_xxy, to_xzzz_xxyyy, to_xzzz_xxyyz, to_xzzz_xxyzz, to_xzzz_xxz, to_xzzz_xxzzz, to_xzzz_xyy, to_xzzz_xyyyy, to_xzzz_xyyyz, to_xzzz_xyyzz, to_xzzz_xyz, to_xzzz_xyzzz, to_xzzz_xzz, to_xzzz_xzzzz, to_xzzz_yyy, to_xzzz_yyz, to_xzzz_yzz, to_xzzz_zzz, to_z_x_xzz_xxxx, to_z_x_xzz_xxxy, to_z_x_xzz_xxxz, to_z_x_xzz_xxyy, to_z_x_xzz_xxyz, to_z_x_xzz_xxzz, to_z_x_xzz_xyyy, to_z_x_xzz_xyyz, to_z_x_xzz_xyzz, to_z_x_xzz_xzzz, to_z_x_xzz_yyyy, to_z_x_xzz_yyyz, to_z_x_xzz_yyzz, to_z_x_xzz_yzzz, to_z_x_xzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xzz_xxxx[k] = 8.0 * to_xz_xxx[k] - 4.0 * to_xz_xxxxx[k] * tke_0 - 8.0 * to_xzzz_xxx[k] * tbe_0 + 4.0 * to_xzzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_xzz_xxxy[k] = 6.0 * to_xz_xxy[k] - 4.0 * to_xz_xxxxy[k] * tke_0 - 6.0 * to_xzzz_xxy[k] * tbe_0 + 4.0 * to_xzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_xzz_xxxz[k] = 6.0 * to_xz_xxz[k] - 4.0 * to_xz_xxxxz[k] * tke_0 - 6.0 * to_xzzz_xxz[k] * tbe_0 + 4.0 * to_xzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_xzz_xxyy[k] = 4.0 * to_xz_xyy[k] - 4.0 * to_xz_xxxyy[k] * tke_0 - 4.0 * to_xzzz_xyy[k] * tbe_0 + 4.0 * to_xzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_xzz_xxyz[k] = 4.0 * to_xz_xyz[k] - 4.0 * to_xz_xxxyz[k] * tke_0 - 4.0 * to_xzzz_xyz[k] * tbe_0 + 4.0 * to_xzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_xzz_xxzz[k] = 4.0 * to_xz_xzz[k] - 4.0 * to_xz_xxxzz[k] * tke_0 - 4.0 * to_xzzz_xzz[k] * tbe_0 + 4.0 * to_xzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_xzz_xyyy[k] = 2.0 * to_xz_yyy[k] - 4.0 * to_xz_xxyyy[k] * tke_0 - 2.0 * to_xzzz_yyy[k] * tbe_0 + 4.0 * to_xzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_xzz_xyyz[k] = 2.0 * to_xz_yyz[k] - 4.0 * to_xz_xxyyz[k] * tke_0 - 2.0 * to_xzzz_yyz[k] * tbe_0 + 4.0 * to_xzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_xzz_xyzz[k] = 2.0 * to_xz_yzz[k] - 4.0 * to_xz_xxyzz[k] * tke_0 - 2.0 * to_xzzz_yzz[k] * tbe_0 + 4.0 * to_xzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_xzz_xzzz[k] = 2.0 * to_xz_zzz[k] - 4.0 * to_xz_xxzzz[k] * tke_0 - 2.0 * to_xzzz_zzz[k] * tbe_0 + 4.0 * to_xzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_xzz_yyyy[k] = -4.0 * to_xz_xyyyy[k] * tke_0 + 4.0 * to_xzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_xzz_yyyz[k] = -4.0 * to_xz_xyyyz[k] * tke_0 + 4.0 * to_xzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_xzz_yyzz[k] = -4.0 * to_xz_xyyzz[k] * tke_0 + 4.0 * to_xzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_xzz_yzzz[k] = -4.0 * to_xz_xyzzz[k] * tke_0 + 4.0 * to_xzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_xzz_zzzz[k] = -4.0 * to_xz_xzzzz[k] * tke_0 + 4.0 * to_xzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 990-1005 components of targeted buffer : FG

        auto to_z_x_yyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 90);

        auto to_z_x_yyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 91);

        auto to_z_x_yyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 92);

        auto to_z_x_yyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 93);

        auto to_z_x_yyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 94);

        auto to_z_x_yyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 95);

        auto to_z_x_yyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 96);

        auto to_z_x_yyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 97);

        auto to_z_x_yyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 98);

        auto to_z_x_yyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 99);

        auto to_z_x_yyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 100);

        auto to_z_x_yyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 101);

        auto to_z_x_yyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 102);

        auto to_z_x_yyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 103);

        auto to_z_x_yyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_yyyz_xxx, to_yyyz_xxxxx, to_yyyz_xxxxy, to_yyyz_xxxxz, to_yyyz_xxxyy, to_yyyz_xxxyz, to_yyyz_xxxzz, to_yyyz_xxy, to_yyyz_xxyyy, to_yyyz_xxyyz, to_yyyz_xxyzz, to_yyyz_xxz, to_yyyz_xxzzz, to_yyyz_xyy, to_yyyz_xyyyy, to_yyyz_xyyyz, to_yyyz_xyyzz, to_yyyz_xyz, to_yyyz_xyzzz, to_yyyz_xzz, to_yyyz_xzzzz, to_yyyz_yyy, to_yyyz_yyz, to_yyyz_yzz, to_yyyz_zzz, to_z_x_yyy_xxxx, to_z_x_yyy_xxxy, to_z_x_yyy_xxxz, to_z_x_yyy_xxyy, to_z_x_yyy_xxyz, to_z_x_yyy_xxzz, to_z_x_yyy_xyyy, to_z_x_yyy_xyyz, to_z_x_yyy_xyzz, to_z_x_yyy_xzzz, to_z_x_yyy_yyyy, to_z_x_yyy_yyyz, to_z_x_yyy_yyzz, to_z_x_yyy_yzzz, to_z_x_yyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyy_xxxx[k] = -8.0 * to_yyyz_xxx[k] * tbe_0 + 4.0 * to_yyyz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_yyy_xxxy[k] = -6.0 * to_yyyz_xxy[k] * tbe_0 + 4.0 * to_yyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_yyy_xxxz[k] = -6.0 * to_yyyz_xxz[k] * tbe_0 + 4.0 * to_yyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_yyy_xxyy[k] = -4.0 * to_yyyz_xyy[k] * tbe_0 + 4.0 * to_yyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_yyy_xxyz[k] = -4.0 * to_yyyz_xyz[k] * tbe_0 + 4.0 * to_yyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_yyy_xxzz[k] = -4.0 * to_yyyz_xzz[k] * tbe_0 + 4.0 * to_yyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_yyy_xyyy[k] = -2.0 * to_yyyz_yyy[k] * tbe_0 + 4.0 * to_yyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_yyy_xyyz[k] = -2.0 * to_yyyz_yyz[k] * tbe_0 + 4.0 * to_yyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_yyy_xyzz[k] = -2.0 * to_yyyz_yzz[k] * tbe_0 + 4.0 * to_yyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_yyy_xzzz[k] = -2.0 * to_yyyz_zzz[k] * tbe_0 + 4.0 * to_yyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_yyy_yyyy[k] = 4.0 * to_yyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_yyy_yyyz[k] = 4.0 * to_yyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_yyy_yyzz[k] = 4.0 * to_yyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_yyy_yzzz[k] = 4.0 * to_yyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_yyy_zzzz[k] = 4.0 * to_yyyz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1005-1020 components of targeted buffer : FG

        auto to_z_x_yyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 105);

        auto to_z_x_yyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 106);

        auto to_z_x_yyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 107);

        auto to_z_x_yyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 108);

        auto to_z_x_yyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 109);

        auto to_z_x_yyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 110);

        auto to_z_x_yyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 111);

        auto to_z_x_yyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 112);

        auto to_z_x_yyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 113);

        auto to_z_x_yyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 114);

        auto to_z_x_yyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 115);

        auto to_z_x_yyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 116);

        auto to_z_x_yyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 117);

        auto to_z_x_yyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 118);

        auto to_z_x_yyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_yy_xxx, to_yy_xxxxx, to_yy_xxxxy, to_yy_xxxxz, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyz, to_yy_yzz, to_yy_zzz, to_yyzz_xxx, to_yyzz_xxxxx, to_yyzz_xxxxy, to_yyzz_xxxxz, to_yyzz_xxxyy, to_yyzz_xxxyz, to_yyzz_xxxzz, to_yyzz_xxy, to_yyzz_xxyyy, to_yyzz_xxyyz, to_yyzz_xxyzz, to_yyzz_xxz, to_yyzz_xxzzz, to_yyzz_xyy, to_yyzz_xyyyy, to_yyzz_xyyyz, to_yyzz_xyyzz, to_yyzz_xyz, to_yyzz_xyzzz, to_yyzz_xzz, to_yyzz_xzzzz, to_yyzz_yyy, to_yyzz_yyz, to_yyzz_yzz, to_yyzz_zzz, to_z_x_yyz_xxxx, to_z_x_yyz_xxxy, to_z_x_yyz_xxxz, to_z_x_yyz_xxyy, to_z_x_yyz_xxyz, to_z_x_yyz_xxzz, to_z_x_yyz_xyyy, to_z_x_yyz_xyyz, to_z_x_yyz_xyzz, to_z_x_yyz_xzzz, to_z_x_yyz_yyyy, to_z_x_yyz_yyyz, to_z_x_yyz_yyzz, to_z_x_yyz_yzzz, to_z_x_yyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyz_xxxx[k] = 4.0 * to_yy_xxx[k] - 2.0 * to_yy_xxxxx[k] * tke_0 - 8.0 * to_yyzz_xxx[k] * tbe_0 + 4.0 * to_yyzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_yyz_xxxy[k] = 3.0 * to_yy_xxy[k] - 2.0 * to_yy_xxxxy[k] * tke_0 - 6.0 * to_yyzz_xxy[k] * tbe_0 + 4.0 * to_yyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_yyz_xxxz[k] = 3.0 * to_yy_xxz[k] - 2.0 * to_yy_xxxxz[k] * tke_0 - 6.0 * to_yyzz_xxz[k] * tbe_0 + 4.0 * to_yyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_yyz_xxyy[k] = 2.0 * to_yy_xyy[k] - 2.0 * to_yy_xxxyy[k] * tke_0 - 4.0 * to_yyzz_xyy[k] * tbe_0 + 4.0 * to_yyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_yyz_xxyz[k] = 2.0 * to_yy_xyz[k] - 2.0 * to_yy_xxxyz[k] * tke_0 - 4.0 * to_yyzz_xyz[k] * tbe_0 + 4.0 * to_yyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_yyz_xxzz[k] = 2.0 * to_yy_xzz[k] - 2.0 * to_yy_xxxzz[k] * tke_0 - 4.0 * to_yyzz_xzz[k] * tbe_0 + 4.0 * to_yyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_yyz_xyyy[k] = to_yy_yyy[k] - 2.0 * to_yy_xxyyy[k] * tke_0 - 2.0 * to_yyzz_yyy[k] * tbe_0 + 4.0 * to_yyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_yyz_xyyz[k] = to_yy_yyz[k] - 2.0 * to_yy_xxyyz[k] * tke_0 - 2.0 * to_yyzz_yyz[k] * tbe_0 + 4.0 * to_yyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_yyz_xyzz[k] = to_yy_yzz[k] - 2.0 * to_yy_xxyzz[k] * tke_0 - 2.0 * to_yyzz_yzz[k] * tbe_0 + 4.0 * to_yyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_yyz_xzzz[k] = to_yy_zzz[k] - 2.0 * to_yy_xxzzz[k] * tke_0 - 2.0 * to_yyzz_zzz[k] * tbe_0 + 4.0 * to_yyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_yyz_yyyy[k] = -2.0 * to_yy_xyyyy[k] * tke_0 + 4.0 * to_yyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_yyz_yyyz[k] = -2.0 * to_yy_xyyyz[k] * tke_0 + 4.0 * to_yyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_yyz_yyzz[k] = -2.0 * to_yy_xyyzz[k] * tke_0 + 4.0 * to_yyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_yyz_yzzz[k] = -2.0 * to_yy_xyzzz[k] * tke_0 + 4.0 * to_yyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_yyz_zzzz[k] = -2.0 * to_yy_xzzzz[k] * tke_0 + 4.0 * to_yyzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1020-1035 components of targeted buffer : FG

        auto to_z_x_yzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 120);

        auto to_z_x_yzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 121);

        auto to_z_x_yzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 122);

        auto to_z_x_yzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 123);

        auto to_z_x_yzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 124);

        auto to_z_x_yzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 125);

        auto to_z_x_yzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 126);

        auto to_z_x_yzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 127);

        auto to_z_x_yzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 128);

        auto to_z_x_yzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 129);

        auto to_z_x_yzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 130);

        auto to_z_x_yzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 131);

        auto to_z_x_yzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 132);

        auto to_z_x_yzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 133);

        auto to_z_x_yzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_yz_xxx, to_yz_xxxxx, to_yz_xxxxy, to_yz_xxxxz, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_zzz, to_yzzz_xxx, to_yzzz_xxxxx, to_yzzz_xxxxy, to_yzzz_xxxxz, to_yzzz_xxxyy, to_yzzz_xxxyz, to_yzzz_xxxzz, to_yzzz_xxy, to_yzzz_xxyyy, to_yzzz_xxyyz, to_yzzz_xxyzz, to_yzzz_xxz, to_yzzz_xxzzz, to_yzzz_xyy, to_yzzz_xyyyy, to_yzzz_xyyyz, to_yzzz_xyyzz, to_yzzz_xyz, to_yzzz_xyzzz, to_yzzz_xzz, to_yzzz_xzzzz, to_yzzz_yyy, to_yzzz_yyz, to_yzzz_yzz, to_yzzz_zzz, to_z_x_yzz_xxxx, to_z_x_yzz_xxxy, to_z_x_yzz_xxxz, to_z_x_yzz_xxyy, to_z_x_yzz_xxyz, to_z_x_yzz_xxzz, to_z_x_yzz_xyyy, to_z_x_yzz_xyyz, to_z_x_yzz_xyzz, to_z_x_yzz_xzzz, to_z_x_yzz_yyyy, to_z_x_yzz_yyyz, to_z_x_yzz_yyzz, to_z_x_yzz_yzzz, to_z_x_yzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yzz_xxxx[k] = 8.0 * to_yz_xxx[k] - 4.0 * to_yz_xxxxx[k] * tke_0 - 8.0 * to_yzzz_xxx[k] * tbe_0 + 4.0 * to_yzzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_yzz_xxxy[k] = 6.0 * to_yz_xxy[k] - 4.0 * to_yz_xxxxy[k] * tke_0 - 6.0 * to_yzzz_xxy[k] * tbe_0 + 4.0 * to_yzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_yzz_xxxz[k] = 6.0 * to_yz_xxz[k] - 4.0 * to_yz_xxxxz[k] * tke_0 - 6.0 * to_yzzz_xxz[k] * tbe_0 + 4.0 * to_yzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_yzz_xxyy[k] = 4.0 * to_yz_xyy[k] - 4.0 * to_yz_xxxyy[k] * tke_0 - 4.0 * to_yzzz_xyy[k] * tbe_0 + 4.0 * to_yzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_yzz_xxyz[k] = 4.0 * to_yz_xyz[k] - 4.0 * to_yz_xxxyz[k] * tke_0 - 4.0 * to_yzzz_xyz[k] * tbe_0 + 4.0 * to_yzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_yzz_xxzz[k] = 4.0 * to_yz_xzz[k] - 4.0 * to_yz_xxxzz[k] * tke_0 - 4.0 * to_yzzz_xzz[k] * tbe_0 + 4.0 * to_yzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_yzz_xyyy[k] = 2.0 * to_yz_yyy[k] - 4.0 * to_yz_xxyyy[k] * tke_0 - 2.0 * to_yzzz_yyy[k] * tbe_0 + 4.0 * to_yzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_yzz_xyyz[k] = 2.0 * to_yz_yyz[k] - 4.0 * to_yz_xxyyz[k] * tke_0 - 2.0 * to_yzzz_yyz[k] * tbe_0 + 4.0 * to_yzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_yzz_xyzz[k] = 2.0 * to_yz_yzz[k] - 4.0 * to_yz_xxyzz[k] * tke_0 - 2.0 * to_yzzz_yzz[k] * tbe_0 + 4.0 * to_yzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_yzz_xzzz[k] = 2.0 * to_yz_zzz[k] - 4.0 * to_yz_xxzzz[k] * tke_0 - 2.0 * to_yzzz_zzz[k] * tbe_0 + 4.0 * to_yzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_yzz_yyyy[k] = -4.0 * to_yz_xyyyy[k] * tke_0 + 4.0 * to_yzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_yzz_yyyz[k] = -4.0 * to_yz_xyyyz[k] * tke_0 + 4.0 * to_yzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_yzz_yyzz[k] = -4.0 * to_yz_xyyzz[k] * tke_0 + 4.0 * to_yzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_yzz_yzzz[k] = -4.0 * to_yz_xyzzz[k] * tke_0 + 4.0 * to_yzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_yzz_zzzz[k] = -4.0 * to_yz_xzzzz[k] * tke_0 + 4.0 * to_yzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1035-1050 components of targeted buffer : FG

        auto to_z_x_zzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 135);

        auto to_z_x_zzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 136);

        auto to_z_x_zzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 137);

        auto to_z_x_zzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 138);

        auto to_z_x_zzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 139);

        auto to_z_x_zzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 140);

        auto to_z_x_zzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 141);

        auto to_z_x_zzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 142);

        auto to_z_x_zzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 143);

        auto to_z_x_zzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 144);

        auto to_z_x_zzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 145);

        auto to_z_x_zzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 146);

        auto to_z_x_zzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 147);

        auto to_z_x_zzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 148);

        auto to_z_x_zzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 6 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_z_x_zzz_xxxx, to_z_x_zzz_xxxy, to_z_x_zzz_xxxz, to_z_x_zzz_xxyy, to_z_x_zzz_xxyz, to_z_x_zzz_xxzz, to_z_x_zzz_xyyy, to_z_x_zzz_xyyz, to_z_x_zzz_xyzz, to_z_x_zzz_xzzz, to_z_x_zzz_yyyy, to_z_x_zzz_yyyz, to_z_x_zzz_yyzz, to_z_x_zzz_yzzz, to_z_x_zzz_zzzz, to_zz_xxx, to_zz_xxxxx, to_zz_xxxxy, to_zz_xxxxz, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyz, to_zz_yzz, to_zz_zzz, to_zzzz_xxx, to_zzzz_xxxxx, to_zzzz_xxxxy, to_zzzz_xxxxz, to_zzzz_xxxyy, to_zzzz_xxxyz, to_zzzz_xxxzz, to_zzzz_xxy, to_zzzz_xxyyy, to_zzzz_xxyyz, to_zzzz_xxyzz, to_zzzz_xxz, to_zzzz_xxzzz, to_zzzz_xyy, to_zzzz_xyyyy, to_zzzz_xyyyz, to_zzzz_xyyzz, to_zzzz_xyz, to_zzzz_xyzzz, to_zzzz_xzz, to_zzzz_xzzzz, to_zzzz_yyy, to_zzzz_yyz, to_zzzz_yzz, to_zzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zzz_xxxx[k] = 12.0 * to_zz_xxx[k] - 6.0 * to_zz_xxxxx[k] * tke_0 - 8.0 * to_zzzz_xxx[k] * tbe_0 + 4.0 * to_zzzz_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_zzz_xxxy[k] = 9.0 * to_zz_xxy[k] - 6.0 * to_zz_xxxxy[k] * tke_0 - 6.0 * to_zzzz_xxy[k] * tbe_0 + 4.0 * to_zzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_zzz_xxxz[k] = 9.0 * to_zz_xxz[k] - 6.0 * to_zz_xxxxz[k] * tke_0 - 6.0 * to_zzzz_xxz[k] * tbe_0 + 4.0 * to_zzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_zzz_xxyy[k] = 6.0 * to_zz_xyy[k] - 6.0 * to_zz_xxxyy[k] * tke_0 - 4.0 * to_zzzz_xyy[k] * tbe_0 + 4.0 * to_zzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_zzz_xxyz[k] = 6.0 * to_zz_xyz[k] - 6.0 * to_zz_xxxyz[k] * tke_0 - 4.0 * to_zzzz_xyz[k] * tbe_0 + 4.0 * to_zzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_zzz_xxzz[k] = 6.0 * to_zz_xzz[k] - 6.0 * to_zz_xxxzz[k] * tke_0 - 4.0 * to_zzzz_xzz[k] * tbe_0 + 4.0 * to_zzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_zzz_xyyy[k] = 3.0 * to_zz_yyy[k] - 6.0 * to_zz_xxyyy[k] * tke_0 - 2.0 * to_zzzz_yyy[k] * tbe_0 + 4.0 * to_zzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_zzz_xyyz[k] = 3.0 * to_zz_yyz[k] - 6.0 * to_zz_xxyyz[k] * tke_0 - 2.0 * to_zzzz_yyz[k] * tbe_0 + 4.0 * to_zzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_zzz_xyzz[k] = 3.0 * to_zz_yzz[k] - 6.0 * to_zz_xxyzz[k] * tke_0 - 2.0 * to_zzzz_yzz[k] * tbe_0 + 4.0 * to_zzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_zzz_xzzz[k] = 3.0 * to_zz_zzz[k] - 6.0 * to_zz_xxzzz[k] * tke_0 - 2.0 * to_zzzz_zzz[k] * tbe_0 + 4.0 * to_zzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_zzz_yyyy[k] = -6.0 * to_zz_xyyyy[k] * tke_0 + 4.0 * to_zzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_zzz_yyyz[k] = -6.0 * to_zz_xyyyz[k] * tke_0 + 4.0 * to_zzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_zzz_yyzz[k] = -6.0 * to_zz_xyyzz[k] * tke_0 + 4.0 * to_zzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_zzz_yzzz[k] = -6.0 * to_zz_xyzzz[k] * tke_0 + 4.0 * to_zzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_zzz_zzzz[k] = -6.0 * to_zz_xzzzz[k] * tke_0 + 4.0 * to_zzzz_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1050-1065 components of targeted buffer : FG

        auto to_z_y_xxx_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 0);

        auto to_z_y_xxx_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 1);

        auto to_z_y_xxx_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 2);

        auto to_z_y_xxx_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 3);

        auto to_z_y_xxx_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 4);

        auto to_z_y_xxx_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 5);

        auto to_z_y_xxx_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 6);

        auto to_z_y_xxx_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 7);

        auto to_z_y_xxx_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 8);

        auto to_z_y_xxx_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 9);

        auto to_z_y_xxx_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 10);

        auto to_z_y_xxx_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 11);

        auto to_z_y_xxx_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 12);

        auto to_z_y_xxx_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 13);

        auto to_z_y_xxx_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_xxxz_xxx, to_xxxz_xxxxy, to_xxxz_xxxyy, to_xxxz_xxxyz, to_xxxz_xxy, to_xxxz_xxyyy, to_xxxz_xxyyz, to_xxxz_xxyzz, to_xxxz_xxz, to_xxxz_xyy, to_xxxz_xyyyy, to_xxxz_xyyyz, to_xxxz_xyyzz, to_xxxz_xyz, to_xxxz_xyzzz, to_xxxz_xzz, to_xxxz_yyy, to_xxxz_yyyyy, to_xxxz_yyyyz, to_xxxz_yyyzz, to_xxxz_yyz, to_xxxz_yyzzz, to_xxxz_yzz, to_xxxz_yzzzz, to_xxxz_zzz, to_z_y_xxx_xxxx, to_z_y_xxx_xxxy, to_z_y_xxx_xxxz, to_z_y_xxx_xxyy, to_z_y_xxx_xxyz, to_z_y_xxx_xxzz, to_z_y_xxx_xyyy, to_z_y_xxx_xyyz, to_z_y_xxx_xyzz, to_z_y_xxx_xzzz, to_z_y_xxx_yyyy, to_z_y_xxx_yyyz, to_z_y_xxx_yyzz, to_z_y_xxx_yzzz, to_z_y_xxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxx_xxxx[k] = 4.0 * to_xxxz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xxx_xxxy[k] = -2.0 * to_xxxz_xxx[k] * tbe_0 + 4.0 * to_xxxz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xxx_xxxz[k] = 4.0 * to_xxxz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xxx_xxyy[k] = -4.0 * to_xxxz_xxy[k] * tbe_0 + 4.0 * to_xxxz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xxx_xxyz[k] = -2.0 * to_xxxz_xxz[k] * tbe_0 + 4.0 * to_xxxz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xxx_xxzz[k] = 4.0 * to_xxxz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xxx_xyyy[k] = -6.0 * to_xxxz_xyy[k] * tbe_0 + 4.0 * to_xxxz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxx_xyyz[k] = -4.0 * to_xxxz_xyz[k] * tbe_0 + 4.0 * to_xxxz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxx_xyzz[k] = -2.0 * to_xxxz_xzz[k] * tbe_0 + 4.0 * to_xxxz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxx_xzzz[k] = 4.0 * to_xxxz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxx_yyyy[k] = -8.0 * to_xxxz_yyy[k] * tbe_0 + 4.0 * to_xxxz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxx_yyyz[k] = -6.0 * to_xxxz_yyz[k] * tbe_0 + 4.0 * to_xxxz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxx_yyzz[k] = -4.0 * to_xxxz_yzz[k] * tbe_0 + 4.0 * to_xxxz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxx_yzzz[k] = -2.0 * to_xxxz_zzz[k] * tbe_0 + 4.0 * to_xxxz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxx_zzzz[k] = 4.0 * to_xxxz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1065-1080 components of targeted buffer : FG

        auto to_z_y_xxy_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 15);

        auto to_z_y_xxy_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 16);

        auto to_z_y_xxy_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 17);

        auto to_z_y_xxy_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 18);

        auto to_z_y_xxy_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 19);

        auto to_z_y_xxy_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 20);

        auto to_z_y_xxy_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 21);

        auto to_z_y_xxy_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 22);

        auto to_z_y_xxy_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 23);

        auto to_z_y_xxy_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 24);

        auto to_z_y_xxy_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 25);

        auto to_z_y_xxy_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 26);

        auto to_z_y_xxy_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 27);

        auto to_z_y_xxy_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 28);

        auto to_z_y_xxy_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xxyz_xxx, to_xxyz_xxxxy, to_xxyz_xxxyy, to_xxyz_xxxyz, to_xxyz_xxy, to_xxyz_xxyyy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xyy, to_xxyz_xyyyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_yyy, to_xxyz_yyyyy, to_xxyz_yyyyz, to_xxyz_yyyzz, to_xxyz_yyz, to_xxyz_yyzzz, to_xxyz_yzz, to_xxyz_yzzzz, to_xxyz_zzz, to_z_y_xxy_xxxx, to_z_y_xxy_xxxy, to_z_y_xxy_xxxz, to_z_y_xxy_xxyy, to_z_y_xxy_xxyz, to_z_y_xxy_xxzz, to_z_y_xxy_xyyy, to_z_y_xxy_xyyz, to_z_y_xxy_xyzz, to_z_y_xxy_xzzz, to_z_y_xxy_yyyy, to_z_y_xxy_yyyz, to_z_y_xxy_yyzz, to_z_y_xxy_yzzz, to_z_y_xxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxy_xxxx[k] = 4.0 * to_xxyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xxy_xxxy[k] = -2.0 * to_xxyz_xxx[k] * tbe_0 + 4.0 * to_xxyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xxy_xxxz[k] = 4.0 * to_xxyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xxy_xxyy[k] = -4.0 * to_xxyz_xxy[k] * tbe_0 + 4.0 * to_xxyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xxy_xxyz[k] = -2.0 * to_xxyz_xxz[k] * tbe_0 + 4.0 * to_xxyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xxy_xxzz[k] = 4.0 * to_xxyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xxy_xyyy[k] = -6.0 * to_xxyz_xyy[k] * tbe_0 + 4.0 * to_xxyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxy_xyyz[k] = -4.0 * to_xxyz_xyz[k] * tbe_0 + 4.0 * to_xxyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxy_xyzz[k] = -2.0 * to_xxyz_xzz[k] * tbe_0 + 4.0 * to_xxyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxy_xzzz[k] = 4.0 * to_xxyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxy_yyyy[k] = -8.0 * to_xxyz_yyy[k] * tbe_0 + 4.0 * to_xxyz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxy_yyyz[k] = -6.0 * to_xxyz_yyz[k] * tbe_0 + 4.0 * to_xxyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxy_yyzz[k] = -4.0 * to_xxyz_yzz[k] * tbe_0 + 4.0 * to_xxyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxy_yzzz[k] = -2.0 * to_xxyz_zzz[k] * tbe_0 + 4.0 * to_xxyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxy_zzzz[k] = 4.0 * to_xxyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1080-1095 components of targeted buffer : FG

        auto to_z_y_xxz_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 30);

        auto to_z_y_xxz_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 31);

        auto to_z_y_xxz_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 32);

        auto to_z_y_xxz_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 33);

        auto to_z_y_xxz_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 34);

        auto to_z_y_xxz_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 35);

        auto to_z_y_xxz_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 36);

        auto to_z_y_xxz_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 37);

        auto to_z_y_xxz_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 38);

        auto to_z_y_xxz_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 39);

        auto to_z_y_xxz_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 40);

        auto to_z_y_xxz_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 41);

        auto to_z_y_xxz_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 42);

        auto to_z_y_xxz_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 43);

        auto to_z_y_xxz_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_xx_xxx, to_xx_xxxxy, to_xx_xxxyy, to_xx_xxxyz, to_xx_xxy, to_xx_xxyyy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xyy, to_xx_xyyyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_yyy, to_xx_yyyyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, to_xxzz_xxx, to_xxzz_xxxxy, to_xxzz_xxxyy, to_xxzz_xxxyz, to_xxzz_xxy, to_xxzz_xxyyy, to_xxzz_xxyyz, to_xxzz_xxyzz, to_xxzz_xxz, to_xxzz_xyy, to_xxzz_xyyyy, to_xxzz_xyyyz, to_xxzz_xyyzz, to_xxzz_xyz, to_xxzz_xyzzz, to_xxzz_xzz, to_xxzz_yyy, to_xxzz_yyyyy, to_xxzz_yyyyz, to_xxzz_yyyzz, to_xxzz_yyz, to_xxzz_yyzzz, to_xxzz_yzz, to_xxzz_yzzzz, to_xxzz_zzz, to_z_y_xxz_xxxx, to_z_y_xxz_xxxy, to_z_y_xxz_xxxz, to_z_y_xxz_xxyy, to_z_y_xxz_xxyz, to_z_y_xxz_xxzz, to_z_y_xxz_xyyy, to_z_y_xxz_xyyz, to_z_y_xxz_xyzz, to_z_y_xxz_xzzz, to_z_y_xxz_yyyy, to_z_y_xxz_yyyz, to_z_y_xxz_yyzz, to_z_y_xxz_yzzz, to_z_y_xxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxz_xxxx[k] = -2.0 * to_xx_xxxxy[k] * tke_0 + 4.0 * to_xxzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xxz_xxxy[k] = to_xx_xxx[k] - 2.0 * to_xx_xxxyy[k] * tke_0 - 2.0 * to_xxzz_xxx[k] * tbe_0 + 4.0 * to_xxzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xxz_xxxz[k] = -2.0 * to_xx_xxxyz[k] * tke_0 + 4.0 * to_xxzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xxz_xxyy[k] = 2.0 * to_xx_xxy[k] - 2.0 * to_xx_xxyyy[k] * tke_0 - 4.0 * to_xxzz_xxy[k] * tbe_0 + 4.0 * to_xxzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xxz_xxyz[k] = to_xx_xxz[k] - 2.0 * to_xx_xxyyz[k] * tke_0 - 2.0 * to_xxzz_xxz[k] * tbe_0 + 4.0 * to_xxzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xxz_xxzz[k] = -2.0 * to_xx_xxyzz[k] * tke_0 + 4.0 * to_xxzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xxz_xyyy[k] = 3.0 * to_xx_xyy[k] - 2.0 * to_xx_xyyyy[k] * tke_0 - 6.0 * to_xxzz_xyy[k] * tbe_0 + 4.0 * to_xxzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxz_xyyz[k] = 2.0 * to_xx_xyz[k] - 2.0 * to_xx_xyyyz[k] * tke_0 - 4.0 * to_xxzz_xyz[k] * tbe_0 + 4.0 * to_xxzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxz_xyzz[k] = to_xx_xzz[k] - 2.0 * to_xx_xyyzz[k] * tke_0 - 2.0 * to_xxzz_xzz[k] * tbe_0 + 4.0 * to_xxzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxz_xzzz[k] = -2.0 * to_xx_xyzzz[k] * tke_0 + 4.0 * to_xxzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxz_yyyy[k] = 4.0 * to_xx_yyy[k] - 2.0 * to_xx_yyyyy[k] * tke_0 - 8.0 * to_xxzz_yyy[k] * tbe_0 + 4.0 * to_xxzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xxz_yyyz[k] = 3.0 * to_xx_yyz[k] - 2.0 * to_xx_yyyyz[k] * tke_0 - 6.0 * to_xxzz_yyz[k] * tbe_0 + 4.0 * to_xxzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xxz_yyzz[k] = 2.0 * to_xx_yzz[k] - 2.0 * to_xx_yyyzz[k] * tke_0 - 4.0 * to_xxzz_yzz[k] * tbe_0 + 4.0 * to_xxzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xxz_yzzz[k] = to_xx_zzz[k] - 2.0 * to_xx_yyzzz[k] * tke_0 - 2.0 * to_xxzz_zzz[k] * tbe_0 + 4.0 * to_xxzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xxz_zzzz[k] = -2.0 * to_xx_yzzzz[k] * tke_0 + 4.0 * to_xxzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1095-1110 components of targeted buffer : FG

        auto to_z_y_xyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 45);

        auto to_z_y_xyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 46);

        auto to_z_y_xyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 47);

        auto to_z_y_xyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 48);

        auto to_z_y_xyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 49);

        auto to_z_y_xyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 50);

        auto to_z_y_xyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 51);

        auto to_z_y_xyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 52);

        auto to_z_y_xyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 53);

        auto to_z_y_xyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 54);

        auto to_z_y_xyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 55);

        auto to_z_y_xyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 56);

        auto to_z_y_xyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 57);

        auto to_z_y_xyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 58);

        auto to_z_y_xyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xyyz_xxx, to_xyyz_xxxxy, to_xyyz_xxxyy, to_xyyz_xxxyz, to_xyyz_xxy, to_xyyz_xxyyy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xyy, to_xyyz_xyyyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_yyy, to_xyyz_yyyyy, to_xyyz_yyyyz, to_xyyz_yyyzz, to_xyyz_yyz, to_xyyz_yyzzz, to_xyyz_yzz, to_xyyz_yzzzz, to_xyyz_zzz, to_z_y_xyy_xxxx, to_z_y_xyy_xxxy, to_z_y_xyy_xxxz, to_z_y_xyy_xxyy, to_z_y_xyy_xxyz, to_z_y_xyy_xxzz, to_z_y_xyy_xyyy, to_z_y_xyy_xyyz, to_z_y_xyy_xyzz, to_z_y_xyy_xzzz, to_z_y_xyy_yyyy, to_z_y_xyy_yyyz, to_z_y_xyy_yyzz, to_z_y_xyy_yzzz, to_z_y_xyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyy_xxxx[k] = 4.0 * to_xyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xyy_xxxy[k] = -2.0 * to_xyyz_xxx[k] * tbe_0 + 4.0 * to_xyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xyy_xxxz[k] = 4.0 * to_xyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xyy_xxyy[k] = -4.0 * to_xyyz_xxy[k] * tbe_0 + 4.0 * to_xyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xyy_xxyz[k] = -2.0 * to_xyyz_xxz[k] * tbe_0 + 4.0 * to_xyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xyy_xxzz[k] = 4.0 * to_xyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xyy_xyyy[k] = -6.0 * to_xyyz_xyy[k] * tbe_0 + 4.0 * to_xyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyy_xyyz[k] = -4.0 * to_xyyz_xyz[k] * tbe_0 + 4.0 * to_xyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyy_xyzz[k] = -2.0 * to_xyyz_xzz[k] * tbe_0 + 4.0 * to_xyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyy_xzzz[k] = 4.0 * to_xyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyy_yyyy[k] = -8.0 * to_xyyz_yyy[k] * tbe_0 + 4.0 * to_xyyz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyy_yyyz[k] = -6.0 * to_xyyz_yyz[k] * tbe_0 + 4.0 * to_xyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyy_yyzz[k] = -4.0 * to_xyyz_yzz[k] * tbe_0 + 4.0 * to_xyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyy_yzzz[k] = -2.0 * to_xyyz_zzz[k] * tbe_0 + 4.0 * to_xyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyy_zzzz[k] = 4.0 * to_xyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1110-1125 components of targeted buffer : FG

        auto to_z_y_xyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 60);

        auto to_z_y_xyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 61);

        auto to_z_y_xyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 62);

        auto to_z_y_xyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 63);

        auto to_z_y_xyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 64);

        auto to_z_y_xyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 65);

        auto to_z_y_xyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 66);

        auto to_z_y_xyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 67);

        auto to_z_y_xyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 68);

        auto to_z_y_xyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 69);

        auto to_z_y_xyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 70);

        auto to_z_y_xyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 71);

        auto to_z_y_xyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 72);

        auto to_z_y_xyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 73);

        auto to_z_y_xyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_xy_xxx, to_xy_xxxxy, to_xy_xxxyy, to_xy_xxxyz, to_xy_xxy, to_xy_xxyyy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xyy, to_xy_xyyyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_yyy, to_xy_yyyyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, to_xyzz_xxx, to_xyzz_xxxxy, to_xyzz_xxxyy, to_xyzz_xxxyz, to_xyzz_xxy, to_xyzz_xxyyy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xyy, to_xyzz_xyyyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_yyy, to_xyzz_yyyyy, to_xyzz_yyyyz, to_xyzz_yyyzz, to_xyzz_yyz, to_xyzz_yyzzz, to_xyzz_yzz, to_xyzz_yzzzz, to_xyzz_zzz, to_z_y_xyz_xxxx, to_z_y_xyz_xxxy, to_z_y_xyz_xxxz, to_z_y_xyz_xxyy, to_z_y_xyz_xxyz, to_z_y_xyz_xxzz, to_z_y_xyz_xyyy, to_z_y_xyz_xyyz, to_z_y_xyz_xyzz, to_z_y_xyz_xzzz, to_z_y_xyz_yyyy, to_z_y_xyz_yyyz, to_z_y_xyz_yyzz, to_z_y_xyz_yzzz, to_z_y_xyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyz_xxxx[k] = -2.0 * to_xy_xxxxy[k] * tke_0 + 4.0 * to_xyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xyz_xxxy[k] = to_xy_xxx[k] - 2.0 * to_xy_xxxyy[k] * tke_0 - 2.0 * to_xyzz_xxx[k] * tbe_0 + 4.0 * to_xyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xyz_xxxz[k] = -2.0 * to_xy_xxxyz[k] * tke_0 + 4.0 * to_xyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xyz_xxyy[k] = 2.0 * to_xy_xxy[k] - 2.0 * to_xy_xxyyy[k] * tke_0 - 4.0 * to_xyzz_xxy[k] * tbe_0 + 4.0 * to_xyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xyz_xxyz[k] = to_xy_xxz[k] - 2.0 * to_xy_xxyyz[k] * tke_0 - 2.0 * to_xyzz_xxz[k] * tbe_0 + 4.0 * to_xyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xyz_xxzz[k] = -2.0 * to_xy_xxyzz[k] * tke_0 + 4.0 * to_xyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xyz_xyyy[k] = 3.0 * to_xy_xyy[k] - 2.0 * to_xy_xyyyy[k] * tke_0 - 6.0 * to_xyzz_xyy[k] * tbe_0 + 4.0 * to_xyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyz_xyyz[k] = 2.0 * to_xy_xyz[k] - 2.0 * to_xy_xyyyz[k] * tke_0 - 4.0 * to_xyzz_xyz[k] * tbe_0 + 4.0 * to_xyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyz_xyzz[k] = to_xy_xzz[k] - 2.0 * to_xy_xyyzz[k] * tke_0 - 2.0 * to_xyzz_xzz[k] * tbe_0 + 4.0 * to_xyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyz_xzzz[k] = -2.0 * to_xy_xyzzz[k] * tke_0 + 4.0 * to_xyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyz_yyyy[k] = 4.0 * to_xy_yyy[k] - 2.0 * to_xy_yyyyy[k] * tke_0 - 8.0 * to_xyzz_yyy[k] * tbe_0 + 4.0 * to_xyzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xyz_yyyz[k] = 3.0 * to_xy_yyz[k] - 2.0 * to_xy_yyyyz[k] * tke_0 - 6.0 * to_xyzz_yyz[k] * tbe_0 + 4.0 * to_xyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xyz_yyzz[k] = 2.0 * to_xy_yzz[k] - 2.0 * to_xy_yyyzz[k] * tke_0 - 4.0 * to_xyzz_yzz[k] * tbe_0 + 4.0 * to_xyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xyz_yzzz[k] = to_xy_zzz[k] - 2.0 * to_xy_yyzzz[k] * tke_0 - 2.0 * to_xyzz_zzz[k] * tbe_0 + 4.0 * to_xyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xyz_zzzz[k] = -2.0 * to_xy_yzzzz[k] * tke_0 + 4.0 * to_xyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1125-1140 components of targeted buffer : FG

        auto to_z_y_xzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 75);

        auto to_z_y_xzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 76);

        auto to_z_y_xzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 77);

        auto to_z_y_xzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 78);

        auto to_z_y_xzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 79);

        auto to_z_y_xzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 80);

        auto to_z_y_xzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 81);

        auto to_z_y_xzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 82);

        auto to_z_y_xzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 83);

        auto to_z_y_xzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 84);

        auto to_z_y_xzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 85);

        auto to_z_y_xzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 86);

        auto to_z_y_xzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 87);

        auto to_z_y_xzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 88);

        auto to_z_y_xzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xz_xxx, to_xz_xxxxy, to_xz_xxxyy, to_xz_xxxyz, to_xz_xxy, to_xz_xxyyy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xyy, to_xz_xyyyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_yyy, to_xz_yyyyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, to_xzzz_xxx, to_xzzz_xxxxy, to_xzzz_xxxyy, to_xzzz_xxxyz, to_xzzz_xxy, to_xzzz_xxyyy, to_xzzz_xxyyz, to_xzzz_xxyzz, to_xzzz_xxz, to_xzzz_xyy, to_xzzz_xyyyy, to_xzzz_xyyyz, to_xzzz_xyyzz, to_xzzz_xyz, to_xzzz_xyzzz, to_xzzz_xzz, to_xzzz_yyy, to_xzzz_yyyyy, to_xzzz_yyyyz, to_xzzz_yyyzz, to_xzzz_yyz, to_xzzz_yyzzz, to_xzzz_yzz, to_xzzz_yzzzz, to_xzzz_zzz, to_z_y_xzz_xxxx, to_z_y_xzz_xxxy, to_z_y_xzz_xxxz, to_z_y_xzz_xxyy, to_z_y_xzz_xxyz, to_z_y_xzz_xxzz, to_z_y_xzz_xyyy, to_z_y_xzz_xyyz, to_z_y_xzz_xyzz, to_z_y_xzz_xzzz, to_z_y_xzz_yyyy, to_z_y_xzz_yyyz, to_z_y_xzz_yyzz, to_z_y_xzz_yzzz, to_z_y_xzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xzz_xxxx[k] = -4.0 * to_xz_xxxxy[k] * tke_0 + 4.0 * to_xzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_xzz_xxxy[k] = 2.0 * to_xz_xxx[k] - 4.0 * to_xz_xxxyy[k] * tke_0 - 2.0 * to_xzzz_xxx[k] * tbe_0 + 4.0 * to_xzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_xzz_xxxz[k] = -4.0 * to_xz_xxxyz[k] * tke_0 + 4.0 * to_xzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_xzz_xxyy[k] = 4.0 * to_xz_xxy[k] - 4.0 * to_xz_xxyyy[k] * tke_0 - 4.0 * to_xzzz_xxy[k] * tbe_0 + 4.0 * to_xzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_xzz_xxyz[k] = 2.0 * to_xz_xxz[k] - 4.0 * to_xz_xxyyz[k] * tke_0 - 2.0 * to_xzzz_xxz[k] * tbe_0 + 4.0 * to_xzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_xzz_xxzz[k] = -4.0 * to_xz_xxyzz[k] * tke_0 + 4.0 * to_xzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_xzz_xyyy[k] = 6.0 * to_xz_xyy[k] - 4.0 * to_xz_xyyyy[k] * tke_0 - 6.0 * to_xzzz_xyy[k] * tbe_0 + 4.0 * to_xzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_xzz_xyyz[k] = 4.0 * to_xz_xyz[k] - 4.0 * to_xz_xyyyz[k] * tke_0 - 4.0 * to_xzzz_xyz[k] * tbe_0 + 4.0 * to_xzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_xzz_xyzz[k] = 2.0 * to_xz_xzz[k] - 4.0 * to_xz_xyyzz[k] * tke_0 - 2.0 * to_xzzz_xzz[k] * tbe_0 + 4.0 * to_xzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_xzz_xzzz[k] = -4.0 * to_xz_xyzzz[k] * tke_0 + 4.0 * to_xzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_xzz_yyyy[k] = 8.0 * to_xz_yyy[k] - 4.0 * to_xz_yyyyy[k] * tke_0 - 8.0 * to_xzzz_yyy[k] * tbe_0 + 4.0 * to_xzzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_xzz_yyyz[k] = 6.0 * to_xz_yyz[k] - 4.0 * to_xz_yyyyz[k] * tke_0 - 6.0 * to_xzzz_yyz[k] * tbe_0 + 4.0 * to_xzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_xzz_yyzz[k] = 4.0 * to_xz_yzz[k] - 4.0 * to_xz_yyyzz[k] * tke_0 - 4.0 * to_xzzz_yzz[k] * tbe_0 + 4.0 * to_xzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_xzz_yzzz[k] = 2.0 * to_xz_zzz[k] - 4.0 * to_xz_yyzzz[k] * tke_0 - 2.0 * to_xzzz_zzz[k] * tbe_0 + 4.0 * to_xzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_xzz_zzzz[k] = -4.0 * to_xz_yzzzz[k] * tke_0 + 4.0 * to_xzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1140-1155 components of targeted buffer : FG

        auto to_z_y_yyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 90);

        auto to_z_y_yyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 91);

        auto to_z_y_yyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 92);

        auto to_z_y_yyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 93);

        auto to_z_y_yyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 94);

        auto to_z_y_yyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 95);

        auto to_z_y_yyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 96);

        auto to_z_y_yyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 97);

        auto to_z_y_yyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 98);

        auto to_z_y_yyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 99);

        auto to_z_y_yyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 100);

        auto to_z_y_yyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 101);

        auto to_z_y_yyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 102);

        auto to_z_y_yyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 103);

        auto to_z_y_yyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_yyyz_xxx, to_yyyz_xxxxy, to_yyyz_xxxyy, to_yyyz_xxxyz, to_yyyz_xxy, to_yyyz_xxyyy, to_yyyz_xxyyz, to_yyyz_xxyzz, to_yyyz_xxz, to_yyyz_xyy, to_yyyz_xyyyy, to_yyyz_xyyyz, to_yyyz_xyyzz, to_yyyz_xyz, to_yyyz_xyzzz, to_yyyz_xzz, to_yyyz_yyy, to_yyyz_yyyyy, to_yyyz_yyyyz, to_yyyz_yyyzz, to_yyyz_yyz, to_yyyz_yyzzz, to_yyyz_yzz, to_yyyz_yzzzz, to_yyyz_zzz, to_z_y_yyy_xxxx, to_z_y_yyy_xxxy, to_z_y_yyy_xxxz, to_z_y_yyy_xxyy, to_z_y_yyy_xxyz, to_z_y_yyy_xxzz, to_z_y_yyy_xyyy, to_z_y_yyy_xyyz, to_z_y_yyy_xyzz, to_z_y_yyy_xzzz, to_z_y_yyy_yyyy, to_z_y_yyy_yyyz, to_z_y_yyy_yyzz, to_z_y_yyy_yzzz, to_z_y_yyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyy_xxxx[k] = 4.0 * to_yyyz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_yyy_xxxy[k] = -2.0 * to_yyyz_xxx[k] * tbe_0 + 4.0 * to_yyyz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_yyy_xxxz[k] = 4.0 * to_yyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_yyy_xxyy[k] = -4.0 * to_yyyz_xxy[k] * tbe_0 + 4.0 * to_yyyz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_yyy_xxyz[k] = -2.0 * to_yyyz_xxz[k] * tbe_0 + 4.0 * to_yyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_yyy_xxzz[k] = 4.0 * to_yyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_yyy_xyyy[k] = -6.0 * to_yyyz_xyy[k] * tbe_0 + 4.0 * to_yyyz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyy_xyyz[k] = -4.0 * to_yyyz_xyz[k] * tbe_0 + 4.0 * to_yyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyy_xyzz[k] = -2.0 * to_yyyz_xzz[k] * tbe_0 + 4.0 * to_yyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyy_xzzz[k] = 4.0 * to_yyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyy_yyyy[k] = -8.0 * to_yyyz_yyy[k] * tbe_0 + 4.0 * to_yyyz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyy_yyyz[k] = -6.0 * to_yyyz_yyz[k] * tbe_0 + 4.0 * to_yyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyy_yyzz[k] = -4.0 * to_yyyz_yzz[k] * tbe_0 + 4.0 * to_yyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyy_yzzz[k] = -2.0 * to_yyyz_zzz[k] * tbe_0 + 4.0 * to_yyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyy_zzzz[k] = 4.0 * to_yyyz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1155-1170 components of targeted buffer : FG

        auto to_z_y_yyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 105);

        auto to_z_y_yyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 106);

        auto to_z_y_yyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 107);

        auto to_z_y_yyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 108);

        auto to_z_y_yyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 109);

        auto to_z_y_yyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 110);

        auto to_z_y_yyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 111);

        auto to_z_y_yyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 112);

        auto to_z_y_yyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 113);

        auto to_z_y_yyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 114);

        auto to_z_y_yyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 115);

        auto to_z_y_yyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 116);

        auto to_z_y_yyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 117);

        auto to_z_y_yyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 118);

        auto to_z_y_yyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_yy_xxx, to_yy_xxxxy, to_yy_xxxyy, to_yy_xxxyz, to_yy_xxy, to_yy_xxyyy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xyy, to_yy_xyyyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_yyy, to_yy_yyyyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, to_yyzz_xxx, to_yyzz_xxxxy, to_yyzz_xxxyy, to_yyzz_xxxyz, to_yyzz_xxy, to_yyzz_xxyyy, to_yyzz_xxyyz, to_yyzz_xxyzz, to_yyzz_xxz, to_yyzz_xyy, to_yyzz_xyyyy, to_yyzz_xyyyz, to_yyzz_xyyzz, to_yyzz_xyz, to_yyzz_xyzzz, to_yyzz_xzz, to_yyzz_yyy, to_yyzz_yyyyy, to_yyzz_yyyyz, to_yyzz_yyyzz, to_yyzz_yyz, to_yyzz_yyzzz, to_yyzz_yzz, to_yyzz_yzzzz, to_yyzz_zzz, to_z_y_yyz_xxxx, to_z_y_yyz_xxxy, to_z_y_yyz_xxxz, to_z_y_yyz_xxyy, to_z_y_yyz_xxyz, to_z_y_yyz_xxzz, to_z_y_yyz_xyyy, to_z_y_yyz_xyyz, to_z_y_yyz_xyzz, to_z_y_yyz_xzzz, to_z_y_yyz_yyyy, to_z_y_yyz_yyyz, to_z_y_yyz_yyzz, to_z_y_yyz_yzzz, to_z_y_yyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyz_xxxx[k] = -2.0 * to_yy_xxxxy[k] * tke_0 + 4.0 * to_yyzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_yyz_xxxy[k] = to_yy_xxx[k] - 2.0 * to_yy_xxxyy[k] * tke_0 - 2.0 * to_yyzz_xxx[k] * tbe_0 + 4.0 * to_yyzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_yyz_xxxz[k] = -2.0 * to_yy_xxxyz[k] * tke_0 + 4.0 * to_yyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_yyz_xxyy[k] = 2.0 * to_yy_xxy[k] - 2.0 * to_yy_xxyyy[k] * tke_0 - 4.0 * to_yyzz_xxy[k] * tbe_0 + 4.0 * to_yyzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_yyz_xxyz[k] = to_yy_xxz[k] - 2.0 * to_yy_xxyyz[k] * tke_0 - 2.0 * to_yyzz_xxz[k] * tbe_0 + 4.0 * to_yyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_yyz_xxzz[k] = -2.0 * to_yy_xxyzz[k] * tke_0 + 4.0 * to_yyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_yyz_xyyy[k] = 3.0 * to_yy_xyy[k] - 2.0 * to_yy_xyyyy[k] * tke_0 - 6.0 * to_yyzz_xyy[k] * tbe_0 + 4.0 * to_yyzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyz_xyyz[k] = 2.0 * to_yy_xyz[k] - 2.0 * to_yy_xyyyz[k] * tke_0 - 4.0 * to_yyzz_xyz[k] * tbe_0 + 4.0 * to_yyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyz_xyzz[k] = to_yy_xzz[k] - 2.0 * to_yy_xyyzz[k] * tke_0 - 2.0 * to_yyzz_xzz[k] * tbe_0 + 4.0 * to_yyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyz_xzzz[k] = -2.0 * to_yy_xyzzz[k] * tke_0 + 4.0 * to_yyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyz_yyyy[k] = 4.0 * to_yy_yyy[k] - 2.0 * to_yy_yyyyy[k] * tke_0 - 8.0 * to_yyzz_yyy[k] * tbe_0 + 4.0 * to_yyzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_yyz_yyyz[k] = 3.0 * to_yy_yyz[k] - 2.0 * to_yy_yyyyz[k] * tke_0 - 6.0 * to_yyzz_yyz[k] * tbe_0 + 4.0 * to_yyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_yyz_yyzz[k] = 2.0 * to_yy_yzz[k] - 2.0 * to_yy_yyyzz[k] * tke_0 - 4.0 * to_yyzz_yzz[k] * tbe_0 + 4.0 * to_yyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_yyz_yzzz[k] = to_yy_zzz[k] - 2.0 * to_yy_yyzzz[k] * tke_0 - 2.0 * to_yyzz_zzz[k] * tbe_0 + 4.0 * to_yyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_yyz_zzzz[k] = -2.0 * to_yy_yzzzz[k] * tke_0 + 4.0 * to_yyzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1170-1185 components of targeted buffer : FG

        auto to_z_y_yzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 120);

        auto to_z_y_yzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 121);

        auto to_z_y_yzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 122);

        auto to_z_y_yzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 123);

        auto to_z_y_yzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 124);

        auto to_z_y_yzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 125);

        auto to_z_y_yzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 126);

        auto to_z_y_yzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 127);

        auto to_z_y_yzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 128);

        auto to_z_y_yzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 129);

        auto to_z_y_yzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 130);

        auto to_z_y_yzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 131);

        auto to_z_y_yzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 132);

        auto to_z_y_yzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 133);

        auto to_z_y_yzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_yz_xxx, to_yz_xxxxy, to_yz_xxxyy, to_yz_xxxyz, to_yz_xxy, to_yz_xxyyy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xyy, to_yz_xyyyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_yyy, to_yz_yyyyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, to_yzzz_xxx, to_yzzz_xxxxy, to_yzzz_xxxyy, to_yzzz_xxxyz, to_yzzz_xxy, to_yzzz_xxyyy, to_yzzz_xxyyz, to_yzzz_xxyzz, to_yzzz_xxz, to_yzzz_xyy, to_yzzz_xyyyy, to_yzzz_xyyyz, to_yzzz_xyyzz, to_yzzz_xyz, to_yzzz_xyzzz, to_yzzz_xzz, to_yzzz_yyy, to_yzzz_yyyyy, to_yzzz_yyyyz, to_yzzz_yyyzz, to_yzzz_yyz, to_yzzz_yyzzz, to_yzzz_yzz, to_yzzz_yzzzz, to_yzzz_zzz, to_z_y_yzz_xxxx, to_z_y_yzz_xxxy, to_z_y_yzz_xxxz, to_z_y_yzz_xxyy, to_z_y_yzz_xxyz, to_z_y_yzz_xxzz, to_z_y_yzz_xyyy, to_z_y_yzz_xyyz, to_z_y_yzz_xyzz, to_z_y_yzz_xzzz, to_z_y_yzz_yyyy, to_z_y_yzz_yyyz, to_z_y_yzz_yyzz, to_z_y_yzz_yzzz, to_z_y_yzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yzz_xxxx[k] = -4.0 * to_yz_xxxxy[k] * tke_0 + 4.0 * to_yzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_yzz_xxxy[k] = 2.0 * to_yz_xxx[k] - 4.0 * to_yz_xxxyy[k] * tke_0 - 2.0 * to_yzzz_xxx[k] * tbe_0 + 4.0 * to_yzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_yzz_xxxz[k] = -4.0 * to_yz_xxxyz[k] * tke_0 + 4.0 * to_yzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_yzz_xxyy[k] = 4.0 * to_yz_xxy[k] - 4.0 * to_yz_xxyyy[k] * tke_0 - 4.0 * to_yzzz_xxy[k] * tbe_0 + 4.0 * to_yzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_yzz_xxyz[k] = 2.0 * to_yz_xxz[k] - 4.0 * to_yz_xxyyz[k] * tke_0 - 2.0 * to_yzzz_xxz[k] * tbe_0 + 4.0 * to_yzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_yzz_xxzz[k] = -4.0 * to_yz_xxyzz[k] * tke_0 + 4.0 * to_yzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_yzz_xyyy[k] = 6.0 * to_yz_xyy[k] - 4.0 * to_yz_xyyyy[k] * tke_0 - 6.0 * to_yzzz_xyy[k] * tbe_0 + 4.0 * to_yzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_yzz_xyyz[k] = 4.0 * to_yz_xyz[k] - 4.0 * to_yz_xyyyz[k] * tke_0 - 4.0 * to_yzzz_xyz[k] * tbe_0 + 4.0 * to_yzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_yzz_xyzz[k] = 2.0 * to_yz_xzz[k] - 4.0 * to_yz_xyyzz[k] * tke_0 - 2.0 * to_yzzz_xzz[k] * tbe_0 + 4.0 * to_yzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_yzz_xzzz[k] = -4.0 * to_yz_xyzzz[k] * tke_0 + 4.0 * to_yzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_yzz_yyyy[k] = 8.0 * to_yz_yyy[k] - 4.0 * to_yz_yyyyy[k] * tke_0 - 8.0 * to_yzzz_yyy[k] * tbe_0 + 4.0 * to_yzzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_yzz_yyyz[k] = 6.0 * to_yz_yyz[k] - 4.0 * to_yz_yyyyz[k] * tke_0 - 6.0 * to_yzzz_yyz[k] * tbe_0 + 4.0 * to_yzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_yzz_yyzz[k] = 4.0 * to_yz_yzz[k] - 4.0 * to_yz_yyyzz[k] * tke_0 - 4.0 * to_yzzz_yzz[k] * tbe_0 + 4.0 * to_yzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_yzz_yzzz[k] = 2.0 * to_yz_zzz[k] - 4.0 * to_yz_yyzzz[k] * tke_0 - 2.0 * to_yzzz_zzz[k] * tbe_0 + 4.0 * to_yzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_yzz_zzzz[k] = -4.0 * to_yz_yzzzz[k] * tke_0 + 4.0 * to_yzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1185-1200 components of targeted buffer : FG

        auto to_z_y_zzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 135);

        auto to_z_y_zzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 136);

        auto to_z_y_zzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 137);

        auto to_z_y_zzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 138);

        auto to_z_y_zzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 139);

        auto to_z_y_zzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 140);

        auto to_z_y_zzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 141);

        auto to_z_y_zzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 142);

        auto to_z_y_zzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 143);

        auto to_z_y_zzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 144);

        auto to_z_y_zzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 145);

        auto to_z_y_zzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 146);

        auto to_z_y_zzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 147);

        auto to_z_y_zzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 148);

        auto to_z_y_zzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 7 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_z_y_zzz_xxxx, to_z_y_zzz_xxxy, to_z_y_zzz_xxxz, to_z_y_zzz_xxyy, to_z_y_zzz_xxyz, to_z_y_zzz_xxzz, to_z_y_zzz_xyyy, to_z_y_zzz_xyyz, to_z_y_zzz_xyzz, to_z_y_zzz_xzzz, to_z_y_zzz_yyyy, to_z_y_zzz_yyyz, to_z_y_zzz_yyzz, to_z_y_zzz_yzzz, to_z_y_zzz_zzzz, to_zz_xxx, to_zz_xxxxy, to_zz_xxxyy, to_zz_xxxyz, to_zz_xxy, to_zz_xxyyy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xyy, to_zz_xyyyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_yyy, to_zz_yyyyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, to_zzzz_xxx, to_zzzz_xxxxy, to_zzzz_xxxyy, to_zzzz_xxxyz, to_zzzz_xxy, to_zzzz_xxyyy, to_zzzz_xxyyz, to_zzzz_xxyzz, to_zzzz_xxz, to_zzzz_xyy, to_zzzz_xyyyy, to_zzzz_xyyyz, to_zzzz_xyyzz, to_zzzz_xyz, to_zzzz_xyzzz, to_zzzz_xzz, to_zzzz_yyy, to_zzzz_yyyyy, to_zzzz_yyyyz, to_zzzz_yyyzz, to_zzzz_yyz, to_zzzz_yyzzz, to_zzzz_yzz, to_zzzz_yzzzz, to_zzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zzz_xxxx[k] = -6.0 * to_zz_xxxxy[k] * tke_0 + 4.0 * to_zzzz_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_zzz_xxxy[k] = 3.0 * to_zz_xxx[k] - 6.0 * to_zz_xxxyy[k] * tke_0 - 2.0 * to_zzzz_xxx[k] * tbe_0 + 4.0 * to_zzzz_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_zzz_xxxz[k] = -6.0 * to_zz_xxxyz[k] * tke_0 + 4.0 * to_zzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_zzz_xxyy[k] = 6.0 * to_zz_xxy[k] - 6.0 * to_zz_xxyyy[k] * tke_0 - 4.0 * to_zzzz_xxy[k] * tbe_0 + 4.0 * to_zzzz_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_zzz_xxyz[k] = 3.0 * to_zz_xxz[k] - 6.0 * to_zz_xxyyz[k] * tke_0 - 2.0 * to_zzzz_xxz[k] * tbe_0 + 4.0 * to_zzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_zzz_xxzz[k] = -6.0 * to_zz_xxyzz[k] * tke_0 + 4.0 * to_zzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_zzz_xyyy[k] = 9.0 * to_zz_xyy[k] - 6.0 * to_zz_xyyyy[k] * tke_0 - 6.0 * to_zzzz_xyy[k] * tbe_0 + 4.0 * to_zzzz_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_zzz_xyyz[k] = 6.0 * to_zz_xyz[k] - 6.0 * to_zz_xyyyz[k] * tke_0 - 4.0 * to_zzzz_xyz[k] * tbe_0 + 4.0 * to_zzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_zzz_xyzz[k] = 3.0 * to_zz_xzz[k] - 6.0 * to_zz_xyyzz[k] * tke_0 - 2.0 * to_zzzz_xzz[k] * tbe_0 + 4.0 * to_zzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_zzz_xzzz[k] = -6.0 * to_zz_xyzzz[k] * tke_0 + 4.0 * to_zzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_zzz_yyyy[k] = 12.0 * to_zz_yyy[k] - 6.0 * to_zz_yyyyy[k] * tke_0 - 8.0 * to_zzzz_yyy[k] * tbe_0 + 4.0 * to_zzzz_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_zzz_yyyz[k] = 9.0 * to_zz_yyz[k] - 6.0 * to_zz_yyyyz[k] * tke_0 - 6.0 * to_zzzz_yyz[k] * tbe_0 + 4.0 * to_zzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_zzz_yyzz[k] = 6.0 * to_zz_yzz[k] - 6.0 * to_zz_yyyzz[k] * tke_0 - 4.0 * to_zzzz_yzz[k] * tbe_0 + 4.0 * to_zzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_zzz_yzzz[k] = 3.0 * to_zz_zzz[k] - 6.0 * to_zz_yyzzz[k] * tke_0 - 2.0 * to_zzzz_zzz[k] * tbe_0 + 4.0 * to_zzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_zzz_zzzz[k] = -6.0 * to_zz_yzzzz[k] * tke_0 + 4.0 * to_zzzz_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1200-1215 components of targeted buffer : FG

        auto to_z_z_xxx_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 0);

        auto to_z_z_xxx_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 1);

        auto to_z_z_xxx_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 2);

        auto to_z_z_xxx_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 3);

        auto to_z_z_xxx_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 4);

        auto to_z_z_xxx_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 5);

        auto to_z_z_xxx_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 6);

        auto to_z_z_xxx_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 7);

        auto to_z_z_xxx_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 8);

        auto to_z_z_xxx_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 9);

        auto to_z_z_xxx_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 10);

        auto to_z_z_xxx_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 11);

        auto to_z_z_xxx_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 12);

        auto to_z_z_xxx_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 13);

        auto to_z_z_xxx_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 14);

        #pragma omp simd aligned(to_xxxz_xxx, to_xxxz_xxxxz, to_xxxz_xxxyz, to_xxxz_xxxzz, to_xxxz_xxy, to_xxxz_xxyyz, to_xxxz_xxyzz, to_xxxz_xxz, to_xxxz_xxzzz, to_xxxz_xyy, to_xxxz_xyyyz, to_xxxz_xyyzz, to_xxxz_xyz, to_xxxz_xyzzz, to_xxxz_xzz, to_xxxz_xzzzz, to_xxxz_yyy, to_xxxz_yyyyz, to_xxxz_yyyzz, to_xxxz_yyz, to_xxxz_yyzzz, to_xxxz_yzz, to_xxxz_yzzzz, to_xxxz_zzz, to_xxxz_zzzzz, to_z_z_xxx_xxxx, to_z_z_xxx_xxxy, to_z_z_xxx_xxxz, to_z_z_xxx_xxyy, to_z_z_xxx_xxyz, to_z_z_xxx_xxzz, to_z_z_xxx_xyyy, to_z_z_xxx_xyyz, to_z_z_xxx_xyzz, to_z_z_xxx_xzzz, to_z_z_xxx_yyyy, to_z_z_xxx_yyyz, to_z_z_xxx_yyzz, to_z_z_xxx_yzzz, to_z_z_xxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxx_xxxx[k] = 4.0 * to_xxxz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xxxy[k] = 4.0 * to_xxxz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xxxz[k] = -2.0 * to_xxxz_xxx[k] * tbe_0 + 4.0 * to_xxxz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xxyy[k] = 4.0 * to_xxxz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xxyz[k] = -2.0 * to_xxxz_xxy[k] * tbe_0 + 4.0 * to_xxxz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xxzz[k] = -4.0 * to_xxxz_xxz[k] * tbe_0 + 4.0 * to_xxxz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xyyy[k] = 4.0 * to_xxxz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xyyz[k] = -2.0 * to_xxxz_xyy[k] * tbe_0 + 4.0 * to_xxxz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xyzz[k] = -4.0 * to_xxxz_xyz[k] * tbe_0 + 4.0 * to_xxxz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xzzz[k] = -6.0 * to_xxxz_xzz[k] * tbe_0 + 4.0 * to_xxxz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_yyyy[k] = 4.0 * to_xxxz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxx_yyyz[k] = -2.0 * to_xxxz_yyy[k] * tbe_0 + 4.0 * to_xxxz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_yyzz[k] = -4.0 * to_xxxz_yyz[k] * tbe_0 + 4.0 * to_xxxz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_yzzz[k] = -6.0 * to_xxxz_yzz[k] * tbe_0 + 4.0 * to_xxxz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_zzzz[k] = -8.0 * to_xxxz_zzz[k] * tbe_0 + 4.0 * to_xxxz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1215-1230 components of targeted buffer : FG

        auto to_z_z_xxy_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 15);

        auto to_z_z_xxy_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 16);

        auto to_z_z_xxy_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 17);

        auto to_z_z_xxy_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 18);

        auto to_z_z_xxy_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 19);

        auto to_z_z_xxy_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 20);

        auto to_z_z_xxy_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 21);

        auto to_z_z_xxy_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 22);

        auto to_z_z_xxy_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 23);

        auto to_z_z_xxy_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 24);

        auto to_z_z_xxy_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 25);

        auto to_z_z_xxy_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 26);

        auto to_z_z_xxy_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 27);

        auto to_z_z_xxy_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 28);

        auto to_z_z_xxy_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 29);

        #pragma omp simd aligned(to_xxyz_xxx, to_xxyz_xxxxz, to_xxyz_xxxyz, to_xxyz_xxxzz, to_xxyz_xxy, to_xxyz_xxyyz, to_xxyz_xxyzz, to_xxyz_xxz, to_xxyz_xxzzz, to_xxyz_xyy, to_xxyz_xyyyz, to_xxyz_xyyzz, to_xxyz_xyz, to_xxyz_xyzzz, to_xxyz_xzz, to_xxyz_xzzzz, to_xxyz_yyy, to_xxyz_yyyyz, to_xxyz_yyyzz, to_xxyz_yyz, to_xxyz_yyzzz, to_xxyz_yzz, to_xxyz_yzzzz, to_xxyz_zzz, to_xxyz_zzzzz, to_z_z_xxy_xxxx, to_z_z_xxy_xxxy, to_z_z_xxy_xxxz, to_z_z_xxy_xxyy, to_z_z_xxy_xxyz, to_z_z_xxy_xxzz, to_z_z_xxy_xyyy, to_z_z_xxy_xyyz, to_z_z_xxy_xyzz, to_z_z_xxy_xzzz, to_z_z_xxy_yyyy, to_z_z_xxy_yyyz, to_z_z_xxy_yyzz, to_z_z_xxy_yzzz, to_z_z_xxy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxy_xxxx[k] = 4.0 * to_xxyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xxxy[k] = 4.0 * to_xxyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xxxz[k] = -2.0 * to_xxyz_xxx[k] * tbe_0 + 4.0 * to_xxyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xxyy[k] = 4.0 * to_xxyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xxyz[k] = -2.0 * to_xxyz_xxy[k] * tbe_0 + 4.0 * to_xxyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xxzz[k] = -4.0 * to_xxyz_xxz[k] * tbe_0 + 4.0 * to_xxyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xyyy[k] = 4.0 * to_xxyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xyyz[k] = -2.0 * to_xxyz_xyy[k] * tbe_0 + 4.0 * to_xxyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xyzz[k] = -4.0 * to_xxyz_xyz[k] * tbe_0 + 4.0 * to_xxyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xzzz[k] = -6.0 * to_xxyz_xzz[k] * tbe_0 + 4.0 * to_xxyz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_yyyy[k] = 4.0 * to_xxyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxy_yyyz[k] = -2.0 * to_xxyz_yyy[k] * tbe_0 + 4.0 * to_xxyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_yyzz[k] = -4.0 * to_xxyz_yyz[k] * tbe_0 + 4.0 * to_xxyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_yzzz[k] = -6.0 * to_xxyz_yzz[k] * tbe_0 + 4.0 * to_xxyz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_zzzz[k] = -8.0 * to_xxyz_zzz[k] * tbe_0 + 4.0 * to_xxyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1230-1245 components of targeted buffer : FG

        auto to_z_z_xxz_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 30);

        auto to_z_z_xxz_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 31);

        auto to_z_z_xxz_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 32);

        auto to_z_z_xxz_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 33);

        auto to_z_z_xxz_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 34);

        auto to_z_z_xxz_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 35);

        auto to_z_z_xxz_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 36);

        auto to_z_z_xxz_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 37);

        auto to_z_z_xxz_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 38);

        auto to_z_z_xxz_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 39);

        auto to_z_z_xxz_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 40);

        auto to_z_z_xxz_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 41);

        auto to_z_z_xxz_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 42);

        auto to_z_z_xxz_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 43);

        auto to_z_z_xxz_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 44);

        #pragma omp simd aligned(to_xx_xxx, to_xx_xxxxz, to_xx_xxxyz, to_xx_xxxzz, to_xx_xxy, to_xx_xxyyz, to_xx_xxyzz, to_xx_xxz, to_xx_xxzzz, to_xx_xyy, to_xx_xyyyz, to_xx_xyyzz, to_xx_xyz, to_xx_xyzzz, to_xx_xzz, to_xx_xzzzz, to_xx_yyy, to_xx_yyyyz, to_xx_yyyzz, to_xx_yyz, to_xx_yyzzz, to_xx_yzz, to_xx_yzzzz, to_xx_zzz, to_xx_zzzzz, to_xxzz_xxx, to_xxzz_xxxxz, to_xxzz_xxxyz, to_xxzz_xxxzz, to_xxzz_xxy, to_xxzz_xxyyz, to_xxzz_xxyzz, to_xxzz_xxz, to_xxzz_xxzzz, to_xxzz_xyy, to_xxzz_xyyyz, to_xxzz_xyyzz, to_xxzz_xyz, to_xxzz_xyzzz, to_xxzz_xzz, to_xxzz_xzzzz, to_xxzz_yyy, to_xxzz_yyyyz, to_xxzz_yyyzz, to_xxzz_yyz, to_xxzz_yyzzz, to_xxzz_yzz, to_xxzz_yzzzz, to_xxzz_zzz, to_xxzz_zzzzz, to_z_z_xxz_xxxx, to_z_z_xxz_xxxy, to_z_z_xxz_xxxz, to_z_z_xxz_xxyy, to_z_z_xxz_xxyz, to_z_z_xxz_xxzz, to_z_z_xxz_xyyy, to_z_z_xxz_xyyz, to_z_z_xxz_xyzz, to_z_z_xxz_xzzz, to_z_z_xxz_yyyy, to_z_z_xxz_yyyz, to_z_z_xxz_yyzz, to_z_z_xxz_yzzz, to_z_z_xxz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxz_xxxx[k] = -2.0 * to_xx_xxxxz[k] * tke_0 + 4.0 * to_xxzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xxxy[k] = -2.0 * to_xx_xxxyz[k] * tke_0 + 4.0 * to_xxzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xxxz[k] = to_xx_xxx[k] - 2.0 * to_xx_xxxzz[k] * tke_0 - 2.0 * to_xxzz_xxx[k] * tbe_0 + 4.0 * to_xxzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xxyy[k] = -2.0 * to_xx_xxyyz[k] * tke_0 + 4.0 * to_xxzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xxyz[k] = to_xx_xxy[k] - 2.0 * to_xx_xxyzz[k] * tke_0 - 2.0 * to_xxzz_xxy[k] * tbe_0 + 4.0 * to_xxzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xxzz[k] = 2.0 * to_xx_xxz[k] - 2.0 * to_xx_xxzzz[k] * tke_0 - 4.0 * to_xxzz_xxz[k] * tbe_0 + 4.0 * to_xxzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xyyy[k] = -2.0 * to_xx_xyyyz[k] * tke_0 + 4.0 * to_xxzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xyyz[k] = to_xx_xyy[k] - 2.0 * to_xx_xyyzz[k] * tke_0 - 2.0 * to_xxzz_xyy[k] * tbe_0 + 4.0 * to_xxzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xyzz[k] = 2.0 * to_xx_xyz[k] - 2.0 * to_xx_xyzzz[k] * tke_0 - 4.0 * to_xxzz_xyz[k] * tbe_0 + 4.0 * to_xxzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xzzz[k] = 3.0 * to_xx_xzz[k] - 2.0 * to_xx_xzzzz[k] * tke_0 - 6.0 * to_xxzz_xzz[k] * tbe_0 + 4.0 * to_xxzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_yyyy[k] = -2.0 * to_xx_yyyyz[k] * tke_0 + 4.0 * to_xxzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xxz_yyyz[k] = to_xx_yyy[k] - 2.0 * to_xx_yyyzz[k] * tke_0 - 2.0 * to_xxzz_yyy[k] * tbe_0 + 4.0 * to_xxzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_yyzz[k] = 2.0 * to_xx_yyz[k] - 2.0 * to_xx_yyzzz[k] * tke_0 - 4.0 * to_xxzz_yyz[k] * tbe_0 + 4.0 * to_xxzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_yzzz[k] = 3.0 * to_xx_yzz[k] - 2.0 * to_xx_yzzzz[k] * tke_0 - 6.0 * to_xxzz_yzz[k] * tbe_0 + 4.0 * to_xxzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_zzzz[k] = 4.0 * to_xx_zzz[k] - 2.0 * to_xx_zzzzz[k] * tke_0 - 8.0 * to_xxzz_zzz[k] * tbe_0 + 4.0 * to_xxzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1245-1260 components of targeted buffer : FG

        auto to_z_z_xyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 45);

        auto to_z_z_xyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 46);

        auto to_z_z_xyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 47);

        auto to_z_z_xyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 48);

        auto to_z_z_xyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 49);

        auto to_z_z_xyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 50);

        auto to_z_z_xyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 51);

        auto to_z_z_xyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 52);

        auto to_z_z_xyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 53);

        auto to_z_z_xyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 54);

        auto to_z_z_xyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 55);

        auto to_z_z_xyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 56);

        auto to_z_z_xyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 57);

        auto to_z_z_xyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 58);

        auto to_z_z_xyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 59);

        #pragma omp simd aligned(to_xyyz_xxx, to_xyyz_xxxxz, to_xyyz_xxxyz, to_xyyz_xxxzz, to_xyyz_xxy, to_xyyz_xxyyz, to_xyyz_xxyzz, to_xyyz_xxz, to_xyyz_xxzzz, to_xyyz_xyy, to_xyyz_xyyyz, to_xyyz_xyyzz, to_xyyz_xyz, to_xyyz_xyzzz, to_xyyz_xzz, to_xyyz_xzzzz, to_xyyz_yyy, to_xyyz_yyyyz, to_xyyz_yyyzz, to_xyyz_yyz, to_xyyz_yyzzz, to_xyyz_yzz, to_xyyz_yzzzz, to_xyyz_zzz, to_xyyz_zzzzz, to_z_z_xyy_xxxx, to_z_z_xyy_xxxy, to_z_z_xyy_xxxz, to_z_z_xyy_xxyy, to_z_z_xyy_xxyz, to_z_z_xyy_xxzz, to_z_z_xyy_xyyy, to_z_z_xyy_xyyz, to_z_z_xyy_xyzz, to_z_z_xyy_xzzz, to_z_z_xyy_yyyy, to_z_z_xyy_yyyz, to_z_z_xyy_yyzz, to_z_z_xyy_yzzz, to_z_z_xyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyy_xxxx[k] = 4.0 * to_xyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xxxy[k] = 4.0 * to_xyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xxxz[k] = -2.0 * to_xyyz_xxx[k] * tbe_0 + 4.0 * to_xyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xxyy[k] = 4.0 * to_xyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xxyz[k] = -2.0 * to_xyyz_xxy[k] * tbe_0 + 4.0 * to_xyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xxzz[k] = -4.0 * to_xyyz_xxz[k] * tbe_0 + 4.0 * to_xyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xyyy[k] = 4.0 * to_xyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xyyz[k] = -2.0 * to_xyyz_xyy[k] * tbe_0 + 4.0 * to_xyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xyzz[k] = -4.0 * to_xyyz_xyz[k] * tbe_0 + 4.0 * to_xyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xzzz[k] = -6.0 * to_xyyz_xzz[k] * tbe_0 + 4.0 * to_xyyz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_yyyy[k] = 4.0 * to_xyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyy_yyyz[k] = -2.0 * to_xyyz_yyy[k] * tbe_0 + 4.0 * to_xyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_yyzz[k] = -4.0 * to_xyyz_yyz[k] * tbe_0 + 4.0 * to_xyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_yzzz[k] = -6.0 * to_xyyz_yzz[k] * tbe_0 + 4.0 * to_xyyz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_zzzz[k] = -8.0 * to_xyyz_zzz[k] * tbe_0 + 4.0 * to_xyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1260-1275 components of targeted buffer : FG

        auto to_z_z_xyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 60);

        auto to_z_z_xyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 61);

        auto to_z_z_xyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 62);

        auto to_z_z_xyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 63);

        auto to_z_z_xyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 64);

        auto to_z_z_xyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 65);

        auto to_z_z_xyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 66);

        auto to_z_z_xyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 67);

        auto to_z_z_xyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 68);

        auto to_z_z_xyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 69);

        auto to_z_z_xyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 70);

        auto to_z_z_xyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 71);

        auto to_z_z_xyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 72);

        auto to_z_z_xyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 73);

        auto to_z_z_xyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 74);

        #pragma omp simd aligned(to_xy_xxx, to_xy_xxxxz, to_xy_xxxyz, to_xy_xxxzz, to_xy_xxy, to_xy_xxyyz, to_xy_xxyzz, to_xy_xxz, to_xy_xxzzz, to_xy_xyy, to_xy_xyyyz, to_xy_xyyzz, to_xy_xyz, to_xy_xyzzz, to_xy_xzz, to_xy_xzzzz, to_xy_yyy, to_xy_yyyyz, to_xy_yyyzz, to_xy_yyz, to_xy_yyzzz, to_xy_yzz, to_xy_yzzzz, to_xy_zzz, to_xy_zzzzz, to_xyzz_xxx, to_xyzz_xxxxz, to_xyzz_xxxyz, to_xyzz_xxxzz, to_xyzz_xxy, to_xyzz_xxyyz, to_xyzz_xxyzz, to_xyzz_xxz, to_xyzz_xxzzz, to_xyzz_xyy, to_xyzz_xyyyz, to_xyzz_xyyzz, to_xyzz_xyz, to_xyzz_xyzzz, to_xyzz_xzz, to_xyzz_xzzzz, to_xyzz_yyy, to_xyzz_yyyyz, to_xyzz_yyyzz, to_xyzz_yyz, to_xyzz_yyzzz, to_xyzz_yzz, to_xyzz_yzzzz, to_xyzz_zzz, to_xyzz_zzzzz, to_z_z_xyz_xxxx, to_z_z_xyz_xxxy, to_z_z_xyz_xxxz, to_z_z_xyz_xxyy, to_z_z_xyz_xxyz, to_z_z_xyz_xxzz, to_z_z_xyz_xyyy, to_z_z_xyz_xyyz, to_z_z_xyz_xyzz, to_z_z_xyz_xzzz, to_z_z_xyz_yyyy, to_z_z_xyz_yyyz, to_z_z_xyz_yyzz, to_z_z_xyz_yzzz, to_z_z_xyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyz_xxxx[k] = -2.0 * to_xy_xxxxz[k] * tke_0 + 4.0 * to_xyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xxxy[k] = -2.0 * to_xy_xxxyz[k] * tke_0 + 4.0 * to_xyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xxxz[k] = to_xy_xxx[k] - 2.0 * to_xy_xxxzz[k] * tke_0 - 2.0 * to_xyzz_xxx[k] * tbe_0 + 4.0 * to_xyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xxyy[k] = -2.0 * to_xy_xxyyz[k] * tke_0 + 4.0 * to_xyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xxyz[k] = to_xy_xxy[k] - 2.0 * to_xy_xxyzz[k] * tke_0 - 2.0 * to_xyzz_xxy[k] * tbe_0 + 4.0 * to_xyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xxzz[k] = 2.0 * to_xy_xxz[k] - 2.0 * to_xy_xxzzz[k] * tke_0 - 4.0 * to_xyzz_xxz[k] * tbe_0 + 4.0 * to_xyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xyyy[k] = -2.0 * to_xy_xyyyz[k] * tke_0 + 4.0 * to_xyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xyyz[k] = to_xy_xyy[k] - 2.0 * to_xy_xyyzz[k] * tke_0 - 2.0 * to_xyzz_xyy[k] * tbe_0 + 4.0 * to_xyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xyzz[k] = 2.0 * to_xy_xyz[k] - 2.0 * to_xy_xyzzz[k] * tke_0 - 4.0 * to_xyzz_xyz[k] * tbe_0 + 4.0 * to_xyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xzzz[k] = 3.0 * to_xy_xzz[k] - 2.0 * to_xy_xzzzz[k] * tke_0 - 6.0 * to_xyzz_xzz[k] * tbe_0 + 4.0 * to_xyzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_yyyy[k] = -2.0 * to_xy_yyyyz[k] * tke_0 + 4.0 * to_xyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xyz_yyyz[k] = to_xy_yyy[k] - 2.0 * to_xy_yyyzz[k] * tke_0 - 2.0 * to_xyzz_yyy[k] * tbe_0 + 4.0 * to_xyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_yyzz[k] = 2.0 * to_xy_yyz[k] - 2.0 * to_xy_yyzzz[k] * tke_0 - 4.0 * to_xyzz_yyz[k] * tbe_0 + 4.0 * to_xyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_yzzz[k] = 3.0 * to_xy_yzz[k] - 2.0 * to_xy_yzzzz[k] * tke_0 - 6.0 * to_xyzz_yzz[k] * tbe_0 + 4.0 * to_xyzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_zzzz[k] = 4.0 * to_xy_zzz[k] - 2.0 * to_xy_zzzzz[k] * tke_0 - 8.0 * to_xyzz_zzz[k] * tbe_0 + 4.0 * to_xyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1275-1290 components of targeted buffer : FG

        auto to_z_z_xzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 75);

        auto to_z_z_xzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 76);

        auto to_z_z_xzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 77);

        auto to_z_z_xzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 78);

        auto to_z_z_xzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 79);

        auto to_z_z_xzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 80);

        auto to_z_z_xzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 81);

        auto to_z_z_xzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 82);

        auto to_z_z_xzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 83);

        auto to_z_z_xzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 84);

        auto to_z_z_xzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 85);

        auto to_z_z_xzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 86);

        auto to_z_z_xzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 87);

        auto to_z_z_xzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 88);

        auto to_z_z_xzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 89);

        #pragma omp simd aligned(to_xz_xxx, to_xz_xxxxz, to_xz_xxxyz, to_xz_xxxzz, to_xz_xxy, to_xz_xxyyz, to_xz_xxyzz, to_xz_xxz, to_xz_xxzzz, to_xz_xyy, to_xz_xyyyz, to_xz_xyyzz, to_xz_xyz, to_xz_xyzzz, to_xz_xzz, to_xz_xzzzz, to_xz_yyy, to_xz_yyyyz, to_xz_yyyzz, to_xz_yyz, to_xz_yyzzz, to_xz_yzz, to_xz_yzzzz, to_xz_zzz, to_xz_zzzzz, to_xzzz_xxx, to_xzzz_xxxxz, to_xzzz_xxxyz, to_xzzz_xxxzz, to_xzzz_xxy, to_xzzz_xxyyz, to_xzzz_xxyzz, to_xzzz_xxz, to_xzzz_xxzzz, to_xzzz_xyy, to_xzzz_xyyyz, to_xzzz_xyyzz, to_xzzz_xyz, to_xzzz_xyzzz, to_xzzz_xzz, to_xzzz_xzzzz, to_xzzz_yyy, to_xzzz_yyyyz, to_xzzz_yyyzz, to_xzzz_yyz, to_xzzz_yyzzz, to_xzzz_yzz, to_xzzz_yzzzz, to_xzzz_zzz, to_xzzz_zzzzz, to_z_z_xzz_xxxx, to_z_z_xzz_xxxy, to_z_z_xzz_xxxz, to_z_z_xzz_xxyy, to_z_z_xzz_xxyz, to_z_z_xzz_xxzz, to_z_z_xzz_xyyy, to_z_z_xzz_xyyz, to_z_z_xzz_xyzz, to_z_z_xzz_xzzz, to_z_z_xzz_yyyy, to_z_z_xzz_yyyz, to_z_z_xzz_yyzz, to_z_z_xzz_yzzz, to_z_z_xzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xzz_xxxx[k] = -4.0 * to_xz_xxxxz[k] * tke_0 + 4.0 * to_xzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xxxy[k] = -4.0 * to_xz_xxxyz[k] * tke_0 + 4.0 * to_xzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xxxz[k] = 2.0 * to_xz_xxx[k] - 4.0 * to_xz_xxxzz[k] * tke_0 - 2.0 * to_xzzz_xxx[k] * tbe_0 + 4.0 * to_xzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xxyy[k] = -4.0 * to_xz_xxyyz[k] * tke_0 + 4.0 * to_xzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xxyz[k] = 2.0 * to_xz_xxy[k] - 4.0 * to_xz_xxyzz[k] * tke_0 - 2.0 * to_xzzz_xxy[k] * tbe_0 + 4.0 * to_xzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xxzz[k] = 4.0 * to_xz_xxz[k] - 4.0 * to_xz_xxzzz[k] * tke_0 - 4.0 * to_xzzz_xxz[k] * tbe_0 + 4.0 * to_xzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xyyy[k] = -4.0 * to_xz_xyyyz[k] * tke_0 + 4.0 * to_xzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xyyz[k] = 2.0 * to_xz_xyy[k] - 4.0 * to_xz_xyyzz[k] * tke_0 - 2.0 * to_xzzz_xyy[k] * tbe_0 + 4.0 * to_xzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xyzz[k] = 4.0 * to_xz_xyz[k] - 4.0 * to_xz_xyzzz[k] * tke_0 - 4.0 * to_xzzz_xyz[k] * tbe_0 + 4.0 * to_xzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xzzz[k] = 6.0 * to_xz_xzz[k] - 4.0 * to_xz_xzzzz[k] * tke_0 - 6.0 * to_xzzz_xzz[k] * tbe_0 + 4.0 * to_xzzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_yyyy[k] = -4.0 * to_xz_yyyyz[k] * tke_0 + 4.0 * to_xzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_xzz_yyyz[k] = 2.0 * to_xz_yyy[k] - 4.0 * to_xz_yyyzz[k] * tke_0 - 2.0 * to_xzzz_yyy[k] * tbe_0 + 4.0 * to_xzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_yyzz[k] = 4.0 * to_xz_yyz[k] - 4.0 * to_xz_yyzzz[k] * tke_0 - 4.0 * to_xzzz_yyz[k] * tbe_0 + 4.0 * to_xzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_yzzz[k] = 6.0 * to_xz_yzz[k] - 4.0 * to_xz_yzzzz[k] * tke_0 - 6.0 * to_xzzz_yzz[k] * tbe_0 + 4.0 * to_xzzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_zzzz[k] = 8.0 * to_xz_zzz[k] - 4.0 * to_xz_zzzzz[k] * tke_0 - 8.0 * to_xzzz_zzz[k] * tbe_0 + 4.0 * to_xzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1290-1305 components of targeted buffer : FG

        auto to_z_z_yyy_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 90);

        auto to_z_z_yyy_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 91);

        auto to_z_z_yyy_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 92);

        auto to_z_z_yyy_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 93);

        auto to_z_z_yyy_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 94);

        auto to_z_z_yyy_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 95);

        auto to_z_z_yyy_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 96);

        auto to_z_z_yyy_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 97);

        auto to_z_z_yyy_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 98);

        auto to_z_z_yyy_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 99);

        auto to_z_z_yyy_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 100);

        auto to_z_z_yyy_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 101);

        auto to_z_z_yyy_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 102);

        auto to_z_z_yyy_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 103);

        auto to_z_z_yyy_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 104);

        #pragma omp simd aligned(to_yyyz_xxx, to_yyyz_xxxxz, to_yyyz_xxxyz, to_yyyz_xxxzz, to_yyyz_xxy, to_yyyz_xxyyz, to_yyyz_xxyzz, to_yyyz_xxz, to_yyyz_xxzzz, to_yyyz_xyy, to_yyyz_xyyyz, to_yyyz_xyyzz, to_yyyz_xyz, to_yyyz_xyzzz, to_yyyz_xzz, to_yyyz_xzzzz, to_yyyz_yyy, to_yyyz_yyyyz, to_yyyz_yyyzz, to_yyyz_yyz, to_yyyz_yyzzz, to_yyyz_yzz, to_yyyz_yzzzz, to_yyyz_zzz, to_yyyz_zzzzz, to_z_z_yyy_xxxx, to_z_z_yyy_xxxy, to_z_z_yyy_xxxz, to_z_z_yyy_xxyy, to_z_z_yyy_xxyz, to_z_z_yyy_xxzz, to_z_z_yyy_xyyy, to_z_z_yyy_xyyz, to_z_z_yyy_xyzz, to_z_z_yyy_xzzz, to_z_z_yyy_yyyy, to_z_z_yyy_yyyz, to_z_z_yyy_yyzz, to_z_z_yyy_yzzz, to_z_z_yyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyy_xxxx[k] = 4.0 * to_yyyz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xxxy[k] = 4.0 * to_yyyz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xxxz[k] = -2.0 * to_yyyz_xxx[k] * tbe_0 + 4.0 * to_yyyz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xxyy[k] = 4.0 * to_yyyz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xxyz[k] = -2.0 * to_yyyz_xxy[k] * tbe_0 + 4.0 * to_yyyz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xxzz[k] = -4.0 * to_yyyz_xxz[k] * tbe_0 + 4.0 * to_yyyz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xyyy[k] = 4.0 * to_yyyz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xyyz[k] = -2.0 * to_yyyz_xyy[k] * tbe_0 + 4.0 * to_yyyz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xyzz[k] = -4.0 * to_yyyz_xyz[k] * tbe_0 + 4.0 * to_yyyz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xzzz[k] = -6.0 * to_yyyz_xzz[k] * tbe_0 + 4.0 * to_yyyz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_yyyy[k] = 4.0 * to_yyyz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyy_yyyz[k] = -2.0 * to_yyyz_yyy[k] * tbe_0 + 4.0 * to_yyyz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_yyzz[k] = -4.0 * to_yyyz_yyz[k] * tbe_0 + 4.0 * to_yyyz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_yzzz[k] = -6.0 * to_yyyz_yzz[k] * tbe_0 + 4.0 * to_yyyz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_zzzz[k] = -8.0 * to_yyyz_zzz[k] * tbe_0 + 4.0 * to_yyyz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1305-1320 components of targeted buffer : FG

        auto to_z_z_yyz_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 105);

        auto to_z_z_yyz_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 106);

        auto to_z_z_yyz_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 107);

        auto to_z_z_yyz_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 108);

        auto to_z_z_yyz_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 109);

        auto to_z_z_yyz_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 110);

        auto to_z_z_yyz_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 111);

        auto to_z_z_yyz_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 112);

        auto to_z_z_yyz_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 113);

        auto to_z_z_yyz_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 114);

        auto to_z_z_yyz_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 115);

        auto to_z_z_yyz_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 116);

        auto to_z_z_yyz_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 117);

        auto to_z_z_yyz_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 118);

        auto to_z_z_yyz_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 119);

        #pragma omp simd aligned(to_yy_xxx, to_yy_xxxxz, to_yy_xxxyz, to_yy_xxxzz, to_yy_xxy, to_yy_xxyyz, to_yy_xxyzz, to_yy_xxz, to_yy_xxzzz, to_yy_xyy, to_yy_xyyyz, to_yy_xyyzz, to_yy_xyz, to_yy_xyzzz, to_yy_xzz, to_yy_xzzzz, to_yy_yyy, to_yy_yyyyz, to_yy_yyyzz, to_yy_yyz, to_yy_yyzzz, to_yy_yzz, to_yy_yzzzz, to_yy_zzz, to_yy_zzzzz, to_yyzz_xxx, to_yyzz_xxxxz, to_yyzz_xxxyz, to_yyzz_xxxzz, to_yyzz_xxy, to_yyzz_xxyyz, to_yyzz_xxyzz, to_yyzz_xxz, to_yyzz_xxzzz, to_yyzz_xyy, to_yyzz_xyyyz, to_yyzz_xyyzz, to_yyzz_xyz, to_yyzz_xyzzz, to_yyzz_xzz, to_yyzz_xzzzz, to_yyzz_yyy, to_yyzz_yyyyz, to_yyzz_yyyzz, to_yyzz_yyz, to_yyzz_yyzzz, to_yyzz_yzz, to_yyzz_yzzzz, to_yyzz_zzz, to_yyzz_zzzzz, to_z_z_yyz_xxxx, to_z_z_yyz_xxxy, to_z_z_yyz_xxxz, to_z_z_yyz_xxyy, to_z_z_yyz_xxyz, to_z_z_yyz_xxzz, to_z_z_yyz_xyyy, to_z_z_yyz_xyyz, to_z_z_yyz_xyzz, to_z_z_yyz_xzzz, to_z_z_yyz_yyyy, to_z_z_yyz_yyyz, to_z_z_yyz_yyzz, to_z_z_yyz_yzzz, to_z_z_yyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyz_xxxx[k] = -2.0 * to_yy_xxxxz[k] * tke_0 + 4.0 * to_yyzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xxxy[k] = -2.0 * to_yy_xxxyz[k] * tke_0 + 4.0 * to_yyzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xxxz[k] = to_yy_xxx[k] - 2.0 * to_yy_xxxzz[k] * tke_0 - 2.0 * to_yyzz_xxx[k] * tbe_0 + 4.0 * to_yyzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xxyy[k] = -2.0 * to_yy_xxyyz[k] * tke_0 + 4.0 * to_yyzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xxyz[k] = to_yy_xxy[k] - 2.0 * to_yy_xxyzz[k] * tke_0 - 2.0 * to_yyzz_xxy[k] * tbe_0 + 4.0 * to_yyzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xxzz[k] = 2.0 * to_yy_xxz[k] - 2.0 * to_yy_xxzzz[k] * tke_0 - 4.0 * to_yyzz_xxz[k] * tbe_0 + 4.0 * to_yyzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xyyy[k] = -2.0 * to_yy_xyyyz[k] * tke_0 + 4.0 * to_yyzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xyyz[k] = to_yy_xyy[k] - 2.0 * to_yy_xyyzz[k] * tke_0 - 2.0 * to_yyzz_xyy[k] * tbe_0 + 4.0 * to_yyzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xyzz[k] = 2.0 * to_yy_xyz[k] - 2.0 * to_yy_xyzzz[k] * tke_0 - 4.0 * to_yyzz_xyz[k] * tbe_0 + 4.0 * to_yyzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xzzz[k] = 3.0 * to_yy_xzz[k] - 2.0 * to_yy_xzzzz[k] * tke_0 - 6.0 * to_yyzz_xzz[k] * tbe_0 + 4.0 * to_yyzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_yyyy[k] = -2.0 * to_yy_yyyyz[k] * tke_0 + 4.0 * to_yyzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_yyz_yyyz[k] = to_yy_yyy[k] - 2.0 * to_yy_yyyzz[k] * tke_0 - 2.0 * to_yyzz_yyy[k] * tbe_0 + 4.0 * to_yyzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_yyzz[k] = 2.0 * to_yy_yyz[k] - 2.0 * to_yy_yyzzz[k] * tke_0 - 4.0 * to_yyzz_yyz[k] * tbe_0 + 4.0 * to_yyzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_yzzz[k] = 3.0 * to_yy_yzz[k] - 2.0 * to_yy_yzzzz[k] * tke_0 - 6.0 * to_yyzz_yzz[k] * tbe_0 + 4.0 * to_yyzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_zzzz[k] = 4.0 * to_yy_zzz[k] - 2.0 * to_yy_zzzzz[k] * tke_0 - 8.0 * to_yyzz_zzz[k] * tbe_0 + 4.0 * to_yyzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1320-1335 components of targeted buffer : FG

        auto to_z_z_yzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 120);

        auto to_z_z_yzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 121);

        auto to_z_z_yzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 122);

        auto to_z_z_yzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 123);

        auto to_z_z_yzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 124);

        auto to_z_z_yzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 125);

        auto to_z_z_yzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 126);

        auto to_z_z_yzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 127);

        auto to_z_z_yzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 128);

        auto to_z_z_yzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 129);

        auto to_z_z_yzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 130);

        auto to_z_z_yzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 131);

        auto to_z_z_yzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 132);

        auto to_z_z_yzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 133);

        auto to_z_z_yzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 134);

        #pragma omp simd aligned(to_yz_xxx, to_yz_xxxxz, to_yz_xxxyz, to_yz_xxxzz, to_yz_xxy, to_yz_xxyyz, to_yz_xxyzz, to_yz_xxz, to_yz_xxzzz, to_yz_xyy, to_yz_xyyyz, to_yz_xyyzz, to_yz_xyz, to_yz_xyzzz, to_yz_xzz, to_yz_xzzzz, to_yz_yyy, to_yz_yyyyz, to_yz_yyyzz, to_yz_yyz, to_yz_yyzzz, to_yz_yzz, to_yz_yzzzz, to_yz_zzz, to_yz_zzzzz, to_yzzz_xxx, to_yzzz_xxxxz, to_yzzz_xxxyz, to_yzzz_xxxzz, to_yzzz_xxy, to_yzzz_xxyyz, to_yzzz_xxyzz, to_yzzz_xxz, to_yzzz_xxzzz, to_yzzz_xyy, to_yzzz_xyyyz, to_yzzz_xyyzz, to_yzzz_xyz, to_yzzz_xyzzz, to_yzzz_xzz, to_yzzz_xzzzz, to_yzzz_yyy, to_yzzz_yyyyz, to_yzzz_yyyzz, to_yzzz_yyz, to_yzzz_yyzzz, to_yzzz_yzz, to_yzzz_yzzzz, to_yzzz_zzz, to_yzzz_zzzzz, to_z_z_yzz_xxxx, to_z_z_yzz_xxxy, to_z_z_yzz_xxxz, to_z_z_yzz_xxyy, to_z_z_yzz_xxyz, to_z_z_yzz_xxzz, to_z_z_yzz_xyyy, to_z_z_yzz_xyyz, to_z_z_yzz_xyzz, to_z_z_yzz_xzzz, to_z_z_yzz_yyyy, to_z_z_yzz_yyyz, to_z_z_yzz_yyzz, to_z_z_yzz_yzzz, to_z_z_yzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yzz_xxxx[k] = -4.0 * to_yz_xxxxz[k] * tke_0 + 4.0 * to_yzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xxxy[k] = -4.0 * to_yz_xxxyz[k] * tke_0 + 4.0 * to_yzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xxxz[k] = 2.0 * to_yz_xxx[k] - 4.0 * to_yz_xxxzz[k] * tke_0 - 2.0 * to_yzzz_xxx[k] * tbe_0 + 4.0 * to_yzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xxyy[k] = -4.0 * to_yz_xxyyz[k] * tke_0 + 4.0 * to_yzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xxyz[k] = 2.0 * to_yz_xxy[k] - 4.0 * to_yz_xxyzz[k] * tke_0 - 2.0 * to_yzzz_xxy[k] * tbe_0 + 4.0 * to_yzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xxzz[k] = 4.0 * to_yz_xxz[k] - 4.0 * to_yz_xxzzz[k] * tke_0 - 4.0 * to_yzzz_xxz[k] * tbe_0 + 4.0 * to_yzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xyyy[k] = -4.0 * to_yz_xyyyz[k] * tke_0 + 4.0 * to_yzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xyyz[k] = 2.0 * to_yz_xyy[k] - 4.0 * to_yz_xyyzz[k] * tke_0 - 2.0 * to_yzzz_xyy[k] * tbe_0 + 4.0 * to_yzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xyzz[k] = 4.0 * to_yz_xyz[k] - 4.0 * to_yz_xyzzz[k] * tke_0 - 4.0 * to_yzzz_xyz[k] * tbe_0 + 4.0 * to_yzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xzzz[k] = 6.0 * to_yz_xzz[k] - 4.0 * to_yz_xzzzz[k] * tke_0 - 6.0 * to_yzzz_xzz[k] * tbe_0 + 4.0 * to_yzzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_yyyy[k] = -4.0 * to_yz_yyyyz[k] * tke_0 + 4.0 * to_yzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_yzz_yyyz[k] = 2.0 * to_yz_yyy[k] - 4.0 * to_yz_yyyzz[k] * tke_0 - 2.0 * to_yzzz_yyy[k] * tbe_0 + 4.0 * to_yzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_yyzz[k] = 4.0 * to_yz_yyz[k] - 4.0 * to_yz_yyzzz[k] * tke_0 - 4.0 * to_yzzz_yyz[k] * tbe_0 + 4.0 * to_yzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_yzzz[k] = 6.0 * to_yz_yzz[k] - 4.0 * to_yz_yzzzz[k] * tke_0 - 6.0 * to_yzzz_yzz[k] * tbe_0 + 4.0 * to_yzzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_zzzz[k] = 8.0 * to_yz_zzz[k] - 4.0 * to_yz_zzzzz[k] * tke_0 - 8.0 * to_yzzz_zzz[k] * tbe_0 + 4.0 * to_yzzz_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 1335-1350 components of targeted buffer : FG

        auto to_z_z_zzz_xxxx = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 135);

        auto to_z_z_zzz_xxxy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 136);

        auto to_z_z_zzz_xxxz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 137);

        auto to_z_z_zzz_xxyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 138);

        auto to_z_z_zzz_xxyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 139);

        auto to_z_z_zzz_xxzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 140);

        auto to_z_z_zzz_xyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 141);

        auto to_z_z_zzz_xyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 142);

        auto to_z_z_zzz_xyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 143);

        auto to_z_z_zzz_xzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 144);

        auto to_z_z_zzz_yyyy = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 145);

        auto to_z_z_zzz_yyyz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 146);

        auto to_z_z_zzz_yyzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 147);

        auto to_z_z_zzz_yzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 148);

        auto to_z_z_zzz_zzzz = pbuffer.data(idx_op_geom_101_fg + 8 * op_comps * 150 + i * 150 + 149);

        #pragma omp simd aligned(to_z_z_zzz_xxxx, to_z_z_zzz_xxxy, to_z_z_zzz_xxxz, to_z_z_zzz_xxyy, to_z_z_zzz_xxyz, to_z_z_zzz_xxzz, to_z_z_zzz_xyyy, to_z_z_zzz_xyyz, to_z_z_zzz_xyzz, to_z_z_zzz_xzzz, to_z_z_zzz_yyyy, to_z_z_zzz_yyyz, to_z_z_zzz_yyzz, to_z_z_zzz_yzzz, to_z_z_zzz_zzzz, to_zz_xxx, to_zz_xxxxz, to_zz_xxxyz, to_zz_xxxzz, to_zz_xxy, to_zz_xxyyz, to_zz_xxyzz, to_zz_xxz, to_zz_xxzzz, to_zz_xyy, to_zz_xyyyz, to_zz_xyyzz, to_zz_xyz, to_zz_xyzzz, to_zz_xzz, to_zz_xzzzz, to_zz_yyy, to_zz_yyyyz, to_zz_yyyzz, to_zz_yyz, to_zz_yyzzz, to_zz_yzz, to_zz_yzzzz, to_zz_zzz, to_zz_zzzzz, to_zzzz_xxx, to_zzzz_xxxxz, to_zzzz_xxxyz, to_zzzz_xxxzz, to_zzzz_xxy, to_zzzz_xxyyz, to_zzzz_xxyzz, to_zzzz_xxz, to_zzzz_xxzzz, to_zzzz_xyy, to_zzzz_xyyyz, to_zzzz_xyyzz, to_zzzz_xyz, to_zzzz_xyzzz, to_zzzz_xzz, to_zzzz_xzzzz, to_zzzz_yyy, to_zzzz_yyyyz, to_zzzz_yyyzz, to_zzzz_yyz, to_zzzz_yyzzz, to_zzzz_yzz, to_zzzz_yzzzz, to_zzzz_zzz, to_zzzz_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zzz_xxxx[k] = -6.0 * to_zz_xxxxz[k] * tke_0 + 4.0 * to_zzzz_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xxxy[k] = -6.0 * to_zz_xxxyz[k] * tke_0 + 4.0 * to_zzzz_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xxxz[k] = 3.0 * to_zz_xxx[k] - 6.0 * to_zz_xxxzz[k] * tke_0 - 2.0 * to_zzzz_xxx[k] * tbe_0 + 4.0 * to_zzzz_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xxyy[k] = -6.0 * to_zz_xxyyz[k] * tke_0 + 4.0 * to_zzzz_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xxyz[k] = 3.0 * to_zz_xxy[k] - 6.0 * to_zz_xxyzz[k] * tke_0 - 2.0 * to_zzzz_xxy[k] * tbe_0 + 4.0 * to_zzzz_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xxzz[k] = 6.0 * to_zz_xxz[k] - 6.0 * to_zz_xxzzz[k] * tke_0 - 4.0 * to_zzzz_xxz[k] * tbe_0 + 4.0 * to_zzzz_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xyyy[k] = -6.0 * to_zz_xyyyz[k] * tke_0 + 4.0 * to_zzzz_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xyyz[k] = 3.0 * to_zz_xyy[k] - 6.0 * to_zz_xyyzz[k] * tke_0 - 2.0 * to_zzzz_xyy[k] * tbe_0 + 4.0 * to_zzzz_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xyzz[k] = 6.0 * to_zz_xyz[k] - 6.0 * to_zz_xyzzz[k] * tke_0 - 4.0 * to_zzzz_xyz[k] * tbe_0 + 4.0 * to_zzzz_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xzzz[k] = 9.0 * to_zz_xzz[k] - 6.0 * to_zz_xzzzz[k] * tke_0 - 6.0 * to_zzzz_xzz[k] * tbe_0 + 4.0 * to_zzzz_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_yyyy[k] = -6.0 * to_zz_yyyyz[k] * tke_0 + 4.0 * to_zzzz_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_zzz_yyyz[k] = 3.0 * to_zz_yyy[k] - 6.0 * to_zz_yyyzz[k] * tke_0 - 2.0 * to_zzzz_yyy[k] * tbe_0 + 4.0 * to_zzzz_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_yyzz[k] = 6.0 * to_zz_yyz[k] - 6.0 * to_zz_yyzzz[k] * tke_0 - 4.0 * to_zzzz_yyz[k] * tbe_0 + 4.0 * to_zzzz_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_yzzz[k] = 9.0 * to_zz_yzz[k] - 6.0 * to_zz_yzzzz[k] * tke_0 - 6.0 * to_zzzz_yzz[k] * tbe_0 + 4.0 * to_zzzz_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_zzzz[k] = 12.0 * to_zz_zzz[k] - 6.0 * to_zz_zzzzz[k] * tke_0 - 8.0 * to_zzzz_zzz[k] * tbe_0 + 4.0 * to_zzzz_zzzzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

