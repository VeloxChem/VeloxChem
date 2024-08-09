#include "GeometricalDerivatives1X1ForFD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_fd(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_fd,
                        const size_t idx_op_dp,
                        const size_t idx_op_df,
                        const size_t idx_op_gp,
                        const size_t idx_op_gf,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : DP

        auto to_xx_x = pbuffer.data(idx_op_dp + i * 18 + 0);

        auto to_xx_y = pbuffer.data(idx_op_dp + i * 18 + 1);

        auto to_xx_z = pbuffer.data(idx_op_dp + i * 18 + 2);

        auto to_xy_x = pbuffer.data(idx_op_dp + i * 18 + 3);

        auto to_xy_y = pbuffer.data(idx_op_dp + i * 18 + 4);

        auto to_xy_z = pbuffer.data(idx_op_dp + i * 18 + 5);

        auto to_xz_x = pbuffer.data(idx_op_dp + i * 18 + 6);

        auto to_xz_y = pbuffer.data(idx_op_dp + i * 18 + 7);

        auto to_xz_z = pbuffer.data(idx_op_dp + i * 18 + 8);

        auto to_yy_x = pbuffer.data(idx_op_dp + i * 18 + 9);

        auto to_yy_y = pbuffer.data(idx_op_dp + i * 18 + 10);

        auto to_yy_z = pbuffer.data(idx_op_dp + i * 18 + 11);

        auto to_yz_x = pbuffer.data(idx_op_dp + i * 18 + 12);

        auto to_yz_y = pbuffer.data(idx_op_dp + i * 18 + 13);

        auto to_yz_z = pbuffer.data(idx_op_dp + i * 18 + 14);

        auto to_zz_x = pbuffer.data(idx_op_dp + i * 18 + 15);

        auto to_zz_y = pbuffer.data(idx_op_dp + i * 18 + 16);

        auto to_zz_z = pbuffer.data(idx_op_dp + i * 18 + 17);

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

        // Set up components of auxiliary buffer : GP

        auto to_xxxx_x = pbuffer.data(idx_op_gp + i * 45 + 0);

        auto to_xxxx_y = pbuffer.data(idx_op_gp + i * 45 + 1);

        auto to_xxxx_z = pbuffer.data(idx_op_gp + i * 45 + 2);

        auto to_xxxy_x = pbuffer.data(idx_op_gp + i * 45 + 3);

        auto to_xxxy_y = pbuffer.data(idx_op_gp + i * 45 + 4);

        auto to_xxxy_z = pbuffer.data(idx_op_gp + i * 45 + 5);

        auto to_xxxz_x = pbuffer.data(idx_op_gp + i * 45 + 6);

        auto to_xxxz_y = pbuffer.data(idx_op_gp + i * 45 + 7);

        auto to_xxxz_z = pbuffer.data(idx_op_gp + i * 45 + 8);

        auto to_xxyy_x = pbuffer.data(idx_op_gp + i * 45 + 9);

        auto to_xxyy_y = pbuffer.data(idx_op_gp + i * 45 + 10);

        auto to_xxyy_z = pbuffer.data(idx_op_gp + i * 45 + 11);

        auto to_xxyz_x = pbuffer.data(idx_op_gp + i * 45 + 12);

        auto to_xxyz_y = pbuffer.data(idx_op_gp + i * 45 + 13);

        auto to_xxyz_z = pbuffer.data(idx_op_gp + i * 45 + 14);

        auto to_xxzz_x = pbuffer.data(idx_op_gp + i * 45 + 15);

        auto to_xxzz_y = pbuffer.data(idx_op_gp + i * 45 + 16);

        auto to_xxzz_z = pbuffer.data(idx_op_gp + i * 45 + 17);

        auto to_xyyy_x = pbuffer.data(idx_op_gp + i * 45 + 18);

        auto to_xyyy_y = pbuffer.data(idx_op_gp + i * 45 + 19);

        auto to_xyyy_z = pbuffer.data(idx_op_gp + i * 45 + 20);

        auto to_xyyz_x = pbuffer.data(idx_op_gp + i * 45 + 21);

        auto to_xyyz_y = pbuffer.data(idx_op_gp + i * 45 + 22);

        auto to_xyyz_z = pbuffer.data(idx_op_gp + i * 45 + 23);

        auto to_xyzz_x = pbuffer.data(idx_op_gp + i * 45 + 24);

        auto to_xyzz_y = pbuffer.data(idx_op_gp + i * 45 + 25);

        auto to_xyzz_z = pbuffer.data(idx_op_gp + i * 45 + 26);

        auto to_xzzz_x = pbuffer.data(idx_op_gp + i * 45 + 27);

        auto to_xzzz_y = pbuffer.data(idx_op_gp + i * 45 + 28);

        auto to_xzzz_z = pbuffer.data(idx_op_gp + i * 45 + 29);

        auto to_yyyy_x = pbuffer.data(idx_op_gp + i * 45 + 30);

        auto to_yyyy_y = pbuffer.data(idx_op_gp + i * 45 + 31);

        auto to_yyyy_z = pbuffer.data(idx_op_gp + i * 45 + 32);

        auto to_yyyz_x = pbuffer.data(idx_op_gp + i * 45 + 33);

        auto to_yyyz_y = pbuffer.data(idx_op_gp + i * 45 + 34);

        auto to_yyyz_z = pbuffer.data(idx_op_gp + i * 45 + 35);

        auto to_yyzz_x = pbuffer.data(idx_op_gp + i * 45 + 36);

        auto to_yyzz_y = pbuffer.data(idx_op_gp + i * 45 + 37);

        auto to_yyzz_z = pbuffer.data(idx_op_gp + i * 45 + 38);

        auto to_yzzz_x = pbuffer.data(idx_op_gp + i * 45 + 39);

        auto to_yzzz_y = pbuffer.data(idx_op_gp + i * 45 + 40);

        auto to_yzzz_z = pbuffer.data(idx_op_gp + i * 45 + 41);

        auto to_zzzz_x = pbuffer.data(idx_op_gp + i * 45 + 42);

        auto to_zzzz_y = pbuffer.data(idx_op_gp + i * 45 + 43);

        auto to_zzzz_z = pbuffer.data(idx_op_gp + i * 45 + 44);

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

        // Set up 0-6 components of targeted buffer : FD

        auto to_x_x_xxx_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 0);

        auto to_x_x_xxx_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 1);

        auto to_x_x_xxx_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 2);

        auto to_x_x_xxx_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 3);

        auto to_x_x_xxx_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 4);

        auto to_x_x_xxx_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_x_x_xxx_xx, to_x_x_xxx_xy, to_x_x_xxx_xz, to_x_x_xxx_yy, to_x_x_xxx_yz, to_x_x_xxx_zz, to_xx_x, to_xx_xxx, to_xx_xxy, to_xx_xxz, to_xx_xyy, to_xx_xyz, to_xx_xzz, to_xx_y, to_xx_z, to_xxxx_x, to_xxxx_xxx, to_xxxx_xxy, to_xxxx_xxz, to_xxxx_xyy, to_xxxx_xyz, to_xxxx_xzz, to_xxxx_y, to_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxx_xx[k] = 6.0 * to_xx_x[k] - 6.0 * to_xx_xxx[k] * tke_0 - 4.0 * to_xxxx_x[k] * tbe_0 + 4.0 * to_xxxx_xxx[k] * tbe_0 * tke_0;

            to_x_x_xxx_xy[k] = 3.0 * to_xx_y[k] - 6.0 * to_xx_xxy[k] * tke_0 - 2.0 * to_xxxx_y[k] * tbe_0 + 4.0 * to_xxxx_xxy[k] * tbe_0 * tke_0;

            to_x_x_xxx_xz[k] = 3.0 * to_xx_z[k] - 6.0 * to_xx_xxz[k] * tke_0 - 2.0 * to_xxxx_z[k] * tbe_0 + 4.0 * to_xxxx_xxz[k] * tbe_0 * tke_0;

            to_x_x_xxx_yy[k] = -6.0 * to_xx_xyy[k] * tke_0 + 4.0 * to_xxxx_xyy[k] * tbe_0 * tke_0;

            to_x_x_xxx_yz[k] = -6.0 * to_xx_xyz[k] * tke_0 + 4.0 * to_xxxx_xyz[k] * tbe_0 * tke_0;

            to_x_x_xxx_zz[k] = -6.0 * to_xx_xzz[k] * tke_0 + 4.0 * to_xxxx_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 6-12 components of targeted buffer : FD

        auto to_x_x_xxy_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 6);

        auto to_x_x_xxy_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 7);

        auto to_x_x_xxy_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 8);

        auto to_x_x_xxy_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 9);

        auto to_x_x_xxy_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 10);

        auto to_x_x_xxy_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_x_x_xxy_xx, to_x_x_xxy_xy, to_x_x_xxy_xz, to_x_x_xxy_yy, to_x_x_xxy_yz, to_x_x_xxy_zz, to_xxxy_x, to_xxxy_xxx, to_xxxy_xxy, to_xxxy_xxz, to_xxxy_xyy, to_xxxy_xyz, to_xxxy_xzz, to_xxxy_y, to_xxxy_z, to_xy_x, to_xy_xxx, to_xy_xxy, to_xy_xxz, to_xy_xyy, to_xy_xyz, to_xy_xzz, to_xy_y, to_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxy_xx[k] = 4.0 * to_xy_x[k] - 4.0 * to_xy_xxx[k] * tke_0 - 4.0 * to_xxxy_x[k] * tbe_0 + 4.0 * to_xxxy_xxx[k] * tbe_0 * tke_0;

            to_x_x_xxy_xy[k] = 2.0 * to_xy_y[k] - 4.0 * to_xy_xxy[k] * tke_0 - 2.0 * to_xxxy_y[k] * tbe_0 + 4.0 * to_xxxy_xxy[k] * tbe_0 * tke_0;

            to_x_x_xxy_xz[k] = 2.0 * to_xy_z[k] - 4.0 * to_xy_xxz[k] * tke_0 - 2.0 * to_xxxy_z[k] * tbe_0 + 4.0 * to_xxxy_xxz[k] * tbe_0 * tke_0;

            to_x_x_xxy_yy[k] = -4.0 * to_xy_xyy[k] * tke_0 + 4.0 * to_xxxy_xyy[k] * tbe_0 * tke_0;

            to_x_x_xxy_yz[k] = -4.0 * to_xy_xyz[k] * tke_0 + 4.0 * to_xxxy_xyz[k] * tbe_0 * tke_0;

            to_x_x_xxy_zz[k] = -4.0 * to_xy_xzz[k] * tke_0 + 4.0 * to_xxxy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 12-18 components of targeted buffer : FD

        auto to_x_x_xxz_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 12);

        auto to_x_x_xxz_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 13);

        auto to_x_x_xxz_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 14);

        auto to_x_x_xxz_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 15);

        auto to_x_x_xxz_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 16);

        auto to_x_x_xxz_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_x_x_xxz_xx, to_x_x_xxz_xy, to_x_x_xxz_xz, to_x_x_xxz_yy, to_x_x_xxz_yz, to_x_x_xxz_zz, to_xxxz_x, to_xxxz_xxx, to_xxxz_xxy, to_xxxz_xxz, to_xxxz_xyy, to_xxxz_xyz, to_xxxz_xzz, to_xxxz_y, to_xxxz_z, to_xz_x, to_xz_xxx, to_xz_xxy, to_xz_xxz, to_xz_xyy, to_xz_xyz, to_xz_xzz, to_xz_y, to_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxz_xx[k] = 4.0 * to_xz_x[k] - 4.0 * to_xz_xxx[k] * tke_0 - 4.0 * to_xxxz_x[k] * tbe_0 + 4.0 * to_xxxz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xxz_xy[k] = 2.0 * to_xz_y[k] - 4.0 * to_xz_xxy[k] * tke_0 - 2.0 * to_xxxz_y[k] * tbe_0 + 4.0 * to_xxxz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xxz_xz[k] = 2.0 * to_xz_z[k] - 4.0 * to_xz_xxz[k] * tke_0 - 2.0 * to_xxxz_z[k] * tbe_0 + 4.0 * to_xxxz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xxz_yy[k] = -4.0 * to_xz_xyy[k] * tke_0 + 4.0 * to_xxxz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xxz_yz[k] = -4.0 * to_xz_xyz[k] * tke_0 + 4.0 * to_xxxz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xxz_zz[k] = -4.0 * to_xz_xzz[k] * tke_0 + 4.0 * to_xxxz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 18-24 components of targeted buffer : FD

        auto to_x_x_xyy_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 18);

        auto to_x_x_xyy_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 19);

        auto to_x_x_xyy_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 20);

        auto to_x_x_xyy_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 21);

        auto to_x_x_xyy_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 22);

        auto to_x_x_xyy_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_x_x_xyy_xx, to_x_x_xyy_xy, to_x_x_xyy_xz, to_x_x_xyy_yy, to_x_x_xyy_yz, to_x_x_xyy_zz, to_xxyy_x, to_xxyy_xxx, to_xxyy_xxy, to_xxyy_xxz, to_xxyy_xyy, to_xxyy_xyz, to_xxyy_xzz, to_xxyy_y, to_xxyy_z, to_yy_x, to_yy_xxx, to_yy_xxy, to_yy_xxz, to_yy_xyy, to_yy_xyz, to_yy_xzz, to_yy_y, to_yy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyy_xx[k] = 2.0 * to_yy_x[k] - 2.0 * to_yy_xxx[k] * tke_0 - 4.0 * to_xxyy_x[k] * tbe_0 + 4.0 * to_xxyy_xxx[k] * tbe_0 * tke_0;

            to_x_x_xyy_xy[k] = to_yy_y[k] - 2.0 * to_yy_xxy[k] * tke_0 - 2.0 * to_xxyy_y[k] * tbe_0 + 4.0 * to_xxyy_xxy[k] * tbe_0 * tke_0;

            to_x_x_xyy_xz[k] = to_yy_z[k] - 2.0 * to_yy_xxz[k] * tke_0 - 2.0 * to_xxyy_z[k] * tbe_0 + 4.0 * to_xxyy_xxz[k] * tbe_0 * tke_0;

            to_x_x_xyy_yy[k] = -2.0 * to_yy_xyy[k] * tke_0 + 4.0 * to_xxyy_xyy[k] * tbe_0 * tke_0;

            to_x_x_xyy_yz[k] = -2.0 * to_yy_xyz[k] * tke_0 + 4.0 * to_xxyy_xyz[k] * tbe_0 * tke_0;

            to_x_x_xyy_zz[k] = -2.0 * to_yy_xzz[k] * tke_0 + 4.0 * to_xxyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 24-30 components of targeted buffer : FD

        auto to_x_x_xyz_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 24);

        auto to_x_x_xyz_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 25);

        auto to_x_x_xyz_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 26);

        auto to_x_x_xyz_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 27);

        auto to_x_x_xyz_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 28);

        auto to_x_x_xyz_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_x_x_xyz_xx, to_x_x_xyz_xy, to_x_x_xyz_xz, to_x_x_xyz_yy, to_x_x_xyz_yz, to_x_x_xyz_zz, to_xxyz_x, to_xxyz_xxx, to_xxyz_xxy, to_xxyz_xxz, to_xxyz_xyy, to_xxyz_xyz, to_xxyz_xzz, to_xxyz_y, to_xxyz_z, to_yz_x, to_yz_xxx, to_yz_xxy, to_yz_xxz, to_yz_xyy, to_yz_xyz, to_yz_xzz, to_yz_y, to_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyz_xx[k] = 2.0 * to_yz_x[k] - 2.0 * to_yz_xxx[k] * tke_0 - 4.0 * to_xxyz_x[k] * tbe_0 + 4.0 * to_xxyz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xyz_xy[k] = to_yz_y[k] - 2.0 * to_yz_xxy[k] * tke_0 - 2.0 * to_xxyz_y[k] * tbe_0 + 4.0 * to_xxyz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xyz_xz[k] = to_yz_z[k] - 2.0 * to_yz_xxz[k] * tke_0 - 2.0 * to_xxyz_z[k] * tbe_0 + 4.0 * to_xxyz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xyz_yy[k] = -2.0 * to_yz_xyy[k] * tke_0 + 4.0 * to_xxyz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xyz_yz[k] = -2.0 * to_yz_xyz[k] * tke_0 + 4.0 * to_xxyz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xyz_zz[k] = -2.0 * to_yz_xzz[k] * tke_0 + 4.0 * to_xxyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-36 components of targeted buffer : FD

        auto to_x_x_xzz_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 30);

        auto to_x_x_xzz_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 31);

        auto to_x_x_xzz_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 32);

        auto to_x_x_xzz_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 33);

        auto to_x_x_xzz_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 34);

        auto to_x_x_xzz_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_x_x_xzz_xx, to_x_x_xzz_xy, to_x_x_xzz_xz, to_x_x_xzz_yy, to_x_x_xzz_yz, to_x_x_xzz_zz, to_xxzz_x, to_xxzz_xxx, to_xxzz_xxy, to_xxzz_xxz, to_xxzz_xyy, to_xxzz_xyz, to_xxzz_xzz, to_xxzz_y, to_xxzz_z, to_zz_x, to_zz_xxx, to_zz_xxy, to_zz_xxz, to_zz_xyy, to_zz_xyz, to_zz_xzz, to_zz_y, to_zz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xzz_xx[k] = 2.0 * to_zz_x[k] - 2.0 * to_zz_xxx[k] * tke_0 - 4.0 * to_xxzz_x[k] * tbe_0 + 4.0 * to_xxzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xzz_xy[k] = to_zz_y[k] - 2.0 * to_zz_xxy[k] * tke_0 - 2.0 * to_xxzz_y[k] * tbe_0 + 4.0 * to_xxzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xzz_xz[k] = to_zz_z[k] - 2.0 * to_zz_xxz[k] * tke_0 - 2.0 * to_xxzz_z[k] * tbe_0 + 4.0 * to_xxzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xzz_yy[k] = -2.0 * to_zz_xyy[k] * tke_0 + 4.0 * to_xxzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xzz_yz[k] = -2.0 * to_zz_xyz[k] * tke_0 + 4.0 * to_xxzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xzz_zz[k] = -2.0 * to_zz_xzz[k] * tke_0 + 4.0 * to_xxzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 36-42 components of targeted buffer : FD

        auto to_x_x_yyy_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 36);

        auto to_x_x_yyy_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 37);

        auto to_x_x_yyy_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 38);

        auto to_x_x_yyy_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 39);

        auto to_x_x_yyy_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 40);

        auto to_x_x_yyy_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_x_x_yyy_xx, to_x_x_yyy_xy, to_x_x_yyy_xz, to_x_x_yyy_yy, to_x_x_yyy_yz, to_x_x_yyy_zz, to_xyyy_x, to_xyyy_xxx, to_xyyy_xxy, to_xyyy_xxz, to_xyyy_xyy, to_xyyy_xyz, to_xyyy_xzz, to_xyyy_y, to_xyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyy_xx[k] = -4.0 * to_xyyy_x[k] * tbe_0 + 4.0 * to_xyyy_xxx[k] * tbe_0 * tke_0;

            to_x_x_yyy_xy[k] = -2.0 * to_xyyy_y[k] * tbe_0 + 4.0 * to_xyyy_xxy[k] * tbe_0 * tke_0;

            to_x_x_yyy_xz[k] = -2.0 * to_xyyy_z[k] * tbe_0 + 4.0 * to_xyyy_xxz[k] * tbe_0 * tke_0;

            to_x_x_yyy_yy[k] = 4.0 * to_xyyy_xyy[k] * tbe_0 * tke_0;

            to_x_x_yyy_yz[k] = 4.0 * to_xyyy_xyz[k] * tbe_0 * tke_0;

            to_x_x_yyy_zz[k] = 4.0 * to_xyyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 42-48 components of targeted buffer : FD

        auto to_x_x_yyz_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 42);

        auto to_x_x_yyz_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 43);

        auto to_x_x_yyz_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 44);

        auto to_x_x_yyz_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 45);

        auto to_x_x_yyz_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 46);

        auto to_x_x_yyz_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_x_x_yyz_xx, to_x_x_yyz_xy, to_x_x_yyz_xz, to_x_x_yyz_yy, to_x_x_yyz_yz, to_x_x_yyz_zz, to_xyyz_x, to_xyyz_xxx, to_xyyz_xxy, to_xyyz_xxz, to_xyyz_xyy, to_xyyz_xyz, to_xyyz_xzz, to_xyyz_y, to_xyyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyz_xx[k] = -4.0 * to_xyyz_x[k] * tbe_0 + 4.0 * to_xyyz_xxx[k] * tbe_0 * tke_0;

            to_x_x_yyz_xy[k] = -2.0 * to_xyyz_y[k] * tbe_0 + 4.0 * to_xyyz_xxy[k] * tbe_0 * tke_0;

            to_x_x_yyz_xz[k] = -2.0 * to_xyyz_z[k] * tbe_0 + 4.0 * to_xyyz_xxz[k] * tbe_0 * tke_0;

            to_x_x_yyz_yy[k] = 4.0 * to_xyyz_xyy[k] * tbe_0 * tke_0;

            to_x_x_yyz_yz[k] = 4.0 * to_xyyz_xyz[k] * tbe_0 * tke_0;

            to_x_x_yyz_zz[k] = 4.0 * to_xyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 48-54 components of targeted buffer : FD

        auto to_x_x_yzz_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 48);

        auto to_x_x_yzz_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 49);

        auto to_x_x_yzz_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 50);

        auto to_x_x_yzz_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 51);

        auto to_x_x_yzz_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 52);

        auto to_x_x_yzz_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_x_x_yzz_xx, to_x_x_yzz_xy, to_x_x_yzz_xz, to_x_x_yzz_yy, to_x_x_yzz_yz, to_x_x_yzz_zz, to_xyzz_x, to_xyzz_xxx, to_xyzz_xxy, to_xyzz_xxz, to_xyzz_xyy, to_xyzz_xyz, to_xyzz_xzz, to_xyzz_y, to_xyzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yzz_xx[k] = -4.0 * to_xyzz_x[k] * tbe_0 + 4.0 * to_xyzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_yzz_xy[k] = -2.0 * to_xyzz_y[k] * tbe_0 + 4.0 * to_xyzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_yzz_xz[k] = -2.0 * to_xyzz_z[k] * tbe_0 + 4.0 * to_xyzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_yzz_yy[k] = 4.0 * to_xyzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_yzz_yz[k] = 4.0 * to_xyzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_yzz_zz[k] = 4.0 * to_xyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 54-60 components of targeted buffer : FD

        auto to_x_x_zzz_xx = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 54);

        auto to_x_x_zzz_xy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 55);

        auto to_x_x_zzz_xz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 56);

        auto to_x_x_zzz_yy = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 57);

        auto to_x_x_zzz_yz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 58);

        auto to_x_x_zzz_zz = pbuffer.data(idx_op_geom_101_fd + 0 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_x_x_zzz_xx, to_x_x_zzz_xy, to_x_x_zzz_xz, to_x_x_zzz_yy, to_x_x_zzz_yz, to_x_x_zzz_zz, to_xzzz_x, to_xzzz_xxx, to_xzzz_xxy, to_xzzz_xxz, to_xzzz_xyy, to_xzzz_xyz, to_xzzz_xzz, to_xzzz_y, to_xzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zzz_xx[k] = -4.0 * to_xzzz_x[k] * tbe_0 + 4.0 * to_xzzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_zzz_xy[k] = -2.0 * to_xzzz_y[k] * tbe_0 + 4.0 * to_xzzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_zzz_xz[k] = -2.0 * to_xzzz_z[k] * tbe_0 + 4.0 * to_xzzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_zzz_yy[k] = 4.0 * to_xzzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_zzz_yz[k] = 4.0 * to_xzzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_zzz_zz[k] = 4.0 * to_xzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-66 components of targeted buffer : FD

        auto to_x_y_xxx_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 0);

        auto to_x_y_xxx_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 1);

        auto to_x_y_xxx_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 2);

        auto to_x_y_xxx_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 3);

        auto to_x_y_xxx_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 4);

        auto to_x_y_xxx_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_x_y_xxx_xx, to_x_y_xxx_xy, to_x_y_xxx_xz, to_x_y_xxx_yy, to_x_y_xxx_yz, to_x_y_xxx_zz, to_xx_x, to_xx_xxy, to_xx_xyy, to_xx_xyz, to_xx_y, to_xx_yyy, to_xx_yyz, to_xx_yzz, to_xx_z, to_xxxx_x, to_xxxx_xxy, to_xxxx_xyy, to_xxxx_xyz, to_xxxx_y, to_xxxx_yyy, to_xxxx_yyz, to_xxxx_yzz, to_xxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxx_xx[k] = -6.0 * to_xx_xxy[k] * tke_0 + 4.0 * to_xxxx_xxy[k] * tbe_0 * tke_0;

            to_x_y_xxx_xy[k] = 3.0 * to_xx_x[k] - 6.0 * to_xx_xyy[k] * tke_0 - 2.0 * to_xxxx_x[k] * tbe_0 + 4.0 * to_xxxx_xyy[k] * tbe_0 * tke_0;

            to_x_y_xxx_xz[k] = -6.0 * to_xx_xyz[k] * tke_0 + 4.0 * to_xxxx_xyz[k] * tbe_0 * tke_0;

            to_x_y_xxx_yy[k] = 6.0 * to_xx_y[k] - 6.0 * to_xx_yyy[k] * tke_0 - 4.0 * to_xxxx_y[k] * tbe_0 + 4.0 * to_xxxx_yyy[k] * tbe_0 * tke_0;

            to_x_y_xxx_yz[k] = 3.0 * to_xx_z[k] - 6.0 * to_xx_yyz[k] * tke_0 - 2.0 * to_xxxx_z[k] * tbe_0 + 4.0 * to_xxxx_yyz[k] * tbe_0 * tke_0;

            to_x_y_xxx_zz[k] = -6.0 * to_xx_yzz[k] * tke_0 + 4.0 * to_xxxx_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 66-72 components of targeted buffer : FD

        auto to_x_y_xxy_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 6);

        auto to_x_y_xxy_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 7);

        auto to_x_y_xxy_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 8);

        auto to_x_y_xxy_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 9);

        auto to_x_y_xxy_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 10);

        auto to_x_y_xxy_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_x_y_xxy_xx, to_x_y_xxy_xy, to_x_y_xxy_xz, to_x_y_xxy_yy, to_x_y_xxy_yz, to_x_y_xxy_zz, to_xxxy_x, to_xxxy_xxy, to_xxxy_xyy, to_xxxy_xyz, to_xxxy_y, to_xxxy_yyy, to_xxxy_yyz, to_xxxy_yzz, to_xxxy_z, to_xy_x, to_xy_xxy, to_xy_xyy, to_xy_xyz, to_xy_y, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxy_xx[k] = -4.0 * to_xy_xxy[k] * tke_0 + 4.0 * to_xxxy_xxy[k] * tbe_0 * tke_0;

            to_x_y_xxy_xy[k] = 2.0 * to_xy_x[k] - 4.0 * to_xy_xyy[k] * tke_0 - 2.0 * to_xxxy_x[k] * tbe_0 + 4.0 * to_xxxy_xyy[k] * tbe_0 * tke_0;

            to_x_y_xxy_xz[k] = -4.0 * to_xy_xyz[k] * tke_0 + 4.0 * to_xxxy_xyz[k] * tbe_0 * tke_0;

            to_x_y_xxy_yy[k] = 4.0 * to_xy_y[k] - 4.0 * to_xy_yyy[k] * tke_0 - 4.0 * to_xxxy_y[k] * tbe_0 + 4.0 * to_xxxy_yyy[k] * tbe_0 * tke_0;

            to_x_y_xxy_yz[k] = 2.0 * to_xy_z[k] - 4.0 * to_xy_yyz[k] * tke_0 - 2.0 * to_xxxy_z[k] * tbe_0 + 4.0 * to_xxxy_yyz[k] * tbe_0 * tke_0;

            to_x_y_xxy_zz[k] = -4.0 * to_xy_yzz[k] * tke_0 + 4.0 * to_xxxy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 72-78 components of targeted buffer : FD

        auto to_x_y_xxz_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 12);

        auto to_x_y_xxz_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 13);

        auto to_x_y_xxz_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 14);

        auto to_x_y_xxz_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 15);

        auto to_x_y_xxz_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 16);

        auto to_x_y_xxz_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_x_y_xxz_xx, to_x_y_xxz_xy, to_x_y_xxz_xz, to_x_y_xxz_yy, to_x_y_xxz_yz, to_x_y_xxz_zz, to_xxxz_x, to_xxxz_xxy, to_xxxz_xyy, to_xxxz_xyz, to_xxxz_y, to_xxxz_yyy, to_xxxz_yyz, to_xxxz_yzz, to_xxxz_z, to_xz_x, to_xz_xxy, to_xz_xyy, to_xz_xyz, to_xz_y, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxz_xx[k] = -4.0 * to_xz_xxy[k] * tke_0 + 4.0 * to_xxxz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xxz_xy[k] = 2.0 * to_xz_x[k] - 4.0 * to_xz_xyy[k] * tke_0 - 2.0 * to_xxxz_x[k] * tbe_0 + 4.0 * to_xxxz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xxz_xz[k] = -4.0 * to_xz_xyz[k] * tke_0 + 4.0 * to_xxxz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xxz_yy[k] = 4.0 * to_xz_y[k] - 4.0 * to_xz_yyy[k] * tke_0 - 4.0 * to_xxxz_y[k] * tbe_0 + 4.0 * to_xxxz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xxz_yz[k] = 2.0 * to_xz_z[k] - 4.0 * to_xz_yyz[k] * tke_0 - 2.0 * to_xxxz_z[k] * tbe_0 + 4.0 * to_xxxz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xxz_zz[k] = -4.0 * to_xz_yzz[k] * tke_0 + 4.0 * to_xxxz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 78-84 components of targeted buffer : FD

        auto to_x_y_xyy_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 18);

        auto to_x_y_xyy_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 19);

        auto to_x_y_xyy_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 20);

        auto to_x_y_xyy_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 21);

        auto to_x_y_xyy_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 22);

        auto to_x_y_xyy_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_x_y_xyy_xx, to_x_y_xyy_xy, to_x_y_xyy_xz, to_x_y_xyy_yy, to_x_y_xyy_yz, to_x_y_xyy_zz, to_xxyy_x, to_xxyy_xxy, to_xxyy_xyy, to_xxyy_xyz, to_xxyy_y, to_xxyy_yyy, to_xxyy_yyz, to_xxyy_yzz, to_xxyy_z, to_yy_x, to_yy_xxy, to_yy_xyy, to_yy_xyz, to_yy_y, to_yy_yyy, to_yy_yyz, to_yy_yzz, to_yy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyy_xx[k] = -2.0 * to_yy_xxy[k] * tke_0 + 4.0 * to_xxyy_xxy[k] * tbe_0 * tke_0;

            to_x_y_xyy_xy[k] = to_yy_x[k] - 2.0 * to_yy_xyy[k] * tke_0 - 2.0 * to_xxyy_x[k] * tbe_0 + 4.0 * to_xxyy_xyy[k] * tbe_0 * tke_0;

            to_x_y_xyy_xz[k] = -2.0 * to_yy_xyz[k] * tke_0 + 4.0 * to_xxyy_xyz[k] * tbe_0 * tke_0;

            to_x_y_xyy_yy[k] = 2.0 * to_yy_y[k] - 2.0 * to_yy_yyy[k] * tke_0 - 4.0 * to_xxyy_y[k] * tbe_0 + 4.0 * to_xxyy_yyy[k] * tbe_0 * tke_0;

            to_x_y_xyy_yz[k] = to_yy_z[k] - 2.0 * to_yy_yyz[k] * tke_0 - 2.0 * to_xxyy_z[k] * tbe_0 + 4.0 * to_xxyy_yyz[k] * tbe_0 * tke_0;

            to_x_y_xyy_zz[k] = -2.0 * to_yy_yzz[k] * tke_0 + 4.0 * to_xxyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 84-90 components of targeted buffer : FD

        auto to_x_y_xyz_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 24);

        auto to_x_y_xyz_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 25);

        auto to_x_y_xyz_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 26);

        auto to_x_y_xyz_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 27);

        auto to_x_y_xyz_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 28);

        auto to_x_y_xyz_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_x_y_xyz_xx, to_x_y_xyz_xy, to_x_y_xyz_xz, to_x_y_xyz_yy, to_x_y_xyz_yz, to_x_y_xyz_zz, to_xxyz_x, to_xxyz_xxy, to_xxyz_xyy, to_xxyz_xyz, to_xxyz_y, to_xxyz_yyy, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_z, to_yz_x, to_yz_xxy, to_yz_xyy, to_yz_xyz, to_yz_y, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyz_xx[k] = -2.0 * to_yz_xxy[k] * tke_0 + 4.0 * to_xxyz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xyz_xy[k] = to_yz_x[k] - 2.0 * to_yz_xyy[k] * tke_0 - 2.0 * to_xxyz_x[k] * tbe_0 + 4.0 * to_xxyz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xyz_xz[k] = -2.0 * to_yz_xyz[k] * tke_0 + 4.0 * to_xxyz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xyz_yy[k] = 2.0 * to_yz_y[k] - 2.0 * to_yz_yyy[k] * tke_0 - 4.0 * to_xxyz_y[k] * tbe_0 + 4.0 * to_xxyz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xyz_yz[k] = to_yz_z[k] - 2.0 * to_yz_yyz[k] * tke_0 - 2.0 * to_xxyz_z[k] * tbe_0 + 4.0 * to_xxyz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xyz_zz[k] = -2.0 * to_yz_yzz[k] * tke_0 + 4.0 * to_xxyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-96 components of targeted buffer : FD

        auto to_x_y_xzz_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 30);

        auto to_x_y_xzz_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 31);

        auto to_x_y_xzz_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 32);

        auto to_x_y_xzz_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 33);

        auto to_x_y_xzz_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 34);

        auto to_x_y_xzz_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_x_y_xzz_xx, to_x_y_xzz_xy, to_x_y_xzz_xz, to_x_y_xzz_yy, to_x_y_xzz_yz, to_x_y_xzz_zz, to_xxzz_x, to_xxzz_xxy, to_xxzz_xyy, to_xxzz_xyz, to_xxzz_y, to_xxzz_yyy, to_xxzz_yyz, to_xxzz_yzz, to_xxzz_z, to_zz_x, to_zz_xxy, to_zz_xyy, to_zz_xyz, to_zz_y, to_zz_yyy, to_zz_yyz, to_zz_yzz, to_zz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xzz_xx[k] = -2.0 * to_zz_xxy[k] * tke_0 + 4.0 * to_xxzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xzz_xy[k] = to_zz_x[k] - 2.0 * to_zz_xyy[k] * tke_0 - 2.0 * to_xxzz_x[k] * tbe_0 + 4.0 * to_xxzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xzz_xz[k] = -2.0 * to_zz_xyz[k] * tke_0 + 4.0 * to_xxzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xzz_yy[k] = 2.0 * to_zz_y[k] - 2.0 * to_zz_yyy[k] * tke_0 - 4.0 * to_xxzz_y[k] * tbe_0 + 4.0 * to_xxzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xzz_yz[k] = to_zz_z[k] - 2.0 * to_zz_yyz[k] * tke_0 - 2.0 * to_xxzz_z[k] * tbe_0 + 4.0 * to_xxzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xzz_zz[k] = -2.0 * to_zz_yzz[k] * tke_0 + 4.0 * to_xxzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 96-102 components of targeted buffer : FD

        auto to_x_y_yyy_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 36);

        auto to_x_y_yyy_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 37);

        auto to_x_y_yyy_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 38);

        auto to_x_y_yyy_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 39);

        auto to_x_y_yyy_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 40);

        auto to_x_y_yyy_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_x_y_yyy_xx, to_x_y_yyy_xy, to_x_y_yyy_xz, to_x_y_yyy_yy, to_x_y_yyy_yz, to_x_y_yyy_zz, to_xyyy_x, to_xyyy_xxy, to_xyyy_xyy, to_xyyy_xyz, to_xyyy_y, to_xyyy_yyy, to_xyyy_yyz, to_xyyy_yzz, to_xyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyy_xx[k] = 4.0 * to_xyyy_xxy[k] * tbe_0 * tke_0;

            to_x_y_yyy_xy[k] = -2.0 * to_xyyy_x[k] * tbe_0 + 4.0 * to_xyyy_xyy[k] * tbe_0 * tke_0;

            to_x_y_yyy_xz[k] = 4.0 * to_xyyy_xyz[k] * tbe_0 * tke_0;

            to_x_y_yyy_yy[k] = -4.0 * to_xyyy_y[k] * tbe_0 + 4.0 * to_xyyy_yyy[k] * tbe_0 * tke_0;

            to_x_y_yyy_yz[k] = -2.0 * to_xyyy_z[k] * tbe_0 + 4.0 * to_xyyy_yyz[k] * tbe_0 * tke_0;

            to_x_y_yyy_zz[k] = 4.0 * to_xyyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 102-108 components of targeted buffer : FD

        auto to_x_y_yyz_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 42);

        auto to_x_y_yyz_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 43);

        auto to_x_y_yyz_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 44);

        auto to_x_y_yyz_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 45);

        auto to_x_y_yyz_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 46);

        auto to_x_y_yyz_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_x_y_yyz_xx, to_x_y_yyz_xy, to_x_y_yyz_xz, to_x_y_yyz_yy, to_x_y_yyz_yz, to_x_y_yyz_zz, to_xyyz_x, to_xyyz_xxy, to_xyyz_xyy, to_xyyz_xyz, to_xyyz_y, to_xyyz_yyy, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyz_xx[k] = 4.0 * to_xyyz_xxy[k] * tbe_0 * tke_0;

            to_x_y_yyz_xy[k] = -2.0 * to_xyyz_x[k] * tbe_0 + 4.0 * to_xyyz_xyy[k] * tbe_0 * tke_0;

            to_x_y_yyz_xz[k] = 4.0 * to_xyyz_xyz[k] * tbe_0 * tke_0;

            to_x_y_yyz_yy[k] = -4.0 * to_xyyz_y[k] * tbe_0 + 4.0 * to_xyyz_yyy[k] * tbe_0 * tke_0;

            to_x_y_yyz_yz[k] = -2.0 * to_xyyz_z[k] * tbe_0 + 4.0 * to_xyyz_yyz[k] * tbe_0 * tke_0;

            to_x_y_yyz_zz[k] = 4.0 * to_xyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 108-114 components of targeted buffer : FD

        auto to_x_y_yzz_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 48);

        auto to_x_y_yzz_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 49);

        auto to_x_y_yzz_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 50);

        auto to_x_y_yzz_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 51);

        auto to_x_y_yzz_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 52);

        auto to_x_y_yzz_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_x_y_yzz_xx, to_x_y_yzz_xy, to_x_y_yzz_xz, to_x_y_yzz_yy, to_x_y_yzz_yz, to_x_y_yzz_zz, to_xyzz_x, to_xyzz_xxy, to_xyzz_xyy, to_xyzz_xyz, to_xyzz_y, to_xyzz_yyy, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yzz_xx[k] = 4.0 * to_xyzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_yzz_xy[k] = -2.0 * to_xyzz_x[k] * tbe_0 + 4.0 * to_xyzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_yzz_xz[k] = 4.0 * to_xyzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_yzz_yy[k] = -4.0 * to_xyzz_y[k] * tbe_0 + 4.0 * to_xyzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_yzz_yz[k] = -2.0 * to_xyzz_z[k] * tbe_0 + 4.0 * to_xyzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_yzz_zz[k] = 4.0 * to_xyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 114-120 components of targeted buffer : FD

        auto to_x_y_zzz_xx = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 54);

        auto to_x_y_zzz_xy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 55);

        auto to_x_y_zzz_xz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 56);

        auto to_x_y_zzz_yy = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 57);

        auto to_x_y_zzz_yz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 58);

        auto to_x_y_zzz_zz = pbuffer.data(idx_op_geom_101_fd + 1 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_x_y_zzz_xx, to_x_y_zzz_xy, to_x_y_zzz_xz, to_x_y_zzz_yy, to_x_y_zzz_yz, to_x_y_zzz_zz, to_xzzz_x, to_xzzz_xxy, to_xzzz_xyy, to_xzzz_xyz, to_xzzz_y, to_xzzz_yyy, to_xzzz_yyz, to_xzzz_yzz, to_xzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zzz_xx[k] = 4.0 * to_xzzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_zzz_xy[k] = -2.0 * to_xzzz_x[k] * tbe_0 + 4.0 * to_xzzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_zzz_xz[k] = 4.0 * to_xzzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_zzz_yy[k] = -4.0 * to_xzzz_y[k] * tbe_0 + 4.0 * to_xzzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_zzz_yz[k] = -2.0 * to_xzzz_z[k] * tbe_0 + 4.0 * to_xzzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_zzz_zz[k] = 4.0 * to_xzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-126 components of targeted buffer : FD

        auto to_x_z_xxx_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 0);

        auto to_x_z_xxx_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 1);

        auto to_x_z_xxx_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 2);

        auto to_x_z_xxx_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 3);

        auto to_x_z_xxx_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 4);

        auto to_x_z_xxx_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_x_z_xxx_xx, to_x_z_xxx_xy, to_x_z_xxx_xz, to_x_z_xxx_yy, to_x_z_xxx_yz, to_x_z_xxx_zz, to_xx_x, to_xx_xxz, to_xx_xyz, to_xx_xzz, to_xx_y, to_xx_yyz, to_xx_yzz, to_xx_z, to_xx_zzz, to_xxxx_x, to_xxxx_xxz, to_xxxx_xyz, to_xxxx_xzz, to_xxxx_y, to_xxxx_yyz, to_xxxx_yzz, to_xxxx_z, to_xxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxx_xx[k] = -6.0 * to_xx_xxz[k] * tke_0 + 4.0 * to_xxxx_xxz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xy[k] = -6.0 * to_xx_xyz[k] * tke_0 + 4.0 * to_xxxx_xyz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xz[k] = 3.0 * to_xx_x[k] - 6.0 * to_xx_xzz[k] * tke_0 - 2.0 * to_xxxx_x[k] * tbe_0 + 4.0 * to_xxxx_xzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_yy[k] = -6.0 * to_xx_yyz[k] * tke_0 + 4.0 * to_xxxx_yyz[k] * tbe_0 * tke_0;

            to_x_z_xxx_yz[k] = 3.0 * to_xx_y[k] - 6.0 * to_xx_yzz[k] * tke_0 - 2.0 * to_xxxx_y[k] * tbe_0 + 4.0 * to_xxxx_yzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_zz[k] = 6.0 * to_xx_z[k] - 6.0 * to_xx_zzz[k] * tke_0 - 4.0 * to_xxxx_z[k] * tbe_0 + 4.0 * to_xxxx_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 126-132 components of targeted buffer : FD

        auto to_x_z_xxy_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 6);

        auto to_x_z_xxy_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 7);

        auto to_x_z_xxy_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 8);

        auto to_x_z_xxy_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 9);

        auto to_x_z_xxy_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 10);

        auto to_x_z_xxy_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_x_z_xxy_xx, to_x_z_xxy_xy, to_x_z_xxy_xz, to_x_z_xxy_yy, to_x_z_xxy_yz, to_x_z_xxy_zz, to_xxxy_x, to_xxxy_xxz, to_xxxy_xyz, to_xxxy_xzz, to_xxxy_y, to_xxxy_yyz, to_xxxy_yzz, to_xxxy_z, to_xxxy_zzz, to_xy_x, to_xy_xxz, to_xy_xyz, to_xy_xzz, to_xy_y, to_xy_yyz, to_xy_yzz, to_xy_z, to_xy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxy_xx[k] = -4.0 * to_xy_xxz[k] * tke_0 + 4.0 * to_xxxy_xxz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xy[k] = -4.0 * to_xy_xyz[k] * tke_0 + 4.0 * to_xxxy_xyz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xz[k] = 2.0 * to_xy_x[k] - 4.0 * to_xy_xzz[k] * tke_0 - 2.0 * to_xxxy_x[k] * tbe_0 + 4.0 * to_xxxy_xzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_yy[k] = -4.0 * to_xy_yyz[k] * tke_0 + 4.0 * to_xxxy_yyz[k] * tbe_0 * tke_0;

            to_x_z_xxy_yz[k] = 2.0 * to_xy_y[k] - 4.0 * to_xy_yzz[k] * tke_0 - 2.0 * to_xxxy_y[k] * tbe_0 + 4.0 * to_xxxy_yzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_zz[k] = 4.0 * to_xy_z[k] - 4.0 * to_xy_zzz[k] * tke_0 - 4.0 * to_xxxy_z[k] * tbe_0 + 4.0 * to_xxxy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 132-138 components of targeted buffer : FD

        auto to_x_z_xxz_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 12);

        auto to_x_z_xxz_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 13);

        auto to_x_z_xxz_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 14);

        auto to_x_z_xxz_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 15);

        auto to_x_z_xxz_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 16);

        auto to_x_z_xxz_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_x_z_xxz_xx, to_x_z_xxz_xy, to_x_z_xxz_xz, to_x_z_xxz_yy, to_x_z_xxz_yz, to_x_z_xxz_zz, to_xxxz_x, to_xxxz_xxz, to_xxxz_xyz, to_xxxz_xzz, to_xxxz_y, to_xxxz_yyz, to_xxxz_yzz, to_xxxz_z, to_xxxz_zzz, to_xz_x, to_xz_xxz, to_xz_xyz, to_xz_xzz, to_xz_y, to_xz_yyz, to_xz_yzz, to_xz_z, to_xz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxz_xx[k] = -4.0 * to_xz_xxz[k] * tke_0 + 4.0 * to_xxxz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xy[k] = -4.0 * to_xz_xyz[k] * tke_0 + 4.0 * to_xxxz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xz[k] = 2.0 * to_xz_x[k] - 4.0 * to_xz_xzz[k] * tke_0 - 2.0 * to_xxxz_x[k] * tbe_0 + 4.0 * to_xxxz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_yy[k] = -4.0 * to_xz_yyz[k] * tke_0 + 4.0 * to_xxxz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xxz_yz[k] = 2.0 * to_xz_y[k] - 4.0 * to_xz_yzz[k] * tke_0 - 2.0 * to_xxxz_y[k] * tbe_0 + 4.0 * to_xxxz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_zz[k] = 4.0 * to_xz_z[k] - 4.0 * to_xz_zzz[k] * tke_0 - 4.0 * to_xxxz_z[k] * tbe_0 + 4.0 * to_xxxz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 138-144 components of targeted buffer : FD

        auto to_x_z_xyy_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 18);

        auto to_x_z_xyy_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 19);

        auto to_x_z_xyy_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 20);

        auto to_x_z_xyy_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 21);

        auto to_x_z_xyy_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 22);

        auto to_x_z_xyy_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_x_z_xyy_xx, to_x_z_xyy_xy, to_x_z_xyy_xz, to_x_z_xyy_yy, to_x_z_xyy_yz, to_x_z_xyy_zz, to_xxyy_x, to_xxyy_xxz, to_xxyy_xyz, to_xxyy_xzz, to_xxyy_y, to_xxyy_yyz, to_xxyy_yzz, to_xxyy_z, to_xxyy_zzz, to_yy_x, to_yy_xxz, to_yy_xyz, to_yy_xzz, to_yy_y, to_yy_yyz, to_yy_yzz, to_yy_z, to_yy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyy_xx[k] = -2.0 * to_yy_xxz[k] * tke_0 + 4.0 * to_xxyy_xxz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xy[k] = -2.0 * to_yy_xyz[k] * tke_0 + 4.0 * to_xxyy_xyz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xz[k] = to_yy_x[k] - 2.0 * to_yy_xzz[k] * tke_0 - 2.0 * to_xxyy_x[k] * tbe_0 + 4.0 * to_xxyy_xzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_yy[k] = -2.0 * to_yy_yyz[k] * tke_0 + 4.0 * to_xxyy_yyz[k] * tbe_0 * tke_0;

            to_x_z_xyy_yz[k] = to_yy_y[k] - 2.0 * to_yy_yzz[k] * tke_0 - 2.0 * to_xxyy_y[k] * tbe_0 + 4.0 * to_xxyy_yzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_zz[k] = 2.0 * to_yy_z[k] - 2.0 * to_yy_zzz[k] * tke_0 - 4.0 * to_xxyy_z[k] * tbe_0 + 4.0 * to_xxyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 144-150 components of targeted buffer : FD

        auto to_x_z_xyz_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 24);

        auto to_x_z_xyz_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 25);

        auto to_x_z_xyz_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 26);

        auto to_x_z_xyz_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 27);

        auto to_x_z_xyz_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 28);

        auto to_x_z_xyz_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_x_z_xyz_xx, to_x_z_xyz_xy, to_x_z_xyz_xz, to_x_z_xyz_yy, to_x_z_xyz_yz, to_x_z_xyz_zz, to_xxyz_x, to_xxyz_xxz, to_xxyz_xyz, to_xxyz_xzz, to_xxyz_y, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_z, to_xxyz_zzz, to_yz_x, to_yz_xxz, to_yz_xyz, to_yz_xzz, to_yz_y, to_yz_yyz, to_yz_yzz, to_yz_z, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyz_xx[k] = -2.0 * to_yz_xxz[k] * tke_0 + 4.0 * to_xxyz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xy[k] = -2.0 * to_yz_xyz[k] * tke_0 + 4.0 * to_xxyz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xz[k] = to_yz_x[k] - 2.0 * to_yz_xzz[k] * tke_0 - 2.0 * to_xxyz_x[k] * tbe_0 + 4.0 * to_xxyz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_yy[k] = -2.0 * to_yz_yyz[k] * tke_0 + 4.0 * to_xxyz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xyz_yz[k] = to_yz_y[k] - 2.0 * to_yz_yzz[k] * tke_0 - 2.0 * to_xxyz_y[k] * tbe_0 + 4.0 * to_xxyz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_zz[k] = 2.0 * to_yz_z[k] - 2.0 * to_yz_zzz[k] * tke_0 - 4.0 * to_xxyz_z[k] * tbe_0 + 4.0 * to_xxyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-156 components of targeted buffer : FD

        auto to_x_z_xzz_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 30);

        auto to_x_z_xzz_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 31);

        auto to_x_z_xzz_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 32);

        auto to_x_z_xzz_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 33);

        auto to_x_z_xzz_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 34);

        auto to_x_z_xzz_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_x_z_xzz_xx, to_x_z_xzz_xy, to_x_z_xzz_xz, to_x_z_xzz_yy, to_x_z_xzz_yz, to_x_z_xzz_zz, to_xxzz_x, to_xxzz_xxz, to_xxzz_xyz, to_xxzz_xzz, to_xxzz_y, to_xxzz_yyz, to_xxzz_yzz, to_xxzz_z, to_xxzz_zzz, to_zz_x, to_zz_xxz, to_zz_xyz, to_zz_xzz, to_zz_y, to_zz_yyz, to_zz_yzz, to_zz_z, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xzz_xx[k] = -2.0 * to_zz_xxz[k] * tke_0 + 4.0 * to_xxzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xy[k] = -2.0 * to_zz_xyz[k] * tke_0 + 4.0 * to_xxzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xz[k] = to_zz_x[k] - 2.0 * to_zz_xzz[k] * tke_0 - 2.0 * to_xxzz_x[k] * tbe_0 + 4.0 * to_xxzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_yy[k] = -2.0 * to_zz_yyz[k] * tke_0 + 4.0 * to_xxzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xzz_yz[k] = to_zz_y[k] - 2.0 * to_zz_yzz[k] * tke_0 - 2.0 * to_xxzz_y[k] * tbe_0 + 4.0 * to_xxzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_zz[k] = 2.0 * to_zz_z[k] - 2.0 * to_zz_zzz[k] * tke_0 - 4.0 * to_xxzz_z[k] * tbe_0 + 4.0 * to_xxzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 156-162 components of targeted buffer : FD

        auto to_x_z_yyy_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 36);

        auto to_x_z_yyy_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 37);

        auto to_x_z_yyy_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 38);

        auto to_x_z_yyy_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 39);

        auto to_x_z_yyy_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 40);

        auto to_x_z_yyy_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_x_z_yyy_xx, to_x_z_yyy_xy, to_x_z_yyy_xz, to_x_z_yyy_yy, to_x_z_yyy_yz, to_x_z_yyy_zz, to_xyyy_x, to_xyyy_xxz, to_xyyy_xyz, to_xyyy_xzz, to_xyyy_y, to_xyyy_yyz, to_xyyy_yzz, to_xyyy_z, to_xyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyy_xx[k] = 4.0 * to_xyyy_xxz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xy[k] = 4.0 * to_xyyy_xyz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xz[k] = -2.0 * to_xyyy_x[k] * tbe_0 + 4.0 * to_xyyy_xzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_yy[k] = 4.0 * to_xyyy_yyz[k] * tbe_0 * tke_0;

            to_x_z_yyy_yz[k] = -2.0 * to_xyyy_y[k] * tbe_0 + 4.0 * to_xyyy_yzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_zz[k] = -4.0 * to_xyyy_z[k] * tbe_0 + 4.0 * to_xyyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 162-168 components of targeted buffer : FD

        auto to_x_z_yyz_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 42);

        auto to_x_z_yyz_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 43);

        auto to_x_z_yyz_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 44);

        auto to_x_z_yyz_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 45);

        auto to_x_z_yyz_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 46);

        auto to_x_z_yyz_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_x_z_yyz_xx, to_x_z_yyz_xy, to_x_z_yyz_xz, to_x_z_yyz_yy, to_x_z_yyz_yz, to_x_z_yyz_zz, to_xyyz_x, to_xyyz_xxz, to_xyyz_xyz, to_xyyz_xzz, to_xyyz_y, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_z, to_xyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyz_xx[k] = 4.0 * to_xyyz_xxz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xy[k] = 4.0 * to_xyyz_xyz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xz[k] = -2.0 * to_xyyz_x[k] * tbe_0 + 4.0 * to_xyyz_xzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_yy[k] = 4.0 * to_xyyz_yyz[k] * tbe_0 * tke_0;

            to_x_z_yyz_yz[k] = -2.0 * to_xyyz_y[k] * tbe_0 + 4.0 * to_xyyz_yzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_zz[k] = -4.0 * to_xyyz_z[k] * tbe_0 + 4.0 * to_xyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 168-174 components of targeted buffer : FD

        auto to_x_z_yzz_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 48);

        auto to_x_z_yzz_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 49);

        auto to_x_z_yzz_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 50);

        auto to_x_z_yzz_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 51);

        auto to_x_z_yzz_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 52);

        auto to_x_z_yzz_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_x_z_yzz_xx, to_x_z_yzz_xy, to_x_z_yzz_xz, to_x_z_yzz_yy, to_x_z_yzz_yz, to_x_z_yzz_zz, to_xyzz_x, to_xyzz_xxz, to_xyzz_xyz, to_xyzz_xzz, to_xyzz_y, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_z, to_xyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yzz_xx[k] = 4.0 * to_xyzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xy[k] = 4.0 * to_xyzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xz[k] = -2.0 * to_xyzz_x[k] * tbe_0 + 4.0 * to_xyzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_yy[k] = 4.0 * to_xyzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_yzz_yz[k] = -2.0 * to_xyzz_y[k] * tbe_0 + 4.0 * to_xyzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_zz[k] = -4.0 * to_xyzz_z[k] * tbe_0 + 4.0 * to_xyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 174-180 components of targeted buffer : FD

        auto to_x_z_zzz_xx = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 54);

        auto to_x_z_zzz_xy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 55);

        auto to_x_z_zzz_xz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 56);

        auto to_x_z_zzz_yy = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 57);

        auto to_x_z_zzz_yz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 58);

        auto to_x_z_zzz_zz = pbuffer.data(idx_op_geom_101_fd + 2 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_x_z_zzz_xx, to_x_z_zzz_xy, to_x_z_zzz_xz, to_x_z_zzz_yy, to_x_z_zzz_yz, to_x_z_zzz_zz, to_xzzz_x, to_xzzz_xxz, to_xzzz_xyz, to_xzzz_xzz, to_xzzz_y, to_xzzz_yyz, to_xzzz_yzz, to_xzzz_z, to_xzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zzz_xx[k] = 4.0 * to_xzzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xy[k] = 4.0 * to_xzzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xz[k] = -2.0 * to_xzzz_x[k] * tbe_0 + 4.0 * to_xzzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_yy[k] = 4.0 * to_xzzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_zzz_yz[k] = -2.0 * to_xzzz_y[k] * tbe_0 + 4.0 * to_xzzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_zz[k] = -4.0 * to_xzzz_z[k] * tbe_0 + 4.0 * to_xzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-186 components of targeted buffer : FD

        auto to_y_x_xxx_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 0);

        auto to_y_x_xxx_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 1);

        auto to_y_x_xxx_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 2);

        auto to_y_x_xxx_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 3);

        auto to_y_x_xxx_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 4);

        auto to_y_x_xxx_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_xxxy_x, to_xxxy_xxx, to_xxxy_xxy, to_xxxy_xxz, to_xxxy_xyy, to_xxxy_xyz, to_xxxy_xzz, to_xxxy_y, to_xxxy_z, to_y_x_xxx_xx, to_y_x_xxx_xy, to_y_x_xxx_xz, to_y_x_xxx_yy, to_y_x_xxx_yz, to_y_x_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxx_xx[k] = -4.0 * to_xxxy_x[k] * tbe_0 + 4.0 * to_xxxy_xxx[k] * tbe_0 * tke_0;

            to_y_x_xxx_xy[k] = -2.0 * to_xxxy_y[k] * tbe_0 + 4.0 * to_xxxy_xxy[k] * tbe_0 * tke_0;

            to_y_x_xxx_xz[k] = -2.0 * to_xxxy_z[k] * tbe_0 + 4.0 * to_xxxy_xxz[k] * tbe_0 * tke_0;

            to_y_x_xxx_yy[k] = 4.0 * to_xxxy_xyy[k] * tbe_0 * tke_0;

            to_y_x_xxx_yz[k] = 4.0 * to_xxxy_xyz[k] * tbe_0 * tke_0;

            to_y_x_xxx_zz[k] = 4.0 * to_xxxy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 186-192 components of targeted buffer : FD

        auto to_y_x_xxy_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 6);

        auto to_y_x_xxy_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 7);

        auto to_y_x_xxy_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 8);

        auto to_y_x_xxy_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 9);

        auto to_y_x_xxy_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 10);

        auto to_y_x_xxy_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_xx_x, to_xx_xxx, to_xx_xxy, to_xx_xxz, to_xx_xyy, to_xx_xyz, to_xx_xzz, to_xx_y, to_xx_z, to_xxyy_x, to_xxyy_xxx, to_xxyy_xxy, to_xxyy_xxz, to_xxyy_xyy, to_xxyy_xyz, to_xxyy_xzz, to_xxyy_y, to_xxyy_z, to_y_x_xxy_xx, to_y_x_xxy_xy, to_y_x_xxy_xz, to_y_x_xxy_yy, to_y_x_xxy_yz, to_y_x_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxy_xx[k] = 2.0 * to_xx_x[k] - 2.0 * to_xx_xxx[k] * tke_0 - 4.0 * to_xxyy_x[k] * tbe_0 + 4.0 * to_xxyy_xxx[k] * tbe_0 * tke_0;

            to_y_x_xxy_xy[k] = to_xx_y[k] - 2.0 * to_xx_xxy[k] * tke_0 - 2.0 * to_xxyy_y[k] * tbe_0 + 4.0 * to_xxyy_xxy[k] * tbe_0 * tke_0;

            to_y_x_xxy_xz[k] = to_xx_z[k] - 2.0 * to_xx_xxz[k] * tke_0 - 2.0 * to_xxyy_z[k] * tbe_0 + 4.0 * to_xxyy_xxz[k] * tbe_0 * tke_0;

            to_y_x_xxy_yy[k] = -2.0 * to_xx_xyy[k] * tke_0 + 4.0 * to_xxyy_xyy[k] * tbe_0 * tke_0;

            to_y_x_xxy_yz[k] = -2.0 * to_xx_xyz[k] * tke_0 + 4.0 * to_xxyy_xyz[k] * tbe_0 * tke_0;

            to_y_x_xxy_zz[k] = -2.0 * to_xx_xzz[k] * tke_0 + 4.0 * to_xxyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 192-198 components of targeted buffer : FD

        auto to_y_x_xxz_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 12);

        auto to_y_x_xxz_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 13);

        auto to_y_x_xxz_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 14);

        auto to_y_x_xxz_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 15);

        auto to_y_x_xxz_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 16);

        auto to_y_x_xxz_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_xxyz_x, to_xxyz_xxx, to_xxyz_xxy, to_xxyz_xxz, to_xxyz_xyy, to_xxyz_xyz, to_xxyz_xzz, to_xxyz_y, to_xxyz_z, to_y_x_xxz_xx, to_y_x_xxz_xy, to_y_x_xxz_xz, to_y_x_xxz_yy, to_y_x_xxz_yz, to_y_x_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxz_xx[k] = -4.0 * to_xxyz_x[k] * tbe_0 + 4.0 * to_xxyz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xxz_xy[k] = -2.0 * to_xxyz_y[k] * tbe_0 + 4.0 * to_xxyz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xxz_xz[k] = -2.0 * to_xxyz_z[k] * tbe_0 + 4.0 * to_xxyz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xxz_yy[k] = 4.0 * to_xxyz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xxz_yz[k] = 4.0 * to_xxyz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xxz_zz[k] = 4.0 * to_xxyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 198-204 components of targeted buffer : FD

        auto to_y_x_xyy_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 18);

        auto to_y_x_xyy_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 19);

        auto to_y_x_xyy_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 20);

        auto to_y_x_xyy_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 21);

        auto to_y_x_xyy_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 22);

        auto to_y_x_xyy_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_xy_x, to_xy_xxx, to_xy_xxy, to_xy_xxz, to_xy_xyy, to_xy_xyz, to_xy_xzz, to_xy_y, to_xy_z, to_xyyy_x, to_xyyy_xxx, to_xyyy_xxy, to_xyyy_xxz, to_xyyy_xyy, to_xyyy_xyz, to_xyyy_xzz, to_xyyy_y, to_xyyy_z, to_y_x_xyy_xx, to_y_x_xyy_xy, to_y_x_xyy_xz, to_y_x_xyy_yy, to_y_x_xyy_yz, to_y_x_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyy_xx[k] = 4.0 * to_xy_x[k] - 4.0 * to_xy_xxx[k] * tke_0 - 4.0 * to_xyyy_x[k] * tbe_0 + 4.0 * to_xyyy_xxx[k] * tbe_0 * tke_0;

            to_y_x_xyy_xy[k] = 2.0 * to_xy_y[k] - 4.0 * to_xy_xxy[k] * tke_0 - 2.0 * to_xyyy_y[k] * tbe_0 + 4.0 * to_xyyy_xxy[k] * tbe_0 * tke_0;

            to_y_x_xyy_xz[k] = 2.0 * to_xy_z[k] - 4.0 * to_xy_xxz[k] * tke_0 - 2.0 * to_xyyy_z[k] * tbe_0 + 4.0 * to_xyyy_xxz[k] * tbe_0 * tke_0;

            to_y_x_xyy_yy[k] = -4.0 * to_xy_xyy[k] * tke_0 + 4.0 * to_xyyy_xyy[k] * tbe_0 * tke_0;

            to_y_x_xyy_yz[k] = -4.0 * to_xy_xyz[k] * tke_0 + 4.0 * to_xyyy_xyz[k] * tbe_0 * tke_0;

            to_y_x_xyy_zz[k] = -4.0 * to_xy_xzz[k] * tke_0 + 4.0 * to_xyyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 204-210 components of targeted buffer : FD

        auto to_y_x_xyz_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 24);

        auto to_y_x_xyz_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 25);

        auto to_y_x_xyz_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 26);

        auto to_y_x_xyz_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 27);

        auto to_y_x_xyz_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 28);

        auto to_y_x_xyz_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_xyyz_x, to_xyyz_xxx, to_xyyz_xxy, to_xyyz_xxz, to_xyyz_xyy, to_xyyz_xyz, to_xyyz_xzz, to_xyyz_y, to_xyyz_z, to_xz_x, to_xz_xxx, to_xz_xxy, to_xz_xxz, to_xz_xyy, to_xz_xyz, to_xz_xzz, to_xz_y, to_xz_z, to_y_x_xyz_xx, to_y_x_xyz_xy, to_y_x_xyz_xz, to_y_x_xyz_yy, to_y_x_xyz_yz, to_y_x_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyz_xx[k] = 2.0 * to_xz_x[k] - 2.0 * to_xz_xxx[k] * tke_0 - 4.0 * to_xyyz_x[k] * tbe_0 + 4.0 * to_xyyz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xyz_xy[k] = to_xz_y[k] - 2.0 * to_xz_xxy[k] * tke_0 - 2.0 * to_xyyz_y[k] * tbe_0 + 4.0 * to_xyyz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xyz_xz[k] = to_xz_z[k] - 2.0 * to_xz_xxz[k] * tke_0 - 2.0 * to_xyyz_z[k] * tbe_0 + 4.0 * to_xyyz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xyz_yy[k] = -2.0 * to_xz_xyy[k] * tke_0 + 4.0 * to_xyyz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xyz_yz[k] = -2.0 * to_xz_xyz[k] * tke_0 + 4.0 * to_xyyz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xyz_zz[k] = -2.0 * to_xz_xzz[k] * tke_0 + 4.0 * to_xyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-216 components of targeted buffer : FD

        auto to_y_x_xzz_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 30);

        auto to_y_x_xzz_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 31);

        auto to_y_x_xzz_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 32);

        auto to_y_x_xzz_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 33);

        auto to_y_x_xzz_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 34);

        auto to_y_x_xzz_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_xyzz_x, to_xyzz_xxx, to_xyzz_xxy, to_xyzz_xxz, to_xyzz_xyy, to_xyzz_xyz, to_xyzz_xzz, to_xyzz_y, to_xyzz_z, to_y_x_xzz_xx, to_y_x_xzz_xy, to_y_x_xzz_xz, to_y_x_xzz_yy, to_y_x_xzz_yz, to_y_x_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xzz_xx[k] = -4.0 * to_xyzz_x[k] * tbe_0 + 4.0 * to_xyzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xzz_xy[k] = -2.0 * to_xyzz_y[k] * tbe_0 + 4.0 * to_xyzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xzz_xz[k] = -2.0 * to_xyzz_z[k] * tbe_0 + 4.0 * to_xyzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xzz_yy[k] = 4.0 * to_xyzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xzz_yz[k] = 4.0 * to_xyzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xzz_zz[k] = 4.0 * to_xyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 216-222 components of targeted buffer : FD

        auto to_y_x_yyy_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 36);

        auto to_y_x_yyy_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 37);

        auto to_y_x_yyy_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 38);

        auto to_y_x_yyy_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 39);

        auto to_y_x_yyy_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 40);

        auto to_y_x_yyy_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_y_x_yyy_xx, to_y_x_yyy_xy, to_y_x_yyy_xz, to_y_x_yyy_yy, to_y_x_yyy_yz, to_y_x_yyy_zz, to_yy_x, to_yy_xxx, to_yy_xxy, to_yy_xxz, to_yy_xyy, to_yy_xyz, to_yy_xzz, to_yy_y, to_yy_z, to_yyyy_x, to_yyyy_xxx, to_yyyy_xxy, to_yyyy_xxz, to_yyyy_xyy, to_yyyy_xyz, to_yyyy_xzz, to_yyyy_y, to_yyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyy_xx[k] = 6.0 * to_yy_x[k] - 6.0 * to_yy_xxx[k] * tke_0 - 4.0 * to_yyyy_x[k] * tbe_0 + 4.0 * to_yyyy_xxx[k] * tbe_0 * tke_0;

            to_y_x_yyy_xy[k] = 3.0 * to_yy_y[k] - 6.0 * to_yy_xxy[k] * tke_0 - 2.0 * to_yyyy_y[k] * tbe_0 + 4.0 * to_yyyy_xxy[k] * tbe_0 * tke_0;

            to_y_x_yyy_xz[k] = 3.0 * to_yy_z[k] - 6.0 * to_yy_xxz[k] * tke_0 - 2.0 * to_yyyy_z[k] * tbe_0 + 4.0 * to_yyyy_xxz[k] * tbe_0 * tke_0;

            to_y_x_yyy_yy[k] = -6.0 * to_yy_xyy[k] * tke_0 + 4.0 * to_yyyy_xyy[k] * tbe_0 * tke_0;

            to_y_x_yyy_yz[k] = -6.0 * to_yy_xyz[k] * tke_0 + 4.0 * to_yyyy_xyz[k] * tbe_0 * tke_0;

            to_y_x_yyy_zz[k] = -6.0 * to_yy_xzz[k] * tke_0 + 4.0 * to_yyyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 222-228 components of targeted buffer : FD

        auto to_y_x_yyz_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 42);

        auto to_y_x_yyz_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 43);

        auto to_y_x_yyz_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 44);

        auto to_y_x_yyz_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 45);

        auto to_y_x_yyz_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 46);

        auto to_y_x_yyz_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_y_x_yyz_xx, to_y_x_yyz_xy, to_y_x_yyz_xz, to_y_x_yyz_yy, to_y_x_yyz_yz, to_y_x_yyz_zz, to_yyyz_x, to_yyyz_xxx, to_yyyz_xxy, to_yyyz_xxz, to_yyyz_xyy, to_yyyz_xyz, to_yyyz_xzz, to_yyyz_y, to_yyyz_z, to_yz_x, to_yz_xxx, to_yz_xxy, to_yz_xxz, to_yz_xyy, to_yz_xyz, to_yz_xzz, to_yz_y, to_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyz_xx[k] = 4.0 * to_yz_x[k] - 4.0 * to_yz_xxx[k] * tke_0 - 4.0 * to_yyyz_x[k] * tbe_0 + 4.0 * to_yyyz_xxx[k] * tbe_0 * tke_0;

            to_y_x_yyz_xy[k] = 2.0 * to_yz_y[k] - 4.0 * to_yz_xxy[k] * tke_0 - 2.0 * to_yyyz_y[k] * tbe_0 + 4.0 * to_yyyz_xxy[k] * tbe_0 * tke_0;

            to_y_x_yyz_xz[k] = 2.0 * to_yz_z[k] - 4.0 * to_yz_xxz[k] * tke_0 - 2.0 * to_yyyz_z[k] * tbe_0 + 4.0 * to_yyyz_xxz[k] * tbe_0 * tke_0;

            to_y_x_yyz_yy[k] = -4.0 * to_yz_xyy[k] * tke_0 + 4.0 * to_yyyz_xyy[k] * tbe_0 * tke_0;

            to_y_x_yyz_yz[k] = -4.0 * to_yz_xyz[k] * tke_0 + 4.0 * to_yyyz_xyz[k] * tbe_0 * tke_0;

            to_y_x_yyz_zz[k] = -4.0 * to_yz_xzz[k] * tke_0 + 4.0 * to_yyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 228-234 components of targeted buffer : FD

        auto to_y_x_yzz_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 48);

        auto to_y_x_yzz_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 49);

        auto to_y_x_yzz_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 50);

        auto to_y_x_yzz_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 51);

        auto to_y_x_yzz_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 52);

        auto to_y_x_yzz_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_y_x_yzz_xx, to_y_x_yzz_xy, to_y_x_yzz_xz, to_y_x_yzz_yy, to_y_x_yzz_yz, to_y_x_yzz_zz, to_yyzz_x, to_yyzz_xxx, to_yyzz_xxy, to_yyzz_xxz, to_yyzz_xyy, to_yyzz_xyz, to_yyzz_xzz, to_yyzz_y, to_yyzz_z, to_zz_x, to_zz_xxx, to_zz_xxy, to_zz_xxz, to_zz_xyy, to_zz_xyz, to_zz_xzz, to_zz_y, to_zz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yzz_xx[k] = 2.0 * to_zz_x[k] - 2.0 * to_zz_xxx[k] * tke_0 - 4.0 * to_yyzz_x[k] * tbe_0 + 4.0 * to_yyzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_yzz_xy[k] = to_zz_y[k] - 2.0 * to_zz_xxy[k] * tke_0 - 2.0 * to_yyzz_y[k] * tbe_0 + 4.0 * to_yyzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_yzz_xz[k] = to_zz_z[k] - 2.0 * to_zz_xxz[k] * tke_0 - 2.0 * to_yyzz_z[k] * tbe_0 + 4.0 * to_yyzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_yzz_yy[k] = -2.0 * to_zz_xyy[k] * tke_0 + 4.0 * to_yyzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_yzz_yz[k] = -2.0 * to_zz_xyz[k] * tke_0 + 4.0 * to_yyzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_yzz_zz[k] = -2.0 * to_zz_xzz[k] * tke_0 + 4.0 * to_yyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 234-240 components of targeted buffer : FD

        auto to_y_x_zzz_xx = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 54);

        auto to_y_x_zzz_xy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 55);

        auto to_y_x_zzz_xz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 56);

        auto to_y_x_zzz_yy = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 57);

        auto to_y_x_zzz_yz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 58);

        auto to_y_x_zzz_zz = pbuffer.data(idx_op_geom_101_fd + 3 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_y_x_zzz_xx, to_y_x_zzz_xy, to_y_x_zzz_xz, to_y_x_zzz_yy, to_y_x_zzz_yz, to_y_x_zzz_zz, to_yzzz_x, to_yzzz_xxx, to_yzzz_xxy, to_yzzz_xxz, to_yzzz_xyy, to_yzzz_xyz, to_yzzz_xzz, to_yzzz_y, to_yzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zzz_xx[k] = -4.0 * to_yzzz_x[k] * tbe_0 + 4.0 * to_yzzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_zzz_xy[k] = -2.0 * to_yzzz_y[k] * tbe_0 + 4.0 * to_yzzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_zzz_xz[k] = -2.0 * to_yzzz_z[k] * tbe_0 + 4.0 * to_yzzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_zzz_yy[k] = 4.0 * to_yzzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_zzz_yz[k] = 4.0 * to_yzzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_zzz_zz[k] = 4.0 * to_yzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-246 components of targeted buffer : FD

        auto to_y_y_xxx_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 0);

        auto to_y_y_xxx_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 1);

        auto to_y_y_xxx_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 2);

        auto to_y_y_xxx_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 3);

        auto to_y_y_xxx_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 4);

        auto to_y_y_xxx_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_xxxy_x, to_xxxy_xxy, to_xxxy_xyy, to_xxxy_xyz, to_xxxy_y, to_xxxy_yyy, to_xxxy_yyz, to_xxxy_yzz, to_xxxy_z, to_y_y_xxx_xx, to_y_y_xxx_xy, to_y_y_xxx_xz, to_y_y_xxx_yy, to_y_y_xxx_yz, to_y_y_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxx_xx[k] = 4.0 * to_xxxy_xxy[k] * tbe_0 * tke_0;

            to_y_y_xxx_xy[k] = -2.0 * to_xxxy_x[k] * tbe_0 + 4.0 * to_xxxy_xyy[k] * tbe_0 * tke_0;

            to_y_y_xxx_xz[k] = 4.0 * to_xxxy_xyz[k] * tbe_0 * tke_0;

            to_y_y_xxx_yy[k] = -4.0 * to_xxxy_y[k] * tbe_0 + 4.0 * to_xxxy_yyy[k] * tbe_0 * tke_0;

            to_y_y_xxx_yz[k] = -2.0 * to_xxxy_z[k] * tbe_0 + 4.0 * to_xxxy_yyz[k] * tbe_0 * tke_0;

            to_y_y_xxx_zz[k] = 4.0 * to_xxxy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 246-252 components of targeted buffer : FD

        auto to_y_y_xxy_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 6);

        auto to_y_y_xxy_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 7);

        auto to_y_y_xxy_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 8);

        auto to_y_y_xxy_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 9);

        auto to_y_y_xxy_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 10);

        auto to_y_y_xxy_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_xx_x, to_xx_xxy, to_xx_xyy, to_xx_xyz, to_xx_y, to_xx_yyy, to_xx_yyz, to_xx_yzz, to_xx_z, to_xxyy_x, to_xxyy_xxy, to_xxyy_xyy, to_xxyy_xyz, to_xxyy_y, to_xxyy_yyy, to_xxyy_yyz, to_xxyy_yzz, to_xxyy_z, to_y_y_xxy_xx, to_y_y_xxy_xy, to_y_y_xxy_xz, to_y_y_xxy_yy, to_y_y_xxy_yz, to_y_y_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxy_xx[k] = -2.0 * to_xx_xxy[k] * tke_0 + 4.0 * to_xxyy_xxy[k] * tbe_0 * tke_0;

            to_y_y_xxy_xy[k] = to_xx_x[k] - 2.0 * to_xx_xyy[k] * tke_0 - 2.0 * to_xxyy_x[k] * tbe_0 + 4.0 * to_xxyy_xyy[k] * tbe_0 * tke_0;

            to_y_y_xxy_xz[k] = -2.0 * to_xx_xyz[k] * tke_0 + 4.0 * to_xxyy_xyz[k] * tbe_0 * tke_0;

            to_y_y_xxy_yy[k] = 2.0 * to_xx_y[k] - 2.0 * to_xx_yyy[k] * tke_0 - 4.0 * to_xxyy_y[k] * tbe_0 + 4.0 * to_xxyy_yyy[k] * tbe_0 * tke_0;

            to_y_y_xxy_yz[k] = to_xx_z[k] - 2.0 * to_xx_yyz[k] * tke_0 - 2.0 * to_xxyy_z[k] * tbe_0 + 4.0 * to_xxyy_yyz[k] * tbe_0 * tke_0;

            to_y_y_xxy_zz[k] = -2.0 * to_xx_yzz[k] * tke_0 + 4.0 * to_xxyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 252-258 components of targeted buffer : FD

        auto to_y_y_xxz_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 12);

        auto to_y_y_xxz_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 13);

        auto to_y_y_xxz_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 14);

        auto to_y_y_xxz_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 15);

        auto to_y_y_xxz_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 16);

        auto to_y_y_xxz_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_xxyz_x, to_xxyz_xxy, to_xxyz_xyy, to_xxyz_xyz, to_xxyz_y, to_xxyz_yyy, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_z, to_y_y_xxz_xx, to_y_y_xxz_xy, to_y_y_xxz_xz, to_y_y_xxz_yy, to_y_y_xxz_yz, to_y_y_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxz_xx[k] = 4.0 * to_xxyz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xxz_xy[k] = -2.0 * to_xxyz_x[k] * tbe_0 + 4.0 * to_xxyz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xxz_xz[k] = 4.0 * to_xxyz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xxz_yy[k] = -4.0 * to_xxyz_y[k] * tbe_0 + 4.0 * to_xxyz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xxz_yz[k] = -2.0 * to_xxyz_z[k] * tbe_0 + 4.0 * to_xxyz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xxz_zz[k] = 4.0 * to_xxyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 258-264 components of targeted buffer : FD

        auto to_y_y_xyy_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 18);

        auto to_y_y_xyy_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 19);

        auto to_y_y_xyy_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 20);

        auto to_y_y_xyy_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 21);

        auto to_y_y_xyy_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 22);

        auto to_y_y_xyy_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_xy_x, to_xy_xxy, to_xy_xyy, to_xy_xyz, to_xy_y, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_z, to_xyyy_x, to_xyyy_xxy, to_xyyy_xyy, to_xyyy_xyz, to_xyyy_y, to_xyyy_yyy, to_xyyy_yyz, to_xyyy_yzz, to_xyyy_z, to_y_y_xyy_xx, to_y_y_xyy_xy, to_y_y_xyy_xz, to_y_y_xyy_yy, to_y_y_xyy_yz, to_y_y_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyy_xx[k] = -4.0 * to_xy_xxy[k] * tke_0 + 4.0 * to_xyyy_xxy[k] * tbe_0 * tke_0;

            to_y_y_xyy_xy[k] = 2.0 * to_xy_x[k] - 4.0 * to_xy_xyy[k] * tke_0 - 2.0 * to_xyyy_x[k] * tbe_0 + 4.0 * to_xyyy_xyy[k] * tbe_0 * tke_0;

            to_y_y_xyy_xz[k] = -4.0 * to_xy_xyz[k] * tke_0 + 4.0 * to_xyyy_xyz[k] * tbe_0 * tke_0;

            to_y_y_xyy_yy[k] = 4.0 * to_xy_y[k] - 4.0 * to_xy_yyy[k] * tke_0 - 4.0 * to_xyyy_y[k] * tbe_0 + 4.0 * to_xyyy_yyy[k] * tbe_0 * tke_0;

            to_y_y_xyy_yz[k] = 2.0 * to_xy_z[k] - 4.0 * to_xy_yyz[k] * tke_0 - 2.0 * to_xyyy_z[k] * tbe_0 + 4.0 * to_xyyy_yyz[k] * tbe_0 * tke_0;

            to_y_y_xyy_zz[k] = -4.0 * to_xy_yzz[k] * tke_0 + 4.0 * to_xyyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 264-270 components of targeted buffer : FD

        auto to_y_y_xyz_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 24);

        auto to_y_y_xyz_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 25);

        auto to_y_y_xyz_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 26);

        auto to_y_y_xyz_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 27);

        auto to_y_y_xyz_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 28);

        auto to_y_y_xyz_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_xyyz_x, to_xyyz_xxy, to_xyyz_xyy, to_xyyz_xyz, to_xyyz_y, to_xyyz_yyy, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_z, to_xz_x, to_xz_xxy, to_xz_xyy, to_xz_xyz, to_xz_y, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_z, to_y_y_xyz_xx, to_y_y_xyz_xy, to_y_y_xyz_xz, to_y_y_xyz_yy, to_y_y_xyz_yz, to_y_y_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyz_xx[k] = -2.0 * to_xz_xxy[k] * tke_0 + 4.0 * to_xyyz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xyz_xy[k] = to_xz_x[k] - 2.0 * to_xz_xyy[k] * tke_0 - 2.0 * to_xyyz_x[k] * tbe_0 + 4.0 * to_xyyz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xyz_xz[k] = -2.0 * to_xz_xyz[k] * tke_0 + 4.0 * to_xyyz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xyz_yy[k] = 2.0 * to_xz_y[k] - 2.0 * to_xz_yyy[k] * tke_0 - 4.0 * to_xyyz_y[k] * tbe_0 + 4.0 * to_xyyz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xyz_yz[k] = to_xz_z[k] - 2.0 * to_xz_yyz[k] * tke_0 - 2.0 * to_xyyz_z[k] * tbe_0 + 4.0 * to_xyyz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xyz_zz[k] = -2.0 * to_xz_yzz[k] * tke_0 + 4.0 * to_xyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-276 components of targeted buffer : FD

        auto to_y_y_xzz_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 30);

        auto to_y_y_xzz_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 31);

        auto to_y_y_xzz_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 32);

        auto to_y_y_xzz_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 33);

        auto to_y_y_xzz_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 34);

        auto to_y_y_xzz_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_xyzz_x, to_xyzz_xxy, to_xyzz_xyy, to_xyzz_xyz, to_xyzz_y, to_xyzz_yyy, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_z, to_y_y_xzz_xx, to_y_y_xzz_xy, to_y_y_xzz_xz, to_y_y_xzz_yy, to_y_y_xzz_yz, to_y_y_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xzz_xx[k] = 4.0 * to_xyzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xzz_xy[k] = -2.0 * to_xyzz_x[k] * tbe_0 + 4.0 * to_xyzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xzz_xz[k] = 4.0 * to_xyzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xzz_yy[k] = -4.0 * to_xyzz_y[k] * tbe_0 + 4.0 * to_xyzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xzz_yz[k] = -2.0 * to_xyzz_z[k] * tbe_0 + 4.0 * to_xyzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xzz_zz[k] = 4.0 * to_xyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 276-282 components of targeted buffer : FD

        auto to_y_y_yyy_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 36);

        auto to_y_y_yyy_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 37);

        auto to_y_y_yyy_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 38);

        auto to_y_y_yyy_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 39);

        auto to_y_y_yyy_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 40);

        auto to_y_y_yyy_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_y_y_yyy_xx, to_y_y_yyy_xy, to_y_y_yyy_xz, to_y_y_yyy_yy, to_y_y_yyy_yz, to_y_y_yyy_zz, to_yy_x, to_yy_xxy, to_yy_xyy, to_yy_xyz, to_yy_y, to_yy_yyy, to_yy_yyz, to_yy_yzz, to_yy_z, to_yyyy_x, to_yyyy_xxy, to_yyyy_xyy, to_yyyy_xyz, to_yyyy_y, to_yyyy_yyy, to_yyyy_yyz, to_yyyy_yzz, to_yyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyy_xx[k] = -6.0 * to_yy_xxy[k] * tke_0 + 4.0 * to_yyyy_xxy[k] * tbe_0 * tke_0;

            to_y_y_yyy_xy[k] = 3.0 * to_yy_x[k] - 6.0 * to_yy_xyy[k] * tke_0 - 2.0 * to_yyyy_x[k] * tbe_0 + 4.0 * to_yyyy_xyy[k] * tbe_0 * tke_0;

            to_y_y_yyy_xz[k] = -6.0 * to_yy_xyz[k] * tke_0 + 4.0 * to_yyyy_xyz[k] * tbe_0 * tke_0;

            to_y_y_yyy_yy[k] = 6.0 * to_yy_y[k] - 6.0 * to_yy_yyy[k] * tke_0 - 4.0 * to_yyyy_y[k] * tbe_0 + 4.0 * to_yyyy_yyy[k] * tbe_0 * tke_0;

            to_y_y_yyy_yz[k] = 3.0 * to_yy_z[k] - 6.0 * to_yy_yyz[k] * tke_0 - 2.0 * to_yyyy_z[k] * tbe_0 + 4.0 * to_yyyy_yyz[k] * tbe_0 * tke_0;

            to_y_y_yyy_zz[k] = -6.0 * to_yy_yzz[k] * tke_0 + 4.0 * to_yyyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 282-288 components of targeted buffer : FD

        auto to_y_y_yyz_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 42);

        auto to_y_y_yyz_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 43);

        auto to_y_y_yyz_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 44);

        auto to_y_y_yyz_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 45);

        auto to_y_y_yyz_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 46);

        auto to_y_y_yyz_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_y_y_yyz_xx, to_y_y_yyz_xy, to_y_y_yyz_xz, to_y_y_yyz_yy, to_y_y_yyz_yz, to_y_y_yyz_zz, to_yyyz_x, to_yyyz_xxy, to_yyyz_xyy, to_yyyz_xyz, to_yyyz_y, to_yyyz_yyy, to_yyyz_yyz, to_yyyz_yzz, to_yyyz_z, to_yz_x, to_yz_xxy, to_yz_xyy, to_yz_xyz, to_yz_y, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyz_xx[k] = -4.0 * to_yz_xxy[k] * tke_0 + 4.0 * to_yyyz_xxy[k] * tbe_0 * tke_0;

            to_y_y_yyz_xy[k] = 2.0 * to_yz_x[k] - 4.0 * to_yz_xyy[k] * tke_0 - 2.0 * to_yyyz_x[k] * tbe_0 + 4.0 * to_yyyz_xyy[k] * tbe_0 * tke_0;

            to_y_y_yyz_xz[k] = -4.0 * to_yz_xyz[k] * tke_0 + 4.0 * to_yyyz_xyz[k] * tbe_0 * tke_0;

            to_y_y_yyz_yy[k] = 4.0 * to_yz_y[k] - 4.0 * to_yz_yyy[k] * tke_0 - 4.0 * to_yyyz_y[k] * tbe_0 + 4.0 * to_yyyz_yyy[k] * tbe_0 * tke_0;

            to_y_y_yyz_yz[k] = 2.0 * to_yz_z[k] - 4.0 * to_yz_yyz[k] * tke_0 - 2.0 * to_yyyz_z[k] * tbe_0 + 4.0 * to_yyyz_yyz[k] * tbe_0 * tke_0;

            to_y_y_yyz_zz[k] = -4.0 * to_yz_yzz[k] * tke_0 + 4.0 * to_yyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 288-294 components of targeted buffer : FD

        auto to_y_y_yzz_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 48);

        auto to_y_y_yzz_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 49);

        auto to_y_y_yzz_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 50);

        auto to_y_y_yzz_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 51);

        auto to_y_y_yzz_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 52);

        auto to_y_y_yzz_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_y_y_yzz_xx, to_y_y_yzz_xy, to_y_y_yzz_xz, to_y_y_yzz_yy, to_y_y_yzz_yz, to_y_y_yzz_zz, to_yyzz_x, to_yyzz_xxy, to_yyzz_xyy, to_yyzz_xyz, to_yyzz_y, to_yyzz_yyy, to_yyzz_yyz, to_yyzz_yzz, to_yyzz_z, to_zz_x, to_zz_xxy, to_zz_xyy, to_zz_xyz, to_zz_y, to_zz_yyy, to_zz_yyz, to_zz_yzz, to_zz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yzz_xx[k] = -2.0 * to_zz_xxy[k] * tke_0 + 4.0 * to_yyzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_yzz_xy[k] = to_zz_x[k] - 2.0 * to_zz_xyy[k] * tke_0 - 2.0 * to_yyzz_x[k] * tbe_0 + 4.0 * to_yyzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_yzz_xz[k] = -2.0 * to_zz_xyz[k] * tke_0 + 4.0 * to_yyzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_yzz_yy[k] = 2.0 * to_zz_y[k] - 2.0 * to_zz_yyy[k] * tke_0 - 4.0 * to_yyzz_y[k] * tbe_0 + 4.0 * to_yyzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_yzz_yz[k] = to_zz_z[k] - 2.0 * to_zz_yyz[k] * tke_0 - 2.0 * to_yyzz_z[k] * tbe_0 + 4.0 * to_yyzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_yzz_zz[k] = -2.0 * to_zz_yzz[k] * tke_0 + 4.0 * to_yyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 294-300 components of targeted buffer : FD

        auto to_y_y_zzz_xx = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 54);

        auto to_y_y_zzz_xy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 55);

        auto to_y_y_zzz_xz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 56);

        auto to_y_y_zzz_yy = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 57);

        auto to_y_y_zzz_yz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 58);

        auto to_y_y_zzz_zz = pbuffer.data(idx_op_geom_101_fd + 4 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_y_y_zzz_xx, to_y_y_zzz_xy, to_y_y_zzz_xz, to_y_y_zzz_yy, to_y_y_zzz_yz, to_y_y_zzz_zz, to_yzzz_x, to_yzzz_xxy, to_yzzz_xyy, to_yzzz_xyz, to_yzzz_y, to_yzzz_yyy, to_yzzz_yyz, to_yzzz_yzz, to_yzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zzz_xx[k] = 4.0 * to_yzzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_zzz_xy[k] = -2.0 * to_yzzz_x[k] * tbe_0 + 4.0 * to_yzzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_zzz_xz[k] = 4.0 * to_yzzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_zzz_yy[k] = -4.0 * to_yzzz_y[k] * tbe_0 + 4.0 * to_yzzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_zzz_yz[k] = -2.0 * to_yzzz_z[k] * tbe_0 + 4.0 * to_yzzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_zzz_zz[k] = 4.0 * to_yzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-306 components of targeted buffer : FD

        auto to_y_z_xxx_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 0);

        auto to_y_z_xxx_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 1);

        auto to_y_z_xxx_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 2);

        auto to_y_z_xxx_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 3);

        auto to_y_z_xxx_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 4);

        auto to_y_z_xxx_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_xxxy_x, to_xxxy_xxz, to_xxxy_xyz, to_xxxy_xzz, to_xxxy_y, to_xxxy_yyz, to_xxxy_yzz, to_xxxy_z, to_xxxy_zzz, to_y_z_xxx_xx, to_y_z_xxx_xy, to_y_z_xxx_xz, to_y_z_xxx_yy, to_y_z_xxx_yz, to_y_z_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxx_xx[k] = 4.0 * to_xxxy_xxz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xy[k] = 4.0 * to_xxxy_xyz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xz[k] = -2.0 * to_xxxy_x[k] * tbe_0 + 4.0 * to_xxxy_xzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_yy[k] = 4.0 * to_xxxy_yyz[k] * tbe_0 * tke_0;

            to_y_z_xxx_yz[k] = -2.0 * to_xxxy_y[k] * tbe_0 + 4.0 * to_xxxy_yzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_zz[k] = -4.0 * to_xxxy_z[k] * tbe_0 + 4.0 * to_xxxy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 306-312 components of targeted buffer : FD

        auto to_y_z_xxy_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 6);

        auto to_y_z_xxy_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 7);

        auto to_y_z_xxy_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 8);

        auto to_y_z_xxy_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 9);

        auto to_y_z_xxy_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 10);

        auto to_y_z_xxy_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_xx_x, to_xx_xxz, to_xx_xyz, to_xx_xzz, to_xx_y, to_xx_yyz, to_xx_yzz, to_xx_z, to_xx_zzz, to_xxyy_x, to_xxyy_xxz, to_xxyy_xyz, to_xxyy_xzz, to_xxyy_y, to_xxyy_yyz, to_xxyy_yzz, to_xxyy_z, to_xxyy_zzz, to_y_z_xxy_xx, to_y_z_xxy_xy, to_y_z_xxy_xz, to_y_z_xxy_yy, to_y_z_xxy_yz, to_y_z_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxy_xx[k] = -2.0 * to_xx_xxz[k] * tke_0 + 4.0 * to_xxyy_xxz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xy[k] = -2.0 * to_xx_xyz[k] * tke_0 + 4.0 * to_xxyy_xyz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xz[k] = to_xx_x[k] - 2.0 * to_xx_xzz[k] * tke_0 - 2.0 * to_xxyy_x[k] * tbe_0 + 4.0 * to_xxyy_xzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_yy[k] = -2.0 * to_xx_yyz[k] * tke_0 + 4.0 * to_xxyy_yyz[k] * tbe_0 * tke_0;

            to_y_z_xxy_yz[k] = to_xx_y[k] - 2.0 * to_xx_yzz[k] * tke_0 - 2.0 * to_xxyy_y[k] * tbe_0 + 4.0 * to_xxyy_yzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_zz[k] = 2.0 * to_xx_z[k] - 2.0 * to_xx_zzz[k] * tke_0 - 4.0 * to_xxyy_z[k] * tbe_0 + 4.0 * to_xxyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 312-318 components of targeted buffer : FD

        auto to_y_z_xxz_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 12);

        auto to_y_z_xxz_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 13);

        auto to_y_z_xxz_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 14);

        auto to_y_z_xxz_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 15);

        auto to_y_z_xxz_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 16);

        auto to_y_z_xxz_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_xxyz_x, to_xxyz_xxz, to_xxyz_xyz, to_xxyz_xzz, to_xxyz_y, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_z, to_xxyz_zzz, to_y_z_xxz_xx, to_y_z_xxz_xy, to_y_z_xxz_xz, to_y_z_xxz_yy, to_y_z_xxz_yz, to_y_z_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxz_xx[k] = 4.0 * to_xxyz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xy[k] = 4.0 * to_xxyz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xz[k] = -2.0 * to_xxyz_x[k] * tbe_0 + 4.0 * to_xxyz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_yy[k] = 4.0 * to_xxyz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xxz_yz[k] = -2.0 * to_xxyz_y[k] * tbe_0 + 4.0 * to_xxyz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_zz[k] = -4.0 * to_xxyz_z[k] * tbe_0 + 4.0 * to_xxyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 318-324 components of targeted buffer : FD

        auto to_y_z_xyy_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 18);

        auto to_y_z_xyy_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 19);

        auto to_y_z_xyy_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 20);

        auto to_y_z_xyy_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 21);

        auto to_y_z_xyy_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 22);

        auto to_y_z_xyy_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_xy_x, to_xy_xxz, to_xy_xyz, to_xy_xzz, to_xy_y, to_xy_yyz, to_xy_yzz, to_xy_z, to_xy_zzz, to_xyyy_x, to_xyyy_xxz, to_xyyy_xyz, to_xyyy_xzz, to_xyyy_y, to_xyyy_yyz, to_xyyy_yzz, to_xyyy_z, to_xyyy_zzz, to_y_z_xyy_xx, to_y_z_xyy_xy, to_y_z_xyy_xz, to_y_z_xyy_yy, to_y_z_xyy_yz, to_y_z_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyy_xx[k] = -4.0 * to_xy_xxz[k] * tke_0 + 4.0 * to_xyyy_xxz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xy[k] = -4.0 * to_xy_xyz[k] * tke_0 + 4.0 * to_xyyy_xyz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xz[k] = 2.0 * to_xy_x[k] - 4.0 * to_xy_xzz[k] * tke_0 - 2.0 * to_xyyy_x[k] * tbe_0 + 4.0 * to_xyyy_xzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_yy[k] = -4.0 * to_xy_yyz[k] * tke_0 + 4.0 * to_xyyy_yyz[k] * tbe_0 * tke_0;

            to_y_z_xyy_yz[k] = 2.0 * to_xy_y[k] - 4.0 * to_xy_yzz[k] * tke_0 - 2.0 * to_xyyy_y[k] * tbe_0 + 4.0 * to_xyyy_yzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_zz[k] = 4.0 * to_xy_z[k] - 4.0 * to_xy_zzz[k] * tke_0 - 4.0 * to_xyyy_z[k] * tbe_0 + 4.0 * to_xyyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 324-330 components of targeted buffer : FD

        auto to_y_z_xyz_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 24);

        auto to_y_z_xyz_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 25);

        auto to_y_z_xyz_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 26);

        auto to_y_z_xyz_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 27);

        auto to_y_z_xyz_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 28);

        auto to_y_z_xyz_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_xyyz_x, to_xyyz_xxz, to_xyyz_xyz, to_xyyz_xzz, to_xyyz_y, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_z, to_xyyz_zzz, to_xz_x, to_xz_xxz, to_xz_xyz, to_xz_xzz, to_xz_y, to_xz_yyz, to_xz_yzz, to_xz_z, to_xz_zzz, to_y_z_xyz_xx, to_y_z_xyz_xy, to_y_z_xyz_xz, to_y_z_xyz_yy, to_y_z_xyz_yz, to_y_z_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyz_xx[k] = -2.0 * to_xz_xxz[k] * tke_0 + 4.0 * to_xyyz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xy[k] = -2.0 * to_xz_xyz[k] * tke_0 + 4.0 * to_xyyz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xz[k] = to_xz_x[k] - 2.0 * to_xz_xzz[k] * tke_0 - 2.0 * to_xyyz_x[k] * tbe_0 + 4.0 * to_xyyz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_yy[k] = -2.0 * to_xz_yyz[k] * tke_0 + 4.0 * to_xyyz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xyz_yz[k] = to_xz_y[k] - 2.0 * to_xz_yzz[k] * tke_0 - 2.0 * to_xyyz_y[k] * tbe_0 + 4.0 * to_xyyz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_zz[k] = 2.0 * to_xz_z[k] - 2.0 * to_xz_zzz[k] * tke_0 - 4.0 * to_xyyz_z[k] * tbe_0 + 4.0 * to_xyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 330-336 components of targeted buffer : FD

        auto to_y_z_xzz_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 30);

        auto to_y_z_xzz_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 31);

        auto to_y_z_xzz_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 32);

        auto to_y_z_xzz_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 33);

        auto to_y_z_xzz_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 34);

        auto to_y_z_xzz_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_xyzz_x, to_xyzz_xxz, to_xyzz_xyz, to_xyzz_xzz, to_xyzz_y, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_z, to_xyzz_zzz, to_y_z_xzz_xx, to_y_z_xzz_xy, to_y_z_xzz_xz, to_y_z_xzz_yy, to_y_z_xzz_yz, to_y_z_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xzz_xx[k] = 4.0 * to_xyzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xy[k] = 4.0 * to_xyzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xz[k] = -2.0 * to_xyzz_x[k] * tbe_0 + 4.0 * to_xyzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_yy[k] = 4.0 * to_xyzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xzz_yz[k] = -2.0 * to_xyzz_y[k] * tbe_0 + 4.0 * to_xyzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_zz[k] = -4.0 * to_xyzz_z[k] * tbe_0 + 4.0 * to_xyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 336-342 components of targeted buffer : FD

        auto to_y_z_yyy_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 36);

        auto to_y_z_yyy_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 37);

        auto to_y_z_yyy_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 38);

        auto to_y_z_yyy_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 39);

        auto to_y_z_yyy_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 40);

        auto to_y_z_yyy_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_y_z_yyy_xx, to_y_z_yyy_xy, to_y_z_yyy_xz, to_y_z_yyy_yy, to_y_z_yyy_yz, to_y_z_yyy_zz, to_yy_x, to_yy_xxz, to_yy_xyz, to_yy_xzz, to_yy_y, to_yy_yyz, to_yy_yzz, to_yy_z, to_yy_zzz, to_yyyy_x, to_yyyy_xxz, to_yyyy_xyz, to_yyyy_xzz, to_yyyy_y, to_yyyy_yyz, to_yyyy_yzz, to_yyyy_z, to_yyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyy_xx[k] = -6.0 * to_yy_xxz[k] * tke_0 + 4.0 * to_yyyy_xxz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xy[k] = -6.0 * to_yy_xyz[k] * tke_0 + 4.0 * to_yyyy_xyz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xz[k] = 3.0 * to_yy_x[k] - 6.0 * to_yy_xzz[k] * tke_0 - 2.0 * to_yyyy_x[k] * tbe_0 + 4.0 * to_yyyy_xzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_yy[k] = -6.0 * to_yy_yyz[k] * tke_0 + 4.0 * to_yyyy_yyz[k] * tbe_0 * tke_0;

            to_y_z_yyy_yz[k] = 3.0 * to_yy_y[k] - 6.0 * to_yy_yzz[k] * tke_0 - 2.0 * to_yyyy_y[k] * tbe_0 + 4.0 * to_yyyy_yzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_zz[k] = 6.0 * to_yy_z[k] - 6.0 * to_yy_zzz[k] * tke_0 - 4.0 * to_yyyy_z[k] * tbe_0 + 4.0 * to_yyyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 342-348 components of targeted buffer : FD

        auto to_y_z_yyz_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 42);

        auto to_y_z_yyz_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 43);

        auto to_y_z_yyz_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 44);

        auto to_y_z_yyz_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 45);

        auto to_y_z_yyz_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 46);

        auto to_y_z_yyz_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_y_z_yyz_xx, to_y_z_yyz_xy, to_y_z_yyz_xz, to_y_z_yyz_yy, to_y_z_yyz_yz, to_y_z_yyz_zz, to_yyyz_x, to_yyyz_xxz, to_yyyz_xyz, to_yyyz_xzz, to_yyyz_y, to_yyyz_yyz, to_yyyz_yzz, to_yyyz_z, to_yyyz_zzz, to_yz_x, to_yz_xxz, to_yz_xyz, to_yz_xzz, to_yz_y, to_yz_yyz, to_yz_yzz, to_yz_z, to_yz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyz_xx[k] = -4.0 * to_yz_xxz[k] * tke_0 + 4.0 * to_yyyz_xxz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xy[k] = -4.0 * to_yz_xyz[k] * tke_0 + 4.0 * to_yyyz_xyz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xz[k] = 2.0 * to_yz_x[k] - 4.0 * to_yz_xzz[k] * tke_0 - 2.0 * to_yyyz_x[k] * tbe_0 + 4.0 * to_yyyz_xzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_yy[k] = -4.0 * to_yz_yyz[k] * tke_0 + 4.0 * to_yyyz_yyz[k] * tbe_0 * tke_0;

            to_y_z_yyz_yz[k] = 2.0 * to_yz_y[k] - 4.0 * to_yz_yzz[k] * tke_0 - 2.0 * to_yyyz_y[k] * tbe_0 + 4.0 * to_yyyz_yzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_zz[k] = 4.0 * to_yz_z[k] - 4.0 * to_yz_zzz[k] * tke_0 - 4.0 * to_yyyz_z[k] * tbe_0 + 4.0 * to_yyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 348-354 components of targeted buffer : FD

        auto to_y_z_yzz_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 48);

        auto to_y_z_yzz_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 49);

        auto to_y_z_yzz_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 50);

        auto to_y_z_yzz_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 51);

        auto to_y_z_yzz_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 52);

        auto to_y_z_yzz_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_y_z_yzz_xx, to_y_z_yzz_xy, to_y_z_yzz_xz, to_y_z_yzz_yy, to_y_z_yzz_yz, to_y_z_yzz_zz, to_yyzz_x, to_yyzz_xxz, to_yyzz_xyz, to_yyzz_xzz, to_yyzz_y, to_yyzz_yyz, to_yyzz_yzz, to_yyzz_z, to_yyzz_zzz, to_zz_x, to_zz_xxz, to_zz_xyz, to_zz_xzz, to_zz_y, to_zz_yyz, to_zz_yzz, to_zz_z, to_zz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yzz_xx[k] = -2.0 * to_zz_xxz[k] * tke_0 + 4.0 * to_yyzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xy[k] = -2.0 * to_zz_xyz[k] * tke_0 + 4.0 * to_yyzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xz[k] = to_zz_x[k] - 2.0 * to_zz_xzz[k] * tke_0 - 2.0 * to_yyzz_x[k] * tbe_0 + 4.0 * to_yyzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_yy[k] = -2.0 * to_zz_yyz[k] * tke_0 + 4.0 * to_yyzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_yzz_yz[k] = to_zz_y[k] - 2.0 * to_zz_yzz[k] * tke_0 - 2.0 * to_yyzz_y[k] * tbe_0 + 4.0 * to_yyzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_zz[k] = 2.0 * to_zz_z[k] - 2.0 * to_zz_zzz[k] * tke_0 - 4.0 * to_yyzz_z[k] * tbe_0 + 4.0 * to_yyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 354-360 components of targeted buffer : FD

        auto to_y_z_zzz_xx = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 54);

        auto to_y_z_zzz_xy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 55);

        auto to_y_z_zzz_xz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 56);

        auto to_y_z_zzz_yy = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 57);

        auto to_y_z_zzz_yz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 58);

        auto to_y_z_zzz_zz = pbuffer.data(idx_op_geom_101_fd + 5 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_y_z_zzz_xx, to_y_z_zzz_xy, to_y_z_zzz_xz, to_y_z_zzz_yy, to_y_z_zzz_yz, to_y_z_zzz_zz, to_yzzz_x, to_yzzz_xxz, to_yzzz_xyz, to_yzzz_xzz, to_yzzz_y, to_yzzz_yyz, to_yzzz_yzz, to_yzzz_z, to_yzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zzz_xx[k] = 4.0 * to_yzzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xy[k] = 4.0 * to_yzzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xz[k] = -2.0 * to_yzzz_x[k] * tbe_0 + 4.0 * to_yzzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_yy[k] = 4.0 * to_yzzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_zzz_yz[k] = -2.0 * to_yzzz_y[k] * tbe_0 + 4.0 * to_yzzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_zz[k] = -4.0 * to_yzzz_z[k] * tbe_0 + 4.0 * to_yzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 360-366 components of targeted buffer : FD

        auto to_z_x_xxx_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 0);

        auto to_z_x_xxx_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 1);

        auto to_z_x_xxx_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 2);

        auto to_z_x_xxx_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 3);

        auto to_z_x_xxx_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 4);

        auto to_z_x_xxx_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_xxxz_x, to_xxxz_xxx, to_xxxz_xxy, to_xxxz_xxz, to_xxxz_xyy, to_xxxz_xyz, to_xxxz_xzz, to_xxxz_y, to_xxxz_z, to_z_x_xxx_xx, to_z_x_xxx_xy, to_z_x_xxx_xz, to_z_x_xxx_yy, to_z_x_xxx_yz, to_z_x_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxx_xx[k] = -4.0 * to_xxxz_x[k] * tbe_0 + 4.0 * to_xxxz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xxx_xy[k] = -2.0 * to_xxxz_y[k] * tbe_0 + 4.0 * to_xxxz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xxx_xz[k] = -2.0 * to_xxxz_z[k] * tbe_0 + 4.0 * to_xxxz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xxx_yy[k] = 4.0 * to_xxxz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xxx_yz[k] = 4.0 * to_xxxz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xxx_zz[k] = 4.0 * to_xxxz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 366-372 components of targeted buffer : FD

        auto to_z_x_xxy_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 6);

        auto to_z_x_xxy_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 7);

        auto to_z_x_xxy_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 8);

        auto to_z_x_xxy_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 9);

        auto to_z_x_xxy_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 10);

        auto to_z_x_xxy_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_xxyz_x, to_xxyz_xxx, to_xxyz_xxy, to_xxyz_xxz, to_xxyz_xyy, to_xxyz_xyz, to_xxyz_xzz, to_xxyz_y, to_xxyz_z, to_z_x_xxy_xx, to_z_x_xxy_xy, to_z_x_xxy_xz, to_z_x_xxy_yy, to_z_x_xxy_yz, to_z_x_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxy_xx[k] = -4.0 * to_xxyz_x[k] * tbe_0 + 4.0 * to_xxyz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xxy_xy[k] = -2.0 * to_xxyz_y[k] * tbe_0 + 4.0 * to_xxyz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xxy_xz[k] = -2.0 * to_xxyz_z[k] * tbe_0 + 4.0 * to_xxyz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xxy_yy[k] = 4.0 * to_xxyz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xxy_yz[k] = 4.0 * to_xxyz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xxy_zz[k] = 4.0 * to_xxyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 372-378 components of targeted buffer : FD

        auto to_z_x_xxz_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 12);

        auto to_z_x_xxz_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 13);

        auto to_z_x_xxz_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 14);

        auto to_z_x_xxz_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 15);

        auto to_z_x_xxz_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 16);

        auto to_z_x_xxz_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_xx_x, to_xx_xxx, to_xx_xxy, to_xx_xxz, to_xx_xyy, to_xx_xyz, to_xx_xzz, to_xx_y, to_xx_z, to_xxzz_x, to_xxzz_xxx, to_xxzz_xxy, to_xxzz_xxz, to_xxzz_xyy, to_xxzz_xyz, to_xxzz_xzz, to_xxzz_y, to_xxzz_z, to_z_x_xxz_xx, to_z_x_xxz_xy, to_z_x_xxz_xz, to_z_x_xxz_yy, to_z_x_xxz_yz, to_z_x_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxz_xx[k] = 2.0 * to_xx_x[k] - 2.0 * to_xx_xxx[k] * tke_0 - 4.0 * to_xxzz_x[k] * tbe_0 + 4.0 * to_xxzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xxz_xy[k] = to_xx_y[k] - 2.0 * to_xx_xxy[k] * tke_0 - 2.0 * to_xxzz_y[k] * tbe_0 + 4.0 * to_xxzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xxz_xz[k] = to_xx_z[k] - 2.0 * to_xx_xxz[k] * tke_0 - 2.0 * to_xxzz_z[k] * tbe_0 + 4.0 * to_xxzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xxz_yy[k] = -2.0 * to_xx_xyy[k] * tke_0 + 4.0 * to_xxzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xxz_yz[k] = -2.0 * to_xx_xyz[k] * tke_0 + 4.0 * to_xxzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xxz_zz[k] = -2.0 * to_xx_xzz[k] * tke_0 + 4.0 * to_xxzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 378-384 components of targeted buffer : FD

        auto to_z_x_xyy_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 18);

        auto to_z_x_xyy_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 19);

        auto to_z_x_xyy_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 20);

        auto to_z_x_xyy_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 21);

        auto to_z_x_xyy_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 22);

        auto to_z_x_xyy_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_xyyz_x, to_xyyz_xxx, to_xyyz_xxy, to_xyyz_xxz, to_xyyz_xyy, to_xyyz_xyz, to_xyyz_xzz, to_xyyz_y, to_xyyz_z, to_z_x_xyy_xx, to_z_x_xyy_xy, to_z_x_xyy_xz, to_z_x_xyy_yy, to_z_x_xyy_yz, to_z_x_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyy_xx[k] = -4.0 * to_xyyz_x[k] * tbe_0 + 4.0 * to_xyyz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xyy_xy[k] = -2.0 * to_xyyz_y[k] * tbe_0 + 4.0 * to_xyyz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xyy_xz[k] = -2.0 * to_xyyz_z[k] * tbe_0 + 4.0 * to_xyyz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xyy_yy[k] = 4.0 * to_xyyz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xyy_yz[k] = 4.0 * to_xyyz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xyy_zz[k] = 4.0 * to_xyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 384-390 components of targeted buffer : FD

        auto to_z_x_xyz_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 24);

        auto to_z_x_xyz_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 25);

        auto to_z_x_xyz_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 26);

        auto to_z_x_xyz_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 27);

        auto to_z_x_xyz_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 28);

        auto to_z_x_xyz_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_xy_x, to_xy_xxx, to_xy_xxy, to_xy_xxz, to_xy_xyy, to_xy_xyz, to_xy_xzz, to_xy_y, to_xy_z, to_xyzz_x, to_xyzz_xxx, to_xyzz_xxy, to_xyzz_xxz, to_xyzz_xyy, to_xyzz_xyz, to_xyzz_xzz, to_xyzz_y, to_xyzz_z, to_z_x_xyz_xx, to_z_x_xyz_xy, to_z_x_xyz_xz, to_z_x_xyz_yy, to_z_x_xyz_yz, to_z_x_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyz_xx[k] = 2.0 * to_xy_x[k] - 2.0 * to_xy_xxx[k] * tke_0 - 4.0 * to_xyzz_x[k] * tbe_0 + 4.0 * to_xyzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xyz_xy[k] = to_xy_y[k] - 2.0 * to_xy_xxy[k] * tke_0 - 2.0 * to_xyzz_y[k] * tbe_0 + 4.0 * to_xyzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xyz_xz[k] = to_xy_z[k] - 2.0 * to_xy_xxz[k] * tke_0 - 2.0 * to_xyzz_z[k] * tbe_0 + 4.0 * to_xyzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xyz_yy[k] = -2.0 * to_xy_xyy[k] * tke_0 + 4.0 * to_xyzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xyz_yz[k] = -2.0 * to_xy_xyz[k] * tke_0 + 4.0 * to_xyzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xyz_zz[k] = -2.0 * to_xy_xzz[k] * tke_0 + 4.0 * to_xyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 390-396 components of targeted buffer : FD

        auto to_z_x_xzz_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 30);

        auto to_z_x_xzz_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 31);

        auto to_z_x_xzz_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 32);

        auto to_z_x_xzz_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 33);

        auto to_z_x_xzz_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 34);

        auto to_z_x_xzz_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_xz_x, to_xz_xxx, to_xz_xxy, to_xz_xxz, to_xz_xyy, to_xz_xyz, to_xz_xzz, to_xz_y, to_xz_z, to_xzzz_x, to_xzzz_xxx, to_xzzz_xxy, to_xzzz_xxz, to_xzzz_xyy, to_xzzz_xyz, to_xzzz_xzz, to_xzzz_y, to_xzzz_z, to_z_x_xzz_xx, to_z_x_xzz_xy, to_z_x_xzz_xz, to_z_x_xzz_yy, to_z_x_xzz_yz, to_z_x_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xzz_xx[k] = 4.0 * to_xz_x[k] - 4.0 * to_xz_xxx[k] * tke_0 - 4.0 * to_xzzz_x[k] * tbe_0 + 4.0 * to_xzzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xzz_xy[k] = 2.0 * to_xz_y[k] - 4.0 * to_xz_xxy[k] * tke_0 - 2.0 * to_xzzz_y[k] * tbe_0 + 4.0 * to_xzzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xzz_xz[k] = 2.0 * to_xz_z[k] - 4.0 * to_xz_xxz[k] * tke_0 - 2.0 * to_xzzz_z[k] * tbe_0 + 4.0 * to_xzzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xzz_yy[k] = -4.0 * to_xz_xyy[k] * tke_0 + 4.0 * to_xzzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xzz_yz[k] = -4.0 * to_xz_xyz[k] * tke_0 + 4.0 * to_xzzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xzz_zz[k] = -4.0 * to_xz_xzz[k] * tke_0 + 4.0 * to_xzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 396-402 components of targeted buffer : FD

        auto to_z_x_yyy_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 36);

        auto to_z_x_yyy_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 37);

        auto to_z_x_yyy_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 38);

        auto to_z_x_yyy_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 39);

        auto to_z_x_yyy_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 40);

        auto to_z_x_yyy_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_yyyz_x, to_yyyz_xxx, to_yyyz_xxy, to_yyyz_xxz, to_yyyz_xyy, to_yyyz_xyz, to_yyyz_xzz, to_yyyz_y, to_yyyz_z, to_z_x_yyy_xx, to_z_x_yyy_xy, to_z_x_yyy_xz, to_z_x_yyy_yy, to_z_x_yyy_yz, to_z_x_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyy_xx[k] = -4.0 * to_yyyz_x[k] * tbe_0 + 4.0 * to_yyyz_xxx[k] * tbe_0 * tke_0;

            to_z_x_yyy_xy[k] = -2.0 * to_yyyz_y[k] * tbe_0 + 4.0 * to_yyyz_xxy[k] * tbe_0 * tke_0;

            to_z_x_yyy_xz[k] = -2.0 * to_yyyz_z[k] * tbe_0 + 4.0 * to_yyyz_xxz[k] * tbe_0 * tke_0;

            to_z_x_yyy_yy[k] = 4.0 * to_yyyz_xyy[k] * tbe_0 * tke_0;

            to_z_x_yyy_yz[k] = 4.0 * to_yyyz_xyz[k] * tbe_0 * tke_0;

            to_z_x_yyy_zz[k] = 4.0 * to_yyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 402-408 components of targeted buffer : FD

        auto to_z_x_yyz_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 42);

        auto to_z_x_yyz_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 43);

        auto to_z_x_yyz_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 44);

        auto to_z_x_yyz_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 45);

        auto to_z_x_yyz_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 46);

        auto to_z_x_yyz_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_yy_x, to_yy_xxx, to_yy_xxy, to_yy_xxz, to_yy_xyy, to_yy_xyz, to_yy_xzz, to_yy_y, to_yy_z, to_yyzz_x, to_yyzz_xxx, to_yyzz_xxy, to_yyzz_xxz, to_yyzz_xyy, to_yyzz_xyz, to_yyzz_xzz, to_yyzz_y, to_yyzz_z, to_z_x_yyz_xx, to_z_x_yyz_xy, to_z_x_yyz_xz, to_z_x_yyz_yy, to_z_x_yyz_yz, to_z_x_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyz_xx[k] = 2.0 * to_yy_x[k] - 2.0 * to_yy_xxx[k] * tke_0 - 4.0 * to_yyzz_x[k] * tbe_0 + 4.0 * to_yyzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_yyz_xy[k] = to_yy_y[k] - 2.0 * to_yy_xxy[k] * tke_0 - 2.0 * to_yyzz_y[k] * tbe_0 + 4.0 * to_yyzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_yyz_xz[k] = to_yy_z[k] - 2.0 * to_yy_xxz[k] * tke_0 - 2.0 * to_yyzz_z[k] * tbe_0 + 4.0 * to_yyzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_yyz_yy[k] = -2.0 * to_yy_xyy[k] * tke_0 + 4.0 * to_yyzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_yyz_yz[k] = -2.0 * to_yy_xyz[k] * tke_0 + 4.0 * to_yyzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_yyz_zz[k] = -2.0 * to_yy_xzz[k] * tke_0 + 4.0 * to_yyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 408-414 components of targeted buffer : FD

        auto to_z_x_yzz_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 48);

        auto to_z_x_yzz_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 49);

        auto to_z_x_yzz_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 50);

        auto to_z_x_yzz_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 51);

        auto to_z_x_yzz_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 52);

        auto to_z_x_yzz_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_yz_x, to_yz_xxx, to_yz_xxy, to_yz_xxz, to_yz_xyy, to_yz_xyz, to_yz_xzz, to_yz_y, to_yz_z, to_yzzz_x, to_yzzz_xxx, to_yzzz_xxy, to_yzzz_xxz, to_yzzz_xyy, to_yzzz_xyz, to_yzzz_xzz, to_yzzz_y, to_yzzz_z, to_z_x_yzz_xx, to_z_x_yzz_xy, to_z_x_yzz_xz, to_z_x_yzz_yy, to_z_x_yzz_yz, to_z_x_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yzz_xx[k] = 4.0 * to_yz_x[k] - 4.0 * to_yz_xxx[k] * tke_0 - 4.0 * to_yzzz_x[k] * tbe_0 + 4.0 * to_yzzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_yzz_xy[k] = 2.0 * to_yz_y[k] - 4.0 * to_yz_xxy[k] * tke_0 - 2.0 * to_yzzz_y[k] * tbe_0 + 4.0 * to_yzzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_yzz_xz[k] = 2.0 * to_yz_z[k] - 4.0 * to_yz_xxz[k] * tke_0 - 2.0 * to_yzzz_z[k] * tbe_0 + 4.0 * to_yzzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_yzz_yy[k] = -4.0 * to_yz_xyy[k] * tke_0 + 4.0 * to_yzzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_yzz_yz[k] = -4.0 * to_yz_xyz[k] * tke_0 + 4.0 * to_yzzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_yzz_zz[k] = -4.0 * to_yz_xzz[k] * tke_0 + 4.0 * to_yzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 414-420 components of targeted buffer : FD

        auto to_z_x_zzz_xx = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 54);

        auto to_z_x_zzz_xy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 55);

        auto to_z_x_zzz_xz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 56);

        auto to_z_x_zzz_yy = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 57);

        auto to_z_x_zzz_yz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 58);

        auto to_z_x_zzz_zz = pbuffer.data(idx_op_geom_101_fd + 6 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_z_x_zzz_xx, to_z_x_zzz_xy, to_z_x_zzz_xz, to_z_x_zzz_yy, to_z_x_zzz_yz, to_z_x_zzz_zz, to_zz_x, to_zz_xxx, to_zz_xxy, to_zz_xxz, to_zz_xyy, to_zz_xyz, to_zz_xzz, to_zz_y, to_zz_z, to_zzzz_x, to_zzzz_xxx, to_zzzz_xxy, to_zzzz_xxz, to_zzzz_xyy, to_zzzz_xyz, to_zzzz_xzz, to_zzzz_y, to_zzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zzz_xx[k] = 6.0 * to_zz_x[k] - 6.0 * to_zz_xxx[k] * tke_0 - 4.0 * to_zzzz_x[k] * tbe_0 + 4.0 * to_zzzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_zzz_xy[k] = 3.0 * to_zz_y[k] - 6.0 * to_zz_xxy[k] * tke_0 - 2.0 * to_zzzz_y[k] * tbe_0 + 4.0 * to_zzzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_zzz_xz[k] = 3.0 * to_zz_z[k] - 6.0 * to_zz_xxz[k] * tke_0 - 2.0 * to_zzzz_z[k] * tbe_0 + 4.0 * to_zzzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_zzz_yy[k] = -6.0 * to_zz_xyy[k] * tke_0 + 4.0 * to_zzzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_zzz_yz[k] = -6.0 * to_zz_xyz[k] * tke_0 + 4.0 * to_zzzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_zzz_zz[k] = -6.0 * to_zz_xzz[k] * tke_0 + 4.0 * to_zzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 420-426 components of targeted buffer : FD

        auto to_z_y_xxx_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 0);

        auto to_z_y_xxx_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 1);

        auto to_z_y_xxx_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 2);

        auto to_z_y_xxx_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 3);

        auto to_z_y_xxx_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 4);

        auto to_z_y_xxx_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_xxxz_x, to_xxxz_xxy, to_xxxz_xyy, to_xxxz_xyz, to_xxxz_y, to_xxxz_yyy, to_xxxz_yyz, to_xxxz_yzz, to_xxxz_z, to_z_y_xxx_xx, to_z_y_xxx_xy, to_z_y_xxx_xz, to_z_y_xxx_yy, to_z_y_xxx_yz, to_z_y_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxx_xx[k] = 4.0 * to_xxxz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xxx_xy[k] = -2.0 * to_xxxz_x[k] * tbe_0 + 4.0 * to_xxxz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xxx_xz[k] = 4.0 * to_xxxz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xxx_yy[k] = -4.0 * to_xxxz_y[k] * tbe_0 + 4.0 * to_xxxz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xxx_yz[k] = -2.0 * to_xxxz_z[k] * tbe_0 + 4.0 * to_xxxz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xxx_zz[k] = 4.0 * to_xxxz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 426-432 components of targeted buffer : FD

        auto to_z_y_xxy_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 6);

        auto to_z_y_xxy_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 7);

        auto to_z_y_xxy_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 8);

        auto to_z_y_xxy_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 9);

        auto to_z_y_xxy_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 10);

        auto to_z_y_xxy_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_xxyz_x, to_xxyz_xxy, to_xxyz_xyy, to_xxyz_xyz, to_xxyz_y, to_xxyz_yyy, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_z, to_z_y_xxy_xx, to_z_y_xxy_xy, to_z_y_xxy_xz, to_z_y_xxy_yy, to_z_y_xxy_yz, to_z_y_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxy_xx[k] = 4.0 * to_xxyz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xxy_xy[k] = -2.0 * to_xxyz_x[k] * tbe_0 + 4.0 * to_xxyz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xxy_xz[k] = 4.0 * to_xxyz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xxy_yy[k] = -4.0 * to_xxyz_y[k] * tbe_0 + 4.0 * to_xxyz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xxy_yz[k] = -2.0 * to_xxyz_z[k] * tbe_0 + 4.0 * to_xxyz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xxy_zz[k] = 4.0 * to_xxyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 432-438 components of targeted buffer : FD

        auto to_z_y_xxz_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 12);

        auto to_z_y_xxz_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 13);

        auto to_z_y_xxz_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 14);

        auto to_z_y_xxz_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 15);

        auto to_z_y_xxz_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 16);

        auto to_z_y_xxz_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_xx_x, to_xx_xxy, to_xx_xyy, to_xx_xyz, to_xx_y, to_xx_yyy, to_xx_yyz, to_xx_yzz, to_xx_z, to_xxzz_x, to_xxzz_xxy, to_xxzz_xyy, to_xxzz_xyz, to_xxzz_y, to_xxzz_yyy, to_xxzz_yyz, to_xxzz_yzz, to_xxzz_z, to_z_y_xxz_xx, to_z_y_xxz_xy, to_z_y_xxz_xz, to_z_y_xxz_yy, to_z_y_xxz_yz, to_z_y_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxz_xx[k] = -2.0 * to_xx_xxy[k] * tke_0 + 4.0 * to_xxzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xxz_xy[k] = to_xx_x[k] - 2.0 * to_xx_xyy[k] * tke_0 - 2.0 * to_xxzz_x[k] * tbe_0 + 4.0 * to_xxzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xxz_xz[k] = -2.0 * to_xx_xyz[k] * tke_0 + 4.0 * to_xxzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xxz_yy[k] = 2.0 * to_xx_y[k] - 2.0 * to_xx_yyy[k] * tke_0 - 4.0 * to_xxzz_y[k] * tbe_0 + 4.0 * to_xxzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xxz_yz[k] = to_xx_z[k] - 2.0 * to_xx_yyz[k] * tke_0 - 2.0 * to_xxzz_z[k] * tbe_0 + 4.0 * to_xxzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xxz_zz[k] = -2.0 * to_xx_yzz[k] * tke_0 + 4.0 * to_xxzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 438-444 components of targeted buffer : FD

        auto to_z_y_xyy_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 18);

        auto to_z_y_xyy_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 19);

        auto to_z_y_xyy_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 20);

        auto to_z_y_xyy_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 21);

        auto to_z_y_xyy_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 22);

        auto to_z_y_xyy_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_xyyz_x, to_xyyz_xxy, to_xyyz_xyy, to_xyyz_xyz, to_xyyz_y, to_xyyz_yyy, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_z, to_z_y_xyy_xx, to_z_y_xyy_xy, to_z_y_xyy_xz, to_z_y_xyy_yy, to_z_y_xyy_yz, to_z_y_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyy_xx[k] = 4.0 * to_xyyz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xyy_xy[k] = -2.0 * to_xyyz_x[k] * tbe_0 + 4.0 * to_xyyz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xyy_xz[k] = 4.0 * to_xyyz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xyy_yy[k] = -4.0 * to_xyyz_y[k] * tbe_0 + 4.0 * to_xyyz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xyy_yz[k] = -2.0 * to_xyyz_z[k] * tbe_0 + 4.0 * to_xyyz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xyy_zz[k] = 4.0 * to_xyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 444-450 components of targeted buffer : FD

        auto to_z_y_xyz_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 24);

        auto to_z_y_xyz_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 25);

        auto to_z_y_xyz_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 26);

        auto to_z_y_xyz_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 27);

        auto to_z_y_xyz_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 28);

        auto to_z_y_xyz_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_xy_x, to_xy_xxy, to_xy_xyy, to_xy_xyz, to_xy_y, to_xy_yyy, to_xy_yyz, to_xy_yzz, to_xy_z, to_xyzz_x, to_xyzz_xxy, to_xyzz_xyy, to_xyzz_xyz, to_xyzz_y, to_xyzz_yyy, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_z, to_z_y_xyz_xx, to_z_y_xyz_xy, to_z_y_xyz_xz, to_z_y_xyz_yy, to_z_y_xyz_yz, to_z_y_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyz_xx[k] = -2.0 * to_xy_xxy[k] * tke_0 + 4.0 * to_xyzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xyz_xy[k] = to_xy_x[k] - 2.0 * to_xy_xyy[k] * tke_0 - 2.0 * to_xyzz_x[k] * tbe_0 + 4.0 * to_xyzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xyz_xz[k] = -2.0 * to_xy_xyz[k] * tke_0 + 4.0 * to_xyzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xyz_yy[k] = 2.0 * to_xy_y[k] - 2.0 * to_xy_yyy[k] * tke_0 - 4.0 * to_xyzz_y[k] * tbe_0 + 4.0 * to_xyzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xyz_yz[k] = to_xy_z[k] - 2.0 * to_xy_yyz[k] * tke_0 - 2.0 * to_xyzz_z[k] * tbe_0 + 4.0 * to_xyzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xyz_zz[k] = -2.0 * to_xy_yzz[k] * tke_0 + 4.0 * to_xyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 450-456 components of targeted buffer : FD

        auto to_z_y_xzz_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 30);

        auto to_z_y_xzz_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 31);

        auto to_z_y_xzz_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 32);

        auto to_z_y_xzz_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 33);

        auto to_z_y_xzz_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 34);

        auto to_z_y_xzz_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_xz_x, to_xz_xxy, to_xz_xyy, to_xz_xyz, to_xz_y, to_xz_yyy, to_xz_yyz, to_xz_yzz, to_xz_z, to_xzzz_x, to_xzzz_xxy, to_xzzz_xyy, to_xzzz_xyz, to_xzzz_y, to_xzzz_yyy, to_xzzz_yyz, to_xzzz_yzz, to_xzzz_z, to_z_y_xzz_xx, to_z_y_xzz_xy, to_z_y_xzz_xz, to_z_y_xzz_yy, to_z_y_xzz_yz, to_z_y_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xzz_xx[k] = -4.0 * to_xz_xxy[k] * tke_0 + 4.0 * to_xzzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xzz_xy[k] = 2.0 * to_xz_x[k] - 4.0 * to_xz_xyy[k] * tke_0 - 2.0 * to_xzzz_x[k] * tbe_0 + 4.0 * to_xzzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xzz_xz[k] = -4.0 * to_xz_xyz[k] * tke_0 + 4.0 * to_xzzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xzz_yy[k] = 4.0 * to_xz_y[k] - 4.0 * to_xz_yyy[k] * tke_0 - 4.0 * to_xzzz_y[k] * tbe_0 + 4.0 * to_xzzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xzz_yz[k] = 2.0 * to_xz_z[k] - 4.0 * to_xz_yyz[k] * tke_0 - 2.0 * to_xzzz_z[k] * tbe_0 + 4.0 * to_xzzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xzz_zz[k] = -4.0 * to_xz_yzz[k] * tke_0 + 4.0 * to_xzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 456-462 components of targeted buffer : FD

        auto to_z_y_yyy_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 36);

        auto to_z_y_yyy_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 37);

        auto to_z_y_yyy_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 38);

        auto to_z_y_yyy_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 39);

        auto to_z_y_yyy_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 40);

        auto to_z_y_yyy_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_yyyz_x, to_yyyz_xxy, to_yyyz_xyy, to_yyyz_xyz, to_yyyz_y, to_yyyz_yyy, to_yyyz_yyz, to_yyyz_yzz, to_yyyz_z, to_z_y_yyy_xx, to_z_y_yyy_xy, to_z_y_yyy_xz, to_z_y_yyy_yy, to_z_y_yyy_yz, to_z_y_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyy_xx[k] = 4.0 * to_yyyz_xxy[k] * tbe_0 * tke_0;

            to_z_y_yyy_xy[k] = -2.0 * to_yyyz_x[k] * tbe_0 + 4.0 * to_yyyz_xyy[k] * tbe_0 * tke_0;

            to_z_y_yyy_xz[k] = 4.0 * to_yyyz_xyz[k] * tbe_0 * tke_0;

            to_z_y_yyy_yy[k] = -4.0 * to_yyyz_y[k] * tbe_0 + 4.0 * to_yyyz_yyy[k] * tbe_0 * tke_0;

            to_z_y_yyy_yz[k] = -2.0 * to_yyyz_z[k] * tbe_0 + 4.0 * to_yyyz_yyz[k] * tbe_0 * tke_0;

            to_z_y_yyy_zz[k] = 4.0 * to_yyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 462-468 components of targeted buffer : FD

        auto to_z_y_yyz_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 42);

        auto to_z_y_yyz_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 43);

        auto to_z_y_yyz_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 44);

        auto to_z_y_yyz_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 45);

        auto to_z_y_yyz_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 46);

        auto to_z_y_yyz_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_yy_x, to_yy_xxy, to_yy_xyy, to_yy_xyz, to_yy_y, to_yy_yyy, to_yy_yyz, to_yy_yzz, to_yy_z, to_yyzz_x, to_yyzz_xxy, to_yyzz_xyy, to_yyzz_xyz, to_yyzz_y, to_yyzz_yyy, to_yyzz_yyz, to_yyzz_yzz, to_yyzz_z, to_z_y_yyz_xx, to_z_y_yyz_xy, to_z_y_yyz_xz, to_z_y_yyz_yy, to_z_y_yyz_yz, to_z_y_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyz_xx[k] = -2.0 * to_yy_xxy[k] * tke_0 + 4.0 * to_yyzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_yyz_xy[k] = to_yy_x[k] - 2.0 * to_yy_xyy[k] * tke_0 - 2.0 * to_yyzz_x[k] * tbe_0 + 4.0 * to_yyzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_yyz_xz[k] = -2.0 * to_yy_xyz[k] * tke_0 + 4.0 * to_yyzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_yyz_yy[k] = 2.0 * to_yy_y[k] - 2.0 * to_yy_yyy[k] * tke_0 - 4.0 * to_yyzz_y[k] * tbe_0 + 4.0 * to_yyzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_yyz_yz[k] = to_yy_z[k] - 2.0 * to_yy_yyz[k] * tke_0 - 2.0 * to_yyzz_z[k] * tbe_0 + 4.0 * to_yyzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_yyz_zz[k] = -2.0 * to_yy_yzz[k] * tke_0 + 4.0 * to_yyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 468-474 components of targeted buffer : FD

        auto to_z_y_yzz_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 48);

        auto to_z_y_yzz_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 49);

        auto to_z_y_yzz_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 50);

        auto to_z_y_yzz_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 51);

        auto to_z_y_yzz_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 52);

        auto to_z_y_yzz_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_yz_x, to_yz_xxy, to_yz_xyy, to_yz_xyz, to_yz_y, to_yz_yyy, to_yz_yyz, to_yz_yzz, to_yz_z, to_yzzz_x, to_yzzz_xxy, to_yzzz_xyy, to_yzzz_xyz, to_yzzz_y, to_yzzz_yyy, to_yzzz_yyz, to_yzzz_yzz, to_yzzz_z, to_z_y_yzz_xx, to_z_y_yzz_xy, to_z_y_yzz_xz, to_z_y_yzz_yy, to_z_y_yzz_yz, to_z_y_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yzz_xx[k] = -4.0 * to_yz_xxy[k] * tke_0 + 4.0 * to_yzzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_yzz_xy[k] = 2.0 * to_yz_x[k] - 4.0 * to_yz_xyy[k] * tke_0 - 2.0 * to_yzzz_x[k] * tbe_0 + 4.0 * to_yzzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_yzz_xz[k] = -4.0 * to_yz_xyz[k] * tke_0 + 4.0 * to_yzzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_yzz_yy[k] = 4.0 * to_yz_y[k] - 4.0 * to_yz_yyy[k] * tke_0 - 4.0 * to_yzzz_y[k] * tbe_0 + 4.0 * to_yzzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_yzz_yz[k] = 2.0 * to_yz_z[k] - 4.0 * to_yz_yyz[k] * tke_0 - 2.0 * to_yzzz_z[k] * tbe_0 + 4.0 * to_yzzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_yzz_zz[k] = -4.0 * to_yz_yzz[k] * tke_0 + 4.0 * to_yzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 474-480 components of targeted buffer : FD

        auto to_z_y_zzz_xx = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 54);

        auto to_z_y_zzz_xy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 55);

        auto to_z_y_zzz_xz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 56);

        auto to_z_y_zzz_yy = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 57);

        auto to_z_y_zzz_yz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 58);

        auto to_z_y_zzz_zz = pbuffer.data(idx_op_geom_101_fd + 7 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_z_y_zzz_xx, to_z_y_zzz_xy, to_z_y_zzz_xz, to_z_y_zzz_yy, to_z_y_zzz_yz, to_z_y_zzz_zz, to_zz_x, to_zz_xxy, to_zz_xyy, to_zz_xyz, to_zz_y, to_zz_yyy, to_zz_yyz, to_zz_yzz, to_zz_z, to_zzzz_x, to_zzzz_xxy, to_zzzz_xyy, to_zzzz_xyz, to_zzzz_y, to_zzzz_yyy, to_zzzz_yyz, to_zzzz_yzz, to_zzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zzz_xx[k] = -6.0 * to_zz_xxy[k] * tke_0 + 4.0 * to_zzzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_zzz_xy[k] = 3.0 * to_zz_x[k] - 6.0 * to_zz_xyy[k] * tke_0 - 2.0 * to_zzzz_x[k] * tbe_0 + 4.0 * to_zzzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_zzz_xz[k] = -6.0 * to_zz_xyz[k] * tke_0 + 4.0 * to_zzzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_zzz_yy[k] = 6.0 * to_zz_y[k] - 6.0 * to_zz_yyy[k] * tke_0 - 4.0 * to_zzzz_y[k] * tbe_0 + 4.0 * to_zzzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_zzz_yz[k] = 3.0 * to_zz_z[k] - 6.0 * to_zz_yyz[k] * tke_0 - 2.0 * to_zzzz_z[k] * tbe_0 + 4.0 * to_zzzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_zzz_zz[k] = -6.0 * to_zz_yzz[k] * tke_0 + 4.0 * to_zzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 480-486 components of targeted buffer : FD

        auto to_z_z_xxx_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 0);

        auto to_z_z_xxx_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 1);

        auto to_z_z_xxx_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 2);

        auto to_z_z_xxx_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 3);

        auto to_z_z_xxx_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 4);

        auto to_z_z_xxx_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 5);

        #pragma omp simd aligned(to_xxxz_x, to_xxxz_xxz, to_xxxz_xyz, to_xxxz_xzz, to_xxxz_y, to_xxxz_yyz, to_xxxz_yzz, to_xxxz_z, to_xxxz_zzz, to_z_z_xxx_xx, to_z_z_xxx_xy, to_z_z_xxx_xz, to_z_z_xxx_yy, to_z_z_xxx_yz, to_z_z_xxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxx_xx[k] = 4.0 * to_xxxz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xy[k] = 4.0 * to_xxxz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xz[k] = -2.0 * to_xxxz_x[k] * tbe_0 + 4.0 * to_xxxz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_yy[k] = 4.0 * to_xxxz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xxx_yz[k] = -2.0 * to_xxxz_y[k] * tbe_0 + 4.0 * to_xxxz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_zz[k] = -4.0 * to_xxxz_z[k] * tbe_0 + 4.0 * to_xxxz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 486-492 components of targeted buffer : FD

        auto to_z_z_xxy_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 6);

        auto to_z_z_xxy_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 7);

        auto to_z_z_xxy_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 8);

        auto to_z_z_xxy_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 9);

        auto to_z_z_xxy_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 10);

        auto to_z_z_xxy_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 11);

        #pragma omp simd aligned(to_xxyz_x, to_xxyz_xxz, to_xxyz_xyz, to_xxyz_xzz, to_xxyz_y, to_xxyz_yyz, to_xxyz_yzz, to_xxyz_z, to_xxyz_zzz, to_z_z_xxy_xx, to_z_z_xxy_xy, to_z_z_xxy_xz, to_z_z_xxy_yy, to_z_z_xxy_yz, to_z_z_xxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxy_xx[k] = 4.0 * to_xxyz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xy[k] = 4.0 * to_xxyz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xz[k] = -2.0 * to_xxyz_x[k] * tbe_0 + 4.0 * to_xxyz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_yy[k] = 4.0 * to_xxyz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xxy_yz[k] = -2.0 * to_xxyz_y[k] * tbe_0 + 4.0 * to_xxyz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_zz[k] = -4.0 * to_xxyz_z[k] * tbe_0 + 4.0 * to_xxyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 492-498 components of targeted buffer : FD

        auto to_z_z_xxz_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 12);

        auto to_z_z_xxz_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 13);

        auto to_z_z_xxz_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 14);

        auto to_z_z_xxz_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 15);

        auto to_z_z_xxz_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 16);

        auto to_z_z_xxz_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 17);

        #pragma omp simd aligned(to_xx_x, to_xx_xxz, to_xx_xyz, to_xx_xzz, to_xx_y, to_xx_yyz, to_xx_yzz, to_xx_z, to_xx_zzz, to_xxzz_x, to_xxzz_xxz, to_xxzz_xyz, to_xxzz_xzz, to_xxzz_y, to_xxzz_yyz, to_xxzz_yzz, to_xxzz_z, to_xxzz_zzz, to_z_z_xxz_xx, to_z_z_xxz_xy, to_z_z_xxz_xz, to_z_z_xxz_yy, to_z_z_xxz_yz, to_z_z_xxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxz_xx[k] = -2.0 * to_xx_xxz[k] * tke_0 + 4.0 * to_xxzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xy[k] = -2.0 * to_xx_xyz[k] * tke_0 + 4.0 * to_xxzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xz[k] = to_xx_x[k] - 2.0 * to_xx_xzz[k] * tke_0 - 2.0 * to_xxzz_x[k] * tbe_0 + 4.0 * to_xxzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_yy[k] = -2.0 * to_xx_yyz[k] * tke_0 + 4.0 * to_xxzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xxz_yz[k] = to_xx_y[k] - 2.0 * to_xx_yzz[k] * tke_0 - 2.0 * to_xxzz_y[k] * tbe_0 + 4.0 * to_xxzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_zz[k] = 2.0 * to_xx_z[k] - 2.0 * to_xx_zzz[k] * tke_0 - 4.0 * to_xxzz_z[k] * tbe_0 + 4.0 * to_xxzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 498-504 components of targeted buffer : FD

        auto to_z_z_xyy_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 18);

        auto to_z_z_xyy_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 19);

        auto to_z_z_xyy_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 20);

        auto to_z_z_xyy_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 21);

        auto to_z_z_xyy_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 22);

        auto to_z_z_xyy_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 23);

        #pragma omp simd aligned(to_xyyz_x, to_xyyz_xxz, to_xyyz_xyz, to_xyyz_xzz, to_xyyz_y, to_xyyz_yyz, to_xyyz_yzz, to_xyyz_z, to_xyyz_zzz, to_z_z_xyy_xx, to_z_z_xyy_xy, to_z_z_xyy_xz, to_z_z_xyy_yy, to_z_z_xyy_yz, to_z_z_xyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyy_xx[k] = 4.0 * to_xyyz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xy[k] = 4.0 * to_xyyz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xz[k] = -2.0 * to_xyyz_x[k] * tbe_0 + 4.0 * to_xyyz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_yy[k] = 4.0 * to_xyyz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xyy_yz[k] = -2.0 * to_xyyz_y[k] * tbe_0 + 4.0 * to_xyyz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_zz[k] = -4.0 * to_xyyz_z[k] * tbe_0 + 4.0 * to_xyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 504-510 components of targeted buffer : FD

        auto to_z_z_xyz_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 24);

        auto to_z_z_xyz_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 25);

        auto to_z_z_xyz_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 26);

        auto to_z_z_xyz_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 27);

        auto to_z_z_xyz_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 28);

        auto to_z_z_xyz_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 29);

        #pragma omp simd aligned(to_xy_x, to_xy_xxz, to_xy_xyz, to_xy_xzz, to_xy_y, to_xy_yyz, to_xy_yzz, to_xy_z, to_xy_zzz, to_xyzz_x, to_xyzz_xxz, to_xyzz_xyz, to_xyzz_xzz, to_xyzz_y, to_xyzz_yyz, to_xyzz_yzz, to_xyzz_z, to_xyzz_zzz, to_z_z_xyz_xx, to_z_z_xyz_xy, to_z_z_xyz_xz, to_z_z_xyz_yy, to_z_z_xyz_yz, to_z_z_xyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyz_xx[k] = -2.0 * to_xy_xxz[k] * tke_0 + 4.0 * to_xyzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xy[k] = -2.0 * to_xy_xyz[k] * tke_0 + 4.0 * to_xyzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xz[k] = to_xy_x[k] - 2.0 * to_xy_xzz[k] * tke_0 - 2.0 * to_xyzz_x[k] * tbe_0 + 4.0 * to_xyzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_yy[k] = -2.0 * to_xy_yyz[k] * tke_0 + 4.0 * to_xyzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xyz_yz[k] = to_xy_y[k] - 2.0 * to_xy_yzz[k] * tke_0 - 2.0 * to_xyzz_y[k] * tbe_0 + 4.0 * to_xyzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_zz[k] = 2.0 * to_xy_z[k] - 2.0 * to_xy_zzz[k] * tke_0 - 4.0 * to_xyzz_z[k] * tbe_0 + 4.0 * to_xyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 510-516 components of targeted buffer : FD

        auto to_z_z_xzz_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 30);

        auto to_z_z_xzz_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 31);

        auto to_z_z_xzz_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 32);

        auto to_z_z_xzz_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 33);

        auto to_z_z_xzz_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 34);

        auto to_z_z_xzz_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 35);

        #pragma omp simd aligned(to_xz_x, to_xz_xxz, to_xz_xyz, to_xz_xzz, to_xz_y, to_xz_yyz, to_xz_yzz, to_xz_z, to_xz_zzz, to_xzzz_x, to_xzzz_xxz, to_xzzz_xyz, to_xzzz_xzz, to_xzzz_y, to_xzzz_yyz, to_xzzz_yzz, to_xzzz_z, to_xzzz_zzz, to_z_z_xzz_xx, to_z_z_xzz_xy, to_z_z_xzz_xz, to_z_z_xzz_yy, to_z_z_xzz_yz, to_z_z_xzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xzz_xx[k] = -4.0 * to_xz_xxz[k] * tke_0 + 4.0 * to_xzzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xy[k] = -4.0 * to_xz_xyz[k] * tke_0 + 4.0 * to_xzzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xz[k] = 2.0 * to_xz_x[k] - 4.0 * to_xz_xzz[k] * tke_0 - 2.0 * to_xzzz_x[k] * tbe_0 + 4.0 * to_xzzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_yy[k] = -4.0 * to_xz_yyz[k] * tke_0 + 4.0 * to_xzzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xzz_yz[k] = 2.0 * to_xz_y[k] - 4.0 * to_xz_yzz[k] * tke_0 - 2.0 * to_xzzz_y[k] * tbe_0 + 4.0 * to_xzzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_zz[k] = 4.0 * to_xz_z[k] - 4.0 * to_xz_zzz[k] * tke_0 - 4.0 * to_xzzz_z[k] * tbe_0 + 4.0 * to_xzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 516-522 components of targeted buffer : FD

        auto to_z_z_yyy_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 36);

        auto to_z_z_yyy_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 37);

        auto to_z_z_yyy_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 38);

        auto to_z_z_yyy_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 39);

        auto to_z_z_yyy_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 40);

        auto to_z_z_yyy_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 41);

        #pragma omp simd aligned(to_yyyz_x, to_yyyz_xxz, to_yyyz_xyz, to_yyyz_xzz, to_yyyz_y, to_yyyz_yyz, to_yyyz_yzz, to_yyyz_z, to_yyyz_zzz, to_z_z_yyy_xx, to_z_z_yyy_xy, to_z_z_yyy_xz, to_z_z_yyy_yy, to_z_z_yyy_yz, to_z_z_yyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyy_xx[k] = 4.0 * to_yyyz_xxz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xy[k] = 4.0 * to_yyyz_xyz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xz[k] = -2.0 * to_yyyz_x[k] * tbe_0 + 4.0 * to_yyyz_xzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_yy[k] = 4.0 * to_yyyz_yyz[k] * tbe_0 * tke_0;

            to_z_z_yyy_yz[k] = -2.0 * to_yyyz_y[k] * tbe_0 + 4.0 * to_yyyz_yzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_zz[k] = -4.0 * to_yyyz_z[k] * tbe_0 + 4.0 * to_yyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 522-528 components of targeted buffer : FD

        auto to_z_z_yyz_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 42);

        auto to_z_z_yyz_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 43);

        auto to_z_z_yyz_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 44);

        auto to_z_z_yyz_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 45);

        auto to_z_z_yyz_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 46);

        auto to_z_z_yyz_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 47);

        #pragma omp simd aligned(to_yy_x, to_yy_xxz, to_yy_xyz, to_yy_xzz, to_yy_y, to_yy_yyz, to_yy_yzz, to_yy_z, to_yy_zzz, to_yyzz_x, to_yyzz_xxz, to_yyzz_xyz, to_yyzz_xzz, to_yyzz_y, to_yyzz_yyz, to_yyzz_yzz, to_yyzz_z, to_yyzz_zzz, to_z_z_yyz_xx, to_z_z_yyz_xy, to_z_z_yyz_xz, to_z_z_yyz_yy, to_z_z_yyz_yz, to_z_z_yyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyz_xx[k] = -2.0 * to_yy_xxz[k] * tke_0 + 4.0 * to_yyzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xy[k] = -2.0 * to_yy_xyz[k] * tke_0 + 4.0 * to_yyzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xz[k] = to_yy_x[k] - 2.0 * to_yy_xzz[k] * tke_0 - 2.0 * to_yyzz_x[k] * tbe_0 + 4.0 * to_yyzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_yy[k] = -2.0 * to_yy_yyz[k] * tke_0 + 4.0 * to_yyzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_yyz_yz[k] = to_yy_y[k] - 2.0 * to_yy_yzz[k] * tke_0 - 2.0 * to_yyzz_y[k] * tbe_0 + 4.0 * to_yyzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_zz[k] = 2.0 * to_yy_z[k] - 2.0 * to_yy_zzz[k] * tke_0 - 4.0 * to_yyzz_z[k] * tbe_0 + 4.0 * to_yyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 528-534 components of targeted buffer : FD

        auto to_z_z_yzz_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 48);

        auto to_z_z_yzz_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 49);

        auto to_z_z_yzz_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 50);

        auto to_z_z_yzz_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 51);

        auto to_z_z_yzz_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 52);

        auto to_z_z_yzz_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 53);

        #pragma omp simd aligned(to_yz_x, to_yz_xxz, to_yz_xyz, to_yz_xzz, to_yz_y, to_yz_yyz, to_yz_yzz, to_yz_z, to_yz_zzz, to_yzzz_x, to_yzzz_xxz, to_yzzz_xyz, to_yzzz_xzz, to_yzzz_y, to_yzzz_yyz, to_yzzz_yzz, to_yzzz_z, to_yzzz_zzz, to_z_z_yzz_xx, to_z_z_yzz_xy, to_z_z_yzz_xz, to_z_z_yzz_yy, to_z_z_yzz_yz, to_z_z_yzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yzz_xx[k] = -4.0 * to_yz_xxz[k] * tke_0 + 4.0 * to_yzzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xy[k] = -4.0 * to_yz_xyz[k] * tke_0 + 4.0 * to_yzzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xz[k] = 2.0 * to_yz_x[k] - 4.0 * to_yz_xzz[k] * tke_0 - 2.0 * to_yzzz_x[k] * tbe_0 + 4.0 * to_yzzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_yy[k] = -4.0 * to_yz_yyz[k] * tke_0 + 4.0 * to_yzzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_yzz_yz[k] = 2.0 * to_yz_y[k] - 4.0 * to_yz_yzz[k] * tke_0 - 2.0 * to_yzzz_y[k] * tbe_0 + 4.0 * to_yzzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_zz[k] = 4.0 * to_yz_z[k] - 4.0 * to_yz_zzz[k] * tke_0 - 4.0 * to_yzzz_z[k] * tbe_0 + 4.0 * to_yzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 534-540 components of targeted buffer : FD

        auto to_z_z_zzz_xx = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 54);

        auto to_z_z_zzz_xy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 55);

        auto to_z_z_zzz_xz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 56);

        auto to_z_z_zzz_yy = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 57);

        auto to_z_z_zzz_yz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 58);

        auto to_z_z_zzz_zz = pbuffer.data(idx_op_geom_101_fd + 8 * op_comps * 60 + i * 60 + 59);

        #pragma omp simd aligned(to_z_z_zzz_xx, to_z_z_zzz_xy, to_z_z_zzz_xz, to_z_z_zzz_yy, to_z_z_zzz_yz, to_z_z_zzz_zz, to_zz_x, to_zz_xxz, to_zz_xyz, to_zz_xzz, to_zz_y, to_zz_yyz, to_zz_yzz, to_zz_z, to_zz_zzz, to_zzzz_x, to_zzzz_xxz, to_zzzz_xyz, to_zzzz_xzz, to_zzzz_y, to_zzzz_yyz, to_zzzz_yzz, to_zzzz_z, to_zzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zzz_xx[k] = -6.0 * to_zz_xxz[k] * tke_0 + 4.0 * to_zzzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xy[k] = -6.0 * to_zz_xyz[k] * tke_0 + 4.0 * to_zzzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xz[k] = 3.0 * to_zz_x[k] - 6.0 * to_zz_xzz[k] * tke_0 - 2.0 * to_zzzz_x[k] * tbe_0 + 4.0 * to_zzzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_yy[k] = -6.0 * to_zz_yyz[k] * tke_0 + 4.0 * to_zzzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_zzz_yz[k] = 3.0 * to_zz_y[k] - 6.0 * to_zz_yzz[k] * tke_0 - 2.0 * to_zzzz_y[k] * tbe_0 + 4.0 * to_zzzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_zz[k] = 6.0 * to_zz_z[k] - 6.0 * to_zz_zzz[k] * tke_0 - 4.0 * to_zzzz_z[k] * tbe_0 + 4.0 * to_zzzz_zzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

