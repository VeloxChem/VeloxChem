#include "GeometricalDerivatives1X1ForFF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_ff(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_ff,
                        const size_t idx_op_dd,
                        const size_t idx_op_dg,
                        const size_t idx_op_gd,
                        const size_t idx_op_gg,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : DD

        auto to_xx_xx = pbuffer.data(idx_op_dd + i * 36 + 0);

        auto to_xx_xy = pbuffer.data(idx_op_dd + i * 36 + 1);

        auto to_xx_xz = pbuffer.data(idx_op_dd + i * 36 + 2);

        auto to_xx_yy = pbuffer.data(idx_op_dd + i * 36 + 3);

        auto to_xx_yz = pbuffer.data(idx_op_dd + i * 36 + 4);

        auto to_xx_zz = pbuffer.data(idx_op_dd + i * 36 + 5);

        auto to_xy_xx = pbuffer.data(idx_op_dd + i * 36 + 6);

        auto to_xy_xy = pbuffer.data(idx_op_dd + i * 36 + 7);

        auto to_xy_xz = pbuffer.data(idx_op_dd + i * 36 + 8);

        auto to_xy_yy = pbuffer.data(idx_op_dd + i * 36 + 9);

        auto to_xy_yz = pbuffer.data(idx_op_dd + i * 36 + 10);

        auto to_xy_zz = pbuffer.data(idx_op_dd + i * 36 + 11);

        auto to_xz_xx = pbuffer.data(idx_op_dd + i * 36 + 12);

        auto to_xz_xy = pbuffer.data(idx_op_dd + i * 36 + 13);

        auto to_xz_xz = pbuffer.data(idx_op_dd + i * 36 + 14);

        auto to_xz_yy = pbuffer.data(idx_op_dd + i * 36 + 15);

        auto to_xz_yz = pbuffer.data(idx_op_dd + i * 36 + 16);

        auto to_xz_zz = pbuffer.data(idx_op_dd + i * 36 + 17);

        auto to_yy_xx = pbuffer.data(idx_op_dd + i * 36 + 18);

        auto to_yy_xy = pbuffer.data(idx_op_dd + i * 36 + 19);

        auto to_yy_xz = pbuffer.data(idx_op_dd + i * 36 + 20);

        auto to_yy_yy = pbuffer.data(idx_op_dd + i * 36 + 21);

        auto to_yy_yz = pbuffer.data(idx_op_dd + i * 36 + 22);

        auto to_yy_zz = pbuffer.data(idx_op_dd + i * 36 + 23);

        auto to_yz_xx = pbuffer.data(idx_op_dd + i * 36 + 24);

        auto to_yz_xy = pbuffer.data(idx_op_dd + i * 36 + 25);

        auto to_yz_xz = pbuffer.data(idx_op_dd + i * 36 + 26);

        auto to_yz_yy = pbuffer.data(idx_op_dd + i * 36 + 27);

        auto to_yz_yz = pbuffer.data(idx_op_dd + i * 36 + 28);

        auto to_yz_zz = pbuffer.data(idx_op_dd + i * 36 + 29);

        auto to_zz_xx = pbuffer.data(idx_op_dd + i * 36 + 30);

        auto to_zz_xy = pbuffer.data(idx_op_dd + i * 36 + 31);

        auto to_zz_xz = pbuffer.data(idx_op_dd + i * 36 + 32);

        auto to_zz_yy = pbuffer.data(idx_op_dd + i * 36 + 33);

        auto to_zz_yz = pbuffer.data(idx_op_dd + i * 36 + 34);

        auto to_zz_zz = pbuffer.data(idx_op_dd + i * 36 + 35);

        // Set up components of auxiliary buffer : DG

        auto to_xx_xxxx = pbuffer.data(idx_op_dg + i * 90 + 0);

        auto to_xx_xxxy = pbuffer.data(idx_op_dg + i * 90 + 1);

        auto to_xx_xxxz = pbuffer.data(idx_op_dg + i * 90 + 2);

        auto to_xx_xxyy = pbuffer.data(idx_op_dg + i * 90 + 3);

        auto to_xx_xxyz = pbuffer.data(idx_op_dg + i * 90 + 4);

        auto to_xx_xxzz = pbuffer.data(idx_op_dg + i * 90 + 5);

        auto to_xx_xyyy = pbuffer.data(idx_op_dg + i * 90 + 6);

        auto to_xx_xyyz = pbuffer.data(idx_op_dg + i * 90 + 7);

        auto to_xx_xyzz = pbuffer.data(idx_op_dg + i * 90 + 8);

        auto to_xx_xzzz = pbuffer.data(idx_op_dg + i * 90 + 9);

        auto to_xx_yyyy = pbuffer.data(idx_op_dg + i * 90 + 10);

        auto to_xx_yyyz = pbuffer.data(idx_op_dg + i * 90 + 11);

        auto to_xx_yyzz = pbuffer.data(idx_op_dg + i * 90 + 12);

        auto to_xx_yzzz = pbuffer.data(idx_op_dg + i * 90 + 13);

        auto to_xx_zzzz = pbuffer.data(idx_op_dg + i * 90 + 14);

        auto to_xy_xxxx = pbuffer.data(idx_op_dg + i * 90 + 15);

        auto to_xy_xxxy = pbuffer.data(idx_op_dg + i * 90 + 16);

        auto to_xy_xxxz = pbuffer.data(idx_op_dg + i * 90 + 17);

        auto to_xy_xxyy = pbuffer.data(idx_op_dg + i * 90 + 18);

        auto to_xy_xxyz = pbuffer.data(idx_op_dg + i * 90 + 19);

        auto to_xy_xxzz = pbuffer.data(idx_op_dg + i * 90 + 20);

        auto to_xy_xyyy = pbuffer.data(idx_op_dg + i * 90 + 21);

        auto to_xy_xyyz = pbuffer.data(idx_op_dg + i * 90 + 22);

        auto to_xy_xyzz = pbuffer.data(idx_op_dg + i * 90 + 23);

        auto to_xy_xzzz = pbuffer.data(idx_op_dg + i * 90 + 24);

        auto to_xy_yyyy = pbuffer.data(idx_op_dg + i * 90 + 25);

        auto to_xy_yyyz = pbuffer.data(idx_op_dg + i * 90 + 26);

        auto to_xy_yyzz = pbuffer.data(idx_op_dg + i * 90 + 27);

        auto to_xy_yzzz = pbuffer.data(idx_op_dg + i * 90 + 28);

        auto to_xy_zzzz = pbuffer.data(idx_op_dg + i * 90 + 29);

        auto to_xz_xxxx = pbuffer.data(idx_op_dg + i * 90 + 30);

        auto to_xz_xxxy = pbuffer.data(idx_op_dg + i * 90 + 31);

        auto to_xz_xxxz = pbuffer.data(idx_op_dg + i * 90 + 32);

        auto to_xz_xxyy = pbuffer.data(idx_op_dg + i * 90 + 33);

        auto to_xz_xxyz = pbuffer.data(idx_op_dg + i * 90 + 34);

        auto to_xz_xxzz = pbuffer.data(idx_op_dg + i * 90 + 35);

        auto to_xz_xyyy = pbuffer.data(idx_op_dg + i * 90 + 36);

        auto to_xz_xyyz = pbuffer.data(idx_op_dg + i * 90 + 37);

        auto to_xz_xyzz = pbuffer.data(idx_op_dg + i * 90 + 38);

        auto to_xz_xzzz = pbuffer.data(idx_op_dg + i * 90 + 39);

        auto to_xz_yyyy = pbuffer.data(idx_op_dg + i * 90 + 40);

        auto to_xz_yyyz = pbuffer.data(idx_op_dg + i * 90 + 41);

        auto to_xz_yyzz = pbuffer.data(idx_op_dg + i * 90 + 42);

        auto to_xz_yzzz = pbuffer.data(idx_op_dg + i * 90 + 43);

        auto to_xz_zzzz = pbuffer.data(idx_op_dg + i * 90 + 44);

        auto to_yy_xxxx = pbuffer.data(idx_op_dg + i * 90 + 45);

        auto to_yy_xxxy = pbuffer.data(idx_op_dg + i * 90 + 46);

        auto to_yy_xxxz = pbuffer.data(idx_op_dg + i * 90 + 47);

        auto to_yy_xxyy = pbuffer.data(idx_op_dg + i * 90 + 48);

        auto to_yy_xxyz = pbuffer.data(idx_op_dg + i * 90 + 49);

        auto to_yy_xxzz = pbuffer.data(idx_op_dg + i * 90 + 50);

        auto to_yy_xyyy = pbuffer.data(idx_op_dg + i * 90 + 51);

        auto to_yy_xyyz = pbuffer.data(idx_op_dg + i * 90 + 52);

        auto to_yy_xyzz = pbuffer.data(idx_op_dg + i * 90 + 53);

        auto to_yy_xzzz = pbuffer.data(idx_op_dg + i * 90 + 54);

        auto to_yy_yyyy = pbuffer.data(idx_op_dg + i * 90 + 55);

        auto to_yy_yyyz = pbuffer.data(idx_op_dg + i * 90 + 56);

        auto to_yy_yyzz = pbuffer.data(idx_op_dg + i * 90 + 57);

        auto to_yy_yzzz = pbuffer.data(idx_op_dg + i * 90 + 58);

        auto to_yy_zzzz = pbuffer.data(idx_op_dg + i * 90 + 59);

        auto to_yz_xxxx = pbuffer.data(idx_op_dg + i * 90 + 60);

        auto to_yz_xxxy = pbuffer.data(idx_op_dg + i * 90 + 61);

        auto to_yz_xxxz = pbuffer.data(idx_op_dg + i * 90 + 62);

        auto to_yz_xxyy = pbuffer.data(idx_op_dg + i * 90 + 63);

        auto to_yz_xxyz = pbuffer.data(idx_op_dg + i * 90 + 64);

        auto to_yz_xxzz = pbuffer.data(idx_op_dg + i * 90 + 65);

        auto to_yz_xyyy = pbuffer.data(idx_op_dg + i * 90 + 66);

        auto to_yz_xyyz = pbuffer.data(idx_op_dg + i * 90 + 67);

        auto to_yz_xyzz = pbuffer.data(idx_op_dg + i * 90 + 68);

        auto to_yz_xzzz = pbuffer.data(idx_op_dg + i * 90 + 69);

        auto to_yz_yyyy = pbuffer.data(idx_op_dg + i * 90 + 70);

        auto to_yz_yyyz = pbuffer.data(idx_op_dg + i * 90 + 71);

        auto to_yz_yyzz = pbuffer.data(idx_op_dg + i * 90 + 72);

        auto to_yz_yzzz = pbuffer.data(idx_op_dg + i * 90 + 73);

        auto to_yz_zzzz = pbuffer.data(idx_op_dg + i * 90 + 74);

        auto to_zz_xxxx = pbuffer.data(idx_op_dg + i * 90 + 75);

        auto to_zz_xxxy = pbuffer.data(idx_op_dg + i * 90 + 76);

        auto to_zz_xxxz = pbuffer.data(idx_op_dg + i * 90 + 77);

        auto to_zz_xxyy = pbuffer.data(idx_op_dg + i * 90 + 78);

        auto to_zz_xxyz = pbuffer.data(idx_op_dg + i * 90 + 79);

        auto to_zz_xxzz = pbuffer.data(idx_op_dg + i * 90 + 80);

        auto to_zz_xyyy = pbuffer.data(idx_op_dg + i * 90 + 81);

        auto to_zz_xyyz = pbuffer.data(idx_op_dg + i * 90 + 82);

        auto to_zz_xyzz = pbuffer.data(idx_op_dg + i * 90 + 83);

        auto to_zz_xzzz = pbuffer.data(idx_op_dg + i * 90 + 84);

        auto to_zz_yyyy = pbuffer.data(idx_op_dg + i * 90 + 85);

        auto to_zz_yyyz = pbuffer.data(idx_op_dg + i * 90 + 86);

        auto to_zz_yyzz = pbuffer.data(idx_op_dg + i * 90 + 87);

        auto to_zz_yzzz = pbuffer.data(idx_op_dg + i * 90 + 88);

        auto to_zz_zzzz = pbuffer.data(idx_op_dg + i * 90 + 89);

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

        // Set up 0-10 components of targeted buffer : FF

        auto to_x_x_xxx_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 0);

        auto to_x_x_xxx_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 1);

        auto to_x_x_xxx_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 2);

        auto to_x_x_xxx_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 3);

        auto to_x_x_xxx_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 4);

        auto to_x_x_xxx_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 5);

        auto to_x_x_xxx_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 6);

        auto to_x_x_xxx_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 7);

        auto to_x_x_xxx_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 8);

        auto to_x_x_xxx_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_x_x_xxx_xxx, to_x_x_xxx_xxy, to_x_x_xxx_xxz, to_x_x_xxx_xyy, to_x_x_xxx_xyz, to_x_x_xxx_xzz, to_x_x_xxx_yyy, to_x_x_xxx_yyz, to_x_x_xxx_yzz, to_x_x_xxx_zzz, to_xx_xx, to_xx_xxxx, to_xx_xxxy, to_xx_xxxz, to_xx_xxyy, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yz, to_xx_zz, to_xxxx_xx, to_xxxx_xxxx, to_xxxx_xxxy, to_xxxx_xxxz, to_xxxx_xxyy, to_xxxx_xxyz, to_xxxx_xxzz, to_xxxx_xy, to_xxxx_xyyy, to_xxxx_xyyz, to_xxxx_xyzz, to_xxxx_xz, to_xxxx_xzzz, to_xxxx_yy, to_xxxx_yz, to_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxx_xxx[k] = 9.0 * to_xx_xx[k] - 6.0 * to_xx_xxxx[k] * tke_0 - 6.0 * to_xxxx_xx[k] * tbe_0 + 4.0 * to_xxxx_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xxx_xxy[k] = 6.0 * to_xx_xy[k] - 6.0 * to_xx_xxxy[k] * tke_0 - 4.0 * to_xxxx_xy[k] * tbe_0 + 4.0 * to_xxxx_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xxx_xxz[k] = 6.0 * to_xx_xz[k] - 6.0 * to_xx_xxxz[k] * tke_0 - 4.0 * to_xxxx_xz[k] * tbe_0 + 4.0 * to_xxxx_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xxx_xyy[k] = 3.0 * to_xx_yy[k] - 6.0 * to_xx_xxyy[k] * tke_0 - 2.0 * to_xxxx_yy[k] * tbe_0 + 4.0 * to_xxxx_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xxx_xyz[k] = 3.0 * to_xx_yz[k] - 6.0 * to_xx_xxyz[k] * tke_0 - 2.0 * to_xxxx_yz[k] * tbe_0 + 4.0 * to_xxxx_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xxx_xzz[k] = 3.0 * to_xx_zz[k] - 6.0 * to_xx_xxzz[k] * tke_0 - 2.0 * to_xxxx_zz[k] * tbe_0 + 4.0 * to_xxxx_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xxx_yyy[k] = -6.0 * to_xx_xyyy[k] * tke_0 + 4.0 * to_xxxx_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xxx_yyz[k] = -6.0 * to_xx_xyyz[k] * tke_0 + 4.0 * to_xxxx_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xxx_yzz[k] = -6.0 * to_xx_xyzz[k] * tke_0 + 4.0 * to_xxxx_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xxx_zzz[k] = -6.0 * to_xx_xzzz[k] * tke_0 + 4.0 * to_xxxx_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 10-20 components of targeted buffer : FF

        auto to_x_x_xxy_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 10);

        auto to_x_x_xxy_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 11);

        auto to_x_x_xxy_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 12);

        auto to_x_x_xxy_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 13);

        auto to_x_x_xxy_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 14);

        auto to_x_x_xxy_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 15);

        auto to_x_x_xxy_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 16);

        auto to_x_x_xxy_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 17);

        auto to_x_x_xxy_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 18);

        auto to_x_x_xxy_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_x_x_xxy_xxx, to_x_x_xxy_xxy, to_x_x_xxy_xxz, to_x_x_xxy_xyy, to_x_x_xxy_xyz, to_x_x_xxy_xzz, to_x_x_xxy_yyy, to_x_x_xxy_yyz, to_x_x_xxy_yzz, to_x_x_xxy_zzz, to_xxxy_xx, to_xxxy_xxxx, to_xxxy_xxxy, to_xxxy_xxxz, to_xxxy_xxyy, to_xxxy_xxyz, to_xxxy_xxzz, to_xxxy_xy, to_xxxy_xyyy, to_xxxy_xyyz, to_xxxy_xyzz, to_xxxy_xz, to_xxxy_xzzz, to_xxxy_yy, to_xxxy_yz, to_xxxy_zz, to_xy_xx, to_xy_xxxx, to_xy_xxxy, to_xy_xxxz, to_xy_xxyy, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yz, to_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxy_xxx[k] = 6.0 * to_xy_xx[k] - 4.0 * to_xy_xxxx[k] * tke_0 - 6.0 * to_xxxy_xx[k] * tbe_0 + 4.0 * to_xxxy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xxy_xxy[k] = 4.0 * to_xy_xy[k] - 4.0 * to_xy_xxxy[k] * tke_0 - 4.0 * to_xxxy_xy[k] * tbe_0 + 4.0 * to_xxxy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xxy_xxz[k] = 4.0 * to_xy_xz[k] - 4.0 * to_xy_xxxz[k] * tke_0 - 4.0 * to_xxxy_xz[k] * tbe_0 + 4.0 * to_xxxy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xxy_xyy[k] = 2.0 * to_xy_yy[k] - 4.0 * to_xy_xxyy[k] * tke_0 - 2.0 * to_xxxy_yy[k] * tbe_0 + 4.0 * to_xxxy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xxy_xyz[k] = 2.0 * to_xy_yz[k] - 4.0 * to_xy_xxyz[k] * tke_0 - 2.0 * to_xxxy_yz[k] * tbe_0 + 4.0 * to_xxxy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xxy_xzz[k] = 2.0 * to_xy_zz[k] - 4.0 * to_xy_xxzz[k] * tke_0 - 2.0 * to_xxxy_zz[k] * tbe_0 + 4.0 * to_xxxy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xxy_yyy[k] = -4.0 * to_xy_xyyy[k] * tke_0 + 4.0 * to_xxxy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xxy_yyz[k] = -4.0 * to_xy_xyyz[k] * tke_0 + 4.0 * to_xxxy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xxy_yzz[k] = -4.0 * to_xy_xyzz[k] * tke_0 + 4.0 * to_xxxy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xxy_zzz[k] = -4.0 * to_xy_xzzz[k] * tke_0 + 4.0 * to_xxxy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 20-30 components of targeted buffer : FF

        auto to_x_x_xxz_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 20);

        auto to_x_x_xxz_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 21);

        auto to_x_x_xxz_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 22);

        auto to_x_x_xxz_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 23);

        auto to_x_x_xxz_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 24);

        auto to_x_x_xxz_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 25);

        auto to_x_x_xxz_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 26);

        auto to_x_x_xxz_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 27);

        auto to_x_x_xxz_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 28);

        auto to_x_x_xxz_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_x_x_xxz_xxx, to_x_x_xxz_xxy, to_x_x_xxz_xxz, to_x_x_xxz_xyy, to_x_x_xxz_xyz, to_x_x_xxz_xzz, to_x_x_xxz_yyy, to_x_x_xxz_yyz, to_x_x_xxz_yzz, to_x_x_xxz_zzz, to_xxxz_xx, to_xxxz_xxxx, to_xxxz_xxxy, to_xxxz_xxxz, to_xxxz_xxyy, to_xxxz_xxyz, to_xxxz_xxzz, to_xxxz_xy, to_xxxz_xyyy, to_xxxz_xyyz, to_xxxz_xyzz, to_xxxz_xz, to_xxxz_xzzz, to_xxxz_yy, to_xxxz_yz, to_xxxz_zz, to_xz_xx, to_xz_xxxx, to_xz_xxxy, to_xz_xxxz, to_xz_xxyy, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yz, to_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxz_xxx[k] = 6.0 * to_xz_xx[k] - 4.0 * to_xz_xxxx[k] * tke_0 - 6.0 * to_xxxz_xx[k] * tbe_0 + 4.0 * to_xxxz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xxz_xxy[k] = 4.0 * to_xz_xy[k] - 4.0 * to_xz_xxxy[k] * tke_0 - 4.0 * to_xxxz_xy[k] * tbe_0 + 4.0 * to_xxxz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xxz_xxz[k] = 4.0 * to_xz_xz[k] - 4.0 * to_xz_xxxz[k] * tke_0 - 4.0 * to_xxxz_xz[k] * tbe_0 + 4.0 * to_xxxz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xxz_xyy[k] = 2.0 * to_xz_yy[k] - 4.0 * to_xz_xxyy[k] * tke_0 - 2.0 * to_xxxz_yy[k] * tbe_0 + 4.0 * to_xxxz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xxz_xyz[k] = 2.0 * to_xz_yz[k] - 4.0 * to_xz_xxyz[k] * tke_0 - 2.0 * to_xxxz_yz[k] * tbe_0 + 4.0 * to_xxxz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xxz_xzz[k] = 2.0 * to_xz_zz[k] - 4.0 * to_xz_xxzz[k] * tke_0 - 2.0 * to_xxxz_zz[k] * tbe_0 + 4.0 * to_xxxz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xxz_yyy[k] = -4.0 * to_xz_xyyy[k] * tke_0 + 4.0 * to_xxxz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xxz_yyz[k] = -4.0 * to_xz_xyyz[k] * tke_0 + 4.0 * to_xxxz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xxz_yzz[k] = -4.0 * to_xz_xyzz[k] * tke_0 + 4.0 * to_xxxz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xxz_zzz[k] = -4.0 * to_xz_xzzz[k] * tke_0 + 4.0 * to_xxxz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-40 components of targeted buffer : FF

        auto to_x_x_xyy_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 30);

        auto to_x_x_xyy_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 31);

        auto to_x_x_xyy_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 32);

        auto to_x_x_xyy_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 33);

        auto to_x_x_xyy_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 34);

        auto to_x_x_xyy_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 35);

        auto to_x_x_xyy_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 36);

        auto to_x_x_xyy_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 37);

        auto to_x_x_xyy_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 38);

        auto to_x_x_xyy_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_x_x_xyy_xxx, to_x_x_xyy_xxy, to_x_x_xyy_xxz, to_x_x_xyy_xyy, to_x_x_xyy_xyz, to_x_x_xyy_xzz, to_x_x_xyy_yyy, to_x_x_xyy_yyz, to_x_x_xyy_yzz, to_x_x_xyy_zzz, to_xxyy_xx, to_xxyy_xxxx, to_xxyy_xxxy, to_xxyy_xxxz, to_xxyy_xxyy, to_xxyy_xxyz, to_xxyy_xxzz, to_xxyy_xy, to_xxyy_xyyy, to_xxyy_xyyz, to_xxyy_xyzz, to_xxyy_xz, to_xxyy_xzzz, to_xxyy_yy, to_xxyy_yz, to_xxyy_zz, to_yy_xx, to_yy_xxxx, to_yy_xxxy, to_yy_xxxz, to_yy_xxyy, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yz, to_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyy_xxx[k] = 3.0 * to_yy_xx[k] - 2.0 * to_yy_xxxx[k] * tke_0 - 6.0 * to_xxyy_xx[k] * tbe_0 + 4.0 * to_xxyy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xyy_xxy[k] = 2.0 * to_yy_xy[k] - 2.0 * to_yy_xxxy[k] * tke_0 - 4.0 * to_xxyy_xy[k] * tbe_0 + 4.0 * to_xxyy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xyy_xxz[k] = 2.0 * to_yy_xz[k] - 2.0 * to_yy_xxxz[k] * tke_0 - 4.0 * to_xxyy_xz[k] * tbe_0 + 4.0 * to_xxyy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xyy_xyy[k] = to_yy_yy[k] - 2.0 * to_yy_xxyy[k] * tke_0 - 2.0 * to_xxyy_yy[k] * tbe_0 + 4.0 * to_xxyy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xyy_xyz[k] = to_yy_yz[k] - 2.0 * to_yy_xxyz[k] * tke_0 - 2.0 * to_xxyy_yz[k] * tbe_0 + 4.0 * to_xxyy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xyy_xzz[k] = to_yy_zz[k] - 2.0 * to_yy_xxzz[k] * tke_0 - 2.0 * to_xxyy_zz[k] * tbe_0 + 4.0 * to_xxyy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xyy_yyy[k] = -2.0 * to_yy_xyyy[k] * tke_0 + 4.0 * to_xxyy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xyy_yyz[k] = -2.0 * to_yy_xyyz[k] * tke_0 + 4.0 * to_xxyy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xyy_yzz[k] = -2.0 * to_yy_xyzz[k] * tke_0 + 4.0 * to_xxyy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xyy_zzz[k] = -2.0 * to_yy_xzzz[k] * tke_0 + 4.0 * to_xxyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 40-50 components of targeted buffer : FF

        auto to_x_x_xyz_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 40);

        auto to_x_x_xyz_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 41);

        auto to_x_x_xyz_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 42);

        auto to_x_x_xyz_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 43);

        auto to_x_x_xyz_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 44);

        auto to_x_x_xyz_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 45);

        auto to_x_x_xyz_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 46);

        auto to_x_x_xyz_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 47);

        auto to_x_x_xyz_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 48);

        auto to_x_x_xyz_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_x_x_xyz_xxx, to_x_x_xyz_xxy, to_x_x_xyz_xxz, to_x_x_xyz_xyy, to_x_x_xyz_xyz, to_x_x_xyz_xzz, to_x_x_xyz_yyy, to_x_x_xyz_yyz, to_x_x_xyz_yzz, to_x_x_xyz_zzz, to_xxyz_xx, to_xxyz_xxxx, to_xxyz_xxxy, to_xxyz_xxxz, to_xxyz_xxyy, to_xxyz_xxyz, to_xxyz_xxzz, to_xxyz_xy, to_xxyz_xyyy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_xzzz, to_xxyz_yy, to_xxyz_yz, to_xxyz_zz, to_yz_xx, to_yz_xxxx, to_yz_xxxy, to_yz_xxxz, to_yz_xxyy, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyz_xxx[k] = 3.0 * to_yz_xx[k] - 2.0 * to_yz_xxxx[k] * tke_0 - 6.0 * to_xxyz_xx[k] * tbe_0 + 4.0 * to_xxyz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xyz_xxy[k] = 2.0 * to_yz_xy[k] - 2.0 * to_yz_xxxy[k] * tke_0 - 4.0 * to_xxyz_xy[k] * tbe_0 + 4.0 * to_xxyz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xyz_xxz[k] = 2.0 * to_yz_xz[k] - 2.0 * to_yz_xxxz[k] * tke_0 - 4.0 * to_xxyz_xz[k] * tbe_0 + 4.0 * to_xxyz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xyz_xyy[k] = to_yz_yy[k] - 2.0 * to_yz_xxyy[k] * tke_0 - 2.0 * to_xxyz_yy[k] * tbe_0 + 4.0 * to_xxyz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xyz_xyz[k] = to_yz_yz[k] - 2.0 * to_yz_xxyz[k] * tke_0 - 2.0 * to_xxyz_yz[k] * tbe_0 + 4.0 * to_xxyz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xyz_xzz[k] = to_yz_zz[k] - 2.0 * to_yz_xxzz[k] * tke_0 - 2.0 * to_xxyz_zz[k] * tbe_0 + 4.0 * to_xxyz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xyz_yyy[k] = -2.0 * to_yz_xyyy[k] * tke_0 + 4.0 * to_xxyz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xyz_yyz[k] = -2.0 * to_yz_xyyz[k] * tke_0 + 4.0 * to_xxyz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xyz_yzz[k] = -2.0 * to_yz_xyzz[k] * tke_0 + 4.0 * to_xxyz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xyz_zzz[k] = -2.0 * to_yz_xzzz[k] * tke_0 + 4.0 * to_xxyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 50-60 components of targeted buffer : FF

        auto to_x_x_xzz_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 50);

        auto to_x_x_xzz_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 51);

        auto to_x_x_xzz_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 52);

        auto to_x_x_xzz_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 53);

        auto to_x_x_xzz_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 54);

        auto to_x_x_xzz_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 55);

        auto to_x_x_xzz_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 56);

        auto to_x_x_xzz_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 57);

        auto to_x_x_xzz_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 58);

        auto to_x_x_xzz_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_x_x_xzz_xxx, to_x_x_xzz_xxy, to_x_x_xzz_xxz, to_x_x_xzz_xyy, to_x_x_xzz_xyz, to_x_x_xzz_xzz, to_x_x_xzz_yyy, to_x_x_xzz_yyz, to_x_x_xzz_yzz, to_x_x_xzz_zzz, to_xxzz_xx, to_xxzz_xxxx, to_xxzz_xxxy, to_xxzz_xxxz, to_xxzz_xxyy, to_xxzz_xxyz, to_xxzz_xxzz, to_xxzz_xy, to_xxzz_xyyy, to_xxzz_xyyz, to_xxzz_xyzz, to_xxzz_xz, to_xxzz_xzzz, to_xxzz_yy, to_xxzz_yz, to_xxzz_zz, to_zz_xx, to_zz_xxxx, to_zz_xxxy, to_zz_xxxz, to_zz_xxyy, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xzz_xxx[k] = 3.0 * to_zz_xx[k] - 2.0 * to_zz_xxxx[k] * tke_0 - 6.0 * to_xxzz_xx[k] * tbe_0 + 4.0 * to_xxzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_xzz_xxy[k] = 2.0 * to_zz_xy[k] - 2.0 * to_zz_xxxy[k] * tke_0 - 4.0 * to_xxzz_xy[k] * tbe_0 + 4.0 * to_xxzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_xzz_xxz[k] = 2.0 * to_zz_xz[k] - 2.0 * to_zz_xxxz[k] * tke_0 - 4.0 * to_xxzz_xz[k] * tbe_0 + 4.0 * to_xxzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_xzz_xyy[k] = to_zz_yy[k] - 2.0 * to_zz_xxyy[k] * tke_0 - 2.0 * to_xxzz_yy[k] * tbe_0 + 4.0 * to_xxzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_xzz_xyz[k] = to_zz_yz[k] - 2.0 * to_zz_xxyz[k] * tke_0 - 2.0 * to_xxzz_yz[k] * tbe_0 + 4.0 * to_xxzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_xzz_xzz[k] = to_zz_zz[k] - 2.0 * to_zz_xxzz[k] * tke_0 - 2.0 * to_xxzz_zz[k] * tbe_0 + 4.0 * to_xxzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_xzz_yyy[k] = -2.0 * to_zz_xyyy[k] * tke_0 + 4.0 * to_xxzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_xzz_yyz[k] = -2.0 * to_zz_xyyz[k] * tke_0 + 4.0 * to_xxzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_xzz_yzz[k] = -2.0 * to_zz_xyzz[k] * tke_0 + 4.0 * to_xxzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_xzz_zzz[k] = -2.0 * to_zz_xzzz[k] * tke_0 + 4.0 * to_xxzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-70 components of targeted buffer : FF

        auto to_x_x_yyy_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 60);

        auto to_x_x_yyy_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 61);

        auto to_x_x_yyy_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 62);

        auto to_x_x_yyy_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 63);

        auto to_x_x_yyy_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 64);

        auto to_x_x_yyy_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 65);

        auto to_x_x_yyy_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 66);

        auto to_x_x_yyy_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 67);

        auto to_x_x_yyy_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 68);

        auto to_x_x_yyy_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_x_x_yyy_xxx, to_x_x_yyy_xxy, to_x_x_yyy_xxz, to_x_x_yyy_xyy, to_x_x_yyy_xyz, to_x_x_yyy_xzz, to_x_x_yyy_yyy, to_x_x_yyy_yyz, to_x_x_yyy_yzz, to_x_x_yyy_zzz, to_xyyy_xx, to_xyyy_xxxx, to_xyyy_xxxy, to_xyyy_xxxz, to_xyyy_xxyy, to_xyyy_xxyz, to_xyyy_xxzz, to_xyyy_xy, to_xyyy_xyyy, to_xyyy_xyyz, to_xyyy_xyzz, to_xyyy_xz, to_xyyy_xzzz, to_xyyy_yy, to_xyyy_yz, to_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyy_xxx[k] = -6.0 * to_xyyy_xx[k] * tbe_0 + 4.0 * to_xyyy_xxxx[k] * tbe_0 * tke_0;

            to_x_x_yyy_xxy[k] = -4.0 * to_xyyy_xy[k] * tbe_0 + 4.0 * to_xyyy_xxxy[k] * tbe_0 * tke_0;

            to_x_x_yyy_xxz[k] = -4.0 * to_xyyy_xz[k] * tbe_0 + 4.0 * to_xyyy_xxxz[k] * tbe_0 * tke_0;

            to_x_x_yyy_xyy[k] = -2.0 * to_xyyy_yy[k] * tbe_0 + 4.0 * to_xyyy_xxyy[k] * tbe_0 * tke_0;

            to_x_x_yyy_xyz[k] = -2.0 * to_xyyy_yz[k] * tbe_0 + 4.0 * to_xyyy_xxyz[k] * tbe_0 * tke_0;

            to_x_x_yyy_xzz[k] = -2.0 * to_xyyy_zz[k] * tbe_0 + 4.0 * to_xyyy_xxzz[k] * tbe_0 * tke_0;

            to_x_x_yyy_yyy[k] = 4.0 * to_xyyy_xyyy[k] * tbe_0 * tke_0;

            to_x_x_yyy_yyz[k] = 4.0 * to_xyyy_xyyz[k] * tbe_0 * tke_0;

            to_x_x_yyy_yzz[k] = 4.0 * to_xyyy_xyzz[k] * tbe_0 * tke_0;

            to_x_x_yyy_zzz[k] = 4.0 * to_xyyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 70-80 components of targeted buffer : FF

        auto to_x_x_yyz_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 70);

        auto to_x_x_yyz_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 71);

        auto to_x_x_yyz_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 72);

        auto to_x_x_yyz_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 73);

        auto to_x_x_yyz_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 74);

        auto to_x_x_yyz_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 75);

        auto to_x_x_yyz_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 76);

        auto to_x_x_yyz_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 77);

        auto to_x_x_yyz_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 78);

        auto to_x_x_yyz_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_x_x_yyz_xxx, to_x_x_yyz_xxy, to_x_x_yyz_xxz, to_x_x_yyz_xyy, to_x_x_yyz_xyz, to_x_x_yyz_xzz, to_x_x_yyz_yyy, to_x_x_yyz_yyz, to_x_x_yyz_yzz, to_x_x_yyz_zzz, to_xyyz_xx, to_xyyz_xxxx, to_xyyz_xxxy, to_xyyz_xxxz, to_xyyz_xxyy, to_xyyz_xxyz, to_xyyz_xxzz, to_xyyz_xy, to_xyyz_xyyy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_xzzz, to_xyyz_yy, to_xyyz_yz, to_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyz_xxx[k] = -6.0 * to_xyyz_xx[k] * tbe_0 + 4.0 * to_xyyz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_yyz_xxy[k] = -4.0 * to_xyyz_xy[k] * tbe_0 + 4.0 * to_xyyz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_yyz_xxz[k] = -4.0 * to_xyyz_xz[k] * tbe_0 + 4.0 * to_xyyz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_yyz_xyy[k] = -2.0 * to_xyyz_yy[k] * tbe_0 + 4.0 * to_xyyz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_yyz_xyz[k] = -2.0 * to_xyyz_yz[k] * tbe_0 + 4.0 * to_xyyz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_yyz_xzz[k] = -2.0 * to_xyyz_zz[k] * tbe_0 + 4.0 * to_xyyz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_yyz_yyy[k] = 4.0 * to_xyyz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_yyz_yyz[k] = 4.0 * to_xyyz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_yyz_yzz[k] = 4.0 * to_xyyz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_yyz_zzz[k] = 4.0 * to_xyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 80-90 components of targeted buffer : FF

        auto to_x_x_yzz_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 80);

        auto to_x_x_yzz_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 81);

        auto to_x_x_yzz_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 82);

        auto to_x_x_yzz_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 83);

        auto to_x_x_yzz_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 84);

        auto to_x_x_yzz_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 85);

        auto to_x_x_yzz_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 86);

        auto to_x_x_yzz_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 87);

        auto to_x_x_yzz_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 88);

        auto to_x_x_yzz_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_x_x_yzz_xxx, to_x_x_yzz_xxy, to_x_x_yzz_xxz, to_x_x_yzz_xyy, to_x_x_yzz_xyz, to_x_x_yzz_xzz, to_x_x_yzz_yyy, to_x_x_yzz_yyz, to_x_x_yzz_yzz, to_x_x_yzz_zzz, to_xyzz_xx, to_xyzz_xxxx, to_xyzz_xxxy, to_xyzz_xxxz, to_xyzz_xxyy, to_xyzz_xxyz, to_xyzz_xxzz, to_xyzz_xy, to_xyzz_xyyy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_xzzz, to_xyzz_yy, to_xyzz_yz, to_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yzz_xxx[k] = -6.0 * to_xyzz_xx[k] * tbe_0 + 4.0 * to_xyzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_yzz_xxy[k] = -4.0 * to_xyzz_xy[k] * tbe_0 + 4.0 * to_xyzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_yzz_xxz[k] = -4.0 * to_xyzz_xz[k] * tbe_0 + 4.0 * to_xyzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_yzz_xyy[k] = -2.0 * to_xyzz_yy[k] * tbe_0 + 4.0 * to_xyzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_yzz_xyz[k] = -2.0 * to_xyzz_yz[k] * tbe_0 + 4.0 * to_xyzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_yzz_xzz[k] = -2.0 * to_xyzz_zz[k] * tbe_0 + 4.0 * to_xyzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_yzz_yyy[k] = 4.0 * to_xyzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_yzz_yyz[k] = 4.0 * to_xyzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_yzz_yzz[k] = 4.0 * to_xyzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_yzz_zzz[k] = 4.0 * to_xyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-100 components of targeted buffer : FF

        auto to_x_x_zzz_xxx = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 90);

        auto to_x_x_zzz_xxy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 91);

        auto to_x_x_zzz_xxz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 92);

        auto to_x_x_zzz_xyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 93);

        auto to_x_x_zzz_xyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 94);

        auto to_x_x_zzz_xzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 95);

        auto to_x_x_zzz_yyy = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 96);

        auto to_x_x_zzz_yyz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 97);

        auto to_x_x_zzz_yzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 98);

        auto to_x_x_zzz_zzz = pbuffer.data(idx_op_geom_101_ff + 0 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_x_x_zzz_xxx, to_x_x_zzz_xxy, to_x_x_zzz_xxz, to_x_x_zzz_xyy, to_x_x_zzz_xyz, to_x_x_zzz_xzz, to_x_x_zzz_yyy, to_x_x_zzz_yyz, to_x_x_zzz_yzz, to_x_x_zzz_zzz, to_xzzz_xx, to_xzzz_xxxx, to_xzzz_xxxy, to_xzzz_xxxz, to_xzzz_xxyy, to_xzzz_xxyz, to_xzzz_xxzz, to_xzzz_xy, to_xzzz_xyyy, to_xzzz_xyyz, to_xzzz_xyzz, to_xzzz_xz, to_xzzz_xzzz, to_xzzz_yy, to_xzzz_yz, to_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zzz_xxx[k] = -6.0 * to_xzzz_xx[k] * tbe_0 + 4.0 * to_xzzz_xxxx[k] * tbe_0 * tke_0;

            to_x_x_zzz_xxy[k] = -4.0 * to_xzzz_xy[k] * tbe_0 + 4.0 * to_xzzz_xxxy[k] * tbe_0 * tke_0;

            to_x_x_zzz_xxz[k] = -4.0 * to_xzzz_xz[k] * tbe_0 + 4.0 * to_xzzz_xxxz[k] * tbe_0 * tke_0;

            to_x_x_zzz_xyy[k] = -2.0 * to_xzzz_yy[k] * tbe_0 + 4.0 * to_xzzz_xxyy[k] * tbe_0 * tke_0;

            to_x_x_zzz_xyz[k] = -2.0 * to_xzzz_yz[k] * tbe_0 + 4.0 * to_xzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_x_zzz_xzz[k] = -2.0 * to_xzzz_zz[k] * tbe_0 + 4.0 * to_xzzz_xxzz[k] * tbe_0 * tke_0;

            to_x_x_zzz_yyy[k] = 4.0 * to_xzzz_xyyy[k] * tbe_0 * tke_0;

            to_x_x_zzz_yyz[k] = 4.0 * to_xzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_x_zzz_yzz[k] = 4.0 * to_xzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_x_zzz_zzz[k] = 4.0 * to_xzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 100-110 components of targeted buffer : FF

        auto to_x_y_xxx_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 0);

        auto to_x_y_xxx_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 1);

        auto to_x_y_xxx_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 2);

        auto to_x_y_xxx_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 3);

        auto to_x_y_xxx_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 4);

        auto to_x_y_xxx_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 5);

        auto to_x_y_xxx_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 6);

        auto to_x_y_xxx_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 7);

        auto to_x_y_xxx_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 8);

        auto to_x_y_xxx_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_x_y_xxx_xxx, to_x_y_xxx_xxy, to_x_y_xxx_xxz, to_x_y_xxx_xyy, to_x_y_xxx_xyz, to_x_y_xxx_xzz, to_x_y_xxx_yyy, to_x_y_xxx_yyz, to_x_y_xxx_yzz, to_x_y_xxx_zzz, to_xx_xx, to_xx_xxxy, to_xx_xxyy, to_xx_xxyz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_yy, to_xx_yyyy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, to_xxxx_xx, to_xxxx_xxxy, to_xxxx_xxyy, to_xxxx_xxyz, to_xxxx_xy, to_xxxx_xyyy, to_xxxx_xyyz, to_xxxx_xyzz, to_xxxx_xz, to_xxxx_yy, to_xxxx_yyyy, to_xxxx_yyyz, to_xxxx_yyzz, to_xxxx_yz, to_xxxx_yzzz, to_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxx_xxx[k] = -6.0 * to_xx_xxxy[k] * tke_0 + 4.0 * to_xxxx_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xxx_xxy[k] = 3.0 * to_xx_xx[k] - 6.0 * to_xx_xxyy[k] * tke_0 - 2.0 * to_xxxx_xx[k] * tbe_0 + 4.0 * to_xxxx_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xxx_xxz[k] = -6.0 * to_xx_xxyz[k] * tke_0 + 4.0 * to_xxxx_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xxx_xyy[k] = 6.0 * to_xx_xy[k] - 6.0 * to_xx_xyyy[k] * tke_0 - 4.0 * to_xxxx_xy[k] * tbe_0 + 4.0 * to_xxxx_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xxx_xyz[k] = 3.0 * to_xx_xz[k] - 6.0 * to_xx_xyyz[k] * tke_0 - 2.0 * to_xxxx_xz[k] * tbe_0 + 4.0 * to_xxxx_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xxx_xzz[k] = -6.0 * to_xx_xyzz[k] * tke_0 + 4.0 * to_xxxx_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xxx_yyy[k] = 9.0 * to_xx_yy[k] - 6.0 * to_xx_yyyy[k] * tke_0 - 6.0 * to_xxxx_yy[k] * tbe_0 + 4.0 * to_xxxx_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xxx_yyz[k] = 6.0 * to_xx_yz[k] - 6.0 * to_xx_yyyz[k] * tke_0 - 4.0 * to_xxxx_yz[k] * tbe_0 + 4.0 * to_xxxx_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xxx_yzz[k] = 3.0 * to_xx_zz[k] - 6.0 * to_xx_yyzz[k] * tke_0 - 2.0 * to_xxxx_zz[k] * tbe_0 + 4.0 * to_xxxx_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xxx_zzz[k] = -6.0 * to_xx_yzzz[k] * tke_0 + 4.0 * to_xxxx_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 110-120 components of targeted buffer : FF

        auto to_x_y_xxy_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 10);

        auto to_x_y_xxy_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 11);

        auto to_x_y_xxy_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 12);

        auto to_x_y_xxy_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 13);

        auto to_x_y_xxy_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 14);

        auto to_x_y_xxy_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 15);

        auto to_x_y_xxy_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 16);

        auto to_x_y_xxy_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 17);

        auto to_x_y_xxy_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 18);

        auto to_x_y_xxy_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_x_y_xxy_xxx, to_x_y_xxy_xxy, to_x_y_xxy_xxz, to_x_y_xxy_xyy, to_x_y_xxy_xyz, to_x_y_xxy_xzz, to_x_y_xxy_yyy, to_x_y_xxy_yyz, to_x_y_xxy_yzz, to_x_y_xxy_zzz, to_xxxy_xx, to_xxxy_xxxy, to_xxxy_xxyy, to_xxxy_xxyz, to_xxxy_xy, to_xxxy_xyyy, to_xxxy_xyyz, to_xxxy_xyzz, to_xxxy_xz, to_xxxy_yy, to_xxxy_yyyy, to_xxxy_yyyz, to_xxxy_yyzz, to_xxxy_yz, to_xxxy_yzzz, to_xxxy_zz, to_xy_xx, to_xy_xxxy, to_xy_xxyy, to_xy_xxyz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_yy, to_xy_yyyy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxy_xxx[k] = -4.0 * to_xy_xxxy[k] * tke_0 + 4.0 * to_xxxy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xxy_xxy[k] = 2.0 * to_xy_xx[k] - 4.0 * to_xy_xxyy[k] * tke_0 - 2.0 * to_xxxy_xx[k] * tbe_0 + 4.0 * to_xxxy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xxy_xxz[k] = -4.0 * to_xy_xxyz[k] * tke_0 + 4.0 * to_xxxy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xxy_xyy[k] = 4.0 * to_xy_xy[k] - 4.0 * to_xy_xyyy[k] * tke_0 - 4.0 * to_xxxy_xy[k] * tbe_0 + 4.0 * to_xxxy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xxy_xyz[k] = 2.0 * to_xy_xz[k] - 4.0 * to_xy_xyyz[k] * tke_0 - 2.0 * to_xxxy_xz[k] * tbe_0 + 4.0 * to_xxxy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xxy_xzz[k] = -4.0 * to_xy_xyzz[k] * tke_0 + 4.0 * to_xxxy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xxy_yyy[k] = 6.0 * to_xy_yy[k] - 4.0 * to_xy_yyyy[k] * tke_0 - 6.0 * to_xxxy_yy[k] * tbe_0 + 4.0 * to_xxxy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xxy_yyz[k] = 4.0 * to_xy_yz[k] - 4.0 * to_xy_yyyz[k] * tke_0 - 4.0 * to_xxxy_yz[k] * tbe_0 + 4.0 * to_xxxy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xxy_yzz[k] = 2.0 * to_xy_zz[k] - 4.0 * to_xy_yyzz[k] * tke_0 - 2.0 * to_xxxy_zz[k] * tbe_0 + 4.0 * to_xxxy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xxy_zzz[k] = -4.0 * to_xy_yzzz[k] * tke_0 + 4.0 * to_xxxy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-130 components of targeted buffer : FF

        auto to_x_y_xxz_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 20);

        auto to_x_y_xxz_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 21);

        auto to_x_y_xxz_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 22);

        auto to_x_y_xxz_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 23);

        auto to_x_y_xxz_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 24);

        auto to_x_y_xxz_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 25);

        auto to_x_y_xxz_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 26);

        auto to_x_y_xxz_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 27);

        auto to_x_y_xxz_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 28);

        auto to_x_y_xxz_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_x_y_xxz_xxx, to_x_y_xxz_xxy, to_x_y_xxz_xxz, to_x_y_xxz_xyy, to_x_y_xxz_xyz, to_x_y_xxz_xzz, to_x_y_xxz_yyy, to_x_y_xxz_yyz, to_x_y_xxz_yzz, to_x_y_xxz_zzz, to_xxxz_xx, to_xxxz_xxxy, to_xxxz_xxyy, to_xxxz_xxyz, to_xxxz_xy, to_xxxz_xyyy, to_xxxz_xyyz, to_xxxz_xyzz, to_xxxz_xz, to_xxxz_yy, to_xxxz_yyyy, to_xxxz_yyyz, to_xxxz_yyzz, to_xxxz_yz, to_xxxz_yzzz, to_xxxz_zz, to_xz_xx, to_xz_xxxy, to_xz_xxyy, to_xz_xxyz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_yy, to_xz_yyyy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxz_xxx[k] = -4.0 * to_xz_xxxy[k] * tke_0 + 4.0 * to_xxxz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xxz_xxy[k] = 2.0 * to_xz_xx[k] - 4.0 * to_xz_xxyy[k] * tke_0 - 2.0 * to_xxxz_xx[k] * tbe_0 + 4.0 * to_xxxz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xxz_xxz[k] = -4.0 * to_xz_xxyz[k] * tke_0 + 4.0 * to_xxxz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xxz_xyy[k] = 4.0 * to_xz_xy[k] - 4.0 * to_xz_xyyy[k] * tke_0 - 4.0 * to_xxxz_xy[k] * tbe_0 + 4.0 * to_xxxz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xxz_xyz[k] = 2.0 * to_xz_xz[k] - 4.0 * to_xz_xyyz[k] * tke_0 - 2.0 * to_xxxz_xz[k] * tbe_0 + 4.0 * to_xxxz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xxz_xzz[k] = -4.0 * to_xz_xyzz[k] * tke_0 + 4.0 * to_xxxz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xxz_yyy[k] = 6.0 * to_xz_yy[k] - 4.0 * to_xz_yyyy[k] * tke_0 - 6.0 * to_xxxz_yy[k] * tbe_0 + 4.0 * to_xxxz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xxz_yyz[k] = 4.0 * to_xz_yz[k] - 4.0 * to_xz_yyyz[k] * tke_0 - 4.0 * to_xxxz_yz[k] * tbe_0 + 4.0 * to_xxxz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xxz_yzz[k] = 2.0 * to_xz_zz[k] - 4.0 * to_xz_yyzz[k] * tke_0 - 2.0 * to_xxxz_zz[k] * tbe_0 + 4.0 * to_xxxz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xxz_zzz[k] = -4.0 * to_xz_yzzz[k] * tke_0 + 4.0 * to_xxxz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 130-140 components of targeted buffer : FF

        auto to_x_y_xyy_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 30);

        auto to_x_y_xyy_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 31);

        auto to_x_y_xyy_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 32);

        auto to_x_y_xyy_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 33);

        auto to_x_y_xyy_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 34);

        auto to_x_y_xyy_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 35);

        auto to_x_y_xyy_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 36);

        auto to_x_y_xyy_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 37);

        auto to_x_y_xyy_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 38);

        auto to_x_y_xyy_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_x_y_xyy_xxx, to_x_y_xyy_xxy, to_x_y_xyy_xxz, to_x_y_xyy_xyy, to_x_y_xyy_xyz, to_x_y_xyy_xzz, to_x_y_xyy_yyy, to_x_y_xyy_yyz, to_x_y_xyy_yzz, to_x_y_xyy_zzz, to_xxyy_xx, to_xxyy_xxxy, to_xxyy_xxyy, to_xxyy_xxyz, to_xxyy_xy, to_xxyy_xyyy, to_xxyy_xyyz, to_xxyy_xyzz, to_xxyy_xz, to_xxyy_yy, to_xxyy_yyyy, to_xxyy_yyyz, to_xxyy_yyzz, to_xxyy_yz, to_xxyy_yzzz, to_xxyy_zz, to_yy_xx, to_yy_xxxy, to_yy_xxyy, to_yy_xxyz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_yy, to_yy_yyyy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyy_xxx[k] = -2.0 * to_yy_xxxy[k] * tke_0 + 4.0 * to_xxyy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xyy_xxy[k] = to_yy_xx[k] - 2.0 * to_yy_xxyy[k] * tke_0 - 2.0 * to_xxyy_xx[k] * tbe_0 + 4.0 * to_xxyy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xyy_xxz[k] = -2.0 * to_yy_xxyz[k] * tke_0 + 4.0 * to_xxyy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xyy_xyy[k] = 2.0 * to_yy_xy[k] - 2.0 * to_yy_xyyy[k] * tke_0 - 4.0 * to_xxyy_xy[k] * tbe_0 + 4.0 * to_xxyy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xyy_xyz[k] = to_yy_xz[k] - 2.0 * to_yy_xyyz[k] * tke_0 - 2.0 * to_xxyy_xz[k] * tbe_0 + 4.0 * to_xxyy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xyy_xzz[k] = -2.0 * to_yy_xyzz[k] * tke_0 + 4.0 * to_xxyy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xyy_yyy[k] = 3.0 * to_yy_yy[k] - 2.0 * to_yy_yyyy[k] * tke_0 - 6.0 * to_xxyy_yy[k] * tbe_0 + 4.0 * to_xxyy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xyy_yyz[k] = 2.0 * to_yy_yz[k] - 2.0 * to_yy_yyyz[k] * tke_0 - 4.0 * to_xxyy_yz[k] * tbe_0 + 4.0 * to_xxyy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xyy_yzz[k] = to_yy_zz[k] - 2.0 * to_yy_yyzz[k] * tke_0 - 2.0 * to_xxyy_zz[k] * tbe_0 + 4.0 * to_xxyy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xyy_zzz[k] = -2.0 * to_yy_yzzz[k] * tke_0 + 4.0 * to_xxyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 140-150 components of targeted buffer : FF

        auto to_x_y_xyz_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 40);

        auto to_x_y_xyz_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 41);

        auto to_x_y_xyz_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 42);

        auto to_x_y_xyz_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 43);

        auto to_x_y_xyz_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 44);

        auto to_x_y_xyz_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 45);

        auto to_x_y_xyz_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 46);

        auto to_x_y_xyz_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 47);

        auto to_x_y_xyz_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 48);

        auto to_x_y_xyz_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_x_y_xyz_xxx, to_x_y_xyz_xxy, to_x_y_xyz_xxz, to_x_y_xyz_xyy, to_x_y_xyz_xyz, to_x_y_xyz_xzz, to_x_y_xyz_yyy, to_x_y_xyz_yyz, to_x_y_xyz_yzz, to_x_y_xyz_zzz, to_xxyz_xx, to_xxyz_xxxy, to_xxyz_xxyy, to_xxyz_xxyz, to_xxyz_xy, to_xxyz_xyyy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_yy, to_xxyz_yyyy, to_xxyz_yyyz, to_xxyz_yyzz, to_xxyz_yz, to_xxyz_yzzz, to_xxyz_zz, to_yz_xx, to_yz_xxxy, to_yz_xxyy, to_yz_xxyz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_yy, to_yz_yyyy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyz_xxx[k] = -2.0 * to_yz_xxxy[k] * tke_0 + 4.0 * to_xxyz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xyz_xxy[k] = to_yz_xx[k] - 2.0 * to_yz_xxyy[k] * tke_0 - 2.0 * to_xxyz_xx[k] * tbe_0 + 4.0 * to_xxyz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xyz_xxz[k] = -2.0 * to_yz_xxyz[k] * tke_0 + 4.0 * to_xxyz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xyz_xyy[k] = 2.0 * to_yz_xy[k] - 2.0 * to_yz_xyyy[k] * tke_0 - 4.0 * to_xxyz_xy[k] * tbe_0 + 4.0 * to_xxyz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xyz_xyz[k] = to_yz_xz[k] - 2.0 * to_yz_xyyz[k] * tke_0 - 2.0 * to_xxyz_xz[k] * tbe_0 + 4.0 * to_xxyz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xyz_xzz[k] = -2.0 * to_yz_xyzz[k] * tke_0 + 4.0 * to_xxyz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xyz_yyy[k] = 3.0 * to_yz_yy[k] - 2.0 * to_yz_yyyy[k] * tke_0 - 6.0 * to_xxyz_yy[k] * tbe_0 + 4.0 * to_xxyz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xyz_yyz[k] = 2.0 * to_yz_yz[k] - 2.0 * to_yz_yyyz[k] * tke_0 - 4.0 * to_xxyz_yz[k] * tbe_0 + 4.0 * to_xxyz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xyz_yzz[k] = to_yz_zz[k] - 2.0 * to_yz_yyzz[k] * tke_0 - 2.0 * to_xxyz_zz[k] * tbe_0 + 4.0 * to_xxyz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xyz_zzz[k] = -2.0 * to_yz_yzzz[k] * tke_0 + 4.0 * to_xxyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-160 components of targeted buffer : FF

        auto to_x_y_xzz_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 50);

        auto to_x_y_xzz_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 51);

        auto to_x_y_xzz_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 52);

        auto to_x_y_xzz_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 53);

        auto to_x_y_xzz_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 54);

        auto to_x_y_xzz_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 55);

        auto to_x_y_xzz_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 56);

        auto to_x_y_xzz_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 57);

        auto to_x_y_xzz_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 58);

        auto to_x_y_xzz_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_x_y_xzz_xxx, to_x_y_xzz_xxy, to_x_y_xzz_xxz, to_x_y_xzz_xyy, to_x_y_xzz_xyz, to_x_y_xzz_xzz, to_x_y_xzz_yyy, to_x_y_xzz_yyz, to_x_y_xzz_yzz, to_x_y_xzz_zzz, to_xxzz_xx, to_xxzz_xxxy, to_xxzz_xxyy, to_xxzz_xxyz, to_xxzz_xy, to_xxzz_xyyy, to_xxzz_xyyz, to_xxzz_xyzz, to_xxzz_xz, to_xxzz_yy, to_xxzz_yyyy, to_xxzz_yyyz, to_xxzz_yyzz, to_xxzz_yz, to_xxzz_yzzz, to_xxzz_zz, to_zz_xx, to_zz_xxxy, to_zz_xxyy, to_zz_xxyz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_yy, to_zz_yyyy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xzz_xxx[k] = -2.0 * to_zz_xxxy[k] * tke_0 + 4.0 * to_xxzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_xzz_xxy[k] = to_zz_xx[k] - 2.0 * to_zz_xxyy[k] * tke_0 - 2.0 * to_xxzz_xx[k] * tbe_0 + 4.0 * to_xxzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_xzz_xxz[k] = -2.0 * to_zz_xxyz[k] * tke_0 + 4.0 * to_xxzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_xzz_xyy[k] = 2.0 * to_zz_xy[k] - 2.0 * to_zz_xyyy[k] * tke_0 - 4.0 * to_xxzz_xy[k] * tbe_0 + 4.0 * to_xxzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_xzz_xyz[k] = to_zz_xz[k] - 2.0 * to_zz_xyyz[k] * tke_0 - 2.0 * to_xxzz_xz[k] * tbe_0 + 4.0 * to_xxzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_xzz_xzz[k] = -2.0 * to_zz_xyzz[k] * tke_0 + 4.0 * to_xxzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_xzz_yyy[k] = 3.0 * to_zz_yy[k] - 2.0 * to_zz_yyyy[k] * tke_0 - 6.0 * to_xxzz_yy[k] * tbe_0 + 4.0 * to_xxzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_xzz_yyz[k] = 2.0 * to_zz_yz[k] - 2.0 * to_zz_yyyz[k] * tke_0 - 4.0 * to_xxzz_yz[k] * tbe_0 + 4.0 * to_xxzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_xzz_yzz[k] = to_zz_zz[k] - 2.0 * to_zz_yyzz[k] * tke_0 - 2.0 * to_xxzz_zz[k] * tbe_0 + 4.0 * to_xxzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_xzz_zzz[k] = -2.0 * to_zz_yzzz[k] * tke_0 + 4.0 * to_xxzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 160-170 components of targeted buffer : FF

        auto to_x_y_yyy_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 60);

        auto to_x_y_yyy_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 61);

        auto to_x_y_yyy_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 62);

        auto to_x_y_yyy_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 63);

        auto to_x_y_yyy_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 64);

        auto to_x_y_yyy_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 65);

        auto to_x_y_yyy_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 66);

        auto to_x_y_yyy_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 67);

        auto to_x_y_yyy_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 68);

        auto to_x_y_yyy_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_x_y_yyy_xxx, to_x_y_yyy_xxy, to_x_y_yyy_xxz, to_x_y_yyy_xyy, to_x_y_yyy_xyz, to_x_y_yyy_xzz, to_x_y_yyy_yyy, to_x_y_yyy_yyz, to_x_y_yyy_yzz, to_x_y_yyy_zzz, to_xyyy_xx, to_xyyy_xxxy, to_xyyy_xxyy, to_xyyy_xxyz, to_xyyy_xy, to_xyyy_xyyy, to_xyyy_xyyz, to_xyyy_xyzz, to_xyyy_xz, to_xyyy_yy, to_xyyy_yyyy, to_xyyy_yyyz, to_xyyy_yyzz, to_xyyy_yz, to_xyyy_yzzz, to_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyy_xxx[k] = 4.0 * to_xyyy_xxxy[k] * tbe_0 * tke_0;

            to_x_y_yyy_xxy[k] = -2.0 * to_xyyy_xx[k] * tbe_0 + 4.0 * to_xyyy_xxyy[k] * tbe_0 * tke_0;

            to_x_y_yyy_xxz[k] = 4.0 * to_xyyy_xxyz[k] * tbe_0 * tke_0;

            to_x_y_yyy_xyy[k] = -4.0 * to_xyyy_xy[k] * tbe_0 + 4.0 * to_xyyy_xyyy[k] * tbe_0 * tke_0;

            to_x_y_yyy_xyz[k] = -2.0 * to_xyyy_xz[k] * tbe_0 + 4.0 * to_xyyy_xyyz[k] * tbe_0 * tke_0;

            to_x_y_yyy_xzz[k] = 4.0 * to_xyyy_xyzz[k] * tbe_0 * tke_0;

            to_x_y_yyy_yyy[k] = -6.0 * to_xyyy_yy[k] * tbe_0 + 4.0 * to_xyyy_yyyy[k] * tbe_0 * tke_0;

            to_x_y_yyy_yyz[k] = -4.0 * to_xyyy_yz[k] * tbe_0 + 4.0 * to_xyyy_yyyz[k] * tbe_0 * tke_0;

            to_x_y_yyy_yzz[k] = -2.0 * to_xyyy_zz[k] * tbe_0 + 4.0 * to_xyyy_yyzz[k] * tbe_0 * tke_0;

            to_x_y_yyy_zzz[k] = 4.0 * to_xyyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 170-180 components of targeted buffer : FF

        auto to_x_y_yyz_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 70);

        auto to_x_y_yyz_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 71);

        auto to_x_y_yyz_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 72);

        auto to_x_y_yyz_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 73);

        auto to_x_y_yyz_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 74);

        auto to_x_y_yyz_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 75);

        auto to_x_y_yyz_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 76);

        auto to_x_y_yyz_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 77);

        auto to_x_y_yyz_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 78);

        auto to_x_y_yyz_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_x_y_yyz_xxx, to_x_y_yyz_xxy, to_x_y_yyz_xxz, to_x_y_yyz_xyy, to_x_y_yyz_xyz, to_x_y_yyz_xzz, to_x_y_yyz_yyy, to_x_y_yyz_yyz, to_x_y_yyz_yzz, to_x_y_yyz_zzz, to_xyyz_xx, to_xyyz_xxxy, to_xyyz_xxyy, to_xyyz_xxyz, to_xyyz_xy, to_xyyz_xyyy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_yy, to_xyyz_yyyy, to_xyyz_yyyz, to_xyyz_yyzz, to_xyyz_yz, to_xyyz_yzzz, to_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyz_xxx[k] = 4.0 * to_xyyz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_yyz_xxy[k] = -2.0 * to_xyyz_xx[k] * tbe_0 + 4.0 * to_xyyz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_yyz_xxz[k] = 4.0 * to_xyyz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_yyz_xyy[k] = -4.0 * to_xyyz_xy[k] * tbe_0 + 4.0 * to_xyyz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_yyz_xyz[k] = -2.0 * to_xyyz_xz[k] * tbe_0 + 4.0 * to_xyyz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_yyz_xzz[k] = 4.0 * to_xyyz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_yyz_yyy[k] = -6.0 * to_xyyz_yy[k] * tbe_0 + 4.0 * to_xyyz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_yyz_yyz[k] = -4.0 * to_xyyz_yz[k] * tbe_0 + 4.0 * to_xyyz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_yyz_yzz[k] = -2.0 * to_xyyz_zz[k] * tbe_0 + 4.0 * to_xyyz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_yyz_zzz[k] = 4.0 * to_xyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-190 components of targeted buffer : FF

        auto to_x_y_yzz_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 80);

        auto to_x_y_yzz_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 81);

        auto to_x_y_yzz_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 82);

        auto to_x_y_yzz_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 83);

        auto to_x_y_yzz_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 84);

        auto to_x_y_yzz_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 85);

        auto to_x_y_yzz_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 86);

        auto to_x_y_yzz_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 87);

        auto to_x_y_yzz_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 88);

        auto to_x_y_yzz_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_x_y_yzz_xxx, to_x_y_yzz_xxy, to_x_y_yzz_xxz, to_x_y_yzz_xyy, to_x_y_yzz_xyz, to_x_y_yzz_xzz, to_x_y_yzz_yyy, to_x_y_yzz_yyz, to_x_y_yzz_yzz, to_x_y_yzz_zzz, to_xyzz_xx, to_xyzz_xxxy, to_xyzz_xxyy, to_xyzz_xxyz, to_xyzz_xy, to_xyzz_xyyy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_yy, to_xyzz_yyyy, to_xyzz_yyyz, to_xyzz_yyzz, to_xyzz_yz, to_xyzz_yzzz, to_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yzz_xxx[k] = 4.0 * to_xyzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_yzz_xxy[k] = -2.0 * to_xyzz_xx[k] * tbe_0 + 4.0 * to_xyzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_yzz_xxz[k] = 4.0 * to_xyzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_yzz_xyy[k] = -4.0 * to_xyzz_xy[k] * tbe_0 + 4.0 * to_xyzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_yzz_xyz[k] = -2.0 * to_xyzz_xz[k] * tbe_0 + 4.0 * to_xyzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_yzz_xzz[k] = 4.0 * to_xyzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_yzz_yyy[k] = -6.0 * to_xyzz_yy[k] * tbe_0 + 4.0 * to_xyzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_yzz_yyz[k] = -4.0 * to_xyzz_yz[k] * tbe_0 + 4.0 * to_xyzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_yzz_yzz[k] = -2.0 * to_xyzz_zz[k] * tbe_0 + 4.0 * to_xyzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_yzz_zzz[k] = 4.0 * to_xyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 190-200 components of targeted buffer : FF

        auto to_x_y_zzz_xxx = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 90);

        auto to_x_y_zzz_xxy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 91);

        auto to_x_y_zzz_xxz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 92);

        auto to_x_y_zzz_xyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 93);

        auto to_x_y_zzz_xyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 94);

        auto to_x_y_zzz_xzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 95);

        auto to_x_y_zzz_yyy = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 96);

        auto to_x_y_zzz_yyz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 97);

        auto to_x_y_zzz_yzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 98);

        auto to_x_y_zzz_zzz = pbuffer.data(idx_op_geom_101_ff + 1 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_x_y_zzz_xxx, to_x_y_zzz_xxy, to_x_y_zzz_xxz, to_x_y_zzz_xyy, to_x_y_zzz_xyz, to_x_y_zzz_xzz, to_x_y_zzz_yyy, to_x_y_zzz_yyz, to_x_y_zzz_yzz, to_x_y_zzz_zzz, to_xzzz_xx, to_xzzz_xxxy, to_xzzz_xxyy, to_xzzz_xxyz, to_xzzz_xy, to_xzzz_xyyy, to_xzzz_xyyz, to_xzzz_xyzz, to_xzzz_xz, to_xzzz_yy, to_xzzz_yyyy, to_xzzz_yyyz, to_xzzz_yyzz, to_xzzz_yz, to_xzzz_yzzz, to_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zzz_xxx[k] = 4.0 * to_xzzz_xxxy[k] * tbe_0 * tke_0;

            to_x_y_zzz_xxy[k] = -2.0 * to_xzzz_xx[k] * tbe_0 + 4.0 * to_xzzz_xxyy[k] * tbe_0 * tke_0;

            to_x_y_zzz_xxz[k] = 4.0 * to_xzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_y_zzz_xyy[k] = -4.0 * to_xzzz_xy[k] * tbe_0 + 4.0 * to_xzzz_xyyy[k] * tbe_0 * tke_0;

            to_x_y_zzz_xyz[k] = -2.0 * to_xzzz_xz[k] * tbe_0 + 4.0 * to_xzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_y_zzz_xzz[k] = 4.0 * to_xzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_y_zzz_yyy[k] = -6.0 * to_xzzz_yy[k] * tbe_0 + 4.0 * to_xzzz_yyyy[k] * tbe_0 * tke_0;

            to_x_y_zzz_yyz[k] = -4.0 * to_xzzz_yz[k] * tbe_0 + 4.0 * to_xzzz_yyyz[k] * tbe_0 * tke_0;

            to_x_y_zzz_yzz[k] = -2.0 * to_xzzz_zz[k] * tbe_0 + 4.0 * to_xzzz_yyzz[k] * tbe_0 * tke_0;

            to_x_y_zzz_zzz[k] = 4.0 * to_xzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 200-210 components of targeted buffer : FF

        auto to_x_z_xxx_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 0);

        auto to_x_z_xxx_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 1);

        auto to_x_z_xxx_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 2);

        auto to_x_z_xxx_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 3);

        auto to_x_z_xxx_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 4);

        auto to_x_z_xxx_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 5);

        auto to_x_z_xxx_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 6);

        auto to_x_z_xxx_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 7);

        auto to_x_z_xxx_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 8);

        auto to_x_z_xxx_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_x_z_xxx_xxx, to_x_z_xxx_xxy, to_x_z_xxx_xxz, to_x_z_xxx_xyy, to_x_z_xxx_xyz, to_x_z_xxx_xzz, to_x_z_xxx_yyy, to_x_z_xxx_yyz, to_x_z_xxx_yzz, to_x_z_xxx_zzz, to_xx_xx, to_xx_xxxz, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, to_xx_zzzz, to_xxxx_xx, to_xxxx_xxxz, to_xxxx_xxyz, to_xxxx_xxzz, to_xxxx_xy, to_xxxx_xyyz, to_xxxx_xyzz, to_xxxx_xz, to_xxxx_xzzz, to_xxxx_yy, to_xxxx_yyyz, to_xxxx_yyzz, to_xxxx_yz, to_xxxx_yzzz, to_xxxx_zz, to_xxxx_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxx_xxx[k] = -6.0 * to_xx_xxxz[k] * tke_0 + 4.0 * to_xxxx_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xxy[k] = -6.0 * to_xx_xxyz[k] * tke_0 + 4.0 * to_xxxx_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xxz[k] = 3.0 * to_xx_xx[k] - 6.0 * to_xx_xxzz[k] * tke_0 - 2.0 * to_xxxx_xx[k] * tbe_0 + 4.0 * to_xxxx_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xyy[k] = -6.0 * to_xx_xyyz[k] * tke_0 + 4.0 * to_xxxx_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xyz[k] = 3.0 * to_xx_xy[k] - 6.0 * to_xx_xyzz[k] * tke_0 - 2.0 * to_xxxx_xy[k] * tbe_0 + 4.0 * to_xxxx_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_xzz[k] = 6.0 * to_xx_xz[k] - 6.0 * to_xx_xzzz[k] * tke_0 - 4.0 * to_xxxx_xz[k] * tbe_0 + 4.0 * to_xxxx_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_yyy[k] = -6.0 * to_xx_yyyz[k] * tke_0 + 4.0 * to_xxxx_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xxx_yyz[k] = 3.0 * to_xx_yy[k] - 6.0 * to_xx_yyzz[k] * tke_0 - 2.0 * to_xxxx_yy[k] * tbe_0 + 4.0 * to_xxxx_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_yzz[k] = 6.0 * to_xx_yz[k] - 6.0 * to_xx_yzzz[k] * tke_0 - 4.0 * to_xxxx_yz[k] * tbe_0 + 4.0 * to_xxxx_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xxx_zzz[k] = 9.0 * to_xx_zz[k] - 6.0 * to_xx_zzzz[k] * tke_0 - 6.0 * to_xxxx_zz[k] * tbe_0 + 4.0 * to_xxxx_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-220 components of targeted buffer : FF

        auto to_x_z_xxy_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 10);

        auto to_x_z_xxy_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 11);

        auto to_x_z_xxy_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 12);

        auto to_x_z_xxy_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 13);

        auto to_x_z_xxy_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 14);

        auto to_x_z_xxy_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 15);

        auto to_x_z_xxy_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 16);

        auto to_x_z_xxy_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 17);

        auto to_x_z_xxy_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 18);

        auto to_x_z_xxy_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_x_z_xxy_xxx, to_x_z_xxy_xxy, to_x_z_xxy_xxz, to_x_z_xxy_xyy, to_x_z_xxy_xyz, to_x_z_xxy_xzz, to_x_z_xxy_yyy, to_x_z_xxy_yyz, to_x_z_xxy_yzz, to_x_z_xxy_zzz, to_xxxy_xx, to_xxxy_xxxz, to_xxxy_xxyz, to_xxxy_xxzz, to_xxxy_xy, to_xxxy_xyyz, to_xxxy_xyzz, to_xxxy_xz, to_xxxy_xzzz, to_xxxy_yy, to_xxxy_yyyz, to_xxxy_yyzz, to_xxxy_yz, to_xxxy_yzzz, to_xxxy_zz, to_xxxy_zzzz, to_xy_xx, to_xy_xxxz, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, to_xy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxy_xxx[k] = -4.0 * to_xy_xxxz[k] * tke_0 + 4.0 * to_xxxy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xxy[k] = -4.0 * to_xy_xxyz[k] * tke_0 + 4.0 * to_xxxy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xxz[k] = 2.0 * to_xy_xx[k] - 4.0 * to_xy_xxzz[k] * tke_0 - 2.0 * to_xxxy_xx[k] * tbe_0 + 4.0 * to_xxxy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xyy[k] = -4.0 * to_xy_xyyz[k] * tke_0 + 4.0 * to_xxxy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xyz[k] = 2.0 * to_xy_xy[k] - 4.0 * to_xy_xyzz[k] * tke_0 - 2.0 * to_xxxy_xy[k] * tbe_0 + 4.0 * to_xxxy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_xzz[k] = 4.0 * to_xy_xz[k] - 4.0 * to_xy_xzzz[k] * tke_0 - 4.0 * to_xxxy_xz[k] * tbe_0 + 4.0 * to_xxxy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_yyy[k] = -4.0 * to_xy_yyyz[k] * tke_0 + 4.0 * to_xxxy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xxy_yyz[k] = 2.0 * to_xy_yy[k] - 4.0 * to_xy_yyzz[k] * tke_0 - 2.0 * to_xxxy_yy[k] * tbe_0 + 4.0 * to_xxxy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_yzz[k] = 4.0 * to_xy_yz[k] - 4.0 * to_xy_yzzz[k] * tke_0 - 4.0 * to_xxxy_yz[k] * tbe_0 + 4.0 * to_xxxy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xxy_zzz[k] = 6.0 * to_xy_zz[k] - 4.0 * to_xy_zzzz[k] * tke_0 - 6.0 * to_xxxy_zz[k] * tbe_0 + 4.0 * to_xxxy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 220-230 components of targeted buffer : FF

        auto to_x_z_xxz_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 20);

        auto to_x_z_xxz_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 21);

        auto to_x_z_xxz_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 22);

        auto to_x_z_xxz_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 23);

        auto to_x_z_xxz_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 24);

        auto to_x_z_xxz_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 25);

        auto to_x_z_xxz_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 26);

        auto to_x_z_xxz_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 27);

        auto to_x_z_xxz_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 28);

        auto to_x_z_xxz_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_x_z_xxz_xxx, to_x_z_xxz_xxy, to_x_z_xxz_xxz, to_x_z_xxz_xyy, to_x_z_xxz_xyz, to_x_z_xxz_xzz, to_x_z_xxz_yyy, to_x_z_xxz_yyz, to_x_z_xxz_yzz, to_x_z_xxz_zzz, to_xxxz_xx, to_xxxz_xxxz, to_xxxz_xxyz, to_xxxz_xxzz, to_xxxz_xy, to_xxxz_xyyz, to_xxxz_xyzz, to_xxxz_xz, to_xxxz_xzzz, to_xxxz_yy, to_xxxz_yyyz, to_xxxz_yyzz, to_xxxz_yz, to_xxxz_yzzz, to_xxxz_zz, to_xxxz_zzzz, to_xz_xx, to_xz_xxxz, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, to_xz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxz_xxx[k] = -4.0 * to_xz_xxxz[k] * tke_0 + 4.0 * to_xxxz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xxy[k] = -4.0 * to_xz_xxyz[k] * tke_0 + 4.0 * to_xxxz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xxz[k] = 2.0 * to_xz_xx[k] - 4.0 * to_xz_xxzz[k] * tke_0 - 2.0 * to_xxxz_xx[k] * tbe_0 + 4.0 * to_xxxz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xyy[k] = -4.0 * to_xz_xyyz[k] * tke_0 + 4.0 * to_xxxz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xyz[k] = 2.0 * to_xz_xy[k] - 4.0 * to_xz_xyzz[k] * tke_0 - 2.0 * to_xxxz_xy[k] * tbe_0 + 4.0 * to_xxxz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_xzz[k] = 4.0 * to_xz_xz[k] - 4.0 * to_xz_xzzz[k] * tke_0 - 4.0 * to_xxxz_xz[k] * tbe_0 + 4.0 * to_xxxz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_yyy[k] = -4.0 * to_xz_yyyz[k] * tke_0 + 4.0 * to_xxxz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xxz_yyz[k] = 2.0 * to_xz_yy[k] - 4.0 * to_xz_yyzz[k] * tke_0 - 2.0 * to_xxxz_yy[k] * tbe_0 + 4.0 * to_xxxz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_yzz[k] = 4.0 * to_xz_yz[k] - 4.0 * to_xz_yzzz[k] * tke_0 - 4.0 * to_xxxz_yz[k] * tbe_0 + 4.0 * to_xxxz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xxz_zzz[k] = 6.0 * to_xz_zz[k] - 4.0 * to_xz_zzzz[k] * tke_0 - 6.0 * to_xxxz_zz[k] * tbe_0 + 4.0 * to_xxxz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 230-240 components of targeted buffer : FF

        auto to_x_z_xyy_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 30);

        auto to_x_z_xyy_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 31);

        auto to_x_z_xyy_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 32);

        auto to_x_z_xyy_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 33);

        auto to_x_z_xyy_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 34);

        auto to_x_z_xyy_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 35);

        auto to_x_z_xyy_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 36);

        auto to_x_z_xyy_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 37);

        auto to_x_z_xyy_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 38);

        auto to_x_z_xyy_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_x_z_xyy_xxx, to_x_z_xyy_xxy, to_x_z_xyy_xxz, to_x_z_xyy_xyy, to_x_z_xyy_xyz, to_x_z_xyy_xzz, to_x_z_xyy_yyy, to_x_z_xyy_yyz, to_x_z_xyy_yzz, to_x_z_xyy_zzz, to_xxyy_xx, to_xxyy_xxxz, to_xxyy_xxyz, to_xxyy_xxzz, to_xxyy_xy, to_xxyy_xyyz, to_xxyy_xyzz, to_xxyy_xz, to_xxyy_xzzz, to_xxyy_yy, to_xxyy_yyyz, to_xxyy_yyzz, to_xxyy_yz, to_xxyy_yzzz, to_xxyy_zz, to_xxyy_zzzz, to_yy_xx, to_yy_xxxz, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, to_yy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyy_xxx[k] = -2.0 * to_yy_xxxz[k] * tke_0 + 4.0 * to_xxyy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xxy[k] = -2.0 * to_yy_xxyz[k] * tke_0 + 4.0 * to_xxyy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xxz[k] = to_yy_xx[k] - 2.0 * to_yy_xxzz[k] * tke_0 - 2.0 * to_xxyy_xx[k] * tbe_0 + 4.0 * to_xxyy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xyy[k] = -2.0 * to_yy_xyyz[k] * tke_0 + 4.0 * to_xxyy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xyz[k] = to_yy_xy[k] - 2.0 * to_yy_xyzz[k] * tke_0 - 2.0 * to_xxyy_xy[k] * tbe_0 + 4.0 * to_xxyy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_xzz[k] = 2.0 * to_yy_xz[k] - 2.0 * to_yy_xzzz[k] * tke_0 - 4.0 * to_xxyy_xz[k] * tbe_0 + 4.0 * to_xxyy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_yyy[k] = -2.0 * to_yy_yyyz[k] * tke_0 + 4.0 * to_xxyy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xyy_yyz[k] = to_yy_yy[k] - 2.0 * to_yy_yyzz[k] * tke_0 - 2.0 * to_xxyy_yy[k] * tbe_0 + 4.0 * to_xxyy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_yzz[k] = 2.0 * to_yy_yz[k] - 2.0 * to_yy_yzzz[k] * tke_0 - 4.0 * to_xxyy_yz[k] * tbe_0 + 4.0 * to_xxyy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xyy_zzz[k] = 3.0 * to_yy_zz[k] - 2.0 * to_yy_zzzz[k] * tke_0 - 6.0 * to_xxyy_zz[k] * tbe_0 + 4.0 * to_xxyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-250 components of targeted buffer : FF

        auto to_x_z_xyz_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 40);

        auto to_x_z_xyz_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 41);

        auto to_x_z_xyz_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 42);

        auto to_x_z_xyz_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 43);

        auto to_x_z_xyz_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 44);

        auto to_x_z_xyz_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 45);

        auto to_x_z_xyz_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 46);

        auto to_x_z_xyz_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 47);

        auto to_x_z_xyz_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 48);

        auto to_x_z_xyz_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_x_z_xyz_xxx, to_x_z_xyz_xxy, to_x_z_xyz_xxz, to_x_z_xyz_xyy, to_x_z_xyz_xyz, to_x_z_xyz_xzz, to_x_z_xyz_yyy, to_x_z_xyz_yyz, to_x_z_xyz_yzz, to_x_z_xyz_zzz, to_xxyz_xx, to_xxyz_xxxz, to_xxyz_xxyz, to_xxyz_xxzz, to_xxyz_xy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_xzzz, to_xxyz_yy, to_xxyz_yyyz, to_xxyz_yyzz, to_xxyz_yz, to_xxyz_yzzz, to_xxyz_zz, to_xxyz_zzzz, to_yz_xx, to_yz_xxxz, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, to_yz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyz_xxx[k] = -2.0 * to_yz_xxxz[k] * tke_0 + 4.0 * to_xxyz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xxy[k] = -2.0 * to_yz_xxyz[k] * tke_0 + 4.0 * to_xxyz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xxz[k] = to_yz_xx[k] - 2.0 * to_yz_xxzz[k] * tke_0 - 2.0 * to_xxyz_xx[k] * tbe_0 + 4.0 * to_xxyz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xyy[k] = -2.0 * to_yz_xyyz[k] * tke_0 + 4.0 * to_xxyz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xyz[k] = to_yz_xy[k] - 2.0 * to_yz_xyzz[k] * tke_0 - 2.0 * to_xxyz_xy[k] * tbe_0 + 4.0 * to_xxyz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_xzz[k] = 2.0 * to_yz_xz[k] - 2.0 * to_yz_xzzz[k] * tke_0 - 4.0 * to_xxyz_xz[k] * tbe_0 + 4.0 * to_xxyz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_yyy[k] = -2.0 * to_yz_yyyz[k] * tke_0 + 4.0 * to_xxyz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xyz_yyz[k] = to_yz_yy[k] - 2.0 * to_yz_yyzz[k] * tke_0 - 2.0 * to_xxyz_yy[k] * tbe_0 + 4.0 * to_xxyz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_yzz[k] = 2.0 * to_yz_yz[k] - 2.0 * to_yz_yzzz[k] * tke_0 - 4.0 * to_xxyz_yz[k] * tbe_0 + 4.0 * to_xxyz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xyz_zzz[k] = 3.0 * to_yz_zz[k] - 2.0 * to_yz_zzzz[k] * tke_0 - 6.0 * to_xxyz_zz[k] * tbe_0 + 4.0 * to_xxyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 250-260 components of targeted buffer : FF

        auto to_x_z_xzz_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 50);

        auto to_x_z_xzz_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 51);

        auto to_x_z_xzz_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 52);

        auto to_x_z_xzz_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 53);

        auto to_x_z_xzz_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 54);

        auto to_x_z_xzz_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 55);

        auto to_x_z_xzz_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 56);

        auto to_x_z_xzz_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 57);

        auto to_x_z_xzz_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 58);

        auto to_x_z_xzz_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_x_z_xzz_xxx, to_x_z_xzz_xxy, to_x_z_xzz_xxz, to_x_z_xzz_xyy, to_x_z_xzz_xyz, to_x_z_xzz_xzz, to_x_z_xzz_yyy, to_x_z_xzz_yyz, to_x_z_xzz_yzz, to_x_z_xzz_zzz, to_xxzz_xx, to_xxzz_xxxz, to_xxzz_xxyz, to_xxzz_xxzz, to_xxzz_xy, to_xxzz_xyyz, to_xxzz_xyzz, to_xxzz_xz, to_xxzz_xzzz, to_xxzz_yy, to_xxzz_yyyz, to_xxzz_yyzz, to_xxzz_yz, to_xxzz_yzzz, to_xxzz_zz, to_xxzz_zzzz, to_zz_xx, to_zz_xxxz, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, to_zz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xzz_xxx[k] = -2.0 * to_zz_xxxz[k] * tke_0 + 4.0 * to_xxzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xxy[k] = -2.0 * to_zz_xxyz[k] * tke_0 + 4.0 * to_xxzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xxz[k] = to_zz_xx[k] - 2.0 * to_zz_xxzz[k] * tke_0 - 2.0 * to_xxzz_xx[k] * tbe_0 + 4.0 * to_xxzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xyy[k] = -2.0 * to_zz_xyyz[k] * tke_0 + 4.0 * to_xxzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xyz[k] = to_zz_xy[k] - 2.0 * to_zz_xyzz[k] * tke_0 - 2.0 * to_xxzz_xy[k] * tbe_0 + 4.0 * to_xxzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_xzz[k] = 2.0 * to_zz_xz[k] - 2.0 * to_zz_xzzz[k] * tke_0 - 4.0 * to_xxzz_xz[k] * tbe_0 + 4.0 * to_xxzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_yyy[k] = -2.0 * to_zz_yyyz[k] * tke_0 + 4.0 * to_xxzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_xzz_yyz[k] = to_zz_yy[k] - 2.0 * to_zz_yyzz[k] * tke_0 - 2.0 * to_xxzz_yy[k] * tbe_0 + 4.0 * to_xxzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_yzz[k] = 2.0 * to_zz_yz[k] - 2.0 * to_zz_yzzz[k] * tke_0 - 4.0 * to_xxzz_yz[k] * tbe_0 + 4.0 * to_xxzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_xzz_zzz[k] = 3.0 * to_zz_zz[k] - 2.0 * to_zz_zzzz[k] * tke_0 - 6.0 * to_xxzz_zz[k] * tbe_0 + 4.0 * to_xxzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 260-270 components of targeted buffer : FF

        auto to_x_z_yyy_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 60);

        auto to_x_z_yyy_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 61);

        auto to_x_z_yyy_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 62);

        auto to_x_z_yyy_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 63);

        auto to_x_z_yyy_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 64);

        auto to_x_z_yyy_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 65);

        auto to_x_z_yyy_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 66);

        auto to_x_z_yyy_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 67);

        auto to_x_z_yyy_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 68);

        auto to_x_z_yyy_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_x_z_yyy_xxx, to_x_z_yyy_xxy, to_x_z_yyy_xxz, to_x_z_yyy_xyy, to_x_z_yyy_xyz, to_x_z_yyy_xzz, to_x_z_yyy_yyy, to_x_z_yyy_yyz, to_x_z_yyy_yzz, to_x_z_yyy_zzz, to_xyyy_xx, to_xyyy_xxxz, to_xyyy_xxyz, to_xyyy_xxzz, to_xyyy_xy, to_xyyy_xyyz, to_xyyy_xyzz, to_xyyy_xz, to_xyyy_xzzz, to_xyyy_yy, to_xyyy_yyyz, to_xyyy_yyzz, to_xyyy_yz, to_xyyy_yzzz, to_xyyy_zz, to_xyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyy_xxx[k] = 4.0 * to_xyyy_xxxz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xxy[k] = 4.0 * to_xyyy_xxyz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xxz[k] = -2.0 * to_xyyy_xx[k] * tbe_0 + 4.0 * to_xyyy_xxzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xyy[k] = 4.0 * to_xyyy_xyyz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xyz[k] = -2.0 * to_xyyy_xy[k] * tbe_0 + 4.0 * to_xyyy_xyzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_xzz[k] = -4.0 * to_xyyy_xz[k] * tbe_0 + 4.0 * to_xyyy_xzzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_yyy[k] = 4.0 * to_xyyy_yyyz[k] * tbe_0 * tke_0;

            to_x_z_yyy_yyz[k] = -2.0 * to_xyyy_yy[k] * tbe_0 + 4.0 * to_xyyy_yyzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_yzz[k] = -4.0 * to_xyyy_yz[k] * tbe_0 + 4.0 * to_xyyy_yzzz[k] * tbe_0 * tke_0;

            to_x_z_yyy_zzz[k] = -6.0 * to_xyyy_zz[k] * tbe_0 + 4.0 * to_xyyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-280 components of targeted buffer : FF

        auto to_x_z_yyz_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 70);

        auto to_x_z_yyz_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 71);

        auto to_x_z_yyz_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 72);

        auto to_x_z_yyz_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 73);

        auto to_x_z_yyz_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 74);

        auto to_x_z_yyz_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 75);

        auto to_x_z_yyz_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 76);

        auto to_x_z_yyz_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 77);

        auto to_x_z_yyz_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 78);

        auto to_x_z_yyz_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_x_z_yyz_xxx, to_x_z_yyz_xxy, to_x_z_yyz_xxz, to_x_z_yyz_xyy, to_x_z_yyz_xyz, to_x_z_yyz_xzz, to_x_z_yyz_yyy, to_x_z_yyz_yyz, to_x_z_yyz_yzz, to_x_z_yyz_zzz, to_xyyz_xx, to_xyyz_xxxz, to_xyyz_xxyz, to_xyyz_xxzz, to_xyyz_xy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_xzzz, to_xyyz_yy, to_xyyz_yyyz, to_xyyz_yyzz, to_xyyz_yz, to_xyyz_yzzz, to_xyyz_zz, to_xyyz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyz_xxx[k] = 4.0 * to_xyyz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xxy[k] = 4.0 * to_xyyz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xxz[k] = -2.0 * to_xyyz_xx[k] * tbe_0 + 4.0 * to_xyyz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xyy[k] = 4.0 * to_xyyz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xyz[k] = -2.0 * to_xyyz_xy[k] * tbe_0 + 4.0 * to_xyyz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_xzz[k] = -4.0 * to_xyyz_xz[k] * tbe_0 + 4.0 * to_xyyz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_yyy[k] = 4.0 * to_xyyz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_yyz_yyz[k] = -2.0 * to_xyyz_yy[k] * tbe_0 + 4.0 * to_xyyz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_yzz[k] = -4.0 * to_xyyz_yz[k] * tbe_0 + 4.0 * to_xyyz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_yyz_zzz[k] = -6.0 * to_xyyz_zz[k] * tbe_0 + 4.0 * to_xyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 280-290 components of targeted buffer : FF

        auto to_x_z_yzz_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 80);

        auto to_x_z_yzz_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 81);

        auto to_x_z_yzz_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 82);

        auto to_x_z_yzz_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 83);

        auto to_x_z_yzz_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 84);

        auto to_x_z_yzz_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 85);

        auto to_x_z_yzz_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 86);

        auto to_x_z_yzz_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 87);

        auto to_x_z_yzz_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 88);

        auto to_x_z_yzz_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_x_z_yzz_xxx, to_x_z_yzz_xxy, to_x_z_yzz_xxz, to_x_z_yzz_xyy, to_x_z_yzz_xyz, to_x_z_yzz_xzz, to_x_z_yzz_yyy, to_x_z_yzz_yyz, to_x_z_yzz_yzz, to_x_z_yzz_zzz, to_xyzz_xx, to_xyzz_xxxz, to_xyzz_xxyz, to_xyzz_xxzz, to_xyzz_xy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_xzzz, to_xyzz_yy, to_xyzz_yyyz, to_xyzz_yyzz, to_xyzz_yz, to_xyzz_yzzz, to_xyzz_zz, to_xyzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yzz_xxx[k] = 4.0 * to_xyzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xxy[k] = 4.0 * to_xyzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xxz[k] = -2.0 * to_xyzz_xx[k] * tbe_0 + 4.0 * to_xyzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xyy[k] = 4.0 * to_xyzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xyz[k] = -2.0 * to_xyzz_xy[k] * tbe_0 + 4.0 * to_xyzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_xzz[k] = -4.0 * to_xyzz_xz[k] * tbe_0 + 4.0 * to_xyzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_yyy[k] = 4.0 * to_xyzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_yzz_yyz[k] = -2.0 * to_xyzz_yy[k] * tbe_0 + 4.0 * to_xyzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_yzz[k] = -4.0 * to_xyzz_yz[k] * tbe_0 + 4.0 * to_xyzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_yzz_zzz[k] = -6.0 * to_xyzz_zz[k] * tbe_0 + 4.0 * to_xyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 290-300 components of targeted buffer : FF

        auto to_x_z_zzz_xxx = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 90);

        auto to_x_z_zzz_xxy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 91);

        auto to_x_z_zzz_xxz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 92);

        auto to_x_z_zzz_xyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 93);

        auto to_x_z_zzz_xyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 94);

        auto to_x_z_zzz_xzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 95);

        auto to_x_z_zzz_yyy = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 96);

        auto to_x_z_zzz_yyz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 97);

        auto to_x_z_zzz_yzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 98);

        auto to_x_z_zzz_zzz = pbuffer.data(idx_op_geom_101_ff + 2 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_x_z_zzz_xxx, to_x_z_zzz_xxy, to_x_z_zzz_xxz, to_x_z_zzz_xyy, to_x_z_zzz_xyz, to_x_z_zzz_xzz, to_x_z_zzz_yyy, to_x_z_zzz_yyz, to_x_z_zzz_yzz, to_x_z_zzz_zzz, to_xzzz_xx, to_xzzz_xxxz, to_xzzz_xxyz, to_xzzz_xxzz, to_xzzz_xy, to_xzzz_xyyz, to_xzzz_xyzz, to_xzzz_xz, to_xzzz_xzzz, to_xzzz_yy, to_xzzz_yyyz, to_xzzz_yyzz, to_xzzz_yz, to_xzzz_yzzz, to_xzzz_zz, to_xzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zzz_xxx[k] = 4.0 * to_xzzz_xxxz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xxy[k] = 4.0 * to_xzzz_xxyz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xxz[k] = -2.0 * to_xzzz_xx[k] * tbe_0 + 4.0 * to_xzzz_xxzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xyy[k] = 4.0 * to_xzzz_xyyz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xyz[k] = -2.0 * to_xzzz_xy[k] * tbe_0 + 4.0 * to_xzzz_xyzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_xzz[k] = -4.0 * to_xzzz_xz[k] * tbe_0 + 4.0 * to_xzzz_xzzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_yyy[k] = 4.0 * to_xzzz_yyyz[k] * tbe_0 * tke_0;

            to_x_z_zzz_yyz[k] = -2.0 * to_xzzz_yy[k] * tbe_0 + 4.0 * to_xzzz_yyzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_yzz[k] = -4.0 * to_xzzz_yz[k] * tbe_0 + 4.0 * to_xzzz_yzzz[k] * tbe_0 * tke_0;

            to_x_z_zzz_zzz[k] = -6.0 * to_xzzz_zz[k] * tbe_0 + 4.0 * to_xzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-310 components of targeted buffer : FF

        auto to_y_x_xxx_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 0);

        auto to_y_x_xxx_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 1);

        auto to_y_x_xxx_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 2);

        auto to_y_x_xxx_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 3);

        auto to_y_x_xxx_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 4);

        auto to_y_x_xxx_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 5);

        auto to_y_x_xxx_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 6);

        auto to_y_x_xxx_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 7);

        auto to_y_x_xxx_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 8);

        auto to_y_x_xxx_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_xxxy_xx, to_xxxy_xxxx, to_xxxy_xxxy, to_xxxy_xxxz, to_xxxy_xxyy, to_xxxy_xxyz, to_xxxy_xxzz, to_xxxy_xy, to_xxxy_xyyy, to_xxxy_xyyz, to_xxxy_xyzz, to_xxxy_xz, to_xxxy_xzzz, to_xxxy_yy, to_xxxy_yz, to_xxxy_zz, to_y_x_xxx_xxx, to_y_x_xxx_xxy, to_y_x_xxx_xxz, to_y_x_xxx_xyy, to_y_x_xxx_xyz, to_y_x_xxx_xzz, to_y_x_xxx_yyy, to_y_x_xxx_yyz, to_y_x_xxx_yzz, to_y_x_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxx_xxx[k] = -6.0 * to_xxxy_xx[k] * tbe_0 + 4.0 * to_xxxy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xxx_xxy[k] = -4.0 * to_xxxy_xy[k] * tbe_0 + 4.0 * to_xxxy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xxx_xxz[k] = -4.0 * to_xxxy_xz[k] * tbe_0 + 4.0 * to_xxxy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xxx_xyy[k] = -2.0 * to_xxxy_yy[k] * tbe_0 + 4.0 * to_xxxy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xxx_xyz[k] = -2.0 * to_xxxy_yz[k] * tbe_0 + 4.0 * to_xxxy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xxx_xzz[k] = -2.0 * to_xxxy_zz[k] * tbe_0 + 4.0 * to_xxxy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xxx_yyy[k] = 4.0 * to_xxxy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xxx_yyz[k] = 4.0 * to_xxxy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xxx_yzz[k] = 4.0 * to_xxxy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xxx_zzz[k] = 4.0 * to_xxxy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 310-320 components of targeted buffer : FF

        auto to_y_x_xxy_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 10);

        auto to_y_x_xxy_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 11);

        auto to_y_x_xxy_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 12);

        auto to_y_x_xxy_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 13);

        auto to_y_x_xxy_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 14);

        auto to_y_x_xxy_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 15);

        auto to_y_x_xxy_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 16);

        auto to_y_x_xxy_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 17);

        auto to_y_x_xxy_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 18);

        auto to_y_x_xxy_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_xx_xx, to_xx_xxxx, to_xx_xxxy, to_xx_xxxz, to_xx_xxyy, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yz, to_xx_zz, to_xxyy_xx, to_xxyy_xxxx, to_xxyy_xxxy, to_xxyy_xxxz, to_xxyy_xxyy, to_xxyy_xxyz, to_xxyy_xxzz, to_xxyy_xy, to_xxyy_xyyy, to_xxyy_xyyz, to_xxyy_xyzz, to_xxyy_xz, to_xxyy_xzzz, to_xxyy_yy, to_xxyy_yz, to_xxyy_zz, to_y_x_xxy_xxx, to_y_x_xxy_xxy, to_y_x_xxy_xxz, to_y_x_xxy_xyy, to_y_x_xxy_xyz, to_y_x_xxy_xzz, to_y_x_xxy_yyy, to_y_x_xxy_yyz, to_y_x_xxy_yzz, to_y_x_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxy_xxx[k] = 3.0 * to_xx_xx[k] - 2.0 * to_xx_xxxx[k] * tke_0 - 6.0 * to_xxyy_xx[k] * tbe_0 + 4.0 * to_xxyy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xxy_xxy[k] = 2.0 * to_xx_xy[k] - 2.0 * to_xx_xxxy[k] * tke_0 - 4.0 * to_xxyy_xy[k] * tbe_0 + 4.0 * to_xxyy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xxy_xxz[k] = 2.0 * to_xx_xz[k] - 2.0 * to_xx_xxxz[k] * tke_0 - 4.0 * to_xxyy_xz[k] * tbe_0 + 4.0 * to_xxyy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xxy_xyy[k] = to_xx_yy[k] - 2.0 * to_xx_xxyy[k] * tke_0 - 2.0 * to_xxyy_yy[k] * tbe_0 + 4.0 * to_xxyy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xxy_xyz[k] = to_xx_yz[k] - 2.0 * to_xx_xxyz[k] * tke_0 - 2.0 * to_xxyy_yz[k] * tbe_0 + 4.0 * to_xxyy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xxy_xzz[k] = to_xx_zz[k] - 2.0 * to_xx_xxzz[k] * tke_0 - 2.0 * to_xxyy_zz[k] * tbe_0 + 4.0 * to_xxyy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xxy_yyy[k] = -2.0 * to_xx_xyyy[k] * tke_0 + 4.0 * to_xxyy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xxy_yyz[k] = -2.0 * to_xx_xyyz[k] * tke_0 + 4.0 * to_xxyy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xxy_yzz[k] = -2.0 * to_xx_xyzz[k] * tke_0 + 4.0 * to_xxyy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xxy_zzz[k] = -2.0 * to_xx_xzzz[k] * tke_0 + 4.0 * to_xxyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 320-330 components of targeted buffer : FF

        auto to_y_x_xxz_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 20);

        auto to_y_x_xxz_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 21);

        auto to_y_x_xxz_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 22);

        auto to_y_x_xxz_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 23);

        auto to_y_x_xxz_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 24);

        auto to_y_x_xxz_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 25);

        auto to_y_x_xxz_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 26);

        auto to_y_x_xxz_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 27);

        auto to_y_x_xxz_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 28);

        auto to_y_x_xxz_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_xxyz_xx, to_xxyz_xxxx, to_xxyz_xxxy, to_xxyz_xxxz, to_xxyz_xxyy, to_xxyz_xxyz, to_xxyz_xxzz, to_xxyz_xy, to_xxyz_xyyy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_xzzz, to_xxyz_yy, to_xxyz_yz, to_xxyz_zz, to_y_x_xxz_xxx, to_y_x_xxz_xxy, to_y_x_xxz_xxz, to_y_x_xxz_xyy, to_y_x_xxz_xyz, to_y_x_xxz_xzz, to_y_x_xxz_yyy, to_y_x_xxz_yyz, to_y_x_xxz_yzz, to_y_x_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxz_xxx[k] = -6.0 * to_xxyz_xx[k] * tbe_0 + 4.0 * to_xxyz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xxz_xxy[k] = -4.0 * to_xxyz_xy[k] * tbe_0 + 4.0 * to_xxyz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xxz_xxz[k] = -4.0 * to_xxyz_xz[k] * tbe_0 + 4.0 * to_xxyz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xxz_xyy[k] = -2.0 * to_xxyz_yy[k] * tbe_0 + 4.0 * to_xxyz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xxz_xyz[k] = -2.0 * to_xxyz_yz[k] * tbe_0 + 4.0 * to_xxyz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xxz_xzz[k] = -2.0 * to_xxyz_zz[k] * tbe_0 + 4.0 * to_xxyz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xxz_yyy[k] = 4.0 * to_xxyz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xxz_yyz[k] = 4.0 * to_xxyz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xxz_yzz[k] = 4.0 * to_xxyz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xxz_zzz[k] = 4.0 * to_xxyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 330-340 components of targeted buffer : FF

        auto to_y_x_xyy_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 30);

        auto to_y_x_xyy_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 31);

        auto to_y_x_xyy_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 32);

        auto to_y_x_xyy_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 33);

        auto to_y_x_xyy_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 34);

        auto to_y_x_xyy_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 35);

        auto to_y_x_xyy_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 36);

        auto to_y_x_xyy_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 37);

        auto to_y_x_xyy_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 38);

        auto to_y_x_xyy_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_xy_xx, to_xy_xxxx, to_xy_xxxy, to_xy_xxxz, to_xy_xxyy, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yz, to_xy_zz, to_xyyy_xx, to_xyyy_xxxx, to_xyyy_xxxy, to_xyyy_xxxz, to_xyyy_xxyy, to_xyyy_xxyz, to_xyyy_xxzz, to_xyyy_xy, to_xyyy_xyyy, to_xyyy_xyyz, to_xyyy_xyzz, to_xyyy_xz, to_xyyy_xzzz, to_xyyy_yy, to_xyyy_yz, to_xyyy_zz, to_y_x_xyy_xxx, to_y_x_xyy_xxy, to_y_x_xyy_xxz, to_y_x_xyy_xyy, to_y_x_xyy_xyz, to_y_x_xyy_xzz, to_y_x_xyy_yyy, to_y_x_xyy_yyz, to_y_x_xyy_yzz, to_y_x_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyy_xxx[k] = 6.0 * to_xy_xx[k] - 4.0 * to_xy_xxxx[k] * tke_0 - 6.0 * to_xyyy_xx[k] * tbe_0 + 4.0 * to_xyyy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xyy_xxy[k] = 4.0 * to_xy_xy[k] - 4.0 * to_xy_xxxy[k] * tke_0 - 4.0 * to_xyyy_xy[k] * tbe_0 + 4.0 * to_xyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xyy_xxz[k] = 4.0 * to_xy_xz[k] - 4.0 * to_xy_xxxz[k] * tke_0 - 4.0 * to_xyyy_xz[k] * tbe_0 + 4.0 * to_xyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xyy_xyy[k] = 2.0 * to_xy_yy[k] - 4.0 * to_xy_xxyy[k] * tke_0 - 2.0 * to_xyyy_yy[k] * tbe_0 + 4.0 * to_xyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xyy_xyz[k] = 2.0 * to_xy_yz[k] - 4.0 * to_xy_xxyz[k] * tke_0 - 2.0 * to_xyyy_yz[k] * tbe_0 + 4.0 * to_xyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xyy_xzz[k] = 2.0 * to_xy_zz[k] - 4.0 * to_xy_xxzz[k] * tke_0 - 2.0 * to_xyyy_zz[k] * tbe_0 + 4.0 * to_xyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xyy_yyy[k] = -4.0 * to_xy_xyyy[k] * tke_0 + 4.0 * to_xyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xyy_yyz[k] = -4.0 * to_xy_xyyz[k] * tke_0 + 4.0 * to_xyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xyy_yzz[k] = -4.0 * to_xy_xyzz[k] * tke_0 + 4.0 * to_xyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xyy_zzz[k] = -4.0 * to_xy_xzzz[k] * tke_0 + 4.0 * to_xyyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 340-350 components of targeted buffer : FF

        auto to_y_x_xyz_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 40);

        auto to_y_x_xyz_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 41);

        auto to_y_x_xyz_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 42);

        auto to_y_x_xyz_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 43);

        auto to_y_x_xyz_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 44);

        auto to_y_x_xyz_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 45);

        auto to_y_x_xyz_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 46);

        auto to_y_x_xyz_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 47);

        auto to_y_x_xyz_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 48);

        auto to_y_x_xyz_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_xyyz_xx, to_xyyz_xxxx, to_xyyz_xxxy, to_xyyz_xxxz, to_xyyz_xxyy, to_xyyz_xxyz, to_xyyz_xxzz, to_xyyz_xy, to_xyyz_xyyy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_xzzz, to_xyyz_yy, to_xyyz_yz, to_xyyz_zz, to_xz_xx, to_xz_xxxx, to_xz_xxxy, to_xz_xxxz, to_xz_xxyy, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yz, to_xz_zz, to_y_x_xyz_xxx, to_y_x_xyz_xxy, to_y_x_xyz_xxz, to_y_x_xyz_xyy, to_y_x_xyz_xyz, to_y_x_xyz_xzz, to_y_x_xyz_yyy, to_y_x_xyz_yyz, to_y_x_xyz_yzz, to_y_x_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyz_xxx[k] = 3.0 * to_xz_xx[k] - 2.0 * to_xz_xxxx[k] * tke_0 - 6.0 * to_xyyz_xx[k] * tbe_0 + 4.0 * to_xyyz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xyz_xxy[k] = 2.0 * to_xz_xy[k] - 2.0 * to_xz_xxxy[k] * tke_0 - 4.0 * to_xyyz_xy[k] * tbe_0 + 4.0 * to_xyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xyz_xxz[k] = 2.0 * to_xz_xz[k] - 2.0 * to_xz_xxxz[k] * tke_0 - 4.0 * to_xyyz_xz[k] * tbe_0 + 4.0 * to_xyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xyz_xyy[k] = to_xz_yy[k] - 2.0 * to_xz_xxyy[k] * tke_0 - 2.0 * to_xyyz_yy[k] * tbe_0 + 4.0 * to_xyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xyz_xyz[k] = to_xz_yz[k] - 2.0 * to_xz_xxyz[k] * tke_0 - 2.0 * to_xyyz_yz[k] * tbe_0 + 4.0 * to_xyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xyz_xzz[k] = to_xz_zz[k] - 2.0 * to_xz_xxzz[k] * tke_0 - 2.0 * to_xyyz_zz[k] * tbe_0 + 4.0 * to_xyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xyz_yyy[k] = -2.0 * to_xz_xyyy[k] * tke_0 + 4.0 * to_xyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xyz_yyz[k] = -2.0 * to_xz_xyyz[k] * tke_0 + 4.0 * to_xyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xyz_yzz[k] = -2.0 * to_xz_xyzz[k] * tke_0 + 4.0 * to_xyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xyz_zzz[k] = -2.0 * to_xz_xzzz[k] * tke_0 + 4.0 * to_xyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 350-360 components of targeted buffer : FF

        auto to_y_x_xzz_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 50);

        auto to_y_x_xzz_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 51);

        auto to_y_x_xzz_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 52);

        auto to_y_x_xzz_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 53);

        auto to_y_x_xzz_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 54);

        auto to_y_x_xzz_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 55);

        auto to_y_x_xzz_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 56);

        auto to_y_x_xzz_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 57);

        auto to_y_x_xzz_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 58);

        auto to_y_x_xzz_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_xyzz_xx, to_xyzz_xxxx, to_xyzz_xxxy, to_xyzz_xxxz, to_xyzz_xxyy, to_xyzz_xxyz, to_xyzz_xxzz, to_xyzz_xy, to_xyzz_xyyy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_xzzz, to_xyzz_yy, to_xyzz_yz, to_xyzz_zz, to_y_x_xzz_xxx, to_y_x_xzz_xxy, to_y_x_xzz_xxz, to_y_x_xzz_xyy, to_y_x_xzz_xyz, to_y_x_xzz_xzz, to_y_x_xzz_yyy, to_y_x_xzz_yyz, to_y_x_xzz_yzz, to_y_x_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xzz_xxx[k] = -6.0 * to_xyzz_xx[k] * tbe_0 + 4.0 * to_xyzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_xzz_xxy[k] = -4.0 * to_xyzz_xy[k] * tbe_0 + 4.0 * to_xyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_xzz_xxz[k] = -4.0 * to_xyzz_xz[k] * tbe_0 + 4.0 * to_xyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_xzz_xyy[k] = -2.0 * to_xyzz_yy[k] * tbe_0 + 4.0 * to_xyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_xzz_xyz[k] = -2.0 * to_xyzz_yz[k] * tbe_0 + 4.0 * to_xyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_xzz_xzz[k] = -2.0 * to_xyzz_zz[k] * tbe_0 + 4.0 * to_xyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_xzz_yyy[k] = 4.0 * to_xyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_xzz_yyz[k] = 4.0 * to_xyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_xzz_yzz[k] = 4.0 * to_xyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_xzz_zzz[k] = 4.0 * to_xyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 360-370 components of targeted buffer : FF

        auto to_y_x_yyy_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 60);

        auto to_y_x_yyy_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 61);

        auto to_y_x_yyy_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 62);

        auto to_y_x_yyy_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 63);

        auto to_y_x_yyy_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 64);

        auto to_y_x_yyy_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 65);

        auto to_y_x_yyy_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 66);

        auto to_y_x_yyy_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 67);

        auto to_y_x_yyy_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 68);

        auto to_y_x_yyy_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_y_x_yyy_xxx, to_y_x_yyy_xxy, to_y_x_yyy_xxz, to_y_x_yyy_xyy, to_y_x_yyy_xyz, to_y_x_yyy_xzz, to_y_x_yyy_yyy, to_y_x_yyy_yyz, to_y_x_yyy_yzz, to_y_x_yyy_zzz, to_yy_xx, to_yy_xxxx, to_yy_xxxy, to_yy_xxxz, to_yy_xxyy, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yz, to_yy_zz, to_yyyy_xx, to_yyyy_xxxx, to_yyyy_xxxy, to_yyyy_xxxz, to_yyyy_xxyy, to_yyyy_xxyz, to_yyyy_xxzz, to_yyyy_xy, to_yyyy_xyyy, to_yyyy_xyyz, to_yyyy_xyzz, to_yyyy_xz, to_yyyy_xzzz, to_yyyy_yy, to_yyyy_yz, to_yyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyy_xxx[k] = 9.0 * to_yy_xx[k] - 6.0 * to_yy_xxxx[k] * tke_0 - 6.0 * to_yyyy_xx[k] * tbe_0 + 4.0 * to_yyyy_xxxx[k] * tbe_0 * tke_0;

            to_y_x_yyy_xxy[k] = 6.0 * to_yy_xy[k] - 6.0 * to_yy_xxxy[k] * tke_0 - 4.0 * to_yyyy_xy[k] * tbe_0 + 4.0 * to_yyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_x_yyy_xxz[k] = 6.0 * to_yy_xz[k] - 6.0 * to_yy_xxxz[k] * tke_0 - 4.0 * to_yyyy_xz[k] * tbe_0 + 4.0 * to_yyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_x_yyy_xyy[k] = 3.0 * to_yy_yy[k] - 6.0 * to_yy_xxyy[k] * tke_0 - 2.0 * to_yyyy_yy[k] * tbe_0 + 4.0 * to_yyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_x_yyy_xyz[k] = 3.0 * to_yy_yz[k] - 6.0 * to_yy_xxyz[k] * tke_0 - 2.0 * to_yyyy_yz[k] * tbe_0 + 4.0 * to_yyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_x_yyy_xzz[k] = 3.0 * to_yy_zz[k] - 6.0 * to_yy_xxzz[k] * tke_0 - 2.0 * to_yyyy_zz[k] * tbe_0 + 4.0 * to_yyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_x_yyy_yyy[k] = -6.0 * to_yy_xyyy[k] * tke_0 + 4.0 * to_yyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_x_yyy_yyz[k] = -6.0 * to_yy_xyyz[k] * tke_0 + 4.0 * to_yyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_x_yyy_yzz[k] = -6.0 * to_yy_xyzz[k] * tke_0 + 4.0 * to_yyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_x_yyy_zzz[k] = -6.0 * to_yy_xzzz[k] * tke_0 + 4.0 * to_yyyy_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 370-380 components of targeted buffer : FF

        auto to_y_x_yyz_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 70);

        auto to_y_x_yyz_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 71);

        auto to_y_x_yyz_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 72);

        auto to_y_x_yyz_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 73);

        auto to_y_x_yyz_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 74);

        auto to_y_x_yyz_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 75);

        auto to_y_x_yyz_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 76);

        auto to_y_x_yyz_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 77);

        auto to_y_x_yyz_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 78);

        auto to_y_x_yyz_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_y_x_yyz_xxx, to_y_x_yyz_xxy, to_y_x_yyz_xxz, to_y_x_yyz_xyy, to_y_x_yyz_xyz, to_y_x_yyz_xzz, to_y_x_yyz_yyy, to_y_x_yyz_yyz, to_y_x_yyz_yzz, to_y_x_yyz_zzz, to_yyyz_xx, to_yyyz_xxxx, to_yyyz_xxxy, to_yyyz_xxxz, to_yyyz_xxyy, to_yyyz_xxyz, to_yyyz_xxzz, to_yyyz_xy, to_yyyz_xyyy, to_yyyz_xyyz, to_yyyz_xyzz, to_yyyz_xz, to_yyyz_xzzz, to_yyyz_yy, to_yyyz_yz, to_yyyz_zz, to_yz_xx, to_yz_xxxx, to_yz_xxxy, to_yz_xxxz, to_yz_xxyy, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyz_xxx[k] = 6.0 * to_yz_xx[k] - 4.0 * to_yz_xxxx[k] * tke_0 - 6.0 * to_yyyz_xx[k] * tbe_0 + 4.0 * to_yyyz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_yyz_xxy[k] = 4.0 * to_yz_xy[k] - 4.0 * to_yz_xxxy[k] * tke_0 - 4.0 * to_yyyz_xy[k] * tbe_0 + 4.0 * to_yyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_yyz_xxz[k] = 4.0 * to_yz_xz[k] - 4.0 * to_yz_xxxz[k] * tke_0 - 4.0 * to_yyyz_xz[k] * tbe_0 + 4.0 * to_yyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_yyz_xyy[k] = 2.0 * to_yz_yy[k] - 4.0 * to_yz_xxyy[k] * tke_0 - 2.0 * to_yyyz_yy[k] * tbe_0 + 4.0 * to_yyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_yyz_xyz[k] = 2.0 * to_yz_yz[k] - 4.0 * to_yz_xxyz[k] * tke_0 - 2.0 * to_yyyz_yz[k] * tbe_0 + 4.0 * to_yyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_yyz_xzz[k] = 2.0 * to_yz_zz[k] - 4.0 * to_yz_xxzz[k] * tke_0 - 2.0 * to_yyyz_zz[k] * tbe_0 + 4.0 * to_yyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_yyz_yyy[k] = -4.0 * to_yz_xyyy[k] * tke_0 + 4.0 * to_yyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_yyz_yyz[k] = -4.0 * to_yz_xyyz[k] * tke_0 + 4.0 * to_yyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_yyz_yzz[k] = -4.0 * to_yz_xyzz[k] * tke_0 + 4.0 * to_yyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_yyz_zzz[k] = -4.0 * to_yz_xzzz[k] * tke_0 + 4.0 * to_yyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 380-390 components of targeted buffer : FF

        auto to_y_x_yzz_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 80);

        auto to_y_x_yzz_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 81);

        auto to_y_x_yzz_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 82);

        auto to_y_x_yzz_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 83);

        auto to_y_x_yzz_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 84);

        auto to_y_x_yzz_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 85);

        auto to_y_x_yzz_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 86);

        auto to_y_x_yzz_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 87);

        auto to_y_x_yzz_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 88);

        auto to_y_x_yzz_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_y_x_yzz_xxx, to_y_x_yzz_xxy, to_y_x_yzz_xxz, to_y_x_yzz_xyy, to_y_x_yzz_xyz, to_y_x_yzz_xzz, to_y_x_yzz_yyy, to_y_x_yzz_yyz, to_y_x_yzz_yzz, to_y_x_yzz_zzz, to_yyzz_xx, to_yyzz_xxxx, to_yyzz_xxxy, to_yyzz_xxxz, to_yyzz_xxyy, to_yyzz_xxyz, to_yyzz_xxzz, to_yyzz_xy, to_yyzz_xyyy, to_yyzz_xyyz, to_yyzz_xyzz, to_yyzz_xz, to_yyzz_xzzz, to_yyzz_yy, to_yyzz_yz, to_yyzz_zz, to_zz_xx, to_zz_xxxx, to_zz_xxxy, to_zz_xxxz, to_zz_xxyy, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yzz_xxx[k] = 3.0 * to_zz_xx[k] - 2.0 * to_zz_xxxx[k] * tke_0 - 6.0 * to_yyzz_xx[k] * tbe_0 + 4.0 * to_yyzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_yzz_xxy[k] = 2.0 * to_zz_xy[k] - 2.0 * to_zz_xxxy[k] * tke_0 - 4.0 * to_yyzz_xy[k] * tbe_0 + 4.0 * to_yyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_yzz_xxz[k] = 2.0 * to_zz_xz[k] - 2.0 * to_zz_xxxz[k] * tke_0 - 4.0 * to_yyzz_xz[k] * tbe_0 + 4.0 * to_yyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_yzz_xyy[k] = to_zz_yy[k] - 2.0 * to_zz_xxyy[k] * tke_0 - 2.0 * to_yyzz_yy[k] * tbe_0 + 4.0 * to_yyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_yzz_xyz[k] = to_zz_yz[k] - 2.0 * to_zz_xxyz[k] * tke_0 - 2.0 * to_yyzz_yz[k] * tbe_0 + 4.0 * to_yyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_yzz_xzz[k] = to_zz_zz[k] - 2.0 * to_zz_xxzz[k] * tke_0 - 2.0 * to_yyzz_zz[k] * tbe_0 + 4.0 * to_yyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_yzz_yyy[k] = -2.0 * to_zz_xyyy[k] * tke_0 + 4.0 * to_yyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_yzz_yyz[k] = -2.0 * to_zz_xyyz[k] * tke_0 + 4.0 * to_yyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_yzz_yzz[k] = -2.0 * to_zz_xyzz[k] * tke_0 + 4.0 * to_yyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_yzz_zzz[k] = -2.0 * to_zz_xzzz[k] * tke_0 + 4.0 * to_yyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 390-400 components of targeted buffer : FF

        auto to_y_x_zzz_xxx = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 90);

        auto to_y_x_zzz_xxy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 91);

        auto to_y_x_zzz_xxz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 92);

        auto to_y_x_zzz_xyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 93);

        auto to_y_x_zzz_xyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 94);

        auto to_y_x_zzz_xzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 95);

        auto to_y_x_zzz_yyy = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 96);

        auto to_y_x_zzz_yyz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 97);

        auto to_y_x_zzz_yzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 98);

        auto to_y_x_zzz_zzz = pbuffer.data(idx_op_geom_101_ff + 3 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_y_x_zzz_xxx, to_y_x_zzz_xxy, to_y_x_zzz_xxz, to_y_x_zzz_xyy, to_y_x_zzz_xyz, to_y_x_zzz_xzz, to_y_x_zzz_yyy, to_y_x_zzz_yyz, to_y_x_zzz_yzz, to_y_x_zzz_zzz, to_yzzz_xx, to_yzzz_xxxx, to_yzzz_xxxy, to_yzzz_xxxz, to_yzzz_xxyy, to_yzzz_xxyz, to_yzzz_xxzz, to_yzzz_xy, to_yzzz_xyyy, to_yzzz_xyyz, to_yzzz_xyzz, to_yzzz_xz, to_yzzz_xzzz, to_yzzz_yy, to_yzzz_yz, to_yzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zzz_xxx[k] = -6.0 * to_yzzz_xx[k] * tbe_0 + 4.0 * to_yzzz_xxxx[k] * tbe_0 * tke_0;

            to_y_x_zzz_xxy[k] = -4.0 * to_yzzz_xy[k] * tbe_0 + 4.0 * to_yzzz_xxxy[k] * tbe_0 * tke_0;

            to_y_x_zzz_xxz[k] = -4.0 * to_yzzz_xz[k] * tbe_0 + 4.0 * to_yzzz_xxxz[k] * tbe_0 * tke_0;

            to_y_x_zzz_xyy[k] = -2.0 * to_yzzz_yy[k] * tbe_0 + 4.0 * to_yzzz_xxyy[k] * tbe_0 * tke_0;

            to_y_x_zzz_xyz[k] = -2.0 * to_yzzz_yz[k] * tbe_0 + 4.0 * to_yzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_x_zzz_xzz[k] = -2.0 * to_yzzz_zz[k] * tbe_0 + 4.0 * to_yzzz_xxzz[k] * tbe_0 * tke_0;

            to_y_x_zzz_yyy[k] = 4.0 * to_yzzz_xyyy[k] * tbe_0 * tke_0;

            to_y_x_zzz_yyz[k] = 4.0 * to_yzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_x_zzz_yzz[k] = 4.0 * to_yzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_x_zzz_zzz[k] = 4.0 * to_yzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 400-410 components of targeted buffer : FF

        auto to_y_y_xxx_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 0);

        auto to_y_y_xxx_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 1);

        auto to_y_y_xxx_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 2);

        auto to_y_y_xxx_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 3);

        auto to_y_y_xxx_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 4);

        auto to_y_y_xxx_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 5);

        auto to_y_y_xxx_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 6);

        auto to_y_y_xxx_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 7);

        auto to_y_y_xxx_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 8);

        auto to_y_y_xxx_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_xxxy_xx, to_xxxy_xxxy, to_xxxy_xxyy, to_xxxy_xxyz, to_xxxy_xy, to_xxxy_xyyy, to_xxxy_xyyz, to_xxxy_xyzz, to_xxxy_xz, to_xxxy_yy, to_xxxy_yyyy, to_xxxy_yyyz, to_xxxy_yyzz, to_xxxy_yz, to_xxxy_yzzz, to_xxxy_zz, to_y_y_xxx_xxx, to_y_y_xxx_xxy, to_y_y_xxx_xxz, to_y_y_xxx_xyy, to_y_y_xxx_xyz, to_y_y_xxx_xzz, to_y_y_xxx_yyy, to_y_y_xxx_yyz, to_y_y_xxx_yzz, to_y_y_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxx_xxx[k] = 4.0 * to_xxxy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xxx_xxy[k] = -2.0 * to_xxxy_xx[k] * tbe_0 + 4.0 * to_xxxy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xxx_xxz[k] = 4.0 * to_xxxy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xxx_xyy[k] = -4.0 * to_xxxy_xy[k] * tbe_0 + 4.0 * to_xxxy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xxx_xyz[k] = -2.0 * to_xxxy_xz[k] * tbe_0 + 4.0 * to_xxxy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xxx_xzz[k] = 4.0 * to_xxxy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xxx_yyy[k] = -6.0 * to_xxxy_yy[k] * tbe_0 + 4.0 * to_xxxy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xxx_yyz[k] = -4.0 * to_xxxy_yz[k] * tbe_0 + 4.0 * to_xxxy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xxx_yzz[k] = -2.0 * to_xxxy_zz[k] * tbe_0 + 4.0 * to_xxxy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xxx_zzz[k] = 4.0 * to_xxxy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 410-420 components of targeted buffer : FF

        auto to_y_y_xxy_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 10);

        auto to_y_y_xxy_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 11);

        auto to_y_y_xxy_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 12);

        auto to_y_y_xxy_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 13);

        auto to_y_y_xxy_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 14);

        auto to_y_y_xxy_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 15);

        auto to_y_y_xxy_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 16);

        auto to_y_y_xxy_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 17);

        auto to_y_y_xxy_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 18);

        auto to_y_y_xxy_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_xx_xx, to_xx_xxxy, to_xx_xxyy, to_xx_xxyz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_yy, to_xx_yyyy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, to_xxyy_xx, to_xxyy_xxxy, to_xxyy_xxyy, to_xxyy_xxyz, to_xxyy_xy, to_xxyy_xyyy, to_xxyy_xyyz, to_xxyy_xyzz, to_xxyy_xz, to_xxyy_yy, to_xxyy_yyyy, to_xxyy_yyyz, to_xxyy_yyzz, to_xxyy_yz, to_xxyy_yzzz, to_xxyy_zz, to_y_y_xxy_xxx, to_y_y_xxy_xxy, to_y_y_xxy_xxz, to_y_y_xxy_xyy, to_y_y_xxy_xyz, to_y_y_xxy_xzz, to_y_y_xxy_yyy, to_y_y_xxy_yyz, to_y_y_xxy_yzz, to_y_y_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxy_xxx[k] = -2.0 * to_xx_xxxy[k] * tke_0 + 4.0 * to_xxyy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xxy_xxy[k] = to_xx_xx[k] - 2.0 * to_xx_xxyy[k] * tke_0 - 2.0 * to_xxyy_xx[k] * tbe_0 + 4.0 * to_xxyy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xxy_xxz[k] = -2.0 * to_xx_xxyz[k] * tke_0 + 4.0 * to_xxyy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xxy_xyy[k] = 2.0 * to_xx_xy[k] - 2.0 * to_xx_xyyy[k] * tke_0 - 4.0 * to_xxyy_xy[k] * tbe_0 + 4.0 * to_xxyy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xxy_xyz[k] = to_xx_xz[k] - 2.0 * to_xx_xyyz[k] * tke_0 - 2.0 * to_xxyy_xz[k] * tbe_0 + 4.0 * to_xxyy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xxy_xzz[k] = -2.0 * to_xx_xyzz[k] * tke_0 + 4.0 * to_xxyy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xxy_yyy[k] = 3.0 * to_xx_yy[k] - 2.0 * to_xx_yyyy[k] * tke_0 - 6.0 * to_xxyy_yy[k] * tbe_0 + 4.0 * to_xxyy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xxy_yyz[k] = 2.0 * to_xx_yz[k] - 2.0 * to_xx_yyyz[k] * tke_0 - 4.0 * to_xxyy_yz[k] * tbe_0 + 4.0 * to_xxyy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xxy_yzz[k] = to_xx_zz[k] - 2.0 * to_xx_yyzz[k] * tke_0 - 2.0 * to_xxyy_zz[k] * tbe_0 + 4.0 * to_xxyy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xxy_zzz[k] = -2.0 * to_xx_yzzz[k] * tke_0 + 4.0 * to_xxyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 420-430 components of targeted buffer : FF

        auto to_y_y_xxz_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 20);

        auto to_y_y_xxz_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 21);

        auto to_y_y_xxz_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 22);

        auto to_y_y_xxz_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 23);

        auto to_y_y_xxz_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 24);

        auto to_y_y_xxz_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 25);

        auto to_y_y_xxz_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 26);

        auto to_y_y_xxz_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 27);

        auto to_y_y_xxz_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 28);

        auto to_y_y_xxz_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_xxyz_xx, to_xxyz_xxxy, to_xxyz_xxyy, to_xxyz_xxyz, to_xxyz_xy, to_xxyz_xyyy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_yy, to_xxyz_yyyy, to_xxyz_yyyz, to_xxyz_yyzz, to_xxyz_yz, to_xxyz_yzzz, to_xxyz_zz, to_y_y_xxz_xxx, to_y_y_xxz_xxy, to_y_y_xxz_xxz, to_y_y_xxz_xyy, to_y_y_xxz_xyz, to_y_y_xxz_xzz, to_y_y_xxz_yyy, to_y_y_xxz_yyz, to_y_y_xxz_yzz, to_y_y_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxz_xxx[k] = 4.0 * to_xxyz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xxz_xxy[k] = -2.0 * to_xxyz_xx[k] * tbe_0 + 4.0 * to_xxyz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xxz_xxz[k] = 4.0 * to_xxyz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xxz_xyy[k] = -4.0 * to_xxyz_xy[k] * tbe_0 + 4.0 * to_xxyz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xxz_xyz[k] = -2.0 * to_xxyz_xz[k] * tbe_0 + 4.0 * to_xxyz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xxz_xzz[k] = 4.0 * to_xxyz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xxz_yyy[k] = -6.0 * to_xxyz_yy[k] * tbe_0 + 4.0 * to_xxyz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xxz_yyz[k] = -4.0 * to_xxyz_yz[k] * tbe_0 + 4.0 * to_xxyz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xxz_yzz[k] = -2.0 * to_xxyz_zz[k] * tbe_0 + 4.0 * to_xxyz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xxz_zzz[k] = 4.0 * to_xxyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 430-440 components of targeted buffer : FF

        auto to_y_y_xyy_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 30);

        auto to_y_y_xyy_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 31);

        auto to_y_y_xyy_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 32);

        auto to_y_y_xyy_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 33);

        auto to_y_y_xyy_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 34);

        auto to_y_y_xyy_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 35);

        auto to_y_y_xyy_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 36);

        auto to_y_y_xyy_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 37);

        auto to_y_y_xyy_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 38);

        auto to_y_y_xyy_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_xy_xx, to_xy_xxxy, to_xy_xxyy, to_xy_xxyz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_yy, to_xy_yyyy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, to_xyyy_xx, to_xyyy_xxxy, to_xyyy_xxyy, to_xyyy_xxyz, to_xyyy_xy, to_xyyy_xyyy, to_xyyy_xyyz, to_xyyy_xyzz, to_xyyy_xz, to_xyyy_yy, to_xyyy_yyyy, to_xyyy_yyyz, to_xyyy_yyzz, to_xyyy_yz, to_xyyy_yzzz, to_xyyy_zz, to_y_y_xyy_xxx, to_y_y_xyy_xxy, to_y_y_xyy_xxz, to_y_y_xyy_xyy, to_y_y_xyy_xyz, to_y_y_xyy_xzz, to_y_y_xyy_yyy, to_y_y_xyy_yyz, to_y_y_xyy_yzz, to_y_y_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyy_xxx[k] = -4.0 * to_xy_xxxy[k] * tke_0 + 4.0 * to_xyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xyy_xxy[k] = 2.0 * to_xy_xx[k] - 4.0 * to_xy_xxyy[k] * tke_0 - 2.0 * to_xyyy_xx[k] * tbe_0 + 4.0 * to_xyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xyy_xxz[k] = -4.0 * to_xy_xxyz[k] * tke_0 + 4.0 * to_xyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xyy_xyy[k] = 4.0 * to_xy_xy[k] - 4.0 * to_xy_xyyy[k] * tke_0 - 4.0 * to_xyyy_xy[k] * tbe_0 + 4.0 * to_xyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xyy_xyz[k] = 2.0 * to_xy_xz[k] - 4.0 * to_xy_xyyz[k] * tke_0 - 2.0 * to_xyyy_xz[k] * tbe_0 + 4.0 * to_xyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xyy_xzz[k] = -4.0 * to_xy_xyzz[k] * tke_0 + 4.0 * to_xyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xyy_yyy[k] = 6.0 * to_xy_yy[k] - 4.0 * to_xy_yyyy[k] * tke_0 - 6.0 * to_xyyy_yy[k] * tbe_0 + 4.0 * to_xyyy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xyy_yyz[k] = 4.0 * to_xy_yz[k] - 4.0 * to_xy_yyyz[k] * tke_0 - 4.0 * to_xyyy_yz[k] * tbe_0 + 4.0 * to_xyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xyy_yzz[k] = 2.0 * to_xy_zz[k] - 4.0 * to_xy_yyzz[k] * tke_0 - 2.0 * to_xyyy_zz[k] * tbe_0 + 4.0 * to_xyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xyy_zzz[k] = -4.0 * to_xy_yzzz[k] * tke_0 + 4.0 * to_xyyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 440-450 components of targeted buffer : FF

        auto to_y_y_xyz_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 40);

        auto to_y_y_xyz_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 41);

        auto to_y_y_xyz_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 42);

        auto to_y_y_xyz_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 43);

        auto to_y_y_xyz_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 44);

        auto to_y_y_xyz_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 45);

        auto to_y_y_xyz_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 46);

        auto to_y_y_xyz_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 47);

        auto to_y_y_xyz_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 48);

        auto to_y_y_xyz_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_xyyz_xx, to_xyyz_xxxy, to_xyyz_xxyy, to_xyyz_xxyz, to_xyyz_xy, to_xyyz_xyyy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_yy, to_xyyz_yyyy, to_xyyz_yyyz, to_xyyz_yyzz, to_xyyz_yz, to_xyyz_yzzz, to_xyyz_zz, to_xz_xx, to_xz_xxxy, to_xz_xxyy, to_xz_xxyz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_yy, to_xz_yyyy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, to_y_y_xyz_xxx, to_y_y_xyz_xxy, to_y_y_xyz_xxz, to_y_y_xyz_xyy, to_y_y_xyz_xyz, to_y_y_xyz_xzz, to_y_y_xyz_yyy, to_y_y_xyz_yyz, to_y_y_xyz_yzz, to_y_y_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyz_xxx[k] = -2.0 * to_xz_xxxy[k] * tke_0 + 4.0 * to_xyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xyz_xxy[k] = to_xz_xx[k] - 2.0 * to_xz_xxyy[k] * tke_0 - 2.0 * to_xyyz_xx[k] * tbe_0 + 4.0 * to_xyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xyz_xxz[k] = -2.0 * to_xz_xxyz[k] * tke_0 + 4.0 * to_xyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xyz_xyy[k] = 2.0 * to_xz_xy[k] - 2.0 * to_xz_xyyy[k] * tke_0 - 4.0 * to_xyyz_xy[k] * tbe_0 + 4.0 * to_xyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xyz_xyz[k] = to_xz_xz[k] - 2.0 * to_xz_xyyz[k] * tke_0 - 2.0 * to_xyyz_xz[k] * tbe_0 + 4.0 * to_xyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xyz_xzz[k] = -2.0 * to_xz_xyzz[k] * tke_0 + 4.0 * to_xyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xyz_yyy[k] = 3.0 * to_xz_yy[k] - 2.0 * to_xz_yyyy[k] * tke_0 - 6.0 * to_xyyz_yy[k] * tbe_0 + 4.0 * to_xyyz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xyz_yyz[k] = 2.0 * to_xz_yz[k] - 2.0 * to_xz_yyyz[k] * tke_0 - 4.0 * to_xyyz_yz[k] * tbe_0 + 4.0 * to_xyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xyz_yzz[k] = to_xz_zz[k] - 2.0 * to_xz_yyzz[k] * tke_0 - 2.0 * to_xyyz_zz[k] * tbe_0 + 4.0 * to_xyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xyz_zzz[k] = -2.0 * to_xz_yzzz[k] * tke_0 + 4.0 * to_xyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 450-460 components of targeted buffer : FF

        auto to_y_y_xzz_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 50);

        auto to_y_y_xzz_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 51);

        auto to_y_y_xzz_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 52);

        auto to_y_y_xzz_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 53);

        auto to_y_y_xzz_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 54);

        auto to_y_y_xzz_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 55);

        auto to_y_y_xzz_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 56);

        auto to_y_y_xzz_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 57);

        auto to_y_y_xzz_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 58);

        auto to_y_y_xzz_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_xyzz_xx, to_xyzz_xxxy, to_xyzz_xxyy, to_xyzz_xxyz, to_xyzz_xy, to_xyzz_xyyy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_yy, to_xyzz_yyyy, to_xyzz_yyyz, to_xyzz_yyzz, to_xyzz_yz, to_xyzz_yzzz, to_xyzz_zz, to_y_y_xzz_xxx, to_y_y_xzz_xxy, to_y_y_xzz_xxz, to_y_y_xzz_xyy, to_y_y_xzz_xyz, to_y_y_xzz_xzz, to_y_y_xzz_yyy, to_y_y_xzz_yyz, to_y_y_xzz_yzz, to_y_y_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xzz_xxx[k] = 4.0 * to_xyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_xzz_xxy[k] = -2.0 * to_xyzz_xx[k] * tbe_0 + 4.0 * to_xyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_xzz_xxz[k] = 4.0 * to_xyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_xzz_xyy[k] = -4.0 * to_xyzz_xy[k] * tbe_0 + 4.0 * to_xyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_xzz_xyz[k] = -2.0 * to_xyzz_xz[k] * tbe_0 + 4.0 * to_xyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_xzz_xzz[k] = 4.0 * to_xyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_xzz_yyy[k] = -6.0 * to_xyzz_yy[k] * tbe_0 + 4.0 * to_xyzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_xzz_yyz[k] = -4.0 * to_xyzz_yz[k] * tbe_0 + 4.0 * to_xyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_xzz_yzz[k] = -2.0 * to_xyzz_zz[k] * tbe_0 + 4.0 * to_xyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_xzz_zzz[k] = 4.0 * to_xyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 460-470 components of targeted buffer : FF

        auto to_y_y_yyy_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 60);

        auto to_y_y_yyy_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 61);

        auto to_y_y_yyy_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 62);

        auto to_y_y_yyy_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 63);

        auto to_y_y_yyy_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 64);

        auto to_y_y_yyy_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 65);

        auto to_y_y_yyy_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 66);

        auto to_y_y_yyy_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 67);

        auto to_y_y_yyy_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 68);

        auto to_y_y_yyy_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_y_y_yyy_xxx, to_y_y_yyy_xxy, to_y_y_yyy_xxz, to_y_y_yyy_xyy, to_y_y_yyy_xyz, to_y_y_yyy_xzz, to_y_y_yyy_yyy, to_y_y_yyy_yyz, to_y_y_yyy_yzz, to_y_y_yyy_zzz, to_yy_xx, to_yy_xxxy, to_yy_xxyy, to_yy_xxyz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_yy, to_yy_yyyy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, to_yyyy_xx, to_yyyy_xxxy, to_yyyy_xxyy, to_yyyy_xxyz, to_yyyy_xy, to_yyyy_xyyy, to_yyyy_xyyz, to_yyyy_xyzz, to_yyyy_xz, to_yyyy_yy, to_yyyy_yyyy, to_yyyy_yyyz, to_yyyy_yyzz, to_yyyy_yz, to_yyyy_yzzz, to_yyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyy_xxx[k] = -6.0 * to_yy_xxxy[k] * tke_0 + 4.0 * to_yyyy_xxxy[k] * tbe_0 * tke_0;

            to_y_y_yyy_xxy[k] = 3.0 * to_yy_xx[k] - 6.0 * to_yy_xxyy[k] * tke_0 - 2.0 * to_yyyy_xx[k] * tbe_0 + 4.0 * to_yyyy_xxyy[k] * tbe_0 * tke_0;

            to_y_y_yyy_xxz[k] = -6.0 * to_yy_xxyz[k] * tke_0 + 4.0 * to_yyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_y_yyy_xyy[k] = 6.0 * to_yy_xy[k] - 6.0 * to_yy_xyyy[k] * tke_0 - 4.0 * to_yyyy_xy[k] * tbe_0 + 4.0 * to_yyyy_xyyy[k] * tbe_0 * tke_0;

            to_y_y_yyy_xyz[k] = 3.0 * to_yy_xz[k] - 6.0 * to_yy_xyyz[k] * tke_0 - 2.0 * to_yyyy_xz[k] * tbe_0 + 4.0 * to_yyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_y_yyy_xzz[k] = -6.0 * to_yy_xyzz[k] * tke_0 + 4.0 * to_yyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_y_yyy_yyy[k] = 9.0 * to_yy_yy[k] - 6.0 * to_yy_yyyy[k] * tke_0 - 6.0 * to_yyyy_yy[k] * tbe_0 + 4.0 * to_yyyy_yyyy[k] * tbe_0 * tke_0;

            to_y_y_yyy_yyz[k] = 6.0 * to_yy_yz[k] - 6.0 * to_yy_yyyz[k] * tke_0 - 4.0 * to_yyyy_yz[k] * tbe_0 + 4.0 * to_yyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_y_yyy_yzz[k] = 3.0 * to_yy_zz[k] - 6.0 * to_yy_yyzz[k] * tke_0 - 2.0 * to_yyyy_zz[k] * tbe_0 + 4.0 * to_yyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_y_yyy_zzz[k] = -6.0 * to_yy_yzzz[k] * tke_0 + 4.0 * to_yyyy_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 470-480 components of targeted buffer : FF

        auto to_y_y_yyz_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 70);

        auto to_y_y_yyz_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 71);

        auto to_y_y_yyz_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 72);

        auto to_y_y_yyz_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 73);

        auto to_y_y_yyz_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 74);

        auto to_y_y_yyz_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 75);

        auto to_y_y_yyz_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 76);

        auto to_y_y_yyz_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 77);

        auto to_y_y_yyz_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 78);

        auto to_y_y_yyz_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_y_y_yyz_xxx, to_y_y_yyz_xxy, to_y_y_yyz_xxz, to_y_y_yyz_xyy, to_y_y_yyz_xyz, to_y_y_yyz_xzz, to_y_y_yyz_yyy, to_y_y_yyz_yyz, to_y_y_yyz_yzz, to_y_y_yyz_zzz, to_yyyz_xx, to_yyyz_xxxy, to_yyyz_xxyy, to_yyyz_xxyz, to_yyyz_xy, to_yyyz_xyyy, to_yyyz_xyyz, to_yyyz_xyzz, to_yyyz_xz, to_yyyz_yy, to_yyyz_yyyy, to_yyyz_yyyz, to_yyyz_yyzz, to_yyyz_yz, to_yyyz_yzzz, to_yyyz_zz, to_yz_xx, to_yz_xxxy, to_yz_xxyy, to_yz_xxyz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_yy, to_yz_yyyy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyz_xxx[k] = -4.0 * to_yz_xxxy[k] * tke_0 + 4.0 * to_yyyz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_yyz_xxy[k] = 2.0 * to_yz_xx[k] - 4.0 * to_yz_xxyy[k] * tke_0 - 2.0 * to_yyyz_xx[k] * tbe_0 + 4.0 * to_yyyz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_yyz_xxz[k] = -4.0 * to_yz_xxyz[k] * tke_0 + 4.0 * to_yyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_yyz_xyy[k] = 4.0 * to_yz_xy[k] - 4.0 * to_yz_xyyy[k] * tke_0 - 4.0 * to_yyyz_xy[k] * tbe_0 + 4.0 * to_yyyz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_yyz_xyz[k] = 2.0 * to_yz_xz[k] - 4.0 * to_yz_xyyz[k] * tke_0 - 2.0 * to_yyyz_xz[k] * tbe_0 + 4.0 * to_yyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_yyz_xzz[k] = -4.0 * to_yz_xyzz[k] * tke_0 + 4.0 * to_yyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_yyz_yyy[k] = 6.0 * to_yz_yy[k] - 4.0 * to_yz_yyyy[k] * tke_0 - 6.0 * to_yyyz_yy[k] * tbe_0 + 4.0 * to_yyyz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_yyz_yyz[k] = 4.0 * to_yz_yz[k] - 4.0 * to_yz_yyyz[k] * tke_0 - 4.0 * to_yyyz_yz[k] * tbe_0 + 4.0 * to_yyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_yyz_yzz[k] = 2.0 * to_yz_zz[k] - 4.0 * to_yz_yyzz[k] * tke_0 - 2.0 * to_yyyz_zz[k] * tbe_0 + 4.0 * to_yyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_yyz_zzz[k] = -4.0 * to_yz_yzzz[k] * tke_0 + 4.0 * to_yyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 480-490 components of targeted buffer : FF

        auto to_y_y_yzz_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 80);

        auto to_y_y_yzz_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 81);

        auto to_y_y_yzz_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 82);

        auto to_y_y_yzz_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 83);

        auto to_y_y_yzz_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 84);

        auto to_y_y_yzz_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 85);

        auto to_y_y_yzz_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 86);

        auto to_y_y_yzz_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 87);

        auto to_y_y_yzz_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 88);

        auto to_y_y_yzz_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_y_y_yzz_xxx, to_y_y_yzz_xxy, to_y_y_yzz_xxz, to_y_y_yzz_xyy, to_y_y_yzz_xyz, to_y_y_yzz_xzz, to_y_y_yzz_yyy, to_y_y_yzz_yyz, to_y_y_yzz_yzz, to_y_y_yzz_zzz, to_yyzz_xx, to_yyzz_xxxy, to_yyzz_xxyy, to_yyzz_xxyz, to_yyzz_xy, to_yyzz_xyyy, to_yyzz_xyyz, to_yyzz_xyzz, to_yyzz_xz, to_yyzz_yy, to_yyzz_yyyy, to_yyzz_yyyz, to_yyzz_yyzz, to_yyzz_yz, to_yyzz_yzzz, to_yyzz_zz, to_zz_xx, to_zz_xxxy, to_zz_xxyy, to_zz_xxyz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_yy, to_zz_yyyy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yzz_xxx[k] = -2.0 * to_zz_xxxy[k] * tke_0 + 4.0 * to_yyzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_yzz_xxy[k] = to_zz_xx[k] - 2.0 * to_zz_xxyy[k] * tke_0 - 2.0 * to_yyzz_xx[k] * tbe_0 + 4.0 * to_yyzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_yzz_xxz[k] = -2.0 * to_zz_xxyz[k] * tke_0 + 4.0 * to_yyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_yzz_xyy[k] = 2.0 * to_zz_xy[k] - 2.0 * to_zz_xyyy[k] * tke_0 - 4.0 * to_yyzz_xy[k] * tbe_0 + 4.0 * to_yyzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_yzz_xyz[k] = to_zz_xz[k] - 2.0 * to_zz_xyyz[k] * tke_0 - 2.0 * to_yyzz_xz[k] * tbe_0 + 4.0 * to_yyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_yzz_xzz[k] = -2.0 * to_zz_xyzz[k] * tke_0 + 4.0 * to_yyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_yzz_yyy[k] = 3.0 * to_zz_yy[k] - 2.0 * to_zz_yyyy[k] * tke_0 - 6.0 * to_yyzz_yy[k] * tbe_0 + 4.0 * to_yyzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_yzz_yyz[k] = 2.0 * to_zz_yz[k] - 2.0 * to_zz_yyyz[k] * tke_0 - 4.0 * to_yyzz_yz[k] * tbe_0 + 4.0 * to_yyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_yzz_yzz[k] = to_zz_zz[k] - 2.0 * to_zz_yyzz[k] * tke_0 - 2.0 * to_yyzz_zz[k] * tbe_0 + 4.0 * to_yyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_yzz_zzz[k] = -2.0 * to_zz_yzzz[k] * tke_0 + 4.0 * to_yyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 490-500 components of targeted buffer : FF

        auto to_y_y_zzz_xxx = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 90);

        auto to_y_y_zzz_xxy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 91);

        auto to_y_y_zzz_xxz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 92);

        auto to_y_y_zzz_xyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 93);

        auto to_y_y_zzz_xyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 94);

        auto to_y_y_zzz_xzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 95);

        auto to_y_y_zzz_yyy = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 96);

        auto to_y_y_zzz_yyz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 97);

        auto to_y_y_zzz_yzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 98);

        auto to_y_y_zzz_zzz = pbuffer.data(idx_op_geom_101_ff + 4 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_y_y_zzz_xxx, to_y_y_zzz_xxy, to_y_y_zzz_xxz, to_y_y_zzz_xyy, to_y_y_zzz_xyz, to_y_y_zzz_xzz, to_y_y_zzz_yyy, to_y_y_zzz_yyz, to_y_y_zzz_yzz, to_y_y_zzz_zzz, to_yzzz_xx, to_yzzz_xxxy, to_yzzz_xxyy, to_yzzz_xxyz, to_yzzz_xy, to_yzzz_xyyy, to_yzzz_xyyz, to_yzzz_xyzz, to_yzzz_xz, to_yzzz_yy, to_yzzz_yyyy, to_yzzz_yyyz, to_yzzz_yyzz, to_yzzz_yz, to_yzzz_yzzz, to_yzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zzz_xxx[k] = 4.0 * to_yzzz_xxxy[k] * tbe_0 * tke_0;

            to_y_y_zzz_xxy[k] = -2.0 * to_yzzz_xx[k] * tbe_0 + 4.0 * to_yzzz_xxyy[k] * tbe_0 * tke_0;

            to_y_y_zzz_xxz[k] = 4.0 * to_yzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_y_zzz_xyy[k] = -4.0 * to_yzzz_xy[k] * tbe_0 + 4.0 * to_yzzz_xyyy[k] * tbe_0 * tke_0;

            to_y_y_zzz_xyz[k] = -2.0 * to_yzzz_xz[k] * tbe_0 + 4.0 * to_yzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_y_zzz_xzz[k] = 4.0 * to_yzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_y_zzz_yyy[k] = -6.0 * to_yzzz_yy[k] * tbe_0 + 4.0 * to_yzzz_yyyy[k] * tbe_0 * tke_0;

            to_y_y_zzz_yyz[k] = -4.0 * to_yzzz_yz[k] * tbe_0 + 4.0 * to_yzzz_yyyz[k] * tbe_0 * tke_0;

            to_y_y_zzz_yzz[k] = -2.0 * to_yzzz_zz[k] * tbe_0 + 4.0 * to_yzzz_yyzz[k] * tbe_0 * tke_0;

            to_y_y_zzz_zzz[k] = 4.0 * to_yzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 500-510 components of targeted buffer : FF

        auto to_y_z_xxx_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 0);

        auto to_y_z_xxx_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 1);

        auto to_y_z_xxx_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 2);

        auto to_y_z_xxx_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 3);

        auto to_y_z_xxx_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 4);

        auto to_y_z_xxx_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 5);

        auto to_y_z_xxx_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 6);

        auto to_y_z_xxx_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 7);

        auto to_y_z_xxx_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 8);

        auto to_y_z_xxx_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_xxxy_xx, to_xxxy_xxxz, to_xxxy_xxyz, to_xxxy_xxzz, to_xxxy_xy, to_xxxy_xyyz, to_xxxy_xyzz, to_xxxy_xz, to_xxxy_xzzz, to_xxxy_yy, to_xxxy_yyyz, to_xxxy_yyzz, to_xxxy_yz, to_xxxy_yzzz, to_xxxy_zz, to_xxxy_zzzz, to_y_z_xxx_xxx, to_y_z_xxx_xxy, to_y_z_xxx_xxz, to_y_z_xxx_xyy, to_y_z_xxx_xyz, to_y_z_xxx_xzz, to_y_z_xxx_yyy, to_y_z_xxx_yyz, to_y_z_xxx_yzz, to_y_z_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxx_xxx[k] = 4.0 * to_xxxy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xxy[k] = 4.0 * to_xxxy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xxz[k] = -2.0 * to_xxxy_xx[k] * tbe_0 + 4.0 * to_xxxy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xyy[k] = 4.0 * to_xxxy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xyz[k] = -2.0 * to_xxxy_xy[k] * tbe_0 + 4.0 * to_xxxy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_xzz[k] = -4.0 * to_xxxy_xz[k] * tbe_0 + 4.0 * to_xxxy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_yyy[k] = 4.0 * to_xxxy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xxx_yyz[k] = -2.0 * to_xxxy_yy[k] * tbe_0 + 4.0 * to_xxxy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_yzz[k] = -4.0 * to_xxxy_yz[k] * tbe_0 + 4.0 * to_xxxy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xxx_zzz[k] = -6.0 * to_xxxy_zz[k] * tbe_0 + 4.0 * to_xxxy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 510-520 components of targeted buffer : FF

        auto to_y_z_xxy_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 10);

        auto to_y_z_xxy_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 11);

        auto to_y_z_xxy_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 12);

        auto to_y_z_xxy_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 13);

        auto to_y_z_xxy_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 14);

        auto to_y_z_xxy_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 15);

        auto to_y_z_xxy_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 16);

        auto to_y_z_xxy_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 17);

        auto to_y_z_xxy_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 18);

        auto to_y_z_xxy_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_xx_xx, to_xx_xxxz, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, to_xx_zzzz, to_xxyy_xx, to_xxyy_xxxz, to_xxyy_xxyz, to_xxyy_xxzz, to_xxyy_xy, to_xxyy_xyyz, to_xxyy_xyzz, to_xxyy_xz, to_xxyy_xzzz, to_xxyy_yy, to_xxyy_yyyz, to_xxyy_yyzz, to_xxyy_yz, to_xxyy_yzzz, to_xxyy_zz, to_xxyy_zzzz, to_y_z_xxy_xxx, to_y_z_xxy_xxy, to_y_z_xxy_xxz, to_y_z_xxy_xyy, to_y_z_xxy_xyz, to_y_z_xxy_xzz, to_y_z_xxy_yyy, to_y_z_xxy_yyz, to_y_z_xxy_yzz, to_y_z_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxy_xxx[k] = -2.0 * to_xx_xxxz[k] * tke_0 + 4.0 * to_xxyy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xxy[k] = -2.0 * to_xx_xxyz[k] * tke_0 + 4.0 * to_xxyy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xxz[k] = to_xx_xx[k] - 2.0 * to_xx_xxzz[k] * tke_0 - 2.0 * to_xxyy_xx[k] * tbe_0 + 4.0 * to_xxyy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xyy[k] = -2.0 * to_xx_xyyz[k] * tke_0 + 4.0 * to_xxyy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xyz[k] = to_xx_xy[k] - 2.0 * to_xx_xyzz[k] * tke_0 - 2.0 * to_xxyy_xy[k] * tbe_0 + 4.0 * to_xxyy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_xzz[k] = 2.0 * to_xx_xz[k] - 2.0 * to_xx_xzzz[k] * tke_0 - 4.0 * to_xxyy_xz[k] * tbe_0 + 4.0 * to_xxyy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_yyy[k] = -2.0 * to_xx_yyyz[k] * tke_0 + 4.0 * to_xxyy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xxy_yyz[k] = to_xx_yy[k] - 2.0 * to_xx_yyzz[k] * tke_0 - 2.0 * to_xxyy_yy[k] * tbe_0 + 4.0 * to_xxyy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_yzz[k] = 2.0 * to_xx_yz[k] - 2.0 * to_xx_yzzz[k] * tke_0 - 4.0 * to_xxyy_yz[k] * tbe_0 + 4.0 * to_xxyy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xxy_zzz[k] = 3.0 * to_xx_zz[k] - 2.0 * to_xx_zzzz[k] * tke_0 - 6.0 * to_xxyy_zz[k] * tbe_0 + 4.0 * to_xxyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 520-530 components of targeted buffer : FF

        auto to_y_z_xxz_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 20);

        auto to_y_z_xxz_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 21);

        auto to_y_z_xxz_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 22);

        auto to_y_z_xxz_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 23);

        auto to_y_z_xxz_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 24);

        auto to_y_z_xxz_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 25);

        auto to_y_z_xxz_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 26);

        auto to_y_z_xxz_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 27);

        auto to_y_z_xxz_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 28);

        auto to_y_z_xxz_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_xxyz_xx, to_xxyz_xxxz, to_xxyz_xxyz, to_xxyz_xxzz, to_xxyz_xy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_xzzz, to_xxyz_yy, to_xxyz_yyyz, to_xxyz_yyzz, to_xxyz_yz, to_xxyz_yzzz, to_xxyz_zz, to_xxyz_zzzz, to_y_z_xxz_xxx, to_y_z_xxz_xxy, to_y_z_xxz_xxz, to_y_z_xxz_xyy, to_y_z_xxz_xyz, to_y_z_xxz_xzz, to_y_z_xxz_yyy, to_y_z_xxz_yyz, to_y_z_xxz_yzz, to_y_z_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxz_xxx[k] = 4.0 * to_xxyz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xxy[k] = 4.0 * to_xxyz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xxz[k] = -2.0 * to_xxyz_xx[k] * tbe_0 + 4.0 * to_xxyz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xyy[k] = 4.0 * to_xxyz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xyz[k] = -2.0 * to_xxyz_xy[k] * tbe_0 + 4.0 * to_xxyz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_xzz[k] = -4.0 * to_xxyz_xz[k] * tbe_0 + 4.0 * to_xxyz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_yyy[k] = 4.0 * to_xxyz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xxz_yyz[k] = -2.0 * to_xxyz_yy[k] * tbe_0 + 4.0 * to_xxyz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_yzz[k] = -4.0 * to_xxyz_yz[k] * tbe_0 + 4.0 * to_xxyz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xxz_zzz[k] = -6.0 * to_xxyz_zz[k] * tbe_0 + 4.0 * to_xxyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 530-540 components of targeted buffer : FF

        auto to_y_z_xyy_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 30);

        auto to_y_z_xyy_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 31);

        auto to_y_z_xyy_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 32);

        auto to_y_z_xyy_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 33);

        auto to_y_z_xyy_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 34);

        auto to_y_z_xyy_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 35);

        auto to_y_z_xyy_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 36);

        auto to_y_z_xyy_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 37);

        auto to_y_z_xyy_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 38);

        auto to_y_z_xyy_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_xy_xx, to_xy_xxxz, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, to_xy_zzzz, to_xyyy_xx, to_xyyy_xxxz, to_xyyy_xxyz, to_xyyy_xxzz, to_xyyy_xy, to_xyyy_xyyz, to_xyyy_xyzz, to_xyyy_xz, to_xyyy_xzzz, to_xyyy_yy, to_xyyy_yyyz, to_xyyy_yyzz, to_xyyy_yz, to_xyyy_yzzz, to_xyyy_zz, to_xyyy_zzzz, to_y_z_xyy_xxx, to_y_z_xyy_xxy, to_y_z_xyy_xxz, to_y_z_xyy_xyy, to_y_z_xyy_xyz, to_y_z_xyy_xzz, to_y_z_xyy_yyy, to_y_z_xyy_yyz, to_y_z_xyy_yzz, to_y_z_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyy_xxx[k] = -4.0 * to_xy_xxxz[k] * tke_0 + 4.0 * to_xyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xxy[k] = -4.0 * to_xy_xxyz[k] * tke_0 + 4.0 * to_xyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xxz[k] = 2.0 * to_xy_xx[k] - 4.0 * to_xy_xxzz[k] * tke_0 - 2.0 * to_xyyy_xx[k] * tbe_0 + 4.0 * to_xyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xyy[k] = -4.0 * to_xy_xyyz[k] * tke_0 + 4.0 * to_xyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xyz[k] = 2.0 * to_xy_xy[k] - 4.0 * to_xy_xyzz[k] * tke_0 - 2.0 * to_xyyy_xy[k] * tbe_0 + 4.0 * to_xyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_xzz[k] = 4.0 * to_xy_xz[k] - 4.0 * to_xy_xzzz[k] * tke_0 - 4.0 * to_xyyy_xz[k] * tbe_0 + 4.0 * to_xyyy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_yyy[k] = -4.0 * to_xy_yyyz[k] * tke_0 + 4.0 * to_xyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xyy_yyz[k] = 2.0 * to_xy_yy[k] - 4.0 * to_xy_yyzz[k] * tke_0 - 2.0 * to_xyyy_yy[k] * tbe_0 + 4.0 * to_xyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_yzz[k] = 4.0 * to_xy_yz[k] - 4.0 * to_xy_yzzz[k] * tke_0 - 4.0 * to_xyyy_yz[k] * tbe_0 + 4.0 * to_xyyy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xyy_zzz[k] = 6.0 * to_xy_zz[k] - 4.0 * to_xy_zzzz[k] * tke_0 - 6.0 * to_xyyy_zz[k] * tbe_0 + 4.0 * to_xyyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 540-550 components of targeted buffer : FF

        auto to_y_z_xyz_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 40);

        auto to_y_z_xyz_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 41);

        auto to_y_z_xyz_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 42);

        auto to_y_z_xyz_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 43);

        auto to_y_z_xyz_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 44);

        auto to_y_z_xyz_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 45);

        auto to_y_z_xyz_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 46);

        auto to_y_z_xyz_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 47);

        auto to_y_z_xyz_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 48);

        auto to_y_z_xyz_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_xyyz_xx, to_xyyz_xxxz, to_xyyz_xxyz, to_xyyz_xxzz, to_xyyz_xy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_xzzz, to_xyyz_yy, to_xyyz_yyyz, to_xyyz_yyzz, to_xyyz_yz, to_xyyz_yzzz, to_xyyz_zz, to_xyyz_zzzz, to_xz_xx, to_xz_xxxz, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, to_xz_zzzz, to_y_z_xyz_xxx, to_y_z_xyz_xxy, to_y_z_xyz_xxz, to_y_z_xyz_xyy, to_y_z_xyz_xyz, to_y_z_xyz_xzz, to_y_z_xyz_yyy, to_y_z_xyz_yyz, to_y_z_xyz_yzz, to_y_z_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyz_xxx[k] = -2.0 * to_xz_xxxz[k] * tke_0 + 4.0 * to_xyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xxy[k] = -2.0 * to_xz_xxyz[k] * tke_0 + 4.0 * to_xyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xxz[k] = to_xz_xx[k] - 2.0 * to_xz_xxzz[k] * tke_0 - 2.0 * to_xyyz_xx[k] * tbe_0 + 4.0 * to_xyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xyy[k] = -2.0 * to_xz_xyyz[k] * tke_0 + 4.0 * to_xyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xyz[k] = to_xz_xy[k] - 2.0 * to_xz_xyzz[k] * tke_0 - 2.0 * to_xyyz_xy[k] * tbe_0 + 4.0 * to_xyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_xzz[k] = 2.0 * to_xz_xz[k] - 2.0 * to_xz_xzzz[k] * tke_0 - 4.0 * to_xyyz_xz[k] * tbe_0 + 4.0 * to_xyyz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_yyy[k] = -2.0 * to_xz_yyyz[k] * tke_0 + 4.0 * to_xyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xyz_yyz[k] = to_xz_yy[k] - 2.0 * to_xz_yyzz[k] * tke_0 - 2.0 * to_xyyz_yy[k] * tbe_0 + 4.0 * to_xyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_yzz[k] = 2.0 * to_xz_yz[k] - 2.0 * to_xz_yzzz[k] * tke_0 - 4.0 * to_xyyz_yz[k] * tbe_0 + 4.0 * to_xyyz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xyz_zzz[k] = 3.0 * to_xz_zz[k] - 2.0 * to_xz_zzzz[k] * tke_0 - 6.0 * to_xyyz_zz[k] * tbe_0 + 4.0 * to_xyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 550-560 components of targeted buffer : FF

        auto to_y_z_xzz_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 50);

        auto to_y_z_xzz_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 51);

        auto to_y_z_xzz_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 52);

        auto to_y_z_xzz_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 53);

        auto to_y_z_xzz_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 54);

        auto to_y_z_xzz_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 55);

        auto to_y_z_xzz_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 56);

        auto to_y_z_xzz_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 57);

        auto to_y_z_xzz_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 58);

        auto to_y_z_xzz_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_xyzz_xx, to_xyzz_xxxz, to_xyzz_xxyz, to_xyzz_xxzz, to_xyzz_xy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_xzzz, to_xyzz_yy, to_xyzz_yyyz, to_xyzz_yyzz, to_xyzz_yz, to_xyzz_yzzz, to_xyzz_zz, to_xyzz_zzzz, to_y_z_xzz_xxx, to_y_z_xzz_xxy, to_y_z_xzz_xxz, to_y_z_xzz_xyy, to_y_z_xzz_xyz, to_y_z_xzz_xzz, to_y_z_xzz_yyy, to_y_z_xzz_yyz, to_y_z_xzz_yzz, to_y_z_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xzz_xxx[k] = 4.0 * to_xyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xxy[k] = 4.0 * to_xyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xxz[k] = -2.0 * to_xyzz_xx[k] * tbe_0 + 4.0 * to_xyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xyy[k] = 4.0 * to_xyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xyz[k] = -2.0 * to_xyzz_xy[k] * tbe_0 + 4.0 * to_xyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_xzz[k] = -4.0 * to_xyzz_xz[k] * tbe_0 + 4.0 * to_xyzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_yyy[k] = 4.0 * to_xyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_xzz_yyz[k] = -2.0 * to_xyzz_yy[k] * tbe_0 + 4.0 * to_xyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_yzz[k] = -4.0 * to_xyzz_yz[k] * tbe_0 + 4.0 * to_xyzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_xzz_zzz[k] = -6.0 * to_xyzz_zz[k] * tbe_0 + 4.0 * to_xyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 560-570 components of targeted buffer : FF

        auto to_y_z_yyy_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 60);

        auto to_y_z_yyy_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 61);

        auto to_y_z_yyy_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 62);

        auto to_y_z_yyy_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 63);

        auto to_y_z_yyy_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 64);

        auto to_y_z_yyy_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 65);

        auto to_y_z_yyy_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 66);

        auto to_y_z_yyy_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 67);

        auto to_y_z_yyy_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 68);

        auto to_y_z_yyy_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_y_z_yyy_xxx, to_y_z_yyy_xxy, to_y_z_yyy_xxz, to_y_z_yyy_xyy, to_y_z_yyy_xyz, to_y_z_yyy_xzz, to_y_z_yyy_yyy, to_y_z_yyy_yyz, to_y_z_yyy_yzz, to_y_z_yyy_zzz, to_yy_xx, to_yy_xxxz, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, to_yy_zzzz, to_yyyy_xx, to_yyyy_xxxz, to_yyyy_xxyz, to_yyyy_xxzz, to_yyyy_xy, to_yyyy_xyyz, to_yyyy_xyzz, to_yyyy_xz, to_yyyy_xzzz, to_yyyy_yy, to_yyyy_yyyz, to_yyyy_yyzz, to_yyyy_yz, to_yyyy_yzzz, to_yyyy_zz, to_yyyy_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyy_xxx[k] = -6.0 * to_yy_xxxz[k] * tke_0 + 4.0 * to_yyyy_xxxz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xxy[k] = -6.0 * to_yy_xxyz[k] * tke_0 + 4.0 * to_yyyy_xxyz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xxz[k] = 3.0 * to_yy_xx[k] - 6.0 * to_yy_xxzz[k] * tke_0 - 2.0 * to_yyyy_xx[k] * tbe_0 + 4.0 * to_yyyy_xxzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xyy[k] = -6.0 * to_yy_xyyz[k] * tke_0 + 4.0 * to_yyyy_xyyz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xyz[k] = 3.0 * to_yy_xy[k] - 6.0 * to_yy_xyzz[k] * tke_0 - 2.0 * to_yyyy_xy[k] * tbe_0 + 4.0 * to_yyyy_xyzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_xzz[k] = 6.0 * to_yy_xz[k] - 6.0 * to_yy_xzzz[k] * tke_0 - 4.0 * to_yyyy_xz[k] * tbe_0 + 4.0 * to_yyyy_xzzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_yyy[k] = -6.0 * to_yy_yyyz[k] * tke_0 + 4.0 * to_yyyy_yyyz[k] * tbe_0 * tke_0;

            to_y_z_yyy_yyz[k] = 3.0 * to_yy_yy[k] - 6.0 * to_yy_yyzz[k] * tke_0 - 2.0 * to_yyyy_yy[k] * tbe_0 + 4.0 * to_yyyy_yyzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_yzz[k] = 6.0 * to_yy_yz[k] - 6.0 * to_yy_yzzz[k] * tke_0 - 4.0 * to_yyyy_yz[k] * tbe_0 + 4.0 * to_yyyy_yzzz[k] * tbe_0 * tke_0;

            to_y_z_yyy_zzz[k] = 9.0 * to_yy_zz[k] - 6.0 * to_yy_zzzz[k] * tke_0 - 6.0 * to_yyyy_zz[k] * tbe_0 + 4.0 * to_yyyy_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 570-580 components of targeted buffer : FF

        auto to_y_z_yyz_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 70);

        auto to_y_z_yyz_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 71);

        auto to_y_z_yyz_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 72);

        auto to_y_z_yyz_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 73);

        auto to_y_z_yyz_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 74);

        auto to_y_z_yyz_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 75);

        auto to_y_z_yyz_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 76);

        auto to_y_z_yyz_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 77);

        auto to_y_z_yyz_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 78);

        auto to_y_z_yyz_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_y_z_yyz_xxx, to_y_z_yyz_xxy, to_y_z_yyz_xxz, to_y_z_yyz_xyy, to_y_z_yyz_xyz, to_y_z_yyz_xzz, to_y_z_yyz_yyy, to_y_z_yyz_yyz, to_y_z_yyz_yzz, to_y_z_yyz_zzz, to_yyyz_xx, to_yyyz_xxxz, to_yyyz_xxyz, to_yyyz_xxzz, to_yyyz_xy, to_yyyz_xyyz, to_yyyz_xyzz, to_yyyz_xz, to_yyyz_xzzz, to_yyyz_yy, to_yyyz_yyyz, to_yyyz_yyzz, to_yyyz_yz, to_yyyz_yzzz, to_yyyz_zz, to_yyyz_zzzz, to_yz_xx, to_yz_xxxz, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, to_yz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyz_xxx[k] = -4.0 * to_yz_xxxz[k] * tke_0 + 4.0 * to_yyyz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xxy[k] = -4.0 * to_yz_xxyz[k] * tke_0 + 4.0 * to_yyyz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xxz[k] = 2.0 * to_yz_xx[k] - 4.0 * to_yz_xxzz[k] * tke_0 - 2.0 * to_yyyz_xx[k] * tbe_0 + 4.0 * to_yyyz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xyy[k] = -4.0 * to_yz_xyyz[k] * tke_0 + 4.0 * to_yyyz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xyz[k] = 2.0 * to_yz_xy[k] - 4.0 * to_yz_xyzz[k] * tke_0 - 2.0 * to_yyyz_xy[k] * tbe_0 + 4.0 * to_yyyz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_xzz[k] = 4.0 * to_yz_xz[k] - 4.0 * to_yz_xzzz[k] * tke_0 - 4.0 * to_yyyz_xz[k] * tbe_0 + 4.0 * to_yyyz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_yyy[k] = -4.0 * to_yz_yyyz[k] * tke_0 + 4.0 * to_yyyz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_yyz_yyz[k] = 2.0 * to_yz_yy[k] - 4.0 * to_yz_yyzz[k] * tke_0 - 2.0 * to_yyyz_yy[k] * tbe_0 + 4.0 * to_yyyz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_yzz[k] = 4.0 * to_yz_yz[k] - 4.0 * to_yz_yzzz[k] * tke_0 - 4.0 * to_yyyz_yz[k] * tbe_0 + 4.0 * to_yyyz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_yyz_zzz[k] = 6.0 * to_yz_zz[k] - 4.0 * to_yz_zzzz[k] * tke_0 - 6.0 * to_yyyz_zz[k] * tbe_0 + 4.0 * to_yyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 580-590 components of targeted buffer : FF

        auto to_y_z_yzz_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 80);

        auto to_y_z_yzz_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 81);

        auto to_y_z_yzz_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 82);

        auto to_y_z_yzz_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 83);

        auto to_y_z_yzz_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 84);

        auto to_y_z_yzz_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 85);

        auto to_y_z_yzz_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 86);

        auto to_y_z_yzz_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 87);

        auto to_y_z_yzz_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 88);

        auto to_y_z_yzz_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_y_z_yzz_xxx, to_y_z_yzz_xxy, to_y_z_yzz_xxz, to_y_z_yzz_xyy, to_y_z_yzz_xyz, to_y_z_yzz_xzz, to_y_z_yzz_yyy, to_y_z_yzz_yyz, to_y_z_yzz_yzz, to_y_z_yzz_zzz, to_yyzz_xx, to_yyzz_xxxz, to_yyzz_xxyz, to_yyzz_xxzz, to_yyzz_xy, to_yyzz_xyyz, to_yyzz_xyzz, to_yyzz_xz, to_yyzz_xzzz, to_yyzz_yy, to_yyzz_yyyz, to_yyzz_yyzz, to_yyzz_yz, to_yyzz_yzzz, to_yyzz_zz, to_yyzz_zzzz, to_zz_xx, to_zz_xxxz, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, to_zz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yzz_xxx[k] = -2.0 * to_zz_xxxz[k] * tke_0 + 4.0 * to_yyzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xxy[k] = -2.0 * to_zz_xxyz[k] * tke_0 + 4.0 * to_yyzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xxz[k] = to_zz_xx[k] - 2.0 * to_zz_xxzz[k] * tke_0 - 2.0 * to_yyzz_xx[k] * tbe_0 + 4.0 * to_yyzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xyy[k] = -2.0 * to_zz_xyyz[k] * tke_0 + 4.0 * to_yyzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xyz[k] = to_zz_xy[k] - 2.0 * to_zz_xyzz[k] * tke_0 - 2.0 * to_yyzz_xy[k] * tbe_0 + 4.0 * to_yyzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_xzz[k] = 2.0 * to_zz_xz[k] - 2.0 * to_zz_xzzz[k] * tke_0 - 4.0 * to_yyzz_xz[k] * tbe_0 + 4.0 * to_yyzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_yyy[k] = -2.0 * to_zz_yyyz[k] * tke_0 + 4.0 * to_yyzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_yzz_yyz[k] = to_zz_yy[k] - 2.0 * to_zz_yyzz[k] * tke_0 - 2.0 * to_yyzz_yy[k] * tbe_0 + 4.0 * to_yyzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_yzz[k] = 2.0 * to_zz_yz[k] - 2.0 * to_zz_yzzz[k] * tke_0 - 4.0 * to_yyzz_yz[k] * tbe_0 + 4.0 * to_yyzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_yzz_zzz[k] = 3.0 * to_zz_zz[k] - 2.0 * to_zz_zzzz[k] * tke_0 - 6.0 * to_yyzz_zz[k] * tbe_0 + 4.0 * to_yyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 590-600 components of targeted buffer : FF

        auto to_y_z_zzz_xxx = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 90);

        auto to_y_z_zzz_xxy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 91);

        auto to_y_z_zzz_xxz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 92);

        auto to_y_z_zzz_xyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 93);

        auto to_y_z_zzz_xyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 94);

        auto to_y_z_zzz_xzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 95);

        auto to_y_z_zzz_yyy = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 96);

        auto to_y_z_zzz_yyz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 97);

        auto to_y_z_zzz_yzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 98);

        auto to_y_z_zzz_zzz = pbuffer.data(idx_op_geom_101_ff + 5 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_y_z_zzz_xxx, to_y_z_zzz_xxy, to_y_z_zzz_xxz, to_y_z_zzz_xyy, to_y_z_zzz_xyz, to_y_z_zzz_xzz, to_y_z_zzz_yyy, to_y_z_zzz_yyz, to_y_z_zzz_yzz, to_y_z_zzz_zzz, to_yzzz_xx, to_yzzz_xxxz, to_yzzz_xxyz, to_yzzz_xxzz, to_yzzz_xy, to_yzzz_xyyz, to_yzzz_xyzz, to_yzzz_xz, to_yzzz_xzzz, to_yzzz_yy, to_yzzz_yyyz, to_yzzz_yyzz, to_yzzz_yz, to_yzzz_yzzz, to_yzzz_zz, to_yzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zzz_xxx[k] = 4.0 * to_yzzz_xxxz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xxy[k] = 4.0 * to_yzzz_xxyz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xxz[k] = -2.0 * to_yzzz_xx[k] * tbe_0 + 4.0 * to_yzzz_xxzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xyy[k] = 4.0 * to_yzzz_xyyz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xyz[k] = -2.0 * to_yzzz_xy[k] * tbe_0 + 4.0 * to_yzzz_xyzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_xzz[k] = -4.0 * to_yzzz_xz[k] * tbe_0 + 4.0 * to_yzzz_xzzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_yyy[k] = 4.0 * to_yzzz_yyyz[k] * tbe_0 * tke_0;

            to_y_z_zzz_yyz[k] = -2.0 * to_yzzz_yy[k] * tbe_0 + 4.0 * to_yzzz_yyzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_yzz[k] = -4.0 * to_yzzz_yz[k] * tbe_0 + 4.0 * to_yzzz_yzzz[k] * tbe_0 * tke_0;

            to_y_z_zzz_zzz[k] = -6.0 * to_yzzz_zz[k] * tbe_0 + 4.0 * to_yzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 600-610 components of targeted buffer : FF

        auto to_z_x_xxx_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 0);

        auto to_z_x_xxx_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 1);

        auto to_z_x_xxx_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 2);

        auto to_z_x_xxx_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 3);

        auto to_z_x_xxx_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 4);

        auto to_z_x_xxx_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 5);

        auto to_z_x_xxx_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 6);

        auto to_z_x_xxx_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 7);

        auto to_z_x_xxx_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 8);

        auto to_z_x_xxx_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_xxxz_xx, to_xxxz_xxxx, to_xxxz_xxxy, to_xxxz_xxxz, to_xxxz_xxyy, to_xxxz_xxyz, to_xxxz_xxzz, to_xxxz_xy, to_xxxz_xyyy, to_xxxz_xyyz, to_xxxz_xyzz, to_xxxz_xz, to_xxxz_xzzz, to_xxxz_yy, to_xxxz_yz, to_xxxz_zz, to_z_x_xxx_xxx, to_z_x_xxx_xxy, to_z_x_xxx_xxz, to_z_x_xxx_xyy, to_z_x_xxx_xyz, to_z_x_xxx_xzz, to_z_x_xxx_yyy, to_z_x_xxx_yyz, to_z_x_xxx_yzz, to_z_x_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxx_xxx[k] = -6.0 * to_xxxz_xx[k] * tbe_0 + 4.0 * to_xxxz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xxx_xxy[k] = -4.0 * to_xxxz_xy[k] * tbe_0 + 4.0 * to_xxxz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xxx_xxz[k] = -4.0 * to_xxxz_xz[k] * tbe_0 + 4.0 * to_xxxz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xxx_xyy[k] = -2.0 * to_xxxz_yy[k] * tbe_0 + 4.0 * to_xxxz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xxx_xyz[k] = -2.0 * to_xxxz_yz[k] * tbe_0 + 4.0 * to_xxxz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xxx_xzz[k] = -2.0 * to_xxxz_zz[k] * tbe_0 + 4.0 * to_xxxz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xxx_yyy[k] = 4.0 * to_xxxz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xxx_yyz[k] = 4.0 * to_xxxz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xxx_yzz[k] = 4.0 * to_xxxz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xxx_zzz[k] = 4.0 * to_xxxz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 610-620 components of targeted buffer : FF

        auto to_z_x_xxy_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 10);

        auto to_z_x_xxy_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 11);

        auto to_z_x_xxy_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 12);

        auto to_z_x_xxy_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 13);

        auto to_z_x_xxy_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 14);

        auto to_z_x_xxy_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 15);

        auto to_z_x_xxy_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 16);

        auto to_z_x_xxy_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 17);

        auto to_z_x_xxy_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 18);

        auto to_z_x_xxy_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_xxyz_xx, to_xxyz_xxxx, to_xxyz_xxxy, to_xxyz_xxxz, to_xxyz_xxyy, to_xxyz_xxyz, to_xxyz_xxzz, to_xxyz_xy, to_xxyz_xyyy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_xzzz, to_xxyz_yy, to_xxyz_yz, to_xxyz_zz, to_z_x_xxy_xxx, to_z_x_xxy_xxy, to_z_x_xxy_xxz, to_z_x_xxy_xyy, to_z_x_xxy_xyz, to_z_x_xxy_xzz, to_z_x_xxy_yyy, to_z_x_xxy_yyz, to_z_x_xxy_yzz, to_z_x_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxy_xxx[k] = -6.0 * to_xxyz_xx[k] * tbe_0 + 4.0 * to_xxyz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xxy_xxy[k] = -4.0 * to_xxyz_xy[k] * tbe_0 + 4.0 * to_xxyz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xxy_xxz[k] = -4.0 * to_xxyz_xz[k] * tbe_0 + 4.0 * to_xxyz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xxy_xyy[k] = -2.0 * to_xxyz_yy[k] * tbe_0 + 4.0 * to_xxyz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xxy_xyz[k] = -2.0 * to_xxyz_yz[k] * tbe_0 + 4.0 * to_xxyz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xxy_xzz[k] = -2.0 * to_xxyz_zz[k] * tbe_0 + 4.0 * to_xxyz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xxy_yyy[k] = 4.0 * to_xxyz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xxy_yyz[k] = 4.0 * to_xxyz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xxy_yzz[k] = 4.0 * to_xxyz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xxy_zzz[k] = 4.0 * to_xxyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 620-630 components of targeted buffer : FF

        auto to_z_x_xxz_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 20);

        auto to_z_x_xxz_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 21);

        auto to_z_x_xxz_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 22);

        auto to_z_x_xxz_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 23);

        auto to_z_x_xxz_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 24);

        auto to_z_x_xxz_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 25);

        auto to_z_x_xxz_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 26);

        auto to_z_x_xxz_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 27);

        auto to_z_x_xxz_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 28);

        auto to_z_x_xxz_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_xx_xx, to_xx_xxxx, to_xx_xxxy, to_xx_xxxz, to_xx_xxyy, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yz, to_xx_zz, to_xxzz_xx, to_xxzz_xxxx, to_xxzz_xxxy, to_xxzz_xxxz, to_xxzz_xxyy, to_xxzz_xxyz, to_xxzz_xxzz, to_xxzz_xy, to_xxzz_xyyy, to_xxzz_xyyz, to_xxzz_xyzz, to_xxzz_xz, to_xxzz_xzzz, to_xxzz_yy, to_xxzz_yz, to_xxzz_zz, to_z_x_xxz_xxx, to_z_x_xxz_xxy, to_z_x_xxz_xxz, to_z_x_xxz_xyy, to_z_x_xxz_xyz, to_z_x_xxz_xzz, to_z_x_xxz_yyy, to_z_x_xxz_yyz, to_z_x_xxz_yzz, to_z_x_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxz_xxx[k] = 3.0 * to_xx_xx[k] - 2.0 * to_xx_xxxx[k] * tke_0 - 6.0 * to_xxzz_xx[k] * tbe_0 + 4.0 * to_xxzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xxz_xxy[k] = 2.0 * to_xx_xy[k] - 2.0 * to_xx_xxxy[k] * tke_0 - 4.0 * to_xxzz_xy[k] * tbe_0 + 4.0 * to_xxzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xxz_xxz[k] = 2.0 * to_xx_xz[k] - 2.0 * to_xx_xxxz[k] * tke_0 - 4.0 * to_xxzz_xz[k] * tbe_0 + 4.0 * to_xxzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xxz_xyy[k] = to_xx_yy[k] - 2.0 * to_xx_xxyy[k] * tke_0 - 2.0 * to_xxzz_yy[k] * tbe_0 + 4.0 * to_xxzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xxz_xyz[k] = to_xx_yz[k] - 2.0 * to_xx_xxyz[k] * tke_0 - 2.0 * to_xxzz_yz[k] * tbe_0 + 4.0 * to_xxzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xxz_xzz[k] = to_xx_zz[k] - 2.0 * to_xx_xxzz[k] * tke_0 - 2.0 * to_xxzz_zz[k] * tbe_0 + 4.0 * to_xxzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xxz_yyy[k] = -2.0 * to_xx_xyyy[k] * tke_0 + 4.0 * to_xxzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xxz_yyz[k] = -2.0 * to_xx_xyyz[k] * tke_0 + 4.0 * to_xxzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xxz_yzz[k] = -2.0 * to_xx_xyzz[k] * tke_0 + 4.0 * to_xxzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xxz_zzz[k] = -2.0 * to_xx_xzzz[k] * tke_0 + 4.0 * to_xxzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 630-640 components of targeted buffer : FF

        auto to_z_x_xyy_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 30);

        auto to_z_x_xyy_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 31);

        auto to_z_x_xyy_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 32);

        auto to_z_x_xyy_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 33);

        auto to_z_x_xyy_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 34);

        auto to_z_x_xyy_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 35);

        auto to_z_x_xyy_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 36);

        auto to_z_x_xyy_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 37);

        auto to_z_x_xyy_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 38);

        auto to_z_x_xyy_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_xyyz_xx, to_xyyz_xxxx, to_xyyz_xxxy, to_xyyz_xxxz, to_xyyz_xxyy, to_xyyz_xxyz, to_xyyz_xxzz, to_xyyz_xy, to_xyyz_xyyy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_xzzz, to_xyyz_yy, to_xyyz_yz, to_xyyz_zz, to_z_x_xyy_xxx, to_z_x_xyy_xxy, to_z_x_xyy_xxz, to_z_x_xyy_xyy, to_z_x_xyy_xyz, to_z_x_xyy_xzz, to_z_x_xyy_yyy, to_z_x_xyy_yyz, to_z_x_xyy_yzz, to_z_x_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyy_xxx[k] = -6.0 * to_xyyz_xx[k] * tbe_0 + 4.0 * to_xyyz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xyy_xxy[k] = -4.0 * to_xyyz_xy[k] * tbe_0 + 4.0 * to_xyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xyy_xxz[k] = -4.0 * to_xyyz_xz[k] * tbe_0 + 4.0 * to_xyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xyy_xyy[k] = -2.0 * to_xyyz_yy[k] * tbe_0 + 4.0 * to_xyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xyy_xyz[k] = -2.0 * to_xyyz_yz[k] * tbe_0 + 4.0 * to_xyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xyy_xzz[k] = -2.0 * to_xyyz_zz[k] * tbe_0 + 4.0 * to_xyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xyy_yyy[k] = 4.0 * to_xyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xyy_yyz[k] = 4.0 * to_xyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xyy_yzz[k] = 4.0 * to_xyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xyy_zzz[k] = 4.0 * to_xyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 640-650 components of targeted buffer : FF

        auto to_z_x_xyz_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 40);

        auto to_z_x_xyz_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 41);

        auto to_z_x_xyz_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 42);

        auto to_z_x_xyz_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 43);

        auto to_z_x_xyz_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 44);

        auto to_z_x_xyz_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 45);

        auto to_z_x_xyz_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 46);

        auto to_z_x_xyz_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 47);

        auto to_z_x_xyz_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 48);

        auto to_z_x_xyz_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_xy_xx, to_xy_xxxx, to_xy_xxxy, to_xy_xxxz, to_xy_xxyy, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yz, to_xy_zz, to_xyzz_xx, to_xyzz_xxxx, to_xyzz_xxxy, to_xyzz_xxxz, to_xyzz_xxyy, to_xyzz_xxyz, to_xyzz_xxzz, to_xyzz_xy, to_xyzz_xyyy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_xzzz, to_xyzz_yy, to_xyzz_yz, to_xyzz_zz, to_z_x_xyz_xxx, to_z_x_xyz_xxy, to_z_x_xyz_xxz, to_z_x_xyz_xyy, to_z_x_xyz_xyz, to_z_x_xyz_xzz, to_z_x_xyz_yyy, to_z_x_xyz_yyz, to_z_x_xyz_yzz, to_z_x_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyz_xxx[k] = 3.0 * to_xy_xx[k] - 2.0 * to_xy_xxxx[k] * tke_0 - 6.0 * to_xyzz_xx[k] * tbe_0 + 4.0 * to_xyzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xyz_xxy[k] = 2.0 * to_xy_xy[k] - 2.0 * to_xy_xxxy[k] * tke_0 - 4.0 * to_xyzz_xy[k] * tbe_0 + 4.0 * to_xyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xyz_xxz[k] = 2.0 * to_xy_xz[k] - 2.0 * to_xy_xxxz[k] * tke_0 - 4.0 * to_xyzz_xz[k] * tbe_0 + 4.0 * to_xyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xyz_xyy[k] = to_xy_yy[k] - 2.0 * to_xy_xxyy[k] * tke_0 - 2.0 * to_xyzz_yy[k] * tbe_0 + 4.0 * to_xyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xyz_xyz[k] = to_xy_yz[k] - 2.0 * to_xy_xxyz[k] * tke_0 - 2.0 * to_xyzz_yz[k] * tbe_0 + 4.0 * to_xyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xyz_xzz[k] = to_xy_zz[k] - 2.0 * to_xy_xxzz[k] * tke_0 - 2.0 * to_xyzz_zz[k] * tbe_0 + 4.0 * to_xyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xyz_yyy[k] = -2.0 * to_xy_xyyy[k] * tke_0 + 4.0 * to_xyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xyz_yyz[k] = -2.0 * to_xy_xyyz[k] * tke_0 + 4.0 * to_xyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xyz_yzz[k] = -2.0 * to_xy_xyzz[k] * tke_0 + 4.0 * to_xyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xyz_zzz[k] = -2.0 * to_xy_xzzz[k] * tke_0 + 4.0 * to_xyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 650-660 components of targeted buffer : FF

        auto to_z_x_xzz_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 50);

        auto to_z_x_xzz_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 51);

        auto to_z_x_xzz_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 52);

        auto to_z_x_xzz_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 53);

        auto to_z_x_xzz_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 54);

        auto to_z_x_xzz_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 55);

        auto to_z_x_xzz_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 56);

        auto to_z_x_xzz_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 57);

        auto to_z_x_xzz_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 58);

        auto to_z_x_xzz_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_xz_xx, to_xz_xxxx, to_xz_xxxy, to_xz_xxxz, to_xz_xxyy, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yz, to_xz_zz, to_xzzz_xx, to_xzzz_xxxx, to_xzzz_xxxy, to_xzzz_xxxz, to_xzzz_xxyy, to_xzzz_xxyz, to_xzzz_xxzz, to_xzzz_xy, to_xzzz_xyyy, to_xzzz_xyyz, to_xzzz_xyzz, to_xzzz_xz, to_xzzz_xzzz, to_xzzz_yy, to_xzzz_yz, to_xzzz_zz, to_z_x_xzz_xxx, to_z_x_xzz_xxy, to_z_x_xzz_xxz, to_z_x_xzz_xyy, to_z_x_xzz_xyz, to_z_x_xzz_xzz, to_z_x_xzz_yyy, to_z_x_xzz_yyz, to_z_x_xzz_yzz, to_z_x_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xzz_xxx[k] = 6.0 * to_xz_xx[k] - 4.0 * to_xz_xxxx[k] * tke_0 - 6.0 * to_xzzz_xx[k] * tbe_0 + 4.0 * to_xzzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_xzz_xxy[k] = 4.0 * to_xz_xy[k] - 4.0 * to_xz_xxxy[k] * tke_0 - 4.0 * to_xzzz_xy[k] * tbe_0 + 4.0 * to_xzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_xzz_xxz[k] = 4.0 * to_xz_xz[k] - 4.0 * to_xz_xxxz[k] * tke_0 - 4.0 * to_xzzz_xz[k] * tbe_0 + 4.0 * to_xzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_xzz_xyy[k] = 2.0 * to_xz_yy[k] - 4.0 * to_xz_xxyy[k] * tke_0 - 2.0 * to_xzzz_yy[k] * tbe_0 + 4.0 * to_xzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_xzz_xyz[k] = 2.0 * to_xz_yz[k] - 4.0 * to_xz_xxyz[k] * tke_0 - 2.0 * to_xzzz_yz[k] * tbe_0 + 4.0 * to_xzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_xzz_xzz[k] = 2.0 * to_xz_zz[k] - 4.0 * to_xz_xxzz[k] * tke_0 - 2.0 * to_xzzz_zz[k] * tbe_0 + 4.0 * to_xzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_xzz_yyy[k] = -4.0 * to_xz_xyyy[k] * tke_0 + 4.0 * to_xzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_xzz_yyz[k] = -4.0 * to_xz_xyyz[k] * tke_0 + 4.0 * to_xzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_xzz_yzz[k] = -4.0 * to_xz_xyzz[k] * tke_0 + 4.0 * to_xzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_xzz_zzz[k] = -4.0 * to_xz_xzzz[k] * tke_0 + 4.0 * to_xzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 660-670 components of targeted buffer : FF

        auto to_z_x_yyy_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 60);

        auto to_z_x_yyy_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 61);

        auto to_z_x_yyy_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 62);

        auto to_z_x_yyy_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 63);

        auto to_z_x_yyy_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 64);

        auto to_z_x_yyy_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 65);

        auto to_z_x_yyy_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 66);

        auto to_z_x_yyy_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 67);

        auto to_z_x_yyy_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 68);

        auto to_z_x_yyy_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_yyyz_xx, to_yyyz_xxxx, to_yyyz_xxxy, to_yyyz_xxxz, to_yyyz_xxyy, to_yyyz_xxyz, to_yyyz_xxzz, to_yyyz_xy, to_yyyz_xyyy, to_yyyz_xyyz, to_yyyz_xyzz, to_yyyz_xz, to_yyyz_xzzz, to_yyyz_yy, to_yyyz_yz, to_yyyz_zz, to_z_x_yyy_xxx, to_z_x_yyy_xxy, to_z_x_yyy_xxz, to_z_x_yyy_xyy, to_z_x_yyy_xyz, to_z_x_yyy_xzz, to_z_x_yyy_yyy, to_z_x_yyy_yyz, to_z_x_yyy_yzz, to_z_x_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyy_xxx[k] = -6.0 * to_yyyz_xx[k] * tbe_0 + 4.0 * to_yyyz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_yyy_xxy[k] = -4.0 * to_yyyz_xy[k] * tbe_0 + 4.0 * to_yyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_yyy_xxz[k] = -4.0 * to_yyyz_xz[k] * tbe_0 + 4.0 * to_yyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_yyy_xyy[k] = -2.0 * to_yyyz_yy[k] * tbe_0 + 4.0 * to_yyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_yyy_xyz[k] = -2.0 * to_yyyz_yz[k] * tbe_0 + 4.0 * to_yyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_yyy_xzz[k] = -2.0 * to_yyyz_zz[k] * tbe_0 + 4.0 * to_yyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_yyy_yyy[k] = 4.0 * to_yyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_yyy_yyz[k] = 4.0 * to_yyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_yyy_yzz[k] = 4.0 * to_yyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_yyy_zzz[k] = 4.0 * to_yyyz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 670-680 components of targeted buffer : FF

        auto to_z_x_yyz_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 70);

        auto to_z_x_yyz_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 71);

        auto to_z_x_yyz_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 72);

        auto to_z_x_yyz_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 73);

        auto to_z_x_yyz_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 74);

        auto to_z_x_yyz_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 75);

        auto to_z_x_yyz_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 76);

        auto to_z_x_yyz_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 77);

        auto to_z_x_yyz_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 78);

        auto to_z_x_yyz_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_yy_xx, to_yy_xxxx, to_yy_xxxy, to_yy_xxxz, to_yy_xxyy, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yz, to_yy_zz, to_yyzz_xx, to_yyzz_xxxx, to_yyzz_xxxy, to_yyzz_xxxz, to_yyzz_xxyy, to_yyzz_xxyz, to_yyzz_xxzz, to_yyzz_xy, to_yyzz_xyyy, to_yyzz_xyyz, to_yyzz_xyzz, to_yyzz_xz, to_yyzz_xzzz, to_yyzz_yy, to_yyzz_yz, to_yyzz_zz, to_z_x_yyz_xxx, to_z_x_yyz_xxy, to_z_x_yyz_xxz, to_z_x_yyz_xyy, to_z_x_yyz_xyz, to_z_x_yyz_xzz, to_z_x_yyz_yyy, to_z_x_yyz_yyz, to_z_x_yyz_yzz, to_z_x_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyz_xxx[k] = 3.0 * to_yy_xx[k] - 2.0 * to_yy_xxxx[k] * tke_0 - 6.0 * to_yyzz_xx[k] * tbe_0 + 4.0 * to_yyzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_yyz_xxy[k] = 2.0 * to_yy_xy[k] - 2.0 * to_yy_xxxy[k] * tke_0 - 4.0 * to_yyzz_xy[k] * tbe_0 + 4.0 * to_yyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_yyz_xxz[k] = 2.0 * to_yy_xz[k] - 2.0 * to_yy_xxxz[k] * tke_0 - 4.0 * to_yyzz_xz[k] * tbe_0 + 4.0 * to_yyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_yyz_xyy[k] = to_yy_yy[k] - 2.0 * to_yy_xxyy[k] * tke_0 - 2.0 * to_yyzz_yy[k] * tbe_0 + 4.0 * to_yyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_yyz_xyz[k] = to_yy_yz[k] - 2.0 * to_yy_xxyz[k] * tke_0 - 2.0 * to_yyzz_yz[k] * tbe_0 + 4.0 * to_yyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_yyz_xzz[k] = to_yy_zz[k] - 2.0 * to_yy_xxzz[k] * tke_0 - 2.0 * to_yyzz_zz[k] * tbe_0 + 4.0 * to_yyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_yyz_yyy[k] = -2.0 * to_yy_xyyy[k] * tke_0 + 4.0 * to_yyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_yyz_yyz[k] = -2.0 * to_yy_xyyz[k] * tke_0 + 4.0 * to_yyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_yyz_yzz[k] = -2.0 * to_yy_xyzz[k] * tke_0 + 4.0 * to_yyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_yyz_zzz[k] = -2.0 * to_yy_xzzz[k] * tke_0 + 4.0 * to_yyzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 680-690 components of targeted buffer : FF

        auto to_z_x_yzz_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 80);

        auto to_z_x_yzz_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 81);

        auto to_z_x_yzz_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 82);

        auto to_z_x_yzz_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 83);

        auto to_z_x_yzz_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 84);

        auto to_z_x_yzz_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 85);

        auto to_z_x_yzz_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 86);

        auto to_z_x_yzz_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 87);

        auto to_z_x_yzz_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 88);

        auto to_z_x_yzz_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_yz_xx, to_yz_xxxx, to_yz_xxxy, to_yz_xxxz, to_yz_xxyy, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yz, to_yz_zz, to_yzzz_xx, to_yzzz_xxxx, to_yzzz_xxxy, to_yzzz_xxxz, to_yzzz_xxyy, to_yzzz_xxyz, to_yzzz_xxzz, to_yzzz_xy, to_yzzz_xyyy, to_yzzz_xyyz, to_yzzz_xyzz, to_yzzz_xz, to_yzzz_xzzz, to_yzzz_yy, to_yzzz_yz, to_yzzz_zz, to_z_x_yzz_xxx, to_z_x_yzz_xxy, to_z_x_yzz_xxz, to_z_x_yzz_xyy, to_z_x_yzz_xyz, to_z_x_yzz_xzz, to_z_x_yzz_yyy, to_z_x_yzz_yyz, to_z_x_yzz_yzz, to_z_x_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yzz_xxx[k] = 6.0 * to_yz_xx[k] - 4.0 * to_yz_xxxx[k] * tke_0 - 6.0 * to_yzzz_xx[k] * tbe_0 + 4.0 * to_yzzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_yzz_xxy[k] = 4.0 * to_yz_xy[k] - 4.0 * to_yz_xxxy[k] * tke_0 - 4.0 * to_yzzz_xy[k] * tbe_0 + 4.0 * to_yzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_yzz_xxz[k] = 4.0 * to_yz_xz[k] - 4.0 * to_yz_xxxz[k] * tke_0 - 4.0 * to_yzzz_xz[k] * tbe_0 + 4.0 * to_yzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_yzz_xyy[k] = 2.0 * to_yz_yy[k] - 4.0 * to_yz_xxyy[k] * tke_0 - 2.0 * to_yzzz_yy[k] * tbe_0 + 4.0 * to_yzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_yzz_xyz[k] = 2.0 * to_yz_yz[k] - 4.0 * to_yz_xxyz[k] * tke_0 - 2.0 * to_yzzz_yz[k] * tbe_0 + 4.0 * to_yzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_yzz_xzz[k] = 2.0 * to_yz_zz[k] - 4.0 * to_yz_xxzz[k] * tke_0 - 2.0 * to_yzzz_zz[k] * tbe_0 + 4.0 * to_yzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_yzz_yyy[k] = -4.0 * to_yz_xyyy[k] * tke_0 + 4.0 * to_yzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_yzz_yyz[k] = -4.0 * to_yz_xyyz[k] * tke_0 + 4.0 * to_yzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_yzz_yzz[k] = -4.0 * to_yz_xyzz[k] * tke_0 + 4.0 * to_yzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_yzz_zzz[k] = -4.0 * to_yz_xzzz[k] * tke_0 + 4.0 * to_yzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 690-700 components of targeted buffer : FF

        auto to_z_x_zzz_xxx = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 90);

        auto to_z_x_zzz_xxy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 91);

        auto to_z_x_zzz_xxz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 92);

        auto to_z_x_zzz_xyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 93);

        auto to_z_x_zzz_xyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 94);

        auto to_z_x_zzz_xzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 95);

        auto to_z_x_zzz_yyy = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 96);

        auto to_z_x_zzz_yyz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 97);

        auto to_z_x_zzz_yzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 98);

        auto to_z_x_zzz_zzz = pbuffer.data(idx_op_geom_101_ff + 6 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_z_x_zzz_xxx, to_z_x_zzz_xxy, to_z_x_zzz_xxz, to_z_x_zzz_xyy, to_z_x_zzz_xyz, to_z_x_zzz_xzz, to_z_x_zzz_yyy, to_z_x_zzz_yyz, to_z_x_zzz_yzz, to_z_x_zzz_zzz, to_zz_xx, to_zz_xxxx, to_zz_xxxy, to_zz_xxxz, to_zz_xxyy, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yz, to_zz_zz, to_zzzz_xx, to_zzzz_xxxx, to_zzzz_xxxy, to_zzzz_xxxz, to_zzzz_xxyy, to_zzzz_xxyz, to_zzzz_xxzz, to_zzzz_xy, to_zzzz_xyyy, to_zzzz_xyyz, to_zzzz_xyzz, to_zzzz_xz, to_zzzz_xzzz, to_zzzz_yy, to_zzzz_yz, to_zzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zzz_xxx[k] = 9.0 * to_zz_xx[k] - 6.0 * to_zz_xxxx[k] * tke_0 - 6.0 * to_zzzz_xx[k] * tbe_0 + 4.0 * to_zzzz_xxxx[k] * tbe_0 * tke_0;

            to_z_x_zzz_xxy[k] = 6.0 * to_zz_xy[k] - 6.0 * to_zz_xxxy[k] * tke_0 - 4.0 * to_zzzz_xy[k] * tbe_0 + 4.0 * to_zzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_x_zzz_xxz[k] = 6.0 * to_zz_xz[k] - 6.0 * to_zz_xxxz[k] * tke_0 - 4.0 * to_zzzz_xz[k] * tbe_0 + 4.0 * to_zzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_x_zzz_xyy[k] = 3.0 * to_zz_yy[k] - 6.0 * to_zz_xxyy[k] * tke_0 - 2.0 * to_zzzz_yy[k] * tbe_0 + 4.0 * to_zzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_x_zzz_xyz[k] = 3.0 * to_zz_yz[k] - 6.0 * to_zz_xxyz[k] * tke_0 - 2.0 * to_zzzz_yz[k] * tbe_0 + 4.0 * to_zzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_x_zzz_xzz[k] = 3.0 * to_zz_zz[k] - 6.0 * to_zz_xxzz[k] * tke_0 - 2.0 * to_zzzz_zz[k] * tbe_0 + 4.0 * to_zzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_x_zzz_yyy[k] = -6.0 * to_zz_xyyy[k] * tke_0 + 4.0 * to_zzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_x_zzz_yyz[k] = -6.0 * to_zz_xyyz[k] * tke_0 + 4.0 * to_zzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_x_zzz_yzz[k] = -6.0 * to_zz_xyzz[k] * tke_0 + 4.0 * to_zzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_x_zzz_zzz[k] = -6.0 * to_zz_xzzz[k] * tke_0 + 4.0 * to_zzzz_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 700-710 components of targeted buffer : FF

        auto to_z_y_xxx_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 0);

        auto to_z_y_xxx_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 1);

        auto to_z_y_xxx_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 2);

        auto to_z_y_xxx_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 3);

        auto to_z_y_xxx_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 4);

        auto to_z_y_xxx_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 5);

        auto to_z_y_xxx_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 6);

        auto to_z_y_xxx_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 7);

        auto to_z_y_xxx_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 8);

        auto to_z_y_xxx_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_xxxz_xx, to_xxxz_xxxy, to_xxxz_xxyy, to_xxxz_xxyz, to_xxxz_xy, to_xxxz_xyyy, to_xxxz_xyyz, to_xxxz_xyzz, to_xxxz_xz, to_xxxz_yy, to_xxxz_yyyy, to_xxxz_yyyz, to_xxxz_yyzz, to_xxxz_yz, to_xxxz_yzzz, to_xxxz_zz, to_z_y_xxx_xxx, to_z_y_xxx_xxy, to_z_y_xxx_xxz, to_z_y_xxx_xyy, to_z_y_xxx_xyz, to_z_y_xxx_xzz, to_z_y_xxx_yyy, to_z_y_xxx_yyz, to_z_y_xxx_yzz, to_z_y_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxx_xxx[k] = 4.0 * to_xxxz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xxx_xxy[k] = -2.0 * to_xxxz_xx[k] * tbe_0 + 4.0 * to_xxxz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xxx_xxz[k] = 4.0 * to_xxxz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xxx_xyy[k] = -4.0 * to_xxxz_xy[k] * tbe_0 + 4.0 * to_xxxz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xxx_xyz[k] = -2.0 * to_xxxz_xz[k] * tbe_0 + 4.0 * to_xxxz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xxx_xzz[k] = 4.0 * to_xxxz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xxx_yyy[k] = -6.0 * to_xxxz_yy[k] * tbe_0 + 4.0 * to_xxxz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xxx_yyz[k] = -4.0 * to_xxxz_yz[k] * tbe_0 + 4.0 * to_xxxz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xxx_yzz[k] = -2.0 * to_xxxz_zz[k] * tbe_0 + 4.0 * to_xxxz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xxx_zzz[k] = 4.0 * to_xxxz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 710-720 components of targeted buffer : FF

        auto to_z_y_xxy_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 10);

        auto to_z_y_xxy_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 11);

        auto to_z_y_xxy_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 12);

        auto to_z_y_xxy_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 13);

        auto to_z_y_xxy_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 14);

        auto to_z_y_xxy_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 15);

        auto to_z_y_xxy_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 16);

        auto to_z_y_xxy_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 17);

        auto to_z_y_xxy_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 18);

        auto to_z_y_xxy_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_xxyz_xx, to_xxyz_xxxy, to_xxyz_xxyy, to_xxyz_xxyz, to_xxyz_xy, to_xxyz_xyyy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_yy, to_xxyz_yyyy, to_xxyz_yyyz, to_xxyz_yyzz, to_xxyz_yz, to_xxyz_yzzz, to_xxyz_zz, to_z_y_xxy_xxx, to_z_y_xxy_xxy, to_z_y_xxy_xxz, to_z_y_xxy_xyy, to_z_y_xxy_xyz, to_z_y_xxy_xzz, to_z_y_xxy_yyy, to_z_y_xxy_yyz, to_z_y_xxy_yzz, to_z_y_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxy_xxx[k] = 4.0 * to_xxyz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xxy_xxy[k] = -2.0 * to_xxyz_xx[k] * tbe_0 + 4.0 * to_xxyz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xxy_xxz[k] = 4.0 * to_xxyz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xxy_xyy[k] = -4.0 * to_xxyz_xy[k] * tbe_0 + 4.0 * to_xxyz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xxy_xyz[k] = -2.0 * to_xxyz_xz[k] * tbe_0 + 4.0 * to_xxyz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xxy_xzz[k] = 4.0 * to_xxyz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xxy_yyy[k] = -6.0 * to_xxyz_yy[k] * tbe_0 + 4.0 * to_xxyz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xxy_yyz[k] = -4.0 * to_xxyz_yz[k] * tbe_0 + 4.0 * to_xxyz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xxy_yzz[k] = -2.0 * to_xxyz_zz[k] * tbe_0 + 4.0 * to_xxyz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xxy_zzz[k] = 4.0 * to_xxyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 720-730 components of targeted buffer : FF

        auto to_z_y_xxz_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 20);

        auto to_z_y_xxz_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 21);

        auto to_z_y_xxz_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 22);

        auto to_z_y_xxz_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 23);

        auto to_z_y_xxz_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 24);

        auto to_z_y_xxz_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 25);

        auto to_z_y_xxz_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 26);

        auto to_z_y_xxz_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 27);

        auto to_z_y_xxz_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 28);

        auto to_z_y_xxz_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_xx_xx, to_xx_xxxy, to_xx_xxyy, to_xx_xxyz, to_xx_xy, to_xx_xyyy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_yy, to_xx_yyyy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, to_xxzz_xx, to_xxzz_xxxy, to_xxzz_xxyy, to_xxzz_xxyz, to_xxzz_xy, to_xxzz_xyyy, to_xxzz_xyyz, to_xxzz_xyzz, to_xxzz_xz, to_xxzz_yy, to_xxzz_yyyy, to_xxzz_yyyz, to_xxzz_yyzz, to_xxzz_yz, to_xxzz_yzzz, to_xxzz_zz, to_z_y_xxz_xxx, to_z_y_xxz_xxy, to_z_y_xxz_xxz, to_z_y_xxz_xyy, to_z_y_xxz_xyz, to_z_y_xxz_xzz, to_z_y_xxz_yyy, to_z_y_xxz_yyz, to_z_y_xxz_yzz, to_z_y_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxz_xxx[k] = -2.0 * to_xx_xxxy[k] * tke_0 + 4.0 * to_xxzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xxz_xxy[k] = to_xx_xx[k] - 2.0 * to_xx_xxyy[k] * tke_0 - 2.0 * to_xxzz_xx[k] * tbe_0 + 4.0 * to_xxzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xxz_xxz[k] = -2.0 * to_xx_xxyz[k] * tke_0 + 4.0 * to_xxzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xxz_xyy[k] = 2.0 * to_xx_xy[k] - 2.0 * to_xx_xyyy[k] * tke_0 - 4.0 * to_xxzz_xy[k] * tbe_0 + 4.0 * to_xxzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xxz_xyz[k] = to_xx_xz[k] - 2.0 * to_xx_xyyz[k] * tke_0 - 2.0 * to_xxzz_xz[k] * tbe_0 + 4.0 * to_xxzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xxz_xzz[k] = -2.0 * to_xx_xyzz[k] * tke_0 + 4.0 * to_xxzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xxz_yyy[k] = 3.0 * to_xx_yy[k] - 2.0 * to_xx_yyyy[k] * tke_0 - 6.0 * to_xxzz_yy[k] * tbe_0 + 4.0 * to_xxzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xxz_yyz[k] = 2.0 * to_xx_yz[k] - 2.0 * to_xx_yyyz[k] * tke_0 - 4.0 * to_xxzz_yz[k] * tbe_0 + 4.0 * to_xxzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xxz_yzz[k] = to_xx_zz[k] - 2.0 * to_xx_yyzz[k] * tke_0 - 2.0 * to_xxzz_zz[k] * tbe_0 + 4.0 * to_xxzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xxz_zzz[k] = -2.0 * to_xx_yzzz[k] * tke_0 + 4.0 * to_xxzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 730-740 components of targeted buffer : FF

        auto to_z_y_xyy_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 30);

        auto to_z_y_xyy_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 31);

        auto to_z_y_xyy_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 32);

        auto to_z_y_xyy_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 33);

        auto to_z_y_xyy_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 34);

        auto to_z_y_xyy_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 35);

        auto to_z_y_xyy_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 36);

        auto to_z_y_xyy_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 37);

        auto to_z_y_xyy_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 38);

        auto to_z_y_xyy_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_xyyz_xx, to_xyyz_xxxy, to_xyyz_xxyy, to_xyyz_xxyz, to_xyyz_xy, to_xyyz_xyyy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_yy, to_xyyz_yyyy, to_xyyz_yyyz, to_xyyz_yyzz, to_xyyz_yz, to_xyyz_yzzz, to_xyyz_zz, to_z_y_xyy_xxx, to_z_y_xyy_xxy, to_z_y_xyy_xxz, to_z_y_xyy_xyy, to_z_y_xyy_xyz, to_z_y_xyy_xzz, to_z_y_xyy_yyy, to_z_y_xyy_yyz, to_z_y_xyy_yzz, to_z_y_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyy_xxx[k] = 4.0 * to_xyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xyy_xxy[k] = -2.0 * to_xyyz_xx[k] * tbe_0 + 4.0 * to_xyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xyy_xxz[k] = 4.0 * to_xyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xyy_xyy[k] = -4.0 * to_xyyz_xy[k] * tbe_0 + 4.0 * to_xyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xyy_xyz[k] = -2.0 * to_xyyz_xz[k] * tbe_0 + 4.0 * to_xyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xyy_xzz[k] = 4.0 * to_xyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xyy_yyy[k] = -6.0 * to_xyyz_yy[k] * tbe_0 + 4.0 * to_xyyz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xyy_yyz[k] = -4.0 * to_xyyz_yz[k] * tbe_0 + 4.0 * to_xyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xyy_yzz[k] = -2.0 * to_xyyz_zz[k] * tbe_0 + 4.0 * to_xyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xyy_zzz[k] = 4.0 * to_xyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 740-750 components of targeted buffer : FF

        auto to_z_y_xyz_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 40);

        auto to_z_y_xyz_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 41);

        auto to_z_y_xyz_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 42);

        auto to_z_y_xyz_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 43);

        auto to_z_y_xyz_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 44);

        auto to_z_y_xyz_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 45);

        auto to_z_y_xyz_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 46);

        auto to_z_y_xyz_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 47);

        auto to_z_y_xyz_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 48);

        auto to_z_y_xyz_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_xy_xx, to_xy_xxxy, to_xy_xxyy, to_xy_xxyz, to_xy_xy, to_xy_xyyy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_yy, to_xy_yyyy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, to_xyzz_xx, to_xyzz_xxxy, to_xyzz_xxyy, to_xyzz_xxyz, to_xyzz_xy, to_xyzz_xyyy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_yy, to_xyzz_yyyy, to_xyzz_yyyz, to_xyzz_yyzz, to_xyzz_yz, to_xyzz_yzzz, to_xyzz_zz, to_z_y_xyz_xxx, to_z_y_xyz_xxy, to_z_y_xyz_xxz, to_z_y_xyz_xyy, to_z_y_xyz_xyz, to_z_y_xyz_xzz, to_z_y_xyz_yyy, to_z_y_xyz_yyz, to_z_y_xyz_yzz, to_z_y_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyz_xxx[k] = -2.0 * to_xy_xxxy[k] * tke_0 + 4.0 * to_xyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xyz_xxy[k] = to_xy_xx[k] - 2.0 * to_xy_xxyy[k] * tke_0 - 2.0 * to_xyzz_xx[k] * tbe_0 + 4.0 * to_xyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xyz_xxz[k] = -2.0 * to_xy_xxyz[k] * tke_0 + 4.0 * to_xyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xyz_xyy[k] = 2.0 * to_xy_xy[k] - 2.0 * to_xy_xyyy[k] * tke_0 - 4.0 * to_xyzz_xy[k] * tbe_0 + 4.0 * to_xyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xyz_xyz[k] = to_xy_xz[k] - 2.0 * to_xy_xyyz[k] * tke_0 - 2.0 * to_xyzz_xz[k] * tbe_0 + 4.0 * to_xyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xyz_xzz[k] = -2.0 * to_xy_xyzz[k] * tke_0 + 4.0 * to_xyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xyz_yyy[k] = 3.0 * to_xy_yy[k] - 2.0 * to_xy_yyyy[k] * tke_0 - 6.0 * to_xyzz_yy[k] * tbe_0 + 4.0 * to_xyzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xyz_yyz[k] = 2.0 * to_xy_yz[k] - 2.0 * to_xy_yyyz[k] * tke_0 - 4.0 * to_xyzz_yz[k] * tbe_0 + 4.0 * to_xyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xyz_yzz[k] = to_xy_zz[k] - 2.0 * to_xy_yyzz[k] * tke_0 - 2.0 * to_xyzz_zz[k] * tbe_0 + 4.0 * to_xyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xyz_zzz[k] = -2.0 * to_xy_yzzz[k] * tke_0 + 4.0 * to_xyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 750-760 components of targeted buffer : FF

        auto to_z_y_xzz_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 50);

        auto to_z_y_xzz_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 51);

        auto to_z_y_xzz_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 52);

        auto to_z_y_xzz_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 53);

        auto to_z_y_xzz_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 54);

        auto to_z_y_xzz_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 55);

        auto to_z_y_xzz_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 56);

        auto to_z_y_xzz_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 57);

        auto to_z_y_xzz_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 58);

        auto to_z_y_xzz_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_xz_xx, to_xz_xxxy, to_xz_xxyy, to_xz_xxyz, to_xz_xy, to_xz_xyyy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_yy, to_xz_yyyy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, to_xzzz_xx, to_xzzz_xxxy, to_xzzz_xxyy, to_xzzz_xxyz, to_xzzz_xy, to_xzzz_xyyy, to_xzzz_xyyz, to_xzzz_xyzz, to_xzzz_xz, to_xzzz_yy, to_xzzz_yyyy, to_xzzz_yyyz, to_xzzz_yyzz, to_xzzz_yz, to_xzzz_yzzz, to_xzzz_zz, to_z_y_xzz_xxx, to_z_y_xzz_xxy, to_z_y_xzz_xxz, to_z_y_xzz_xyy, to_z_y_xzz_xyz, to_z_y_xzz_xzz, to_z_y_xzz_yyy, to_z_y_xzz_yyz, to_z_y_xzz_yzz, to_z_y_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xzz_xxx[k] = -4.0 * to_xz_xxxy[k] * tke_0 + 4.0 * to_xzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_xzz_xxy[k] = 2.0 * to_xz_xx[k] - 4.0 * to_xz_xxyy[k] * tke_0 - 2.0 * to_xzzz_xx[k] * tbe_0 + 4.0 * to_xzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_xzz_xxz[k] = -4.0 * to_xz_xxyz[k] * tke_0 + 4.0 * to_xzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_xzz_xyy[k] = 4.0 * to_xz_xy[k] - 4.0 * to_xz_xyyy[k] * tke_0 - 4.0 * to_xzzz_xy[k] * tbe_0 + 4.0 * to_xzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_xzz_xyz[k] = 2.0 * to_xz_xz[k] - 4.0 * to_xz_xyyz[k] * tke_0 - 2.0 * to_xzzz_xz[k] * tbe_0 + 4.0 * to_xzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_xzz_xzz[k] = -4.0 * to_xz_xyzz[k] * tke_0 + 4.0 * to_xzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_xzz_yyy[k] = 6.0 * to_xz_yy[k] - 4.0 * to_xz_yyyy[k] * tke_0 - 6.0 * to_xzzz_yy[k] * tbe_0 + 4.0 * to_xzzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_xzz_yyz[k] = 4.0 * to_xz_yz[k] - 4.0 * to_xz_yyyz[k] * tke_0 - 4.0 * to_xzzz_yz[k] * tbe_0 + 4.0 * to_xzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_xzz_yzz[k] = 2.0 * to_xz_zz[k] - 4.0 * to_xz_yyzz[k] * tke_0 - 2.0 * to_xzzz_zz[k] * tbe_0 + 4.0 * to_xzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_xzz_zzz[k] = -4.0 * to_xz_yzzz[k] * tke_0 + 4.0 * to_xzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 760-770 components of targeted buffer : FF

        auto to_z_y_yyy_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 60);

        auto to_z_y_yyy_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 61);

        auto to_z_y_yyy_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 62);

        auto to_z_y_yyy_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 63);

        auto to_z_y_yyy_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 64);

        auto to_z_y_yyy_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 65);

        auto to_z_y_yyy_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 66);

        auto to_z_y_yyy_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 67);

        auto to_z_y_yyy_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 68);

        auto to_z_y_yyy_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_yyyz_xx, to_yyyz_xxxy, to_yyyz_xxyy, to_yyyz_xxyz, to_yyyz_xy, to_yyyz_xyyy, to_yyyz_xyyz, to_yyyz_xyzz, to_yyyz_xz, to_yyyz_yy, to_yyyz_yyyy, to_yyyz_yyyz, to_yyyz_yyzz, to_yyyz_yz, to_yyyz_yzzz, to_yyyz_zz, to_z_y_yyy_xxx, to_z_y_yyy_xxy, to_z_y_yyy_xxz, to_z_y_yyy_xyy, to_z_y_yyy_xyz, to_z_y_yyy_xzz, to_z_y_yyy_yyy, to_z_y_yyy_yyz, to_z_y_yyy_yzz, to_z_y_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyy_xxx[k] = 4.0 * to_yyyz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_yyy_xxy[k] = -2.0 * to_yyyz_xx[k] * tbe_0 + 4.0 * to_yyyz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_yyy_xxz[k] = 4.0 * to_yyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_yyy_xyy[k] = -4.0 * to_yyyz_xy[k] * tbe_0 + 4.0 * to_yyyz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_yyy_xyz[k] = -2.0 * to_yyyz_xz[k] * tbe_0 + 4.0 * to_yyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_yyy_xzz[k] = 4.0 * to_yyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_yyy_yyy[k] = -6.0 * to_yyyz_yy[k] * tbe_0 + 4.0 * to_yyyz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_yyy_yyz[k] = -4.0 * to_yyyz_yz[k] * tbe_0 + 4.0 * to_yyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_yyy_yzz[k] = -2.0 * to_yyyz_zz[k] * tbe_0 + 4.0 * to_yyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_yyy_zzz[k] = 4.0 * to_yyyz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 770-780 components of targeted buffer : FF

        auto to_z_y_yyz_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 70);

        auto to_z_y_yyz_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 71);

        auto to_z_y_yyz_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 72);

        auto to_z_y_yyz_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 73);

        auto to_z_y_yyz_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 74);

        auto to_z_y_yyz_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 75);

        auto to_z_y_yyz_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 76);

        auto to_z_y_yyz_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 77);

        auto to_z_y_yyz_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 78);

        auto to_z_y_yyz_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_yy_xx, to_yy_xxxy, to_yy_xxyy, to_yy_xxyz, to_yy_xy, to_yy_xyyy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_yy, to_yy_yyyy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, to_yyzz_xx, to_yyzz_xxxy, to_yyzz_xxyy, to_yyzz_xxyz, to_yyzz_xy, to_yyzz_xyyy, to_yyzz_xyyz, to_yyzz_xyzz, to_yyzz_xz, to_yyzz_yy, to_yyzz_yyyy, to_yyzz_yyyz, to_yyzz_yyzz, to_yyzz_yz, to_yyzz_yzzz, to_yyzz_zz, to_z_y_yyz_xxx, to_z_y_yyz_xxy, to_z_y_yyz_xxz, to_z_y_yyz_xyy, to_z_y_yyz_xyz, to_z_y_yyz_xzz, to_z_y_yyz_yyy, to_z_y_yyz_yyz, to_z_y_yyz_yzz, to_z_y_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyz_xxx[k] = -2.0 * to_yy_xxxy[k] * tke_0 + 4.0 * to_yyzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_yyz_xxy[k] = to_yy_xx[k] - 2.0 * to_yy_xxyy[k] * tke_0 - 2.0 * to_yyzz_xx[k] * tbe_0 + 4.0 * to_yyzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_yyz_xxz[k] = -2.0 * to_yy_xxyz[k] * tke_0 + 4.0 * to_yyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_yyz_xyy[k] = 2.0 * to_yy_xy[k] - 2.0 * to_yy_xyyy[k] * tke_0 - 4.0 * to_yyzz_xy[k] * tbe_0 + 4.0 * to_yyzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_yyz_xyz[k] = to_yy_xz[k] - 2.0 * to_yy_xyyz[k] * tke_0 - 2.0 * to_yyzz_xz[k] * tbe_0 + 4.0 * to_yyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_yyz_xzz[k] = -2.0 * to_yy_xyzz[k] * tke_0 + 4.0 * to_yyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_yyz_yyy[k] = 3.0 * to_yy_yy[k] - 2.0 * to_yy_yyyy[k] * tke_0 - 6.0 * to_yyzz_yy[k] * tbe_0 + 4.0 * to_yyzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_yyz_yyz[k] = 2.0 * to_yy_yz[k] - 2.0 * to_yy_yyyz[k] * tke_0 - 4.0 * to_yyzz_yz[k] * tbe_0 + 4.0 * to_yyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_yyz_yzz[k] = to_yy_zz[k] - 2.0 * to_yy_yyzz[k] * tke_0 - 2.0 * to_yyzz_zz[k] * tbe_0 + 4.0 * to_yyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_yyz_zzz[k] = -2.0 * to_yy_yzzz[k] * tke_0 + 4.0 * to_yyzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 780-790 components of targeted buffer : FF

        auto to_z_y_yzz_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 80);

        auto to_z_y_yzz_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 81);

        auto to_z_y_yzz_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 82);

        auto to_z_y_yzz_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 83);

        auto to_z_y_yzz_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 84);

        auto to_z_y_yzz_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 85);

        auto to_z_y_yzz_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 86);

        auto to_z_y_yzz_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 87);

        auto to_z_y_yzz_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 88);

        auto to_z_y_yzz_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_yz_xx, to_yz_xxxy, to_yz_xxyy, to_yz_xxyz, to_yz_xy, to_yz_xyyy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_yy, to_yz_yyyy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, to_yzzz_xx, to_yzzz_xxxy, to_yzzz_xxyy, to_yzzz_xxyz, to_yzzz_xy, to_yzzz_xyyy, to_yzzz_xyyz, to_yzzz_xyzz, to_yzzz_xz, to_yzzz_yy, to_yzzz_yyyy, to_yzzz_yyyz, to_yzzz_yyzz, to_yzzz_yz, to_yzzz_yzzz, to_yzzz_zz, to_z_y_yzz_xxx, to_z_y_yzz_xxy, to_z_y_yzz_xxz, to_z_y_yzz_xyy, to_z_y_yzz_xyz, to_z_y_yzz_xzz, to_z_y_yzz_yyy, to_z_y_yzz_yyz, to_z_y_yzz_yzz, to_z_y_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yzz_xxx[k] = -4.0 * to_yz_xxxy[k] * tke_0 + 4.0 * to_yzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_yzz_xxy[k] = 2.0 * to_yz_xx[k] - 4.0 * to_yz_xxyy[k] * tke_0 - 2.0 * to_yzzz_xx[k] * tbe_0 + 4.0 * to_yzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_yzz_xxz[k] = -4.0 * to_yz_xxyz[k] * tke_0 + 4.0 * to_yzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_yzz_xyy[k] = 4.0 * to_yz_xy[k] - 4.0 * to_yz_xyyy[k] * tke_0 - 4.0 * to_yzzz_xy[k] * tbe_0 + 4.0 * to_yzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_yzz_xyz[k] = 2.0 * to_yz_xz[k] - 4.0 * to_yz_xyyz[k] * tke_0 - 2.0 * to_yzzz_xz[k] * tbe_0 + 4.0 * to_yzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_yzz_xzz[k] = -4.0 * to_yz_xyzz[k] * tke_0 + 4.0 * to_yzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_yzz_yyy[k] = 6.0 * to_yz_yy[k] - 4.0 * to_yz_yyyy[k] * tke_0 - 6.0 * to_yzzz_yy[k] * tbe_0 + 4.0 * to_yzzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_yzz_yyz[k] = 4.0 * to_yz_yz[k] - 4.0 * to_yz_yyyz[k] * tke_0 - 4.0 * to_yzzz_yz[k] * tbe_0 + 4.0 * to_yzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_yzz_yzz[k] = 2.0 * to_yz_zz[k] - 4.0 * to_yz_yyzz[k] * tke_0 - 2.0 * to_yzzz_zz[k] * tbe_0 + 4.0 * to_yzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_yzz_zzz[k] = -4.0 * to_yz_yzzz[k] * tke_0 + 4.0 * to_yzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 790-800 components of targeted buffer : FF

        auto to_z_y_zzz_xxx = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 90);

        auto to_z_y_zzz_xxy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 91);

        auto to_z_y_zzz_xxz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 92);

        auto to_z_y_zzz_xyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 93);

        auto to_z_y_zzz_xyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 94);

        auto to_z_y_zzz_xzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 95);

        auto to_z_y_zzz_yyy = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 96);

        auto to_z_y_zzz_yyz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 97);

        auto to_z_y_zzz_yzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 98);

        auto to_z_y_zzz_zzz = pbuffer.data(idx_op_geom_101_ff + 7 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_z_y_zzz_xxx, to_z_y_zzz_xxy, to_z_y_zzz_xxz, to_z_y_zzz_xyy, to_z_y_zzz_xyz, to_z_y_zzz_xzz, to_z_y_zzz_yyy, to_z_y_zzz_yyz, to_z_y_zzz_yzz, to_z_y_zzz_zzz, to_zz_xx, to_zz_xxxy, to_zz_xxyy, to_zz_xxyz, to_zz_xy, to_zz_xyyy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_yy, to_zz_yyyy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, to_zzzz_xx, to_zzzz_xxxy, to_zzzz_xxyy, to_zzzz_xxyz, to_zzzz_xy, to_zzzz_xyyy, to_zzzz_xyyz, to_zzzz_xyzz, to_zzzz_xz, to_zzzz_yy, to_zzzz_yyyy, to_zzzz_yyyz, to_zzzz_yyzz, to_zzzz_yz, to_zzzz_yzzz, to_zzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zzz_xxx[k] = -6.0 * to_zz_xxxy[k] * tke_0 + 4.0 * to_zzzz_xxxy[k] * tbe_0 * tke_0;

            to_z_y_zzz_xxy[k] = 3.0 * to_zz_xx[k] - 6.0 * to_zz_xxyy[k] * tke_0 - 2.0 * to_zzzz_xx[k] * tbe_0 + 4.0 * to_zzzz_xxyy[k] * tbe_0 * tke_0;

            to_z_y_zzz_xxz[k] = -6.0 * to_zz_xxyz[k] * tke_0 + 4.0 * to_zzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_y_zzz_xyy[k] = 6.0 * to_zz_xy[k] - 6.0 * to_zz_xyyy[k] * tke_0 - 4.0 * to_zzzz_xy[k] * tbe_0 + 4.0 * to_zzzz_xyyy[k] * tbe_0 * tke_0;

            to_z_y_zzz_xyz[k] = 3.0 * to_zz_xz[k] - 6.0 * to_zz_xyyz[k] * tke_0 - 2.0 * to_zzzz_xz[k] * tbe_0 + 4.0 * to_zzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_y_zzz_xzz[k] = -6.0 * to_zz_xyzz[k] * tke_0 + 4.0 * to_zzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_y_zzz_yyy[k] = 9.0 * to_zz_yy[k] - 6.0 * to_zz_yyyy[k] * tke_0 - 6.0 * to_zzzz_yy[k] * tbe_0 + 4.0 * to_zzzz_yyyy[k] * tbe_0 * tke_0;

            to_z_y_zzz_yyz[k] = 6.0 * to_zz_yz[k] - 6.0 * to_zz_yyyz[k] * tke_0 - 4.0 * to_zzzz_yz[k] * tbe_0 + 4.0 * to_zzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_y_zzz_yzz[k] = 3.0 * to_zz_zz[k] - 6.0 * to_zz_yyzz[k] * tke_0 - 2.0 * to_zzzz_zz[k] * tbe_0 + 4.0 * to_zzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_y_zzz_zzz[k] = -6.0 * to_zz_yzzz[k] * tke_0 + 4.0 * to_zzzz_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 800-810 components of targeted buffer : FF

        auto to_z_z_xxx_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 0);

        auto to_z_z_xxx_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 1);

        auto to_z_z_xxx_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 2);

        auto to_z_z_xxx_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 3);

        auto to_z_z_xxx_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 4);

        auto to_z_z_xxx_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 5);

        auto to_z_z_xxx_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 6);

        auto to_z_z_xxx_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 7);

        auto to_z_z_xxx_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 8);

        auto to_z_z_xxx_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 9);

        #pragma omp simd aligned(to_xxxz_xx, to_xxxz_xxxz, to_xxxz_xxyz, to_xxxz_xxzz, to_xxxz_xy, to_xxxz_xyyz, to_xxxz_xyzz, to_xxxz_xz, to_xxxz_xzzz, to_xxxz_yy, to_xxxz_yyyz, to_xxxz_yyzz, to_xxxz_yz, to_xxxz_yzzz, to_xxxz_zz, to_xxxz_zzzz, to_z_z_xxx_xxx, to_z_z_xxx_xxy, to_z_z_xxx_xxz, to_z_z_xxx_xyy, to_z_z_xxx_xyz, to_z_z_xxx_xzz, to_z_z_xxx_yyy, to_z_z_xxx_yyz, to_z_z_xxx_yzz, to_z_z_xxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxx_xxx[k] = 4.0 * to_xxxz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xxy[k] = 4.0 * to_xxxz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xxz[k] = -2.0 * to_xxxz_xx[k] * tbe_0 + 4.0 * to_xxxz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xyy[k] = 4.0 * to_xxxz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xyz[k] = -2.0 * to_xxxz_xy[k] * tbe_0 + 4.0 * to_xxxz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_xzz[k] = -4.0 * to_xxxz_xz[k] * tbe_0 + 4.0 * to_xxxz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_yyy[k] = 4.0 * to_xxxz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xxx_yyz[k] = -2.0 * to_xxxz_yy[k] * tbe_0 + 4.0 * to_xxxz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_yzz[k] = -4.0 * to_xxxz_yz[k] * tbe_0 + 4.0 * to_xxxz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xxx_zzz[k] = -6.0 * to_xxxz_zz[k] * tbe_0 + 4.0 * to_xxxz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 810-820 components of targeted buffer : FF

        auto to_z_z_xxy_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 10);

        auto to_z_z_xxy_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 11);

        auto to_z_z_xxy_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 12);

        auto to_z_z_xxy_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 13);

        auto to_z_z_xxy_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 14);

        auto to_z_z_xxy_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 15);

        auto to_z_z_xxy_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 16);

        auto to_z_z_xxy_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 17);

        auto to_z_z_xxy_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 18);

        auto to_z_z_xxy_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 19);

        #pragma omp simd aligned(to_xxyz_xx, to_xxyz_xxxz, to_xxyz_xxyz, to_xxyz_xxzz, to_xxyz_xy, to_xxyz_xyyz, to_xxyz_xyzz, to_xxyz_xz, to_xxyz_xzzz, to_xxyz_yy, to_xxyz_yyyz, to_xxyz_yyzz, to_xxyz_yz, to_xxyz_yzzz, to_xxyz_zz, to_xxyz_zzzz, to_z_z_xxy_xxx, to_z_z_xxy_xxy, to_z_z_xxy_xxz, to_z_z_xxy_xyy, to_z_z_xxy_xyz, to_z_z_xxy_xzz, to_z_z_xxy_yyy, to_z_z_xxy_yyz, to_z_z_xxy_yzz, to_z_z_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxy_xxx[k] = 4.0 * to_xxyz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xxy[k] = 4.0 * to_xxyz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xxz[k] = -2.0 * to_xxyz_xx[k] * tbe_0 + 4.0 * to_xxyz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xyy[k] = 4.0 * to_xxyz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xyz[k] = -2.0 * to_xxyz_xy[k] * tbe_0 + 4.0 * to_xxyz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_xzz[k] = -4.0 * to_xxyz_xz[k] * tbe_0 + 4.0 * to_xxyz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_yyy[k] = 4.0 * to_xxyz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xxy_yyz[k] = -2.0 * to_xxyz_yy[k] * tbe_0 + 4.0 * to_xxyz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_yzz[k] = -4.0 * to_xxyz_yz[k] * tbe_0 + 4.0 * to_xxyz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xxy_zzz[k] = -6.0 * to_xxyz_zz[k] * tbe_0 + 4.0 * to_xxyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 820-830 components of targeted buffer : FF

        auto to_z_z_xxz_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 20);

        auto to_z_z_xxz_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 21);

        auto to_z_z_xxz_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 22);

        auto to_z_z_xxz_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 23);

        auto to_z_z_xxz_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 24);

        auto to_z_z_xxz_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 25);

        auto to_z_z_xxz_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 26);

        auto to_z_z_xxz_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 27);

        auto to_z_z_xxz_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 28);

        auto to_z_z_xxz_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 29);

        #pragma omp simd aligned(to_xx_xx, to_xx_xxxz, to_xx_xxyz, to_xx_xxzz, to_xx_xy, to_xx_xyyz, to_xx_xyzz, to_xx_xz, to_xx_xzzz, to_xx_yy, to_xx_yyyz, to_xx_yyzz, to_xx_yz, to_xx_yzzz, to_xx_zz, to_xx_zzzz, to_xxzz_xx, to_xxzz_xxxz, to_xxzz_xxyz, to_xxzz_xxzz, to_xxzz_xy, to_xxzz_xyyz, to_xxzz_xyzz, to_xxzz_xz, to_xxzz_xzzz, to_xxzz_yy, to_xxzz_yyyz, to_xxzz_yyzz, to_xxzz_yz, to_xxzz_yzzz, to_xxzz_zz, to_xxzz_zzzz, to_z_z_xxz_xxx, to_z_z_xxz_xxy, to_z_z_xxz_xxz, to_z_z_xxz_xyy, to_z_z_xxz_xyz, to_z_z_xxz_xzz, to_z_z_xxz_yyy, to_z_z_xxz_yyz, to_z_z_xxz_yzz, to_z_z_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxz_xxx[k] = -2.0 * to_xx_xxxz[k] * tke_0 + 4.0 * to_xxzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xxy[k] = -2.0 * to_xx_xxyz[k] * tke_0 + 4.0 * to_xxzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xxz[k] = to_xx_xx[k] - 2.0 * to_xx_xxzz[k] * tke_0 - 2.0 * to_xxzz_xx[k] * tbe_0 + 4.0 * to_xxzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xyy[k] = -2.0 * to_xx_xyyz[k] * tke_0 + 4.0 * to_xxzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xyz[k] = to_xx_xy[k] - 2.0 * to_xx_xyzz[k] * tke_0 - 2.0 * to_xxzz_xy[k] * tbe_0 + 4.0 * to_xxzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_xzz[k] = 2.0 * to_xx_xz[k] - 2.0 * to_xx_xzzz[k] * tke_0 - 4.0 * to_xxzz_xz[k] * tbe_0 + 4.0 * to_xxzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_yyy[k] = -2.0 * to_xx_yyyz[k] * tke_0 + 4.0 * to_xxzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xxz_yyz[k] = to_xx_yy[k] - 2.0 * to_xx_yyzz[k] * tke_0 - 2.0 * to_xxzz_yy[k] * tbe_0 + 4.0 * to_xxzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_yzz[k] = 2.0 * to_xx_yz[k] - 2.0 * to_xx_yzzz[k] * tke_0 - 4.0 * to_xxzz_yz[k] * tbe_0 + 4.0 * to_xxzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xxz_zzz[k] = 3.0 * to_xx_zz[k] - 2.0 * to_xx_zzzz[k] * tke_0 - 6.0 * to_xxzz_zz[k] * tbe_0 + 4.0 * to_xxzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 830-840 components of targeted buffer : FF

        auto to_z_z_xyy_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 30);

        auto to_z_z_xyy_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 31);

        auto to_z_z_xyy_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 32);

        auto to_z_z_xyy_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 33);

        auto to_z_z_xyy_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 34);

        auto to_z_z_xyy_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 35);

        auto to_z_z_xyy_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 36);

        auto to_z_z_xyy_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 37);

        auto to_z_z_xyy_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 38);

        auto to_z_z_xyy_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 39);

        #pragma omp simd aligned(to_xyyz_xx, to_xyyz_xxxz, to_xyyz_xxyz, to_xyyz_xxzz, to_xyyz_xy, to_xyyz_xyyz, to_xyyz_xyzz, to_xyyz_xz, to_xyyz_xzzz, to_xyyz_yy, to_xyyz_yyyz, to_xyyz_yyzz, to_xyyz_yz, to_xyyz_yzzz, to_xyyz_zz, to_xyyz_zzzz, to_z_z_xyy_xxx, to_z_z_xyy_xxy, to_z_z_xyy_xxz, to_z_z_xyy_xyy, to_z_z_xyy_xyz, to_z_z_xyy_xzz, to_z_z_xyy_yyy, to_z_z_xyy_yyz, to_z_z_xyy_yzz, to_z_z_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyy_xxx[k] = 4.0 * to_xyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xxy[k] = 4.0 * to_xyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xxz[k] = -2.0 * to_xyyz_xx[k] * tbe_0 + 4.0 * to_xyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xyy[k] = 4.0 * to_xyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xyz[k] = -2.0 * to_xyyz_xy[k] * tbe_0 + 4.0 * to_xyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_xzz[k] = -4.0 * to_xyyz_xz[k] * tbe_0 + 4.0 * to_xyyz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_yyy[k] = 4.0 * to_xyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xyy_yyz[k] = -2.0 * to_xyyz_yy[k] * tbe_0 + 4.0 * to_xyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_yzz[k] = -4.0 * to_xyyz_yz[k] * tbe_0 + 4.0 * to_xyyz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xyy_zzz[k] = -6.0 * to_xyyz_zz[k] * tbe_0 + 4.0 * to_xyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 840-850 components of targeted buffer : FF

        auto to_z_z_xyz_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 40);

        auto to_z_z_xyz_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 41);

        auto to_z_z_xyz_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 42);

        auto to_z_z_xyz_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 43);

        auto to_z_z_xyz_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 44);

        auto to_z_z_xyz_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 45);

        auto to_z_z_xyz_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 46);

        auto to_z_z_xyz_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 47);

        auto to_z_z_xyz_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 48);

        auto to_z_z_xyz_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 49);

        #pragma omp simd aligned(to_xy_xx, to_xy_xxxz, to_xy_xxyz, to_xy_xxzz, to_xy_xy, to_xy_xyyz, to_xy_xyzz, to_xy_xz, to_xy_xzzz, to_xy_yy, to_xy_yyyz, to_xy_yyzz, to_xy_yz, to_xy_yzzz, to_xy_zz, to_xy_zzzz, to_xyzz_xx, to_xyzz_xxxz, to_xyzz_xxyz, to_xyzz_xxzz, to_xyzz_xy, to_xyzz_xyyz, to_xyzz_xyzz, to_xyzz_xz, to_xyzz_xzzz, to_xyzz_yy, to_xyzz_yyyz, to_xyzz_yyzz, to_xyzz_yz, to_xyzz_yzzz, to_xyzz_zz, to_xyzz_zzzz, to_z_z_xyz_xxx, to_z_z_xyz_xxy, to_z_z_xyz_xxz, to_z_z_xyz_xyy, to_z_z_xyz_xyz, to_z_z_xyz_xzz, to_z_z_xyz_yyy, to_z_z_xyz_yyz, to_z_z_xyz_yzz, to_z_z_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyz_xxx[k] = -2.0 * to_xy_xxxz[k] * tke_0 + 4.0 * to_xyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xxy[k] = -2.0 * to_xy_xxyz[k] * tke_0 + 4.0 * to_xyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xxz[k] = to_xy_xx[k] - 2.0 * to_xy_xxzz[k] * tke_0 - 2.0 * to_xyzz_xx[k] * tbe_0 + 4.0 * to_xyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xyy[k] = -2.0 * to_xy_xyyz[k] * tke_0 + 4.0 * to_xyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xyz[k] = to_xy_xy[k] - 2.0 * to_xy_xyzz[k] * tke_0 - 2.0 * to_xyzz_xy[k] * tbe_0 + 4.0 * to_xyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_xzz[k] = 2.0 * to_xy_xz[k] - 2.0 * to_xy_xzzz[k] * tke_0 - 4.0 * to_xyzz_xz[k] * tbe_0 + 4.0 * to_xyzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_yyy[k] = -2.0 * to_xy_yyyz[k] * tke_0 + 4.0 * to_xyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xyz_yyz[k] = to_xy_yy[k] - 2.0 * to_xy_yyzz[k] * tke_0 - 2.0 * to_xyzz_yy[k] * tbe_0 + 4.0 * to_xyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_yzz[k] = 2.0 * to_xy_yz[k] - 2.0 * to_xy_yzzz[k] * tke_0 - 4.0 * to_xyzz_yz[k] * tbe_0 + 4.0 * to_xyzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xyz_zzz[k] = 3.0 * to_xy_zz[k] - 2.0 * to_xy_zzzz[k] * tke_0 - 6.0 * to_xyzz_zz[k] * tbe_0 + 4.0 * to_xyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 850-860 components of targeted buffer : FF

        auto to_z_z_xzz_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 50);

        auto to_z_z_xzz_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 51);

        auto to_z_z_xzz_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 52);

        auto to_z_z_xzz_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 53);

        auto to_z_z_xzz_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 54);

        auto to_z_z_xzz_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 55);

        auto to_z_z_xzz_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 56);

        auto to_z_z_xzz_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 57);

        auto to_z_z_xzz_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 58);

        auto to_z_z_xzz_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 59);

        #pragma omp simd aligned(to_xz_xx, to_xz_xxxz, to_xz_xxyz, to_xz_xxzz, to_xz_xy, to_xz_xyyz, to_xz_xyzz, to_xz_xz, to_xz_xzzz, to_xz_yy, to_xz_yyyz, to_xz_yyzz, to_xz_yz, to_xz_yzzz, to_xz_zz, to_xz_zzzz, to_xzzz_xx, to_xzzz_xxxz, to_xzzz_xxyz, to_xzzz_xxzz, to_xzzz_xy, to_xzzz_xyyz, to_xzzz_xyzz, to_xzzz_xz, to_xzzz_xzzz, to_xzzz_yy, to_xzzz_yyyz, to_xzzz_yyzz, to_xzzz_yz, to_xzzz_yzzz, to_xzzz_zz, to_xzzz_zzzz, to_z_z_xzz_xxx, to_z_z_xzz_xxy, to_z_z_xzz_xxz, to_z_z_xzz_xyy, to_z_z_xzz_xyz, to_z_z_xzz_xzz, to_z_z_xzz_yyy, to_z_z_xzz_yyz, to_z_z_xzz_yzz, to_z_z_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xzz_xxx[k] = -4.0 * to_xz_xxxz[k] * tke_0 + 4.0 * to_xzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xxy[k] = -4.0 * to_xz_xxyz[k] * tke_0 + 4.0 * to_xzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xxz[k] = 2.0 * to_xz_xx[k] - 4.0 * to_xz_xxzz[k] * tke_0 - 2.0 * to_xzzz_xx[k] * tbe_0 + 4.0 * to_xzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xyy[k] = -4.0 * to_xz_xyyz[k] * tke_0 + 4.0 * to_xzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xyz[k] = 2.0 * to_xz_xy[k] - 4.0 * to_xz_xyzz[k] * tke_0 - 2.0 * to_xzzz_xy[k] * tbe_0 + 4.0 * to_xzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_xzz[k] = 4.0 * to_xz_xz[k] - 4.0 * to_xz_xzzz[k] * tke_0 - 4.0 * to_xzzz_xz[k] * tbe_0 + 4.0 * to_xzzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_yyy[k] = -4.0 * to_xz_yyyz[k] * tke_0 + 4.0 * to_xzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_xzz_yyz[k] = 2.0 * to_xz_yy[k] - 4.0 * to_xz_yyzz[k] * tke_0 - 2.0 * to_xzzz_yy[k] * tbe_0 + 4.0 * to_xzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_yzz[k] = 4.0 * to_xz_yz[k] - 4.0 * to_xz_yzzz[k] * tke_0 - 4.0 * to_xzzz_yz[k] * tbe_0 + 4.0 * to_xzzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_xzz_zzz[k] = 6.0 * to_xz_zz[k] - 4.0 * to_xz_zzzz[k] * tke_0 - 6.0 * to_xzzz_zz[k] * tbe_0 + 4.0 * to_xzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 860-870 components of targeted buffer : FF

        auto to_z_z_yyy_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 60);

        auto to_z_z_yyy_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 61);

        auto to_z_z_yyy_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 62);

        auto to_z_z_yyy_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 63);

        auto to_z_z_yyy_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 64);

        auto to_z_z_yyy_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 65);

        auto to_z_z_yyy_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 66);

        auto to_z_z_yyy_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 67);

        auto to_z_z_yyy_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 68);

        auto to_z_z_yyy_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 69);

        #pragma omp simd aligned(to_yyyz_xx, to_yyyz_xxxz, to_yyyz_xxyz, to_yyyz_xxzz, to_yyyz_xy, to_yyyz_xyyz, to_yyyz_xyzz, to_yyyz_xz, to_yyyz_xzzz, to_yyyz_yy, to_yyyz_yyyz, to_yyyz_yyzz, to_yyyz_yz, to_yyyz_yzzz, to_yyyz_zz, to_yyyz_zzzz, to_z_z_yyy_xxx, to_z_z_yyy_xxy, to_z_z_yyy_xxz, to_z_z_yyy_xyy, to_z_z_yyy_xyz, to_z_z_yyy_xzz, to_z_z_yyy_yyy, to_z_z_yyy_yyz, to_z_z_yyy_yzz, to_z_z_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyy_xxx[k] = 4.0 * to_yyyz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xxy[k] = 4.0 * to_yyyz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xxz[k] = -2.0 * to_yyyz_xx[k] * tbe_0 + 4.0 * to_yyyz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xyy[k] = 4.0 * to_yyyz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xyz[k] = -2.0 * to_yyyz_xy[k] * tbe_0 + 4.0 * to_yyyz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_xzz[k] = -4.0 * to_yyyz_xz[k] * tbe_0 + 4.0 * to_yyyz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_yyy[k] = 4.0 * to_yyyz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_yyy_yyz[k] = -2.0 * to_yyyz_yy[k] * tbe_0 + 4.0 * to_yyyz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_yzz[k] = -4.0 * to_yyyz_yz[k] * tbe_0 + 4.0 * to_yyyz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_yyy_zzz[k] = -6.0 * to_yyyz_zz[k] * tbe_0 + 4.0 * to_yyyz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 870-880 components of targeted buffer : FF

        auto to_z_z_yyz_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 70);

        auto to_z_z_yyz_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 71);

        auto to_z_z_yyz_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 72);

        auto to_z_z_yyz_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 73);

        auto to_z_z_yyz_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 74);

        auto to_z_z_yyz_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 75);

        auto to_z_z_yyz_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 76);

        auto to_z_z_yyz_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 77);

        auto to_z_z_yyz_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 78);

        auto to_z_z_yyz_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 79);

        #pragma omp simd aligned(to_yy_xx, to_yy_xxxz, to_yy_xxyz, to_yy_xxzz, to_yy_xy, to_yy_xyyz, to_yy_xyzz, to_yy_xz, to_yy_xzzz, to_yy_yy, to_yy_yyyz, to_yy_yyzz, to_yy_yz, to_yy_yzzz, to_yy_zz, to_yy_zzzz, to_yyzz_xx, to_yyzz_xxxz, to_yyzz_xxyz, to_yyzz_xxzz, to_yyzz_xy, to_yyzz_xyyz, to_yyzz_xyzz, to_yyzz_xz, to_yyzz_xzzz, to_yyzz_yy, to_yyzz_yyyz, to_yyzz_yyzz, to_yyzz_yz, to_yyzz_yzzz, to_yyzz_zz, to_yyzz_zzzz, to_z_z_yyz_xxx, to_z_z_yyz_xxy, to_z_z_yyz_xxz, to_z_z_yyz_xyy, to_z_z_yyz_xyz, to_z_z_yyz_xzz, to_z_z_yyz_yyy, to_z_z_yyz_yyz, to_z_z_yyz_yzz, to_z_z_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyz_xxx[k] = -2.0 * to_yy_xxxz[k] * tke_0 + 4.0 * to_yyzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xxy[k] = -2.0 * to_yy_xxyz[k] * tke_0 + 4.0 * to_yyzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xxz[k] = to_yy_xx[k] - 2.0 * to_yy_xxzz[k] * tke_0 - 2.0 * to_yyzz_xx[k] * tbe_0 + 4.0 * to_yyzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xyy[k] = -2.0 * to_yy_xyyz[k] * tke_0 + 4.0 * to_yyzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xyz[k] = to_yy_xy[k] - 2.0 * to_yy_xyzz[k] * tke_0 - 2.0 * to_yyzz_xy[k] * tbe_0 + 4.0 * to_yyzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_xzz[k] = 2.0 * to_yy_xz[k] - 2.0 * to_yy_xzzz[k] * tke_0 - 4.0 * to_yyzz_xz[k] * tbe_0 + 4.0 * to_yyzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_yyy[k] = -2.0 * to_yy_yyyz[k] * tke_0 + 4.0 * to_yyzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_yyz_yyz[k] = to_yy_yy[k] - 2.0 * to_yy_yyzz[k] * tke_0 - 2.0 * to_yyzz_yy[k] * tbe_0 + 4.0 * to_yyzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_yzz[k] = 2.0 * to_yy_yz[k] - 2.0 * to_yy_yzzz[k] * tke_0 - 4.0 * to_yyzz_yz[k] * tbe_0 + 4.0 * to_yyzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_yyz_zzz[k] = 3.0 * to_yy_zz[k] - 2.0 * to_yy_zzzz[k] * tke_0 - 6.0 * to_yyzz_zz[k] * tbe_0 + 4.0 * to_yyzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 880-890 components of targeted buffer : FF

        auto to_z_z_yzz_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 80);

        auto to_z_z_yzz_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 81);

        auto to_z_z_yzz_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 82);

        auto to_z_z_yzz_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 83);

        auto to_z_z_yzz_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 84);

        auto to_z_z_yzz_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 85);

        auto to_z_z_yzz_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 86);

        auto to_z_z_yzz_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 87);

        auto to_z_z_yzz_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 88);

        auto to_z_z_yzz_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 89);

        #pragma omp simd aligned(to_yz_xx, to_yz_xxxz, to_yz_xxyz, to_yz_xxzz, to_yz_xy, to_yz_xyyz, to_yz_xyzz, to_yz_xz, to_yz_xzzz, to_yz_yy, to_yz_yyyz, to_yz_yyzz, to_yz_yz, to_yz_yzzz, to_yz_zz, to_yz_zzzz, to_yzzz_xx, to_yzzz_xxxz, to_yzzz_xxyz, to_yzzz_xxzz, to_yzzz_xy, to_yzzz_xyyz, to_yzzz_xyzz, to_yzzz_xz, to_yzzz_xzzz, to_yzzz_yy, to_yzzz_yyyz, to_yzzz_yyzz, to_yzzz_yz, to_yzzz_yzzz, to_yzzz_zz, to_yzzz_zzzz, to_z_z_yzz_xxx, to_z_z_yzz_xxy, to_z_z_yzz_xxz, to_z_z_yzz_xyy, to_z_z_yzz_xyz, to_z_z_yzz_xzz, to_z_z_yzz_yyy, to_z_z_yzz_yyz, to_z_z_yzz_yzz, to_z_z_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yzz_xxx[k] = -4.0 * to_yz_xxxz[k] * tke_0 + 4.0 * to_yzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xxy[k] = -4.0 * to_yz_xxyz[k] * tke_0 + 4.0 * to_yzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xxz[k] = 2.0 * to_yz_xx[k] - 4.0 * to_yz_xxzz[k] * tke_0 - 2.0 * to_yzzz_xx[k] * tbe_0 + 4.0 * to_yzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xyy[k] = -4.0 * to_yz_xyyz[k] * tke_0 + 4.0 * to_yzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xyz[k] = 2.0 * to_yz_xy[k] - 4.0 * to_yz_xyzz[k] * tke_0 - 2.0 * to_yzzz_xy[k] * tbe_0 + 4.0 * to_yzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_xzz[k] = 4.0 * to_yz_xz[k] - 4.0 * to_yz_xzzz[k] * tke_0 - 4.0 * to_yzzz_xz[k] * tbe_0 + 4.0 * to_yzzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_yyy[k] = -4.0 * to_yz_yyyz[k] * tke_0 + 4.0 * to_yzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_yzz_yyz[k] = 2.0 * to_yz_yy[k] - 4.0 * to_yz_yyzz[k] * tke_0 - 2.0 * to_yzzz_yy[k] * tbe_0 + 4.0 * to_yzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_yzz[k] = 4.0 * to_yz_yz[k] - 4.0 * to_yz_yzzz[k] * tke_0 - 4.0 * to_yzzz_yz[k] * tbe_0 + 4.0 * to_yzzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_yzz_zzz[k] = 6.0 * to_yz_zz[k] - 4.0 * to_yz_zzzz[k] * tke_0 - 6.0 * to_yzzz_zz[k] * tbe_0 + 4.0 * to_yzzz_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 890-900 components of targeted buffer : FF

        auto to_z_z_zzz_xxx = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 90);

        auto to_z_z_zzz_xxy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 91);

        auto to_z_z_zzz_xxz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 92);

        auto to_z_z_zzz_xyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 93);

        auto to_z_z_zzz_xyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 94);

        auto to_z_z_zzz_xzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 95);

        auto to_z_z_zzz_yyy = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 96);

        auto to_z_z_zzz_yyz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 97);

        auto to_z_z_zzz_yzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 98);

        auto to_z_z_zzz_zzz = pbuffer.data(idx_op_geom_101_ff + 8 * op_comps * 100 + i * 100 + 99);

        #pragma omp simd aligned(to_z_z_zzz_xxx, to_z_z_zzz_xxy, to_z_z_zzz_xxz, to_z_z_zzz_xyy, to_z_z_zzz_xyz, to_z_z_zzz_xzz, to_z_z_zzz_yyy, to_z_z_zzz_yyz, to_z_z_zzz_yzz, to_z_z_zzz_zzz, to_zz_xx, to_zz_xxxz, to_zz_xxyz, to_zz_xxzz, to_zz_xy, to_zz_xyyz, to_zz_xyzz, to_zz_xz, to_zz_xzzz, to_zz_yy, to_zz_yyyz, to_zz_yyzz, to_zz_yz, to_zz_yzzz, to_zz_zz, to_zz_zzzz, to_zzzz_xx, to_zzzz_xxxz, to_zzzz_xxyz, to_zzzz_xxzz, to_zzzz_xy, to_zzzz_xyyz, to_zzzz_xyzz, to_zzzz_xz, to_zzzz_xzzz, to_zzzz_yy, to_zzzz_yyyz, to_zzzz_yyzz, to_zzzz_yz, to_zzzz_yzzz, to_zzzz_zz, to_zzzz_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zzz_xxx[k] = -6.0 * to_zz_xxxz[k] * tke_0 + 4.0 * to_zzzz_xxxz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xxy[k] = -6.0 * to_zz_xxyz[k] * tke_0 + 4.0 * to_zzzz_xxyz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xxz[k] = 3.0 * to_zz_xx[k] - 6.0 * to_zz_xxzz[k] * tke_0 - 2.0 * to_zzzz_xx[k] * tbe_0 + 4.0 * to_zzzz_xxzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xyy[k] = -6.0 * to_zz_xyyz[k] * tke_0 + 4.0 * to_zzzz_xyyz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xyz[k] = 3.0 * to_zz_xy[k] - 6.0 * to_zz_xyzz[k] * tke_0 - 2.0 * to_zzzz_xy[k] * tbe_0 + 4.0 * to_zzzz_xyzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_xzz[k] = 6.0 * to_zz_xz[k] - 6.0 * to_zz_xzzz[k] * tke_0 - 4.0 * to_zzzz_xz[k] * tbe_0 + 4.0 * to_zzzz_xzzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_yyy[k] = -6.0 * to_zz_yyyz[k] * tke_0 + 4.0 * to_zzzz_yyyz[k] * tbe_0 * tke_0;

            to_z_z_zzz_yyz[k] = 3.0 * to_zz_yy[k] - 6.0 * to_zz_yyzz[k] * tke_0 - 2.0 * to_zzzz_yy[k] * tbe_0 + 4.0 * to_zzzz_yyzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_yzz[k] = 6.0 * to_zz_yz[k] - 6.0 * to_zz_yzzz[k] * tke_0 - 4.0 * to_zzzz_yz[k] * tbe_0 + 4.0 * to_zzzz_yzzz[k] * tbe_0 * tke_0;

            to_z_z_zzz_zzz[k] = 9.0 * to_zz_zz[k] - 6.0 * to_zz_zzzz[k] * tke_0 - 6.0 * to_zzzz_zz[k] * tbe_0 + 4.0 * to_zzzz_zzzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

