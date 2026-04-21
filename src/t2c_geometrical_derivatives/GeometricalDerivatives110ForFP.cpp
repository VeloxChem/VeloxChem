#include "GeometricalDerivatives110ForFP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_fp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_fp,
                         const int idx_op_pp,
                         const int idx_op_ds,
                         const int idx_op_dd,
                         const int idx_op_fp,
                         const int idx_op_gs,
                         const int idx_op_gd,
                         const int idx_op_hp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : PP

    auto tr_x_x = pbuffer.data(idx_op_pp);

    auto tr_x_y = pbuffer.data(idx_op_pp + 1);

    auto tr_x_z = pbuffer.data(idx_op_pp + 2);

    auto tr_y_x = pbuffer.data(idx_op_pp + 3);

    auto tr_y_y = pbuffer.data(idx_op_pp + 4);

    auto tr_y_z = pbuffer.data(idx_op_pp + 5);

    auto tr_z_x = pbuffer.data(idx_op_pp + 6);

    auto tr_z_y = pbuffer.data(idx_op_pp + 7);

    auto tr_z_z = pbuffer.data(idx_op_pp + 8);

    // Set up components of auxiliary buffer : DS

    auto tr_xx_0 = pbuffer.data(idx_op_ds);

    auto tr_xy_0 = pbuffer.data(idx_op_ds + 1);

    auto tr_xz_0 = pbuffer.data(idx_op_ds + 2);

    auto tr_yy_0 = pbuffer.data(idx_op_ds + 3);

    auto tr_yz_0 = pbuffer.data(idx_op_ds + 4);

    auto tr_zz_0 = pbuffer.data(idx_op_ds + 5);

    // Set up components of auxiliary buffer : DD

    auto tr_xx_xx = pbuffer.data(idx_op_dd);

    auto tr_xx_xy = pbuffer.data(idx_op_dd + 1);

    auto tr_xx_xz = pbuffer.data(idx_op_dd + 2);

    auto tr_xx_yy = pbuffer.data(idx_op_dd + 3);

    auto tr_xx_yz = pbuffer.data(idx_op_dd + 4);

    auto tr_xx_zz = pbuffer.data(idx_op_dd + 5);

    auto tr_xy_xx = pbuffer.data(idx_op_dd + 6);

    auto tr_xy_xy = pbuffer.data(idx_op_dd + 7);

    auto tr_xy_xz = pbuffer.data(idx_op_dd + 8);

    auto tr_xy_yy = pbuffer.data(idx_op_dd + 9);

    auto tr_xy_yz = pbuffer.data(idx_op_dd + 10);

    auto tr_xy_zz = pbuffer.data(idx_op_dd + 11);

    auto tr_xz_xx = pbuffer.data(idx_op_dd + 12);

    auto tr_xz_xy = pbuffer.data(idx_op_dd + 13);

    auto tr_xz_xz = pbuffer.data(idx_op_dd + 14);

    auto tr_xz_yy = pbuffer.data(idx_op_dd + 15);

    auto tr_xz_yz = pbuffer.data(idx_op_dd + 16);

    auto tr_xz_zz = pbuffer.data(idx_op_dd + 17);

    auto tr_yy_xx = pbuffer.data(idx_op_dd + 18);

    auto tr_yy_xy = pbuffer.data(idx_op_dd + 19);

    auto tr_yy_xz = pbuffer.data(idx_op_dd + 20);

    auto tr_yy_yy = pbuffer.data(idx_op_dd + 21);

    auto tr_yy_yz = pbuffer.data(idx_op_dd + 22);

    auto tr_yy_zz = pbuffer.data(idx_op_dd + 23);

    auto tr_yz_xx = pbuffer.data(idx_op_dd + 24);

    auto tr_yz_xy = pbuffer.data(idx_op_dd + 25);

    auto tr_yz_xz = pbuffer.data(idx_op_dd + 26);

    auto tr_yz_yy = pbuffer.data(idx_op_dd + 27);

    auto tr_yz_yz = pbuffer.data(idx_op_dd + 28);

    auto tr_yz_zz = pbuffer.data(idx_op_dd + 29);

    auto tr_zz_xx = pbuffer.data(idx_op_dd + 30);

    auto tr_zz_xy = pbuffer.data(idx_op_dd + 31);

    auto tr_zz_xz = pbuffer.data(idx_op_dd + 32);

    auto tr_zz_yy = pbuffer.data(idx_op_dd + 33);

    auto tr_zz_yz = pbuffer.data(idx_op_dd + 34);

    auto tr_zz_zz = pbuffer.data(idx_op_dd + 35);

    // Set up components of auxiliary buffer : FP

    auto tr_xxx_x = pbuffer.data(idx_op_fp);

    auto tr_xxx_y = pbuffer.data(idx_op_fp + 1);

    auto tr_xxx_z = pbuffer.data(idx_op_fp + 2);

    auto tr_xxy_x = pbuffer.data(idx_op_fp + 3);

    auto tr_xxy_y = pbuffer.data(idx_op_fp + 4);

    auto tr_xxy_z = pbuffer.data(idx_op_fp + 5);

    auto tr_xxz_x = pbuffer.data(idx_op_fp + 6);

    auto tr_xxz_y = pbuffer.data(idx_op_fp + 7);

    auto tr_xxz_z = pbuffer.data(idx_op_fp + 8);

    auto tr_xyy_x = pbuffer.data(idx_op_fp + 9);

    auto tr_xyy_y = pbuffer.data(idx_op_fp + 10);

    auto tr_xyy_z = pbuffer.data(idx_op_fp + 11);

    auto tr_xyz_x = pbuffer.data(idx_op_fp + 12);

    auto tr_xyz_y = pbuffer.data(idx_op_fp + 13);

    auto tr_xyz_z = pbuffer.data(idx_op_fp + 14);

    auto tr_xzz_x = pbuffer.data(idx_op_fp + 15);

    auto tr_xzz_y = pbuffer.data(idx_op_fp + 16);

    auto tr_xzz_z = pbuffer.data(idx_op_fp + 17);

    auto tr_yyy_x = pbuffer.data(idx_op_fp + 18);

    auto tr_yyy_y = pbuffer.data(idx_op_fp + 19);

    auto tr_yyy_z = pbuffer.data(idx_op_fp + 20);

    auto tr_yyz_x = pbuffer.data(idx_op_fp + 21);

    auto tr_yyz_y = pbuffer.data(idx_op_fp + 22);

    auto tr_yyz_z = pbuffer.data(idx_op_fp + 23);

    auto tr_yzz_x = pbuffer.data(idx_op_fp + 24);

    auto tr_yzz_y = pbuffer.data(idx_op_fp + 25);

    auto tr_yzz_z = pbuffer.data(idx_op_fp + 26);

    auto tr_zzz_x = pbuffer.data(idx_op_fp + 27);

    auto tr_zzz_y = pbuffer.data(idx_op_fp + 28);

    auto tr_zzz_z = pbuffer.data(idx_op_fp + 29);

    // Set up components of auxiliary buffer : GS

    auto tr_xxxx_0 = pbuffer.data(idx_op_gs);

    auto tr_xxxy_0 = pbuffer.data(idx_op_gs + 1);

    auto tr_xxxz_0 = pbuffer.data(idx_op_gs + 2);

    auto tr_xxyy_0 = pbuffer.data(idx_op_gs + 3);

    auto tr_xxyz_0 = pbuffer.data(idx_op_gs + 4);

    auto tr_xxzz_0 = pbuffer.data(idx_op_gs + 5);

    auto tr_xyyy_0 = pbuffer.data(idx_op_gs + 6);

    auto tr_xyyz_0 = pbuffer.data(idx_op_gs + 7);

    auto tr_xyzz_0 = pbuffer.data(idx_op_gs + 8);

    auto tr_xzzz_0 = pbuffer.data(idx_op_gs + 9);

    auto tr_yyyy_0 = pbuffer.data(idx_op_gs + 10);

    auto tr_yyyz_0 = pbuffer.data(idx_op_gs + 11);

    auto tr_yyzz_0 = pbuffer.data(idx_op_gs + 12);

    auto tr_yzzz_0 = pbuffer.data(idx_op_gs + 13);

    auto tr_zzzz_0 = pbuffer.data(idx_op_gs + 14);

    // Set up components of auxiliary buffer : GD

    auto tr_xxxx_xx = pbuffer.data(idx_op_gd);

    auto tr_xxxx_xy = pbuffer.data(idx_op_gd + 1);

    auto tr_xxxx_xz = pbuffer.data(idx_op_gd + 2);

    auto tr_xxxx_yy = pbuffer.data(idx_op_gd + 3);

    auto tr_xxxx_yz = pbuffer.data(idx_op_gd + 4);

    auto tr_xxxx_zz = pbuffer.data(idx_op_gd + 5);

    auto tr_xxxy_xx = pbuffer.data(idx_op_gd + 6);

    auto tr_xxxy_xy = pbuffer.data(idx_op_gd + 7);

    auto tr_xxxy_xz = pbuffer.data(idx_op_gd + 8);

    auto tr_xxxy_yy = pbuffer.data(idx_op_gd + 9);

    auto tr_xxxy_yz = pbuffer.data(idx_op_gd + 10);

    auto tr_xxxy_zz = pbuffer.data(idx_op_gd + 11);

    auto tr_xxxz_xx = pbuffer.data(idx_op_gd + 12);

    auto tr_xxxz_xy = pbuffer.data(idx_op_gd + 13);

    auto tr_xxxz_xz = pbuffer.data(idx_op_gd + 14);

    auto tr_xxxz_yy = pbuffer.data(idx_op_gd + 15);

    auto tr_xxxz_yz = pbuffer.data(idx_op_gd + 16);

    auto tr_xxxz_zz = pbuffer.data(idx_op_gd + 17);

    auto tr_xxyy_xx = pbuffer.data(idx_op_gd + 18);

    auto tr_xxyy_xy = pbuffer.data(idx_op_gd + 19);

    auto tr_xxyy_xz = pbuffer.data(idx_op_gd + 20);

    auto tr_xxyy_yy = pbuffer.data(idx_op_gd + 21);

    auto tr_xxyy_yz = pbuffer.data(idx_op_gd + 22);

    auto tr_xxyy_zz = pbuffer.data(idx_op_gd + 23);

    auto tr_xxyz_xx = pbuffer.data(idx_op_gd + 24);

    auto tr_xxyz_xy = pbuffer.data(idx_op_gd + 25);

    auto tr_xxyz_xz = pbuffer.data(idx_op_gd + 26);

    auto tr_xxyz_yy = pbuffer.data(idx_op_gd + 27);

    auto tr_xxyz_yz = pbuffer.data(idx_op_gd + 28);

    auto tr_xxyz_zz = pbuffer.data(idx_op_gd + 29);

    auto tr_xxzz_xx = pbuffer.data(idx_op_gd + 30);

    auto tr_xxzz_xy = pbuffer.data(idx_op_gd + 31);

    auto tr_xxzz_xz = pbuffer.data(idx_op_gd + 32);

    auto tr_xxzz_yy = pbuffer.data(idx_op_gd + 33);

    auto tr_xxzz_yz = pbuffer.data(idx_op_gd + 34);

    auto tr_xxzz_zz = pbuffer.data(idx_op_gd + 35);

    auto tr_xyyy_xx = pbuffer.data(idx_op_gd + 36);

    auto tr_xyyy_xy = pbuffer.data(idx_op_gd + 37);

    auto tr_xyyy_xz = pbuffer.data(idx_op_gd + 38);

    auto tr_xyyy_yy = pbuffer.data(idx_op_gd + 39);

    auto tr_xyyy_yz = pbuffer.data(idx_op_gd + 40);

    auto tr_xyyy_zz = pbuffer.data(idx_op_gd + 41);

    auto tr_xyyz_xx = pbuffer.data(idx_op_gd + 42);

    auto tr_xyyz_xy = pbuffer.data(idx_op_gd + 43);

    auto tr_xyyz_xz = pbuffer.data(idx_op_gd + 44);

    auto tr_xyyz_yy = pbuffer.data(idx_op_gd + 45);

    auto tr_xyyz_yz = pbuffer.data(idx_op_gd + 46);

    auto tr_xyyz_zz = pbuffer.data(idx_op_gd + 47);

    auto tr_xyzz_xx = pbuffer.data(idx_op_gd + 48);

    auto tr_xyzz_xy = pbuffer.data(idx_op_gd + 49);

    auto tr_xyzz_xz = pbuffer.data(idx_op_gd + 50);

    auto tr_xyzz_yy = pbuffer.data(idx_op_gd + 51);

    auto tr_xyzz_yz = pbuffer.data(idx_op_gd + 52);

    auto tr_xyzz_zz = pbuffer.data(idx_op_gd + 53);

    auto tr_xzzz_xx = pbuffer.data(idx_op_gd + 54);

    auto tr_xzzz_xy = pbuffer.data(idx_op_gd + 55);

    auto tr_xzzz_xz = pbuffer.data(idx_op_gd + 56);

    auto tr_xzzz_yy = pbuffer.data(idx_op_gd + 57);

    auto tr_xzzz_yz = pbuffer.data(idx_op_gd + 58);

    auto tr_xzzz_zz = pbuffer.data(idx_op_gd + 59);

    auto tr_yyyy_xx = pbuffer.data(idx_op_gd + 60);

    auto tr_yyyy_xy = pbuffer.data(idx_op_gd + 61);

    auto tr_yyyy_xz = pbuffer.data(idx_op_gd + 62);

    auto tr_yyyy_yy = pbuffer.data(idx_op_gd + 63);

    auto tr_yyyy_yz = pbuffer.data(idx_op_gd + 64);

    auto tr_yyyy_zz = pbuffer.data(idx_op_gd + 65);

    auto tr_yyyz_xx = pbuffer.data(idx_op_gd + 66);

    auto tr_yyyz_xy = pbuffer.data(idx_op_gd + 67);

    auto tr_yyyz_xz = pbuffer.data(idx_op_gd + 68);

    auto tr_yyyz_yy = pbuffer.data(idx_op_gd + 69);

    auto tr_yyyz_yz = pbuffer.data(idx_op_gd + 70);

    auto tr_yyyz_zz = pbuffer.data(idx_op_gd + 71);

    auto tr_yyzz_xx = pbuffer.data(idx_op_gd + 72);

    auto tr_yyzz_xy = pbuffer.data(idx_op_gd + 73);

    auto tr_yyzz_xz = pbuffer.data(idx_op_gd + 74);

    auto tr_yyzz_yy = pbuffer.data(idx_op_gd + 75);

    auto tr_yyzz_yz = pbuffer.data(idx_op_gd + 76);

    auto tr_yyzz_zz = pbuffer.data(idx_op_gd + 77);

    auto tr_yzzz_xx = pbuffer.data(idx_op_gd + 78);

    auto tr_yzzz_xy = pbuffer.data(idx_op_gd + 79);

    auto tr_yzzz_xz = pbuffer.data(idx_op_gd + 80);

    auto tr_yzzz_yy = pbuffer.data(idx_op_gd + 81);

    auto tr_yzzz_yz = pbuffer.data(idx_op_gd + 82);

    auto tr_yzzz_zz = pbuffer.data(idx_op_gd + 83);

    auto tr_zzzz_xx = pbuffer.data(idx_op_gd + 84);

    auto tr_zzzz_xy = pbuffer.data(idx_op_gd + 85);

    auto tr_zzzz_xz = pbuffer.data(idx_op_gd + 86);

    auto tr_zzzz_yy = pbuffer.data(idx_op_gd + 87);

    auto tr_zzzz_yz = pbuffer.data(idx_op_gd + 88);

    auto tr_zzzz_zz = pbuffer.data(idx_op_gd + 89);

    // Set up components of auxiliary buffer : HP

    auto tr_xxxxx_x = pbuffer.data(idx_op_hp);

    auto tr_xxxxx_y = pbuffer.data(idx_op_hp + 1);

    auto tr_xxxxx_z = pbuffer.data(idx_op_hp + 2);

    auto tr_xxxxy_x = pbuffer.data(idx_op_hp + 3);

    auto tr_xxxxy_y = pbuffer.data(idx_op_hp + 4);

    auto tr_xxxxy_z = pbuffer.data(idx_op_hp + 5);

    auto tr_xxxxz_x = pbuffer.data(idx_op_hp + 6);

    auto tr_xxxxz_y = pbuffer.data(idx_op_hp + 7);

    auto tr_xxxxz_z = pbuffer.data(idx_op_hp + 8);

    auto tr_xxxyy_x = pbuffer.data(idx_op_hp + 9);

    auto tr_xxxyy_y = pbuffer.data(idx_op_hp + 10);

    auto tr_xxxyy_z = pbuffer.data(idx_op_hp + 11);

    auto tr_xxxyz_x = pbuffer.data(idx_op_hp + 12);

    auto tr_xxxyz_y = pbuffer.data(idx_op_hp + 13);

    auto tr_xxxyz_z = pbuffer.data(idx_op_hp + 14);

    auto tr_xxxzz_x = pbuffer.data(idx_op_hp + 15);

    auto tr_xxxzz_y = pbuffer.data(idx_op_hp + 16);

    auto tr_xxxzz_z = pbuffer.data(idx_op_hp + 17);

    auto tr_xxyyy_x = pbuffer.data(idx_op_hp + 18);

    auto tr_xxyyy_y = pbuffer.data(idx_op_hp + 19);

    auto tr_xxyyy_z = pbuffer.data(idx_op_hp + 20);

    auto tr_xxyyz_x = pbuffer.data(idx_op_hp + 21);

    auto tr_xxyyz_y = pbuffer.data(idx_op_hp + 22);

    auto tr_xxyyz_z = pbuffer.data(idx_op_hp + 23);

    auto tr_xxyzz_x = pbuffer.data(idx_op_hp + 24);

    auto tr_xxyzz_y = pbuffer.data(idx_op_hp + 25);

    auto tr_xxyzz_z = pbuffer.data(idx_op_hp + 26);

    auto tr_xxzzz_x = pbuffer.data(idx_op_hp + 27);

    auto tr_xxzzz_y = pbuffer.data(idx_op_hp + 28);

    auto tr_xxzzz_z = pbuffer.data(idx_op_hp + 29);

    auto tr_xyyyy_x = pbuffer.data(idx_op_hp + 30);

    auto tr_xyyyy_y = pbuffer.data(idx_op_hp + 31);

    auto tr_xyyyy_z = pbuffer.data(idx_op_hp + 32);

    auto tr_xyyyz_x = pbuffer.data(idx_op_hp + 33);

    auto tr_xyyyz_y = pbuffer.data(idx_op_hp + 34);

    auto tr_xyyyz_z = pbuffer.data(idx_op_hp + 35);

    auto tr_xyyzz_x = pbuffer.data(idx_op_hp + 36);

    auto tr_xyyzz_y = pbuffer.data(idx_op_hp + 37);

    auto tr_xyyzz_z = pbuffer.data(idx_op_hp + 38);

    auto tr_xyzzz_x = pbuffer.data(idx_op_hp + 39);

    auto tr_xyzzz_y = pbuffer.data(idx_op_hp + 40);

    auto tr_xyzzz_z = pbuffer.data(idx_op_hp + 41);

    auto tr_xzzzz_x = pbuffer.data(idx_op_hp + 42);

    auto tr_xzzzz_y = pbuffer.data(idx_op_hp + 43);

    auto tr_xzzzz_z = pbuffer.data(idx_op_hp + 44);

    auto tr_yyyyy_x = pbuffer.data(idx_op_hp + 45);

    auto tr_yyyyy_y = pbuffer.data(idx_op_hp + 46);

    auto tr_yyyyy_z = pbuffer.data(idx_op_hp + 47);

    auto tr_yyyyz_x = pbuffer.data(idx_op_hp + 48);

    auto tr_yyyyz_y = pbuffer.data(idx_op_hp + 49);

    auto tr_yyyyz_z = pbuffer.data(idx_op_hp + 50);

    auto tr_yyyzz_x = pbuffer.data(idx_op_hp + 51);

    auto tr_yyyzz_y = pbuffer.data(idx_op_hp + 52);

    auto tr_yyyzz_z = pbuffer.data(idx_op_hp + 53);

    auto tr_yyzzz_x = pbuffer.data(idx_op_hp + 54);

    auto tr_yyzzz_y = pbuffer.data(idx_op_hp + 55);

    auto tr_yyzzz_z = pbuffer.data(idx_op_hp + 56);

    auto tr_yzzzz_x = pbuffer.data(idx_op_hp + 57);

    auto tr_yzzzz_y = pbuffer.data(idx_op_hp + 58);

    auto tr_yzzzz_z = pbuffer.data(idx_op_hp + 59);

    auto tr_zzzzz_x = pbuffer.data(idx_op_hp + 60);

    auto tr_zzzzz_y = pbuffer.data(idx_op_hp + 61);

    auto tr_zzzzz_z = pbuffer.data(idx_op_hp + 62);

    // Set up 0-3 components of targeted buffer : FP

    auto tr_x_0_x_xxx_x = pbuffer.data(idx_op_geom_110_fp);

    auto tr_x_0_x_xxx_y = pbuffer.data(idx_op_geom_110_fp + 1);

    auto tr_x_0_x_xxx_z = pbuffer.data(idx_op_geom_110_fp + 2);

    #pragma omp simd aligned(tr_x_0_x_xxx_x, tr_x_0_x_xxx_y, tr_x_0_x_xxx_z, tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxx_0, tr_xxxx_xx, tr_xxxx_xy, tr_xxxx_xz, tr_xxxxx_x, tr_xxxxx_y, tr_xxxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxx_x[i] = 6.0 * tr_x_x[i] + 3.0 * tr_xx_0[i] - 6.0 * tr_xx_xx[i] * tke_0 - 14.0 * tr_xxx_x[i] * tbe_0 - 2.0 * tr_xxxx_0[i] * tbe_0 + 4.0 * tr_xxxx_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_y[i] = 6.0 * tr_x_y[i] - 6.0 * tr_xx_xy[i] * tke_0 - 14.0 * tr_xxx_y[i] * tbe_0 + 4.0 * tr_xxxx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_z[i] = 6.0 * tr_x_z[i] - 6.0 * tr_xx_xz[i] * tke_0 - 14.0 * tr_xxx_z[i] * tbe_0 + 4.0 * tr_xxxx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_z[i] * tbe_0 * tbe_0;
    }

    // Set up 3-6 components of targeted buffer : FP

    auto tr_x_0_x_xxy_x = pbuffer.data(idx_op_geom_110_fp + 3);

    auto tr_x_0_x_xxy_y = pbuffer.data(idx_op_geom_110_fp + 4);

    auto tr_x_0_x_xxy_z = pbuffer.data(idx_op_geom_110_fp + 5);

    #pragma omp simd aligned(tr_x_0_x_xxy_x, tr_x_0_x_xxy_y, tr_x_0_x_xxy_z, tr_xxxxy_x, tr_xxxxy_y, tr_xxxxy_z, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xy_0, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_y_x, tr_y_y, tr_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxy_x[i] = 2.0 * tr_y_x[i] + 2.0 * tr_xy_0[i] - 4.0 * tr_xy_xx[i] * tke_0 - 10.0 * tr_xxy_x[i] * tbe_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_y[i] = 2.0 * tr_y_y[i] - 4.0 * tr_xy_xy[i] * tke_0 - 10.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_z[i] = 2.0 * tr_y_z[i] - 4.0 * tr_xy_xz[i] * tke_0 - 10.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 6-9 components of targeted buffer : FP

    auto tr_x_0_x_xxz_x = pbuffer.data(idx_op_geom_110_fp + 6);

    auto tr_x_0_x_xxz_y = pbuffer.data(idx_op_geom_110_fp + 7);

    auto tr_x_0_x_xxz_z = pbuffer.data(idx_op_geom_110_fp + 8);

    #pragma omp simd aligned(tr_x_0_x_xxz_x, tr_x_0_x_xxz_y, tr_x_0_x_xxz_z, tr_xxxxz_x, tr_xxxxz_y, tr_xxxxz_z, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xz_0, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxz_x[i] = 2.0 * tr_z_x[i] + 2.0 * tr_xz_0[i] - 4.0 * tr_xz_xx[i] * tke_0 - 10.0 * tr_xxz_x[i] * tbe_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_y[i] = 2.0 * tr_z_y[i] - 4.0 * tr_xz_xy[i] * tke_0 - 10.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_z[i] = 2.0 * tr_z_z[i] - 4.0 * tr_xz_xz[i] * tke_0 - 10.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 9-12 components of targeted buffer : FP

    auto tr_x_0_x_xyy_x = pbuffer.data(idx_op_geom_110_fp + 9);

    auto tr_x_0_x_xyy_y = pbuffer.data(idx_op_geom_110_fp + 10);

    auto tr_x_0_x_xyy_z = pbuffer.data(idx_op_geom_110_fp + 11);

    #pragma omp simd aligned(tr_x_0_x_xyy_x, tr_x_0_x_xyy_y, tr_x_0_x_xyy_z, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_yy_0, tr_yy_xx, tr_yy_xy, tr_yy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyy_x[i] = tr_yy_0[i] - 2.0 * tr_yy_xx[i] * tke_0 - 6.0 * tr_xyy_x[i] * tbe_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_y[i] = -2.0 * tr_yy_xy[i] * tke_0 - 6.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xxyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_z[i] = -2.0 * tr_yy_xz[i] * tke_0 - 6.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 12-15 components of targeted buffer : FP

    auto tr_x_0_x_xyz_x = pbuffer.data(idx_op_geom_110_fp + 12);

    auto tr_x_0_x_xyz_y = pbuffer.data(idx_op_geom_110_fp + 13);

    auto tr_x_0_x_xyz_z = pbuffer.data(idx_op_geom_110_fp + 14);

    #pragma omp simd aligned(tr_x_0_x_xyz_x, tr_x_0_x_xyz_y, tr_x_0_x_xyz_z, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_yz_0, tr_yz_xx, tr_yz_xy, tr_yz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyz_x[i] = tr_yz_0[i] - 2.0 * tr_yz_xx[i] * tke_0 - 6.0 * tr_xyz_x[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_y[i] = -2.0 * tr_yz_xy[i] * tke_0 - 6.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_z[i] = -2.0 * tr_yz_xz[i] * tke_0 - 6.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 15-18 components of targeted buffer : FP

    auto tr_x_0_x_xzz_x = pbuffer.data(idx_op_geom_110_fp + 15);

    auto tr_x_0_x_xzz_y = pbuffer.data(idx_op_geom_110_fp + 16);

    auto tr_x_0_x_xzz_z = pbuffer.data(idx_op_geom_110_fp + 17);

    #pragma omp simd aligned(tr_x_0_x_xzz_x, tr_x_0_x_xzz_y, tr_x_0_x_xzz_z, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_zz_0, tr_zz_xx, tr_zz_xy, tr_zz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xzz_x[i] = tr_zz_0[i] - 2.0 * tr_zz_xx[i] * tke_0 - 6.0 * tr_xzz_x[i] * tbe_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_y[i] = -2.0 * tr_zz_xy[i] * tke_0 - 6.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_z[i] = -2.0 * tr_zz_xz[i] * tke_0 - 6.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xxzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 18-21 components of targeted buffer : FP

    auto tr_x_0_x_yyy_x = pbuffer.data(idx_op_geom_110_fp + 18);

    auto tr_x_0_x_yyy_y = pbuffer.data(idx_op_geom_110_fp + 19);

    auto tr_x_0_x_yyy_z = pbuffer.data(idx_op_geom_110_fp + 20);

    #pragma omp simd aligned(tr_x_0_x_yyy_x, tr_x_0_x_yyy_y, tr_x_0_x_yyy_z, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_yyy_x, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyy_x[i] = -2.0 * tr_yyy_x[i] * tbe_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_y[i] = -2.0 * tr_yyy_y[i] * tbe_0 + 4.0 * tr_xyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_z[i] = -2.0 * tr_yyy_z[i] * tbe_0 + 4.0 * tr_xyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 21-24 components of targeted buffer : FP

    auto tr_x_0_x_yyz_x = pbuffer.data(idx_op_geom_110_fp + 21);

    auto tr_x_0_x_yyz_y = pbuffer.data(idx_op_geom_110_fp + 22);

    auto tr_x_0_x_yyz_z = pbuffer.data(idx_op_geom_110_fp + 23);

    #pragma omp simd aligned(tr_x_0_x_yyz_x, tr_x_0_x_yyz_y, tr_x_0_x_yyz_z, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_yyz_x, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyz_x[i] = -2.0 * tr_yyz_x[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_y[i] = -2.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_z[i] = -2.0 * tr_yyz_z[i] * tbe_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 24-27 components of targeted buffer : FP

    auto tr_x_0_x_yzz_x = pbuffer.data(idx_op_geom_110_fp + 24);

    auto tr_x_0_x_yzz_y = pbuffer.data(idx_op_geom_110_fp + 25);

    auto tr_x_0_x_yzz_z = pbuffer.data(idx_op_geom_110_fp + 26);

    #pragma omp simd aligned(tr_x_0_x_yzz_x, tr_x_0_x_yzz_y, tr_x_0_x_yzz_z, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yzz_x[i] = -2.0 * tr_yzz_x[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_y[i] = -2.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_z[i] = -2.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 27-30 components of targeted buffer : FP

    auto tr_x_0_x_zzz_x = pbuffer.data(idx_op_geom_110_fp + 27);

    auto tr_x_0_x_zzz_y = pbuffer.data(idx_op_geom_110_fp + 28);

    auto tr_x_0_x_zzz_z = pbuffer.data(idx_op_geom_110_fp + 29);

    #pragma omp simd aligned(tr_x_0_x_zzz_x, tr_x_0_x_zzz_y, tr_x_0_x_zzz_z, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_zzz_x[i] = -2.0 * tr_zzz_x[i] * tbe_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_y[i] = -2.0 * tr_zzz_y[i] * tbe_0 + 4.0 * tr_xzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_z[i] = -2.0 * tr_zzz_z[i] * tbe_0 + 4.0 * tr_xzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 30-33 components of targeted buffer : FP

    auto tr_x_0_y_xxx_x = pbuffer.data(idx_op_geom_110_fp + 30);

    auto tr_x_0_y_xxx_y = pbuffer.data(idx_op_geom_110_fp + 31);

    auto tr_x_0_y_xxx_z = pbuffer.data(idx_op_geom_110_fp + 32);

    #pragma omp simd aligned(tr_x_0_y_xxx_x, tr_x_0_y_xxx_y, tr_x_0_y_xxx_z, tr_xx_0, tr_xx_xy, tr_xx_yy, tr_xx_yz, tr_xxxx_0, tr_xxxx_xy, tr_xxxx_yy, tr_xxxx_yz, tr_xxxxy_x, tr_xxxxy_y, tr_xxxxy_z, tr_xxy_x, tr_xxy_y, tr_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxx_x[i] = -6.0 * tr_xx_xy[i] * tke_0 - 6.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxxx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_y[i] = 3.0 * tr_xx_0[i] - 6.0 * tr_xx_yy[i] * tke_0 - 6.0 * tr_xxy_y[i] * tbe_0 - 2.0 * tr_xxxx_0[i] * tbe_0 + 4.0 * tr_xxxx_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_z[i] = -6.0 * tr_xx_yz[i] * tke_0 - 6.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxxx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 33-36 components of targeted buffer : FP

    auto tr_x_0_y_xxy_x = pbuffer.data(idx_op_geom_110_fp + 33);

    auto tr_x_0_y_xxy_y = pbuffer.data(idx_op_geom_110_fp + 34);

    auto tr_x_0_y_xxy_z = pbuffer.data(idx_op_geom_110_fp + 35);

    #pragma omp simd aligned(tr_x_0_y_xxy_x, tr_x_0_y_xxy_y, tr_x_0_y_xxy_z, tr_x_x, tr_x_y, tr_x_z, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxy_0, tr_xxxy_xy, tr_xxxy_yy, tr_xxxy_yz, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, tr_xy_0, tr_xy_xy, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_y, tr_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxy_x[i] = 2.0 * tr_x_x[i] - 4.0 * tr_xy_xy[i] * tke_0 - 4.0 * tr_xyy_x[i] * tbe_0 - 2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_y[i] = 2.0 * tr_x_y[i] + 2.0 * tr_xy_0[i] - 4.0 * tr_xy_yy[i] * tke_0 - 4.0 * tr_xyy_y[i] * tbe_0 - 2.0 * tr_xxx_y[i] * tbe_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_z[i] = 2.0 * tr_x_z[i] - 4.0 * tr_xy_yz[i] * tke_0 - 4.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xxx_z[i] * tbe_0 + 4.0 * tr_xxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 36-39 components of targeted buffer : FP

    auto tr_x_0_y_xxz_x = pbuffer.data(idx_op_geom_110_fp + 36);

    auto tr_x_0_y_xxz_y = pbuffer.data(idx_op_geom_110_fp + 37);

    auto tr_x_0_y_xxz_z = pbuffer.data(idx_op_geom_110_fp + 38);

    #pragma omp simd aligned(tr_x_0_y_xxz_x, tr_x_0_y_xxz_y, tr_x_0_y_xxz_z, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxxz_0, tr_xxxz_xy, tr_xxxz_yy, tr_xxxz_yz, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xz_0, tr_xz_xy, tr_xz_yy, tr_xz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxz_x[i] = -4.0 * tr_xz_xy[i] * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_y[i] = 2.0 * tr_xz_0[i] - 4.0 * tr_xz_yy[i] * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_z[i] = -4.0 * tr_xz_yz[i] * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 39-42 components of targeted buffer : FP

    auto tr_x_0_y_xyy_x = pbuffer.data(idx_op_geom_110_fp + 39);

    auto tr_x_0_y_xyy_y = pbuffer.data(idx_op_geom_110_fp + 40);

    auto tr_x_0_y_xyy_z = pbuffer.data(idx_op_geom_110_fp + 41);

    #pragma omp simd aligned(tr_x_0_y_xyy_x, tr_x_0_y_xyy_y, tr_x_0_y_xyy_z, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyy_0, tr_xxyy_xy, tr_xxyy_yy, tr_xxyy_yz, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yy_xy, tr_yy_yy, tr_yy_yz, tr_yyy_x, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyy_x[i] = 2.0 * tr_y_x[i] - 2.0 * tr_yy_xy[i] * tke_0 - 2.0 * tr_yyy_x[i] * tbe_0 - 4.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_y[i] = 2.0 * tr_y_y[i] + tr_yy_0[i] - 2.0 * tr_yy_yy[i] * tke_0 - 2.0 * tr_yyy_y[i] * tbe_0 - 4.0 * tr_xxy_y[i] * tbe_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_z[i] = 2.0 * tr_y_z[i] - 2.0 * tr_yy_yz[i] * tke_0 - 2.0 * tr_yyy_z[i] * tbe_0 - 4.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 42-45 components of targeted buffer : FP

    auto tr_x_0_y_xyz_x = pbuffer.data(idx_op_geom_110_fp + 42);

    auto tr_x_0_y_xyz_y = pbuffer.data(idx_op_geom_110_fp + 43);

    auto tr_x_0_y_xyz_z = pbuffer.data(idx_op_geom_110_fp + 44);

    #pragma omp simd aligned(tr_x_0_y_xyz_x, tr_x_0_y_xyz_y, tr_x_0_y_xyz_z, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxyz_0, tr_xxyz_xy, tr_xxyz_yy, tr_xxyz_yz, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yz_0, tr_yz_xy, tr_yz_yy, tr_yz_yz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyz_x[i] = tr_z_x[i] - 2.0 * tr_yz_xy[i] * tke_0 - 2.0 * tr_yyz_x[i] * tbe_0 - 2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_y[i] = tr_z_y[i] + tr_yz_0[i] - 2.0 * tr_yz_yy[i] * tke_0 - 2.0 * tr_yyz_y[i] * tbe_0 - 2.0 * tr_xxz_y[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_z[i] = tr_z_z[i] - 2.0 * tr_yz_yz[i] * tke_0 - 2.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 45-48 components of targeted buffer : FP

    auto tr_x_0_y_xzz_x = pbuffer.data(idx_op_geom_110_fp + 45);

    auto tr_x_0_y_xzz_y = pbuffer.data(idx_op_geom_110_fp + 46);

    auto tr_x_0_y_xzz_z = pbuffer.data(idx_op_geom_110_fp + 47);

    #pragma omp simd aligned(tr_x_0_y_xzz_x, tr_x_0_y_xzz_y, tr_x_0_y_xzz_z, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xxzz_0, tr_xxzz_xy, tr_xxzz_yy, tr_xxzz_yz, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_zz_0, tr_zz_xy, tr_zz_yy, tr_zz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xzz_x[i] = -2.0 * tr_zz_xy[i] * tke_0 - 2.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_y[i] = tr_zz_0[i] - 2.0 * tr_zz_yy[i] * tke_0 - 2.0 * tr_yzz_y[i] * tbe_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_z[i] = -2.0 * tr_zz_yz[i] * tke_0 - 2.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 48-51 components of targeted buffer : FP

    auto tr_x_0_y_yyy_x = pbuffer.data(idx_op_geom_110_fp + 48);

    auto tr_x_0_y_yyy_y = pbuffer.data(idx_op_geom_110_fp + 49);

    auto tr_x_0_y_yyy_z = pbuffer.data(idx_op_geom_110_fp + 50);

    #pragma omp simd aligned(tr_x_0_y_yyy_x, tr_x_0_y_yyy_y, tr_x_0_y_yyy_z, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyy_0, tr_xyyy_xy, tr_xyyy_yy, tr_xyyy_yz, tr_xyyyy_x, tr_xyyyy_y, tr_xyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyy_x[i] = -6.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_y[i] = -6.0 * tr_xyy_y[i] * tbe_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_z[i] = -6.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 51-54 components of targeted buffer : FP

    auto tr_x_0_y_yyz_x = pbuffer.data(idx_op_geom_110_fp + 51);

    auto tr_x_0_y_yyz_y = pbuffer.data(idx_op_geom_110_fp + 52);

    auto tr_x_0_y_yyz_z = pbuffer.data(idx_op_geom_110_fp + 53);

    #pragma omp simd aligned(tr_x_0_y_yyz_x, tr_x_0_y_yyz_y, tr_x_0_y_yyz_z, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyyz_0, tr_xyyz_xy, tr_xyyz_yy, tr_xyyz_yz, tr_xyz_x, tr_xyz_y, tr_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyz_x[i] = -4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_y[i] = -4.0 * tr_xyz_y[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_z[i] = -4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 54-57 components of targeted buffer : FP

    auto tr_x_0_y_yzz_x = pbuffer.data(idx_op_geom_110_fp + 54);

    auto tr_x_0_y_yzz_y = pbuffer.data(idx_op_geom_110_fp + 55);

    auto tr_x_0_y_yzz_z = pbuffer.data(idx_op_geom_110_fp + 56);

    #pragma omp simd aligned(tr_x_0_y_yzz_x, tr_x_0_y_yzz_y, tr_x_0_y_yzz_z, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xyzz_0, tr_xyzz_xy, tr_xyzz_yy, tr_xyzz_yz, tr_xzz_x, tr_xzz_y, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yzz_x[i] = -2.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_y[i] = -2.0 * tr_xzz_y[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_z[i] = -2.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 57-60 components of targeted buffer : FP

    auto tr_x_0_y_zzz_x = pbuffer.data(idx_op_geom_110_fp + 57);

    auto tr_x_0_y_zzz_y = pbuffer.data(idx_op_geom_110_fp + 58);

    auto tr_x_0_y_zzz_z = pbuffer.data(idx_op_geom_110_fp + 59);

    #pragma omp simd aligned(tr_x_0_y_zzz_x, tr_x_0_y_zzz_y, tr_x_0_y_zzz_z, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_xzzz_0, tr_xzzz_xy, tr_xzzz_yy, tr_xzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_zzz_x[i] = 4.0 * tr_xzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_y[i] = -2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_z[i] = 4.0 * tr_xzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 60-63 components of targeted buffer : FP

    auto tr_x_0_z_xxx_x = pbuffer.data(idx_op_geom_110_fp + 60);

    auto tr_x_0_z_xxx_y = pbuffer.data(idx_op_geom_110_fp + 61);

    auto tr_x_0_z_xxx_z = pbuffer.data(idx_op_geom_110_fp + 62);

    #pragma omp simd aligned(tr_x_0_z_xxx_x, tr_x_0_z_xxx_y, tr_x_0_z_xxx_z, tr_xx_0, tr_xx_xz, tr_xx_yz, tr_xx_zz, tr_xxxx_0, tr_xxxx_xz, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxz_x, tr_xxxxz_y, tr_xxxxz_z, tr_xxz_x, tr_xxz_y, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxx_x[i] = -6.0 * tr_xx_xz[i] * tke_0 - 6.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxxx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_y[i] = -6.0 * tr_xx_yz[i] * tke_0 - 6.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxxx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_z[i] = 3.0 * tr_xx_0[i] - 6.0 * tr_xx_zz[i] * tke_0 - 6.0 * tr_xxz_z[i] * tbe_0 - 2.0 * tr_xxxx_0[i] * tbe_0 + 4.0 * tr_xxxx_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 63-66 components of targeted buffer : FP

    auto tr_x_0_z_xxy_x = pbuffer.data(idx_op_geom_110_fp + 63);

    auto tr_x_0_z_xxy_y = pbuffer.data(idx_op_geom_110_fp + 64);

    auto tr_x_0_z_xxy_z = pbuffer.data(idx_op_geom_110_fp + 65);

    #pragma omp simd aligned(tr_x_0_z_xxy_x, tr_x_0_z_xxy_y, tr_x_0_z_xxy_z, tr_xxxy_0, tr_xxxy_xz, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xy_0, tr_xy_xz, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_y, tr_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxy_x[i] = -4.0 * tr_xy_xz[i] * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_y[i] = -4.0 * tr_xy_yz[i] * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_z[i] = 2.0 * tr_xy_0[i] - 4.0 * tr_xy_zz[i] * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 66-69 components of targeted buffer : FP

    auto tr_x_0_z_xxz_x = pbuffer.data(idx_op_geom_110_fp + 66);

    auto tr_x_0_z_xxz_y = pbuffer.data(idx_op_geom_110_fp + 67);

    auto tr_x_0_z_xxz_z = pbuffer.data(idx_op_geom_110_fp + 68);

    #pragma omp simd aligned(tr_x_0_z_xxz_x, tr_x_0_z_xxz_y, tr_x_0_z_xxz_z, tr_x_x, tr_x_y, tr_x_z, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxz_0, tr_xxxz_xz, tr_xxxz_yz, tr_xxxz_zz, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, tr_xz_0, tr_xz_xz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_y, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxz_x[i] = 2.0 * tr_x_x[i] - 4.0 * tr_xz_xz[i] * tke_0 - 4.0 * tr_xzz_x[i] * tbe_0 - 2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_y[i] = 2.0 * tr_x_y[i] - 4.0 * tr_xz_yz[i] * tke_0 - 4.0 * tr_xzz_y[i] * tbe_0 - 2.0 * tr_xxx_y[i] * tbe_0 + 4.0 * tr_xxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_z[i] = 2.0 * tr_x_z[i] + 2.0 * tr_xz_0[i] - 4.0 * tr_xz_zz[i] * tke_0 - 4.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xxx_z[i] * tbe_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 69-72 components of targeted buffer : FP

    auto tr_x_0_z_xyy_x = pbuffer.data(idx_op_geom_110_fp + 69);

    auto tr_x_0_z_xyy_y = pbuffer.data(idx_op_geom_110_fp + 70);

    auto tr_x_0_z_xyy_z = pbuffer.data(idx_op_geom_110_fp + 71);

    #pragma omp simd aligned(tr_x_0_z_xyy_x, tr_x_0_z_xyy_y, tr_x_0_z_xyy_z, tr_xxyy_0, tr_xxyy_xz, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_yy_0, tr_yy_xz, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyy_x[i] = -2.0 * tr_yy_xz[i] * tke_0 - 2.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_y[i] = -2.0 * tr_yy_yz[i] * tke_0 - 2.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_z[i] = tr_yy_0[i] - 2.0 * tr_yy_zz[i] * tke_0 - 2.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 72-75 components of targeted buffer : FP

    auto tr_x_0_z_xyz_x = pbuffer.data(idx_op_geom_110_fp + 72);

    auto tr_x_0_z_xyz_y = pbuffer.data(idx_op_geom_110_fp + 73);

    auto tr_x_0_z_xyz_z = pbuffer.data(idx_op_geom_110_fp + 74);

    #pragma omp simd aligned(tr_x_0_z_xyz_x, tr_x_0_z_xyz_y, tr_x_0_z_xyz_z, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyz_0, tr_xxyz_xz, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_y_x, tr_y_y, tr_y_z, tr_yz_0, tr_yz_xz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyz_x[i] = tr_y_x[i] - 2.0 * tr_yz_xz[i] * tke_0 - 2.0 * tr_yzz_x[i] * tbe_0 - 2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_y[i] = tr_y_y[i] - 2.0 * tr_yz_yz[i] * tke_0 - 2.0 * tr_yzz_y[i] * tbe_0 - 2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_z[i] = tr_y_z[i] + tr_yz_0[i] - 2.0 * tr_yz_zz[i] * tke_0 - 2.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_xxy_z[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 75-78 components of targeted buffer : FP

    auto tr_x_0_z_xzz_x = pbuffer.data(idx_op_geom_110_fp + 75);

    auto tr_x_0_z_xzz_y = pbuffer.data(idx_op_geom_110_fp + 76);

    auto tr_x_0_z_xzz_z = pbuffer.data(idx_op_geom_110_fp + 77);

    #pragma omp simd aligned(tr_x_0_z_xzz_x, tr_x_0_z_xzz_y, tr_x_0_z_xzz_z, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xxzz_0, tr_xxzz_xz, tr_xxzz_yz, tr_xxzz_zz, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, tr_zz_xz, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xzz_x[i] = 2.0 * tr_z_x[i] - 2.0 * tr_zz_xz[i] * tke_0 - 2.0 * tr_zzz_x[i] * tbe_0 - 4.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_y[i] = 2.0 * tr_z_y[i] - 2.0 * tr_zz_yz[i] * tke_0 - 2.0 * tr_zzz_y[i] * tbe_0 - 4.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_z[i] = 2.0 * tr_z_z[i] + tr_zz_0[i] - 2.0 * tr_zz_zz[i] * tke_0 - 2.0 * tr_zzz_z[i] * tbe_0 - 4.0 * tr_xxz_z[i] * tbe_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 78-81 components of targeted buffer : FP

    auto tr_x_0_z_yyy_x = pbuffer.data(idx_op_geom_110_fp + 78);

    auto tr_x_0_z_yyy_y = pbuffer.data(idx_op_geom_110_fp + 79);

    auto tr_x_0_z_yyy_z = pbuffer.data(idx_op_geom_110_fp + 80);

    #pragma omp simd aligned(tr_x_0_z_yyy_x, tr_x_0_z_yyy_y, tr_x_0_z_yyy_z, tr_xyyy_0, tr_xyyy_xz, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyy_x[i] = 4.0 * tr_xyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_y[i] = 4.0 * tr_xyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_z[i] = -2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 81-84 components of targeted buffer : FP

    auto tr_x_0_z_yyz_x = pbuffer.data(idx_op_geom_110_fp + 81);

    auto tr_x_0_z_yyz_y = pbuffer.data(idx_op_geom_110_fp + 82);

    auto tr_x_0_z_yyz_z = pbuffer.data(idx_op_geom_110_fp + 83);

    #pragma omp simd aligned(tr_x_0_z_yyz_x, tr_x_0_z_yyz_y, tr_x_0_z_yyz_z, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyz_0, tr_xyyz_xz, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyz_x[i] = -2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_y[i] = -2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_z[i] = -2.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 84-87 components of targeted buffer : FP

    auto tr_x_0_z_yzz_x = pbuffer.data(idx_op_geom_110_fp + 84);

    auto tr_x_0_z_yzz_y = pbuffer.data(idx_op_geom_110_fp + 85);

    auto tr_x_0_z_yzz_z = pbuffer.data(idx_op_geom_110_fp + 86);

    #pragma omp simd aligned(tr_x_0_z_yzz_x, tr_x_0_z_yzz_y, tr_x_0_z_yzz_z, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xyzz_xz, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yzz_x[i] = -4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_y[i] = -4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_z[i] = -4.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 87-90 components of targeted buffer : FP

    auto tr_x_0_z_zzz_x = pbuffer.data(idx_op_geom_110_fp + 87);

    auto tr_x_0_z_zzz_y = pbuffer.data(idx_op_geom_110_fp + 88);

    auto tr_x_0_z_zzz_z = pbuffer.data(idx_op_geom_110_fp + 89);

    #pragma omp simd aligned(tr_x_0_z_zzz_x, tr_x_0_z_zzz_y, tr_x_0_z_zzz_z, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_xzzz_0, tr_xzzz_xz, tr_xzzz_yz, tr_xzzz_zz, tr_xzzzz_x, tr_xzzzz_y, tr_xzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_zzz_x[i] = -6.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_y[i] = -6.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_z[i] = -6.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 90-93 components of targeted buffer : FP

    auto tr_y_0_x_xxx_x = pbuffer.data(idx_op_geom_110_fp + 90);

    auto tr_y_0_x_xxx_y = pbuffer.data(idx_op_geom_110_fp + 91);

    auto tr_y_0_x_xxx_z = pbuffer.data(idx_op_geom_110_fp + 92);

    #pragma omp simd aligned(tr_xxxxy_x, tr_xxxxy_y, tr_xxxxy_z, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_y_0_x_xxx_x, tr_y_0_x_xxx_y, tr_y_0_x_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxx_x[i] = -6.0 * tr_xxy_x[i] * tbe_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_y[i] = -6.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_z[i] = -6.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 93-96 components of targeted buffer : FP

    auto tr_y_0_x_xxy_x = pbuffer.data(idx_op_geom_110_fp + 93);

    auto tr_y_0_x_xxy_y = pbuffer.data(idx_op_geom_110_fp + 94);

    auto tr_y_0_x_xxy_z = pbuffer.data(idx_op_geom_110_fp + 95);

    #pragma omp simd aligned(tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_y_0_x_xxy_x, tr_y_0_x_xxy_y, tr_y_0_x_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxy_x[i] = 2.0 * tr_x_x[i] - 4.0 * tr_xyy_x[i] * tbe_0 + tr_xx_0[i] - 2.0 * tr_xx_xx[i] * tke_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxxyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_y[i] = 2.0 * tr_x_y[i] - 4.0 * tr_xyy_y[i] * tbe_0 - 2.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xxyy_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_y[i] * tbe_0 + 4.0 * tr_xxxyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_z[i] = 2.0 * tr_x_z[i] - 4.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_z[i] * tbe_0 + 4.0 * tr_xxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 96-99 components of targeted buffer : FP

    auto tr_y_0_x_xxz_x = pbuffer.data(idx_op_geom_110_fp + 96);

    auto tr_y_0_x_xxz_y = pbuffer.data(idx_op_geom_110_fp + 97);

    auto tr_y_0_x_xxz_z = pbuffer.data(idx_op_geom_110_fp + 98);

    #pragma omp simd aligned(tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_y_0_x_xxz_x, tr_y_0_x_xxz_y, tr_y_0_x_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxz_x[i] = -4.0 * tr_xyz_x[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_y[i] = -4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_z[i] = -4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 99-102 components of targeted buffer : FP

    auto tr_y_0_x_xyy_x = pbuffer.data(idx_op_geom_110_fp + 99);

    auto tr_y_0_x_xyy_y = pbuffer.data(idx_op_geom_110_fp + 100);

    auto tr_y_0_x_xyy_z = pbuffer.data(idx_op_geom_110_fp + 101);

    #pragma omp simd aligned(tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, tr_xy_0, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_y_0_x_xyy_x, tr_y_0_x_xyy_y, tr_y_0_x_xyy_z, tr_y_x, tr_y_y, tr_y_z, tr_yyy_x, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyy_x[i] = 2.0 * tr_y_x[i] - 2.0 * tr_yyy_x[i] * tbe_0 + 2.0 * tr_xy_0[i] - 4.0 * tr_xy_xx[i] * tke_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_y[i] = 2.0 * tr_y_y[i] - 2.0 * tr_yyy_y[i] * tbe_0 - 4.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xyyy_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_z[i] = 2.0 * tr_y_z[i] - 2.0 * tr_yyy_z[i] * tbe_0 - 4.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xyyy_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 102-105 components of targeted buffer : FP

    auto tr_y_0_x_xyz_x = pbuffer.data(idx_op_geom_110_fp + 102);

    auto tr_y_0_x_xyz_y = pbuffer.data(idx_op_geom_110_fp + 103);

    auto tr_y_0_x_xyz_z = pbuffer.data(idx_op_geom_110_fp + 104);

    #pragma omp simd aligned(tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xz_0, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_y_0_x_xyz_x, tr_y_0_x_xyz_y, tr_y_0_x_xyz_z, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyz_x[i] = tr_z_x[i] - 2.0 * tr_yyz_x[i] * tbe_0 + tr_xz_0[i] - 2.0 * tr_xz_xx[i] * tke_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_y[i] = tr_z_y[i] - 2.0 * tr_yyz_y[i] * tbe_0 - 2.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_z[i] = tr_z_z[i] - 2.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 105-108 components of targeted buffer : FP

    auto tr_y_0_x_xzz_x = pbuffer.data(idx_op_geom_110_fp + 105);

    auto tr_y_0_x_xzz_y = pbuffer.data(idx_op_geom_110_fp + 106);

    auto tr_y_0_x_xzz_z = pbuffer.data(idx_op_geom_110_fp + 107);

    #pragma omp simd aligned(tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_y_0_x_xzz_x, tr_y_0_x_xzz_y, tr_y_0_x_xzz_z, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xzz_x[i] = -2.0 * tr_yzz_x[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_y[i] = -2.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_z[i] = -2.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 108-111 components of targeted buffer : FP

    auto tr_y_0_x_yyy_x = pbuffer.data(idx_op_geom_110_fp + 108);

    auto tr_y_0_x_yyy_y = pbuffer.data(idx_op_geom_110_fp + 109);

    auto tr_y_0_x_yyy_z = pbuffer.data(idx_op_geom_110_fp + 110);

    #pragma omp simd aligned(tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyyy_x, tr_xyyyy_y, tr_xyyyy_z, tr_y_0_x_yyy_x, tr_y_0_x_yyy_y, tr_y_0_x_yyy_z, tr_yy_0, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yyyy_0, tr_yyyy_xx, tr_yyyy_xy, tr_yyyy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyy_x[i] = 3.0 * tr_yy_0[i] - 6.0 * tr_yy_xx[i] * tke_0 - 2.0 * tr_yyyy_0[i] * tbe_0 + 4.0 * tr_yyyy_xx[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_y[i] = -6.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yyyy_xy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_z[i] = -6.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yyyy_xz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 111-114 components of targeted buffer : FP

    auto tr_y_0_x_yyz_x = pbuffer.data(idx_op_geom_110_fp + 111);

    auto tr_y_0_x_yyz_y = pbuffer.data(idx_op_geom_110_fp + 112);

    auto tr_y_0_x_yyz_z = pbuffer.data(idx_op_geom_110_fp + 113);

    #pragma omp simd aligned(tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_y_0_x_yyz_x, tr_y_0_x_yyz_y, tr_y_0_x_yyz_z, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yz_0, tr_yz_xx, tr_yz_xy, tr_yz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyz_x[i] = 2.0 * tr_yz_0[i] - 4.0 * tr_yz_xx[i] * tke_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_y[i] = -4.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yyyz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_z[i] = -4.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yyyz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 114-117 components of targeted buffer : FP

    auto tr_y_0_x_yzz_x = pbuffer.data(idx_op_geom_110_fp + 114);

    auto tr_y_0_x_yzz_y = pbuffer.data(idx_op_geom_110_fp + 115);

    auto tr_y_0_x_yzz_z = pbuffer.data(idx_op_geom_110_fp + 116);

    #pragma omp simd aligned(tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_y_0_x_yzz_x, tr_y_0_x_yzz_y, tr_y_0_x_yzz_z, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_zz_0, tr_zz_xx, tr_zz_xy, tr_zz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yzz_x[i] = tr_zz_0[i] - 2.0 * tr_zz_xx[i] * tke_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_y[i] = -2.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_z[i] = -2.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 117-120 components of targeted buffer : FP

    auto tr_y_0_x_zzz_x = pbuffer.data(idx_op_geom_110_fp + 117);

    auto tr_y_0_x_zzz_y = pbuffer.data(idx_op_geom_110_fp + 118);

    auto tr_y_0_x_zzz_z = pbuffer.data(idx_op_geom_110_fp + 119);

    #pragma omp simd aligned(tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_y_0_x_zzz_x, tr_y_0_x_zzz_y, tr_y_0_x_zzz_z, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_zzz_x[i] = -2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_y[i] = 4.0 * tr_yzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_z[i] = 4.0 * tr_yzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 120-123 components of targeted buffer : FP

    auto tr_y_0_y_xxx_x = pbuffer.data(idx_op_geom_110_fp + 120);

    auto tr_y_0_y_xxx_y = pbuffer.data(idx_op_geom_110_fp + 121);

    auto tr_y_0_y_xxx_z = pbuffer.data(idx_op_geom_110_fp + 122);

    #pragma omp simd aligned(tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxy_0, tr_xxxy_xy, tr_xxxy_yy, tr_xxxy_yz, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, tr_y_0_y_xxx_x, tr_y_0_y_xxx_y, tr_y_0_y_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxx_x[i] = -2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_y[i] = -2.0 * tr_xxx_y[i] * tbe_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_z[i] = -2.0 * tr_xxx_z[i] * tbe_0 + 4.0 * tr_xxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 123-126 components of targeted buffer : FP

    auto tr_y_0_y_xxy_x = pbuffer.data(idx_op_geom_110_fp + 123);

    auto tr_y_0_y_xxy_y = pbuffer.data(idx_op_geom_110_fp + 124);

    auto tr_y_0_y_xxy_z = pbuffer.data(idx_op_geom_110_fp + 125);

    #pragma omp simd aligned(tr_xx_0, tr_xx_xy, tr_xx_yy, tr_xx_yz, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyy_0, tr_xxyy_xy, tr_xxyy_yy, tr_xxyy_yz, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, tr_y_0_y_xxy_x, tr_y_0_y_xxy_y, tr_y_0_y_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxy_x[i] = -2.0 * tr_xx_xy[i] * tke_0 - 6.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_y[i] = tr_xx_0[i] - 2.0 * tr_xx_yy[i] * tke_0 - 6.0 * tr_xxy_y[i] * tbe_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_z[i] = -2.0 * tr_xx_yz[i] * tke_0 - 6.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 126-129 components of targeted buffer : FP

    auto tr_y_0_y_xxz_x = pbuffer.data(idx_op_geom_110_fp + 126);

    auto tr_y_0_y_xxz_y = pbuffer.data(idx_op_geom_110_fp + 127);

    auto tr_y_0_y_xxz_z = pbuffer.data(idx_op_geom_110_fp + 128);

    #pragma omp simd aligned(tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxyz_0, tr_xxyz_xy, tr_xxyz_yy, tr_xxyz_yz, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_y_0_y_xxz_x, tr_y_0_y_xxz_y, tr_y_0_y_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxz_x[i] = -2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_y[i] = -2.0 * tr_xxz_y[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_z[i] = -2.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 129-132 components of targeted buffer : FP

    auto tr_y_0_y_xyy_x = pbuffer.data(idx_op_geom_110_fp + 129);

    auto tr_y_0_y_xyy_y = pbuffer.data(idx_op_geom_110_fp + 130);

    auto tr_y_0_y_xyy_z = pbuffer.data(idx_op_geom_110_fp + 131);

    #pragma omp simd aligned(tr_x_x, tr_x_y, tr_x_z, tr_xy_0, tr_xy_xy, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyy_0, tr_xyyy_xy, tr_xyyy_yy, tr_xyyy_yz, tr_xyyyy_x, tr_xyyyy_y, tr_xyyyy_z, tr_y_0_y_xyy_x, tr_y_0_y_xyy_y, tr_y_0_y_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyy_x[i] = 2.0 * tr_x_x[i] - 4.0 * tr_xy_xy[i] * tke_0 - 10.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_y[i] = 2.0 * tr_x_y[i] + 2.0 * tr_xy_0[i] - 4.0 * tr_xy_yy[i] * tke_0 - 10.0 * tr_xyy_y[i] * tbe_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_z[i] = 2.0 * tr_x_z[i] - 4.0 * tr_xy_yz[i] * tke_0 - 10.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 132-135 components of targeted buffer : FP

    auto tr_y_0_y_xyz_x = pbuffer.data(idx_op_geom_110_fp + 132);

    auto tr_y_0_y_xyz_y = pbuffer.data(idx_op_geom_110_fp + 133);

    auto tr_y_0_y_xyz_z = pbuffer.data(idx_op_geom_110_fp + 134);

    #pragma omp simd aligned(tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyyz_0, tr_xyyz_xy, tr_xyyz_yy, tr_xyyz_yz, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xz_0, tr_xz_xy, tr_xz_yy, tr_xz_yz, tr_y_0_y_xyz_x, tr_y_0_y_xyz_y, tr_y_0_y_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyz_x[i] = -2.0 * tr_xz_xy[i] * tke_0 - 6.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_y[i] = tr_xz_0[i] - 2.0 * tr_xz_yy[i] * tke_0 - 6.0 * tr_xyz_y[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_z[i] = -2.0 * tr_xz_yz[i] * tke_0 - 6.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 135-138 components of targeted buffer : FP

    auto tr_y_0_y_xzz_x = pbuffer.data(idx_op_geom_110_fp + 135);

    auto tr_y_0_y_xzz_y = pbuffer.data(idx_op_geom_110_fp + 136);

    auto tr_y_0_y_xzz_z = pbuffer.data(idx_op_geom_110_fp + 137);

    #pragma omp simd aligned(tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xyzz_0, tr_xyzz_xy, tr_xyzz_yy, tr_xyzz_yz, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_y_0_y_xzz_x, tr_y_0_y_xzz_y, tr_y_0_y_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xzz_x[i] = -2.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_y[i] = -2.0 * tr_xzz_y[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_z[i] = -2.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 138-141 components of targeted buffer : FP

    auto tr_y_0_y_yyy_x = pbuffer.data(idx_op_geom_110_fp + 138);

    auto tr_y_0_y_yyy_y = pbuffer.data(idx_op_geom_110_fp + 139);

    auto tr_y_0_y_yyy_z = pbuffer.data(idx_op_geom_110_fp + 140);

    #pragma omp simd aligned(tr_y_0_y_yyy_x, tr_y_0_y_yyy_y, tr_y_0_y_yyy_z, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yy_xy, tr_yy_yy, tr_yy_yz, tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyy_0, tr_yyyy_xy, tr_yyyy_yy, tr_yyyy_yz, tr_yyyyy_x, tr_yyyyy_y, tr_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyy_x[i] = 6.0 * tr_y_x[i] - 6.0 * tr_yy_xy[i] * tke_0 - 14.0 * tr_yyy_x[i] * tbe_0 + 4.0 * tr_yyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_y[i] = 6.0 * tr_y_y[i] + 3.0 * tr_yy_0[i] - 6.0 * tr_yy_yy[i] * tke_0 - 14.0 * tr_yyy_y[i] * tbe_0 - 2.0 * tr_yyyy_0[i] * tbe_0 + 4.0 * tr_yyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_z[i] = 6.0 * tr_y_z[i] - 6.0 * tr_yy_yz[i] * tke_0 - 14.0 * tr_yyy_z[i] * tbe_0 + 4.0 * tr_yyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 141-144 components of targeted buffer : FP

    auto tr_y_0_y_yyz_x = pbuffer.data(idx_op_geom_110_fp + 141);

    auto tr_y_0_y_yyz_y = pbuffer.data(idx_op_geom_110_fp + 142);

    auto tr_y_0_y_yyz_z = pbuffer.data(idx_op_geom_110_fp + 143);

    #pragma omp simd aligned(tr_y_0_y_yyz_x, tr_y_0_y_yyz_y, tr_y_0_y_yyz_z, tr_yyyyz_x, tr_yyyyz_y, tr_yyyyz_z, tr_yyyz_0, tr_yyyz_xy, tr_yyyz_yy, tr_yyyz_yz, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yz_0, tr_yz_xy, tr_yz_yy, tr_yz_yz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyz_x[i] = 2.0 * tr_z_x[i] - 4.0 * tr_yz_xy[i] * tke_0 - 10.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_y[i] = 2.0 * tr_z_y[i] + 2.0 * tr_yz_0[i] - 4.0 * tr_yz_yy[i] * tke_0 - 10.0 * tr_yyz_y[i] * tbe_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_z[i] = 2.0 * tr_z_z[i] - 4.0 * tr_yz_yz[i] * tke_0 - 10.0 * tr_yyz_z[i] * tbe_0 + 4.0 * tr_yyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 144-147 components of targeted buffer : FP

    auto tr_y_0_y_yzz_x = pbuffer.data(idx_op_geom_110_fp + 144);

    auto tr_y_0_y_yzz_y = pbuffer.data(idx_op_geom_110_fp + 145);

    auto tr_y_0_y_yzz_z = pbuffer.data(idx_op_geom_110_fp + 146);

    #pragma omp simd aligned(tr_y_0_y_yzz_x, tr_y_0_y_yzz_y, tr_y_0_y_yzz_z, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, tr_yyzz_0, tr_yyzz_xy, tr_yyzz_yy, tr_yyzz_yz, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_zz_0, tr_zz_xy, tr_zz_yy, tr_zz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yzz_x[i] = -2.0 * tr_zz_xy[i] * tke_0 - 6.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_y[i] = tr_zz_0[i] - 2.0 * tr_zz_yy[i] * tke_0 - 6.0 * tr_yzz_y[i] * tbe_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_z[i] = -2.0 * tr_zz_yz[i] * tke_0 - 6.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_yyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 147-150 components of targeted buffer : FP

    auto tr_y_0_y_zzz_x = pbuffer.data(idx_op_geom_110_fp + 147);

    auto tr_y_0_y_zzz_y = pbuffer.data(idx_op_geom_110_fp + 148);

    auto tr_y_0_y_zzz_z = pbuffer.data(idx_op_geom_110_fp + 149);

    #pragma omp simd aligned(tr_y_0_y_zzz_x, tr_y_0_y_zzz_y, tr_y_0_y_zzz_z, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, tr_yzzz_0, tr_yzzz_xy, tr_yzzz_yy, tr_yzzz_yz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_zzz_x[i] = -2.0 * tr_zzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_y[i] = -2.0 * tr_zzz_y[i] * tbe_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_z[i] = -2.0 * tr_zzz_z[i] * tbe_0 + 4.0 * tr_yzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 150-153 components of targeted buffer : FP

    auto tr_y_0_z_xxx_x = pbuffer.data(idx_op_geom_110_fp + 150);

    auto tr_y_0_z_xxx_y = pbuffer.data(idx_op_geom_110_fp + 151);

    auto tr_y_0_z_xxx_z = pbuffer.data(idx_op_geom_110_fp + 152);

    #pragma omp simd aligned(tr_xxxy_0, tr_xxxy_xz, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_y_0_z_xxx_x, tr_y_0_z_xxx_y, tr_y_0_z_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxx_x[i] = 4.0 * tr_xxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_y[i] = 4.0 * tr_xxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_z[i] = -2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 153-156 components of targeted buffer : FP

    auto tr_y_0_z_xxy_x = pbuffer.data(idx_op_geom_110_fp + 153);

    auto tr_y_0_z_xxy_y = pbuffer.data(idx_op_geom_110_fp + 154);

    auto tr_y_0_z_xxy_z = pbuffer.data(idx_op_geom_110_fp + 155);

    #pragma omp simd aligned(tr_xx_0, tr_xx_xz, tr_xx_yz, tr_xx_zz, tr_xxyy_0, tr_xxyy_xz, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_y_0_z_xxy_x, tr_y_0_z_xxy_y, tr_y_0_z_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxy_x[i] = -2.0 * tr_xx_xz[i] * tke_0 - 2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_y[i] = -2.0 * tr_xx_yz[i] * tke_0 - 2.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_z[i] = tr_xx_0[i] - 2.0 * tr_xx_zz[i] * tke_0 - 2.0 * tr_xxz_z[i] * tbe_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 156-159 components of targeted buffer : FP

    auto tr_y_0_z_xxz_x = pbuffer.data(idx_op_geom_110_fp + 156);

    auto tr_y_0_z_xxz_y = pbuffer.data(idx_op_geom_110_fp + 157);

    auto tr_y_0_z_xxz_z = pbuffer.data(idx_op_geom_110_fp + 158);

    #pragma omp simd aligned(tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyz_0, tr_xxyz_xz, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_y_0_z_xxz_x, tr_y_0_z_xxz_y, tr_y_0_z_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxz_x[i] = -2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_y[i] = -2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_z[i] = -2.0 * tr_xxy_z[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 159-162 components of targeted buffer : FP

    auto tr_y_0_z_xyy_x = pbuffer.data(idx_op_geom_110_fp + 159);

    auto tr_y_0_z_xyy_y = pbuffer.data(idx_op_geom_110_fp + 160);

    auto tr_y_0_z_xyy_z = pbuffer.data(idx_op_geom_110_fp + 161);

    #pragma omp simd aligned(tr_xy_0, tr_xy_xz, tr_xy_yz, tr_xy_zz, tr_xyyy_0, tr_xyyy_xz, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_y_0_z_xyy_x, tr_y_0_z_xyy_y, tr_y_0_z_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyy_x[i] = -4.0 * tr_xy_xz[i] * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_y[i] = -4.0 * tr_xy_yz[i] * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_z[i] = 2.0 * tr_xy_0[i] - 4.0 * tr_xy_zz[i] * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 162-165 components of targeted buffer : FP

    auto tr_y_0_z_xyz_x = pbuffer.data(idx_op_geom_110_fp + 162);

    auto tr_y_0_z_xyz_y = pbuffer.data(idx_op_geom_110_fp + 163);

    auto tr_y_0_z_xyz_z = pbuffer.data(idx_op_geom_110_fp + 164);

    #pragma omp simd aligned(tr_x_x, tr_x_y, tr_x_z, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyz_0, tr_xyyz_xz, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xz_0, tr_xz_xz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_y_0_z_xyz_x, tr_y_0_z_xyz_y, tr_y_0_z_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyz_x[i] = tr_x_x[i] - 2.0 * tr_xz_xz[i] * tke_0 - 2.0 * tr_xzz_x[i] * tbe_0 - 2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_y[i] = tr_x_y[i] - 2.0 * tr_xz_yz[i] * tke_0 - 2.0 * tr_xzz_y[i] * tbe_0 - 2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_z[i] = tr_x_z[i] + tr_xz_0[i] - 2.0 * tr_xz_zz[i] * tke_0 - 2.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 165-168 components of targeted buffer : FP

    auto tr_y_0_z_xzz_x = pbuffer.data(idx_op_geom_110_fp + 165);

    auto tr_y_0_z_xzz_y = pbuffer.data(idx_op_geom_110_fp + 166);

    auto tr_y_0_z_xzz_z = pbuffer.data(idx_op_geom_110_fp + 167);

    #pragma omp simd aligned(tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xyzz_xz, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_y_0_z_xzz_x, tr_y_0_z_xzz_y, tr_y_0_z_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xzz_x[i] = -4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_y[i] = -4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_z[i] = -4.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 168-171 components of targeted buffer : FP

    auto tr_y_0_z_yyy_x = pbuffer.data(idx_op_geom_110_fp + 168);

    auto tr_y_0_z_yyy_y = pbuffer.data(idx_op_geom_110_fp + 169);

    auto tr_y_0_z_yyy_z = pbuffer.data(idx_op_geom_110_fp + 170);

    #pragma omp simd aligned(tr_y_0_z_yyy_x, tr_y_0_z_yyy_y, tr_y_0_z_yyy_z, tr_yy_0, tr_yy_xz, tr_yy_yz, tr_yy_zz, tr_yyyy_0, tr_yyyy_xz, tr_yyyy_yz, tr_yyyy_zz, tr_yyyyz_x, tr_yyyyz_y, tr_yyyyz_z, tr_yyz_x, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyy_x[i] = -6.0 * tr_yy_xz[i] * tke_0 - 6.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_y[i] = -6.0 * tr_yy_yz[i] * tke_0 - 6.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_yyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_z[i] = 3.0 * tr_yy_0[i] - 6.0 * tr_yy_zz[i] * tke_0 - 6.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_yyyy_0[i] * tbe_0 + 4.0 * tr_yyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 171-174 components of targeted buffer : FP

    auto tr_y_0_z_yyz_x = pbuffer.data(idx_op_geom_110_fp + 171);

    auto tr_y_0_z_yyz_y = pbuffer.data(idx_op_geom_110_fp + 172);

    auto tr_y_0_z_yyz_z = pbuffer.data(idx_op_geom_110_fp + 173);

    #pragma omp simd aligned(tr_y_0_z_yyz_x, tr_y_0_z_yyz_y, tr_y_0_z_yyz_z, tr_y_x, tr_y_y, tr_y_z, tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyz_0, tr_yyyz_xz, tr_yyyz_yz, tr_yyyz_zz, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, tr_yz_0, tr_yz_xz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyz_x[i] = 2.0 * tr_y_x[i] - 4.0 * tr_yz_xz[i] * tke_0 - 4.0 * tr_yzz_x[i] * tbe_0 - 2.0 * tr_yyy_x[i] * tbe_0 + 4.0 * tr_yyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_y[i] = 2.0 * tr_y_y[i] - 4.0 * tr_yz_yz[i] * tke_0 - 4.0 * tr_yzz_y[i] * tbe_0 - 2.0 * tr_yyy_y[i] * tbe_0 + 4.0 * tr_yyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_z[i] = 2.0 * tr_y_z[i] + 2.0 * tr_yz_0[i] - 4.0 * tr_yz_zz[i] * tke_0 - 4.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_yyy_z[i] * tbe_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 174-177 components of targeted buffer : FP

    auto tr_y_0_z_yzz_x = pbuffer.data(idx_op_geom_110_fp + 174);

    auto tr_y_0_z_yzz_y = pbuffer.data(idx_op_geom_110_fp + 175);

    auto tr_y_0_z_yzz_z = pbuffer.data(idx_op_geom_110_fp + 176);

    #pragma omp simd aligned(tr_y_0_z_yzz_x, tr_y_0_z_yzz_y, tr_y_0_z_yzz_z, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yyzz_0, tr_yyzz_xz, tr_yyzz_yz, tr_yyzz_zz, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, tr_zz_xz, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yzz_x[i] = 2.0 * tr_z_x[i] - 2.0 * tr_zz_xz[i] * tke_0 - 2.0 * tr_zzz_x[i] * tbe_0 - 4.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_y[i] = 2.0 * tr_z_y[i] - 2.0 * tr_zz_yz[i] * tke_0 - 2.0 * tr_zzz_y[i] * tbe_0 - 4.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_yyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_z[i] = 2.0 * tr_z_z[i] + tr_zz_0[i] - 2.0 * tr_zz_zz[i] * tke_0 - 2.0 * tr_zzz_z[i] * tbe_0 - 4.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 177-180 components of targeted buffer : FP

    auto tr_y_0_z_zzz_x = pbuffer.data(idx_op_geom_110_fp + 177);

    auto tr_y_0_z_zzz_y = pbuffer.data(idx_op_geom_110_fp + 178);

    auto tr_y_0_z_zzz_z = pbuffer.data(idx_op_geom_110_fp + 179);

    #pragma omp simd aligned(tr_y_0_z_zzz_x, tr_y_0_z_zzz_y, tr_y_0_z_zzz_z, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_yzzz_0, tr_yzzz_xz, tr_yzzz_yz, tr_yzzz_zz, tr_yzzzz_x, tr_yzzzz_y, tr_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_zzz_x[i] = -6.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_y[i] = -6.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_yzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_z[i] = -6.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 180-183 components of targeted buffer : FP

    auto tr_z_0_x_xxx_x = pbuffer.data(idx_op_geom_110_fp + 180);

    auto tr_z_0_x_xxx_y = pbuffer.data(idx_op_geom_110_fp + 181);

    auto tr_z_0_x_xxx_z = pbuffer.data(idx_op_geom_110_fp + 182);

    #pragma omp simd aligned(tr_xxxxz_x, tr_xxxxz_y, tr_xxxxz_z, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_z_0_x_xxx_x, tr_z_0_x_xxx_y, tr_z_0_x_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxx_x[i] = -6.0 * tr_xxz_x[i] * tbe_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_y[i] = -6.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_z[i] = -6.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 183-186 components of targeted buffer : FP

    auto tr_z_0_x_xxy_x = pbuffer.data(idx_op_geom_110_fp + 183);

    auto tr_z_0_x_xxy_y = pbuffer.data(idx_op_geom_110_fp + 184);

    auto tr_z_0_x_xxy_z = pbuffer.data(idx_op_geom_110_fp + 185);

    #pragma omp simd aligned(tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_z_0_x_xxy_x, tr_z_0_x_xxy_y, tr_z_0_x_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxy_x[i] = -4.0 * tr_xyz_x[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_y[i] = -4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_z[i] = -4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 186-189 components of targeted buffer : FP

    auto tr_z_0_x_xxz_x = pbuffer.data(idx_op_geom_110_fp + 186);

    auto tr_z_0_x_xxz_y = pbuffer.data(idx_op_geom_110_fp + 187);

    auto tr_z_0_x_xxz_z = pbuffer.data(idx_op_geom_110_fp + 188);

    #pragma omp simd aligned(tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_z_0_x_xxz_x, tr_z_0_x_xxz_y, tr_z_0_x_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxz_x[i] = 2.0 * tr_x_x[i] - 4.0 * tr_xzz_x[i] * tbe_0 + tr_xx_0[i] - 2.0 * tr_xx_xx[i] * tke_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxxzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_y[i] = 2.0 * tr_x_y[i] - 4.0 * tr_xzz_y[i] * tbe_0 - 2.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_y[i] * tbe_0 + 4.0 * tr_xxxzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_z[i] = 2.0 * tr_x_z[i] - 4.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xxzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_z[i] * tbe_0 + 4.0 * tr_xxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 189-192 components of targeted buffer : FP

    auto tr_z_0_x_xyy_x = pbuffer.data(idx_op_geom_110_fp + 189);

    auto tr_z_0_x_xyy_y = pbuffer.data(idx_op_geom_110_fp + 190);

    auto tr_z_0_x_xyy_z = pbuffer.data(idx_op_geom_110_fp + 191);

    #pragma omp simd aligned(tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_z_0_x_xyy_x, tr_z_0_x_xyy_y, tr_z_0_x_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyy_x[i] = -2.0 * tr_yyz_x[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_y[i] = -2.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_z[i] = -2.0 * tr_yyz_z[i] * tbe_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 192-195 components of targeted buffer : FP

    auto tr_z_0_x_xyz_x = pbuffer.data(idx_op_geom_110_fp + 192);

    auto tr_z_0_x_xyz_y = pbuffer.data(idx_op_geom_110_fp + 193);

    auto tr_z_0_x_xyz_z = pbuffer.data(idx_op_geom_110_fp + 194);

    #pragma omp simd aligned(tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xy_0, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_y_x, tr_y_y, tr_y_z, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_z_0_x_xyz_x, tr_z_0_x_xyz_y, tr_z_0_x_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyz_x[i] = tr_y_x[i] - 2.0 * tr_yzz_x[i] * tbe_0 + tr_xy_0[i] - 2.0 * tr_xy_xx[i] * tke_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_y[i] = tr_y_y[i] - 2.0 * tr_yzz_y[i] * tbe_0 - 2.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_z[i] = tr_y_z[i] - 2.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 195-198 components of targeted buffer : FP

    auto tr_z_0_x_xzz_x = pbuffer.data(idx_op_geom_110_fp + 195);

    auto tr_z_0_x_xzz_y = pbuffer.data(idx_op_geom_110_fp + 196);

    auto tr_z_0_x_xzz_z = pbuffer.data(idx_op_geom_110_fp + 197);

    #pragma omp simd aligned(tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, tr_xz_0, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_z_0_x_xzz_x, tr_z_0_x_xzz_y, tr_z_0_x_xzz_z, tr_z_x, tr_z_y, tr_z_z, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xzz_x[i] = 2.0 * tr_z_x[i] - 2.0 * tr_zzz_x[i] * tbe_0 + 2.0 * tr_xz_0[i] - 4.0 * tr_xz_xx[i] * tke_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_y[i] = 2.0 * tr_z_y[i] - 2.0 * tr_zzz_y[i] * tbe_0 - 4.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_z[i] = 2.0 * tr_z_z[i] - 2.0 * tr_zzz_z[i] * tbe_0 - 4.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xzzz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 198-201 components of targeted buffer : FP

    auto tr_z_0_x_yyy_x = pbuffer.data(idx_op_geom_110_fp + 198);

    auto tr_z_0_x_yyy_y = pbuffer.data(idx_op_geom_110_fp + 199);

    auto tr_z_0_x_yyy_z = pbuffer.data(idx_op_geom_110_fp + 200);

    #pragma omp simd aligned(tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_z_0_x_yyy_x, tr_z_0_x_yyy_y, tr_z_0_x_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyy_x[i] = -2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_y[i] = 4.0 * tr_yyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_z[i] = 4.0 * tr_yyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 201-204 components of targeted buffer : FP

    auto tr_z_0_x_yyz_x = pbuffer.data(idx_op_geom_110_fp + 201);

    auto tr_z_0_x_yyz_y = pbuffer.data(idx_op_geom_110_fp + 202);

    auto tr_z_0_x_yyz_z = pbuffer.data(idx_op_geom_110_fp + 203);

    #pragma omp simd aligned(tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_yy_0, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_z_0_x_yyz_x, tr_z_0_x_yyz_y, tr_z_0_x_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyz_x[i] = tr_yy_0[i] - 2.0 * tr_yy_xx[i] * tke_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_y[i] = -2.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_z[i] = -2.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 204-207 components of targeted buffer : FP

    auto tr_z_0_x_yzz_x = pbuffer.data(idx_op_geom_110_fp + 204);

    auto tr_z_0_x_yzz_y = pbuffer.data(idx_op_geom_110_fp + 205);

    auto tr_z_0_x_yzz_z = pbuffer.data(idx_op_geom_110_fp + 206);

    #pragma omp simd aligned(tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_yz_0, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_z_0_x_yzz_x, tr_z_0_x_yzz_y, tr_z_0_x_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yzz_x[i] = 2.0 * tr_yz_0[i] - 4.0 * tr_yz_xx[i] * tke_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_y[i] = -4.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_z[i] = -4.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yzzz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 207-210 components of targeted buffer : FP

    auto tr_z_0_x_zzz_x = pbuffer.data(idx_op_geom_110_fp + 207);

    auto tr_z_0_x_zzz_y = pbuffer.data(idx_op_geom_110_fp + 208);

    auto tr_z_0_x_zzz_z = pbuffer.data(idx_op_geom_110_fp + 209);

    #pragma omp simd aligned(tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_xzzzz_x, tr_xzzzz_y, tr_xzzzz_z, tr_z_0_x_zzz_x, tr_z_0_x_zzz_y, tr_z_0_x_zzz_z, tr_zz_0, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zzzz_0, tr_zzzz_xx, tr_zzzz_xy, tr_zzzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_zzz_x[i] = 3.0 * tr_zz_0[i] - 6.0 * tr_zz_xx[i] * tke_0 - 2.0 * tr_zzzz_0[i] * tbe_0 + 4.0 * tr_zzzz_xx[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_y[i] = -6.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zzzz_xy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_z[i] = -6.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_zzzz_xz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 210-213 components of targeted buffer : FP

    auto tr_z_0_y_xxx_x = pbuffer.data(idx_op_geom_110_fp + 210);

    auto tr_z_0_y_xxx_y = pbuffer.data(idx_op_geom_110_fp + 211);

    auto tr_z_0_y_xxx_z = pbuffer.data(idx_op_geom_110_fp + 212);

    #pragma omp simd aligned(tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxxz_0, tr_xxxz_xy, tr_xxxz_yy, tr_xxxz_yz, tr_z_0_y_xxx_x, tr_z_0_y_xxx_y, tr_z_0_y_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxx_x[i] = 4.0 * tr_xxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_y[i] = -2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_z[i] = 4.0 * tr_xxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 213-216 components of targeted buffer : FP

    auto tr_z_0_y_xxy_x = pbuffer.data(idx_op_geom_110_fp + 213);

    auto tr_z_0_y_xxy_y = pbuffer.data(idx_op_geom_110_fp + 214);

    auto tr_z_0_y_xxy_z = pbuffer.data(idx_op_geom_110_fp + 215);

    #pragma omp simd aligned(tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxyz_0, tr_xxyz_xy, tr_xxyz_yy, tr_xxyz_yz, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_z_0_y_xxy_x, tr_z_0_y_xxy_y, tr_z_0_y_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxy_x[i] = -2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_y[i] = -2.0 * tr_xxz_y[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_z[i] = -2.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 216-219 components of targeted buffer : FP

    auto tr_z_0_y_xxz_x = pbuffer.data(idx_op_geom_110_fp + 216);

    auto tr_z_0_y_xxz_y = pbuffer.data(idx_op_geom_110_fp + 217);

    auto tr_z_0_y_xxz_z = pbuffer.data(idx_op_geom_110_fp + 218);

    #pragma omp simd aligned(tr_xx_0, tr_xx_xy, tr_xx_yy, tr_xx_yz, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xxzz_0, tr_xxzz_xy, tr_xxzz_yy, tr_xxzz_yz, tr_z_0_y_xxz_x, tr_z_0_y_xxz_y, tr_z_0_y_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxz_x[i] = -2.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_y[i] = tr_xx_0[i] - 2.0 * tr_xx_yy[i] * tke_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_z[i] = -2.0 * tr_xx_yz[i] * tke_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 219-222 components of targeted buffer : FP

    auto tr_z_0_y_xyy_x = pbuffer.data(idx_op_geom_110_fp + 219);

    auto tr_z_0_y_xyy_y = pbuffer.data(idx_op_geom_110_fp + 220);

    auto tr_z_0_y_xyy_z = pbuffer.data(idx_op_geom_110_fp + 221);

    #pragma omp simd aligned(tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyyz_0, tr_xyyz_xy, tr_xyyz_yy, tr_xyyz_yz, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_z_0_y_xyy_x, tr_z_0_y_xyy_y, tr_z_0_y_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyy_x[i] = -4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_y[i] = -4.0 * tr_xyz_y[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_z[i] = -4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 222-225 components of targeted buffer : FP

    auto tr_z_0_y_xyz_x = pbuffer.data(idx_op_geom_110_fp + 222);

    auto tr_z_0_y_xyz_y = pbuffer.data(idx_op_geom_110_fp + 223);

    auto tr_z_0_y_xyz_z = pbuffer.data(idx_op_geom_110_fp + 224);

    #pragma omp simd aligned(tr_x_x, tr_x_y, tr_x_z, tr_xy_0, tr_xy_xy, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xyzz_0, tr_xyzz_xy, tr_xyzz_yy, tr_xyzz_yz, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_z_0_y_xyz_x, tr_z_0_y_xyz_y, tr_z_0_y_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyz_x[i] = tr_x_x[i] - 2.0 * tr_xzz_x[i] * tbe_0 - 2.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_y[i] = tr_x_y[i] - 2.0 * tr_xzz_y[i] * tbe_0 + tr_xy_0[i] - 2.0 * tr_xy_yy[i] * tke_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_z[i] = tr_x_z[i] - 2.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xy_yz[i] * tke_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 225-228 components of targeted buffer : FP

    auto tr_z_0_y_xzz_x = pbuffer.data(idx_op_geom_110_fp + 225);

    auto tr_z_0_y_xzz_y = pbuffer.data(idx_op_geom_110_fp + 226);

    auto tr_z_0_y_xzz_z = pbuffer.data(idx_op_geom_110_fp + 227);

    #pragma omp simd aligned(tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_xz_0, tr_xz_xy, tr_xz_yy, tr_xz_yz, tr_xzzz_0, tr_xzzz_xy, tr_xzzz_yy, tr_xzzz_yz, tr_z_0_y_xzz_x, tr_z_0_y_xzz_y, tr_z_0_y_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xzz_x[i] = -4.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_y[i] = 2.0 * tr_xz_0[i] - 4.0 * tr_xz_yy[i] * tke_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_yy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_z[i] = -4.0 * tr_xz_yz[i] * tke_0 + 4.0 * tr_xzzz_yz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 228-231 components of targeted buffer : FP

    auto tr_z_0_y_yyy_x = pbuffer.data(idx_op_geom_110_fp + 228);

    auto tr_z_0_y_yyy_y = pbuffer.data(idx_op_geom_110_fp + 229);

    auto tr_z_0_y_yyy_z = pbuffer.data(idx_op_geom_110_fp + 230);

    #pragma omp simd aligned(tr_yyyyz_x, tr_yyyyz_y, tr_yyyyz_z, tr_yyyz_0, tr_yyyz_xy, tr_yyyz_yy, tr_yyyz_yz, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_z_0_y_yyy_x, tr_z_0_y_yyy_y, tr_z_0_y_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyy_x[i] = -6.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_y[i] = -6.0 * tr_yyz_y[i] * tbe_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_z[i] = -6.0 * tr_yyz_z[i] * tbe_0 + 4.0 * tr_yyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 231-234 components of targeted buffer : FP

    auto tr_z_0_y_yyz_x = pbuffer.data(idx_op_geom_110_fp + 231);

    auto tr_z_0_y_yyz_y = pbuffer.data(idx_op_geom_110_fp + 232);

    auto tr_z_0_y_yyz_z = pbuffer.data(idx_op_geom_110_fp + 233);

    #pragma omp simd aligned(tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yy_xy, tr_yy_yy, tr_yy_yz, tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, tr_yyzz_0, tr_yyzz_xy, tr_yyzz_yy, tr_yyzz_yz, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_z_0_y_yyz_x, tr_z_0_y_yyz_y, tr_z_0_y_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyz_x[i] = 2.0 * tr_y_x[i] - 4.0 * tr_yzz_x[i] * tbe_0 - 2.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_x[i] * tbe_0 + 4.0 * tr_yyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_y[i] = 2.0 * tr_y_y[i] - 4.0 * tr_yzz_y[i] * tbe_0 + tr_yy_0[i] - 2.0 * tr_yy_yy[i] * tke_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_y[i] * tbe_0 + 4.0 * tr_yyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_z[i] = 2.0 * tr_y_z[i] - 4.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_yy_yz[i] * tke_0 + 4.0 * tr_yyzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_z[i] * tbe_0 + 4.0 * tr_yyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 234-237 components of targeted buffer : FP

    auto tr_z_0_y_yzz_x = pbuffer.data(idx_op_geom_110_fp + 234);

    auto tr_z_0_y_yzz_y = pbuffer.data(idx_op_geom_110_fp + 235);

    auto tr_z_0_y_yzz_z = pbuffer.data(idx_op_geom_110_fp + 236);

    #pragma omp simd aligned(tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, tr_yz_0, tr_yz_xy, tr_yz_yy, tr_yz_yz, tr_yzzz_0, tr_yzzz_xy, tr_yzzz_yy, tr_yzzz_yz, tr_z_0_y_yzz_x, tr_z_0_y_yzz_y, tr_z_0_y_yzz_z, tr_z_x, tr_z_y, tr_z_z, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yzz_x[i] = 2.0 * tr_z_x[i] - 2.0 * tr_zzz_x[i] * tbe_0 - 4.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_y[i] = 2.0 * tr_z_y[i] - 2.0 * tr_zzz_y[i] * tbe_0 + 2.0 * tr_yz_0[i] - 4.0 * tr_yz_yy[i] * tke_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_yy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_yyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_z[i] = 2.0 * tr_z_z[i] - 2.0 * tr_zzz_z[i] * tbe_0 - 4.0 * tr_yz_yz[i] * tke_0 + 4.0 * tr_yzzz_yz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_z[i] * tbe_0 + 4.0 * tr_yyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 237-240 components of targeted buffer : FP

    auto tr_z_0_y_zzz_x = pbuffer.data(idx_op_geom_110_fp + 237);

    auto tr_z_0_y_zzz_y = pbuffer.data(idx_op_geom_110_fp + 238);

    auto tr_z_0_y_zzz_z = pbuffer.data(idx_op_geom_110_fp + 239);

    #pragma omp simd aligned(tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_yzzzz_x, tr_yzzzz_y, tr_yzzzz_z, tr_z_0_y_zzz_x, tr_z_0_y_zzz_y, tr_z_0_y_zzz_z, tr_zz_0, tr_zz_xy, tr_zz_yy, tr_zz_yz, tr_zzzz_0, tr_zzzz_xy, tr_zzzz_yy, tr_zzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_zzz_x[i] = -6.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zzzz_xy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_yzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_y[i] = 3.0 * tr_zz_0[i] - 6.0 * tr_zz_yy[i] * tke_0 - 2.0 * tr_zzzz_0[i] * tbe_0 + 4.0 * tr_zzzz_yy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_yzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_z[i] = -6.0 * tr_zz_yz[i] * tke_0 + 4.0 * tr_zzzz_yz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_yzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 240-243 components of targeted buffer : FP

    auto tr_z_0_z_xxx_x = pbuffer.data(idx_op_geom_110_fp + 240);

    auto tr_z_0_z_xxx_y = pbuffer.data(idx_op_geom_110_fp + 241);

    auto tr_z_0_z_xxx_z = pbuffer.data(idx_op_geom_110_fp + 242);

    #pragma omp simd aligned(tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxz_0, tr_xxxz_xz, tr_xxxz_yz, tr_xxxz_zz, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, tr_z_0_z_xxx_x, tr_z_0_z_xxx_y, tr_z_0_z_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxx_x[i] = -2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_y[i] = -2.0 * tr_xxx_y[i] * tbe_0 + 4.0 * tr_xxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_z[i] = -2.0 * tr_xxx_z[i] * tbe_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 243-246 components of targeted buffer : FP

    auto tr_z_0_z_xxy_x = pbuffer.data(idx_op_geom_110_fp + 243);

    auto tr_z_0_z_xxy_y = pbuffer.data(idx_op_geom_110_fp + 244);

    auto tr_z_0_z_xxy_z = pbuffer.data(idx_op_geom_110_fp + 245);

    #pragma omp simd aligned(tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyz_0, tr_xxyz_xz, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_z_0_z_xxy_x, tr_z_0_z_xxy_y, tr_z_0_z_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxy_x[i] = -2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_y[i] = -2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_z[i] = -2.0 * tr_xxy_z[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 246-249 components of targeted buffer : FP

    auto tr_z_0_z_xxz_x = pbuffer.data(idx_op_geom_110_fp + 246);

    auto tr_z_0_z_xxz_y = pbuffer.data(idx_op_geom_110_fp + 247);

    auto tr_z_0_z_xxz_z = pbuffer.data(idx_op_geom_110_fp + 248);

    #pragma omp simd aligned(tr_xx_0, tr_xx_xz, tr_xx_yz, tr_xx_zz, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xxzz_0, tr_xxzz_xz, tr_xxzz_yz, tr_xxzz_zz, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, tr_z_0_z_xxz_x, tr_z_0_z_xxz_y, tr_z_0_z_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxz_x[i] = -2.0 * tr_xx_xz[i] * tke_0 - 6.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_y[i] = -2.0 * tr_xx_yz[i] * tke_0 - 6.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_z[i] = tr_xx_0[i] - 2.0 * tr_xx_zz[i] * tke_0 - 6.0 * tr_xxz_z[i] * tbe_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 249-252 components of targeted buffer : FP

    auto tr_z_0_z_xyy_x = pbuffer.data(idx_op_geom_110_fp + 249);

    auto tr_z_0_z_xyy_y = pbuffer.data(idx_op_geom_110_fp + 250);

    auto tr_z_0_z_xyy_z = pbuffer.data(idx_op_geom_110_fp + 251);

    #pragma omp simd aligned(tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyz_0, tr_xyyz_xz, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_z_0_z_xyy_x, tr_z_0_z_xyy_y, tr_z_0_z_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyy_x[i] = -2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_y[i] = -2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_z[i] = -2.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 252-255 components of targeted buffer : FP

    auto tr_z_0_z_xyz_x = pbuffer.data(idx_op_geom_110_fp + 252);

    auto tr_z_0_z_xyz_y = pbuffer.data(idx_op_geom_110_fp + 253);

    auto tr_z_0_z_xyz_z = pbuffer.data(idx_op_geom_110_fp + 254);

    #pragma omp simd aligned(tr_xy_0, tr_xy_xz, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xyzz_xz, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_z_0_z_xyz_x, tr_z_0_z_xyz_y, tr_z_0_z_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyz_x[i] = -2.0 * tr_xy_xz[i] * tke_0 - 6.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_y[i] = -2.0 * tr_xy_yz[i] * tke_0 - 6.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_z[i] = tr_xy_0[i] - 2.0 * tr_xy_zz[i] * tke_0 - 6.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 255-258 components of targeted buffer : FP

    auto tr_z_0_z_xzz_x = pbuffer.data(idx_op_geom_110_fp + 255);

    auto tr_z_0_z_xzz_y = pbuffer.data(idx_op_geom_110_fp + 256);

    auto tr_z_0_z_xzz_z = pbuffer.data(idx_op_geom_110_fp + 257);

    #pragma omp simd aligned(tr_x_x, tr_x_y, tr_x_z, tr_xz_0, tr_xz_xz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_xzzz_0, tr_xzzz_xz, tr_xzzz_yz, tr_xzzz_zz, tr_xzzzz_x, tr_xzzzz_y, tr_xzzzz_z, tr_z_0_z_xzz_x, tr_z_0_z_xzz_y, tr_z_0_z_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xzz_x[i] = 2.0 * tr_x_x[i] - 4.0 * tr_xz_xz[i] * tke_0 - 10.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_y[i] = 2.0 * tr_x_y[i] - 4.0 * tr_xz_yz[i] * tke_0 - 10.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_z[i] = 2.0 * tr_x_z[i] + 2.0 * tr_xz_0[i] - 4.0 * tr_xz_zz[i] * tke_0 - 10.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 258-261 components of targeted buffer : FP

    auto tr_z_0_z_yyy_x = pbuffer.data(idx_op_geom_110_fp + 258);

    auto tr_z_0_z_yyy_y = pbuffer.data(idx_op_geom_110_fp + 259);

    auto tr_z_0_z_yyy_z = pbuffer.data(idx_op_geom_110_fp + 260);

    #pragma omp simd aligned(tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyz_0, tr_yyyz_xz, tr_yyyz_yz, tr_yyyz_zz, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, tr_z_0_z_yyy_x, tr_z_0_z_yyy_y, tr_z_0_z_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyy_x[i] = -2.0 * tr_yyy_x[i] * tbe_0 + 4.0 * tr_yyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_y[i] = -2.0 * tr_yyy_y[i] * tbe_0 + 4.0 * tr_yyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_z[i] = -2.0 * tr_yyy_z[i] * tbe_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 261-264 components of targeted buffer : FP

    auto tr_z_0_z_yyz_x = pbuffer.data(idx_op_geom_110_fp + 261);

    auto tr_z_0_z_yyz_y = pbuffer.data(idx_op_geom_110_fp + 262);

    auto tr_z_0_z_yyz_z = pbuffer.data(idx_op_geom_110_fp + 263);

    #pragma omp simd aligned(tr_yy_0, tr_yy_xz, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yyzz_0, tr_yyzz_xz, tr_yyzz_yz, tr_yyzz_zz, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, tr_z_0_z_yyz_x, tr_z_0_z_yyz_y, tr_z_0_z_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyz_x[i] = -2.0 * tr_yy_xz[i] * tke_0 - 6.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_y[i] = -2.0 * tr_yy_yz[i] * tke_0 - 6.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_yyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_z[i] = tr_yy_0[i] - 2.0 * tr_yy_zz[i] * tke_0 - 6.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 264-267 components of targeted buffer : FP

    auto tr_z_0_z_yzz_x = pbuffer.data(idx_op_geom_110_fp + 264);

    auto tr_z_0_z_yzz_y = pbuffer.data(idx_op_geom_110_fp + 265);

    auto tr_z_0_z_yzz_z = pbuffer.data(idx_op_geom_110_fp + 266);

    #pragma omp simd aligned(tr_y_x, tr_y_y, tr_y_z, tr_yz_0, tr_yz_xz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_yzzz_0, tr_yzzz_xz, tr_yzzz_yz, tr_yzzz_zz, tr_yzzzz_x, tr_yzzzz_y, tr_yzzzz_z, tr_z_0_z_yzz_x, tr_z_0_z_yzz_y, tr_z_0_z_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yzz_x[i] = 2.0 * tr_y_x[i] - 4.0 * tr_yz_xz[i] * tke_0 - 10.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_y[i] = 2.0 * tr_y_y[i] - 4.0 * tr_yz_yz[i] * tke_0 - 10.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_yzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_z[i] = 2.0 * tr_y_z[i] + 2.0 * tr_yz_0[i] - 4.0 * tr_yz_zz[i] * tke_0 - 10.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 267-270 components of targeted buffer : FP

    auto tr_z_0_z_zzz_x = pbuffer.data(idx_op_geom_110_fp + 267);

    auto tr_z_0_z_zzz_y = pbuffer.data(idx_op_geom_110_fp + 268);

    auto tr_z_0_z_zzz_z = pbuffer.data(idx_op_geom_110_fp + 269);

    #pragma omp simd aligned(tr_z_0_z_zzz_x, tr_z_0_z_zzz_y, tr_z_0_z_zzz_z, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, tr_zz_xz, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_y, tr_zzz_z, tr_zzzz_0, tr_zzzz_xz, tr_zzzz_yz, tr_zzzz_zz, tr_zzzzz_x, tr_zzzzz_y, tr_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_zzz_x[i] = 6.0 * tr_z_x[i] - 6.0 * tr_zz_xz[i] * tke_0 - 14.0 * tr_zzz_x[i] * tbe_0 + 4.0 * tr_zzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_y[i] = 6.0 * tr_z_y[i] - 6.0 * tr_zz_yz[i] * tke_0 - 14.0 * tr_zzz_y[i] * tbe_0 + 4.0 * tr_zzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_z[i] = 6.0 * tr_z_z[i] + 3.0 * tr_zz_0[i] - 6.0 * tr_zz_zz[i] * tke_0 - 14.0 * tr_zzz_z[i] * tbe_0 - 2.0 * tr_zzzz_0[i] * tbe_0 + 4.0 * tr_zzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_z[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

