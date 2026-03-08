#include "GeometricalDerivatives010ForGP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_gp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_gp,
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

    // Set up 0-3 components of targeted buffer : GP

    auto tr_0_0_x_xxxx_x = pbuffer.data(idx_op_geom_010_gp);

    auto tr_0_0_x_xxxx_y = pbuffer.data(idx_op_geom_010_gp + 1);

    auto tr_0_0_x_xxxx_z = pbuffer.data(idx_op_geom_010_gp + 2);

    #pragma omp simd aligned(tr_0_0_x_xxxx_x, tr_0_0_x_xxxx_y, tr_0_0_x_xxxx_z, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxx_0, tr_xxxx_xx, tr_xxxx_xy, tr_xxxx_xz, tr_xxxxx_x, tr_xxxxx_y, tr_xxxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxx_x[i] = 2.0 * tr_xxxxx_x[i] * tbe_0 + 2.0 * tr_xxxx_xx[i] * tke_0 - 4.0 * tr_xxx_x[i] - tr_xxxx_0[i];

        tr_0_0_x_xxxx_y[i] = 2.0 * tr_xxxxx_y[i] * tbe_0 + 2.0 * tr_xxxx_xy[i] * tke_0 - 4.0 * tr_xxx_y[i];

        tr_0_0_x_xxxx_z[i] = 2.0 * tr_xxxxx_z[i] * tbe_0 + 2.0 * tr_xxxx_xz[i] * tke_0 - 4.0 * tr_xxx_z[i];
    }

    // Set up 3-6 components of targeted buffer : GP

    auto tr_0_0_x_xxxy_x = pbuffer.data(idx_op_geom_010_gp + 3);

    auto tr_0_0_x_xxxy_y = pbuffer.data(idx_op_geom_010_gp + 4);

    auto tr_0_0_x_xxxy_z = pbuffer.data(idx_op_geom_010_gp + 5);

    #pragma omp simd aligned(tr_0_0_x_xxxy_x, tr_0_0_x_xxxy_y, tr_0_0_x_xxxy_z, tr_xxxxy_x, tr_xxxxy_y, tr_xxxxy_z, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxy_x, tr_xxy_y, tr_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxy_x[i] = 2.0 * tr_xxxxy_x[i] * tbe_0 + 2.0 * tr_xxxy_xx[i] * tke_0 - 3.0 * tr_xxy_x[i] - tr_xxxy_0[i];

        tr_0_0_x_xxxy_y[i] = 2.0 * tr_xxxxy_y[i] * tbe_0 + 2.0 * tr_xxxy_xy[i] * tke_0 - 3.0 * tr_xxy_y[i];

        tr_0_0_x_xxxy_z[i] = 2.0 * tr_xxxxy_z[i] * tbe_0 + 2.0 * tr_xxxy_xz[i] * tke_0 - 3.0 * tr_xxy_z[i];
    }

    // Set up 6-9 components of targeted buffer : GP

    auto tr_0_0_x_xxxz_x = pbuffer.data(idx_op_geom_010_gp + 6);

    auto tr_0_0_x_xxxz_y = pbuffer.data(idx_op_geom_010_gp + 7);

    auto tr_0_0_x_xxxz_z = pbuffer.data(idx_op_geom_010_gp + 8);

    #pragma omp simd aligned(tr_0_0_x_xxxz_x, tr_0_0_x_xxxz_y, tr_0_0_x_xxxz_z, tr_xxxxz_x, tr_xxxxz_y, tr_xxxxz_z, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxz_x, tr_xxz_y, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxz_x[i] = 2.0 * tr_xxxxz_x[i] * tbe_0 + 2.0 * tr_xxxz_xx[i] * tke_0 - 3.0 * tr_xxz_x[i] - tr_xxxz_0[i];

        tr_0_0_x_xxxz_y[i] = 2.0 * tr_xxxxz_y[i] * tbe_0 + 2.0 * tr_xxxz_xy[i] * tke_0 - 3.0 * tr_xxz_y[i];

        tr_0_0_x_xxxz_z[i] = 2.0 * tr_xxxxz_z[i] * tbe_0 + 2.0 * tr_xxxz_xz[i] * tke_0 - 3.0 * tr_xxz_z[i];
    }

    // Set up 9-12 components of targeted buffer : GP

    auto tr_0_0_x_xxyy_x = pbuffer.data(idx_op_geom_010_gp + 9);

    auto tr_0_0_x_xxyy_y = pbuffer.data(idx_op_geom_010_gp + 10);

    auto tr_0_0_x_xxyy_z = pbuffer.data(idx_op_geom_010_gp + 11);

    #pragma omp simd aligned(tr_0_0_x_xxyy_x, tr_0_0_x_xxyy_y, tr_0_0_x_xxyy_z, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xyy_x, tr_xyy_y, tr_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxyy_x[i] = 2.0 * tr_xxxyy_x[i] * tbe_0 + 2.0 * tr_xxyy_xx[i] * tke_0 - 2.0 * tr_xyy_x[i] - tr_xxyy_0[i];

        tr_0_0_x_xxyy_y[i] = 2.0 * tr_xxxyy_y[i] * tbe_0 + 2.0 * tr_xxyy_xy[i] * tke_0 - 2.0 * tr_xyy_y[i];

        tr_0_0_x_xxyy_z[i] = 2.0 * tr_xxxyy_z[i] * tbe_0 + 2.0 * tr_xxyy_xz[i] * tke_0 - 2.0 * tr_xyy_z[i];
    }

    // Set up 12-15 components of targeted buffer : GP

    auto tr_0_0_x_xxyz_x = pbuffer.data(idx_op_geom_010_gp + 12);

    auto tr_0_0_x_xxyz_y = pbuffer.data(idx_op_geom_010_gp + 13);

    auto tr_0_0_x_xxyz_z = pbuffer.data(idx_op_geom_010_gp + 14);

    #pragma omp simd aligned(tr_0_0_x_xxyz_x, tr_0_0_x_xxyz_y, tr_0_0_x_xxyz_z, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xyz_x, tr_xyz_y, tr_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxyz_x[i] = 2.0 * tr_xxxyz_x[i] * tbe_0 + 2.0 * tr_xxyz_xx[i] * tke_0 - 2.0 * tr_xyz_x[i] - tr_xxyz_0[i];

        tr_0_0_x_xxyz_y[i] = 2.0 * tr_xxxyz_y[i] * tbe_0 + 2.0 * tr_xxyz_xy[i] * tke_0 - 2.0 * tr_xyz_y[i];

        tr_0_0_x_xxyz_z[i] = 2.0 * tr_xxxyz_z[i] * tbe_0 + 2.0 * tr_xxyz_xz[i] * tke_0 - 2.0 * tr_xyz_z[i];
    }

    // Set up 15-18 components of targeted buffer : GP

    auto tr_0_0_x_xxzz_x = pbuffer.data(idx_op_geom_010_gp + 15);

    auto tr_0_0_x_xxzz_y = pbuffer.data(idx_op_geom_010_gp + 16);

    auto tr_0_0_x_xxzz_z = pbuffer.data(idx_op_geom_010_gp + 17);

    #pragma omp simd aligned(tr_0_0_x_xxzz_x, tr_0_0_x_xxzz_y, tr_0_0_x_xxzz_z, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xzz_x, tr_xzz_y, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxzz_x[i] = 2.0 * tr_xxxzz_x[i] * tbe_0 + 2.0 * tr_xxzz_xx[i] * tke_0 - 2.0 * tr_xzz_x[i] - tr_xxzz_0[i];

        tr_0_0_x_xxzz_y[i] = 2.0 * tr_xxxzz_y[i] * tbe_0 + 2.0 * tr_xxzz_xy[i] * tke_0 - 2.0 * tr_xzz_y[i];

        tr_0_0_x_xxzz_z[i] = 2.0 * tr_xxxzz_z[i] * tbe_0 + 2.0 * tr_xxzz_xz[i] * tke_0 - 2.0 * tr_xzz_z[i];
    }

    // Set up 18-21 components of targeted buffer : GP

    auto tr_0_0_x_xyyy_x = pbuffer.data(idx_op_geom_010_gp + 18);

    auto tr_0_0_x_xyyy_y = pbuffer.data(idx_op_geom_010_gp + 19);

    auto tr_0_0_x_xyyy_z = pbuffer.data(idx_op_geom_010_gp + 20);

    #pragma omp simd aligned(tr_0_0_x_xyyy_x, tr_0_0_x_xyyy_y, tr_0_0_x_xyyy_z, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_yyy_x, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyyy_x[i] = 2.0 * tr_xxyyy_x[i] * tbe_0 + 2.0 * tr_xyyy_xx[i] * tke_0 - tr_yyy_x[i] - tr_xyyy_0[i];

        tr_0_0_x_xyyy_y[i] = 2.0 * tr_xxyyy_y[i] * tbe_0 + 2.0 * tr_xyyy_xy[i] * tke_0 - tr_yyy_y[i];

        tr_0_0_x_xyyy_z[i] = 2.0 * tr_xxyyy_z[i] * tbe_0 + 2.0 * tr_xyyy_xz[i] * tke_0 - tr_yyy_z[i];
    }

    // Set up 21-24 components of targeted buffer : GP

    auto tr_0_0_x_xyyz_x = pbuffer.data(idx_op_geom_010_gp + 21);

    auto tr_0_0_x_xyyz_y = pbuffer.data(idx_op_geom_010_gp + 22);

    auto tr_0_0_x_xyyz_z = pbuffer.data(idx_op_geom_010_gp + 23);

    #pragma omp simd aligned(tr_0_0_x_xyyz_x, tr_0_0_x_xyyz_y, tr_0_0_x_xyyz_z, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_yyz_x, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyyz_x[i] = 2.0 * tr_xxyyz_x[i] * tbe_0 + 2.0 * tr_xyyz_xx[i] * tke_0 - tr_yyz_x[i] - tr_xyyz_0[i];

        tr_0_0_x_xyyz_y[i] = 2.0 * tr_xxyyz_y[i] * tbe_0 + 2.0 * tr_xyyz_xy[i] * tke_0 - tr_yyz_y[i];

        tr_0_0_x_xyyz_z[i] = 2.0 * tr_xxyyz_z[i] * tbe_0 + 2.0 * tr_xyyz_xz[i] * tke_0 - tr_yyz_z[i];
    }

    // Set up 24-27 components of targeted buffer : GP

    auto tr_0_0_x_xyzz_x = pbuffer.data(idx_op_geom_010_gp + 24);

    auto tr_0_0_x_xyzz_y = pbuffer.data(idx_op_geom_010_gp + 25);

    auto tr_0_0_x_xyzz_z = pbuffer.data(idx_op_geom_010_gp + 26);

    #pragma omp simd aligned(tr_0_0_x_xyzz_x, tr_0_0_x_xyzz_y, tr_0_0_x_xyzz_z, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyzz_x[i] = 2.0 * tr_xxyzz_x[i] * tbe_0 + 2.0 * tr_xyzz_xx[i] * tke_0 - tr_yzz_x[i] - tr_xyzz_0[i];

        tr_0_0_x_xyzz_y[i] = 2.0 * tr_xxyzz_y[i] * tbe_0 + 2.0 * tr_xyzz_xy[i] * tke_0 - tr_yzz_y[i];

        tr_0_0_x_xyzz_z[i] = 2.0 * tr_xxyzz_z[i] * tbe_0 + 2.0 * tr_xyzz_xz[i] * tke_0 - tr_yzz_z[i];
    }

    // Set up 27-30 components of targeted buffer : GP

    auto tr_0_0_x_xzzz_x = pbuffer.data(idx_op_geom_010_gp + 27);

    auto tr_0_0_x_xzzz_y = pbuffer.data(idx_op_geom_010_gp + 28);

    auto tr_0_0_x_xzzz_z = pbuffer.data(idx_op_geom_010_gp + 29);

    #pragma omp simd aligned(tr_0_0_x_xzzz_x, tr_0_0_x_xzzz_y, tr_0_0_x_xzzz_z, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xzzz_x[i] = 2.0 * tr_xxzzz_x[i] * tbe_0 + 2.0 * tr_xzzz_xx[i] * tke_0 - tr_zzz_x[i] - tr_xzzz_0[i];

        tr_0_0_x_xzzz_y[i] = 2.0 * tr_xxzzz_y[i] * tbe_0 + 2.0 * tr_xzzz_xy[i] * tke_0 - tr_zzz_y[i];

        tr_0_0_x_xzzz_z[i] = 2.0 * tr_xxzzz_z[i] * tbe_0 + 2.0 * tr_xzzz_xz[i] * tke_0 - tr_zzz_z[i];
    }

    // Set up 30-33 components of targeted buffer : GP

    auto tr_0_0_x_yyyy_x = pbuffer.data(idx_op_geom_010_gp + 30);

    auto tr_0_0_x_yyyy_y = pbuffer.data(idx_op_geom_010_gp + 31);

    auto tr_0_0_x_yyyy_z = pbuffer.data(idx_op_geom_010_gp + 32);

    #pragma omp simd aligned(tr_0_0_x_yyyy_x, tr_0_0_x_yyyy_y, tr_0_0_x_yyyy_z, tr_xyyyy_x, tr_xyyyy_y, tr_xyyyy_z, tr_yyyy_0, tr_yyyy_xx, tr_yyyy_xy, tr_yyyy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyyy_x[i] = 2.0 * tr_xyyyy_x[i] * tbe_0 + 2.0 * tr_yyyy_xx[i] * tke_0 - tr_yyyy_0[i];

        tr_0_0_x_yyyy_y[i] = 2.0 * tr_xyyyy_y[i] * tbe_0 + 2.0 * tr_yyyy_xy[i] * tke_0;

        tr_0_0_x_yyyy_z[i] = 2.0 * tr_xyyyy_z[i] * tbe_0 + 2.0 * tr_yyyy_xz[i] * tke_0;
    }

    // Set up 33-36 components of targeted buffer : GP

    auto tr_0_0_x_yyyz_x = pbuffer.data(idx_op_geom_010_gp + 33);

    auto tr_0_0_x_yyyz_y = pbuffer.data(idx_op_geom_010_gp + 34);

    auto tr_0_0_x_yyyz_z = pbuffer.data(idx_op_geom_010_gp + 35);

    #pragma omp simd aligned(tr_0_0_x_yyyz_x, tr_0_0_x_yyyz_y, tr_0_0_x_yyyz_z, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyyz_x[i] = 2.0 * tr_xyyyz_x[i] * tbe_0 + 2.0 * tr_yyyz_xx[i] * tke_0 - tr_yyyz_0[i];

        tr_0_0_x_yyyz_y[i] = 2.0 * tr_xyyyz_y[i] * tbe_0 + 2.0 * tr_yyyz_xy[i] * tke_0;

        tr_0_0_x_yyyz_z[i] = 2.0 * tr_xyyyz_z[i] * tbe_0 + 2.0 * tr_yyyz_xz[i] * tke_0;
    }

    // Set up 36-39 components of targeted buffer : GP

    auto tr_0_0_x_yyzz_x = pbuffer.data(idx_op_geom_010_gp + 36);

    auto tr_0_0_x_yyzz_y = pbuffer.data(idx_op_geom_010_gp + 37);

    auto tr_0_0_x_yyzz_z = pbuffer.data(idx_op_geom_010_gp + 38);

    #pragma omp simd aligned(tr_0_0_x_yyzz_x, tr_0_0_x_yyzz_y, tr_0_0_x_yyzz_z, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyzz_x[i] = 2.0 * tr_xyyzz_x[i] * tbe_0 + 2.0 * tr_yyzz_xx[i] * tke_0 - tr_yyzz_0[i];

        tr_0_0_x_yyzz_y[i] = 2.0 * tr_xyyzz_y[i] * tbe_0 + 2.0 * tr_yyzz_xy[i] * tke_0;

        tr_0_0_x_yyzz_z[i] = 2.0 * tr_xyyzz_z[i] * tbe_0 + 2.0 * tr_yyzz_xz[i] * tke_0;
    }

    // Set up 39-42 components of targeted buffer : GP

    auto tr_0_0_x_yzzz_x = pbuffer.data(idx_op_geom_010_gp + 39);

    auto tr_0_0_x_yzzz_y = pbuffer.data(idx_op_geom_010_gp + 40);

    auto tr_0_0_x_yzzz_z = pbuffer.data(idx_op_geom_010_gp + 41);

    #pragma omp simd aligned(tr_0_0_x_yzzz_x, tr_0_0_x_yzzz_y, tr_0_0_x_yzzz_z, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yzzz_x[i] = 2.0 * tr_xyzzz_x[i] * tbe_0 + 2.0 * tr_yzzz_xx[i] * tke_0 - tr_yzzz_0[i];

        tr_0_0_x_yzzz_y[i] = 2.0 * tr_xyzzz_y[i] * tbe_0 + 2.0 * tr_yzzz_xy[i] * tke_0;

        tr_0_0_x_yzzz_z[i] = 2.0 * tr_xyzzz_z[i] * tbe_0 + 2.0 * tr_yzzz_xz[i] * tke_0;
    }

    // Set up 42-45 components of targeted buffer : GP

    auto tr_0_0_x_zzzz_x = pbuffer.data(idx_op_geom_010_gp + 42);

    auto tr_0_0_x_zzzz_y = pbuffer.data(idx_op_geom_010_gp + 43);

    auto tr_0_0_x_zzzz_z = pbuffer.data(idx_op_geom_010_gp + 44);

    #pragma omp simd aligned(tr_0_0_x_zzzz_x, tr_0_0_x_zzzz_y, tr_0_0_x_zzzz_z, tr_xzzzz_x, tr_xzzzz_y, tr_xzzzz_z, tr_zzzz_0, tr_zzzz_xx, tr_zzzz_xy, tr_zzzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zzzz_x[i] = 2.0 * tr_xzzzz_x[i] * tbe_0 + 2.0 * tr_zzzz_xx[i] * tke_0 - tr_zzzz_0[i];

        tr_0_0_x_zzzz_y[i] = 2.0 * tr_xzzzz_y[i] * tbe_0 + 2.0 * tr_zzzz_xy[i] * tke_0;

        tr_0_0_x_zzzz_z[i] = 2.0 * tr_xzzzz_z[i] * tbe_0 + 2.0 * tr_zzzz_xz[i] * tke_0;
    }

    // Set up 45-48 components of targeted buffer : GP

    auto tr_0_0_y_xxxx_x = pbuffer.data(idx_op_geom_010_gp + 45);

    auto tr_0_0_y_xxxx_y = pbuffer.data(idx_op_geom_010_gp + 46);

    auto tr_0_0_y_xxxx_z = pbuffer.data(idx_op_geom_010_gp + 47);

    #pragma omp simd aligned(tr_0_0_y_xxxx_x, tr_0_0_y_xxxx_y, tr_0_0_y_xxxx_z, tr_xxxx_0, tr_xxxx_xy, tr_xxxx_yy, tr_xxxx_yz, tr_xxxxy_x, tr_xxxxy_y, tr_xxxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxx_x[i] = 2.0 * tr_xxxxy_x[i] * tbe_0 + 2.0 * tr_xxxx_xy[i] * tke_0;

        tr_0_0_y_xxxx_y[i] = 2.0 * tr_xxxxy_y[i] * tbe_0 + 2.0 * tr_xxxx_yy[i] * tke_0 - tr_xxxx_0[i];

        tr_0_0_y_xxxx_z[i] = 2.0 * tr_xxxxy_z[i] * tbe_0 + 2.0 * tr_xxxx_yz[i] * tke_0;
    }

    // Set up 48-51 components of targeted buffer : GP

    auto tr_0_0_y_xxxy_x = pbuffer.data(idx_op_geom_010_gp + 48);

    auto tr_0_0_y_xxxy_y = pbuffer.data(idx_op_geom_010_gp + 49);

    auto tr_0_0_y_xxxy_z = pbuffer.data(idx_op_geom_010_gp + 50);

    #pragma omp simd aligned(tr_0_0_y_xxxy_x, tr_0_0_y_xxxy_y, tr_0_0_y_xxxy_z, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxy_0, tr_xxxy_xy, tr_xxxy_yy, tr_xxxy_yz, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxy_x[i] = 2.0 * tr_xxxyy_x[i] * tbe_0 + 2.0 * tr_xxxy_xy[i] * tke_0 - tr_xxx_x[i];

        tr_0_0_y_xxxy_y[i] = 2.0 * tr_xxxyy_y[i] * tbe_0 + 2.0 * tr_xxxy_yy[i] * tke_0 - tr_xxx_y[i] - tr_xxxy_0[i];

        tr_0_0_y_xxxy_z[i] = 2.0 * tr_xxxyy_z[i] * tbe_0 + 2.0 * tr_xxxy_yz[i] * tke_0 - tr_xxx_z[i];
    }

    // Set up 51-54 components of targeted buffer : GP

    auto tr_0_0_y_xxxz_x = pbuffer.data(idx_op_geom_010_gp + 51);

    auto tr_0_0_y_xxxz_y = pbuffer.data(idx_op_geom_010_gp + 52);

    auto tr_0_0_y_xxxz_z = pbuffer.data(idx_op_geom_010_gp + 53);

    #pragma omp simd aligned(tr_0_0_y_xxxz_x, tr_0_0_y_xxxz_y, tr_0_0_y_xxxz_z, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxxz_0, tr_xxxz_xy, tr_xxxz_yy, tr_xxxz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxz_x[i] = 2.0 * tr_xxxyz_x[i] * tbe_0 + 2.0 * tr_xxxz_xy[i] * tke_0;

        tr_0_0_y_xxxz_y[i] = 2.0 * tr_xxxyz_y[i] * tbe_0 + 2.0 * tr_xxxz_yy[i] * tke_0 - tr_xxxz_0[i];

        tr_0_0_y_xxxz_z[i] = 2.0 * tr_xxxyz_z[i] * tbe_0 + 2.0 * tr_xxxz_yz[i] * tke_0;
    }

    // Set up 54-57 components of targeted buffer : GP

    auto tr_0_0_y_xxyy_x = pbuffer.data(idx_op_geom_010_gp + 54);

    auto tr_0_0_y_xxyy_y = pbuffer.data(idx_op_geom_010_gp + 55);

    auto tr_0_0_y_xxyy_z = pbuffer.data(idx_op_geom_010_gp + 56);

    #pragma omp simd aligned(tr_0_0_y_xxyy_x, tr_0_0_y_xxyy_y, tr_0_0_y_xxyy_z, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyy_0, tr_xxyy_xy, tr_xxyy_yy, tr_xxyy_yz, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxyy_x[i] = 2.0 * tr_xxyyy_x[i] * tbe_0 + 2.0 * tr_xxyy_xy[i] * tke_0 - 2.0 * tr_xxy_x[i];

        tr_0_0_y_xxyy_y[i] = 2.0 * tr_xxyyy_y[i] * tbe_0 + 2.0 * tr_xxyy_yy[i] * tke_0 - 2.0 * tr_xxy_y[i] - tr_xxyy_0[i];

        tr_0_0_y_xxyy_z[i] = 2.0 * tr_xxyyy_z[i] * tbe_0 + 2.0 * tr_xxyy_yz[i] * tke_0 - 2.0 * tr_xxy_z[i];
    }

    // Set up 57-60 components of targeted buffer : GP

    auto tr_0_0_y_xxyz_x = pbuffer.data(idx_op_geom_010_gp + 57);

    auto tr_0_0_y_xxyz_y = pbuffer.data(idx_op_geom_010_gp + 58);

    auto tr_0_0_y_xxyz_z = pbuffer.data(idx_op_geom_010_gp + 59);

    #pragma omp simd aligned(tr_0_0_y_xxyz_x, tr_0_0_y_xxyz_y, tr_0_0_y_xxyz_z, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxyz_0, tr_xxyz_xy, tr_xxyz_yy, tr_xxyz_yz, tr_xxz_x, tr_xxz_y, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxyz_x[i] = 2.0 * tr_xxyyz_x[i] * tbe_0 + 2.0 * tr_xxyz_xy[i] * tke_0 - tr_xxz_x[i];

        tr_0_0_y_xxyz_y[i] = 2.0 * tr_xxyyz_y[i] * tbe_0 + 2.0 * tr_xxyz_yy[i] * tke_0 - tr_xxz_y[i] - tr_xxyz_0[i];

        tr_0_0_y_xxyz_z[i] = 2.0 * tr_xxyyz_z[i] * tbe_0 + 2.0 * tr_xxyz_yz[i] * tke_0 - tr_xxz_z[i];
    }

    // Set up 60-63 components of targeted buffer : GP

    auto tr_0_0_y_xxzz_x = pbuffer.data(idx_op_geom_010_gp + 60);

    auto tr_0_0_y_xxzz_y = pbuffer.data(idx_op_geom_010_gp + 61);

    auto tr_0_0_y_xxzz_z = pbuffer.data(idx_op_geom_010_gp + 62);

    #pragma omp simd aligned(tr_0_0_y_xxzz_x, tr_0_0_y_xxzz_y, tr_0_0_y_xxzz_z, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xxzz_0, tr_xxzz_xy, tr_xxzz_yy, tr_xxzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxzz_x[i] = 2.0 * tr_xxyzz_x[i] * tbe_0 + 2.0 * tr_xxzz_xy[i] * tke_0;

        tr_0_0_y_xxzz_y[i] = 2.0 * tr_xxyzz_y[i] * tbe_0 + 2.0 * tr_xxzz_yy[i] * tke_0 - tr_xxzz_0[i];

        tr_0_0_y_xxzz_z[i] = 2.0 * tr_xxyzz_z[i] * tbe_0 + 2.0 * tr_xxzz_yz[i] * tke_0;
    }

    // Set up 63-66 components of targeted buffer : GP

    auto tr_0_0_y_xyyy_x = pbuffer.data(idx_op_geom_010_gp + 63);

    auto tr_0_0_y_xyyy_y = pbuffer.data(idx_op_geom_010_gp + 64);

    auto tr_0_0_y_xyyy_z = pbuffer.data(idx_op_geom_010_gp + 65);

    #pragma omp simd aligned(tr_0_0_y_xyyy_x, tr_0_0_y_xyyy_y, tr_0_0_y_xyyy_z, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyy_0, tr_xyyy_xy, tr_xyyy_yy, tr_xyyy_yz, tr_xyyyy_x, tr_xyyyy_y, tr_xyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyyy_x[i] = 2.0 * tr_xyyyy_x[i] * tbe_0 + 2.0 * tr_xyyy_xy[i] * tke_0 - 3.0 * tr_xyy_x[i];

        tr_0_0_y_xyyy_y[i] = 2.0 * tr_xyyyy_y[i] * tbe_0 + 2.0 * tr_xyyy_yy[i] * tke_0 - 3.0 * tr_xyy_y[i] - tr_xyyy_0[i];

        tr_0_0_y_xyyy_z[i] = 2.0 * tr_xyyyy_z[i] * tbe_0 + 2.0 * tr_xyyy_yz[i] * tke_0 - 3.0 * tr_xyy_z[i];
    }

    // Set up 66-69 components of targeted buffer : GP

    auto tr_0_0_y_xyyz_x = pbuffer.data(idx_op_geom_010_gp + 66);

    auto tr_0_0_y_xyyz_y = pbuffer.data(idx_op_geom_010_gp + 67);

    auto tr_0_0_y_xyyz_z = pbuffer.data(idx_op_geom_010_gp + 68);

    #pragma omp simd aligned(tr_0_0_y_xyyz_x, tr_0_0_y_xyyz_y, tr_0_0_y_xyyz_z, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyyz_0, tr_xyyz_xy, tr_xyyz_yy, tr_xyyz_yz, tr_xyz_x, tr_xyz_y, tr_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyyz_x[i] = 2.0 * tr_xyyyz_x[i] * tbe_0 + 2.0 * tr_xyyz_xy[i] * tke_0 - 2.0 * tr_xyz_x[i];

        tr_0_0_y_xyyz_y[i] = 2.0 * tr_xyyyz_y[i] * tbe_0 + 2.0 * tr_xyyz_yy[i] * tke_0 - 2.0 * tr_xyz_y[i] - tr_xyyz_0[i];

        tr_0_0_y_xyyz_z[i] = 2.0 * tr_xyyyz_z[i] * tbe_0 + 2.0 * tr_xyyz_yz[i] * tke_0 - 2.0 * tr_xyz_z[i];
    }

    // Set up 69-72 components of targeted buffer : GP

    auto tr_0_0_y_xyzz_x = pbuffer.data(idx_op_geom_010_gp + 69);

    auto tr_0_0_y_xyzz_y = pbuffer.data(idx_op_geom_010_gp + 70);

    auto tr_0_0_y_xyzz_z = pbuffer.data(idx_op_geom_010_gp + 71);

    #pragma omp simd aligned(tr_0_0_y_xyzz_x, tr_0_0_y_xyzz_y, tr_0_0_y_xyzz_z, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xyzz_0, tr_xyzz_xy, tr_xyzz_yy, tr_xyzz_yz, tr_xzz_x, tr_xzz_y, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyzz_x[i] = 2.0 * tr_xyyzz_x[i] * tbe_0 + 2.0 * tr_xyzz_xy[i] * tke_0 - tr_xzz_x[i];

        tr_0_0_y_xyzz_y[i] = 2.0 * tr_xyyzz_y[i] * tbe_0 + 2.0 * tr_xyzz_yy[i] * tke_0 - tr_xzz_y[i] - tr_xyzz_0[i];

        tr_0_0_y_xyzz_z[i] = 2.0 * tr_xyyzz_z[i] * tbe_0 + 2.0 * tr_xyzz_yz[i] * tke_0 - tr_xzz_z[i];
    }

    // Set up 72-75 components of targeted buffer : GP

    auto tr_0_0_y_xzzz_x = pbuffer.data(idx_op_geom_010_gp + 72);

    auto tr_0_0_y_xzzz_y = pbuffer.data(idx_op_geom_010_gp + 73);

    auto tr_0_0_y_xzzz_z = pbuffer.data(idx_op_geom_010_gp + 74);

    #pragma omp simd aligned(tr_0_0_y_xzzz_x, tr_0_0_y_xzzz_y, tr_0_0_y_xzzz_z, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_xzzz_0, tr_xzzz_xy, tr_xzzz_yy, tr_xzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xzzz_x[i] = 2.0 * tr_xyzzz_x[i] * tbe_0 + 2.0 * tr_xzzz_xy[i] * tke_0;

        tr_0_0_y_xzzz_y[i] = 2.0 * tr_xyzzz_y[i] * tbe_0 + 2.0 * tr_xzzz_yy[i] * tke_0 - tr_xzzz_0[i];

        tr_0_0_y_xzzz_z[i] = 2.0 * tr_xyzzz_z[i] * tbe_0 + 2.0 * tr_xzzz_yz[i] * tke_0;
    }

    // Set up 75-78 components of targeted buffer : GP

    auto tr_0_0_y_yyyy_x = pbuffer.data(idx_op_geom_010_gp + 75);

    auto tr_0_0_y_yyyy_y = pbuffer.data(idx_op_geom_010_gp + 76);

    auto tr_0_0_y_yyyy_z = pbuffer.data(idx_op_geom_010_gp + 77);

    #pragma omp simd aligned(tr_0_0_y_yyyy_x, tr_0_0_y_yyyy_y, tr_0_0_y_yyyy_z, tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyy_0, tr_yyyy_xy, tr_yyyy_yy, tr_yyyy_yz, tr_yyyyy_x, tr_yyyyy_y, tr_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyyy_x[i] = 2.0 * tr_yyyyy_x[i] * tbe_0 + 2.0 * tr_yyyy_xy[i] * tke_0 - 4.0 * tr_yyy_x[i];

        tr_0_0_y_yyyy_y[i] = 2.0 * tr_yyyyy_y[i] * tbe_0 + 2.0 * tr_yyyy_yy[i] * tke_0 - 4.0 * tr_yyy_y[i] - tr_yyyy_0[i];

        tr_0_0_y_yyyy_z[i] = 2.0 * tr_yyyyy_z[i] * tbe_0 + 2.0 * tr_yyyy_yz[i] * tke_0 - 4.0 * tr_yyy_z[i];
    }

    // Set up 78-81 components of targeted buffer : GP

    auto tr_0_0_y_yyyz_x = pbuffer.data(idx_op_geom_010_gp + 78);

    auto tr_0_0_y_yyyz_y = pbuffer.data(idx_op_geom_010_gp + 79);

    auto tr_0_0_y_yyyz_z = pbuffer.data(idx_op_geom_010_gp + 80);

    #pragma omp simd aligned(tr_0_0_y_yyyz_x, tr_0_0_y_yyyz_y, tr_0_0_y_yyyz_z, tr_yyyyz_x, tr_yyyyz_y, tr_yyyyz_z, tr_yyyz_0, tr_yyyz_xy, tr_yyyz_yy, tr_yyyz_yz, tr_yyz_x, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyyz_x[i] = 2.0 * tr_yyyyz_x[i] * tbe_0 + 2.0 * tr_yyyz_xy[i] * tke_0 - 3.0 * tr_yyz_x[i];

        tr_0_0_y_yyyz_y[i] = 2.0 * tr_yyyyz_y[i] * tbe_0 + 2.0 * tr_yyyz_yy[i] * tke_0 - 3.0 * tr_yyz_y[i] - tr_yyyz_0[i];

        tr_0_0_y_yyyz_z[i] = 2.0 * tr_yyyyz_z[i] * tbe_0 + 2.0 * tr_yyyz_yz[i] * tke_0 - 3.0 * tr_yyz_z[i];
    }

    // Set up 81-84 components of targeted buffer : GP

    auto tr_0_0_y_yyzz_x = pbuffer.data(idx_op_geom_010_gp + 81);

    auto tr_0_0_y_yyzz_y = pbuffer.data(idx_op_geom_010_gp + 82);

    auto tr_0_0_y_yyzz_z = pbuffer.data(idx_op_geom_010_gp + 83);

    #pragma omp simd aligned(tr_0_0_y_yyzz_x, tr_0_0_y_yyzz_y, tr_0_0_y_yyzz_z, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, tr_yyzz_0, tr_yyzz_xy, tr_yyzz_yy, tr_yyzz_yz, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyzz_x[i] = 2.0 * tr_yyyzz_x[i] * tbe_0 + 2.0 * tr_yyzz_xy[i] * tke_0 - 2.0 * tr_yzz_x[i];

        tr_0_0_y_yyzz_y[i] = 2.0 * tr_yyyzz_y[i] * tbe_0 + 2.0 * tr_yyzz_yy[i] * tke_0 - 2.0 * tr_yzz_y[i] - tr_yyzz_0[i];

        tr_0_0_y_yyzz_z[i] = 2.0 * tr_yyyzz_z[i] * tbe_0 + 2.0 * tr_yyzz_yz[i] * tke_0 - 2.0 * tr_yzz_z[i];
    }

    // Set up 84-87 components of targeted buffer : GP

    auto tr_0_0_y_yzzz_x = pbuffer.data(idx_op_geom_010_gp + 84);

    auto tr_0_0_y_yzzz_y = pbuffer.data(idx_op_geom_010_gp + 85);

    auto tr_0_0_y_yzzz_z = pbuffer.data(idx_op_geom_010_gp + 86);

    #pragma omp simd aligned(tr_0_0_y_yzzz_x, tr_0_0_y_yzzz_y, tr_0_0_y_yzzz_z, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, tr_yzzz_0, tr_yzzz_xy, tr_yzzz_yy, tr_yzzz_yz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yzzz_x[i] = 2.0 * tr_yyzzz_x[i] * tbe_0 + 2.0 * tr_yzzz_xy[i] * tke_0 - tr_zzz_x[i];

        tr_0_0_y_yzzz_y[i] = 2.0 * tr_yyzzz_y[i] * tbe_0 + 2.0 * tr_yzzz_yy[i] * tke_0 - tr_zzz_y[i] - tr_yzzz_0[i];

        tr_0_0_y_yzzz_z[i] = 2.0 * tr_yyzzz_z[i] * tbe_0 + 2.0 * tr_yzzz_yz[i] * tke_0 - tr_zzz_z[i];
    }

    // Set up 87-90 components of targeted buffer : GP

    auto tr_0_0_y_zzzz_x = pbuffer.data(idx_op_geom_010_gp + 87);

    auto tr_0_0_y_zzzz_y = pbuffer.data(idx_op_geom_010_gp + 88);

    auto tr_0_0_y_zzzz_z = pbuffer.data(idx_op_geom_010_gp + 89);

    #pragma omp simd aligned(tr_0_0_y_zzzz_x, tr_0_0_y_zzzz_y, tr_0_0_y_zzzz_z, tr_yzzzz_x, tr_yzzzz_y, tr_yzzzz_z, tr_zzzz_0, tr_zzzz_xy, tr_zzzz_yy, tr_zzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zzzz_x[i] = 2.0 * tr_yzzzz_x[i] * tbe_0 + 2.0 * tr_zzzz_xy[i] * tke_0;

        tr_0_0_y_zzzz_y[i] = 2.0 * tr_yzzzz_y[i] * tbe_0 + 2.0 * tr_zzzz_yy[i] * tke_0 - tr_zzzz_0[i];

        tr_0_0_y_zzzz_z[i] = 2.0 * tr_yzzzz_z[i] * tbe_0 + 2.0 * tr_zzzz_yz[i] * tke_0;
    }

    // Set up 90-93 components of targeted buffer : GP

    auto tr_0_0_z_xxxx_x = pbuffer.data(idx_op_geom_010_gp + 90);

    auto tr_0_0_z_xxxx_y = pbuffer.data(idx_op_geom_010_gp + 91);

    auto tr_0_0_z_xxxx_z = pbuffer.data(idx_op_geom_010_gp + 92);

    #pragma omp simd aligned(tr_0_0_z_xxxx_x, tr_0_0_z_xxxx_y, tr_0_0_z_xxxx_z, tr_xxxx_0, tr_xxxx_xz, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxz_x, tr_xxxxz_y, tr_xxxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxx_x[i] = 2.0 * tr_xxxxz_x[i] * tbe_0 + 2.0 * tr_xxxx_xz[i] * tke_0;

        tr_0_0_z_xxxx_y[i] = 2.0 * tr_xxxxz_y[i] * tbe_0 + 2.0 * tr_xxxx_yz[i] * tke_0;

        tr_0_0_z_xxxx_z[i] = 2.0 * tr_xxxxz_z[i] * tbe_0 + 2.0 * tr_xxxx_zz[i] * tke_0 - tr_xxxx_0[i];
    }

    // Set up 93-96 components of targeted buffer : GP

    auto tr_0_0_z_xxxy_x = pbuffer.data(idx_op_geom_010_gp + 93);

    auto tr_0_0_z_xxxy_y = pbuffer.data(idx_op_geom_010_gp + 94);

    auto tr_0_0_z_xxxy_z = pbuffer.data(idx_op_geom_010_gp + 95);

    #pragma omp simd aligned(tr_0_0_z_xxxy_x, tr_0_0_z_xxxy_y, tr_0_0_z_xxxy_z, tr_xxxy_0, tr_xxxy_xz, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxy_x[i] = 2.0 * tr_xxxyz_x[i] * tbe_0 + 2.0 * tr_xxxy_xz[i] * tke_0;

        tr_0_0_z_xxxy_y[i] = 2.0 * tr_xxxyz_y[i] * tbe_0 + 2.0 * tr_xxxy_yz[i] * tke_0;

        tr_0_0_z_xxxy_z[i] = 2.0 * tr_xxxyz_z[i] * tbe_0 + 2.0 * tr_xxxy_zz[i] * tke_0 - tr_xxxy_0[i];
    }

    // Set up 96-99 components of targeted buffer : GP

    auto tr_0_0_z_xxxz_x = pbuffer.data(idx_op_geom_010_gp + 96);

    auto tr_0_0_z_xxxz_y = pbuffer.data(idx_op_geom_010_gp + 97);

    auto tr_0_0_z_xxxz_z = pbuffer.data(idx_op_geom_010_gp + 98);

    #pragma omp simd aligned(tr_0_0_z_xxxz_x, tr_0_0_z_xxxz_y, tr_0_0_z_xxxz_z, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxz_0, tr_xxxz_xz, tr_xxxz_yz, tr_xxxz_zz, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxz_x[i] = 2.0 * tr_xxxzz_x[i] * tbe_0 + 2.0 * tr_xxxz_xz[i] * tke_0 - tr_xxx_x[i];

        tr_0_0_z_xxxz_y[i] = 2.0 * tr_xxxzz_y[i] * tbe_0 + 2.0 * tr_xxxz_yz[i] * tke_0 - tr_xxx_y[i];

        tr_0_0_z_xxxz_z[i] = 2.0 * tr_xxxzz_z[i] * tbe_0 + 2.0 * tr_xxxz_zz[i] * tke_0 - tr_xxx_z[i] - tr_xxxz_0[i];
    }

    // Set up 99-102 components of targeted buffer : GP

    auto tr_0_0_z_xxyy_x = pbuffer.data(idx_op_geom_010_gp + 99);

    auto tr_0_0_z_xxyy_y = pbuffer.data(idx_op_geom_010_gp + 100);

    auto tr_0_0_z_xxyy_z = pbuffer.data(idx_op_geom_010_gp + 101);

    #pragma omp simd aligned(tr_0_0_z_xxyy_x, tr_0_0_z_xxyy_y, tr_0_0_z_xxyy_z, tr_xxyy_0, tr_xxyy_xz, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxyy_x[i] = 2.0 * tr_xxyyz_x[i] * tbe_0 + 2.0 * tr_xxyy_xz[i] * tke_0;

        tr_0_0_z_xxyy_y[i] = 2.0 * tr_xxyyz_y[i] * tbe_0 + 2.0 * tr_xxyy_yz[i] * tke_0;

        tr_0_0_z_xxyy_z[i] = 2.0 * tr_xxyyz_z[i] * tbe_0 + 2.0 * tr_xxyy_zz[i] * tke_0 - tr_xxyy_0[i];
    }

    // Set up 102-105 components of targeted buffer : GP

    auto tr_0_0_z_xxyz_x = pbuffer.data(idx_op_geom_010_gp + 102);

    auto tr_0_0_z_xxyz_y = pbuffer.data(idx_op_geom_010_gp + 103);

    auto tr_0_0_z_xxyz_z = pbuffer.data(idx_op_geom_010_gp + 104);

    #pragma omp simd aligned(tr_0_0_z_xxyz_x, tr_0_0_z_xxyz_y, tr_0_0_z_xxyz_z, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyz_0, tr_xxyz_xz, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxyz_x[i] = 2.0 * tr_xxyzz_x[i] * tbe_0 + 2.0 * tr_xxyz_xz[i] * tke_0 - tr_xxy_x[i];

        tr_0_0_z_xxyz_y[i] = 2.0 * tr_xxyzz_y[i] * tbe_0 + 2.0 * tr_xxyz_yz[i] * tke_0 - tr_xxy_y[i];

        tr_0_0_z_xxyz_z[i] = 2.0 * tr_xxyzz_z[i] * tbe_0 + 2.0 * tr_xxyz_zz[i] * tke_0 - tr_xxy_z[i] - tr_xxyz_0[i];
    }

    // Set up 105-108 components of targeted buffer : GP

    auto tr_0_0_z_xxzz_x = pbuffer.data(idx_op_geom_010_gp + 105);

    auto tr_0_0_z_xxzz_y = pbuffer.data(idx_op_geom_010_gp + 106);

    auto tr_0_0_z_xxzz_z = pbuffer.data(idx_op_geom_010_gp + 107);

    #pragma omp simd aligned(tr_0_0_z_xxzz_x, tr_0_0_z_xxzz_y, tr_0_0_z_xxzz_z, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xxzz_0, tr_xxzz_xz, tr_xxzz_yz, tr_xxzz_zz, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxzz_x[i] = 2.0 * tr_xxzzz_x[i] * tbe_0 + 2.0 * tr_xxzz_xz[i] * tke_0 - 2.0 * tr_xxz_x[i];

        tr_0_0_z_xxzz_y[i] = 2.0 * tr_xxzzz_y[i] * tbe_0 + 2.0 * tr_xxzz_yz[i] * tke_0 - 2.0 * tr_xxz_y[i];

        tr_0_0_z_xxzz_z[i] = 2.0 * tr_xxzzz_z[i] * tbe_0 + 2.0 * tr_xxzz_zz[i] * tke_0 - 2.0 * tr_xxz_z[i] - tr_xxzz_0[i];
    }

    // Set up 108-111 components of targeted buffer : GP

    auto tr_0_0_z_xyyy_x = pbuffer.data(idx_op_geom_010_gp + 108);

    auto tr_0_0_z_xyyy_y = pbuffer.data(idx_op_geom_010_gp + 109);

    auto tr_0_0_z_xyyy_z = pbuffer.data(idx_op_geom_010_gp + 110);

    #pragma omp simd aligned(tr_0_0_z_xyyy_x, tr_0_0_z_xyyy_y, tr_0_0_z_xyyy_z, tr_xyyy_0, tr_xyyy_xz, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyyy_x[i] = 2.0 * tr_xyyyz_x[i] * tbe_0 + 2.0 * tr_xyyy_xz[i] * tke_0;

        tr_0_0_z_xyyy_y[i] = 2.0 * tr_xyyyz_y[i] * tbe_0 + 2.0 * tr_xyyy_yz[i] * tke_0;

        tr_0_0_z_xyyy_z[i] = 2.0 * tr_xyyyz_z[i] * tbe_0 + 2.0 * tr_xyyy_zz[i] * tke_0 - tr_xyyy_0[i];
    }

    // Set up 111-114 components of targeted buffer : GP

    auto tr_0_0_z_xyyz_x = pbuffer.data(idx_op_geom_010_gp + 111);

    auto tr_0_0_z_xyyz_y = pbuffer.data(idx_op_geom_010_gp + 112);

    auto tr_0_0_z_xyyz_z = pbuffer.data(idx_op_geom_010_gp + 113);

    #pragma omp simd aligned(tr_0_0_z_xyyz_x, tr_0_0_z_xyyz_y, tr_0_0_z_xyyz_z, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyz_0, tr_xyyz_xz, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyyz_x[i] = 2.0 * tr_xyyzz_x[i] * tbe_0 + 2.0 * tr_xyyz_xz[i] * tke_0 - tr_xyy_x[i];

        tr_0_0_z_xyyz_y[i] = 2.0 * tr_xyyzz_y[i] * tbe_0 + 2.0 * tr_xyyz_yz[i] * tke_0 - tr_xyy_y[i];

        tr_0_0_z_xyyz_z[i] = 2.0 * tr_xyyzz_z[i] * tbe_0 + 2.0 * tr_xyyz_zz[i] * tke_0 - tr_xyy_z[i] - tr_xyyz_0[i];
    }

    // Set up 114-117 components of targeted buffer : GP

    auto tr_0_0_z_xyzz_x = pbuffer.data(idx_op_geom_010_gp + 114);

    auto tr_0_0_z_xyzz_y = pbuffer.data(idx_op_geom_010_gp + 115);

    auto tr_0_0_z_xyzz_z = pbuffer.data(idx_op_geom_010_gp + 116);

    #pragma omp simd aligned(tr_0_0_z_xyzz_x, tr_0_0_z_xyzz_y, tr_0_0_z_xyzz_z, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xyzz_xz, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyzz_x[i] = 2.0 * tr_xyzzz_x[i] * tbe_0 + 2.0 * tr_xyzz_xz[i] * tke_0 - 2.0 * tr_xyz_x[i];

        tr_0_0_z_xyzz_y[i] = 2.0 * tr_xyzzz_y[i] * tbe_0 + 2.0 * tr_xyzz_yz[i] * tke_0 - 2.0 * tr_xyz_y[i];

        tr_0_0_z_xyzz_z[i] = 2.0 * tr_xyzzz_z[i] * tbe_0 + 2.0 * tr_xyzz_zz[i] * tke_0 - 2.0 * tr_xyz_z[i] - tr_xyzz_0[i];
    }

    // Set up 117-120 components of targeted buffer : GP

    auto tr_0_0_z_xzzz_x = pbuffer.data(idx_op_geom_010_gp + 117);

    auto tr_0_0_z_xzzz_y = pbuffer.data(idx_op_geom_010_gp + 118);

    auto tr_0_0_z_xzzz_z = pbuffer.data(idx_op_geom_010_gp + 119);

    #pragma omp simd aligned(tr_0_0_z_xzzz_x, tr_0_0_z_xzzz_y, tr_0_0_z_xzzz_z, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_xzzz_0, tr_xzzz_xz, tr_xzzz_yz, tr_xzzz_zz, tr_xzzzz_x, tr_xzzzz_y, tr_xzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xzzz_x[i] = 2.0 * tr_xzzzz_x[i] * tbe_0 + 2.0 * tr_xzzz_xz[i] * tke_0 - 3.0 * tr_xzz_x[i];

        tr_0_0_z_xzzz_y[i] = 2.0 * tr_xzzzz_y[i] * tbe_0 + 2.0 * tr_xzzz_yz[i] * tke_0 - 3.0 * tr_xzz_y[i];

        tr_0_0_z_xzzz_z[i] = 2.0 * tr_xzzzz_z[i] * tbe_0 + 2.0 * tr_xzzz_zz[i] * tke_0 - 3.0 * tr_xzz_z[i] - tr_xzzz_0[i];
    }

    // Set up 120-123 components of targeted buffer : GP

    auto tr_0_0_z_yyyy_x = pbuffer.data(idx_op_geom_010_gp + 120);

    auto tr_0_0_z_yyyy_y = pbuffer.data(idx_op_geom_010_gp + 121);

    auto tr_0_0_z_yyyy_z = pbuffer.data(idx_op_geom_010_gp + 122);

    #pragma omp simd aligned(tr_0_0_z_yyyy_x, tr_0_0_z_yyyy_y, tr_0_0_z_yyyy_z, tr_yyyy_0, tr_yyyy_xz, tr_yyyy_yz, tr_yyyy_zz, tr_yyyyz_x, tr_yyyyz_y, tr_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyyy_x[i] = 2.0 * tr_yyyyz_x[i] * tbe_0 + 2.0 * tr_yyyy_xz[i] * tke_0;

        tr_0_0_z_yyyy_y[i] = 2.0 * tr_yyyyz_y[i] * tbe_0 + 2.0 * tr_yyyy_yz[i] * tke_0;

        tr_0_0_z_yyyy_z[i] = 2.0 * tr_yyyyz_z[i] * tbe_0 + 2.0 * tr_yyyy_zz[i] * tke_0 - tr_yyyy_0[i];
    }

    // Set up 123-126 components of targeted buffer : GP

    auto tr_0_0_z_yyyz_x = pbuffer.data(idx_op_geom_010_gp + 123);

    auto tr_0_0_z_yyyz_y = pbuffer.data(idx_op_geom_010_gp + 124);

    auto tr_0_0_z_yyyz_z = pbuffer.data(idx_op_geom_010_gp + 125);

    #pragma omp simd aligned(tr_0_0_z_yyyz_x, tr_0_0_z_yyyz_y, tr_0_0_z_yyyz_z, tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyz_0, tr_yyyz_xz, tr_yyyz_yz, tr_yyyz_zz, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyyz_x[i] = 2.0 * tr_yyyzz_x[i] * tbe_0 + 2.0 * tr_yyyz_xz[i] * tke_0 - tr_yyy_x[i];

        tr_0_0_z_yyyz_y[i] = 2.0 * tr_yyyzz_y[i] * tbe_0 + 2.0 * tr_yyyz_yz[i] * tke_0 - tr_yyy_y[i];

        tr_0_0_z_yyyz_z[i] = 2.0 * tr_yyyzz_z[i] * tbe_0 + 2.0 * tr_yyyz_zz[i] * tke_0 - tr_yyy_z[i] - tr_yyyz_0[i];
    }

    // Set up 126-129 components of targeted buffer : GP

    auto tr_0_0_z_yyzz_x = pbuffer.data(idx_op_geom_010_gp + 126);

    auto tr_0_0_z_yyzz_y = pbuffer.data(idx_op_geom_010_gp + 127);

    auto tr_0_0_z_yyzz_z = pbuffer.data(idx_op_geom_010_gp + 128);

    #pragma omp simd aligned(tr_0_0_z_yyzz_x, tr_0_0_z_yyzz_y, tr_0_0_z_yyzz_z, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yyzz_0, tr_yyzz_xz, tr_yyzz_yz, tr_yyzz_zz, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyzz_x[i] = 2.0 * tr_yyzzz_x[i] * tbe_0 + 2.0 * tr_yyzz_xz[i] * tke_0 - 2.0 * tr_yyz_x[i];

        tr_0_0_z_yyzz_y[i] = 2.0 * tr_yyzzz_y[i] * tbe_0 + 2.0 * tr_yyzz_yz[i] * tke_0 - 2.0 * tr_yyz_y[i];

        tr_0_0_z_yyzz_z[i] = 2.0 * tr_yyzzz_z[i] * tbe_0 + 2.0 * tr_yyzz_zz[i] * tke_0 - 2.0 * tr_yyz_z[i] - tr_yyzz_0[i];
    }

    // Set up 129-132 components of targeted buffer : GP

    auto tr_0_0_z_yzzz_x = pbuffer.data(idx_op_geom_010_gp + 129);

    auto tr_0_0_z_yzzz_y = pbuffer.data(idx_op_geom_010_gp + 130);

    auto tr_0_0_z_yzzz_z = pbuffer.data(idx_op_geom_010_gp + 131);

    #pragma omp simd aligned(tr_0_0_z_yzzz_x, tr_0_0_z_yzzz_y, tr_0_0_z_yzzz_z, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_yzzz_0, tr_yzzz_xz, tr_yzzz_yz, tr_yzzz_zz, tr_yzzzz_x, tr_yzzzz_y, tr_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yzzz_x[i] = 2.0 * tr_yzzzz_x[i] * tbe_0 + 2.0 * tr_yzzz_xz[i] * tke_0 - 3.0 * tr_yzz_x[i];

        tr_0_0_z_yzzz_y[i] = 2.0 * tr_yzzzz_y[i] * tbe_0 + 2.0 * tr_yzzz_yz[i] * tke_0 - 3.0 * tr_yzz_y[i];

        tr_0_0_z_yzzz_z[i] = 2.0 * tr_yzzzz_z[i] * tbe_0 + 2.0 * tr_yzzz_zz[i] * tke_0 - 3.0 * tr_yzz_z[i] - tr_yzzz_0[i];
    }

    // Set up 132-135 components of targeted buffer : GP

    auto tr_0_0_z_zzzz_x = pbuffer.data(idx_op_geom_010_gp + 132);

    auto tr_0_0_z_zzzz_y = pbuffer.data(idx_op_geom_010_gp + 133);

    auto tr_0_0_z_zzzz_z = pbuffer.data(idx_op_geom_010_gp + 134);

    #pragma omp simd aligned(tr_0_0_z_zzzz_x, tr_0_0_z_zzzz_y, tr_0_0_z_zzzz_z, tr_zzz_x, tr_zzz_y, tr_zzz_z, tr_zzzz_0, tr_zzzz_xz, tr_zzzz_yz, tr_zzzz_zz, tr_zzzzz_x, tr_zzzzz_y, tr_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zzzz_x[i] = 2.0 * tr_zzzzz_x[i] * tbe_0 + 2.0 * tr_zzzz_xz[i] * tke_0 - 4.0 * tr_zzz_x[i];

        tr_0_0_z_zzzz_y[i] = 2.0 * tr_zzzzz_y[i] * tbe_0 + 2.0 * tr_zzzz_yz[i] * tke_0 - 4.0 * tr_zzz_y[i];

        tr_0_0_z_zzzz_z[i] = 2.0 * tr_zzzzz_z[i] * tbe_0 + 2.0 * tr_zzzz_zz[i] * tke_0 - 4.0 * tr_zzz_z[i] - tr_zzzz_0[i];
    }

}

} // t2cgeom namespace

