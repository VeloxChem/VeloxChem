#include "GeometricalDerivatives010ForPG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_pg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_pg,
                         const int idx_op_sg,
                         const int idx_op_pf,
                         const int idx_op_ph,
                         const int idx_op_dg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SG

    auto tr_0_xxxx = pbuffer.data(idx_op_sg);

    auto tr_0_xxxy = pbuffer.data(idx_op_sg + 1);

    auto tr_0_xxxz = pbuffer.data(idx_op_sg + 2);

    auto tr_0_xxyy = pbuffer.data(idx_op_sg + 3);

    auto tr_0_xxyz = pbuffer.data(idx_op_sg + 4);

    auto tr_0_xxzz = pbuffer.data(idx_op_sg + 5);

    auto tr_0_xyyy = pbuffer.data(idx_op_sg + 6);

    auto tr_0_xyyz = pbuffer.data(idx_op_sg + 7);

    auto tr_0_xyzz = pbuffer.data(idx_op_sg + 8);

    auto tr_0_xzzz = pbuffer.data(idx_op_sg + 9);

    auto tr_0_yyyy = pbuffer.data(idx_op_sg + 10);

    auto tr_0_yyyz = pbuffer.data(idx_op_sg + 11);

    auto tr_0_yyzz = pbuffer.data(idx_op_sg + 12);

    auto tr_0_yzzz = pbuffer.data(idx_op_sg + 13);

    auto tr_0_zzzz = pbuffer.data(idx_op_sg + 14);

    // Set up components of auxiliary buffer : PF

    auto tr_x_xxx = pbuffer.data(idx_op_pf);

    auto tr_x_xxy = pbuffer.data(idx_op_pf + 1);

    auto tr_x_xxz = pbuffer.data(idx_op_pf + 2);

    auto tr_x_xyy = pbuffer.data(idx_op_pf + 3);

    auto tr_x_xyz = pbuffer.data(idx_op_pf + 4);

    auto tr_x_xzz = pbuffer.data(idx_op_pf + 5);

    auto tr_x_yyy = pbuffer.data(idx_op_pf + 6);

    auto tr_x_yyz = pbuffer.data(idx_op_pf + 7);

    auto tr_x_yzz = pbuffer.data(idx_op_pf + 8);

    auto tr_x_zzz = pbuffer.data(idx_op_pf + 9);

    auto tr_y_xxx = pbuffer.data(idx_op_pf + 10);

    auto tr_y_xxy = pbuffer.data(idx_op_pf + 11);

    auto tr_y_xxz = pbuffer.data(idx_op_pf + 12);

    auto tr_y_xyy = pbuffer.data(idx_op_pf + 13);

    auto tr_y_xyz = pbuffer.data(idx_op_pf + 14);

    auto tr_y_xzz = pbuffer.data(idx_op_pf + 15);

    auto tr_y_yyy = pbuffer.data(idx_op_pf + 16);

    auto tr_y_yyz = pbuffer.data(idx_op_pf + 17);

    auto tr_y_yzz = pbuffer.data(idx_op_pf + 18);

    auto tr_y_zzz = pbuffer.data(idx_op_pf + 19);

    auto tr_z_xxx = pbuffer.data(idx_op_pf + 20);

    auto tr_z_xxy = pbuffer.data(idx_op_pf + 21);

    auto tr_z_xxz = pbuffer.data(idx_op_pf + 22);

    auto tr_z_xyy = pbuffer.data(idx_op_pf + 23);

    auto tr_z_xyz = pbuffer.data(idx_op_pf + 24);

    auto tr_z_xzz = pbuffer.data(idx_op_pf + 25);

    auto tr_z_yyy = pbuffer.data(idx_op_pf + 26);

    auto tr_z_yyz = pbuffer.data(idx_op_pf + 27);

    auto tr_z_yzz = pbuffer.data(idx_op_pf + 28);

    auto tr_z_zzz = pbuffer.data(idx_op_pf + 29);

    // Set up components of auxiliary buffer : PH

    auto tr_x_xxxxx = pbuffer.data(idx_op_ph);

    auto tr_x_xxxxy = pbuffer.data(idx_op_ph + 1);

    auto tr_x_xxxxz = pbuffer.data(idx_op_ph + 2);

    auto tr_x_xxxyy = pbuffer.data(idx_op_ph + 3);

    auto tr_x_xxxyz = pbuffer.data(idx_op_ph + 4);

    auto tr_x_xxxzz = pbuffer.data(idx_op_ph + 5);

    auto tr_x_xxyyy = pbuffer.data(idx_op_ph + 6);

    auto tr_x_xxyyz = pbuffer.data(idx_op_ph + 7);

    auto tr_x_xxyzz = pbuffer.data(idx_op_ph + 8);

    auto tr_x_xxzzz = pbuffer.data(idx_op_ph + 9);

    auto tr_x_xyyyy = pbuffer.data(idx_op_ph + 10);

    auto tr_x_xyyyz = pbuffer.data(idx_op_ph + 11);

    auto tr_x_xyyzz = pbuffer.data(idx_op_ph + 12);

    auto tr_x_xyzzz = pbuffer.data(idx_op_ph + 13);

    auto tr_x_xzzzz = pbuffer.data(idx_op_ph + 14);

    auto tr_x_yyyyy = pbuffer.data(idx_op_ph + 15);

    auto tr_x_yyyyz = pbuffer.data(idx_op_ph + 16);

    auto tr_x_yyyzz = pbuffer.data(idx_op_ph + 17);

    auto tr_x_yyzzz = pbuffer.data(idx_op_ph + 18);

    auto tr_x_yzzzz = pbuffer.data(idx_op_ph + 19);

    auto tr_x_zzzzz = pbuffer.data(idx_op_ph + 20);

    auto tr_y_xxxxx = pbuffer.data(idx_op_ph + 21);

    auto tr_y_xxxxy = pbuffer.data(idx_op_ph + 22);

    auto tr_y_xxxxz = pbuffer.data(idx_op_ph + 23);

    auto tr_y_xxxyy = pbuffer.data(idx_op_ph + 24);

    auto tr_y_xxxyz = pbuffer.data(idx_op_ph + 25);

    auto tr_y_xxxzz = pbuffer.data(idx_op_ph + 26);

    auto tr_y_xxyyy = pbuffer.data(idx_op_ph + 27);

    auto tr_y_xxyyz = pbuffer.data(idx_op_ph + 28);

    auto tr_y_xxyzz = pbuffer.data(idx_op_ph + 29);

    auto tr_y_xxzzz = pbuffer.data(idx_op_ph + 30);

    auto tr_y_xyyyy = pbuffer.data(idx_op_ph + 31);

    auto tr_y_xyyyz = pbuffer.data(idx_op_ph + 32);

    auto tr_y_xyyzz = pbuffer.data(idx_op_ph + 33);

    auto tr_y_xyzzz = pbuffer.data(idx_op_ph + 34);

    auto tr_y_xzzzz = pbuffer.data(idx_op_ph + 35);

    auto tr_y_yyyyy = pbuffer.data(idx_op_ph + 36);

    auto tr_y_yyyyz = pbuffer.data(idx_op_ph + 37);

    auto tr_y_yyyzz = pbuffer.data(idx_op_ph + 38);

    auto tr_y_yyzzz = pbuffer.data(idx_op_ph + 39);

    auto tr_y_yzzzz = pbuffer.data(idx_op_ph + 40);

    auto tr_y_zzzzz = pbuffer.data(idx_op_ph + 41);

    auto tr_z_xxxxx = pbuffer.data(idx_op_ph + 42);

    auto tr_z_xxxxy = pbuffer.data(idx_op_ph + 43);

    auto tr_z_xxxxz = pbuffer.data(idx_op_ph + 44);

    auto tr_z_xxxyy = pbuffer.data(idx_op_ph + 45);

    auto tr_z_xxxyz = pbuffer.data(idx_op_ph + 46);

    auto tr_z_xxxzz = pbuffer.data(idx_op_ph + 47);

    auto tr_z_xxyyy = pbuffer.data(idx_op_ph + 48);

    auto tr_z_xxyyz = pbuffer.data(idx_op_ph + 49);

    auto tr_z_xxyzz = pbuffer.data(idx_op_ph + 50);

    auto tr_z_xxzzz = pbuffer.data(idx_op_ph + 51);

    auto tr_z_xyyyy = pbuffer.data(idx_op_ph + 52);

    auto tr_z_xyyyz = pbuffer.data(idx_op_ph + 53);

    auto tr_z_xyyzz = pbuffer.data(idx_op_ph + 54);

    auto tr_z_xyzzz = pbuffer.data(idx_op_ph + 55);

    auto tr_z_xzzzz = pbuffer.data(idx_op_ph + 56);

    auto tr_z_yyyyy = pbuffer.data(idx_op_ph + 57);

    auto tr_z_yyyyz = pbuffer.data(idx_op_ph + 58);

    auto tr_z_yyyzz = pbuffer.data(idx_op_ph + 59);

    auto tr_z_yyzzz = pbuffer.data(idx_op_ph + 60);

    auto tr_z_yzzzz = pbuffer.data(idx_op_ph + 61);

    auto tr_z_zzzzz = pbuffer.data(idx_op_ph + 62);

    // Set up components of auxiliary buffer : DG

    auto tr_xx_xxxx = pbuffer.data(idx_op_dg);

    auto tr_xx_xxxy = pbuffer.data(idx_op_dg + 1);

    auto tr_xx_xxxz = pbuffer.data(idx_op_dg + 2);

    auto tr_xx_xxyy = pbuffer.data(idx_op_dg + 3);

    auto tr_xx_xxyz = pbuffer.data(idx_op_dg + 4);

    auto tr_xx_xxzz = pbuffer.data(idx_op_dg + 5);

    auto tr_xx_xyyy = pbuffer.data(idx_op_dg + 6);

    auto tr_xx_xyyz = pbuffer.data(idx_op_dg + 7);

    auto tr_xx_xyzz = pbuffer.data(idx_op_dg + 8);

    auto tr_xx_xzzz = pbuffer.data(idx_op_dg + 9);

    auto tr_xx_yyyy = pbuffer.data(idx_op_dg + 10);

    auto tr_xx_yyyz = pbuffer.data(idx_op_dg + 11);

    auto tr_xx_yyzz = pbuffer.data(idx_op_dg + 12);

    auto tr_xx_yzzz = pbuffer.data(idx_op_dg + 13);

    auto tr_xx_zzzz = pbuffer.data(idx_op_dg + 14);

    auto tr_xy_xxxx = pbuffer.data(idx_op_dg + 15);

    auto tr_xy_xxxy = pbuffer.data(idx_op_dg + 16);

    auto tr_xy_xxxz = pbuffer.data(idx_op_dg + 17);

    auto tr_xy_xxyy = pbuffer.data(idx_op_dg + 18);

    auto tr_xy_xxyz = pbuffer.data(idx_op_dg + 19);

    auto tr_xy_xxzz = pbuffer.data(idx_op_dg + 20);

    auto tr_xy_xyyy = pbuffer.data(idx_op_dg + 21);

    auto tr_xy_xyyz = pbuffer.data(idx_op_dg + 22);

    auto tr_xy_xyzz = pbuffer.data(idx_op_dg + 23);

    auto tr_xy_xzzz = pbuffer.data(idx_op_dg + 24);

    auto tr_xy_yyyy = pbuffer.data(idx_op_dg + 25);

    auto tr_xy_yyyz = pbuffer.data(idx_op_dg + 26);

    auto tr_xy_yyzz = pbuffer.data(idx_op_dg + 27);

    auto tr_xy_yzzz = pbuffer.data(idx_op_dg + 28);

    auto tr_xy_zzzz = pbuffer.data(idx_op_dg + 29);

    auto tr_xz_xxxx = pbuffer.data(idx_op_dg + 30);

    auto tr_xz_xxxy = pbuffer.data(idx_op_dg + 31);

    auto tr_xz_xxxz = pbuffer.data(idx_op_dg + 32);

    auto tr_xz_xxyy = pbuffer.data(idx_op_dg + 33);

    auto tr_xz_xxyz = pbuffer.data(idx_op_dg + 34);

    auto tr_xz_xxzz = pbuffer.data(idx_op_dg + 35);

    auto tr_xz_xyyy = pbuffer.data(idx_op_dg + 36);

    auto tr_xz_xyyz = pbuffer.data(idx_op_dg + 37);

    auto tr_xz_xyzz = pbuffer.data(idx_op_dg + 38);

    auto tr_xz_xzzz = pbuffer.data(idx_op_dg + 39);

    auto tr_xz_yyyy = pbuffer.data(idx_op_dg + 40);

    auto tr_xz_yyyz = pbuffer.data(idx_op_dg + 41);

    auto tr_xz_yyzz = pbuffer.data(idx_op_dg + 42);

    auto tr_xz_yzzz = pbuffer.data(idx_op_dg + 43);

    auto tr_xz_zzzz = pbuffer.data(idx_op_dg + 44);

    auto tr_yy_xxxx = pbuffer.data(idx_op_dg + 45);

    auto tr_yy_xxxy = pbuffer.data(idx_op_dg + 46);

    auto tr_yy_xxxz = pbuffer.data(idx_op_dg + 47);

    auto tr_yy_xxyy = pbuffer.data(idx_op_dg + 48);

    auto tr_yy_xxyz = pbuffer.data(idx_op_dg + 49);

    auto tr_yy_xxzz = pbuffer.data(idx_op_dg + 50);

    auto tr_yy_xyyy = pbuffer.data(idx_op_dg + 51);

    auto tr_yy_xyyz = pbuffer.data(idx_op_dg + 52);

    auto tr_yy_xyzz = pbuffer.data(idx_op_dg + 53);

    auto tr_yy_xzzz = pbuffer.data(idx_op_dg + 54);

    auto tr_yy_yyyy = pbuffer.data(idx_op_dg + 55);

    auto tr_yy_yyyz = pbuffer.data(idx_op_dg + 56);

    auto tr_yy_yyzz = pbuffer.data(idx_op_dg + 57);

    auto tr_yy_yzzz = pbuffer.data(idx_op_dg + 58);

    auto tr_yy_zzzz = pbuffer.data(idx_op_dg + 59);

    auto tr_yz_xxxx = pbuffer.data(idx_op_dg + 60);

    auto tr_yz_xxxy = pbuffer.data(idx_op_dg + 61);

    auto tr_yz_xxxz = pbuffer.data(idx_op_dg + 62);

    auto tr_yz_xxyy = pbuffer.data(idx_op_dg + 63);

    auto tr_yz_xxyz = pbuffer.data(idx_op_dg + 64);

    auto tr_yz_xxzz = pbuffer.data(idx_op_dg + 65);

    auto tr_yz_xyyy = pbuffer.data(idx_op_dg + 66);

    auto tr_yz_xyyz = pbuffer.data(idx_op_dg + 67);

    auto tr_yz_xyzz = pbuffer.data(idx_op_dg + 68);

    auto tr_yz_xzzz = pbuffer.data(idx_op_dg + 69);

    auto tr_yz_yyyy = pbuffer.data(idx_op_dg + 70);

    auto tr_yz_yyyz = pbuffer.data(idx_op_dg + 71);

    auto tr_yz_yyzz = pbuffer.data(idx_op_dg + 72);

    auto tr_yz_yzzz = pbuffer.data(idx_op_dg + 73);

    auto tr_yz_zzzz = pbuffer.data(idx_op_dg + 74);

    auto tr_zz_xxxx = pbuffer.data(idx_op_dg + 75);

    auto tr_zz_xxxy = pbuffer.data(idx_op_dg + 76);

    auto tr_zz_xxxz = pbuffer.data(idx_op_dg + 77);

    auto tr_zz_xxyy = pbuffer.data(idx_op_dg + 78);

    auto tr_zz_xxyz = pbuffer.data(idx_op_dg + 79);

    auto tr_zz_xxzz = pbuffer.data(idx_op_dg + 80);

    auto tr_zz_xyyy = pbuffer.data(idx_op_dg + 81);

    auto tr_zz_xyyz = pbuffer.data(idx_op_dg + 82);

    auto tr_zz_xyzz = pbuffer.data(idx_op_dg + 83);

    auto tr_zz_xzzz = pbuffer.data(idx_op_dg + 84);

    auto tr_zz_yyyy = pbuffer.data(idx_op_dg + 85);

    auto tr_zz_yyyz = pbuffer.data(idx_op_dg + 86);

    auto tr_zz_yyzz = pbuffer.data(idx_op_dg + 87);

    auto tr_zz_yzzz = pbuffer.data(idx_op_dg + 88);

    auto tr_zz_zzzz = pbuffer.data(idx_op_dg + 89);

    // Set up 0-15 components of targeted buffer : PG

    auto tr_0_0_x_x_xxxx = pbuffer.data(idx_op_geom_010_pg);

    auto tr_0_0_x_x_xxxy = pbuffer.data(idx_op_geom_010_pg + 1);

    auto tr_0_0_x_x_xxxz = pbuffer.data(idx_op_geom_010_pg + 2);

    auto tr_0_0_x_x_xxyy = pbuffer.data(idx_op_geom_010_pg + 3);

    auto tr_0_0_x_x_xxyz = pbuffer.data(idx_op_geom_010_pg + 4);

    auto tr_0_0_x_x_xxzz = pbuffer.data(idx_op_geom_010_pg + 5);

    auto tr_0_0_x_x_xyyy = pbuffer.data(idx_op_geom_010_pg + 6);

    auto tr_0_0_x_x_xyyz = pbuffer.data(idx_op_geom_010_pg + 7);

    auto tr_0_0_x_x_xyzz = pbuffer.data(idx_op_geom_010_pg + 8);

    auto tr_0_0_x_x_xzzz = pbuffer.data(idx_op_geom_010_pg + 9);

    auto tr_0_0_x_x_yyyy = pbuffer.data(idx_op_geom_010_pg + 10);

    auto tr_0_0_x_x_yyyz = pbuffer.data(idx_op_geom_010_pg + 11);

    auto tr_0_0_x_x_yyzz = pbuffer.data(idx_op_geom_010_pg + 12);

    auto tr_0_0_x_x_yzzz = pbuffer.data(idx_op_geom_010_pg + 13);

    auto tr_0_0_x_x_zzzz = pbuffer.data(idx_op_geom_010_pg + 14);

    #pragma omp simd aligned(tr_0_0_x_x_xxxx, tr_0_0_x_x_xxxy, tr_0_0_x_x_xxxz, tr_0_0_x_x_xxyy, tr_0_0_x_x_xxyz, tr_0_0_x_x_xxzz, tr_0_0_x_x_xyyy, tr_0_0_x_x_xyyz, tr_0_0_x_x_xyzz, tr_0_0_x_x_xzzz, tr_0_0_x_x_yyyy, tr_0_0_x_x_yyyz, tr_0_0_x_x_yyzz, tr_0_0_x_x_yzzz, tr_0_0_x_x_zzzz, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_x_xxx, tr_x_xxxxx, tr_x_xxxxy, tr_x_xxxxz, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_x_xxxx[i] = 2.0 * tr_xx_xxxx[i] * tbe_0 + 2.0 * tr_x_xxxxx[i] * tke_0 - tr_0_xxxx[i] - 4.0 * tr_x_xxx[i];

        tr_0_0_x_x_xxxy[i] = 2.0 * tr_xx_xxxy[i] * tbe_0 + 2.0 * tr_x_xxxxy[i] * tke_0 - tr_0_xxxy[i] - 3.0 * tr_x_xxy[i];

        tr_0_0_x_x_xxxz[i] = 2.0 * tr_xx_xxxz[i] * tbe_0 + 2.0 * tr_x_xxxxz[i] * tke_0 - tr_0_xxxz[i] - 3.0 * tr_x_xxz[i];

        tr_0_0_x_x_xxyy[i] = 2.0 * tr_xx_xxyy[i] * tbe_0 + 2.0 * tr_x_xxxyy[i] * tke_0 - tr_0_xxyy[i] - 2.0 * tr_x_xyy[i];

        tr_0_0_x_x_xxyz[i] = 2.0 * tr_xx_xxyz[i] * tbe_0 + 2.0 * tr_x_xxxyz[i] * tke_0 - tr_0_xxyz[i] - 2.0 * tr_x_xyz[i];

        tr_0_0_x_x_xxzz[i] = 2.0 * tr_xx_xxzz[i] * tbe_0 + 2.0 * tr_x_xxxzz[i] * tke_0 - tr_0_xxzz[i] - 2.0 * tr_x_xzz[i];

        tr_0_0_x_x_xyyy[i] = 2.0 * tr_xx_xyyy[i] * tbe_0 + 2.0 * tr_x_xxyyy[i] * tke_0 - tr_0_xyyy[i] - tr_x_yyy[i];

        tr_0_0_x_x_xyyz[i] = 2.0 * tr_xx_xyyz[i] * tbe_0 + 2.0 * tr_x_xxyyz[i] * tke_0 - tr_0_xyyz[i] - tr_x_yyz[i];

        tr_0_0_x_x_xyzz[i] = 2.0 * tr_xx_xyzz[i] * tbe_0 + 2.0 * tr_x_xxyzz[i] * tke_0 - tr_0_xyzz[i] - tr_x_yzz[i];

        tr_0_0_x_x_xzzz[i] = 2.0 * tr_xx_xzzz[i] * tbe_0 + 2.0 * tr_x_xxzzz[i] * tke_0 - tr_0_xzzz[i] - tr_x_zzz[i];

        tr_0_0_x_x_yyyy[i] = 2.0 * tr_xx_yyyy[i] * tbe_0 + 2.0 * tr_x_xyyyy[i] * tke_0 - tr_0_yyyy[i];

        tr_0_0_x_x_yyyz[i] = 2.0 * tr_xx_yyyz[i] * tbe_0 + 2.0 * tr_x_xyyyz[i] * tke_0 - tr_0_yyyz[i];

        tr_0_0_x_x_yyzz[i] = 2.0 * tr_xx_yyzz[i] * tbe_0 + 2.0 * tr_x_xyyzz[i] * tke_0 - tr_0_yyzz[i];

        tr_0_0_x_x_yzzz[i] = 2.0 * tr_xx_yzzz[i] * tbe_0 + 2.0 * tr_x_xyzzz[i] * tke_0 - tr_0_yzzz[i];

        tr_0_0_x_x_zzzz[i] = 2.0 * tr_xx_zzzz[i] * tbe_0 + 2.0 * tr_x_xzzzz[i] * tke_0 - tr_0_zzzz[i];
    }

    // Set up 15-30 components of targeted buffer : PG

    auto tr_0_0_x_y_xxxx = pbuffer.data(idx_op_geom_010_pg + 15);

    auto tr_0_0_x_y_xxxy = pbuffer.data(idx_op_geom_010_pg + 16);

    auto tr_0_0_x_y_xxxz = pbuffer.data(idx_op_geom_010_pg + 17);

    auto tr_0_0_x_y_xxyy = pbuffer.data(idx_op_geom_010_pg + 18);

    auto tr_0_0_x_y_xxyz = pbuffer.data(idx_op_geom_010_pg + 19);

    auto tr_0_0_x_y_xxzz = pbuffer.data(idx_op_geom_010_pg + 20);

    auto tr_0_0_x_y_xyyy = pbuffer.data(idx_op_geom_010_pg + 21);

    auto tr_0_0_x_y_xyyz = pbuffer.data(idx_op_geom_010_pg + 22);

    auto tr_0_0_x_y_xyzz = pbuffer.data(idx_op_geom_010_pg + 23);

    auto tr_0_0_x_y_xzzz = pbuffer.data(idx_op_geom_010_pg + 24);

    auto tr_0_0_x_y_yyyy = pbuffer.data(idx_op_geom_010_pg + 25);

    auto tr_0_0_x_y_yyyz = pbuffer.data(idx_op_geom_010_pg + 26);

    auto tr_0_0_x_y_yyzz = pbuffer.data(idx_op_geom_010_pg + 27);

    auto tr_0_0_x_y_yzzz = pbuffer.data(idx_op_geom_010_pg + 28);

    auto tr_0_0_x_y_zzzz = pbuffer.data(idx_op_geom_010_pg + 29);

    #pragma omp simd aligned(tr_0_0_x_y_xxxx, tr_0_0_x_y_xxxy, tr_0_0_x_y_xxxz, tr_0_0_x_y_xxyy, tr_0_0_x_y_xxyz, tr_0_0_x_y_xxzz, tr_0_0_x_y_xyyy, tr_0_0_x_y_xyyz, tr_0_0_x_y_xyzz, tr_0_0_x_y_xzzz, tr_0_0_x_y_yyyy, tr_0_0_x_y_yyyz, tr_0_0_x_y_yyzz, tr_0_0_x_y_yzzz, tr_0_0_x_y_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_y_xxx, tr_y_xxxxx, tr_y_xxxxy, tr_y_xxxxz, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_y_xxxx[i] = 2.0 * tr_xy_xxxx[i] * tbe_0 + 2.0 * tr_y_xxxxx[i] * tke_0 - 4.0 * tr_y_xxx[i];

        tr_0_0_x_y_xxxy[i] = 2.0 * tr_xy_xxxy[i] * tbe_0 + 2.0 * tr_y_xxxxy[i] * tke_0 - 3.0 * tr_y_xxy[i];

        tr_0_0_x_y_xxxz[i] = 2.0 * tr_xy_xxxz[i] * tbe_0 + 2.0 * tr_y_xxxxz[i] * tke_0 - 3.0 * tr_y_xxz[i];

        tr_0_0_x_y_xxyy[i] = 2.0 * tr_xy_xxyy[i] * tbe_0 + 2.0 * tr_y_xxxyy[i] * tke_0 - 2.0 * tr_y_xyy[i];

        tr_0_0_x_y_xxyz[i] = 2.0 * tr_xy_xxyz[i] * tbe_0 + 2.0 * tr_y_xxxyz[i] * tke_0 - 2.0 * tr_y_xyz[i];

        tr_0_0_x_y_xxzz[i] = 2.0 * tr_xy_xxzz[i] * tbe_0 + 2.0 * tr_y_xxxzz[i] * tke_0 - 2.0 * tr_y_xzz[i];

        tr_0_0_x_y_xyyy[i] = 2.0 * tr_xy_xyyy[i] * tbe_0 + 2.0 * tr_y_xxyyy[i] * tke_0 - tr_y_yyy[i];

        tr_0_0_x_y_xyyz[i] = 2.0 * tr_xy_xyyz[i] * tbe_0 + 2.0 * tr_y_xxyyz[i] * tke_0 - tr_y_yyz[i];

        tr_0_0_x_y_xyzz[i] = 2.0 * tr_xy_xyzz[i] * tbe_0 + 2.0 * tr_y_xxyzz[i] * tke_0 - tr_y_yzz[i];

        tr_0_0_x_y_xzzz[i] = 2.0 * tr_xy_xzzz[i] * tbe_0 + 2.0 * tr_y_xxzzz[i] * tke_0 - tr_y_zzz[i];

        tr_0_0_x_y_yyyy[i] = 2.0 * tr_xy_yyyy[i] * tbe_0 + 2.0 * tr_y_xyyyy[i] * tke_0;

        tr_0_0_x_y_yyyz[i] = 2.0 * tr_xy_yyyz[i] * tbe_0 + 2.0 * tr_y_xyyyz[i] * tke_0;

        tr_0_0_x_y_yyzz[i] = 2.0 * tr_xy_yyzz[i] * tbe_0 + 2.0 * tr_y_xyyzz[i] * tke_0;

        tr_0_0_x_y_yzzz[i] = 2.0 * tr_xy_yzzz[i] * tbe_0 + 2.0 * tr_y_xyzzz[i] * tke_0;

        tr_0_0_x_y_zzzz[i] = 2.0 * tr_xy_zzzz[i] * tbe_0 + 2.0 * tr_y_xzzzz[i] * tke_0;
    }

    // Set up 30-45 components of targeted buffer : PG

    auto tr_0_0_x_z_xxxx = pbuffer.data(idx_op_geom_010_pg + 30);

    auto tr_0_0_x_z_xxxy = pbuffer.data(idx_op_geom_010_pg + 31);

    auto tr_0_0_x_z_xxxz = pbuffer.data(idx_op_geom_010_pg + 32);

    auto tr_0_0_x_z_xxyy = pbuffer.data(idx_op_geom_010_pg + 33);

    auto tr_0_0_x_z_xxyz = pbuffer.data(idx_op_geom_010_pg + 34);

    auto tr_0_0_x_z_xxzz = pbuffer.data(idx_op_geom_010_pg + 35);

    auto tr_0_0_x_z_xyyy = pbuffer.data(idx_op_geom_010_pg + 36);

    auto tr_0_0_x_z_xyyz = pbuffer.data(idx_op_geom_010_pg + 37);

    auto tr_0_0_x_z_xyzz = pbuffer.data(idx_op_geom_010_pg + 38);

    auto tr_0_0_x_z_xzzz = pbuffer.data(idx_op_geom_010_pg + 39);

    auto tr_0_0_x_z_yyyy = pbuffer.data(idx_op_geom_010_pg + 40);

    auto tr_0_0_x_z_yyyz = pbuffer.data(idx_op_geom_010_pg + 41);

    auto tr_0_0_x_z_yyzz = pbuffer.data(idx_op_geom_010_pg + 42);

    auto tr_0_0_x_z_yzzz = pbuffer.data(idx_op_geom_010_pg + 43);

    auto tr_0_0_x_z_zzzz = pbuffer.data(idx_op_geom_010_pg + 44);

    #pragma omp simd aligned(tr_0_0_x_z_xxxx, tr_0_0_x_z_xxxy, tr_0_0_x_z_xxxz, tr_0_0_x_z_xxyy, tr_0_0_x_z_xxyz, tr_0_0_x_z_xxzz, tr_0_0_x_z_xyyy, tr_0_0_x_z_xyyz, tr_0_0_x_z_xyzz, tr_0_0_x_z_xzzz, tr_0_0_x_z_yyyy, tr_0_0_x_z_yyyz, tr_0_0_x_z_yyzz, tr_0_0_x_z_yzzz, tr_0_0_x_z_zzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_z_xxx, tr_z_xxxxx, tr_z_xxxxy, tr_z_xxxxz, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_z_xxxx[i] = 2.0 * tr_xz_xxxx[i] * tbe_0 + 2.0 * tr_z_xxxxx[i] * tke_0 - 4.0 * tr_z_xxx[i];

        tr_0_0_x_z_xxxy[i] = 2.0 * tr_xz_xxxy[i] * tbe_0 + 2.0 * tr_z_xxxxy[i] * tke_0 - 3.0 * tr_z_xxy[i];

        tr_0_0_x_z_xxxz[i] = 2.0 * tr_xz_xxxz[i] * tbe_0 + 2.0 * tr_z_xxxxz[i] * tke_0 - 3.0 * tr_z_xxz[i];

        tr_0_0_x_z_xxyy[i] = 2.0 * tr_xz_xxyy[i] * tbe_0 + 2.0 * tr_z_xxxyy[i] * tke_0 - 2.0 * tr_z_xyy[i];

        tr_0_0_x_z_xxyz[i] = 2.0 * tr_xz_xxyz[i] * tbe_0 + 2.0 * tr_z_xxxyz[i] * tke_0 - 2.0 * tr_z_xyz[i];

        tr_0_0_x_z_xxzz[i] = 2.0 * tr_xz_xxzz[i] * tbe_0 + 2.0 * tr_z_xxxzz[i] * tke_0 - 2.0 * tr_z_xzz[i];

        tr_0_0_x_z_xyyy[i] = 2.0 * tr_xz_xyyy[i] * tbe_0 + 2.0 * tr_z_xxyyy[i] * tke_0 - tr_z_yyy[i];

        tr_0_0_x_z_xyyz[i] = 2.0 * tr_xz_xyyz[i] * tbe_0 + 2.0 * tr_z_xxyyz[i] * tke_0 - tr_z_yyz[i];

        tr_0_0_x_z_xyzz[i] = 2.0 * tr_xz_xyzz[i] * tbe_0 + 2.0 * tr_z_xxyzz[i] * tke_0 - tr_z_yzz[i];

        tr_0_0_x_z_xzzz[i] = 2.0 * tr_xz_xzzz[i] * tbe_0 + 2.0 * tr_z_xxzzz[i] * tke_0 - tr_z_zzz[i];

        tr_0_0_x_z_yyyy[i] = 2.0 * tr_xz_yyyy[i] * tbe_0 + 2.0 * tr_z_xyyyy[i] * tke_0;

        tr_0_0_x_z_yyyz[i] = 2.0 * tr_xz_yyyz[i] * tbe_0 + 2.0 * tr_z_xyyyz[i] * tke_0;

        tr_0_0_x_z_yyzz[i] = 2.0 * tr_xz_yyzz[i] * tbe_0 + 2.0 * tr_z_xyyzz[i] * tke_0;

        tr_0_0_x_z_yzzz[i] = 2.0 * tr_xz_yzzz[i] * tbe_0 + 2.0 * tr_z_xyzzz[i] * tke_0;

        tr_0_0_x_z_zzzz[i] = 2.0 * tr_xz_zzzz[i] * tbe_0 + 2.0 * tr_z_xzzzz[i] * tke_0;
    }

    // Set up 45-60 components of targeted buffer : PG

    auto tr_0_0_y_x_xxxx = pbuffer.data(idx_op_geom_010_pg + 45);

    auto tr_0_0_y_x_xxxy = pbuffer.data(idx_op_geom_010_pg + 46);

    auto tr_0_0_y_x_xxxz = pbuffer.data(idx_op_geom_010_pg + 47);

    auto tr_0_0_y_x_xxyy = pbuffer.data(idx_op_geom_010_pg + 48);

    auto tr_0_0_y_x_xxyz = pbuffer.data(idx_op_geom_010_pg + 49);

    auto tr_0_0_y_x_xxzz = pbuffer.data(idx_op_geom_010_pg + 50);

    auto tr_0_0_y_x_xyyy = pbuffer.data(idx_op_geom_010_pg + 51);

    auto tr_0_0_y_x_xyyz = pbuffer.data(idx_op_geom_010_pg + 52);

    auto tr_0_0_y_x_xyzz = pbuffer.data(idx_op_geom_010_pg + 53);

    auto tr_0_0_y_x_xzzz = pbuffer.data(idx_op_geom_010_pg + 54);

    auto tr_0_0_y_x_yyyy = pbuffer.data(idx_op_geom_010_pg + 55);

    auto tr_0_0_y_x_yyyz = pbuffer.data(idx_op_geom_010_pg + 56);

    auto tr_0_0_y_x_yyzz = pbuffer.data(idx_op_geom_010_pg + 57);

    auto tr_0_0_y_x_yzzz = pbuffer.data(idx_op_geom_010_pg + 58);

    auto tr_0_0_y_x_zzzz = pbuffer.data(idx_op_geom_010_pg + 59);

    #pragma omp simd aligned(tr_0_0_y_x_xxxx, tr_0_0_y_x_xxxy, tr_0_0_y_x_xxxz, tr_0_0_y_x_xxyy, tr_0_0_y_x_xxyz, tr_0_0_y_x_xxzz, tr_0_0_y_x_xyyy, tr_0_0_y_x_xyyz, tr_0_0_y_x_xyzz, tr_0_0_y_x_xzzz, tr_0_0_y_x_yyyy, tr_0_0_y_x_yyyz, tr_0_0_y_x_yyzz, tr_0_0_y_x_yzzz, tr_0_0_y_x_zzzz, tr_x_xxx, tr_x_xxxxy, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_yyy, tr_x_yyyyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_zzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_x_xxxx[i] = 2.0 * tr_xy_xxxx[i] * tbe_0 + 2.0 * tr_x_xxxxy[i] * tke_0;

        tr_0_0_y_x_xxxy[i] = 2.0 * tr_xy_xxxy[i] * tbe_0 + 2.0 * tr_x_xxxyy[i] * tke_0 - tr_x_xxx[i];

        tr_0_0_y_x_xxxz[i] = 2.0 * tr_xy_xxxz[i] * tbe_0 + 2.0 * tr_x_xxxyz[i] * tke_0;

        tr_0_0_y_x_xxyy[i] = 2.0 * tr_xy_xxyy[i] * tbe_0 + 2.0 * tr_x_xxyyy[i] * tke_0 - 2.0 * tr_x_xxy[i];

        tr_0_0_y_x_xxyz[i] = 2.0 * tr_xy_xxyz[i] * tbe_0 + 2.0 * tr_x_xxyyz[i] * tke_0 - tr_x_xxz[i];

        tr_0_0_y_x_xxzz[i] = 2.0 * tr_xy_xxzz[i] * tbe_0 + 2.0 * tr_x_xxyzz[i] * tke_0;

        tr_0_0_y_x_xyyy[i] = 2.0 * tr_xy_xyyy[i] * tbe_0 + 2.0 * tr_x_xyyyy[i] * tke_0 - 3.0 * tr_x_xyy[i];

        tr_0_0_y_x_xyyz[i] = 2.0 * tr_xy_xyyz[i] * tbe_0 + 2.0 * tr_x_xyyyz[i] * tke_0 - 2.0 * tr_x_xyz[i];

        tr_0_0_y_x_xyzz[i] = 2.0 * tr_xy_xyzz[i] * tbe_0 + 2.0 * tr_x_xyyzz[i] * tke_0 - tr_x_xzz[i];

        tr_0_0_y_x_xzzz[i] = 2.0 * tr_xy_xzzz[i] * tbe_0 + 2.0 * tr_x_xyzzz[i] * tke_0;

        tr_0_0_y_x_yyyy[i] = 2.0 * tr_xy_yyyy[i] * tbe_0 + 2.0 * tr_x_yyyyy[i] * tke_0 - 4.0 * tr_x_yyy[i];

        tr_0_0_y_x_yyyz[i] = 2.0 * tr_xy_yyyz[i] * tbe_0 + 2.0 * tr_x_yyyyz[i] * tke_0 - 3.0 * tr_x_yyz[i];

        tr_0_0_y_x_yyzz[i] = 2.0 * tr_xy_yyzz[i] * tbe_0 + 2.0 * tr_x_yyyzz[i] * tke_0 - 2.0 * tr_x_yzz[i];

        tr_0_0_y_x_yzzz[i] = 2.0 * tr_xy_yzzz[i] * tbe_0 + 2.0 * tr_x_yyzzz[i] * tke_0 - tr_x_zzz[i];

        tr_0_0_y_x_zzzz[i] = 2.0 * tr_xy_zzzz[i] * tbe_0 + 2.0 * tr_x_yzzzz[i] * tke_0;
    }

    // Set up 60-75 components of targeted buffer : PG

    auto tr_0_0_y_y_xxxx = pbuffer.data(idx_op_geom_010_pg + 60);

    auto tr_0_0_y_y_xxxy = pbuffer.data(idx_op_geom_010_pg + 61);

    auto tr_0_0_y_y_xxxz = pbuffer.data(idx_op_geom_010_pg + 62);

    auto tr_0_0_y_y_xxyy = pbuffer.data(idx_op_geom_010_pg + 63);

    auto tr_0_0_y_y_xxyz = pbuffer.data(idx_op_geom_010_pg + 64);

    auto tr_0_0_y_y_xxzz = pbuffer.data(idx_op_geom_010_pg + 65);

    auto tr_0_0_y_y_xyyy = pbuffer.data(idx_op_geom_010_pg + 66);

    auto tr_0_0_y_y_xyyz = pbuffer.data(idx_op_geom_010_pg + 67);

    auto tr_0_0_y_y_xyzz = pbuffer.data(idx_op_geom_010_pg + 68);

    auto tr_0_0_y_y_xzzz = pbuffer.data(idx_op_geom_010_pg + 69);

    auto tr_0_0_y_y_yyyy = pbuffer.data(idx_op_geom_010_pg + 70);

    auto tr_0_0_y_y_yyyz = pbuffer.data(idx_op_geom_010_pg + 71);

    auto tr_0_0_y_y_yyzz = pbuffer.data(idx_op_geom_010_pg + 72);

    auto tr_0_0_y_y_yzzz = pbuffer.data(idx_op_geom_010_pg + 73);

    auto tr_0_0_y_y_zzzz = pbuffer.data(idx_op_geom_010_pg + 74);

    #pragma omp simd aligned(tr_0_0_y_y_xxxx, tr_0_0_y_y_xxxy, tr_0_0_y_y_xxxz, tr_0_0_y_y_xxyy, tr_0_0_y_y_xxyz, tr_0_0_y_y_xxzz, tr_0_0_y_y_xyyy, tr_0_0_y_y_xyyz, tr_0_0_y_y_xyzz, tr_0_0_y_y_xzzz, tr_0_0_y_y_yyyy, tr_0_0_y_y_yyyz, tr_0_0_y_y_yyzz, tr_0_0_y_y_yzzz, tr_0_0_y_y_zzzz, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_y_xxx, tr_y_xxxxy, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_yyy, tr_y_yyyyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_zzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_y_xxxx[i] = 2.0 * tr_yy_xxxx[i] * tbe_0 + 2.0 * tr_y_xxxxy[i] * tke_0 - tr_0_xxxx[i];

        tr_0_0_y_y_xxxy[i] = 2.0 * tr_yy_xxxy[i] * tbe_0 + 2.0 * tr_y_xxxyy[i] * tke_0 - tr_0_xxxy[i] - tr_y_xxx[i];

        tr_0_0_y_y_xxxz[i] = 2.0 * tr_yy_xxxz[i] * tbe_0 + 2.0 * tr_y_xxxyz[i] * tke_0 - tr_0_xxxz[i];

        tr_0_0_y_y_xxyy[i] = 2.0 * tr_yy_xxyy[i] * tbe_0 + 2.0 * tr_y_xxyyy[i] * tke_0 - tr_0_xxyy[i] - 2.0 * tr_y_xxy[i];

        tr_0_0_y_y_xxyz[i] = 2.0 * tr_yy_xxyz[i] * tbe_0 + 2.0 * tr_y_xxyyz[i] * tke_0 - tr_0_xxyz[i] - tr_y_xxz[i];

        tr_0_0_y_y_xxzz[i] = 2.0 * tr_yy_xxzz[i] * tbe_0 + 2.0 * tr_y_xxyzz[i] * tke_0 - tr_0_xxzz[i];

        tr_0_0_y_y_xyyy[i] = 2.0 * tr_yy_xyyy[i] * tbe_0 + 2.0 * tr_y_xyyyy[i] * tke_0 - tr_0_xyyy[i] - 3.0 * tr_y_xyy[i];

        tr_0_0_y_y_xyyz[i] = 2.0 * tr_yy_xyyz[i] * tbe_0 + 2.0 * tr_y_xyyyz[i] * tke_0 - tr_0_xyyz[i] - 2.0 * tr_y_xyz[i];

        tr_0_0_y_y_xyzz[i] = 2.0 * tr_yy_xyzz[i] * tbe_0 + 2.0 * tr_y_xyyzz[i] * tke_0 - tr_0_xyzz[i] - tr_y_xzz[i];

        tr_0_0_y_y_xzzz[i] = 2.0 * tr_yy_xzzz[i] * tbe_0 + 2.0 * tr_y_xyzzz[i] * tke_0 - tr_0_xzzz[i];

        tr_0_0_y_y_yyyy[i] = 2.0 * tr_yy_yyyy[i] * tbe_0 + 2.0 * tr_y_yyyyy[i] * tke_0 - tr_0_yyyy[i] - 4.0 * tr_y_yyy[i];

        tr_0_0_y_y_yyyz[i] = 2.0 * tr_yy_yyyz[i] * tbe_0 + 2.0 * tr_y_yyyyz[i] * tke_0 - tr_0_yyyz[i] - 3.0 * tr_y_yyz[i];

        tr_0_0_y_y_yyzz[i] = 2.0 * tr_yy_yyzz[i] * tbe_0 + 2.0 * tr_y_yyyzz[i] * tke_0 - tr_0_yyzz[i] - 2.0 * tr_y_yzz[i];

        tr_0_0_y_y_yzzz[i] = 2.0 * tr_yy_yzzz[i] * tbe_0 + 2.0 * tr_y_yyzzz[i] * tke_0 - tr_0_yzzz[i] - tr_y_zzz[i];

        tr_0_0_y_y_zzzz[i] = 2.0 * tr_yy_zzzz[i] * tbe_0 + 2.0 * tr_y_yzzzz[i] * tke_0 - tr_0_zzzz[i];
    }

    // Set up 75-90 components of targeted buffer : PG

    auto tr_0_0_y_z_xxxx = pbuffer.data(idx_op_geom_010_pg + 75);

    auto tr_0_0_y_z_xxxy = pbuffer.data(idx_op_geom_010_pg + 76);

    auto tr_0_0_y_z_xxxz = pbuffer.data(idx_op_geom_010_pg + 77);

    auto tr_0_0_y_z_xxyy = pbuffer.data(idx_op_geom_010_pg + 78);

    auto tr_0_0_y_z_xxyz = pbuffer.data(idx_op_geom_010_pg + 79);

    auto tr_0_0_y_z_xxzz = pbuffer.data(idx_op_geom_010_pg + 80);

    auto tr_0_0_y_z_xyyy = pbuffer.data(idx_op_geom_010_pg + 81);

    auto tr_0_0_y_z_xyyz = pbuffer.data(idx_op_geom_010_pg + 82);

    auto tr_0_0_y_z_xyzz = pbuffer.data(idx_op_geom_010_pg + 83);

    auto tr_0_0_y_z_xzzz = pbuffer.data(idx_op_geom_010_pg + 84);

    auto tr_0_0_y_z_yyyy = pbuffer.data(idx_op_geom_010_pg + 85);

    auto tr_0_0_y_z_yyyz = pbuffer.data(idx_op_geom_010_pg + 86);

    auto tr_0_0_y_z_yyzz = pbuffer.data(idx_op_geom_010_pg + 87);

    auto tr_0_0_y_z_yzzz = pbuffer.data(idx_op_geom_010_pg + 88);

    auto tr_0_0_y_z_zzzz = pbuffer.data(idx_op_geom_010_pg + 89);

    #pragma omp simd aligned(tr_0_0_y_z_xxxx, tr_0_0_y_z_xxxy, tr_0_0_y_z_xxxz, tr_0_0_y_z_xxyy, tr_0_0_y_z_xxyz, tr_0_0_y_z_xxzz, tr_0_0_y_z_xyyy, tr_0_0_y_z_xyyz, tr_0_0_y_z_xyzz, tr_0_0_y_z_xzzz, tr_0_0_y_z_yyyy, tr_0_0_y_z_yyyz, tr_0_0_y_z_yyzz, tr_0_0_y_z_yzzz, tr_0_0_y_z_zzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_z_xxx, tr_z_xxxxy, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_yyy, tr_z_yyyyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_z_xxxx[i] = 2.0 * tr_yz_xxxx[i] * tbe_0 + 2.0 * tr_z_xxxxy[i] * tke_0;

        tr_0_0_y_z_xxxy[i] = 2.0 * tr_yz_xxxy[i] * tbe_0 + 2.0 * tr_z_xxxyy[i] * tke_0 - tr_z_xxx[i];

        tr_0_0_y_z_xxxz[i] = 2.0 * tr_yz_xxxz[i] * tbe_0 + 2.0 * tr_z_xxxyz[i] * tke_0;

        tr_0_0_y_z_xxyy[i] = 2.0 * tr_yz_xxyy[i] * tbe_0 + 2.0 * tr_z_xxyyy[i] * tke_0 - 2.0 * tr_z_xxy[i];

        tr_0_0_y_z_xxyz[i] = 2.0 * tr_yz_xxyz[i] * tbe_0 + 2.0 * tr_z_xxyyz[i] * tke_0 - tr_z_xxz[i];

        tr_0_0_y_z_xxzz[i] = 2.0 * tr_yz_xxzz[i] * tbe_0 + 2.0 * tr_z_xxyzz[i] * tke_0;

        tr_0_0_y_z_xyyy[i] = 2.0 * tr_yz_xyyy[i] * tbe_0 + 2.0 * tr_z_xyyyy[i] * tke_0 - 3.0 * tr_z_xyy[i];

        tr_0_0_y_z_xyyz[i] = 2.0 * tr_yz_xyyz[i] * tbe_0 + 2.0 * tr_z_xyyyz[i] * tke_0 - 2.0 * tr_z_xyz[i];

        tr_0_0_y_z_xyzz[i] = 2.0 * tr_yz_xyzz[i] * tbe_0 + 2.0 * tr_z_xyyzz[i] * tke_0 - tr_z_xzz[i];

        tr_0_0_y_z_xzzz[i] = 2.0 * tr_yz_xzzz[i] * tbe_0 + 2.0 * tr_z_xyzzz[i] * tke_0;

        tr_0_0_y_z_yyyy[i] = 2.0 * tr_yz_yyyy[i] * tbe_0 + 2.0 * tr_z_yyyyy[i] * tke_0 - 4.0 * tr_z_yyy[i];

        tr_0_0_y_z_yyyz[i] = 2.0 * tr_yz_yyyz[i] * tbe_0 + 2.0 * tr_z_yyyyz[i] * tke_0 - 3.0 * tr_z_yyz[i];

        tr_0_0_y_z_yyzz[i] = 2.0 * tr_yz_yyzz[i] * tbe_0 + 2.0 * tr_z_yyyzz[i] * tke_0 - 2.0 * tr_z_yzz[i];

        tr_0_0_y_z_yzzz[i] = 2.0 * tr_yz_yzzz[i] * tbe_0 + 2.0 * tr_z_yyzzz[i] * tke_0 - tr_z_zzz[i];

        tr_0_0_y_z_zzzz[i] = 2.0 * tr_yz_zzzz[i] * tbe_0 + 2.0 * tr_z_yzzzz[i] * tke_0;
    }

    // Set up 90-105 components of targeted buffer : PG

    auto tr_0_0_z_x_xxxx = pbuffer.data(idx_op_geom_010_pg + 90);

    auto tr_0_0_z_x_xxxy = pbuffer.data(idx_op_geom_010_pg + 91);

    auto tr_0_0_z_x_xxxz = pbuffer.data(idx_op_geom_010_pg + 92);

    auto tr_0_0_z_x_xxyy = pbuffer.data(idx_op_geom_010_pg + 93);

    auto tr_0_0_z_x_xxyz = pbuffer.data(idx_op_geom_010_pg + 94);

    auto tr_0_0_z_x_xxzz = pbuffer.data(idx_op_geom_010_pg + 95);

    auto tr_0_0_z_x_xyyy = pbuffer.data(idx_op_geom_010_pg + 96);

    auto tr_0_0_z_x_xyyz = pbuffer.data(idx_op_geom_010_pg + 97);

    auto tr_0_0_z_x_xyzz = pbuffer.data(idx_op_geom_010_pg + 98);

    auto tr_0_0_z_x_xzzz = pbuffer.data(idx_op_geom_010_pg + 99);

    auto tr_0_0_z_x_yyyy = pbuffer.data(idx_op_geom_010_pg + 100);

    auto tr_0_0_z_x_yyyz = pbuffer.data(idx_op_geom_010_pg + 101);

    auto tr_0_0_z_x_yyzz = pbuffer.data(idx_op_geom_010_pg + 102);

    auto tr_0_0_z_x_yzzz = pbuffer.data(idx_op_geom_010_pg + 103);

    auto tr_0_0_z_x_zzzz = pbuffer.data(idx_op_geom_010_pg + 104);

    #pragma omp simd aligned(tr_0_0_z_x_xxxx, tr_0_0_z_x_xxxy, tr_0_0_z_x_xxxz, tr_0_0_z_x_xxyy, tr_0_0_z_x_xxyz, tr_0_0_z_x_xxzz, tr_0_0_z_x_xyyy, tr_0_0_z_x_xyyz, tr_0_0_z_x_xyzz, tr_0_0_z_x_xzzz, tr_0_0_z_x_yyyy, tr_0_0_z_x_yyyz, tr_0_0_z_x_yyzz, tr_0_0_z_x_yzzz, tr_0_0_z_x_zzzz, tr_x_xxx, tr_x_xxxxz, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_yyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_zzz, tr_x_zzzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_x_xxxx[i] = 2.0 * tr_xz_xxxx[i] * tbe_0 + 2.0 * tr_x_xxxxz[i] * tke_0;

        tr_0_0_z_x_xxxy[i] = 2.0 * tr_xz_xxxy[i] * tbe_0 + 2.0 * tr_x_xxxyz[i] * tke_0;

        tr_0_0_z_x_xxxz[i] = 2.0 * tr_xz_xxxz[i] * tbe_0 + 2.0 * tr_x_xxxzz[i] * tke_0 - tr_x_xxx[i];

        tr_0_0_z_x_xxyy[i] = 2.0 * tr_xz_xxyy[i] * tbe_0 + 2.0 * tr_x_xxyyz[i] * tke_0;

        tr_0_0_z_x_xxyz[i] = 2.0 * tr_xz_xxyz[i] * tbe_0 + 2.0 * tr_x_xxyzz[i] * tke_0 - tr_x_xxy[i];

        tr_0_0_z_x_xxzz[i] = 2.0 * tr_xz_xxzz[i] * tbe_0 + 2.0 * tr_x_xxzzz[i] * tke_0 - 2.0 * tr_x_xxz[i];

        tr_0_0_z_x_xyyy[i] = 2.0 * tr_xz_xyyy[i] * tbe_0 + 2.0 * tr_x_xyyyz[i] * tke_0;

        tr_0_0_z_x_xyyz[i] = 2.0 * tr_xz_xyyz[i] * tbe_0 + 2.0 * tr_x_xyyzz[i] * tke_0 - tr_x_xyy[i];

        tr_0_0_z_x_xyzz[i] = 2.0 * tr_xz_xyzz[i] * tbe_0 + 2.0 * tr_x_xyzzz[i] * tke_0 - 2.0 * tr_x_xyz[i];

        tr_0_0_z_x_xzzz[i] = 2.0 * tr_xz_xzzz[i] * tbe_0 + 2.0 * tr_x_xzzzz[i] * tke_0 - 3.0 * tr_x_xzz[i];

        tr_0_0_z_x_yyyy[i] = 2.0 * tr_xz_yyyy[i] * tbe_0 + 2.0 * tr_x_yyyyz[i] * tke_0;

        tr_0_0_z_x_yyyz[i] = 2.0 * tr_xz_yyyz[i] * tbe_0 + 2.0 * tr_x_yyyzz[i] * tke_0 - tr_x_yyy[i];

        tr_0_0_z_x_yyzz[i] = 2.0 * tr_xz_yyzz[i] * tbe_0 + 2.0 * tr_x_yyzzz[i] * tke_0 - 2.0 * tr_x_yyz[i];

        tr_0_0_z_x_yzzz[i] = 2.0 * tr_xz_yzzz[i] * tbe_0 + 2.0 * tr_x_yzzzz[i] * tke_0 - 3.0 * tr_x_yzz[i];

        tr_0_0_z_x_zzzz[i] = 2.0 * tr_xz_zzzz[i] * tbe_0 + 2.0 * tr_x_zzzzz[i] * tke_0 - 4.0 * tr_x_zzz[i];
    }

    // Set up 105-120 components of targeted buffer : PG

    auto tr_0_0_z_y_xxxx = pbuffer.data(idx_op_geom_010_pg + 105);

    auto tr_0_0_z_y_xxxy = pbuffer.data(idx_op_geom_010_pg + 106);

    auto tr_0_0_z_y_xxxz = pbuffer.data(idx_op_geom_010_pg + 107);

    auto tr_0_0_z_y_xxyy = pbuffer.data(idx_op_geom_010_pg + 108);

    auto tr_0_0_z_y_xxyz = pbuffer.data(idx_op_geom_010_pg + 109);

    auto tr_0_0_z_y_xxzz = pbuffer.data(idx_op_geom_010_pg + 110);

    auto tr_0_0_z_y_xyyy = pbuffer.data(idx_op_geom_010_pg + 111);

    auto tr_0_0_z_y_xyyz = pbuffer.data(idx_op_geom_010_pg + 112);

    auto tr_0_0_z_y_xyzz = pbuffer.data(idx_op_geom_010_pg + 113);

    auto tr_0_0_z_y_xzzz = pbuffer.data(idx_op_geom_010_pg + 114);

    auto tr_0_0_z_y_yyyy = pbuffer.data(idx_op_geom_010_pg + 115);

    auto tr_0_0_z_y_yyyz = pbuffer.data(idx_op_geom_010_pg + 116);

    auto tr_0_0_z_y_yyzz = pbuffer.data(idx_op_geom_010_pg + 117);

    auto tr_0_0_z_y_yzzz = pbuffer.data(idx_op_geom_010_pg + 118);

    auto tr_0_0_z_y_zzzz = pbuffer.data(idx_op_geom_010_pg + 119);

    #pragma omp simd aligned(tr_0_0_z_y_xxxx, tr_0_0_z_y_xxxy, tr_0_0_z_y_xxxz, tr_0_0_z_y_xxyy, tr_0_0_z_y_xxyz, tr_0_0_z_y_xxzz, tr_0_0_z_y_xyyy, tr_0_0_z_y_xyyz, tr_0_0_z_y_xyzz, tr_0_0_z_y_xzzz, tr_0_0_z_y_yyyy, tr_0_0_z_y_yyyz, tr_0_0_z_y_yyzz, tr_0_0_z_y_yzzz, tr_0_0_z_y_zzzz, tr_y_xxx, tr_y_xxxxz, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_yyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_zzz, tr_y_zzzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_y_xxxx[i] = 2.0 * tr_yz_xxxx[i] * tbe_0 + 2.0 * tr_y_xxxxz[i] * tke_0;

        tr_0_0_z_y_xxxy[i] = 2.0 * tr_yz_xxxy[i] * tbe_0 + 2.0 * tr_y_xxxyz[i] * tke_0;

        tr_0_0_z_y_xxxz[i] = 2.0 * tr_yz_xxxz[i] * tbe_0 + 2.0 * tr_y_xxxzz[i] * tke_0 - tr_y_xxx[i];

        tr_0_0_z_y_xxyy[i] = 2.0 * tr_yz_xxyy[i] * tbe_0 + 2.0 * tr_y_xxyyz[i] * tke_0;

        tr_0_0_z_y_xxyz[i] = 2.0 * tr_yz_xxyz[i] * tbe_0 + 2.0 * tr_y_xxyzz[i] * tke_0 - tr_y_xxy[i];

        tr_0_0_z_y_xxzz[i] = 2.0 * tr_yz_xxzz[i] * tbe_0 + 2.0 * tr_y_xxzzz[i] * tke_0 - 2.0 * tr_y_xxz[i];

        tr_0_0_z_y_xyyy[i] = 2.0 * tr_yz_xyyy[i] * tbe_0 + 2.0 * tr_y_xyyyz[i] * tke_0;

        tr_0_0_z_y_xyyz[i] = 2.0 * tr_yz_xyyz[i] * tbe_0 + 2.0 * tr_y_xyyzz[i] * tke_0 - tr_y_xyy[i];

        tr_0_0_z_y_xyzz[i] = 2.0 * tr_yz_xyzz[i] * tbe_0 + 2.0 * tr_y_xyzzz[i] * tke_0 - 2.0 * tr_y_xyz[i];

        tr_0_0_z_y_xzzz[i] = 2.0 * tr_yz_xzzz[i] * tbe_0 + 2.0 * tr_y_xzzzz[i] * tke_0 - 3.0 * tr_y_xzz[i];

        tr_0_0_z_y_yyyy[i] = 2.0 * tr_yz_yyyy[i] * tbe_0 + 2.0 * tr_y_yyyyz[i] * tke_0;

        tr_0_0_z_y_yyyz[i] = 2.0 * tr_yz_yyyz[i] * tbe_0 + 2.0 * tr_y_yyyzz[i] * tke_0 - tr_y_yyy[i];

        tr_0_0_z_y_yyzz[i] = 2.0 * tr_yz_yyzz[i] * tbe_0 + 2.0 * tr_y_yyzzz[i] * tke_0 - 2.0 * tr_y_yyz[i];

        tr_0_0_z_y_yzzz[i] = 2.0 * tr_yz_yzzz[i] * tbe_0 + 2.0 * tr_y_yzzzz[i] * tke_0 - 3.0 * tr_y_yzz[i];

        tr_0_0_z_y_zzzz[i] = 2.0 * tr_yz_zzzz[i] * tbe_0 + 2.0 * tr_y_zzzzz[i] * tke_0 - 4.0 * tr_y_zzz[i];
    }

    // Set up 120-135 components of targeted buffer : PG

    auto tr_0_0_z_z_xxxx = pbuffer.data(idx_op_geom_010_pg + 120);

    auto tr_0_0_z_z_xxxy = pbuffer.data(idx_op_geom_010_pg + 121);

    auto tr_0_0_z_z_xxxz = pbuffer.data(idx_op_geom_010_pg + 122);

    auto tr_0_0_z_z_xxyy = pbuffer.data(idx_op_geom_010_pg + 123);

    auto tr_0_0_z_z_xxyz = pbuffer.data(idx_op_geom_010_pg + 124);

    auto tr_0_0_z_z_xxzz = pbuffer.data(idx_op_geom_010_pg + 125);

    auto tr_0_0_z_z_xyyy = pbuffer.data(idx_op_geom_010_pg + 126);

    auto tr_0_0_z_z_xyyz = pbuffer.data(idx_op_geom_010_pg + 127);

    auto tr_0_0_z_z_xyzz = pbuffer.data(idx_op_geom_010_pg + 128);

    auto tr_0_0_z_z_xzzz = pbuffer.data(idx_op_geom_010_pg + 129);

    auto tr_0_0_z_z_yyyy = pbuffer.data(idx_op_geom_010_pg + 130);

    auto tr_0_0_z_z_yyyz = pbuffer.data(idx_op_geom_010_pg + 131);

    auto tr_0_0_z_z_yyzz = pbuffer.data(idx_op_geom_010_pg + 132);

    auto tr_0_0_z_z_yzzz = pbuffer.data(idx_op_geom_010_pg + 133);

    auto tr_0_0_z_z_zzzz = pbuffer.data(idx_op_geom_010_pg + 134);

    #pragma omp simd aligned(tr_0_0_z_z_xxxx, tr_0_0_z_z_xxxy, tr_0_0_z_z_xxxz, tr_0_0_z_z_xxyy, tr_0_0_z_z_xxyz, tr_0_0_z_z_xxzz, tr_0_0_z_z_xyyy, tr_0_0_z_z_xyyz, tr_0_0_z_z_xyzz, tr_0_0_z_z_xzzz, tr_0_0_z_z_yyyy, tr_0_0_z_z_yyyz, tr_0_0_z_z_yyzz, tr_0_0_z_z_yzzz, tr_0_0_z_z_zzzz, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_z_xxx, tr_z_xxxxz, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_yyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_zzz, tr_z_zzzzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_z_xxxx[i] = 2.0 * tr_zz_xxxx[i] * tbe_0 + 2.0 * tr_z_xxxxz[i] * tke_0 - tr_0_xxxx[i];

        tr_0_0_z_z_xxxy[i] = 2.0 * tr_zz_xxxy[i] * tbe_0 + 2.0 * tr_z_xxxyz[i] * tke_0 - tr_0_xxxy[i];

        tr_0_0_z_z_xxxz[i] = 2.0 * tr_zz_xxxz[i] * tbe_0 + 2.0 * tr_z_xxxzz[i] * tke_0 - tr_0_xxxz[i] - tr_z_xxx[i];

        tr_0_0_z_z_xxyy[i] = 2.0 * tr_zz_xxyy[i] * tbe_0 + 2.0 * tr_z_xxyyz[i] * tke_0 - tr_0_xxyy[i];

        tr_0_0_z_z_xxyz[i] = 2.0 * tr_zz_xxyz[i] * tbe_0 + 2.0 * tr_z_xxyzz[i] * tke_0 - tr_0_xxyz[i] - tr_z_xxy[i];

        tr_0_0_z_z_xxzz[i] = 2.0 * tr_zz_xxzz[i] * tbe_0 + 2.0 * tr_z_xxzzz[i] * tke_0 - tr_0_xxzz[i] - 2.0 * tr_z_xxz[i];

        tr_0_0_z_z_xyyy[i] = 2.0 * tr_zz_xyyy[i] * tbe_0 + 2.0 * tr_z_xyyyz[i] * tke_0 - tr_0_xyyy[i];

        tr_0_0_z_z_xyyz[i] = 2.0 * tr_zz_xyyz[i] * tbe_0 + 2.0 * tr_z_xyyzz[i] * tke_0 - tr_0_xyyz[i] - tr_z_xyy[i];

        tr_0_0_z_z_xyzz[i] = 2.0 * tr_zz_xyzz[i] * tbe_0 + 2.0 * tr_z_xyzzz[i] * tke_0 - tr_0_xyzz[i] - 2.0 * tr_z_xyz[i];

        tr_0_0_z_z_xzzz[i] = 2.0 * tr_zz_xzzz[i] * tbe_0 + 2.0 * tr_z_xzzzz[i] * tke_0 - tr_0_xzzz[i] - 3.0 * tr_z_xzz[i];

        tr_0_0_z_z_yyyy[i] = 2.0 * tr_zz_yyyy[i] * tbe_0 + 2.0 * tr_z_yyyyz[i] * tke_0 - tr_0_yyyy[i];

        tr_0_0_z_z_yyyz[i] = 2.0 * tr_zz_yyyz[i] * tbe_0 + 2.0 * tr_z_yyyzz[i] * tke_0 - tr_0_yyyz[i] - tr_z_yyy[i];

        tr_0_0_z_z_yyzz[i] = 2.0 * tr_zz_yyzz[i] * tbe_0 + 2.0 * tr_z_yyzzz[i] * tke_0 - tr_0_yyzz[i] - 2.0 * tr_z_yyz[i];

        tr_0_0_z_z_yzzz[i] = 2.0 * tr_zz_yzzz[i] * tbe_0 + 2.0 * tr_z_yzzzz[i] * tke_0 - tr_0_yzzz[i] - 3.0 * tr_z_yzz[i];

        tr_0_0_z_z_zzzz[i] = 2.0 * tr_zz_zzzz[i] * tbe_0 + 2.0 * tr_z_zzzzz[i] * tke_0 - tr_0_zzzz[i] - 4.0 * tr_z_zzz[i];
    }

}

} // t2cgeom namespace

