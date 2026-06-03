#include "GeometricalDerivatives010ForSG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_sg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_sg,
                         const int idx_op_sf,
                         const int idx_op_sh,
                         const int idx_op_pg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SF

    auto tr_0_xxx = pbuffer.data(idx_op_sf);

    auto tr_0_xxy = pbuffer.data(idx_op_sf + 1);

    auto tr_0_xxz = pbuffer.data(idx_op_sf + 2);

    auto tr_0_xyy = pbuffer.data(idx_op_sf + 3);

    auto tr_0_xyz = pbuffer.data(idx_op_sf + 4);

    auto tr_0_xzz = pbuffer.data(idx_op_sf + 5);

    auto tr_0_yyy = pbuffer.data(idx_op_sf + 6);

    auto tr_0_yyz = pbuffer.data(idx_op_sf + 7);

    auto tr_0_yzz = pbuffer.data(idx_op_sf + 8);

    auto tr_0_zzz = pbuffer.data(idx_op_sf + 9);

    // Set up components of auxiliary buffer : SH

    auto tr_0_xxxxx = pbuffer.data(idx_op_sh);

    auto tr_0_xxxxy = pbuffer.data(idx_op_sh + 1);

    auto tr_0_xxxxz = pbuffer.data(idx_op_sh + 2);

    auto tr_0_xxxyy = pbuffer.data(idx_op_sh + 3);

    auto tr_0_xxxyz = pbuffer.data(idx_op_sh + 4);

    auto tr_0_xxxzz = pbuffer.data(idx_op_sh + 5);

    auto tr_0_xxyyy = pbuffer.data(idx_op_sh + 6);

    auto tr_0_xxyyz = pbuffer.data(idx_op_sh + 7);

    auto tr_0_xxyzz = pbuffer.data(idx_op_sh + 8);

    auto tr_0_xxzzz = pbuffer.data(idx_op_sh + 9);

    auto tr_0_xyyyy = pbuffer.data(idx_op_sh + 10);

    auto tr_0_xyyyz = pbuffer.data(idx_op_sh + 11);

    auto tr_0_xyyzz = pbuffer.data(idx_op_sh + 12);

    auto tr_0_xyzzz = pbuffer.data(idx_op_sh + 13);

    auto tr_0_xzzzz = pbuffer.data(idx_op_sh + 14);

    auto tr_0_yyyyy = pbuffer.data(idx_op_sh + 15);

    auto tr_0_yyyyz = pbuffer.data(idx_op_sh + 16);

    auto tr_0_yyyzz = pbuffer.data(idx_op_sh + 17);

    auto tr_0_yyzzz = pbuffer.data(idx_op_sh + 18);

    auto tr_0_yzzzz = pbuffer.data(idx_op_sh + 19);

    auto tr_0_zzzzz = pbuffer.data(idx_op_sh + 20);

    // Set up components of auxiliary buffer : PG

    auto tr_x_xxxx = pbuffer.data(idx_op_pg);

    auto tr_x_xxxy = pbuffer.data(idx_op_pg + 1);

    auto tr_x_xxxz = pbuffer.data(idx_op_pg + 2);

    auto tr_x_xxyy = pbuffer.data(idx_op_pg + 3);

    auto tr_x_xxyz = pbuffer.data(idx_op_pg + 4);

    auto tr_x_xxzz = pbuffer.data(idx_op_pg + 5);

    auto tr_x_xyyy = pbuffer.data(idx_op_pg + 6);

    auto tr_x_xyyz = pbuffer.data(idx_op_pg + 7);

    auto tr_x_xyzz = pbuffer.data(idx_op_pg + 8);

    auto tr_x_xzzz = pbuffer.data(idx_op_pg + 9);

    auto tr_x_yyyy = pbuffer.data(idx_op_pg + 10);

    auto tr_x_yyyz = pbuffer.data(idx_op_pg + 11);

    auto tr_x_yyzz = pbuffer.data(idx_op_pg + 12);

    auto tr_x_yzzz = pbuffer.data(idx_op_pg + 13);

    auto tr_x_zzzz = pbuffer.data(idx_op_pg + 14);

    auto tr_y_xxxx = pbuffer.data(idx_op_pg + 15);

    auto tr_y_xxxy = pbuffer.data(idx_op_pg + 16);

    auto tr_y_xxxz = pbuffer.data(idx_op_pg + 17);

    auto tr_y_xxyy = pbuffer.data(idx_op_pg + 18);

    auto tr_y_xxyz = pbuffer.data(idx_op_pg + 19);

    auto tr_y_xxzz = pbuffer.data(idx_op_pg + 20);

    auto tr_y_xyyy = pbuffer.data(idx_op_pg + 21);

    auto tr_y_xyyz = pbuffer.data(idx_op_pg + 22);

    auto tr_y_xyzz = pbuffer.data(idx_op_pg + 23);

    auto tr_y_xzzz = pbuffer.data(idx_op_pg + 24);

    auto tr_y_yyyy = pbuffer.data(idx_op_pg + 25);

    auto tr_y_yyyz = pbuffer.data(idx_op_pg + 26);

    auto tr_y_yyzz = pbuffer.data(idx_op_pg + 27);

    auto tr_y_yzzz = pbuffer.data(idx_op_pg + 28);

    auto tr_y_zzzz = pbuffer.data(idx_op_pg + 29);

    auto tr_z_xxxx = pbuffer.data(idx_op_pg + 30);

    auto tr_z_xxxy = pbuffer.data(idx_op_pg + 31);

    auto tr_z_xxxz = pbuffer.data(idx_op_pg + 32);

    auto tr_z_xxyy = pbuffer.data(idx_op_pg + 33);

    auto tr_z_xxyz = pbuffer.data(idx_op_pg + 34);

    auto tr_z_xxzz = pbuffer.data(idx_op_pg + 35);

    auto tr_z_xyyy = pbuffer.data(idx_op_pg + 36);

    auto tr_z_xyyz = pbuffer.data(idx_op_pg + 37);

    auto tr_z_xyzz = pbuffer.data(idx_op_pg + 38);

    auto tr_z_xzzz = pbuffer.data(idx_op_pg + 39);

    auto tr_z_yyyy = pbuffer.data(idx_op_pg + 40);

    auto tr_z_yyyz = pbuffer.data(idx_op_pg + 41);

    auto tr_z_yyzz = pbuffer.data(idx_op_pg + 42);

    auto tr_z_yzzz = pbuffer.data(idx_op_pg + 43);

    auto tr_z_zzzz = pbuffer.data(idx_op_pg + 44);

    // Set up components of targeted buffer : SG

    auto tr_0_0_x_0_xxxx = pbuffer.data(idx_op_geom_010_sg);

    auto tr_0_0_x_0_xxxy = pbuffer.data(idx_op_geom_010_sg + 1);

    auto tr_0_0_x_0_xxxz = pbuffer.data(idx_op_geom_010_sg + 2);

    auto tr_0_0_x_0_xxyy = pbuffer.data(idx_op_geom_010_sg + 3);

    auto tr_0_0_x_0_xxyz = pbuffer.data(idx_op_geom_010_sg + 4);

    auto tr_0_0_x_0_xxzz = pbuffer.data(idx_op_geom_010_sg + 5);

    auto tr_0_0_x_0_xyyy = pbuffer.data(idx_op_geom_010_sg + 6);

    auto tr_0_0_x_0_xyyz = pbuffer.data(idx_op_geom_010_sg + 7);

    auto tr_0_0_x_0_xyzz = pbuffer.data(idx_op_geom_010_sg + 8);

    auto tr_0_0_x_0_xzzz = pbuffer.data(idx_op_geom_010_sg + 9);

    auto tr_0_0_x_0_yyyy = pbuffer.data(idx_op_geom_010_sg + 10);

    auto tr_0_0_x_0_yyyz = pbuffer.data(idx_op_geom_010_sg + 11);

    auto tr_0_0_x_0_yyzz = pbuffer.data(idx_op_geom_010_sg + 12);

    auto tr_0_0_x_0_yzzz = pbuffer.data(idx_op_geom_010_sg + 13);

    auto tr_0_0_x_0_zzzz = pbuffer.data(idx_op_geom_010_sg + 14);

    auto tr_0_0_y_0_xxxx = pbuffer.data(idx_op_geom_010_sg + 15);

    auto tr_0_0_y_0_xxxy = pbuffer.data(idx_op_geom_010_sg + 16);

    auto tr_0_0_y_0_xxxz = pbuffer.data(idx_op_geom_010_sg + 17);

    auto tr_0_0_y_0_xxyy = pbuffer.data(idx_op_geom_010_sg + 18);

    auto tr_0_0_y_0_xxyz = pbuffer.data(idx_op_geom_010_sg + 19);

    auto tr_0_0_y_0_xxzz = pbuffer.data(idx_op_geom_010_sg + 20);

    auto tr_0_0_y_0_xyyy = pbuffer.data(idx_op_geom_010_sg + 21);

    auto tr_0_0_y_0_xyyz = pbuffer.data(idx_op_geom_010_sg + 22);

    auto tr_0_0_y_0_xyzz = pbuffer.data(idx_op_geom_010_sg + 23);

    auto tr_0_0_y_0_xzzz = pbuffer.data(idx_op_geom_010_sg + 24);

    auto tr_0_0_y_0_yyyy = pbuffer.data(idx_op_geom_010_sg + 25);

    auto tr_0_0_y_0_yyyz = pbuffer.data(idx_op_geom_010_sg + 26);

    auto tr_0_0_y_0_yyzz = pbuffer.data(idx_op_geom_010_sg + 27);

    auto tr_0_0_y_0_yzzz = pbuffer.data(idx_op_geom_010_sg + 28);

    auto tr_0_0_y_0_zzzz = pbuffer.data(idx_op_geom_010_sg + 29);

    auto tr_0_0_z_0_xxxx = pbuffer.data(idx_op_geom_010_sg + 30);

    auto tr_0_0_z_0_xxxy = pbuffer.data(idx_op_geom_010_sg + 31);

    auto tr_0_0_z_0_xxxz = pbuffer.data(idx_op_geom_010_sg + 32);

    auto tr_0_0_z_0_xxyy = pbuffer.data(idx_op_geom_010_sg + 33);

    auto tr_0_0_z_0_xxyz = pbuffer.data(idx_op_geom_010_sg + 34);

    auto tr_0_0_z_0_xxzz = pbuffer.data(idx_op_geom_010_sg + 35);

    auto tr_0_0_z_0_xyyy = pbuffer.data(idx_op_geom_010_sg + 36);

    auto tr_0_0_z_0_xyyz = pbuffer.data(idx_op_geom_010_sg + 37);

    auto tr_0_0_z_0_xyzz = pbuffer.data(idx_op_geom_010_sg + 38);

    auto tr_0_0_z_0_xzzz = pbuffer.data(idx_op_geom_010_sg + 39);

    auto tr_0_0_z_0_yyyy = pbuffer.data(idx_op_geom_010_sg + 40);

    auto tr_0_0_z_0_yyyz = pbuffer.data(idx_op_geom_010_sg + 41);

    auto tr_0_0_z_0_yyzz = pbuffer.data(idx_op_geom_010_sg + 42);

    auto tr_0_0_z_0_yzzz = pbuffer.data(idx_op_geom_010_sg + 43);

    auto tr_0_0_z_0_zzzz = pbuffer.data(idx_op_geom_010_sg + 44);

    #pragma omp simd aligned(tr_0_0_x_0_xxxx, tr_0_0_x_0_xxxy, tr_0_0_x_0_xxxz, tr_0_0_x_0_xxyy, tr_0_0_x_0_xxyz, tr_0_0_x_0_xxzz, tr_0_0_x_0_xyyy, tr_0_0_x_0_xyyz, tr_0_0_x_0_xyzz, tr_0_0_x_0_xzzz, tr_0_0_x_0_yyyy, tr_0_0_x_0_yyyz, tr_0_0_x_0_yyzz, tr_0_0_x_0_yzzz, tr_0_0_x_0_zzzz, tr_0_0_y_0_xxxx, tr_0_0_y_0_xxxy, tr_0_0_y_0_xxxz, tr_0_0_y_0_xxyy, tr_0_0_y_0_xxyz, tr_0_0_y_0_xxzz, tr_0_0_y_0_xyyy, tr_0_0_y_0_xyyz, tr_0_0_y_0_xyzz, tr_0_0_y_0_xzzz, tr_0_0_y_0_yyyy, tr_0_0_y_0_yyyz, tr_0_0_y_0_yyzz, tr_0_0_y_0_yzzz, tr_0_0_y_0_zzzz, tr_0_0_z_0_xxxx, tr_0_0_z_0_xxxy, tr_0_0_z_0_xxxz, tr_0_0_z_0_xxyy, tr_0_0_z_0_xxyz, tr_0_0_z_0_xxzz, tr_0_0_z_0_xyyy, tr_0_0_z_0_xyyz, tr_0_0_z_0_xyzz, tr_0_0_z_0_xzzz, tr_0_0_z_0_yyyy, tr_0_0_z_0_yyyz, tr_0_0_z_0_yyzz, tr_0_0_z_0_yzzz, tr_0_0_z_0_zzzz, tr_0_xxx, tr_0_xxxxx, tr_0_xxxxy, tr_0_xxxxz, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyyyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_0_zzzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_0_xxxx[i] = 2.0 * tr_x_xxxx[i] * tbe_0 + 2.0 * tr_0_xxxxx[i] * tke_0 - 4.0 * tr_0_xxx[i];

        tr_0_0_x_0_xxxy[i] = 2.0 * tr_x_xxxy[i] * tbe_0 + 2.0 * tr_0_xxxxy[i] * tke_0 - 3.0 * tr_0_xxy[i];

        tr_0_0_x_0_xxxz[i] = 2.0 * tr_x_xxxz[i] * tbe_0 + 2.0 * tr_0_xxxxz[i] * tke_0 - 3.0 * tr_0_xxz[i];

        tr_0_0_x_0_xxyy[i] = 2.0 * tr_x_xxyy[i] * tbe_0 + 2.0 * tr_0_xxxyy[i] * tke_0 - 2.0 * tr_0_xyy[i];

        tr_0_0_x_0_xxyz[i] = 2.0 * tr_x_xxyz[i] * tbe_0 + 2.0 * tr_0_xxxyz[i] * tke_0 - 2.0 * tr_0_xyz[i];

        tr_0_0_x_0_xxzz[i] = 2.0 * tr_x_xxzz[i] * tbe_0 + 2.0 * tr_0_xxxzz[i] * tke_0 - 2.0 * tr_0_xzz[i];

        tr_0_0_x_0_xyyy[i] = 2.0 * tr_x_xyyy[i] * tbe_0 + 2.0 * tr_0_xxyyy[i] * tke_0 - tr_0_yyy[i];

        tr_0_0_x_0_xyyz[i] = 2.0 * tr_x_xyyz[i] * tbe_0 + 2.0 * tr_0_xxyyz[i] * tke_0 - tr_0_yyz[i];

        tr_0_0_x_0_xyzz[i] = 2.0 * tr_x_xyzz[i] * tbe_0 + 2.0 * tr_0_xxyzz[i] * tke_0 - tr_0_yzz[i];

        tr_0_0_x_0_xzzz[i] = 2.0 * tr_x_xzzz[i] * tbe_0 + 2.0 * tr_0_xxzzz[i] * tke_0 - tr_0_zzz[i];

        tr_0_0_x_0_yyyy[i] = 2.0 * tr_x_yyyy[i] * tbe_0 + 2.0 * tr_0_xyyyy[i] * tke_0;

        tr_0_0_x_0_yyyz[i] = 2.0 * tr_x_yyyz[i] * tbe_0 + 2.0 * tr_0_xyyyz[i] * tke_0;

        tr_0_0_x_0_yyzz[i] = 2.0 * tr_x_yyzz[i] * tbe_0 + 2.0 * tr_0_xyyzz[i] * tke_0;

        tr_0_0_x_0_yzzz[i] = 2.0 * tr_x_yzzz[i] * tbe_0 + 2.0 * tr_0_xyzzz[i] * tke_0;

        tr_0_0_x_0_zzzz[i] = 2.0 * tr_x_zzzz[i] * tbe_0 + 2.0 * tr_0_xzzzz[i] * tke_0;

        tr_0_0_y_0_xxxx[i] = 2.0 * tr_y_xxxx[i] * tbe_0 + 2.0 * tr_0_xxxxy[i] * tke_0;

        tr_0_0_y_0_xxxy[i] = 2.0 * tr_y_xxxy[i] * tbe_0 + 2.0 * tr_0_xxxyy[i] * tke_0 - tr_0_xxx[i];

        tr_0_0_y_0_xxxz[i] = 2.0 * tr_y_xxxz[i] * tbe_0 + 2.0 * tr_0_xxxyz[i] * tke_0;

        tr_0_0_y_0_xxyy[i] = 2.0 * tr_y_xxyy[i] * tbe_0 + 2.0 * tr_0_xxyyy[i] * tke_0 - 2.0 * tr_0_xxy[i];

        tr_0_0_y_0_xxyz[i] = 2.0 * tr_y_xxyz[i] * tbe_0 + 2.0 * tr_0_xxyyz[i] * tke_0 - tr_0_xxz[i];

        tr_0_0_y_0_xxzz[i] = 2.0 * tr_y_xxzz[i] * tbe_0 + 2.0 * tr_0_xxyzz[i] * tke_0;

        tr_0_0_y_0_xyyy[i] = 2.0 * tr_y_xyyy[i] * tbe_0 + 2.0 * tr_0_xyyyy[i] * tke_0 - 3.0 * tr_0_xyy[i];

        tr_0_0_y_0_xyyz[i] = 2.0 * tr_y_xyyz[i] * tbe_0 + 2.0 * tr_0_xyyyz[i] * tke_0 - 2.0 * tr_0_xyz[i];

        tr_0_0_y_0_xyzz[i] = 2.0 * tr_y_xyzz[i] * tbe_0 + 2.0 * tr_0_xyyzz[i] * tke_0 - tr_0_xzz[i];

        tr_0_0_y_0_xzzz[i] = 2.0 * tr_y_xzzz[i] * tbe_0 + 2.0 * tr_0_xyzzz[i] * tke_0;

        tr_0_0_y_0_yyyy[i] = 2.0 * tr_y_yyyy[i] * tbe_0 + 2.0 * tr_0_yyyyy[i] * tke_0 - 4.0 * tr_0_yyy[i];

        tr_0_0_y_0_yyyz[i] = 2.0 * tr_y_yyyz[i] * tbe_0 + 2.0 * tr_0_yyyyz[i] * tke_0 - 3.0 * tr_0_yyz[i];

        tr_0_0_y_0_yyzz[i] = 2.0 * tr_y_yyzz[i] * tbe_0 + 2.0 * tr_0_yyyzz[i] * tke_0 - 2.0 * tr_0_yzz[i];

        tr_0_0_y_0_yzzz[i] = 2.0 * tr_y_yzzz[i] * tbe_0 + 2.0 * tr_0_yyzzz[i] * tke_0 - tr_0_zzz[i];

        tr_0_0_y_0_zzzz[i] = 2.0 * tr_y_zzzz[i] * tbe_0 + 2.0 * tr_0_yzzzz[i] * tke_0;

        tr_0_0_z_0_xxxx[i] = 2.0 * tr_z_xxxx[i] * tbe_0 + 2.0 * tr_0_xxxxz[i] * tke_0;

        tr_0_0_z_0_xxxy[i] = 2.0 * tr_z_xxxy[i] * tbe_0 + 2.0 * tr_0_xxxyz[i] * tke_0;

        tr_0_0_z_0_xxxz[i] = 2.0 * tr_z_xxxz[i] * tbe_0 + 2.0 * tr_0_xxxzz[i] * tke_0 - tr_0_xxx[i];

        tr_0_0_z_0_xxyy[i] = 2.0 * tr_z_xxyy[i] * tbe_0 + 2.0 * tr_0_xxyyz[i] * tke_0;

        tr_0_0_z_0_xxyz[i] = 2.0 * tr_z_xxyz[i] * tbe_0 + 2.0 * tr_0_xxyzz[i] * tke_0 - tr_0_xxy[i];

        tr_0_0_z_0_xxzz[i] = 2.0 * tr_z_xxzz[i] * tbe_0 + 2.0 * tr_0_xxzzz[i] * tke_0 - 2.0 * tr_0_xxz[i];

        tr_0_0_z_0_xyyy[i] = 2.0 * tr_z_xyyy[i] * tbe_0 + 2.0 * tr_0_xyyyz[i] * tke_0;

        tr_0_0_z_0_xyyz[i] = 2.0 * tr_z_xyyz[i] * tbe_0 + 2.0 * tr_0_xyyzz[i] * tke_0 - tr_0_xyy[i];

        tr_0_0_z_0_xyzz[i] = 2.0 * tr_z_xyzz[i] * tbe_0 + 2.0 * tr_0_xyzzz[i] * tke_0 - 2.0 * tr_0_xyz[i];

        tr_0_0_z_0_xzzz[i] = 2.0 * tr_z_xzzz[i] * tbe_0 + 2.0 * tr_0_xzzzz[i] * tke_0 - 3.0 * tr_0_xzz[i];

        tr_0_0_z_0_yyyy[i] = 2.0 * tr_z_yyyy[i] * tbe_0 + 2.0 * tr_0_yyyyz[i] * tke_0;

        tr_0_0_z_0_yyyz[i] = 2.0 * tr_z_yyyz[i] * tbe_0 + 2.0 * tr_0_yyyzz[i] * tke_0 - tr_0_yyy[i];

        tr_0_0_z_0_yyzz[i] = 2.0 * tr_z_yyzz[i] * tbe_0 + 2.0 * tr_0_yyzzz[i] * tke_0 - 2.0 * tr_0_yyz[i];

        tr_0_0_z_0_yzzz[i] = 2.0 * tr_z_yzzz[i] * tbe_0 + 2.0 * tr_0_yzzzz[i] * tke_0 - 3.0 * tr_0_yzz[i];

        tr_0_0_z_0_zzzz[i] = 2.0 * tr_z_zzzz[i] * tbe_0 + 2.0 * tr_0_zzzzz[i] * tke_0 - 4.0 * tr_0_zzz[i];
    }
}

} // t2cgeom namespace

