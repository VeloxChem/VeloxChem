#include "GeometricalDerivatives0X1ForPG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_pg(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_pg,
                       const int idx_op_pf,
                       const int idx_op_ph,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : PF

        auto to_x_xxx = pbuffer.data(idx_op_pf + i * 30 + 0);

        auto to_x_xxy = pbuffer.data(idx_op_pf + i * 30 + 1);

        auto to_x_xxz = pbuffer.data(idx_op_pf + i * 30 + 2);

        auto to_x_xyy = pbuffer.data(idx_op_pf + i * 30 + 3);

        auto to_x_xyz = pbuffer.data(idx_op_pf + i * 30 + 4);

        auto to_x_xzz = pbuffer.data(idx_op_pf + i * 30 + 5);

        auto to_x_yyy = pbuffer.data(idx_op_pf + i * 30 + 6);

        auto to_x_yyz = pbuffer.data(idx_op_pf + i * 30 + 7);

        auto to_x_yzz = pbuffer.data(idx_op_pf + i * 30 + 8);

        auto to_x_zzz = pbuffer.data(idx_op_pf + i * 30 + 9);

        auto to_y_xxx = pbuffer.data(idx_op_pf + i * 30 + 10);

        auto to_y_xxy = pbuffer.data(idx_op_pf + i * 30 + 11);

        auto to_y_xxz = pbuffer.data(idx_op_pf + i * 30 + 12);

        auto to_y_xyy = pbuffer.data(idx_op_pf + i * 30 + 13);

        auto to_y_xyz = pbuffer.data(idx_op_pf + i * 30 + 14);

        auto to_y_xzz = pbuffer.data(idx_op_pf + i * 30 + 15);

        auto to_y_yyy = pbuffer.data(idx_op_pf + i * 30 + 16);

        auto to_y_yyz = pbuffer.data(idx_op_pf + i * 30 + 17);

        auto to_y_yzz = pbuffer.data(idx_op_pf + i * 30 + 18);

        auto to_y_zzz = pbuffer.data(idx_op_pf + i * 30 + 19);

        auto to_z_xxx = pbuffer.data(idx_op_pf + i * 30 + 20);

        auto to_z_xxy = pbuffer.data(idx_op_pf + i * 30 + 21);

        auto to_z_xxz = pbuffer.data(idx_op_pf + i * 30 + 22);

        auto to_z_xyy = pbuffer.data(idx_op_pf + i * 30 + 23);

        auto to_z_xyz = pbuffer.data(idx_op_pf + i * 30 + 24);

        auto to_z_xzz = pbuffer.data(idx_op_pf + i * 30 + 25);

        auto to_z_yyy = pbuffer.data(idx_op_pf + i * 30 + 26);

        auto to_z_yyz = pbuffer.data(idx_op_pf + i * 30 + 27);

        auto to_z_yzz = pbuffer.data(idx_op_pf + i * 30 + 28);

        auto to_z_zzz = pbuffer.data(idx_op_pf + i * 30 + 29);

        // Set up components of auxiliary buffer : PH

        auto to_x_xxxxx = pbuffer.data(idx_op_ph + i * 63 + 0);

        auto to_x_xxxxy = pbuffer.data(idx_op_ph + i * 63 + 1);

        auto to_x_xxxxz = pbuffer.data(idx_op_ph + i * 63 + 2);

        auto to_x_xxxyy = pbuffer.data(idx_op_ph + i * 63 + 3);

        auto to_x_xxxyz = pbuffer.data(idx_op_ph + i * 63 + 4);

        auto to_x_xxxzz = pbuffer.data(idx_op_ph + i * 63 + 5);

        auto to_x_xxyyy = pbuffer.data(idx_op_ph + i * 63 + 6);

        auto to_x_xxyyz = pbuffer.data(idx_op_ph + i * 63 + 7);

        auto to_x_xxyzz = pbuffer.data(idx_op_ph + i * 63 + 8);

        auto to_x_xxzzz = pbuffer.data(idx_op_ph + i * 63 + 9);

        auto to_x_xyyyy = pbuffer.data(idx_op_ph + i * 63 + 10);

        auto to_x_xyyyz = pbuffer.data(idx_op_ph + i * 63 + 11);

        auto to_x_xyyzz = pbuffer.data(idx_op_ph + i * 63 + 12);

        auto to_x_xyzzz = pbuffer.data(idx_op_ph + i * 63 + 13);

        auto to_x_xzzzz = pbuffer.data(idx_op_ph + i * 63 + 14);

        auto to_x_yyyyy = pbuffer.data(idx_op_ph + i * 63 + 15);

        auto to_x_yyyyz = pbuffer.data(idx_op_ph + i * 63 + 16);

        auto to_x_yyyzz = pbuffer.data(idx_op_ph + i * 63 + 17);

        auto to_x_yyzzz = pbuffer.data(idx_op_ph + i * 63 + 18);

        auto to_x_yzzzz = pbuffer.data(idx_op_ph + i * 63 + 19);

        auto to_x_zzzzz = pbuffer.data(idx_op_ph + i * 63 + 20);

        auto to_y_xxxxx = pbuffer.data(idx_op_ph + i * 63 + 21);

        auto to_y_xxxxy = pbuffer.data(idx_op_ph + i * 63 + 22);

        auto to_y_xxxxz = pbuffer.data(idx_op_ph + i * 63 + 23);

        auto to_y_xxxyy = pbuffer.data(idx_op_ph + i * 63 + 24);

        auto to_y_xxxyz = pbuffer.data(idx_op_ph + i * 63 + 25);

        auto to_y_xxxzz = pbuffer.data(idx_op_ph + i * 63 + 26);

        auto to_y_xxyyy = pbuffer.data(idx_op_ph + i * 63 + 27);

        auto to_y_xxyyz = pbuffer.data(idx_op_ph + i * 63 + 28);

        auto to_y_xxyzz = pbuffer.data(idx_op_ph + i * 63 + 29);

        auto to_y_xxzzz = pbuffer.data(idx_op_ph + i * 63 + 30);

        auto to_y_xyyyy = pbuffer.data(idx_op_ph + i * 63 + 31);

        auto to_y_xyyyz = pbuffer.data(idx_op_ph + i * 63 + 32);

        auto to_y_xyyzz = pbuffer.data(idx_op_ph + i * 63 + 33);

        auto to_y_xyzzz = pbuffer.data(idx_op_ph + i * 63 + 34);

        auto to_y_xzzzz = pbuffer.data(idx_op_ph + i * 63 + 35);

        auto to_y_yyyyy = pbuffer.data(idx_op_ph + i * 63 + 36);

        auto to_y_yyyyz = pbuffer.data(idx_op_ph + i * 63 + 37);

        auto to_y_yyyzz = pbuffer.data(idx_op_ph + i * 63 + 38);

        auto to_y_yyzzz = pbuffer.data(idx_op_ph + i * 63 + 39);

        auto to_y_yzzzz = pbuffer.data(idx_op_ph + i * 63 + 40);

        auto to_y_zzzzz = pbuffer.data(idx_op_ph + i * 63 + 41);

        auto to_z_xxxxx = pbuffer.data(idx_op_ph + i * 63 + 42);

        auto to_z_xxxxy = pbuffer.data(idx_op_ph + i * 63 + 43);

        auto to_z_xxxxz = pbuffer.data(idx_op_ph + i * 63 + 44);

        auto to_z_xxxyy = pbuffer.data(idx_op_ph + i * 63 + 45);

        auto to_z_xxxyz = pbuffer.data(idx_op_ph + i * 63 + 46);

        auto to_z_xxxzz = pbuffer.data(idx_op_ph + i * 63 + 47);

        auto to_z_xxyyy = pbuffer.data(idx_op_ph + i * 63 + 48);

        auto to_z_xxyyz = pbuffer.data(idx_op_ph + i * 63 + 49);

        auto to_z_xxyzz = pbuffer.data(idx_op_ph + i * 63 + 50);

        auto to_z_xxzzz = pbuffer.data(idx_op_ph + i * 63 + 51);

        auto to_z_xyyyy = pbuffer.data(idx_op_ph + i * 63 + 52);

        auto to_z_xyyyz = pbuffer.data(idx_op_ph + i * 63 + 53);

        auto to_z_xyyzz = pbuffer.data(idx_op_ph + i * 63 + 54);

        auto to_z_xyzzz = pbuffer.data(idx_op_ph + i * 63 + 55);

        auto to_z_xzzzz = pbuffer.data(idx_op_ph + i * 63 + 56);

        auto to_z_yyyyy = pbuffer.data(idx_op_ph + i * 63 + 57);

        auto to_z_yyyyz = pbuffer.data(idx_op_ph + i * 63 + 58);

        auto to_z_yyyzz = pbuffer.data(idx_op_ph + i * 63 + 59);

        auto to_z_yyzzz = pbuffer.data(idx_op_ph + i * 63 + 60);

        auto to_z_yzzzz = pbuffer.data(idx_op_ph + i * 63 + 61);

        auto to_z_zzzzz = pbuffer.data(idx_op_ph + i * 63 + 62);

        // Set up 0-15 components of targeted buffer : PG

        auto to_0_x_x_xxxx = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 0);

        auto to_0_x_x_xxxy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 1);

        auto to_0_x_x_xxxz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 2);

        auto to_0_x_x_xxyy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 3);

        auto to_0_x_x_xxyz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 4);

        auto to_0_x_x_xxzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 5);

        auto to_0_x_x_xyyy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 6);

        auto to_0_x_x_xyyz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 7);

        auto to_0_x_x_xyzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 8);

        auto to_0_x_x_xzzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 9);

        auto to_0_x_x_yyyy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 10);

        auto to_0_x_x_yyyz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 11);

        auto to_0_x_x_yyzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 12);

        auto to_0_x_x_yzzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 13);

        auto to_0_x_x_zzzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_0_x_x_xxxx, to_0_x_x_xxxy, to_0_x_x_xxxz, to_0_x_x_xxyy, to_0_x_x_xxyz, to_0_x_x_xxzz, to_0_x_x_xyyy, to_0_x_x_xyyz, to_0_x_x_xyzz, to_0_x_x_xzzz, to_0_x_x_yyyy, to_0_x_x_yyyz, to_0_x_x_yyzz, to_0_x_x_yzzz, to_0_x_x_zzzz, to_x_xxx, to_x_xxxxx, to_x_xxxxy, to_x_xxxxz, to_x_xxxyy, to_x_xxxyz, to_x_xxxzz, to_x_xxy, to_x_xxyyy, to_x_xxyyz, to_x_xxyzz, to_x_xxz, to_x_xxzzz, to_x_xyy, to_x_xyyyy, to_x_xyyyz, to_x_xyyzz, to_x_xyz, to_x_xyzzz, to_x_xzz, to_x_xzzzz, to_x_yyy, to_x_yyz, to_x_yzz, to_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_x_xxxx[k] = -4.0 * to_x_xxx[k] + 2.0 * to_x_xxxxx[k] * tke_0;

            to_0_x_x_xxxy[k] = -3.0 * to_x_xxy[k] + 2.0 * to_x_xxxxy[k] * tke_0;

            to_0_x_x_xxxz[k] = -3.0 * to_x_xxz[k] + 2.0 * to_x_xxxxz[k] * tke_0;

            to_0_x_x_xxyy[k] = -2.0 * to_x_xyy[k] + 2.0 * to_x_xxxyy[k] * tke_0;

            to_0_x_x_xxyz[k] = -2.0 * to_x_xyz[k] + 2.0 * to_x_xxxyz[k] * tke_0;

            to_0_x_x_xxzz[k] = -2.0 * to_x_xzz[k] + 2.0 * to_x_xxxzz[k] * tke_0;

            to_0_x_x_xyyy[k] = -to_x_yyy[k] + 2.0 * to_x_xxyyy[k] * tke_0;

            to_0_x_x_xyyz[k] = -to_x_yyz[k] + 2.0 * to_x_xxyyz[k] * tke_0;

            to_0_x_x_xyzz[k] = -to_x_yzz[k] + 2.0 * to_x_xxyzz[k] * tke_0;

            to_0_x_x_xzzz[k] = -to_x_zzz[k] + 2.0 * to_x_xxzzz[k] * tke_0;

            to_0_x_x_yyyy[k] = 2.0 * to_x_xyyyy[k] * tke_0;

            to_0_x_x_yyyz[k] = 2.0 * to_x_xyyyz[k] * tke_0;

            to_0_x_x_yyzz[k] = 2.0 * to_x_xyyzz[k] * tke_0;

            to_0_x_x_yzzz[k] = 2.0 * to_x_xyzzz[k] * tke_0;

            to_0_x_x_zzzz[k] = 2.0 * to_x_xzzzz[k] * tke_0;
        }

        // Set up 15-30 components of targeted buffer : PG

        auto to_0_x_y_xxxx = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 15);

        auto to_0_x_y_xxxy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 16);

        auto to_0_x_y_xxxz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 17);

        auto to_0_x_y_xxyy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 18);

        auto to_0_x_y_xxyz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 19);

        auto to_0_x_y_xxzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 20);

        auto to_0_x_y_xyyy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 21);

        auto to_0_x_y_xyyz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 22);

        auto to_0_x_y_xyzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 23);

        auto to_0_x_y_xzzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 24);

        auto to_0_x_y_yyyy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 25);

        auto to_0_x_y_yyyz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 26);

        auto to_0_x_y_yyzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 27);

        auto to_0_x_y_yzzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 28);

        auto to_0_x_y_zzzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_0_x_y_xxxx, to_0_x_y_xxxy, to_0_x_y_xxxz, to_0_x_y_xxyy, to_0_x_y_xxyz, to_0_x_y_xxzz, to_0_x_y_xyyy, to_0_x_y_xyyz, to_0_x_y_xyzz, to_0_x_y_xzzz, to_0_x_y_yyyy, to_0_x_y_yyyz, to_0_x_y_yyzz, to_0_x_y_yzzz, to_0_x_y_zzzz, to_y_xxx, to_y_xxxxx, to_y_xxxxy, to_y_xxxxz, to_y_xxxyy, to_y_xxxyz, to_y_xxxzz, to_y_xxy, to_y_xxyyy, to_y_xxyyz, to_y_xxyzz, to_y_xxz, to_y_xxzzz, to_y_xyy, to_y_xyyyy, to_y_xyyyz, to_y_xyyzz, to_y_xyz, to_y_xyzzz, to_y_xzz, to_y_xzzzz, to_y_yyy, to_y_yyz, to_y_yzz, to_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_y_xxxx[k] = -4.0 * to_y_xxx[k] + 2.0 * to_y_xxxxx[k] * tke_0;

            to_0_x_y_xxxy[k] = -3.0 * to_y_xxy[k] + 2.0 * to_y_xxxxy[k] * tke_0;

            to_0_x_y_xxxz[k] = -3.0 * to_y_xxz[k] + 2.0 * to_y_xxxxz[k] * tke_0;

            to_0_x_y_xxyy[k] = -2.0 * to_y_xyy[k] + 2.0 * to_y_xxxyy[k] * tke_0;

            to_0_x_y_xxyz[k] = -2.0 * to_y_xyz[k] + 2.0 * to_y_xxxyz[k] * tke_0;

            to_0_x_y_xxzz[k] = -2.0 * to_y_xzz[k] + 2.0 * to_y_xxxzz[k] * tke_0;

            to_0_x_y_xyyy[k] = -to_y_yyy[k] + 2.0 * to_y_xxyyy[k] * tke_0;

            to_0_x_y_xyyz[k] = -to_y_yyz[k] + 2.0 * to_y_xxyyz[k] * tke_0;

            to_0_x_y_xyzz[k] = -to_y_yzz[k] + 2.0 * to_y_xxyzz[k] * tke_0;

            to_0_x_y_xzzz[k] = -to_y_zzz[k] + 2.0 * to_y_xxzzz[k] * tke_0;

            to_0_x_y_yyyy[k] = 2.0 * to_y_xyyyy[k] * tke_0;

            to_0_x_y_yyyz[k] = 2.0 * to_y_xyyyz[k] * tke_0;

            to_0_x_y_yyzz[k] = 2.0 * to_y_xyyzz[k] * tke_0;

            to_0_x_y_yzzz[k] = 2.0 * to_y_xyzzz[k] * tke_0;

            to_0_x_y_zzzz[k] = 2.0 * to_y_xzzzz[k] * tke_0;
        }

        // Set up 30-45 components of targeted buffer : PG

        auto to_0_x_z_xxxx = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 30);

        auto to_0_x_z_xxxy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 31);

        auto to_0_x_z_xxxz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 32);

        auto to_0_x_z_xxyy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 33);

        auto to_0_x_z_xxyz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 34);

        auto to_0_x_z_xxzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 35);

        auto to_0_x_z_xyyy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 36);

        auto to_0_x_z_xyyz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 37);

        auto to_0_x_z_xyzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 38);

        auto to_0_x_z_xzzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 39);

        auto to_0_x_z_yyyy = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 40);

        auto to_0_x_z_yyyz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 41);

        auto to_0_x_z_yyzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 42);

        auto to_0_x_z_yzzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 43);

        auto to_0_x_z_zzzz = pbuffer.data(idx_op_geom_001_pg + 0 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_0_x_z_xxxx, to_0_x_z_xxxy, to_0_x_z_xxxz, to_0_x_z_xxyy, to_0_x_z_xxyz, to_0_x_z_xxzz, to_0_x_z_xyyy, to_0_x_z_xyyz, to_0_x_z_xyzz, to_0_x_z_xzzz, to_0_x_z_yyyy, to_0_x_z_yyyz, to_0_x_z_yyzz, to_0_x_z_yzzz, to_0_x_z_zzzz, to_z_xxx, to_z_xxxxx, to_z_xxxxy, to_z_xxxxz, to_z_xxxyy, to_z_xxxyz, to_z_xxxzz, to_z_xxy, to_z_xxyyy, to_z_xxyyz, to_z_xxyzz, to_z_xxz, to_z_xxzzz, to_z_xyy, to_z_xyyyy, to_z_xyyyz, to_z_xyyzz, to_z_xyz, to_z_xyzzz, to_z_xzz, to_z_xzzzz, to_z_yyy, to_z_yyz, to_z_yzz, to_z_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_z_xxxx[k] = -4.0 * to_z_xxx[k] + 2.0 * to_z_xxxxx[k] * tke_0;

            to_0_x_z_xxxy[k] = -3.0 * to_z_xxy[k] + 2.0 * to_z_xxxxy[k] * tke_0;

            to_0_x_z_xxxz[k] = -3.0 * to_z_xxz[k] + 2.0 * to_z_xxxxz[k] * tke_0;

            to_0_x_z_xxyy[k] = -2.0 * to_z_xyy[k] + 2.0 * to_z_xxxyy[k] * tke_0;

            to_0_x_z_xxyz[k] = -2.0 * to_z_xyz[k] + 2.0 * to_z_xxxyz[k] * tke_0;

            to_0_x_z_xxzz[k] = -2.0 * to_z_xzz[k] + 2.0 * to_z_xxxzz[k] * tke_0;

            to_0_x_z_xyyy[k] = -to_z_yyy[k] + 2.0 * to_z_xxyyy[k] * tke_0;

            to_0_x_z_xyyz[k] = -to_z_yyz[k] + 2.0 * to_z_xxyyz[k] * tke_0;

            to_0_x_z_xyzz[k] = -to_z_yzz[k] + 2.0 * to_z_xxyzz[k] * tke_0;

            to_0_x_z_xzzz[k] = -to_z_zzz[k] + 2.0 * to_z_xxzzz[k] * tke_0;

            to_0_x_z_yyyy[k] = 2.0 * to_z_xyyyy[k] * tke_0;

            to_0_x_z_yyyz[k] = 2.0 * to_z_xyyyz[k] * tke_0;

            to_0_x_z_yyzz[k] = 2.0 * to_z_xyyzz[k] * tke_0;

            to_0_x_z_yzzz[k] = 2.0 * to_z_xyzzz[k] * tke_0;

            to_0_x_z_zzzz[k] = 2.0 * to_z_xzzzz[k] * tke_0;
        }

        // Set up 45-60 components of targeted buffer : PG

        auto to_0_y_x_xxxx = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 0);

        auto to_0_y_x_xxxy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 1);

        auto to_0_y_x_xxxz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 2);

        auto to_0_y_x_xxyy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 3);

        auto to_0_y_x_xxyz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 4);

        auto to_0_y_x_xxzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 5);

        auto to_0_y_x_xyyy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 6);

        auto to_0_y_x_xyyz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 7);

        auto to_0_y_x_xyzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 8);

        auto to_0_y_x_xzzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 9);

        auto to_0_y_x_yyyy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 10);

        auto to_0_y_x_yyyz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 11);

        auto to_0_y_x_yyzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 12);

        auto to_0_y_x_yzzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 13);

        auto to_0_y_x_zzzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_0_y_x_xxxx, to_0_y_x_xxxy, to_0_y_x_xxxz, to_0_y_x_xxyy, to_0_y_x_xxyz, to_0_y_x_xxzz, to_0_y_x_xyyy, to_0_y_x_xyyz, to_0_y_x_xyzz, to_0_y_x_xzzz, to_0_y_x_yyyy, to_0_y_x_yyyz, to_0_y_x_yyzz, to_0_y_x_yzzz, to_0_y_x_zzzz, to_x_xxx, to_x_xxxxy, to_x_xxxyy, to_x_xxxyz, to_x_xxy, to_x_xxyyy, to_x_xxyyz, to_x_xxyzz, to_x_xxz, to_x_xyy, to_x_xyyyy, to_x_xyyyz, to_x_xyyzz, to_x_xyz, to_x_xyzzz, to_x_xzz, to_x_yyy, to_x_yyyyy, to_x_yyyyz, to_x_yyyzz, to_x_yyz, to_x_yyzzz, to_x_yzz, to_x_yzzzz, to_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_x_xxxx[k] = 2.0 * to_x_xxxxy[k] * tke_0;

            to_0_y_x_xxxy[k] = -to_x_xxx[k] + 2.0 * to_x_xxxyy[k] * tke_0;

            to_0_y_x_xxxz[k] = 2.0 * to_x_xxxyz[k] * tke_0;

            to_0_y_x_xxyy[k] = -2.0 * to_x_xxy[k] + 2.0 * to_x_xxyyy[k] * tke_0;

            to_0_y_x_xxyz[k] = -to_x_xxz[k] + 2.0 * to_x_xxyyz[k] * tke_0;

            to_0_y_x_xxzz[k] = 2.0 * to_x_xxyzz[k] * tke_0;

            to_0_y_x_xyyy[k] = -3.0 * to_x_xyy[k] + 2.0 * to_x_xyyyy[k] * tke_0;

            to_0_y_x_xyyz[k] = -2.0 * to_x_xyz[k] + 2.0 * to_x_xyyyz[k] * tke_0;

            to_0_y_x_xyzz[k] = -to_x_xzz[k] + 2.0 * to_x_xyyzz[k] * tke_0;

            to_0_y_x_xzzz[k] = 2.0 * to_x_xyzzz[k] * tke_0;

            to_0_y_x_yyyy[k] = -4.0 * to_x_yyy[k] + 2.0 * to_x_yyyyy[k] * tke_0;

            to_0_y_x_yyyz[k] = -3.0 * to_x_yyz[k] + 2.0 * to_x_yyyyz[k] * tke_0;

            to_0_y_x_yyzz[k] = -2.0 * to_x_yzz[k] + 2.0 * to_x_yyyzz[k] * tke_0;

            to_0_y_x_yzzz[k] = -to_x_zzz[k] + 2.0 * to_x_yyzzz[k] * tke_0;

            to_0_y_x_zzzz[k] = 2.0 * to_x_yzzzz[k] * tke_0;
        }

        // Set up 60-75 components of targeted buffer : PG

        auto to_0_y_y_xxxx = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 15);

        auto to_0_y_y_xxxy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 16);

        auto to_0_y_y_xxxz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 17);

        auto to_0_y_y_xxyy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 18);

        auto to_0_y_y_xxyz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 19);

        auto to_0_y_y_xxzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 20);

        auto to_0_y_y_xyyy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 21);

        auto to_0_y_y_xyyz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 22);

        auto to_0_y_y_xyzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 23);

        auto to_0_y_y_xzzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 24);

        auto to_0_y_y_yyyy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 25);

        auto to_0_y_y_yyyz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 26);

        auto to_0_y_y_yyzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 27);

        auto to_0_y_y_yzzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 28);

        auto to_0_y_y_zzzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_0_y_y_xxxx, to_0_y_y_xxxy, to_0_y_y_xxxz, to_0_y_y_xxyy, to_0_y_y_xxyz, to_0_y_y_xxzz, to_0_y_y_xyyy, to_0_y_y_xyyz, to_0_y_y_xyzz, to_0_y_y_xzzz, to_0_y_y_yyyy, to_0_y_y_yyyz, to_0_y_y_yyzz, to_0_y_y_yzzz, to_0_y_y_zzzz, to_y_xxx, to_y_xxxxy, to_y_xxxyy, to_y_xxxyz, to_y_xxy, to_y_xxyyy, to_y_xxyyz, to_y_xxyzz, to_y_xxz, to_y_xyy, to_y_xyyyy, to_y_xyyyz, to_y_xyyzz, to_y_xyz, to_y_xyzzz, to_y_xzz, to_y_yyy, to_y_yyyyy, to_y_yyyyz, to_y_yyyzz, to_y_yyz, to_y_yyzzz, to_y_yzz, to_y_yzzzz, to_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_y_xxxx[k] = 2.0 * to_y_xxxxy[k] * tke_0;

            to_0_y_y_xxxy[k] = -to_y_xxx[k] + 2.0 * to_y_xxxyy[k] * tke_0;

            to_0_y_y_xxxz[k] = 2.0 * to_y_xxxyz[k] * tke_0;

            to_0_y_y_xxyy[k] = -2.0 * to_y_xxy[k] + 2.0 * to_y_xxyyy[k] * tke_0;

            to_0_y_y_xxyz[k] = -to_y_xxz[k] + 2.0 * to_y_xxyyz[k] * tke_0;

            to_0_y_y_xxzz[k] = 2.0 * to_y_xxyzz[k] * tke_0;

            to_0_y_y_xyyy[k] = -3.0 * to_y_xyy[k] + 2.0 * to_y_xyyyy[k] * tke_0;

            to_0_y_y_xyyz[k] = -2.0 * to_y_xyz[k] + 2.0 * to_y_xyyyz[k] * tke_0;

            to_0_y_y_xyzz[k] = -to_y_xzz[k] + 2.0 * to_y_xyyzz[k] * tke_0;

            to_0_y_y_xzzz[k] = 2.0 * to_y_xyzzz[k] * tke_0;

            to_0_y_y_yyyy[k] = -4.0 * to_y_yyy[k] + 2.0 * to_y_yyyyy[k] * tke_0;

            to_0_y_y_yyyz[k] = -3.0 * to_y_yyz[k] + 2.0 * to_y_yyyyz[k] * tke_0;

            to_0_y_y_yyzz[k] = -2.0 * to_y_yzz[k] + 2.0 * to_y_yyyzz[k] * tke_0;

            to_0_y_y_yzzz[k] = -to_y_zzz[k] + 2.0 * to_y_yyzzz[k] * tke_0;

            to_0_y_y_zzzz[k] = 2.0 * to_y_yzzzz[k] * tke_0;
        }

        // Set up 75-90 components of targeted buffer : PG

        auto to_0_y_z_xxxx = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 30);

        auto to_0_y_z_xxxy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 31);

        auto to_0_y_z_xxxz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 32);

        auto to_0_y_z_xxyy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 33);

        auto to_0_y_z_xxyz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 34);

        auto to_0_y_z_xxzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 35);

        auto to_0_y_z_xyyy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 36);

        auto to_0_y_z_xyyz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 37);

        auto to_0_y_z_xyzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 38);

        auto to_0_y_z_xzzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 39);

        auto to_0_y_z_yyyy = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 40);

        auto to_0_y_z_yyyz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 41);

        auto to_0_y_z_yyzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 42);

        auto to_0_y_z_yzzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 43);

        auto to_0_y_z_zzzz = pbuffer.data(idx_op_geom_001_pg + 1 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_0_y_z_xxxx, to_0_y_z_xxxy, to_0_y_z_xxxz, to_0_y_z_xxyy, to_0_y_z_xxyz, to_0_y_z_xxzz, to_0_y_z_xyyy, to_0_y_z_xyyz, to_0_y_z_xyzz, to_0_y_z_xzzz, to_0_y_z_yyyy, to_0_y_z_yyyz, to_0_y_z_yyzz, to_0_y_z_yzzz, to_0_y_z_zzzz, to_z_xxx, to_z_xxxxy, to_z_xxxyy, to_z_xxxyz, to_z_xxy, to_z_xxyyy, to_z_xxyyz, to_z_xxyzz, to_z_xxz, to_z_xyy, to_z_xyyyy, to_z_xyyyz, to_z_xyyzz, to_z_xyz, to_z_xyzzz, to_z_xzz, to_z_yyy, to_z_yyyyy, to_z_yyyyz, to_z_yyyzz, to_z_yyz, to_z_yyzzz, to_z_yzz, to_z_yzzzz, to_z_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_z_xxxx[k] = 2.0 * to_z_xxxxy[k] * tke_0;

            to_0_y_z_xxxy[k] = -to_z_xxx[k] + 2.0 * to_z_xxxyy[k] * tke_0;

            to_0_y_z_xxxz[k] = 2.0 * to_z_xxxyz[k] * tke_0;

            to_0_y_z_xxyy[k] = -2.0 * to_z_xxy[k] + 2.0 * to_z_xxyyy[k] * tke_0;

            to_0_y_z_xxyz[k] = -to_z_xxz[k] + 2.0 * to_z_xxyyz[k] * tke_0;

            to_0_y_z_xxzz[k] = 2.0 * to_z_xxyzz[k] * tke_0;

            to_0_y_z_xyyy[k] = -3.0 * to_z_xyy[k] + 2.0 * to_z_xyyyy[k] * tke_0;

            to_0_y_z_xyyz[k] = -2.0 * to_z_xyz[k] + 2.0 * to_z_xyyyz[k] * tke_0;

            to_0_y_z_xyzz[k] = -to_z_xzz[k] + 2.0 * to_z_xyyzz[k] * tke_0;

            to_0_y_z_xzzz[k] = 2.0 * to_z_xyzzz[k] * tke_0;

            to_0_y_z_yyyy[k] = -4.0 * to_z_yyy[k] + 2.0 * to_z_yyyyy[k] * tke_0;

            to_0_y_z_yyyz[k] = -3.0 * to_z_yyz[k] + 2.0 * to_z_yyyyz[k] * tke_0;

            to_0_y_z_yyzz[k] = -2.0 * to_z_yzz[k] + 2.0 * to_z_yyyzz[k] * tke_0;

            to_0_y_z_yzzz[k] = -to_z_zzz[k] + 2.0 * to_z_yyzzz[k] * tke_0;

            to_0_y_z_zzzz[k] = 2.0 * to_z_yzzzz[k] * tke_0;
        }

        // Set up 90-105 components of targeted buffer : PG

        auto to_0_z_x_xxxx = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 0);

        auto to_0_z_x_xxxy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 1);

        auto to_0_z_x_xxxz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 2);

        auto to_0_z_x_xxyy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 3);

        auto to_0_z_x_xxyz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 4);

        auto to_0_z_x_xxzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 5);

        auto to_0_z_x_xyyy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 6);

        auto to_0_z_x_xyyz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 7);

        auto to_0_z_x_xyzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 8);

        auto to_0_z_x_xzzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 9);

        auto to_0_z_x_yyyy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 10);

        auto to_0_z_x_yyyz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 11);

        auto to_0_z_x_yyzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 12);

        auto to_0_z_x_yzzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 13);

        auto to_0_z_x_zzzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_0_z_x_xxxx, to_0_z_x_xxxy, to_0_z_x_xxxz, to_0_z_x_xxyy, to_0_z_x_xxyz, to_0_z_x_xxzz, to_0_z_x_xyyy, to_0_z_x_xyyz, to_0_z_x_xyzz, to_0_z_x_xzzz, to_0_z_x_yyyy, to_0_z_x_yyyz, to_0_z_x_yyzz, to_0_z_x_yzzz, to_0_z_x_zzzz, to_x_xxx, to_x_xxxxz, to_x_xxxyz, to_x_xxxzz, to_x_xxy, to_x_xxyyz, to_x_xxyzz, to_x_xxz, to_x_xxzzz, to_x_xyy, to_x_xyyyz, to_x_xyyzz, to_x_xyz, to_x_xyzzz, to_x_xzz, to_x_xzzzz, to_x_yyy, to_x_yyyyz, to_x_yyyzz, to_x_yyz, to_x_yyzzz, to_x_yzz, to_x_yzzzz, to_x_zzz, to_x_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_x_xxxx[k] = 2.0 * to_x_xxxxz[k] * tke_0;

            to_0_z_x_xxxy[k] = 2.0 * to_x_xxxyz[k] * tke_0;

            to_0_z_x_xxxz[k] = -to_x_xxx[k] + 2.0 * to_x_xxxzz[k] * tke_0;

            to_0_z_x_xxyy[k] = 2.0 * to_x_xxyyz[k] * tke_0;

            to_0_z_x_xxyz[k] = -to_x_xxy[k] + 2.0 * to_x_xxyzz[k] * tke_0;

            to_0_z_x_xxzz[k] = -2.0 * to_x_xxz[k] + 2.0 * to_x_xxzzz[k] * tke_0;

            to_0_z_x_xyyy[k] = 2.0 * to_x_xyyyz[k] * tke_0;

            to_0_z_x_xyyz[k] = -to_x_xyy[k] + 2.0 * to_x_xyyzz[k] * tke_0;

            to_0_z_x_xyzz[k] = -2.0 * to_x_xyz[k] + 2.0 * to_x_xyzzz[k] * tke_0;

            to_0_z_x_xzzz[k] = -3.0 * to_x_xzz[k] + 2.0 * to_x_xzzzz[k] * tke_0;

            to_0_z_x_yyyy[k] = 2.0 * to_x_yyyyz[k] * tke_0;

            to_0_z_x_yyyz[k] = -to_x_yyy[k] + 2.0 * to_x_yyyzz[k] * tke_0;

            to_0_z_x_yyzz[k] = -2.0 * to_x_yyz[k] + 2.0 * to_x_yyzzz[k] * tke_0;

            to_0_z_x_yzzz[k] = -3.0 * to_x_yzz[k] + 2.0 * to_x_yzzzz[k] * tke_0;

            to_0_z_x_zzzz[k] = -4.0 * to_x_zzz[k] + 2.0 * to_x_zzzzz[k] * tke_0;
        }

        // Set up 105-120 components of targeted buffer : PG

        auto to_0_z_y_xxxx = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 15);

        auto to_0_z_y_xxxy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 16);

        auto to_0_z_y_xxxz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 17);

        auto to_0_z_y_xxyy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 18);

        auto to_0_z_y_xxyz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 19);

        auto to_0_z_y_xxzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 20);

        auto to_0_z_y_xyyy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 21);

        auto to_0_z_y_xyyz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 22);

        auto to_0_z_y_xyzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 23);

        auto to_0_z_y_xzzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 24);

        auto to_0_z_y_yyyy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 25);

        auto to_0_z_y_yyyz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 26);

        auto to_0_z_y_yyzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 27);

        auto to_0_z_y_yzzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 28);

        auto to_0_z_y_zzzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_0_z_y_xxxx, to_0_z_y_xxxy, to_0_z_y_xxxz, to_0_z_y_xxyy, to_0_z_y_xxyz, to_0_z_y_xxzz, to_0_z_y_xyyy, to_0_z_y_xyyz, to_0_z_y_xyzz, to_0_z_y_xzzz, to_0_z_y_yyyy, to_0_z_y_yyyz, to_0_z_y_yyzz, to_0_z_y_yzzz, to_0_z_y_zzzz, to_y_xxx, to_y_xxxxz, to_y_xxxyz, to_y_xxxzz, to_y_xxy, to_y_xxyyz, to_y_xxyzz, to_y_xxz, to_y_xxzzz, to_y_xyy, to_y_xyyyz, to_y_xyyzz, to_y_xyz, to_y_xyzzz, to_y_xzz, to_y_xzzzz, to_y_yyy, to_y_yyyyz, to_y_yyyzz, to_y_yyz, to_y_yyzzz, to_y_yzz, to_y_yzzzz, to_y_zzz, to_y_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_y_xxxx[k] = 2.0 * to_y_xxxxz[k] * tke_0;

            to_0_z_y_xxxy[k] = 2.0 * to_y_xxxyz[k] * tke_0;

            to_0_z_y_xxxz[k] = -to_y_xxx[k] + 2.0 * to_y_xxxzz[k] * tke_0;

            to_0_z_y_xxyy[k] = 2.0 * to_y_xxyyz[k] * tke_0;

            to_0_z_y_xxyz[k] = -to_y_xxy[k] + 2.0 * to_y_xxyzz[k] * tke_0;

            to_0_z_y_xxzz[k] = -2.0 * to_y_xxz[k] + 2.0 * to_y_xxzzz[k] * tke_0;

            to_0_z_y_xyyy[k] = 2.0 * to_y_xyyyz[k] * tke_0;

            to_0_z_y_xyyz[k] = -to_y_xyy[k] + 2.0 * to_y_xyyzz[k] * tke_0;

            to_0_z_y_xyzz[k] = -2.0 * to_y_xyz[k] + 2.0 * to_y_xyzzz[k] * tke_0;

            to_0_z_y_xzzz[k] = -3.0 * to_y_xzz[k] + 2.0 * to_y_xzzzz[k] * tke_0;

            to_0_z_y_yyyy[k] = 2.0 * to_y_yyyyz[k] * tke_0;

            to_0_z_y_yyyz[k] = -to_y_yyy[k] + 2.0 * to_y_yyyzz[k] * tke_0;

            to_0_z_y_yyzz[k] = -2.0 * to_y_yyz[k] + 2.0 * to_y_yyzzz[k] * tke_0;

            to_0_z_y_yzzz[k] = -3.0 * to_y_yzz[k] + 2.0 * to_y_yzzzz[k] * tke_0;

            to_0_z_y_zzzz[k] = -4.0 * to_y_zzz[k] + 2.0 * to_y_zzzzz[k] * tke_0;
        }

        // Set up 120-135 components of targeted buffer : PG

        auto to_0_z_z_xxxx = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 30);

        auto to_0_z_z_xxxy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 31);

        auto to_0_z_z_xxxz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 32);

        auto to_0_z_z_xxyy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 33);

        auto to_0_z_z_xxyz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 34);

        auto to_0_z_z_xxzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 35);

        auto to_0_z_z_xyyy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 36);

        auto to_0_z_z_xyyz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 37);

        auto to_0_z_z_xyzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 38);

        auto to_0_z_z_xzzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 39);

        auto to_0_z_z_yyyy = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 40);

        auto to_0_z_z_yyyz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 41);

        auto to_0_z_z_yyzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 42);

        auto to_0_z_z_yzzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 43);

        auto to_0_z_z_zzzz = pbuffer.data(idx_op_geom_001_pg + 2 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_0_z_z_xxxx, to_0_z_z_xxxy, to_0_z_z_xxxz, to_0_z_z_xxyy, to_0_z_z_xxyz, to_0_z_z_xxzz, to_0_z_z_xyyy, to_0_z_z_xyyz, to_0_z_z_xyzz, to_0_z_z_xzzz, to_0_z_z_yyyy, to_0_z_z_yyyz, to_0_z_z_yyzz, to_0_z_z_yzzz, to_0_z_z_zzzz, to_z_xxx, to_z_xxxxz, to_z_xxxyz, to_z_xxxzz, to_z_xxy, to_z_xxyyz, to_z_xxyzz, to_z_xxz, to_z_xxzzz, to_z_xyy, to_z_xyyyz, to_z_xyyzz, to_z_xyz, to_z_xyzzz, to_z_xzz, to_z_xzzzz, to_z_yyy, to_z_yyyyz, to_z_yyyzz, to_z_yyz, to_z_yyzzz, to_z_yzz, to_z_yzzzz, to_z_zzz, to_z_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_z_xxxx[k] = 2.0 * to_z_xxxxz[k] * tke_0;

            to_0_z_z_xxxy[k] = 2.0 * to_z_xxxyz[k] * tke_0;

            to_0_z_z_xxxz[k] = -to_z_xxx[k] + 2.0 * to_z_xxxzz[k] * tke_0;

            to_0_z_z_xxyy[k] = 2.0 * to_z_xxyyz[k] * tke_0;

            to_0_z_z_xxyz[k] = -to_z_xxy[k] + 2.0 * to_z_xxyzz[k] * tke_0;

            to_0_z_z_xxzz[k] = -2.0 * to_z_xxz[k] + 2.0 * to_z_xxzzz[k] * tke_0;

            to_0_z_z_xyyy[k] = 2.0 * to_z_xyyyz[k] * tke_0;

            to_0_z_z_xyyz[k] = -to_z_xyy[k] + 2.0 * to_z_xyyzz[k] * tke_0;

            to_0_z_z_xyzz[k] = -2.0 * to_z_xyz[k] + 2.0 * to_z_xyzzz[k] * tke_0;

            to_0_z_z_xzzz[k] = -3.0 * to_z_xzz[k] + 2.0 * to_z_xzzzz[k] * tke_0;

            to_0_z_z_yyyy[k] = 2.0 * to_z_yyyyz[k] * tke_0;

            to_0_z_z_yyyz[k] = -to_z_yyy[k] + 2.0 * to_z_yyyzz[k] * tke_0;

            to_0_z_z_yyzz[k] = -2.0 * to_z_yyz[k] + 2.0 * to_z_yyzzz[k] * tke_0;

            to_0_z_z_yzzz[k] = -3.0 * to_z_yzz[k] + 2.0 * to_z_yzzzz[k] * tke_0;

            to_0_z_z_zzzz[k] = -4.0 * to_z_zzz[k] + 2.0 * to_z_zzzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

