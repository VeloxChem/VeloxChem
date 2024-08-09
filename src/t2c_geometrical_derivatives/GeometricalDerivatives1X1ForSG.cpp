#include "GeometricalDerivatives1X1ForSG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_sg(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_sg,
                        const size_t idx_op_pf,
                        const size_t idx_op_ph,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
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

        // Set up 0-15 components of targeted buffer : SG

        auto to_x_x_0_xxxx = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 0);

        auto to_x_x_0_xxxy = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 1);

        auto to_x_x_0_xxxz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 2);

        auto to_x_x_0_xxyy = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 3);

        auto to_x_x_0_xxyz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 4);

        auto to_x_x_0_xxzz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 5);

        auto to_x_x_0_xyyy = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 6);

        auto to_x_x_0_xyyz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 7);

        auto to_x_x_0_xyzz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 8);

        auto to_x_x_0_xzzz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 9);

        auto to_x_x_0_yyyy = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 10);

        auto to_x_x_0_yyyz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 11);

        auto to_x_x_0_yyzz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 12);

        auto to_x_x_0_yzzz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 13);

        auto to_x_x_0_zzzz = pbuffer.data(idx_op_geom_101_sg + 0 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_x_x_0_xxxx, to_x_x_0_xxxy, to_x_x_0_xxxz, to_x_x_0_xxyy, to_x_x_0_xxyz, to_x_x_0_xxzz, to_x_x_0_xyyy, to_x_x_0_xyyz, to_x_x_0_xyzz, to_x_x_0_xzzz, to_x_x_0_yyyy, to_x_x_0_yyyz, to_x_x_0_yyzz, to_x_x_0_yzzz, to_x_x_0_zzzz, to_x_xxx, to_x_xxxxx, to_x_xxxxy, to_x_xxxxz, to_x_xxxyy, to_x_xxxyz, to_x_xxxzz, to_x_xxy, to_x_xxyyy, to_x_xxyyz, to_x_xxyzz, to_x_xxz, to_x_xxzzz, to_x_xyy, to_x_xyyyy, to_x_xyyyz, to_x_xyyzz, to_x_xyz, to_x_xyzzz, to_x_xzz, to_x_xzzzz, to_x_yyy, to_x_yyz, to_x_yzz, to_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_0_xxxx[k] = -8.0 * to_x_xxx[k] * tbe_0 + 4.0 * to_x_xxxxx[k] * tbe_0 * tke_0;

            to_x_x_0_xxxy[k] = -6.0 * to_x_xxy[k] * tbe_0 + 4.0 * to_x_xxxxy[k] * tbe_0 * tke_0;

            to_x_x_0_xxxz[k] = -6.0 * to_x_xxz[k] * tbe_0 + 4.0 * to_x_xxxxz[k] * tbe_0 * tke_0;

            to_x_x_0_xxyy[k] = -4.0 * to_x_xyy[k] * tbe_0 + 4.0 * to_x_xxxyy[k] * tbe_0 * tke_0;

            to_x_x_0_xxyz[k] = -4.0 * to_x_xyz[k] * tbe_0 + 4.0 * to_x_xxxyz[k] * tbe_0 * tke_0;

            to_x_x_0_xxzz[k] = -4.0 * to_x_xzz[k] * tbe_0 + 4.0 * to_x_xxxzz[k] * tbe_0 * tke_0;

            to_x_x_0_xyyy[k] = -2.0 * to_x_yyy[k] * tbe_0 + 4.0 * to_x_xxyyy[k] * tbe_0 * tke_0;

            to_x_x_0_xyyz[k] = -2.0 * to_x_yyz[k] * tbe_0 + 4.0 * to_x_xxyyz[k] * tbe_0 * tke_0;

            to_x_x_0_xyzz[k] = -2.0 * to_x_yzz[k] * tbe_0 + 4.0 * to_x_xxyzz[k] * tbe_0 * tke_0;

            to_x_x_0_xzzz[k] = -2.0 * to_x_zzz[k] * tbe_0 + 4.0 * to_x_xxzzz[k] * tbe_0 * tke_0;

            to_x_x_0_yyyy[k] = 4.0 * to_x_xyyyy[k] * tbe_0 * tke_0;

            to_x_x_0_yyyz[k] = 4.0 * to_x_xyyyz[k] * tbe_0 * tke_0;

            to_x_x_0_yyzz[k] = 4.0 * to_x_xyyzz[k] * tbe_0 * tke_0;

            to_x_x_0_yzzz[k] = 4.0 * to_x_xyzzz[k] * tbe_0 * tke_0;

            to_x_x_0_zzzz[k] = 4.0 * to_x_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 15-30 components of targeted buffer : SG

        auto to_x_y_0_xxxx = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 0);

        auto to_x_y_0_xxxy = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 1);

        auto to_x_y_0_xxxz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 2);

        auto to_x_y_0_xxyy = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 3);

        auto to_x_y_0_xxyz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 4);

        auto to_x_y_0_xxzz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 5);

        auto to_x_y_0_xyyy = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 6);

        auto to_x_y_0_xyyz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 7);

        auto to_x_y_0_xyzz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 8);

        auto to_x_y_0_xzzz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 9);

        auto to_x_y_0_yyyy = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 10);

        auto to_x_y_0_yyyz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 11);

        auto to_x_y_0_yyzz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 12);

        auto to_x_y_0_yzzz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 13);

        auto to_x_y_0_zzzz = pbuffer.data(idx_op_geom_101_sg + 1 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_x_xxx, to_x_xxxxy, to_x_xxxyy, to_x_xxxyz, to_x_xxy, to_x_xxyyy, to_x_xxyyz, to_x_xxyzz, to_x_xxz, to_x_xyy, to_x_xyyyy, to_x_xyyyz, to_x_xyyzz, to_x_xyz, to_x_xyzzz, to_x_xzz, to_x_y_0_xxxx, to_x_y_0_xxxy, to_x_y_0_xxxz, to_x_y_0_xxyy, to_x_y_0_xxyz, to_x_y_0_xxzz, to_x_y_0_xyyy, to_x_y_0_xyyz, to_x_y_0_xyzz, to_x_y_0_xzzz, to_x_y_0_yyyy, to_x_y_0_yyyz, to_x_y_0_yyzz, to_x_y_0_yzzz, to_x_y_0_zzzz, to_x_yyy, to_x_yyyyy, to_x_yyyyz, to_x_yyyzz, to_x_yyz, to_x_yyzzz, to_x_yzz, to_x_yzzzz, to_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_0_xxxx[k] = 4.0 * to_x_xxxxy[k] * tbe_0 * tke_0;

            to_x_y_0_xxxy[k] = -2.0 * to_x_xxx[k] * tbe_0 + 4.0 * to_x_xxxyy[k] * tbe_0 * tke_0;

            to_x_y_0_xxxz[k] = 4.0 * to_x_xxxyz[k] * tbe_0 * tke_0;

            to_x_y_0_xxyy[k] = -4.0 * to_x_xxy[k] * tbe_0 + 4.0 * to_x_xxyyy[k] * tbe_0 * tke_0;

            to_x_y_0_xxyz[k] = -2.0 * to_x_xxz[k] * tbe_0 + 4.0 * to_x_xxyyz[k] * tbe_0 * tke_0;

            to_x_y_0_xxzz[k] = 4.0 * to_x_xxyzz[k] * tbe_0 * tke_0;

            to_x_y_0_xyyy[k] = -6.0 * to_x_xyy[k] * tbe_0 + 4.0 * to_x_xyyyy[k] * tbe_0 * tke_0;

            to_x_y_0_xyyz[k] = -4.0 * to_x_xyz[k] * tbe_0 + 4.0 * to_x_xyyyz[k] * tbe_0 * tke_0;

            to_x_y_0_xyzz[k] = -2.0 * to_x_xzz[k] * tbe_0 + 4.0 * to_x_xyyzz[k] * tbe_0 * tke_0;

            to_x_y_0_xzzz[k] = 4.0 * to_x_xyzzz[k] * tbe_0 * tke_0;

            to_x_y_0_yyyy[k] = -8.0 * to_x_yyy[k] * tbe_0 + 4.0 * to_x_yyyyy[k] * tbe_0 * tke_0;

            to_x_y_0_yyyz[k] = -6.0 * to_x_yyz[k] * tbe_0 + 4.0 * to_x_yyyyz[k] * tbe_0 * tke_0;

            to_x_y_0_yyzz[k] = -4.0 * to_x_yzz[k] * tbe_0 + 4.0 * to_x_yyyzz[k] * tbe_0 * tke_0;

            to_x_y_0_yzzz[k] = -2.0 * to_x_zzz[k] * tbe_0 + 4.0 * to_x_yyzzz[k] * tbe_0 * tke_0;

            to_x_y_0_zzzz[k] = 4.0 * to_x_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-45 components of targeted buffer : SG

        auto to_x_z_0_xxxx = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 0);

        auto to_x_z_0_xxxy = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 1);

        auto to_x_z_0_xxxz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 2);

        auto to_x_z_0_xxyy = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 3);

        auto to_x_z_0_xxyz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 4);

        auto to_x_z_0_xxzz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 5);

        auto to_x_z_0_xyyy = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 6);

        auto to_x_z_0_xyyz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 7);

        auto to_x_z_0_xyzz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 8);

        auto to_x_z_0_xzzz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 9);

        auto to_x_z_0_yyyy = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 10);

        auto to_x_z_0_yyyz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 11);

        auto to_x_z_0_yyzz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 12);

        auto to_x_z_0_yzzz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 13);

        auto to_x_z_0_zzzz = pbuffer.data(idx_op_geom_101_sg + 2 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_x_xxx, to_x_xxxxz, to_x_xxxyz, to_x_xxxzz, to_x_xxy, to_x_xxyyz, to_x_xxyzz, to_x_xxz, to_x_xxzzz, to_x_xyy, to_x_xyyyz, to_x_xyyzz, to_x_xyz, to_x_xyzzz, to_x_xzz, to_x_xzzzz, to_x_yyy, to_x_yyyyz, to_x_yyyzz, to_x_yyz, to_x_yyzzz, to_x_yzz, to_x_yzzzz, to_x_z_0_xxxx, to_x_z_0_xxxy, to_x_z_0_xxxz, to_x_z_0_xxyy, to_x_z_0_xxyz, to_x_z_0_xxzz, to_x_z_0_xyyy, to_x_z_0_xyyz, to_x_z_0_xyzz, to_x_z_0_xzzz, to_x_z_0_yyyy, to_x_z_0_yyyz, to_x_z_0_yyzz, to_x_z_0_yzzz, to_x_z_0_zzzz, to_x_zzz, to_x_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_0_xxxx[k] = 4.0 * to_x_xxxxz[k] * tbe_0 * tke_0;

            to_x_z_0_xxxy[k] = 4.0 * to_x_xxxyz[k] * tbe_0 * tke_0;

            to_x_z_0_xxxz[k] = -2.0 * to_x_xxx[k] * tbe_0 + 4.0 * to_x_xxxzz[k] * tbe_0 * tke_0;

            to_x_z_0_xxyy[k] = 4.0 * to_x_xxyyz[k] * tbe_0 * tke_0;

            to_x_z_0_xxyz[k] = -2.0 * to_x_xxy[k] * tbe_0 + 4.0 * to_x_xxyzz[k] * tbe_0 * tke_0;

            to_x_z_0_xxzz[k] = -4.0 * to_x_xxz[k] * tbe_0 + 4.0 * to_x_xxzzz[k] * tbe_0 * tke_0;

            to_x_z_0_xyyy[k] = 4.0 * to_x_xyyyz[k] * tbe_0 * tke_0;

            to_x_z_0_xyyz[k] = -2.0 * to_x_xyy[k] * tbe_0 + 4.0 * to_x_xyyzz[k] * tbe_0 * tke_0;

            to_x_z_0_xyzz[k] = -4.0 * to_x_xyz[k] * tbe_0 + 4.0 * to_x_xyzzz[k] * tbe_0 * tke_0;

            to_x_z_0_xzzz[k] = -6.0 * to_x_xzz[k] * tbe_0 + 4.0 * to_x_xzzzz[k] * tbe_0 * tke_0;

            to_x_z_0_yyyy[k] = 4.0 * to_x_yyyyz[k] * tbe_0 * tke_0;

            to_x_z_0_yyyz[k] = -2.0 * to_x_yyy[k] * tbe_0 + 4.0 * to_x_yyyzz[k] * tbe_0 * tke_0;

            to_x_z_0_yyzz[k] = -4.0 * to_x_yyz[k] * tbe_0 + 4.0 * to_x_yyzzz[k] * tbe_0 * tke_0;

            to_x_z_0_yzzz[k] = -6.0 * to_x_yzz[k] * tbe_0 + 4.0 * to_x_yzzzz[k] * tbe_0 * tke_0;

            to_x_z_0_zzzz[k] = -8.0 * to_x_zzz[k] * tbe_0 + 4.0 * to_x_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 45-60 components of targeted buffer : SG

        auto to_y_x_0_xxxx = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 0);

        auto to_y_x_0_xxxy = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 1);

        auto to_y_x_0_xxxz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 2);

        auto to_y_x_0_xxyy = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 3);

        auto to_y_x_0_xxyz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 4);

        auto to_y_x_0_xxzz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 5);

        auto to_y_x_0_xyyy = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 6);

        auto to_y_x_0_xyyz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 7);

        auto to_y_x_0_xyzz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 8);

        auto to_y_x_0_xzzz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 9);

        auto to_y_x_0_yyyy = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 10);

        auto to_y_x_0_yyyz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 11);

        auto to_y_x_0_yyzz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 12);

        auto to_y_x_0_yzzz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 13);

        auto to_y_x_0_zzzz = pbuffer.data(idx_op_geom_101_sg + 3 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_y_x_0_xxxx, to_y_x_0_xxxy, to_y_x_0_xxxz, to_y_x_0_xxyy, to_y_x_0_xxyz, to_y_x_0_xxzz, to_y_x_0_xyyy, to_y_x_0_xyyz, to_y_x_0_xyzz, to_y_x_0_xzzz, to_y_x_0_yyyy, to_y_x_0_yyyz, to_y_x_0_yyzz, to_y_x_0_yzzz, to_y_x_0_zzzz, to_y_xxx, to_y_xxxxx, to_y_xxxxy, to_y_xxxxz, to_y_xxxyy, to_y_xxxyz, to_y_xxxzz, to_y_xxy, to_y_xxyyy, to_y_xxyyz, to_y_xxyzz, to_y_xxz, to_y_xxzzz, to_y_xyy, to_y_xyyyy, to_y_xyyyz, to_y_xyyzz, to_y_xyz, to_y_xyzzz, to_y_xzz, to_y_xzzzz, to_y_yyy, to_y_yyz, to_y_yzz, to_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_0_xxxx[k] = -8.0 * to_y_xxx[k] * tbe_0 + 4.0 * to_y_xxxxx[k] * tbe_0 * tke_0;

            to_y_x_0_xxxy[k] = -6.0 * to_y_xxy[k] * tbe_0 + 4.0 * to_y_xxxxy[k] * tbe_0 * tke_0;

            to_y_x_0_xxxz[k] = -6.0 * to_y_xxz[k] * tbe_0 + 4.0 * to_y_xxxxz[k] * tbe_0 * tke_0;

            to_y_x_0_xxyy[k] = -4.0 * to_y_xyy[k] * tbe_0 + 4.0 * to_y_xxxyy[k] * tbe_0 * tke_0;

            to_y_x_0_xxyz[k] = -4.0 * to_y_xyz[k] * tbe_0 + 4.0 * to_y_xxxyz[k] * tbe_0 * tke_0;

            to_y_x_0_xxzz[k] = -4.0 * to_y_xzz[k] * tbe_0 + 4.0 * to_y_xxxzz[k] * tbe_0 * tke_0;

            to_y_x_0_xyyy[k] = -2.0 * to_y_yyy[k] * tbe_0 + 4.0 * to_y_xxyyy[k] * tbe_0 * tke_0;

            to_y_x_0_xyyz[k] = -2.0 * to_y_yyz[k] * tbe_0 + 4.0 * to_y_xxyyz[k] * tbe_0 * tke_0;

            to_y_x_0_xyzz[k] = -2.0 * to_y_yzz[k] * tbe_0 + 4.0 * to_y_xxyzz[k] * tbe_0 * tke_0;

            to_y_x_0_xzzz[k] = -2.0 * to_y_zzz[k] * tbe_0 + 4.0 * to_y_xxzzz[k] * tbe_0 * tke_0;

            to_y_x_0_yyyy[k] = 4.0 * to_y_xyyyy[k] * tbe_0 * tke_0;

            to_y_x_0_yyyz[k] = 4.0 * to_y_xyyyz[k] * tbe_0 * tke_0;

            to_y_x_0_yyzz[k] = 4.0 * to_y_xyyzz[k] * tbe_0 * tke_0;

            to_y_x_0_yzzz[k] = 4.0 * to_y_xyzzz[k] * tbe_0 * tke_0;

            to_y_x_0_zzzz[k] = 4.0 * to_y_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-75 components of targeted buffer : SG

        auto to_y_y_0_xxxx = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 0);

        auto to_y_y_0_xxxy = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 1);

        auto to_y_y_0_xxxz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 2);

        auto to_y_y_0_xxyy = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 3);

        auto to_y_y_0_xxyz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 4);

        auto to_y_y_0_xxzz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 5);

        auto to_y_y_0_xyyy = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 6);

        auto to_y_y_0_xyyz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 7);

        auto to_y_y_0_xyzz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 8);

        auto to_y_y_0_xzzz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 9);

        auto to_y_y_0_yyyy = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 10);

        auto to_y_y_0_yyyz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 11);

        auto to_y_y_0_yyzz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 12);

        auto to_y_y_0_yzzz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 13);

        auto to_y_y_0_zzzz = pbuffer.data(idx_op_geom_101_sg + 4 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_y_xxx, to_y_xxxxy, to_y_xxxyy, to_y_xxxyz, to_y_xxy, to_y_xxyyy, to_y_xxyyz, to_y_xxyzz, to_y_xxz, to_y_xyy, to_y_xyyyy, to_y_xyyyz, to_y_xyyzz, to_y_xyz, to_y_xyzzz, to_y_xzz, to_y_y_0_xxxx, to_y_y_0_xxxy, to_y_y_0_xxxz, to_y_y_0_xxyy, to_y_y_0_xxyz, to_y_y_0_xxzz, to_y_y_0_xyyy, to_y_y_0_xyyz, to_y_y_0_xyzz, to_y_y_0_xzzz, to_y_y_0_yyyy, to_y_y_0_yyyz, to_y_y_0_yyzz, to_y_y_0_yzzz, to_y_y_0_zzzz, to_y_yyy, to_y_yyyyy, to_y_yyyyz, to_y_yyyzz, to_y_yyz, to_y_yyzzz, to_y_yzz, to_y_yzzzz, to_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_0_xxxx[k] = 4.0 * to_y_xxxxy[k] * tbe_0 * tke_0;

            to_y_y_0_xxxy[k] = -2.0 * to_y_xxx[k] * tbe_0 + 4.0 * to_y_xxxyy[k] * tbe_0 * tke_0;

            to_y_y_0_xxxz[k] = 4.0 * to_y_xxxyz[k] * tbe_0 * tke_0;

            to_y_y_0_xxyy[k] = -4.0 * to_y_xxy[k] * tbe_0 + 4.0 * to_y_xxyyy[k] * tbe_0 * tke_0;

            to_y_y_0_xxyz[k] = -2.0 * to_y_xxz[k] * tbe_0 + 4.0 * to_y_xxyyz[k] * tbe_0 * tke_0;

            to_y_y_0_xxzz[k] = 4.0 * to_y_xxyzz[k] * tbe_0 * tke_0;

            to_y_y_0_xyyy[k] = -6.0 * to_y_xyy[k] * tbe_0 + 4.0 * to_y_xyyyy[k] * tbe_0 * tke_0;

            to_y_y_0_xyyz[k] = -4.0 * to_y_xyz[k] * tbe_0 + 4.0 * to_y_xyyyz[k] * tbe_0 * tke_0;

            to_y_y_0_xyzz[k] = -2.0 * to_y_xzz[k] * tbe_0 + 4.0 * to_y_xyyzz[k] * tbe_0 * tke_0;

            to_y_y_0_xzzz[k] = 4.0 * to_y_xyzzz[k] * tbe_0 * tke_0;

            to_y_y_0_yyyy[k] = -8.0 * to_y_yyy[k] * tbe_0 + 4.0 * to_y_yyyyy[k] * tbe_0 * tke_0;

            to_y_y_0_yyyz[k] = -6.0 * to_y_yyz[k] * tbe_0 + 4.0 * to_y_yyyyz[k] * tbe_0 * tke_0;

            to_y_y_0_yyzz[k] = -4.0 * to_y_yzz[k] * tbe_0 + 4.0 * to_y_yyyzz[k] * tbe_0 * tke_0;

            to_y_y_0_yzzz[k] = -2.0 * to_y_zzz[k] * tbe_0 + 4.0 * to_y_yyzzz[k] * tbe_0 * tke_0;

            to_y_y_0_zzzz[k] = 4.0 * to_y_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 75-90 components of targeted buffer : SG

        auto to_y_z_0_xxxx = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 0);

        auto to_y_z_0_xxxy = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 1);

        auto to_y_z_0_xxxz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 2);

        auto to_y_z_0_xxyy = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 3);

        auto to_y_z_0_xxyz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 4);

        auto to_y_z_0_xxzz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 5);

        auto to_y_z_0_xyyy = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 6);

        auto to_y_z_0_xyyz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 7);

        auto to_y_z_0_xyzz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 8);

        auto to_y_z_0_xzzz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 9);

        auto to_y_z_0_yyyy = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 10);

        auto to_y_z_0_yyyz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 11);

        auto to_y_z_0_yyzz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 12);

        auto to_y_z_0_yzzz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 13);

        auto to_y_z_0_zzzz = pbuffer.data(idx_op_geom_101_sg + 5 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_y_xxx, to_y_xxxxz, to_y_xxxyz, to_y_xxxzz, to_y_xxy, to_y_xxyyz, to_y_xxyzz, to_y_xxz, to_y_xxzzz, to_y_xyy, to_y_xyyyz, to_y_xyyzz, to_y_xyz, to_y_xyzzz, to_y_xzz, to_y_xzzzz, to_y_yyy, to_y_yyyyz, to_y_yyyzz, to_y_yyz, to_y_yyzzz, to_y_yzz, to_y_yzzzz, to_y_z_0_xxxx, to_y_z_0_xxxy, to_y_z_0_xxxz, to_y_z_0_xxyy, to_y_z_0_xxyz, to_y_z_0_xxzz, to_y_z_0_xyyy, to_y_z_0_xyyz, to_y_z_0_xyzz, to_y_z_0_xzzz, to_y_z_0_yyyy, to_y_z_0_yyyz, to_y_z_0_yyzz, to_y_z_0_yzzz, to_y_z_0_zzzz, to_y_zzz, to_y_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_0_xxxx[k] = 4.0 * to_y_xxxxz[k] * tbe_0 * tke_0;

            to_y_z_0_xxxy[k] = 4.0 * to_y_xxxyz[k] * tbe_0 * tke_0;

            to_y_z_0_xxxz[k] = -2.0 * to_y_xxx[k] * tbe_0 + 4.0 * to_y_xxxzz[k] * tbe_0 * tke_0;

            to_y_z_0_xxyy[k] = 4.0 * to_y_xxyyz[k] * tbe_0 * tke_0;

            to_y_z_0_xxyz[k] = -2.0 * to_y_xxy[k] * tbe_0 + 4.0 * to_y_xxyzz[k] * tbe_0 * tke_0;

            to_y_z_0_xxzz[k] = -4.0 * to_y_xxz[k] * tbe_0 + 4.0 * to_y_xxzzz[k] * tbe_0 * tke_0;

            to_y_z_0_xyyy[k] = 4.0 * to_y_xyyyz[k] * tbe_0 * tke_0;

            to_y_z_0_xyyz[k] = -2.0 * to_y_xyy[k] * tbe_0 + 4.0 * to_y_xyyzz[k] * tbe_0 * tke_0;

            to_y_z_0_xyzz[k] = -4.0 * to_y_xyz[k] * tbe_0 + 4.0 * to_y_xyzzz[k] * tbe_0 * tke_0;

            to_y_z_0_xzzz[k] = -6.0 * to_y_xzz[k] * tbe_0 + 4.0 * to_y_xzzzz[k] * tbe_0 * tke_0;

            to_y_z_0_yyyy[k] = 4.0 * to_y_yyyyz[k] * tbe_0 * tke_0;

            to_y_z_0_yyyz[k] = -2.0 * to_y_yyy[k] * tbe_0 + 4.0 * to_y_yyyzz[k] * tbe_0 * tke_0;

            to_y_z_0_yyzz[k] = -4.0 * to_y_yyz[k] * tbe_0 + 4.0 * to_y_yyzzz[k] * tbe_0 * tke_0;

            to_y_z_0_yzzz[k] = -6.0 * to_y_yzz[k] * tbe_0 + 4.0 * to_y_yzzzz[k] * tbe_0 * tke_0;

            to_y_z_0_zzzz[k] = -8.0 * to_y_zzz[k] * tbe_0 + 4.0 * to_y_zzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-105 components of targeted buffer : SG

        auto to_z_x_0_xxxx = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 0);

        auto to_z_x_0_xxxy = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 1);

        auto to_z_x_0_xxxz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 2);

        auto to_z_x_0_xxyy = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 3);

        auto to_z_x_0_xxyz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 4);

        auto to_z_x_0_xxzz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 5);

        auto to_z_x_0_xyyy = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 6);

        auto to_z_x_0_xyyz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 7);

        auto to_z_x_0_xyzz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 8);

        auto to_z_x_0_xzzz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 9);

        auto to_z_x_0_yyyy = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 10);

        auto to_z_x_0_yyyz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 11);

        auto to_z_x_0_yyzz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 12);

        auto to_z_x_0_yzzz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 13);

        auto to_z_x_0_zzzz = pbuffer.data(idx_op_geom_101_sg + 6 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_z_x_0_xxxx, to_z_x_0_xxxy, to_z_x_0_xxxz, to_z_x_0_xxyy, to_z_x_0_xxyz, to_z_x_0_xxzz, to_z_x_0_xyyy, to_z_x_0_xyyz, to_z_x_0_xyzz, to_z_x_0_xzzz, to_z_x_0_yyyy, to_z_x_0_yyyz, to_z_x_0_yyzz, to_z_x_0_yzzz, to_z_x_0_zzzz, to_z_xxx, to_z_xxxxx, to_z_xxxxy, to_z_xxxxz, to_z_xxxyy, to_z_xxxyz, to_z_xxxzz, to_z_xxy, to_z_xxyyy, to_z_xxyyz, to_z_xxyzz, to_z_xxz, to_z_xxzzz, to_z_xyy, to_z_xyyyy, to_z_xyyyz, to_z_xyyzz, to_z_xyz, to_z_xyzzz, to_z_xzz, to_z_xzzzz, to_z_yyy, to_z_yyz, to_z_yzz, to_z_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_0_xxxx[k] = -8.0 * to_z_xxx[k] * tbe_0 + 4.0 * to_z_xxxxx[k] * tbe_0 * tke_0;

            to_z_x_0_xxxy[k] = -6.0 * to_z_xxy[k] * tbe_0 + 4.0 * to_z_xxxxy[k] * tbe_0 * tke_0;

            to_z_x_0_xxxz[k] = -6.0 * to_z_xxz[k] * tbe_0 + 4.0 * to_z_xxxxz[k] * tbe_0 * tke_0;

            to_z_x_0_xxyy[k] = -4.0 * to_z_xyy[k] * tbe_0 + 4.0 * to_z_xxxyy[k] * tbe_0 * tke_0;

            to_z_x_0_xxyz[k] = -4.0 * to_z_xyz[k] * tbe_0 + 4.0 * to_z_xxxyz[k] * tbe_0 * tke_0;

            to_z_x_0_xxzz[k] = -4.0 * to_z_xzz[k] * tbe_0 + 4.0 * to_z_xxxzz[k] * tbe_0 * tke_0;

            to_z_x_0_xyyy[k] = -2.0 * to_z_yyy[k] * tbe_0 + 4.0 * to_z_xxyyy[k] * tbe_0 * tke_0;

            to_z_x_0_xyyz[k] = -2.0 * to_z_yyz[k] * tbe_0 + 4.0 * to_z_xxyyz[k] * tbe_0 * tke_0;

            to_z_x_0_xyzz[k] = -2.0 * to_z_yzz[k] * tbe_0 + 4.0 * to_z_xxyzz[k] * tbe_0 * tke_0;

            to_z_x_0_xzzz[k] = -2.0 * to_z_zzz[k] * tbe_0 + 4.0 * to_z_xxzzz[k] * tbe_0 * tke_0;

            to_z_x_0_yyyy[k] = 4.0 * to_z_xyyyy[k] * tbe_0 * tke_0;

            to_z_x_0_yyyz[k] = 4.0 * to_z_xyyyz[k] * tbe_0 * tke_0;

            to_z_x_0_yyzz[k] = 4.0 * to_z_xyyzz[k] * tbe_0 * tke_0;

            to_z_x_0_yzzz[k] = 4.0 * to_z_xyzzz[k] * tbe_0 * tke_0;

            to_z_x_0_zzzz[k] = 4.0 * to_z_xzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 105-120 components of targeted buffer : SG

        auto to_z_y_0_xxxx = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 0);

        auto to_z_y_0_xxxy = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 1);

        auto to_z_y_0_xxxz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 2);

        auto to_z_y_0_xxyy = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 3);

        auto to_z_y_0_xxyz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 4);

        auto to_z_y_0_xxzz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 5);

        auto to_z_y_0_xyyy = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 6);

        auto to_z_y_0_xyyz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 7);

        auto to_z_y_0_xyzz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 8);

        auto to_z_y_0_xzzz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 9);

        auto to_z_y_0_yyyy = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 10);

        auto to_z_y_0_yyyz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 11);

        auto to_z_y_0_yyzz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 12);

        auto to_z_y_0_yzzz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 13);

        auto to_z_y_0_zzzz = pbuffer.data(idx_op_geom_101_sg + 7 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_z_xxx, to_z_xxxxy, to_z_xxxyy, to_z_xxxyz, to_z_xxy, to_z_xxyyy, to_z_xxyyz, to_z_xxyzz, to_z_xxz, to_z_xyy, to_z_xyyyy, to_z_xyyyz, to_z_xyyzz, to_z_xyz, to_z_xyzzz, to_z_xzz, to_z_y_0_xxxx, to_z_y_0_xxxy, to_z_y_0_xxxz, to_z_y_0_xxyy, to_z_y_0_xxyz, to_z_y_0_xxzz, to_z_y_0_xyyy, to_z_y_0_xyyz, to_z_y_0_xyzz, to_z_y_0_xzzz, to_z_y_0_yyyy, to_z_y_0_yyyz, to_z_y_0_yyzz, to_z_y_0_yzzz, to_z_y_0_zzzz, to_z_yyy, to_z_yyyyy, to_z_yyyyz, to_z_yyyzz, to_z_yyz, to_z_yyzzz, to_z_yzz, to_z_yzzzz, to_z_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_0_xxxx[k] = 4.0 * to_z_xxxxy[k] * tbe_0 * tke_0;

            to_z_y_0_xxxy[k] = -2.0 * to_z_xxx[k] * tbe_0 + 4.0 * to_z_xxxyy[k] * tbe_0 * tke_0;

            to_z_y_0_xxxz[k] = 4.0 * to_z_xxxyz[k] * tbe_0 * tke_0;

            to_z_y_0_xxyy[k] = -4.0 * to_z_xxy[k] * tbe_0 + 4.0 * to_z_xxyyy[k] * tbe_0 * tke_0;

            to_z_y_0_xxyz[k] = -2.0 * to_z_xxz[k] * tbe_0 + 4.0 * to_z_xxyyz[k] * tbe_0 * tke_0;

            to_z_y_0_xxzz[k] = 4.0 * to_z_xxyzz[k] * tbe_0 * tke_0;

            to_z_y_0_xyyy[k] = -6.0 * to_z_xyy[k] * tbe_0 + 4.0 * to_z_xyyyy[k] * tbe_0 * tke_0;

            to_z_y_0_xyyz[k] = -4.0 * to_z_xyz[k] * tbe_0 + 4.0 * to_z_xyyyz[k] * tbe_0 * tke_0;

            to_z_y_0_xyzz[k] = -2.0 * to_z_xzz[k] * tbe_0 + 4.0 * to_z_xyyzz[k] * tbe_0 * tke_0;

            to_z_y_0_xzzz[k] = 4.0 * to_z_xyzzz[k] * tbe_0 * tke_0;

            to_z_y_0_yyyy[k] = -8.0 * to_z_yyy[k] * tbe_0 + 4.0 * to_z_yyyyy[k] * tbe_0 * tke_0;

            to_z_y_0_yyyz[k] = -6.0 * to_z_yyz[k] * tbe_0 + 4.0 * to_z_yyyyz[k] * tbe_0 * tke_0;

            to_z_y_0_yyzz[k] = -4.0 * to_z_yzz[k] * tbe_0 + 4.0 * to_z_yyyzz[k] * tbe_0 * tke_0;

            to_z_y_0_yzzz[k] = -2.0 * to_z_zzz[k] * tbe_0 + 4.0 * to_z_yyzzz[k] * tbe_0 * tke_0;

            to_z_y_0_zzzz[k] = 4.0 * to_z_yzzzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-135 components of targeted buffer : SG

        auto to_z_z_0_xxxx = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 0);

        auto to_z_z_0_xxxy = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 1);

        auto to_z_z_0_xxxz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 2);

        auto to_z_z_0_xxyy = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 3);

        auto to_z_z_0_xxyz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 4);

        auto to_z_z_0_xxzz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 5);

        auto to_z_z_0_xyyy = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 6);

        auto to_z_z_0_xyyz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 7);

        auto to_z_z_0_xyzz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 8);

        auto to_z_z_0_xzzz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 9);

        auto to_z_z_0_yyyy = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 10);

        auto to_z_z_0_yyyz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 11);

        auto to_z_z_0_yyzz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 12);

        auto to_z_z_0_yzzz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 13);

        auto to_z_z_0_zzzz = pbuffer.data(idx_op_geom_101_sg + 8 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_z_xxx, to_z_xxxxz, to_z_xxxyz, to_z_xxxzz, to_z_xxy, to_z_xxyyz, to_z_xxyzz, to_z_xxz, to_z_xxzzz, to_z_xyy, to_z_xyyyz, to_z_xyyzz, to_z_xyz, to_z_xyzzz, to_z_xzz, to_z_xzzzz, to_z_yyy, to_z_yyyyz, to_z_yyyzz, to_z_yyz, to_z_yyzzz, to_z_yzz, to_z_yzzzz, to_z_z_0_xxxx, to_z_z_0_xxxy, to_z_z_0_xxxz, to_z_z_0_xxyy, to_z_z_0_xxyz, to_z_z_0_xxzz, to_z_z_0_xyyy, to_z_z_0_xyyz, to_z_z_0_xyzz, to_z_z_0_xzzz, to_z_z_0_yyyy, to_z_z_0_yyyz, to_z_z_0_yyzz, to_z_z_0_yzzz, to_z_z_0_zzzz, to_z_zzz, to_z_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_0_xxxx[k] = 4.0 * to_z_xxxxz[k] * tbe_0 * tke_0;

            to_z_z_0_xxxy[k] = 4.0 * to_z_xxxyz[k] * tbe_0 * tke_0;

            to_z_z_0_xxxz[k] = -2.0 * to_z_xxx[k] * tbe_0 + 4.0 * to_z_xxxzz[k] * tbe_0 * tke_0;

            to_z_z_0_xxyy[k] = 4.0 * to_z_xxyyz[k] * tbe_0 * tke_0;

            to_z_z_0_xxyz[k] = -2.0 * to_z_xxy[k] * tbe_0 + 4.0 * to_z_xxyzz[k] * tbe_0 * tke_0;

            to_z_z_0_xxzz[k] = -4.0 * to_z_xxz[k] * tbe_0 + 4.0 * to_z_xxzzz[k] * tbe_0 * tke_0;

            to_z_z_0_xyyy[k] = 4.0 * to_z_xyyyz[k] * tbe_0 * tke_0;

            to_z_z_0_xyyz[k] = -2.0 * to_z_xyy[k] * tbe_0 + 4.0 * to_z_xyyzz[k] * tbe_0 * tke_0;

            to_z_z_0_xyzz[k] = -4.0 * to_z_xyz[k] * tbe_0 + 4.0 * to_z_xyzzz[k] * tbe_0 * tke_0;

            to_z_z_0_xzzz[k] = -6.0 * to_z_xzz[k] * tbe_0 + 4.0 * to_z_xzzzz[k] * tbe_0 * tke_0;

            to_z_z_0_yyyy[k] = 4.0 * to_z_yyyyz[k] * tbe_0 * tke_0;

            to_z_z_0_yyyz[k] = -2.0 * to_z_yyy[k] * tbe_0 + 4.0 * to_z_yyyzz[k] * tbe_0 * tke_0;

            to_z_z_0_yyzz[k] = -4.0 * to_z_yyz[k] * tbe_0 + 4.0 * to_z_yyzzz[k] * tbe_0 * tke_0;

            to_z_z_0_yzzz[k] = -6.0 * to_z_yzz[k] * tbe_0 + 4.0 * to_z_yzzzz[k] * tbe_0 * tke_0;

            to_z_z_0_zzzz[k] = -8.0 * to_z_zzz[k] * tbe_0 + 4.0 * to_z_zzzzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

