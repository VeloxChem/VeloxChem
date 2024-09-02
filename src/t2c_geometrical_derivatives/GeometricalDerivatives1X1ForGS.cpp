#include "GeometricalDerivatives1X1ForGS.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_11_gs(CSimdArray<double>&       pbuffer,
                        const size_t              idx_op_geom_101_gs,
                        const size_t              idx_op_fp,
                        const size_t              idx_op_hp,
                        const size_t              op_comps,
                        const CSimdArray<double>& factors,
                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : FP

        auto to_xxx_x = pbuffer.data(idx_op_fp + i * 30 + 0);

        auto to_xxx_y = pbuffer.data(idx_op_fp + i * 30 + 1);

        auto to_xxx_z = pbuffer.data(idx_op_fp + i * 30 + 2);

        auto to_xxy_x = pbuffer.data(idx_op_fp + i * 30 + 3);

        auto to_xxy_y = pbuffer.data(idx_op_fp + i * 30 + 4);

        auto to_xxy_z = pbuffer.data(idx_op_fp + i * 30 + 5);

        auto to_xxz_x = pbuffer.data(idx_op_fp + i * 30 + 6);

        auto to_xxz_y = pbuffer.data(idx_op_fp + i * 30 + 7);

        auto to_xxz_z = pbuffer.data(idx_op_fp + i * 30 + 8);

        auto to_xyy_x = pbuffer.data(idx_op_fp + i * 30 + 9);

        auto to_xyy_y = pbuffer.data(idx_op_fp + i * 30 + 10);

        auto to_xyy_z = pbuffer.data(idx_op_fp + i * 30 + 11);

        auto to_xyz_x = pbuffer.data(idx_op_fp + i * 30 + 12);

        auto to_xyz_y = pbuffer.data(idx_op_fp + i * 30 + 13);

        auto to_xyz_z = pbuffer.data(idx_op_fp + i * 30 + 14);

        auto to_xzz_x = pbuffer.data(idx_op_fp + i * 30 + 15);

        auto to_xzz_y = pbuffer.data(idx_op_fp + i * 30 + 16);

        auto to_xzz_z = pbuffer.data(idx_op_fp + i * 30 + 17);

        auto to_yyy_x = pbuffer.data(idx_op_fp + i * 30 + 18);

        auto to_yyy_y = pbuffer.data(idx_op_fp + i * 30 + 19);

        auto to_yyy_z = pbuffer.data(idx_op_fp + i * 30 + 20);

        auto to_yyz_x = pbuffer.data(idx_op_fp + i * 30 + 21);

        auto to_yyz_y = pbuffer.data(idx_op_fp + i * 30 + 22);

        auto to_yyz_z = pbuffer.data(idx_op_fp + i * 30 + 23);

        auto to_yzz_x = pbuffer.data(idx_op_fp + i * 30 + 24);

        auto to_yzz_y = pbuffer.data(idx_op_fp + i * 30 + 25);

        auto to_yzz_z = pbuffer.data(idx_op_fp + i * 30 + 26);

        auto to_zzz_x = pbuffer.data(idx_op_fp + i * 30 + 27);

        auto to_zzz_y = pbuffer.data(idx_op_fp + i * 30 + 28);

        auto to_zzz_z = pbuffer.data(idx_op_fp + i * 30 + 29);

        // Set up components of auxiliary buffer : HP

        auto to_xxxxx_x = pbuffer.data(idx_op_hp + i * 63 + 0);

        auto to_xxxxx_y = pbuffer.data(idx_op_hp + i * 63 + 1);

        auto to_xxxxx_z = pbuffer.data(idx_op_hp + i * 63 + 2);

        auto to_xxxxy_x = pbuffer.data(idx_op_hp + i * 63 + 3);

        auto to_xxxxy_y = pbuffer.data(idx_op_hp + i * 63 + 4);

        auto to_xxxxy_z = pbuffer.data(idx_op_hp + i * 63 + 5);

        auto to_xxxxz_x = pbuffer.data(idx_op_hp + i * 63 + 6);

        auto to_xxxxz_y = pbuffer.data(idx_op_hp + i * 63 + 7);

        auto to_xxxxz_z = pbuffer.data(idx_op_hp + i * 63 + 8);

        auto to_xxxyy_x = pbuffer.data(idx_op_hp + i * 63 + 9);

        auto to_xxxyy_y = pbuffer.data(idx_op_hp + i * 63 + 10);

        auto to_xxxyy_z = pbuffer.data(idx_op_hp + i * 63 + 11);

        auto to_xxxyz_x = pbuffer.data(idx_op_hp + i * 63 + 12);

        auto to_xxxyz_y = pbuffer.data(idx_op_hp + i * 63 + 13);

        auto to_xxxyz_z = pbuffer.data(idx_op_hp + i * 63 + 14);

        auto to_xxxzz_x = pbuffer.data(idx_op_hp + i * 63 + 15);

        auto to_xxxzz_y = pbuffer.data(idx_op_hp + i * 63 + 16);

        auto to_xxxzz_z = pbuffer.data(idx_op_hp + i * 63 + 17);

        auto to_xxyyy_x = pbuffer.data(idx_op_hp + i * 63 + 18);

        auto to_xxyyy_y = pbuffer.data(idx_op_hp + i * 63 + 19);

        auto to_xxyyy_z = pbuffer.data(idx_op_hp + i * 63 + 20);

        auto to_xxyyz_x = pbuffer.data(idx_op_hp + i * 63 + 21);

        auto to_xxyyz_y = pbuffer.data(idx_op_hp + i * 63 + 22);

        auto to_xxyyz_z = pbuffer.data(idx_op_hp + i * 63 + 23);

        auto to_xxyzz_x = pbuffer.data(idx_op_hp + i * 63 + 24);

        auto to_xxyzz_y = pbuffer.data(idx_op_hp + i * 63 + 25);

        auto to_xxyzz_z = pbuffer.data(idx_op_hp + i * 63 + 26);

        auto to_xxzzz_x = pbuffer.data(idx_op_hp + i * 63 + 27);

        auto to_xxzzz_y = pbuffer.data(idx_op_hp + i * 63 + 28);

        auto to_xxzzz_z = pbuffer.data(idx_op_hp + i * 63 + 29);

        auto to_xyyyy_x = pbuffer.data(idx_op_hp + i * 63 + 30);

        auto to_xyyyy_y = pbuffer.data(idx_op_hp + i * 63 + 31);

        auto to_xyyyy_z = pbuffer.data(idx_op_hp + i * 63 + 32);

        auto to_xyyyz_x = pbuffer.data(idx_op_hp + i * 63 + 33);

        auto to_xyyyz_y = pbuffer.data(idx_op_hp + i * 63 + 34);

        auto to_xyyyz_z = pbuffer.data(idx_op_hp + i * 63 + 35);

        auto to_xyyzz_x = pbuffer.data(idx_op_hp + i * 63 + 36);

        auto to_xyyzz_y = pbuffer.data(idx_op_hp + i * 63 + 37);

        auto to_xyyzz_z = pbuffer.data(idx_op_hp + i * 63 + 38);

        auto to_xyzzz_x = pbuffer.data(idx_op_hp + i * 63 + 39);

        auto to_xyzzz_y = pbuffer.data(idx_op_hp + i * 63 + 40);

        auto to_xyzzz_z = pbuffer.data(idx_op_hp + i * 63 + 41);

        auto to_xzzzz_x = pbuffer.data(idx_op_hp + i * 63 + 42);

        auto to_xzzzz_y = pbuffer.data(idx_op_hp + i * 63 + 43);

        auto to_xzzzz_z = pbuffer.data(idx_op_hp + i * 63 + 44);

        auto to_yyyyy_x = pbuffer.data(idx_op_hp + i * 63 + 45);

        auto to_yyyyy_y = pbuffer.data(idx_op_hp + i * 63 + 46);

        auto to_yyyyy_z = pbuffer.data(idx_op_hp + i * 63 + 47);

        auto to_yyyyz_x = pbuffer.data(idx_op_hp + i * 63 + 48);

        auto to_yyyyz_y = pbuffer.data(idx_op_hp + i * 63 + 49);

        auto to_yyyyz_z = pbuffer.data(idx_op_hp + i * 63 + 50);

        auto to_yyyzz_x = pbuffer.data(idx_op_hp + i * 63 + 51);

        auto to_yyyzz_y = pbuffer.data(idx_op_hp + i * 63 + 52);

        auto to_yyyzz_z = pbuffer.data(idx_op_hp + i * 63 + 53);

        auto to_yyzzz_x = pbuffer.data(idx_op_hp + i * 63 + 54);

        auto to_yyzzz_y = pbuffer.data(idx_op_hp + i * 63 + 55);

        auto to_yyzzz_z = pbuffer.data(idx_op_hp + i * 63 + 56);

        auto to_yzzzz_x = pbuffer.data(idx_op_hp + i * 63 + 57);

        auto to_yzzzz_y = pbuffer.data(idx_op_hp + i * 63 + 58);

        auto to_yzzzz_z = pbuffer.data(idx_op_hp + i * 63 + 59);

        auto to_zzzzz_x = pbuffer.data(idx_op_hp + i * 63 + 60);

        auto to_zzzzz_y = pbuffer.data(idx_op_hp + i * 63 + 61);

        auto to_zzzzz_z = pbuffer.data(idx_op_hp + i * 63 + 62);

        // Set up 0-15 components of targeted buffer : GS

        auto to_x_x_xxxx_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 0);

        auto to_x_x_xxxy_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 1);

        auto to_x_x_xxxz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 2);

        auto to_x_x_xxyy_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 3);

        auto to_x_x_xxyz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 4);

        auto to_x_x_xxzz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 5);

        auto to_x_x_xyyy_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 6);

        auto to_x_x_xyyz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 7);

        auto to_x_x_xyzz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 8);

        auto to_x_x_xzzz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 9);

        auto to_x_x_yyyy_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 10);

        auto to_x_x_yyyz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 11);

        auto to_x_x_yyzz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 12);

        auto to_x_x_yzzz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 13);

        auto to_x_x_zzzz_0 = pbuffer.data(idx_op_geom_101_gs + 0 * op_comps * 15 + i * 15 + 14);

#pragma omp simd aligned(to_x_x_xxxx_0,     \
                             to_x_x_xxxy_0, \
                             to_x_x_xxxz_0, \
                             to_x_x_xxyy_0, \
                             to_x_x_xxyz_0, \
                             to_x_x_xxzz_0, \
                             to_x_x_xyyy_0, \
                             to_x_x_xyyz_0, \
                             to_x_x_xyzz_0, \
                             to_x_x_xzzz_0, \
                             to_x_x_yyyy_0, \
                             to_x_x_yyyz_0, \
                             to_x_x_yyzz_0, \
                             to_x_x_yzzz_0, \
                             to_x_x_zzzz_0, \
                             to_xxx_x,      \
                             to_xxxxx_x,    \
                             to_xxxxy_x,    \
                             to_xxxxz_x,    \
                             to_xxxyy_x,    \
                             to_xxxyz_x,    \
                             to_xxxzz_x,    \
                             to_xxy_x,      \
                             to_xxyyy_x,    \
                             to_xxyyz_x,    \
                             to_xxyzz_x,    \
                             to_xxz_x,      \
                             to_xxzzz_x,    \
                             to_xyy_x,      \
                             to_xyyyy_x,    \
                             to_xyyyz_x,    \
                             to_xyyzz_x,    \
                             to_xyz_x,      \
                             to_xyzzz_x,    \
                             to_xzz_x,      \
                             to_xzzzz_x,    \
                             to_yyy_x,      \
                             to_yyz_x,      \
                             to_yzz_x,      \
                             to_zzz_x,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxx_0[k] = -8.0 * to_xxx_x[k] * tke_0 + 4.0 * to_xxxxx_x[k] * tbe_0 * tke_0;

            to_x_x_xxxy_0[k] = -6.0 * to_xxy_x[k] * tke_0 + 4.0 * to_xxxxy_x[k] * tbe_0 * tke_0;

            to_x_x_xxxz_0[k] = -6.0 * to_xxz_x[k] * tke_0 + 4.0 * to_xxxxz_x[k] * tbe_0 * tke_0;

            to_x_x_xxyy_0[k] = -4.0 * to_xyy_x[k] * tke_0 + 4.0 * to_xxxyy_x[k] * tbe_0 * tke_0;

            to_x_x_xxyz_0[k] = -4.0 * to_xyz_x[k] * tke_0 + 4.0 * to_xxxyz_x[k] * tbe_0 * tke_0;

            to_x_x_xxzz_0[k] = -4.0 * to_xzz_x[k] * tke_0 + 4.0 * to_xxxzz_x[k] * tbe_0 * tke_0;

            to_x_x_xyyy_0[k] = -2.0 * to_yyy_x[k] * tke_0 + 4.0 * to_xxyyy_x[k] * tbe_0 * tke_0;

            to_x_x_xyyz_0[k] = -2.0 * to_yyz_x[k] * tke_0 + 4.0 * to_xxyyz_x[k] * tbe_0 * tke_0;

            to_x_x_xyzz_0[k] = -2.0 * to_yzz_x[k] * tke_0 + 4.0 * to_xxyzz_x[k] * tbe_0 * tke_0;

            to_x_x_xzzz_0[k] = -2.0 * to_zzz_x[k] * tke_0 + 4.0 * to_xxzzz_x[k] * tbe_0 * tke_0;

            to_x_x_yyyy_0[k] = 4.0 * to_xyyyy_x[k] * tbe_0 * tke_0;

            to_x_x_yyyz_0[k] = 4.0 * to_xyyyz_x[k] * tbe_0 * tke_0;

            to_x_x_yyzz_0[k] = 4.0 * to_xyyzz_x[k] * tbe_0 * tke_0;

            to_x_x_yzzz_0[k] = 4.0 * to_xyzzz_x[k] * tbe_0 * tke_0;

            to_x_x_zzzz_0[k] = 4.0 * to_xzzzz_x[k] * tbe_0 * tke_0;
        }

        // Set up 15-30 components of targeted buffer : GS

        auto to_x_y_xxxx_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 0);

        auto to_x_y_xxxy_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 1);

        auto to_x_y_xxxz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 2);

        auto to_x_y_xxyy_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 3);

        auto to_x_y_xxyz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 4);

        auto to_x_y_xxzz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 5);

        auto to_x_y_xyyy_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 6);

        auto to_x_y_xyyz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 7);

        auto to_x_y_xyzz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 8);

        auto to_x_y_xzzz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 9);

        auto to_x_y_yyyy_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 10);

        auto to_x_y_yyyz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 11);

        auto to_x_y_yyzz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 12);

        auto to_x_y_yzzz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 13);

        auto to_x_y_zzzz_0 = pbuffer.data(idx_op_geom_101_gs + 1 * op_comps * 15 + i * 15 + 14);

#pragma omp simd aligned(to_x_y_xxxx_0,     \
                             to_x_y_xxxy_0, \
                             to_x_y_xxxz_0, \
                             to_x_y_xxyy_0, \
                             to_x_y_xxyz_0, \
                             to_x_y_xxzz_0, \
                             to_x_y_xyyy_0, \
                             to_x_y_xyyz_0, \
                             to_x_y_xyzz_0, \
                             to_x_y_xzzz_0, \
                             to_x_y_yyyy_0, \
                             to_x_y_yyyz_0, \
                             to_x_y_yyzz_0, \
                             to_x_y_yzzz_0, \
                             to_x_y_zzzz_0, \
                             to_xxx_y,      \
                             to_xxxxx_y,    \
                             to_xxxxy_y,    \
                             to_xxxxz_y,    \
                             to_xxxyy_y,    \
                             to_xxxyz_y,    \
                             to_xxxzz_y,    \
                             to_xxy_y,      \
                             to_xxyyy_y,    \
                             to_xxyyz_y,    \
                             to_xxyzz_y,    \
                             to_xxz_y,      \
                             to_xxzzz_y,    \
                             to_xyy_y,      \
                             to_xyyyy_y,    \
                             to_xyyyz_y,    \
                             to_xyyzz_y,    \
                             to_xyz_y,      \
                             to_xyzzz_y,    \
                             to_xzz_y,      \
                             to_xzzzz_y,    \
                             to_yyy_y,      \
                             to_yyz_y,      \
                             to_yzz_y,      \
                             to_zzz_y,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxx_0[k] = -8.0 * to_xxx_y[k] * tke_0 + 4.0 * to_xxxxx_y[k] * tbe_0 * tke_0;

            to_x_y_xxxy_0[k] = -6.0 * to_xxy_y[k] * tke_0 + 4.0 * to_xxxxy_y[k] * tbe_0 * tke_0;

            to_x_y_xxxz_0[k] = -6.0 * to_xxz_y[k] * tke_0 + 4.0 * to_xxxxz_y[k] * tbe_0 * tke_0;

            to_x_y_xxyy_0[k] = -4.0 * to_xyy_y[k] * tke_0 + 4.0 * to_xxxyy_y[k] * tbe_0 * tke_0;

            to_x_y_xxyz_0[k] = -4.0 * to_xyz_y[k] * tke_0 + 4.0 * to_xxxyz_y[k] * tbe_0 * tke_0;

            to_x_y_xxzz_0[k] = -4.0 * to_xzz_y[k] * tke_0 + 4.0 * to_xxxzz_y[k] * tbe_0 * tke_0;

            to_x_y_xyyy_0[k] = -2.0 * to_yyy_y[k] * tke_0 + 4.0 * to_xxyyy_y[k] * tbe_0 * tke_0;

            to_x_y_xyyz_0[k] = -2.0 * to_yyz_y[k] * tke_0 + 4.0 * to_xxyyz_y[k] * tbe_0 * tke_0;

            to_x_y_xyzz_0[k] = -2.0 * to_yzz_y[k] * tke_0 + 4.0 * to_xxyzz_y[k] * tbe_0 * tke_0;

            to_x_y_xzzz_0[k] = -2.0 * to_zzz_y[k] * tke_0 + 4.0 * to_xxzzz_y[k] * tbe_0 * tke_0;

            to_x_y_yyyy_0[k] = 4.0 * to_xyyyy_y[k] * tbe_0 * tke_0;

            to_x_y_yyyz_0[k] = 4.0 * to_xyyyz_y[k] * tbe_0 * tke_0;

            to_x_y_yyzz_0[k] = 4.0 * to_xyyzz_y[k] * tbe_0 * tke_0;

            to_x_y_yzzz_0[k] = 4.0 * to_xyzzz_y[k] * tbe_0 * tke_0;

            to_x_y_zzzz_0[k] = 4.0 * to_xzzzz_y[k] * tbe_0 * tke_0;
        }

        // Set up 30-45 components of targeted buffer : GS

        auto to_x_z_xxxx_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 0);

        auto to_x_z_xxxy_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 1);

        auto to_x_z_xxxz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 2);

        auto to_x_z_xxyy_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 3);

        auto to_x_z_xxyz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 4);

        auto to_x_z_xxzz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 5);

        auto to_x_z_xyyy_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 6);

        auto to_x_z_xyyz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 7);

        auto to_x_z_xyzz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 8);

        auto to_x_z_xzzz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 9);

        auto to_x_z_yyyy_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 10);

        auto to_x_z_yyyz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 11);

        auto to_x_z_yyzz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 12);

        auto to_x_z_yzzz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 13);

        auto to_x_z_zzzz_0 = pbuffer.data(idx_op_geom_101_gs + 2 * op_comps * 15 + i * 15 + 14);

#pragma omp simd aligned(to_x_z_xxxx_0,     \
                             to_x_z_xxxy_0, \
                             to_x_z_xxxz_0, \
                             to_x_z_xxyy_0, \
                             to_x_z_xxyz_0, \
                             to_x_z_xxzz_0, \
                             to_x_z_xyyy_0, \
                             to_x_z_xyyz_0, \
                             to_x_z_xyzz_0, \
                             to_x_z_xzzz_0, \
                             to_x_z_yyyy_0, \
                             to_x_z_yyyz_0, \
                             to_x_z_yyzz_0, \
                             to_x_z_yzzz_0, \
                             to_x_z_zzzz_0, \
                             to_xxx_z,      \
                             to_xxxxx_z,    \
                             to_xxxxy_z,    \
                             to_xxxxz_z,    \
                             to_xxxyy_z,    \
                             to_xxxyz_z,    \
                             to_xxxzz_z,    \
                             to_xxy_z,      \
                             to_xxyyy_z,    \
                             to_xxyyz_z,    \
                             to_xxyzz_z,    \
                             to_xxz_z,      \
                             to_xxzzz_z,    \
                             to_xyy_z,      \
                             to_xyyyy_z,    \
                             to_xyyyz_z,    \
                             to_xyyzz_z,    \
                             to_xyz_z,      \
                             to_xyzzz_z,    \
                             to_xzz_z,      \
                             to_xzzzz_z,    \
                             to_yyy_z,      \
                             to_yyz_z,      \
                             to_yzz_z,      \
                             to_zzz_z,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxx_0[k] = -8.0 * to_xxx_z[k] * tke_0 + 4.0 * to_xxxxx_z[k] * tbe_0 * tke_0;

            to_x_z_xxxy_0[k] = -6.0 * to_xxy_z[k] * tke_0 + 4.0 * to_xxxxy_z[k] * tbe_0 * tke_0;

            to_x_z_xxxz_0[k] = -6.0 * to_xxz_z[k] * tke_0 + 4.0 * to_xxxxz_z[k] * tbe_0 * tke_0;

            to_x_z_xxyy_0[k] = -4.0 * to_xyy_z[k] * tke_0 + 4.0 * to_xxxyy_z[k] * tbe_0 * tke_0;

            to_x_z_xxyz_0[k] = -4.0 * to_xyz_z[k] * tke_0 + 4.0 * to_xxxyz_z[k] * tbe_0 * tke_0;

            to_x_z_xxzz_0[k] = -4.0 * to_xzz_z[k] * tke_0 + 4.0 * to_xxxzz_z[k] * tbe_0 * tke_0;

            to_x_z_xyyy_0[k] = -2.0 * to_yyy_z[k] * tke_0 + 4.0 * to_xxyyy_z[k] * tbe_0 * tke_0;

            to_x_z_xyyz_0[k] = -2.0 * to_yyz_z[k] * tke_0 + 4.0 * to_xxyyz_z[k] * tbe_0 * tke_0;

            to_x_z_xyzz_0[k] = -2.0 * to_yzz_z[k] * tke_0 + 4.0 * to_xxyzz_z[k] * tbe_0 * tke_0;

            to_x_z_xzzz_0[k] = -2.0 * to_zzz_z[k] * tke_0 + 4.0 * to_xxzzz_z[k] * tbe_0 * tke_0;

            to_x_z_yyyy_0[k] = 4.0 * to_xyyyy_z[k] * tbe_0 * tke_0;

            to_x_z_yyyz_0[k] = 4.0 * to_xyyyz_z[k] * tbe_0 * tke_0;

            to_x_z_yyzz_0[k] = 4.0 * to_xyyzz_z[k] * tbe_0 * tke_0;

            to_x_z_yzzz_0[k] = 4.0 * to_xyzzz_z[k] * tbe_0 * tke_0;

            to_x_z_zzzz_0[k] = 4.0 * to_xzzzz_z[k] * tbe_0 * tke_0;
        }

        // Set up 45-60 components of targeted buffer : GS

        auto to_y_x_xxxx_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 0);

        auto to_y_x_xxxy_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 1);

        auto to_y_x_xxxz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 2);

        auto to_y_x_xxyy_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 3);

        auto to_y_x_xxyz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 4);

        auto to_y_x_xxzz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 5);

        auto to_y_x_xyyy_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 6);

        auto to_y_x_xyyz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 7);

        auto to_y_x_xyzz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 8);

        auto to_y_x_xzzz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 9);

        auto to_y_x_yyyy_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 10);

        auto to_y_x_yyyz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 11);

        auto to_y_x_yyzz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 12);

        auto to_y_x_yzzz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 13);

        auto to_y_x_zzzz_0 = pbuffer.data(idx_op_geom_101_gs + 3 * op_comps * 15 + i * 15 + 14);

#pragma omp simd aligned(to_xxx_x,          \
                             to_xxxxy_x,    \
                             to_xxxyy_x,    \
                             to_xxxyz_x,    \
                             to_xxy_x,      \
                             to_xxyyy_x,    \
                             to_xxyyz_x,    \
                             to_xxyzz_x,    \
                             to_xxz_x,      \
                             to_xyy_x,      \
                             to_xyyyy_x,    \
                             to_xyyyz_x,    \
                             to_xyyzz_x,    \
                             to_xyz_x,      \
                             to_xyzzz_x,    \
                             to_xzz_x,      \
                             to_y_x_xxxx_0, \
                             to_y_x_xxxy_0, \
                             to_y_x_xxxz_0, \
                             to_y_x_xxyy_0, \
                             to_y_x_xxyz_0, \
                             to_y_x_xxzz_0, \
                             to_y_x_xyyy_0, \
                             to_y_x_xyyz_0, \
                             to_y_x_xyzz_0, \
                             to_y_x_xzzz_0, \
                             to_y_x_yyyy_0, \
                             to_y_x_yyyz_0, \
                             to_y_x_yyzz_0, \
                             to_y_x_yzzz_0, \
                             to_y_x_zzzz_0, \
                             to_yyy_x,      \
                             to_yyyyy_x,    \
                             to_yyyyz_x,    \
                             to_yyyzz_x,    \
                             to_yyz_x,      \
                             to_yyzzz_x,    \
                             to_yzz_x,      \
                             to_yzzzz_x,    \
                             to_zzz_x,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxx_0[k] = 4.0 * to_xxxxy_x[k] * tbe_0 * tke_0;

            to_y_x_xxxy_0[k] = -2.0 * to_xxx_x[k] * tke_0 + 4.0 * to_xxxyy_x[k] * tbe_0 * tke_0;

            to_y_x_xxxz_0[k] = 4.0 * to_xxxyz_x[k] * tbe_0 * tke_0;

            to_y_x_xxyy_0[k] = -4.0 * to_xxy_x[k] * tke_0 + 4.0 * to_xxyyy_x[k] * tbe_0 * tke_0;

            to_y_x_xxyz_0[k] = -2.0 * to_xxz_x[k] * tke_0 + 4.0 * to_xxyyz_x[k] * tbe_0 * tke_0;

            to_y_x_xxzz_0[k] = 4.0 * to_xxyzz_x[k] * tbe_0 * tke_0;

            to_y_x_xyyy_0[k] = -6.0 * to_xyy_x[k] * tke_0 + 4.0 * to_xyyyy_x[k] * tbe_0 * tke_0;

            to_y_x_xyyz_0[k] = -4.0 * to_xyz_x[k] * tke_0 + 4.0 * to_xyyyz_x[k] * tbe_0 * tke_0;

            to_y_x_xyzz_0[k] = -2.0 * to_xzz_x[k] * tke_0 + 4.0 * to_xyyzz_x[k] * tbe_0 * tke_0;

            to_y_x_xzzz_0[k] = 4.0 * to_xyzzz_x[k] * tbe_0 * tke_0;

            to_y_x_yyyy_0[k] = -8.0 * to_yyy_x[k] * tke_0 + 4.0 * to_yyyyy_x[k] * tbe_0 * tke_0;

            to_y_x_yyyz_0[k] = -6.0 * to_yyz_x[k] * tke_0 + 4.0 * to_yyyyz_x[k] * tbe_0 * tke_0;

            to_y_x_yyzz_0[k] = -4.0 * to_yzz_x[k] * tke_0 + 4.0 * to_yyyzz_x[k] * tbe_0 * tke_0;

            to_y_x_yzzz_0[k] = -2.0 * to_zzz_x[k] * tke_0 + 4.0 * to_yyzzz_x[k] * tbe_0 * tke_0;

            to_y_x_zzzz_0[k] = 4.0 * to_yzzzz_x[k] * tbe_0 * tke_0;
        }

        // Set up 60-75 components of targeted buffer : GS

        auto to_y_y_xxxx_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 0);

        auto to_y_y_xxxy_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 1);

        auto to_y_y_xxxz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 2);

        auto to_y_y_xxyy_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 3);

        auto to_y_y_xxyz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 4);

        auto to_y_y_xxzz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 5);

        auto to_y_y_xyyy_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 6);

        auto to_y_y_xyyz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 7);

        auto to_y_y_xyzz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 8);

        auto to_y_y_xzzz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 9);

        auto to_y_y_yyyy_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 10);

        auto to_y_y_yyyz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 11);

        auto to_y_y_yyzz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 12);

        auto to_y_y_yzzz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 13);

        auto to_y_y_zzzz_0 = pbuffer.data(idx_op_geom_101_gs + 4 * op_comps * 15 + i * 15 + 14);

#pragma omp simd aligned(to_xxx_y,          \
                             to_xxxxy_y,    \
                             to_xxxyy_y,    \
                             to_xxxyz_y,    \
                             to_xxy_y,      \
                             to_xxyyy_y,    \
                             to_xxyyz_y,    \
                             to_xxyzz_y,    \
                             to_xxz_y,      \
                             to_xyy_y,      \
                             to_xyyyy_y,    \
                             to_xyyyz_y,    \
                             to_xyyzz_y,    \
                             to_xyz_y,      \
                             to_xyzzz_y,    \
                             to_xzz_y,      \
                             to_y_y_xxxx_0, \
                             to_y_y_xxxy_0, \
                             to_y_y_xxxz_0, \
                             to_y_y_xxyy_0, \
                             to_y_y_xxyz_0, \
                             to_y_y_xxzz_0, \
                             to_y_y_xyyy_0, \
                             to_y_y_xyyz_0, \
                             to_y_y_xyzz_0, \
                             to_y_y_xzzz_0, \
                             to_y_y_yyyy_0, \
                             to_y_y_yyyz_0, \
                             to_y_y_yyzz_0, \
                             to_y_y_yzzz_0, \
                             to_y_y_zzzz_0, \
                             to_yyy_y,      \
                             to_yyyyy_y,    \
                             to_yyyyz_y,    \
                             to_yyyzz_y,    \
                             to_yyz_y,      \
                             to_yyzzz_y,    \
                             to_yzz_y,      \
                             to_yzzzz_y,    \
                             to_zzz_y,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxx_0[k] = 4.0 * to_xxxxy_y[k] * tbe_0 * tke_0;

            to_y_y_xxxy_0[k] = -2.0 * to_xxx_y[k] * tke_0 + 4.0 * to_xxxyy_y[k] * tbe_0 * tke_0;

            to_y_y_xxxz_0[k] = 4.0 * to_xxxyz_y[k] * tbe_0 * tke_0;

            to_y_y_xxyy_0[k] = -4.0 * to_xxy_y[k] * tke_0 + 4.0 * to_xxyyy_y[k] * tbe_0 * tke_0;

            to_y_y_xxyz_0[k] = -2.0 * to_xxz_y[k] * tke_0 + 4.0 * to_xxyyz_y[k] * tbe_0 * tke_0;

            to_y_y_xxzz_0[k] = 4.0 * to_xxyzz_y[k] * tbe_0 * tke_0;

            to_y_y_xyyy_0[k] = -6.0 * to_xyy_y[k] * tke_0 + 4.0 * to_xyyyy_y[k] * tbe_0 * tke_0;

            to_y_y_xyyz_0[k] = -4.0 * to_xyz_y[k] * tke_0 + 4.0 * to_xyyyz_y[k] * tbe_0 * tke_0;

            to_y_y_xyzz_0[k] = -2.0 * to_xzz_y[k] * tke_0 + 4.0 * to_xyyzz_y[k] * tbe_0 * tke_0;

            to_y_y_xzzz_0[k] = 4.0 * to_xyzzz_y[k] * tbe_0 * tke_0;

            to_y_y_yyyy_0[k] = -8.0 * to_yyy_y[k] * tke_0 + 4.0 * to_yyyyy_y[k] * tbe_0 * tke_0;

            to_y_y_yyyz_0[k] = -6.0 * to_yyz_y[k] * tke_0 + 4.0 * to_yyyyz_y[k] * tbe_0 * tke_0;

            to_y_y_yyzz_0[k] = -4.0 * to_yzz_y[k] * tke_0 + 4.0 * to_yyyzz_y[k] * tbe_0 * tke_0;

            to_y_y_yzzz_0[k] = -2.0 * to_zzz_y[k] * tke_0 + 4.0 * to_yyzzz_y[k] * tbe_0 * tke_0;

            to_y_y_zzzz_0[k] = 4.0 * to_yzzzz_y[k] * tbe_0 * tke_0;
        }

        // Set up 75-90 components of targeted buffer : GS

        auto to_y_z_xxxx_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 0);

        auto to_y_z_xxxy_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 1);

        auto to_y_z_xxxz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 2);

        auto to_y_z_xxyy_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 3);

        auto to_y_z_xxyz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 4);

        auto to_y_z_xxzz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 5);

        auto to_y_z_xyyy_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 6);

        auto to_y_z_xyyz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 7);

        auto to_y_z_xyzz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 8);

        auto to_y_z_xzzz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 9);

        auto to_y_z_yyyy_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 10);

        auto to_y_z_yyyz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 11);

        auto to_y_z_yyzz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 12);

        auto to_y_z_yzzz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 13);

        auto to_y_z_zzzz_0 = pbuffer.data(idx_op_geom_101_gs + 5 * op_comps * 15 + i * 15 + 14);

#pragma omp simd aligned(to_xxx_z,          \
                             to_xxxxy_z,    \
                             to_xxxyy_z,    \
                             to_xxxyz_z,    \
                             to_xxy_z,      \
                             to_xxyyy_z,    \
                             to_xxyyz_z,    \
                             to_xxyzz_z,    \
                             to_xxz_z,      \
                             to_xyy_z,      \
                             to_xyyyy_z,    \
                             to_xyyyz_z,    \
                             to_xyyzz_z,    \
                             to_xyz_z,      \
                             to_xyzzz_z,    \
                             to_xzz_z,      \
                             to_y_z_xxxx_0, \
                             to_y_z_xxxy_0, \
                             to_y_z_xxxz_0, \
                             to_y_z_xxyy_0, \
                             to_y_z_xxyz_0, \
                             to_y_z_xxzz_0, \
                             to_y_z_xyyy_0, \
                             to_y_z_xyyz_0, \
                             to_y_z_xyzz_0, \
                             to_y_z_xzzz_0, \
                             to_y_z_yyyy_0, \
                             to_y_z_yyyz_0, \
                             to_y_z_yyzz_0, \
                             to_y_z_yzzz_0, \
                             to_y_z_zzzz_0, \
                             to_yyy_z,      \
                             to_yyyyy_z,    \
                             to_yyyyz_z,    \
                             to_yyyzz_z,    \
                             to_yyz_z,      \
                             to_yyzzz_z,    \
                             to_yzz_z,      \
                             to_yzzzz_z,    \
                             to_zzz_z,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxx_0[k] = 4.0 * to_xxxxy_z[k] * tbe_0 * tke_0;

            to_y_z_xxxy_0[k] = -2.0 * to_xxx_z[k] * tke_0 + 4.0 * to_xxxyy_z[k] * tbe_0 * tke_0;

            to_y_z_xxxz_0[k] = 4.0 * to_xxxyz_z[k] * tbe_0 * tke_0;

            to_y_z_xxyy_0[k] = -4.0 * to_xxy_z[k] * tke_0 + 4.0 * to_xxyyy_z[k] * tbe_0 * tke_0;

            to_y_z_xxyz_0[k] = -2.0 * to_xxz_z[k] * tke_0 + 4.0 * to_xxyyz_z[k] * tbe_0 * tke_0;

            to_y_z_xxzz_0[k] = 4.0 * to_xxyzz_z[k] * tbe_0 * tke_0;

            to_y_z_xyyy_0[k] = -6.0 * to_xyy_z[k] * tke_0 + 4.0 * to_xyyyy_z[k] * tbe_0 * tke_0;

            to_y_z_xyyz_0[k] = -4.0 * to_xyz_z[k] * tke_0 + 4.0 * to_xyyyz_z[k] * tbe_0 * tke_0;

            to_y_z_xyzz_0[k] = -2.0 * to_xzz_z[k] * tke_0 + 4.0 * to_xyyzz_z[k] * tbe_0 * tke_0;

            to_y_z_xzzz_0[k] = 4.0 * to_xyzzz_z[k] * tbe_0 * tke_0;

            to_y_z_yyyy_0[k] = -8.0 * to_yyy_z[k] * tke_0 + 4.0 * to_yyyyy_z[k] * tbe_0 * tke_0;

            to_y_z_yyyz_0[k] = -6.0 * to_yyz_z[k] * tke_0 + 4.0 * to_yyyyz_z[k] * tbe_0 * tke_0;

            to_y_z_yyzz_0[k] = -4.0 * to_yzz_z[k] * tke_0 + 4.0 * to_yyyzz_z[k] * tbe_0 * tke_0;

            to_y_z_yzzz_0[k] = -2.0 * to_zzz_z[k] * tke_0 + 4.0 * to_yyzzz_z[k] * tbe_0 * tke_0;

            to_y_z_zzzz_0[k] = 4.0 * to_yzzzz_z[k] * tbe_0 * tke_0;
        }

        // Set up 90-105 components of targeted buffer : GS

        auto to_z_x_xxxx_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 0);

        auto to_z_x_xxxy_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 1);

        auto to_z_x_xxxz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 2);

        auto to_z_x_xxyy_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 3);

        auto to_z_x_xxyz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 4);

        auto to_z_x_xxzz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 5);

        auto to_z_x_xyyy_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 6);

        auto to_z_x_xyyz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 7);

        auto to_z_x_xyzz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 8);

        auto to_z_x_xzzz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 9);

        auto to_z_x_yyyy_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 10);

        auto to_z_x_yyyz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 11);

        auto to_z_x_yyzz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 12);

        auto to_z_x_yzzz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 13);

        auto to_z_x_zzzz_0 = pbuffer.data(idx_op_geom_101_gs + 6 * op_comps * 15 + i * 15 + 14);

#pragma omp simd aligned(to_xxx_x,          \
                             to_xxxxz_x,    \
                             to_xxxyz_x,    \
                             to_xxxzz_x,    \
                             to_xxy_x,      \
                             to_xxyyz_x,    \
                             to_xxyzz_x,    \
                             to_xxz_x,      \
                             to_xxzzz_x,    \
                             to_xyy_x,      \
                             to_xyyyz_x,    \
                             to_xyyzz_x,    \
                             to_xyz_x,      \
                             to_xyzzz_x,    \
                             to_xzz_x,      \
                             to_xzzzz_x,    \
                             to_yyy_x,      \
                             to_yyyyz_x,    \
                             to_yyyzz_x,    \
                             to_yyz_x,      \
                             to_yyzzz_x,    \
                             to_yzz_x,      \
                             to_yzzzz_x,    \
                             to_z_x_xxxx_0, \
                             to_z_x_xxxy_0, \
                             to_z_x_xxxz_0, \
                             to_z_x_xxyy_0, \
                             to_z_x_xxyz_0, \
                             to_z_x_xxzz_0, \
                             to_z_x_xyyy_0, \
                             to_z_x_xyyz_0, \
                             to_z_x_xyzz_0, \
                             to_z_x_xzzz_0, \
                             to_z_x_yyyy_0, \
                             to_z_x_yyyz_0, \
                             to_z_x_yyzz_0, \
                             to_z_x_yzzz_0, \
                             to_z_x_zzzz_0, \
                             to_zzz_x,      \
                             to_zzzzz_x,    \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxx_0[k] = 4.0 * to_xxxxz_x[k] * tbe_0 * tke_0;

            to_z_x_xxxy_0[k] = 4.0 * to_xxxyz_x[k] * tbe_0 * tke_0;

            to_z_x_xxxz_0[k] = -2.0 * to_xxx_x[k] * tke_0 + 4.0 * to_xxxzz_x[k] * tbe_0 * tke_0;

            to_z_x_xxyy_0[k] = 4.0 * to_xxyyz_x[k] * tbe_0 * tke_0;

            to_z_x_xxyz_0[k] = -2.0 * to_xxy_x[k] * tke_0 + 4.0 * to_xxyzz_x[k] * tbe_0 * tke_0;

            to_z_x_xxzz_0[k] = -4.0 * to_xxz_x[k] * tke_0 + 4.0 * to_xxzzz_x[k] * tbe_0 * tke_0;

            to_z_x_xyyy_0[k] = 4.0 * to_xyyyz_x[k] * tbe_0 * tke_0;

            to_z_x_xyyz_0[k] = -2.0 * to_xyy_x[k] * tke_0 + 4.0 * to_xyyzz_x[k] * tbe_0 * tke_0;

            to_z_x_xyzz_0[k] = -4.0 * to_xyz_x[k] * tke_0 + 4.0 * to_xyzzz_x[k] * tbe_0 * tke_0;

            to_z_x_xzzz_0[k] = -6.0 * to_xzz_x[k] * tke_0 + 4.0 * to_xzzzz_x[k] * tbe_0 * tke_0;

            to_z_x_yyyy_0[k] = 4.0 * to_yyyyz_x[k] * tbe_0 * tke_0;

            to_z_x_yyyz_0[k] = -2.0 * to_yyy_x[k] * tke_0 + 4.0 * to_yyyzz_x[k] * tbe_0 * tke_0;

            to_z_x_yyzz_0[k] = -4.0 * to_yyz_x[k] * tke_0 + 4.0 * to_yyzzz_x[k] * tbe_0 * tke_0;

            to_z_x_yzzz_0[k] = -6.0 * to_yzz_x[k] * tke_0 + 4.0 * to_yzzzz_x[k] * tbe_0 * tke_0;

            to_z_x_zzzz_0[k] = -8.0 * to_zzz_x[k] * tke_0 + 4.0 * to_zzzzz_x[k] * tbe_0 * tke_0;
        }

        // Set up 105-120 components of targeted buffer : GS

        auto to_z_y_xxxx_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 0);

        auto to_z_y_xxxy_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 1);

        auto to_z_y_xxxz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 2);

        auto to_z_y_xxyy_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 3);

        auto to_z_y_xxyz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 4);

        auto to_z_y_xxzz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 5);

        auto to_z_y_xyyy_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 6);

        auto to_z_y_xyyz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 7);

        auto to_z_y_xyzz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 8);

        auto to_z_y_xzzz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 9);

        auto to_z_y_yyyy_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 10);

        auto to_z_y_yyyz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 11);

        auto to_z_y_yyzz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 12);

        auto to_z_y_yzzz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 13);

        auto to_z_y_zzzz_0 = pbuffer.data(idx_op_geom_101_gs + 7 * op_comps * 15 + i * 15 + 14);

#pragma omp simd aligned(to_xxx_y,          \
                             to_xxxxz_y,    \
                             to_xxxyz_y,    \
                             to_xxxzz_y,    \
                             to_xxy_y,      \
                             to_xxyyz_y,    \
                             to_xxyzz_y,    \
                             to_xxz_y,      \
                             to_xxzzz_y,    \
                             to_xyy_y,      \
                             to_xyyyz_y,    \
                             to_xyyzz_y,    \
                             to_xyz_y,      \
                             to_xyzzz_y,    \
                             to_xzz_y,      \
                             to_xzzzz_y,    \
                             to_yyy_y,      \
                             to_yyyyz_y,    \
                             to_yyyzz_y,    \
                             to_yyz_y,      \
                             to_yyzzz_y,    \
                             to_yzz_y,      \
                             to_yzzzz_y,    \
                             to_z_y_xxxx_0, \
                             to_z_y_xxxy_0, \
                             to_z_y_xxxz_0, \
                             to_z_y_xxyy_0, \
                             to_z_y_xxyz_0, \
                             to_z_y_xxzz_0, \
                             to_z_y_xyyy_0, \
                             to_z_y_xyyz_0, \
                             to_z_y_xyzz_0, \
                             to_z_y_xzzz_0, \
                             to_z_y_yyyy_0, \
                             to_z_y_yyyz_0, \
                             to_z_y_yyzz_0, \
                             to_z_y_yzzz_0, \
                             to_z_y_zzzz_0, \
                             to_zzz_y,      \
                             to_zzzzz_y,    \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxx_0[k] = 4.0 * to_xxxxz_y[k] * tbe_0 * tke_0;

            to_z_y_xxxy_0[k] = 4.0 * to_xxxyz_y[k] * tbe_0 * tke_0;

            to_z_y_xxxz_0[k] = -2.0 * to_xxx_y[k] * tke_0 + 4.0 * to_xxxzz_y[k] * tbe_0 * tke_0;

            to_z_y_xxyy_0[k] = 4.0 * to_xxyyz_y[k] * tbe_0 * tke_0;

            to_z_y_xxyz_0[k] = -2.0 * to_xxy_y[k] * tke_0 + 4.0 * to_xxyzz_y[k] * tbe_0 * tke_0;

            to_z_y_xxzz_0[k] = -4.0 * to_xxz_y[k] * tke_0 + 4.0 * to_xxzzz_y[k] * tbe_0 * tke_0;

            to_z_y_xyyy_0[k] = 4.0 * to_xyyyz_y[k] * tbe_0 * tke_0;

            to_z_y_xyyz_0[k] = -2.0 * to_xyy_y[k] * tke_0 + 4.0 * to_xyyzz_y[k] * tbe_0 * tke_0;

            to_z_y_xyzz_0[k] = -4.0 * to_xyz_y[k] * tke_0 + 4.0 * to_xyzzz_y[k] * tbe_0 * tke_0;

            to_z_y_xzzz_0[k] = -6.0 * to_xzz_y[k] * tke_0 + 4.0 * to_xzzzz_y[k] * tbe_0 * tke_0;

            to_z_y_yyyy_0[k] = 4.0 * to_yyyyz_y[k] * tbe_0 * tke_0;

            to_z_y_yyyz_0[k] = -2.0 * to_yyy_y[k] * tke_0 + 4.0 * to_yyyzz_y[k] * tbe_0 * tke_0;

            to_z_y_yyzz_0[k] = -4.0 * to_yyz_y[k] * tke_0 + 4.0 * to_yyzzz_y[k] * tbe_0 * tke_0;

            to_z_y_yzzz_0[k] = -6.0 * to_yzz_y[k] * tke_0 + 4.0 * to_yzzzz_y[k] * tbe_0 * tke_0;

            to_z_y_zzzz_0[k] = -8.0 * to_zzz_y[k] * tke_0 + 4.0 * to_zzzzz_y[k] * tbe_0 * tke_0;
        }

        // Set up 120-135 components of targeted buffer : GS

        auto to_z_z_xxxx_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 0);

        auto to_z_z_xxxy_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 1);

        auto to_z_z_xxxz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 2);

        auto to_z_z_xxyy_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 3);

        auto to_z_z_xxyz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 4);

        auto to_z_z_xxzz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 5);

        auto to_z_z_xyyy_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 6);

        auto to_z_z_xyyz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 7);

        auto to_z_z_xyzz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 8);

        auto to_z_z_xzzz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 9);

        auto to_z_z_yyyy_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 10);

        auto to_z_z_yyyz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 11);

        auto to_z_z_yyzz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 12);

        auto to_z_z_yzzz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 13);

        auto to_z_z_zzzz_0 = pbuffer.data(idx_op_geom_101_gs + 8 * op_comps * 15 + i * 15 + 14);

#pragma omp simd aligned(to_xxx_z,          \
                             to_xxxxz_z,    \
                             to_xxxyz_z,    \
                             to_xxxzz_z,    \
                             to_xxy_z,      \
                             to_xxyyz_z,    \
                             to_xxyzz_z,    \
                             to_xxz_z,      \
                             to_xxzzz_z,    \
                             to_xyy_z,      \
                             to_xyyyz_z,    \
                             to_xyyzz_z,    \
                             to_xyz_z,      \
                             to_xyzzz_z,    \
                             to_xzz_z,      \
                             to_xzzzz_z,    \
                             to_yyy_z,      \
                             to_yyyyz_z,    \
                             to_yyyzz_z,    \
                             to_yyz_z,      \
                             to_yyzzz_z,    \
                             to_yzz_z,      \
                             to_yzzzz_z,    \
                             to_z_z_xxxx_0, \
                             to_z_z_xxxy_0, \
                             to_z_z_xxxz_0, \
                             to_z_z_xxyy_0, \
                             to_z_z_xxyz_0, \
                             to_z_z_xxzz_0, \
                             to_z_z_xyyy_0, \
                             to_z_z_xyyz_0, \
                             to_z_z_xyzz_0, \
                             to_z_z_xzzz_0, \
                             to_z_z_yyyy_0, \
                             to_z_z_yyyz_0, \
                             to_z_z_yyzz_0, \
                             to_z_z_yzzz_0, \
                             to_z_z_zzzz_0, \
                             to_zzz_z,      \
                             to_zzzzz_z,    \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxx_0[k] = 4.0 * to_xxxxz_z[k] * tbe_0 * tke_0;

            to_z_z_xxxy_0[k] = 4.0 * to_xxxyz_z[k] * tbe_0 * tke_0;

            to_z_z_xxxz_0[k] = -2.0 * to_xxx_z[k] * tke_0 + 4.0 * to_xxxzz_z[k] * tbe_0 * tke_0;

            to_z_z_xxyy_0[k] = 4.0 * to_xxyyz_z[k] * tbe_0 * tke_0;

            to_z_z_xxyz_0[k] = -2.0 * to_xxy_z[k] * tke_0 + 4.0 * to_xxyzz_z[k] * tbe_0 * tke_0;

            to_z_z_xxzz_0[k] = -4.0 * to_xxz_z[k] * tke_0 + 4.0 * to_xxzzz_z[k] * tbe_0 * tke_0;

            to_z_z_xyyy_0[k] = 4.0 * to_xyyyz_z[k] * tbe_0 * tke_0;

            to_z_z_xyyz_0[k] = -2.0 * to_xyy_z[k] * tke_0 + 4.0 * to_xyyzz_z[k] * tbe_0 * tke_0;

            to_z_z_xyzz_0[k] = -4.0 * to_xyz_z[k] * tke_0 + 4.0 * to_xyzzz_z[k] * tbe_0 * tke_0;

            to_z_z_xzzz_0[k] = -6.0 * to_xzz_z[k] * tke_0 + 4.0 * to_xzzzz_z[k] * tbe_0 * tke_0;

            to_z_z_yyyy_0[k] = 4.0 * to_yyyyz_z[k] * tbe_0 * tke_0;

            to_z_z_yyyz_0[k] = -2.0 * to_yyy_z[k] * tke_0 + 4.0 * to_yyyzz_z[k] * tbe_0 * tke_0;

            to_z_z_yyzz_0[k] = -4.0 * to_yyz_z[k] * tke_0 + 4.0 * to_yyzzz_z[k] * tbe_0 * tke_0;

            to_z_z_yzzz_0[k] = -6.0 * to_yzz_z[k] * tke_0 + 4.0 * to_yzzzz_z[k] * tbe_0 * tke_0;

            to_z_z_zzzz_0[k] = -8.0 * to_zzz_z[k] * tke_0 + 4.0 * to_zzzzz_z[k] * tbe_0 * tke_0;
        }
    }
}

}  // namespace t2cgeom
