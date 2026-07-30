#include "GeometricalDerivatives010ForGS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_gs(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_gs,
                         const int idx_op_fs,
                         const int idx_op_gp,
                         const int idx_op_hs,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : FS

    auto tr_xxx_0 = pbuffer.data(idx_op_fs);

    auto tr_xxy_0 = pbuffer.data(idx_op_fs + 1);

    auto tr_xxz_0 = pbuffer.data(idx_op_fs + 2);

    auto tr_xyy_0 = pbuffer.data(idx_op_fs + 3);

    auto tr_xyz_0 = pbuffer.data(idx_op_fs + 4);

    auto tr_xzz_0 = pbuffer.data(idx_op_fs + 5);

    auto tr_yyy_0 = pbuffer.data(idx_op_fs + 6);

    auto tr_yyz_0 = pbuffer.data(idx_op_fs + 7);

    auto tr_yzz_0 = pbuffer.data(idx_op_fs + 8);

    auto tr_zzz_0 = pbuffer.data(idx_op_fs + 9);

    // Set up components of auxiliary buffer : GP

    auto tr_xxxx_x = pbuffer.data(idx_op_gp);

    auto tr_xxxx_y = pbuffer.data(idx_op_gp + 1);

    auto tr_xxxx_z = pbuffer.data(idx_op_gp + 2);

    auto tr_xxxy_x = pbuffer.data(idx_op_gp + 3);

    auto tr_xxxy_y = pbuffer.data(idx_op_gp + 4);

    auto tr_xxxy_z = pbuffer.data(idx_op_gp + 5);

    auto tr_xxxz_x = pbuffer.data(idx_op_gp + 6);

    auto tr_xxxz_y = pbuffer.data(idx_op_gp + 7);

    auto tr_xxxz_z = pbuffer.data(idx_op_gp + 8);

    auto tr_xxyy_x = pbuffer.data(idx_op_gp + 9);

    auto tr_xxyy_y = pbuffer.data(idx_op_gp + 10);

    auto tr_xxyy_z = pbuffer.data(idx_op_gp + 11);

    auto tr_xxyz_x = pbuffer.data(idx_op_gp + 12);

    auto tr_xxyz_y = pbuffer.data(idx_op_gp + 13);

    auto tr_xxyz_z = pbuffer.data(idx_op_gp + 14);

    auto tr_xxzz_x = pbuffer.data(idx_op_gp + 15);

    auto tr_xxzz_y = pbuffer.data(idx_op_gp + 16);

    auto tr_xxzz_z = pbuffer.data(idx_op_gp + 17);

    auto tr_xyyy_x = pbuffer.data(idx_op_gp + 18);

    auto tr_xyyy_y = pbuffer.data(idx_op_gp + 19);

    auto tr_xyyy_z = pbuffer.data(idx_op_gp + 20);

    auto tr_xyyz_x = pbuffer.data(idx_op_gp + 21);

    auto tr_xyyz_y = pbuffer.data(idx_op_gp + 22);

    auto tr_xyyz_z = pbuffer.data(idx_op_gp + 23);

    auto tr_xyzz_x = pbuffer.data(idx_op_gp + 24);

    auto tr_xyzz_y = pbuffer.data(idx_op_gp + 25);

    auto tr_xyzz_z = pbuffer.data(idx_op_gp + 26);

    auto tr_xzzz_x = pbuffer.data(idx_op_gp + 27);

    auto tr_xzzz_y = pbuffer.data(idx_op_gp + 28);

    auto tr_xzzz_z = pbuffer.data(idx_op_gp + 29);

    auto tr_yyyy_x = pbuffer.data(idx_op_gp + 30);

    auto tr_yyyy_y = pbuffer.data(idx_op_gp + 31);

    auto tr_yyyy_z = pbuffer.data(idx_op_gp + 32);

    auto tr_yyyz_x = pbuffer.data(idx_op_gp + 33);

    auto tr_yyyz_y = pbuffer.data(idx_op_gp + 34);

    auto tr_yyyz_z = pbuffer.data(idx_op_gp + 35);

    auto tr_yyzz_x = pbuffer.data(idx_op_gp + 36);

    auto tr_yyzz_y = pbuffer.data(idx_op_gp + 37);

    auto tr_yyzz_z = pbuffer.data(idx_op_gp + 38);

    auto tr_yzzz_x = pbuffer.data(idx_op_gp + 39);

    auto tr_yzzz_y = pbuffer.data(idx_op_gp + 40);

    auto tr_yzzz_z = pbuffer.data(idx_op_gp + 41);

    auto tr_zzzz_x = pbuffer.data(idx_op_gp + 42);

    auto tr_zzzz_y = pbuffer.data(idx_op_gp + 43);

    auto tr_zzzz_z = pbuffer.data(idx_op_gp + 44);

    // Set up components of auxiliary buffer : HS

    auto tr_xxxxx_0 = pbuffer.data(idx_op_hs);

    auto tr_xxxxy_0 = pbuffer.data(idx_op_hs + 1);

    auto tr_xxxxz_0 = pbuffer.data(idx_op_hs + 2);

    auto tr_xxxyy_0 = pbuffer.data(idx_op_hs + 3);

    auto tr_xxxyz_0 = pbuffer.data(idx_op_hs + 4);

    auto tr_xxxzz_0 = pbuffer.data(idx_op_hs + 5);

    auto tr_xxyyy_0 = pbuffer.data(idx_op_hs + 6);

    auto tr_xxyyz_0 = pbuffer.data(idx_op_hs + 7);

    auto tr_xxyzz_0 = pbuffer.data(idx_op_hs + 8);

    auto tr_xxzzz_0 = pbuffer.data(idx_op_hs + 9);

    auto tr_xyyyy_0 = pbuffer.data(idx_op_hs + 10);

    auto tr_xyyyz_0 = pbuffer.data(idx_op_hs + 11);

    auto tr_xyyzz_0 = pbuffer.data(idx_op_hs + 12);

    auto tr_xyzzz_0 = pbuffer.data(idx_op_hs + 13);

    auto tr_xzzzz_0 = pbuffer.data(idx_op_hs + 14);

    auto tr_yyyyy_0 = pbuffer.data(idx_op_hs + 15);

    auto tr_yyyyz_0 = pbuffer.data(idx_op_hs + 16);

    auto tr_yyyzz_0 = pbuffer.data(idx_op_hs + 17);

    auto tr_yyzzz_0 = pbuffer.data(idx_op_hs + 18);

    auto tr_yzzzz_0 = pbuffer.data(idx_op_hs + 19);

    auto tr_zzzzz_0 = pbuffer.data(idx_op_hs + 20);

    // Set up components of targeted buffer : GS

    auto tr_0_0_x_xxxx_0 = pbuffer.data(idx_op_geom_010_gs);

    auto tr_0_0_x_xxxy_0 = pbuffer.data(idx_op_geom_010_gs + 1);

    auto tr_0_0_x_xxxz_0 = pbuffer.data(idx_op_geom_010_gs + 2);

    auto tr_0_0_x_xxyy_0 = pbuffer.data(idx_op_geom_010_gs + 3);

    auto tr_0_0_x_xxyz_0 = pbuffer.data(idx_op_geom_010_gs + 4);

    auto tr_0_0_x_xxzz_0 = pbuffer.data(idx_op_geom_010_gs + 5);

    auto tr_0_0_x_xyyy_0 = pbuffer.data(idx_op_geom_010_gs + 6);

    auto tr_0_0_x_xyyz_0 = pbuffer.data(idx_op_geom_010_gs + 7);

    auto tr_0_0_x_xyzz_0 = pbuffer.data(idx_op_geom_010_gs + 8);

    auto tr_0_0_x_xzzz_0 = pbuffer.data(idx_op_geom_010_gs + 9);

    auto tr_0_0_x_yyyy_0 = pbuffer.data(idx_op_geom_010_gs + 10);

    auto tr_0_0_x_yyyz_0 = pbuffer.data(idx_op_geom_010_gs + 11);

    auto tr_0_0_x_yyzz_0 = pbuffer.data(idx_op_geom_010_gs + 12);

    auto tr_0_0_x_yzzz_0 = pbuffer.data(idx_op_geom_010_gs + 13);

    auto tr_0_0_x_zzzz_0 = pbuffer.data(idx_op_geom_010_gs + 14);

    auto tr_0_0_y_xxxx_0 = pbuffer.data(idx_op_geom_010_gs + 15);

    auto tr_0_0_y_xxxy_0 = pbuffer.data(idx_op_geom_010_gs + 16);

    auto tr_0_0_y_xxxz_0 = pbuffer.data(idx_op_geom_010_gs + 17);

    auto tr_0_0_y_xxyy_0 = pbuffer.data(idx_op_geom_010_gs + 18);

    auto tr_0_0_y_xxyz_0 = pbuffer.data(idx_op_geom_010_gs + 19);

    auto tr_0_0_y_xxzz_0 = pbuffer.data(idx_op_geom_010_gs + 20);

    auto tr_0_0_y_xyyy_0 = pbuffer.data(idx_op_geom_010_gs + 21);

    auto tr_0_0_y_xyyz_0 = pbuffer.data(idx_op_geom_010_gs + 22);

    auto tr_0_0_y_xyzz_0 = pbuffer.data(idx_op_geom_010_gs + 23);

    auto tr_0_0_y_xzzz_0 = pbuffer.data(idx_op_geom_010_gs + 24);

    auto tr_0_0_y_yyyy_0 = pbuffer.data(idx_op_geom_010_gs + 25);

    auto tr_0_0_y_yyyz_0 = pbuffer.data(idx_op_geom_010_gs + 26);

    auto tr_0_0_y_yyzz_0 = pbuffer.data(idx_op_geom_010_gs + 27);

    auto tr_0_0_y_yzzz_0 = pbuffer.data(idx_op_geom_010_gs + 28);

    auto tr_0_0_y_zzzz_0 = pbuffer.data(idx_op_geom_010_gs + 29);

    auto tr_0_0_z_xxxx_0 = pbuffer.data(idx_op_geom_010_gs + 30);

    auto tr_0_0_z_xxxy_0 = pbuffer.data(idx_op_geom_010_gs + 31);

    auto tr_0_0_z_xxxz_0 = pbuffer.data(idx_op_geom_010_gs + 32);

    auto tr_0_0_z_xxyy_0 = pbuffer.data(idx_op_geom_010_gs + 33);

    auto tr_0_0_z_xxyz_0 = pbuffer.data(idx_op_geom_010_gs + 34);

    auto tr_0_0_z_xxzz_0 = pbuffer.data(idx_op_geom_010_gs + 35);

    auto tr_0_0_z_xyyy_0 = pbuffer.data(idx_op_geom_010_gs + 36);

    auto tr_0_0_z_xyyz_0 = pbuffer.data(idx_op_geom_010_gs + 37);

    auto tr_0_0_z_xyzz_0 = pbuffer.data(idx_op_geom_010_gs + 38);

    auto tr_0_0_z_xzzz_0 = pbuffer.data(idx_op_geom_010_gs + 39);

    auto tr_0_0_z_yyyy_0 = pbuffer.data(idx_op_geom_010_gs + 40);

    auto tr_0_0_z_yyyz_0 = pbuffer.data(idx_op_geom_010_gs + 41);

    auto tr_0_0_z_yyzz_0 = pbuffer.data(idx_op_geom_010_gs + 42);

    auto tr_0_0_z_yzzz_0 = pbuffer.data(idx_op_geom_010_gs + 43);

    auto tr_0_0_z_zzzz_0 = pbuffer.data(idx_op_geom_010_gs + 44);

    #pragma omp simd aligned(tr_0_0_x_xxxx_0, tr_0_0_x_xxxy_0, tr_0_0_x_xxxz_0, tr_0_0_x_xxyy_0, tr_0_0_x_xxyz_0, tr_0_0_x_xxzz_0, tr_0_0_x_xyyy_0, tr_0_0_x_xyyz_0, tr_0_0_x_xyzz_0, tr_0_0_x_xzzz_0, tr_0_0_x_yyyy_0, tr_0_0_x_yyyz_0, tr_0_0_x_yyzz_0, tr_0_0_x_yzzz_0, tr_0_0_x_zzzz_0, tr_0_0_y_xxxx_0, tr_0_0_y_xxxy_0, tr_0_0_y_xxxz_0, tr_0_0_y_xxyy_0, tr_0_0_y_xxyz_0, tr_0_0_y_xxzz_0, tr_0_0_y_xyyy_0, tr_0_0_y_xyyz_0, tr_0_0_y_xyzz_0, tr_0_0_y_xzzz_0, tr_0_0_y_yyyy_0, tr_0_0_y_yyyz_0, tr_0_0_y_yyzz_0, tr_0_0_y_yzzz_0, tr_0_0_y_zzzz_0, tr_0_0_z_xxxx_0, tr_0_0_z_xxxy_0, tr_0_0_z_xxxz_0, tr_0_0_z_xxyy_0, tr_0_0_z_xxyz_0, tr_0_0_z_xxzz_0, tr_0_0_z_xyyy_0, tr_0_0_z_xyyz_0, tr_0_0_z_xyzz_0, tr_0_0_z_xzzz_0, tr_0_0_z_yyyy_0, tr_0_0_z_yyyz_0, tr_0_0_z_yyzz_0, tr_0_0_z_yzzz_0, tr_0_0_z_zzzz_0, tr_xxx_0, tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, tr_xxxxx_0, tr_xxxxy_0, tr_xxxxz_0, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxxyy_0, tr_xxxyz_0, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxxzz_0, tr_xxy_0, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xxyyy_0, tr_xxyyz_0, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxyzz_0, tr_xxz_0, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xxzzz_0, tr_xyy_0, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_xyyyy_0, tr_xyyyz_0, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyyzz_0, tr_xyz_0, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xyzzz_0, tr_xzz_0, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_xzzzz_0, tr_yyy_0, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, tr_yyyyy_0, tr_yyyyz_0, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyyzz_0, tr_yyz_0, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yyzzz_0, tr_yzz_0, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_yzzzz_0, tr_zzz_0, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, tr_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxx_0[i] = 2.0 * tr_xxxxx_0[i] * tbe_0 + 2.0 * tr_xxxx_x[i] * tke_0 - 4.0 * tr_xxx_0[i];

        tr_0_0_x_xxxy_0[i] = 2.0 * tr_xxxxy_0[i] * tbe_0 + 2.0 * tr_xxxy_x[i] * tke_0 - 3.0 * tr_xxy_0[i];

        tr_0_0_x_xxxz_0[i] = 2.0 * tr_xxxxz_0[i] * tbe_0 + 2.0 * tr_xxxz_x[i] * tke_0 - 3.0 * tr_xxz_0[i];

        tr_0_0_x_xxyy_0[i] = 2.0 * tr_xxxyy_0[i] * tbe_0 + 2.0 * tr_xxyy_x[i] * tke_0 - 2.0 * tr_xyy_0[i];

        tr_0_0_x_xxyz_0[i] = 2.0 * tr_xxxyz_0[i] * tbe_0 + 2.0 * tr_xxyz_x[i] * tke_0 - 2.0 * tr_xyz_0[i];

        tr_0_0_x_xxzz_0[i] = 2.0 * tr_xxxzz_0[i] * tbe_0 + 2.0 * tr_xxzz_x[i] * tke_0 - 2.0 * tr_xzz_0[i];

        tr_0_0_x_xyyy_0[i] = 2.0 * tr_xxyyy_0[i] * tbe_0 + 2.0 * tr_xyyy_x[i] * tke_0 - tr_yyy_0[i];

        tr_0_0_x_xyyz_0[i] = 2.0 * tr_xxyyz_0[i] * tbe_0 + 2.0 * tr_xyyz_x[i] * tke_0 - tr_yyz_0[i];

        tr_0_0_x_xyzz_0[i] = 2.0 * tr_xxyzz_0[i] * tbe_0 + 2.0 * tr_xyzz_x[i] * tke_0 - tr_yzz_0[i];

        tr_0_0_x_xzzz_0[i] = 2.0 * tr_xxzzz_0[i] * tbe_0 + 2.0 * tr_xzzz_x[i] * tke_0 - tr_zzz_0[i];

        tr_0_0_x_yyyy_0[i] = 2.0 * tr_xyyyy_0[i] * tbe_0 + 2.0 * tr_yyyy_x[i] * tke_0;

        tr_0_0_x_yyyz_0[i] = 2.0 * tr_xyyyz_0[i] * tbe_0 + 2.0 * tr_yyyz_x[i] * tke_0;

        tr_0_0_x_yyzz_0[i] = 2.0 * tr_xyyzz_0[i] * tbe_0 + 2.0 * tr_yyzz_x[i] * tke_0;

        tr_0_0_x_yzzz_0[i] = 2.0 * tr_xyzzz_0[i] * tbe_0 + 2.0 * tr_yzzz_x[i] * tke_0;

        tr_0_0_x_zzzz_0[i] = 2.0 * tr_xzzzz_0[i] * tbe_0 + 2.0 * tr_zzzz_x[i] * tke_0;

        tr_0_0_y_xxxx_0[i] = 2.0 * tr_xxxxy_0[i] * tbe_0 + 2.0 * tr_xxxx_y[i] * tke_0;

        tr_0_0_y_xxxy_0[i] = 2.0 * tr_xxxyy_0[i] * tbe_0 + 2.0 * tr_xxxy_y[i] * tke_0 - tr_xxx_0[i];

        tr_0_0_y_xxxz_0[i] = 2.0 * tr_xxxyz_0[i] * tbe_0 + 2.0 * tr_xxxz_y[i] * tke_0;

        tr_0_0_y_xxyy_0[i] = 2.0 * tr_xxyyy_0[i] * tbe_0 + 2.0 * tr_xxyy_y[i] * tke_0 - 2.0 * tr_xxy_0[i];

        tr_0_0_y_xxyz_0[i] = 2.0 * tr_xxyyz_0[i] * tbe_0 + 2.0 * tr_xxyz_y[i] * tke_0 - tr_xxz_0[i];

        tr_0_0_y_xxzz_0[i] = 2.0 * tr_xxyzz_0[i] * tbe_0 + 2.0 * tr_xxzz_y[i] * tke_0;

        tr_0_0_y_xyyy_0[i] = 2.0 * tr_xyyyy_0[i] * tbe_0 + 2.0 * tr_xyyy_y[i] * tke_0 - 3.0 * tr_xyy_0[i];

        tr_0_0_y_xyyz_0[i] = 2.0 * tr_xyyyz_0[i] * tbe_0 + 2.0 * tr_xyyz_y[i] * tke_0 - 2.0 * tr_xyz_0[i];

        tr_0_0_y_xyzz_0[i] = 2.0 * tr_xyyzz_0[i] * tbe_0 + 2.0 * tr_xyzz_y[i] * tke_0 - tr_xzz_0[i];

        tr_0_0_y_xzzz_0[i] = 2.0 * tr_xyzzz_0[i] * tbe_0 + 2.0 * tr_xzzz_y[i] * tke_0;

        tr_0_0_y_yyyy_0[i] = 2.0 * tr_yyyyy_0[i] * tbe_0 + 2.0 * tr_yyyy_y[i] * tke_0 - 4.0 * tr_yyy_0[i];

        tr_0_0_y_yyyz_0[i] = 2.0 * tr_yyyyz_0[i] * tbe_0 + 2.0 * tr_yyyz_y[i] * tke_0 - 3.0 * tr_yyz_0[i];

        tr_0_0_y_yyzz_0[i] = 2.0 * tr_yyyzz_0[i] * tbe_0 + 2.0 * tr_yyzz_y[i] * tke_0 - 2.0 * tr_yzz_0[i];

        tr_0_0_y_yzzz_0[i] = 2.0 * tr_yyzzz_0[i] * tbe_0 + 2.0 * tr_yzzz_y[i] * tke_0 - tr_zzz_0[i];

        tr_0_0_y_zzzz_0[i] = 2.0 * tr_yzzzz_0[i] * tbe_0 + 2.0 * tr_zzzz_y[i] * tke_0;

        tr_0_0_z_xxxx_0[i] = 2.0 * tr_xxxxz_0[i] * tbe_0 + 2.0 * tr_xxxx_z[i] * tke_0;

        tr_0_0_z_xxxy_0[i] = 2.0 * tr_xxxyz_0[i] * tbe_0 + 2.0 * tr_xxxy_z[i] * tke_0;

        tr_0_0_z_xxxz_0[i] = 2.0 * tr_xxxzz_0[i] * tbe_0 + 2.0 * tr_xxxz_z[i] * tke_0 - tr_xxx_0[i];

        tr_0_0_z_xxyy_0[i] = 2.0 * tr_xxyyz_0[i] * tbe_0 + 2.0 * tr_xxyy_z[i] * tke_0;

        tr_0_0_z_xxyz_0[i] = 2.0 * tr_xxyzz_0[i] * tbe_0 + 2.0 * tr_xxyz_z[i] * tke_0 - tr_xxy_0[i];

        tr_0_0_z_xxzz_0[i] = 2.0 * tr_xxzzz_0[i] * tbe_0 + 2.0 * tr_xxzz_z[i] * tke_0 - 2.0 * tr_xxz_0[i];

        tr_0_0_z_xyyy_0[i] = 2.0 * tr_xyyyz_0[i] * tbe_0 + 2.0 * tr_xyyy_z[i] * tke_0;

        tr_0_0_z_xyyz_0[i] = 2.0 * tr_xyyzz_0[i] * tbe_0 + 2.0 * tr_xyyz_z[i] * tke_0 - tr_xyy_0[i];

        tr_0_0_z_xyzz_0[i] = 2.0 * tr_xyzzz_0[i] * tbe_0 + 2.0 * tr_xyzz_z[i] * tke_0 - 2.0 * tr_xyz_0[i];

        tr_0_0_z_xzzz_0[i] = 2.0 * tr_xzzzz_0[i] * tbe_0 + 2.0 * tr_xzzz_z[i] * tke_0 - 3.0 * tr_xzz_0[i];

        tr_0_0_z_yyyy_0[i] = 2.0 * tr_yyyyz_0[i] * tbe_0 + 2.0 * tr_yyyy_z[i] * tke_0;

        tr_0_0_z_yyyz_0[i] = 2.0 * tr_yyyzz_0[i] * tbe_0 + 2.0 * tr_yyyz_z[i] * tke_0 - tr_yyy_0[i];

        tr_0_0_z_yyzz_0[i] = 2.0 * tr_yyzzz_0[i] * tbe_0 + 2.0 * tr_yyzz_z[i] * tke_0 - 2.0 * tr_yyz_0[i];

        tr_0_0_z_yzzz_0[i] = 2.0 * tr_yzzzz_0[i] * tbe_0 + 2.0 * tr_yzzz_z[i] * tke_0 - 3.0 * tr_yzz_0[i];

        tr_0_0_z_zzzz_0[i] = 2.0 * tr_zzzzz_0[i] * tbe_0 + 2.0 * tr_zzzz_z[i] * tke_0 - 4.0 * tr_zzz_0[i];
    }
}

} // t2cgeom namespace

