#include "TwoCenterElectronRepulsionPrimRecPG.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_pg(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_pg,
                                const size_t idx_eri_1_sf,
                                const size_t idx_eri_1_sg,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SF

    auto g_0_xxx_1 = pbuffer.data(idx_eri_1_sf);

    auto g_0_xxy_1 = pbuffer.data(idx_eri_1_sf + 1);

    auto g_0_xxz_1 = pbuffer.data(idx_eri_1_sf + 2);

    auto g_0_xyy_1 = pbuffer.data(idx_eri_1_sf + 3);

    auto g_0_xyz_1 = pbuffer.data(idx_eri_1_sf + 4);

    auto g_0_xzz_1 = pbuffer.data(idx_eri_1_sf + 5);

    auto g_0_yyy_1 = pbuffer.data(idx_eri_1_sf + 6);

    auto g_0_yyz_1 = pbuffer.data(idx_eri_1_sf + 7);

    auto g_0_yzz_1 = pbuffer.data(idx_eri_1_sf + 8);

    auto g_0_zzz_1 = pbuffer.data(idx_eri_1_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto g_0_xxxx_1 = pbuffer.data(idx_eri_1_sg);

    auto g_0_xxxy_1 = pbuffer.data(idx_eri_1_sg + 1);

    auto g_0_xxxz_1 = pbuffer.data(idx_eri_1_sg + 2);

    auto g_0_xxyy_1 = pbuffer.data(idx_eri_1_sg + 3);

    auto g_0_xxyz_1 = pbuffer.data(idx_eri_1_sg + 4);

    auto g_0_xxzz_1 = pbuffer.data(idx_eri_1_sg + 5);

    auto g_0_xyyy_1 = pbuffer.data(idx_eri_1_sg + 6);

    auto g_0_xyyz_1 = pbuffer.data(idx_eri_1_sg + 7);

    auto g_0_xyzz_1 = pbuffer.data(idx_eri_1_sg + 8);

    auto g_0_xzzz_1 = pbuffer.data(idx_eri_1_sg + 9);

    auto g_0_yyyy_1 = pbuffer.data(idx_eri_1_sg + 10);

    auto g_0_yyyz_1 = pbuffer.data(idx_eri_1_sg + 11);

    auto g_0_yyzz_1 = pbuffer.data(idx_eri_1_sg + 12);

    auto g_0_yzzz_1 = pbuffer.data(idx_eri_1_sg + 13);

    auto g_0_zzzz_1 = pbuffer.data(idx_eri_1_sg + 14);

    // Set up 0-15 components of targeted buffer : PG

    auto g_x_xxxx_0 = pbuffer.data(idx_eri_0_pg);

    auto g_x_xxxy_0 = pbuffer.data(idx_eri_0_pg + 1);

    auto g_x_xxxz_0 = pbuffer.data(idx_eri_0_pg + 2);

    auto g_x_xxyy_0 = pbuffer.data(idx_eri_0_pg + 3);

    auto g_x_xxyz_0 = pbuffer.data(idx_eri_0_pg + 4);

    auto g_x_xxzz_0 = pbuffer.data(idx_eri_0_pg + 5);

    auto g_x_xyyy_0 = pbuffer.data(idx_eri_0_pg + 6);

    auto g_x_xyyz_0 = pbuffer.data(idx_eri_0_pg + 7);

    auto g_x_xyzz_0 = pbuffer.data(idx_eri_0_pg + 8);

    auto g_x_xzzz_0 = pbuffer.data(idx_eri_0_pg + 9);

    auto g_x_yyyy_0 = pbuffer.data(idx_eri_0_pg + 10);

    auto g_x_yyyz_0 = pbuffer.data(idx_eri_0_pg + 11);

    auto g_x_yyzz_0 = pbuffer.data(idx_eri_0_pg + 12);

    auto g_x_yzzz_0 = pbuffer.data(idx_eri_0_pg + 13);

    auto g_x_zzzz_0 = pbuffer.data(idx_eri_0_pg + 14);

    #pragma omp simd aligned(g_0_xxx_1, g_0_xxxx_1, g_0_xxxy_1, g_0_xxxz_1, g_0_xxy_1, g_0_xxyy_1, g_0_xxyz_1, g_0_xxz_1, g_0_xxzz_1, g_0_xyy_1, g_0_xyyy_1, g_0_xyyz_1, g_0_xyz_1, g_0_xyzz_1, g_0_xzz_1, g_0_xzzz_1, g_0_yyy_1, g_0_yyyy_1, g_0_yyyz_1, g_0_yyz_1, g_0_yyzz_1, g_0_yzz_1, g_0_yzzz_1, g_0_zzz_1, g_0_zzzz_1, g_x_xxxx_0, g_x_xxxy_0, g_x_xxxz_0, g_x_xxyy_0, g_x_xxyz_0, g_x_xxzz_0, g_x_xyyy_0, g_x_xyyz_0, g_x_xyzz_0, g_x_xzzz_0, g_x_yyyy_0, g_x_yyyz_0, g_x_yyzz_0, g_x_yzzz_0, g_x_zzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_x_xxxx_0[i] = 4.0 * g_0_xxx_1[i] * fe_0 + g_0_xxxx_1[i] * pa_x[i];

        g_x_xxxy_0[i] = 3.0 * g_0_xxy_1[i] * fe_0 + g_0_xxxy_1[i] * pa_x[i];

        g_x_xxxz_0[i] = 3.0 * g_0_xxz_1[i] * fe_0 + g_0_xxxz_1[i] * pa_x[i];

        g_x_xxyy_0[i] = 2.0 * g_0_xyy_1[i] * fe_0 + g_0_xxyy_1[i] * pa_x[i];

        g_x_xxyz_0[i] = 2.0 * g_0_xyz_1[i] * fe_0 + g_0_xxyz_1[i] * pa_x[i];

        g_x_xxzz_0[i] = 2.0 * g_0_xzz_1[i] * fe_0 + g_0_xxzz_1[i] * pa_x[i];

        g_x_xyyy_0[i] = g_0_yyy_1[i] * fe_0 + g_0_xyyy_1[i] * pa_x[i];

        g_x_xyyz_0[i] = g_0_yyz_1[i] * fe_0 + g_0_xyyz_1[i] * pa_x[i];

        g_x_xyzz_0[i] = g_0_yzz_1[i] * fe_0 + g_0_xyzz_1[i] * pa_x[i];

        g_x_xzzz_0[i] = g_0_zzz_1[i] * fe_0 + g_0_xzzz_1[i] * pa_x[i];

        g_x_yyyy_0[i] = g_0_yyyy_1[i] * pa_x[i];

        g_x_yyyz_0[i] = g_0_yyyz_1[i] * pa_x[i];

        g_x_yyzz_0[i] = g_0_yyzz_1[i] * pa_x[i];

        g_x_yzzz_0[i] = g_0_yzzz_1[i] * pa_x[i];

        g_x_zzzz_0[i] = g_0_zzzz_1[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : PG

    auto g_y_xxxx_0 = pbuffer.data(idx_eri_0_pg + 15);

    auto g_y_xxxy_0 = pbuffer.data(idx_eri_0_pg + 16);

    auto g_y_xxxz_0 = pbuffer.data(idx_eri_0_pg + 17);

    auto g_y_xxyy_0 = pbuffer.data(idx_eri_0_pg + 18);

    auto g_y_xxyz_0 = pbuffer.data(idx_eri_0_pg + 19);

    auto g_y_xxzz_0 = pbuffer.data(idx_eri_0_pg + 20);

    auto g_y_xyyy_0 = pbuffer.data(idx_eri_0_pg + 21);

    auto g_y_xyyz_0 = pbuffer.data(idx_eri_0_pg + 22);

    auto g_y_xyzz_0 = pbuffer.data(idx_eri_0_pg + 23);

    auto g_y_xzzz_0 = pbuffer.data(idx_eri_0_pg + 24);

    auto g_y_yyyy_0 = pbuffer.data(idx_eri_0_pg + 25);

    auto g_y_yyyz_0 = pbuffer.data(idx_eri_0_pg + 26);

    auto g_y_yyzz_0 = pbuffer.data(idx_eri_0_pg + 27);

    auto g_y_yzzz_0 = pbuffer.data(idx_eri_0_pg + 28);

    auto g_y_zzzz_0 = pbuffer.data(idx_eri_0_pg + 29);

    #pragma omp simd aligned(g_0_xxx_1, g_0_xxxx_1, g_0_xxxy_1, g_0_xxxz_1, g_0_xxy_1, g_0_xxyy_1, g_0_xxyz_1, g_0_xxz_1, g_0_xxzz_1, g_0_xyy_1, g_0_xyyy_1, g_0_xyyz_1, g_0_xyz_1, g_0_xyzz_1, g_0_xzz_1, g_0_xzzz_1, g_0_yyy_1, g_0_yyyy_1, g_0_yyyz_1, g_0_yyz_1, g_0_yyzz_1, g_0_yzz_1, g_0_yzzz_1, g_0_zzz_1, g_0_zzzz_1, g_y_xxxx_0, g_y_xxxy_0, g_y_xxxz_0, g_y_xxyy_0, g_y_xxyz_0, g_y_xxzz_0, g_y_xyyy_0, g_y_xyyz_0, g_y_xyzz_0, g_y_xzzz_0, g_y_yyyy_0, g_y_yyyz_0, g_y_yyzz_0, g_y_yzzz_0, g_y_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_y_xxxx_0[i] = g_0_xxxx_1[i] * pa_y[i];

        g_y_xxxy_0[i] = g_0_xxx_1[i] * fe_0 + g_0_xxxy_1[i] * pa_y[i];

        g_y_xxxz_0[i] = g_0_xxxz_1[i] * pa_y[i];

        g_y_xxyy_0[i] = 2.0 * g_0_xxy_1[i] * fe_0 + g_0_xxyy_1[i] * pa_y[i];

        g_y_xxyz_0[i] = g_0_xxz_1[i] * fe_0 + g_0_xxyz_1[i] * pa_y[i];

        g_y_xxzz_0[i] = g_0_xxzz_1[i] * pa_y[i];

        g_y_xyyy_0[i] = 3.0 * g_0_xyy_1[i] * fe_0 + g_0_xyyy_1[i] * pa_y[i];

        g_y_xyyz_0[i] = 2.0 * g_0_xyz_1[i] * fe_0 + g_0_xyyz_1[i] * pa_y[i];

        g_y_xyzz_0[i] = g_0_xzz_1[i] * fe_0 + g_0_xyzz_1[i] * pa_y[i];

        g_y_xzzz_0[i] = g_0_xzzz_1[i] * pa_y[i];

        g_y_yyyy_0[i] = 4.0 * g_0_yyy_1[i] * fe_0 + g_0_yyyy_1[i] * pa_y[i];

        g_y_yyyz_0[i] = 3.0 * g_0_yyz_1[i] * fe_0 + g_0_yyyz_1[i] * pa_y[i];

        g_y_yyzz_0[i] = 2.0 * g_0_yzz_1[i] * fe_0 + g_0_yyzz_1[i] * pa_y[i];

        g_y_yzzz_0[i] = g_0_zzz_1[i] * fe_0 + g_0_yzzz_1[i] * pa_y[i];

        g_y_zzzz_0[i] = g_0_zzzz_1[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : PG

    auto g_z_xxxx_0 = pbuffer.data(idx_eri_0_pg + 30);

    auto g_z_xxxy_0 = pbuffer.data(idx_eri_0_pg + 31);

    auto g_z_xxxz_0 = pbuffer.data(idx_eri_0_pg + 32);

    auto g_z_xxyy_0 = pbuffer.data(idx_eri_0_pg + 33);

    auto g_z_xxyz_0 = pbuffer.data(idx_eri_0_pg + 34);

    auto g_z_xxzz_0 = pbuffer.data(idx_eri_0_pg + 35);

    auto g_z_xyyy_0 = pbuffer.data(idx_eri_0_pg + 36);

    auto g_z_xyyz_0 = pbuffer.data(idx_eri_0_pg + 37);

    auto g_z_xyzz_0 = pbuffer.data(idx_eri_0_pg + 38);

    auto g_z_xzzz_0 = pbuffer.data(idx_eri_0_pg + 39);

    auto g_z_yyyy_0 = pbuffer.data(idx_eri_0_pg + 40);

    auto g_z_yyyz_0 = pbuffer.data(idx_eri_0_pg + 41);

    auto g_z_yyzz_0 = pbuffer.data(idx_eri_0_pg + 42);

    auto g_z_yzzz_0 = pbuffer.data(idx_eri_0_pg + 43);

    auto g_z_zzzz_0 = pbuffer.data(idx_eri_0_pg + 44);

    #pragma omp simd aligned(g_0_xxx_1, g_0_xxxx_1, g_0_xxxy_1, g_0_xxxz_1, g_0_xxy_1, g_0_xxyy_1, g_0_xxyz_1, g_0_xxz_1, g_0_xxzz_1, g_0_xyy_1, g_0_xyyy_1, g_0_xyyz_1, g_0_xyz_1, g_0_xyzz_1, g_0_xzz_1, g_0_xzzz_1, g_0_yyy_1, g_0_yyyy_1, g_0_yyyz_1, g_0_yyz_1, g_0_yyzz_1, g_0_yzz_1, g_0_yzzz_1, g_0_zzz_1, g_0_zzzz_1, g_z_xxxx_0, g_z_xxxy_0, g_z_xxxz_0, g_z_xxyy_0, g_z_xxyz_0, g_z_xxzz_0, g_z_xyyy_0, g_z_xyyz_0, g_z_xyzz_0, g_z_xzzz_0, g_z_yyyy_0, g_z_yyyz_0, g_z_yyzz_0, g_z_yzzz_0, g_z_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_z_xxxx_0[i] = g_0_xxxx_1[i] * pa_z[i];

        g_z_xxxy_0[i] = g_0_xxxy_1[i] * pa_z[i];

        g_z_xxxz_0[i] = g_0_xxx_1[i] * fe_0 + g_0_xxxz_1[i] * pa_z[i];

        g_z_xxyy_0[i] = g_0_xxyy_1[i] * pa_z[i];

        g_z_xxyz_0[i] = g_0_xxy_1[i] * fe_0 + g_0_xxyz_1[i] * pa_z[i];

        g_z_xxzz_0[i] = 2.0 * g_0_xxz_1[i] * fe_0 + g_0_xxzz_1[i] * pa_z[i];

        g_z_xyyy_0[i] = g_0_xyyy_1[i] * pa_z[i];

        g_z_xyyz_0[i] = g_0_xyy_1[i] * fe_0 + g_0_xyyz_1[i] * pa_z[i];

        g_z_xyzz_0[i] = 2.0 * g_0_xyz_1[i] * fe_0 + g_0_xyzz_1[i] * pa_z[i];

        g_z_xzzz_0[i] = 3.0 * g_0_xzz_1[i] * fe_0 + g_0_xzzz_1[i] * pa_z[i];

        g_z_yyyy_0[i] = g_0_yyyy_1[i] * pa_z[i];

        g_z_yyyz_0[i] = g_0_yyy_1[i] * fe_0 + g_0_yyyz_1[i] * pa_z[i];

        g_z_yyzz_0[i] = 2.0 * g_0_yyz_1[i] * fe_0 + g_0_yyzz_1[i] * pa_z[i];

        g_z_yzzz_0[i] = 3.0 * g_0_yzz_1[i] * fe_0 + g_0_yzzz_1[i] * pa_z[i];

        g_z_zzzz_0[i] = 4.0 * g_0_zzz_1[i] * fe_0 + g_0_zzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

