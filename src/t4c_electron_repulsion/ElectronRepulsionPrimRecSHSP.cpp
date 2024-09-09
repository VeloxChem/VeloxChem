#include "ElectronRepulsionPrimRecSHSP.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_shsp(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_shsp,
                                  size_t                idx_eri_0_sfsp,
                                  size_t                idx_eri_1_sfsp,
                                  size_t                idx_eri_1_sgss,
                                  size_t                idx_eri_0_sgsp,
                                  size_t                idx_eri_1_sgsp,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SFSP

    auto g_0_xxx_0_x_0 = pbuffer.data(idx_eri_0_sfsp);

    auto g_0_xxx_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 1);

    auto g_0_xxx_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 2);

    auto g_0_xxy_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 3);

    auto g_0_xxz_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 6);

    auto g_0_xyy_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 10);

    auto g_0_xyy_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 11);

    auto g_0_xzz_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 16);

    auto g_0_xzz_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 17);

    auto g_0_yyy_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 18);

    auto g_0_yyy_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 19);

    auto g_0_yyy_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 20);

    auto g_0_yyz_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 22);

    auto g_0_yzz_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 24);

    auto g_0_yzz_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 26);

    auto g_0_zzz_0_x_0 = pbuffer.data(idx_eri_0_sfsp + 27);

    auto g_0_zzz_0_y_0 = pbuffer.data(idx_eri_0_sfsp + 28);

    auto g_0_zzz_0_z_0 = pbuffer.data(idx_eri_0_sfsp + 29);

    /// Set up components of auxilary buffer : SFSP

    auto g_0_xxx_0_x_1 = pbuffer.data(idx_eri_1_sfsp);

    auto g_0_xxx_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 1);

    auto g_0_xxx_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 2);

    auto g_0_xxy_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 3);

    auto g_0_xxz_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 6);

    auto g_0_xyy_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 10);

    auto g_0_xyy_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 11);

    auto g_0_xzz_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 16);

    auto g_0_xzz_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 17);

    auto g_0_yyy_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 18);

    auto g_0_yyy_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 19);

    auto g_0_yyy_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 20);

    auto g_0_yyz_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 22);

    auto g_0_yzz_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 24);

    auto g_0_yzz_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 26);

    auto g_0_zzz_0_x_1 = pbuffer.data(idx_eri_1_sfsp + 27);

    auto g_0_zzz_0_y_1 = pbuffer.data(idx_eri_1_sfsp + 28);

    auto g_0_zzz_0_z_1 = pbuffer.data(idx_eri_1_sfsp + 29);

    /// Set up components of auxilary buffer : SGSS

    auto g_0_xxxx_0_0_1 = pbuffer.data(idx_eri_1_sgss);

    auto g_0_xxyy_0_0_1 = pbuffer.data(idx_eri_1_sgss + 3);

    auto g_0_xxzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 5);

    auto g_0_yyyy_0_0_1 = pbuffer.data(idx_eri_1_sgss + 10);

    auto g_0_yyzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 12);

    auto g_0_zzzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 14);

    /// Set up components of auxilary buffer : SGSP

    auto g_0_xxxx_0_x_0 = pbuffer.data(idx_eri_0_sgsp);

    auto g_0_xxxx_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 1);

    auto g_0_xxxx_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 2);

    auto g_0_xxxy_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 3);

    auto g_0_xxxy_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 4);

    auto g_0_xxxz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 6);

    auto g_0_xxxz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 8);

    auto g_0_xxyy_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 9);

    auto g_0_xxyy_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 10);

    auto g_0_xxyy_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 11);

    auto g_0_xxzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 15);

    auto g_0_xxzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 16);

    auto g_0_xxzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 17);

    auto g_0_xyyy_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 18);

    auto g_0_xyyy_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 19);

    auto g_0_xyyy_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 20);

    auto g_0_xzzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 27);

    auto g_0_xzzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 28);

    auto g_0_xzzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 29);

    auto g_0_yyyy_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 30);

    auto g_0_yyyy_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 31);

    auto g_0_yyyy_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 32);

    auto g_0_yyyz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 34);

    auto g_0_yyyz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 35);

    auto g_0_yyzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 36);

    auto g_0_yyzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 37);

    auto g_0_yyzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 38);

    auto g_0_yzzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 39);

    auto g_0_yzzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 40);

    auto g_0_yzzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 41);

    auto g_0_zzzz_0_x_0 = pbuffer.data(idx_eri_0_sgsp + 42);

    auto g_0_zzzz_0_y_0 = pbuffer.data(idx_eri_0_sgsp + 43);

    auto g_0_zzzz_0_z_0 = pbuffer.data(idx_eri_0_sgsp + 44);

    /// Set up components of auxilary buffer : SGSP

    auto g_0_xxxx_0_x_1 = pbuffer.data(idx_eri_1_sgsp);

    auto g_0_xxxx_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 1);

    auto g_0_xxxx_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 2);

    auto g_0_xxxy_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 3);

    auto g_0_xxxy_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 4);

    auto g_0_xxxz_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 6);

    auto g_0_xxxz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 8);

    auto g_0_xxyy_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 9);

    auto g_0_xxyy_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 10);

    auto g_0_xxyy_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 11);

    auto g_0_xxzz_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 15);

    auto g_0_xxzz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 16);

    auto g_0_xxzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 17);

    auto g_0_xyyy_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 18);

    auto g_0_xyyy_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 19);

    auto g_0_xyyy_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 20);

    auto g_0_xzzz_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 27);

    auto g_0_xzzz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 28);

    auto g_0_xzzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 29);

    auto g_0_yyyy_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 30);

    auto g_0_yyyy_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 31);

    auto g_0_yyyy_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 32);

    auto g_0_yyyz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 34);

    auto g_0_yyyz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 35);

    auto g_0_yyzz_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 36);

    auto g_0_yyzz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 37);

    auto g_0_yyzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 38);

    auto g_0_yzzz_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 39);

    auto g_0_yzzz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 40);

    auto g_0_yzzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 41);

    auto g_0_zzzz_0_x_1 = pbuffer.data(idx_eri_1_sgsp + 42);

    auto g_0_zzzz_0_y_1 = pbuffer.data(idx_eri_1_sgsp + 43);

    auto g_0_zzzz_0_z_1 = pbuffer.data(idx_eri_1_sgsp + 44);

    /// Set up 0-3 components of targeted buffer : SHSP

    auto g_0_xxxxx_0_x_0 = pbuffer.data(idx_eri_0_shsp);

    auto g_0_xxxxx_0_y_0 = pbuffer.data(idx_eri_0_shsp + 1);

    auto g_0_xxxxx_0_z_0 = pbuffer.data(idx_eri_0_shsp + 2);

#pragma omp simd aligned(g_0_xxx_0_x_0,       \
                             g_0_xxx_0_x_1,   \
                             g_0_xxx_0_y_0,   \
                             g_0_xxx_0_y_1,   \
                             g_0_xxx_0_z_0,   \
                             g_0_xxx_0_z_1,   \
                             g_0_xxxx_0_0_1,  \
                             g_0_xxxx_0_x_0,  \
                             g_0_xxxx_0_x_1,  \
                             g_0_xxxx_0_y_0,  \
                             g_0_xxxx_0_y_1,  \
                             g_0_xxxx_0_z_0,  \
                             g_0_xxxx_0_z_1,  \
                             g_0_xxxxx_0_x_0, \
                             g_0_xxxxx_0_y_0, \
                             g_0_xxxxx_0_z_0, \
                             wp_x,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxx_0_x_0[i] = 4.0 * g_0_xxx_0_x_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_x_1[i] * fti_ab_0 + g_0_xxxx_0_0_1[i] * fi_abcd_0 +
                             g_0_xxxx_0_x_0[i] * pb_x + g_0_xxxx_0_x_1[i] * wp_x[i];

        g_0_xxxxx_0_y_0[i] =
            4.0 * g_0_xxx_0_y_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_y_1[i] * fti_ab_0 + g_0_xxxx_0_y_0[i] * pb_x + g_0_xxxx_0_y_1[i] * wp_x[i];

        g_0_xxxxx_0_z_0[i] =
            4.0 * g_0_xxx_0_z_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_z_1[i] * fti_ab_0 + g_0_xxxx_0_z_0[i] * pb_x + g_0_xxxx_0_z_1[i] * wp_x[i];
    }

    /// Set up 3-6 components of targeted buffer : SHSP

    auto g_0_xxxxy_0_x_0 = pbuffer.data(idx_eri_0_shsp + 3);

    auto g_0_xxxxy_0_y_0 = pbuffer.data(idx_eri_0_shsp + 4);

    auto g_0_xxxxy_0_z_0 = pbuffer.data(idx_eri_0_shsp + 5);

#pragma omp simd aligned(g_0_xxxx_0_0_1,      \
                             g_0_xxxx_0_x_0,  \
                             g_0_xxxx_0_x_1,  \
                             g_0_xxxx_0_y_0,  \
                             g_0_xxxx_0_y_1,  \
                             g_0_xxxx_0_z_0,  \
                             g_0_xxxx_0_z_1,  \
                             g_0_xxxxy_0_x_0, \
                             g_0_xxxxy_0_y_0, \
                             g_0_xxxxy_0_z_0, \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxy_0_x_0[i] = g_0_xxxx_0_x_0[i] * pb_y + g_0_xxxx_0_x_1[i] * wp_y[i];

        g_0_xxxxy_0_y_0[i] = g_0_xxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxx_0_y_0[i] * pb_y + g_0_xxxx_0_y_1[i] * wp_y[i];

        g_0_xxxxy_0_z_0[i] = g_0_xxxx_0_z_0[i] * pb_y + g_0_xxxx_0_z_1[i] * wp_y[i];
    }

    /// Set up 6-9 components of targeted buffer : SHSP

    auto g_0_xxxxz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 6);

    auto g_0_xxxxz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 7);

    auto g_0_xxxxz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 8);

#pragma omp simd aligned(g_0_xxxx_0_0_1,      \
                             g_0_xxxx_0_x_0,  \
                             g_0_xxxx_0_x_1,  \
                             g_0_xxxx_0_y_0,  \
                             g_0_xxxx_0_y_1,  \
                             g_0_xxxx_0_z_0,  \
                             g_0_xxxx_0_z_1,  \
                             g_0_xxxxz_0_x_0, \
                             g_0_xxxxz_0_y_0, \
                             g_0_xxxxz_0_z_0, \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxz_0_x_0[i] = g_0_xxxx_0_x_0[i] * pb_z + g_0_xxxx_0_x_1[i] * wp_z[i];

        g_0_xxxxz_0_y_0[i] = g_0_xxxx_0_y_0[i] * pb_z + g_0_xxxx_0_y_1[i] * wp_z[i];

        g_0_xxxxz_0_z_0[i] = g_0_xxxx_0_0_1[i] * fi_abcd_0 + g_0_xxxx_0_z_0[i] * pb_z + g_0_xxxx_0_z_1[i] * wp_z[i];
    }

    /// Set up 9-12 components of targeted buffer : SHSP

    auto g_0_xxxyy_0_x_0 = pbuffer.data(idx_eri_0_shsp + 9);

    auto g_0_xxxyy_0_y_0 = pbuffer.data(idx_eri_0_shsp + 10);

    auto g_0_xxxyy_0_z_0 = pbuffer.data(idx_eri_0_shsp + 11);

#pragma omp simd aligned(g_0_xxx_0_x_0,       \
                             g_0_xxx_0_x_1,   \
                             g_0_xxxy_0_x_0,  \
                             g_0_xxxy_0_x_1,  \
                             g_0_xxxyy_0_x_0, \
                             g_0_xxxyy_0_y_0, \
                             g_0_xxxyy_0_z_0, \
                             g_0_xxyy_0_y_0,  \
                             g_0_xxyy_0_y_1,  \
                             g_0_xxyy_0_z_0,  \
                             g_0_xxyy_0_z_1,  \
                             g_0_xyy_0_y_0,   \
                             g_0_xyy_0_y_1,   \
                             g_0_xyy_0_z_0,   \
                             g_0_xyy_0_z_1,   \
                             wp_x,            \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyy_0_x_0[i] = g_0_xxx_0_x_0[i] * fi_ab_0 - g_0_xxx_0_x_1[i] * fti_ab_0 + g_0_xxxy_0_x_0[i] * pb_y + g_0_xxxy_0_x_1[i] * wp_y[i];

        g_0_xxxyy_0_y_0[i] =
            2.0 * g_0_xyy_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_y_1[i] * fti_ab_0 + g_0_xxyy_0_y_0[i] * pb_x + g_0_xxyy_0_y_1[i] * wp_x[i];

        g_0_xxxyy_0_z_0[i] =
            2.0 * g_0_xyy_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_z_1[i] * fti_ab_0 + g_0_xxyy_0_z_0[i] * pb_x + g_0_xxyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 12-15 components of targeted buffer : SHSP

    auto g_0_xxxyz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 12);

    auto g_0_xxxyz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 13);

    auto g_0_xxxyz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 14);

#pragma omp simd aligned(g_0_xxxy_0_y_0,      \
                             g_0_xxxy_0_y_1,  \
                             g_0_xxxyz_0_x_0, \
                             g_0_xxxyz_0_y_0, \
                             g_0_xxxyz_0_z_0, \
                             g_0_xxxz_0_x_0,  \
                             g_0_xxxz_0_x_1,  \
                             g_0_xxxz_0_z_0,  \
                             g_0_xxxz_0_z_1,  \
                             wp_y,            \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xxxyz_0_x_0[i] = g_0_xxxz_0_x_0[i] * pb_y + g_0_xxxz_0_x_1[i] * wp_y[i];

        g_0_xxxyz_0_y_0[i] = g_0_xxxy_0_y_0[i] * pb_z + g_0_xxxy_0_y_1[i] * wp_z[i];

        g_0_xxxyz_0_z_0[i] = g_0_xxxz_0_z_0[i] * pb_y + g_0_xxxz_0_z_1[i] * wp_y[i];
    }

    /// Set up 15-18 components of targeted buffer : SHSP

    auto g_0_xxxzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 15);

    auto g_0_xxxzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 16);

    auto g_0_xxxzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 17);

#pragma omp simd aligned(g_0_xxx_0_x_0,       \
                             g_0_xxx_0_x_1,   \
                             g_0_xxxz_0_x_0,  \
                             g_0_xxxz_0_x_1,  \
                             g_0_xxxzz_0_x_0, \
                             g_0_xxxzz_0_y_0, \
                             g_0_xxxzz_0_z_0, \
                             g_0_xxzz_0_y_0,  \
                             g_0_xxzz_0_y_1,  \
                             g_0_xxzz_0_z_0,  \
                             g_0_xxzz_0_z_1,  \
                             g_0_xzz_0_y_0,   \
                             g_0_xzz_0_y_1,   \
                             g_0_xzz_0_z_0,   \
                             g_0_xzz_0_z_1,   \
                             wp_x,            \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxzz_0_x_0[i] = g_0_xxx_0_x_0[i] * fi_ab_0 - g_0_xxx_0_x_1[i] * fti_ab_0 + g_0_xxxz_0_x_0[i] * pb_z + g_0_xxxz_0_x_1[i] * wp_z[i];

        g_0_xxxzz_0_y_0[i] =
            2.0 * g_0_xzz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_y_1[i] * fti_ab_0 + g_0_xxzz_0_y_0[i] * pb_x + g_0_xxzz_0_y_1[i] * wp_x[i];

        g_0_xxxzz_0_z_0[i] =
            2.0 * g_0_xzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_z_1[i] * fti_ab_0 + g_0_xxzz_0_z_0[i] * pb_x + g_0_xxzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 18-21 components of targeted buffer : SHSP

    auto g_0_xxyyy_0_x_0 = pbuffer.data(idx_eri_0_shsp + 18);

    auto g_0_xxyyy_0_y_0 = pbuffer.data(idx_eri_0_shsp + 19);

    auto g_0_xxyyy_0_z_0 = pbuffer.data(idx_eri_0_shsp + 20);

#pragma omp simd aligned(g_0_xxy_0_x_0,       \
                             g_0_xxy_0_x_1,   \
                             g_0_xxyy_0_x_0,  \
                             g_0_xxyy_0_x_1,  \
                             g_0_xxyyy_0_x_0, \
                             g_0_xxyyy_0_y_0, \
                             g_0_xxyyy_0_z_0, \
                             g_0_xyyy_0_y_0,  \
                             g_0_xyyy_0_y_1,  \
                             g_0_xyyy_0_z_0,  \
                             g_0_xyyy_0_z_1,  \
                             g_0_yyy_0_y_0,   \
                             g_0_yyy_0_y_1,   \
                             g_0_yyy_0_z_0,   \
                             g_0_yyy_0_z_1,   \
                             wp_x,            \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyy_0_x_0[i] =
            2.0 * g_0_xxy_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_x_1[i] * fti_ab_0 + g_0_xxyy_0_x_0[i] * pb_y + g_0_xxyy_0_x_1[i] * wp_y[i];

        g_0_xxyyy_0_y_0[i] = g_0_yyy_0_y_0[i] * fi_ab_0 - g_0_yyy_0_y_1[i] * fti_ab_0 + g_0_xyyy_0_y_0[i] * pb_x + g_0_xyyy_0_y_1[i] * wp_x[i];

        g_0_xxyyy_0_z_0[i] = g_0_yyy_0_z_0[i] * fi_ab_0 - g_0_yyy_0_z_1[i] * fti_ab_0 + g_0_xyyy_0_z_0[i] * pb_x + g_0_xyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 21-24 components of targeted buffer : SHSP

    auto g_0_xxyyz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 21);

    auto g_0_xxyyz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 22);

    auto g_0_xxyyz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 23);

#pragma omp simd aligned(g_0_xxyy_0_0_1,      \
                             g_0_xxyy_0_x_0,  \
                             g_0_xxyy_0_x_1,  \
                             g_0_xxyy_0_y_0,  \
                             g_0_xxyy_0_y_1,  \
                             g_0_xxyy_0_z_0,  \
                             g_0_xxyy_0_z_1,  \
                             g_0_xxyyz_0_x_0, \
                             g_0_xxyyz_0_y_0, \
                             g_0_xxyyz_0_z_0, \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyz_0_x_0[i] = g_0_xxyy_0_x_0[i] * pb_z + g_0_xxyy_0_x_1[i] * wp_z[i];

        g_0_xxyyz_0_y_0[i] = g_0_xxyy_0_y_0[i] * pb_z + g_0_xxyy_0_y_1[i] * wp_z[i];

        g_0_xxyyz_0_z_0[i] = g_0_xxyy_0_0_1[i] * fi_abcd_0 + g_0_xxyy_0_z_0[i] * pb_z + g_0_xxyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 24-27 components of targeted buffer : SHSP

    auto g_0_xxyzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 24);

    auto g_0_xxyzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 25);

    auto g_0_xxyzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 26);

#pragma omp simd aligned(g_0_xxyzz_0_x_0,     \
                             g_0_xxyzz_0_y_0, \
                             g_0_xxyzz_0_z_0, \
                             g_0_xxzz_0_0_1,  \
                             g_0_xxzz_0_x_0,  \
                             g_0_xxzz_0_x_1,  \
                             g_0_xxzz_0_y_0,  \
                             g_0_xxzz_0_y_1,  \
                             g_0_xxzz_0_z_0,  \
                             g_0_xxzz_0_z_1,  \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzz_0_x_0[i] = g_0_xxzz_0_x_0[i] * pb_y + g_0_xxzz_0_x_1[i] * wp_y[i];

        g_0_xxyzz_0_y_0[i] = g_0_xxzz_0_0_1[i] * fi_abcd_0 + g_0_xxzz_0_y_0[i] * pb_y + g_0_xxzz_0_y_1[i] * wp_y[i];

        g_0_xxyzz_0_z_0[i] = g_0_xxzz_0_z_0[i] * pb_y + g_0_xxzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 27-30 components of targeted buffer : SHSP

    auto g_0_xxzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 27);

    auto g_0_xxzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 28);

    auto g_0_xxzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 29);

#pragma omp simd aligned(g_0_xxz_0_x_0,       \
                             g_0_xxz_0_x_1,   \
                             g_0_xxzz_0_x_0,  \
                             g_0_xxzz_0_x_1,  \
                             g_0_xxzzz_0_x_0, \
                             g_0_xxzzz_0_y_0, \
                             g_0_xxzzz_0_z_0, \
                             g_0_xzzz_0_y_0,  \
                             g_0_xzzz_0_y_1,  \
                             g_0_xzzz_0_z_0,  \
                             g_0_xzzz_0_z_1,  \
                             g_0_zzz_0_y_0,   \
                             g_0_zzz_0_y_1,   \
                             g_0_zzz_0_z_0,   \
                             g_0_zzz_0_z_1,   \
                             wp_x,            \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxzzz_0_x_0[i] =
            2.0 * g_0_xxz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_x_1[i] * fti_ab_0 + g_0_xxzz_0_x_0[i] * pb_z + g_0_xxzz_0_x_1[i] * wp_z[i];

        g_0_xxzzz_0_y_0[i] = g_0_zzz_0_y_0[i] * fi_ab_0 - g_0_zzz_0_y_1[i] * fti_ab_0 + g_0_xzzz_0_y_0[i] * pb_x + g_0_xzzz_0_y_1[i] * wp_x[i];

        g_0_xxzzz_0_z_0[i] = g_0_zzz_0_z_0[i] * fi_ab_0 - g_0_zzz_0_z_1[i] * fti_ab_0 + g_0_xzzz_0_z_0[i] * pb_x + g_0_xzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 30-33 components of targeted buffer : SHSP

    auto g_0_xyyyy_0_x_0 = pbuffer.data(idx_eri_0_shsp + 30);

    auto g_0_xyyyy_0_y_0 = pbuffer.data(idx_eri_0_shsp + 31);

    auto g_0_xyyyy_0_z_0 = pbuffer.data(idx_eri_0_shsp + 32);

#pragma omp simd aligned(g_0_xyyyy_0_x_0,     \
                             g_0_xyyyy_0_y_0, \
                             g_0_xyyyy_0_z_0, \
                             g_0_yyyy_0_0_1,  \
                             g_0_yyyy_0_x_0,  \
                             g_0_yyyy_0_x_1,  \
                             g_0_yyyy_0_y_0,  \
                             g_0_yyyy_0_y_1,  \
                             g_0_yyyy_0_z_0,  \
                             g_0_yyyy_0_z_1,  \
                             wp_x,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyy_0_x_0[i] = g_0_yyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyy_0_x_0[i] * pb_x + g_0_yyyy_0_x_1[i] * wp_x[i];

        g_0_xyyyy_0_y_0[i] = g_0_yyyy_0_y_0[i] * pb_x + g_0_yyyy_0_y_1[i] * wp_x[i];

        g_0_xyyyy_0_z_0[i] = g_0_yyyy_0_z_0[i] * pb_x + g_0_yyyy_0_z_1[i] * wp_x[i];
    }

    /// Set up 33-36 components of targeted buffer : SHSP

    auto g_0_xyyyz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 33);

    auto g_0_xyyyz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 34);

    auto g_0_xyyyz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 35);

#pragma omp simd aligned(g_0_xyyy_0_x_0,      \
                             g_0_xyyy_0_x_1,  \
                             g_0_xyyyz_0_x_0, \
                             g_0_xyyyz_0_y_0, \
                             g_0_xyyyz_0_z_0, \
                             g_0_yyyz_0_y_0,  \
                             g_0_yyyz_0_y_1,  \
                             g_0_yyyz_0_z_0,  \
                             g_0_yyyz_0_z_1,  \
                             wp_x,            \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xyyyz_0_x_0[i] = g_0_xyyy_0_x_0[i] * pb_z + g_0_xyyy_0_x_1[i] * wp_z[i];

        g_0_xyyyz_0_y_0[i] = g_0_yyyz_0_y_0[i] * pb_x + g_0_yyyz_0_y_1[i] * wp_x[i];

        g_0_xyyyz_0_z_0[i] = g_0_yyyz_0_z_0[i] * pb_x + g_0_yyyz_0_z_1[i] * wp_x[i];
    }

    /// Set up 36-39 components of targeted buffer : SHSP

    auto g_0_xyyzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 36);

    auto g_0_xyyzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 37);

    auto g_0_xyyzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 38);

#pragma omp simd aligned(g_0_xyyzz_0_x_0,     \
                             g_0_xyyzz_0_y_0, \
                             g_0_xyyzz_0_z_0, \
                             g_0_yyzz_0_0_1,  \
                             g_0_yyzz_0_x_0,  \
                             g_0_yyzz_0_x_1,  \
                             g_0_yyzz_0_y_0,  \
                             g_0_yyzz_0_y_1,  \
                             g_0_yyzz_0_z_0,  \
                             g_0_yyzz_0_z_1,  \
                             wp_x,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzz_0_x_0[i] = g_0_yyzz_0_0_1[i] * fi_abcd_0 + g_0_yyzz_0_x_0[i] * pb_x + g_0_yyzz_0_x_1[i] * wp_x[i];

        g_0_xyyzz_0_y_0[i] = g_0_yyzz_0_y_0[i] * pb_x + g_0_yyzz_0_y_1[i] * wp_x[i];

        g_0_xyyzz_0_z_0[i] = g_0_yyzz_0_z_0[i] * pb_x + g_0_yyzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 39-42 components of targeted buffer : SHSP

    auto g_0_xyzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 39);

    auto g_0_xyzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 40);

    auto g_0_xyzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 41);

#pragma omp simd aligned(g_0_xyzzz_0_x_0,     \
                             g_0_xyzzz_0_y_0, \
                             g_0_xyzzz_0_z_0, \
                             g_0_xzzz_0_x_0,  \
                             g_0_xzzz_0_x_1,  \
                             g_0_yzzz_0_y_0,  \
                             g_0_yzzz_0_y_1,  \
                             g_0_yzzz_0_z_0,  \
                             g_0_yzzz_0_z_1,  \
                             wp_x,            \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xyzzz_0_x_0[i] = g_0_xzzz_0_x_0[i] * pb_y + g_0_xzzz_0_x_1[i] * wp_y[i];

        g_0_xyzzz_0_y_0[i] = g_0_yzzz_0_y_0[i] * pb_x + g_0_yzzz_0_y_1[i] * wp_x[i];

        g_0_xyzzz_0_z_0[i] = g_0_yzzz_0_z_0[i] * pb_x + g_0_yzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 42-45 components of targeted buffer : SHSP

    auto g_0_xzzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 42);

    auto g_0_xzzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 43);

    auto g_0_xzzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 44);

#pragma omp simd aligned(g_0_xzzzz_0_x_0,     \
                             g_0_xzzzz_0_y_0, \
                             g_0_xzzzz_0_z_0, \
                             g_0_zzzz_0_0_1,  \
                             g_0_zzzz_0_x_0,  \
                             g_0_zzzz_0_x_1,  \
                             g_0_zzzz_0_y_0,  \
                             g_0_zzzz_0_y_1,  \
                             g_0_zzzz_0_z_0,  \
                             g_0_zzzz_0_z_1,  \
                             wp_x,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzz_0_x_0[i] = g_0_zzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzz_0_x_0[i] * pb_x + g_0_zzzz_0_x_1[i] * wp_x[i];

        g_0_xzzzz_0_y_0[i] = g_0_zzzz_0_y_0[i] * pb_x + g_0_zzzz_0_y_1[i] * wp_x[i];

        g_0_xzzzz_0_z_0[i] = g_0_zzzz_0_z_0[i] * pb_x + g_0_zzzz_0_z_1[i] * wp_x[i];
    }

    /// Set up 45-48 components of targeted buffer : SHSP

    auto g_0_yyyyy_0_x_0 = pbuffer.data(idx_eri_0_shsp + 45);

    auto g_0_yyyyy_0_y_0 = pbuffer.data(idx_eri_0_shsp + 46);

    auto g_0_yyyyy_0_z_0 = pbuffer.data(idx_eri_0_shsp + 47);

#pragma omp simd aligned(g_0_yyy_0_x_0,       \
                             g_0_yyy_0_x_1,   \
                             g_0_yyy_0_y_0,   \
                             g_0_yyy_0_y_1,   \
                             g_0_yyy_0_z_0,   \
                             g_0_yyy_0_z_1,   \
                             g_0_yyyy_0_0_1,  \
                             g_0_yyyy_0_x_0,  \
                             g_0_yyyy_0_x_1,  \
                             g_0_yyyy_0_y_0,  \
                             g_0_yyyy_0_y_1,  \
                             g_0_yyyy_0_z_0,  \
                             g_0_yyyy_0_z_1,  \
                             g_0_yyyyy_0_x_0, \
                             g_0_yyyyy_0_y_0, \
                             g_0_yyyyy_0_z_0, \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyy_0_x_0[i] =
            4.0 * g_0_yyy_0_x_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_x_1[i] * fti_ab_0 + g_0_yyyy_0_x_0[i] * pb_y + g_0_yyyy_0_x_1[i] * wp_y[i];

        g_0_yyyyy_0_y_0[i] = 4.0 * g_0_yyy_0_y_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_y_1[i] * fti_ab_0 + g_0_yyyy_0_0_1[i] * fi_abcd_0 +
                             g_0_yyyy_0_y_0[i] * pb_y + g_0_yyyy_0_y_1[i] * wp_y[i];

        g_0_yyyyy_0_z_0[i] =
            4.0 * g_0_yyy_0_z_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_z_1[i] * fti_ab_0 + g_0_yyyy_0_z_0[i] * pb_y + g_0_yyyy_0_z_1[i] * wp_y[i];
    }

    /// Set up 48-51 components of targeted buffer : SHSP

    auto g_0_yyyyz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 48);

    auto g_0_yyyyz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 49);

    auto g_0_yyyyz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 50);

#pragma omp simd aligned(g_0_yyyy_0_0_1,      \
                             g_0_yyyy_0_x_0,  \
                             g_0_yyyy_0_x_1,  \
                             g_0_yyyy_0_y_0,  \
                             g_0_yyyy_0_y_1,  \
                             g_0_yyyy_0_z_0,  \
                             g_0_yyyy_0_z_1,  \
                             g_0_yyyyz_0_x_0, \
                             g_0_yyyyz_0_y_0, \
                             g_0_yyyyz_0_z_0, \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyz_0_x_0[i] = g_0_yyyy_0_x_0[i] * pb_z + g_0_yyyy_0_x_1[i] * wp_z[i];

        g_0_yyyyz_0_y_0[i] = g_0_yyyy_0_y_0[i] * pb_z + g_0_yyyy_0_y_1[i] * wp_z[i];

        g_0_yyyyz_0_z_0[i] = g_0_yyyy_0_0_1[i] * fi_abcd_0 + g_0_yyyy_0_z_0[i] * pb_z + g_0_yyyy_0_z_1[i] * wp_z[i];
    }

    /// Set up 51-54 components of targeted buffer : SHSP

    auto g_0_yyyzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 51);

    auto g_0_yyyzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 52);

    auto g_0_yyyzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 53);

#pragma omp simd aligned(g_0_yyy_0_y_0,       \
                             g_0_yyy_0_y_1,   \
                             g_0_yyyz_0_y_0,  \
                             g_0_yyyz_0_y_1,  \
                             g_0_yyyzz_0_x_0, \
                             g_0_yyyzz_0_y_0, \
                             g_0_yyyzz_0_z_0, \
                             g_0_yyzz_0_x_0,  \
                             g_0_yyzz_0_x_1,  \
                             g_0_yyzz_0_z_0,  \
                             g_0_yyzz_0_z_1,  \
                             g_0_yzz_0_x_0,   \
                             g_0_yzz_0_x_1,   \
                             g_0_yzz_0_z_0,   \
                             g_0_yzz_0_z_1,   \
                             wp_y,            \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyzz_0_x_0[i] =
            2.0 * g_0_yzz_0_x_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_x_1[i] * fti_ab_0 + g_0_yyzz_0_x_0[i] * pb_y + g_0_yyzz_0_x_1[i] * wp_y[i];

        g_0_yyyzz_0_y_0[i] = g_0_yyy_0_y_0[i] * fi_ab_0 - g_0_yyy_0_y_1[i] * fti_ab_0 + g_0_yyyz_0_y_0[i] * pb_z + g_0_yyyz_0_y_1[i] * wp_z[i];

        g_0_yyyzz_0_z_0[i] =
            2.0 * g_0_yzz_0_z_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_z_1[i] * fti_ab_0 + g_0_yyzz_0_z_0[i] * pb_y + g_0_yyzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 54-57 components of targeted buffer : SHSP

    auto g_0_yyzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 54);

    auto g_0_yyzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 55);

    auto g_0_yyzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 56);

#pragma omp simd aligned(g_0_yyz_0_y_0,       \
                             g_0_yyz_0_y_1,   \
                             g_0_yyzz_0_y_0,  \
                             g_0_yyzz_0_y_1,  \
                             g_0_yyzzz_0_x_0, \
                             g_0_yyzzz_0_y_0, \
                             g_0_yyzzz_0_z_0, \
                             g_0_yzzz_0_x_0,  \
                             g_0_yzzz_0_x_1,  \
                             g_0_yzzz_0_z_0,  \
                             g_0_yzzz_0_z_1,  \
                             g_0_zzz_0_x_0,   \
                             g_0_zzz_0_x_1,   \
                             g_0_zzz_0_z_0,   \
                             g_0_zzz_0_z_1,   \
                             wp_y,            \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyzzz_0_x_0[i] = g_0_zzz_0_x_0[i] * fi_ab_0 - g_0_zzz_0_x_1[i] * fti_ab_0 + g_0_yzzz_0_x_0[i] * pb_y + g_0_yzzz_0_x_1[i] * wp_y[i];

        g_0_yyzzz_0_y_0[i] =
            2.0 * g_0_yyz_0_y_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_y_1[i] * fti_ab_0 + g_0_yyzz_0_y_0[i] * pb_z + g_0_yyzz_0_y_1[i] * wp_z[i];

        g_0_yyzzz_0_z_0[i] = g_0_zzz_0_z_0[i] * fi_ab_0 - g_0_zzz_0_z_1[i] * fti_ab_0 + g_0_yzzz_0_z_0[i] * pb_y + g_0_yzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 57-60 components of targeted buffer : SHSP

    auto g_0_yzzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 57);

    auto g_0_yzzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 58);

    auto g_0_yzzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 59);

#pragma omp simd aligned(g_0_yzzzz_0_x_0,     \
                             g_0_yzzzz_0_y_0, \
                             g_0_yzzzz_0_z_0, \
                             g_0_zzzz_0_0_1,  \
                             g_0_zzzz_0_x_0,  \
                             g_0_zzzz_0_x_1,  \
                             g_0_zzzz_0_y_0,  \
                             g_0_zzzz_0_y_1,  \
                             g_0_zzzz_0_z_0,  \
                             g_0_zzzz_0_z_1,  \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzz_0_x_0[i] = g_0_zzzz_0_x_0[i] * pb_y + g_0_zzzz_0_x_1[i] * wp_y[i];

        g_0_yzzzz_0_y_0[i] = g_0_zzzz_0_0_1[i] * fi_abcd_0 + g_0_zzzz_0_y_0[i] * pb_y + g_0_zzzz_0_y_1[i] * wp_y[i];

        g_0_yzzzz_0_z_0[i] = g_0_zzzz_0_z_0[i] * pb_y + g_0_zzzz_0_z_1[i] * wp_y[i];
    }

    /// Set up 60-63 components of targeted buffer : SHSP

    auto g_0_zzzzz_0_x_0 = pbuffer.data(idx_eri_0_shsp + 60);

    auto g_0_zzzzz_0_y_0 = pbuffer.data(idx_eri_0_shsp + 61);

    auto g_0_zzzzz_0_z_0 = pbuffer.data(idx_eri_0_shsp + 62);

#pragma omp simd aligned(g_0_zzz_0_x_0,       \
                             g_0_zzz_0_x_1,   \
                             g_0_zzz_0_y_0,   \
                             g_0_zzz_0_y_1,   \
                             g_0_zzz_0_z_0,   \
                             g_0_zzz_0_z_1,   \
                             g_0_zzzz_0_0_1,  \
                             g_0_zzzz_0_x_0,  \
                             g_0_zzzz_0_x_1,  \
                             g_0_zzzz_0_y_0,  \
                             g_0_zzzz_0_y_1,  \
                             g_0_zzzz_0_z_0,  \
                             g_0_zzzz_0_z_1,  \
                             g_0_zzzzz_0_x_0, \
                             g_0_zzzzz_0_y_0, \
                             g_0_zzzzz_0_z_0, \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzz_0_x_0[i] =
            4.0 * g_0_zzz_0_x_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_x_1[i] * fti_ab_0 + g_0_zzzz_0_x_0[i] * pb_z + g_0_zzzz_0_x_1[i] * wp_z[i];

        g_0_zzzzz_0_y_0[i] =
            4.0 * g_0_zzz_0_y_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_y_1[i] * fti_ab_0 + g_0_zzzz_0_y_0[i] * pb_z + g_0_zzzz_0_y_1[i] * wp_z[i];

        g_0_zzzzz_0_z_0[i] = 4.0 * g_0_zzz_0_z_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_z_1[i] * fti_ab_0 + g_0_zzzz_0_0_1[i] * fi_abcd_0 +
                             g_0_zzzz_0_z_0[i] * pb_z + g_0_zzzz_0_z_1[i] * wp_z[i];
    }
}

}  // namespace erirec
