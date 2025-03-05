#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_psg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psg,
                                 size_t idx_eri_1_ssf,
                                 size_t idx_eri_1_ssg,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : SSF

    auto g_0_0_xxx_1 = pbuffer.data(idx_eri_1_ssf);

    auto g_0_0_xxy_1 = pbuffer.data(idx_eri_1_ssf + 1);

    auto g_0_0_xxz_1 = pbuffer.data(idx_eri_1_ssf + 2);

    auto g_0_0_xyy_1 = pbuffer.data(idx_eri_1_ssf + 3);

    auto g_0_0_xyz_1 = pbuffer.data(idx_eri_1_ssf + 4);

    auto g_0_0_xzz_1 = pbuffer.data(idx_eri_1_ssf + 5);

    auto g_0_0_yyy_1 = pbuffer.data(idx_eri_1_ssf + 6);

    auto g_0_0_yyz_1 = pbuffer.data(idx_eri_1_ssf + 7);

    auto g_0_0_yzz_1 = pbuffer.data(idx_eri_1_ssf + 8);

    auto g_0_0_zzz_1 = pbuffer.data(idx_eri_1_ssf + 9);

    /// Set up components of auxilary buffer : SSG

    auto g_0_0_xxxx_1 = pbuffer.data(idx_eri_1_ssg);

    auto g_0_0_xxxy_1 = pbuffer.data(idx_eri_1_ssg + 1);

    auto g_0_0_xxxz_1 = pbuffer.data(idx_eri_1_ssg + 2);

    auto g_0_0_xxyy_1 = pbuffer.data(idx_eri_1_ssg + 3);

    auto g_0_0_xxyz_1 = pbuffer.data(idx_eri_1_ssg + 4);

    auto g_0_0_xxzz_1 = pbuffer.data(idx_eri_1_ssg + 5);

    auto g_0_0_xyyy_1 = pbuffer.data(idx_eri_1_ssg + 6);

    auto g_0_0_xyyz_1 = pbuffer.data(idx_eri_1_ssg + 7);

    auto g_0_0_xyzz_1 = pbuffer.data(idx_eri_1_ssg + 8);

    auto g_0_0_xzzz_1 = pbuffer.data(idx_eri_1_ssg + 9);

    auto g_0_0_yyyy_1 = pbuffer.data(idx_eri_1_ssg + 10);

    auto g_0_0_yyyz_1 = pbuffer.data(idx_eri_1_ssg + 11);

    auto g_0_0_yyzz_1 = pbuffer.data(idx_eri_1_ssg + 12);

    auto g_0_0_yzzz_1 = pbuffer.data(idx_eri_1_ssg + 13);

    auto g_0_0_zzzz_1 = pbuffer.data(idx_eri_1_ssg + 14);

    /// Set up 0-15 components of targeted buffer : PSG

    auto g_x_0_xxxx_0 = pbuffer.data(idx_eri_0_psg);

    auto g_x_0_xxxy_0 = pbuffer.data(idx_eri_0_psg + 1);

    auto g_x_0_xxxz_0 = pbuffer.data(idx_eri_0_psg + 2);

    auto g_x_0_xxyy_0 = pbuffer.data(idx_eri_0_psg + 3);

    auto g_x_0_xxyz_0 = pbuffer.data(idx_eri_0_psg + 4);

    auto g_x_0_xxzz_0 = pbuffer.data(idx_eri_0_psg + 5);

    auto g_x_0_xyyy_0 = pbuffer.data(idx_eri_0_psg + 6);

    auto g_x_0_xyyz_0 = pbuffer.data(idx_eri_0_psg + 7);

    auto g_x_0_xyzz_0 = pbuffer.data(idx_eri_0_psg + 8);

    auto g_x_0_xzzz_0 = pbuffer.data(idx_eri_0_psg + 9);

    auto g_x_0_yyyy_0 = pbuffer.data(idx_eri_0_psg + 10);

    auto g_x_0_yyyz_0 = pbuffer.data(idx_eri_0_psg + 11);

    auto g_x_0_yyzz_0 = pbuffer.data(idx_eri_0_psg + 12);

    auto g_x_0_yzzz_0 = pbuffer.data(idx_eri_0_psg + 13);

    auto g_x_0_zzzz_0 = pbuffer.data(idx_eri_0_psg + 14);

    #pragma omp simd aligned(g_0_0_xxx_1, g_0_0_xxxx_1, g_0_0_xxxy_1, g_0_0_xxxz_1, g_0_0_xxy_1, g_0_0_xxyy_1, g_0_0_xxyz_1, g_0_0_xxz_1, g_0_0_xxzz_1, g_0_0_xyy_1, g_0_0_xyyy_1, g_0_0_xyyz_1, g_0_0_xyz_1, g_0_0_xyzz_1, g_0_0_xzz_1, g_0_0_xzzz_1, g_0_0_yyy_1, g_0_0_yyyy_1, g_0_0_yyyz_1, g_0_0_yyz_1, g_0_0_yyzz_1, g_0_0_yzz_1, g_0_0_yzzz_1, g_0_0_zzz_1, g_0_0_zzzz_1, g_x_0_xxxx_0, g_x_0_xxxy_0, g_x_0_xxxz_0, g_x_0_xxyy_0, g_x_0_xxyz_0, g_x_0_xxzz_0, g_x_0_xyyy_0, g_x_0_xyyz_0, g_x_0_xyzz_0, g_x_0_xzzz_0, g_x_0_yyyy_0, g_x_0_yyyz_0, g_x_0_yyzz_0, g_x_0_yzzz_0, g_x_0_zzzz_0, wa_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_x_0_xxxx_0[i] = 4.0 * g_0_0_xxx_1[i] * fi_acd_0 + g_0_0_xxxx_1[i] * wa_x[i];

        g_x_0_xxxy_0[i] = 3.0 * g_0_0_xxy_1[i] * fi_acd_0 + g_0_0_xxxy_1[i] * wa_x[i];

        g_x_0_xxxz_0[i] = 3.0 * g_0_0_xxz_1[i] * fi_acd_0 + g_0_0_xxxz_1[i] * wa_x[i];

        g_x_0_xxyy_0[i] = 2.0 * g_0_0_xyy_1[i] * fi_acd_0 + g_0_0_xxyy_1[i] * wa_x[i];

        g_x_0_xxyz_0[i] = 2.0 * g_0_0_xyz_1[i] * fi_acd_0 + g_0_0_xxyz_1[i] * wa_x[i];

        g_x_0_xxzz_0[i] = 2.0 * g_0_0_xzz_1[i] * fi_acd_0 + g_0_0_xxzz_1[i] * wa_x[i];

        g_x_0_xyyy_0[i] = g_0_0_yyy_1[i] * fi_acd_0 + g_0_0_xyyy_1[i] * wa_x[i];

        g_x_0_xyyz_0[i] = g_0_0_yyz_1[i] * fi_acd_0 + g_0_0_xyyz_1[i] * wa_x[i];

        g_x_0_xyzz_0[i] = g_0_0_yzz_1[i] * fi_acd_0 + g_0_0_xyzz_1[i] * wa_x[i];

        g_x_0_xzzz_0[i] = g_0_0_zzz_1[i] * fi_acd_0 + g_0_0_xzzz_1[i] * wa_x[i];

        g_x_0_yyyy_0[i] = g_0_0_yyyy_1[i] * wa_x[i];

        g_x_0_yyyz_0[i] = g_0_0_yyyz_1[i] * wa_x[i];

        g_x_0_yyzz_0[i] = g_0_0_yyzz_1[i] * wa_x[i];

        g_x_0_yzzz_0[i] = g_0_0_yzzz_1[i] * wa_x[i];

        g_x_0_zzzz_0[i] = g_0_0_zzzz_1[i] * wa_x[i];
    }

    /// Set up 15-30 components of targeted buffer : PSG

    auto g_y_0_xxxx_0 = pbuffer.data(idx_eri_0_psg + 15);

    auto g_y_0_xxxy_0 = pbuffer.data(idx_eri_0_psg + 16);

    auto g_y_0_xxxz_0 = pbuffer.data(idx_eri_0_psg + 17);

    auto g_y_0_xxyy_0 = pbuffer.data(idx_eri_0_psg + 18);

    auto g_y_0_xxyz_0 = pbuffer.data(idx_eri_0_psg + 19);

    auto g_y_0_xxzz_0 = pbuffer.data(idx_eri_0_psg + 20);

    auto g_y_0_xyyy_0 = pbuffer.data(idx_eri_0_psg + 21);

    auto g_y_0_xyyz_0 = pbuffer.data(idx_eri_0_psg + 22);

    auto g_y_0_xyzz_0 = pbuffer.data(idx_eri_0_psg + 23);

    auto g_y_0_xzzz_0 = pbuffer.data(idx_eri_0_psg + 24);

    auto g_y_0_yyyy_0 = pbuffer.data(idx_eri_0_psg + 25);

    auto g_y_0_yyyz_0 = pbuffer.data(idx_eri_0_psg + 26);

    auto g_y_0_yyzz_0 = pbuffer.data(idx_eri_0_psg + 27);

    auto g_y_0_yzzz_0 = pbuffer.data(idx_eri_0_psg + 28);

    auto g_y_0_zzzz_0 = pbuffer.data(idx_eri_0_psg + 29);

    #pragma omp simd aligned(g_0_0_xxx_1, g_0_0_xxxx_1, g_0_0_xxxy_1, g_0_0_xxxz_1, g_0_0_xxy_1, g_0_0_xxyy_1, g_0_0_xxyz_1, g_0_0_xxz_1, g_0_0_xxzz_1, g_0_0_xyy_1, g_0_0_xyyy_1, g_0_0_xyyz_1, g_0_0_xyz_1, g_0_0_xyzz_1, g_0_0_xzz_1, g_0_0_xzzz_1, g_0_0_yyy_1, g_0_0_yyyy_1, g_0_0_yyyz_1, g_0_0_yyz_1, g_0_0_yyzz_1, g_0_0_yzz_1, g_0_0_yzzz_1, g_0_0_zzz_1, g_0_0_zzzz_1, g_y_0_xxxx_0, g_y_0_xxxy_0, g_y_0_xxxz_0, g_y_0_xxyy_0, g_y_0_xxyz_0, g_y_0_xxzz_0, g_y_0_xyyy_0, g_y_0_xyyz_0, g_y_0_xyzz_0, g_y_0_xzzz_0, g_y_0_yyyy_0, g_y_0_yyyz_0, g_y_0_yyzz_0, g_y_0_yzzz_0, g_y_0_zzzz_0, wa_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_y_0_xxxx_0[i] = g_0_0_xxxx_1[i] * wa_y[i];

        g_y_0_xxxy_0[i] = g_0_0_xxx_1[i] * fi_acd_0 + g_0_0_xxxy_1[i] * wa_y[i];

        g_y_0_xxxz_0[i] = g_0_0_xxxz_1[i] * wa_y[i];

        g_y_0_xxyy_0[i] = 2.0 * g_0_0_xxy_1[i] * fi_acd_0 + g_0_0_xxyy_1[i] * wa_y[i];

        g_y_0_xxyz_0[i] = g_0_0_xxz_1[i] * fi_acd_0 + g_0_0_xxyz_1[i] * wa_y[i];

        g_y_0_xxzz_0[i] = g_0_0_xxzz_1[i] * wa_y[i];

        g_y_0_xyyy_0[i] = 3.0 * g_0_0_xyy_1[i] * fi_acd_0 + g_0_0_xyyy_1[i] * wa_y[i];

        g_y_0_xyyz_0[i] = 2.0 * g_0_0_xyz_1[i] * fi_acd_0 + g_0_0_xyyz_1[i] * wa_y[i];

        g_y_0_xyzz_0[i] = g_0_0_xzz_1[i] * fi_acd_0 + g_0_0_xyzz_1[i] * wa_y[i];

        g_y_0_xzzz_0[i] = g_0_0_xzzz_1[i] * wa_y[i];

        g_y_0_yyyy_0[i] = 4.0 * g_0_0_yyy_1[i] * fi_acd_0 + g_0_0_yyyy_1[i] * wa_y[i];

        g_y_0_yyyz_0[i] = 3.0 * g_0_0_yyz_1[i] * fi_acd_0 + g_0_0_yyyz_1[i] * wa_y[i];

        g_y_0_yyzz_0[i] = 2.0 * g_0_0_yzz_1[i] * fi_acd_0 + g_0_0_yyzz_1[i] * wa_y[i];

        g_y_0_yzzz_0[i] = g_0_0_zzz_1[i] * fi_acd_0 + g_0_0_yzzz_1[i] * wa_y[i];

        g_y_0_zzzz_0[i] = g_0_0_zzzz_1[i] * wa_y[i];
    }

    /// Set up 30-45 components of targeted buffer : PSG

    auto g_z_0_xxxx_0 = pbuffer.data(idx_eri_0_psg + 30);

    auto g_z_0_xxxy_0 = pbuffer.data(idx_eri_0_psg + 31);

    auto g_z_0_xxxz_0 = pbuffer.data(idx_eri_0_psg + 32);

    auto g_z_0_xxyy_0 = pbuffer.data(idx_eri_0_psg + 33);

    auto g_z_0_xxyz_0 = pbuffer.data(idx_eri_0_psg + 34);

    auto g_z_0_xxzz_0 = pbuffer.data(idx_eri_0_psg + 35);

    auto g_z_0_xyyy_0 = pbuffer.data(idx_eri_0_psg + 36);

    auto g_z_0_xyyz_0 = pbuffer.data(idx_eri_0_psg + 37);

    auto g_z_0_xyzz_0 = pbuffer.data(idx_eri_0_psg + 38);

    auto g_z_0_xzzz_0 = pbuffer.data(idx_eri_0_psg + 39);

    auto g_z_0_yyyy_0 = pbuffer.data(idx_eri_0_psg + 40);

    auto g_z_0_yyyz_0 = pbuffer.data(idx_eri_0_psg + 41);

    auto g_z_0_yyzz_0 = pbuffer.data(idx_eri_0_psg + 42);

    auto g_z_0_yzzz_0 = pbuffer.data(idx_eri_0_psg + 43);

    auto g_z_0_zzzz_0 = pbuffer.data(idx_eri_0_psg + 44);

    #pragma omp simd aligned(g_0_0_xxx_1, g_0_0_xxxx_1, g_0_0_xxxy_1, g_0_0_xxxz_1, g_0_0_xxy_1, g_0_0_xxyy_1, g_0_0_xxyz_1, g_0_0_xxz_1, g_0_0_xxzz_1, g_0_0_xyy_1, g_0_0_xyyy_1, g_0_0_xyyz_1, g_0_0_xyz_1, g_0_0_xyzz_1, g_0_0_xzz_1, g_0_0_xzzz_1, g_0_0_yyy_1, g_0_0_yyyy_1, g_0_0_yyyz_1, g_0_0_yyz_1, g_0_0_yyzz_1, g_0_0_yzz_1, g_0_0_yzzz_1, g_0_0_zzz_1, g_0_0_zzzz_1, g_z_0_xxxx_0, g_z_0_xxxy_0, g_z_0_xxxz_0, g_z_0_xxyy_0, g_z_0_xxyz_0, g_z_0_xxzz_0, g_z_0_xyyy_0, g_z_0_xyyz_0, g_z_0_xyzz_0, g_z_0_xzzz_0, g_z_0_yyyy_0, g_z_0_yyyz_0, g_z_0_yyzz_0, g_z_0_yzzz_0, g_z_0_zzzz_0, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_z_0_xxxx_0[i] = g_0_0_xxxx_1[i] * wa_z[i];

        g_z_0_xxxy_0[i] = g_0_0_xxxy_1[i] * wa_z[i];

        g_z_0_xxxz_0[i] = g_0_0_xxx_1[i] * fi_acd_0 + g_0_0_xxxz_1[i] * wa_z[i];

        g_z_0_xxyy_0[i] = g_0_0_xxyy_1[i] * wa_z[i];

        g_z_0_xxyz_0[i] = g_0_0_xxy_1[i] * fi_acd_0 + g_0_0_xxyz_1[i] * wa_z[i];

        g_z_0_xxzz_0[i] = 2.0 * g_0_0_xxz_1[i] * fi_acd_0 + g_0_0_xxzz_1[i] * wa_z[i];

        g_z_0_xyyy_0[i] = g_0_0_xyyy_1[i] * wa_z[i];

        g_z_0_xyyz_0[i] = g_0_0_xyy_1[i] * fi_acd_0 + g_0_0_xyyz_1[i] * wa_z[i];

        g_z_0_xyzz_0[i] = 2.0 * g_0_0_xyz_1[i] * fi_acd_0 + g_0_0_xyzz_1[i] * wa_z[i];

        g_z_0_xzzz_0[i] = 3.0 * g_0_0_xzz_1[i] * fi_acd_0 + g_0_0_xzzz_1[i] * wa_z[i];

        g_z_0_yyyy_0[i] = g_0_0_yyyy_1[i] * wa_z[i];

        g_z_0_yyyz_0[i] = g_0_0_yyy_1[i] * fi_acd_0 + g_0_0_yyyz_1[i] * wa_z[i];

        g_z_0_yyzz_0[i] = 2.0 * g_0_0_yyz_1[i] * fi_acd_0 + g_0_0_yyzz_1[i] * wa_z[i];

        g_z_0_yzzz_0[i] = 3.0 * g_0_0_yzz_1[i] * fi_acd_0 + g_0_0_yzzz_1[i] * wa_z[i];

        g_z_0_zzzz_0[i] = 4.0 * g_0_0_zzz_1[i] * fi_acd_0 + g_0_0_zzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

