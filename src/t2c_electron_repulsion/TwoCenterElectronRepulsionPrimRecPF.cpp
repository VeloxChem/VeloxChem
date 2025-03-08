#include "TwoCenterElectronRepulsionPrimRecPF.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_pf(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_pf,
                                const size_t idx_eri_1_sd,
                                const size_t idx_eri_1_sf,
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

    // Set up components of auxiliary buffer : SD

    auto g_0_xx_1 = pbuffer.data(idx_eri_1_sd);

    auto g_0_xy_1 = pbuffer.data(idx_eri_1_sd + 1);

    auto g_0_xz_1 = pbuffer.data(idx_eri_1_sd + 2);

    auto g_0_yy_1 = pbuffer.data(idx_eri_1_sd + 3);

    auto g_0_yz_1 = pbuffer.data(idx_eri_1_sd + 4);

    auto g_0_zz_1 = pbuffer.data(idx_eri_1_sd + 5);

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

    // Set up 0-10 components of targeted buffer : PF

    auto g_x_xxx_0 = pbuffer.data(idx_eri_0_pf);

    auto g_x_xxy_0 = pbuffer.data(idx_eri_0_pf + 1);

    auto g_x_xxz_0 = pbuffer.data(idx_eri_0_pf + 2);

    auto g_x_xyy_0 = pbuffer.data(idx_eri_0_pf + 3);

    auto g_x_xyz_0 = pbuffer.data(idx_eri_0_pf + 4);

    auto g_x_xzz_0 = pbuffer.data(idx_eri_0_pf + 5);

    auto g_x_yyy_0 = pbuffer.data(idx_eri_0_pf + 6);

    auto g_x_yyz_0 = pbuffer.data(idx_eri_0_pf + 7);

    auto g_x_yzz_0 = pbuffer.data(idx_eri_0_pf + 8);

    auto g_x_zzz_0 = pbuffer.data(idx_eri_0_pf + 9);

    #pragma omp simd aligned(g_0_xx_1, g_0_xxx_1, g_0_xxy_1, g_0_xxz_1, g_0_xy_1, g_0_xyy_1, g_0_xyz_1, g_0_xz_1, g_0_xzz_1, g_0_yy_1, g_0_yyy_1, g_0_yyz_1, g_0_yz_1, g_0_yzz_1, g_0_zz_1, g_0_zzz_1, g_x_xxx_0, g_x_xxy_0, g_x_xxz_0, g_x_xyy_0, g_x_xyz_0, g_x_xzz_0, g_x_yyy_0, g_x_yyz_0, g_x_yzz_0, g_x_zzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_x_xxx_0[i] = 3.0 * g_0_xx_1[i] * fe_0 + g_0_xxx_1[i] * pa_x[i];

        g_x_xxy_0[i] = 2.0 * g_0_xy_1[i] * fe_0 + g_0_xxy_1[i] * pa_x[i];

        g_x_xxz_0[i] = 2.0 * g_0_xz_1[i] * fe_0 + g_0_xxz_1[i] * pa_x[i];

        g_x_xyy_0[i] = g_0_yy_1[i] * fe_0 + g_0_xyy_1[i] * pa_x[i];

        g_x_xyz_0[i] = g_0_yz_1[i] * fe_0 + g_0_xyz_1[i] * pa_x[i];

        g_x_xzz_0[i] = g_0_zz_1[i] * fe_0 + g_0_xzz_1[i] * pa_x[i];

        g_x_yyy_0[i] = g_0_yyy_1[i] * pa_x[i];

        g_x_yyz_0[i] = g_0_yyz_1[i] * pa_x[i];

        g_x_yzz_0[i] = g_0_yzz_1[i] * pa_x[i];

        g_x_zzz_0[i] = g_0_zzz_1[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : PF

    auto g_y_xxx_0 = pbuffer.data(idx_eri_0_pf + 10);

    auto g_y_xxy_0 = pbuffer.data(idx_eri_0_pf + 11);

    auto g_y_xxz_0 = pbuffer.data(idx_eri_0_pf + 12);

    auto g_y_xyy_0 = pbuffer.data(idx_eri_0_pf + 13);

    auto g_y_xyz_0 = pbuffer.data(idx_eri_0_pf + 14);

    auto g_y_xzz_0 = pbuffer.data(idx_eri_0_pf + 15);

    auto g_y_yyy_0 = pbuffer.data(idx_eri_0_pf + 16);

    auto g_y_yyz_0 = pbuffer.data(idx_eri_0_pf + 17);

    auto g_y_yzz_0 = pbuffer.data(idx_eri_0_pf + 18);

    auto g_y_zzz_0 = pbuffer.data(idx_eri_0_pf + 19);

    #pragma omp simd aligned(g_0_xx_1, g_0_xxx_1, g_0_xxy_1, g_0_xxz_1, g_0_xy_1, g_0_xyy_1, g_0_xyz_1, g_0_xz_1, g_0_xzz_1, g_0_yy_1, g_0_yyy_1, g_0_yyz_1, g_0_yz_1, g_0_yzz_1, g_0_zz_1, g_0_zzz_1, g_y_xxx_0, g_y_xxy_0, g_y_xxz_0, g_y_xyy_0, g_y_xyz_0, g_y_xzz_0, g_y_yyy_0, g_y_yyz_0, g_y_yzz_0, g_y_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_y_xxx_0[i] = g_0_xxx_1[i] * pa_y[i];

        g_y_xxy_0[i] = g_0_xx_1[i] * fe_0 + g_0_xxy_1[i] * pa_y[i];

        g_y_xxz_0[i] = g_0_xxz_1[i] * pa_y[i];

        g_y_xyy_0[i] = 2.0 * g_0_xy_1[i] * fe_0 + g_0_xyy_1[i] * pa_y[i];

        g_y_xyz_0[i] = g_0_xz_1[i] * fe_0 + g_0_xyz_1[i] * pa_y[i];

        g_y_xzz_0[i] = g_0_xzz_1[i] * pa_y[i];

        g_y_yyy_0[i] = 3.0 * g_0_yy_1[i] * fe_0 + g_0_yyy_1[i] * pa_y[i];

        g_y_yyz_0[i] = 2.0 * g_0_yz_1[i] * fe_0 + g_0_yyz_1[i] * pa_y[i];

        g_y_yzz_0[i] = g_0_zz_1[i] * fe_0 + g_0_yzz_1[i] * pa_y[i];

        g_y_zzz_0[i] = g_0_zzz_1[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : PF

    auto g_z_xxx_0 = pbuffer.data(idx_eri_0_pf + 20);

    auto g_z_xxy_0 = pbuffer.data(idx_eri_0_pf + 21);

    auto g_z_xxz_0 = pbuffer.data(idx_eri_0_pf + 22);

    auto g_z_xyy_0 = pbuffer.data(idx_eri_0_pf + 23);

    auto g_z_xyz_0 = pbuffer.data(idx_eri_0_pf + 24);

    auto g_z_xzz_0 = pbuffer.data(idx_eri_0_pf + 25);

    auto g_z_yyy_0 = pbuffer.data(idx_eri_0_pf + 26);

    auto g_z_yyz_0 = pbuffer.data(idx_eri_0_pf + 27);

    auto g_z_yzz_0 = pbuffer.data(idx_eri_0_pf + 28);

    auto g_z_zzz_0 = pbuffer.data(idx_eri_0_pf + 29);

    #pragma omp simd aligned(g_0_xx_1, g_0_xxx_1, g_0_xxy_1, g_0_xxz_1, g_0_xy_1, g_0_xyy_1, g_0_xyz_1, g_0_xz_1, g_0_xzz_1, g_0_yy_1, g_0_yyy_1, g_0_yyz_1, g_0_yz_1, g_0_yzz_1, g_0_zz_1, g_0_zzz_1, g_z_xxx_0, g_z_xxy_0, g_z_xxz_0, g_z_xyy_0, g_z_xyz_0, g_z_xzz_0, g_z_yyy_0, g_z_yyz_0, g_z_yzz_0, g_z_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_z_xxx_0[i] = g_0_xxx_1[i] * pa_z[i];

        g_z_xxy_0[i] = g_0_xxy_1[i] * pa_z[i];

        g_z_xxz_0[i] = g_0_xx_1[i] * fe_0 + g_0_xxz_1[i] * pa_z[i];

        g_z_xyy_0[i] = g_0_xyy_1[i] * pa_z[i];

        g_z_xyz_0[i] = g_0_xy_1[i] * fe_0 + g_0_xyz_1[i] * pa_z[i];

        g_z_xzz_0[i] = 2.0 * g_0_xz_1[i] * fe_0 + g_0_xzz_1[i] * pa_z[i];

        g_z_yyy_0[i] = g_0_yyy_1[i] * pa_z[i];

        g_z_yyz_0[i] = g_0_yy_1[i] * fe_0 + g_0_yyz_1[i] * pa_z[i];

        g_z_yzz_0[i] = 2.0 * g_0_yz_1[i] * fe_0 + g_0_yzz_1[i] * pa_z[i];

        g_z_zzz_0[i] = 3.0 * g_0_zz_1[i] * fe_0 + g_0_zzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

