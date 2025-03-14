#include "TwoCenterElectronRepulsionPrimRecFD.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_fd(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_fd,
                                const size_t idx_eri_0_pd,
                                const size_t idx_eri_1_pd,
                                const size_t idx_eri_1_dp,
                                const size_t idx_eri_1_dd,
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

    // Set up components of auxiliary buffer : PD

    auto g_x_xx_0 = pbuffer.data(idx_eri_0_pd);

    auto g_x_xy_0 = pbuffer.data(idx_eri_0_pd + 1);

    auto g_x_xz_0 = pbuffer.data(idx_eri_0_pd + 2);

    auto g_x_yy_0 = pbuffer.data(idx_eri_0_pd + 3);

    auto g_x_yz_0 = pbuffer.data(idx_eri_0_pd + 4);

    auto g_x_zz_0 = pbuffer.data(idx_eri_0_pd + 5);

    auto g_y_xx_0 = pbuffer.data(idx_eri_0_pd + 6);

    auto g_y_xy_0 = pbuffer.data(idx_eri_0_pd + 7);

    auto g_y_xz_0 = pbuffer.data(idx_eri_0_pd + 8);

    auto g_y_yy_0 = pbuffer.data(idx_eri_0_pd + 9);

    auto g_y_yz_0 = pbuffer.data(idx_eri_0_pd + 10);

    auto g_y_zz_0 = pbuffer.data(idx_eri_0_pd + 11);

    auto g_z_xx_0 = pbuffer.data(idx_eri_0_pd + 12);

    auto g_z_xy_0 = pbuffer.data(idx_eri_0_pd + 13);

    auto g_z_xz_0 = pbuffer.data(idx_eri_0_pd + 14);

    auto g_z_yy_0 = pbuffer.data(idx_eri_0_pd + 15);

    auto g_z_yz_0 = pbuffer.data(idx_eri_0_pd + 16);

    auto g_z_zz_0 = pbuffer.data(idx_eri_0_pd + 17);

    // Set up components of auxiliary buffer : PD

    auto g_x_xx_1 = pbuffer.data(idx_eri_1_pd);

    auto g_x_xy_1 = pbuffer.data(idx_eri_1_pd + 1);

    auto g_x_xz_1 = pbuffer.data(idx_eri_1_pd + 2);

    auto g_x_yy_1 = pbuffer.data(idx_eri_1_pd + 3);

    auto g_x_yz_1 = pbuffer.data(idx_eri_1_pd + 4);

    auto g_x_zz_1 = pbuffer.data(idx_eri_1_pd + 5);

    auto g_y_xx_1 = pbuffer.data(idx_eri_1_pd + 6);

    auto g_y_xy_1 = pbuffer.data(idx_eri_1_pd + 7);

    auto g_y_xz_1 = pbuffer.data(idx_eri_1_pd + 8);

    auto g_y_yy_1 = pbuffer.data(idx_eri_1_pd + 9);

    auto g_y_yz_1 = pbuffer.data(idx_eri_1_pd + 10);

    auto g_y_zz_1 = pbuffer.data(idx_eri_1_pd + 11);

    auto g_z_xx_1 = pbuffer.data(idx_eri_1_pd + 12);

    auto g_z_xy_1 = pbuffer.data(idx_eri_1_pd + 13);

    auto g_z_xz_1 = pbuffer.data(idx_eri_1_pd + 14);

    auto g_z_yy_1 = pbuffer.data(idx_eri_1_pd + 15);

    auto g_z_yz_1 = pbuffer.data(idx_eri_1_pd + 16);

    auto g_z_zz_1 = pbuffer.data(idx_eri_1_pd + 17);

    // Set up components of auxiliary buffer : DP

    auto g_xx_x_1 = pbuffer.data(idx_eri_1_dp);

    auto g_xx_y_1 = pbuffer.data(idx_eri_1_dp + 1);

    auto g_xx_z_1 = pbuffer.data(idx_eri_1_dp + 2);

    auto g_yy_x_1 = pbuffer.data(idx_eri_1_dp + 9);

    auto g_yy_y_1 = pbuffer.data(idx_eri_1_dp + 10);

    auto g_yy_z_1 = pbuffer.data(idx_eri_1_dp + 11);

    auto g_zz_x_1 = pbuffer.data(idx_eri_1_dp + 15);

    auto g_zz_y_1 = pbuffer.data(idx_eri_1_dp + 16);

    auto g_zz_z_1 = pbuffer.data(idx_eri_1_dp + 17);

    // Set up components of auxiliary buffer : DD

    auto g_xx_xx_1 = pbuffer.data(idx_eri_1_dd);

    auto g_xx_xy_1 = pbuffer.data(idx_eri_1_dd + 1);

    auto g_xx_xz_1 = pbuffer.data(idx_eri_1_dd + 2);

    auto g_xx_yy_1 = pbuffer.data(idx_eri_1_dd + 3);

    auto g_xx_yz_1 = pbuffer.data(idx_eri_1_dd + 4);

    auto g_xx_zz_1 = pbuffer.data(idx_eri_1_dd + 5);

    auto g_xy_xy_1 = pbuffer.data(idx_eri_1_dd + 7);

    auto g_xz_xx_1 = pbuffer.data(idx_eri_1_dd + 12);

    auto g_xz_xz_1 = pbuffer.data(idx_eri_1_dd + 14);

    auto g_yy_xx_1 = pbuffer.data(idx_eri_1_dd + 18);

    auto g_yy_xy_1 = pbuffer.data(idx_eri_1_dd + 19);

    auto g_yy_xz_1 = pbuffer.data(idx_eri_1_dd + 20);

    auto g_yy_yy_1 = pbuffer.data(idx_eri_1_dd + 21);

    auto g_yy_yz_1 = pbuffer.data(idx_eri_1_dd + 22);

    auto g_yy_zz_1 = pbuffer.data(idx_eri_1_dd + 23);

    auto g_yz_yy_1 = pbuffer.data(idx_eri_1_dd + 27);

    auto g_yz_yz_1 = pbuffer.data(idx_eri_1_dd + 28);

    auto g_yz_zz_1 = pbuffer.data(idx_eri_1_dd + 29);

    auto g_zz_xx_1 = pbuffer.data(idx_eri_1_dd + 30);

    auto g_zz_xy_1 = pbuffer.data(idx_eri_1_dd + 31);

    auto g_zz_xz_1 = pbuffer.data(idx_eri_1_dd + 32);

    auto g_zz_yy_1 = pbuffer.data(idx_eri_1_dd + 33);

    auto g_zz_yz_1 = pbuffer.data(idx_eri_1_dd + 34);

    auto g_zz_zz_1 = pbuffer.data(idx_eri_1_dd + 35);

    // Set up 0-6 components of targeted buffer : FD

    auto g_xxx_xx_0 = pbuffer.data(idx_eri_0_fd);

    auto g_xxx_xy_0 = pbuffer.data(idx_eri_0_fd + 1);

    auto g_xxx_xz_0 = pbuffer.data(idx_eri_0_fd + 2);

    auto g_xxx_yy_0 = pbuffer.data(idx_eri_0_fd + 3);

    auto g_xxx_yz_0 = pbuffer.data(idx_eri_0_fd + 4);

    auto g_xxx_zz_0 = pbuffer.data(idx_eri_0_fd + 5);

    #pragma omp simd aligned(g_x_xx_0, g_x_xx_1, g_x_xy_0, g_x_xy_1, g_x_xz_0, g_x_xz_1, g_x_yy_0, g_x_yy_1, g_x_yz_0, g_x_yz_1, g_x_zz_0, g_x_zz_1, g_xx_x_1, g_xx_xx_1, g_xx_xy_1, g_xx_xz_1, g_xx_y_1, g_xx_yy_1, g_xx_yz_1, g_xx_z_1, g_xx_zz_1, g_xxx_xx_0, g_xxx_xy_0, g_xxx_xz_0, g_xxx_yy_0, g_xxx_yz_0, g_xxx_zz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxx_xx_0[i] = 2.0 * g_x_xx_0[i] * fbe_0 - 2.0 * g_x_xx_1[i] * fz_be_0 + 2.0 * g_xx_x_1[i] * fe_0 + g_xx_xx_1[i] * pa_x[i];

        g_xxx_xy_0[i] = 2.0 * g_x_xy_0[i] * fbe_0 - 2.0 * g_x_xy_1[i] * fz_be_0 + g_xx_y_1[i] * fe_0 + g_xx_xy_1[i] * pa_x[i];

        g_xxx_xz_0[i] = 2.0 * g_x_xz_0[i] * fbe_0 - 2.0 * g_x_xz_1[i] * fz_be_0 + g_xx_z_1[i] * fe_0 + g_xx_xz_1[i] * pa_x[i];

        g_xxx_yy_0[i] = 2.0 * g_x_yy_0[i] * fbe_0 - 2.0 * g_x_yy_1[i] * fz_be_0 + g_xx_yy_1[i] * pa_x[i];

        g_xxx_yz_0[i] = 2.0 * g_x_yz_0[i] * fbe_0 - 2.0 * g_x_yz_1[i] * fz_be_0 + g_xx_yz_1[i] * pa_x[i];

        g_xxx_zz_0[i] = 2.0 * g_x_zz_0[i] * fbe_0 - 2.0 * g_x_zz_1[i] * fz_be_0 + g_xx_zz_1[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : FD

    auto g_xxy_xx_0 = pbuffer.data(idx_eri_0_fd + 6);

    auto g_xxy_xy_0 = pbuffer.data(idx_eri_0_fd + 7);

    auto g_xxy_xz_0 = pbuffer.data(idx_eri_0_fd + 8);

    auto g_xxy_yy_0 = pbuffer.data(idx_eri_0_fd + 9);

    auto g_xxy_yz_0 = pbuffer.data(idx_eri_0_fd + 10);

    auto g_xxy_zz_0 = pbuffer.data(idx_eri_0_fd + 11);

    #pragma omp simd aligned(g_xx_x_1, g_xx_xx_1, g_xx_xy_1, g_xx_xz_1, g_xx_y_1, g_xx_yy_1, g_xx_yz_1, g_xx_z_1, g_xx_zz_1, g_xxy_xx_0, g_xxy_xy_0, g_xxy_xz_0, g_xxy_yy_0, g_xxy_yz_0, g_xxy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxy_xx_0[i] = g_xx_xx_1[i] * pa_y[i];

        g_xxy_xy_0[i] = g_xx_x_1[i] * fe_0 + g_xx_xy_1[i] * pa_y[i];

        g_xxy_xz_0[i] = g_xx_xz_1[i] * pa_y[i];

        g_xxy_yy_0[i] = 2.0 * g_xx_y_1[i] * fe_0 + g_xx_yy_1[i] * pa_y[i];

        g_xxy_yz_0[i] = g_xx_z_1[i] * fe_0 + g_xx_yz_1[i] * pa_y[i];

        g_xxy_zz_0[i] = g_xx_zz_1[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : FD

    auto g_xxz_xx_0 = pbuffer.data(idx_eri_0_fd + 12);

    auto g_xxz_xy_0 = pbuffer.data(idx_eri_0_fd + 13);

    auto g_xxz_xz_0 = pbuffer.data(idx_eri_0_fd + 14);

    auto g_xxz_yy_0 = pbuffer.data(idx_eri_0_fd + 15);

    auto g_xxz_yz_0 = pbuffer.data(idx_eri_0_fd + 16);

    auto g_xxz_zz_0 = pbuffer.data(idx_eri_0_fd + 17);

    #pragma omp simd aligned(g_xx_x_1, g_xx_xx_1, g_xx_xy_1, g_xx_xz_1, g_xx_y_1, g_xx_yy_1, g_xx_yz_1, g_xx_z_1, g_xx_zz_1, g_xxz_xx_0, g_xxz_xy_0, g_xxz_xz_0, g_xxz_yy_0, g_xxz_yz_0, g_xxz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxz_xx_0[i] = g_xx_xx_1[i] * pa_z[i];

        g_xxz_xy_0[i] = g_xx_xy_1[i] * pa_z[i];

        g_xxz_xz_0[i] = g_xx_x_1[i] * fe_0 + g_xx_xz_1[i] * pa_z[i];

        g_xxz_yy_0[i] = g_xx_yy_1[i] * pa_z[i];

        g_xxz_yz_0[i] = g_xx_y_1[i] * fe_0 + g_xx_yz_1[i] * pa_z[i];

        g_xxz_zz_0[i] = 2.0 * g_xx_z_1[i] * fe_0 + g_xx_zz_1[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : FD

    auto g_xyy_xx_0 = pbuffer.data(idx_eri_0_fd + 18);

    auto g_xyy_xy_0 = pbuffer.data(idx_eri_0_fd + 19);

    auto g_xyy_xz_0 = pbuffer.data(idx_eri_0_fd + 20);

    auto g_xyy_yy_0 = pbuffer.data(idx_eri_0_fd + 21);

    auto g_xyy_yz_0 = pbuffer.data(idx_eri_0_fd + 22);

    auto g_xyy_zz_0 = pbuffer.data(idx_eri_0_fd + 23);

    #pragma omp simd aligned(g_xyy_xx_0, g_xyy_xy_0, g_xyy_xz_0, g_xyy_yy_0, g_xyy_yz_0, g_xyy_zz_0, g_yy_x_1, g_yy_xx_1, g_yy_xy_1, g_yy_xz_1, g_yy_y_1, g_yy_yy_1, g_yy_yz_1, g_yy_z_1, g_yy_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyy_xx_0[i] = 2.0 * g_yy_x_1[i] * fe_0 + g_yy_xx_1[i] * pa_x[i];

        g_xyy_xy_0[i] = g_yy_y_1[i] * fe_0 + g_yy_xy_1[i] * pa_x[i];

        g_xyy_xz_0[i] = g_yy_z_1[i] * fe_0 + g_yy_xz_1[i] * pa_x[i];

        g_xyy_yy_0[i] = g_yy_yy_1[i] * pa_x[i];

        g_xyy_yz_0[i] = g_yy_yz_1[i] * pa_x[i];

        g_xyy_zz_0[i] = g_yy_zz_1[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : FD

    auto g_xyz_xx_0 = pbuffer.data(idx_eri_0_fd + 24);

    auto g_xyz_xy_0 = pbuffer.data(idx_eri_0_fd + 25);

    auto g_xyz_xz_0 = pbuffer.data(idx_eri_0_fd + 26);

    auto g_xyz_yy_0 = pbuffer.data(idx_eri_0_fd + 27);

    auto g_xyz_yz_0 = pbuffer.data(idx_eri_0_fd + 28);

    auto g_xyz_zz_0 = pbuffer.data(idx_eri_0_fd + 29);

    #pragma omp simd aligned(g_xy_xy_1, g_xyz_xx_0, g_xyz_xy_0, g_xyz_xz_0, g_xyz_yy_0, g_xyz_yz_0, g_xyz_zz_0, g_xz_xx_1, g_xz_xz_1, g_yz_yy_1, g_yz_yz_1, g_yz_zz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyz_xx_0[i] = g_xz_xx_1[i] * pa_y[i];

        g_xyz_xy_0[i] = g_xy_xy_1[i] * pa_z[i];

        g_xyz_xz_0[i] = g_xz_xz_1[i] * pa_y[i];

        g_xyz_yy_0[i] = g_yz_yy_1[i] * pa_x[i];

        g_xyz_yz_0[i] = g_yz_yz_1[i] * pa_x[i];

        g_xyz_zz_0[i] = g_yz_zz_1[i] * pa_x[i];
    }

    // Set up 30-36 components of targeted buffer : FD

    auto g_xzz_xx_0 = pbuffer.data(idx_eri_0_fd + 30);

    auto g_xzz_xy_0 = pbuffer.data(idx_eri_0_fd + 31);

    auto g_xzz_xz_0 = pbuffer.data(idx_eri_0_fd + 32);

    auto g_xzz_yy_0 = pbuffer.data(idx_eri_0_fd + 33);

    auto g_xzz_yz_0 = pbuffer.data(idx_eri_0_fd + 34);

    auto g_xzz_zz_0 = pbuffer.data(idx_eri_0_fd + 35);

    #pragma omp simd aligned(g_xzz_xx_0, g_xzz_xy_0, g_xzz_xz_0, g_xzz_yy_0, g_xzz_yz_0, g_xzz_zz_0, g_zz_x_1, g_zz_xx_1, g_zz_xy_1, g_zz_xz_1, g_zz_y_1, g_zz_yy_1, g_zz_yz_1, g_zz_z_1, g_zz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzz_xx_0[i] = 2.0 * g_zz_x_1[i] * fe_0 + g_zz_xx_1[i] * pa_x[i];

        g_xzz_xy_0[i] = g_zz_y_1[i] * fe_0 + g_zz_xy_1[i] * pa_x[i];

        g_xzz_xz_0[i] = g_zz_z_1[i] * fe_0 + g_zz_xz_1[i] * pa_x[i];

        g_xzz_yy_0[i] = g_zz_yy_1[i] * pa_x[i];

        g_xzz_yz_0[i] = g_zz_yz_1[i] * pa_x[i];

        g_xzz_zz_0[i] = g_zz_zz_1[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : FD

    auto g_yyy_xx_0 = pbuffer.data(idx_eri_0_fd + 36);

    auto g_yyy_xy_0 = pbuffer.data(idx_eri_0_fd + 37);

    auto g_yyy_xz_0 = pbuffer.data(idx_eri_0_fd + 38);

    auto g_yyy_yy_0 = pbuffer.data(idx_eri_0_fd + 39);

    auto g_yyy_yz_0 = pbuffer.data(idx_eri_0_fd + 40);

    auto g_yyy_zz_0 = pbuffer.data(idx_eri_0_fd + 41);

    #pragma omp simd aligned(g_y_xx_0, g_y_xx_1, g_y_xy_0, g_y_xy_1, g_y_xz_0, g_y_xz_1, g_y_yy_0, g_y_yy_1, g_y_yz_0, g_y_yz_1, g_y_zz_0, g_y_zz_1, g_yy_x_1, g_yy_xx_1, g_yy_xy_1, g_yy_xz_1, g_yy_y_1, g_yy_yy_1, g_yy_yz_1, g_yy_z_1, g_yy_zz_1, g_yyy_xx_0, g_yyy_xy_0, g_yyy_xz_0, g_yyy_yy_0, g_yyy_yz_0, g_yyy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyy_xx_0[i] = 2.0 * g_y_xx_0[i] * fbe_0 - 2.0 * g_y_xx_1[i] * fz_be_0 + g_yy_xx_1[i] * pa_y[i];

        g_yyy_xy_0[i] = 2.0 * g_y_xy_0[i] * fbe_0 - 2.0 * g_y_xy_1[i] * fz_be_0 + g_yy_x_1[i] * fe_0 + g_yy_xy_1[i] * pa_y[i];

        g_yyy_xz_0[i] = 2.0 * g_y_xz_0[i] * fbe_0 - 2.0 * g_y_xz_1[i] * fz_be_0 + g_yy_xz_1[i] * pa_y[i];

        g_yyy_yy_0[i] = 2.0 * g_y_yy_0[i] * fbe_0 - 2.0 * g_y_yy_1[i] * fz_be_0 + 2.0 * g_yy_y_1[i] * fe_0 + g_yy_yy_1[i] * pa_y[i];

        g_yyy_yz_0[i] = 2.0 * g_y_yz_0[i] * fbe_0 - 2.0 * g_y_yz_1[i] * fz_be_0 + g_yy_z_1[i] * fe_0 + g_yy_yz_1[i] * pa_y[i];

        g_yyy_zz_0[i] = 2.0 * g_y_zz_0[i] * fbe_0 - 2.0 * g_y_zz_1[i] * fz_be_0 + g_yy_zz_1[i] * pa_y[i];
    }

    // Set up 42-48 components of targeted buffer : FD

    auto g_yyz_xx_0 = pbuffer.data(idx_eri_0_fd + 42);

    auto g_yyz_xy_0 = pbuffer.data(idx_eri_0_fd + 43);

    auto g_yyz_xz_0 = pbuffer.data(idx_eri_0_fd + 44);

    auto g_yyz_yy_0 = pbuffer.data(idx_eri_0_fd + 45);

    auto g_yyz_yz_0 = pbuffer.data(idx_eri_0_fd + 46);

    auto g_yyz_zz_0 = pbuffer.data(idx_eri_0_fd + 47);

    #pragma omp simd aligned(g_yy_x_1, g_yy_xx_1, g_yy_xy_1, g_yy_xz_1, g_yy_y_1, g_yy_yy_1, g_yy_yz_1, g_yy_z_1, g_yy_zz_1, g_yyz_xx_0, g_yyz_xy_0, g_yyz_xz_0, g_yyz_yy_0, g_yyz_yz_0, g_yyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyz_xx_0[i] = g_yy_xx_1[i] * pa_z[i];

        g_yyz_xy_0[i] = g_yy_xy_1[i] * pa_z[i];

        g_yyz_xz_0[i] = g_yy_x_1[i] * fe_0 + g_yy_xz_1[i] * pa_z[i];

        g_yyz_yy_0[i] = g_yy_yy_1[i] * pa_z[i];

        g_yyz_yz_0[i] = g_yy_y_1[i] * fe_0 + g_yy_yz_1[i] * pa_z[i];

        g_yyz_zz_0[i] = 2.0 * g_yy_z_1[i] * fe_0 + g_yy_zz_1[i] * pa_z[i];
    }

    // Set up 48-54 components of targeted buffer : FD

    auto g_yzz_xx_0 = pbuffer.data(idx_eri_0_fd + 48);

    auto g_yzz_xy_0 = pbuffer.data(idx_eri_0_fd + 49);

    auto g_yzz_xz_0 = pbuffer.data(idx_eri_0_fd + 50);

    auto g_yzz_yy_0 = pbuffer.data(idx_eri_0_fd + 51);

    auto g_yzz_yz_0 = pbuffer.data(idx_eri_0_fd + 52);

    auto g_yzz_zz_0 = pbuffer.data(idx_eri_0_fd + 53);

    #pragma omp simd aligned(g_yzz_xx_0, g_yzz_xy_0, g_yzz_xz_0, g_yzz_yy_0, g_yzz_yz_0, g_yzz_zz_0, g_zz_x_1, g_zz_xx_1, g_zz_xy_1, g_zz_xz_1, g_zz_y_1, g_zz_yy_1, g_zz_yz_1, g_zz_z_1, g_zz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzz_xx_0[i] = g_zz_xx_1[i] * pa_y[i];

        g_yzz_xy_0[i] = g_zz_x_1[i] * fe_0 + g_zz_xy_1[i] * pa_y[i];

        g_yzz_xz_0[i] = g_zz_xz_1[i] * pa_y[i];

        g_yzz_yy_0[i] = 2.0 * g_zz_y_1[i] * fe_0 + g_zz_yy_1[i] * pa_y[i];

        g_yzz_yz_0[i] = g_zz_z_1[i] * fe_0 + g_zz_yz_1[i] * pa_y[i];

        g_yzz_zz_0[i] = g_zz_zz_1[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : FD

    auto g_zzz_xx_0 = pbuffer.data(idx_eri_0_fd + 54);

    auto g_zzz_xy_0 = pbuffer.data(idx_eri_0_fd + 55);

    auto g_zzz_xz_0 = pbuffer.data(idx_eri_0_fd + 56);

    auto g_zzz_yy_0 = pbuffer.data(idx_eri_0_fd + 57);

    auto g_zzz_yz_0 = pbuffer.data(idx_eri_0_fd + 58);

    auto g_zzz_zz_0 = pbuffer.data(idx_eri_0_fd + 59);

    #pragma omp simd aligned(g_z_xx_0, g_z_xx_1, g_z_xy_0, g_z_xy_1, g_z_xz_0, g_z_xz_1, g_z_yy_0, g_z_yy_1, g_z_yz_0, g_z_yz_1, g_z_zz_0, g_z_zz_1, g_zz_x_1, g_zz_xx_1, g_zz_xy_1, g_zz_xz_1, g_zz_y_1, g_zz_yy_1, g_zz_yz_1, g_zz_z_1, g_zz_zz_1, g_zzz_xx_0, g_zzz_xy_0, g_zzz_xz_0, g_zzz_yy_0, g_zzz_yz_0, g_zzz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzz_xx_0[i] = 2.0 * g_z_xx_0[i] * fbe_0 - 2.0 * g_z_xx_1[i] * fz_be_0 + g_zz_xx_1[i] * pa_z[i];

        g_zzz_xy_0[i] = 2.0 * g_z_xy_0[i] * fbe_0 - 2.0 * g_z_xy_1[i] * fz_be_0 + g_zz_xy_1[i] * pa_z[i];

        g_zzz_xz_0[i] = 2.0 * g_z_xz_0[i] * fbe_0 - 2.0 * g_z_xz_1[i] * fz_be_0 + g_zz_x_1[i] * fe_0 + g_zz_xz_1[i] * pa_z[i];

        g_zzz_yy_0[i] = 2.0 * g_z_yy_0[i] * fbe_0 - 2.0 * g_z_yy_1[i] * fz_be_0 + g_zz_yy_1[i] * pa_z[i];

        g_zzz_yz_0[i] = 2.0 * g_z_yz_0[i] * fbe_0 - 2.0 * g_z_yz_1[i] * fz_be_0 + g_zz_y_1[i] * fe_0 + g_zz_yz_1[i] * pa_z[i];

        g_zzz_zz_0[i] = 2.0 * g_z_zz_0[i] * fbe_0 - 2.0 * g_z_zz_1[i] * fz_be_0 + 2.0 * g_zz_z_1[i] * fe_0 + g_zz_zz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

