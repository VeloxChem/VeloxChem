#include "TwoCenterElectronRepulsionPrimRecGD.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_gd(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_gd,
                                const size_t idx_eri_0_dd,
                                const size_t idx_eri_1_dd,
                                const size_t idx_eri_1_fp,
                                const size_t idx_eri_1_fd,
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

    // Set up components of auxiliary buffer : DD

    auto g_xx_xx_0 = pbuffer.data(idx_eri_0_dd);

    auto g_xx_xy_0 = pbuffer.data(idx_eri_0_dd + 1);

    auto g_xx_xz_0 = pbuffer.data(idx_eri_0_dd + 2);

    auto g_xx_yy_0 = pbuffer.data(idx_eri_0_dd + 3);

    auto g_xx_yz_0 = pbuffer.data(idx_eri_0_dd + 4);

    auto g_xx_zz_0 = pbuffer.data(idx_eri_0_dd + 5);

    auto g_yy_xx_0 = pbuffer.data(idx_eri_0_dd + 18);

    auto g_yy_xy_0 = pbuffer.data(idx_eri_0_dd + 19);

    auto g_yy_xz_0 = pbuffer.data(idx_eri_0_dd + 20);

    auto g_yy_yy_0 = pbuffer.data(idx_eri_0_dd + 21);

    auto g_yy_yz_0 = pbuffer.data(idx_eri_0_dd + 22);

    auto g_yy_zz_0 = pbuffer.data(idx_eri_0_dd + 23);

    auto g_zz_xx_0 = pbuffer.data(idx_eri_0_dd + 30);

    auto g_zz_xy_0 = pbuffer.data(idx_eri_0_dd + 31);

    auto g_zz_xz_0 = pbuffer.data(idx_eri_0_dd + 32);

    auto g_zz_yy_0 = pbuffer.data(idx_eri_0_dd + 33);

    auto g_zz_yz_0 = pbuffer.data(idx_eri_0_dd + 34);

    auto g_zz_zz_0 = pbuffer.data(idx_eri_0_dd + 35);

    // Set up components of auxiliary buffer : DD

    auto g_xx_xx_1 = pbuffer.data(idx_eri_1_dd);

    auto g_xx_xy_1 = pbuffer.data(idx_eri_1_dd + 1);

    auto g_xx_xz_1 = pbuffer.data(idx_eri_1_dd + 2);

    auto g_xx_yy_1 = pbuffer.data(idx_eri_1_dd + 3);

    auto g_xx_yz_1 = pbuffer.data(idx_eri_1_dd + 4);

    auto g_xx_zz_1 = pbuffer.data(idx_eri_1_dd + 5);

    auto g_yy_xx_1 = pbuffer.data(idx_eri_1_dd + 18);

    auto g_yy_xy_1 = pbuffer.data(idx_eri_1_dd + 19);

    auto g_yy_xz_1 = pbuffer.data(idx_eri_1_dd + 20);

    auto g_yy_yy_1 = pbuffer.data(idx_eri_1_dd + 21);

    auto g_yy_yz_1 = pbuffer.data(idx_eri_1_dd + 22);

    auto g_yy_zz_1 = pbuffer.data(idx_eri_1_dd + 23);

    auto g_zz_xx_1 = pbuffer.data(idx_eri_1_dd + 30);

    auto g_zz_xy_1 = pbuffer.data(idx_eri_1_dd + 31);

    auto g_zz_xz_1 = pbuffer.data(idx_eri_1_dd + 32);

    auto g_zz_yy_1 = pbuffer.data(idx_eri_1_dd + 33);

    auto g_zz_yz_1 = pbuffer.data(idx_eri_1_dd + 34);

    auto g_zz_zz_1 = pbuffer.data(idx_eri_1_dd + 35);

    // Set up components of auxiliary buffer : FP

    auto g_xxx_x_1 = pbuffer.data(idx_eri_1_fp);

    auto g_xxx_y_1 = pbuffer.data(idx_eri_1_fp + 1);

    auto g_xxx_z_1 = pbuffer.data(idx_eri_1_fp + 2);

    auto g_xxz_z_1 = pbuffer.data(idx_eri_1_fp + 8);

    auto g_xyy_y_1 = pbuffer.data(idx_eri_1_fp + 10);

    auto g_xzz_z_1 = pbuffer.data(idx_eri_1_fp + 17);

    auto g_yyy_x_1 = pbuffer.data(idx_eri_1_fp + 18);

    auto g_yyy_y_1 = pbuffer.data(idx_eri_1_fp + 19);

    auto g_yyy_z_1 = pbuffer.data(idx_eri_1_fp + 20);

    auto g_yyz_z_1 = pbuffer.data(idx_eri_1_fp + 23);

    auto g_yzz_y_1 = pbuffer.data(idx_eri_1_fp + 25);

    auto g_yzz_z_1 = pbuffer.data(idx_eri_1_fp + 26);

    auto g_zzz_x_1 = pbuffer.data(idx_eri_1_fp + 27);

    auto g_zzz_y_1 = pbuffer.data(idx_eri_1_fp + 28);

    auto g_zzz_z_1 = pbuffer.data(idx_eri_1_fp + 29);

    // Set up components of auxiliary buffer : FD

    auto g_xxx_xx_1 = pbuffer.data(idx_eri_1_fd);

    auto g_xxx_xy_1 = pbuffer.data(idx_eri_1_fd + 1);

    auto g_xxx_xz_1 = pbuffer.data(idx_eri_1_fd + 2);

    auto g_xxx_yy_1 = pbuffer.data(idx_eri_1_fd + 3);

    auto g_xxx_yz_1 = pbuffer.data(idx_eri_1_fd + 4);

    auto g_xxx_zz_1 = pbuffer.data(idx_eri_1_fd + 5);

    auto g_xxy_xx_1 = pbuffer.data(idx_eri_1_fd + 6);

    auto g_xxy_xy_1 = pbuffer.data(idx_eri_1_fd + 7);

    auto g_xxy_xz_1 = pbuffer.data(idx_eri_1_fd + 8);

    auto g_xxy_yy_1 = pbuffer.data(idx_eri_1_fd + 9);

    auto g_xxz_xx_1 = pbuffer.data(idx_eri_1_fd + 12);

    auto g_xxz_xy_1 = pbuffer.data(idx_eri_1_fd + 13);

    auto g_xxz_xz_1 = pbuffer.data(idx_eri_1_fd + 14);

    auto g_xxz_yz_1 = pbuffer.data(idx_eri_1_fd + 16);

    auto g_xxz_zz_1 = pbuffer.data(idx_eri_1_fd + 17);

    auto g_xyy_xx_1 = pbuffer.data(idx_eri_1_fd + 18);

    auto g_xyy_xy_1 = pbuffer.data(idx_eri_1_fd + 19);

    auto g_xyy_yy_1 = pbuffer.data(idx_eri_1_fd + 21);

    auto g_xyy_yz_1 = pbuffer.data(idx_eri_1_fd + 22);

    auto g_xyy_zz_1 = pbuffer.data(idx_eri_1_fd + 23);

    auto g_xzz_xx_1 = pbuffer.data(idx_eri_1_fd + 30);

    auto g_xzz_xz_1 = pbuffer.data(idx_eri_1_fd + 32);

    auto g_xzz_yy_1 = pbuffer.data(idx_eri_1_fd + 33);

    auto g_xzz_yz_1 = pbuffer.data(idx_eri_1_fd + 34);

    auto g_xzz_zz_1 = pbuffer.data(idx_eri_1_fd + 35);

    auto g_yyy_xx_1 = pbuffer.data(idx_eri_1_fd + 36);

    auto g_yyy_xy_1 = pbuffer.data(idx_eri_1_fd + 37);

    auto g_yyy_xz_1 = pbuffer.data(idx_eri_1_fd + 38);

    auto g_yyy_yy_1 = pbuffer.data(idx_eri_1_fd + 39);

    auto g_yyy_yz_1 = pbuffer.data(idx_eri_1_fd + 40);

    auto g_yyy_zz_1 = pbuffer.data(idx_eri_1_fd + 41);

    auto g_yyz_xy_1 = pbuffer.data(idx_eri_1_fd + 43);

    auto g_yyz_xz_1 = pbuffer.data(idx_eri_1_fd + 44);

    auto g_yyz_yy_1 = pbuffer.data(idx_eri_1_fd + 45);

    auto g_yyz_yz_1 = pbuffer.data(idx_eri_1_fd + 46);

    auto g_yyz_zz_1 = pbuffer.data(idx_eri_1_fd + 47);

    auto g_yzz_xx_1 = pbuffer.data(idx_eri_1_fd + 48);

    auto g_yzz_xy_1 = pbuffer.data(idx_eri_1_fd + 49);

    auto g_yzz_xz_1 = pbuffer.data(idx_eri_1_fd + 50);

    auto g_yzz_yy_1 = pbuffer.data(idx_eri_1_fd + 51);

    auto g_yzz_yz_1 = pbuffer.data(idx_eri_1_fd + 52);

    auto g_yzz_zz_1 = pbuffer.data(idx_eri_1_fd + 53);

    auto g_zzz_xx_1 = pbuffer.data(idx_eri_1_fd + 54);

    auto g_zzz_xy_1 = pbuffer.data(idx_eri_1_fd + 55);

    auto g_zzz_xz_1 = pbuffer.data(idx_eri_1_fd + 56);

    auto g_zzz_yy_1 = pbuffer.data(idx_eri_1_fd + 57);

    auto g_zzz_yz_1 = pbuffer.data(idx_eri_1_fd + 58);

    auto g_zzz_zz_1 = pbuffer.data(idx_eri_1_fd + 59);

    // Set up 0-6 components of targeted buffer : GD

    auto g_xxxx_xx_0 = pbuffer.data(idx_eri_0_gd);

    auto g_xxxx_xy_0 = pbuffer.data(idx_eri_0_gd + 1);

    auto g_xxxx_xz_0 = pbuffer.data(idx_eri_0_gd + 2);

    auto g_xxxx_yy_0 = pbuffer.data(idx_eri_0_gd + 3);

    auto g_xxxx_yz_0 = pbuffer.data(idx_eri_0_gd + 4);

    auto g_xxxx_zz_0 = pbuffer.data(idx_eri_0_gd + 5);

    #pragma omp simd aligned(g_xx_xx_0, g_xx_xx_1, g_xx_xy_0, g_xx_xy_1, g_xx_xz_0, g_xx_xz_1, g_xx_yy_0, g_xx_yy_1, g_xx_yz_0, g_xx_yz_1, g_xx_zz_0, g_xx_zz_1, g_xxx_x_1, g_xxx_xx_1, g_xxx_xy_1, g_xxx_xz_1, g_xxx_y_1, g_xxx_yy_1, g_xxx_yz_1, g_xxx_z_1, g_xxx_zz_1, g_xxxx_xx_0, g_xxxx_xy_0, g_xxxx_xz_0, g_xxxx_yy_0, g_xxxx_yz_0, g_xxxx_zz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxx_xx_0[i] = 3.0 * g_xx_xx_0[i] * fbe_0 - 3.0 * g_xx_xx_1[i] * fz_be_0 + 2.0 * g_xxx_x_1[i] * fe_0 + g_xxx_xx_1[i] * pa_x[i];

        g_xxxx_xy_0[i] = 3.0 * g_xx_xy_0[i] * fbe_0 - 3.0 * g_xx_xy_1[i] * fz_be_0 + g_xxx_y_1[i] * fe_0 + g_xxx_xy_1[i] * pa_x[i];

        g_xxxx_xz_0[i] = 3.0 * g_xx_xz_0[i] * fbe_0 - 3.0 * g_xx_xz_1[i] * fz_be_0 + g_xxx_z_1[i] * fe_0 + g_xxx_xz_1[i] * pa_x[i];

        g_xxxx_yy_0[i] = 3.0 * g_xx_yy_0[i] * fbe_0 - 3.0 * g_xx_yy_1[i] * fz_be_0 + g_xxx_yy_1[i] * pa_x[i];

        g_xxxx_yz_0[i] = 3.0 * g_xx_yz_0[i] * fbe_0 - 3.0 * g_xx_yz_1[i] * fz_be_0 + g_xxx_yz_1[i] * pa_x[i];

        g_xxxx_zz_0[i] = 3.0 * g_xx_zz_0[i] * fbe_0 - 3.0 * g_xx_zz_1[i] * fz_be_0 + g_xxx_zz_1[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : GD

    auto g_xxxy_xx_0 = pbuffer.data(idx_eri_0_gd + 6);

    auto g_xxxy_xy_0 = pbuffer.data(idx_eri_0_gd + 7);

    auto g_xxxy_xz_0 = pbuffer.data(idx_eri_0_gd + 8);

    auto g_xxxy_yy_0 = pbuffer.data(idx_eri_0_gd + 9);

    auto g_xxxy_yz_0 = pbuffer.data(idx_eri_0_gd + 10);

    auto g_xxxy_zz_0 = pbuffer.data(idx_eri_0_gd + 11);

    #pragma omp simd aligned(g_xxx_x_1, g_xxx_xx_1, g_xxx_xy_1, g_xxx_xz_1, g_xxx_y_1, g_xxx_yy_1, g_xxx_yz_1, g_xxx_z_1, g_xxx_zz_1, g_xxxy_xx_0, g_xxxy_xy_0, g_xxxy_xz_0, g_xxxy_yy_0, g_xxxy_yz_0, g_xxxy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxy_xx_0[i] = g_xxx_xx_1[i] * pa_y[i];

        g_xxxy_xy_0[i] = g_xxx_x_1[i] * fe_0 + g_xxx_xy_1[i] * pa_y[i];

        g_xxxy_xz_0[i] = g_xxx_xz_1[i] * pa_y[i];

        g_xxxy_yy_0[i] = 2.0 * g_xxx_y_1[i] * fe_0 + g_xxx_yy_1[i] * pa_y[i];

        g_xxxy_yz_0[i] = g_xxx_z_1[i] * fe_0 + g_xxx_yz_1[i] * pa_y[i];

        g_xxxy_zz_0[i] = g_xxx_zz_1[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : GD

    auto g_xxxz_xx_0 = pbuffer.data(idx_eri_0_gd + 12);

    auto g_xxxz_xy_0 = pbuffer.data(idx_eri_0_gd + 13);

    auto g_xxxz_xz_0 = pbuffer.data(idx_eri_0_gd + 14);

    auto g_xxxz_yy_0 = pbuffer.data(idx_eri_0_gd + 15);

    auto g_xxxz_yz_0 = pbuffer.data(idx_eri_0_gd + 16);

    auto g_xxxz_zz_0 = pbuffer.data(idx_eri_0_gd + 17);

    #pragma omp simd aligned(g_xxx_x_1, g_xxx_xx_1, g_xxx_xy_1, g_xxx_xz_1, g_xxx_y_1, g_xxx_yy_1, g_xxx_yz_1, g_xxx_z_1, g_xxx_zz_1, g_xxxz_xx_0, g_xxxz_xy_0, g_xxxz_xz_0, g_xxxz_yy_0, g_xxxz_yz_0, g_xxxz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxz_xx_0[i] = g_xxx_xx_1[i] * pa_z[i];

        g_xxxz_xy_0[i] = g_xxx_xy_1[i] * pa_z[i];

        g_xxxz_xz_0[i] = g_xxx_x_1[i] * fe_0 + g_xxx_xz_1[i] * pa_z[i];

        g_xxxz_yy_0[i] = g_xxx_yy_1[i] * pa_z[i];

        g_xxxz_yz_0[i] = g_xxx_y_1[i] * fe_0 + g_xxx_yz_1[i] * pa_z[i];

        g_xxxz_zz_0[i] = 2.0 * g_xxx_z_1[i] * fe_0 + g_xxx_zz_1[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : GD

    auto g_xxyy_xx_0 = pbuffer.data(idx_eri_0_gd + 18);

    auto g_xxyy_xy_0 = pbuffer.data(idx_eri_0_gd + 19);

    auto g_xxyy_xz_0 = pbuffer.data(idx_eri_0_gd + 20);

    auto g_xxyy_yy_0 = pbuffer.data(idx_eri_0_gd + 21);

    auto g_xxyy_yz_0 = pbuffer.data(idx_eri_0_gd + 22);

    auto g_xxyy_zz_0 = pbuffer.data(idx_eri_0_gd + 23);

    #pragma omp simd aligned(g_xx_xx_0, g_xx_xx_1, g_xx_xz_0, g_xx_xz_1, g_xxy_xx_1, g_xxy_xz_1, g_xxyy_xx_0, g_xxyy_xy_0, g_xxyy_xz_0, g_xxyy_yy_0, g_xxyy_yz_0, g_xxyy_zz_0, g_xyy_xy_1, g_xyy_y_1, g_xyy_yy_1, g_xyy_yz_1, g_xyy_zz_1, g_yy_xy_0, g_yy_xy_1, g_yy_yy_0, g_yy_yy_1, g_yy_yz_0, g_yy_yz_1, g_yy_zz_0, g_yy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyy_xx_0[i] = g_xx_xx_0[i] * fbe_0 - g_xx_xx_1[i] * fz_be_0 + g_xxy_xx_1[i] * pa_y[i];

        g_xxyy_xy_0[i] = g_yy_xy_0[i] * fbe_0 - g_yy_xy_1[i] * fz_be_0 + g_xyy_y_1[i] * fe_0 + g_xyy_xy_1[i] * pa_x[i];

        g_xxyy_xz_0[i] = g_xx_xz_0[i] * fbe_0 - g_xx_xz_1[i] * fz_be_0 + g_xxy_xz_1[i] * pa_y[i];

        g_xxyy_yy_0[i] = g_yy_yy_0[i] * fbe_0 - g_yy_yy_1[i] * fz_be_0 + g_xyy_yy_1[i] * pa_x[i];

        g_xxyy_yz_0[i] = g_yy_yz_0[i] * fbe_0 - g_yy_yz_1[i] * fz_be_0 + g_xyy_yz_1[i] * pa_x[i];

        g_xxyy_zz_0[i] = g_yy_zz_0[i] * fbe_0 - g_yy_zz_1[i] * fz_be_0 + g_xyy_zz_1[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : GD

    auto g_xxyz_xx_0 = pbuffer.data(idx_eri_0_gd + 24);

    auto g_xxyz_xy_0 = pbuffer.data(idx_eri_0_gd + 25);

    auto g_xxyz_xz_0 = pbuffer.data(idx_eri_0_gd + 26);

    auto g_xxyz_yy_0 = pbuffer.data(idx_eri_0_gd + 27);

    auto g_xxyz_yz_0 = pbuffer.data(idx_eri_0_gd + 28);

    auto g_xxyz_zz_0 = pbuffer.data(idx_eri_0_gd + 29);

    #pragma omp simd aligned(g_xxy_xy_1, g_xxy_yy_1, g_xxyz_xx_0, g_xxyz_xy_0, g_xxyz_xz_0, g_xxyz_yy_0, g_xxyz_yz_0, g_xxyz_zz_0, g_xxz_xx_1, g_xxz_xz_1, g_xxz_yz_1, g_xxz_z_1, g_xxz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyz_xx_0[i] = g_xxz_xx_1[i] * pa_y[i];

        g_xxyz_xy_0[i] = g_xxy_xy_1[i] * pa_z[i];

        g_xxyz_xz_0[i] = g_xxz_xz_1[i] * pa_y[i];

        g_xxyz_yy_0[i] = g_xxy_yy_1[i] * pa_z[i];

        g_xxyz_yz_0[i] = g_xxz_z_1[i] * fe_0 + g_xxz_yz_1[i] * pa_y[i];

        g_xxyz_zz_0[i] = g_xxz_zz_1[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : GD

    auto g_xxzz_xx_0 = pbuffer.data(idx_eri_0_gd + 30);

    auto g_xxzz_xy_0 = pbuffer.data(idx_eri_0_gd + 31);

    auto g_xxzz_xz_0 = pbuffer.data(idx_eri_0_gd + 32);

    auto g_xxzz_yy_0 = pbuffer.data(idx_eri_0_gd + 33);

    auto g_xxzz_yz_0 = pbuffer.data(idx_eri_0_gd + 34);

    auto g_xxzz_zz_0 = pbuffer.data(idx_eri_0_gd + 35);

    #pragma omp simd aligned(g_xx_xx_0, g_xx_xx_1, g_xx_xy_0, g_xx_xy_1, g_xxz_xx_1, g_xxz_xy_1, g_xxzz_xx_0, g_xxzz_xy_0, g_xxzz_xz_0, g_xxzz_yy_0, g_xxzz_yz_0, g_xxzz_zz_0, g_xzz_xz_1, g_xzz_yy_1, g_xzz_yz_1, g_xzz_z_1, g_xzz_zz_1, g_zz_xz_0, g_zz_xz_1, g_zz_yy_0, g_zz_yy_1, g_zz_yz_0, g_zz_yz_1, g_zz_zz_0, g_zz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzz_xx_0[i] = g_xx_xx_0[i] * fbe_0 - g_xx_xx_1[i] * fz_be_0 + g_xxz_xx_1[i] * pa_z[i];

        g_xxzz_xy_0[i] = g_xx_xy_0[i] * fbe_0 - g_xx_xy_1[i] * fz_be_0 + g_xxz_xy_1[i] * pa_z[i];

        g_xxzz_xz_0[i] = g_zz_xz_0[i] * fbe_0 - g_zz_xz_1[i] * fz_be_0 + g_xzz_z_1[i] * fe_0 + g_xzz_xz_1[i] * pa_x[i];

        g_xxzz_yy_0[i] = g_zz_yy_0[i] * fbe_0 - g_zz_yy_1[i] * fz_be_0 + g_xzz_yy_1[i] * pa_x[i];

        g_xxzz_yz_0[i] = g_zz_yz_0[i] * fbe_0 - g_zz_yz_1[i] * fz_be_0 + g_xzz_yz_1[i] * pa_x[i];

        g_xxzz_zz_0[i] = g_zz_zz_0[i] * fbe_0 - g_zz_zz_1[i] * fz_be_0 + g_xzz_zz_1[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : GD

    auto g_xyyy_xx_0 = pbuffer.data(idx_eri_0_gd + 36);

    auto g_xyyy_xy_0 = pbuffer.data(idx_eri_0_gd + 37);

    auto g_xyyy_xz_0 = pbuffer.data(idx_eri_0_gd + 38);

    auto g_xyyy_yy_0 = pbuffer.data(idx_eri_0_gd + 39);

    auto g_xyyy_yz_0 = pbuffer.data(idx_eri_0_gd + 40);

    auto g_xyyy_zz_0 = pbuffer.data(idx_eri_0_gd + 41);

    #pragma omp simd aligned(g_xyyy_xx_0, g_xyyy_xy_0, g_xyyy_xz_0, g_xyyy_yy_0, g_xyyy_yz_0, g_xyyy_zz_0, g_yyy_x_1, g_yyy_xx_1, g_yyy_xy_1, g_yyy_xz_1, g_yyy_y_1, g_yyy_yy_1, g_yyy_yz_1, g_yyy_z_1, g_yyy_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyy_xx_0[i] = 2.0 * g_yyy_x_1[i] * fe_0 + g_yyy_xx_1[i] * pa_x[i];

        g_xyyy_xy_0[i] = g_yyy_y_1[i] * fe_0 + g_yyy_xy_1[i] * pa_x[i];

        g_xyyy_xz_0[i] = g_yyy_z_1[i] * fe_0 + g_yyy_xz_1[i] * pa_x[i];

        g_xyyy_yy_0[i] = g_yyy_yy_1[i] * pa_x[i];

        g_xyyy_yz_0[i] = g_yyy_yz_1[i] * pa_x[i];

        g_xyyy_zz_0[i] = g_yyy_zz_1[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : GD

    auto g_xyyz_xx_0 = pbuffer.data(idx_eri_0_gd + 42);

    auto g_xyyz_xy_0 = pbuffer.data(idx_eri_0_gd + 43);

    auto g_xyyz_xz_0 = pbuffer.data(idx_eri_0_gd + 44);

    auto g_xyyz_yy_0 = pbuffer.data(idx_eri_0_gd + 45);

    auto g_xyyz_yz_0 = pbuffer.data(idx_eri_0_gd + 46);

    auto g_xyyz_zz_0 = pbuffer.data(idx_eri_0_gd + 47);

    #pragma omp simd aligned(g_xyy_xx_1, g_xyy_xy_1, g_xyyz_xx_0, g_xyyz_xy_0, g_xyyz_xz_0, g_xyyz_yy_0, g_xyyz_yz_0, g_xyyz_zz_0, g_yyz_xz_1, g_yyz_yy_1, g_yyz_yz_1, g_yyz_z_1, g_yyz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyz_xx_0[i] = g_xyy_xx_1[i] * pa_z[i];

        g_xyyz_xy_0[i] = g_xyy_xy_1[i] * pa_z[i];

        g_xyyz_xz_0[i] = g_yyz_z_1[i] * fe_0 + g_yyz_xz_1[i] * pa_x[i];

        g_xyyz_yy_0[i] = g_yyz_yy_1[i] * pa_x[i];

        g_xyyz_yz_0[i] = g_yyz_yz_1[i] * pa_x[i];

        g_xyyz_zz_0[i] = g_yyz_zz_1[i] * pa_x[i];
    }

    // Set up 48-54 components of targeted buffer : GD

    auto g_xyzz_xx_0 = pbuffer.data(idx_eri_0_gd + 48);

    auto g_xyzz_xy_0 = pbuffer.data(idx_eri_0_gd + 49);

    auto g_xyzz_xz_0 = pbuffer.data(idx_eri_0_gd + 50);

    auto g_xyzz_yy_0 = pbuffer.data(idx_eri_0_gd + 51);

    auto g_xyzz_yz_0 = pbuffer.data(idx_eri_0_gd + 52);

    auto g_xyzz_zz_0 = pbuffer.data(idx_eri_0_gd + 53);

    #pragma omp simd aligned(g_xyzz_xx_0, g_xyzz_xy_0, g_xyzz_xz_0, g_xyzz_yy_0, g_xyzz_yz_0, g_xyzz_zz_0, g_xzz_xx_1, g_xzz_xz_1, g_yzz_xy_1, g_yzz_y_1, g_yzz_yy_1, g_yzz_yz_1, g_yzz_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzz_xx_0[i] = g_xzz_xx_1[i] * pa_y[i];

        g_xyzz_xy_0[i] = g_yzz_y_1[i] * fe_0 + g_yzz_xy_1[i] * pa_x[i];

        g_xyzz_xz_0[i] = g_xzz_xz_1[i] * pa_y[i];

        g_xyzz_yy_0[i] = g_yzz_yy_1[i] * pa_x[i];

        g_xyzz_yz_0[i] = g_yzz_yz_1[i] * pa_x[i];

        g_xyzz_zz_0[i] = g_yzz_zz_1[i] * pa_x[i];
    }

    // Set up 54-60 components of targeted buffer : GD

    auto g_xzzz_xx_0 = pbuffer.data(idx_eri_0_gd + 54);

    auto g_xzzz_xy_0 = pbuffer.data(idx_eri_0_gd + 55);

    auto g_xzzz_xz_0 = pbuffer.data(idx_eri_0_gd + 56);

    auto g_xzzz_yy_0 = pbuffer.data(idx_eri_0_gd + 57);

    auto g_xzzz_yz_0 = pbuffer.data(idx_eri_0_gd + 58);

    auto g_xzzz_zz_0 = pbuffer.data(idx_eri_0_gd + 59);

    #pragma omp simd aligned(g_xzzz_xx_0, g_xzzz_xy_0, g_xzzz_xz_0, g_xzzz_yy_0, g_xzzz_yz_0, g_xzzz_zz_0, g_zzz_x_1, g_zzz_xx_1, g_zzz_xy_1, g_zzz_xz_1, g_zzz_y_1, g_zzz_yy_1, g_zzz_yz_1, g_zzz_z_1, g_zzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzz_xx_0[i] = 2.0 * g_zzz_x_1[i] * fe_0 + g_zzz_xx_1[i] * pa_x[i];

        g_xzzz_xy_0[i] = g_zzz_y_1[i] * fe_0 + g_zzz_xy_1[i] * pa_x[i];

        g_xzzz_xz_0[i] = g_zzz_z_1[i] * fe_0 + g_zzz_xz_1[i] * pa_x[i];

        g_xzzz_yy_0[i] = g_zzz_yy_1[i] * pa_x[i];

        g_xzzz_yz_0[i] = g_zzz_yz_1[i] * pa_x[i];

        g_xzzz_zz_0[i] = g_zzz_zz_1[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : GD

    auto g_yyyy_xx_0 = pbuffer.data(idx_eri_0_gd + 60);

    auto g_yyyy_xy_0 = pbuffer.data(idx_eri_0_gd + 61);

    auto g_yyyy_xz_0 = pbuffer.data(idx_eri_0_gd + 62);

    auto g_yyyy_yy_0 = pbuffer.data(idx_eri_0_gd + 63);

    auto g_yyyy_yz_0 = pbuffer.data(idx_eri_0_gd + 64);

    auto g_yyyy_zz_0 = pbuffer.data(idx_eri_0_gd + 65);

    #pragma omp simd aligned(g_yy_xx_0, g_yy_xx_1, g_yy_xy_0, g_yy_xy_1, g_yy_xz_0, g_yy_xz_1, g_yy_yy_0, g_yy_yy_1, g_yy_yz_0, g_yy_yz_1, g_yy_zz_0, g_yy_zz_1, g_yyy_x_1, g_yyy_xx_1, g_yyy_xy_1, g_yyy_xz_1, g_yyy_y_1, g_yyy_yy_1, g_yyy_yz_1, g_yyy_z_1, g_yyy_zz_1, g_yyyy_xx_0, g_yyyy_xy_0, g_yyyy_xz_0, g_yyyy_yy_0, g_yyyy_yz_0, g_yyyy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyy_xx_0[i] = 3.0 * g_yy_xx_0[i] * fbe_0 - 3.0 * g_yy_xx_1[i] * fz_be_0 + g_yyy_xx_1[i] * pa_y[i];

        g_yyyy_xy_0[i] = 3.0 * g_yy_xy_0[i] * fbe_0 - 3.0 * g_yy_xy_1[i] * fz_be_0 + g_yyy_x_1[i] * fe_0 + g_yyy_xy_1[i] * pa_y[i];

        g_yyyy_xz_0[i] = 3.0 * g_yy_xz_0[i] * fbe_0 - 3.0 * g_yy_xz_1[i] * fz_be_0 + g_yyy_xz_1[i] * pa_y[i];

        g_yyyy_yy_0[i] = 3.0 * g_yy_yy_0[i] * fbe_0 - 3.0 * g_yy_yy_1[i] * fz_be_0 + 2.0 * g_yyy_y_1[i] * fe_0 + g_yyy_yy_1[i] * pa_y[i];

        g_yyyy_yz_0[i] = 3.0 * g_yy_yz_0[i] * fbe_0 - 3.0 * g_yy_yz_1[i] * fz_be_0 + g_yyy_z_1[i] * fe_0 + g_yyy_yz_1[i] * pa_y[i];

        g_yyyy_zz_0[i] = 3.0 * g_yy_zz_0[i] * fbe_0 - 3.0 * g_yy_zz_1[i] * fz_be_0 + g_yyy_zz_1[i] * pa_y[i];
    }

    // Set up 66-72 components of targeted buffer : GD

    auto g_yyyz_xx_0 = pbuffer.data(idx_eri_0_gd + 66);

    auto g_yyyz_xy_0 = pbuffer.data(idx_eri_0_gd + 67);

    auto g_yyyz_xz_0 = pbuffer.data(idx_eri_0_gd + 68);

    auto g_yyyz_yy_0 = pbuffer.data(idx_eri_0_gd + 69);

    auto g_yyyz_yz_0 = pbuffer.data(idx_eri_0_gd + 70);

    auto g_yyyz_zz_0 = pbuffer.data(idx_eri_0_gd + 71);

    #pragma omp simd aligned(g_yyy_x_1, g_yyy_xx_1, g_yyy_xy_1, g_yyy_xz_1, g_yyy_y_1, g_yyy_yy_1, g_yyy_yz_1, g_yyy_z_1, g_yyy_zz_1, g_yyyz_xx_0, g_yyyz_xy_0, g_yyyz_xz_0, g_yyyz_yy_0, g_yyyz_yz_0, g_yyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyz_xx_0[i] = g_yyy_xx_1[i] * pa_z[i];

        g_yyyz_xy_0[i] = g_yyy_xy_1[i] * pa_z[i];

        g_yyyz_xz_0[i] = g_yyy_x_1[i] * fe_0 + g_yyy_xz_1[i] * pa_z[i];

        g_yyyz_yy_0[i] = g_yyy_yy_1[i] * pa_z[i];

        g_yyyz_yz_0[i] = g_yyy_y_1[i] * fe_0 + g_yyy_yz_1[i] * pa_z[i];

        g_yyyz_zz_0[i] = 2.0 * g_yyy_z_1[i] * fe_0 + g_yyy_zz_1[i] * pa_z[i];
    }

    // Set up 72-78 components of targeted buffer : GD

    auto g_yyzz_xx_0 = pbuffer.data(idx_eri_0_gd + 72);

    auto g_yyzz_xy_0 = pbuffer.data(idx_eri_0_gd + 73);

    auto g_yyzz_xz_0 = pbuffer.data(idx_eri_0_gd + 74);

    auto g_yyzz_yy_0 = pbuffer.data(idx_eri_0_gd + 75);

    auto g_yyzz_yz_0 = pbuffer.data(idx_eri_0_gd + 76);

    auto g_yyzz_zz_0 = pbuffer.data(idx_eri_0_gd + 77);

    #pragma omp simd aligned(g_yy_xy_0, g_yy_xy_1, g_yy_yy_0, g_yy_yy_1, g_yyz_xy_1, g_yyz_yy_1, g_yyzz_xx_0, g_yyzz_xy_0, g_yyzz_xz_0, g_yyzz_yy_0, g_yyzz_yz_0, g_yyzz_zz_0, g_yzz_xx_1, g_yzz_xz_1, g_yzz_yz_1, g_yzz_z_1, g_yzz_zz_1, g_zz_xx_0, g_zz_xx_1, g_zz_xz_0, g_zz_xz_1, g_zz_yz_0, g_zz_yz_1, g_zz_zz_0, g_zz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzz_xx_0[i] = g_zz_xx_0[i] * fbe_0 - g_zz_xx_1[i] * fz_be_0 + g_yzz_xx_1[i] * pa_y[i];

        g_yyzz_xy_0[i] = g_yy_xy_0[i] * fbe_0 - g_yy_xy_1[i] * fz_be_0 + g_yyz_xy_1[i] * pa_z[i];

        g_yyzz_xz_0[i] = g_zz_xz_0[i] * fbe_0 - g_zz_xz_1[i] * fz_be_0 + g_yzz_xz_1[i] * pa_y[i];

        g_yyzz_yy_0[i] = g_yy_yy_0[i] * fbe_0 - g_yy_yy_1[i] * fz_be_0 + g_yyz_yy_1[i] * pa_z[i];

        g_yyzz_yz_0[i] = g_zz_yz_0[i] * fbe_0 - g_zz_yz_1[i] * fz_be_0 + g_yzz_z_1[i] * fe_0 + g_yzz_yz_1[i] * pa_y[i];

        g_yyzz_zz_0[i] = g_zz_zz_0[i] * fbe_0 - g_zz_zz_1[i] * fz_be_0 + g_yzz_zz_1[i] * pa_y[i];
    }

    // Set up 78-84 components of targeted buffer : GD

    auto g_yzzz_xx_0 = pbuffer.data(idx_eri_0_gd + 78);

    auto g_yzzz_xy_0 = pbuffer.data(idx_eri_0_gd + 79);

    auto g_yzzz_xz_0 = pbuffer.data(idx_eri_0_gd + 80);

    auto g_yzzz_yy_0 = pbuffer.data(idx_eri_0_gd + 81);

    auto g_yzzz_yz_0 = pbuffer.data(idx_eri_0_gd + 82);

    auto g_yzzz_zz_0 = pbuffer.data(idx_eri_0_gd + 83);

    #pragma omp simd aligned(g_yzzz_xx_0, g_yzzz_xy_0, g_yzzz_xz_0, g_yzzz_yy_0, g_yzzz_yz_0, g_yzzz_zz_0, g_zzz_x_1, g_zzz_xx_1, g_zzz_xy_1, g_zzz_xz_1, g_zzz_y_1, g_zzz_yy_1, g_zzz_yz_1, g_zzz_z_1, g_zzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzz_xx_0[i] = g_zzz_xx_1[i] * pa_y[i];

        g_yzzz_xy_0[i] = g_zzz_x_1[i] * fe_0 + g_zzz_xy_1[i] * pa_y[i];

        g_yzzz_xz_0[i] = g_zzz_xz_1[i] * pa_y[i];

        g_yzzz_yy_0[i] = 2.0 * g_zzz_y_1[i] * fe_0 + g_zzz_yy_1[i] * pa_y[i];

        g_yzzz_yz_0[i] = g_zzz_z_1[i] * fe_0 + g_zzz_yz_1[i] * pa_y[i];

        g_yzzz_zz_0[i] = g_zzz_zz_1[i] * pa_y[i];
    }

    // Set up 84-90 components of targeted buffer : GD

    auto g_zzzz_xx_0 = pbuffer.data(idx_eri_0_gd + 84);

    auto g_zzzz_xy_0 = pbuffer.data(idx_eri_0_gd + 85);

    auto g_zzzz_xz_0 = pbuffer.data(idx_eri_0_gd + 86);

    auto g_zzzz_yy_0 = pbuffer.data(idx_eri_0_gd + 87);

    auto g_zzzz_yz_0 = pbuffer.data(idx_eri_0_gd + 88);

    auto g_zzzz_zz_0 = pbuffer.data(idx_eri_0_gd + 89);

    #pragma omp simd aligned(g_zz_xx_0, g_zz_xx_1, g_zz_xy_0, g_zz_xy_1, g_zz_xz_0, g_zz_xz_1, g_zz_yy_0, g_zz_yy_1, g_zz_yz_0, g_zz_yz_1, g_zz_zz_0, g_zz_zz_1, g_zzz_x_1, g_zzz_xx_1, g_zzz_xy_1, g_zzz_xz_1, g_zzz_y_1, g_zzz_yy_1, g_zzz_yz_1, g_zzz_z_1, g_zzz_zz_1, g_zzzz_xx_0, g_zzzz_xy_0, g_zzzz_xz_0, g_zzzz_yy_0, g_zzzz_yz_0, g_zzzz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzz_xx_0[i] = 3.0 * g_zz_xx_0[i] * fbe_0 - 3.0 * g_zz_xx_1[i] * fz_be_0 + g_zzz_xx_1[i] * pa_z[i];

        g_zzzz_xy_0[i] = 3.0 * g_zz_xy_0[i] * fbe_0 - 3.0 * g_zz_xy_1[i] * fz_be_0 + g_zzz_xy_1[i] * pa_z[i];

        g_zzzz_xz_0[i] = 3.0 * g_zz_xz_0[i] * fbe_0 - 3.0 * g_zz_xz_1[i] * fz_be_0 + g_zzz_x_1[i] * fe_0 + g_zzz_xz_1[i] * pa_z[i];

        g_zzzz_yy_0[i] = 3.0 * g_zz_yy_0[i] * fbe_0 - 3.0 * g_zz_yy_1[i] * fz_be_0 + g_zzz_yy_1[i] * pa_z[i];

        g_zzzz_yz_0[i] = 3.0 * g_zz_yz_0[i] * fbe_0 - 3.0 * g_zz_yz_1[i] * fz_be_0 + g_zzz_y_1[i] * fe_0 + g_zzz_yz_1[i] * pa_z[i];

        g_zzzz_zz_0[i] = 3.0 * g_zz_zz_0[i] * fbe_0 - 3.0 * g_zz_zz_1[i] * fz_be_0 + 2.0 * g_zzz_z_1[i] * fe_0 + g_zzz_zz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

