#include "TwoCenterElectronRepulsionPrimRecFS.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_fs(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_fs,
                                const size_t idx_eri_0_ps,
                                const size_t idx_eri_1_ps,
                                const size_t idx_eri_1_ds,
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

    // Set up components of auxiliary buffer : PS

    auto g_x_0_0 = pbuffer.data(idx_eri_0_ps);

    auto g_y_0_0 = pbuffer.data(idx_eri_0_ps + 1);

    auto g_z_0_0 = pbuffer.data(idx_eri_0_ps + 2);

    // Set up components of auxiliary buffer : PS

    auto g_x_0_1 = pbuffer.data(idx_eri_1_ps);

    auto g_y_0_1 = pbuffer.data(idx_eri_1_ps + 1);

    auto g_z_0_1 = pbuffer.data(idx_eri_1_ps + 2);

    // Set up components of auxiliary buffer : DS

    auto g_xx_0_1 = pbuffer.data(idx_eri_1_ds);

    auto g_yy_0_1 = pbuffer.data(idx_eri_1_ds + 3);

    auto g_yz_0_1 = pbuffer.data(idx_eri_1_ds + 4);

    auto g_zz_0_1 = pbuffer.data(idx_eri_1_ds + 5);

    // Set up components of targeted buffer : FS

    auto g_xxx_0_0 = pbuffer.data(idx_eri_0_fs);

    auto g_xxy_0_0 = pbuffer.data(idx_eri_0_fs + 1);

    auto g_xxz_0_0 = pbuffer.data(idx_eri_0_fs + 2);

    auto g_xyy_0_0 = pbuffer.data(idx_eri_0_fs + 3);

    auto g_xyz_0_0 = pbuffer.data(idx_eri_0_fs + 4);

    auto g_xzz_0_0 = pbuffer.data(idx_eri_0_fs + 5);

    auto g_yyy_0_0 = pbuffer.data(idx_eri_0_fs + 6);

    auto g_yyz_0_0 = pbuffer.data(idx_eri_0_fs + 7);

    auto g_yzz_0_0 = pbuffer.data(idx_eri_0_fs + 8);

    auto g_zzz_0_0 = pbuffer.data(idx_eri_0_fs + 9);

    #pragma omp simd aligned(g_x_0_0, g_x_0_1, g_xx_0_1, g_xxx_0_0, g_xxy_0_0, g_xxz_0_0, g_xyy_0_0, g_xyz_0_0, g_xzz_0_0, g_y_0_0, g_y_0_1, g_yy_0_1, g_yyy_0_0, g_yyz_0_0, g_yz_0_1, g_yzz_0_0, g_z_0_0, g_z_0_1, g_zz_0_1, g_zzz_0_0, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxx_0_0[i] = 2.0 * g_x_0_0[i] * fbe_0 - 2.0 * g_x_0_1[i] * fz_be_0 + g_xx_0_1[i] * pa_x[i];

        g_xxy_0_0[i] = g_xx_0_1[i] * pa_y[i];

        g_xxz_0_0[i] = g_xx_0_1[i] * pa_z[i];

        g_xyy_0_0[i] = g_yy_0_1[i] * pa_x[i];

        g_xyz_0_0[i] = g_yz_0_1[i] * pa_x[i];

        g_xzz_0_0[i] = g_zz_0_1[i] * pa_x[i];

        g_yyy_0_0[i] = 2.0 * g_y_0_0[i] * fbe_0 - 2.0 * g_y_0_1[i] * fz_be_0 + g_yy_0_1[i] * pa_y[i];

        g_yyz_0_0[i] = g_yy_0_1[i] * pa_z[i];

        g_yzz_0_0[i] = g_zz_0_1[i] * pa_y[i];

        g_zzz_0_0[i] = 2.0 * g_z_0_0[i] * fbe_0 - 2.0 * g_z_0_1[i] * fz_be_0 + g_zz_0_1[i] * pa_z[i];
    }
}

} // t2ceri namespace
