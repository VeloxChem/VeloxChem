#include "TwoCenterElectronRepulsionPrimRecSF.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_sf(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_sf,
                                const size_t idx_eri_0_sp,
                                const size_t idx_eri_1_sp,
                                const size_t idx_eri_1_sd,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpb,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SP

    auto g_0_x_0 = pbuffer.data(idx_eri_0_sp);

    auto g_0_y_0 = pbuffer.data(idx_eri_0_sp + 1);

    auto g_0_z_0 = pbuffer.data(idx_eri_0_sp + 2);

    // Set up components of auxiliary buffer : SP

    auto g_0_x_1 = pbuffer.data(idx_eri_1_sp);

    auto g_0_y_1 = pbuffer.data(idx_eri_1_sp + 1);

    auto g_0_z_1 = pbuffer.data(idx_eri_1_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto g_0_xx_1 = pbuffer.data(idx_eri_1_sd);

    auto g_0_yy_1 = pbuffer.data(idx_eri_1_sd + 3);

    auto g_0_yz_1 = pbuffer.data(idx_eri_1_sd + 4);

    auto g_0_zz_1 = pbuffer.data(idx_eri_1_sd + 5);

    // Set up components of targeted buffer : SF

    auto g_0_xxx_0 = pbuffer.data(idx_eri_0_sf);

    auto g_0_xxy_0 = pbuffer.data(idx_eri_0_sf + 1);

    auto g_0_xxz_0 = pbuffer.data(idx_eri_0_sf + 2);

    auto g_0_xyy_0 = pbuffer.data(idx_eri_0_sf + 3);

    auto g_0_xyz_0 = pbuffer.data(idx_eri_0_sf + 4);

    auto g_0_xzz_0 = pbuffer.data(idx_eri_0_sf + 5);

    auto g_0_yyy_0 = pbuffer.data(idx_eri_0_sf + 6);

    auto g_0_yyz_0 = pbuffer.data(idx_eri_0_sf + 7);

    auto g_0_yzz_0 = pbuffer.data(idx_eri_0_sf + 8);

    auto g_0_zzz_0 = pbuffer.data(idx_eri_0_sf + 9);

    #pragma omp simd aligned(g_0_x_0, g_0_x_1, g_0_xx_1, g_0_xxx_0, g_0_xxy_0, g_0_xxz_0, g_0_xyy_0, g_0_xyz_0, g_0_xzz_0, g_0_y_0, g_0_y_1, g_0_yy_1, g_0_yyy_0, g_0_yyz_0, g_0_yz_1, g_0_yzz_0, g_0_z_0, g_0_z_1, g_0_zz_1, g_0_zzz_0, pb_x, pb_y, pb_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fke_0 = 0.5 / b_exps[i];

        const double fz_ke_0 = a_exp * fke_0 / (a_exp + b_exps[i]);

        g_0_xxx_0[i] = 2.0 * g_0_x_0[i] * fke_0 - 2.0 * g_0_x_1[i] * fz_ke_0 + g_0_xx_1[i] * pb_x[i];

        g_0_xxy_0[i] = g_0_xx_1[i] * pb_y[i];

        g_0_xxz_0[i] = g_0_xx_1[i] * pb_z[i];

        g_0_xyy_0[i] = g_0_yy_1[i] * pb_x[i];

        g_0_xyz_0[i] = g_0_yz_1[i] * pb_x[i];

        g_0_xzz_0[i] = g_0_zz_1[i] * pb_x[i];

        g_0_yyy_0[i] = 2.0 * g_0_y_0[i] * fke_0 - 2.0 * g_0_y_1[i] * fz_ke_0 + g_0_yy_1[i] * pb_y[i];

        g_0_yyz_0[i] = g_0_yy_1[i] * pb_z[i];

        g_0_yzz_0[i] = g_0_zz_1[i] * pb_y[i];

        g_0_zzz_0[i] = 2.0 * g_0_z_0[i] * fke_0 - 2.0 * g_0_z_1[i] * fz_ke_0 + g_0_zz_1[i] * pb_z[i];
    }
}

} // t2ceri namespace

