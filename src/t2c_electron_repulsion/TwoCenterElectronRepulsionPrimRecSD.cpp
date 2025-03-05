#include "TwoCenterElectronRepulsionPrimRecSD.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_sd(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_sd,
                                const size_t idx_eri_0_ss,
                                const size_t idx_eri_1_ss,
                                const size_t idx_eri_1_sp,
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

    // Set up components of auxiliary buffer : SS

    auto g_0_0_0 = pbuffer.data(idx_eri_0_ss);

    // Set up components of auxiliary buffer : SS

    auto g_0_0_1 = pbuffer.data(idx_eri_1_ss);

    // Set up components of auxiliary buffer : SP

    auto g_0_x_1 = pbuffer.data(idx_eri_1_sp);

    auto g_0_y_1 = pbuffer.data(idx_eri_1_sp + 1);

    auto g_0_z_1 = pbuffer.data(idx_eri_1_sp + 2);

    // Set up components of targeted buffer : SD

    auto g_0_xx_0 = pbuffer.data(idx_eri_0_sd);

    auto g_0_xy_0 = pbuffer.data(idx_eri_0_sd + 1);

    auto g_0_xz_0 = pbuffer.data(idx_eri_0_sd + 2);

    auto g_0_yy_0 = pbuffer.data(idx_eri_0_sd + 3);

    auto g_0_yz_0 = pbuffer.data(idx_eri_0_sd + 4);

    auto g_0_zz_0 = pbuffer.data(idx_eri_0_sd + 5);

    #pragma omp simd aligned(g_0_0_0, g_0_0_1, g_0_x_1, g_0_xx_0, g_0_xy_0, g_0_xz_0, g_0_y_1, g_0_yy_0, g_0_yz_0, g_0_z_1, g_0_zz_0, pb_x, pb_y, pb_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fke_0 = 0.5 / b_exps[i];

        const double fz_ke_0 = a_exp * fke_0 / (a_exp + b_exps[i]);

        g_0_xx_0[i] = g_0_0_0[i] * fke_0 - g_0_0_1[i] * fz_ke_0 + g_0_x_1[i] * pb_x[i];

        g_0_xy_0[i] = g_0_y_1[i] * pb_x[i];

        g_0_xz_0[i] = g_0_z_1[i] * pb_x[i];

        g_0_yy_0[i] = g_0_0_0[i] * fke_0 - g_0_0_1[i] * fz_ke_0 + g_0_y_1[i] * pb_y[i];

        g_0_yz_0[i] = g_0_z_1[i] * pb_y[i];

        g_0_zz_0[i] = g_0_0_0[i] * fke_0 - g_0_0_1[i] * fz_ke_0 + g_0_z_1[i] * pb_z[i];
    }
}

} // t2ceri namespace

