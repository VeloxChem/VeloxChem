#include "TwoCenterElectronRepulsionPrimRecPP.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_pp(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_pp,
                                const size_t idx_eri_1_ss,
                                const size_t idx_eri_1_sp,
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

    // Set up components of auxiliary buffer : SS

    auto g_0_0_1 = pbuffer.data(idx_eri_1_ss);

    // Set up components of auxiliary buffer : SP

    auto g_0_x_1 = pbuffer.data(idx_eri_1_sp);

    auto g_0_y_1 = pbuffer.data(idx_eri_1_sp + 1);

    auto g_0_z_1 = pbuffer.data(idx_eri_1_sp + 2);

    // Set up 0-3 components of targeted buffer : PP

    auto g_x_x_0 = pbuffer.data(idx_eri_0_pp);

    auto g_x_y_0 = pbuffer.data(idx_eri_0_pp + 1);

    auto g_x_z_0 = pbuffer.data(idx_eri_0_pp + 2);

    #pragma omp simd aligned(g_0_0_1, g_0_x_1, g_0_y_1, g_0_z_1, g_x_x_0, g_x_y_0, g_x_z_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_x_x_0[i] = g_0_0_1[i] * fe_0 + g_0_x_1[i] * pa_x[i];

        g_x_y_0[i] = g_0_y_1[i] * pa_x[i];

        g_x_z_0[i] = g_0_z_1[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : PP

    auto g_y_x_0 = pbuffer.data(idx_eri_0_pp + 3);

    auto g_y_y_0 = pbuffer.data(idx_eri_0_pp + 4);

    auto g_y_z_0 = pbuffer.data(idx_eri_0_pp + 5);

    #pragma omp simd aligned(g_0_0_1, g_0_x_1, g_0_y_1, g_0_z_1, g_y_x_0, g_y_y_0, g_y_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_y_x_0[i] = g_0_x_1[i] * pa_y[i];

        g_y_y_0[i] = g_0_0_1[i] * fe_0 + g_0_y_1[i] * pa_y[i];

        g_y_z_0[i] = g_0_z_1[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : PP

    auto g_z_x_0 = pbuffer.data(idx_eri_0_pp + 6);

    auto g_z_y_0 = pbuffer.data(idx_eri_0_pp + 7);

    auto g_z_z_0 = pbuffer.data(idx_eri_0_pp + 8);

    #pragma omp simd aligned(g_0_0_1, g_0_x_1, g_0_y_1, g_0_z_1, g_z_x_0, g_z_y_0, g_z_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_z_x_0[i] = g_0_x_1[i] * pa_z[i];

        g_z_y_0[i] = g_0_y_1[i] * pa_z[i];

        g_z_z_0[i] = g_0_0_1[i] * fe_0 + g_0_z_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

