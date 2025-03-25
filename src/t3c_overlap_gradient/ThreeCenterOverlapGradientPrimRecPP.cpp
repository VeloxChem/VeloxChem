#include "ThreeCenterOverlapGradientPrimRecPP.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_pp(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_pp,
                              const size_t idx_sp,
                              const size_t idx_ps,
                              const size_t idx_pp,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_sp);

    auto ts_0_y = pbuffer.data(idx_sp + 1);

    auto ts_0_z = pbuffer.data(idx_sp + 2);

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ps);

    auto ts_y_0 = pbuffer.data(idx_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ps + 2);

    // Set up components of auxiliary buffer : PP

    auto ts_x_x = pbuffer.data(idx_pp);

    auto ts_x_y = pbuffer.data(idx_pp + 1);

    auto ts_x_z = pbuffer.data(idx_pp + 2);

    auto ts_y_x = pbuffer.data(idx_pp + 3);

    auto ts_y_y = pbuffer.data(idx_pp + 4);

    auto ts_y_z = pbuffer.data(idx_pp + 5);

    auto ts_z_x = pbuffer.data(idx_pp + 6);

    auto ts_z_y = pbuffer.data(idx_pp + 7);

    auto ts_z_z = pbuffer.data(idx_pp + 8);

    // Set up 0-3 components of targeted buffer : PP

    auto gs_x_x_x = pbuffer.data(idx_g_pp);

    auto gs_x_x_y = pbuffer.data(idx_g_pp + 1);

    auto gs_x_x_z = pbuffer.data(idx_g_pp + 2);

    #pragma omp simd aligned(gc_x, gs_x_x_x, gs_x_x_y, gs_x_x_z, ts_0_x, ts_0_y, ts_0_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_x_x[i] = 2.0 * ts_0_x[i] * gfe_0 * tce_0 + 2.0 * ts_x_0[i] * gfe_0 * tce_0 + 2.0 * ts_x_x[i] * gc_x[i] * tce_0;

        gs_x_x_y[i] = 2.0 * ts_0_y[i] * gfe_0 * tce_0 + 2.0 * ts_x_y[i] * gc_x[i] * tce_0;

        gs_x_x_z[i] = 2.0 * ts_0_z[i] * gfe_0 * tce_0 + 2.0 * ts_x_z[i] * gc_x[i] * tce_0;
    }

    // Set up 3-6 components of targeted buffer : PP

    auto gs_x_y_x = pbuffer.data(idx_g_pp + 3);

    auto gs_x_y_y = pbuffer.data(idx_g_pp + 4);

    auto gs_x_y_z = pbuffer.data(idx_g_pp + 5);

    #pragma omp simd aligned(gc_x, gs_x_y_x, gs_x_y_y, gs_x_y_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_y_x[i] = 2.0 * ts_y_0[i] * gfe_0 * tce_0 + 2.0 * ts_y_x[i] * gc_x[i] * tce_0;

        gs_x_y_y[i] = 2.0 * ts_y_y[i] * gc_x[i] * tce_0;

        gs_x_y_z[i] = 2.0 * ts_y_z[i] * gc_x[i] * tce_0;
    }

    // Set up 6-9 components of targeted buffer : PP

    auto gs_x_z_x = pbuffer.data(idx_g_pp + 6);

    auto gs_x_z_y = pbuffer.data(idx_g_pp + 7);

    auto gs_x_z_z = pbuffer.data(idx_g_pp + 8);

    #pragma omp simd aligned(gc_x, gs_x_z_x, gs_x_z_y, gs_x_z_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_z_x[i] = 2.0 * ts_z_0[i] * gfe_0 * tce_0 + 2.0 * ts_z_x[i] * gc_x[i] * tce_0;

        gs_x_z_y[i] = 2.0 * ts_z_y[i] * gc_x[i] * tce_0;

        gs_x_z_z[i] = 2.0 * ts_z_z[i] * gc_x[i] * tce_0;
    }

    // Set up 9-12 components of targeted buffer : PP

    auto gs_y_x_x = pbuffer.data(idx_g_pp + 9);

    auto gs_y_x_y = pbuffer.data(idx_g_pp + 10);

    auto gs_y_x_z = pbuffer.data(idx_g_pp + 11);

    #pragma omp simd aligned(gc_y, gs_y_x_x, gs_y_x_y, gs_y_x_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_x_x[i] = 2.0 * ts_x_x[i] * gc_y[i] * tce_0;

        gs_y_x_y[i] = 2.0 * ts_x_0[i] * gfe_0 * tce_0 + 2.0 * ts_x_y[i] * gc_y[i] * tce_0;

        gs_y_x_z[i] = 2.0 * ts_x_z[i] * gc_y[i] * tce_0;
    }

    // Set up 12-15 components of targeted buffer : PP

    auto gs_y_y_x = pbuffer.data(idx_g_pp + 12);

    auto gs_y_y_y = pbuffer.data(idx_g_pp + 13);

    auto gs_y_y_z = pbuffer.data(idx_g_pp + 14);

    #pragma omp simd aligned(gc_y, gs_y_y_x, gs_y_y_y, gs_y_y_z, ts_0_x, ts_0_y, ts_0_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_y_x[i] = 2.0 * ts_0_x[i] * gfe_0 * tce_0 + 2.0 * ts_y_x[i] * gc_y[i] * tce_0;

        gs_y_y_y[i] = 2.0 * ts_0_y[i] * gfe_0 * tce_0 + 2.0 * ts_y_0[i] * gfe_0 * tce_0 + 2.0 * ts_y_y[i] * gc_y[i] * tce_0;

        gs_y_y_z[i] = 2.0 * ts_0_z[i] * gfe_0 * tce_0 + 2.0 * ts_y_z[i] * gc_y[i] * tce_0;
    }

    // Set up 15-18 components of targeted buffer : PP

    auto gs_y_z_x = pbuffer.data(idx_g_pp + 15);

    auto gs_y_z_y = pbuffer.data(idx_g_pp + 16);

    auto gs_y_z_z = pbuffer.data(idx_g_pp + 17);

    #pragma omp simd aligned(gc_y, gs_y_z_x, gs_y_z_y, gs_y_z_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_z_x[i] = 2.0 * ts_z_x[i] * gc_y[i] * tce_0;

        gs_y_z_y[i] = 2.0 * ts_z_0[i] * gfe_0 * tce_0 + 2.0 * ts_z_y[i] * gc_y[i] * tce_0;

        gs_y_z_z[i] = 2.0 * ts_z_z[i] * gc_y[i] * tce_0;
    }

    // Set up 18-21 components of targeted buffer : PP

    auto gs_z_x_x = pbuffer.data(idx_g_pp + 18);

    auto gs_z_x_y = pbuffer.data(idx_g_pp + 19);

    auto gs_z_x_z = pbuffer.data(idx_g_pp + 20);

    #pragma omp simd aligned(gc_z, gs_z_x_x, gs_z_x_y, gs_z_x_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_x_x[i] = 2.0 * ts_x_x[i] * gc_z[i] * tce_0;

        gs_z_x_y[i] = 2.0 * ts_x_y[i] * gc_z[i] * tce_0;

        gs_z_x_z[i] = 2.0 * ts_x_0[i] * gfe_0 * tce_0 + 2.0 * ts_x_z[i] * gc_z[i] * tce_0;
    }

    // Set up 21-24 components of targeted buffer : PP

    auto gs_z_y_x = pbuffer.data(idx_g_pp + 21);

    auto gs_z_y_y = pbuffer.data(idx_g_pp + 22);

    auto gs_z_y_z = pbuffer.data(idx_g_pp + 23);

    #pragma omp simd aligned(gc_z, gs_z_y_x, gs_z_y_y, gs_z_y_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_y_x[i] = 2.0 * ts_y_x[i] * gc_z[i] * tce_0;

        gs_z_y_y[i] = 2.0 * ts_y_y[i] * gc_z[i] * tce_0;

        gs_z_y_z[i] = 2.0 * ts_y_0[i] * gfe_0 * tce_0 + 2.0 * ts_y_z[i] * gc_z[i] * tce_0;
    }

    // Set up 24-27 components of targeted buffer : PP

    auto gs_z_z_x = pbuffer.data(idx_g_pp + 24);

    auto gs_z_z_y = pbuffer.data(idx_g_pp + 25);

    auto gs_z_z_z = pbuffer.data(idx_g_pp + 26);

    #pragma omp simd aligned(gc_z, gs_z_z_x, gs_z_z_y, gs_z_z_z, ts_0_x, ts_0_y, ts_0_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_z_x[i] = 2.0 * ts_0_x[i] * gfe_0 * tce_0 + 2.0 * ts_z_x[i] * gc_z[i] * tce_0;

        gs_z_z_y[i] = 2.0 * ts_0_y[i] * gfe_0 * tce_0 + 2.0 * ts_z_y[i] * gc_z[i] * tce_0;

        gs_z_z_z[i] = 2.0 * ts_0_z[i] * gfe_0 * tce_0 + 2.0 * ts_z_0[i] * gfe_0 * tce_0 + 2.0 * ts_z_z[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

