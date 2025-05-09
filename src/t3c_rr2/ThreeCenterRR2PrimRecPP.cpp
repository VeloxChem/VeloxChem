#include "ThreeCenterRR2PrimRecPP.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_pp(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_pp,
                  const size_t idx_sp,
                  const size_t idx_g_sp,
                  const size_t idx_ps,
                  const size_t idx_g_ps,
                  const size_t idx_pp,
                  const size_t idx_g_pp,
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

    // Set up components of auxiliary buffer : SP

    auto gr_0_x = pbuffer.data(idx_g_sp);

    auto gr_0_y = pbuffer.data(idx_g_sp + 1);

    auto gr_0_z = pbuffer.data(idx_g_sp + 2);

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ps);

    auto ts_y_0 = pbuffer.data(idx_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ps + 2);

    // Set up components of auxiliary buffer : PS

    auto gr_x_0 = pbuffer.data(idx_g_ps);

    auto gr_y_0 = pbuffer.data(idx_g_ps + 1);

    auto gr_z_0 = pbuffer.data(idx_g_ps + 2);

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

    // Set up components of auxiliary buffer : PP

    auto gr_x_x = pbuffer.data(idx_g_pp);

    auto gr_x_y = pbuffer.data(idx_g_pp + 1);

    auto gr_x_z = pbuffer.data(idx_g_pp + 2);

    auto gr_y_x = pbuffer.data(idx_g_pp + 3);

    auto gr_y_y = pbuffer.data(idx_g_pp + 4);

    auto gr_y_z = pbuffer.data(idx_g_pp + 5);

    auto gr_z_x = pbuffer.data(idx_g_pp + 6);

    auto gr_z_y = pbuffer.data(idx_g_pp + 7);

    auto gr_z_z = pbuffer.data(idx_g_pp + 8);

    // Set up 0-3 components of targeted buffer : PP

    auto grr_x_x_x = pbuffer.data(idx_gr_pp);

    auto grr_x_x_y = pbuffer.data(idx_gr_pp + 1);

    auto grr_x_x_z = pbuffer.data(idx_gr_pp + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_x, gr_0_y, gr_0_z, gr_x_0, gr_x_x, gr_x_y, gr_x_z, grr_x_x_x, grr_x_x_y, grr_x_x_z, ts_0_x, ts_0_y, ts_0_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_x_x[i] = ts_0_x[i] * gfe2_0 + gr_0_x[i] * gfe_0 + ts_x_0[i] * gfe2_0 + gr_x_0[i] * gfe_0 + ts_x_x[i] * gfe_0 * gc_x[i] + gr_x_x[i] * gc_x[i];

        grr_x_x_y[i] = ts_0_y[i] * gfe2_0 + gr_0_y[i] * gfe_0 + ts_x_y[i] * gfe_0 * gc_x[i] + gr_x_y[i] * gc_x[i];

        grr_x_x_z[i] = ts_0_z[i] * gfe2_0 + gr_0_z[i] * gfe_0 + ts_x_z[i] * gfe_0 * gc_x[i] + gr_x_z[i] * gc_x[i];
    }

    // Set up 3-6 components of targeted buffer : PP

    auto grr_x_y_x = pbuffer.data(idx_gr_pp + 3);

    auto grr_x_y_y = pbuffer.data(idx_gr_pp + 4);

    auto grr_x_y_z = pbuffer.data(idx_gr_pp + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_0, gr_y_x, gr_y_y, gr_y_z, grr_x_y_x, grr_x_y_y, grr_x_y_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_y_x[i] = ts_y_0[i] * gfe2_0 + gr_y_0[i] * gfe_0 + ts_y_x[i] * gfe_0 * gc_x[i] + gr_y_x[i] * gc_x[i];

        grr_x_y_y[i] = ts_y_y[i] * gfe_0 * gc_x[i] + gr_y_y[i] * gc_x[i];

        grr_x_y_z[i] = ts_y_z[i] * gfe_0 * gc_x[i] + gr_y_z[i] * gc_x[i];
    }

    // Set up 6-9 components of targeted buffer : PP

    auto grr_x_z_x = pbuffer.data(idx_gr_pp + 6);

    auto grr_x_z_y = pbuffer.data(idx_gr_pp + 7);

    auto grr_x_z_z = pbuffer.data(idx_gr_pp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_0, gr_z_x, gr_z_y, gr_z_z, grr_x_z_x, grr_x_z_y, grr_x_z_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_z_x[i] = ts_z_0[i] * gfe2_0 + gr_z_0[i] * gfe_0 + ts_z_x[i] * gfe_0 * gc_x[i] + gr_z_x[i] * gc_x[i];

        grr_x_z_y[i] = ts_z_y[i] * gfe_0 * gc_x[i] + gr_z_y[i] * gc_x[i];

        grr_x_z_z[i] = ts_z_z[i] * gfe_0 * gc_x[i] + gr_z_z[i] * gc_x[i];
    }

    // Set up 9-12 components of targeted buffer : PP

    auto grr_y_x_x = pbuffer.data(idx_gr_pp + 9);

    auto grr_y_x_y = pbuffer.data(idx_gr_pp + 10);

    auto grr_y_x_z = pbuffer.data(idx_gr_pp + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_0, gr_x_x, gr_x_y, gr_x_z, grr_y_x_x, grr_y_x_y, grr_y_x_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_x_x[i] = ts_x_x[i] * gfe_0 * gc_y[i] + gr_x_x[i] * gc_y[i];

        grr_y_x_y[i] = ts_x_0[i] * gfe2_0 + gr_x_0[i] * gfe_0 + ts_x_y[i] * gfe_0 * gc_y[i] + gr_x_y[i] * gc_y[i];

        grr_y_x_z[i] = ts_x_z[i] * gfe_0 * gc_y[i] + gr_x_z[i] * gc_y[i];
    }

    // Set up 12-15 components of targeted buffer : PP

    auto grr_y_y_x = pbuffer.data(idx_gr_pp + 12);

    auto grr_y_y_y = pbuffer.data(idx_gr_pp + 13);

    auto grr_y_y_z = pbuffer.data(idx_gr_pp + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_x, gr_0_y, gr_0_z, gr_y_0, gr_y_x, gr_y_y, gr_y_z, grr_y_y_x, grr_y_y_y, grr_y_y_z, ts_0_x, ts_0_y, ts_0_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_y_x[i] = ts_0_x[i] * gfe2_0 + gr_0_x[i] * gfe_0 + ts_y_x[i] * gfe_0 * gc_y[i] + gr_y_x[i] * gc_y[i];

        grr_y_y_y[i] = ts_0_y[i] * gfe2_0 + gr_0_y[i] * gfe_0 + ts_y_0[i] * gfe2_0 + gr_y_0[i] * gfe_0 + ts_y_y[i] * gfe_0 * gc_y[i] + gr_y_y[i] * gc_y[i];

        grr_y_y_z[i] = ts_0_z[i] * gfe2_0 + gr_0_z[i] * gfe_0 + ts_y_z[i] * gfe_0 * gc_y[i] + gr_y_z[i] * gc_y[i];
    }

    // Set up 15-18 components of targeted buffer : PP

    auto grr_y_z_x = pbuffer.data(idx_gr_pp + 15);

    auto grr_y_z_y = pbuffer.data(idx_gr_pp + 16);

    auto grr_y_z_z = pbuffer.data(idx_gr_pp + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_0, gr_z_x, gr_z_y, gr_z_z, grr_y_z_x, grr_y_z_y, grr_y_z_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_z_x[i] = ts_z_x[i] * gfe_0 * gc_y[i] + gr_z_x[i] * gc_y[i];

        grr_y_z_y[i] = ts_z_0[i] * gfe2_0 + gr_z_0[i] * gfe_0 + ts_z_y[i] * gfe_0 * gc_y[i] + gr_z_y[i] * gc_y[i];

        grr_y_z_z[i] = ts_z_z[i] * gfe_0 * gc_y[i] + gr_z_z[i] * gc_y[i];
    }

    // Set up 18-21 components of targeted buffer : PP

    auto grr_z_x_x = pbuffer.data(idx_gr_pp + 18);

    auto grr_z_x_y = pbuffer.data(idx_gr_pp + 19);

    auto grr_z_x_z = pbuffer.data(idx_gr_pp + 20);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_0, gr_x_x, gr_x_y, gr_x_z, grr_z_x_x, grr_z_x_y, grr_z_x_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_x_x[i] = ts_x_x[i] * gfe_0 * gc_z[i] + gr_x_x[i] * gc_z[i];

        grr_z_x_y[i] = ts_x_y[i] * gfe_0 * gc_z[i] + gr_x_y[i] * gc_z[i];

        grr_z_x_z[i] = ts_x_0[i] * gfe2_0 + gr_x_0[i] * gfe_0 + ts_x_z[i] * gfe_0 * gc_z[i] + gr_x_z[i] * gc_z[i];
    }

    // Set up 21-24 components of targeted buffer : PP

    auto grr_z_y_x = pbuffer.data(idx_gr_pp + 21);

    auto grr_z_y_y = pbuffer.data(idx_gr_pp + 22);

    auto grr_z_y_z = pbuffer.data(idx_gr_pp + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_0, gr_y_x, gr_y_y, gr_y_z, grr_z_y_x, grr_z_y_y, grr_z_y_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_y_x[i] = ts_y_x[i] * gfe_0 * gc_z[i] + gr_y_x[i] * gc_z[i];

        grr_z_y_y[i] = ts_y_y[i] * gfe_0 * gc_z[i] + gr_y_y[i] * gc_z[i];

        grr_z_y_z[i] = ts_y_0[i] * gfe2_0 + gr_y_0[i] * gfe_0 + ts_y_z[i] * gfe_0 * gc_z[i] + gr_y_z[i] * gc_z[i];
    }

    // Set up 24-27 components of targeted buffer : PP

    auto grr_z_z_x = pbuffer.data(idx_gr_pp + 24);

    auto grr_z_z_y = pbuffer.data(idx_gr_pp + 25);

    auto grr_z_z_z = pbuffer.data(idx_gr_pp + 26);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_x, gr_0_y, gr_0_z, gr_z_0, gr_z_x, gr_z_y, gr_z_z, grr_z_z_x, grr_z_z_y, grr_z_z_z, ts_0_x, ts_0_y, ts_0_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_z_x[i] = ts_0_x[i] * gfe2_0 + gr_0_x[i] * gfe_0 + ts_z_x[i] * gfe_0 * gc_z[i] + gr_z_x[i] * gc_z[i];

        grr_z_z_y[i] = ts_0_y[i] * gfe2_0 + gr_0_y[i] * gfe_0 + ts_z_y[i] * gfe_0 * gc_z[i] + gr_z_y[i] * gc_z[i];

        grr_z_z_z[i] = ts_0_z[i] * gfe2_0 + gr_0_z[i] * gfe_0 + ts_z_0[i] * gfe2_0 + gr_z_0[i] * gfe_0 + ts_z_z[i] * gfe_0 * gc_z[i] + gr_z_z[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

