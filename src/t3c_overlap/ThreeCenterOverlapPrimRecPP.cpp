#include "ThreeCenterOverlapPrimRecPP.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_pp(CSimdArray<double>& pbuffer, 
                     const size_t idx_pp,
                     const size_t idx_ss,
                     const size_t idx_sp,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ss);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_sp);

    auto ts_0_y = pbuffer.data(idx_sp + 1);

    auto ts_0_z = pbuffer.data(idx_sp + 2);

    // Set up 0-3 components of targeted buffer : PP

    auto ts_x_x = pbuffer.data(idx_pp);

    auto ts_x_y = pbuffer.data(idx_pp + 1);

    auto ts_x_z = pbuffer.data(idx_pp + 2);

    #pragma omp simd aligned(ga_x, ts_0_0, ts_0_x, ts_0_y, ts_0_z, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_x_x[i] = ts_0_0[i] * gfe_0 + ts_0_x[i] * ga_x[i];

        ts_x_y[i] = ts_0_y[i] * ga_x[i];

        ts_x_z[i] = ts_0_z[i] * ga_x[i];
    }

    // Set up 3-6 components of targeted buffer : PP

    auto ts_y_x = pbuffer.data(idx_pp + 3);

    auto ts_y_y = pbuffer.data(idx_pp + 4);

    auto ts_y_z = pbuffer.data(idx_pp + 5);

    #pragma omp simd aligned(ga_y, ts_0_0, ts_0_x, ts_0_y, ts_0_z, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_y_x[i] = ts_0_x[i] * ga_y[i];

        ts_y_y[i] = ts_0_0[i] * gfe_0 + ts_0_y[i] * ga_y[i];

        ts_y_z[i] = ts_0_z[i] * ga_y[i];
    }

    // Set up 6-9 components of targeted buffer : PP

    auto ts_z_x = pbuffer.data(idx_pp + 6);

    auto ts_z_y = pbuffer.data(idx_pp + 7);

    auto ts_z_z = pbuffer.data(idx_pp + 8);

    #pragma omp simd aligned(ga_z, ts_0_0, ts_0_x, ts_0_y, ts_0_z, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_z_x[i] = ts_0_x[i] * ga_z[i];

        ts_z_y[i] = ts_0_y[i] * ga_z[i];

        ts_z_z[i] = ts_0_0[i] * gfe_0 + ts_0_z[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

