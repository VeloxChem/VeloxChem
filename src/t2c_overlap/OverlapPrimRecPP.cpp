#include "OverlapPrimRecPP.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_pp(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_pp,
                     const size_t              idx_ovl_ss,
                     const size_t              idx_ovl_sp,
                     const CSimdArray<double>& factors,
                     const size_t              idx_rpa,
                     const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

    // Set up 0-3 components of targeted buffer : PP

    auto ts_x_x = pbuffer.data(idx_ovl_pp);

    auto ts_x_y = pbuffer.data(idx_ovl_pp + 1);

    auto ts_x_z = pbuffer.data(idx_ovl_pp + 2);

#pragma omp simd aligned(pa_x, ts_0_0, ts_0_x, ts_0_y, ts_0_z, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_x_x[i] = ts_0_0[i] * fe_0 + ts_0_x[i] * pa_x[i];

        ts_x_y[i] = ts_0_y[i] * pa_x[i];

        ts_x_z[i] = ts_0_z[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : PP

    auto ts_y_x = pbuffer.data(idx_ovl_pp + 3);

    auto ts_y_y = pbuffer.data(idx_ovl_pp + 4);

    auto ts_y_z = pbuffer.data(idx_ovl_pp + 5);

#pragma omp simd aligned(pa_y, ts_0_0, ts_0_x, ts_0_y, ts_0_z, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_y_x[i] = ts_0_x[i] * pa_y[i];

        ts_y_y[i] = ts_0_0[i] * fe_0 + ts_0_y[i] * pa_y[i];

        ts_y_z[i] = ts_0_z[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : PP

    auto ts_z_x = pbuffer.data(idx_ovl_pp + 6);

    auto ts_z_y = pbuffer.data(idx_ovl_pp + 7);

    auto ts_z_z = pbuffer.data(idx_ovl_pp + 8);

#pragma omp simd aligned(pa_z, ts_0_0, ts_0_x, ts_0_y, ts_0_z, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_z_x[i] = ts_0_x[i] * pa_z[i];

        ts_z_y[i] = ts_0_y[i] * pa_z[i];

        ts_z_z[i] = ts_0_0[i] * fe_0 + ts_0_z[i] * pa_z[i];
    }
}

}  // namespace ovlrec