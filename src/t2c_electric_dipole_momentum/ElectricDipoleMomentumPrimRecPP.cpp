#include "ElectricDipoleMomentumPrimRecPP.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_pp(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_pp,
                                      const size_t              idx_dip_ss,
                                      const size_t              idx_ovl_sp,
                                      const size_t              idx_dip_sp,
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

    auto tr_x_0_0 = pbuffer.data(idx_dip_ss);

    auto tr_y_0_0 = pbuffer.data(idx_dip_ss + 1);

    auto tr_z_0_0 = pbuffer.data(idx_dip_ss + 2);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

    // Set up components of auxiliary buffer : SP

    auto tr_x_0_x = pbuffer.data(idx_dip_sp);

    auto tr_x_0_y = pbuffer.data(idx_dip_sp + 1);

    auto tr_x_0_z = pbuffer.data(idx_dip_sp + 2);

    auto tr_y_0_x = pbuffer.data(idx_dip_sp + 3);

    auto tr_y_0_y = pbuffer.data(idx_dip_sp + 4);

    auto tr_y_0_z = pbuffer.data(idx_dip_sp + 5);

    auto tr_z_0_x = pbuffer.data(idx_dip_sp + 6);

    auto tr_z_0_y = pbuffer.data(idx_dip_sp + 7);

    auto tr_z_0_z = pbuffer.data(idx_dip_sp + 8);

    // Set up 0-3 components of targeted buffer : PP

    auto tr_x_x_x = pbuffer.data(idx_dip_pp);

    auto tr_x_x_y = pbuffer.data(idx_dip_pp + 1);

    auto tr_x_x_z = pbuffer.data(idx_dip_pp + 2);

#pragma omp simd aligned(pa_x, tr_x_0_0, tr_x_0_x, tr_x_0_y, tr_x_0_z, tr_x_x_x, tr_x_x_y, tr_x_x_z, ts_0_x, ts_0_y, ts_0_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_x_x[i] = tr_x_0_0[i] * fe_0 + ts_0_x[i] * fe_0 + tr_x_0_x[i] * pa_x[i];

        tr_x_x_y[i] = ts_0_y[i] * fe_0 + tr_x_0_y[i] * pa_x[i];

        tr_x_x_z[i] = ts_0_z[i] * fe_0 + tr_x_0_z[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : PP

    auto tr_x_y_x = pbuffer.data(idx_dip_pp + 3);

    auto tr_x_y_y = pbuffer.data(idx_dip_pp + 4);

    auto tr_x_y_z = pbuffer.data(idx_dip_pp + 5);

#pragma omp simd aligned(pa_y, tr_x_0_0, tr_x_0_x, tr_x_0_y, tr_x_0_z, tr_x_y_x, tr_x_y_y, tr_x_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_y_x[i] = tr_x_0_x[i] * pa_y[i];

        tr_x_y_y[i] = tr_x_0_0[i] * fe_0 + tr_x_0_y[i] * pa_y[i];

        tr_x_y_z[i] = tr_x_0_z[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : PP

    auto tr_x_z_x = pbuffer.data(idx_dip_pp + 6);

    auto tr_x_z_y = pbuffer.data(idx_dip_pp + 7);

    auto tr_x_z_z = pbuffer.data(idx_dip_pp + 8);

#pragma omp simd aligned(pa_z, tr_x_0_0, tr_x_0_x, tr_x_0_y, tr_x_0_z, tr_x_z_x, tr_x_z_y, tr_x_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_z_x[i] = tr_x_0_x[i] * pa_z[i];

        tr_x_z_y[i] = tr_x_0_y[i] * pa_z[i];

        tr_x_z_z[i] = tr_x_0_0[i] * fe_0 + tr_x_0_z[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : PP

    auto tr_y_x_x = pbuffer.data(idx_dip_pp + 9);

    auto tr_y_x_y = pbuffer.data(idx_dip_pp + 10);

    auto tr_y_x_z = pbuffer.data(idx_dip_pp + 11);

#pragma omp simd aligned(pa_x, tr_y_0_0, tr_y_0_x, tr_y_0_y, tr_y_0_z, tr_y_x_x, tr_y_x_y, tr_y_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_x_x[i] = tr_y_0_0[i] * fe_0 + tr_y_0_x[i] * pa_x[i];

        tr_y_x_y[i] = tr_y_0_y[i] * pa_x[i];

        tr_y_x_z[i] = tr_y_0_z[i] * pa_x[i];
    }

    // Set up 12-15 components of targeted buffer : PP

    auto tr_y_y_x = pbuffer.data(idx_dip_pp + 12);

    auto tr_y_y_y = pbuffer.data(idx_dip_pp + 13);

    auto tr_y_y_z = pbuffer.data(idx_dip_pp + 14);

#pragma omp simd aligned(pa_y, tr_y_0_0, tr_y_0_x, tr_y_0_y, tr_y_0_z, tr_y_y_x, tr_y_y_y, tr_y_y_z, ts_0_x, ts_0_y, ts_0_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_y_x[i] = ts_0_x[i] * fe_0 + tr_y_0_x[i] * pa_y[i];

        tr_y_y_y[i] = tr_y_0_0[i] * fe_0 + ts_0_y[i] * fe_0 + tr_y_0_y[i] * pa_y[i];

        tr_y_y_z[i] = ts_0_z[i] * fe_0 + tr_y_0_z[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : PP

    auto tr_y_z_x = pbuffer.data(idx_dip_pp + 15);

    auto tr_y_z_y = pbuffer.data(idx_dip_pp + 16);

    auto tr_y_z_z = pbuffer.data(idx_dip_pp + 17);

#pragma omp simd aligned(pa_z, tr_y_0_0, tr_y_0_x, tr_y_0_y, tr_y_0_z, tr_y_z_x, tr_y_z_y, tr_y_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_z_x[i] = tr_y_0_x[i] * pa_z[i];

        tr_y_z_y[i] = tr_y_0_y[i] * pa_z[i];

        tr_y_z_z[i] = tr_y_0_0[i] * fe_0 + tr_y_0_z[i] * pa_z[i];
    }

    // Set up 18-21 components of targeted buffer : PP

    auto tr_z_x_x = pbuffer.data(idx_dip_pp + 18);

    auto tr_z_x_y = pbuffer.data(idx_dip_pp + 19);

    auto tr_z_x_z = pbuffer.data(idx_dip_pp + 20);

#pragma omp simd aligned(pa_x, tr_z_0_0, tr_z_0_x, tr_z_0_y, tr_z_0_z, tr_z_x_x, tr_z_x_y, tr_z_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_x_x[i] = tr_z_0_0[i] * fe_0 + tr_z_0_x[i] * pa_x[i];

        tr_z_x_y[i] = tr_z_0_y[i] * pa_x[i];

        tr_z_x_z[i] = tr_z_0_z[i] * pa_x[i];
    }

    // Set up 21-24 components of targeted buffer : PP

    auto tr_z_y_x = pbuffer.data(idx_dip_pp + 21);

    auto tr_z_y_y = pbuffer.data(idx_dip_pp + 22);

    auto tr_z_y_z = pbuffer.data(idx_dip_pp + 23);

#pragma omp simd aligned(pa_y, tr_z_0_0, tr_z_0_x, tr_z_0_y, tr_z_0_z, tr_z_y_x, tr_z_y_y, tr_z_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_y_x[i] = tr_z_0_x[i] * pa_y[i];

        tr_z_y_y[i] = tr_z_0_0[i] * fe_0 + tr_z_0_y[i] * pa_y[i];

        tr_z_y_z[i] = tr_z_0_z[i] * pa_y[i];
    }

    // Set up 24-27 components of targeted buffer : PP

    auto tr_z_z_x = pbuffer.data(idx_dip_pp + 24);

    auto tr_z_z_y = pbuffer.data(idx_dip_pp + 25);

    auto tr_z_z_z = pbuffer.data(idx_dip_pp + 26);

#pragma omp simd aligned(pa_z, tr_z_0_0, tr_z_0_x, tr_z_0_y, tr_z_0_z, tr_z_z_x, tr_z_z_y, tr_z_z_z, ts_0_x, ts_0_y, ts_0_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_z_x[i] = ts_0_x[i] * fe_0 + tr_z_0_x[i] * pa_z[i];

        tr_z_z_y[i] = ts_0_y[i] * fe_0 + tr_z_0_y[i] * pa_z[i];

        tr_z_z_z[i] = tr_z_0_0[i] * fe_0 + ts_0_z[i] * fe_0 + tr_z_0_z[i] * pa_z[i];
    }
}

}  // namespace diprec
