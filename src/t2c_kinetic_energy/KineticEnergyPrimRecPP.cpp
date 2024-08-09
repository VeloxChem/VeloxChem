#include "KineticEnergyPrimRecPP.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_pp(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_pp,
                            const size_t              idx_kin_ss,
                            const size_t              idx_kin_sp,
                            const size_t              idx_ovl_pp,
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

    auto tk_0_0 = pbuffer.data(idx_kin_ss);

    // Set up components of auxiliary buffer : SP

    auto tk_0_x = pbuffer.data(idx_kin_sp);

    auto tk_0_y = pbuffer.data(idx_kin_sp + 1);

    auto tk_0_z = pbuffer.data(idx_kin_sp + 2);

    // Set up components of auxiliary buffer : PP

    auto ts_x_x = pbuffer.data(idx_ovl_pp);

    auto ts_x_y = pbuffer.data(idx_ovl_pp + 1);

    auto ts_x_z = pbuffer.data(idx_ovl_pp + 2);

    auto ts_y_x = pbuffer.data(idx_ovl_pp + 3);

    auto ts_y_y = pbuffer.data(idx_ovl_pp + 4);

    auto ts_y_z = pbuffer.data(idx_ovl_pp + 5);

    auto ts_z_x = pbuffer.data(idx_ovl_pp + 6);

    auto ts_z_y = pbuffer.data(idx_ovl_pp + 7);

    auto ts_z_z = pbuffer.data(idx_ovl_pp + 8);

    // Set up 0-3 components of targeted buffer : PP

    auto tk_x_x = pbuffer.data(idx_kin_pp);

    auto tk_x_y = pbuffer.data(idx_kin_pp + 1);

    auto tk_x_z = pbuffer.data(idx_kin_pp + 2);

#pragma omp simd aligned( \
        pa_x, tk_0_0, tk_0_x, tk_0_y, tk_0_z, tk_x_x, tk_x_y, tk_x_z, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_x_x[i] = tk_0_0[i] * fe_0 + tk_0_x[i] * pa_x[i] + 2.0 * ts_x_x[i] * fz_0;

        tk_x_y[i] = tk_0_y[i] * pa_x[i] + 2.0 * ts_x_y[i] * fz_0;

        tk_x_z[i] = tk_0_z[i] * pa_x[i] + 2.0 * ts_x_z[i] * fz_0;
    }

    // Set up 3-6 components of targeted buffer : PP

    auto tk_y_x = pbuffer.data(idx_kin_pp + 3);

    auto tk_y_y = pbuffer.data(idx_kin_pp + 4);

    auto tk_y_z = pbuffer.data(idx_kin_pp + 5);

#pragma omp simd aligned( \
        pa_y, tk_0_0, tk_0_x, tk_0_y, tk_0_z, tk_y_x, tk_y_y, tk_y_z, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_y_x[i] = tk_0_x[i] * pa_y[i] + 2.0 * ts_y_x[i] * fz_0;

        tk_y_y[i] = tk_0_0[i] * fe_0 + tk_0_y[i] * pa_y[i] + 2.0 * ts_y_y[i] * fz_0;

        tk_y_z[i] = tk_0_z[i] * pa_y[i] + 2.0 * ts_y_z[i] * fz_0;
    }

    // Set up 6-9 components of targeted buffer : PP

    auto tk_z_x = pbuffer.data(idx_kin_pp + 6);

    auto tk_z_y = pbuffer.data(idx_kin_pp + 7);

    auto tk_z_z = pbuffer.data(idx_kin_pp + 8);

#pragma omp simd aligned( \
        pa_z, tk_0_0, tk_0_x, tk_0_y, tk_0_z, tk_z_x, tk_z_y, tk_z_z, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_z_x[i] = tk_0_x[i] * pa_z[i] + 2.0 * ts_z_x[i] * fz_0;

        tk_z_y[i] = tk_0_y[i] * pa_z[i] + 2.0 * ts_z_y[i] * fz_0;

        tk_z_z[i] = tk_0_0[i] * fe_0 + tk_0_z[i] * pa_z[i] + 2.0 * ts_z_z[i] * fz_0;
    }
}

}  // namespace kinrec
