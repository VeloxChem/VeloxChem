#include "KineticEnergyPrimRecPS.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_ps(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_ps,
                            const size_t              idx_kin_ss,
                            const size_t              idx_ovl_ps,
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

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ovl_ps);

    auto ts_y_0 = pbuffer.data(idx_ovl_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ovl_ps + 2);

    // Set up components of targeted buffer : PS

    auto tk_x_0 = pbuffer.data(idx_kin_ps);

    auto tk_y_0 = pbuffer.data(idx_kin_ps + 1);

    auto tk_z_0 = pbuffer.data(idx_kin_ps + 2);

#pragma omp simd aligned(pa_x, pa_y, pa_z, tk_0_0, tk_x_0, tk_y_0, tk_z_0, ts_x_0, ts_y_0, ts_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_x_0[i] = tk_0_0[i] * pa_x[i] + 2.0 * ts_x_0[i] * fz_0;

        tk_y_0[i] = tk_0_0[i] * pa_y[i] + 2.0 * ts_y_0[i] * fz_0;

        tk_z_0[i] = tk_0_0[i] * pa_z[i] + 2.0 * ts_z_0[i] * fz_0;
    }
}

}  // namespace kinrec
