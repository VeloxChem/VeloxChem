#include "KineticEnergyPrimRecSP.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_sp(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_sp,
                            const size_t              idx_kin_ss,
                            const size_t              idx_ovl_sp,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpb,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SS

    auto tk_0_0 = pbuffer.data(idx_kin_ss);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

    // Set up components of targeted buffer : SP

    auto tk_0_x = pbuffer.data(idx_kin_sp);

    auto tk_0_y = pbuffer.data(idx_kin_sp + 1);

    auto tk_0_z = pbuffer.data(idx_kin_sp + 2);

#pragma omp simd aligned(pb_x, pb_y, pb_z, tk_0_0, tk_0_x, tk_0_y, tk_0_z, ts_0_x, ts_0_y, ts_0_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_0_x[i] = tk_0_0[i] * pb_x[i] + 2.0 * ts_0_x[i] * fz_0;

        tk_0_y[i] = tk_0_0[i] * pb_y[i] + 2.0 * ts_0_y[i] * fz_0;

        tk_0_z[i] = tk_0_0[i] * pb_z[i] + 2.0 * ts_0_z[i] * fz_0;
    }
}

}  // namespace kinrec
