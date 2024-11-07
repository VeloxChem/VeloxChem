#include "KineticEnergyPrimRecSS.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_ss(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_ss,
                            const size_t              idx_ovl_ss,
                            const CSimdArray<double>& factors,
                            const double              a_exp) -> void
{
    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(AB) distances

    auto ab_x = factors.data(5);

    auto ab_y = factors.data(6);

    auto ab_z = factors.data(7);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of targeted buffer : SS

    auto tk_0_0 = pbuffer.data(idx_kin_ss);

    /// compute integrals

    const auto nelems = pbuffer.number_of_active_elements();

#pragma omp simd aligned(ts_0_0, tk_0_0, ab_x, ab_y, ab_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        const double ab_2 = ab_x[i] * ab_x[i] + ab_y[i] * ab_y[i] + ab_z[i] * ab_z[i];

        tk_0_0[i] = fz_0 * (3.0 - 2.0 * fz_0 * ab_2) * ts_0_0[i];
    }
}

}  // namespace kinrec
