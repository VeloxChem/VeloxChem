#include "OverlapPrimRecSS.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_ss(CSimdArray<double>& pbuffer,
                     const size_t        idx_ovl_ss,
                     CSimdArray<double>& factors,
                     const double        a_exp,
                     const double        a_norm) -> void
{
    const double fpi = mathconst::pi_value();

    // Set up exponents, normalization factors

    auto b_exps = factors.data(0);

    auto b_norms = factors.data(1);

    // Set up R(AB) distances

    auto ab_x = factors.data(5);

    auto ab_y = factors.data(6);

    auto ab_z = factors.data(7);

    // Set up components of targeted buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    /// compute primitive integrals

    const auto nelems = pbuffer.number_of_active_elements();

#pragma omp simd aligned(ts_0_0, ab_x, ab_y, ab_z, b_exps, b_norms : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fe_0 = 1.0 / (a_exp + b_exps[i]);

        double fz_0 = a_exp * b_exps[i] * fe_0 * (ab_x[i] * ab_x[i] + ab_y[i] * ab_y[i] + ab_z[i] * ab_z[i]);

        fe_0 *= fpi;

        ts_0_0[i] = b_norms[i] * a_norm * fe_0 * std::sqrt(fe_0) * std::exp(-fz_0);
    }
}

}  // namespace ovlrec
