#include "NuclearPotentialPrimRecSS.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_ss(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_ss,
                               const size_t              idx_ovl_ss,
                               const CSimdArray<double>& bf_data,
                               const size_t              idx_vals,
                               CSimdArray<double>&       factors,
                               const double              a_exp) -> void
{
    const auto fpi = 2.0 / std::sqrt(mathconst::pi_value());

    // Set up exponents

    auto b_exps = factors.data(0);

    /// Boys function values

    auto bvals = bf_data.data(idx_vals);

    /// Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    /// Set up components of targeted buffer : SS

    auto ta_0_0 = pbuffer.data(idx_npot_ss);

    /// compute primitive integrals

    const auto nelems = pbuffer.number_of_active_elements();

#pragma omp simd aligned(ta_0_0, ts_0_0, bvals : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_0_0[i] = fpi * std::sqrt(a_exp + b_exps[i]) * ts_0_0[i] * bvals[i];
    }
}

}  // namespace npotrec
