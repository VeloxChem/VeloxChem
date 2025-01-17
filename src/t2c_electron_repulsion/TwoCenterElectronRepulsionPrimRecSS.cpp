#include "TwoCenterElectronRepulsionPrimRecSS.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ss(CSimdArray<double>& pbuffer,
                                const size_t idx_eri_ss,
                                const CSimdArray<double>& bf_data,
                                const size_t              idx_vals,
                                CSimdArray<double>&       factors,
                                const double              a_exp,
                                const double              a_norm) -> void
{
    const auto fpi = mathconst::pi_value();
    
    const auto fact = 2.0 * fpi * fpi * std::sqrt(fpi);

    // Set up exponents, normalization factors

    auto b_exps = factors.data(0);

    auto b_norms = factors.data(1);

    /// Boys function values

    auto bvals = bf_data.data(idx_vals);

    /// Set up components of auxiliary buffer : SS

    auto tg_0_0 = pbuffer.data(idx_eri_ss);

    /// compute primitive integrals

    const auto nelems = pbuffer.number_of_active_elements();

#pragma omp simd aligned(tg_0_0, bvals, b_exps, b_norms : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_0[i] = fact * a_norm * b_norms[i] * bvals[i] / (a_exp * b_exps[i] * std::sqrt(a_exp + b_exps[i])) ;
    }
}

} // t2ceri namespace
