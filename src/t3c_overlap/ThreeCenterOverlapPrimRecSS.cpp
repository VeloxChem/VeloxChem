#include "ThreeCenterOverlapPrimRecSS.hpp"

#include <cmath>

namespace t3ovlrec {  // t3ovlrec namespace

auto
comp_prim_overlap_ss(CSimdArray<double>&       pbuffer,
                     const size_t              idx_t3ovl_ss,
                     const size_t              idx_ovl_ss,
                     CSimdArray<double>&       factors,
                     const size_t              idx_rpc,
                     const double              a_exp,
                     const double              g_exp,
                     const double              g_norm) -> void
{
    // Set up exponents

    auto b_exps = factors.data(0);
    
    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    /// Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    /// Set up components of targeted buffer : SS

    auto ta_0_0 = pbuffer.data(idx_t3ovl_ss);

    /// compute primitive integrals

    const auto nelems = pbuffer.number_of_active_elements();

#pragma omp simd aligned(ta_0_0, ts_0_0, pc_x, pc_y, pc_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double frho = a_exp + b_exps[i];
        
        const double fnu = frho / (frho + g_exp);
        
        ta_0_0[i] = g_norm * std::sqrt(fnu) * fnu * ts_0_0[i]
                  * std::exp(-fnu * g_exp * (pc_x[i] * pc_x[i] + pc_y[i] * pc_y[i] + pc_z[i] * pc_z[i]));
    }
}

}  // namespace t3ovlrec
