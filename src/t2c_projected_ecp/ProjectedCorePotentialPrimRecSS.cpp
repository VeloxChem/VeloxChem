#include "ProjectedCorePotentialPrimRecSS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_ss(const int l,
                                      const int m,
                                      const int p,
                                      const int q,
                                      CSimdArray<double>& pbuffer,
                                      const size_t idx_ss,
                                      const CSimdArray<double>& i_values,
                                      const CSimdArray<double>& l_values,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_gamma,
                                      const size_t idx_mb,
                                      const TPoint<double>& r_a,
                                      const double a_norm,
                                      const double c_norm) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up normalization factors

    auto b_norms = factors.data(1);
    
    // Set up gamma factors
    
    auto f_gamma = factors.data(idx_gamma);
    
    // Set up I_n values
    
    auto i_vals = i_values.data((size_t)(l + m + p + q));
    
    // Set up L_n values
    
    auto l_vals = l_values.data((size_t)(l));
    
    // Set up components of auxiliary buffer : SS

    auto t_0_0 = pbuffer.data(idx_ss);

    #pragma omp simd aligned(f_gamma, i_vals, l_vals, b_norms, t_0_0 : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_0_0[i] = a_norm * b_norms[i] * c_norm * f_gamma[i] * i_vals[i] * l_vals[i];
    }
}

} // t2pecp namespace

