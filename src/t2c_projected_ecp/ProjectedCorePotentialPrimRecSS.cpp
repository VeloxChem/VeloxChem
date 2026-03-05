#include "ProjectedCorePotentialPrimRecSS.hpp"

#include <iostream>

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
    
    // Set up normalization factors

    auto t_args = factors.data(8);
    
    // Set up gamma factors
    
    auto f_gamma = factors.data(idx_gamma);
    
    // Set up I_n values
    
    auto i_vals = i_values.data((size_t)(l + m + p + q));
    
    // Set up L_n values
    
    auto l_vals = l_values.data((size_t)(l));
    
    // Set up components of auxiliary buffer : SS

    auto t_0_0 = pbuffer.data(idx_ss);

    //#pragma omp simd aligned(f_gamma, i_vals, l_vals, b_norms, t_0_0 : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        if (t_args[i] > 149.99)
        {
            t_0_0[i] = 0.0; 
        }
        else
        {
            t_0_0[i] = a_norm * b_norms[i] * c_norm * f_gamma[i] * i_vals[i] * l_vals[i];
        }
    }
    
    if (m > 0)
    {
        // Set up normalization factors

        auto mb = factors.data(idx_mb);
        
        for (int k = 0; k < 2 * m; k++)
        {
            #pragma omp simd aligned(mb, t_0_0 : 64)
            for (size_t i = 0; i < nelems; i++)
            {
                t_0_0[i] *= mb[i];
            }
        }
    }
    
    if (p > 0)
    {
        // set up Cartesian A coordinates

        const auto xyz = r_a.coordinates();

        const auto a_x = xyz[0];

        const auto a_y = xyz[1];

        const auto a_z = xyz[2];
        
        // compute |A|
        
        const auto ma = std::sqrt(a_x * a_x + a_y * a_y + a_z * a_z);
        
        for (int k = 0; k < 2 * p; k++)
        {
            #pragma omp simd aligned(t_0_0 : 64)
            for (size_t i = 0; i < nelems; i++)
            {
                t_0_0[i] *= ma;
            }
        }
    }
}

} // t2pecp namespace

