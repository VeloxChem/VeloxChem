#include "LocalCorePotentialPrimRecSS.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_ss(CSimdArray<double>& pbuffer,
                                  const size_t idx_ss,
                                  const CSimdArray<double>& factors,
                                  const TPoint<double>& r_a,
                                  const double a_exp,
                                  const double c_exp,
                                  const double a_norm,
                                  const double c_norm) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up pi constant
    
    const double fpi = mathconst::pi_value();
    
    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up normalization factors

    auto b_norms = factors.data(1);
    
    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // set up Cartesian B coordinates

    auto b_x = factors.data(2);

    auto b_y = factors.data(3);

    auto b_z = factors.data(4);
    
    // Set up R coordinates

    auto r_x = factors.data(5);

    auto r_y = factors.data(6);

    auto r_z = factors.data(7);

    // Set up components of auxiliary buffer : SS

    auto t_0_0 = pbuffer.data(idx_ss);

    #pragma omp simd aligned(r_x, r_y, r_z, b_x, b_y, b_z, b_exps, b_norms, t_0_0 : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fzeta = b_exps[i] + a_exp + c_exp;
        
        double fact = fpi / fzeta;
        
        double f2rab = fzeta * (r_x[i] * r_x[i] + r_y[i] * r_y[i] + r_z[i] * r_z[i])
              
                     - a_exp * (a_x * a_x + a_y * a_y + a_z * a_z)
        
                     - b_exps[i] * (b_x[i] * b_x[i] + b_y[i] * b_y[i] + b_z[i] * b_z[i]);
        
        
        t_0_0[i] = a_norm * b_norms[i] * c_norm * fact * std::sqrt(fact) * std::exp(f2rab);
    }
}

} // t2lecp namespace

