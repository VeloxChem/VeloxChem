#include "NuclearPotentialGridPrimRecSS.hpp"

#include "MathConst.hpp"

namespace npotrec {  // npotrec namespace

auto comp_on_grid_prim_nuclear_potential_ss(CSubMatrix&  buffer,
                                            const size_t idx_npot_ss,
                                            const size_t idx_bf_vals,
                                            const double foverlap,
                                            const double factor) -> void
{
    const auto fpi = 2.0 / std::sqrt(mathconst::pi_value());

    // set up number of grid points
    
    const auto nelems = buffer.number_of_columns();
    
    /// Boys function values

    auto bvals = &(buffer.data()[idx_bf_vals * nelems]);

    /// Set up components of auxiliary buffer : SS

    auto ta_0_0 = &(buffer.data()[idx_npot_ss * nelems]);

    /// compute primitive integrals

    const auto fact = fpi * std::sqrt(factor) * foverlap;
    
    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        ta_0_0[i] = fact * bvals[i];
    }
}

}  // namespace npotrec
