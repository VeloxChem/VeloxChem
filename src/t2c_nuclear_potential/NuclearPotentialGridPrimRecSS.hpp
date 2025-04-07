#ifndef NuclearPotentialGridPrimRecSS_hpp
#define NuclearPotentialGridPrimRecSS_hpp

#include "SubMatrix.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes primitive [S|A|S]  integrals on grid for set of data buffers.
/// @param buffer The primitive integrals buffer.
/// @param idx_npot_ss The index of integral in primitive integrals buffer.
/// @param idx_bf_vals The primary row index of values in Boys function data.
/// @param foverlap The overlap factor.
/// @param factor The scaling factor.
auto comp_on_grid_prim_nuclear_potential_ss(CSubMatrix&  buffer,
                                            const size_t idx_npot_ss,
                                            const size_t idx_bf_vals,
                                            const double foverlap,
                                            const double factor) -> void;
}  // namespace npotrec

#endif /* NuclearPotentialGridPrimRecSS_hpp */
