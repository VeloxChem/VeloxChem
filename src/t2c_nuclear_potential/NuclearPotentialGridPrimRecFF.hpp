#ifndef NuclearPotentialGridPrimRecFF
#define NuclearPotentialGridPrimRecFF

#include "SubMatrix.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [F|A|F]  integrals for set of data buffers on given grid.
/// @param buffer The primitive integrals buffer.
/// @param idx_npot_0_ff The index of integral in primitive integrals buffer.
/// @param idx_npot_0_pf The index of integral in primitive integrals buffer.
/// @param idx_npot_1_pf The index of integral in primitive integrals buffer.
/// @param idx_npot_0_dd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_dd The index of integral in primitive integrals buffer.
/// @param idx_npot_0_df The index of integral in primitive integrals buffer.
/// @param idx_npot_1_df The index of integral in primitive integrals buffer.
/// @param rpa_x The Cartesian X distance of R(PA) = P - A.
/// @param rpa_y The Cartesian Y distance of R(PA) = P - A.
/// @param rpa_z The Cartesian Z distance of R(PA) = P - A.
/// @param factor The combined exponential factor.
auto
comp_on_grid_prim_nuclear_potential_ff(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_ff,
                                       const size_t idx_npot_0_pf,
                                       const size_t idx_npot_1_pf,
                                       const size_t idx_npot_0_dd,
                                       const size_t idx_npot_1_dd,
                                       const size_t idx_npot_0_df,
                                       const size_t idx_npot_1_df,
                                       const double rpa_x,
                                       const double rpa_y,
                                       const double rpa_z,
                                       const double factor) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGridPrimRecFF */
