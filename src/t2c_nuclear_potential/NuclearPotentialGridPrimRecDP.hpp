#ifndef NuclearPotentialGridPrimRecDP
#define NuclearPotentialGridPrimRecDP

#include "SubMatrix.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [D|A|P]  integrals for set of data buffers on given grid.
/// @param buffer The primitive integrals buffer.
/// @param idx_npot_0_dp The index of integral in primitive integrals buffer.
/// @param idx_npot_0_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_1_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_0_ps The index of integral in primitive integrals buffer.
/// @param idx_npot_1_ps The index of integral in primitive integrals buffer.
/// @param idx_npot_0_pp The index of integral in primitive integrals buffer.
/// @param idx_npot_1_pp The index of integral in primitive integrals buffer.
/// @param rpa_x The Cartesian X distance of R(PA) = P - A.
/// @param rpa_y The Cartesian Y distance of R(PA) = P - A.
/// @param rpa_z The Cartesian Z distance of R(PA) = P - A.
/// @param factor The combined exponential factor.
auto
comp_on_grid_prim_nuclear_potential_dp(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_dp,
                                       const size_t idx_npot_0_sp,
                                       const size_t idx_npot_1_sp,
                                       const size_t idx_npot_0_ps,
                                       const size_t idx_npot_1_ps,
                                       const size_t idx_npot_0_pp,
                                       const size_t idx_npot_1_pp,
                                       const double rpa_x,
                                       const double rpa_y,
                                       const double rpa_z,
                                       const double factor) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGridPrimRecDP */
