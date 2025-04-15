#ifndef NuclearPotentialGridPrimRecPD
#define NuclearPotentialGridPrimRecPD

#include "SubMatrix.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [P|A|D]  integrals for set of data buffers on given grid.
/// @param buffer The primitive integrals buffer.
/// @param idx_npot_0_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_0_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_1_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_0_sd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_sd The index of integral in primitive integrals buffer.
/// @param rpa_x The Cartesian X distance of R(PA) = P - A.
/// @param rpa_y The Cartesian Y distance of R(PA) = P - A.
/// @param rpa_z The Cartesian Z distance of R(PA) = P - A.
/// @param factor The combined exponential factor.
auto
comp_on_grid_prim_nuclear_potential_pd(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_pd,
                                       const size_t idx_npot_0_sp,
                                       const size_t idx_npot_1_sp,
                                       const size_t idx_npot_0_sd,
                                       const size_t idx_npot_1_sd,
                                       const double rpa_x,
                                       const double rpa_y,
                                       const double rpa_z,
                                       const double factor) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGridPrimRecPD */
