#ifndef NuclearPotentialGridPrimRecDF
#define NuclearPotentialGridPrimRecDF

#include "SubMatrix.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [D|A|F]  integrals for set of data buffers on given grid.
/// @param buffer The primitive integrals buffer.
/// @param idx_npot_0_df The index of integral in primitive integrals buffer.
/// @param idx_npot_0_sf The index of integral in primitive integrals buffer.
/// @param idx_npot_1_sf The index of integral in primitive integrals buffer.
/// @param idx_npot_0_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_0_pf The index of integral in primitive integrals buffer.
/// @param idx_npot_1_pf The index of integral in primitive integrals buffer.
/// @param rpa_x The Cartesian X distance of R(PA) = P - A.
/// @param rpa_y The Cartesian Y distance of R(PA) = P - A.
/// @param rpa_z The Cartesian Z distance of R(PA) = P - A.
/// @param factor The combined exponential factor.
auto
comp_on_grid_prim_nuclear_potential_df(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_df,
                                       const size_t idx_npot_0_sf,
                                       const size_t idx_npot_1_sf,
                                       const size_t idx_npot_0_pd,
                                       const size_t idx_npot_1_pd,
                                       const size_t idx_npot_0_pf,
                                       const size_t idx_npot_1_pf,
                                       const double rpa_x,
                                       const double rpa_y,
                                       const double rpa_z,
                                       const double factor) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGridPrimRecDF */
