#ifndef NuclearPotentialGridPrimRecFS
#define NuclearPotentialGridPrimRecFS

#include "SubMatrix.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [F|A|S]  integrals for set of data buffers on given grid.
/// @param buffer The primitive integrals buffer.
/// @param idx_npot_0_fs The index of integral in primitive integrals buffer.
/// @param idx_npot_0_ps The index of integral in primitive integrals buffer.
/// @param idx_npot_1_ps The index of integral in primitive integrals buffer.
/// @param idx_npot_0_ds The index of integral in primitive integrals buffer.
/// @param idx_npot_1_ds The index of integral in primitive integrals buffer.
/// @param rpa_x The Cartesian X distance of R(PA) = P - A.
/// @param rpa_y The Cartesian Y distance of R(PA) = P - A.
/// @param rpa_z The Cartesian Z distance of R(PA) = P - A.
/// @param factor The combined exponential factor.
auto
comp_on_grid_prim_nuclear_potential_fs(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_fs,
                                       const size_t idx_npot_0_ps,
                                       const size_t idx_npot_1_ps,
                                       const size_t idx_npot_0_ds,
                                       const size_t idx_npot_1_ds,
                                       const double rpa_x,
                                       const double rpa_y,
                                       const double rpa_z,
                                       const double factor) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGridPrimRecFS */
