#ifndef NuclearPotentialGridPrimRecFD
#define NuclearPotentialGridPrimRecFD

#include "SubMatrix.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [F|A|D]  integrals for set of data buffers on given grid.
/// @param buffer The primitive integrals buffer.
/// @param idx_npot_0_fd The index of integral in primitive integrals buffer.
/// @param idx_npot_0_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_pd The index of integral in primitive integrals buffer.
/// @param idx_npot_0_dp The index of integral in primitive integrals buffer.
/// @param idx_npot_1_dp The index of integral in primitive integrals buffer.
/// @param idx_npot_0_dd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_dd The index of integral in primitive integrals buffer.
/// @param rpa_x The Cartesian X distance of R(PA) = P - A.
/// @param rpa_y The Cartesian Y distance of R(PA) = P - A.
/// @param rpa_z The Cartesian Z distance of R(PA) = P - A.
/// @param factor The combined exponential factor.
auto
comp_on_grid_prim_nuclear_potential_fd(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_fd,
                                       const size_t idx_npot_0_pd,
                                       const size_t idx_npot_1_pd,
                                       const size_t idx_npot_0_dp,
                                       const size_t idx_npot_1_dp,
                                       const size_t idx_npot_0_dd,
                                       const size_t idx_npot_1_dd,
                                       const double rpa_x,
                                       const double rpa_y,
                                       const double rpa_z,
                                       const double factor) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGridPrimRecFD */
