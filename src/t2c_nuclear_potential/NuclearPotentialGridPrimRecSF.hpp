#ifndef NuclearPotentialGridPrimRecSF
#define NuclearPotentialGridPrimRecSF

#include "SubMatrix.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [S|A|F]  integrals for set of data buffers on given grid.
/// @param buffer The primitive integrals buffer.
/// @param idx_npot_0_sf The index of integral in primitive integrals buffer.
/// @param idx_npot_0_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_1_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_0_sd The index of integral in primitive integrals buffer.
/// @param idx_npot_1_sd The index of integral in primitive integrals buffer.
/// @param rpb_x The Cartesian X distance of R(PB) = P - B.
/// @param rpb_y The Cartesian Y distance of R(PB) = P - B.
/// @param rpb_z The Cartesian Z distance of R(PB) = P - B.
/// @param factor The combined exponential factor.
auto
comp_on_grid_prim_nuclear_potential_sf(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_sf,
                                       const size_t idx_npot_0_sp,
                                       const size_t idx_npot_1_sp,
                                       const size_t idx_npot_0_sd,
                                       const size_t idx_npot_1_sd,
                                       const double rpb_x,
                                       const double rpb_y,
                                       const double rpb_z,
                                       const double factor) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGridPrimRecSF */
