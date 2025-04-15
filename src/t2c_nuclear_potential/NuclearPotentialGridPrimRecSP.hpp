#ifndef NuclearPotentialGridPrimRecSP
#define NuclearPotentialGridPrimRecSP

#include "SubMatrix.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes primitive [S|A|P]  integrals for set of data buffers on given grid.
/// @param buffer The primitive integrals buffer.
/// @param idx_npot_0_sp The index of integral in primitive integrals buffer.
/// @param idx_npot_0_ss The index of integral in primitive integrals buffer.
/// @param idx_npot_1_ss The index of integral in primitive integrals buffer.
/// @param rpb_x The Cartesian X distance of R(PB) = P - B.
/// @param rpb_y The Cartesian Y distance of R(PB) = P - B.
/// @param rpb_z The Cartesian Z distance of R(PB) = P - B.
auto
comp_on_grid_prim_nuclear_potential_sp(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_sp,
                                       const size_t idx_npot_0_ss,
                                       const size_t idx_npot_1_ss,
                                       const double rpb_x,
                                       const double rpb_y,
                                       const double rpb_z) -> void;
} // npotrec namespace

#endif /* NuclearPotentialGridPrimRecSP */
