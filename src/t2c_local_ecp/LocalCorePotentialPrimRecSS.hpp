#ifndef LocalCorePotentialPrimRecSS
#define LocalCorePotentialPrimRecSS

#include "SimdArray.hpp"

namespace t2lecp { // t2lecp namespace

/// @brief Computes primitive [S|U_L|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The core potential function exponent on center C.
/// @param a_norm The primitive basis function normalization factor on center A.
/// @param c_norm The core potential function normalization factor on center C.
auto
comp_prim_local_core_potential_ss(CSimdArray<double>& pbuffer,
                                  const size_t idx_ss,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_r,
                                  const size_t idx_zeta,
                                  const TPoint<double>& r_a,
                                  const double a_exp,
                                  const double c_exp,
                                  const double a_norm,
                                  const double c_norm) -> void;
} // t2lecp namespace

#endif /* LocalCorePotentialPrimRecSS */
