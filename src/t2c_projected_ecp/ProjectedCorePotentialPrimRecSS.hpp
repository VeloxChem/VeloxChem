#ifndef ProjectedCorePotentialPrimRecSS
#define ProjectedCorePotentialPrimRecSS

#include "SimdArray.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes primitive [S|U_l|S]  integrals for set of data buffers.
/// @param l The momentum of projector.
/// @param m The m order of MD  scheme.
/// @param p The p order of MD  scheme.
/// @param q The q order of MD  scheme.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_ss The index of integral in primitive integrals buffer.
/// @param i_values The In buffer.
/// @param l_values The In buffer.
/// @param factors The primitive factors buffer.
/// @param idx_gamma The index of gamma factors.
/// @param idx_mb The index of |B| coordinates.
/// @param r_a The Cartesian A point coordinates.
/// @param a_norm The primitive basis function normalization factor on center A.
/// @param c_norm The core potential function normalization factor on center C.
auto
comp_prim_projected_core_potential_ss(const int l,
                                      const int m,
                                      const int p,
                                      const int q,
                                      CSimdArray<double>& pbuffer,
                                      const size_t idx_ss,
                                      const CSimdArray<double>& i_values,
                                      const CSimdArray<double>& l_values,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_gamma,
                                      const size_t idx_b,
                                      const TPoint<double>& r_a,
                                      const double a_norm,
                                      const double c_norm) -> void;


} // t2pecp namespace

#endif /* ProjectedCorePotentialPrimRecSS */
