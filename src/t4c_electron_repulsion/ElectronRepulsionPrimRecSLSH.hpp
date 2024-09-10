#ifndef ElectronRepulsionPrimRecSLSH_hpp
#define ElectronRepulsionPrimRecSLSH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SL|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_slsh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sish The primitive integrals buffer.
/// @param idx_eri_1_sish The primitive integrals buffer.
/// @param idx_eri_1_sksg The primitive integrals buffer.
/// @param idx_eri_0_sksh The primitive integrals buffer.
/// @param idx_eri_1_sksh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_slsh(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_slsh,
                                  size_t idx_eri_0_sish,
                                  size_t idx_eri_1_sish,
                                  size_t idx_eri_1_sksg,
                                  size_t idx_eri_0_sksh,
                                  size_t idx_eri_1_sksh,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSLSH_hpp */
