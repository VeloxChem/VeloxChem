#ifndef ElectronRepulsionPrimRecSLSI_hpp
#define ElectronRepulsionPrimRecSLSI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SL|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_slsi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sisi The primitive integrals buffer.
/// @param idx_eri_1_sisi The primitive integrals buffer.
/// @param idx_eri_1_sksh The primitive integrals buffer.
/// @param idx_eri_0_sksi The primitive integrals buffer.
/// @param idx_eri_1_sksi The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_slsi(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_slsi,
                                  size_t idx_eri_0_sisi,
                                  size_t idx_eri_1_sisi,
                                  size_t idx_eri_1_sksh,
                                  size_t idx_eri_0_sksi,
                                  size_t idx_eri_1_sksi,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSLSI_hpp */
