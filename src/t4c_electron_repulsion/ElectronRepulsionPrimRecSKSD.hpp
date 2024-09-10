#ifndef ElectronRepulsionPrimRecSKSD_hpp
#define ElectronRepulsionPrimRecSKSD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SK|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sksd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_shsd The primitive integrals buffer.
/// @param idx_eri_1_shsd The primitive integrals buffer.
/// @param idx_eri_1_sisp The primitive integrals buffer.
/// @param idx_eri_0_sisd The primitive integrals buffer.
/// @param idx_eri_1_sisd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_sksd(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sksd,
                                  size_t idx_eri_0_shsd,
                                  size_t idx_eri_1_shsd,
                                  size_t idx_eri_1_sisp,
                                  size_t idx_eri_0_sisd,
                                  size_t idx_eri_1_sisd,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSKSD_hpp */
