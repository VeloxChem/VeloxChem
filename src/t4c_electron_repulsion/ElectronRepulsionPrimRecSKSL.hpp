#ifndef ElectronRepulsionPrimRecSKSL_hpp
#define ElectronRepulsionPrimRecSKSL_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SK|1/|r-r'||SL]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sksl The index of integral in primitive integrals buffer.
/// @param idx_eri_0_shsl The primitive integrals buffer.
/// @param idx_eri_1_shsl The primitive integrals buffer.
/// @param idx_eri_1_sisk The primitive integrals buffer.
/// @param idx_eri_0_sisl The primitive integrals buffer.
/// @param idx_eri_1_sisl The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_sksl(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sksl,
                                  size_t idx_eri_0_shsl,
                                  size_t idx_eri_1_shsl,
                                  size_t idx_eri_1_sisk,
                                  size_t idx_eri_0_sisl,
                                  size_t idx_eri_1_sisl,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSKSL_hpp */
