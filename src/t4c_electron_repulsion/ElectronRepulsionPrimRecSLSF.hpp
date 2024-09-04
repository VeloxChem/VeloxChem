#ifndef ElectronRepulsionPrimRecSLSF_hpp
#define ElectronRepulsionPrimRecSLSF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SL|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_slsf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sisf The primitive integrals buffer.
/// @param idx_eri_1_sisf The primitive integrals buffer.
/// @param idx_eri_1_sksd The primitive integrals buffer.
/// @param idx_eri_0_sksf The primitive integrals buffer.
/// @param idx_eri_1_sksf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_slsf(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_slsf,
                                  size_t idx_eri_0_sisf,
                                  size_t idx_eri_1_sisf,
                                  size_t idx_eri_1_sksd,
                                  size_t idx_eri_0_sksf,
                                  size_t idx_eri_1_sksf,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSLSF_hpp */
