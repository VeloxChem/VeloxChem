#ifndef ElectronRepulsionPrimRecSISG_hpp
#define ElectronRepulsionPrimRecSISG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SI|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sisg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sgsg The primitive integrals buffer.
/// @param idx_eri_1_sgsg The primitive integrals buffer.
/// @param idx_eri_1_shsf The primitive integrals buffer.
/// @param idx_eri_0_shsg The primitive integrals buffer.
/// @param idx_eri_1_shsg The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_sisg(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sisg,
                                  size_t idx_eri_0_sgsg,
                                  size_t idx_eri_1_sgsg,
                                  size_t idx_eri_1_shsf,
                                  size_t idx_eri_0_shsg,
                                  size_t idx_eri_1_shsg,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSISG_hpp */
