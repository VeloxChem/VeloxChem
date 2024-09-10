#ifndef ElectronRepulsionPrimRecSKSG_hpp
#define ElectronRepulsionPrimRecSKSG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SK|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sksg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_shsg The primitive integrals buffer.
/// @param idx_eri_1_shsg The primitive integrals buffer.
/// @param idx_eri_1_sisf The primitive integrals buffer.
/// @param idx_eri_0_sisg The primitive integrals buffer.
/// @param idx_eri_1_sisg The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_sksg(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sksg,
                                  size_t idx_eri_0_shsg,
                                  size_t idx_eri_1_shsg,
                                  size_t idx_eri_1_sisf,
                                  size_t idx_eri_0_sisg,
                                  size_t idx_eri_1_sisg,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSKSG_hpp */
