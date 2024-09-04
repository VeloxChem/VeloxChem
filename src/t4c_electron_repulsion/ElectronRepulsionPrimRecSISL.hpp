#ifndef ElectronRepulsionPrimRecSISL_hpp
#define ElectronRepulsionPrimRecSISL_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SI|1/|r-r'||SL]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sisl The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sgsl The primitive integrals buffer.
/// @param idx_eri_1_sgsl The primitive integrals buffer.
/// @param idx_eri_1_shsk The primitive integrals buffer.
/// @param idx_eri_0_shsl The primitive integrals buffer.
/// @param idx_eri_1_shsl The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_sisl(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sisl,
                                  size_t idx_eri_0_sgsl,
                                  size_t idx_eri_1_sgsl,
                                  size_t idx_eri_1_shsk,
                                  size_t idx_eri_0_shsl,
                                  size_t idx_eri_1_shsl,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSISL_hpp */
