#ifndef ElectronRepulsionPrimRecSGSL_hpp
#define ElectronRepulsionPrimRecSGSL_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SG|1/|r-r'||SL]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sgsl The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sdsl The primitive integrals buffer.
/// @param idx_eri_1_sdsl The primitive integrals buffer.
/// @param idx_eri_1_sfsk The primitive integrals buffer.
/// @param idx_eri_0_sfsl The primitive integrals buffer.
/// @param idx_eri_1_sfsl The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_sgsl(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sgsl,
                                  size_t idx_eri_0_sdsl,
                                  size_t idx_eri_1_sdsl,
                                  size_t idx_eri_1_sfsk,
                                  size_t idx_eri_0_sfsl,
                                  size_t idx_eri_1_sfsl,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSGSL_hpp */
