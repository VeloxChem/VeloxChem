#ifndef ElectronRepulsionPrimRecSDSL_hpp
#define ElectronRepulsionPrimRecSDSL_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SD|1/|r-r'||SL]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sdsl The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sssl The primitive integrals buffer.
/// @param idx_eri_1_sssl The primitive integrals buffer.
/// @param idx_eri_1_spsk The primitive integrals buffer.
/// @param idx_eri_0_spsl The primitive integrals buffer.
/// @param idx_eri_1_spsl The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_sdsl(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sdsl,
                                  size_t idx_eri_0_sssl,
                                  size_t idx_eri_1_sssl,
                                  size_t idx_eri_1_spsk,
                                  size_t idx_eri_0_spsl,
                                  size_t idx_eri_1_spsl,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSDSL_hpp */
