#ifndef ElectronRepulsionPrimRecSPSS_hpp
#define ElectronRepulsionPrimRecSPSS_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SP|1/|r-r'||SS]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_spss The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssss The primitive integrals buffer.
/// @param idx_eri_1_ssss The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
auto
comp_prim_electron_repulsion_spss(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_spss,
                                  size_t idx_eri_0_ssss,
                                  size_t idx_eri_1_ssss,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSPSS_hpp */
