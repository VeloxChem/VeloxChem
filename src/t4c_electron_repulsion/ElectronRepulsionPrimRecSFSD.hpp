#ifndef ElectronRepulsionPrimRecSFSD_hpp
#define ElectronRepulsionPrimRecSFSD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SF|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sfsd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_spsd The primitive integrals buffer.
/// @param idx_eri_1_spsd The primitive integrals buffer.
/// @param idx_eri_1_sdsp The primitive integrals buffer.
/// @param idx_eri_0_sdsd The primitive integrals buffer.
/// @param idx_eri_1_sdsd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sfsd(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sfsd,
                                       size_t                idx_eri_0_spsd,
                                       size_t                idx_eri_1_spsd,
                                       size_t                idx_eri_1_sdsp,
                                       size_t                idx_eri_0_sdsd,
                                       size_t                idx_eri_1_sdsd,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSFSD_hpp */
