#ifndef ElectronRepulsionPrimRecSDSD_hpp
#define ElectronRepulsionPrimRecSDSD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SD|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sdsd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sssd The primitive integrals buffer.
/// @param idx_eri_1_sssd The primitive integrals buffer.
/// @param idx_eri_1_spsp The primitive integrals buffer.
/// @param idx_eri_0_spsd The primitive integrals buffer.
/// @param idx_eri_1_spsd The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sdsd(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sdsd,
                                       size_t                idx_eri_0_sssd,
                                       size_t                idx_eri_1_sssd,
                                       size_t                idx_eri_1_spsp,
                                       size_t                idx_eri_0_spsd,
                                       size_t                idx_eri_1_spsd,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSDSD_hpp */
