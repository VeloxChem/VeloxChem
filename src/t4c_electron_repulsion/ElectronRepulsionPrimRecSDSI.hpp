#ifndef ElectronRepulsionPrimRecSDSI_hpp
#define ElectronRepulsionPrimRecSDSI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SD|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sdsi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sssi The primitive integrals buffer.
/// @param idx_eri_1_sssi The primitive integrals buffer.
/// @param idx_eri_1_spsh The primitive integrals buffer.
/// @param idx_eri_0_spsi The primitive integrals buffer.
/// @param idx_eri_1_spsi The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sdsi(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sdsi,
                                       size_t                idx_eri_0_sssi,
                                       size_t                idx_eri_1_sssi,
                                       size_t                idx_eri_1_spsh,
                                       size_t                idx_eri_0_spsi,
                                       size_t                idx_eri_1_spsi,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSDSI_hpp */
