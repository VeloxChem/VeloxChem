#ifndef ElectronRepulsionPrimRecSISH_hpp
#define ElectronRepulsionPrimRecSISH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SI|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sish The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sgsh The primitive integrals buffer.
/// @param idx_eri_1_sgsh The primitive integrals buffer.
/// @param idx_eri_1_shsg The primitive integrals buffer.
/// @param idx_eri_0_shsh The primitive integrals buffer.
/// @param idx_eri_1_shsh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sish(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sish,
                                       size_t                idx_eri_0_sgsh,
                                       size_t                idx_eri_1_sgsh,
                                       size_t                idx_eri_1_shsg,
                                       size_t                idx_eri_0_shsh,
                                       size_t                idx_eri_1_shsh,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSISH_hpp */
