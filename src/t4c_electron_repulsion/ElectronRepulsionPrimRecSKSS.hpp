#ifndef ElectronRepulsionPrimRecSKSS_hpp
#define ElectronRepulsionPrimRecSKSS_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SK|1/|r-r'||SS]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_skss The index of integral in primitive integrals buffer.
/// @param idx_eri_0_shss The primitive integrals buffer.
/// @param idx_eri_1_shss The primitive integrals buffer.
/// @param idx_eri_0_siss The primitive integrals buffer.
/// @param idx_eri_1_siss The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_skss(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_skss,
                                       size_t                idx_eri_0_shss,
                                       size_t                idx_eri_1_shss,
                                       size_t                idx_eri_0_siss,
                                       size_t                idx_eri_1_siss,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSKSS_hpp */