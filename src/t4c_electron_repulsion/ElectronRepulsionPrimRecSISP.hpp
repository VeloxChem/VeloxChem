#ifndef ElectronRepulsionPrimRecSISP_hpp
#define ElectronRepulsionPrimRecSISP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SI|1/|r-r'||SP]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sisp The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sgsp The primitive integrals buffer.
/// @param idx_eri_1_sgsp The primitive integrals buffer.
/// @param idx_eri_1_shss The primitive integrals buffer.
/// @param idx_eri_0_shsp The primitive integrals buffer.
/// @param idx_eri_1_shsp The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sisp(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sisp,
                                       size_t                idx_eri_0_sgsp,
                                       size_t                idx_eri_1_sgsp,
                                       size_t                idx_eri_1_shss,
                                       size_t                idx_eri_0_shsp,
                                       size_t                idx_eri_1_shsp,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSISP_hpp */