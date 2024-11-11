#ifndef ElectronRepulsionPrimRecSHSP_hpp
#define ElectronRepulsionPrimRecSHSP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SH|1/|r-r'||SP]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_shsp The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sfsp The primitive integrals buffer.
/// @param idx_eri_1_sfsp The primitive integrals buffer.
/// @param idx_eri_1_sgss The primitive integrals buffer.
/// @param idx_eri_0_sgsp The primitive integrals buffer.
/// @param idx_eri_1_sgsp The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_shsp(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_shsp,
                                       size_t                idx_eri_0_sfsp,
                                       size_t                idx_eri_1_sfsp,
                                       size_t                idx_eri_1_sgss,
                                       size_t                idx_eri_0_sgsp,
                                       size_t                idx_eri_1_sgsp,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSHSP_hpp */
