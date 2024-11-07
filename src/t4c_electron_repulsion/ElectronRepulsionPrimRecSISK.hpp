#ifndef ElectronRepulsionPrimRecSISK_hpp
#define ElectronRepulsionPrimRecSISK_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SI|1/|r-r'||SK]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sisk The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sgsk The primitive integrals buffer.
/// @param idx_eri_1_sgsk The primitive integrals buffer.
/// @param idx_eri_1_shsi The primitive integrals buffer.
/// @param idx_eri_0_shsk The primitive integrals buffer.
/// @param idx_eri_1_shsk The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sisk(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sisk,
                                       size_t                idx_eri_0_sgsk,
                                       size_t                idx_eri_1_sgsk,
                                       size_t                idx_eri_1_shsi,
                                       size_t                idx_eri_0_shsk,
                                       size_t                idx_eri_1_shsk,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSISK_hpp */
