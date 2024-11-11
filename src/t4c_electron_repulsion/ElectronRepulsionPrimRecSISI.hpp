#ifndef ElectronRepulsionPrimRecSISI_hpp
#define ElectronRepulsionPrimRecSISI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SI|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sisi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sgsi The primitive integrals buffer.
/// @param idx_eri_1_sgsi The primitive integrals buffer.
/// @param idx_eri_1_shsh The primitive integrals buffer.
/// @param idx_eri_0_shsi The primitive integrals buffer.
/// @param idx_eri_1_shsi The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sisi(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sisi,
                                       size_t                idx_eri_0_sgsi,
                                       size_t                idx_eri_1_sgsi,
                                       size_t                idx_eri_1_shsh,
                                       size_t                idx_eri_0_shsi,
                                       size_t                idx_eri_1_shsi,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSISI_hpp */
