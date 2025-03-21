#ifndef ElectronRepulsionPrimRecSHSF_hpp
#define ElectronRepulsionPrimRecSHSF_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SH|1/|r-r'||SF]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_shsf The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sfsf The primitive integrals buffer.
/// @param idx_eri_1_sfsf The primitive integrals buffer.
/// @param idx_eri_1_sgsd The primitive integrals buffer.
/// @param idx_eri_0_sgsf The primitive integrals buffer.
/// @param idx_eri_1_sgsf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_shsf(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_shsf,
                                       size_t                idx_eri_0_sfsf,
                                       size_t                idx_eri_1_sfsf,
                                       size_t                idx_eri_1_sgsd,
                                       size_t                idx_eri_0_sgsf,
                                       size_t                idx_eri_1_sgsf,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSHSF_hpp */
