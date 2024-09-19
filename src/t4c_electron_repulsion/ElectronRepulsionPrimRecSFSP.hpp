#ifndef ElectronRepulsionPrimRecSFSP_hpp
#define ElectronRepulsionPrimRecSFSP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SF|1/|r-r'||SP]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sfsp The index of integral in primitive integrals buffer.
/// @param idx_eri_0_spsp The primitive integrals buffer.
/// @param idx_eri_1_spsp The primitive integrals buffer.
/// @param idx_eri_1_sdss The primitive integrals buffer.
/// @param idx_eri_0_sdsp The primitive integrals buffer.
/// @param idx_eri_1_sdsp The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sfsp(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sfsp,
                                       size_t                idx_eri_0_spsp,
                                       size_t                idx_eri_1_spsp,
                                       size_t                idx_eri_1_sdss,
                                       size_t                idx_eri_0_sdsp,
                                       size_t                idx_eri_1_sdsp,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSFSP_hpp */
