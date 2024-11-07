#ifndef ElectronRepulsionPrimRecSFSK_hpp
#define ElectronRepulsionPrimRecSFSK_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SF|1/|r-r'||SK]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sfsk The index of integral in primitive integrals buffer.
/// @param idx_eri_0_spsk The primitive integrals buffer.
/// @param idx_eri_1_spsk The primitive integrals buffer.
/// @param idx_eri_1_sdsi The primitive integrals buffer.
/// @param idx_eri_0_sdsk The primitive integrals buffer.
/// @param idx_eri_1_sdsk The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sfsk(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sfsk,
                                       size_t                idx_eri_0_spsk,
                                       size_t                idx_eri_1_spsk,
                                       size_t                idx_eri_1_sdsi,
                                       size_t                idx_eri_0_sdsk,
                                       size_t                idx_eri_1_sdsk,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSFSK_hpp */
