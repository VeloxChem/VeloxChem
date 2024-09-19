#ifndef ElectronRepulsionPrimRecSFSG_hpp
#define ElectronRepulsionPrimRecSFSG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SF|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sfsg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_spsg The primitive integrals buffer.
/// @param idx_eri_1_spsg The primitive integrals buffer.
/// @param idx_eri_1_sdsf The primitive integrals buffer.
/// @param idx_eri_0_sdsg The primitive integrals buffer.
/// @param idx_eri_1_sdsg The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sfsg(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sfsg,
                                       size_t                idx_eri_0_spsg,
                                       size_t                idx_eri_1_spsg,
                                       size_t                idx_eri_1_sdsf,
                                       size_t                idx_eri_0_sdsg,
                                       size_t                idx_eri_1_sdsg,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSFSG_hpp */
