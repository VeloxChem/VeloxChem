#ifndef ElectronRepulsionPrimRecSDSH_hpp
#define ElectronRepulsionPrimRecSDSH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SD|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sdsh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sssh The primitive integrals buffer.
/// @param idx_eri_1_sssh The primitive integrals buffer.
/// @param idx_eri_1_spsg The primitive integrals buffer.
/// @param idx_eri_0_spsh The primitive integrals buffer.
/// @param idx_eri_1_spsh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wp The vector of distances R(WP) = W - P.
/// @param r_pb The Cartesiandistances R(PB) = P - B.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sdsh(CSimdArray<double>&   pbuffer,
                                       const size_t          idx_eri_0_sdsh,
                                       size_t                idx_eri_0_sssh,
                                       size_t                idx_eri_1_sssh,
                                       size_t                idx_eri_1_spsg,
                                       size_t                idx_eri_0_spsh,
                                       size_t                idx_eri_1_spsh,
                                       CSimdArray<double>&   factors,
                                       const size_t          idx_wp,
                                       const TPoint<double>& r_pb,
                                       const double          a_exp,
                                       const double          b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSDSH_hpp */
