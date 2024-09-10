#ifndef ElectronRepulsionPrimRecSSSG_hpp
#define ElectronRepulsionPrimRecSSSG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SS|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sssg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_sssd The primitive integrals buffer.
/// @param idx_eri_1_sssd The primitive integrals buffer.
/// @param idx_eri_0_sssf The primitive integrals buffer.
/// @param idx_eri_1_sssf The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_qd The vector of distances R(QD) = Q - D.
/// @param idx_wq The vector of distances R(WQ) = W - Q.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto
comp_prim_electron_repulsion_sssg(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sssg,
                                  size_t idx_eri_0_sssd,
                                  size_t idx_eri_1_sssd,
                                  size_t idx_eri_0_sssf,
                                  size_t idx_eri_1_sssf,
                                  CSimdArray<double>& factors,
                                  const size_t idx_qd,
                                  const size_t idx_wq,
                                  const double a_exp,
                                  const double b_exp) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSSSG_hpp */
