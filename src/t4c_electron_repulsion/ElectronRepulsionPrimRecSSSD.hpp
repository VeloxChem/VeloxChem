#ifndef ElectronRepulsionPrimRecSSSD_hpp
#define ElectronRepulsionPrimRecSSSD_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SS|1/|r-r'||SD]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sssd The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssss The primitive integrals buffer.
/// @param idx_eri_1_ssss The primitive integrals buffer.
/// @param idx_eri_0_sssp The primitive integrals buffer.
/// @param idx_eri_1_sssp The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_qd The vector of distances R(QD) = Q - D.
/// @param idx_wq The vector of distances R(WQ) = W - Q.
/// @param a_exp The exponent on center A.
/// @param b_exp The exponent on center B.
auto comp_prim_electron_repulsion_sssd(CSimdArray<double>& pbuffer,
                                       const size_t        idx_eri_0_sssd,
                                       size_t              idx_eri_0_ssss,
                                       size_t              idx_eri_1_ssss,
                                       size_t              idx_eri_0_sssp,
                                       size_t              idx_eri_1_sssp,
                                       CSimdArray<double>& factors,
                                       const size_t        idx_qd,
                                       const size_t        idx_wq,
                                       const double        a_exp,
                                       const double        b_exp) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSSSD_hpp */
