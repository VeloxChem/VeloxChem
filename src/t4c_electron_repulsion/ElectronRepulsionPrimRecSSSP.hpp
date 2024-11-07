#ifndef ElectronRepulsionPrimRecSSSP_hpp
#define ElectronRepulsionPrimRecSSSP_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes [SS|1/|r-r'||SP]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_sssp The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssss The primitive integrals buffer.
/// @param idx_eri_1_ssss The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_qd The vector of distances R(QD) = Q - D.
/// @param idx_wq The vector of distances R(WQ) = W - Q.
auto comp_prim_electron_repulsion_sssp(CSimdArray<double>& pbuffer,
                                       const size_t        idx_eri_0_sssp,
                                       size_t              idx_eri_0_ssss,
                                       size_t              idx_eri_1_ssss,
                                       CSimdArray<double>& factors,
                                       const size_t        idx_qd,
                                       const size_t        idx_wq) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionPrimRecSSSP_hpp */
