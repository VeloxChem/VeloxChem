#ifndef ElectricDipoleMomentumPrimRecSS_hpp
#define ElectricDipoleMomentumPrimRecSS_hpp

#include "SimdArray.hpp"

namespace diprec {  // diprec namespace

/// @brief Computes primitive [S|r|S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_dip_ss The index of integral in primitive integrals buffer.
/// @param idx_ovl_ss The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpc The vector of distances R(PC) = P - C.
auto comp_prim_electric_dipole_momentum_ss(CSimdArray<double>&       pbuffer,
                                           const size_t              idx_dip_ss,
                                           const size_t              idx_ovl_ss,
                                           const CSimdArray<double>& factors,
                                           const size_t              idx_rpc) -> void;
}  // namespace diprec

#endif /* ElectricDipoleMomentumPrimRecSS_hpp */
