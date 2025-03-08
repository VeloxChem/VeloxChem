#ifndef ThreeCenterElectronRepulsionPrimRecFSH_hpp
#define ThreeCenterElectronRepulsionPrimRecFSH_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [F|1/|r-r'||SH]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fsh The index of integral in primitive integrals buffer.
/// @param idx_eri_0_psh The primitive integrals buffer.
/// @param idx_eri_1_psh The primitive integrals buffer.
/// @param idx_eri_1_dsg The primitive integrals buffer.
/// @param idx_eri_1_dsh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_fsh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsh,
                                 size_t idx_eri_0_psh,
                                 size_t idx_eri_1_psh,
                                 size_t idx_eri_1_dsg,
                                 size_t idx_eri_1_dsh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecFSH_hpp */
