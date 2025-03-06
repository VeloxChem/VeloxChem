#ifndef ThreeCenterElectronRepulsionPrimRecPSI_hpp
#define ThreeCenterElectronRepulsionPrimRecPSI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [P|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_psi The index of integral in primitive integrals buffer.
/// @param idx_eri_1_ssh The primitive integrals buffer.
/// @param idx_eri_1_ssi The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_psi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psi,
                                 size_t idx_eri_1_ssh,
                                 size_t idx_eri_1_ssi,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecPSI_hpp */
