#ifndef ThreeCenterElectronRepulsionContrRecXPG_hpp
#define ThreeCenterElectronRepulsionContrRecXPG_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||PG)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xpg The contracted integrals buffer.
/// @param idx_xsg The contracted integrals buffer.
/// @param idx_xsh The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xpg(CSimdArray<double>& cbuffer,
                                const size_t idx_xpg,
                                const size_t idx_xsg,
                                const size_t idx_xsh,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXPG_hpp */
