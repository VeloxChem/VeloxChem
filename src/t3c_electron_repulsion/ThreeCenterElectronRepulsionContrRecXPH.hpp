#ifndef ThreeCenterElectronRepulsionContrRecXPH_hpp
#define ThreeCenterElectronRepulsionContrRecXPH_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||PH)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xph The contracted integrals buffer.
/// @param idx_xsh The contracted integrals buffer.
/// @param idx_xsi The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xph(CSimdArray<double>& cbuffer,
                                const size_t idx_xph,
                                const size_t idx_xsh,
                                const size_t idx_xsi,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXPH_hpp */
