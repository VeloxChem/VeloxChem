#ifndef ThreeCenterElectronRepulsionContrRecXPD_hpp
#define ThreeCenterElectronRepulsionContrRecXPD_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||PD)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xpd The contracted integrals buffer.
/// @param idx_xsd The contracted integrals buffer.
/// @param idx_xsf The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xpd(CSimdArray<double>& cbuffer,
                                const size_t idx_xpd,
                                const size_t idx_xsd,
                                const size_t idx_xsf,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXPD_hpp */
