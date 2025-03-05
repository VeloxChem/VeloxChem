#ifndef ThreeCenterElectronRepulsionContrRecXDI_hpp
#define ThreeCenterElectronRepulsionContrRecXDI_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||DI)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xdi The contracted integrals buffer.
/// @param idx_xpi The contracted integrals buffer.
/// @param idx_xpk The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xdi(CSimdArray<double>& cbuffer,
                                const size_t idx_xdi,
                                const size_t idx_xpi,
                                const size_t idx_xpk,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXDI_hpp */
