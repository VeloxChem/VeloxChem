#ifndef ThreeCenterElectronRepulsionContrRecXDF_hpp
#define ThreeCenterElectronRepulsionContrRecXDF_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||DF)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xdf The contracted integrals buffer.
/// @param idx_xpf The contracted integrals buffer.
/// @param idx_xpg The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xdf(CSimdArray<double>& cbuffer,
                                const size_t idx_xdf,
                                const size_t idx_xpf,
                                const size_t idx_xpg,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXDF_hpp */
