#ifndef ThreeCenterElectronRepulsionContrRecXPS_hpp
#define ThreeCenterElectronRepulsionContrRecXPS_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||PS)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xps The contracted integrals buffer.
/// @param idx_xss The contracted integrals buffer.
/// @param idx_xsp The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xps(CSimdArray<double>& cbuffer,
                                const size_t idx_xps,
                                const size_t idx_xss,
                                const size_t idx_xsp,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXPS_hpp */
