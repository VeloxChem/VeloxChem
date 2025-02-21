#ifndef ThreeCenterElectronRepulsionContrRecXGP_hpp
#define ThreeCenterElectronRepulsionContrRecXGP_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||GP)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xgp The contracted integrals buffer.
/// @param idx_xfp The contracted integrals buffer.
/// @param idx_xfd The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xgp(CSimdArray<double>& cbuffer,
                                const size_t idx_xgp,
                                const size_t idx_xfp,
                                const size_t idx_xfd,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXGP_hpp */
