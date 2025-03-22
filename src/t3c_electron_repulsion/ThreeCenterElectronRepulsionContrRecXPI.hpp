#ifndef ThreeCenterElectronRepulsionContrRecXPI_hpp
#define ThreeCenterElectronRepulsionContrRecXPI_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes (X|1/|r-r'||PI)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xpi The contracted integrals buffer.
/// @param idx_xsi The contracted integrals buffer.
/// @param idx_xsk The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
auto
comp_hrr_electron_repulsion_xpi(CSimdArray<double>& cbuffer,
                                const size_t idx_xpi,
                                const size_t idx_xsi,
                                const size_t idx_xsk,
                                const CSimdArray<double>& factors,
                                const size_t idx_cd,
                                const int a_angmom) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionContrRecXPI_hpp */
