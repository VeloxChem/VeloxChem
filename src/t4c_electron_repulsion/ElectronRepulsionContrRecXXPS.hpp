#ifndef ElectronRepulsionContrRecXXPS_hpp
#define ElectronRepulsionContrRecXXPS_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||PS)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxps The contracted integrals buffer.
/// @param pbuffer The Cartesian integrals buffer.
/// @param idx_xxss The contracted integrals buffer.
/// @param idx_xxsp The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_hrr_electron_repulsion_xxps(CSimdArray<double>& cbuffer,
                                     const size_t idx_xxps,
                                     CSimdArray<double>& pbuffer,
                                     const size_t idx_xxss,
                                     const size_t idx_xxsp,
                                     const CSimdArray<double>& factors,
                                     const size_t idx_cd,
                                     const int a_angmom,
                                     const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecXXPS_hpp */
