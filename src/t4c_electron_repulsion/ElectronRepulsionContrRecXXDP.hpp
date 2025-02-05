#ifndef ElectronRepulsionContrRecXXDP_hpp
#define ElectronRepulsionContrRecXXDP_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||DP)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxdp The contracted integrals buffer.
/// @param idx_xxpp The contracted integrals buffer.
/// @param idx_xxpd The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_hrr_electron_repulsion_xxdp(CSimdArray<double>& cbuffer,
                                     const size_t idx_xxdp,
                                     const size_t idx_xxpp,
                                     const size_t idx_xxpd,
                                     const CSimdArray<double>& factors,
                                     const size_t idx_cd,
                                     const int a_angmom,
                                     const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecXXDP_hpp */
