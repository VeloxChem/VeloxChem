#ifndef ElectronRepulsionContrRecXXDS_hpp
#define ElectronRepulsionContrRecXXDS_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||DS)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxds The contracted integrals buffer.
/// @param idx_xxps The contracted integrals buffer.
/// @param idx_xxpp The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_hrr_electron_repulsion_xxds(CSimdArray<double>& cbuffer,
                                     const size_t idx_xxds,
                                     const size_t idx_xxps,
                                     const size_t idx_xxpp,
                                     const CSimdArray<double>& factors,
                                     const size_t idx_cd,
                                     const int a_angmom,
                                     const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecXXDS_hpp */
