#ifndef ElectronRepulsionContrRecXXFH_hpp
#define ElectronRepulsionContrRecXXFH_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||FH)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxfh The contracted integrals buffer.
/// @param pbuffer The Cartesian integrals buffer.
/// @param idx_xxdh The contracted integrals buffer.
/// @param idx_xxdi The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_hrr_electron_repulsion_xxfh(CSimdArray<double>& cbuffer,
                                     const size_t idx_xxfh,
                                     CSimdArray<double>& pbuffer,
                                     const size_t idx_xxdh,
                                     const size_t idx_xxdi,
                                     const CSimdArray<double>& factors,
                                     const size_t idx_cd,
                                     const int a_angmom,
                                     const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecXXFH_hpp */
