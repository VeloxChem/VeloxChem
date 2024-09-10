#ifndef ElectronRepulsionContrRecXXFG_hpp
#define ElectronRepulsionContrRecXXFG_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||FG)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxfg The contracted integrals buffer.
/// @param idx_xxdg The contracted integrals buffer.
/// @param idx_xxdh The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_hrr_electron_repulsion_xxfg(CSimdArray<double>& cbuffer,
                                     const size_t idx_xxfg,
                                     const size_t idx_xxdg,
                                     const size_t idx_xxdh,
                                     const CSimdArray<double>& factors,
                                     const size_t idx_cd,
                                     const int a_angmom,
                                     const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecXXFG_hpp */
