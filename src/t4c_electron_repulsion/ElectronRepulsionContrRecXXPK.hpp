#ifndef ElectronRepulsionContrRecXXPK_hpp
#define ElectronRepulsionContrRecXXPK_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||PK)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxpk The contracted integrals buffer.
/// @param pbuffer The Cartesian integrals buffer.
/// @param idx_xxsk The contracted integrals buffer.
/// @param idx_xxsl The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto
comp_ket_hrr_electron_repulsion_xxpk(CSimdArray<double>& cbuffer,
                                     const size_t idx_xxpk,
                                     CSimdArray<double>& pbuffer,
                                     const size_t idx_xxsk,
                                     const size_t idx_xxsl,
                                     const CSimdArray<double>& factors,
                                     const size_t idx_cd,
                                     const int a_angmom,
                                     const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecXXPK_hpp */
