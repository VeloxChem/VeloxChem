#ifndef ElectronRepulsionContrRecXXPG_hpp
#define ElectronRepulsionContrRecXXPG_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (XX|1/|r-r'||PG)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxpg The contracted integrals buffer.
/// @param pbuffer The Cartesian integrals buffer.
/// @param idx_xxsg The contracted integrals buffer.
/// @param idx_xxsh The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto comp_ket_hrr_electron_repulsion_xxpg(CSimdArray<double>&       cbuffer,
                                          const size_t              idx_xxpg,
                                          CSimdArray<double>&       pbuffer,
                                          const size_t              idx_xxsg,
                                          const size_t              idx_xxsh,
                                          const CSimdArray<double>& factors,
                                          const size_t              idx_cd,
                                          const int                 a_angmom,
                                          const int                 b_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionContrRecXXPG_hpp */
