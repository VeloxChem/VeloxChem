#ifndef ElectronRepulsionContrRecXXPD_hpp
#define ElectronRepulsionContrRecXXPD_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (XX|1/|r-r'||PD)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxpd The contracted integrals buffer.
/// @param pbuffer The Cartesian integrals buffer.
/// @param idx_xxsd The contracted integrals buffer.
/// @param idx_xxsf The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto comp_ket_hrr_electron_repulsion_xxpd(CSimdArray<double>&       cbuffer,
                                          const size_t              idx_xxpd,
                                          CSimdArray<double>&       pbuffer,
                                          const size_t              idx_xxsd,
                                          const size_t              idx_xxsf,
                                          const CSimdArray<double>& factors,
                                          const size_t              idx_cd,
                                          const int                 a_angmom,
                                          const int                 b_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionContrRecXXPD_hpp */
