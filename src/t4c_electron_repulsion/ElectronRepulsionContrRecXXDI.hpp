#ifndef ElectronRepulsionContrRecXXDI_hpp
#define ElectronRepulsionContrRecXXDI_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (XX|1/|r-r'||DI)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxdi The contracted integrals buffer.
/// @param idx_xxpi The contracted integrals buffer.
/// @param idx_xxpk The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto comp_ket_hrr_electron_repulsion_xxdi(CSimdArray<double>&       cbuffer,
                                          const size_t              idx_xxdi,
                                          const size_t              idx_xxpi,
                                          const size_t              idx_xxpk,
                                          const CSimdArray<double>& factors,
                                          const size_t              idx_cd,
                                          const int                 a_angmom,
                                          const int                 b_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionContrRecXXDI_hpp */
