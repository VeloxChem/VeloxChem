#ifndef ElectronRepulsionContrRecXXDF_hpp
#define ElectronRepulsionContrRecXXDF_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (XX|1/|r-r'||DF)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxdf The contracted integrals buffer.
/// @param idx_xxpf The contracted integrals buffer.
/// @param idx_xxpg The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto comp_ket_hrr_electron_repulsion_xxdf(CSimdArray<double>&       cbuffer,
                                          const size_t              idx_xxdf,
                                          const size_t              idx_xxpf,
                                          const size_t              idx_xxpg,
                                          const CSimdArray<double>& factors,
                                          const size_t              idx_cd,
                                          const int                 a_angmom,
                                          const int                 b_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionContrRecXXDF_hpp */
