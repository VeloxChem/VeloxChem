#ifndef ElectronRepulsionContrRecXXPF_hpp
#define ElectronRepulsionContrRecXXPF_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (XX|1/|r-r'||PF)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxpf The contracted integrals buffer.
/// @param pbuffer The Cartesian integrals buffer.
/// @param idx_xxsf The contracted integrals buffer.
/// @param idx_xxsg The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto comp_ket_hrr_electron_repulsion_xxpf(CSimdArray<double>&       cbuffer,
                                          const size_t              idx_xxpf,
                                          CSimdArray<double>&       pbuffer,
                                          const size_t              idx_xxsf,
                                          const size_t              idx_xxsg,
                                          const CSimdArray<double>& factors,
                                          const size_t              idx_cd,
                                          const int                 a_angmom,
                                          const int                 b_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionContrRecXXPF_hpp */
