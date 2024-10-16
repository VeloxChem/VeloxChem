#ifndef ElectronRepulsionContrRecXXFF_hpp
#define ElectronRepulsionContrRecXXFF_hpp

#include <cstddef>

#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (XX|1/|r-r'||FF)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_xxff The contracted integrals buffer.
/// @param idx_xxdf The contracted integrals buffer.
/// @param idx_xxdg The contracted integrals buffer.
/// @param factors The factors buffer.
/// @param idx_cd The vector of distances R(CD) = C - D.
/// @param a_angmom The angular momentum on center A.
/// @param b_angmom The angular momentum on center B.
auto comp_ket_hrr_electron_repulsion_xxff(CSimdArray<double>&       cbuffer,
                                          const size_t              idx_xxff,
                                          const size_t              idx_xxdf,
                                          const size_t              idx_xxdg,
                                          const CSimdArray<double>& factors,
                                          const size_t              idx_cd,
                                          const int                 a_angmom,
                                          const int                 b_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionContrRecXXFF_hpp */
