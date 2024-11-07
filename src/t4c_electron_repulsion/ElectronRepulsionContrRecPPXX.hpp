#ifndef ElectronRepulsionContrRecPPXX_hpp
#define ElectronRepulsionContrRecPPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (PP|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_ppxx The contracted integrals buffer.
/// @param idx_spxx The contracted integrals buffer.
/// @param idx_sdxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto comp_bra_hrr_electron_repulsion_ppxx(CSimdArray<double>&   cbuffer,
                                          const size_t          idx_ppxx,
                                          const size_t          idx_spxx,
                                          const size_t          idx_sdxx,
                                          const TPoint<double>& r_ab,
                                          const int             c_angmom,
                                          const int             d_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionContrRecPPXX_hpp */
