#ifndef ElectronRepulsionContrRecFHXX_hpp
#define ElectronRepulsionContrRecFHXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (FH|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_fhxx The contracted integrals buffer.
/// @param idx_dhxx The contracted integrals buffer.
/// @param idx_dixx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto comp_bra_hrr_electron_repulsion_fhxx(CSimdArray<double>&   cbuffer,
                                          const size_t          idx_fhxx,
                                          const size_t          idx_dhxx,
                                          const size_t          idx_dixx,
                                          const TPoint<double>& r_ab,
                                          const int             c_angmom,
                                          const int             d_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionContrRecFHXX_hpp */
