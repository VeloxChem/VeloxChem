#ifndef ElectronRepulsionContrRecPFXX_hpp
#define ElectronRepulsionContrRecPFXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec {  // erirec namespace

/// Computes (PF|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_pfxx The contracted integrals buffer.
/// @param idx_sfxx The contracted integrals buffer.
/// @param idx_sgxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto comp_bra_hrr_electron_repulsion_pfxx(CSimdArray<double>&   cbuffer,
                                          const size_t          idx_pfxx,
                                          const size_t          idx_sfxx,
                                          const size_t          idx_sgxx,
                                          const TPoint<double>& r_ab,
                                          const int             c_angmom,
                                          const int             d_angmom) -> void;
}  // namespace erirec

#endif /* ElectronRepulsionContrRecPFXX_hpp */
