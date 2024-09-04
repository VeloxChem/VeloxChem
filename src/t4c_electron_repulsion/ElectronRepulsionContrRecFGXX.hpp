#ifndef ElectronRepulsionContrRecFGXX_hpp
#define ElectronRepulsionContrRecFGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (FG|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_fgxx The contracted integrals buffer.
/// @param idx_dgxx The contracted integrals buffer.
/// @param idx_dhxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_fgxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_fgxx,
                                     const size_t idx_dgxx,
                                     const size_t idx_dhxx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecFGXX_hpp */
