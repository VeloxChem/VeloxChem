#ifndef ElectronRepulsionContrRecPHXX_hpp
#define ElectronRepulsionContrRecPHXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (PH|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_phxx The contracted integrals buffer.
/// @param idx_shxx The contracted integrals buffer.
/// @param idx_sixx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_phxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_phxx,
                                     const size_t idx_shxx,
                                     const size_t idx_sixx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecPHXX_hpp */
