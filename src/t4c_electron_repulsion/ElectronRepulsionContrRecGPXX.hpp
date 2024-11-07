#ifndef ElectronRepulsionContrRecGPXX_hpp
#define ElectronRepulsionContrRecGPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GP|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_gpxx The contracted integrals buffer.
/// @param idx_fpxx The contracted integrals buffer.
/// @param idx_fdxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_gpxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_gpxx,
                                     const size_t idx_fpxx,
                                     const size_t idx_fdxx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecGPXX_hpp */
