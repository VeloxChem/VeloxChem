#ifndef ElectronRepulsionContrRecFSXX_hpp
#define ElectronRepulsionContrRecFSXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (FS|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_fsxx The contracted integrals buffer.
/// @param idx_dsxx The contracted integrals buffer.
/// @param idx_dpxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_fsxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_fsxx,
                                     const size_t idx_dsxx,
                                     const size_t idx_dpxx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecFSXX_hpp */
