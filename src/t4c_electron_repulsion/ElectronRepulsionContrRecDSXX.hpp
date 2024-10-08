#ifndef ElectronRepulsionContrRecDSXX_hpp
#define ElectronRepulsionContrRecDSXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (DS|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_dsxx The contracted integrals buffer.
/// @param idx_psxx The contracted integrals buffer.
/// @param idx_ppxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_dsxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_dsxx,
                                     const size_t idx_psxx,
                                     const size_t idx_ppxx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecDSXX_hpp */
