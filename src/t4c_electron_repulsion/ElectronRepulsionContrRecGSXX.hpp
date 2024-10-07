#ifndef ElectronRepulsionContrRecGSXX_hpp
#define ElectronRepulsionContrRecGSXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GS|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_gsxx The contracted integrals buffer.
/// @param idx_fsxx The contracted integrals buffer.
/// @param idx_fpxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_gsxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_gsxx,
                                     const size_t idx_fsxx,
                                     const size_t idx_fpxx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecGSXX_hpp */
