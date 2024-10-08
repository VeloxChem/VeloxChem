#ifndef ElectronRepulsionContrRecFPXX_hpp
#define ElectronRepulsionContrRecFPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (FP|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_fpxx The contracted integrals buffer.
/// @param idx_dpxx The contracted integrals buffer.
/// @param idx_ddxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_fpxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_fpxx,
                                     const size_t idx_dpxx,
                                     const size_t idx_ddxx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecFPXX_hpp */
