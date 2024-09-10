#ifndef ElectronRepulsionContrRecDFXX_hpp
#define ElectronRepulsionContrRecDFXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (DF|1/|r-r'|XX)  integrals for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_dfxx The contracted integrals buffer.
/// @param idx_pfxx The contracted integrals buffer.
/// @param idx_pgxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_dfxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_dfxx,
                                     const size_t idx_pfxx,
                                     const size_t idx_pgxx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecDFXX_hpp */
