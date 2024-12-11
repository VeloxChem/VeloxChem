#ifndef ElectronRepulsionGeom2000ContrRecSXXX_hpp
#define ElectronRepulsionGeom2000ContrRecSXXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (SX|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_20_sxxx The contracted integrals buffer.
/// @param idx_sxxx The contracted integrals buffer.
/// @param idx_dxxx The contracted integrals buffer.
/// @param b_angmom The angular momentum on center B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom20_hrr_electron_repulsion_sxxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_sxxx,
                                            const size_t idx_sxxx,
                                            const size_t idx_dxxx,
                                            const int b_angmom,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom2000ContrRecSXXX_hpp */
