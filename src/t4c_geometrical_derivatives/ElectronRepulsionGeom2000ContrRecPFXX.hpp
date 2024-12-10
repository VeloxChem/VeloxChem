#ifndef ElectronRepulsionGeom2000ContrRecPFXX_hpp
#define ElectronRepulsionGeom2000ContrRecPFXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (PF|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_20_pfxx The contracted integrals buffer.
/// @param idx_geom_10_sfxx The contracted integrals buffer.
/// @param idx_geom_20_sfxx The contracted integrals buffer.
/// @param idx_geom_20_sgxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom20_hrr_electron_repulsion_pfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_pfxx,
                                            const size_t idx_geom_10_sfxx,
                                            const size_t idx_geom_20_sfxx,
                                            const size_t idx_geom_20_sgxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom2000ContrRecPFXX_hpp */
