#ifndef ElectronRepulsionGeom1000ContrRecGFXX_hpp
#define ElectronRepulsionGeom1000ContrRecGFXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GF|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_gfxx The contracted integrals buffer.
/// @param idx_ffxx The contracted integrals buffer.
/// @param idx_geom_10_ffxx The contracted integrals buffer.
/// @param idx_geom_10_fgxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom10_hrr_electron_repulsion_gfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_gfxx,
                                            const size_t idx_ffxx,
                                            const size_t idx_geom_10_ffxx,
                                            const size_t idx_geom_10_fgxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1000ContrRecGFXX_hpp */
