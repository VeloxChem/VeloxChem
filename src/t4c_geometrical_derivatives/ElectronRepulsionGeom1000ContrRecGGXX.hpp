#ifndef ElectronRepulsionGeom1000ContrRecGGXX_hpp
#define ElectronRepulsionGeom1000ContrRecGGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_ggxx The contracted integrals buffer.
/// @param idx_fgxx The contracted integrals buffer.
/// @param idx_geom_10_fgxx The contracted integrals buffer.
/// @param idx_geom_10_fhxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom10_hrr_electron_repulsion_ggxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_ggxx,
                                            const size_t idx_fgxx,
                                            const size_t idx_geom_10_fgxx,
                                            const size_t idx_geom_10_fhxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1000ContrRecGGXX_hpp */
