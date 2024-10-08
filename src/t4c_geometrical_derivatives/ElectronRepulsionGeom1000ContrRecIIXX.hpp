#ifndef ElectronRepulsionGeom1000ContrRecIIXX_hpp
#define ElectronRepulsionGeom1000ContrRecIIXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (II|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_iixx The contracted integrals buffer.
/// @param idx_hixx The contracted integrals buffer.
/// @param idx_geom_10_hixx The contracted integrals buffer.
/// @param idx_geom_10_hkxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom10_hrr_electron_repulsion_iixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_iixx,
                                            const size_t idx_hixx,
                                            const size_t idx_geom_10_hixx,
                                            const size_t idx_geom_10_hkxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1000ContrRecIIXX_hpp */
