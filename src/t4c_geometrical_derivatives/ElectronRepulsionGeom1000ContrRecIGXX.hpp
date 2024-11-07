#ifndef ElectronRepulsionGeom1000ContrRecIGXX_hpp
#define ElectronRepulsionGeom1000ContrRecIGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (IG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_igxx The contracted integrals buffer.
/// @param idx_hgxx The contracted integrals buffer.
/// @param idx_geom_10_hgxx The contracted integrals buffer.
/// @param idx_geom_10_hhxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom10_hrr_electron_repulsion_igxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_igxx,
                                            const size_t idx_hgxx,
                                            const size_t idx_geom_10_hgxx,
                                            const size_t idx_geom_10_hhxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1000ContrRecIGXX_hpp */
