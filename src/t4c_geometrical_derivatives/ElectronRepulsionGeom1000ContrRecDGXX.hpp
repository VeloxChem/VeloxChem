#ifndef ElectronRepulsionGeom1000ContrRecDGXX_hpp
#define ElectronRepulsionGeom1000ContrRecDGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (DG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_dgxx The contracted integrals buffer.
/// @param idx_pgxx The contracted integrals buffer.
/// @param idx_geom_10_pgxx The contracted integrals buffer.
/// @param idx_geom_10_phxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom10_hrr_electron_repulsion_dgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_dgxx,
                                            const size_t idx_pgxx,
                                            const size_t idx_geom_10_pgxx,
                                            const size_t idx_geom_10_phxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1000ContrRecDGXX_hpp */
