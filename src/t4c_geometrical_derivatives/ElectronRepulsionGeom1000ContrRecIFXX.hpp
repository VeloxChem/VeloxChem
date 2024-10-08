#ifndef ElectronRepulsionGeom1000ContrRecIFXX_hpp
#define ElectronRepulsionGeom1000ContrRecIFXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (IF|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_10_ifxx The contracted integrals buffer.
/// @param idx_hfxx The contracted integrals buffer.
/// @param idx_geom_10_hfxx The contracted integrals buffer.
/// @param idx_geom_10_hgxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom10_hrr_electron_repulsion_ifxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_ifxx,
                                            const size_t idx_hfxx,
                                            const size_t idx_geom_10_hfxx,
                                            const size_t idx_geom_10_hgxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1000ContrRecIFXX_hpp */
