#ifndef ElectronRepulsionGeom1100ContrRecIIXX_hpp
#define ElectronRepulsionGeom1100ContrRecIIXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (II|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_11_iixx The contracted integrals buffer.
/// @param idx_geom_01_hixx The contracted integrals buffer.
/// @param idx_geom_10_hixx The contracted integrals buffer.
/// @param idx_geom_11_hixx The contracted integrals buffer.
/// @param idx_geom_11_hkxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom11_hrr_electron_repulsion_iixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_iixx,
                                            const size_t idx_geom_01_hixx,
                                            const size_t idx_geom_10_hixx,
                                            const size_t idx_geom_11_hixx,
                                            const size_t idx_geom_11_hkxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1100ContrRecIIXX_hpp */
