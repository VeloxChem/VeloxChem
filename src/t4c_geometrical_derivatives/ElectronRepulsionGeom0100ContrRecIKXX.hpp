#ifndef ElectronRepulsionGeom0100ContrRecIKXX_hpp
#define ElectronRepulsionGeom0100ContrRecIKXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (IK|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_ikxx The contracted integrals buffer.
/// @param idx_hkxx The contracted integrals buffer.
/// @param idx_geom_01_hkxx The contracted integrals buffer.
/// @param idx_geom_01_hlxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_ikxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_ikxx,
                                            const size_t idx_hkxx,
                                            const size_t idx_geom_01_hkxx,
                                            const size_t idx_geom_01_hlxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecIKXX_hpp */
