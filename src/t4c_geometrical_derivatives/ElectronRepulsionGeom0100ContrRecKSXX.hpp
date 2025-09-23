#ifndef ElectronRepulsionGeom0100ContrRecKSXX_hpp
#define ElectronRepulsionGeom0100ContrRecKSXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (KS|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_ksxx The contracted integrals buffer.
/// @param idx_isxx The contracted integrals buffer.
/// @param idx_geom_01_isxx The contracted integrals buffer.
/// @param idx_geom_01_ipxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_ksxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_ksxx,
                                            const size_t idx_isxx,
                                            const size_t idx_geom_01_isxx,
                                            const size_t idx_geom_01_ipxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecKSXX_hpp */
