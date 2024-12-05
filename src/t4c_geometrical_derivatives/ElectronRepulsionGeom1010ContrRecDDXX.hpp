#ifndef ElectronRepulsionGeom1010ContrRecDDXX_hpp
#define ElectronRepulsionGeom1010ContrRecDDXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (DD|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_1010_ddxx The contracted integrals buffer.
/// @param idx_geom_0010_pdxx The contracted integrals buffer.
/// @param idx_geom_1010_pdxx The contracted integrals buffer.
/// @param idx_geom_1010_pfxx The contracted integrals buffer.
/// @param r_ab The Cartesian distance R(AB) = A - B.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom1010_hrr_electron_repulsion_ddxx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_ddxx,
                                              const size_t idx_geom_0010_pdxx,
                                              const size_t idx_geom_1010_pdxx,
                                              const size_t idx_geom_1010_pfxx,
                                              const TPoint<double>& r_ab,
                                              const int c_angmom,
                                              const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom1010ContrRecDDXX_hpp */
