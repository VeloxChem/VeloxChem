#ifndef ElectronRepulsionGeom0100ContrRecSKXX_hpp
#define ElectronRepulsionGeom0100ContrRecSKXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (SK|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_skxx The contracted integrals buffer.
/// @param idx_sixx The contracted integrals buffer.
/// @param idx_slxx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_skxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_skxx,
                                            const size_t idx_sixx,
                                            const size_t idx_slxx,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecSKXX_hpp */
