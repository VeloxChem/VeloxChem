#ifndef ElectronRepulsionGeom0100ContrRecSIXX_hpp
#define ElectronRepulsionGeom0100ContrRecSIXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (SI|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_sixx The contracted integrals buffer.
/// @param idx_shxx The contracted integrals buffer.
/// @param idx_skxx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_sixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_sixx,
                                            const size_t idx_shxx,
                                            const size_t idx_skxx,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecSIXX_hpp */
