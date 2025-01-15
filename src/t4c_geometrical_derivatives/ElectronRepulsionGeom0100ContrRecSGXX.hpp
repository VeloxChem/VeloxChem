#ifndef ElectronRepulsionGeom0100ContrRecSGXX_hpp
#define ElectronRepulsionGeom0100ContrRecSGXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (SG|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_sgxx The contracted integrals buffer.
/// @param idx_sfxx The contracted integrals buffer.
/// @param idx_shxx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_sgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_sgxx,
                                            const size_t idx_sfxx,
                                            const size_t idx_shxx,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecSGXX_hpp */
