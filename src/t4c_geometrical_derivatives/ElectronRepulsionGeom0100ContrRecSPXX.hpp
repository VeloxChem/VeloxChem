#ifndef ElectronRepulsionGeom0100ContrRecSPXX_hpp
#define ElectronRepulsionGeom0100ContrRecSPXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (SP|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_spxx The contracted integrals buffer.
/// @param idx_ssxx The contracted integrals buffer.
/// @param idx_sdxx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_spxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_spxx,
                                            const size_t idx_ssxx,
                                            const size_t idx_sdxx,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecSPXX_hpp */
