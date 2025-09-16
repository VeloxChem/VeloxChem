#ifndef ElectronRepulsionGeom0100ContrRecSFXX_hpp
#define ElectronRepulsionGeom0100ContrRecSFXX_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (SF|1/|r-r'|XX)  integral derivatives for set of data buffers.
/// @param cbuffer The contracted integrals buffer.
/// @param idx_geom_01_sfxx The contracted integrals buffer.
/// @param idx_sdxx The contracted integrals buffer.
/// @param idx_sgxx The contracted integrals buffer.
/// @param c_angmom The angular momentum on center C.
/// @param d_angmom The angular momentum on center D.
auto
comp_bra_geom01_hrr_electron_repulsion_sfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_sfxx,
                                            const size_t idx_sdxx,
                                            const size_t idx_sgxx,
                                            const int c_angmom,
                                            const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionGeom0100ContrRecSFXX_hpp */
