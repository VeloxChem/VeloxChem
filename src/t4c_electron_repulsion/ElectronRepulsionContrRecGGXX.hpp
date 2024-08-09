#ifndef ElectronRepulsionContrRecGGXX_hpp
#define ElectronRepulsionContrRecGGXX_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (GG|1/|r-r'|XX)  integrals for set of data buffers.
/// - Parameter contr_buffer_ggxx: the contracted integrals buffer.
/// - Parameter contr_buffer_fgxx: the contracted integrals buffer.
/// - Parameter contr_buffer_fhxx: the contracted integrals buffer.
/// - Parameter ab_x: the Cartesian X distance R(AB) = A - B.
/// - Parameter ab_y: the Cartesian Y distance R(AB) = A - B.
/// - Parameter ab_z: the Cartesian Z distance R(AB) = A - B.
/// - Parameter c_angmom: the angular momentum on center C.
/// - Parameter d_angmom: the angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_ggxx(CSimdArray<double>& contr_buffer_ggxx,
                                     const CSimdArray<double>& contr_buffer_fgxx,
                                     const CSimdArray<double>& contr_buffer_fhxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecGGXX_hpp */
