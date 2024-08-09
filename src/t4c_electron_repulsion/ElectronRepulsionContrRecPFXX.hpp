#ifndef ElectronRepulsionContrRecPFXX_hpp
#define ElectronRepulsionContrRecPFXX_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (PF|1/|r-r'|XX)  integrals for set of data buffers.
/// - Parameter contr_buffer_pfxx: the contracted integrals buffer.
/// - Parameter contr_buffer_sfxx: the contracted integrals buffer.
/// - Parameter contr_buffer_sgxx: the contracted integrals buffer.
/// - Parameter ab_x: the Cartesian X distance R(AB) = A - B.
/// - Parameter ab_y: the Cartesian Y distance R(AB) = A - B.
/// - Parameter ab_z: the Cartesian Z distance R(AB) = A - B.
/// - Parameter c_angmom: the angular momentum on center C.
/// - Parameter d_angmom: the angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_pfxx(CSimdArray<double>& contr_buffer_pfxx,
                                     const CSimdArray<double>& contr_buffer_sfxx,
                                     const CSimdArray<double>& contr_buffer_sgxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecPFXX_hpp */
