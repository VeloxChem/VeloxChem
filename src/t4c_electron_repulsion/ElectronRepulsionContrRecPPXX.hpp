#ifndef ElectronRepulsionContrRecPPXX_hpp
#define ElectronRepulsionContrRecPPXX_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (PP|1/|r-r'|XX)  integrals for set of data buffers.
/// - Parameter contr_buffer_ppxx: the contracted integrals buffer.
/// - Parameter contr_buffer_spxx: the contracted integrals buffer.
/// - Parameter contr_buffer_sdxx: the contracted integrals buffer.
/// - Parameter ab_x: the Cartesian X distance R(AB) = A - B.
/// - Parameter ab_y: the Cartesian Y distance R(AB) = A - B.
/// - Parameter ab_z: the Cartesian Z distance R(AB) = A - B.
/// - Parameter c_angmom: the angular momentum on center C.
/// - Parameter d_angmom: the angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_ppxx(CSimdArray<double>& contr_buffer_ppxx,
                                     const CSimdArray<double>& contr_buffer_spxx,
                                     const CSimdArray<double>& contr_buffer_sdxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecPPXX_hpp */
