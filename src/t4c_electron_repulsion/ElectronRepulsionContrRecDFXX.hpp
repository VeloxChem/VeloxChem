#ifndef ElectronRepulsionContrRecDFXX_hpp
#define ElectronRepulsionContrRecDFXX_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (DF|1/|r-r'|XX)  integrals for set of data buffers.
/// - Parameter contr_buffer_dfxx: the contracted integrals buffer.
/// - Parameter contr_buffer_pfxx: the contracted integrals buffer.
/// - Parameter contr_buffer_pgxx: the contracted integrals buffer.
/// - Parameter ab_x: the Cartesian X distance R(AB) = A - B.
/// - Parameter ab_y: the Cartesian Y distance R(AB) = A - B.
/// - Parameter ab_z: the Cartesian Z distance R(AB) = A - B.
/// - Parameter c_angmom: the angular momentum on center C.
/// - Parameter d_angmom: the angular momentum on center D.
auto
comp_bra_hrr_electron_repulsion_dfxx(CSimdArray<double>& contr_buffer_dfxx,
                                     const CSimdArray<double>& contr_buffer_pfxx,
                                     const CSimdArray<double>& contr_buffer_pgxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecDFXX_hpp */
