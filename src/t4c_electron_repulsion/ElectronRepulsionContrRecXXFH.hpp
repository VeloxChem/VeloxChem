#ifndef ElectronRepulsionContrRecXXFH_hpp
#define ElectronRepulsionContrRecXXFH_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||FH)  integrals for set of data buffers.
/// - Parameter contr_buffer_xxfh: the contracted integrals buffer.
/// - Parameter contr_buffer_xxdh: the contracted integrals buffer.
/// - Parameter contr_buffer_xxdi: the contracted integrals buffer.
/// - Parameter cd_x: the vector of Cartesian X distances R(CD) = C - D.
/// - Parameter cd_y: the vector of Cartesian Y distances R(CD) = C - D.
/// - Parameter cd_z: the vector of Cartesian Z distances R(CD) = C - D.
/// - Parameter a_angmom: the angular momentum on center A.
/// - Parameter b_angmom: the angular momentum on center B.
auto
comp_ket_hrr_electron_repulsion_xxfh(CSimdArray<double>& contr_buffer_xxfh,
                                     const CSimdArray<double>& contr_buffer_xxdh,
                                     const CSimdArray<double>& contr_buffer_xxdi,
                                     const double* cd_x,
                                     const double* cd_y,
                                     const double* cd_z,
                                     const int a_angmom,
                                     const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecXXFH_hpp */
