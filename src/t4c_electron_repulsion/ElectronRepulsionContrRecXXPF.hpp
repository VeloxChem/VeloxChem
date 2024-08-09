#ifndef ElectronRepulsionContrRecXXPF_hpp
#define ElectronRepulsionContrRecXXPF_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes (XX|1/|r-r'||PF)  integrals for set of data buffers.
/// - Parameter contr_buffer_xxpf: the contracted integrals buffer.
/// - Parameter contr_buffer_xxsf: the contracted integrals buffer.
/// - Parameter contr_buffer_xxsg: the contracted integrals buffer.
/// - Parameter cd_x: the vector of Cartesian X distances R(CD) = C - D.
/// - Parameter cd_y: the vector of Cartesian Y distances R(CD) = C - D.
/// - Parameter cd_z: the vector of Cartesian Z distances R(CD) = C - D.
/// - Parameter a_angmom: the angular momentum on center A.
/// - Parameter b_angmom: the angular momentum on center B.
auto
comp_ket_hrr_electron_repulsion_xxpf(CSimdArray<double>& contr_buffer_xxpf,
                                     const CSimdArray<double>& contr_buffer_xxsf,
                                     const CSimdArray<double>& contr_buffer_xxsg,
                                     const double* cd_x,
                                     const double* cd_y,
                                     const double* cd_z,
                                     const int a_angmom,
                                     const int b_angmom) -> void;
} // erirec namespace

#endif /* ElectronRepulsionContrRecXXPF_hpp */
