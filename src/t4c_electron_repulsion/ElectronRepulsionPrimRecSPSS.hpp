#ifndef ElectronRepulsionPrimRecSPSS_hpp
#define ElectronRepulsionPrimRecSPSS_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SP|1/|r-r'||SS]  integrals for set of data buffers.
/// - Parameter prim_buffer_0_spss: the primitive integrals buffer.
/// - Parameter prim_buffer_0_ssss: the primitive integrals buffer.
/// - Parameter prim_buffer_1_ssss: the primitive integrals buffer.
/// - Parameter pb_x: the Cartesian X distances R(PB) = P - B.
/// - Parameter pb_y: the Cartesian Y distances R(PB) = P - B.
/// - Parameter pb_z: the Cartesian Z distances R(PB) = P - B.
/// - Parameter wp_x: the vector of Cartesian X distances R(WP) = W - P.
/// - Parameter wp_y: the vector of Cartesian Y distances R(WP) = W - P.
/// - Parameter wp_z: the vector of Cartesian Z distances R(WP) = W - P.
auto
comp_prim_electron_repulsion_spss(CSimdArray<double>& prim_buffer_0_spss,
                                  const CSimdArray<double>& prim_buffer_0_ssss,
                                  const CSimdArray<double>& prim_buffer_1_ssss,
                                  const double pb_x,
                                  const double pb_y,
                                  const double pb_z,
                                  const double* wp_x,
                                  const double* wp_y,
                                  const double* wp_z) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSPSS_hpp */
