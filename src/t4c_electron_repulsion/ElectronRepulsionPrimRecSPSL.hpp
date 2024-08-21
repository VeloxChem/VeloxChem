#ifndef ElectronRepulsionPrimRecSPSL_hpp
#define ElectronRepulsionPrimRecSPSL_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SP|1/|r-r'||SL]  integrals for set of data buffers.
/// - Parameter prim_buffer_0_spsl: the primitive integrals buffer.
/// - Parameter prim_buffer_1_sssk: the primitive integrals buffer.
/// - Parameter prim_buffer_0_sssl: the primitive integrals buffer.
/// - Parameter prim_buffer_1_sssl: the primitive integrals buffer.
/// - Parameter pb_x: the Cartesian X distances R(PB) = P - B.
/// - Parameter pb_y: the Cartesian Y distances R(PB) = P - B.
/// - Parameter pb_z: the Cartesian Z distances R(PB) = P - B.
/// - Parameter wp_x: the vector of Cartesian X distances R(WP) = W - P.
/// - Parameter wp_y: the vector of Cartesian Y distances R(WP) = W - P.
/// - Parameter wp_z: the vector of Cartesian Z distances R(WP) = W - P.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
/// - Parameter c_exps: the vector of exponents on center C.
/// - Parameter d_exps: the vector of exponents on center D.
auto
comp_prim_electron_repulsion_spsl(CSimdArray<double>& prim_buffer_0_spsl,
                                  const CSimdArray<double>& prim_buffer_1_sssk,
                                  const CSimdArray<double>& prim_buffer_0_sssl,
                                  const CSimdArray<double>& prim_buffer_1_sssl,
                                  const double pb_x,
                                  const double pb_y,
                                  const double pb_z,
                                  const double* wp_x,
                                  const double* wp_y,
                                  const double* wp_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSPSL_hpp */