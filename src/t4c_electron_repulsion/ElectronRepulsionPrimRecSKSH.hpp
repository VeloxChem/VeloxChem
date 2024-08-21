#ifndef ElectronRepulsionPrimRecSKSH_hpp
#define ElectronRepulsionPrimRecSKSH_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SK|1/|r-r'||SH]  integrals for set of data buffers.
/// - Parameter prim_buffer_0_sksh: the primitive integrals buffer.
/// - Parameter prim_buffer_0_shsh: the primitive integrals buffer.
/// - Parameter prim_buffer_1_shsh: the primitive integrals buffer.
/// - Parameter prim_buffer_1_sisg: the primitive integrals buffer.
/// - Parameter prim_buffer_0_sish: the primitive integrals buffer.
/// - Parameter prim_buffer_1_sish: the primitive integrals buffer.
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
comp_prim_electron_repulsion_sksh(CSimdArray<double>& prim_buffer_0_sksh,
                                  const CSimdArray<double>& prim_buffer_0_shsh,
                                  const CSimdArray<double>& prim_buffer_1_shsh,
                                  const CSimdArray<double>& prim_buffer_1_sisg,
                                  const CSimdArray<double>& prim_buffer_0_sish,
                                  const CSimdArray<double>& prim_buffer_1_sish,
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

#endif /* ElectronRepulsionPrimRecSKSH_hpp */