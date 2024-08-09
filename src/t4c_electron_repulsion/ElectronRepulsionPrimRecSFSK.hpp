#ifndef ElectronRepulsionPrimRecSFSK_hpp
#define ElectronRepulsionPrimRecSFSK_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SF|1/|r-r'||SK]  integrals for set of data buffers.
/// - Parameter prim_buffer_0_sfsk: the primitive integrals buffer.
/// - Parameter prim_buffer_0_spsk: the primitive integrals buffer.
/// - Parameter prim_buffer_1_spsk: the primitive integrals buffer.
/// - Parameter prim_buffer_1_sdsi: the primitive integrals buffer.
/// - Parameter prim_buffer_0_sdsk: the primitive integrals buffer.
/// - Parameter prim_buffer_1_sdsk: the primitive integrals buffer.
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
comp_prim_electron_repulsion_sfsk(CSimdArray<double>& prim_buffer_0_sfsk,
                                  const CSimdArray<double>& prim_buffer_0_spsk,
                                  const CSimdArray<double>& prim_buffer_1_spsk,
                                  const CSimdArray<double>& prim_buffer_1_sdsi,
                                  const CSimdArray<double>& prim_buffer_0_sdsk,
                                  const CSimdArray<double>& prim_buffer_1_sdsk,
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

#endif /* ElectronRepulsionPrimRecSFSK_hpp */
