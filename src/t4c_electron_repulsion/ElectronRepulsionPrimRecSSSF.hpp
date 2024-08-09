#ifndef ElectronRepulsionPrimRecSSSF_hpp
#define ElectronRepulsionPrimRecSSSF_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SS|1/|r-r'||SF]  integrals for set of data buffers.
/// - Parameter prim_buffer_0_sssf: the primitive integrals buffer.
/// - Parameter prim_buffer_0_sssp: the primitive integrals buffer.
/// - Parameter prim_buffer_1_sssp: the primitive integrals buffer.
/// - Parameter prim_buffer_0_sssd: the primitive integrals buffer.
/// - Parameter prim_buffer_1_sssd: the primitive integrals buffer.
/// - Parameter qd_x: the vector of Cartesian X distances R(QD) = Q - D.
/// - Parameter qd_y: the vector of Cartesian Y distances R(QD) = Q - D.
/// - Parameter qd_z: the vector of Cartesian Z distances R(QD) = Q - D.
/// - Parameter wq_x: the vector of Cartesian X distances R(WQ) = W - Q.
/// - Parameter wq_y: the vector of Cartesian Y distances R(WQ) = W - Q.
/// - Parameter wq_z: the vector of Cartesian Z distances R(WQ) = W - Q.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
/// - Parameter c_exps: the vector of exponents on center C.
/// - Parameter d_exps: the vector of exponents on center D.
auto
comp_prim_electron_repulsion_sssf(CSimdArray<double>& prim_buffer_0_sssf,
                                  const CSimdArray<double>& prim_buffer_0_sssp,
                                  const CSimdArray<double>& prim_buffer_1_sssp,
                                  const CSimdArray<double>& prim_buffer_0_sssd,
                                  const CSimdArray<double>& prim_buffer_1_sssd,
                                  const double* qd_x,
                                  const double* qd_y,
                                  const double* qd_z,
                                  const double* wq_x,
                                  const double* wq_y,
                                  const double* wq_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSSSF_hpp */
