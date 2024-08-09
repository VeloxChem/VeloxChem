#ifndef ElectronRepulsionPrimRecSSSP_hpp
#define ElectronRepulsionPrimRecSSSP_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SS|1/|r-r'||SP]  integrals for set of data buffers.
/// - Parameter prim_buffer_0_sssp: the primitive integrals buffer.
/// - Parameter prim_buffer_0_ssss: the primitive integrals buffer.
/// - Parameter prim_buffer_1_ssss: the primitive integrals buffer.
/// - Parameter qd_x: the vector of Cartesian X distances R(QD) = Q - D.
/// - Parameter qd_y: the vector of Cartesian Y distances R(QD) = Q - D.
/// - Parameter qd_z: the vector of Cartesian Z distances R(QD) = Q - D.
/// - Parameter wq_x: the vector of Cartesian X distances R(WQ) = W - Q.
/// - Parameter wq_y: the vector of Cartesian Y distances R(WQ) = W - Q.
/// - Parameter wq_z: the vector of Cartesian Z distances R(WQ) = W - Q.
auto
comp_prim_electron_repulsion_sssp(CSimdArray<double>& prim_buffer_0_sssp,
                                  const CSimdArray<double>& prim_buffer_0_ssss,
                                  const CSimdArray<double>& prim_buffer_1_ssss,
                                  const double* qd_x,
                                  const double* qd_y,
                                  const double* qd_z,
                                  const double* wq_x,
                                  const double* wq_y,
                                  const double* wq_z) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSSSP_hpp */
