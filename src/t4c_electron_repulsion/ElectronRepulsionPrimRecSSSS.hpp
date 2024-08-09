#ifndef ElectronRepulsionPrimRecSSSS_hpp
#define ElectronRepulsionPrimRecSSSS_hpp

#include "SimdArray.hpp"

namespace erirec { // erirec namespace

/// Computes [SS|1/|r-r'||SS]  integrals for set of data buffers.
/// - Parameter prim_buffer_0_ssss: the primitive integrals buffer.
/// - Parameter fovl_abcd: the vector of combined overlap factors.
/// - Parameter bf_values: the vector of Boys function values.
auto
comp_prim_electron_repulsion_ssss(CSimdArray<double>& prim_buffer_0_ssss,
                                  const double* fovl_abcd,
                                  const double* bf_values) -> void;
} // erirec namespace

#endif /* ElectronRepulsionPrimRecSSSS_hpp */
