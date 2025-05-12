#ifndef GeometricalDerivatives0X1ForSS_hpp
#define GeometricalDerivatives0X1ForSS_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [S|R|d^(1)/dB^(1)S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_001_ss The index of integral in primitive integrals buffer.
/// @param idx_op_sp The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
auto
comp_geom_deriv_0x1_ss(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_ss,
                       const int idx_op_sp,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives0X1ForSS_hpp */
