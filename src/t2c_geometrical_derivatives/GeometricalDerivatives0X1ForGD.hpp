#ifndef GeometricalDerivatives0X1ForGD_hpp
#define GeometricalDerivatives0X1ForGD_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [G|R|d^(1)/dB^(1)D]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_001_gd The index of integral in primitive integrals buffer.
/// @param idx_op_gp The index of integral in primitive integrals buffer.
/// @param idx_op_gf The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
auto
comp_geom_deriv_0x1_gd(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_gd,
                       const int idx_op_gp,
                       const int idx_op_gf,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives0X1ForGD_hpp */
