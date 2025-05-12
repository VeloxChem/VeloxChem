#ifndef GeometricalDerivatives0X1ForGP_hpp
#define GeometricalDerivatives0X1ForGP_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [G|R|d^(1)/dB^(1)P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_001_gp The index of integral in primitive integrals buffer.
/// @param idx_op_gs The index of integral in primitive integrals buffer.
/// @param idx_op_gd The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
auto
comp_geom_deriv_0x1_gp(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_gp,
                       const int idx_op_gs,
                       const int idx_op_gd,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives0X1ForGP_hpp */
