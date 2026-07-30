#ifndef GeometricalDerivatives010ForGD_hpp
#define GeometricalDerivatives010ForGD_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [G|d^(1)R/dX^(1)|D]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_gd The index of integral in primitive integrals buffer.
/// @param idx_op_fd The index of integral in primitive integrals buffer.
/// @param idx_op_gp The index of integral in primitive integrals buffer.
/// @param idx_op_gf The index of integral in primitive integrals buffer.
/// @param idx_op_hd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_gd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_gd,
                         const int idx_op_fd,
                         const int idx_op_gp,
                         const int idx_op_gf,
                         const int idx_op_hd,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForGD_hpp */
