#ifndef GeometricalDerivatives010ForGG_hpp
#define GeometricalDerivatives010ForGG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [G|d^(1)R/dX^(1)|G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_gg The index of integral in primitive integrals buffer.
/// @param idx_op_fg The index of integral in primitive integrals buffer.
/// @param idx_op_gf The index of integral in primitive integrals buffer.
/// @param idx_op_gh The index of integral in primitive integrals buffer.
/// @param idx_op_hg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_gg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_gg,
                         const int idx_op_fg,
                         const int idx_op_gf,
                         const int idx_op_gh,
                         const int idx_op_hg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForGG_hpp */
