#ifndef GeometricalDerivatives010ForGP_hpp
#define GeometricalDerivatives010ForGP_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [G|d^(1)R/dX^(1)|P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_gp The index of integral in primitive integrals buffer.
/// @param idx_op_fp The index of integral in primitive integrals buffer.
/// @param idx_op_gs The index of integral in primitive integrals buffer.
/// @param idx_op_gd The index of integral in primitive integrals buffer.
/// @param idx_op_hp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_gp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_gp,
                         const int idx_op_fp,
                         const int idx_op_gs,
                         const int idx_op_gd,
                         const int idx_op_hp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForGP_hpp */
