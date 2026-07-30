#ifndef GeometricalDerivatives010ForGS_hpp
#define GeometricalDerivatives010ForGS_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [G|d^(1)R/dX^(1)|S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_gs The index of integral in primitive integrals buffer.
/// @param idx_op_fs The index of integral in primitive integrals buffer.
/// @param idx_op_gp The index of integral in primitive integrals buffer.
/// @param idx_op_hs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_gs(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_gs,
                         const int idx_op_fs,
                         const int idx_op_gp,
                         const int idx_op_hs,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForGS_hpp */
