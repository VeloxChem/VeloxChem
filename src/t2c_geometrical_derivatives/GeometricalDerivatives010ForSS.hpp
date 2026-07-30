#ifndef GeometricalDerivatives010ForSS_hpp
#define GeometricalDerivatives010ForSS_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [S|d^(1)R/dX^(1)|S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_ss The index of integral in primitive integrals buffer.
/// @param idx_op_sp The index of integral in primitive integrals buffer.
/// @param idx_op_ps The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_ss(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_ss,
                         const int idx_op_sp,
                         const int idx_op_ps,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForSS_hpp */
