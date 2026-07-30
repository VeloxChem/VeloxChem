#ifndef GeometricalDerivatives010ForDP_hpp
#define GeometricalDerivatives010ForDP_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [D|d^(1)R/dX^(1)|P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_dp The index of integral in primitive integrals buffer.
/// @param idx_op_pp The index of integral in primitive integrals buffer.
/// @param idx_op_ds The index of integral in primitive integrals buffer.
/// @param idx_op_dd The index of integral in primitive integrals buffer.
/// @param idx_op_fp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_dp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_dp,
                         const int idx_op_pp,
                         const int idx_op_ds,
                         const int idx_op_dd,
                         const int idx_op_fp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForDP_hpp */
