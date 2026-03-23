#ifndef GeometricalDerivatives010ForPP_hpp
#define GeometricalDerivatives010ForPP_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [P|d^(1)R/dX^(1)|P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_pp The index of integral in primitive integrals buffer.
/// @param idx_op_sp The index of integral in primitive integrals buffer.
/// @param idx_op_ps The index of integral in primitive integrals buffer.
/// @param idx_op_pd The index of integral in primitive integrals buffer.
/// @param idx_op_dp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_pp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_pp,
                         const int idx_op_sp,
                         const int idx_op_ps,
                         const int idx_op_pd,
                         const int idx_op_dp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForPP_hpp */
