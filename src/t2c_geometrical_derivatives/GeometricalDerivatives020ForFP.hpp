#ifndef GeometricalDerivatives020ForFP_hpp
#define GeometricalDerivatives020ForFP_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [F|d^(2)R/dX^(2)|P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_020_fp The index of integral in primitive integrals buffer.
/// @param idx_op_pp The index of integral in primitive integrals buffer.
/// @param idx_op_ds The index of integral in primitive integrals buffer.
/// @param idx_op_dd The index of integral in primitive integrals buffer.
/// @param idx_op_fp The index of integral in primitive integrals buffer.
/// @param idx_op_ff The index of integral in primitive integrals buffer.
/// @param idx_op_gs The index of integral in primitive integrals buffer.
/// @param idx_op_gd The index of integral in primitive integrals buffer.
/// @param idx_op_hp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_020_fp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_fp,
                         const int idx_op_pp,
                         const int idx_op_ds,
                         const int idx_op_dd,
                         const int idx_op_fp,
                         const int idx_op_ff,
                         const int idx_op_gs,
                         const int idx_op_gd,
                         const int idx_op_hp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives020ForFP_hpp */
