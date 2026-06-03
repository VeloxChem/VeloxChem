#ifndef GeometricalDerivatives110ForFP_hpp
#define GeometricalDerivatives110ForFP_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)F|d^(1)R/dX^(1)|P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_110_fp The index of integral in primitive integrals buffer.
/// @param idx_op_pp The index of integral in primitive integrals buffer.
/// @param idx_op_ds The index of integral in primitive integrals buffer.
/// @param idx_op_dd The index of integral in primitive integrals buffer.
/// @param idx_op_fp The index of integral in primitive integrals buffer.
/// @param idx_op_gs The index of integral in primitive integrals buffer.
/// @param idx_op_gd The index of integral in primitive integrals buffer.
/// @param idx_op_hp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_110_fp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_fp,
                         const int idx_op_pp,
                         const int idx_op_ds,
                         const int idx_op_dd,
                         const int idx_op_fp,
                         const int idx_op_gs,
                         const int idx_op_gd,
                         const int idx_op_hp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives110ForFP_hpp */
