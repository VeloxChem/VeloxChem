#ifndef GeometricalDerivatives110ForSS_hpp
#define GeometricalDerivatives110ForSS_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)S|d^(1)R/dX^(1)|S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_110_ss The index of integral in primitive integrals buffer.
/// @param idx_op_ss The index of integral in primitive integrals buffer.
/// @param idx_op_pp The index of integral in primitive integrals buffer.
/// @param idx_op_ds The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_110_ss(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_ss,
                         const int idx_op_ss,
                         const int idx_op_pp,
                         const int idx_op_ds,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives110ForSS_hpp */
