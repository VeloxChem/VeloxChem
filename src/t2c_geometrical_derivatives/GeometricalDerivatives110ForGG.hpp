#ifndef GeometricalDerivatives110ForGG_hpp
#define GeometricalDerivatives110ForGG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)G|d^(1)R/dX^(1)|G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_110_gg The index of integral in primitive integrals buffer.
/// @param idx_op_dg The index of integral in primitive integrals buffer.
/// @param idx_op_ff The index of integral in primitive integrals buffer.
/// @param idx_op_fh The index of integral in primitive integrals buffer.
/// @param idx_op_gg The index of integral in primitive integrals buffer.
/// @param idx_op_hf The index of integral in primitive integrals buffer.
/// @param idx_op_hh The index of integral in primitive integrals buffer.
/// @param idx_op_ig The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_110_gg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_gg,
                         const int idx_op_dg,
                         const int idx_op_ff,
                         const int idx_op_fh,
                         const int idx_op_gg,
                         const int idx_op_hf,
                         const int idx_op_hh,
                         const int idx_op_ig,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives110ForGG_hpp */
