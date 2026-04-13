#ifndef GeometricalDerivatives020ForGG_hpp
#define GeometricalDerivatives020ForGG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [G|d^(2)R/dX^(2)|G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_020_gg The index of integral in primitive integrals buffer.
/// @param idx_op_dg The index of integral in primitive integrals buffer.
/// @param idx_op_ff The index of integral in primitive integrals buffer.
/// @param idx_op_fh The index of integral in primitive integrals buffer.
/// @param idx_op_gd The index of integral in primitive integrals buffer.
/// @param idx_op_gg The index of integral in primitive integrals buffer.
/// @param idx_op_gi The index of integral in primitive integrals buffer.
/// @param idx_op_hf The index of integral in primitive integrals buffer.
/// @param idx_op_hh The index of integral in primitive integrals buffer.
/// @param idx_op_ig The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_020_gg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_gg,
                         const int idx_op_dg,
                         const int idx_op_ff,
                         const int idx_op_fh,
                         const int idx_op_gd,
                         const int idx_op_gg,
                         const int idx_op_gi,
                         const int idx_op_hf,
                         const int idx_op_hh,
                         const int idx_op_ig,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives020ForGG_hpp */
