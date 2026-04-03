#ifndef GeometricalDerivatives110ForGF_hpp
#define GeometricalDerivatives110ForGF_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)G|d^(1)R/dX^(1)|F]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_110_gf The index of integral in primitive integrals buffer.
/// @param idx_op_df The index of integral in primitive integrals buffer.
/// @param idx_op_fd The index of integral in primitive integrals buffer.
/// @param idx_op_fg The index of integral in primitive integrals buffer.
/// @param idx_op_gf The index of integral in primitive integrals buffer.
/// @param idx_op_hd The index of integral in primitive integrals buffer.
/// @param idx_op_hg The index of integral in primitive integrals buffer.
/// @param idx_op_if The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_110_gf(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_gf,
                         const int idx_op_df,
                         const int idx_op_fd,
                         const int idx_op_fg,
                         const int idx_op_gf,
                         const int idx_op_hd,
                         const int idx_op_hg,
                         const int idx_op_if,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives110ForGF_hpp */
