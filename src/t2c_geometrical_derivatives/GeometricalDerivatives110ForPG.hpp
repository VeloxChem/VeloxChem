#ifndef GeometricalDerivatives110ForPG_hpp
#define GeometricalDerivatives110ForPG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)P|d^(1)R/dX^(1)|G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_110_pg The index of integral in primitive integrals buffer.
/// @param idx_op_sf The index of integral in primitive integrals buffer.
/// @param idx_op_sh The index of integral in primitive integrals buffer.
/// @param idx_op_pg The index of integral in primitive integrals buffer.
/// @param idx_op_df The index of integral in primitive integrals buffer.
/// @param idx_op_dh The index of integral in primitive integrals buffer.
/// @param idx_op_fg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_110_pg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_pg,
                         const int idx_op_sf,
                         const int idx_op_sh,
                         const int idx_op_pg,
                         const int idx_op_df,
                         const int idx_op_dh,
                         const int idx_op_fg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives110ForPG_hpp */
