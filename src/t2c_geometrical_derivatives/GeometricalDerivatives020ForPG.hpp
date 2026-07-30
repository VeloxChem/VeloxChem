#ifndef GeometricalDerivatives020ForPG_hpp
#define GeometricalDerivatives020ForPG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [P|d^(2)R/dX^(2)|G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_020_pg The index of integral in primitive integrals buffer.
/// @param idx_op_sf The index of integral in primitive integrals buffer.
/// @param idx_op_sh The index of integral in primitive integrals buffer.
/// @param idx_op_pd The index of integral in primitive integrals buffer.
/// @param idx_op_pg The index of integral in primitive integrals buffer.
/// @param idx_op_pi The index of integral in primitive integrals buffer.
/// @param idx_op_df The index of integral in primitive integrals buffer.
/// @param idx_op_dh The index of integral in primitive integrals buffer.
/// @param idx_op_fg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_020_pg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_pg,
                         const int idx_op_sf,
                         const int idx_op_sh,
                         const int idx_op_pd,
                         const int idx_op_pg,
                         const int idx_op_pi,
                         const int idx_op_df,
                         const int idx_op_dh,
                         const int idx_op_fg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives020ForPG_hpp */
