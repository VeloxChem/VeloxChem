#ifndef GeometricalDerivatives020ForDF_hpp
#define GeometricalDerivatives020ForDF_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [D|d^(2)R/dX^(2)|F]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_020_df The index of integral in primitive integrals buffer.
/// @param idx_op_sf The index of integral in primitive integrals buffer.
/// @param idx_op_pd The index of integral in primitive integrals buffer.
/// @param idx_op_pg The index of integral in primitive integrals buffer.
/// @param idx_op_dp The index of integral in primitive integrals buffer.
/// @param idx_op_df The index of integral in primitive integrals buffer.
/// @param idx_op_dh The index of integral in primitive integrals buffer.
/// @param idx_op_fd The index of integral in primitive integrals buffer.
/// @param idx_op_fg The index of integral in primitive integrals buffer.
/// @param idx_op_gf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_020_df(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_df,
                         const int idx_op_sf,
                         const int idx_op_pd,
                         const int idx_op_pg,
                         const int idx_op_dp,
                         const int idx_op_df,
                         const int idx_op_dh,
                         const int idx_op_fd,
                         const int idx_op_fg,
                         const int idx_op_gf,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives020ForDF_hpp */
