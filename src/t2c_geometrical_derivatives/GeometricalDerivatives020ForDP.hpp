#ifndef GeometricalDerivatives020ForDP_hpp
#define GeometricalDerivatives020ForDP_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [D|d^(2)R/dX^(2)|P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_020_dp The index of integral in primitive integrals buffer.
/// @param idx_op_sp The index of integral in primitive integrals buffer.
/// @param idx_op_ps The index of integral in primitive integrals buffer.
/// @param idx_op_pd The index of integral in primitive integrals buffer.
/// @param idx_op_dp The index of integral in primitive integrals buffer.
/// @param idx_op_df The index of integral in primitive integrals buffer.
/// @param idx_op_fs The index of integral in primitive integrals buffer.
/// @param idx_op_fd The index of integral in primitive integrals buffer.
/// @param idx_op_gp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_020_dp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_dp,
                         const int idx_op_sp,
                         const int idx_op_ps,
                         const int idx_op_pd,
                         const int idx_op_dp,
                         const int idx_op_df,
                         const int idx_op_fs,
                         const int idx_op_fd,
                         const int idx_op_gp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives020ForDP_hpp */
