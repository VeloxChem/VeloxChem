#ifndef GeometricalDerivatives110ForFD_hpp
#define GeometricalDerivatives110ForFD_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)F|d^(1)R/dX^(1)|D]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_110_fd The index of integral in primitive integrals buffer.
/// @param idx_op_pd The index of integral in primitive integrals buffer.
/// @param idx_op_dp The index of integral in primitive integrals buffer.
/// @param idx_op_df The index of integral in primitive integrals buffer.
/// @param idx_op_fd The index of integral in primitive integrals buffer.
/// @param idx_op_gp The index of integral in primitive integrals buffer.
/// @param idx_op_gf The index of integral in primitive integrals buffer.
/// @param idx_op_hd The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_110_fd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_fd,
                         const int idx_op_pd,
                         const int idx_op_dp,
                         const int idx_op_df,
                         const int idx_op_fd,
                         const int idx_op_gp,
                         const int idx_op_gf,
                         const int idx_op_hd,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives110ForFD_hpp */
