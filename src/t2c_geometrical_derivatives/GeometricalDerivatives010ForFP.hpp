#ifndef GeometricalDerivatives010ForFP_hpp
#define GeometricalDerivatives010ForFP_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [F|d^(1)R/dX^(1)|P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_fp The index of integral in primitive integrals buffer.
/// @param idx_op_dp The index of integral in primitive integrals buffer.
/// @param idx_op_fs The index of integral in primitive integrals buffer.
/// @param idx_op_fd The index of integral in primitive integrals buffer.
/// @param idx_op_gp The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_fp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_fp,
                         const int idx_op_dp,
                         const int idx_op_fs,
                         const int idx_op_fd,
                         const int idx_op_gp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForFP_hpp */
