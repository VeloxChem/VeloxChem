#ifndef GeometricalDerivatives020ForFF_hpp
#define GeometricalDerivatives020ForFF_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [F|d^(2)R/dX^(2)|F]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_020_ff The index of integral in primitive integrals buffer.
/// @param idx_op_pf The index of integral in primitive integrals buffer.
/// @param idx_op_dd The index of integral in primitive integrals buffer.
/// @param idx_op_dg The index of integral in primitive integrals buffer.
/// @param idx_op_fp The index of integral in primitive integrals buffer.
/// @param idx_op_ff The index of integral in primitive integrals buffer.
/// @param idx_op_fh The index of integral in primitive integrals buffer.
/// @param idx_op_gd The index of integral in primitive integrals buffer.
/// @param idx_op_gg The index of integral in primitive integrals buffer.
/// @param idx_op_hf The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_020_ff(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_ff,
                         const int idx_op_pf,
                         const int idx_op_dd,
                         const int idx_op_dg,
                         const int idx_op_fp,
                         const int idx_op_ff,
                         const int idx_op_fh,
                         const int idx_op_gd,
                         const int idx_op_gg,
                         const int idx_op_hf,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives020ForFF_hpp */
