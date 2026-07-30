#ifndef GeometricalDerivatives010ForFS_hpp
#define GeometricalDerivatives010ForFS_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [F|d^(1)R/dX^(1)|S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_fs The index of integral in primitive integrals buffer.
/// @param idx_op_ds The index of integral in primitive integrals buffer.
/// @param idx_op_fp The index of integral in primitive integrals buffer.
/// @param idx_op_gs The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_fs(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_fs,
                         const int idx_op_ds,
                         const int idx_op_fp,
                         const int idx_op_gs,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForFS_hpp */
