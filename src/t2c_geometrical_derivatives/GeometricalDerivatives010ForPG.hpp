#ifndef GeometricalDerivatives010ForPG_hpp
#define GeometricalDerivatives010ForPG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [P|d^(1)R/dX^(1)|G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_pg The index of integral in primitive integrals buffer.
/// @param idx_op_sg The index of integral in primitive integrals buffer.
/// @param idx_op_pf The index of integral in primitive integrals buffer.
/// @param idx_op_ph The index of integral in primitive integrals buffer.
/// @param idx_op_dg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_pg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_pg,
                         const int idx_op_sg,
                         const int idx_op_pf,
                         const int idx_op_ph,
                         const int idx_op_dg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForPG_hpp */
