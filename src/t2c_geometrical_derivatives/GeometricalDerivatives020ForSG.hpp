#ifndef GeometricalDerivatives020ForSG_hpp
#define GeometricalDerivatives020ForSG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [S|d^(2)R/dX^(2)|G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_020_sg The index of integral in primitive integrals buffer.
/// @param idx_op_sd The index of integral in primitive integrals buffer.
/// @param idx_op_sg The index of integral in primitive integrals buffer.
/// @param idx_op_si The index of integral in primitive integrals buffer.
/// @param idx_op_pf The index of integral in primitive integrals buffer.
/// @param idx_op_ph The index of integral in primitive integrals buffer.
/// @param idx_op_dg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_020_sg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_sg,
                         const int idx_op_sd,
                         const int idx_op_sg,
                         const int idx_op_si,
                         const int idx_op_pf,
                         const int idx_op_ph,
                         const int idx_op_dg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives020ForSG_hpp */
