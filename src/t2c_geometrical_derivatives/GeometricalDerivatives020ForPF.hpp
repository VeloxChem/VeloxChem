#ifndef GeometricalDerivatives020ForPF_hpp
#define GeometricalDerivatives020ForPF_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [P|d^(2)R/dX^(2)|F]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_020_pf The index of integral in primitive integrals buffer.
/// @param idx_op_sd The index of integral in primitive integrals buffer.
/// @param idx_op_sg The index of integral in primitive integrals buffer.
/// @param idx_op_pp The index of integral in primitive integrals buffer.
/// @param idx_op_pf The index of integral in primitive integrals buffer.
/// @param idx_op_ph The index of integral in primitive integrals buffer.
/// @param idx_op_dd The index of integral in primitive integrals buffer.
/// @param idx_op_dg The index of integral in primitive integrals buffer.
/// @param idx_op_ff The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_020_pf(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_pf,
                         const int idx_op_sd,
                         const int idx_op_sg,
                         const int idx_op_pp,
                         const int idx_op_pf,
                         const int idx_op_ph,
                         const int idx_op_dd,
                         const int idx_op_dg,
                         const int idx_op_ff,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives020ForPF_hpp */
