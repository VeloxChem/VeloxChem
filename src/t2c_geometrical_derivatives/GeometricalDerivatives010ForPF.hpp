#ifndef GeometricalDerivatives010ForPF_hpp
#define GeometricalDerivatives010ForPF_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [P|d^(1)R/dX^(1)|F]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_010_pf The index of integral in primitive integrals buffer.
/// @param idx_op_sf The index of integral in primitive integrals buffer.
/// @param idx_op_pd The index of integral in primitive integrals buffer.
/// @param idx_op_pg The index of integral in primitive integrals buffer.
/// @param idx_op_df The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_op_geom_010_pf(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_pf,
                         const int idx_op_sf,
                         const int idx_op_pd,
                         const int idx_op_pg,
                         const int idx_op_df,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives010ForPF_hpp */
