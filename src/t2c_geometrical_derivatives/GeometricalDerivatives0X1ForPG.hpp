#ifndef GeometricalDerivatives0X1ForPG_hpp
#define GeometricalDerivatives0X1ForPG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [P|R|d^(1)/dB^(1)G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_001_pg The index of integral in primitive integrals buffer.
/// @param idx_op_pf The index of integral in primitive integrals buffer.
/// @param idx_op_ph The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
auto
comp_geom_deriv_0x1_pg(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_pg,
                       const int idx_op_pf,
                       const int idx_op_ph,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives0X1ForPG_hpp */
