#ifndef GeometricalDerivatives0X1ForDG_hpp
#define GeometricalDerivatives0X1ForDG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [D|R|d^(1)/dB^(1)G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_001_dg The index of integral in primitive integrals buffer.
/// @param idx_op_df The index of integral in primitive integrals buffer.
/// @param idx_op_dh The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
auto
comp_geom_deriv_0x1_dg(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_dg,
                       const int idx_op_df,
                       const int idx_op_dh,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives0X1ForDG_hpp */
