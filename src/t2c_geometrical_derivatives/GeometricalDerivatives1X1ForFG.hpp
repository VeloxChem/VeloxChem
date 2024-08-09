#ifndef GeometricalDerivatives1X1ForFG_hpp
#define GeometricalDerivatives1X1ForFG_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)F|R|d^(1)/dB^(1)G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_101_fgThe index of integral in primitive integrals buffer.
/// @param idx_op_geom_001_dgThe index of integral in primitive integrals buffer.
/// @param idx_op_geom_001_ggThe index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
/// @param a_exp The exponent on center A.
auto
comp_prim_op_geom_11_fg(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_fg,
                        const size_t idx_op_df,
                        const size_t idx_op_dh,
                        const size_t idx_op_gf,
                        const size_t idx_op_gh,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives1X1ForFG_hpp */
