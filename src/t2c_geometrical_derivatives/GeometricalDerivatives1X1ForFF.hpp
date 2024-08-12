#ifndef GeometricalDerivatives1X1ForFF_hpp
#define GeometricalDerivatives1X1ForFF_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)F|R|d^(1)/dB^(1)F]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_101_ff The index of integral in primitive integrals buffer.
/// @param idx_op_dd The index of integral in primitive integrals buffer.
/// @param idx_op_dg The index of integral in primitive integrals buffer.
/// @param idx_op_gd The index of integral in primitive integrals buffer.
/// @param idx_op_gg The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
/// @param a_exp The exponent on center A.
auto
comp_prim_op_geom_11_ff(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_ff,
                        const size_t idx_op_dd,
                        const size_t idx_op_dg,
                        const size_t idx_op_gd,
                        const size_t idx_op_gg,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives1X1ForFF_hpp */
