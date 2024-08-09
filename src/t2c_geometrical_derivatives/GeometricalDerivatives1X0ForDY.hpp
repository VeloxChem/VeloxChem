#ifndef GeometricalDerivatives1X0ForDY_hpp
#define GeometricalDerivatives1X0ForDY_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)D|R|S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_100_ds The index of integral in primitive integrals buffer.
/// @param idx_op_ps The index of integral in primitive integrals buffer.
/// @param idx_op_fs The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param ket_comps The number of ket components.
/// @param a_exp The exponent on center A.
auto
comp_prim_op_geom_10_dx(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_100_ds,
                        const size_t idx_op_ps,
                        const size_t idx_op_fs,
                        const size_t op_comps,
                        const size_t ket_comps,
                        const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives1X0ForDY_hpp */
