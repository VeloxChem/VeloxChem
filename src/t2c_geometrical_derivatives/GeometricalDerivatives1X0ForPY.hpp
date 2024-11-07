#ifndef GeometricalDerivatives1X0ForPY_hpp
#define GeometricalDerivatives1X0ForPY_hpp

#include "SimdArray.hpp"

namespace t2cgeom {  // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)P|R|S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_100_ps The index of integral in primitive integrals buffer.
/// @param idx_op_ss The index of integral in primitive integrals buffer.
/// @param idx_op_ds The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param ket_comps The number of ket components.
/// @param a_exp The exponent on center A.
auto comp_prim_op_geom_10_px(CSimdArray<double>& pbuffer,
                             const size_t        idx_op_geom_100_ps,
                             const size_t        idx_op_ss,
                             const size_t        idx_op_ds,
                             const size_t        op_comps,
                             const size_t        ket_comps,
                             const double        a_exp) -> void;

}  // namespace t2cgeom

#endif /* GeometricalDerivatives1X0ForPY_hpp */
