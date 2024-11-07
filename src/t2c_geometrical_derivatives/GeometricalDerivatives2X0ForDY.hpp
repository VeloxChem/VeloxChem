#ifndef GeometricalDerivatives2X0ForDY_hpp
#define GeometricalDerivatives2X0ForDY_hpp

#include "SimdArray.hpp"

namespace t2cgeom {  // t2cgeom namespace

/// @brief Computes [d^(2)/dA^(2)D|R|S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_200_ds The index of integral in primitive integrals buffer.
/// @param idx_op_ss The index of integral in primitive integrals buffer.
/// @param idx_op_ds The index of integral in primitive integrals buffer.
/// @param idx_op_gs The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param ket_comps The number of ket components.
/// @param a_exp The exponent on center A.
auto comp_prim_op_geom_20_dx(CSimdArray<double>& pbuffer,
                             const size_t        idx_op_geom_200_ds,
                             const size_t        idx_op_ss,
                             const size_t        idx_op_ds,
                             const size_t        idx_op_gs,
                             const size_t        op_comps,
                             const size_t        ket_comps,
                             const double        a_exp) -> void;

}  // namespace t2cgeom

#endif /* GeometricalDerivatives2X0ForDY_hpp */
