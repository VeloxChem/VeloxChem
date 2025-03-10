#ifndef GeometricalDerivatives1X0ForIY_hpp
#define GeometricalDerivatives1X0ForIY_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)I|R|S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_100_is The index of integral in primitive integrals buffer.
/// @param idx_op_hs The index of integral in primitive integrals buffer.
/// @param idx_op_ks The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param ket_comps The number of ket components.
/// @param a_exp The exponent on center A.
auto
comp_prim_op_geom_10_ix(CSimdArray<double>& prim_buffer,
                        const int idx_op_geom_100_is,
                        const int idx_op_hs,
                        const int idx_op_ks,
                        const int op_comps,
                        const int ket_comps,
                        const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives1X0ForIY_hpp */
