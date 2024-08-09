#ifndef GeometricalDerivatives2X0ForFY_hpp
#define GeometricalDerivatives2X0ForFY_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(2)/dA^(2)F|R|S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_200_fsThe index of integral in primitive integrals buffer.
/// @param idx_op_geom_100_dsThe index of integral in primitive integrals buffer.
/// @param idx_op_geom_100_gsThe index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param ket_comps The number of ket components.
/// @param a_exp The exponent on center A.
auto
comp_prim_op_geom_20_fx(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_200_fs,
                        const size_t idx_op_ps,
                        const size_t idx_op_fs,
                        const size_t idx_op_hs,
                        const size_t op_comps,
                        const size_t ket_comps,
                        const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives2X0ForFY_hpp */
