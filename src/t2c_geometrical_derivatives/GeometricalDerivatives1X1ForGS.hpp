#ifndef GeometricalDerivatives1X1ForGS_hpp
#define GeometricalDerivatives1X1ForGS_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)G|R|d^(1)/dB^(1)S]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_101_gsThe index of integral in primitive integrals buffer.
/// @param idx_op_geom_001_fsThe index of integral in primitive integrals buffer.
/// @param idx_op_geom_001_hsThe index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
/// @param a_exp The exponent on center A.
auto
comp_prim_op_geom_11_gs(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_gs,
                        const size_t idx_op_fp,
                        const size_t idx_op_hp,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives1X1ForGS_hpp */
