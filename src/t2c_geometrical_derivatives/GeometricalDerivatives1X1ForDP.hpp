#ifndef GeometricalDerivatives1X1ForDP_hpp
#define GeometricalDerivatives1X1ForDP_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)D|R|d^(1)/dB^(1)P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_101_dpThe index of integral in primitive integrals buffer.
/// @param idx_op_geom_001_ppThe index of integral in primitive integrals buffer.
/// @param idx_op_geom_001_fpThe index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
/// @param a_exp The exponent on center A.
auto
comp_prim_op_geom_11_dp(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_dp,
                        const size_t idx_op_ps,
                        const size_t idx_op_pd,
                        const size_t idx_op_fs,
                        const size_t idx_op_fd,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives1X1ForDP_hpp */
