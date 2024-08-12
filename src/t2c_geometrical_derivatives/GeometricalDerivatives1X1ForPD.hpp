#ifndef GeometricalDerivatives1X1ForPD_hpp
#define GeometricalDerivatives1X1ForPD_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)P|R|d^(1)/dB^(1)D]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_101_pd The index of integral in primitive integrals buffer.
/// @param idx_op_sp The index of integral in primitive integrals buffer.
/// @param idx_op_sf The index of integral in primitive integrals buffer.
/// @param idx_op_dp The index of integral in primitive integrals buffer.
/// @param idx_op_df The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
/// @param a_exp The exponent on center A.
auto
comp_prim_op_geom_11_pd(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_pd,
                        const size_t idx_op_sp,
                        const size_t idx_op_sf,
                        const size_t idx_op_dp,
                        const size_t idx_op_df,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives1X1ForPD_hpp */
