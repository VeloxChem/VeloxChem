#ifndef GeometricalDerivatives1X1ForDF_hpp
#define GeometricalDerivatives1X1ForDF_hpp

#include "SimdArray.hpp"

namespace t2cgeom { // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)D|R|d^(1)/dB^(1)F]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_101_df The index of integral in primitive integrals buffer.
/// @param idx_op_geom_001_pf The index of integral in primitive integrals buffer.
/// @param idx_op_geom_001_ff The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
/// @param a_exp The exponent on center A.
auto
comp_prim_op_geom_11_df(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_df,
                        const size_t idx_op_pd,
                        const size_t idx_op_pg,
                        const size_t idx_op_fd,
                        const size_t idx_op_fg,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void;

} // t2cgeom namespace

#endif /* GeometricalDerivatives1X1ForDF_hpp */
