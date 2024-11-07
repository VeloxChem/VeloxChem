#ifndef GeometricalDerivatives1X1ForGP_hpp
#define GeometricalDerivatives1X1ForGP_hpp

#include "SimdArray.hpp"

namespace t2cgeom {  // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)G|R|d^(1)/dB^(1)P]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_101_gp The index of integral in primitive integrals buffer.
/// @param idx_op_fs The index of integral in primitive integrals buffer.
/// @param idx_op_fd The index of integral in primitive integrals buffer.
/// @param idx_op_hs The index of integral in primitive integrals buffer.
/// @param idx_op_hd The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
/// @param a_exp The exponent on center A.
auto comp_prim_op_geom_11_gp(CSimdArray<double>&       pbuffer,
                             const size_t              idx_op_geom_101_gp,
                             const size_t              idx_op_fs,
                             const size_t              idx_op_fd,
                             const size_t              idx_op_hs,
                             const size_t              idx_op_hd,
                             const size_t              op_comps,
                             const CSimdArray<double>& factors,
                             const double              a_exp) -> void;

}  // namespace t2cgeom

#endif /* GeometricalDerivatives1X1ForGP_hpp */
