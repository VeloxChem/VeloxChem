#ifndef GeometricalDerivatives1X1ForGG_hpp
#define GeometricalDerivatives1X1ForGG_hpp

#include "SimdArray.hpp"

namespace t2cgeom {  // t2cgeom namespace

/// @brief Computes [d^(1)/dA^(1)G|R|d^(1)/dB^(1)G]  integrals for arbitrary operator R.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_op_geom_101_gg The index of integral in primitive integrals buffer.
/// @param idx_op_ff The index of integral in primitive integrals buffer.
/// @param idx_op_fh The index of integral in primitive integrals buffer.
/// @param idx_op_hf The index of integral in primitive integrals buffer.
/// @param idx_op_hh The index of integral in primitive integrals buffer.
/// @param op_comps The number of operator components.
/// @param factors The primitive factors buffer.
/// @param a_exp The exponent on center A.
auto comp_prim_op_geom_11_gg(CSimdArray<double>&       pbuffer,
                             const size_t              idx_op_geom_101_gg,
                             const size_t              idx_op_ff,
                             const size_t              idx_op_fh,
                             const size_t              idx_op_hf,
                             const size_t              idx_op_hh,
                             const size_t              op_comps,
                             const CSimdArray<double>& factors,
                             const double              a_exp) -> void;

}  // namespace t2cgeom

#endif /* GeometricalDerivatives1X1ForGG_hpp */
