#ifndef PrimitiveOverlapGeom300PF_Y_XXZ
#define PrimitiveOverlapGeom300PF_Y_XXZ

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates block of primitive <d^(3)/dA^(3)P_Y|1|F_XXZ> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveOverlapGeom300PF_Y_XXZ(TDoubleArray&       buffer_xxx,
                                         TDoubleArray&       buffer_xxy,
                                         TDoubleArray&       buffer_xxz,
                                         TDoubleArray&       buffer_xyy,
                                         TDoubleArray&       buffer_xyz,
                                         TDoubleArray&       buffer_xzz,
                                         TDoubleArray&       buffer_yyy,
                                         TDoubleArray&       buffer_yyz,
                                         TDoubleArray&       buffer_yzz,
                                         TDoubleArray&       buffer_zzz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

}  // namespace ovlrec

#endif /* PrimitiveOverlapGeom300PF_Y_XXZ */
