#ifndef PrimitiveOverlapGeom301FD_XZZ_XZ
#define PrimitiveOverlapGeom301FD_XZZ_XZ

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates block of primitive <d^(3)/dA^(3)F_XZZ|1|d^(1)/dB^(1)D_XZ> integrals.

 @param buffer_xxx_x the partial integrals buffer.
 @param buffer_xxx_y the partial integrals buffer.
 @param buffer_xxx_z the partial integrals buffer.
 @param buffer_xxy_x the partial integrals buffer.
 @param buffer_xxy_y the partial integrals buffer.
 @param buffer_xxy_z the partial integrals buffer.
 @param buffer_xxz_x the partial integrals buffer.
 @param buffer_xxz_y the partial integrals buffer.
 @param buffer_xxz_z the partial integrals buffer.
 @param buffer_xyy_x the partial integrals buffer.
 @param buffer_xyy_y the partial integrals buffer.
 @param buffer_xyy_z the partial integrals buffer.
 @param buffer_xyz_x the partial integrals buffer.
 @param buffer_xyz_y the partial integrals buffer.
 @param buffer_xyz_z the partial integrals buffer.
 @param buffer_xzz_x the partial integrals buffer.
 @param buffer_xzz_y the partial integrals buffer.
 @param buffer_xzz_z the partial integrals buffer.
 @param buffer_yyy_x the partial integrals buffer.
 @param buffer_yyy_y the partial integrals buffer.
 @param buffer_yyy_z the partial integrals buffer.
 @param buffer_yyz_x the partial integrals buffer.
 @param buffer_yyz_y the partial integrals buffer.
 @param buffer_yyz_z the partial integrals buffer.
 @param buffer_yzz_x the partial integrals buffer.
 @param buffer_yzz_y the partial integrals buffer.
 @param buffer_yzz_z the partial integrals buffer.
 @param buffer_zzz_x the partial integrals buffer.
 @param buffer_zzz_y the partial integrals buffer.
 @param buffer_zzz_z the partial integrals buffer.
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
auto compPrimitiveOverlapGeom301FD_XZZ_XZ(TDoubleArray&       buffer_xxx_x,
                                          TDoubleArray&       buffer_xxx_y,
                                          TDoubleArray&       buffer_xxx_z,
                                          TDoubleArray&       buffer_xxy_x,
                                          TDoubleArray&       buffer_xxy_y,
                                          TDoubleArray&       buffer_xxy_z,
                                          TDoubleArray&       buffer_xxz_x,
                                          TDoubleArray&       buffer_xxz_y,
                                          TDoubleArray&       buffer_xxz_z,
                                          TDoubleArray&       buffer_xyy_x,
                                          TDoubleArray&       buffer_xyy_y,
                                          TDoubleArray&       buffer_xyy_z,
                                          TDoubleArray&       buffer_xyz_x,
                                          TDoubleArray&       buffer_xyz_y,
                                          TDoubleArray&       buffer_xyz_z,
                                          TDoubleArray&       buffer_xzz_x,
                                          TDoubleArray&       buffer_xzz_y,
                                          TDoubleArray&       buffer_xzz_z,
                                          TDoubleArray&       buffer_yyy_x,
                                          TDoubleArray&       buffer_yyy_y,
                                          TDoubleArray&       buffer_yyy_z,
                                          TDoubleArray&       buffer_yyz_x,
                                          TDoubleArray&       buffer_yyz_y,
                                          TDoubleArray&       buffer_yyz_z,
                                          TDoubleArray&       buffer_yzz_x,
                                          TDoubleArray&       buffer_yzz_y,
                                          TDoubleArray&       buffer_yzz_z,
                                          TDoubleArray&       buffer_zzz_x,
                                          TDoubleArray&       buffer_zzz_y,
                                          TDoubleArray&       buffer_zzz_z,
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

#endif /* PrimitiveOverlapGeom301FD_XZZ_XZ */
