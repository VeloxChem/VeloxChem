#ifndef PrimitiveOverlapGeom202PP_Z_Z
#define PrimitiveOverlapGeom202PP_Z_Z

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates block of primitive <d^(2)/dA^(2)P_Z|1|d^(2)/dB^(2)P_Z> integrals.

 @param buffer_xx_xx the partial integrals buffer.
 @param buffer_xx_xy the partial integrals buffer.
 @param buffer_xx_xz the partial integrals buffer.
 @param buffer_xx_yy the partial integrals buffer.
 @param buffer_xx_yz the partial integrals buffer.
 @param buffer_xx_zz the partial integrals buffer.
 @param buffer_xy_xx the partial integrals buffer.
 @param buffer_xy_xy the partial integrals buffer.
 @param buffer_xy_xz the partial integrals buffer.
 @param buffer_xy_yy the partial integrals buffer.
 @param buffer_xy_yz the partial integrals buffer.
 @param buffer_xy_zz the partial integrals buffer.
 @param buffer_xz_xx the partial integrals buffer.
 @param buffer_xz_xy the partial integrals buffer.
 @param buffer_xz_xz the partial integrals buffer.
 @param buffer_xz_yy the partial integrals buffer.
 @param buffer_xz_yz the partial integrals buffer.
 @param buffer_xz_zz the partial integrals buffer.
 @param buffer_yy_xx the partial integrals buffer.
 @param buffer_yy_xy the partial integrals buffer.
 @param buffer_yy_xz the partial integrals buffer.
 @param buffer_yy_yy the partial integrals buffer.
 @param buffer_yy_yz the partial integrals buffer.
 @param buffer_yy_zz the partial integrals buffer.
 @param buffer_yz_xx the partial integrals buffer.
 @param buffer_yz_xy the partial integrals buffer.
 @param buffer_yz_xz the partial integrals buffer.
 @param buffer_yz_yy the partial integrals buffer.
 @param buffer_yz_yz the partial integrals buffer.
 @param buffer_yz_zz the partial integrals buffer.
 @param buffer_zz_xx the partial integrals buffer.
 @param buffer_zz_xy the partial integrals buffer.
 @param buffer_zz_xz the partial integrals buffer.
 @param buffer_zz_yy the partial integrals buffer.
 @param buffer_zz_yz the partial integrals buffer.
 @param buffer_zz_zz the partial integrals buffer.
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
auto compPrimitiveOverlapGeom202PP_Z_Z(TDoubleArray&       buffer_xx_xx,
                                       TDoubleArray&       buffer_xx_xy,
                                       TDoubleArray&       buffer_xx_xz,
                                       TDoubleArray&       buffer_xx_yy,
                                       TDoubleArray&       buffer_xx_yz,
                                       TDoubleArray&       buffer_xx_zz,
                                       TDoubleArray&       buffer_xy_xx,
                                       TDoubleArray&       buffer_xy_xy,
                                       TDoubleArray&       buffer_xy_xz,
                                       TDoubleArray&       buffer_xy_yy,
                                       TDoubleArray&       buffer_xy_yz,
                                       TDoubleArray&       buffer_xy_zz,
                                       TDoubleArray&       buffer_xz_xx,
                                       TDoubleArray&       buffer_xz_xy,
                                       TDoubleArray&       buffer_xz_xz,
                                       TDoubleArray&       buffer_xz_yy,
                                       TDoubleArray&       buffer_xz_yz,
                                       TDoubleArray&       buffer_xz_zz,
                                       TDoubleArray&       buffer_yy_xx,
                                       TDoubleArray&       buffer_yy_xy,
                                       TDoubleArray&       buffer_yy_xz,
                                       TDoubleArray&       buffer_yy_yy,
                                       TDoubleArray&       buffer_yy_yz,
                                       TDoubleArray&       buffer_yy_zz,
                                       TDoubleArray&       buffer_yz_xx,
                                       TDoubleArray&       buffer_yz_xy,
                                       TDoubleArray&       buffer_yz_xz,
                                       TDoubleArray&       buffer_yz_yy,
                                       TDoubleArray&       buffer_yz_yz,
                                       TDoubleArray&       buffer_yz_zz,
                                       TDoubleArray&       buffer_zz_xx,
                                       TDoubleArray&       buffer_zz_xy,
                                       TDoubleArray&       buffer_zz_xz,
                                       TDoubleArray&       buffer_zz_yy,
                                       TDoubleArray&       buffer_zz_yz,
                                       TDoubleArray&       buffer_zz_zz,
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

#endif /* PrimitiveOverlapGeom202PP_Z_Z */
