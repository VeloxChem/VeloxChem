#ifndef PrimitiveOverlapGeom201FF_XYY_ZZZ
#define PrimitiveOverlapGeom201FF_XYY_ZZZ

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates block of primitive <d^(2)/dA^(2)F_XYY|1|d^(1)/dB^(1)F_ZZZ> integrals.

 @param buffer_xx_x the partial integrals buffer.
 @param buffer_xx_y the partial integrals buffer.
 @param buffer_xx_z the partial integrals buffer.
 @param buffer_xy_x the partial integrals buffer.
 @param buffer_xy_y the partial integrals buffer.
 @param buffer_xy_z the partial integrals buffer.
 @param buffer_xz_x the partial integrals buffer.
 @param buffer_xz_y the partial integrals buffer.
 @param buffer_xz_z the partial integrals buffer.
 @param buffer_yy_x the partial integrals buffer.
 @param buffer_yy_y the partial integrals buffer.
 @param buffer_yy_z the partial integrals buffer.
 @param buffer_yz_x the partial integrals buffer.
 @param buffer_yz_y the partial integrals buffer.
 @param buffer_yz_z the partial integrals buffer.
 @param buffer_zz_x the partial integrals buffer.
 @param buffer_zz_y the partial integrals buffer.
 @param buffer_zz_z the partial integrals buffer.
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
auto compPrimitiveOverlapGeom201FF_XYY_ZZZ(TDoubleArray&       buffer_xx_x,
                                           TDoubleArray&       buffer_xx_y,
                                           TDoubleArray&       buffer_xx_z,
                                           TDoubleArray&       buffer_xy_x,
                                           TDoubleArray&       buffer_xy_y,
                                           TDoubleArray&       buffer_xy_z,
                                           TDoubleArray&       buffer_xz_x,
                                           TDoubleArray&       buffer_xz_y,
                                           TDoubleArray&       buffer_xz_z,
                                           TDoubleArray&       buffer_yy_x,
                                           TDoubleArray&       buffer_yy_y,
                                           TDoubleArray&       buffer_yy_z,
                                           TDoubleArray&       buffer_yz_x,
                                           TDoubleArray&       buffer_yz_y,
                                           TDoubleArray&       buffer_yz_z,
                                           TDoubleArray&       buffer_zz_x,
                                           TDoubleArray&       buffer_zz_y,
                                           TDoubleArray&       buffer_zz_z,
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

#endif /* PrimitiveOverlapGeom201FF_XYY_ZZZ */
