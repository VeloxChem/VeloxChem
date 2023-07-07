#ifndef PrimitiveOverlapGeom400DF_XY_ZZZ
#define PrimitiveOverlapGeom400DF_XY_ZZZ

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates block of primitive <d^(4)/dA^(4)D_XY|1|F_ZZZ> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
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
auto compPrimitiveOverlapGeom400DF_XY_ZZZ(TDoubleArray&       buffer_xxxx,
                                          TDoubleArray&       buffer_xxxy,
                                          TDoubleArray&       buffer_xxxz,
                                          TDoubleArray&       buffer_xxyy,
                                          TDoubleArray&       buffer_xxyz,
                                          TDoubleArray&       buffer_xxzz,
                                          TDoubleArray&       buffer_xyyy,
                                          TDoubleArray&       buffer_xyyz,
                                          TDoubleArray&       buffer_xyzz,
                                          TDoubleArray&       buffer_xzzz,
                                          TDoubleArray&       buffer_yyyy,
                                          TDoubleArray&       buffer_yyyz,
                                          TDoubleArray&       buffer_yyzz,
                                          TDoubleArray&       buffer_yzzz,
                                          TDoubleArray&       buffer_zzzz,
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

#endif /* PrimitiveOverlapGeom400DF_XY_ZZZ */
