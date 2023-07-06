#ifndef PrimitiveOverlapGeom100FS_YZZ_0
#define PrimitiveOverlapGeom100FS_YZZ_0

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates block of primitive <d^(1)/dA^(1)F_YZZ|1|S_0> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
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
auto compPrimitiveOverlapGeom100FS_YZZ_0(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
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

#endif /* PrimitiveOverlapGeom100FS_YZZ_0 */
