#ifndef PrimitiveQuadrupoleFD_XXZ_XZ
#define PrimitiveQuadrupoleFD_XXZ_XZ

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace mpol {  // mpol namespace

/**
 Evaluates block of primitive <F_XXZ|r^2|D_XZ> integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param point the coordinates of external point.
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
auto compPrimitiveQuadrupoleFD_XXZ_XZ(TDoubleArray&       buffer_xx,
                                      TDoubleArray&       buffer_xy,
                                      TDoubleArray&       buffer_xz,
                                      TDoubleArray&       buffer_yy,
                                      TDoubleArray&       buffer_yz,
                                      TDoubleArray&       buffer_zz,
                                      const TPoint3D&     point,
                                      const double        bra_exp,
                                      const double        bra_norm,
                                      const TPoint3D&     bra_coord,
                                      const TDoubleArray& ket_exps,
                                      const TDoubleArray& ket_norms,
                                      const TDoubleArray& ket_coords_x,
                                      const TDoubleArray& ket_coords_y,
                                      const TDoubleArray& ket_coords_z,
                                      const int64_t       ket_dim) -> void;

}  // namespace mpol

#endif /* PrimitiveQuadrupoleFD_XXZ_XZ */
