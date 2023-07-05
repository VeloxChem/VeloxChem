#ifndef PrimitiveNuclearPotentialGeom020FF_YZZ_XYY
#define PrimitiveNuclearPotentialGeom020FF_YZZ_XYY

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"
#include "TensorTypes.hpp"

namespace npotg020rec {  // npotg020rec namespace

/**
 Evaluates block of primitive <F_YZZ|AG(2)|F_XYY> integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param quadrupole the quadrupole of external point.
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
auto compPrimitiveNuclearPotentialGeom020FF_YZZ_XYY(TDoubleArray&       buffer_xx,
                                                    TDoubleArray&       buffer_xy,
                                                    TDoubleArray&       buffer_xz,
                                                    TDoubleArray&       buffer_yy,
                                                    TDoubleArray&       buffer_yz,
                                                    TDoubleArray&       buffer_zz,
                                                    const T2Tensor&     quadrupole,
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

}  // namespace npotg020rec

#endif /* PrimitiveNuclearPotentialGeom020FF_YZZ_XYY */
