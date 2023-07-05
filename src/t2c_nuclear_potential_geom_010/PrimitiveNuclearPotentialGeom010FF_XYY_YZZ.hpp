#ifndef PrimitiveNuclearPotentialGeom010FF_XYY_YZZ
#define PrimitiveNuclearPotentialGeom010FF_XYY_YZZ

#include <cstdint>

#include "Point.hpp"
#include "SimdTypes.hpp"

namespace npotg010rec {  // npotg010rec namespace

/**
 Evaluates block of primitive <F_XYY|AG(1)|F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
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
auto compPrimitiveNuclearPotentialGeom010FF_XYY_YZZ(TDoubleArray&       buffer_x,
                                                    TDoubleArray&       buffer_y,
                                                    TDoubleArray&       buffer_z,
                                                    const TPoint3D&     dipole,
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

}  // namespace npotg010rec

#endif /* PrimitiveNuclearPotentialGeom010FF_XYY_YZZ */
