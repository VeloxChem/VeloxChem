#ifndef PrimitiveOctupoleDD_XZ_XY
#define PrimitiveOctupoleDD_XZ_XY

#include <cstdint>

#include "SimdTypes.hpp"
#include "Point.hpp"

namespace octurec { // octurec namespace

/**
 Evaluates block of primitive <D_XZ|r^3|D_XY> integrals.

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
auto
compPrimitiveOctupoleDD_XZ_XY(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
               const TPoint3D& point,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

} // octurec namespace

#endif /* PrimitiveOctupoleDD_XZ_XY */
