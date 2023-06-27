#ifndef NuclearPotentialGeom010RecGS_hpp
#define NuclearPotentialGeom010RecGS_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <G||S>  integrals for given pair of GTOs blocks.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compNuclearPotentialGeom010GS(      CSubMatrix* matrix_x,
                                    CSubMatrix* matrix_y,
                                    CSubMatrix* matrix_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                              const CGtoBlock&  bra_gto_block,
                              const CGtoBlock&  ket_gto_block,
                              const bool        ang_order,
                              const int64_t     bra_first,
                              const int64_t     bra_last) -> void;

/**
 Evaluates block of primitive <G_XXXX||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XXXX_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_XXXY||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XXXY_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_XXXZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XXXZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_XXYY||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XXYY_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_XXYZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XXYZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_XXZZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XXZZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_XYYY||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XYYY_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_XYYZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XYYZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_XYZZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XYZZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_XZZZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_XZZZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_YYYY||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_YYYY_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_YYYZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_YYYZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_YYZZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_YYZZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_YZZZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_YZZZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

/**
 Evaluates block of primitive <G_ZZZZ||S_0> integrals.

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
auto
compPrimitiveNuclearPotentialGeom010GS_ZZZZ_0(      TDoubleArray& buffer_x,
                                                    TDoubleArray& buffer_y,
                                                    TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

} // geom_npotrec namespace

#endif /* NuclearPotentialGeom010RecGS_hpp */
