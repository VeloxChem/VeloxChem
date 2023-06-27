#ifndef NuclearPotentialGeom010RecFF_hpp
#define NuclearPotentialGeom010RecFF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <F||F>  integrals for given GTOs block.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compNuclearPotentialGeom010FF(      CSubMatrix* matrix_x,
                                    CSubMatrix* matrix_y,
                                    CSubMatrix* matrix_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                              const CGtoBlock&  gto_block,
                              const int64_t     bra_first,
                              const int64_t     bra_last) -> void;

/**
 Evaluates <F||F>  integrals for given pair of GTOs blocks.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto
compNuclearPotentialGeom010FF(      CSubMatrix* matrix_x,
                                    CSubMatrix* matrix_y,
                                    CSubMatrix* matrix_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                              const CGtoBlock&  bra_gto_block,
                              const CGtoBlock&  ket_gto_block,
                              const int64_t     bra_first,
                              const int64_t     bra_last,
                              const mat_t       mat_type) -> void;

/**
 Evaluates block of primitive <F_XXX||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXX||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXX||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXX||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXX||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXX||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXX||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXX||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXX||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXX||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXX_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXY||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXY_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XXZ||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XXZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYY||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYY_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XYZ||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XYZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_XZZ||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_XZZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYY||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYY_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YYZ||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YYZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_YZZ||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_YZZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_XXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_XXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_XXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_XYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_XYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_XZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_YYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_YYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_YZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <F_ZZZ||F_ZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010FF_ZZZ_ZZZ(      TDoubleArray& buffer_x,
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

#endif /* NuclearPotentialGeom010RecFF_hpp */
