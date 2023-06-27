#ifndef NuclearPotentialGeom010RecGF_hpp
#define NuclearPotentialGeom010RecGF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <G||F>  integrals for given pair of GTOs blocks.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compNuclearPotentialGeom010GF(      CSubMatrix* matrix_x,
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
 Evaluates block of primitive <G_XXXX||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXX_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXY_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXXZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYY_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXYZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XXZZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYY_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYYZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XYZZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_XZZZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYY_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYYZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YYZZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_YZZZ_ZZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_XXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_XXX(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_XXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_XXY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_XXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_XXZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_XYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_XYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_XYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_XYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_XZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_XZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_YYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_YYY(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_YYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_YYZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_YZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_YZZ(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||F_ZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GF_ZZZZ_ZZZ(      TDoubleArray& buffer_x,
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

#endif /* NuclearPotentialGeom010RecGF_hpp */
