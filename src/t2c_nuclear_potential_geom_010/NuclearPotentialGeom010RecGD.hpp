#ifndef NuclearPotentialGeom010RecGD_hpp
#define NuclearPotentialGeom010RecGD_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <G||D>  integrals for given pair of GTOs blocks.

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
compNuclearPotentialGeom010GD(      CSubMatrix* matrix_x,
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
 Evaluates block of primitive <G_XXXX||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXX_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXX_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXX_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXX_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXX_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXX_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXY_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXY_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXY_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXY_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXY_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXY_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXXZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYY_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYY_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYY_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYY_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYY_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYY_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXYZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXZZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXZZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXZZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXZZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXZZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XXZZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYY_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYY_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYY_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYY_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYY_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYY_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYYZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYZZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYZZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYZZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYZZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYZZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XYZZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XZZZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XZZZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XZZZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XZZZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XZZZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_XZZZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYY_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYY_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYY_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYY_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYY_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYY_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYYZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYZZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYZZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYZZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYZZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYZZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YYZZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YZZZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YZZZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YZZZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YZZZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YZZZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_YZZZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||D_XX> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_ZZZZ_XX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||D_XY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_ZZZZ_XY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||D_XZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_ZZZZ_XZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||D_YY> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_ZZZZ_YY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||D_YZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_ZZZZ_YZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||D_ZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GD_ZZZZ_ZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

#endif /* NuclearPotentialGeom010RecGD_hpp */
