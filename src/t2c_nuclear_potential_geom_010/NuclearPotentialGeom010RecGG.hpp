#ifndef NuclearPotentialGeom010RecGG_hpp
#define NuclearPotentialGeom010RecGG_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <G||G>  integrals for given GTOs block.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compNuclearPotentialGeom010GG(      CSubMatrix* matrix_x,
                                    CSubMatrix* matrix_y,
                                    CSubMatrix* matrix_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                              const CGtoBlock&  gto_block,
                              const int64_t     bra_first,
                              const int64_t     bra_last) -> void;

/**
 Evaluates <G||G>  integrals for given pair of GTOs blocks.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto
compNuclearPotentialGeom010GG(      CSubMatrix* matrix_x,
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
 Evaluates block of primitive <G_XXXX||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXX||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXX_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXY_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXXZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYY_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXYZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XXZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYY_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYYZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XYZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_XZZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYY_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYYZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YYZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_YZZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XXXX> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XXXX(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XXXY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XXXY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XXXZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XXXZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XXYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XXYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XXYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XXYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XXZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XXZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_XZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_XZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_YYYY> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_YYYY(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_YYYZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_YYYZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_YYZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_YYZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_YZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_YZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                 const double        bra_exp,
                                                 const double        bra_norm,
                                                 const TPoint3D&     bra_coord,
                                                 const TDoubleArray& ket_exps,
                                                 const TDoubleArray& ket_norms,
                                                 const TDoubleArray& ket_coords_x,
                                                 const TDoubleArray& ket_coords_y,
                                                 const TDoubleArray& ket_coords_z,
                                                 const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G_ZZZZ> integrals.

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
compPrimitiveNuclearPotentialGeom010GG_ZZZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                       TDoubleArray& buffer_y,
                                                       TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

#endif /* NuclearPotentialGeom010RecGG_hpp */
