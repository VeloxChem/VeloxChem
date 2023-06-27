#ifndef NuclearPotentialGeom010RecFG_hpp
#define NuclearPotentialGeom010RecFG_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <F||G>  integrals for given pair of GTOs blocks.

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
compNuclearPotentialGeom010FG(      CSubMatrix* matrix_x,
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
 Evaluates block of primitive <F_XXX||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXX||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXX_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXY_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XXZ_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYY_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XYZ_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_XZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYY_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YYZ_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_YZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XXXX(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XXXY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XXXZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XXYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XXYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XXZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_XZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_YYYY(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_YYYZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_YYZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_YZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                                const double        bra_exp,
                                                const double        bra_norm,
                                                const TPoint3D&     bra_coord,
                                                const TDoubleArray& ket_exps,
                                                const TDoubleArray& ket_norms,
                                                const TDoubleArray& ket_coords_x,
                                                const TDoubleArray& ket_coords_y,
                                                const TDoubleArray& ket_coords_z,
                                                const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010FG_ZZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                      TDoubleArray& buffer_y,
                                                      TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

#endif /* NuclearPotentialGeom010RecFG_hpp */
