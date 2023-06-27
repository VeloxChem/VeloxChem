#ifndef NuclearPotentialGeom010RecDG_hpp
#define NuclearPotentialGeom010RecDG_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <D||G>  integrals for given pair of GTOs blocks.

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
compNuclearPotentialGeom010DG(      CSubMatrix* matrix_x,
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
 Evaluates block of primitive <D_XX||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XXXX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XXXY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XXXZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XXYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XXYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XXZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_XZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_YYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_YYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_YYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_YZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XX_ZZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XXXX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XXXY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XXXZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XXYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XXYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XXZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_XZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_YYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_YYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_YYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_YZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XY_ZZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XXXX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XXXY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XXXZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XXYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XXYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XXZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_XZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_YYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_YYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_YYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_YZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_XZ_ZZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XXXX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XXXY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XXXZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XXYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XXYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XXZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_XZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_YYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_YYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_YYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_YZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YY_ZZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XXXX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XXXY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XXXZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XXYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XXYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XXZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_XZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_YYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_YYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_YYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_YZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_YZ_ZZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XXXX> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XXXX(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XXXY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XXXY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XXXZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XXXZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XXYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XXYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XXYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XXYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XXZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XXZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_XZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_XZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_YYYY> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_YYYY(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_YYYZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_YYYZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_YYZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_YYZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_YZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_YZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                               const double        bra_exp,
                                               const double        bra_norm,
                                               const TPoint3D&     bra_coord,
                                               const TDoubleArray& ket_exps,
                                               const TDoubleArray& ket_norms,
                                               const TDoubleArray& ket_coords_x,
                                               const TDoubleArray& ket_coords_y,
                                               const TDoubleArray& ket_coords_z,
                                               const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||G_ZZZZ> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DG_ZZ_ZZZZ(      TDoubleArray& buffer_x,
                                                     TDoubleArray& buffer_y,
                                                     TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
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

#endif /* NuclearPotentialGeom010RecDG_hpp */
