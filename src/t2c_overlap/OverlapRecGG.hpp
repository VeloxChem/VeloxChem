#ifndef OverlapRecGG_hpp
#define OverlapRecGG_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace ovlrec { // ovlrec namespace

/**
 Evaluates <G||G>  integrals for given GTOs block.

 @param matrix the pointer to matrix for storage of integrals.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compOverlapGG(      CSubMatrix* matrix,
              const CGtoBlock&  gto_block,
              const int64_t     bra_first,
              const int64_t     bra_last) -> void;

/**
 Evaluates <G||G>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto
compOverlapGG(      CSubMatrix* matrix,
              const CGtoBlock&  bra_gto_block,
              const CGtoBlock&  ket_gto_block,
              const int64_t     bra_first,
              const int64_t     bra_last,
              const mat_t       mat_type) -> void;

/**
 Evaluates block of primitive <G_XXXX||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XXXX_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XXXY_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XXXZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XXYY_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XXYZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XXZZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XYYY_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XYYZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XYZZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_XZZZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_YYYY_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_YYYZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_YYZZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_YZZZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||G> integrals.

 @param buffer_xxxx the partial integrals buffer.
 @param buffer_xxxy the partial integrals buffer.
 @param buffer_xxxz the partial integrals buffer.
 @param buffer_xxyy the partial integrals buffer.
 @param buffer_xxyz the partial integrals buffer.
 @param buffer_xxzz the partial integrals buffer.
 @param buffer_xyyy the partial integrals buffer.
 @param buffer_xyyz the partial integrals buffer.
 @param buffer_xyzz the partial integrals buffer.
 @param buffer_xzzz the partial integrals buffer.
 @param buffer_yyyy the partial integrals buffer.
 @param buffer_yyyz the partial integrals buffer.
 @param buffer_yyzz the partial integrals buffer.
 @param buffer_yzzz the partial integrals buffer.
 @param buffer_zzzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGG_ZZZZ_T(      TDoubleArray& buffer_xxxx,
                                    TDoubleArray& buffer_xxxy,
                                    TDoubleArray& buffer_xxxz,
                                    TDoubleArray& buffer_xxyy,
                                    TDoubleArray& buffer_xxyz,
                                    TDoubleArray& buffer_xxzz,
                                    TDoubleArray& buffer_xyyy,
                                    TDoubleArray& buffer_xyyz,
                                    TDoubleArray& buffer_xyzz,
                                    TDoubleArray& buffer_xzzz,
                                    TDoubleArray& buffer_yyyy,
                                    TDoubleArray& buffer_yyyz,
                                    TDoubleArray& buffer_yyzz,
                                    TDoubleArray& buffer_yzzz,
                                    TDoubleArray& buffer_zzzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

} // ovlrec namespace

#endif /* OverlapRecGG_hpp */
