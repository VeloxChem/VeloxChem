#ifndef NuclearPotentialRecGG_hpp
#define NuclearPotentialRecGG_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace npotrec { // npotrec namespace

/**
 Evaluates <G|A|G>  integrals for given GTOs block.

 @param matrix the pointer to matrix for storage of integrals.
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compNuclearPotentialGG(      CSubMatrix* matrix,
                       const double charge,
                       const TPoint3D& point,
                       const CGtoBlock&  gto_block,
                       const int64_t     bra_first,
                       const int64_t     bra_last) -> void;

/**
 Evaluates <G|A|G>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto
compNuclearPotentialGG(      CSubMatrix* matrix,
                       const double charge,
                       const TPoint3D& point,
                       const CGtoBlock&  bra_gto_block,
                       const CGtoBlock&  ket_gto_block,
                       const int64_t     bra_first,
                       const int64_t     bra_last,
                       const mat_t       mat_type) -> void;

/**
 Evaluates block of primitive <G_XXXX|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XXXX_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XXXY_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XXXZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XXYY_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XXYZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XXZZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XYYY_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XYYZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XYZZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_XZZZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_YYYY_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_YYYZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_YYZZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_YZZZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ|A|G>  integrals.

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
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGG_ZZZZ_T(      TDoubleArray& buffer_xxxx,
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
                       const double charge,
                       const TPoint3D& point,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

} // npotrec namespace

#endif /* NuclearPotentialRecGG_hpp */
