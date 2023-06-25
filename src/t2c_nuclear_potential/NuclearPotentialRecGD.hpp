#ifndef NuclearPotentialRecGD_hpp
#define NuclearPotentialRecGD_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace npotrec { // npotrec namespace

/**
 Evaluates <G|A|D>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
*/
auto
compNuclearPotentialGD(      CSubMatrix* matrix,
                       const CGtoBlock&  bra_gto_block,
                       const CGtoBlock&  ket_gto_block,
                       const bool        ang_order,
                       const int64_t     bra_first,
                       const int64_t     bra_last) -> void;

/**
 Evaluates block of primitive <G_XXXX|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XXXX_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XXXY_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XXXZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XXYY_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XXYZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XXZZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XYYY_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XYYZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XYZZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_XZZZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_YYYY_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_YYYZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_YYZZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_YZZZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
                                       const double        bra_exp,
                                       const double        bra_norm,
                                       const TPoint3D&     bra_coord,
                                       const TDoubleArray& ket_exps,
                                       const TDoubleArray& ket_norms,
                                       const TDoubleArray& ket_coords_x,
                                       const TDoubleArray& ket_coords_y,
                                       const TDoubleArray& ket_coords_z,
                                       const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ|A|D>  integrals.

 @param buffer_xx the partial integrals buffer.
 @param buffer_xy the partial integrals buffer.
 @param buffer_xz the partial integrals buffer.
 @param buffer_yy the partial integrals buffer.
 @param buffer_yz the partial integrals buffer.
 @param buffer_zz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGD_ZZZZ_T(      TDoubleArray& buffer_xx,
                                             TDoubleArray& buffer_xy,
                                             TDoubleArray& buffer_xz,
                                             TDoubleArray& buffer_yy,
                                             TDoubleArray& buffer_yz,
                                             TDoubleArray& buffer_zz,
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

#endif /* NuclearPotentialRecGD_hpp */
