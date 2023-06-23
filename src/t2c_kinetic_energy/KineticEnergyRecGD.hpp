#ifndef KineticEnergyRecGD_hpp
#define KineticEnergyRecGD_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SimdTypes.hpp"
#include "SubMatrix.hpp"

namespace kinrec {  // kinrec namespace

/**
 Evaluates <G|T|D>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
*/
auto compKineticEnergyGD(CSubMatrix*      matrix,
                         const CGtoBlock& bra_gto_block,
                         const CGtoBlock& ket_gto_block,
                         const bool       ang_order,
                         const int64_t    bra_first,
                         const int64_t    bra_last) -> void;

/**
 Evaluates block of primitive <G_XXXX|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XXXX_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XXXY_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XXXZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XXYY_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XXYZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XXZZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XYYY_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XYYZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XYZZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_XZZZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_YYYY_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_YYYZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_YYZZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_YZZZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyGD_ZZZZ_T(TDoubleArray&       buffer_xx,
                                         TDoubleArray&       buffer_xy,
                                         TDoubleArray&       buffer_xz,
                                         TDoubleArray&       buffer_yy,
                                         TDoubleArray&       buffer_yz,
                                         TDoubleArray&       buffer_zz,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

}  // namespace kinrec

#endif /* KineticEnergyRecGD_hpp */
