#ifndef OverlapRecGF_hpp
#define OverlapRecGF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace ovlrec { // ovlrec namespace

/**
 Evaluates <G||F>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compOverlapGF(      CSubMatrix* matrix,
              const CGtoBlock&  bra_gto_block,
              const CGtoBlock&  ket_gto_block,
              const bool        ang_order,
              const int64_t     bra_first,
              const int64_t     bra_last) -> void;

/**
 Evaluates block of primitive <G_XXXX||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XXXX_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XXXY_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XXXZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XXYY_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XXYZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XXZZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XYYY_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XYYZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XYZZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_XZZZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_YYYY_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_YYYZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_YYZZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_YZZZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ||F> integrals.

 @param buffer_xxx the partial integrals buffer.
 @param buffer_xxy the partial integrals buffer.
 @param buffer_xxz the partial integrals buffer.
 @param buffer_xyy the partial integrals buffer.
 @param buffer_xyz the partial integrals buffer.
 @param buffer_xzz the partial integrals buffer.
 @param buffer_yyy the partial integrals buffer.
 @param buffer_yyz the partial integrals buffer.
 @param buffer_yzz the partial integrals buffer.
 @param buffer_zzz the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveOverlapGF_ZZZZ_T(      TDoubleArray& buffer_xxx,
                                    TDoubleArray& buffer_xxy,
                                    TDoubleArray& buffer_xxz,
                                    TDoubleArray& buffer_xyy,
                                    TDoubleArray& buffer_xyz,
                                    TDoubleArray& buffer_xzz,
                                    TDoubleArray& buffer_yyy,
                                    TDoubleArray& buffer_yyz,
                                    TDoubleArray& buffer_yzz,
                                    TDoubleArray& buffer_zzz,
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

#endif /* OverlapRecGF_hpp */
