#ifndef KineticEnergyRecFF_hpp
#define KineticEnergyRecFF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SimdTypes.hpp"
#include "SubMatrix.hpp"

namespace kinrec {  // kinrec namespace

/**
 Evaluates <F|T|F>  integrals for given GTOs block.

 @param matrix the pointer to matrix for storage of integrals.
 @param gto_block the GTOs block.
*/
auto compKineticEnergyFF(CSubMatrix* matrix, const CGtoBlock& gto_block, const int64_t bra_first, const int64_t bra_last) -> void;

/**
 Evaluates <F|T|F>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param mat_type the matrix type.
*/
auto compKineticEnergyFF(CSubMatrix*      matrix,
                         const CGtoBlock& bra_gto_block,
                         const CGtoBlock& ket_gto_block,
                         const int64_t    bra_first,
                         const int64_t    bra_last,
                         const mat_t      mat_type) -> void;

/**
 Evaluates block of primitive <F_XXX|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_XXX_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
                                        const double        bra_exp,
                                        const double        bra_norm,
                                        const TPoint3D&     bra_coord,
                                        const TDoubleArray& ket_exps,
                                        const TDoubleArray& ket_norms,
                                        const TDoubleArray& ket_coords_x,
                                        const TDoubleArray& ket_coords_y,
                                        const TDoubleArray& ket_coords_z,
                                        const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXY|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_XXY_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
                                        const double        bra_exp,
                                        const double        bra_norm,
                                        const TPoint3D&     bra_coord,
                                        const TDoubleArray& ket_exps,
                                        const TDoubleArray& ket_norms,
                                        const TDoubleArray& ket_coords_x,
                                        const TDoubleArray& ket_coords_y,
                                        const TDoubleArray& ket_coords_z,
                                        const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XXZ|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_XXZ_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
                                        const double        bra_exp,
                                        const double        bra_norm,
                                        const TPoint3D&     bra_coord,
                                        const TDoubleArray& ket_exps,
                                        const TDoubleArray& ket_norms,
                                        const TDoubleArray& ket_coords_x,
                                        const TDoubleArray& ket_coords_y,
                                        const TDoubleArray& ket_coords_z,
                                        const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYY|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_XYY_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
                                        const double        bra_exp,
                                        const double        bra_norm,
                                        const TPoint3D&     bra_coord,
                                        const TDoubleArray& ket_exps,
                                        const TDoubleArray& ket_norms,
                                        const TDoubleArray& ket_coords_x,
                                        const TDoubleArray& ket_coords_y,
                                        const TDoubleArray& ket_coords_z,
                                        const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XYZ|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_XYZ_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
                                        const double        bra_exp,
                                        const double        bra_norm,
                                        const TPoint3D&     bra_coord,
                                        const TDoubleArray& ket_exps,
                                        const TDoubleArray& ket_norms,
                                        const TDoubleArray& ket_coords_x,
                                        const TDoubleArray& ket_coords_y,
                                        const TDoubleArray& ket_coords_z,
                                        const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_XZZ|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_XZZ_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
                                        const double        bra_exp,
                                        const double        bra_norm,
                                        const TPoint3D&     bra_coord,
                                        const TDoubleArray& ket_exps,
                                        const TDoubleArray& ket_norms,
                                        const TDoubleArray& ket_coords_x,
                                        const TDoubleArray& ket_coords_y,
                                        const TDoubleArray& ket_coords_z,
                                        const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYY|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_YYY_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
                                        const double        bra_exp,
                                        const double        bra_norm,
                                        const TPoint3D&     bra_coord,
                                        const TDoubleArray& ket_exps,
                                        const TDoubleArray& ket_norms,
                                        const TDoubleArray& ket_coords_x,
                                        const TDoubleArray& ket_coords_y,
                                        const TDoubleArray& ket_coords_z,
                                        const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YYZ|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_YYZ_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
                                        const double        bra_exp,
                                        const double        bra_norm,
                                        const TPoint3D&     bra_coord,
                                        const TDoubleArray& ket_exps,
                                        const TDoubleArray& ket_norms,
                                        const TDoubleArray& ket_coords_x,
                                        const TDoubleArray& ket_coords_y,
                                        const TDoubleArray& ket_coords_z,
                                        const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_YZZ|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_YZZ_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
                                        const double        bra_exp,
                                        const double        bra_norm,
                                        const TPoint3D&     bra_coord,
                                        const TDoubleArray& ket_exps,
                                        const TDoubleArray& ket_norms,
                                        const TDoubleArray& ket_coords_x,
                                        const TDoubleArray& ket_coords_y,
                                        const TDoubleArray& ket_coords_z,
                                        const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <F_ZZZ|T|F>  integrals.

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
auto compPrimitiveKineticEnergyFF_ZZZ_T(TDoubleArray&       buffer_xxx,
                                        TDoubleArray&       buffer_xxy,
                                        TDoubleArray&       buffer_xxz,
                                        TDoubleArray&       buffer_xyy,
                                        TDoubleArray&       buffer_xyz,
                                        TDoubleArray&       buffer_xzz,
                                        TDoubleArray&       buffer_yyy,
                                        TDoubleArray&       buffer_yyz,
                                        TDoubleArray&       buffer_yzz,
                                        TDoubleArray&       buffer_zzz,
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

#endif /* KineticEnergyRecFF_hpp */
