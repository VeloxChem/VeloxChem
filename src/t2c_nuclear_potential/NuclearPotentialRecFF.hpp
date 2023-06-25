#ifndef NuclearPotentialRecFF_hpp
#define NuclearPotentialRecFF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace npotrec { // npotrec namespace

/**
 Evaluates <F|A|F>  integrals for given GTOs block.

 @param matrix the pointer to matrix for storage of integrals.
 @param gto_block the GTOs block.
*/
auto
compNuclearPotentialFF(      CSubMatrix* matrix,
                       const CGtoBlock&  gto_block,
                       const int64_t     bra_first,
                       const int64_t     bra_last) -> void;

/**
 Evaluates <F|A|F>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param mat_type the matrix type.
*/
auto
compNuclearPotentialFF(      CSubMatrix* matrix,
                       const CGtoBlock&  bra_gto_block,
                       const CGtoBlock&  ket_gto_block,
                       const int64_t     bra_first,
                       const int64_t     bra_last,
                       const mat_t       mat_type) -> void;

/**
 Evaluates block of primitive <F_XXX|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_XXX_T(      TDoubleArray& buffer_xxx,
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
 Evaluates block of primitive <F_XXY|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_XXY_T(      TDoubleArray& buffer_xxx,
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
 Evaluates block of primitive <F_XXZ|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_XXZ_T(      TDoubleArray& buffer_xxx,
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
 Evaluates block of primitive <F_XYY|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_XYY_T(      TDoubleArray& buffer_xxx,
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
 Evaluates block of primitive <F_XYZ|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_XYZ_T(      TDoubleArray& buffer_xxx,
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
 Evaluates block of primitive <F_XZZ|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_XZZ_T(      TDoubleArray& buffer_xxx,
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
 Evaluates block of primitive <F_YYY|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_YYY_T(      TDoubleArray& buffer_xxx,
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
 Evaluates block of primitive <F_YYZ|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_YYZ_T(      TDoubleArray& buffer_xxx,
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
 Evaluates block of primitive <F_YZZ|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_YZZ_T(      TDoubleArray& buffer_xxx,
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
 Evaluates block of primitive <F_ZZZ|A|F>  integrals.

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
compPrimitiveNuclearPotentialFF_ZZZ_T(      TDoubleArray& buffer_xxx,
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

} // npotrec namespace

#endif /* NuclearPotentialRecFF_hpp */
