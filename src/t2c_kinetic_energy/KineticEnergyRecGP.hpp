#ifndef KineticEnergyRecGP_hpp
#define KineticEnergyRecGP_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SimdTypes.hpp"
#include "SubMatrix.hpp"

namespace kinrec {  // kinrec namespace

/**
 Evaluates <G|T|P>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
*/
auto compKineticEnergyGP(CSubMatrix*      matrix,
                         const CGtoBlock& bra_gto_block,
                         const CGtoBlock& ket_gto_block,
                         const bool       ang_order,
                         const int64_t    bra_first,
                         const int64_t    bra_last) -> void;

/**
 Evaluates block of primitive <G_XXXX|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XXXX_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXY|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XXXY_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXXZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XXXZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYY|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XXYY_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXYZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XXYZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XXZZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XXZZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYY|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XYYY_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYYZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XYYZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XYZZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XYZZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_XZZZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_XZZZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYY|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_YYYY_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYYZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_YYYZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YYZZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_YYZZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_YZZZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_YZZZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
                                         const double        bra_exp,
                                         const double        bra_norm,
                                         const TPoint3D&     bra_coord,
                                         const TDoubleArray& ket_exps,
                                         const TDoubleArray& ket_norms,
                                         const TDoubleArray& ket_coords_x,
                                         const TDoubleArray& ket_coords_y,
                                         const TDoubleArray& ket_coords_z,
                                         const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <G_ZZZZ|T|P>  integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto compPrimitiveKineticEnergyGP_ZZZZ_T(TDoubleArray&       buffer_x,
                                         TDoubleArray&       buffer_y,
                                         TDoubleArray&       buffer_z,
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

#endif /* KineticEnergyRecGP_hpp */
