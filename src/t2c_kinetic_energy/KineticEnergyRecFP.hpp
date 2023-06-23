#ifndef KineticEnergyRecFP_hpp
#define KineticEnergyRecFP_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SimdTypes.hpp"
#include "SubMatrix.hpp"

namespace kinrec {  // kinrec namespace

/**
 Evaluates <F|T|P>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
*/
auto compKineticEnergyFP(CSubMatrix*      matrix,
                         const CGtoBlock& bra_gto_block,
                         const CGtoBlock& ket_gto_block,
                         const bool       ang_order,
                         const int64_t    bra_first,
                         const int64_t    bra_last) -> void;

/**
 Evaluates block of primitive <F_XXX|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_XXX_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <F_XXY|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_XXY_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <F_XXZ|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_XXZ_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <F_XYY|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_XYY_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <F_XYZ|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_XYZ_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <F_XZZ|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_XZZ_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <F_YYY|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_YYY_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <F_YYZ|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_YYZ_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <F_YZZ|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_YZZ_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <F_ZZZ|T|P>  integrals.

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
auto compPrimitiveKineticEnergyFP_ZZZ_T(TDoubleArray&       buffer_x,
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

#endif /* KineticEnergyRecFP_hpp */
