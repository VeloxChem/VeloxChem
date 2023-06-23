#ifndef OverlapRecPP_hpp
#define OverlapRecPP_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SimdTypes.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <P||P>  integrals for given GTOs block.

 @param matrix the pointer to matrix for storage of integrals.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compOverlapPP(CSubMatrix* matrix, const CGtoBlock& gto_block, const int64_t bra_first, const int64_t bra_last) -> void;

/**
 Evaluates <P||P>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto compOverlapPP(CSubMatrix*      matrix,
                   const CGtoBlock& bra_gto_block,
                   const CGtoBlock& ket_gto_block,
                   const int64_t    bra_first,
                   const int64_t    bra_last,
                   const mat_t      mat_type) -> void;

/**
 Evaluates block of primitive <P_X||P> integrals.

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
auto compPrimitiveOverlapPP_X_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <P_Y||P> integrals.

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
auto compPrimitiveOverlapPP_Y_T(TDoubleArray&       buffer_x,
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
 Evaluates block of primitive <P_Z||P> integrals.

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
auto compPrimitiveOverlapPP_Z_T(TDoubleArray&       buffer_x,
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

}  // namespace ovlrec

#endif /* OverlapRecPP_hpp */
