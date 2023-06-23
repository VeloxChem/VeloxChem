#ifndef KineticEnergyRecDD_hpp
#define KineticEnergyRecDD_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SimdTypes.hpp"
#include "SubMatrix.hpp"

namespace kinrec {  // kinrec namespace

/**
 Evaluates <D|T|D>  integrals for given GTOs block.

 @param matrix the pointer to matrix for storage of integrals.
 @param gto_block the GTOs block.
*/
auto compKineticEnergyDD(CSubMatrix* matrix, const CGtoBlock& gto_block, const int64_t bra_first, const int64_t bra_last) -> void;

/**
 Evaluates <D|T|D>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param mat_type the matrix type.
*/
auto compKineticEnergyDD(CSubMatrix*      matrix,
                         const CGtoBlock& bra_gto_block,
                         const CGtoBlock& ket_gto_block,
                         const int64_t    bra_first,
                         const int64_t    bra_last,
                         const mat_t      mat_type) -> void;

/**
 Evaluates block of primitive <D_XX|T|D>  integrals.

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
auto compPrimitiveKineticEnergyDD_XX_T(TDoubleArray&       buffer_xx,
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
 Evaluates block of primitive <D_XY|T|D>  integrals.

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
auto compPrimitiveKineticEnergyDD_XY_T(TDoubleArray&       buffer_xx,
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
 Evaluates block of primitive <D_XZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyDD_XZ_T(TDoubleArray&       buffer_xx,
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
 Evaluates block of primitive <D_YY|T|D>  integrals.

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
auto compPrimitiveKineticEnergyDD_YY_T(TDoubleArray&       buffer_xx,
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
 Evaluates block of primitive <D_YZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyDD_YZ_T(TDoubleArray&       buffer_xx,
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
 Evaluates block of primitive <D_ZZ|T|D>  integrals.

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
auto compPrimitiveKineticEnergyDD_ZZ_T(TDoubleArray&       buffer_xx,
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

#endif /* KineticEnergyRecDD_hpp */
