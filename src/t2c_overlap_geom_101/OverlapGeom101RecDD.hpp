#ifndef OverlapGeom101RecDD_hpp
#define OverlapGeom101RecDD_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <d^(1)/dA^(1)D|1|d^(1)/dB^(1)D>  integrals for given pair of GTOs blocks.

 @param matrix_x_x the pointer to matrix for storage of Cartesian integral component X_X.
 @param matrix_x_y the pointer to matrix for storage of Cartesian integral component X_Y.
 @param matrix_x_z the pointer to matrix for storage of Cartesian integral component X_Z.
 @param matrix_y_x the pointer to matrix for storage of Cartesian integral component Y_X.
 @param matrix_y_y the pointer to matrix for storage of Cartesian integral component Y_Y.
 @param matrix_y_z the pointer to matrix for storage of Cartesian integral component Y_Z.
 @param matrix_z_x the pointer to matrix for storage of Cartesian integral component Z_X.
 @param matrix_z_y the pointer to matrix for storage of Cartesian integral component Z_Y.
 @param matrix_z_z the pointer to matrix for storage of Cartesian integral component Z_Z.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto compOverlapGeom101DD(CSubMatrix*      matrix_x_x,
                          CSubMatrix*      matrix_x_y,
                          CSubMatrix*      matrix_x_z,
                          CSubMatrix*      matrix_y_x,
                          CSubMatrix*      matrix_y_y,
                          CSubMatrix*      matrix_y_z,
                          CSubMatrix*      matrix_z_x,
                          CSubMatrix*      matrix_z_y,
                          CSubMatrix*      matrix_z_z,
                          const CGtoBlock& bra_gto_block,
                          const CGtoBlock& ket_gto_block,
                          const int64_t    bra_first,
                          const int64_t    bra_last,
                          const mat_t      mat_type) -> void;

}  // namespace ovlrec

#endif /* OverlapGeom101RecDD_hpp */
