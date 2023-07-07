#ifndef OverlapGeom201RecDS_hpp
#define OverlapGeom201RecDS_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <d^(2)/dA^(2)D|1|d^(1)/dB^(1)S>  integrals for given pair of GTOs blocks.

 @param matrix_xx_x the pointer to matrix for storage of Cartesian integral component XX_X.
 @param matrix_xx_y the pointer to matrix for storage of Cartesian integral component XX_Y.
 @param matrix_xx_z the pointer to matrix for storage of Cartesian integral component XX_Z.
 @param matrix_xy_x the pointer to matrix for storage of Cartesian integral component XY_X.
 @param matrix_xy_y the pointer to matrix for storage of Cartesian integral component XY_Y.
 @param matrix_xy_z the pointer to matrix for storage of Cartesian integral component XY_Z.
 @param matrix_xz_x the pointer to matrix for storage of Cartesian integral component XZ_X.
 @param matrix_xz_y the pointer to matrix for storage of Cartesian integral component XZ_Y.
 @param matrix_xz_z the pointer to matrix for storage of Cartesian integral component XZ_Z.
 @param matrix_yy_x the pointer to matrix for storage of Cartesian integral component YY_X.
 @param matrix_yy_y the pointer to matrix for storage of Cartesian integral component YY_Y.
 @param matrix_yy_z the pointer to matrix for storage of Cartesian integral component YY_Z.
 @param matrix_yz_x the pointer to matrix for storage of Cartesian integral component YZ_X.
 @param matrix_yz_y the pointer to matrix for storage of Cartesian integral component YZ_Y.
 @param matrix_yz_z the pointer to matrix for storage of Cartesian integral component YZ_Z.
 @param matrix_zz_x the pointer to matrix for storage of Cartesian integral component ZZ_X.
 @param matrix_zz_y the pointer to matrix for storage of Cartesian integral component ZZ_Y.
 @param matrix_zz_z the pointer to matrix for storage of Cartesian integral component ZZ_Z.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compOverlapGeom201DS(CSubMatrix*      matrix_xx_x,
                          CSubMatrix*      matrix_xx_y,
                          CSubMatrix*      matrix_xx_z,
                          CSubMatrix*      matrix_xy_x,
                          CSubMatrix*      matrix_xy_y,
                          CSubMatrix*      matrix_xy_z,
                          CSubMatrix*      matrix_xz_x,
                          CSubMatrix*      matrix_xz_y,
                          CSubMatrix*      matrix_xz_z,
                          CSubMatrix*      matrix_yy_x,
                          CSubMatrix*      matrix_yy_y,
                          CSubMatrix*      matrix_yy_z,
                          CSubMatrix*      matrix_yz_x,
                          CSubMatrix*      matrix_yz_y,
                          CSubMatrix*      matrix_yz_z,
                          CSubMatrix*      matrix_zz_x,
                          CSubMatrix*      matrix_zz_y,
                          CSubMatrix*      matrix_zz_z,
                          const CGtoBlock& bra_gto_block,
                          const CGtoBlock& ket_gto_block,
                          const bool       ang_order,
                          const int64_t    bra_first,
                          const int64_t    bra_last) -> void;

}  // namespace ovlrec

#endif /* OverlapGeom201RecDS_hpp */
