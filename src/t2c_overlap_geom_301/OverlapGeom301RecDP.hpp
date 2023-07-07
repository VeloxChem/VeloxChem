#ifndef OverlapGeom301RecDP_hpp
#define OverlapGeom301RecDP_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <d^(3)/dA^(3)D|1|d^(1)/dB^(1)P>  integrals for given pair of GTOs blocks.

 @param matrix_xxx_x the pointer to matrix for storage of Cartesian integral component XXX_X.
 @param matrix_xxx_y the pointer to matrix for storage of Cartesian integral component XXX_Y.
 @param matrix_xxx_z the pointer to matrix for storage of Cartesian integral component XXX_Z.
 @param matrix_xxy_x the pointer to matrix for storage of Cartesian integral component XXY_X.
 @param matrix_xxy_y the pointer to matrix for storage of Cartesian integral component XXY_Y.
 @param matrix_xxy_z the pointer to matrix for storage of Cartesian integral component XXY_Z.
 @param matrix_xxz_x the pointer to matrix for storage of Cartesian integral component XXZ_X.
 @param matrix_xxz_y the pointer to matrix for storage of Cartesian integral component XXZ_Y.
 @param matrix_xxz_z the pointer to matrix for storage of Cartesian integral component XXZ_Z.
 @param matrix_xyy_x the pointer to matrix for storage of Cartesian integral component XYY_X.
 @param matrix_xyy_y the pointer to matrix for storage of Cartesian integral component XYY_Y.
 @param matrix_xyy_z the pointer to matrix for storage of Cartesian integral component XYY_Z.
 @param matrix_xyz_x the pointer to matrix for storage of Cartesian integral component XYZ_X.
 @param matrix_xyz_y the pointer to matrix for storage of Cartesian integral component XYZ_Y.
 @param matrix_xyz_z the pointer to matrix for storage of Cartesian integral component XYZ_Z.
 @param matrix_xzz_x the pointer to matrix for storage of Cartesian integral component XZZ_X.
 @param matrix_xzz_y the pointer to matrix for storage of Cartesian integral component XZZ_Y.
 @param matrix_xzz_z the pointer to matrix for storage of Cartesian integral component XZZ_Z.
 @param matrix_yyy_x the pointer to matrix for storage of Cartesian integral component YYY_X.
 @param matrix_yyy_y the pointer to matrix for storage of Cartesian integral component YYY_Y.
 @param matrix_yyy_z the pointer to matrix for storage of Cartesian integral component YYY_Z.
 @param matrix_yyz_x the pointer to matrix for storage of Cartesian integral component YYZ_X.
 @param matrix_yyz_y the pointer to matrix for storage of Cartesian integral component YYZ_Y.
 @param matrix_yyz_z the pointer to matrix for storage of Cartesian integral component YYZ_Z.
 @param matrix_yzz_x the pointer to matrix for storage of Cartesian integral component YZZ_X.
 @param matrix_yzz_y the pointer to matrix for storage of Cartesian integral component YZZ_Y.
 @param matrix_yzz_z the pointer to matrix for storage of Cartesian integral component YZZ_Z.
 @param matrix_zzz_x the pointer to matrix for storage of Cartesian integral component ZZZ_X.
 @param matrix_zzz_y the pointer to matrix for storage of Cartesian integral component ZZZ_Y.
 @param matrix_zzz_z the pointer to matrix for storage of Cartesian integral component ZZZ_Z.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compOverlapGeom301DP(CSubMatrix*      matrix_xxx_x,
                          CSubMatrix*      matrix_xxx_y,
                          CSubMatrix*      matrix_xxx_z,
                          CSubMatrix*      matrix_xxy_x,
                          CSubMatrix*      matrix_xxy_y,
                          CSubMatrix*      matrix_xxy_z,
                          CSubMatrix*      matrix_xxz_x,
                          CSubMatrix*      matrix_xxz_y,
                          CSubMatrix*      matrix_xxz_z,
                          CSubMatrix*      matrix_xyy_x,
                          CSubMatrix*      matrix_xyy_y,
                          CSubMatrix*      matrix_xyy_z,
                          CSubMatrix*      matrix_xyz_x,
                          CSubMatrix*      matrix_xyz_y,
                          CSubMatrix*      matrix_xyz_z,
                          CSubMatrix*      matrix_xzz_x,
                          CSubMatrix*      matrix_xzz_y,
                          CSubMatrix*      matrix_xzz_z,
                          CSubMatrix*      matrix_yyy_x,
                          CSubMatrix*      matrix_yyy_y,
                          CSubMatrix*      matrix_yyy_z,
                          CSubMatrix*      matrix_yyz_x,
                          CSubMatrix*      matrix_yyz_y,
                          CSubMatrix*      matrix_yyz_z,
                          CSubMatrix*      matrix_yzz_x,
                          CSubMatrix*      matrix_yzz_y,
                          CSubMatrix*      matrix_yzz_z,
                          CSubMatrix*      matrix_zzz_x,
                          CSubMatrix*      matrix_zzz_y,
                          CSubMatrix*      matrix_zzz_z,
                          const CGtoBlock& bra_gto_block,
                          const CGtoBlock& ket_gto_block,
                          const bool       ang_order,
                          const int64_t    bra_first,
                          const int64_t    bra_last) -> void;

}  // namespace ovlrec

#endif /* OverlapGeom301RecDP_hpp */
