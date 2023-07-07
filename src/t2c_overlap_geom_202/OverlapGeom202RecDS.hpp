#ifndef OverlapGeom202RecDS_hpp
#define OverlapGeom202RecDS_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <d^(2)/dA^(2)D|1|d^(2)/dB^(2)S>  integrals for given pair of GTOs blocks.

 @param matrix_xx_xx the pointer to matrix for storage of Cartesian integral component XX_XX.
 @param matrix_xx_xy the pointer to matrix for storage of Cartesian integral component XX_XY.
 @param matrix_xx_xz the pointer to matrix for storage of Cartesian integral component XX_XZ.
 @param matrix_xx_yy the pointer to matrix for storage of Cartesian integral component XX_YY.
 @param matrix_xx_yz the pointer to matrix for storage of Cartesian integral component XX_YZ.
 @param matrix_xx_zz the pointer to matrix for storage of Cartesian integral component XX_ZZ.
 @param matrix_xy_xx the pointer to matrix for storage of Cartesian integral component XY_XX.
 @param matrix_xy_xy the pointer to matrix for storage of Cartesian integral component XY_XY.
 @param matrix_xy_xz the pointer to matrix for storage of Cartesian integral component XY_XZ.
 @param matrix_xy_yy the pointer to matrix for storage of Cartesian integral component XY_YY.
 @param matrix_xy_yz the pointer to matrix for storage of Cartesian integral component XY_YZ.
 @param matrix_xy_zz the pointer to matrix for storage of Cartesian integral component XY_ZZ.
 @param matrix_xz_xx the pointer to matrix for storage of Cartesian integral component XZ_XX.
 @param matrix_xz_xy the pointer to matrix for storage of Cartesian integral component XZ_XY.
 @param matrix_xz_xz the pointer to matrix for storage of Cartesian integral component XZ_XZ.
 @param matrix_xz_yy the pointer to matrix for storage of Cartesian integral component XZ_YY.
 @param matrix_xz_yz the pointer to matrix for storage of Cartesian integral component XZ_YZ.
 @param matrix_xz_zz the pointer to matrix for storage of Cartesian integral component XZ_ZZ.
 @param matrix_yy_xx the pointer to matrix for storage of Cartesian integral component YY_XX.
 @param matrix_yy_xy the pointer to matrix for storage of Cartesian integral component YY_XY.
 @param matrix_yy_xz the pointer to matrix for storage of Cartesian integral component YY_XZ.
 @param matrix_yy_yy the pointer to matrix for storage of Cartesian integral component YY_YY.
 @param matrix_yy_yz the pointer to matrix for storage of Cartesian integral component YY_YZ.
 @param matrix_yy_zz the pointer to matrix for storage of Cartesian integral component YY_ZZ.
 @param matrix_yz_xx the pointer to matrix for storage of Cartesian integral component YZ_XX.
 @param matrix_yz_xy the pointer to matrix for storage of Cartesian integral component YZ_XY.
 @param matrix_yz_xz the pointer to matrix for storage of Cartesian integral component YZ_XZ.
 @param matrix_yz_yy the pointer to matrix for storage of Cartesian integral component YZ_YY.
 @param matrix_yz_yz the pointer to matrix for storage of Cartesian integral component YZ_YZ.
 @param matrix_yz_zz the pointer to matrix for storage of Cartesian integral component YZ_ZZ.
 @param matrix_zz_xx the pointer to matrix for storage of Cartesian integral component ZZ_XX.
 @param matrix_zz_xy the pointer to matrix for storage of Cartesian integral component ZZ_XY.
 @param matrix_zz_xz the pointer to matrix for storage of Cartesian integral component ZZ_XZ.
 @param matrix_zz_yy the pointer to matrix for storage of Cartesian integral component ZZ_YY.
 @param matrix_zz_yz the pointer to matrix for storage of Cartesian integral component ZZ_YZ.
 @param matrix_zz_zz the pointer to matrix for storage of Cartesian integral component ZZ_ZZ.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compOverlapGeom202DS(CSubMatrix*      matrix_xx_xx,
                          CSubMatrix*      matrix_xx_xy,
                          CSubMatrix*      matrix_xx_xz,
                          CSubMatrix*      matrix_xx_yy,
                          CSubMatrix*      matrix_xx_yz,
                          CSubMatrix*      matrix_xx_zz,
                          CSubMatrix*      matrix_xy_xx,
                          CSubMatrix*      matrix_xy_xy,
                          CSubMatrix*      matrix_xy_xz,
                          CSubMatrix*      matrix_xy_yy,
                          CSubMatrix*      matrix_xy_yz,
                          CSubMatrix*      matrix_xy_zz,
                          CSubMatrix*      matrix_xz_xx,
                          CSubMatrix*      matrix_xz_xy,
                          CSubMatrix*      matrix_xz_xz,
                          CSubMatrix*      matrix_xz_yy,
                          CSubMatrix*      matrix_xz_yz,
                          CSubMatrix*      matrix_xz_zz,
                          CSubMatrix*      matrix_yy_xx,
                          CSubMatrix*      matrix_yy_xy,
                          CSubMatrix*      matrix_yy_xz,
                          CSubMatrix*      matrix_yy_yy,
                          CSubMatrix*      matrix_yy_yz,
                          CSubMatrix*      matrix_yy_zz,
                          CSubMatrix*      matrix_yz_xx,
                          CSubMatrix*      matrix_yz_xy,
                          CSubMatrix*      matrix_yz_xz,
                          CSubMatrix*      matrix_yz_yy,
                          CSubMatrix*      matrix_yz_yz,
                          CSubMatrix*      matrix_yz_zz,
                          CSubMatrix*      matrix_zz_xx,
                          CSubMatrix*      matrix_zz_xy,
                          CSubMatrix*      matrix_zz_xz,
                          CSubMatrix*      matrix_zz_yy,
                          CSubMatrix*      matrix_zz_yz,
                          CSubMatrix*      matrix_zz_zz,
                          const CGtoBlock& bra_gto_block,
                          const CGtoBlock& ket_gto_block,
                          const bool       ang_order,
                          const int64_t    bra_first,
                          const int64_t    bra_last) -> void;

}  // namespace ovlrec

#endif /* OverlapGeom202RecDS_hpp */
