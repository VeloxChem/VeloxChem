#ifndef QuadrupoleRecFF_hpp
#define QuadrupoleRecFF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SubMatrix.hpp"

namespace mpol {  // mpol namespace

/**
 Evaluates <F|r^2|F>  integrals for given GTOs block.

 @param matrix_xx the pointer to matrix for storage of Cartesian integral component XX.
 @param matrix_xy the pointer to matrix for storage of Cartesian integral component XY.
 @param matrix_xz the pointer to matrix for storage of Cartesian integral component XZ.
 @param matrix_yy the pointer to matrix for storage of Cartesian integral component YY.
 @param matrix_yz the pointer to matrix for storage of Cartesian integral component YZ.
 @param matrix_zz the pointer to matrix for storage of Cartesian integral component ZZ.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compQuadrupoleFF(CSubMatrix*      matrix_xx,
                      CSubMatrix*      matrix_xy,
                      CSubMatrix*      matrix_xz,
                      CSubMatrix*      matrix_yy,
                      CSubMatrix*      matrix_yz,
                      CSubMatrix*      matrix_zz,
                      const TPoint3D&  point,
                      const CGtoBlock& gto_block,
                      const int64_t    bra_first,
                      const int64_t    bra_last) -> void;

/**
 Evaluates <F|r^2|F>  integrals for given pair of GTOs blocks.

 @param matrix_xx the pointer to matrix for storage of Cartesian integral component XX.
 @param matrix_xy the pointer to matrix for storage of Cartesian integral component XY.
 @param matrix_xz the pointer to matrix for storage of Cartesian integral component XZ.
 @param matrix_yy the pointer to matrix for storage of Cartesian integral component YY.
 @param matrix_yz the pointer to matrix for storage of Cartesian integral component YZ.
 @param matrix_zz the pointer to matrix for storage of Cartesian integral component ZZ.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto compQuadrupoleFF(CSubMatrix*      matrix_xx,
                      CSubMatrix*      matrix_xy,
                      CSubMatrix*      matrix_xz,
                      CSubMatrix*      matrix_yy,
                      CSubMatrix*      matrix_yz,
                      CSubMatrix*      matrix_zz,
                      const TPoint3D&  point,
                      const CGtoBlock& bra_gto_block,
                      const CGtoBlock& ket_gto_block,
                      const int64_t    bra_first,
                      const int64_t    bra_last,
                      const mat_t      mat_type) -> void;

}  // namespace mpol

#endif /* QuadrupoleRecFF_hpp */
