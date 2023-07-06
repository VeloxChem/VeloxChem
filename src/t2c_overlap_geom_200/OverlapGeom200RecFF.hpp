#ifndef OverlapGeom200RecFF_hpp
#define OverlapGeom200RecFF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <d^(2)/dA^(2)F|1|F>  integrals for given pair of GTOs blocks.

 @param matrix_xx the pointer to matrix for storage of Cartesian integral component XX.
 @param matrix_xy the pointer to matrix for storage of Cartesian integral component XY.
 @param matrix_xz the pointer to matrix for storage of Cartesian integral component XZ.
 @param matrix_yy the pointer to matrix for storage of Cartesian integral component YY.
 @param matrix_yz the pointer to matrix for storage of Cartesian integral component YZ.
 @param matrix_zz the pointer to matrix for storage of Cartesian integral component ZZ.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto compOverlapGeom200FF(CSubMatrix*      matrix_xx,
                          CSubMatrix*      matrix_xy,
                          CSubMatrix*      matrix_xz,
                          CSubMatrix*      matrix_yy,
                          CSubMatrix*      matrix_yz,
                          CSubMatrix*      matrix_zz,
                          const CGtoBlock& bra_gto_block,
                          const CGtoBlock& ket_gto_block,
                          const int64_t    bra_first,
                          const int64_t    bra_last,
                          const mat_t      mat_type) -> void;

}  // namespace ovlrec

#endif /* OverlapGeom200RecFF_hpp */
