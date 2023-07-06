#ifndef OverlapGeom100RecFF_hpp
#define OverlapGeom100RecFF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <d^(1)/dA^(1)F|1|F>  integrals for given pair of GTOs blocks.

 @param matrix_x the pointer to matrix for storage of Cartesian integral component X.
 @param matrix_y the pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_z the pointer to matrix for storage of Cartesian integral component Z.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto compOverlapGeom100FF(CSubMatrix*      matrix_x,
                          CSubMatrix*      matrix_y,
                          CSubMatrix*      matrix_z,
                          const CGtoBlock& bra_gto_block,
                          const CGtoBlock& ket_gto_block,
                          const int64_t    bra_first,
                          const int64_t    bra_last,
                          const mat_t      mat_type) -> void;

}  // namespace ovlrec

#endif /* OverlapGeom100RecFF_hpp */
