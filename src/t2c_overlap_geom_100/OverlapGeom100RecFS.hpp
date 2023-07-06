#ifndef OverlapGeom100RecFS_hpp
#define OverlapGeom100RecFS_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <d^(1)/dA^(1)F|1|S>  integrals for given pair of GTOs blocks.

 @param matrix_x the pointer to matrix for storage of Cartesian integral component X.
 @param matrix_y the pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_z the pointer to matrix for storage of Cartesian integral component Z.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compOverlapGeom100FS(CSubMatrix*      matrix_x,
                          CSubMatrix*      matrix_y,
                          CSubMatrix*      matrix_z,
                          const CGtoBlock& bra_gto_block,
                          const CGtoBlock& ket_gto_block,
                          const bool       ang_order,
                          const int64_t    bra_first,
                          const int64_t    bra_last) -> void;

}  // namespace ovlrec

#endif /* OverlapGeom100RecFS_hpp */
