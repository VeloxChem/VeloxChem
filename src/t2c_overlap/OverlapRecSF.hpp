#ifndef OverlapRecSF_hpp
#define OverlapRecSF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <S||F>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compOverlapSF(CSubMatrix*      matrix,
                   const CGtoBlock& bra_gto_block,
                   const CGtoBlock& ket_gto_block,
                   const bool       ang_order,
                   const int64_t    bra_first,
                   const int64_t    bra_last) -> void;

}  // namespace ovlrec

#endif /* OverlapRecSF_hpp */
