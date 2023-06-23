#ifndef OverlapFunc_hpp
#define OverlapFunc_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SubMatrix.hpp"

namespace ovlfunc {  // ovlfunc namespace

/**
 Computes overlap integrals for given basis functions block.

 @param matrix the pointer to matrix for storage of integrals.
 @param gto_block the GTOs block.
 @param angmom the angular momentum of bra and ket sides of integrals.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 */
auto compute(CSubMatrix* matrix, const CGtoBlock& gto_block, const int64_t angmom, const int64_t bra_first, const int64_t bra_last) -> void;

/**
 Computes overlap integrals for given of pair basis functions blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_angmom the angular momentum of bra side of integrals.
 @param ket_angmom the angular momentum of bra side of integrals.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the supermatrix type.
 */
auto compute(CSubMatrix*      matrix,
             const CGtoBlock& bra_gto_block,
             const CGtoBlock& ket_gto_block,
             const int64_t    bra_angmom,
             const int64_t    ket_angmom,
             const bool       ang_order,
             const int64_t    bra_first,
             const int64_t    bra_last,
             const mat_t      mat_type) -> void;

}  // namespace ovlfunc

#endif /* OverlapFunc_hpp */
