#ifndef NuclearPotentialRecDD_hpp
#define NuclearPotentialRecDD_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "Point.hpp"
#include "SubMatrix.hpp"

namespace npotrec {  // npotrec namespace

/**
 Evaluates <D|A|D>  integrals for given GTOs block.

 @param matrix the pointer to matrix for storage of integrals.
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compNuclearPotentialDD(CSubMatrix*      matrix,
                            const double     charge,
                            const TPoint3D&  point,
                            const CGtoBlock& gto_block,
                            const int64_t    bra_first,
                            const int64_t    bra_last) -> void;

/**
 Evaluates <D|A|D>  integrals for given pair of GTOs blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto compNuclearPotentialDD(CSubMatrix*      matrix,
                            const double     charge,
                            const TPoint3D&  point,
                            const CGtoBlock& bra_gto_block,
                            const CGtoBlock& ket_gto_block,
                            const int64_t    bra_first,
                            const int64_t    bra_last,
                            const mat_t      mat_type) -> void;

}  // namespace npotrec

#endif /* NuclearPotentialRecDD_hpp */
