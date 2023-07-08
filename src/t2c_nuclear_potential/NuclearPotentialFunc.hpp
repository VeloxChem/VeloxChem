#ifndef NuclearPotentialFunc_hpp
#define NuclearPotentialFunc_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "Point.hpp"
#include "SubMatrix.hpp"

namespace npotfunc {  // npotfunc namespace

/**
 Computes nuclear potential for single external point integrals for given basis functions block.

 @param matrix the pointer to matrix for storage of integrals.
 @param charge the charge of external point.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param angmom the angular momentum of bra and ket sides of integrals.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 */
auto compute(CSubMatrix*      matrix,
             const double     charge,
             const TPoint3D&  point,
             const CGtoBlock& gto_block,
             const int64_t    angmom,
             const int64_t    bra_first,
             const int64_t    bra_last) -> void;

/**
 Computes nuclear potential integrals for given of pair basis functions blocks.

 @param matrix the pointer to matrix for storage of integrals.
 @param charge the charge of external point.
 @param point the coordinates of external point.
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
             const double     charge,
             const TPoint3D&  point,
             const CGtoBlock& bra_gto_block,
             const CGtoBlock& ket_gto_block,
             const int64_t    bra_angmom,
             const int64_t    ket_angmom,
             const bool       ang_order,
             const int64_t    bra_first,
             const int64_t    bra_last,
             const mat_t      mat_type) -> void;

}  // namespace npotfunc

#endif /* NuclearPotentialFunc_hpp */
