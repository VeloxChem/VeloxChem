#ifndef OctupoleRecSF_hpp
#define OctupoleRecSF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"

namespace mpol {  // mpol namespace

/**
 Evaluates <S|r^3|F>  integrals for given pair of GTOs blocks.

 @param matrix_xxx the pointer to matrix for storage of Cartesian integral component XXX.
 @param matrix_xxy the pointer to matrix for storage of Cartesian integral component XXY.
 @param matrix_xxz the pointer to matrix for storage of Cartesian integral component XXZ.
 @param matrix_xyy the pointer to matrix for storage of Cartesian integral component XYY.
 @param matrix_xyz the pointer to matrix for storage of Cartesian integral component XYZ.
 @param matrix_xzz the pointer to matrix for storage of Cartesian integral component XZZ.
 @param matrix_yyy the pointer to matrix for storage of Cartesian integral component YYY.
 @param matrix_yyz the pointer to matrix for storage of Cartesian integral component YYZ.
 @param matrix_yzz the pointer to matrix for storage of Cartesian integral component YZZ.
 @param matrix_zzz the pointer to matrix for storage of Cartesian integral component ZZZ.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compOctupoleSF(CSubMatrix*      matrix_xxx,
                    CSubMatrix*      matrix_xxy,
                    CSubMatrix*      matrix_xxz,
                    CSubMatrix*      matrix_xyy,
                    CSubMatrix*      matrix_xyz,
                    CSubMatrix*      matrix_xzz,
                    CSubMatrix*      matrix_yyy,
                    CSubMatrix*      matrix_yyz,
                    CSubMatrix*      matrix_yzz,
                    CSubMatrix*      matrix_zzz,
                    const TPoint3D&  point,
                    const CGtoBlock& bra_gto_block,
                    const CGtoBlock& ket_gto_block,
                    const bool       ang_order,
                    const int64_t    bra_first,
                    const int64_t    bra_last) -> void;

}  // namespace mpol

#endif /* OctupoleRecSF_hpp */
