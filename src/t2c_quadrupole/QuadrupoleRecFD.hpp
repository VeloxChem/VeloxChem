#ifndef QuadrupoleRecFD_hpp
#define QuadrupoleRecFD_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"

namespace mpol {  // mpol namespace

/**
 Evaluates <F|r^2|D>  integrals for given pair of GTOs blocks.

 @param matrix_xx the pointer to matrix for storage of Cartesian integral component XX.
 @param matrix_xy the pointer to matrix for storage of Cartesian integral component XY.
 @param matrix_xz the pointer to matrix for storage of Cartesian integral component XZ.
 @param matrix_yy the pointer to matrix for storage of Cartesian integral component YY.
 @param matrix_yz the pointer to matrix for storage of Cartesian integral component YZ.
 @param matrix_zz the pointer to matrix for storage of Cartesian integral component ZZ.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compQuadrupoleFD(CSubMatrix*      matrix_xx,
                      CSubMatrix*      matrix_xy,
                      CSubMatrix*      matrix_xz,
                      CSubMatrix*      matrix_yy,
                      CSubMatrix*      matrix_yz,
                      CSubMatrix*      matrix_zz,
                      const TPoint3D&  point,
                      const CGtoBlock& bra_gto_block,
                      const CGtoBlock& ket_gto_block,
                      const bool       ang_order,
                      const int64_t    bra_first,
                      const int64_t    bra_last) -> void;

}  // namespace mpol

#endif /* QuadrupoleRecFD_hpp */
