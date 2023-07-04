#ifndef DipoleRecFF_hpp
#define DipoleRecFF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "SubMatrix.hpp"

namespace mpol {  // mpol namespace

/**
 Evaluates <F|r|F>  integrals for given GTOs block.

 @param matrix_x the pointer to matrix for storage of Cartesian integral component X.
 @param matrix_y the pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_z the pointer to matrix for storage of Cartesian integral component Z.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compDipoleFF(CSubMatrix*      matrix_x,
                  CSubMatrix*      matrix_y,
                  CSubMatrix*      matrix_z,
                  const TPoint3D&  point,
                  const CGtoBlock& gto_block,
                  const int64_t    bra_first,
                  const int64_t    bra_last) -> void;

/**
 Evaluates <F|r|F>  integrals for given pair of GTOs blocks.

 @param matrix_x the pointer to matrix for storage of Cartesian integral component X.
 @param matrix_y the pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_z the pointer to matrix for storage of Cartesian integral component Z.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto compDipoleFF(CSubMatrix*      matrix_x,
                  CSubMatrix*      matrix_y,
                  CSubMatrix*      matrix_z,
                  const TPoint3D&  point,
                  const CGtoBlock& bra_gto_block,
                  const CGtoBlock& ket_gto_block,
                  const int64_t    bra_first,
                  const int64_t    bra_last,
                  const mat_t      mat_type) -> void;

}  // namespace mpol

#endif /* DipoleRecFF_hpp */
