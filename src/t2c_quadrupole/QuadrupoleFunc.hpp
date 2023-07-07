#ifndef QuadrupoleFunc_hpp
#define QuadrupoleFunc_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "Point.hpp"
#include "SubMatrix.hpp"

namespace quadfunc {  // quadfunc namespace

/**
 Computes quadrupole integrals for given basis functions block.

 @param matrix_xx the pointer to matrix for storage of Cartesian integral component XX.
 @param matrix_xy the pointer to matrix for storage of Cartesian integral component XY.
 @param matrix_xz the pointer to matrix for storage of Cartesian integral component XZ.
 @param matrix_yy the pointer to matrix for storage of Cartesian integral component YY.
 @param matrix_yz the pointer to matrix for storage of Cartesian integral component YZ.
 @param matrix_zz the pointer to matrix for storage of Cartesian integral component ZZ.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param angmom the angular momentum of bra and ket sides of integrals.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 */
auto compute(CSubMatrix*      matrix_xx,
             CSubMatrix*      matrix_xy,
             CSubMatrix*      matrix_xz,
             CSubMatrix*      matrix_yy,
             CSubMatrix*      matrix_yz,
             CSubMatrix*      matrix_zz,
             const TPoint3D&  point,
             const CGtoBlock& gto_block,
             const int64_t    angmom,
             const int64_t    bra_first,
             const int64_t    bra_last) -> void;

/**
 Computes quadrupole integrals for given of pair basis functions blocks.

 @param matrix_xx the pointer to matrix for storage of Cartesian integral component XX.
 @param matrix_xy the pointer to matrix for storage of Cartesian integral component XY.
 @param matrix_xz the pointer to matrix for storage of Cartesian integral component XZ.
 @param matrix_yy the pointer to matrix for storage of Cartesian integral component YY.
 @param matrix_yz the pointer to matrix for storage of Cartesian integral component YZ.
 @param matrix_zz the pointer to matrix for storage of Cartesian integral component ZZ.
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
auto compute(CSubMatrix*      matrix_xx,
             CSubMatrix*      matrix_xy,
             CSubMatrix*      matrix_xz,
             CSubMatrix*      matrix_yy,
             CSubMatrix*      matrix_yz,
             CSubMatrix*      matrix_zz,
             const TPoint3D&  point,
             const CGtoBlock& bra_gto_block,
             const CGtoBlock& ket_gto_block,
             const int64_t    bra_angmom,
             const int64_t    ket_angmom,
             const bool       ang_order,
             const int64_t    bra_first,
             const int64_t    bra_last,
             const mat_t      mat_type) -> void;

}  // namespace quadfunc

#endif /* QuadrupoleFunc_hpp */
