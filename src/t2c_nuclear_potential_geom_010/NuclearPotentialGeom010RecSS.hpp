#ifndef NuclearPotentialGeom010RecSS_hpp
#define NuclearPotentialGeom010RecSS_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "Point.hpp"
#include "SubMatrix.hpp"

namespace geom_npotrec {  // geom_npotrec namespace

/**
 Evaluates <S||S>  integrals for given GTOs block.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compNuclearPotentialGeom010SS(CSubMatrix*      matrix_x,
                                   CSubMatrix*      matrix_y,
                                   CSubMatrix*      matrix_z,
                                   const TPoint3D&  dipole,
                                   const TPoint3D&  point,
                                   const CGtoBlock& gto_block,
                                   const int64_t    bra_first,
                                   const int64_t    bra_last) -> void;

/**
 Evaluates <S||S>  integrals for given pair of GTOs blocks.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto compNuclearPotentialGeom010SS(CSubMatrix*      matrix_x,
                                   CSubMatrix*      matrix_y,
                                   CSubMatrix*      matrix_z,
                                   const TPoint3D&  dipole,
                                   const TPoint3D&  point,
                                   const CGtoBlock& bra_gto_block,
                                   const CGtoBlock& ket_gto_block,
                                   const int64_t    bra_first,
                                   const int64_t    bra_last,
                                   const mat_t      mat_type) -> void;

}  // namespace geom_npotrec

#endif /* NuclearPotentialGeom010RecSS_hpp */
