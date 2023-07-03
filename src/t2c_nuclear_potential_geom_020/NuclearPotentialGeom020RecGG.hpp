#ifndef NuclearPotentialGeom020RecGG_hpp
#define NuclearPotentialGeom020RecGG_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MatrixType.hpp"
#include "Point.hpp"
#include "TensorTypes.hpp"
#include "SubMatrix.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <G||G>  integrals for given GTOs block.

 @param matrix_xxthe pointer to matrix for storage of Cartesian integral component XX.
 @param matrix_xythe pointer to matrix for storage of Cartesian integral component XY.
 @param matrix_xzthe pointer to matrix for storage of Cartesian integral component XZ.
 @param matrix_yythe pointer to matrix for storage of Cartesian integral component YY.
 @param matrix_yzthe pointer to matrix for storage of Cartesian integral component YZ.
 @param matrix_zzthe pointer to matrix for storage of Cartesian integral component ZZ.
 @param quadrupole the quadrupole of external point.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compNuclearPotentialGeom020GG(      CSubMatrix* matrix_xx,
                                    CSubMatrix* matrix_xy,
                                    CSubMatrix* matrix_xz,
                                    CSubMatrix* matrix_yy,
                                    CSubMatrix* matrix_yz,
                                    CSubMatrix* matrix_zz,
                              const T2Tensor& quadrupole,
                              const TPoint3D& point,
                              const CGtoBlock&  gto_block,
                              const int64_t     bra_first,
                              const int64_t     bra_last) -> void;

/**
 Evaluates <G||G>  integrals for given pair of GTOs blocks.

 @param matrix_xxthe pointer to matrix for storage of Cartesian integral component XX.
 @param matrix_xythe pointer to matrix for storage of Cartesian integral component XY.
 @param matrix_xzthe pointer to matrix for storage of Cartesian integral component XZ.
 @param matrix_yythe pointer to matrix for storage of Cartesian integral component YY.
 @param matrix_yzthe pointer to matrix for storage of Cartesian integral component YZ.
 @param matrix_zzthe pointer to matrix for storage of Cartesian integral component ZZ.
 @param quadrupole the quadrupole of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto
compNuclearPotentialGeom020GG(      CSubMatrix* matrix_xx,
                                    CSubMatrix* matrix_xy,
                                    CSubMatrix* matrix_xz,
                                    CSubMatrix* matrix_yy,
                                    CSubMatrix* matrix_yz,
                                    CSubMatrix* matrix_zz,
                              const T2Tensor& quadrupole,
                              const TPoint3D& point,
                              const CGtoBlock&  bra_gto_block,
                              const CGtoBlock&  ket_gto_block,
                              const int64_t     bra_first,
                              const int64_t     bra_last,
                              const mat_t       mat_type) -> void;

} // geom_npotrec namespace

#endif /* NuclearPotentialGeom020RecGG_hpp */
