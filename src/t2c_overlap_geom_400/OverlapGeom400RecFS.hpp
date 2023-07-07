#ifndef OverlapGeom400RecFS_hpp
#define OverlapGeom400RecFS_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"

namespace ovlrec {  // ovlrec namespace

/**
 Evaluates <d^(4)/dA^(4)F|1|S>  integrals for given pair of GTOs blocks.

 @param matrix_xxxx the pointer to matrix for storage of Cartesian integral component XXXX.
 @param matrix_xxxy the pointer to matrix for storage of Cartesian integral component XXXY.
 @param matrix_xxxz the pointer to matrix for storage of Cartesian integral component XXXZ.
 @param matrix_xxyy the pointer to matrix for storage of Cartesian integral component XXYY.
 @param matrix_xxyz the pointer to matrix for storage of Cartesian integral component XXYZ.
 @param matrix_xxzz the pointer to matrix for storage of Cartesian integral component XXZZ.
 @param matrix_xyyy the pointer to matrix for storage of Cartesian integral component XYYY.
 @param matrix_xyyz the pointer to matrix for storage of Cartesian integral component XYYZ.
 @param matrix_xyzz the pointer to matrix for storage of Cartesian integral component XYZZ.
 @param matrix_xzzz the pointer to matrix for storage of Cartesian integral component XZZZ.
 @param matrix_yyyy the pointer to matrix for storage of Cartesian integral component YYYY.
 @param matrix_yyyz the pointer to matrix for storage of Cartesian integral component YYYZ.
 @param matrix_yyzz the pointer to matrix for storage of Cartesian integral component YYZZ.
 @param matrix_yzzz the pointer to matrix for storage of Cartesian integral component YZZZ.
 @param matrix_zzzz the pointer to matrix for storage of Cartesian integral component ZZZZ.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compOverlapGeom400FS(CSubMatrix*      matrix_xxxx,
                          CSubMatrix*      matrix_xxxy,
                          CSubMatrix*      matrix_xxxz,
                          CSubMatrix*      matrix_xxyy,
                          CSubMatrix*      matrix_xxyz,
                          CSubMatrix*      matrix_xxzz,
                          CSubMatrix*      matrix_xyyy,
                          CSubMatrix*      matrix_xyyz,
                          CSubMatrix*      matrix_xyzz,
                          CSubMatrix*      matrix_xzzz,
                          CSubMatrix*      matrix_yyyy,
                          CSubMatrix*      matrix_yyyz,
                          CSubMatrix*      matrix_yyzz,
                          CSubMatrix*      matrix_yzzz,
                          CSubMatrix*      matrix_zzzz,
                          const CGtoBlock& bra_gto_block,
                          const CGtoBlock& ket_gto_block,
                          const bool       ang_order,
                          const int64_t    bra_first,
                          const int64_t    bra_last) -> void;

}  // namespace ovlrec

#endif /* OverlapGeom400RecFS_hpp */
