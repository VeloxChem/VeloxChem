#ifndef NuclearPotentialGeom010RecPF_hpp
#define NuclearPotentialGeom010RecPF_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "Point.hpp"
#include "SubMatrix.hpp"

namespace npotg010rec {  // npotg010rec namespace

/**
 Evaluates <P|AG(1)|F>  integrals for given pair of GTOs blocks.

 @param matrix_x the pointer to matrix for storage of Cartesian integral component X.
 @param matrix_y the pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_z the pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto compNuclearPotentialGeom010PF(CSubMatrix*      matrix_x,
                                   CSubMatrix*      matrix_y,
                                   CSubMatrix*      matrix_z,
                                   const TPoint3D&  dipole,
                                   const TPoint3D&  point,
                                   const CGtoBlock& bra_gto_block,
                                   const CGtoBlock& ket_gto_block,
                                   const bool       ang_order,
                                   const int64_t    bra_first,
                                   const int64_t    bra_last) -> void;

}  // namespace npotg010rec

#endif /* NuclearPotentialGeom010RecPF_hpp */