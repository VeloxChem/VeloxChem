#ifndef MatrixFunc_hpp
#define MatrixFunc_hpp

#include "MolecularBasis.hpp"
#include "Matrix.hpp"

namespace matfunc {  // matfunc namespace

/**
 Creates matrix.

 @param basis the molecular basis.
 @param mtype the matrix type.
 @return the matrix.
 */
auto
makeMatrix(const CMolecularBasis& basis,
           const mat_t            mtype) -> CMatrix;

/**
 Creates matrix.

 @param bra_basis the molecular basis on bra side.
 @param ket_basis the molecular basis on ket side.
 @return the matrix.
 */
auto
makeMatrix(const CMolecularBasis& bra_basis,
           const CMolecularBasis& ket_basis) -> CMatrix;

}  // namespace matfunc

#endif /* MatrixFunc_hpp */
