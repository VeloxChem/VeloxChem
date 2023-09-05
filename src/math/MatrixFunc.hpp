#ifndef MatrixFunc_hpp
#define MatrixFunc_hpp

#include <cstdint>
#include <string>

#include "Matrix.hpp"
#include "MolecularBasis.hpp"

namespace matfunc {  // matfunc namespace

/**
 Creates matrix.

 @param basis the molecular basis.
 @param mtype the matrix type.
 @return the matrix.
 */
auto makeMatrix(const CMolecularBasis& basis, const mat_t mtype) -> CMatrix;

/**
 Creates matrix.

 @param bra_basis the molecular basis on bra side.
 @param ket_basis the molecular basis on ket side.
 @return the matrix.
 */
auto makeMatrix(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrix;


/**
 Creates matrix.

 @param label the exchange-correlation functional type label.
 @param nrows the number of rows in submatrix.
 @param ncols the number of columns in submatrix.
 @return the matrix.
 */
auto makeMatrix(const std::string& label, const int64_t nrows, const int64_t ncols) -> CMatrix;

}  // namespace matfunc

#endif /* MatrixFunc_hpp */
