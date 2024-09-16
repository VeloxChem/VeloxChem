#ifndef MatrixFunc_hpp
#define MatrixFunc_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"

namespace matfunc {

/// @brief Creates a matrix.
/// @param basis The molecular basis to create matrix.
/// @param mtype  The matrix type of created matrix.
/// @return The matrix.
auto make_matrix(const CMolecularBasis& basis, const mat_t mtype) -> CMatrix;

/// @brief Creates a matrix.
/// @param bra_basis The molecular basis to create matrix (rows data).
/// @param ket_basis The molecular basis to create matrix (columns data).
/// @return The matrix.
auto make_matrix(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrix;

/// @brief Creates matrix.
/// @param label The exchange-correlation functional type label.
/// @param nrows The number of rows in submatrix.
/// @param ncols The number of columns in submatrix.
/// @return The matrix.
auto make_matrix(const std::string& label, const size_t nrows, const size_t ncols) -> CMatrix;

}  // namespace matfunc

#endif /* MatrixFunc_hpp */
