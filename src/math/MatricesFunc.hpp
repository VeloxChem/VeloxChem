#ifndef MatricesFunc_hpp
#define MatricesFunc_hpp

#include <cstdint>
#include <vector>

#include "Matrices.hpp"
#include "MatrixType.hpp"
#include "MolecularBasis.hpp"

namespace matfunc {  // matfunc namespace

/**
 Creates matrices for given tensor order.

 @param torder the tensor order.
 @param basis the molecular basis.
 @param mtype the matrix type.
 @return the matrices.
 */
auto makeMatrices(const int64_t torder, const CMolecularBasis& basis, const mat_t mtype) -> CMatrices;

/**
 Creates matrices for given tensor order.

 @param torder the tensor order.
 @param bra_basis the molecular basis on bra side.
 @param ket_basis the molecular basis on ket side.
 @return the matrices.
 */
auto makeMatrices(const int64_t torder, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrices;

/**
 Creates matrices for given vector of atoms.

 @param atoms the vector of atoms.
 @param basis the molecular basis.
 @param mtype the matrix type.
 @return the matrices.
 */
auto makeMatrices(const std::vector<int64_t>& atoms, const CMolecularBasis& basis, const mat_t mtype) -> CMatrices;

/**
 Creates matrices for given tensor order.

 @param atoms the vector of atoms.
 @param bra_basis the molecular basis on bra side.
 @param ket_basis the molecular basis on ket side.
 @return the matrices.
 */
auto makeMatrices(const std::vector<int64_t>& atoms, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrices;

/**
 Creates matrices for given vector of atoms.

 @param atoms the vector of atoms.
 @param torder the tensor order.
 @param basis the molecular basis.
 @param mtype the matrix type.
 @return the matrices.
 */
auto makeMatrices(const std::vector<int64_t>& atoms, const int64_t torder, const CMolecularBasis& basis, const mat_t mtype) -> CMatrices;

/**
 Creates matrices for given tensor order.

 @param atoms the vector of atoms.
 @param torder the tensor order.
 @param bra_basis the molecular basis on bra side.
 @param ket_basis the molecular basis on ket side.
 @return the matrices.
 */
auto makeMatrices(const std::vector<int64_t>& atoms, const int64_t torder, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis)
    -> CMatrices;

}  // namespace matfunc

#endif /* MatricesFunc_hpp */
