#include "MatricesFunc.hpp"

#include <string>
#include <vector>

#include "MatrixFunc.hpp"
#include "StringFormat.hpp"

namespace matfunc {  // matfunc namespace

auto
makeMatrices(const int64_t torder, const CMolecularBasis& basis, const mat_t mtype) -> CMatrices
{
    CMatrices matrices;

    for (const auto& tlabel : fstr::to_TensorComponents(torder))
    {
        matrices.add(matfunc::makeMatrix(basis, mtype), tlabel);
    }

    return matrices;
}

auto
makeMatrices(const int64_t torder, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrices
{
    CMatrices matrices;

    for (const auto& tlabel : fstr::to_TensorComponents(torder))
    {
        matrices.add(matfunc::makeMatrix(bra_basis, ket_basis), tlabel);
    }

    return matrices;
}

auto
makeMatrices(const std::vector<int64_t>& atoms, const CMolecularBasis& basis, const mat_t mtype) -> CMatrices
{
    CMatrices matrices;

    for (const auto atom : atoms)
    {
        matrices.add(matfunc::makeMatrix(basis, mtype), atom);
    }

    return matrices;
}

auto
makeMatrices(const std::vector<int64_t>& atoms, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrices
{
    CMatrices matrices;

    for (const auto atom : atoms)
    {
        matrices.add(matfunc::makeMatrix(bra_basis, ket_basis), atom);
    }

    return matrices;
}

auto
makeMatrices(const std::vector<int64_t>& atoms, const int64_t torder, const CMolecularBasis& basis, const mat_t mtype) -> CMatrices
{
    CMatrices matrices;

    for (const auto atom : atoms)
    {
        for (const auto& tlabel : fstr::to_TensorComponents(torder))
        {
            matrices.add(matfunc::makeMatrix(basis, mtype), atom, tlabel);
        }
    }

    return matrices;
}

auto
makeMatrices(const std::vector<int64_t>& atoms, const int64_t torder, const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrices
{
    CMatrices matrices;

    for (const auto atom : atoms)
    {
        for (const auto& tlabel : fstr::to_TensorComponents(torder))
        {
            matrices.add(matfunc::makeMatrix(bra_basis, ket_basis), atom, tlabel);
        }
    }

    return matrices;
}

}  // namespace matfunc
