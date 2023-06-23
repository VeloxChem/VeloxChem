#include "MatrixFunc.hpp"

#include "AngularMomentum.hpp"

namespace matfunc {  // matfunc namespace

auto
makeMatrix(const CMolecularBasis& basis, const mat_t mtype) -> CMatrix
{
    // set up matrix

    auto matrix = CMatrix();

    matrix.setType(mtype);

    // add submatrices

    const auto mang = basis.getMaxAngularMomentum();

    for (int64_t i = 0; i <= mang; i++)
    {
        const auto bra_ncgtos = basis.getNumberOfBasisFunctions(i)

                                * angmom::to_SphericalComponents(i);

        const auto bra_off = basis.getDimensionsOfBasis(i);

        const auto joff = ((mtype == mat_t::symm) || (mtype == mat_t::antisymm)) ? i : 0;

        for (int64_t j = joff; j <= mang; j++)
        {
            const auto ket_ncgtos = basis.getNumberOfBasisFunctions(j)

                                    * angmom::to_SphericalComponents(j);

            const auto ket_off = basis.getDimensionsOfBasis(j);

            matrix.add(T4Index({bra_off, ket_off, bra_ncgtos, ket_ncgtos}), {i, j});
        }
    }

    return matrix;
}

auto
makeMatrix(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrix
{
    // set up matrix

    auto matrix = CMatrix();

    matrix.setType(mat_t::gen);

    // add submatrices

    const auto bra_ang = bra_basis.getMaxAngularMomentum();

    const auto ket_ang = ket_basis.getMaxAngularMomentum();

    for (int64_t i = 0; i <= bra_ang; i++)
    {
        const auto bra_ncgtos = bra_basis.getNumberOfBasisFunctions(i)

                                * angmom::to_SphericalComponents(i);

        const auto bra_off = bra_basis.getDimensionsOfBasis(i);

        for (int64_t j = 0; j <= ket_ang; j++)
        {
            const auto ket_ncgtos = ket_basis.getNumberOfBasisFunctions(j)

                                    * angmom::to_SphericalComponents(j);

            const auto ket_off = ket_basis.getDimensionsOfBasis(j);

            matrix.add(T4Index({bra_off, ket_off, bra_ncgtos, ket_ncgtos}), {i, j});
        }
    }

    return matrix;
}

}  // namespace matfunc
