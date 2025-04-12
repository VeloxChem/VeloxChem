//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "MatrixFunc.hpp"

#include "AngularMomentum.hpp"
#include "StringFormat.hpp"

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

auto
makeMatrix(const std::string& label, const int64_t nrows, const int64_t ncols) -> CMatrix
{
    // set up matrix

    auto matrix = CMatrix();

    matrix.setType(mat_t::gen);

    // add submatrices for LDA case

    if (fstr::upcase(label) == "LDA")
    {
        matrix.add({0, 0, nrows, ncols}, {0, 0});
    }

    // add submatrices for GGA case

    if (fstr::upcase(label) == "GGA")
    {
        matrix.add({0, 0, nrows, ncols}, {0, 0});

        matrix.add({0, 0, nrows, ncols}, {1, 0});

        matrix.add({0, 0, nrows, ncols}, {1, 1});

        matrix.add({0, 0, nrows, ncols}, {1, 2});
    }

    // add submatrices for GGA case

    if (fstr::upcase(label) == "MGGA")
    {
        matrix.add({0, 0, nrows, ncols}, {0, 0});

        matrix.add({0, 0, nrows, ncols}, {1, 0});

        matrix.add({0, 0, nrows, ncols}, {1, 1});

        matrix.add({0, 0, nrows, ncols}, {1, 2});

        matrix.add({0, 0, nrows, ncols}, {2, 0});

        matrix.add({0, 0, nrows, ncols}, {2, 1});

        matrix.add({0, 0, nrows, ncols}, {2, 2});

        matrix.add({0, 0, nrows, ncols}, {2, 3});

        matrix.add({0, 0, nrows, ncols}, {2, 4});

        matrix.add({0, 0, nrows, ncols}, {2, 5});
    }

    return matrix;
}

}  // namespace matfunc
