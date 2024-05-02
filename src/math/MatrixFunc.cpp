//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
