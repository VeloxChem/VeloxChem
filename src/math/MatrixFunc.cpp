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

#include <algorithm>

#include "CustomViews.hpp"
#include "TensorComponents.hpp"
#include "StringFormat.hpp"

namespace matfunc {  // matfunc namespace

auto
make_matrix(const CMolecularBasis& basis, const mat_t mtype) -> CMatrix
{
    if ((mtype == mat_t::symmetric) || (mtype == mat_t::antisymmetric))
    {
        auto matrix = CMatrix();

        matrix.set_type(mtype);

        std::ranges::for_each(views::triangular(basis.max_angular_momentum() + 1), [&](const auto& index) {
            const auto [i, j]       = index;
            const size_t bra_off    = basis.dimensions_of_basis(i);
            const size_t bra_ncgtos = basis.number_of_basis_functions(i) * tensor::number_of_spherical_components(std::array<int, 1>{i});
            const size_t ket_off    = basis.dimensions_of_basis(j);
            const size_t ket_ncgtos = basis.number_of_basis_functions(j) * tensor::number_of_spherical_components(std::array<int, 1>{j});
            matrix.add({bra_off, ket_off, bra_ncgtos, ket_ncgtos}, {i, j});
        });

        return matrix;
    }
    else
    {
        return matfunc::make_matrix(basis, basis);
    }
}

auto
make_matrix(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrix
{
    auto matrix = CMatrix();

    matrix.set_type(mat_t::general);

    std::ranges::for_each(views::rectangular(bra_basis.max_angular_momentum() + 1, ket_basis.max_angular_momentum() + 1), [&](const auto& index) {
        const auto [i, j]       = index;
        const size_t bra_off    = bra_basis.dimensions_of_basis(i);
        const size_t bra_ncgtos = bra_basis.number_of_basis_functions(i) * tensor::number_of_spherical_components(std::array<int, 1>{i});
        const size_t ket_off    = ket_basis.dimensions_of_basis(j);
        const size_t ket_ncgtos = ket_basis.number_of_basis_functions(j) * tensor::number_of_spherical_components(std::array<int, 1>{j});
        matrix.add({bra_off, ket_off, bra_ncgtos, ket_ncgtos}, {i, j});
    });

    return matrix;
}

auto
make_matrix(const std::string& label, const size_t nrows, const size_t ncols) -> CMatrix
{
    // set up matrix

    auto matrix = CMatrix();

    matrix.set_type(mat_t::general);
    
    // add submatrices for LDA case
    
    if (format::upper_case(label) == "LDA")
    {
        matrix.add({0, 0, nrows, ncols}, {0, 0});
    }
    
    // add submatrices for GGA case
    
    if (format::upper_case(label) == "GGA")
    {
        matrix.add({0, 0, nrows, ncols}, {0, 0});
        
        matrix.add({0, 0, nrows, ncols}, {1, 0});
        
        matrix.add({0, 0, nrows, ncols}, {1, 1});
        
        matrix.add({0, 0, nrows, ncols}, {1, 2});
    }
    
    // add submatrices for GGA case
    
    if (format::upper_case(label) == "MGGA")
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
