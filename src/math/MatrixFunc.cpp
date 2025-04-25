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

    // add submatrices for 3rd-order case
    
    if (format::upper_case(label) == "3RD_ORDER")
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

        matrix.add({0, 0, nrows, ncols}, {3, 0});
        matrix.add({0, 0, nrows, ncols}, {3, 1});
        matrix.add({0, 0, nrows, ncols}, {3, 2});
        matrix.add({0, 0, nrows, ncols}, {3, 3});
        matrix.add({0, 0, nrows, ncols}, {3, 4});
        matrix.add({0, 0, nrows, ncols}, {3, 5});
        matrix.add({0, 0, nrows, ncols}, {3, 6});
        matrix.add({0, 0, nrows, ncols}, {3, 7});
        matrix.add({0, 0, nrows, ncols}, {3, 8});
        matrix.add({0, 0, nrows, ncols}, {3, 9});
    }
    
    return matrix;
}

}  // namespace matfunc
