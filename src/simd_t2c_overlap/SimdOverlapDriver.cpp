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


#include "SimdOverlapDriver.hpp"

#include <string>

#include "ErrorHandler.hpp"
#include "Matrix.hpp"
#include "ScreeningFunc.hpp"

auto
CSimdOverlapDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix
{
    // NOTE: the overlap operator is spherically symmetric about an atom, so the
    // diagonal atom pair blocks hold one value for each pair of basis functions
    // with the same angular momentum.

    auto matrix = CSparseMatrix(molecule, basis, screener::overlap, _threshold, mat_t::symmetric, diagstor::scalar);

    matrix.allocate();

    matrix.zero();

    _compute_pair_blocks(matrix, molecule, basis, basis);

    _compute_diagonal_blocks(matrix, molecule, basis, basis);

    return matrix;
}

auto
CSimdOverlapDriver::compute(const CMolecule &molecule, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) const
    -> CSparseMatrix
{
    auto matrix = CSparseMatrix(molecule, bra_basis, ket_basis, screener::overlap, _threshold, mat_t::general, diagstor::scalar);

    matrix.allocate();

    matrix.zero();

    _compute_pair_blocks(matrix, molecule, bra_basis, ket_basis);

    _compute_diagonal_blocks(matrix, molecule, bra_basis, ket_basis);

    return matrix;
}

auto
CSimdOverlapDriver::_compute_pair_blocks(CSparseMatrix         &matrix,
                                         const CMolecule       &molecule,
                                         const CMolecularBasis &bra_basis,
                                         const CMolecularBasis &ket_basis) const -> void
{
    // TODO: for each off-diagonal block of the matrix
    //         create the coordinates of its atom pairs
    //         for each combination of basis functions on bra and ket sides
    //           determine the number of surviving atom pairs of each pair of
    //             primitives
    //           set up the primitive factors of the pairs of primitives
    //           compute the primitive integrals with the overlap recursion
    //           contract the primitive integrals
    //           transform the contracted integrals to the spherical basis
    //           add the integrals to the values of the block

    errors::assertMsgCritical(false, std::string("SimdOverlapDriver: Computation of off-diagonal blocks is not implemented"));
}

auto
CSimdOverlapDriver::_compute_diagonal_blocks(CSparseMatrix         &matrix,
                                             const CMolecule       &molecule,
                                             const CMolecularBasis &bra_basis,
                                             const CMolecularBasis &ket_basis) const -> void
{
    // TODO: for each diagonal block of the matrix
    //         for each combination of basis functions with the same angular
    //           momentum on bra and ket sides
    //           compute the one-center overlap of the pair of basis functions
    //           add the value to the values of the block

    errors::assertMsgCritical(false, std::string("SimdOverlapDriver: Computation of diagonal blocks is not implemented"));
}
