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

auto
CSimdOverlapDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix
{
    // NOTE: the overlap operator is spherically symmetric about an atom, so the
    // diagonal atom pair blocks hold one value for each pair of basis functions
    // with the same angular momentum.

    const auto pattern = make_pattern(molecule, basis, basis, mat_t::symmetric);

    // NOTE: the values blocks are not set to zero after they are allocated, as
    // every value of every block is written below. A combination of basis
    // functions reaching no atom pair holds no values, and a kernel writes the
    // integrals of the atom pairs it reaches and zeros of the remaining ones, so
    // no value is left with the undefined content of the allocation.

    auto matrix = CSparseMatrix(pattern);

    matrix.allocate();

    auto distributor = CSimdT2CDistributor<CSparseMatrix>(&matrix);

    compute(pattern, molecule, basis, basis, distributor);

    return matrix;
}

auto
CSimdOverlapDriver::compute(const CMolecule &molecule, const CMolecularBasis &bra_basis, const CMolecularBasis &ket_basis) const
    -> CSparseMatrix
{
    const auto pattern = make_pattern(molecule, bra_basis, ket_basis, mat_t::general);

    auto matrix = CSparseMatrix(pattern);

    matrix.allocate();

    auto distributor = CSimdT2CDistributor<CSparseMatrix>(&matrix);

    compute(pattern, molecule, bra_basis, ket_basis, distributor);

    return matrix;
}
