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


#include "SimdNuclearPotentialDriver.hpp"

auto
CSimdNuclearPotentialDriver::nuclei_of(const CMolecule &molecule) -> std::pair<std::vector<double>, std::vector<double>>
{
    const auto identifiers = molecule.identifiers();

    const auto &coords = molecule.coordinates("au");

    std::vector<double> charges;

    std::vector<double> points;

    charges.reserve(identifiers.size());

    points.reserve(3 * identifiers.size());

    for (size_t i = 0; i < identifiers.size(); i++)
    {
        charges.push_back(static_cast<double>(identifiers[i]));

        const auto r = coords[i].coordinates();

        points.push_back(r[0]);
        points.push_back(r[1]);
        points.push_back(r[2]);
    }

    return {std::move(charges), std::move(points)};
}

auto
CSimdNuclearPotentialDriver::compute_matrix(const CSparsityPattern    &pattern,
                                            const CMolecule           &molecule,
                                            const CMolecularBasis     &bra_basis,
                                            const CMolecularBasis     &ket_basis,
                                            const std::vector<double> &charges,
                                            const std::vector<double> &points) const -> CSparseMatrix
{
    // NOTE: the values blocks are not set to zero after they are allocated, as every
    // value of every block is written below: a kernel writes the integrals of the
    // atom pairs it reaches and zeros of the remaining ones.

    auto matrix = CSparseMatrix(pattern);

    matrix.allocate();

    auto distributor = CSimdT2CDistributor<CSparseMatrix>(&matrix);

    compute(pattern, molecule, bra_basis, ket_basis, charges, points, distributor);

    return matrix;
}

auto
CSimdNuclearPotentialDriver::compute_matrix(const CMolecule           &molecule,
                                            const CMolecularBasis     &basis,
                                            const std::vector<double> &charges,
                                            const std::vector<double> &points) const -> CSparseMatrix
{
    return compute_matrix(make_pattern(molecule, basis, basis, mat_t::symmetric), molecule, basis, basis, charges, points);
}

auto
CSimdNuclearPotentialDriver::compute_matrix(const CMolecule &molecule, const CMolecularBasis &basis) const -> CSparseMatrix
{
    const auto [charges, points] = nuclei_of(molecule);

    return compute_matrix(molecule, basis, charges, points);
}
