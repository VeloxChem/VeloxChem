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

#include "OverlapDriver.hpp"

#include <cstddef>
#include <vector>

#include "AtomBasis.hpp"
#include "BasisFunction.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisOutline.hpp"
#include "Molecule.hpp"
#include "Point.hpp"

namespace newints {

namespace {

/// @brief STUB: screening test for an off-diagonal (two-center) shell pair.
/// Returns true if the shell-pair block must be computed. Its form will be
/// defined in a later step; for now nothing is screened out.
auto
overlap_screener(const CBasisFunction &bra,
                 const CBasisFunction &ket,
                 const TPoint<double> &bra_center,
                 const TPoint<double> &ket_center,
                 const double          threshold) -> bool
{
    // TODO: implement screening based on exponents, separation and threshold.
    return true;
}

/// @brief STUB: same-center overlap block for two shells of equal angular momentum.
/// Returns a zero block of the correct dimensions until the kernel is implemented.
auto
overlap_kernel_diagonal(const CBasisFunction &bra, const CBasisFunction &ket) -> Block
{
    const std::size_t nrows = 2 * bra.get_angular_momentum() + 1;

    const std::size_t ncols = 2 * ket.get_angular_momentum() + 1;

    // TODO: evaluate same-center overlap integrals for this shell pair.
    return Block{nrows, ncols, std::vector<double>(nrows * ncols, 0.0)};
}

/// @brief STUB: two-center overlap block for two shells on different atoms.
/// Returns a zero block of the correct dimensions until the kernel is implemented.
auto
overlap_kernel(const CBasisFunction &bra,
               const CBasisFunction &ket,
               const TPoint<double> &bra_center,
               const TPoint<double> &ket_center) -> Block
{
    const std::size_t nrows = 2 * bra.get_angular_momentum() + 1;

    const std::size_t ncols = 2 * ket.get_angular_momentum() + 1;

    // TODO: evaluate two-center overlap integrals for this shell pair.
    return Block{nrows, ncols, std::vector<double>(nrows * ncols, 0.0)};
}

}  // namespace

auto
OverlapDriver::compute(const CMolecule &molecule, const CMolecularBasis &basis, const double threshold) const -> SparseMatrix
{
    // the two-center overlap matrix is symmetric
    SparseMatrix matrix(SymmetryType::symmetric);

    const auto outline = MolecularBasisOutline(basis);

    const auto coords = molecule.coordinates("au");

    const auto basis_indices = basis.basis_sets_indices();

    const auto atom_bases = basis.basis_sets();

    const auto natoms = static_cast<int>(outline.number_of_atoms());

    // triangular loop over atom pairs with bra atom a <= ket atom b
    for (int a = 0; a < natoms; a++)
    {
        const auto bra_shells = atom_bases[basis_indices[a]].basis_functions();

        const auto bra_idx = outline.basis_function_indices(a);

        for (int b = a; b < natoms; b++)
        {
            if (a == b)
            {
                // same atom: only l == l' shell pairs contribute; triangular over shells
                for (std::size_t p = 0; p < bra_shells.size(); p++)
                {
                    for (std::size_t q = p; q < bra_shells.size(); q++)
                    {
                        if (bra_shells[p].get_angular_momentum() != bra_shells[q].get_angular_momentum()) continue;

                        matrix.add(bra_idx[p], bra_idx[q], overlap_kernel_diagonal(bra_shells[p], bra_shells[q]));
                    }
                }
            }
            else
            {
                // different atoms: all (l, l') shell pairs, screened
                const auto ket_shells = atom_bases[basis_indices[b]].basis_functions();

                const auto ket_idx = outline.basis_function_indices(b);

                for (std::size_t p = 0; p < bra_shells.size(); p++)
                {
                    for (std::size_t q = 0; q < ket_shells.size(); q++)
                    {
                        if (!overlap_screener(bra_shells[p], ket_shells[q], coords[a], coords[b], threshold)) continue;

                        matrix.add(bra_idx[p], ket_idx[q], overlap_kernel(bra_shells[p], ket_shells[q], coords[a], coords[b]));
                    }
                }
            }
        }
    }

    return matrix;
}

}  // namespace newints
