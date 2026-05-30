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

#include "MolecularBasisOutline.hpp"

#include "AtomBasis.hpp"
#include "BasisFunction.hpp"
#include "MolecularBasis.hpp"

namespace newints {

MolecularBasisOutline::MolecularBasisOutline(const CMolecularBasis &basis)
{
    const auto basis_indices = basis.basis_sets_indices();

    const auto atom_bases = basis.basis_sets();

    const auto natoms = basis_indices.size();

    _atom_shell_offsets.reserve(natoms + 1);

    _atom_shell_offsets.push_back(0);

    int offset = 0;  // running angular-expanded offset = contracted-GTO index

    for (size_t atom = 0; atom < natoms; atom++)
    {
        const auto &atom_basis = atom_bases[basis_indices[atom]];

        for (const auto &bf : atom_basis.basis_functions())
        {
            const auto angmom = bf.get_angular_momentum();

            _indices.push_back(offset);

            _angular_momenta.push_back(angmom);

            _atom_indices.push_back(static_cast<int>(atom));

            offset += 2 * angmom + 1;
        }

        _atom_shell_offsets.push_back(static_cast<int>(_indices.size()));
    }
}

auto
MolecularBasisOutline::operator==(const MolecularBasisOutline &other) const -> bool
{
    return (_indices == other._indices) && (_angular_momenta == other._angular_momenta) && (_atom_indices == other._atom_indices) &&
           (_atom_shell_offsets == other._atom_shell_offsets);
}

auto
MolecularBasisOutline::number_of_atoms() const -> std::size_t
{
    return _atom_shell_offsets.empty() ? 0 : _atom_shell_offsets.size() - 1;
}

auto
MolecularBasisOutline::number_of_basis_functions() const -> std::size_t
{
    return _indices.size();
}

auto
MolecularBasisOutline::number_of_atomic_orbitals() const -> std::size_t
{
    if (_indices.empty()) return 0;

    return static_cast<std::size_t>(_indices.back() + 2 * _angular_momenta.back() + 1);
}

auto
MolecularBasisOutline::basis_function_indices(const int atom) const -> std::vector<int>
{
    const auto start = _atom_shell_offsets[atom];

    const auto end = _atom_shell_offsets[atom + 1];

    return std::vector<int>(_indices.begin() + start, _indices.begin() + end);
}

auto
MolecularBasisOutline::indices() const -> std::vector<int>
{
    return _indices;
}

auto
MolecularBasisOutline::angular_momenta() const -> std::vector<int>
{
    return _angular_momenta;
}

auto
MolecularBasisOutline::atom_indices() const -> std::vector<int>
{
    return _atom_indices;
}

}  // namespace newints
