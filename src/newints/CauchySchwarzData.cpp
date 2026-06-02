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

#include "CauchySchwarzData.hpp"

#include <cmath>

#include "AtomBasis.hpp"
#include "BasisFunction.hpp"
#include "CoulombDiagonal.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisOutline.hpp"

namespace newints {

CauchySchwarzData::CauchySchwarzData(const CMolecularBasis &basis)
{
    const auto outline = MolecularBasisOutline(basis);

    const auto basis_indices = basis.basis_sets_indices();

    const auto atom_bases = basis.basis_sets();

    const auto natoms = static_cast<int>(outline.number_of_atoms());

    q_.reserve(outline.number_of_basis_functions());

    // atom-major shell order, matching MolecularBasisOutline.indices(): for each atom,
    // its shells in natural angular-momentum order. Q_P = sqrt((P|P)) from the concentric
    // two-center Coulomb self-integral (m-independent, so this single value is exact).
    for (int a = 0; a < natoms; a++)
    {
        const auto shells = atom_bases[basis_indices[a]].basis_functions();

        for (const auto &shell : shells)
        {
            q_.push_back(std::sqrt(coulomb_concentric_value(shell, shell)));
        }
    }
}

}  // namespace newints
