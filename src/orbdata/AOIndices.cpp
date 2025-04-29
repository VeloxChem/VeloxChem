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

#include "AOIndices.hpp"

#include <algorithm>

namespace aoindices {  // aoindices namespace

std::vector<std::vector<int>>
getDimerAOIndices(const CMolecule& mol_1, const CMolecule& mol_2, const CMolecularBasis& basis_1, const CMolecularBasis& basis_2)
{
    // get indices of atomic orbitals located on each molecule

    // angular momentum based AO ordering
    // S, P-1, P0, P+1, D-2, D-1, D0, D+1, D+2, ...

    // AO_type:  S S S S S S S S P-1 P-1 P0 P0 P+1 P+1 ...
    // Molecule: A A A A B B B B A   B   A  B  A   B   ...

    int max_angl_1 = basis_1.max_angular_momentum();

    int max_angl_2 = basis_2.max_angular_momentum();

    int max_angl = std::max(max_angl_1, max_angl_2);

    std::vector<std::vector<int>> aoinds_mol;

    aoinds_mol.push_back(std::vector<int>());

    aoinds_mol.push_back(std::vector<int>());

    for (int aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        int nao_1 = basis_1.number_of_basis_functions(angl);

        int nao_2 = basis_2.number_of_basis_functions(angl);

        for (int s = -angl; s <= angl; s++)
        {
            for (int i = 0; i < nao_1; i++, aoidx++)
            {
                aoinds_mol[0].push_back(aoidx);
            }

            for (int i = 0; i < nao_2; i++, aoidx++)
            {
                aoinds_mol[1].push_back(aoidx);
            }
        }
    }

    return aoinds_mol;
}

void
computeAOtoAtomMapping(std::vector<int>& ao_to_atom_ids, const CMolecule& molecule, const CMolecularBasis& basis)
{
    auto natoms = molecule.number_of_atoms();

    auto max_angl = basis.max_angular_momentum();

    // azimuthal quantum number: s,p,d,f,...

    for (int angl = 0, aoidx = 0; angl <= max_angl; angl++)
    {
        auto nsph = angl * 2 + 1;

        // magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...

        for (int isph = 0; isph < nsph; isph++)
        {
            // atoms

            for (int atomidx = 0; atomidx < natoms; atomidx++)
            {
                auto nao = basis.number_of_basis_functions(std::vector<int>({atomidx}), angl);

                // atomic orbitals

                for (int iao = 0; iao < nao; iao++, aoidx++)
                {
                    ao_to_atom_ids[aoidx] = atomidx;
                }
            }
        }
    }
}

}  // namespace aoindices
