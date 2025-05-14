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

void
computeAOtoAtomMapping(std::vector<int64_t>& ao_to_atom_ids, const CMolecule& molecule, const CMolecularBasis& basis)
{
    auto natoms = molecule.getNumberOfAtoms();

    auto max_angl = basis.getMaxAngularMomentum();

    // azimuthal quantum number: s,p,d,f,...

    for (int64_t angl = 0, aoidx = 0; angl <= max_angl; angl++)
    {
        auto nsph = angl * 2 + 1;

        // magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...

        for (int64_t isph = 0; isph < nsph; isph++)
        {
            // atoms

            for (int64_t atomidx = 0; atomidx < natoms; atomidx++)
            {
                auto nao = basis.getNumberOfBasisFunctions(std::vector<int64_t>({atomidx}), angl);

                // atomic orbitals

                for (int64_t iao = 0; iao < nao; iao++, aoidx++)
                {
                    ao_to_atom_ids[aoidx] = atomidx;
                }
            }
        }
    }
}

}  // namespace aoindices
