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



#include "SimdThreeCenterElectronRepulsionDriver.hpp"

#include <string>

#include "ErrorHandler.hpp"
#include "ScreeningFunc.hpp"

auto
CSimdThreeCenterElectronRepulsionDriver::compute(const CMolecule       &molecule,
                                                 const CMolecularBasis &basis,
                                                 const CMolecularBasis &aux_basis,
                                                 const double           threshold) const -> CSparseTensor
{
    // NOTE: the tensor is symmetric in the two sides of the atom pair, so only the
    // upper triangle of the atom basis pair groups is described.

    const auto pattern = make_pattern(molecule, basis, aux_basis, threshold);

    auto tensor = CSparseTensor(pattern);

    tensor.allocate();

    auto distributor = CSimdT3CDistributor<CSparseTensor>(&tensor);

    compute(pattern, molecule, basis, aux_basis, distributor);

    return tensor;
}

auto
CSimdThreeCenterElectronRepulsionDriver::compute(const CMolecule        &molecule,
                                                 const CMolecularBasis  &basis,
                                                 const CMolecularBasis  &aux_basis,
                                                 const double            threshold,
                                                 const std::vector<int> &atoms) const -> CSparseTensor
{
    const auto pattern = make_pattern(molecule, basis, aux_basis, threshold, atoms);

    auto tensor = CSparseTensor(pattern);

    tensor.allocate();

    auto distributor = CSimdT3CDistributor<CSparseTensor>(&tensor);

    compute(pattern, molecule, basis, aux_basis, distributor);

    return tensor;
}

auto
CSimdThreeCenterElectronRepulsionDriver::_make_c_coordinates(const CAtomBasisTripleSparsity &block, const CMolecule &molecule)
    -> CSimdMatrix
{
    const auto &coords = molecule.coordinates();

    const auto &atoms = block.c_atoms();

    const auto natoms = atoms.size();

    auto matrix = CSimdMatrix(3, natoms);

    // NOTE: the coordinates of an axis are stored contiguously over the atoms, as the
    // coordinates of the atom pairs are, so a kernel loads them the same way.

    auto *c_x = matrix.data(0);
    auto *c_y = matrix.data(1);
    auto *c_z = matrix.data(2);

    const auto *rxyz = coords.data();

    for (size_t i = 0; i < natoms; i++)
    {
        errors::assertMsgCritical(atoms[i] < static_cast<int>(coords.size()),
                                  std::string("SimdThreeCenterElectronRepulsionDriver: Atomic index out of range of molecule"));

        const auto r_c = rxyz[atoms[i]].coordinates();

        c_x[i] = r_c[0];

        c_y[i] = r_c[1];

        c_z[i] = r_c[2];
    }

    return matrix;
}
