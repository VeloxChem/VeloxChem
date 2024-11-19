//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef AOIndices_hpp
#define AOIndices_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace aoindices {  // aoindices namespace

/**
 Gets AO indices of the two molecules in a molecular dimer.

 @param mol_1 the first molecule.
 @param mol_2 the second molecule.
 @param basis_1 the basis set for the first molecule.
 @param basis_2 the basis set for the first molecule.
 @return a vector of vector containing the AO indices of the two molecules.
 */
std::vector<std::vector<int>> getDimerAOIndices(const CMolecule&       mol_1,
                                                const CMolecule&       mol_2,
                                                const CMolecularBasis& basis_1,
                                                const CMolecularBasis& basis_2);

/**
 Computes AO-to-atom mapping.

 @param ao_to_atom_ids the vector for storing the mapping.
 @param molecule the molecule.
 @param basis the molecular basis.
 */
void computeAOtoAtomMapping(std::vector<int>& ao_to_atom_ids, const CMolecule& molecule, const CMolecularBasis& basis);

}  // namespace aoindices

#endif /* AOIndices_hpp */
