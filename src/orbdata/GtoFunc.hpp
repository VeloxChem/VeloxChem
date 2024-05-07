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

#ifndef GtoFunc_hpp
#define GtoFunc_hpp

#include <cstdint>
#include <vector>

#include "GtoBlock.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace gtofunc {  // gtofunc namespace

/**
 Creates vector of contracted GTOs blocks.

 @param basis the molecular basis.
 @param molecule the molecule.
 @return the vector of contracted GTOs blocks.
 */
auto makeGtoBlocks(const CMolecularBasis& basis, const CMolecule& molecule) -> std::vector<CGtoBlock>;

/**
 Creates vector of contracted GTOs blocks.

 @param basis the molecular basis.
 @param molecule the molecule.
 @param atoms the vector of atoms to select contracted GTOs.
 @return the vector of contracted GTOs blocks.
 */
auto makeGtoBlocks(const CMolecularBasis& basis, const CMolecule& molecule, const std::vector<int64_t>& atoms) -> std::vector<CGtoBlock>;

/**
 Gets number of atomic orbitals from vector of contracted GTOs blocks.

 @param gto_blocks the vector of contracted GTOs blocks.
 @return the number of atomic orbitals.
 */
auto getNumberOfAtomicOrbitals(const std::vector<CGtoBlock>& gto_blocks) -> int64_t;

}  // namespace gtofunc

#endif /* GtoFunc_hpp */
