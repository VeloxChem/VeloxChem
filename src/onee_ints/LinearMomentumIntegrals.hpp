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

#ifndef LinearMomentumIntegrals_hpp
#define LinearMomentumIntegrals_hpp

#include <vector>

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace onee {  // onee namespace

auto
computeLinearMomentumIntegrals(const CMolecule&           molecule,
                               const CMolecularBasis&     basis) -> std::vector<CDenseMatrix>;

}  // namespace onee

#endif /* LinearMomentumIntegrals_hpp */
