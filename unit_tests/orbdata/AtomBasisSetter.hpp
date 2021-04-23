//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#ifndef AtomBasisSetter_hpp
#define AtomBasisSetter_hpp

#include <cstdint>

#include "AtomBasis.hpp"

namespace vlxbas {  // vlxbas namespace

CAtomBasis getAtomBasisEmpty();

CAtomBasis getAtomBasisForH();

CAtomBasis getAtomBasisForLi();

CAtomBasis getAtomBasisForLiX();

CAtomBasis getNormalizedAtomBasisForH();

CAtomBasis getNormalizedAtomBasisForHe();

CAtomBasis getAtomBasisSPDForHe();

CAtomBasis getNormalizedAtomBasisForO();

CAtomBasis getNormalizedAtomBasisForSe();

CAtomBasis getMinimalBasisForH();

CAtomBasis getMinimalBasisForHe();

CAtomBasis getMinimalBasisForC();

CAtomBasis getMinimalBasisForN();

CAtomBasis getMinimalBasisForO();

CAtomBasis getTestBasisForH();

CAtomBasis getTestBasisForLi();

}  // namespace vlxbas

#endif /* AtomBasisSetter_hpp */
