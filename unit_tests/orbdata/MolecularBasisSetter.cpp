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

#include "MolecularBasisSetter.hpp"

#include "AtomBasisSetter.hpp"

namespace vlxbas {  // vlxbas namespace

CMolecularBasis
getMolecularBasisEmpty()
{
    CMolecularBasis mbas;

    mbas.setMaxAngularMomentum(-1);

    return mbas;
}

CMolecularBasis
getMolecularBasisForLiH()
{
    CMolecularBasis mbas;

    mbas.setLabel({"def2-SVP"});

    mbas.addAtomBasis(getAtomBasisForLi());

    mbas.addAtomBasis(getAtomBasisForH());

    return mbas;
}

CMolecularBasis
getMolecularBasisForLiHX()
{
    CMolecularBasis mbas;

    mbas.setLabel({"def2-SVP-X"});

    mbas.addAtomBasis(getAtomBasisForLiX());

    mbas.addAtomBasis(getAtomBasisForH());

    return mbas;
}

CMolecularBasis
getMolecularBasisForHeAtom()
{
    CMolecularBasis mbas;

    mbas.setLabel({"def2-SVP"});

    mbas.addAtomBasis(getNormalizedAtomBasisForHe());

    return mbas;
}

CMolecularBasis
getMolecularBasisSPDForHeAtom()
{
    CMolecularBasis mbas;

    mbas.setLabel({"XTEST-SPD"});

    mbas.addAtomBasis(getAtomBasisSPDForHe());

    return mbas;
}

CMolecularBasis
getMolecularBasisForH2O()
{
    CMolecularBasis mbas;

    mbas.setLabel({"def2-SVP"});

    mbas.addAtomBasis(getNormalizedAtomBasisForO());

    mbas.addAtomBasis(getNormalizedAtomBasisForH());

    return mbas;
}

CMolecularBasis
getMolecularBasisForH2Se()
{
    CMolecularBasis mbas;

    mbas.setLabel({"def2-SVP"});

    mbas.addAtomBasis(getNormalizedAtomBasisForSe());

    mbas.addAtomBasis(getNormalizedAtomBasisForH());

    return mbas;
}

CMolecularBasis
getMinimalBasisForHeAtom()
{
    CMolecularBasis mbas;

    mbas.setLabel({"AO-START-GUESS"});

    mbas.addAtomBasis(getMinimalBasisForHe());

    return mbas;
}

CMolecularBasis
getMinimalBasisForH2O()
{
    CMolecularBasis mbas;

    mbas.setLabel({"AO-START-GUESS"});

    mbas.addAtomBasis(getMinimalBasisForO());

    mbas.addAtomBasis(getMinimalBasisForH());

    return mbas;
}

CMolecularBasis
getMinimalBasisForNH3CH4()
{
    CMolecularBasis mbas;

    mbas.setLabel({"AO-START-GUESS"});

    mbas.addAtomBasis(getMinimalBasisForC());

    mbas.addAtomBasis(getMinimalBasisForN());

    mbas.addAtomBasis(getMinimalBasisForH());

    return mbas;
}

CMolecularBasis
getTestBasisForLiH()
{
    CMolecularBasis mbas;

    mbas.setLabel({"Test-Basis"});

    mbas.addAtomBasis(getTestBasisForLi());

    mbas.addAtomBasis(getTestBasisForH());

    return mbas;
}

CMolecularBasis
getReducedTestBasisForLiH()
{
    CMolecularBasis mbas;

    mbas.setLabel({"Reduced-Basis"});

    mbas.addAtomBasis(getAtomBasisForLi());

    mbas.addAtomBasis(getTestBasisForH());

    return mbas;
}

}  // namespace vlxbas
