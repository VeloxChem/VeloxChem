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

#include "GtoTransform.hpp"

#include "AngularMomentum.hpp"

namespace gtotra {  // gtotra namespace

CDenseMatrix
to_veloxchem(const CDenseMatrix& matrix, const CMolecularBasis& basis, const CMolecule& molecule)
{
    auto vlxmat = matrix;

    // set up pointers to matrices data

    auto vlxdat = vlxmat.values();

    auto daldat = matrix.values();

    // set up ordering matrix

    auto idsmap = basis.getIndexMapForDalton(molecule);

    auto vlxidx = idsmap.data(0);

    auto dalidx = idsmap.data(1);

    // reorder square matrix

    auto ndim = matrix.getNumberOfColumns();

    for (int32_t i = 0; i < ndim; i++)
    {
        for (int32_t j = 0; j < ndim; j++)
        {
            vlxdat[vlxidx[i] * ndim + vlxidx[j]] = daldat[dalidx[i] * ndim + dalidx[j]];
        }
    }

    return vlxmat;
}

CDenseMatrix
to_dalton(const CDenseMatrix& matrix, const CMolecularBasis& basis, const CMolecule& molecule)
{
    auto dalmat = matrix;

    // set up pointers to matrices data

    auto vlxdat = matrix.values();

    auto daldat = dalmat.values();

    // set up ordering matrix

    auto idsmap = basis.getIndexMapForDalton(molecule);

    auto vlxidx = idsmap.data(0);

    auto dalidx = idsmap.data(1);

    // reorder square matrix

    auto ndim = matrix.getNumberOfColumns();

    for (int32_t i = 0; i < ndim; i++)
    {
        for (int32_t j = 0; j < ndim; j++)
        {
            daldat[dalidx[i] * ndim + dalidx[j]] = vlxdat[vlxidx[i] * ndim + vlxidx[j]];
        }
    }

    return dalmat;
}

std::vector<std::vector<int32_t>>
getBasisFunctionIndicesForAtom(const CMolecule& molecule, const CMolecularBasis& basis, const int32_t atomIdx)
{
    auto idsmap = basis.getIndexMapForDalton(molecule);

    auto vlxidx = idsmap.data(0);

    auto idselm = molecule.getIdsElemental();

    int32_t basisFunctionStart = 0;

    int32_t basisFunctionEnd = 0;

    std::vector<int32_t> angularMomentum;

    for (int32_t i = 0; i <= atomIdx; i++)
    {
        auto atmbas = basis.getAtomBasis(idselm[i]);

        for (int32_t j = 0; j <= atmbas.getMaxAngularMomentum(); j++)
        {
            for (int32_t k = 0; k < atmbas.getNumberOfBasisFunctions(j); k++)
            {
                if (i < atomIdx) basisFunctionStart += angmom::to_SphericalComponents(j);

                basisFunctionEnd += angmom::to_SphericalComponents(j);

                for (int32_t l = 0; l < angmom::to_SphericalComponents(j); l++)
                {
                    angularMomentum.push_back(j);
                }
            }
        }
    }

    std::vector<int32_t> basisFunctionIndices(basisFunctionEnd - basisFunctionStart);

    std::vector<int32_t> basisFunctionAngularMomentum(basisFunctionEnd - basisFunctionStart);

    for (int32_t j = basisFunctionStart; j < basisFunctionEnd; j++)
    {
        basisFunctionIndices[j - basisFunctionStart] = vlxidx[j];

        basisFunctionAngularMomentum[j - basisFunctionStart] = angularMomentum[j];
    }

    return std::vector<std::vector<int32_t>>({basisFunctionIndices, basisFunctionAngularMomentum});
}

}  // namespace gtotra
