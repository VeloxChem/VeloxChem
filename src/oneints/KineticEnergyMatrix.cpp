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

#include "KineticEnergyMatrix.hpp"

#include "DenseLinearAlgebra.hpp"

CKineticEnergyMatrix::CKineticEnergyMatrix()
{
}

CKineticEnergyMatrix::CKineticEnergyMatrix(const CDenseMatrix& matrix)

    : _matrix(matrix)
{
}

CKineticEnergyMatrix::CKineticEnergyMatrix(const CKineticEnergyMatrix& source)

    : _matrix(source._matrix)
{
}

CKineticEnergyMatrix::CKineticEnergyMatrix(CKineticEnergyMatrix&& source) noexcept

    : _matrix(std::move(source._matrix))
{
}

CKineticEnergyMatrix::~CKineticEnergyMatrix()
{
}

CKineticEnergyMatrix&
CKineticEnergyMatrix::operator=(const CKineticEnergyMatrix& source)
{
    if (this == &source) return *this;

    _matrix = source._matrix;

    return *this;
}

CKineticEnergyMatrix&
CKineticEnergyMatrix::operator=(CKineticEnergyMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _matrix = std::move(source._matrix);

    return *this;
}

bool
CKineticEnergyMatrix::operator==(const CKineticEnergyMatrix& other) const
{
    if (_matrix != other._matrix) return false;

    return true;
}

bool
CKineticEnergyMatrix::operator!=(const CKineticEnergyMatrix& other) const
{
    return !(*this == other);
}

std::string
CKineticEnergyMatrix::getString() const
{
    return _matrix.getString();
}

int32_t
CKineticEnergyMatrix::getNumberOfRows() const
{
    return _matrix.getNumberOfRows();
}

int32_t
CKineticEnergyMatrix::getNumberOfColumns() const
{
    return _matrix.getNumberOfColumns();
}

int32_t
CKineticEnergyMatrix::getNumberOfElements() const
{
    return _matrix.getNumberOfElements();
}

const double*
CKineticEnergyMatrix::values() const
{
    return _matrix.values();
}

double
CKineticEnergyMatrix::getKineticEnergy(const CAODensityMatrix& aoDensityMatrix, const int32_t iDensityMatrix) const
{
    if (iDensityMatrix < aoDensityMatrix.getNumberOfMatrices())
    {
        if (aoDensityMatrix.getDensityType() == denmat::rest)
        {
            return denblas::trace(_matrix, aoDensityMatrix.getReferenceToDensity(iDensityMatrix));
        }
        else
        {
            auto idensity_a = 2 * iDensityMatrix;

            auto idensity_b = 2 * iDensityMatrix + 1;

            auto e_a = 0.5 * denblas::trace(_matrix, aoDensityMatrix.getReferenceToDensity(idensity_a));

            auto e_b = 0.5 * denblas::trace(_matrix, aoDensityMatrix.getReferenceToDensity(idensity_b));

            return e_a + e_b;
        }
    }

    return 0.0;
}

std::ostream&
operator<<(std::ostream& output, const CKineticEnergyMatrix& source)
{
    output << std::endl;

    output << "[CKineticEnergyMatrix (Object):" << &source << "]" << std::endl;

    output << "_matrix: " << source._matrix << std::endl;

    return output;
}
