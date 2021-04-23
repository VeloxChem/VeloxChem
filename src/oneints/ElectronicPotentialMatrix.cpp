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

#include "ElectronicPotentialMatrix.hpp"

CElectronicPotentialMatrix::CElectronicPotentialMatrix()
{
}

CElectronicPotentialMatrix::CElectronicPotentialMatrix(const CDenseMatrix& matrix)

    : _matrix(matrix)
{
}

CElectronicPotentialMatrix::CElectronicPotentialMatrix(const CElectronicPotentialMatrix& source)

    : _matrix(source._matrix)
{
}

CElectronicPotentialMatrix::CElectronicPotentialMatrix(CElectronicPotentialMatrix&& source) noexcept

    : _matrix(std::move(source._matrix))
{
}

CElectronicPotentialMatrix::~CElectronicPotentialMatrix()
{
}

CElectronicPotentialMatrix&
CElectronicPotentialMatrix::operator=(const CElectronicPotentialMatrix& source)
{
    if (this == &source) return *this;

    _matrix = source._matrix;

    return *this;
}

CElectronicPotentialMatrix&
CElectronicPotentialMatrix::operator=(CElectronicPotentialMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _matrix = std::move(source._matrix);

    return *this;
}

bool
CElectronicPotentialMatrix::operator==(const CElectronicPotentialMatrix& other) const
{
    if (_matrix != other._matrix) return false;

    return true;
}

bool
CElectronicPotentialMatrix::operator!=(const CElectronicPotentialMatrix& other) const
{
    return !(*this == other);
}

std::ostream&
operator<<(std::ostream& output, const CElectronicPotentialMatrix& source)
{
    output << std::endl;

    output << "[CElectronicPotentialMatrix (Object):" << &source << "]" << std::endl;

    output << "_matrix: " << source._matrix << std::endl;

    return output;
}
