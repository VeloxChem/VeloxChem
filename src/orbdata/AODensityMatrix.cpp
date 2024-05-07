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

#include "AODensityMatrix.hpp"

#include "DenseLinearAlgebra.hpp"
#include "ErrorHandler.hpp"
#include "StringFormat.hpp"

CAODensityMatrix::CAODensityMatrix()

    : _denType(denmat::rest)
{
}

CAODensityMatrix::CAODensityMatrix(const std::vector<CDenseMatrix>& denMatrices, const denmat denType)

    : _denMatrices(denMatrices)

    , _denType(denType)
{
    if (denType == denmat::unrest)
    {
        errors::assertMsgCritical(denMatrices.size() % 2 == 0,
                                  "AODensityMatrix: Odd number of matrices for unrestricted or restricted open-shell density");
    }
}

CAODensityMatrix::CAODensityMatrix(const CAODensityMatrix& source)

    : _denMatrices(source._denMatrices)

    , _denType(source._denType)
{
}

CAODensityMatrix::CAODensityMatrix(CAODensityMatrix&& source) noexcept

    : _denMatrices(std::move(source._denMatrices))

    , _denType(std::move(source._denType))
{
}

CAODensityMatrix::~CAODensityMatrix()
{
}

CAODensityMatrix&
CAODensityMatrix::operator=(const CAODensityMatrix& source)
{
    if (this == &source) return *this;

    _denMatrices = source._denMatrices;

    _denType = source._denType;

    return *this;
}

CAODensityMatrix&
CAODensityMatrix::operator=(CAODensityMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _denMatrices = std::move(source._denMatrices);

    _denType = std::move(source._denType);

    return *this;
}

bool
CAODensityMatrix::operator==(const CAODensityMatrix& other) const
{
    if (_denType != other._denType) return false;

    if (_denMatrices.size() != other._denMatrices.size()) return false;

    for (size_t i = 0; i < _denMatrices.size(); i++)
    {
        if (_denMatrices[i] != other._denMatrices[i]) return false;
    }

    return true;
}

bool
CAODensityMatrix::operator!=(const CAODensityMatrix& other) const
{
    return !(*this == other);
}

bool
CAODensityMatrix::isClosedShell() const
{
    return (_denType == denmat::rest);
}

int64_t
CAODensityMatrix::getNumberOfDensityMatrices() const
{
    if (isClosedShell())
    {
        return static_cast<int64_t>(_denMatrices.size());
    }
    else
    {
        return static_cast<int64_t>(_denMatrices.size()) / 2;
    }
}

denmat
CAODensityMatrix::getDensityType() const
{
    return _denType;
}

int64_t
CAODensityMatrix::getNumberOfRows(const int64_t iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].getNumberOfRows();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix].getNumberOfRows();
        }
    }

    return 0;
}

int64_t
CAODensityMatrix::getNumberOfColumns(const int64_t iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].getNumberOfColumns();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix].getNumberOfColumns();
        }
    }

    return 0;
}

int64_t
CAODensityMatrix::getNumberOfElements(const int64_t iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].getNumberOfElements();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix].getNumberOfElements();
        }
    }

    return 0;
}

const double*
CAODensityMatrix::alphaDensity(const int64_t iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].values();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix].values();
        }
    }

    return nullptr;
}

const double*
CAODensityMatrix::betaDensity(const int64_t iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].values();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix + 1].values();
        }
    }

    return nullptr;
}

const CDenseMatrix&
CAODensityMatrix::getReferenceToAlphaDensity(const int64_t iDensityMatrix) const
{
    if (isClosedShell())
    {
        return _denMatrices[iDensityMatrix];
    }
    else
    {
        return _denMatrices[2 * iDensityMatrix];
    }
}

const CDenseMatrix&
CAODensityMatrix::getReferenceToBetaDensity(const int64_t iDensityMatrix) const
{
    if (isClosedShell())
    {
        return _denMatrices[iDensityMatrix];
    }
    else
    {
        return _denMatrices[2 * iDensityMatrix + 1];
    }
}
