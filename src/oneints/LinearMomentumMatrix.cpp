//
//                           VELOXCHEM 1.0-RC
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

#include "LinearMomentumMatrix.hpp"

#include <cmath>

CLinearMomentumMatrix::CLinearMomentumMatrix()
{
}

CLinearMomentumMatrix::CLinearMomentumMatrix(const CDenseMatrix& xMatrix, const CDenseMatrix& yMatrix, const CDenseMatrix& zMatrix)

    : _xMatrix(xMatrix)

    , _yMatrix(yMatrix)

    , _zMatrix(zMatrix)
{
}

CLinearMomentumMatrix::CLinearMomentumMatrix(const CLinearMomentumMatrix& source)

    : _xMatrix(source._xMatrix)

    , _yMatrix(source._yMatrix)

    , _zMatrix(source._zMatrix)
{
}

CLinearMomentumMatrix::CLinearMomentumMatrix(CLinearMomentumMatrix&& source) noexcept

    : _xMatrix(std::move(source._xMatrix))

    , _yMatrix(std::move(source._yMatrix))

    , _zMatrix(std::move(source._zMatrix))
{
}

CLinearMomentumMatrix::~CLinearMomentumMatrix()
{
}

CLinearMomentumMatrix&
CLinearMomentumMatrix::operator=(const CLinearMomentumMatrix& source)
{
    if (this == &source) return *this;

    _xMatrix = source._xMatrix;

    _yMatrix = source._yMatrix;

    _zMatrix = source._zMatrix;

    return *this;
}

CLinearMomentumMatrix&
CLinearMomentumMatrix::operator=(CLinearMomentumMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _xMatrix = std::move(source._xMatrix);

    _yMatrix = std::move(source._yMatrix);

    _zMatrix = std::move(source._zMatrix);

    return *this;
}

bool
CLinearMomentumMatrix::operator==(const CLinearMomentumMatrix& other) const
{
    if (_xMatrix != other._xMatrix) return false;

    if (_yMatrix != other._yMatrix) return false;

    if (_zMatrix != other._zMatrix) return false;

    return true;
}

bool
CLinearMomentumMatrix::operator!=(const CLinearMomentumMatrix& other) const
{
    return !(*this == other);
}

std::string
CLinearMomentumMatrix::getStringForComponentX() const
{
    return _xMatrix.getString();
}

std::string
CLinearMomentumMatrix::getStringForComponentY() const
{
    return _yMatrix.getString();
}

std::string
CLinearMomentumMatrix::getStringForComponentZ() const
{
    return _zMatrix.getString();
}

int32_t
CLinearMomentumMatrix::getNumberOfRows() const
{
    return _xMatrix.getNumberOfRows();
}

int32_t
CLinearMomentumMatrix::getNumberOfColumns() const
{
    return _xMatrix.getNumberOfColumns();
}

int32_t
CLinearMomentumMatrix::getNumberOfElements() const
{
    return _xMatrix.getNumberOfElements();
}

const double*
CLinearMomentumMatrix::xvalues() const
{
    return _xMatrix.values();
}

const double*
CLinearMomentumMatrix::yvalues() const
{
    return _yMatrix.values();
}

const double*
CLinearMomentumMatrix::zvalues() const
{
    return _zMatrix.values();
}

std::ostream&
operator<<(std::ostream& output, const CLinearMomentumMatrix& source)
{
    output << std::endl;

    output << "[CLinearMomentumMatrix (Object):" << &source << "]" << std::endl;

    output << "_xMatrix: " << source._xMatrix << std::endl;

    output << "_yMatrix: " << source._yMatrix << std::endl;

    output << "_zMatrix: " << source._zMatrix << std::endl;

    return output;
}
