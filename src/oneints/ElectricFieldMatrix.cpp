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

#include "ElectricFieldMatrix.hpp"

#include <array>

CElectricFieldMatrix::CElectricFieldMatrix()
{
}

CElectricFieldMatrix::CElectricFieldMatrix(const CDenseMatrix& xMatrix, const CDenseMatrix& yMatrix, const CDenseMatrix& zMatrix)

    : _xMatrix(xMatrix)

    , _yMatrix(yMatrix)

    , _zMatrix(zMatrix)
{
}

CElectricFieldMatrix::CElectricFieldMatrix(const std::array<CDenseMatrix, 3>& matrices)

    : _xMatrix(matrices[0])

    , _yMatrix(matrices[1])

    , _zMatrix(matrices[2])
{
}

CElectricFieldMatrix::CElectricFieldMatrix(const CElectricFieldMatrix& source)

    : _xMatrix(source._xMatrix)

    , _yMatrix(source._yMatrix)

    , _zMatrix(source._zMatrix)
{
}

CElectricFieldMatrix::CElectricFieldMatrix(CElectricFieldMatrix&& source) noexcept

    : _xMatrix(std::move(source._xMatrix))

    , _yMatrix(std::move(source._yMatrix))

    , _zMatrix(std::move(source._zMatrix))
{
}

CElectricFieldMatrix::~CElectricFieldMatrix()
{
}

CElectricFieldMatrix&
CElectricFieldMatrix::operator=(const CElectricFieldMatrix& source)
{
    if (this == &source) return *this;

    _xMatrix = source._xMatrix;

    _yMatrix = source._yMatrix;

    _zMatrix = source._zMatrix;

    return *this;
}

CElectricFieldMatrix&
CElectricFieldMatrix::operator=(CElectricFieldMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _xMatrix = std::move(source._xMatrix);

    _yMatrix = std::move(source._yMatrix);

    _zMatrix = std::move(source._zMatrix);

    return *this;
}

bool
CElectricFieldMatrix::operator==(const CElectricFieldMatrix& other) const
{
    if (_xMatrix != other._xMatrix) return false;

    if (_yMatrix != other._yMatrix) return false;

    if (_zMatrix != other._zMatrix) return false;

    return true;
}

bool
CElectricFieldMatrix::operator!=(const CElectricFieldMatrix& other) const
{
    return !(*this == other);
}

std::string
CElectricFieldMatrix::getStringForComponentX() const
{
    return _xMatrix.getString();
}

std::string
CElectricFieldMatrix::getStringForComponentY() const
{
    return _yMatrix.getString();
}

std::string
CElectricFieldMatrix::getStringForComponentZ() const
{
    return _zMatrix.getString();
}

std::string
CElectricFieldMatrix::getString() const
{
    return getStringForComponentX() + getStringForComponentY() + getStringForComponentZ();
}

int32_t
CElectricFieldMatrix::getNumberOfRows() const
{
    return _xMatrix.getNumberOfRows();
}

int32_t
CElectricFieldMatrix::getNumberOfColumns() const
{
    return _xMatrix.getNumberOfColumns();
}

int32_t
CElectricFieldMatrix::getNumberOfElements() const
{
    return _xMatrix.getNumberOfElements();
}

const double*
CElectricFieldMatrix::xvalues() const
{
    return _xMatrix.values();
}

const double*
CElectricFieldMatrix::yvalues() const
{
    return _yMatrix.values();
}

const double*
CElectricFieldMatrix::zvalues() const
{
    return _zMatrix.values();
}

const double*
CElectricFieldMatrix::values(cartesians cart) const
{
    const double* retval = nullptr;
    switch (cart)
    {
        case cartesians::X:
            retval = _xMatrix.values();
            break;
        case cartesians::Y:
            retval = _yMatrix.values();
            break;
        case cartesians::Z:
            retval = _zMatrix.values();
            break;
    }
    return retval;
}

const double*
CElectricFieldMatrix::values(int32_t cart) const
{
    auto d = static_cast<cartesians>(cart);
    return values(d);
}

std::string
CElectricFieldMatrix::repr() const
{
    std::ostringstream os;

    os << std::endl;

    os << "[CElectricFieldMatrix (Object):" << this << "]" << std::endl;

    os << "_xMatrix: " << _xMatrix << std::endl;

    os << "_yMatrix: " << _yMatrix << std::endl;

    os << "_zMatrix: " << _zMatrix << std::endl;

    return os.str();
}

std::ostream&
operator<<(std::ostream& output, const CElectricFieldMatrix& source)
{
    return (output << source.repr());
}
