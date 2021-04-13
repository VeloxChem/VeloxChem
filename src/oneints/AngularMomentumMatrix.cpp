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

#include "AngularMomentumMatrix.hpp"

#include <array>
#include <cmath>

CAngularMomentumMatrix::CAngularMomentumMatrix()

    : _xOrigin(0.0)

    , _yOrigin(0.0)

    , _zOrigin(0.0)
{
}

CAngularMomentumMatrix::CAngularMomentumMatrix(const CDenseMatrix& xMatrix,
                                               const CDenseMatrix& yMatrix,
                                               const CDenseMatrix& zMatrix,
                                               const double        xOrigin,
                                               const double        yOrigin,
                                               const double        zOrigin)

    : _xOrigin(xOrigin)

    , _yOrigin(yOrigin)

    , _zOrigin(zOrigin)

    , _xMatrix(xMatrix)

    , _yMatrix(yMatrix)

    , _zMatrix(zMatrix)
{
}

CAngularMomentumMatrix::CAngularMomentumMatrix(const std::array<CDenseMatrix, 3>& matrices, const std::array<double, 3>& origin)
    : _xOrigin(origin[0])

    , _yOrigin(origin[1])

    , _zOrigin(origin[2])

    , _xMatrix(matrices[0])

    , _yMatrix(matrices[1])

    , _zMatrix(matrices[2])
{
}

CAngularMomentumMatrix::CAngularMomentumMatrix(const CAngularMomentumMatrix& source)

    : _xOrigin(source._xOrigin)

    , _yOrigin(source._yOrigin)

    , _zOrigin(source._zOrigin)

    , _xMatrix(source._xMatrix)

    , _yMatrix(source._yMatrix)

    , _zMatrix(source._zMatrix)
{
}

CAngularMomentumMatrix::CAngularMomentumMatrix(CAngularMomentumMatrix&& source) noexcept

    : _xOrigin(std::move(source._xOrigin))

    , _yOrigin(std::move(source._yOrigin))

    , _zOrigin(std::move(source._zOrigin))

    , _xMatrix(std::move(source._xMatrix))

    , _yMatrix(std::move(source._yMatrix))

    , _zMatrix(std::move(source._zMatrix))
{
}

CAngularMomentumMatrix::~CAngularMomentumMatrix()
{
}

CAngularMomentumMatrix&
CAngularMomentumMatrix::operator=(const CAngularMomentumMatrix& source)
{
    if (this == &source) return *this;

    _xMatrix = source._xMatrix;

    _yMatrix = source._yMatrix;

    _zMatrix = source._zMatrix;

    _xOrigin = source._xOrigin;

    _yOrigin = source._yOrigin;

    _zOrigin = source._zOrigin;

    return *this;
}

CAngularMomentumMatrix&
CAngularMomentumMatrix::operator=(CAngularMomentumMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _xMatrix = std::move(source._xMatrix);

    _yMatrix = std::move(source._yMatrix);

    _zMatrix = std::move(source._zMatrix);

    _xOrigin = std::move(source._xOrigin);

    _yOrigin = std::move(source._yOrigin);

    _zOrigin = std::move(source._zOrigin);

    return *this;
}

bool
CAngularMomentumMatrix::operator==(const CAngularMomentumMatrix& other) const
{
    if (_xMatrix != other._xMatrix) return false;

    if (_yMatrix != other._yMatrix) return false;

    if (_zMatrix != other._zMatrix) return false;

    if (std::fabs(_xOrigin - other._xOrigin) > 1.0e-13) return false;

    if (std::fabs(_yOrigin - other._yOrigin) > 1.0e-13) return false;

    if (std::fabs(_zOrigin - other._zOrigin) > 1.0e-13) return false;

    return true;
}

bool
CAngularMomentumMatrix::operator!=(const CAngularMomentumMatrix& other) const
{
    return !(*this == other);
}

void
CAngularMomentumMatrix::setOriginCoordinates(const double xOrigin, const double yOrigin, const double zOrigin)
{
    _xOrigin = xOrigin;

    _yOrigin = yOrigin;

    _zOrigin = zOrigin;
}

void
CAngularMomentumMatrix::setOriginCoordinates(const std::array<double, 3>& origin)
{
    _xOrigin = origin[0];

    _yOrigin = origin[1];

    _zOrigin = origin[2];
}

std::string
CAngularMomentumMatrix::getStringForComponentX() const
{
    return _xMatrix.getString();
}

std::string
CAngularMomentumMatrix::getStringForComponentY() const
{
    return _yMatrix.getString();
}

std::string
CAngularMomentumMatrix::getStringForComponentZ() const
{
    return _zMatrix.getString();
}

std::string
CAngularMomentumMatrix::getString() const
{
    return getStringForComponentX() + getStringForComponentY() + getStringForComponentZ();
}

int32_t
CAngularMomentumMatrix::getNumberOfRows() const
{
    return _xMatrix.getNumberOfRows();
}

int32_t
CAngularMomentumMatrix::getNumberOfColumns() const
{
    return _xMatrix.getNumberOfColumns();
}

int32_t
CAngularMomentumMatrix::getNumberOfElements() const
{
    return _xMatrix.getNumberOfElements();
}

const double*
CAngularMomentumMatrix::xvalues() const
{
    return _xMatrix.values();
}

const double*
CAngularMomentumMatrix::yvalues() const
{
    return _yMatrix.values();
}

const double*
CAngularMomentumMatrix::zvalues() const
{
    return _zMatrix.values();
}

const double*
CAngularMomentumMatrix::values(cartesians cart) const noexcept
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
CAngularMomentumMatrix::values(int32_t cart) const noexcept
{
    auto d = static_cast<cartesians>(cart);
    return values(d);
}

double
CAngularMomentumMatrix::getOriginCoordinateX() const
{
    return _xOrigin;
}

double
CAngularMomentumMatrix::getOriginCoordinateY() const
{
    return _yOrigin;
}

double
CAngularMomentumMatrix::getOriginCoordinateZ() const
{
    return _zOrigin;
}

std::array<double, 3>
CAngularMomentumMatrix::getOriginCoordinates() const
{
    return std::array<double, 3>{{_xOrigin, _yOrigin, _zOrigin}};
}

std::string
CAngularMomentumMatrix::repr() const noexcept
{
    std::ostringstream os;

    os << "[CAngularMomentumMatrix (Object):" << this << "]" << std::endl;

    os << "_xMatrix: " << _xMatrix << std::endl;

    os << "_yMatrix: " << _yMatrix << std::endl;

    os << "_zMatrix: " << _zMatrix << std::endl;

    os << "_xOrigin: " << _xOrigin << std::endl;

    os << "_yOrigin: " << _yOrigin << std::endl;

    os << "_zOrigin: " << _zOrigin << std::endl;

    return os.str();
}

std::ostream&
operator<<(std::ostream& output, const CAngularMomentumMatrix& source)
{
    return (output << source.repr());
}
