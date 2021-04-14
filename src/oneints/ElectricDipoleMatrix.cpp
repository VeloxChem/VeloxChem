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

#include "ElectricDipoleMatrix.hpp"

#include <array>
#include <cmath>
#include <sstream>

#include "CartesianComponents.hpp"

CElectricDipoleMatrix::CElectricDipoleMatrix()

    : _xOrigin(0.0)

    , _yOrigin(0.0)

    , _zOrigin(0.0)
{
}

CElectricDipoleMatrix::CElectricDipoleMatrix(const CDenseMatrix& xMatrix,
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

CElectricDipoleMatrix::CElectricDipoleMatrix(const std::array<CDenseMatrix, 3>& matrices,
                                             const std::array<double, 3>& origin)

    : CElectricDipoleMatrix(matrices[0], matrices[1], matrices[2], origin[0], origin[1], origin[2])
{
}

CElectricDipoleMatrix::CElectricDipoleMatrix(const CElectricDipoleMatrix& source)

    : _xOrigin(source._xOrigin)

    , _yOrigin(source._yOrigin)

    , _zOrigin(source._zOrigin)

    , _xMatrix(source._xMatrix)

    , _yMatrix(source._yMatrix)

    , _zMatrix(source._zMatrix)
{
}

CElectricDipoleMatrix::CElectricDipoleMatrix(CElectricDipoleMatrix&& source) noexcept

    : _xOrigin(std::move(source._xOrigin))

    , _yOrigin(std::move(source._yOrigin))

    , _zOrigin(std::move(source._zOrigin))

    , _xMatrix(std::move(source._xMatrix))

    , _yMatrix(std::move(source._yMatrix))

    , _zMatrix(std::move(source._zMatrix))
{
}

CElectricDipoleMatrix::~CElectricDipoleMatrix()
{
}

CElectricDipoleMatrix&
CElectricDipoleMatrix::operator=(const CElectricDipoleMatrix& source)
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

CElectricDipoleMatrix&
CElectricDipoleMatrix::operator=(CElectricDipoleMatrix&& source) noexcept
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
CElectricDipoleMatrix::operator==(const CElectricDipoleMatrix& other) const
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
CElectricDipoleMatrix::operator!=(const CElectricDipoleMatrix& other) const
{
    return !(*this == other);
}

void
CElectricDipoleMatrix::setOriginCoordinates(const double xOrigin, const double yOrigin, const double zOrigin)
{
    _xOrigin = xOrigin;

    _yOrigin = yOrigin;

    _zOrigin = zOrigin;
}

void
CElectricDipoleMatrix::setOriginCoordinates(const std::array<double, 3>& origin)
{
    _xOrigin = origin[0];

    _yOrigin = origin[1];

    _zOrigin = origin[2];
}

std::string
CElectricDipoleMatrix::getString() const
{
    return getStringForComponentX() + getStringForComponentY() + getStringForComponentZ();
}

std::string
CElectricDipoleMatrix::getStringForComponentX() const
{
    return _xMatrix.getString();
}

std::string
CElectricDipoleMatrix::getStringForComponentY() const
{
    return _yMatrix.getString();
}

std::string
CElectricDipoleMatrix::getStringForComponentZ() const
{
    return _zMatrix.getString();
}

int32_t
CElectricDipoleMatrix::getNumberOfRows() const
{
    return _xMatrix.getNumberOfRows();
}

int32_t
CElectricDipoleMatrix::getNumberOfColumns() const
{
    return _xMatrix.getNumberOfColumns();
}

int32_t
CElectricDipoleMatrix::getNumberOfElements() const
{
    return _xMatrix.getNumberOfElements();
}

const double*
CElectricDipoleMatrix::xvalues() const
{
    return _xMatrix.values();
}

const double*
CElectricDipoleMatrix::yvalues() const
{
    return _yMatrix.values();
}

const double*
CElectricDipoleMatrix::zvalues() const
{
    return _zMatrix.values();
}

const double*
CElectricDipoleMatrix::values(int32_t cart) const
{
    auto d = static_cast<cartesians>(cart);
    return values(d);
}

const double*
CElectricDipoleMatrix::values(cartesians cart) const
{
    const double* retval = nullptr;
    switch (cart) {
        case cartesians::X:
            retval =  _xMatrix.values();
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

double
CElectricDipoleMatrix::getOriginCoordinateX() const
{
    return _xOrigin;
}

double
CElectricDipoleMatrix::getOriginCoordinateY() const
{
    return _yOrigin;
}

double
CElectricDipoleMatrix::getOriginCoordinateZ() const
{
    return _zOrigin;
}

std::array<double, 3>
CElectricDipoleMatrix::getOriginCoordinates() const
{
    return std::array<double, 3>{{_xOrigin, _yOrigin, _zOrigin}};
}

std::string
CElectricDipoleMatrix::repr() const
{
    std::ostringstream os;

    os << std::endl;

    os << "[CElectricDipoleMatrix (Object):" << this << "]" << std::endl;

    os << "_xMatrix: " << _xMatrix << std::endl;

    os << "_yMatrix: " << _yMatrix << std::endl;

    os << "_zMatrix: " << _zMatrix << std::endl;

    os << "_xOrigin: " << _xOrigin << std::endl;

    os << "_yOrigin: " << _yOrigin << std::endl;

    os << "_zOrigin: " << _zOrigin << std::endl;

    return os.str();
}

std::ostream&
operator<<(std::ostream& output, const CElectricDipoleMatrix& source)
{
    return (output << source.repr());
}
