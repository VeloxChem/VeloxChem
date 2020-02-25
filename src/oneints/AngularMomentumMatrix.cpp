//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "AngularMomentumMatrix.hpp"

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

std::ostream&
operator<<(std::ostream& output, const CAngularMomentumMatrix& source)
{
    output << std::endl;

    output << "[CAngularMomentumMatrix (Object):" << &source << "]" << std::endl;

    output << "_xMatrix: " << source._xMatrix << std::endl;

    output << "_yMatrix: " << source._yMatrix << std::endl;

    output << "_zMatrix: " << source._zMatrix << std::endl;

    output << "_xOrigin: " << source._xOrigin << std::endl;

    output << "_yOrigin: " << source._yOrigin << std::endl;

    output << "_zOrigin: " << source._zOrigin << std::endl;

    return output;
}
