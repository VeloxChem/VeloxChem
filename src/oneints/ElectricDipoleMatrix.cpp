//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectricDipoleMatrix.hpp"

#include <cmath>

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

std::ostream&
operator<<(std::ostream& output, const CElectricDipoleMatrix& source)
{
    output << std::endl;

    output << "[CElectricDipoleMatrix (Object):" << &source << "]" << std::endl;

    output << "_xMatrix: " << source._xMatrix << std::endl;

    output << "_yMatrix: " << source._yMatrix << std::endl;

    output << "_zMatrix: " << source._zMatrix << std::endl;

    output << "_xOrigin: " << source._xOrigin << std::endl;

    output << "_yOrigin: " << source._yOrigin << std::endl;

    output << "_zOrigin: " << source._zOrigin << std::endl;

    return output;
}
