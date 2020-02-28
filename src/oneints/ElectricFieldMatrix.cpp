//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectricFieldMatrix.hpp"

CElectricFieldMatrix::CElectricFieldMatrix()
{
}

CElectricFieldMatrix::CElectricFieldMatrix(const CDenseMatrix& xMatrix, const CDenseMatrix& yMatrix, const CDenseMatrix& zMatrix)

    : _xMatrix(xMatrix)

    , _yMatrix(yMatrix)

    , _zMatrix(zMatrix)
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

std::ostream&
operator<<(std::ostream& output, const CElectricFieldMatrix& source)
{
    output << std::endl;

    output << "[CElectricFieldMatrix (Object):" << &source << "]" << std::endl;

    output << "_xMatrix: " << source._xMatrix << std::endl;

    output << "_yMatrix: " << source._yMatrix << std::endl;

    output << "_zMatrix: " << source._zMatrix << std::endl;

    return output;
}
