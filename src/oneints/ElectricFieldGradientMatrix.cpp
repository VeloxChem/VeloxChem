//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectricFieldGradientMatrix.hpp"

CElectricFieldGradientMatrix::CElectricFieldGradientMatrix()
{
}

CElectricFieldGradientMatrix::CElectricFieldGradientMatrix(const CDenseMatrix& xxMatrix,
                                                           const CDenseMatrix& xyMatrix,
                                                           const CDenseMatrix& xzMatrix,
                                                           const CDenseMatrix& yyMatrix,
                                                           const CDenseMatrix& yzMatrix,
                                                           const CDenseMatrix& zzMatrix)

    : _xxMatrix(xxMatrix)

    , _xyMatrix(xyMatrix)

    , _xzMatrix(xzMatrix)

    , _yyMatrix(yyMatrix)

    , _yzMatrix(yzMatrix)

    , _zzMatrix(zzMatrix)
{
}

CElectricFieldGradientMatrix::CElectricFieldGradientMatrix(const CElectricFieldGradientMatrix& source)

    : _xxMatrix(source._xxMatrix)

    , _xyMatrix(source._xyMatrix)

    , _xzMatrix(source._xzMatrix)

    , _yyMatrix(source._yyMatrix)

    , _yzMatrix(source._yzMatrix)

    , _zzMatrix(source._zzMatrix)
{
}

CElectricFieldGradientMatrix::CElectricFieldGradientMatrix(CElectricFieldGradientMatrix&& source) noexcept

    : _xxMatrix(std::move(source._xxMatrix))

    , _xyMatrix(std::move(source._xyMatrix))

    , _xzMatrix(std::move(source._xzMatrix))

    , _yyMatrix(std::move(source._yyMatrix))

    , _yzMatrix(std::move(source._yzMatrix))

    , _zzMatrix(std::move(source._zzMatrix))
{
}

CElectricFieldGradientMatrix::~CElectricFieldGradientMatrix()
{
}

CElectricFieldGradientMatrix&
CElectricFieldGradientMatrix::operator=(const CElectricFieldGradientMatrix& source)
{
    if (this == &source) return *this;

    _xxMatrix = source._xxMatrix;

    _xyMatrix = source._xyMatrix;

    _xzMatrix = source._xzMatrix;

    _yyMatrix = source._yyMatrix;

    _yzMatrix = source._yzMatrix;

    _zzMatrix = source._zzMatrix;

    return *this;
}

CElectricFieldGradientMatrix&
CElectricFieldGradientMatrix::operator=(CElectricFieldGradientMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _xxMatrix = std::move(source._xxMatrix);

    _xyMatrix = std::move(source._xyMatrix);

    _xzMatrix = std::move(source._xzMatrix);

    _yyMatrix = std::move(source._yyMatrix);

    _yzMatrix = std::move(source._yzMatrix);

    _zzMatrix = std::move(source._zzMatrix);

    return *this;
}

bool
CElectricFieldGradientMatrix::operator==(const CElectricFieldGradientMatrix& other) const
{
    if (_xxMatrix != other._xxMatrix) return false;

    if (_xyMatrix != other._xyMatrix) return false;

    if (_xzMatrix != other._xzMatrix) return false;

    if (_yyMatrix != other._yyMatrix) return false;

    if (_yzMatrix != other._yzMatrix) return false;

    if (_zzMatrix != other._zzMatrix) return false;

    return true;
}

bool
CElectricFieldGradientMatrix::operator!=(const CElectricFieldGradientMatrix& other) const
{
    return !(*this == other);
}

std::string
CElectricFieldGradientMatrix::getStringForComponentXX() const
{
    return _xxMatrix.getString();
}

std::string
CElectricFieldGradientMatrix::getStringForComponentXY() const
{
    return _xyMatrix.getString();
}

std::string
CElectricFieldGradientMatrix::getStringForComponentXZ() const
{
    return _xzMatrix.getString();
}

std::string
CElectricFieldGradientMatrix::getStringForComponentYY() const
{
    return _yyMatrix.getString();
}

std::string
CElectricFieldGradientMatrix::getStringForComponentYZ() const
{
    return _yzMatrix.getString();
}

std::string
CElectricFieldGradientMatrix::getStringForComponentZZ() const
{
    return _zzMatrix.getString();
}

int32_t
CElectricFieldGradientMatrix::getNumberOfRows() const
{
    return _xxMatrix.getNumberOfRows();
}

int32_t
CElectricFieldGradientMatrix::getNumberOfColumns() const
{
    return _xxMatrix.getNumberOfColumns();
}

int32_t
CElectricFieldGradientMatrix::getNumberOfElements() const
{
    return _xxMatrix.getNumberOfElements();
}

const double*
CElectricFieldGradientMatrix::xxvalues() const
{
    return _xxMatrix.values();
}

const double*
CElectricFieldGradientMatrix::xyvalues() const
{
    return _xyMatrix.values();
}

const double*
CElectricFieldGradientMatrix::xzvalues() const
{
    return _xzMatrix.values();
}

const double*
CElectricFieldGradientMatrix::yyvalues() const
{
    return _yyMatrix.values();
}

const double*
CElectricFieldGradientMatrix::yzvalues() const
{
    return _yzMatrix.values();
}

const double*
CElectricFieldGradientMatrix::zzvalues() const
{
    return _zzMatrix.values();
}

std::ostream&
operator<<(std::ostream& output, const CElectricFieldGradientMatrix& source)
{
    output << std::endl;

    output << "[CElectricFieldGradientMatrix (Object):" << &source << "]" << std::endl;

    output << "_xxMatrix: " << source._xxMatrix << std::endl;

    output << "_xyMatrix: " << source._xyMatrix << std::endl;

    output << "_xzMatrix: " << source._xzMatrix << std::endl;

    output << "_yyMatrix: " << source._yyMatrix << std::endl;

    output << "_yzMatrix: " << source._yzMatrix << std::endl;

    output << "_zzMatrix: " << source._zzMatrix << std::endl;

    return output;
}
