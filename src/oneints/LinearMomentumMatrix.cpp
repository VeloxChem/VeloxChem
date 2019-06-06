//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "LinearMomentumMatrix.hpp"

#include <cmath>

CLinearMomentumMatrix::CLinearMomentumMatrix()
{
    
}

CLinearMomentumMatrix::CLinearMomentumMatrix(const CDenseMatrix& xMatrix,
                                             const CDenseMatrix& yMatrix,
                                             const CDenseMatrix& zMatrix)

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
operator<<(      std::ostream&  output,
           const CLinearMomentumMatrix& source)
{
    output << std::endl;
    
    output << "[CLinearMomentumMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_xMatrix: " << source._xMatrix << std::endl;
    
    output << "_yMatrix: " << source._yMatrix << std::endl;
    
    output << "_zMatrix: " << source._zMatrix << std::endl;
    
    return output;
}
