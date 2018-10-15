//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapMatrix.hpp"

#include "DenseDiagonalizer.hpp"

COverlapMatrix::COverlapMatrix()
{
    
}

COverlapMatrix::COverlapMatrix(const CDenseMatrix& matrix)

    : _matrix(matrix)
{
    
}

COverlapMatrix::COverlapMatrix(const COverlapMatrix& source)

    : _matrix(source._matrix)
{
    
}

COverlapMatrix::COverlapMatrix(COverlapMatrix&& source) noexcept

    : _matrix(std::move(source._matrix))
{
    
}

COverlapMatrix::~COverlapMatrix()
{
    
}

COverlapMatrix&
COverlapMatrix::operator=(const COverlapMatrix& source)
{
    if (this == &source) return *this;
    
    _matrix = source._matrix;
    
    return *this;
}

COverlapMatrix&
COverlapMatrix::operator=(COverlapMatrix&& source) noexcept
{
    if (this == &source) return *this;
    
    _matrix = std::move(source._matrix);
  
    return *this;
}

bool
COverlapMatrix::operator==(const COverlapMatrix& other) const
{
    if (_matrix != other._matrix) return false;
    
    return true;
}

bool
COverlapMatrix::operator!=(const COverlapMatrix& other) const
{
    return !(*this == other);
}

std::string
COverlapMatrix::getString() const
{
    return _matrix.getString();
}

int32_t
COverlapMatrix::getNumberOfRows() const
{
    return _matrix.getNumberOfRows();
}

int32_t
COverlapMatrix::getNumberOfColumns() const
{
    return _matrix.getNumberOfColumns();
}

int32_t
COverlapMatrix::getNumberOfElements() const
{
    return _matrix.getNumberOfElements();
}

const double*
COverlapMatrix::values() const
{
    return _matrix.values();
}

CDenseMatrix
COverlapMatrix::getOrthogonalizationMatrix(const double threshold) const
{
    CDenseDiagonalizer diagdrv;
    
    diagdrv.diagonalize(_matrix);
    
    return diagdrv.getInvertedSqrtMatrix(threshold);
}

std::ostream&
operator<<(      std::ostream&  output,
           const COverlapMatrix& source)
{
    output << std::endl;
    
    output << "[COverlapMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_matrix: " << source._matrix << std::endl;
    
    return output;
}
