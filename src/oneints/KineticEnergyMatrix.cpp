//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "KineticEnergyMatrix.hpp"

CKineticEnergyMatrix::CKineticEnergyMatrix()
{
    
}

CKineticEnergyMatrix::CKineticEnergyMatrix(const CDenseMatrix& matrix)

    : _matrix(matrix)
{
    
}

CKineticEnergyMatrix::CKineticEnergyMatrix(const CKineticEnergyMatrix& source)

    : _matrix(source._matrix)
{
    
}

CKineticEnergyMatrix::CKineticEnergyMatrix(CKineticEnergyMatrix&& source) noexcept

    : _matrix(std::move(source._matrix))
{
    
}

CKineticEnergyMatrix::~CKineticEnergyMatrix()
{
    
}

CKineticEnergyMatrix&
CKineticEnergyMatrix::operator=(const CKineticEnergyMatrix& source)
{
    if (this == &source) return *this;
    
    _matrix = source._matrix;
    
    return *this;
}

CKineticEnergyMatrix&
CKineticEnergyMatrix::operator=(CKineticEnergyMatrix&& source) noexcept
{
    if (this == &source) return *this;
    
    _matrix = std::move(source._matrix);
    
    return *this;
}

bool
CKineticEnergyMatrix::operator==(const CKineticEnergyMatrix& other) const
{
    if (_matrix != other._matrix) return false;
    
    return true;
}

bool
CKineticEnergyMatrix::operator!=(const CKineticEnergyMatrix& other) const
{
    return !(*this == other);
}

std::string
CKineticEnergyMatrix::getString() const
{
    return _matrix.getString();
}

int32_t
CKineticEnergyMatrix::getNumberOfRows() const
{
    return _matrix.getNumberOfRows();
}

int32_t
CKineticEnergyMatrix::getNumberOfColumns() const
{
    return _matrix.getNumberOfColumns();
}

int32_t
CKineticEnergyMatrix::getNumberOfElements() const
{
    return _matrix.getNumberOfElements();
}

const double*
CKineticEnergyMatrix::values() const
{
    return _matrix.values();
}

std::ostream&
operator<<(      std::ostream&         output,
           const CKineticEnergyMatrix& source)
{
    output << std::endl;
    
    output << "[CKineticEnergyMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_matrix: " << source._matrix << std::endl;
    
    return output;
}
