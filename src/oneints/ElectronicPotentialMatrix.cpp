//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronicPotentialMatrix.hpp"

CElectronicPotentialMatrix::CElectronicPotentialMatrix()
{
    
}

CElectronicPotentialMatrix::CElectronicPotentialMatrix(const CDenseMatrix& matrix)

    : _matrix(matrix)
{
    
}

CElectronicPotentialMatrix::CElectronicPotentialMatrix(const CElectronicPotentialMatrix& source)

    : _matrix(source._matrix)
{
    
}

CElectronicPotentialMatrix::CElectronicPotentialMatrix(CElectronicPotentialMatrix&& source) noexcept

    : _matrix(std::move(source._matrix))
{
    
}

CElectronicPotentialMatrix::~CElectronicPotentialMatrix()
{
    
}

CElectronicPotentialMatrix&
CElectronicPotentialMatrix::operator=(const CElectronicPotentialMatrix& source)
{
    if (this == &source) return *this;
    
    _matrix = source._matrix;
    
    return *this;
}

CElectronicPotentialMatrix&
CElectronicPotentialMatrix::operator=(CElectronicPotentialMatrix&& source) noexcept
{
    if (this == &source) return *this;
    
    _matrix = std::move(source._matrix);
    
    return *this;
}

bool
CElectronicPotentialMatrix::operator==(const CElectronicPotentialMatrix& other) const
{
    if (_matrix != other._matrix) return false;
    
    return true;
}

bool
CElectronicPotentialMatrix::operator!=(const CElectronicPotentialMatrix& other) const
{
    return !(*this == other);
}

std::ostream&
operator<<(      std::ostream&         output,
           const CElectronicPotentialMatrix& source)
{
    output << std::endl;
    
    output << "[CElectronicPotentialMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_matrix: " << source._matrix << std::endl;
    
    return output;
}
