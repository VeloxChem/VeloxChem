//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "NuclearPotentialMatrix.hpp"

CNuclearPotentialMatrix::CNuclearPotentialMatrix()
{
    
}

CNuclearPotentialMatrix::CNuclearPotentialMatrix(const CSparseMatrix& matrix)

    : _matrix(matrix)
{
    
}

CNuclearPotentialMatrix::CNuclearPotentialMatrix(const CNuclearPotentialMatrix& source)

    : _matrix(source._matrix)
{
    
}

CNuclearPotentialMatrix::CNuclearPotentialMatrix(CNuclearPotentialMatrix&& source) noexcept

    : _matrix(std::move(source._matrix))
{
    
}

CNuclearPotentialMatrix::~CNuclearPotentialMatrix()
{
    
}

CNuclearPotentialMatrix&
CNuclearPotentialMatrix::operator=(const CNuclearPotentialMatrix& source)
{
    if (this == &source) return *this;
    
    _matrix = source._matrix;
    
    return *this;
}

CNuclearPotentialMatrix&
CNuclearPotentialMatrix::operator=(CNuclearPotentialMatrix&& source) noexcept
{
    if (this == &source) return *this;
    
    _matrix = std::move(source._matrix);
    
    return *this;
}

bool
CNuclearPotentialMatrix::operator==(const CNuclearPotentialMatrix& other) const
{
    if (_matrix != other._matrix) return false;
    
    return true;
}

bool
CNuclearPotentialMatrix::operator!=(const CNuclearPotentialMatrix& other) const
{
    return !(*this == other);
}

std::ostream&
operator<<(      std::ostream&         output,
           const CNuclearPotentialMatrix& source)
{
    output << std::endl;
    
    output << "[CNuclearPotentialMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_matrix: " << source._matrix << std::endl;
    
    return output;
}
