//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "OverlapMatrix.hpp"

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

std::string
COverlapMatrix::getString() const
{
    return _matrix.getString(); 
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

std::ostream&
operator<<(      std::ostream&  output,
           const COverlapMatrix& source)
{
    output << std::endl;
    
    output << "[COverlapMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_matrix: " << source._matrix << std::endl;
    
    return output;
}
