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
    std::stringstream ss ("");
    ss << std::fixed;

    int32_t nrows = _matrix.getNumberOfRows();
    int32_t ncols = _matrix.getNumberOfColumns();

    const double* data = _matrix.values();

    ss << "[Dimension " << nrows << " x " << ncols << "]" << std::endl;

    for (int row = 0; row < nrows; row++) {
        for (int col = 0; col < ncols; col++) {
            ss << std::setw(12) << std::setprecision(7) << data[row * ncols + col];
        }
        ss << std::endl;
    }

    return ss.str();
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
