//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DenseMatrix.hpp"
#include "StringFormat.hpp"

#include <utility>
#include <cmath>
#include <sstream>

CDenseMatrix::CDenseMatrix()

    : _nRows(0)

    , _nColumns(0)
{
    
}

CDenseMatrix::CDenseMatrix(const std::vector<double>& values,
                           const int32_t              nRows,
                           const int32_t              nColumns)

    : _nRows(nRows)

    , _nColumns(nColumns)

    , _values(CMemBlock<double>(values))
{
    
}

CDenseMatrix::CDenseMatrix(const int32_t nRows,
                           const int32_t nColumns)

    : _nRows(nRows)

    , _nColumns(nColumns)

    , _values(CMemBlock<double>(nRows * nColumns))
{
    
}

CDenseMatrix::CDenseMatrix(const int32_t nRows)

    : _nRows(nRows)

    , _nColumns(nRows)

    , _values(CMemBlock<double>(nRows * nRows))
{
    
}

CDenseMatrix::CDenseMatrix(const CDenseMatrix& source)

    : _nRows(source._nRows)

    , _nColumns(source._nColumns)

    , _values(source._values)
{
    
}

CDenseMatrix::CDenseMatrix(CDenseMatrix&& source) noexcept

    : _nRows(std::move(source._nRows))

    , _nColumns(std::move(source._nColumns))

    , _values(std::move(source._values))
{
    
}

CDenseMatrix::~CDenseMatrix()
{
    
}

CDenseMatrix&
CDenseMatrix::operator=(const CDenseMatrix& source)
{
    if (this == &source) return *this;
    
    _nRows = source._nRows;
    
    _nColumns = source._nColumns;
    
    _values = source._values;
    
    return *this;
}

CDenseMatrix&
CDenseMatrix::operator=(CDenseMatrix&& source) noexcept
{
    if (this == &source) return *this;
    
    _nRows = std::move(source._nRows);
    
    _nColumns = std::move(source._nColumns);
    
    _values = std::move(source._values);
    
    return *this;
}

bool
CDenseMatrix::operator==(const CDenseMatrix& other) const
{
    if (_nRows != other._nRows) return false;
    
    if (_nColumns != other._nColumns) return false;
    
    if (_values != other._values) return false;
    
    return true;
}

bool
CDenseMatrix::operator!=(const CDenseMatrix& other) const
{
    return !(*this == other);
}

int32_t
CDenseMatrix::getNumberOfRows() const
{
    return _nRows;
}

int32_t
CDenseMatrix::getNumberOfColumns() const
{
    return _nColumns;
}

int32_t
CDenseMatrix::getNumberOfElements() const
{
    return _nRows * _nColumns;
}

const double*
CDenseMatrix::values() const
{
    return _values.data();
}

double*
CDenseMatrix::values()
{
    return _values.data();
}

const double*
CDenseMatrix::row(const int32_t iRow) const
{
    if (iRow < _nRows)
    {
        return _values.data(iRow * _nColumns);
    }
    
    return nullptr;
}

double*
CDenseMatrix::row(const int32_t iRow)
{
    if (iRow < _nRows)
    {
        return _values.data(iRow * _nColumns);
    }
    
    return nullptr;
}

std::string
CDenseMatrix::getString() const
{
    std::stringstream sst("");
    
    auto vals = _values.data();
    
    sst << "[Dimension " << _nRows << " x " << _nColumns << "]" << std::endl;
    
    for (int32_t i = 0; i < _nRows; i++)
    {
        for (int32_t j = 0; j < _nColumns; j++)
        {
            sst << fstr::to_string(vals[i * _nColumns + j], 8, 15, fmt::right);
        }
        
        sst << std::endl;
    }
    
    return sst.str();
}

void
CDenseMatrix::zero()
{
    mathfunc::zero(_values.data(), _nRows * _nColumns); 
}

void
CDenseMatrix::symmetrize()
{
    if (_nRows == _nColumns)
    {
        auto fmat = _values.data();
        
        for (int32_t i = 0; i < _nRows; i++)
        {
            for (int32_t j = i; j < _nRows; j++)
            {
                auto ijoff = i * _nColumns + j;
            
                auto jioff = j * _nColumns + i;
            
                auto fval = fmat[ijoff] + fmat[jioff];
                
                fmat[ijoff] = fval;
                
                fmat[jioff] = fval; 
            }
        }
    }
}

std::ostream&
operator<<(      std::ostream&  output,
           const CDenseMatrix& source)
{
    output << std::endl;
    
    output << "[CDenseMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_nRows: " << source._nRows << std::endl;
    
    output << "_nColumns: " << source._nColumns << std::endl;
    
    output << "_values: " << source._values <<  std::endl;
    
    return output;
}

