//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SparseMatrix.hpp"

#include <utility>
#include <cmath>

#include "MathFunc.hpp"

CSparseMatrix::CSparseMatrix()

    : _nRows(0)

    , _nColumns(0)

    , _nMaxElements(0)

    , _nElements(0)

    , _threshold(1.0e-13)
{
    
}

CSparseMatrix::CSparseMatrix(const std::vector<double>&  values,
                             const std::vector<int32_t>& rows,
                             const std::vector<int32_t>& columns,
                             const int32_t               nRows,
                             const int32_t               nColumns,
                             const double                threshold)

    : _values(CMemBlock<double>(values))

    , _rows(CMemBlock<int32_t>(rows))

    , _columns(CMemBlock<int32_t>(columns))

    , _nRows(nRows)

    , _nColumns(nColumns)

    , _rowPositions(CMemBlock<int32_t>(_nRows))

    , _rowSizes(CMemBlock<int32_t>(_nRows))

    , _nElements(_values.size())

    , _nMaxElements(_values.size())

    , _threshold(threshold)
{
    _setAccessPattern(); 
}

CSparseMatrix::CSparseMatrix(const CSparseMatrix& source)

    : _nRows(source._nRows)

    , _nColumns(source._nColumns)

    , _values(source._values)

    , _rows(source._rows)

    , _columns(source._columns)

    , _nMaxElements(source._nMaxElements)

    , _nElements(source._nElements)

    , _rowPositions(source._rowPositions)

    , _rowSizes(source._rowSizes)

    , _threshold(source._threshold)
{
    
}

CSparseMatrix::CSparseMatrix(CSparseMatrix&& source) noexcept

    : _nRows(std::move(source._nRows))

    , _nColumns(std::move(source._nColumns))

    , _values(std::move(source._values))

    , _rows(std::move(source._rows))

    , _columns(std::move(source._columns))

    , _nMaxElements(std::move(source._nMaxElements))

    , _nElements(std::move(source._nElements))

    , _rowPositions(std::move(source._rowPositions))

    , _rowSizes(std::move(source._rowSizes))

    , _threshold(std::move(source._threshold))
{
    
}

CSparseMatrix::~CSparseMatrix()
{
    
}

CSparseMatrix&
CSparseMatrix::operator=(const CSparseMatrix& source)
{
    if (this == &source) return *this;
    
    _nRows = source._nRows;
    
    _nColumns = source._nColumns;
    
    _values = source._values;
    
    _rows = source._rows;
    
    _columns = source._columns;
    
    _nMaxElements = source._nMaxElements;
    
    _nElements = source._nElements;
    
    _rowPositions = source._rowPositions;
    
    _rowSizes = source._rowSizes;
    
    _threshold= source._threshold;
    
    return *this;
}

CSparseMatrix&
CSparseMatrix::operator=(CSparseMatrix&& source) noexcept
{
    if (this == &source) return *this;
    
    _nRows = std::move(source._nRows);
    
    _nColumns = std::move(source._nColumns);
    
    _values = std::move(source._values);
    
    _rows = std::move(source._rows);
    
    _columns = std::move(source._columns);
    
    _nMaxElements = std::move(source._nMaxElements);
    
    _nElements = std::move(source._nElements);
    
    _rowPositions = std::move(source._rowPositions);
    
    _rowSizes = std::move(source._rowSizes);
    
    _threshold = std::move(source._threshold);
   
    return *this;
}

bool
CSparseMatrix::operator==(const CSparseMatrix& other) const
{
    if (_nRows != other._nRows) return false;
    
    if (_nColumns != other._nColumns) return false;
    
    if (_values != other._values) return false;
    
    if (_rows != other._rows) return false;
    
    if (_columns != other._columns) return false;
    
    if (_nMaxElements != other._nMaxElements) return false;
    
    if (_nElements != other._nElements) return false;
    
    if (_rowPositions != other._rowPositions) return false;
    
    if (_rowSizes != other._rowSizes) return false;
    
    if (std::fabs(_threshold - other._threshold) > 10e-13) return false;
    
    return true;
}

bool
CSparseMatrix::operator!=(const CSparseMatrix& other) const
{
    return !(*this == other);
}

void
CSparseMatrix::_setAccessPattern()
{
    mathfunc::zero(_rowSizes.data(), _nRows);
    
    mathfunc::set_to(_rowPositions.data(), -1, _nRows);
    
    if (_nElements > 0)
    {
        int32_t nelem = 0;
        
        int32_t npos = 0;
        
        auto idrow = _rows.at(0);
        
        for (int32_t i = 0; i < _nElements; i++)
        {
            // set up row pattern 
            
            if (idrow != _rows.at(i))
            {
                _rowPositions.at(idrow) = npos;
                
                _rowSizes.at(idrow) = nelem;
                
                idrow = _rows.at(i);
                
                npos += nelem;
                
                nelem = 0;
            }
            
            nelem++;
        }
    }
}

std::ostream&
operator<<(      std::ostream&  output,
           const CSparseMatrix& source)
{
    output << std::endl;
    
    output << "[CSparseMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_nRows: " << source._nRows << std::endl;
    
    output << "_nColumns: " << source._nColumns << std::endl;
    
    output << "_values: " << source._values <<  std::endl;
    
    output << "_rows: " << source._rows <<  std::endl;
    
    output << "_columns: " << source._columns <<  std::endl;
    
    output << "_nMaxElements: " << source._nMaxElements  <<  std::endl;
    
    output << "_nElements: " << source._nElements  <<  std::endl;
    
    output << "_rowPositions: " << source._rowPositions  <<  std::endl;
    
    output << "_rowSizes: " << source._rowSizes  <<  std::endl;
    
    output << "_threshold: " << source._threshold  <<  std::endl;
    
    return output;
}

