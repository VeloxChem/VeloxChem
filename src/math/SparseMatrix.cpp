//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "SparseMatrix.hpp"

#include <cmath>
#include <utility>

#include "MathFunc.hpp"

CSparseMatrix::CSparseMatrix()

    : _nRows(0)

    , _nColumns(0)

    , _nMaxElements(0)

    , _nElements(0)

    , _threshold(1.0e-15)
{
}

CSparseMatrix::CSparseMatrix(const std::vector<double>&  values,
                             const std::vector<int32_t>& rows,
                             const std::vector<int32_t>& columns,
                             const int32_t               nRows,
                             const int32_t               nColumns,
                             const double                threshold)

    : _nRows(nRows)

    , _nColumns(nColumns)

    , _values(CMemBlock<double>(values))

    , _rows(CMemBlock<int32_t>(rows))

    , _columns(CMemBlock<int32_t>(columns))

    , _nMaxElements(static_cast<int32_t>(values.size()))

    , _nElements(static_cast<int32_t>(values.size()))

    , _rowPositions(CMemBlock<int32_t>(nRows))

    , _rowSizes(CMemBlock<int32_t>(nRows))

    , _threshold(threshold)
{
    _setAccessPattern();
}

CSparseMatrix::CSparseMatrix(const int32_t nRows, const int32_t nColumns, const double threshold)

    : _nRows(nRows)

    , _nColumns(nColumns)

    , _nElements(0)

    , _rowPositions(CMemBlock<int32_t>(nRows))

    , _rowSizes(CMemBlock<int32_t>(nRows))

    , _threshold(threshold)
{
    _nMaxElements = _setMaxNumberOfElements();

    _values = CMemBlock<double>(_nMaxElements);

    _rows = CMemBlock<int32_t>(_nMaxElements);

    _columns = CMemBlock<int32_t>(_nMaxElements);

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

    _threshold = source._threshold;

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

    // NOTE: max size of buffer is not uniquely defined and depends on
    // constructor used to initialize sparse matrix.

    if (_nElements != other._nElements) return false;

    // compare relevant elements

    for (int32_t i = 0; i < _nElements; i++)
    {
        if (std::fabs(_values.at(i) - other._values.at(i)) > 1.0e-13) return false;

        if (_rows.at(i) != other._rows.at(i)) return false;

        if (_columns.at(i) != other._columns.at(i)) return false;
    }

    if (_rowPositions != other._rowPositions) return false;

    if (_rowSizes != other._rowSizes) return false;

    if (std::fabs(_threshold - other._threshold) > 1.0e-13) return false;

    return true;
}

bool
CSparseMatrix::operator!=(const CSparseMatrix& other) const
{
    return !(*this == other);
}

void
CSparseMatrix::append(const CMemBlock<double>& rowValues, const CMemBlock<int32_t>& rowColumns, const int32_t nElementsInRow, const int32_t iRow)
{
    if (_isLastRow(iRow))
    {
        CMemBlock<int32_t> idxrow(nElementsInRow);

        mathfunc::set_to(idxrow.data(), iRow, nElementsInRow);

        auto ndim = _nElements + nElementsInRow;

        if (ndim > _nMaxElements)
        {
            _nMaxElements += _getAdditionalRows(iRow) * _nColumns;

            // allocate new buffers

            CMemBlock<double> mvalues(_nMaxElements);

            CMemBlock<int32_t> mrows(_nMaxElements);

            CMemBlock<int32_t> mcolumns(_nMaxElements);

            // copy data to new buffers

            mathfunc::copy(mvalues.data(), 0, _values.data(), 0, _nElements);

            mathfunc::copy(mvalues.data(), _nElements, rowValues.data(), 0, nElementsInRow);

            mathfunc::copy(mrows.data(), 0, _rows.data(), 0, _nElements);

            mathfunc::copy(mrows.data(), _nElements, idxrow.data(), 0, nElementsInRow);

            mathfunc::copy(mcolumns.data(), 0, _columns.data(), 0, _nElements);

            mathfunc::copy(mcolumns.data(), _nElements, rowColumns.data(), 0, nElementsInRow);

            // assign new buffers

            _values = mvalues;

            _rows = mrows;

            _columns = mcolumns;
        }
        else
        {
            // copy data to current buffers

            mathfunc::copy(_values.data(), _nElements, rowValues.data(), 0, nElementsInRow);

            mathfunc::copy(_rows.data(), _nElements, idxrow.data(), 0, nElementsInRow);

            mathfunc::copy(_columns.data(), _nElements, rowColumns.data(), 0, nElementsInRow);
        }

        _nElements = ndim;

        _setAccessPattern();
    }
}

void
CSparseMatrix::append(const CMemBlock<double>& rowValues, const CMemBlock<int32_t>& rowColumns, const int32_t iRow)
{
    append(rowValues, rowColumns, rowValues.size(), iRow);
}

void
CSparseMatrix::optimize_storage()
{
    _values.shrink(_nElements);

    _rows.shrink(_nElements);

    _columns.shrink(_nElements);

    _nMaxElements = _nElements;
}

bool
CSparseMatrix::isOptimizedStorage() const
{
    return (_nElements == _nMaxElements);
}

int32_t
CSparseMatrix::getNumberOfRows() const
{
    return _nRows;
}

int32_t
CSparseMatrix::getNumberOfColumns() const
{
    return _nColumns;
}

int32_t
CSparseMatrix::getNumberOfElements() const
{
    return _nElements;
}

int32_t
CSparseMatrix::getNumberOfElements(const int32_t iRow) const
{
    if (iRow < _nRows) return _rowSizes.at(iRow);

    return 0;
}

double
CSparseMatrix::getThreshold() const
{
    return _threshold;
}

double
CSparseMatrix::getSparsity() const
{
    if ((_nRows > 0) && (_nColumns > 0))
    {
        return static_cast<double>(_nElements) / static_cast<double>(_nRows * _nColumns);
    }

    return 0.0;
}

const double*
CSparseMatrix::row(const int32_t iRow) const
{
    if (iRow < _nRows)
    {
        if (_rowSizes.at(iRow) > 0)
        {
            return _values.data(_rowPositions.at(iRow));
        }

        return nullptr;
    }

    return nullptr;
}

double*
CSparseMatrix::row(const int32_t iRow)
{
    if (iRow < _nRows)
    {
        if (_rowSizes.at(iRow) > 0)
        {
            return _values.data(_rowPositions.at(iRow));
        }

        return nullptr;
    }

    return nullptr;
}

const int32_t*
CSparseMatrix::indexes(const int32_t iRow) const
{
    if (iRow < _nRows)
    {
        if (_rowSizes.at(iRow) > 0)
        {
            return _columns.data(_rowPositions.at(iRow));
        }

        return nullptr;
    }

    return nullptr;
}

int32_t*
CSparseMatrix::indexes(const int32_t iRow)
{
    if (iRow < _nRows)
    {
        if (_rowSizes.at(iRow) > 0)
        {
            return _columns.data(_rowPositions.at(iRow));
        }

        return nullptr;
    }

    return nullptr;
}

const double*
CSparseMatrix::values() const
{
    return _values.data();
}

double*
CSparseMatrix::values()
{
    return _values.data();
}

const int32_t*
CSparseMatrix::rows() const
{
    return _rows.data();
}

const int32_t*
CSparseMatrix::columns() const
{
    return _columns.data();
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
            // set up access pattern for row

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

        // set up access pattern for last non-empty row

        _rowPositions.at(idrow) = npos;

        _rowSizes.at(idrow) = nelem;
    }
}

bool
CSparseMatrix::_isLastRow(const int32_t iRow) const
{
    for (int32_t i = iRow; i < _nRows; i++)
    {
        if (_rowPositions.at(i) != -1) return false;
    }

    return true;
}

int32_t
CSparseMatrix::_getAdditionalRows(const int32_t iRow) const
{
    int32_t ndim = _nRows / 6;

    if ((iRow + ndim) > _nRows) return (_nRows - iRow);

    return ndim;
}

int32_t
CSparseMatrix::_setMaxNumberOfElements() const
{
    if ((_nRows > 20000) && (_nColumns > 20000))
    {
        return _nRows * _nColumns / 3;
    }

    if ((_nRows > 10000) && (_nColumns > 10000))
    {
        return _nRows * _nColumns / 2;
    }

    if ((_nRows > 5000) && (_nColumns > 5000))
    {
        return 8 * _nRows * _nColumns / 10;
    }

    return _nRows * _nColumns;
}

std::ostream&
operator<<(std::ostream& output, const CSparseMatrix& source)
{
    output << std::endl;

    output << "[CSparseMatrix (Object):" << &source << "]" << std::endl;

    output << "_nRows: " << source._nRows << std::endl;

    output << "_nColumns: " << source._nColumns << std::endl;

    output << "_values: " << source._values << std::endl;

    output << "_rows: " << source._rows << std::endl;

    output << "_columns: " << source._columns << std::endl;

    output << "_nMaxElements: " << source._nMaxElements << std::endl;

    output << "_nElements: " << source._nElements << std::endl;

    output << "_rowPositions: " << source._rowPositions << std::endl;

    output << "_rowSizes: " << source._rowSizes << std::endl;

    output << "_threshold: " << source._threshold << std::endl;

    return output;
}
