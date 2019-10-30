//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DenseMatrix.hpp"
#include "StringFormat.hpp"

#include <cmath>
#include <sstream>
#include <utility>

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

CDenseMatrix::CDenseMatrix(const int32_t nRows,
                           const int32_t nColumns,
                           const numa    numaPolicy)

    : _nRows(nRows)

    , _nColumns(nColumns)

    , _values(CMemBlock<double>(nRows * nColumns, numaPolicy))
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

void
CDenseMatrix::zero()
{
    mathfunc::zero(_values.data(), _nRows * _nColumns);
}

CDenseMatrix
CDenseMatrix::transpose() const
{
    CDenseMatrix tmat(_nColumns, _nRows);

    auto cvals = _values.data();

    auto tvals = tmat.values();

    for (int32_t i = 0; i < _nRows; i++)
    {
        for (int32_t j = 0; j < _nColumns; j++)
        {
            tvals[j * _nRows + i] = cvals[i * _nColumns + j];
        }
    }

    return tmat;
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

CDenseMatrix
CDenseMatrix::slice(const int32_t iRow, const int32_t iColumn, const int32_t nRows, const int32_t nColumns) const
{
    if (((iRow + nRows) <= _nRows) && ((iColumn + nColumns) <= _nColumns))
    {
        CDenseMatrix mat(nRows, nColumns);

        for (int32_t i = 0; i < nRows; i++)
        {
            // set up pointers to data

            auto srcrow = row(iRow + i);

            auto dstrow = mat.row(i);

            // copy dense matrix values

            #pragma omp simd
            for (int32_t j = 0; j < nColumns; j++)
            {
                dstrow[j] = srcrow[iColumn + j];
            }
        }

        return mat;
    }

    return CDenseMatrix();
}

CDenseMatrix
CDenseMatrix::selectByColumn(const std::vector<int32_t>& iColumns) const
{
    auto ncol = static_cast<int32_t>(iColumns.size());

    if ((ncol > 0) && (ncol <= _nColumns))
    {
        CDenseMatrix mat(_nRows, ncol);

        for (int32_t i = 0; i < _nRows; i++)
        {
            auto sdat = row(i);

            auto ddat = mat.row(i);

            for (int32_t j = 0; j < ncol; j++)
            {
                ddat[j] = sdat[iColumns[j]];
            }
        }

        return mat;
    }

    return CDenseMatrix();
}

CDenseMatrix
CDenseMatrix::selectByRow(const std::vector<int32_t>& iRows) const
{
    auto nrow = static_cast<int32_t>(iRows.size());

    if ((nrow > 0) && (nrow <= _nRows))
    {
        CDenseMatrix mat(nrow, _nColumns);

        for (int32_t i = 0; i < nrow; i++)
        {
            auto sdat = row(iRows[i]);

            auto ddat = mat.row(i);

            #pragma omp simd
            for (int32_t j = 0; j < _nColumns; j++)
            {
                ddat[j] = sdat[j];
            }
        }

        return mat;
    }

    return CDenseMatrix();
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

    sst << "[Dimension " << _nRows << " x " << _nColumns << "]\n";

    for (int32_t i = 0; i < _nRows; i++)
    {
        for (int32_t j = 0; j < _nColumns; j++)
        {
            sst << fstr::to_string(vals[i * _nColumns + j], 8, 15, fmt::right);
        }

        sst << "\n";
    }

    return sst.str();
}

void
CDenseMatrix::broadcast(int32_t rank, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_nRows, comm);

        mpi::bcast(_nColumns, comm);

        _values.broadcast(rank, comm);
    }
}

void
CDenseMatrix::reduce_sum(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        _values.reduce_sum(rank, nodes, comm);
    }
}

std::ostream&
operator<<(std::ostream& output, const CDenseMatrix& source)
{
    output << std::endl;

    output << "[CDenseMatrix (Object):" << &source << "]" << std::endl;

    output << "_nRows: " << source._nRows << std::endl;

    output << "_nColumns: " << source._nColumns << std::endl;

    output << "_values: " << source._values << std::endl;

    return output;
}
