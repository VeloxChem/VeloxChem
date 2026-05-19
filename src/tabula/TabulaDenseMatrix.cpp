//
//  Tabula — custom-recursion molecular-integral machinery.
//  Matrix storage, independent of VeloxChem's own CMatrix.
//

#include "TabulaDenseMatrix.hpp"

#include <algorithm>
#include <string>

#include "ErrorHandler.hpp"

namespace tabula {  // tabula namespace

DenseMatrix::DenseMatrix()

    : _rows{0}

    , _columns{0}

    , _symmetry{Symmetry::general}

    , _values{}
{
}

DenseMatrix::DenseMatrix(const std::size_t rows, const std::size_t columns, const Symmetry symmetry)

    : _rows{rows}

    , _columns{columns}

    , _symmetry{symmetry}

    , _values(rows * columns, 0.0)
{
    if (symmetry != Symmetry::general)
    {
        errors::assertMsgCritical(rows == columns, std::string("DenseMatrix: a symmetric or antisymmetric matrix must be square"));
    }
}

auto
DenseMatrix::operator()(const std::size_t row, const std::size_t column) -> double&
{
    return _values[row * _columns + column];
}

auto
DenseMatrix::operator()(const std::size_t row, const std::size_t column) const -> double
{
    return _values[row * _columns + column];
}

auto
DenseMatrix::rows() const -> std::size_t
{
    return _rows;
}

auto
DenseMatrix::columns() const -> std::size_t
{
    return _columns;
}

auto
DenseMatrix::symmetry() const -> Symmetry
{
    return _symmetry;
}

auto
DenseMatrix::values() -> double*
{
    return _values.data();
}

auto
DenseMatrix::values() const -> const double*
{
    return _values.data();
}

auto
DenseMatrix::zero() -> void
{
    std::fill(_values.begin(), _values.end(), 0.0);
}

auto
DenseMatrix::scale(const double factor) -> void
{
    for (auto& value : _values)
    {
        value *= factor;
    }
}

auto
DenseMatrix::symmetrize() -> void
{
    if (_symmetry == Symmetry::general) return;

    errors::assertMsgCritical(_rows == _columns, std::string("DenseMatrix: symmetrize requires a square matrix"));

    const auto sign = (_symmetry == Symmetry::antisymmetric) ? -1.0 : 1.0;

    for (std::size_t i = 0; i < _rows; i++)
    {
        if (_symmetry == Symmetry::antisymmetric)
        {
            _values[i * _columns + i] = 0.0;
        }

        for (std::size_t j = i + 1; j < _columns; j++)
        {
            _values[j * _columns + i] = sign * _values[i * _columns + j];
        }
    }
}

}  // namespace tabula
