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

    // cache-blocked, parallel mirror of the upper triangle into the lower:
    // tiling keeps each source/destination tile pair resident in cache rather
    // than striding the whole matrix per element. Each tile row `ti` writes a
    // disjoint set of matrix columns, so the parallel writes do not race.
    const std::size_t     dimension = _rows;
    constexpr std::size_t tile      = 64;

    auto *values = _values.data();

    const auto tiles = (dimension + tile - 1) / tile;

#pragma omp parallel for schedule(dynamic)
    for (std::size_t ti = 0; ti < tiles; ti++)
    {
        const auto i_begin = ti * tile;
        const auto i_end   = std::min(i_begin + tile, dimension);

        for (std::size_t tj = ti; tj < tiles; tj++)
        {
            const auto j_begin = tj * tile;
            const auto j_end   = std::min(j_begin + tile, dimension);

            for (std::size_t i = i_begin; i < i_end; i++)
            {
                for (std::size_t j = std::max(j_begin, i + 1); j < j_end; j++)
                {
                    values[j * dimension + i] = sign * values[i * dimension + j];
                }
            }
        }
    }
}

}  // namespace tabula
