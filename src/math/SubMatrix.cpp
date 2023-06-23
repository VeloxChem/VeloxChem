#include "SubMatrix.hpp"

#include <cstring>

#include "AngularMomentum.hpp"

CSubMatrix::CSubMatrix()

    : _values(nullptr)

    , _dimensions({0, 0, 0, 0})
{
}

CSubMatrix::CSubMatrix(const T4Index& dimensions)

    : _values(nullptr)

    , _dimensions(dimensions)
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
        _values = (double*)std::malloc(nelements * sizeof(double));
    }
}

CSubMatrix::CSubMatrix(const std::vector<double>& values, const T4Index& dimensions)

    : _values(nullptr)

    , _dimensions(dimensions)
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
        _values = (double*)std::malloc(nelements * sizeof(double));

        std::memcpy(_values, values.data(), nelements * sizeof(double));
    }
}

CSubMatrix::CSubMatrix(const CSubMatrix& other)

    : _values(nullptr)

    , _dimensions(other._dimensions)

{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
        _values = (double*)std::malloc(nelements * sizeof(double));

        std::memcpy(_values, other._values, nelements * sizeof(double));
    }
}

CSubMatrix::~CSubMatrix()
{
    if (_values != nullptr)
    {
        std::free(_values);
    }
}

auto
CSubMatrix::setOffsets(const int64_t row_offset, const int64_t col_offset) -> void
{
    _dimensions[0] = row_offset;

    _dimensions[1] = col_offset;
}

auto
CSubMatrix::setValues(const std::vector<double>& values) -> void
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements == values.size())
    {
        std::memcpy(_values, values.data(), nelements * sizeof(double));
    }
}

auto
CSubMatrix::getDimensions() const -> T4Index
{
    return _dimensions;
}

auto
CSubMatrix::getValues() const -> std::vector<double>
{
    if (const auto nelements = _dimensions[2] * _dimensions[3]; nelements > 0)
    {
        std::vector<double> values(nelements, 0.0);

        std::memcpy(values.data(), _values, nelements * sizeof(double));

        return values;
    }
    else
    {
        return std::vector<double>();
    }
}

auto
CSubMatrix::getData() -> double*
{
    return _values;
}

auto
CSubMatrix::getData() const -> const double*
{
    return _values;
}
