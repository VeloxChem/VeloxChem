#include "FockMatrices.hpp"

CFockMatrices::CFockMatrices()

    : _matrices(std::vector<CFockMatrix*>())
{
}

CFockMatrices::CFockMatrices(const std::vector<CFockMatrix>& matrices)

    : _matrices(std::vector<CFockMatrix*>())
{
    for (const auto& matrix : matrices)
    {
        _matrices.push_back(new CFockMatrix(matrix));
    }
}

CFockMatrices::CFockMatrices(const CFockMatrices& other)

    : _matrices(std::vector<CFockMatrix*>())
{
    for (const auto matrix : other._matrices)
    {
        _matrices.push_back(new CFockMatrix(*matrix));
    }
}

CFockMatrices::~CFockMatrices()
{
    for (auto matrix : _matrices)
    {
        delete matrix;
    }
}

auto
CFockMatrices::add(const CFockMatrix& matrix) -> void
{
    _matrices.push_back(new CFockMatrix(matrix));
}

auto
CFockMatrices::zero() -> void
{
    for (auto& matrix : _matrices)
    {
        matrix->zero();
    }
}

auto
CFockMatrices::getMatrix(const int64_t index) -> CFockMatrix*
{
    return _matrices[index];
}

auto
CFockMatrices::getMatrix(const int64_t index) const -> const CFockMatrix*
{
    return _matrices[index];
}

auto
CFockMatrices::getNumberOfMatrices() const -> int64_t
{
    return static_cast<int64_t>(_matrices.size());
}
