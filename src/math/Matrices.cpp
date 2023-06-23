#include "Matrices.hpp"

CMatrices::CMatrices()

    : _matrices(std::map<int64_t, CMatrix*>())
{
}

CMatrices::CMatrices(const std::map<int64_t, CMatrix>& matrices)

    : _matrices(std::map<int64_t, CMatrix*>())
{
    for (const auto& mvalue : matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(mvalue.second)});
    }
}

CMatrices::CMatrices(const CMatrix& matrix, const std::vector<int64_t>& keys)

    : _matrices(std::map<int64_t, CMatrix*>())
{
    for (const auto key : keys)
    {
        _matrices.insert({key, new CMatrix(matrix)});
    }
}

CMatrices::CMatrices(const CMatrices& other)

    : _matrices(std::map<int64_t, CMatrix*>())
{
    for (const auto& mvalue : other._matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(*mvalue.second)});
    }
}

CMatrices::~CMatrices()
{
    for (auto& mvalue : _matrices)
    {
        delete mvalue.second;
    }
}

auto
CMatrices::add(const CMatrix& matrix, const int64_t key) -> void
{
    _matrices.insert({key, new CMatrix(matrix)});
}

auto
CMatrices::getKeys() const -> std::vector<int64_t>
{
    if (!_matrices.empty())
    {
        std::vector<int64_t> keys;

        for (const auto& mvalue : _matrices)
        {
            keys.push_back(mvalue.first);
        }

        return keys;
    }
    else
    {
        return std::vector<int64_t>();
    }
}

auto
CMatrices::getMatrix(const int64_t key) -> CMatrix*
{
    for (const auto& mvalue : _matrices)
    {
        if (mvalue.first == key)
        {
            return mvalue.second;
        }
    }

    return nullptr;
}

auto
CMatrices::getMatrix(const int64_t key) const -> const CMatrix*
{
    for (const auto& mvalue : _matrices)
    {
        if (mvalue.first == key)
        {
            return mvalue.second;
        }
    }

    return nullptr;
}
