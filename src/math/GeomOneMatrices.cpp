#include "GeomOneMatrices.hpp"

#include <string>

CGeomOneMatrices::CGeomOneMatrices()

    : _matrices(std::map<TGeomPair, CMatrix*>())
{
}

CGeomOneMatrices::CGeomOneMatrices(const std::map<TGeomPair, CMatrix>& matrices)

    : _matrices(std::map<TGeomPair, CMatrix*>())
{
    for (const auto& mvalue : matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(mvalue.second)});
    }
}

CGeomOneMatrices::CGeomOneMatrices(const CMatrix& matrix, const std::vector<int64_t>& atoms)

    : _matrices(std::map<TGeomPair, CMatrix*>())
{
    const auto axis = std::string("xyz");

    for (const auto atom : atoms)
    {
        for (size_t i = 0; i < 3; i++)
        {
            _matrices.insert({{atom, axis[i]}, new CMatrix(matrix)});
        }
    }
}

CGeomOneMatrices::CGeomOneMatrices(const CMatrix& matrix, const std::vector<TGeomPair>& keys)

    : _matrices(std::map<TGeomPair, CMatrix*>())
{
    for (const auto& key : keys)
    {
        _matrices.insert({key, new CMatrix(matrix)});
    }
}

CGeomOneMatrices::CGeomOneMatrices(const CGeomOneMatrices& other)

    : _matrices(std::map<TGeomPair, CMatrix*>())
{
    for (const auto& mvalue : other._matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(*mvalue.second)});
    }
}

CGeomOneMatrices::~CGeomOneMatrices()
{
    for (auto& mvalue : _matrices)
    {
        delete mvalue.second;
    }
}

auto
CGeomOneMatrices::add(const CMatrix& matrix, const TGeomPair& key) -> void
{
    _matrices.insert({key, new CMatrix(matrix)});
}

auto
CGeomOneMatrices::getKeys() const -> std::vector<TGeomPair>
{
    if (!_matrices.empty())
    {
        std::vector<TGeomPair> keys;

        for (const auto& mvalue : _matrices)
        {
            keys.push_back(mvalue.first);
        }

        return keys;
    }
    else
    {
        return std::vector<TGeomPair>();
    }
}

auto
CGeomOneMatrices::getMatrix(const TGeomPair& key) -> CMatrix*
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
CGeomOneMatrices::getMatrix(const TGeomPair& key) const -> const CMatrix*
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
