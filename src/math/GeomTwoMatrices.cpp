#include "GeomTwoMatrices.hpp"

#include <string>

CGeomTwoMatrices::CGeomTwoMatrices()

    : _matrices(std::map<T2GeomKey, CMatrix*>())
{
    
}

CGeomTwoMatrices::CGeomTwoMatrices(const std::map<T2GeomKey, CMatrix>& matrices)

    : _matrices(std::map<T2GeomKey, CMatrix*>())
{
    for (const auto& mvalue : matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(mvalue.second)});
    }
}

CGeomTwoMatrices::CGeomTwoMatrices(const CMatrix& matrix, const std::vector<int64_t>& atoms)

    : _matrices(std::map<T2GeomKey, CMatrix*>())
{
    const auto axis = std::string("xyz");
    
    for (const auto atom : atoms)
    {
        for (size_t i = 0; i < 3; i++)
        {
            const auto ipair = T2Pair({atom, axis[i]});
            
            for (size_t j = i; j < 3; j++)
            {
                const auto jpair = T2Pair({atom, axis[j]});
                
                _matrices.insert({{ipair, jpair}, new CMatrix(matrix)});
            }
        }
    }
}

CGeomTwoMatrices::CGeomTwoMatrices(const CMatrix& matrix, const std::vector<T2GeomKey>& keys)

    : _matrices(std::map<T2GeomKey, CMatrix*>())
{
    for (const auto& key : keys)
    {
        _matrices.insert({key, new CMatrix(matrix)});
    }
}

CGeomTwoMatrices::CGeomTwoMatrices(const CGeomTwoMatrices& other)

    : _matrices(std::map<T2GeomKey, CMatrix*>())
{
    for (const auto& mvalue : other._matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(*mvalue.second)});
    }
}

CGeomTwoMatrices::~CGeomTwoMatrices()
{
    for (auto& mvalue : _matrices)
    {
        delete mvalue.second;
    }
}

auto
CGeomTwoMatrices::add(const CMatrix& matrix, const T2GeomKey& key) -> void
{
    _matrices.insert({key, new CMatrix(matrix)});
}

auto
CGeomTwoMatrices::getKeys() const -> std::vector<T2GeomKey>
{
    if (!_matrices.empty())
    {
        std::vector<T2GeomKey> keys;

        for (const auto& mvalue : _matrices)
        {
            keys.push_back(mvalue.first);
        }

        return keys;
    }
    else
    {
        return std::vector<T2GeomKey>();
    }
}

auto
CGeomTwoMatrices::getMatrix(const T2GeomKey& key) -> CMatrix*
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
CGeomTwoMatrices::getMatrix(const T2GeomKey& key) const -> const CMatrix*
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
