#include "GeomThreeMatrices.hpp"

#include <string>

CGeomThreeMatrices::CGeomThreeMatrices()

    : _matrices(std::map<T3GeomKey, CMatrix*>())
{
    
}

CGeomThreeMatrices::CGeomThreeMatrices(const std::map<T3GeomKey, CMatrix>& matrices)

    : _matrices(std::map<T3GeomKey, CMatrix*>())
{
    for (const auto& mvalue : matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(mvalue.second)});
    }
}

CGeomThreeMatrices::CGeomThreeMatrices(const CMatrix& matrix, const std::vector<int64_t>& atoms)

    : _matrices(std::map<T3GeomKey, CMatrix*>())
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
                
                for (size_t k = j; k < 3; k++)
                {
                    const auto kpair = T2Pair({atom, axis[k]});
                 
                    _matrices.insert({{ipair, jpair, kpair}, new CMatrix(matrix)});
                }
            }
        }
    }
}

CGeomThreeMatrices::CGeomThreeMatrices(const CMatrix& matrix, const std::vector<T3GeomKey>& keys)

    : _matrices(std::map<T3GeomKey, CMatrix*>())
{
    for (const auto& key : keys)
    {
        _matrices.insert({key, new CMatrix(matrix)});
    }
}

CGeomThreeMatrices::CGeomThreeMatrices(const CGeomThreeMatrices& other)

    : _matrices(std::map<T3GeomKey, CMatrix*>())
{
    for (const auto& mvalue : other._matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(*mvalue.second)});
    }
}

CGeomThreeMatrices::~CGeomThreeMatrices()
{
    for (auto& mvalue : _matrices)
    {
        delete mvalue.second;
    }
}

auto
CGeomThreeMatrices::add(const CMatrix& matrix, const T3GeomKey& key) -> void
{
    _matrices.insert({key, new CMatrix(matrix)});
}

auto
CGeomThreeMatrices::getKeys() const -> std::vector<T3GeomKey>
{
    if (!_matrices.empty())
    {
        std::vector<T3GeomKey> keys;

        for (const auto& mvalue : _matrices)
        {
            keys.push_back(mvalue.first);
        }

        return keys;
    }
    else
    {
        return std::vector<T3GeomKey>();
    }
}

auto
CGeomThreeMatrices::getMatrix(const T3GeomKey& key) -> CMatrix*
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
CGeomThreeMatrices::getMatrix(const T3GeomKey& key) const -> const CMatrix*
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
