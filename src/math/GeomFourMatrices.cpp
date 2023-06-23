#include "GeomFourMatrices.hpp"

#include <string>

CGeomFourMatrices::CGeomFourMatrices()

    : _matrices(std::map<T4GeomKey, CMatrix*>())
{
    
}

CGeomFourMatrices::CGeomFourMatrices(const std::map<T4GeomKey, CMatrix>& matrices)

    : _matrices(std::map<T4GeomKey, CMatrix*>())
{
    for (const auto& mvalue : matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(mvalue.second)});
    }
}

CGeomFourMatrices::CGeomFourMatrices(const CMatrix& matrix, const std::vector<int64_t>& atoms)

    : _matrices(std::map<T4GeomKey, CMatrix*>())
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
                    
                    for (size_t l = k; l < 3; l++)
                    {
                        const auto lpair = T2Pair({atom, axis[l]});
                        
                        _matrices.insert({{ipair, jpair, kpair, lpair}, new CMatrix(matrix)});
                    }
                }
            }
        }
    }
}

CGeomFourMatrices::CGeomFourMatrices(const CMatrix& matrix, const std::vector<T4GeomKey>& keys)

    : _matrices(std::map<T4GeomKey, CMatrix*>())
{
    for (const auto& key : keys)
    {
        _matrices.insert({key, new CMatrix(matrix)});
    }
}

CGeomFourMatrices::CGeomFourMatrices(const CGeomFourMatrices& other)

    : _matrices(std::map<T4GeomKey, CMatrix*>())
{
    for (const auto& mvalue : other._matrices)
    {
        _matrices.insert({mvalue.first, new CMatrix(*mvalue.second)});
    }
}

CGeomFourMatrices::~CGeomFourMatrices()
{
    for (auto& mvalue : _matrices)
    {
        delete mvalue.second;
    }
}

auto
CGeomFourMatrices::add(const CMatrix& matrix, const T4GeomKey& key) -> void
{
    _matrices.insert({key, new CMatrix(matrix)});
}

auto
CGeomFourMatrices::getKeys() const -> std::vector<T4GeomKey>
{
    if (!_matrices.empty())
    {
        std::vector<T4GeomKey> keys;

        for (const auto& mvalue : _matrices)
        {
            keys.push_back(mvalue.first);
        }

        return keys;
    }
    else
    {
        return std::vector<T4GeomKey>();
    }
}

auto
CGeomFourMatrices::getMatrix(const T4GeomKey& key) -> CMatrix*
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
CGeomFourMatrices::getMatrix(const T4GeomKey& key) const -> const CMatrix*
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
