#include "Matrices.hpp"

#include "StringFormat.hpp"

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
CMatrices::add(const CMatrix& matrix, const std::string& label) -> void
{
    _matrices.insert({_to_key(label), new CMatrix(matrix)});
}

auto
CMatrices::zero() -> void
{
    for (auto& mvalue : _matrices)
    {
        mvalue.second->zero();
    }
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

auto
CMatrices::getMatrix(const std::string& label) -> CMatrix*
{
    const auto key = _to_key(label);
    
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
CMatrices::getMatrix(const std::string& label) const -> const CMatrix*
{
    const auto key = _to_key(label);
    
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
CMatrices::_to_key(const std::string& label) const -> int64_t
{
    const auto index = fstr::upcase(label);
    
    if (index.size() == 1)
    {
        if (index == "X") return 0;
        
        if (index == "Y") return 1;
        
        if (index == "Z") return 2;
    }
    
    if (index.size() == 2)
    {
        if (index == "XX") return 0;
        
        if (index == "XY") return 1;
        
        if (index == "XZ") return 2;
        
        if (index == "YY") return 3;
        
        if (index == "YZ") return 4;
        
        if (index == "ZZ") return 5;
    }
    
    if (index.size() == 3)
    {
        if (index == "XXX") return 0;
        
        if (index == "XXY") return 1;
        
        if (index == "XXZ") return 2;
        
        if (index == "XYY") return 3;
        
        if (index == "XYZ") return 4;
        
        if (index == "XZZ") return 5;
        
        if (index == "YYY") return 6;
        
        if (index == "YYZ") return 7;
        
        if (index == "YZZ") return 8;
        
        if (index == "ZZZ") return 9;
    }
    
    return -1;
}
