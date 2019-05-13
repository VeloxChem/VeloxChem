//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "RecursionMap.hpp"

CRecursionMap::CRecursionMap()
{
    
}

CRecursionMap::CRecursionMap(const std::vector<CRecursionTerm>& recursionTerms,
                             const std::vector<int32_t>&        recursionIndexes)

    : _recursionTerms(recursionTerms)

    , _recursionIndexes(recursionIndexes)
{
    
}

CRecursionMap::CRecursionMap(const CRecursionMap& source)

    : _recursionTerms(source._recursionTerms)

    , _recursionIndexes(source._recursionIndexes)
{
    
}

CRecursionMap::CRecursionMap(CRecursionMap&& source) noexcept

    : _recursionTerms(std::move(source._recursionTerms))

    , _recursionIndexes(std::move(source._recursionIndexes))
{
    
}

CRecursionMap::~CRecursionMap()
{
    
}

CRecursionMap&
CRecursionMap::operator=(const CRecursionMap& source)
{
    if (this == &source) return *this;
    
    _recursionTerms = source._recursionTerms;
    
    _recursionIndexes = source._recursionIndexes;
    
    return *this;
}

CRecursionMap&
CRecursionMap::operator=(CRecursionMap&& source) noexcept
{
    if (this == &source) return *this;
    
    _recursionTerms = std::move(source._recursionTerms);
    
    _recursionIndexes = std::move(source._recursionIndexes);
    
    return *this;
}

bool
CRecursionMap::operator==(const CRecursionMap& other) const
{
    if (_recursionTerms.size() != other._recursionTerms.size()) return false;
    
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (_recursionTerms[i] != other._recursionTerms[i]) return false;
    }
    
    if (_recursionIndexes.size() != other._recursionIndexes.size()) return false;
    
    for (size_t i = 0; i < _recursionIndexes.size(); i++)
    {
        if (_recursionIndexes[i] != other._recursionIndexes[i]) return false;
    }
    
    return true;
}

bool
CRecursionMap::operator!=(const CRecursionMap& other) const
{
    return !(*this == other);
}

std::ostream&
operator<<(      std::ostream&  output,
           const CRecursionMap& source)
{
    output << std::endl;
    
    output << "[CRecursionMap (Object):" << &source << "]" << std::endl;
    
    output << "_recursionTerms: " << std::endl;
    
    for (size_t i = 0; i < source._recursionTerms.size(); i++)
    {
        output << "_recursionTerms[" << i << "]: "<< std::endl;
        
        output << source._recursionTerms[i] << std::endl;
    }

    output << "_recursionIndexes: " << std::endl;
    
    for (size_t i = 0; i < source._recursionIndexes.size(); i++)
    {
        output << "_recursionIndexes[" << i << "]: "<< std::endl;
        
        output << source._recursionIndexes[i] << std::endl;
    }
    
    return output;
}
